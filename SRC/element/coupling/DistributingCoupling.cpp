/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: N. Mora Bowen, P. Palacios, J.A. Abell
// Created: 2026
//
// Description: DistributingCoupling, a distributing coupling (RBE3-type).
// See DistributingCoupling.h for the formulation and references
// (MSC Nastran RBE3; Abaqus *COUPLING, *DISTRIBUTING).

#include <DistributingCoupling.h>
#include <classTags.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <ElementResponse.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

namespace {

// largest |M(i,i)|: stiffness scale of the -host element for -k auto
double maxAbsDiag(const Matrix &M)
{
  double s = 0.0;
  int n = M.noRows() < M.noCols() ? M.noRows() : M.noCols();
  for (int i = 0; i < n; i++)
    if (fabs(M(i, i)) > s) s = fabs(M(i, i));
  return s;
}

bool tokenIs(const char *arg, const char *const *names)
{
  for (int i = 0; names[i] != 0; i++)
    if (strcmp(arg, names[i]) == 0) return true;
  return false;
}

// Cyclic Jacobi eigensolver for a symmetric 3x3 matrix (W.H. Press et al.,
// Numerical Recipes, section 11.1). a is destroyed; d holds the
// eigenvalues and the columns of v the eigenvectors.
void jacobi3(double a[3][3], double d[3], double v[3][3])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      v[i][j] = (i == j) ? 1.0 : 0.0;
  double b[3], z[3];
  for (int i = 0; i < 3; i++) {
    d[i] = a[i][i];
    b[i] = d[i];
    z[i] = 0.0;
  }
#define JROT(M, i, j, k, l) { double g_ = M[i][j]; double h_ = M[k][l]; \
    M[i][j] = g_ - s * (h_ + g_ * tau); M[k][l] = h_ + s * (g_ - h_ * tau); }
  for (int sweep = 0; sweep < 50; sweep++) {
    double sm = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
    if (sm == 0.0) return;
    double tresh = (sweep < 3) ? (0.2 * sm / 9.0) : 0.0;
    for (int p = 0; p < 2; p++) {
      for (int q = p + 1; q < 3; q++) {
        double g = 100.0 * fabs(a[p][q]);
        if (sweep > 3 && (fabs(d[p]) + g) == fabs(d[p]) && (fabs(d[q]) + g) == fabs(d[q])) {
          a[p][q] = 0.0;
        } else if (fabs(a[p][q]) > tresh) {
          double h = d[q] - d[p], t;
          if ((fabs(h) + g) == fabs(h)) {
            t = a[p][q] / h;
          } else {
            double theta = 0.5 * h / a[p][q];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0) t = -t;
          }
          double c = 1.0 / sqrt(1.0 + t * t), s = t * c, tau = s / (1.0 + c);
          double hh = t * a[p][q];
          z[p] -= hh; z[q] += hh; d[p] -= hh; d[q] += hh; a[p][q] = 0.0;
          for (int j = 0; j < p; j++) JROT(a, j, p, j, q)
          for (int j = p + 1; j < q; j++) JROT(a, p, j, j, q)
          for (int j = q + 1; j < 3; j++) JROT(a, p, j, q, j)
          for (int j = 0; j < 3; j++) JROT(v, j, p, j, q)
        }
      }
    }
    for (int i = 0; i < 3; i++) {
      b[i] += z[i];
      d[i] = b[i];
      z[i] = 0.0;
    }
  }
#undef JROT
}

} // namespace

void *OPS_DistributingCoupling(void)
{
  int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) {
    opserr << "WARNING DistributingCoupling: model ndm must be 2 or 3\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient args\n"
           << "Want: element DistributingCoupling tag refNode N i1 ... iN "
           << "<-w w1 ... wN> <-k Kt|auto> <-kAlpha a> <-host eleTag> <-kr Kr> "
           << "<-enforce penalty|al> <-absolute>\n";
    return 0;
  }

  int idata[3];                    // tag, refNode, N
  int n = 3;
  if (OPS_GetIntInput(&n, idata) < 0) {
    opserr << "WARNING DistributingCoupling: invalid tag, refNode or N\n";
    return 0;
  }
  int tag = idata[0], refNode = idata[1], N = idata[2];
  if (N < 1) {
    opserr << "WARNING DistributingCoupling " << tag
           << ": N (number of independent nodes) must be >= 1\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < N) {
    opserr << "WARNING DistributingCoupling " << tag << ": need " << N
           << " independent node tags\n";
    return 0;
  }
  ID indep(N);
  for (int i = 0; i < N; i++) {
    int h;
    n = 1;
    if (OPS_GetIntInput(&n, &h) < 0) {
      opserr << "WARNING DistributingCoupling " << tag << ": invalid independent node "
             << i + 1 << "\n";
      return 0;
    }
    indep(i) = h;
  }

  Vector weights(N);
  for (int i = 0; i < N; i++) weights(i) = 1.0;

  double Kt = 1.0e12;
  bool ktAuto = false;
  double kAlpha = 1.0e3;
  int hostEleTag = -1;
  double Kr = 0.0;
  bool krUser = false;
  int enforce = 0;
  bool initGapCapture = true;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-w") == 0) {
      n = N;
      if (OPS_GetNumRemainingInputArgs() < N || OPS_GetDoubleInput(&n, &weights(0)) < 0) {
        opserr << "WARNING DistributingCoupling " << tag << ": -w wants " << N
               << " weights\n";
        return 0;
      }
      for (int i = 0; i < N; i++)
        if (weights(i) < 0.0) {
          opserr << "WARNING DistributingCoupling " << tag
                 << ": weights must be non-negative\n";
          return 0;
        }
    }
    else if (strcmp(opt, "-k") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING DistributingCoupling " << tag << ": -k wants a value or 'auto'\n";
        return 0;
      }
      char kTok[64];
      OPS_GetStringFromAll(kTok, sizeof(kTok));
      if (strcmp(kTok, "auto") == 0) {
        ktAuto = true;
      } else {
        char *endp = 0;
        double kv = strtod(kTok, &endp);
        if (endp == kTok || *endp != '\0' || kv <= 0.0) {
          opserr << "WARNING DistributingCoupling " << tag
                 << ": -k wants a positive number or 'auto', got '" << kTok << "'\n";
          return 0;
        }
        Kt = kv;
      }
    }
    else if (strcmp(opt, "-kAlpha") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &kAlpha) < 0 || kAlpha <= 0.0) {
        opserr << "WARNING DistributingCoupling " << tag << ": -kAlpha wants a positive value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-host") == 0) {
      n = 1;
      if (OPS_GetIntInput(&n, &hostEleTag) < 0) {
        opserr << "WARNING DistributingCoupling " << tag << ": -host wants an element tag\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-kr") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &Kr) < 0 || Kr < 0.0) {
        opserr << "WARNING DistributingCoupling " << tag << ": -kr wants a non-negative value\n";
        return 0;
      }
      krUser = true;
    }
    else if (strcmp(opt, "-enforce") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING DistributingCoupling " << tag << ": -enforce wants penalty or al\n";
        return 0;
      }
      const char *mode = OPS_GetString();
      if (strcmp(mode, "penalty") == 0)
        enforce = 0;
      else if (strcmp(mode, "al") == 0)
        enforce = 1;
      else {
        opserr << "WARNING DistributingCoupling " << tag << ": unknown -enforce '" << mode
               << "' (want penalty or al)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-absolute") == 0) {
      initGapCapture = false;
    }
    else {
      opserr << "WARNING DistributingCoupling " << tag << ": unknown option '" << opt << "'\n";
      return 0;
    }
  }

  if (ktAuto && hostEleTag < 0) {
    opserr << "WARNING DistributingCoupling " << tag
           << ": -k auto requires a representative -host element\n";
    return 0;
  }

  return new DistributingCoupling(tag, ndm, refNode, indep, weights, Kt, Kr, krUser,
                                  enforce, kAlpha, hostEleTag, ktAuto, initGapCapture);
}

// ===========================================================================
//  construction
// ===========================================================================
DistributingCoupling::DistributingCoupling(int tag, int ndm_, int refNode,
                                           const ID &indepNodes, const Vector &weights_,
                                           double kt, double kr, bool krUser_,
                                           int enforce_, double kAlpha_, int hostEleTag_,
                                           bool ktAuto_, bool initGapCapture_)
  : Element(tag, ELE_TAG_DistributingCoupling),
    ndm(ndm_), nrot((ndm_ == 3) ? 3 : 1), nIndep(indepNodes.Size()),
    connectedNodes(1 + indepNodes.Size()), weights(weights_),
    Kt(kt), Kr(kr), krUser(krUser_), ktAuto(ktAuto_), kAlpha(kAlpha_),
    hostEleTag(hostEleTag_), ktResolved(false), hostMissWarned(false),
    enforce(enforce_), lambda(ndm_), lambda_r((ndm_ == 3) ? 3 : 1),
    valid(false), W(0.0), ell2(0.0), dRef(ndm_), rvec(),
    Icplus((ndm_ == 3) ? 3 : 1, (ndm_ == 3) ? 3 : 1),
    Pproj((ndm_ == 3) ? 3 : 1, (ndm_ == 3) ? 3 : 1), nKept(0),
    nDOF(0), nodeNdf(1 + indepNodes.Size()), dofOffset(1 + indepNodes.Size()),
    B(0), initGapCapture(initGapCapture_), g0Computed(false),
    g0(ndm_), gr0((ndm_ == 3) ? 3 : 1),
    theNodes(0), K(0), P(0), M0(0)
{
  connectedNodes(0) = refNode;
  for (int i = 0; i < nIndep; i++)
    connectedNodes(1 + i) = indepNodes(i);
  theNodes = new Node *[1 + nIndep];
  for (int i = 0; i < 1 + nIndep; i++)
    theNodes[i] = 0;
}

DistributingCoupling::DistributingCoupling()
  : Element(0, ELE_TAG_DistributingCoupling),
    ndm(0), nrot(0), nIndep(0), connectedNodes(), weights(),
    Kt(0.0), Kr(0.0), krUser(false), ktAuto(false), kAlpha(0.0),
    hostEleTag(-1), ktResolved(false), hostMissWarned(false),
    enforce(0), lambda(), lambda_r(),
    valid(false), W(0.0), ell2(0.0), dRef(), rvec(),
    Icplus(), Pproj(), nKept(0),
    nDOF(0), nodeNdf(), dofOffset(),
    B(0), initGapCapture(true), g0Computed(false), g0(), gr0(),
    theNodes(0), K(0), P(0), M0(0)
{
}

DistributingCoupling::~DistributingCoupling()
{
  if (theNodes != 0) delete[] theNodes;
  if (B != 0) delete B;
  if (K != 0) delete K;
  if (P != 0) delete P;
  if (M0 != 0) delete M0;
}

// ===========================================================================
//  domain
// ===========================================================================
int DistributingCoupling::getNumExternalNodes(void) const { return 1 + nIndep; }
const ID &DistributingCoupling::getExternalNodes(void) { return connectedNodes; }
Node **DistributingCoupling::getNodePtrs(void) { return theNodes; }
int DistributingCoupling::getNumDOF(void) { return nDOF; }

void DistributingCoupling::allocate(void)
{
  if (K != 0) delete K;
  K = new Matrix(nDOF, nDOF);
  if (P != 0) delete P;
  P = new Vector(nDOF);
  if (M0 != 0) delete M0;
  M0 = new Matrix(nDOF, nDOF);
}

void DistributingCoupling::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    if (theNodes != 0)
      for (int i = 0; i < 1 + nIndep; i++) theNodes[i] = 0;
    return;
  }

  // the reference node needs translations and rotations, the independent
  // nodes translations
  nodeNdf.resize(1 + nIndep);
  dofOffset.resize(1 + nIndep);
  int pos = 0;
  bool refused = false;
  for (int i = 0; i < 1 + nIndep; i++) {
    theNodes[i] = theDomain->getNode(connectedNodes(i));
    if (theNodes[i] == 0) {
      opserr << "WARNING DistributingCoupling " << this->getTag() << ": node "
             << connectedNodes(i) << " not found\n";
      refused = true;
      nodeNdf(i) = 0;
      dofOffset(i) = pos;
      continue;
    }
    int ndf = theNodes[i]->getNumberDOF();
    int need = (i == 0) ? (ndm + nrot) : ndm;
    if (ndf < need) {
      opserr << "WARNING DistributingCoupling " << this->getTag() << ": "
             << (i == 0 ? "reference" : "independent") << " node " << connectedNodes(i)
             << " has " << ndf << " DOFs, needs at least " << need << "\n";
      refused = true;
    }
    nodeNdf(i) = ndf;
    dofOffset(i) = pos;
    pos += ndf;
  }
  nDOF = pos;

  if (refused)
    valid = false;
  else
    this->resolveGeometry();

  // An ill-posed element is kept inert (zero B, K, P) rather than leaving
  // unallocated matrices for the assembler.
  if (B != 0) delete B;
  B = new Matrix(ndm + nrot, (nDOF > 0 ? nDOF : 1));
  B->Zero();
  this->allocate();
  if (valid) {
    this->buildB();
    this->captureInitialGap();
  }
  this->DomainComponent::setDomain(theDomain);
}

// weighted centroid, r_i, transport lever x_c - x_R, l^2, and the
// pseudo-inverse and range projector of I_c
void DistributingCoupling::resolveGeometry(void)
{
  valid = false;
  W = 0.0;
  for (int i = 0; i < nIndep; i++)
    W += weights(i);
  if (W <= 0.0) {
    opserr << "WARNING DistributingCoupling " << this->getTag()
           << ": total weight is not positive; element ignored\n";
    return;
  }

  Vector xc(ndm);
  for (int i = 0; i < nIndep; i++) {
    const Vector &xi = theNodes[1 + i]->getCrds();
    for (int k = 0; k < ndm; k++) xc(k) += weights(i) * xi(k);
  }
  for (int k = 0; k < ndm; k++) xc(k) /= W;
  const Vector &xR = theNodes[0]->getCrds();
  dRef.resize(ndm);
  for (int k = 0; k < ndm; k++) dRef(k) = xc(k) - xR(k);

  rvec.resize(nIndep, ndm);
  double Lc2 = 0.0, swr2 = 0.0;
  for (int i = 0; i < nIndep; i++) {
    const Vector &xi = theNodes[1 + i]->getCrds();
    double r2 = 0.0;
    for (int k = 0; k < ndm; k++) {
      double rk = xi(k) - xc(k);
      rvec(i, k) = rk;
      r2 += rk * rk;
    }
    swr2 += weights(i) * r2;
    if (r2 > Lc2) Lc2 = r2;
  }
  ell2 = swr2 / W;

  Icplus.resize(nrot, nrot);
  Icplus.Zero();
  Pproj.resize(nrot, nrot);
  Pproj.Zero();
  nKept = 0;

  if (ndm == 2) {
    // polar inertia; one in-plane rotation axis
    double Ic = swr2;
    double floor = 1.0e-8 * Ic;
    double absFloor = 1.0e-10 * W * Lc2;
    if (absFloor > floor) floor = absFloor;
    if (Ic > floor && Ic > 1.0e-300) {
      Icplus(0, 0) = 1.0 / Ic;
      Pproj(0, 0) = 1.0;
      nKept = 1;
    } else {
      opserr << "WARNING DistributingCoupling " << this->getTag()
             << ": the independent nodes have no in-plane spread; the reference "
             << "rotation is unconstrained\n";
    }
    valid = true;
    return;
  }

  double Ic[3][3];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++) Ic[a][b] = 0.0;
  for (int i = 0; i < nIndep; i++) {
    double rx = rvec(i, 0), ry = rvec(i, 1), rz = rvec(i, 2);
    double r2 = rx * rx + ry * ry + rz * rz;
    double w = weights(i);
    Ic[0][0] += w * (r2 - rx * rx);
    Ic[0][1] += w * (-rx * ry);
    Ic[0][2] += w * (-rx * rz);
    Ic[1][1] += w * (r2 - ry * ry);
    Ic[1][2] += w * (-ry * rz);
    Ic[2][2] += w * (r2 - rz * rz);
  }
  Ic[1][0] = Ic[0][1];
  Ic[2][0] = Ic[0][2];
  Ic[2][1] = Ic[1][2];

  double eval[3], evec[3][3];
  jacobi3(Ic, eval, evec);
  double lmax = eval[0];
  for (int k = 1; k < 3; k++)
    if (eval[k] > lmax) lmax = eval[k];
  double floor = 1.0e-8 * lmax;
  double absFloor = 1.0e-10 * W * Lc2;
  if (absFloor > floor) floor = absFloor;

  // drop axes about which the set has no spread (collinear or coincident)
  for (int k = 0; k < 3; k++) {
    if (eval[k] > floor && eval[k] > 1.0e-300) {
      nKept++;
      double inv = 1.0 / eval[k];
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
          Icplus(a, b) += inv * evec[a][k] * evec[b][k];
          Pproj(a, b) += evec[a][k] * evec[b][k];
        }
    } else {
      opserr << "WARNING DistributingCoupling " << this->getTag()
             << ": the reference rotation about (" << evec[0][k] << ", " << evec[1][k]
             << ", " << evec[2][k] << ") is unconstrained (degenerate independent set)\n";
    }
  }
  valid = true;
}

// (r_i x u)_r = sum_j crossOp(i, r, j) u_j
double DistributingCoupling::crossOp(int i, int r, int j) const
{
  if (ndm == 3) {
    double rx = rvec(i, 0), ry = rvec(i, 1), rz = rvec(i, 2);
    switch (r) {
    case 0: return (j == 1) ? -rz : (j == 2) ?  ry : 0.0;   // ry uz - rz uy
    case 1: return (j == 0) ?  rz : (j == 2) ? -rx : 0.0;   // rz ux - rx uz
    case 2: return (j == 0) ? -ry : (j == 1) ?  rx : 0.0;   // rx uy - ry ux
    default: return 0.0;
    }
  }
  double rx = rvec(i, 0), ry = rvec(i, 1);
  return (j == 0) ? -ry : (j == 1) ? rx : 0.0;
}

// (theta x (x_c - x_R))_k = sum_r transOp(k, r) theta_r
double DistributingCoupling::transOp(int k, int r) const
{
  if (ndm == 3) {
    double dx = dRef(0), dy = dRef(1), dz = dRef(2);
    switch (k) {
    case 0: return (r == 1) ?  dz : (r == 2) ? -dy : 0.0;
    case 1: return (r == 0) ? -dz : (r == 2) ?  dx : 0.0;
    case 2: return (r == 0) ?  dy : (r == 1) ? -dx : 0.0;
    default: return 0.0;
    }
  }
  return (k == 0) ? -dRef(1) : dRef(0);
}

// B rows (RBE3 kinematics, Nastran QRG RBE3; Abaqus distributing coupling):
//   0..ndm-1:      u_R + T theta_R - sum (w_i/W) u_i
//   ndm..ndm+nrot: P theta_R - I_c^+ sum w_i C_i u_i
void DistributingCoupling::buildB(void)
{
  B->Zero();
  int oR = dofOffset(0);
  for (int k = 0; k < ndm; k++) {
    (*B)(k, oR + k) += 1.0;
    for (int r = 0; r < nrot; r++)
      (*B)(k, oR + ndm + r) += this->transOp(k, r);
    for (int i = 0; i < nIndep; i++)
      (*B)(k, dofOffset(1 + i) + k) += -(weights(i) / W);
  }
  for (int r = 0; r < nrot; r++) {
    int rowR = ndm + r;
    for (int r2 = 0; r2 < nrot; r2++)
      (*B)(rowR, oR + ndm + r2) += Pproj(r, r2);
    for (int i = 0; i < nIndep; i++)
      for (int j = 0; j < ndm; j++) {
        double v = 0.0;
        for (int s = 0; s < nrot; s++) v += Icplus(r, s) * this->crossOp(i, s, j);
        (*B)(rowR, dofOffset(1 + i) + j) += -weights(i) * v;
      }
  }
}

// K_t from -k auto, default K_r = K_t l^2, and a conditioning check of a
// numeric K_t against the -host element. Deferred until the -host element
// exists in the domain.
void DistributingCoupling::resolvePenalty(void)
{
  if (ktResolved) return;

  Element *host = 0;
  if (hostEleTag >= 0) {
    Domain *theDomain = this->getDomain();
    if (theDomain != 0) host = theDomain->getElement(hostEleTag);
    if (host == 0) {
      if (!krUser) Kr = Kt * ell2;    // provisional
      return;
    }
  }

  if (ktAuto && host != 0) {
    double scale = maxAbsDiag(host->getInitialStiff());
    if (scale > 0.0)
      Kt = kAlpha * scale;
    else
      opserr << "WARNING DistributingCoupling " << this->getTag()
             << ": -host element has a zero stiffness scale; K_t = " << Kt << " kept\n";
  }
  if (!krUser) Kr = Kt * ell2;

  if (!ktAuto && host != 0 && Kt > 0.0) {
    double scale = maxAbsDiag(host->getInitialStiff());
    if (scale > 0.0 && Kt > 1.0e6 * scale)
      opserr << "WARNING DistributingCoupling " << this->getTag() << ": -k " << Kt
             << " is " << Kt / scale << " times the -host element stiffness scale ("
             << scale << "); the system may be ill-conditioned. Consider 1e2 to 1e4 "
             << "times the host stiffness, or -enforce al with a moderate K_t\n";
  }
  ktResolved = true;
}

// ===========================================================================
//  state
// ===========================================================================
int DistributingCoupling::commitState(void)
{
  this->resolvePenalty();
  if (hostEleTag >= 0 && !ktResolved && !hostMissWarned) {
    hostMissWarned = true;
    opserr << "WARNING DistributingCoupling " << this->getTag() << ": -host element "
           << hostEleTag << " is not in the domain; K_t = " << Kt << " is used\n";
  }
  // augmented Lagrangian: one Uzawa update per committed step
  if (enforce == 1 && valid) {
    Vector g(ndm + nrot);
    this->computeGap(g);
    for (int k = 0; k < ndm; k++) lambda(k) += Kt * g(k);
    for (int r = 0; r < nrot; r++) lambda_r(r) += Kr * g(ndm + r);
  }
  return this->Element::commitState();
}

int DistributingCoupling::revertToLastCommit(void)
{
  // the multipliers change only in commitState
  return 0;
}

int DistributingCoupling::revertToStart(void)
{
  lambda.Zero();
  lambda_r.Zero();
  return 0;
}

int DistributingCoupling::update(void)
{
  this->resolvePenalty();
  return 0;
}

void DistributingCoupling::computeGap(Vector &g)
{
  int nGap = ndm + nrot;
  if (g.Size() != nGap) g.resize(nGap);
  g.Zero();
  if (!valid) return;
  Vector uFull(nDOF);
  for (int p = 0; p < 1 + nIndep; p++) {
    const Vector &u = theNodes[p]->getTrialDisp();
    for (int d = 0; d < nodeNdf(p); d++)
      uFull(dofOffset(p) + d) = u(d);
  }
  for (int row = 0; row < nGap; row++) {
    double s = 0.0;
    for (int c = 0; c < nDOF; c++) s += (*B)(row, c) * uFull(c);
    g(row) = s;
  }
  if (g0Computed) {
    for (int k = 0; k < ndm; k++) g(k) -= g0(k);
    for (int r = 0; r < nrot; r++) g(ndm + r) -= gr0(r);
  }
}

// the element is created stress free in the current configuration unless
// -absolute is given
void DistributingCoupling::captureInitialGap(void)
{
  if (!initGapCapture || g0Computed) return;
  Vector g(ndm + nrot);
  this->computeGap(g);
  for (int k = 0; k < ndm; k++) g0(k) = g(k);
  for (int r = 0; r < nrot; r++) gr0(r) = g(ndm + r);
  g0Computed = true;
}

// ===========================================================================
//  matrices and forces
// ===========================================================================
const Vector &DistributingCoupling::getResistingForce(void)
{
  this->resolvePenalty();
  P->Zero();
  if (!valid) return *P;
  int nGap = ndm + nrot;
  Vector g(nGap);
  this->computeGap(g);
  for (int row = 0; row < nGap; row++) {
    double t;
    if (row < ndm)
      t = Kt * g(row) + (enforce == 1 ? lambda(row) : 0.0);
    else
      t = Kr * g(row) + (enforce == 1 ? lambda_r(row - ndm) : 0.0);
    if (t == 0.0) continue;
    for (int c = 0; c < nDOF; c++) {
      double b = (*B)(row, c);
      if (b != 0.0) (*P)(c) += b * t;
    }
  }
  return *P;
}

const Vector &DistributingCoupling::getResistingForceIncInertia(void)
{
  return this->getResistingForce();    // no mass, no damping
}

const Matrix &DistributingCoupling::getTangentStiff(void)
{
  this->resolvePenalty();
  K->Zero();
  if (!valid) return *K;
  int nGap = ndm + nrot;
  for (int row = 0; row < nGap; row++) {
    double dval = (row < ndm) ? Kt : Kr;
    if (dval == 0.0) continue;
    for (int c1 = 0; c1 < nDOF; c1++) {
      double b1 = (*B)(row, c1);
      if (b1 == 0.0) continue;
      double db = dval * b1;
      for (int c2 = 0; c2 < nDOF; c2++) {
        double b2 = (*B)(row, c2);
        if (b2 != 0.0) (*K)(c1, c2) += db * b2;
      }
    }
  }
  return *K;
}

const Matrix &DistributingCoupling::getInitialStiff(void)
{
  return this->getTangentStiff();      // B, K_t and K_r do not depend on the state
}

const Matrix &DistributingCoupling::getMass(void)
{
  M0->Zero();
  return *M0;
}

// Element::getDamp indexes a damping-matrix pool that is only allocated by
// Element::setRayleighDampingFactors; the zero matrix is returned directly.
const Matrix &DistributingCoupling::getDamp(void)
{
  M0->Zero();
  return *M0;
}

int DistributingCoupling::setRayleighDampingFactors(double, double, double, double)
{
  return this->Element::setRayleighDampingFactors(0.0, 0.0, 0.0, 0.0);
}

// ===========================================================================
//  serialization (geometry and B are rebuilt in setDomain)
// ===========================================================================
int DistributingCoupling::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector hdr(12);
  hdr(0) = this->getTag();
  hdr(1) = ndm;
  hdr(2) = nIndep;
  hdr(3) = Kt;
  hdr(4) = Kr;
  hdr(5) = krUser ? 1.0 : 0.0;
  hdr(6) = ktAuto ? 1.0 : 0.0;
  hdr(7) = kAlpha;
  hdr(8) = hostEleTag;
  hdr(9) = enforce;
  hdr(10) = initGapCapture ? 1.0 : 0.0;
  hdr(11) = g0Computed ? 1.0 : 0.0;
  if (theChannel.sendVector(dbTag, commitTag, hdr) < 0) {
    opserr << "DistributingCoupling::sendSelf - failed to send header\n";
    return -1;
  }
  if (theChannel.sendID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "DistributingCoupling::sendSelf - failed to send node tags\n";
    return -1;
  }
  // weights(N), lambda(ndm), lambda_r(nrot), g0(ndm), gr0(nrot)
  Vector payload(nIndep + 2 * ndm + 2 * nrot);
  int off = 0;
  for (int i = 0; i < nIndep; i++) payload(off + i) = weights(i);
  off += nIndep;
  for (int k = 0; k < ndm; k++) payload(off + k) = (lambda.Size() == ndm) ? lambda(k) : 0.0;
  off += ndm;
  for (int r = 0; r < nrot; r++) payload(off + r) = (lambda_r.Size() == nrot) ? lambda_r(r) : 0.0;
  off += nrot;
  for (int k = 0; k < ndm; k++) payload(off + k) = (g0.Size() == ndm) ? g0(k) : 0.0;
  off += ndm;
  for (int r = 0; r < nrot; r++) payload(off + r) = (gr0.Size() == nrot) ? gr0(r) : 0.0;
  if (theChannel.sendVector(dbTag, commitTag, payload) < 0) {
    opserr << "DistributingCoupling::sendSelf - failed to send state\n";
    return -1;
  }
  return 0;
}

int DistributingCoupling::recvSelf(int commitTag, Channel &theChannel,
                                   FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector hdr(12);
  if (theChannel.recvVector(dbTag, commitTag, hdr) < 0) {
    opserr << "DistributingCoupling::recvSelf - failed to receive header\n";
    return -1;
  }
  this->setTag((int)hdr(0));
  ndm = (int)hdr(1);
  nIndep = (int)hdr(2);
  Kt = hdr(3);
  Kr = hdr(4);
  krUser = (hdr(5) != 0.0);
  ktAuto = (hdr(6) != 0.0);
  kAlpha = hdr(7);
  hostEleTag = (int)hdr(8);
  enforce = (int)hdr(9);
  initGapCapture = (hdr(10) != 0.0);
  g0Computed = (hdr(11) != 0.0);
  nrot = (ndm == 3) ? 3 : 1;
  ktResolved = false;
  hostMissWarned = false;
  nDOF = 0;
  valid = false;

  connectedNodes.resize(1 + nIndep);
  if (theChannel.recvID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "DistributingCoupling::recvSelf - failed to receive node tags\n";
    return -1;
  }
  Vector payload(nIndep + 2 * ndm + 2 * nrot);
  if (theChannel.recvVector(dbTag, commitTag, payload) < 0) {
    opserr << "DistributingCoupling::recvSelf - failed to receive state\n";
    return -1;
  }
  int off = 0;
  weights.resize(nIndep);
  for (int i = 0; i < nIndep; i++) weights(i) = payload(off + i);
  off += nIndep;
  lambda.resize(ndm);
  for (int k = 0; k < ndm; k++) lambda(k) = payload(off + k);
  off += ndm;
  lambda_r.resize(nrot);
  for (int r = 0; r < nrot; r++) lambda_r(r) = payload(off + r);
  off += nrot;
  g0.resize(ndm);
  for (int k = 0; k < ndm; k++) g0(k) = payload(off + k);
  off += ndm;
  gr0.resize(nrot);
  for (int r = 0; r < nrot; r++) gr0(r) = payload(off + r);

  dRef.resize(ndm);
  Icplus.resize(nrot, nrot);
  Pproj.resize(nrot, nrot);
  nodeNdf.resize(1 + nIndep);
  dofOffset.resize(1 + nIndep);
  if (theNodes != 0) delete[] theNodes;
  theNodes = new Node *[1 + nIndep];
  for (int i = 0; i < 1 + nIndep; i++) theNodes[i] = 0;
  return 0;
}

void DistributingCoupling::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"DistributingCoupling\", ";
    s << "\"nodes\": [";
    for (int i = 0; i < 1 + nIndep; i++) {
      if (i > 0) s << ", ";
      s << connectedNodes(i);
    }
    s << "], \"kt\": " << Kt << ", \"kr\": " << Kr << "}";
    return;
  }
  s << "DistributingCoupling, tag: " << this->getTag() << "\n";
  s << "  reference node: " << connectedNodes(0) << "  independent nodes:";
  for (int i = 0; i < nIndep; i++) s << " " << connectedNodes(1 + i);
  s << "\n  weights:";
  for (int i = 0; i < nIndep; i++) s << " " << weights(i);
  s << "\n  K_t: " << Kt << "  K_r: " << Kr
    << "  enforce: " << (enforce == 1 ? "al" : "penalty")
    << "  constrained rotation axes: " << nKept << "/" << nrot
    << (initGapCapture ? "" : "  (absolute)") << "\n";
}

// ===========================================================================
//  responses
// ===========================================================================
Response *DistributingCoupling::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  if (argc < 1) return 0;
  static const char *const tForce[] = {"couplingForce", "couplingForces", 0};
  static const char *const tGap[] = {"gap", 0};
  static const char *const tKt[] = {"penalty", "kt", "k", 0};
  static const char *const tKr[] = {"rotPenalty", "kr", 0};
  static const char *const tLam[] = {"lambda", 0};
  static const char *const tLamR[] = {"lambdaR", 0};
  static const char *const tAx[] = {"rotationAxes", "nKept", 0};
  static const char *const tE[] = {"penaltyEnergy", 0};
  static const char *const tV[] = {"constraintViolation", 0};
  int nGap = ndm + nrot;
  if (tokenIs(argv[0], tForce)) return new ElementResponse(this, 1, Vector(nGap));
  if (tokenIs(argv[0], tGap))   return new ElementResponse(this, 2, Vector(nGap));
  if (tokenIs(argv[0], tKt))    return new ElementResponse(this, 3, 0.0);
  if (tokenIs(argv[0], tKr))    return new ElementResponse(this, 4, 0.0);
  if (tokenIs(argv[0], tLam))   return new ElementResponse(this, 5, Vector(ndm));
  if (tokenIs(argv[0], tLamR))  return new ElementResponse(this, 6, Vector(nrot));
  if (tokenIs(argv[0], tAx))    return new ElementResponse(this, 7, 0.0);
  if (tokenIs(argv[0], tE))     return new ElementResponse(this, 8, 0.0);
  if (tokenIs(argv[0], tV))     return new ElementResponse(this, 9, 0.0);
  return this->Element::setResponse(argv, argc, s);
}

int DistributingCoupling::getResponse(int responseID, Information &eleInfo)
{
  this->resolvePenalty();
  int nGap = ndm + nrot;
  switch (responseID) {
  case 1: {                     // tie force D g (+ lambda), translation then rotation
    Vector g(nGap);
    this->computeGap(g);
    for (int k = 0; k < ndm; k++)
      g(k) = Kt * g(k) + (enforce == 1 ? lambda(k) : 0.0);
    for (int r = 0; r < nrot; r++)
      g(ndm + r) = Kr * g(ndm + r) + (enforce == 1 ? lambda_r(r) : 0.0);
    return eleInfo.setVector(g);
  }
  case 2: {
    Vector g(nGap);
    this->computeGap(g);
    return eleInfo.setVector(g);
  }
  case 3: return eleInfo.setDouble(Kt);
  case 4: return eleInfo.setDouble(Kr);
  case 5: return eleInfo.setVector(lambda);
  case 6: return eleInfo.setVector(lambda_r);
  case 7: return eleInfo.setDouble((double)nKept);
  case 8: {                     // 1/2 (K_t |g_t|^2 + K_r |g_r|^2)
    Vector g(nGap);
    this->computeGap(g);
    double Et = 0.0, Er = 0.0;
    for (int k = 0; k < ndm; k++) Et += g(k) * g(k);
    for (int r = 0; r < nrot; r++) Er += g(ndm + r) * g(ndm + r);
    return eleInfo.setDouble(0.5 * (Kt * Et + Kr * Er));
  }
  case 9: {                     // ||g||
    Vector g(nGap);
    this->computeGap(g);
    return eleInfo.setDouble(g.Norm());
  }
  default:
    return this->Element::getResponse(responseID, eleInfo);
  }
}
