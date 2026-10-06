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
// Description: KinematicCoupling, a rigid kinematic coupling (RBE2-type).
// See KinematicCoupling.h for the formulation and references
// (MSC Nastran RBE2; Abaqus *COUPLING, *KINEMATIC).

#include <KinematicCoupling.h>
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

} // namespace

void *OPS_KinematicCoupling(void)
{
  int ndm = OPS_GetNDM();
  if (ndm != 2 && ndm != 3) {
    opserr << "WARNING KinematicCoupling: model ndm must be 2 or 3\n";
    return 0;
  }
  int nrot = (ndm == 3) ? 3 : 1;
  int cmax = ndm + nrot;

  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING insufficient args\n"
           << "Want: element KinematicCoupling tag refNode N s1 ... sN "
           << "<-dof c1 ... cK> <-k Kt|auto> <-kAlpha a> <-host eleTag> <-kr Kr> "
           << "<-enforce penalty|al> <-absolute>\n";
    return 0;
  }

  int idata[3];                    // tag, refNode, N
  int n = 3;
  if (OPS_GetIntInput(&n, idata) < 0) {
    opserr << "WARNING KinematicCoupling: invalid tag, refNode or N\n";
    return 0;
  }
  int tag = idata[0], refNode = idata[1], N = idata[2];
  if (N < 1) {
    opserr << "WARNING KinematicCoupling " << tag
           << ": N (number of slave nodes) must be >= 1\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < N) {
    opserr << "WARNING KinematicCoupling " << tag << ": need " << N
           << " slave node tags\n";
    return 0;
  }
  ID slaves(N);
  for (int i = 0; i < N; i++) {
    int h;
    n = 1;
    if (OPS_GetIntInput(&n, &h) < 0) {
      opserr << "WARNING KinematicCoupling " << tag << ": invalid slave node "
             << i + 1 << "\n";
      return 0;
    }
    slaves(i) = h;
  }

  ID dofSel(0);
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
    if (strcmp(opt, "-dof") == 0) {
      // read integer components until the next flag
      int nsel = 0;
      while (OPS_GetNumRemainingInputArgs() > 0) {
        char tok[64];
        OPS_GetStringFromAll(tok, sizeof(tok));
        char *endp = 0;
        long cv = strtol(tok, &endp, 10);
        if (endp == tok || *endp != '\0') {
          OPS_ResetCurrentInputArg(-1);
          break;
        }
        if (cv < 1 || cv > cmax) {
          opserr << "WARNING KinematicCoupling " << tag << ": -dof component " << (int)cv
                 << " out of range; valid 1.." << cmax << " in " << ndm << "D\n";
          return 0;
        }
        dofSel[nsel++] = (int)cv;
      }
      if (nsel == 0) {
        opserr << "WARNING KinematicCoupling " << tag
               << ": -dof needs at least one component\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-k") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING KinematicCoupling " << tag << ": -k wants a value or 'auto'\n";
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
          opserr << "WARNING KinematicCoupling " << tag
                 << ": -k wants a positive number or 'auto', got '" << kTok << "'\n";
          return 0;
        }
        Kt = kv;
      }
    }
    else if (strcmp(opt, "-kAlpha") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &kAlpha) < 0 || kAlpha <= 0.0) {
        opserr << "WARNING KinematicCoupling " << tag << ": -kAlpha wants a positive value\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-host") == 0) {
      n = 1;
      if (OPS_GetIntInput(&n, &hostEleTag) < 0) {
        opserr << "WARNING KinematicCoupling " << tag << ": -host wants an element tag\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-kr") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &Kr) < 0 || Kr < 0.0) {
        opserr << "WARNING KinematicCoupling " << tag << ": -kr wants a non-negative value\n";
        return 0;
      }
      krUser = true;
    }
    else if (strcmp(opt, "-enforce") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING KinematicCoupling " << tag << ": -enforce wants penalty or al\n";
        return 0;
      }
      const char *mode = OPS_GetString();
      if (strcmp(mode, "penalty") == 0)
        enforce = 0;
      else if (strcmp(mode, "al") == 0)
        enforce = 1;
      else {
        opserr << "WARNING KinematicCoupling " << tag << ": unknown -enforce '" << mode
               << "' (want penalty or al)\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-absolute") == 0) {
      initGapCapture = false;
    }
    else {
      opserr << "WARNING KinematicCoupling " << tag << ": unknown option '" << opt << "'\n";
      return 0;
    }
  }

  if (ktAuto && hostEleTag < 0) {
    opserr << "WARNING KinematicCoupling " << tag
           << ": -k auto requires a representative -host element\n";
    return 0;
  }

  // With the default component list every component 1..ndm+nrot present on
  // a slave is tied. A slave whose ndf is neither ndm nor ndm+nrot (for
  // example a u-p node carrying a pressure DOF) would have that extra DOF
  // tied to a reference rotation, so an explicit -dof list is required.
  if (dofSel.Size() == 0) {
    Domain *dom = OPS_GetDomain();
    for (int i = 0; dom != 0 && i < N; i++) {
      Node *sn = dom->getNode(slaves(i));
      if (sn == 0) continue;
      int sndf = sn->getNumberDOF();
      if (sndf != ndm && sndf != ndm + nrot) {
        opserr << "WARNING KinematicCoupling " << tag << ": slave node " << slaves(i)
               << " has ndf = " << sndf << ", which is neither " << ndm << " nor "
               << ndm + nrot << "; give the tied components with -dof\n";
        return 0;
      }
    }
  }

  return new KinematicCoupling(tag, ndm, refNode, slaves, dofSel, Kt, Kr, krUser,
                               enforce, kAlpha, hostEleTag, ktAuto, initGapCapture);
}

// ===========================================================================
//  construction
// ===========================================================================
KinematicCoupling::KinematicCoupling(int tag, int ndm_, int refNode,
                                     const ID &slaveNodes, const ID &dofSel_,
                                     double kt, double kr, bool krUser_, int enforce_,
                                     double kAlpha_, int hostEleTag_, bool ktAuto_,
                                     bool initGapCapture_)
  : Element(tag, ELE_TAG_KinematicCoupling),
    ndm(ndm_), nrot((ndm_ == 3) ? 3 : 1), nSlave(slaveNodes.Size()),
    connectedNodes(1 + slaveNodes.Size()), dofSel(dofSel_),
    Kt(kt), Kr(kr), krUser(krUser_), ktAuto(ktAuto_), kAlpha(kAlpha_),
    hostEleTag(hostEleTag_), ktResolved(false), hostMissWarned(false), ell2(0.0),
    enforce(enforce_), lambdaAL(), hasRefRot(false),
    valid(false), dvec(), nGap(0), gapNode(), gapDof(), gapIsRot(),
    nDOF(0), nodeNdf(1 + slaveNodes.Size()), dofOffset(1 + slaveNodes.Size()),
    B(0), initGapCapture(initGapCapture_), g0Computed(false), g0(),
    theNodes(0), K(0), P(0), M0(0)
{
  connectedNodes(0) = refNode;
  for (int i = 0; i < nSlave; i++)
    connectedNodes(1 + i) = slaveNodes(i);
  theNodes = new Node *[1 + nSlave];
  for (int i = 0; i < 1 + nSlave; i++)
    theNodes[i] = 0;
}

KinematicCoupling::KinematicCoupling()
  : Element(0, ELE_TAG_KinematicCoupling),
    ndm(0), nrot(0), nSlave(0), connectedNodes(), dofSel(),
    Kt(0.0), Kr(0.0), krUser(false), ktAuto(false), kAlpha(0.0),
    hostEleTag(-1), ktResolved(false), hostMissWarned(false), ell2(0.0),
    enforce(0), lambdaAL(), hasRefRot(false),
    valid(false), dvec(), nGap(0), gapNode(), gapDof(), gapIsRot(),
    nDOF(0), nodeNdf(), dofOffset(),
    B(0), initGapCapture(true), g0Computed(false), g0(),
    theNodes(0), K(0), P(0), M0(0)
{
}

KinematicCoupling::~KinematicCoupling()
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
int KinematicCoupling::getNumExternalNodes(void) const { return 1 + nSlave; }
const ID &KinematicCoupling::getExternalNodes(void) { return connectedNodes; }
Node **KinematicCoupling::getNodePtrs(void) { return theNodes; }
int KinematicCoupling::getNumDOF(void) { return nDOF; }

void KinematicCoupling::allocate(void)
{
  if (K != 0) delete K;
  K = new Matrix(nDOF, nDOF);
  if (P != 0) delete P;
  P = new Vector(nDOF);
  if (M0 != 0) delete M0;
  M0 = new Matrix(nDOF, nDOF);
}

void KinematicCoupling::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    if (theNodes != 0)
      for (int i = 0; i < 1 + nSlave; i++) theNodes[i] = 0;
    return;
  }

  nodeNdf.resize(1 + nSlave);
  dofOffset.resize(1 + nSlave);
  int pos = 0;
  bool refused = false;
  for (int i = 0; i < 1 + nSlave; i++) {
    theNodes[i] = theDomain->getNode(connectedNodes(i));
    if (theNodes[i] == 0) {
      opserr << "WARNING KinematicCoupling " << this->getTag() << ": node "
             << connectedNodes(i) << " not found\n";
      refused = true;
      nodeNdf(i) = 0;
      dofOffset(i) = pos;
      continue;
    }
    int ndf = theNodes[i]->getNumberDOF();
    if (ndf < ndm) {
      opserr << "WARNING KinematicCoupling " << this->getTag() << ": "
             << (i == 0 ? "reference" : "slave") << " node " << connectedNodes(i)
             << " has " << ndf << " DOFs, needs at least " << ndm << "\n";
      refused = true;
    }
    nodeNdf(i) = ndf;
    dofOffset(i) = pos;
    pos += ndf;
  }
  nDOF = pos;

  // a slave equal to the reference or listed twice is refused
  for (int i = 0; i < nSlave && !refused; i++) {
    if (connectedNodes(1 + i) == connectedNodes(0)) {
      opserr << "WARNING KinematicCoupling " << this->getTag() << ": slave node "
             << connectedNodes(1 + i) << " is the reference node; element ignored\n";
      refused = true;
      break;
    }
    for (int j = i + 1; j < nSlave; j++)
      if (connectedNodes(1 + j) == connectedNodes(1 + i)) {
        opserr << "WARNING KinematicCoupling " << this->getTag() << ": slave node "
               << connectedNodes(1 + i) << " is listed more than once; element ignored\n";
        refused = true;
        break;
      }
  }

  if (refused) {
    valid = false;
    nGap = 0;
  } else {
    this->resolveGeometry();
  }

  // An ill-posed element is kept inert (zero B, K, P) rather than leaving
  // unallocated matrices for the assembler.
  if (B != 0) delete B;
  B = new Matrix((nGap > 0 ? nGap : 1), (nDOF > 0 ? nDOF : 1));
  B->Zero();
  this->allocate();

  // keep multipliers and initial gap restored by recvSelf (sizes match)
  bool lamRestored = (nGap > 0 && lambdaAL.Size() == nGap);
  bool g0Restored = (nGap > 0 && g0.Size() == nGap && g0Computed);
  lambdaAL.resize(nGap > 0 ? nGap : 1);
  if (!lamRestored || !valid) lambdaAL.Zero();
  g0.resize(nGap > 0 ? nGap : 1);
  if (!g0Restored) {
    g0.Zero();
    g0Computed = false;
  }
  if (valid) {
    this->buildB();
    this->captureInitialGap();
  }
  this->DomainComponent::setDomain(theDomain);
}

// d_i, the per-slave tied-component layout and the rotation length scale.
void KinematicCoupling::resolveGeometry(void)
{
  valid = false;
  nGap = 0;

  hasRefRot = (nodeNdf(0) >= ndm + nrot);

  dvec.resize(nSlave, ndm);
  const Vector &xR = theNodes[0]->getCrds();
  double sumD2 = 0.0, maxD2 = 0.0;
  for (int i = 0; i < nSlave; i++) {
    const Vector &xi = theNodes[1 + i]->getCrds();
    double d2 = 0.0;
    for (int k = 0; k < ndm; k++) {
      double dk = xi(k) - xR(k);
      dvec(i, k) = dk;
      d2 += dk * dk;
    }
    sumD2 += d2;
    if (d2 > maxD2) maxD2 = d2;
  }

  bool useDefault = (dofSel.Size() == 0);
  bool reqRot = false;
  if (!useDefault)
    for (int s = 0; s < dofSel.Size(); s++)
      if (dofSel(s) > ndm) reqRot = true;

  if (reqRot && !hasRefRot) {
    opserr << "WARNING KinematicCoupling " << this->getTag()
           << ": -dof requests a slave rotation but the reference node has no rotation "
           << "DOFs; element ignored\n";
    return;
  }
  if (!hasRefRot)
    opserr << "KinematicCoupling " << this->getTag()
           << ": reference node has no rotation DOFs; only translations are tied "
           << "(no moment transfer)\n";

  // one row per (slave, component); component c: 1..ndm translation,
  // ndm+1..ndm+nrot rotation, node DOF index c-1
  int count = 0;
  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      nGap = count;
      gapNode.resize(nGap > 0 ? nGap : 1);
      gapDof.resize(nGap > 0 ? nGap : 1);
      gapIsRot.resize(nGap > 0 ? nGap : 1);
      count = 0;
    }
    for (int i = 0; i < nSlave; i++) {
      int sndf = nodeNdf(1 + i);
      int cLast = ndm + (hasRefRot ? nrot : 0);
      for (int c = 1; c <= cLast; c++) {
        bool wanted = useDefault;
        if (!useDefault)
          for (int s = 0; s < dofSel.Size(); s++)
            if (dofSel(s) == c) { wanted = true; break; }
        if (!wanted) continue;
        if ((c - 1) >= sndf) {
          if (pass == 0 && !useDefault)
            opserr << "KinematicCoupling " << this->getTag() << ": slave node "
                   << connectedNodes(1 + i) << " has no component " << c
                   << "; skipped\n";
          continue;
        }
        if (pass == 1) {
          gapNode(count) = 1 + i;
          gapDof(count) = c - 1;
          gapIsRot(count) = (c > ndm) ? 1 : 0;
        }
        count++;
      }
    }
  }

  if (nGap == 0) {
    opserr << "WARNING KinematicCoupling " << this->getTag()
           << ": no slave component is tied; element ignored\n";
    return;
  }

  // l^2 = max(mean |d_i|^2, max |d_i|^2); unit length if all slaves coincide
  // with R and a rotation is tied, so the rotational penalty never vanishes
  ell2 = (nSlave > 0) ? sumD2 / nSlave : 0.0;
  if (ell2 < maxD2) ell2 = maxD2;
  if (ell2 <= 0.0) {
    bool anyRot = false;
    for (int r = 0; r < nGap; r++)
      if (gapIsRot(r)) { anyRot = true; break; }
    if (anyRot) {
      ell2 = 1.0;
      opserr << "KinematicCoupling " << this->getTag()
             << ": all slaves coincide with the reference node; unit length used for "
             << "the default rotational penalty (set -kr to override)\n";
    }
  }
  valid = true;
}

// d (theta x d_i)_k / d theta_r
double KinematicCoupling::transOp(int i, int k, int r) const
{
  if (ndm == 3) {
    double dx = dvec(i, 0), dy = dvec(i, 1), dz = dvec(i, 2);
    switch (k) {
    case 0: return (r == 1) ?  dz : (r == 2) ? -dy : 0.0;   // ty dz - tz dy
    case 1: return (r == 0) ? -dz : (r == 2) ?  dx : 0.0;   // tz dx - tx dz
    case 2: return (r == 0) ?  dy : (r == 1) ? -dx : 0.0;   // tx dy - ty dx
    default: return 0.0;
    }
  }
  // 2D: theta x d = theta (-d_y, d_x)
  return (k == 0) ? -dvec(i, 1) : dvec(i, 0);
}

// B rows (RBE2 kinematics, Nastran QRG RBE2):
//   translation: g = u_i[j] - u_R[j] - sum_r d(theta x d_i)_j/d theta_r theta_R[r]
//   rotation:    g = theta_i[j] - theta_R[j]
void KinematicCoupling::buildB(void)
{
  B->Zero();
  int oR = dofOffset(0);
  for (int row = 0; row < nGap; row++) {
    int i = gapNode(row) - 1;
    int oS = dofOffset(gapNode(row));
    int j = gapDof(row);
    (*B)(row, oS + j) += 1.0;
    (*B)(row, oR + j) += -1.0;
    if (gapIsRot(row) == 0 && hasRefRot)
      for (int r = 0; r < nrot; r++)
        (*B)(row, oR + ndm + r) += -this->transOp(i, j, r);
  }
}

// K_t from -k auto, default K_r = K_t l^2, and a conditioning check of a
// numeric K_t against the -host element. Deferred until the -host element
// exists in the domain (it may be defined after the coupling).
void KinematicCoupling::resolvePenalty(void)
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
      opserr << "WARNING KinematicCoupling " << this->getTag()
             << ": -host element has a zero stiffness scale; K_t = " << Kt << " kept\n";
  }
  if (!krUser) Kr = Kt * ell2;

  // The tie error of a penalty constraint falls as 1/K_t while the condition
  // number of the system grows with K_t.
  if (!ktAuto && host != 0 && Kt > 0.0) {
    double scale = maxAbsDiag(host->getInitialStiff());
    if (scale > 0.0 && Kt > 1.0e6 * scale)
      opserr << "WARNING KinematicCoupling " << this->getTag() << ": -k " << Kt
             << " is " << Kt / scale << " times the -host element stiffness scale ("
             << scale << "); the system may be ill-conditioned. Consider 1e2 to 1e4 "
             << "times the host stiffness, or -enforce al with a moderate K_t\n";
  }
  ktResolved = true;
}

// ===========================================================================
//  state
// ===========================================================================
int KinematicCoupling::commitState(void)
{
  this->resolvePenalty();
  if (hostEleTag >= 0 && !ktResolved && !hostMissWarned) {
    hostMissWarned = true;
    opserr << "WARNING KinematicCoupling " << this->getTag() << ": -host element "
           << hostEleTag << " is not in the domain; K_t = " << Kt << " is used\n";
  }
  // augmented Lagrangian: one Uzawa update per committed step
  if (enforce == 1 && valid) {
    Vector g(nGap);
    this->computeGap(g);
    for (int row = 0; row < nGap; row++)
      lambdaAL(row) += this->rowPenalty(row) * g(row);
  }
  return this->Element::commitState();
}

int KinematicCoupling::revertToLastCommit(void)
{
  // the multipliers change only in commitState
  return 0;
}

int KinematicCoupling::revertToStart(void)
{
  lambdaAL.Zero();
  return 0;
}

int KinematicCoupling::update(void)
{
  this->resolvePenalty();
  return 0;
}

void KinematicCoupling::computeGap(Vector &g)
{
  if (g.Size() != nGap) g.resize(nGap);
  g.Zero();
  if (!valid) return;
  Vector uFull(nDOF);
  for (int p = 0; p < 1 + nSlave; p++) {
    const Vector &u = theNodes[p]->getTrialDisp();
    for (int d = 0; d < nodeNdf(p); d++)
      uFull(dofOffset(p) + d) = u(d);
  }
  for (int row = 0; row < nGap; row++) {
    double s = 0.0;
    for (int c = 0; c < nDOF; c++)
      s += (*B)(row, c) * uFull(c);
    g(row) = s;
  }
  if (g0Computed)
    for (int row = 0; row < nGap; row++) g(row) -= g0(row);
}

// the element is created stress free in the current configuration unless
// -absolute is given
void KinematicCoupling::captureInitialGap(void)
{
  if (!initGapCapture || g0Computed) return;
  Vector g(nGap);
  this->computeGap(g);
  for (int row = 0; row < nGap; row++) g0(row) = g(row);
  g0Computed = true;
}

// ===========================================================================
//  matrices and forces
// ===========================================================================
const Vector &KinematicCoupling::getResistingForce(void)
{
  this->resolvePenalty();
  P->Zero();
  if (!valid) return *P;
  Vector g(nGap);
  this->computeGap(g);
  for (int row = 0; row < nGap; row++) {
    double t = this->rowPenalty(row) * g(row);
    if (enforce == 1) t += lambdaAL(row);
    if (t == 0.0) continue;
    for (int c = 0; c < nDOF; c++) {
      double b = (*B)(row, c);
      if (b != 0.0) (*P)(c) += b * t;
    }
  }
  return *P;
}

const Vector &KinematicCoupling::getResistingForceIncInertia(void)
{
  return this->getResistingForce();    // no mass, no damping
}

const Matrix &KinematicCoupling::getTangentStiff(void)
{
  this->resolvePenalty();
  K->Zero();
  if (!valid) return *K;
  for (int row = 0; row < nGap; row++) {
    double dval = this->rowPenalty(row);
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

const Matrix &KinematicCoupling::getInitialStiff(void)
{
  return this->getTangentStiff();      // B, K_t and K_r do not depend on the state
}

const Matrix &KinematicCoupling::getMass(void)
{
  M0->Zero();
  return *M0;
}

// Element::getDamp indexes a damping-matrix pool that is only allocated by
// Element::setRayleighDampingFactors; the zero matrix is returned directly.
const Matrix &KinematicCoupling::getDamp(void)
{
  M0->Zero();
  return *M0;
}

int KinematicCoupling::setRayleighDampingFactors(double, double, double, double)
{
  return this->Element::setRayleighDampingFactors(0.0, 0.0, 0.0, 0.0);
}

// ===========================================================================
//  serialization (geometry, layout and B are rebuilt in setDomain)
// ===========================================================================
int KinematicCoupling::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector hdr(15);
  hdr(0) = this->getTag();
  hdr(1) = ndm;
  hdr(2) = nSlave;
  hdr(3) = Kt;
  hdr(4) = Kr;
  hdr(5) = krUser ? 1.0 : 0.0;
  hdr(6) = ktAuto ? 1.0 : 0.0;
  hdr(7) = kAlpha;
  hdr(8) = hostEleTag;
  hdr(9) = enforce;
  hdr(10) = initGapCapture ? 1.0 : 0.0;
  hdr(11) = g0Computed ? 1.0 : 0.0;
  hdr(12) = nGap;
  hdr(13) = dofSel.Size();
  hdr(14) = ell2;
  if (theChannel.sendVector(dbTag, commitTag, hdr) < 0) {
    opserr << "KinematicCoupling::sendSelf - failed to send header\n";
    return -1;
  }
  if (theChannel.sendID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "KinematicCoupling::sendSelf - failed to send node tags\n";
    return -1;
  }
  if (dofSel.Size() > 0 && theChannel.sendID(dbTag, commitTag, dofSel) < 0) {
    opserr << "KinematicCoupling::sendSelf - failed to send -dof list\n";
    return -1;
  }
  if (nGap > 0) {
    Vector payload(2 * nGap);
    for (int r = 0; r < nGap; r++) {
      payload(r) = (lambdaAL.Size() == nGap) ? lambdaAL(r) : 0.0;
      payload(nGap + r) = (g0.Size() == nGap) ? g0(r) : 0.0;
    }
    if (theChannel.sendVector(dbTag, commitTag, payload) < 0) {
      opserr << "KinematicCoupling::sendSelf - failed to send state\n";
      return -1;
    }
  }
  return 0;
}

int KinematicCoupling::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector hdr(15);
  if (theChannel.recvVector(dbTag, commitTag, hdr) < 0) {
    opserr << "KinematicCoupling::recvSelf - failed to receive header\n";
    return -1;
  }
  this->setTag((int)hdr(0));
  ndm = (int)hdr(1);
  nSlave = (int)hdr(2);
  Kt = hdr(3);
  Kr = hdr(4);
  krUser = (hdr(5) != 0.0);
  ktAuto = (hdr(6) != 0.0);
  kAlpha = hdr(7);
  hostEleTag = (int)hdr(8);
  enforce = (int)hdr(9);
  initGapCapture = (hdr(10) != 0.0);
  g0Computed = (hdr(11) != 0.0);
  nGap = (int)hdr(12);
  int nDofSel = (int)hdr(13);
  ell2 = hdr(14);
  nrot = (ndm == 3) ? 3 : 1;
  ktResolved = false;
  hostMissWarned = false;
  nDOF = 0;
  valid = false;

  connectedNodes.resize(1 + nSlave);
  if (theChannel.recvID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "KinematicCoupling::recvSelf - failed to receive node tags\n";
    return -1;
  }
  if (nDofSel > 0) {
    dofSel.resize(nDofSel);
    if (theChannel.recvID(dbTag, commitTag, dofSel) < 0) {
      opserr << "KinematicCoupling::recvSelf - failed to receive -dof list\n";
      return -1;
    }
  } else {
    dofSel = ID();
  }
  lambdaAL.resize(nGap > 0 ? nGap : 1);
  lambdaAL.Zero();
  g0.resize(nGap > 0 ? nGap : 1);
  g0.Zero();
  if (nGap > 0) {
    Vector payload(2 * nGap);
    if (theChannel.recvVector(dbTag, commitTag, payload) < 0) {
      opserr << "KinematicCoupling::recvSelf - failed to receive state\n";
      return -1;
    }
    for (int r = 0; r < nGap; r++) {
      lambdaAL(r) = payload(r);
      g0(r) = payload(nGap + r);
    }
  }

  nodeNdf.resize(1 + nSlave);
  dofOffset.resize(1 + nSlave);
  if (theNodes != 0) delete[] theNodes;
  theNodes = new Node *[1 + nSlave];
  for (int i = 0; i < 1 + nSlave; i++) theNodes[i] = 0;
  return 0;
}

void KinematicCoupling::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"KinematicCoupling\", ";
    s << "\"nodes\": [";
    for (int i = 0; i < 1 + nSlave; i++) {
      if (i > 0) s << ", ";
      s << connectedNodes(i);
    }
    s << "], \"kt\": " << Kt << ", \"kr\": " << Kr << "}";
    return;
  }
  s << "KinematicCoupling, tag: " << this->getTag() << "\n";
  s << "  reference node: " << connectedNodes(0) << "  slave nodes:";
  for (int i = 0; i < nSlave; i++) s << " " << connectedNodes(1 + i);
  s << "\n  K_t: " << Kt << "  K_r: " << Kr
    << "  enforce: " << (enforce == 1 ? "al" : "penalty")
    << "  tied components: " << nGap
    << (hasRefRot ? "" : "  (translations only)")
    << (initGapCapture ? "" : "  (absolute)") << "\n";
}

// ===========================================================================
//  responses
// ===========================================================================
Response *KinematicCoupling::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  if (argc < 1) return 0;
  static const char *const tForce[] = {"couplingForce", "couplingForces", 0};
  static const char *const tGap[] = {"gap", 0};
  static const char *const tKt[] = {"penalty", "kt", "k", 0};
  static const char *const tKr[] = {"rotPenalty", "kr", 0};
  static const char *const tLam[] = {"lambda", 0};
  static const char *const tN[] = {"tiedDOFs", "nGap", 0};
  static const char *const tE[] = {"penaltyEnergy", 0};
  static const char *const tV[] = {"constraintViolation", 0};
  int n = (nGap > 0) ? nGap : 1;
  if (tokenIs(argv[0], tForce)) return new ElementResponse(this, 1, Vector(n));
  if (tokenIs(argv[0], tGap))   return new ElementResponse(this, 2, Vector(n));
  if (tokenIs(argv[0], tKt))    return new ElementResponse(this, 3, 0.0);
  if (tokenIs(argv[0], tKr))    return new ElementResponse(this, 4, 0.0);
  if (tokenIs(argv[0], tLam))   return new ElementResponse(this, 5, Vector(n));
  if (tokenIs(argv[0], tN))     return new ElementResponse(this, 6, 0.0);
  if (tokenIs(argv[0], tE))     return new ElementResponse(this, 7, 0.0);
  if (tokenIs(argv[0], tV))     return new ElementResponse(this, 8, 0.0);
  return this->Element::setResponse(argv, argc, s);
}

int KinematicCoupling::getResponse(int responseID, Information &eleInfo)
{
  this->resolvePenalty();
  switch (responseID) {
  case 1: {                     // per-row tie force D g (+ lambda)
    Vector g(nGap);
    this->computeGap(g);
    for (int row = 0; row < nGap; row++) {
      double t = this->rowPenalty(row) * g(row);
      if (enforce == 1) t += lambdaAL(row);
      g(row) = t;
    }
    return eleInfo.setVector(g);
  }
  case 2: {
    Vector g(nGap);
    this->computeGap(g);
    return eleInfo.setVector(g);
  }
  case 3: return eleInfo.setDouble(Kt);
  case 4: return eleInfo.setDouble(Kr);
  case 5: return eleInfo.setVector(lambdaAL);
  case 6: return eleInfo.setDouble((double)nGap);
  case 7: {                     // 1/2 sum k_row g_row^2
    Vector g(nGap);
    this->computeGap(g);
    double E = 0.0;
    for (int row = 0; row < nGap; row++)
      E += this->rowPenalty(row) * g(row) * g(row);
    return eleInfo.setDouble(0.5 * E);
  }
  case 8: {                     // ||g||
    Vector g(nGap);
    this->computeGap(g);
    return eleInfo.setDouble(g.Norm());
  }
  default:
    return this->Element::getResponse(responseID, eleInfo);
  }
}
