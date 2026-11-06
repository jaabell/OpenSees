// SPDX-License-Identifier: BSD-3-Clause
/*
Copyright (c) 2026, Jose Antonio Abell. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <SplitElastic.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <ElementResponse.h>
#include <Information.h>
#include <elementAPI.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <limits>

void *OPS_SplitElastic()
{
    const int d = OPS_GetNDM(), n = d == 2 ? 4 : 8;
    if ((d != 2 && d != 3) || OPS_GetNumRemainingInputArgs() < n + 4) {
        opserr << "SplitElastic: ndm must be 2 or 3; expected tag nodes E nu rho "
               << "<-thickness t> <-planeStress|-planeStrain> <-lMass|-cMass> "
               << "<-coupling eta> <-couplingFactors f1 ... fn>\n";
        return nullptr;
    }
    int tags[9], count = n + 1;
    double prop[3];
    if (OPS_GetIntInput(&count, tags) < 0) return nullptr;
    count = 3;
    if (OPS_GetDoubleInput(&count, prop) < 0) return nullptr;
    double t = 1., eta = 0.;
    bool ps = false, lm = true;
    unsigned seen = 0;
    Vector factors(n);
    for (int i = 0; i < n; ++i) factors(i) = 1.;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        unsigned bit = 0;
        if (!strcmp(opt, "-thickness")) {
            bit = 1; count = 1;
            if (d != 2) { opserr << "SplitElastic: -thickness requires ndm=2\n"; return nullptr; }
            if (OPS_GetDoubleInput(&count, &t) < 0) { opserr << "SplitElastic: expected thickness\n"; return nullptr; }
        } else if (!strcmp(opt, "-planeStress") || !strcmp(opt, "-planeStrain")) {
            bit = 2;
            if (d != 2) { opserr << "SplitElastic: plane options require ndm=2\n"; return nullptr; }
            ps = !strcmp(opt, "-planeStress");
        } else if (!strcmp(opt, "-lMass") || !strcmp(opt, "-cMass")) {
            bit = 4; lm = !strcmp(opt, "-lMass");
        } else if (!strcmp(opt, "-coupling")) {
            bit = 8; count = 1;
            if (OPS_GetDoubleInput(&count, &eta) < 0) { opserr << "SplitElastic: expected coupling coefficient\n"; return nullptr; }
        } else if (!strcmp(opt, "-couplingFactors")) {
            bit = 16; count = n;
            if (OPS_GetDoubleInput(&count, &factors(0)) < 0) { opserr << "SplitElastic: expected one coupling factor per node\n"; return nullptr; }
        } else {
            opserr << "SplitElastic: unknown option " << opt << endln; return nullptr;
        }
        if (seen & bit) { opserr << "SplitElastic: duplicate/conflicting option " << opt << endln; return nullptr; }
        seen |= bit;
    }
    ID ids(n);
    for (int i = 0; i < n; ++i) ids(i) = tags[i + 1];
    auto *element = new SplitElastic(tags[0], d, ids, prop[0], prop[1], prop[2], t, ps, lm, eta, factors);
    element->setDomain(OPS_GetDomain());
    if (!element->isValid()) { delete element; return nullptr; }
    return element;
}

SplitElastic::SplitElastic() : Element(0, ELE_TAG_SplitElastic) {}

SplitElastic::SplitElastic(int tag, int d, const ID &ids, double e, double v,
    double density, double t, bool ps, bool lm, double eta, const Vector &f)
    : Element(tag, ELE_TAG_SplitElastic), ndm(d), nodes(ids), E(e), nu(v),
      rho(density), thickness(t), coupling(eta), planeStress(ps), lumped(lm), factors(f) {}

bool SplitElastic::propertiesValid() const
{
    if ((ndm != 2 && ndm != 3) || nodes.Size() != (ndm == 2 ? 4 : 8) ||
        factors.Size() != nodes.Size()) return false;
    if (!std::isfinite(E) || E <= 0. || !std::isfinite(nu) || nu <= -1. || nu >= .5 ||
        !std::isfinite(rho) || rho < 0. || !std::isfinite(thickness) || thickness <= 0. ||
        !std::isfinite(coupling) || coupling < 0. || (ndm == 3 && planeStress)) return false;
    for (int i = 0; i < nodes.Size(); ++i) {
        if (!std::isfinite(factors(i)) || factors(i) < 0.) return false;
        for (int j = 0; j < i; ++j) if (nodes(i) == nodes(j)) return false;
    }
    return true;
}

void SplitElastic::setDomain(Domain *domain)
{
    valid = false;
    nodePtrs.fill(nullptr);
    if (!domain) { Element::setDomain(nullptr); return; }
    if (!propertiesValid()) {
        opserr << "SplitElastic " << getTag() << ": invalid properties/connectivity\n"; return;
    }
    int offset = 0;
    allShared = true;
    for (int i = 0; i < nodes.Size(); ++i) {
        Node *node = domain->getNode(nodes(i));
        if (!node || node->getCrds().Size() != ndm ||
            (node->getNumberDOF() != ndm && node->getNumberDOF() != 2 * ndm)) {
            opserr << "SplitElastic " << getTag() << ": node " << nodes(i)
                   << " must have ndm coordinates and ndm or 2*ndm DOFs\n"; return;
        }
        nodePtrs[i] = node;
        const bool paired = node->getNumberDOF() == 2 * ndm;
        allShared = allShared && !paired;
        for (int c = 0; c < ndm; ++c) {
            mapA[i * ndm + c] = offset + c;
            mapB[i * ndm + c] = offset + c + (paired ? ndm : 0);
        }
        offset += node->getNumberDOF();
    }
    // Rebinding after deserialization must preserve applied element loads.
    if (numDOF != offset) { load.resize(offset); load.Zero(); }
    numDOF = offset;
    resisting.resize(numDOF);
    stiffness.resize(nodes.Size() * ndm, nodes.Size() * ndm);
    scalarMass.resize(nodes.Size(), nodes.Size());
    volumes.resize(nodes.Size());
    valid = integrate();
    if (!valid) { opserr << "SplitElastic " << getTag() << ": nonpositive/degenerate Jacobian or nonfinite operator\n"; return; }
    Element::setDomain(domain);
}

void SplitElastic::constitutive(double D[6][6]) const
{
    const double mu = E / (2. * (1. + nu));
    const double lambda = planeStress ? E * nu / (1. - nu * nu) : E * nu / ((1. + nu) * (1. - 2. * nu));
    const int ns = ndm == 2 ? 3 : 6;
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) D[i][j] = 0.;
    for (int i = 0; i < ndm; ++i) for (int j = 0; j < ndm; ++j)
        D[i][j] = lambda + (i == j ? 2. * mu : 0.);
    for (int i = ndm; i < ns; ++i) D[i][i] = mu;
}

bool SplitElastic::shape(int gp, double N[8], double grad[8][3], double &weight) const
{
    static const int signs[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
                                   {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    const double q = 1. / std::sqrt(3.);
    double xi[3] = {signs[gp][0] * q, signs[gp][1] * q, signs[gp][2] * q};
    double natural[8][3] = {}, J[3][3] = {}, inv[3][3] = {};
    for (int i = 0; i < nodes.Size(); ++i) {
        N[i] = ndm == 2 ? .25 : .125;
        for (int c = 0; c < ndm; ++c) N[i] *= 1. + signs[i][c] * xi[c];
        for (int c = 0; c < ndm; ++c) {
            natural[i][c] = (ndm == 2 ? .25 : .125) * signs[i][c];
            for (int k = 0; k < ndm; ++k) if (k != c) natural[i][c] *= 1. + signs[i][k] * xi[k];
        }
        const Vector &x = nodePtrs[i]->getCrds();
        // Subtract a reference point to avoid translation-dependent cancellation.
        const Vector &origin = nodePtrs[0]->getCrds();
        for (int r = 0; r < ndm; ++r) for (int c = 0; c < ndm; ++c)
            J[r][c] += (x(r) - origin(r)) * natural[i][c];
    }
    double det;
    if (ndm == 2) {
        det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        inv[0][0] = J[1][1]; inv[0][1] = -J[0][1];
        inv[1][0] = -J[1][0]; inv[1][1] = J[0][0];
    } else {
        inv[0][0] = J[1][1]*J[2][2]-J[1][2]*J[2][1];
        inv[0][1] = J[0][2]*J[2][1]-J[0][1]*J[2][2];
        inv[0][2] = J[0][1]*J[1][2]-J[0][2]*J[1][1];
        inv[1][0] = J[1][2]*J[2][0]-J[1][0]*J[2][2];
        inv[1][1] = J[0][0]*J[2][2]-J[0][2]*J[2][0];
        inv[1][2] = J[0][2]*J[1][0]-J[0][0]*J[1][2];
        inv[2][0] = J[1][0]*J[2][1]-J[1][1]*J[2][0];
        inv[2][1] = J[0][1]*J[2][0]-J[0][0]*J[2][1];
        inv[2][2] = J[0][0]*J[1][1]-J[0][1]*J[1][0];
        det = J[0][0]*inv[0][0]+J[0][1]*inv[1][0]+J[0][2]*inv[2][0];
    }
    double scale = 0.;
    for (int r = 0; r < ndm; ++r) for (int c = 0; c < ndm; ++c) scale = std::max(scale, std::abs(J[r][c]));
    if (!std::isfinite(det) || det <= 64. * std::numeric_limits<double>::epsilon() * std::pow(scale, ndm)) return false;
    for (int i = 0; i < nodes.Size(); ++i) for (int c = 0; c < ndm; ++c) {
        grad[i][c] = 0.;
        for (int k = 0; k < ndm; ++k) grad[i][c] += natural[i][k] * inv[k][c] / det;
    }
    weight = det * (ndm == 2 ? thickness : 1.);
    return std::isfinite(weight);
}

void SplitElastic::strainMatrix(const double grad[8][3], double B[6][24]) const
{
    for (int r = 0; r < 6; ++r) for (int c = 0; c < 24; ++c) B[r][c] = 0.;
    for (int i = 0; i < nodes.Size(); ++i) {
        int a = i * ndm;
        for (int c = 0; c < ndm; ++c) B[c][a+c] = grad[i][c];
        int xy = ndm == 2 ? 2 : 3;
        B[xy][a] = grad[i][1]; B[xy][a+1] = grad[i][0];
        if (ndm == 3) {
            B[4][a+1] = grad[i][2]; B[4][a+2] = grad[i][1];
            B[5][a] = grad[i][2]; B[5][a+2] = grad[i][0];
        }
    }
}

bool SplitElastic::integrate()
{
    stiffness.Zero(); scalarMass.Zero(); volumes.Zero();
    double D[6][6]; constitutive(D);
    const int n = nodes.Size(), p = n * ndm, ns = ndm == 2 ? 3 : 6;
    for (int g = 0; g < n; ++g) {
        double N[8], grad[8][3], w, B[6][24], DB[6][24] = {};
        if (!shape(g, N, grad, w)) return false;
        strainMatrix(grad, B);
        for (int r = 0; r < ns; ++r) for (int j = 0; j < p; ++j)
            for (int c = 0; c < ns; ++c) DB[r][j] += D[r][c] * B[c][j];
        for (int i = 0; i < p; ++i) for (int j = 0; j < p; ++j)
            for (int r = 0; r < ns; ++r) stiffness(i,j) += B[r][i] * DB[r][j] * w;
        for (int i = 0; i < n; ++i) {
            volumes(i) += N[i] * w;
            for (int j = 0; j < n; ++j) scalarMass(i,j) += rho * N[i] * N[j] * w;
        }
    }
    if (lumped) {
        scalarMass.Zero();
        for (int i = 0; i < n; ++i) scalarMass(i,i) = rho * volumes(i);
    }
    for (int i = 0; i < p; ++i) for (int j = 0; j < p; ++j)
        if (!std::isfinite(stiffness(i,j))) return false;
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(coupling * factors(i) * volumes(i))) return false;
        for (int j = 0; j < n; ++j) if (!std::isfinite(scalarMass(i,j))) return false;
    }
    return true;
}

void SplitElastic::assemble(Matrix &out, double kScale, double mScale, bool relative) const
{
    out.resize(numDOF, numDOF); out.Zero();
    if (!valid) return;
    const int p = nodes.Size() * ndm;
    for (int i = 0; i < p; ++i) for (int j = 0; j < p; ++j) {
        const double value = .5 * (kScale * stiffness(i,j) +
            (i % ndm == j % ndm ? mScale * scalarMass(i/ndm,j/ndm) : 0.));
        out(mapA[i],mapA[j]) += value;
        out(mapB[i],mapB[j]) += value;
    }
    if (relative) for (int i = 0; i < nodes.Size(); ++i) {
        const double c = coupling * factors(i) * volumes(i);
        for (int k = 0; k < ndm; ++k) {
            const int a = mapA[i*ndm+k], b = mapB[i*ndm+k];
            if (a == b) continue;
            out(a,a) += c; out(b,b) += c; out(a,b) -= c; out(b,a) -= c;
        }
    }
}

const Matrix &SplitElastic::getTangentStiff()
{
    static thread_local Matrix result;
    assemble(result, 1., 0., false); return result;
}
const Matrix &SplitElastic::getMass()
{
    static thread_local Matrix result;
    assemble(result, 0., 1., false); return result;
}
const Matrix &SplitElastic::getDamp()
{
    static thread_local Matrix result;
    assemble(result, betaK + betaK0 + betaKc, alphaM, true); return result;
}
int SplitElastic::setRayleighDampingFactors(double a, double b, double b0, double bc)
{
    // The tangent is constant: no committed-stiffness copy or reintegration.
    alphaM = a; betaK = b; betaK0 = b0; betaKc = bc;
    return 0;
}

void SplitElastic::gather(int branch, int quantity, double *out) const
{
    int k = 0;
    for (int i = 0; i < nodes.Size(); ++i) {
        const Vector &v = quantity == 0 ? nodePtrs[i]->getTrialDisp() :
                          quantity == 1 ? nodePtrs[i]->getTrialVel() : nodePtrs[i]->getTrialAccel();
        const int start = branch && nodePtrs[i]->getNumberDOF() == 2*ndm ? ndm : 0;
        for (int c = 0; c < ndm; ++c) out[k++] = v(start+c);
    }
}
void SplitElastic::scatter(int branch, const double *in, double scale, Vector &out) const
{
    const auto &map = branch ? mapB : mapA;
    for (int i = 0; i < nodes.Size()*ndm; ++i) out(map[i]) += scale * in[i];
}
void SplitElastic::massAction(const double *in, double *out) const
{
    const int n = nodes.Size();
    std::fill(out, out + n*ndm, 0.);
    if (lumped) {
        for (int i = 0; i < n; ++i) for (int c = 0; c < ndm; ++c) out[i*ndm+c] = scalarMass(i,i) * in[i*ndm+c];
    } else {
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
            for (int c = 0; c < ndm; ++c) out[i*ndm+c] += scalarMass(i,j) * in[j*ndm+c];
    }
}
void SplitElastic::force(bool inertia)
{
    resisting.Zero();
    if (!valid) return;
    const int p = nodes.Size()*ndm;
    const double beta = betaK + betaK0 + betaKc;
    for (int b = 0; b < (allShared ? 1 : 2); ++b) {
        double u[24], v[24] = {}, a[24] = {}, f[24] = {}, ma[24];
        gather(b, 0, u);
        if (inertia) {
            gather(b, 1, v); gather(b, 2, a);
            for (int i = 0; i < p; ++i) { u[i] += beta * v[i]; a[i] += alphaM * v[i]; }
            massAction(a, ma);
            for (int i = 0; i < p; ++i) f[i] = ma[i];
        }
        // Matrix is column-major; visit contiguous rows in the inner loop.
        for (int j = 0; j < p; ++j) for (int i = 0; i < p; ++i) f[i] += stiffness(i,j) * u[j];
        scatter(b, f, allShared ? 1. : .5, resisting);
    }
    if (inertia && coupling != 0. && !allShared) {
        double va[24], vb[24]; gather(0,1,va); gather(1,1,vb);
        for (int i = 0; i < nodes.Size(); ++i) for (int k = 0; k < ndm; ++k) {
            int j = i*ndm+k;
            if (mapA[j] == mapB[j]) continue;
            double f = coupling * factors(i) * volumes(i) * (va[j]-vb[j]);
            resisting(mapA[j]) += f; resisting(mapB[j]) -= f;
        }
    }
    resisting -= load;
}
const Vector &SplitElastic::getResistingForce() { force(false); return resisting; }
const Vector &SplitElastic::getResistingForceIncInertia() { force(true); return resisting; }

int SplitElastic::addLoad(ElementalLoad *, double)
{
    opserr << "SplitElastic " << getTag() << ": element loads unsupported; use nodal loads\n";
    return -1;
}
int SplitElastic::addInertiaLoadToUnbalance(const Vector &accel)
{
    if (!valid) return -1;
    double rv[2][24] = {}, f[24];
    for (int i = 0; i < nodes.Size(); ++i) {
        const Vector &v = nodePtrs[i]->getRV(accel);
        if (v.Size() != nodePtrs[i]->getNumberDOF()) return -1;
        for (int c = 0; c < ndm; ++c) {
            rv[0][i*ndm+c] = v(c);
            // Ground acceleration is a physical translation, common to both
            // branches. UniformExcitation only fills the first ndm R rows;
            // the additional node DOFs are not rotations in this element.
            rv[1][i*ndm+c] = v(c);
        }
    }
    for (int b = 0; b < 2; ++b) { massAction(rv[b], f); scatter(b, f, -.5, load); }
    return 0;
}

double SplitElastic::relativePower() const
{
    if (!valid || coupling == 0. || allShared) return 0.;
    double va[24], vb[24], power = 0.; gather(0,1,va); gather(1,1,vb);
    for (int i = 0; i < nodes.Size(); ++i) for (int c = 0; c < ndm; ++c) {
        double dv = va[i*ndm+c] - vb[i*ndm+c];
        power += coupling * factors(i) * volumes(i) * dv * dv;
    }
    return power;
}

int SplitElastic::sendSelf(int commitTag, Channel &channel)
{
    int db = getDbTag();
    if (!db) { db = channel.getDbTag(); setDbTag(db); }
    ID header(16); header.Zero();
    header(0) = 1; header(1) = getTag(); header(2) = ndm; header(3) = planeStress;
    header(4) = lumped; header(5) = numDOF; header(6) = nodes.Size();
    for (int i = 0; i < nodes.Size(); ++i) header(8+i) = nodes(i);
    Vector data(9 + factors.Size() + load.Size());
    data(0)=E; data(1)=nu; data(2)=rho; data(3)=thickness; data(4)=coupling;
    data(5)=alphaM; data(6)=betaK; data(7)=betaK0; data(8)=betaKc;
    for (int i = 0; i < factors.Size(); ++i) data(9+i) = factors(i);
    for (int i = 0; i < load.Size(); ++i) data(9+factors.Size()+i) = load(i);
    if (channel.sendID(db, commitTag, header) < 0 || channel.sendVector(db, commitTag, data) < 0) return -1;
    return 0;
}
int SplitElastic::recvSelf(int commitTag, Channel &channel, FEM_ObjectBroker &)
{
    ID header(16);
    if (channel.recvID(getDbTag(), commitTag, header) < 0) return -1;
    const int d = header(2), n = header(6), dof = header(5);
    if (header(0) != 1 || (d != 2 && d != 3) || n != (d == 2 ? 4 : 8) ||
        dof < n*d || dof > 2*n*d || (header(3) != 0 && header(3) != 1) ||
        (header(4) != 0 && header(4) != 1)) return -1;
    Vector data(9+n+dof);
    if (channel.recvVector(getDbTag(), commitTag, data) < 0) return -1;
    Element::setDomain(nullptr); nodePtrs.fill(nullptr); valid = false;
    setTag(header(1)); ndm = d; planeStress = header(3); lumped = header(4); numDOF = dof;
    nodes.resize(n); factors.resize(n); load.resize(dof); resisting.resize(dof);
    for (int i = 0; i < n; ++i) { nodes(i) = header(8+i); factors(i) = data(9+i); }
    E=data(0); nu=data(1); rho=data(2); thickness=data(3); coupling=data(4);
    alphaM=data(5); betaK=data(6); betaK0=data(7); betaKc=data(8);
    for (int i = 0; i < dof; ++i) load(i) = data(9+n+i);
    return propertiesValid() ? 0 : -1;
}

void SplitElastic::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "{\"name\": " << getTag() << ", \"type\": \"SplitElastic\", \"nodes\": [";
        for (int i = 0; i < nodes.Size(); ++i) { if (i) s << ", "; s << nodes(i); }
        s << "], \"E\": " << E << ", \"nu\": " << nu << ", \"rho\": " << rho
          << ", \"ndm\": " << ndm << ", \"thickness\": " << thickness
          << ", \"planeStress\": " << (planeStress ? "true" : "false")
          << ", \"coupling\": " << coupling << ", \"lumped\": " << (lumped ? "true" : "false")
          << ", \"couplingFactors\": [";
        for (int i = 0; i < factors.Size(); ++i) { if (i) s << ", "; s << factors(i); }
        s << "]}";
    } else {
        s << "SplitElastic " << getTag() << " ndm=" << ndm << " nodes=" << nodes
          << " E=" << E << " nu=" << nu << " rho=" << rho << " thickness=" << thickness
          << " coupling=" << coupling << " lumped=" << lumped << endln;
    }
}

Response *SplitElastic::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc < 1) return nullptr;
    int id = 0, size = 0;
    if (!strcmp(argv[0],"forces") || !strcmp(argv[0],"force")) { id=1; size=numDOF; }
    else if (!strcmp(argv[0],"relativePower")) { id=2; size=1; }
    else if (!strcmp(argv[0],"energy")) { id=3; size=3; }
    else if (!strcmp(argv[0],"stresses")) { id=4; }
    else if (!strcmp(argv[0],"stresses1")) { id=5; }
    else if (!strcmp(argv[0],"stresses2")) { id=6; }
    else if (!strcmp(argv[0],"strains")) { id=7; }
    else if (!strcmp(argv[0],"strains1")) { id=8; }
    else if (!strcmp(argv[0],"strains2")) { id=9; }
    else if (!strcmp(argv[0],"stiffness")) { id=11; size=numDOF*numDOF; }
    else if (!strcmp(argv[0],"mass")) { id=12; size=numDOF*numDOF; }
    else if (!strcmp(argv[0],"damping")) { id=13; size=numDOF*numDOF; }
    else return Element::setResponse(argv, argc, output);
    if (!size) size = nodes.Size() * (ndm == 2 ? 3 : 6);
    output.tag("ElementOutput"); output.attr("eleType","SplitElastic"); output.attr("eleTag",getTag());
    // Expanded node DOFs are branch translations, never rotational quantities.
    const char *components[6] = {"xx", "yy", "zz", "xy", "yz", "zx"};
    const int ns = ndm == 2 ? 3 : 6;
    for (int i = 0; i < size; ++i) {
        char label[80];
        if (id == 3) {
            const char *names[3] = {"strainEnergy", "kineticEnergy", "relativePower"};
            std::snprintf(label, sizeof(label), "%s", names[i]);
        } else if (id >= 4 && id <= 9) {
            int c = i % ns;
            const char *component = ndm == 2 && c == 2 ? "xy" : components[c];
            std::snprintf(label, sizeof(label), "%s_gp%d_%s", argv[0], i/ns+1, component);
        } else {
            std::snprintf(label, sizeof(label), "%s_%d", argv[0], i+1);
        }
        output.tag("ResponseType", label);
    }
    auto *response = new ElementResponse(this, id, Vector(size));
    output.endTag(); return response;
}
int SplitElastic::getResponse(int id, Information &info)
{
    if (!valid) return -1;
    if (id == 1) return info.setVector(getResistingForce());
    if (id == 2) { Vector value(1); value(0)=relativePower(); return info.setVector(value); }
    if (id == 3) {
        Vector energy(3); energy.Zero();
        const int p = nodes.Size()*ndm;
        for (int b = 0; b < 2; ++b) {
            double u[24], v[24], mv[24]; gather(b,0,u); gather(b,1,v); massAction(v,mv);
            for (int i = 0; i < p; ++i) {
                energy(1) += .25 * v[i] * mv[i];
                for (int j = 0; j < p; ++j) energy(0) += .25 * u[i] * stiffness(i,j) * u[j];
            }
        }
        energy(2)=relativePower(); return info.setVector(energy);
    }
    if (id >= 11 && id <= 13) {
        const Matrix &mat = id == 11 ? getTangentStiff() : id == 12 ? getMass() : getDamp();
        Vector values(numDOF*numDOF);
        for (int i = 0; i < numDOF; ++i) for (int j = 0; j < numDOF; ++j) values(i*numDOF+j)=mat(i,j);
        return info.setVector(values);
    }
    if (id >= 4 && id <= 9) {
        const int ns = ndm == 2 ? 3 : 6, p = nodes.Size()*ndm;
        const bool stress = id <= 6;
        const int which = stress ? id - 4 : id - 7; // 0 common, 1 branch1, 2 branch2
        Vector values(nodes.Size()*ns); values.Zero();
        double D[6][6]; constitutive(D);
        for (int b = 0; b < 2; ++b) {
            if (which && which != b+1) continue;
            double u[24]; gather(b,0,u);
            for (int gp = 0; gp < nodes.Size(); ++gp) {
                double N[8], grad[8][3], w, B[6][24], strain[6] = {};
                if (!shape(gp,N,grad,w)) return -1;
                strainMatrix(grad,B);
                for (int r = 0; r < ns; ++r) for (int j = 0; j < p; ++j) strain[r] += B[r][j]*u[j];
                for (int r = 0; r < ns; ++r) {
                    if (stress) for (int c = 0; c < ns; ++c) values(gp*ns+r) += .5*D[r][c]*strain[c];
                    else values(gp*ns+r) += (which ? 1. : .5)*strain[r];
                }
            }
        }
        return info.setVector(values);
    }
    return -1;
}
