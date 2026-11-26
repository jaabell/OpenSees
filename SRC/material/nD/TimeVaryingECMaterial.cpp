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

// Written: J.A. Abell (Universidad de los Andes), with I. Calderon.
// See TimeVaryingECMaterial.h for the formulation, scope and limits.

#include <TimeVaryingECMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <Domain.h>
#include <ID.h>
#include <elementAPI.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <string.h>
#include <math.h>

Vector TimeVaryingECMaterial::workVector(6);
Matrix TimeVaryingECMaterial::workMatrix(6, 6);

// sendSelf layout version
static const int TVEC_VERSION = 1;
static const int TVEC_STATE_SIZE = 6 * 7 + 4 + 36;   // see packState

static bool allFinite(const Vector &v)
{
    for (int i = 0; i < v.Size(); i++)
        if (!std::isfinite(v(i)))
            return false;
    return true;
}

static bool allFinite(const Matrix &m)
{
    for (int i = 0; i < m.noRows(); i++)
        for (int j = 0; j < m.noCols(); j++)
            if (!std::isfinite(m(i, j)))
                return false;
    return true;
}

void *OPS_TimeVaryingECMaterial(void)
{
    const char *usage =
        "nDMaterial TimeVaryingEC $tag $refMatTag $N t1 .. tN E1 .. EN K1 .. KN A1 .. AN\n";

    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "nDMaterial TimeVaryingEC Error: too few arguments.\n" << usage;
        return 0;
    }
    int iData[2];
    int numData = 2;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "nDMaterial TimeVaryingEC Error: invalid nDMaterial tags.\n";
        return 0;
    }
    int N = 0;
    numData = 1;
    if (OPS_GetInt(&numData, &N) != 0 || N < 1) {
        opserr << "nDMaterial TimeVaryingEC Error: invalid number of data points N\n";
        return 0;
    }
    if (OPS_GetNumRemainingInputArgs() < 4 * N) {
        opserr << "nDMaterial TimeVaryingEC Error: expected " << 4 * N
               << " values (t, E, K, A) for tag " << iData[0] << ".\n" << usage;
        return 0;
    }
    Vector t(N), E(N), K(N), A(N);
    Vector *h[4] = {&t, &E, &K, &A};
    for (int k = 0; k < 4; k++) {
        numData = N;
        if (OPS_GetDouble(&numData, &(*h[k])(0)) != 0) {
            opserr << "nDMaterial TimeVaryingEC Error: reading the property law (tag "
                   << iData[0] << ").\n";
            return 0;
        }
    }
    for (int k = 0; k < N; k++) {
        if (!std::isfinite(t(k)) || !std::isfinite(E(k)) || !std::isfinite(K(k)) ||
            !std::isfinite(A(k)) || (k > 0 && t(k) <= t(k - 1))) {
            opserr << "nDMaterial TimeVaryingEC Error: driver values must be finite and "
                      "strictly increasing (tag " << iData[0] << ").\n";
            return 0;
        }
        if (A(k) <= 0.0 || K(k) <= 0.0 || E(k) <= 0.0 || E(k) >= 9.0 * K(k)) {
            opserr << "nDMaterial TimeVaryingEC Error: need A > 0, K > 0 and 0 < E < 9K "
                      "(so that G > 0) at every point (tag " << iData[0] << ", point "
                   << k + 1 << ").\n";
            return 0;
        }
    }
    NDMaterial *ref = OPS_getNDMaterial(iData[1]);
    if (ref == 0) {
        opserr << "nDMaterial TimeVaryingEC Error: no nDMaterial with tag " << iData[1] << ".\n";
        return 0;
    }
    TimeVaryingECMaterial *mat = new TimeVaryingECMaterial(iData[0], *ref, t, E, K, A);
    if (mat->getOrder() != 6) {   // constructor failed to get a usable reference
        delete mat;
        return 0;
    }
    return mat;
}

TimeVaryingECMaterial::State::State()
    : eps(6), epsStar(6), epsM(6), driver(0.0), K(0.0), G(0.0), A(1.0),
      epsP(6), sigP(6), sig(6), Calg(6, 6), eta(6)
{
}

TimeVaryingECMaterial::TimeVaryingECMaterial(int tag, NDMaterial &ref,
                                             const Vector &t, const Vector &E,
                                             const Vector &K, const Vector &A)
    : NDMaterial(tag, ND_TAG_TimeVaryingECMaterial),
      theReference(0), C0(6, 6), C0inv(6, 6),
      trialSinceCommit(false), referenceInSync(true), pristine(true),
      inputEpsStar(6), committedEpsStar(6),
      inputDriver(0.0), committedDriver(0.0),
      useLocalDriver(false), committedUseLocalDriver(false),
      referenceDbTag(0)
{
    std::shared_ptr<Law> l = std::make_shared<Law>();
    l->t = t; l->E = E; l->K = K; l->A = A;
    law = l;

    theReference = ref.getCopy("ThreeDimensional");
    if (theReference == 0 || theReference->getOrder() != 6) {
        opserr << "nDMaterial TimeVaryingEC Error: the reference material (tag "
               << ref.getTag() << ") has no 3D copy\n";
        if (theReference != 0) delete theReference;
        theReference = 0;
        return;
    }
    // reference elastic tensor of the virgin reference material
    C0 = theReference->getInitialTangent();
    if (!allFinite(C0) || C0.Invert(C0inv) < 0) {
        opserr << "nDMaterial TimeVaryingEC Error: singular reference initial tangent\n";
        delete theReference;
        theReference = 0;
        return;
    }
    resetToStart();
}

TimeVaryingECMaterial::TimeVaryingECMaterial()
    : NDMaterial(0, ND_TAG_TimeVaryingECMaterial),
      theReference(0), C0(6, 6), C0inv(6, 6),
      trialSinceCommit(false), referenceInSync(true), pristine(true),
      inputEpsStar(6), committedEpsStar(6),
      inputDriver(0.0), committedDriver(0.0),
      useLocalDriver(false), committedUseLocalDriver(false),
      referenceDbTag(0)
{
}

TimeVaryingECMaterial::~TimeVaryingECMaterial()
{
    if (theReference != 0)
        delete theReference;
}

double TimeVaryingECMaterial::getRho(void)
{
    return theReference != 0 ? theReference->getRho() : 0.0;
}

double TimeVaryingECMaterial::currentDriver(void) const
{
    if (useLocalDriver)
        return inputDriver;
    // the domain being updated: in OpenSeesSP a remote partition is an
    // ActorSubdomain whose time is not that of the process's global domain
    // (same choice as the TDConcrete family)
    Domain *dom = ops_TheActiveDomain != 0 ? ops_TheActiveDomain : OPS_GetDomain();
    return dom != 0 ? dom->getCurrentTime() : 0.0;
}

// same interpolation as TimeVaryingMaterial::getParameters
int TimeVaryingECMaterial::evaluateLaw(double x, double &K, double &G, double &A) const
{
    if (!law || law->t.Size() < 1 || !std::isfinite(x))
        return -1;
    const Law &L = *law;
    int n = L.t.Size();
    int index = 0;
    for (; index < n; ++index)
        if (L.t(index) >= x)
            break;
    double E;
    if (index == 0) {
        E = L.E(0); A = L.A(0); K = L.K(0);
    } else if (index < n) {
        double t1 = L.t(index - 1), t2 = L.t(index);
        double a = (x - t1) / (t2 - t1);
        E = (1.0 - a) * L.E(index - 1) + a * L.E(index);
        A = (1.0 - a) * L.A(index - 1) + a * L.A(index);
        K = (1.0 - a) * L.K(index - 1) + a * L.K(index);
    } else {
        E = L.E(n - 1); A = L.A(n - 1); K = L.K(n - 1);
    }
    G = (3.0 * K * E) / (9.0 * K - E);
    if (!(K > 0.0) || !(G > 0.0) || !(A > 0.0) || !std::isfinite(G))
        return -1;
    return 0;
}

// isotropic C in Voigt, engineering shear strains
void TimeVaryingECMaterial::physicalStiffness(double K, double G, Matrix &C) const
{
    C.Zero();
    double l = K - 2.0 * G / 3.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            C(i, j) = l;
        C(i, i) = l + 2.0 * G;
        C(i + 3, i + 3) = G;
    }
}

void TimeVaryingECMaterial::mapS(double K, double G, double A, Matrix &S) const
{
    static Matrix C(6, 6);
    physicalStiffness(K, G, C);
    S.addMatrixProduct(0.0, C0inv, C, A);
}

// eta = C^-1 sig for the isotropic C of the state (closed form)
void TimeVaryingECMaterial::recoverableStrain(State &s) const
{
    double p = (s.sig(0) + s.sig(1) + s.sig(2)) / 3.0;
    for (int i = 0; i < 3; i++) {
        s.eta(i) = p / (3.0 * s.K) + (s.sig(i) - p) / (2.0 * s.G);
        s.eta(i + 3) = s.sig(i + 3) / s.G;
    }
}

void TimeVaryingECMaterial::resetToStart(void)
{
    inputEpsStar.Zero();
    committedEpsStar.Zero();
    inputDriver = 0.0;
    committedDriver = 0.0;
    useLocalDriver = false;
    committedUseLocalDriver = false;

    State s;
    s.driver = currentDriver();
    if (evaluateLaw(s.driver, s.K, s.G, s.A) != 0) {
        s.K = C0(0, 0) - 4.0 * C0(3, 3) / 3.0;
        s.G = C0(3, 3);
        s.A = 1.0;
    }
    physicalStiffness(s.K, s.G, s.Calg);
    committed = s;
    trial = s;
    trialSinceCommit = false;
    referenceInSync = true;
    pristine = true;
}

int TimeVaryingECMaterial::setTrialStrain(const Vector &strain, const Vector &rate)
{
    return this->setTrialStrain(strain);
}

int TimeVaryingECMaterial::setTrialStrain(const Vector &strain)
{
    if (theReference == 0)
        return -1;
    if (strain.Size() != 6 || !allFinite(strain)) {
        opserr << "nDMaterial TimeVaryingEC (tag " << this->getTag()
               << ") Error: invalid trial strain; state unchanged\n";
        return -1;
    }

    State s;
    s.eps = strain;
    s.epsStar = inputEpsStar;
    s.epsM = strain;
    s.epsM.addVector(1.0, inputEpsStar, -1.0);
    s.driver = currentDriver();
    if (evaluateLaw(s.driver, s.K, s.G, s.A) != 0) {
        opserr << "nDMaterial TimeVaryingEC (tag " << this->getTag()
               << ") Error: invalid properties at driver " << s.driver << "; state unchanged\n";
        return -1;
    }

    const State &n = committed;
    static Matrix S(6, 6), Sn(6, 6), dS(6, 6), C(6, 6), Cn(6, 6);
    static Vector dEpsM(6), dSeta(6);
    mapS(s.K, s.G, s.A, S);
    mapS(n.K, n.G, n.A, Sn);
    dS = S;
    dS.addMatrix(1.0, Sn, -1.0);
    dEpsM = s.epsM;
    dEpsM.addVector(1.0, n.epsM, -1.0);
    dSeta.addMatrixVector(0.0, dS, n.eta, 1.0);           // (S - S_n) eta_n

    s.epsP = n.epsP;
    s.epsP.addMatrixVector(1.0, S, dEpsM, 1.0);
    s.epsP.addVector(1.0, dSeta, 1.0);

    int res = theReference->setTrialStrain(s.epsP);
    if (res == 0) {
        s.sigP = theReference->getStress();
        static Matrix Ct(6, 6);
        Ct = theReference->getTangent();
        // sig = sig_n + [sigP - sigP_n - C0 (S - S_n) eta_n] / A + (C - C_n) eta_n
        static Vector work(6);
        work = s.sigP;
        work.addVector(1.0, n.sigP, -1.0);
        work.addMatrixVector(1.0, C0, dSeta, -1.0);
        s.sig = n.sig;
        s.sig.addVector(1.0, work, 1.0 / s.A);
        physicalStiffness(s.K, s.G, C);
        physicalStiffness(n.K, n.G, Cn);
        C.addMatrix(1.0, Cn, -1.0);
        s.sig.addMatrixVector(1.0, C, n.eta, 1.0);
        s.Calg.addMatrixProduct(0.0, Ct, S, 1.0 / s.A);
        if (!allFinite(s.sigP) || !allFinite(s.sig) || !allFinite(s.Calg))
            res = -2;
    }
    if (res != 0) {
        // the reference trial no longer matches the kept trial; it is redone
        // from the reference's committed state before a commit
        referenceInSync = false;
        opserr << "nDMaterial TimeVaryingEC (tag " << this->getTag()
               << ") Error: the reference material failed; previous trial kept\n";
        return res;
    }
    trial = s;
    trialSinceCommit = true;
    referenceInSync = true;
    return 0;
}

const Vector &TimeVaryingECMaterial::getStrain(void)
{
    return trial.eps;
}

const Vector &TimeVaryingECMaterial::getStress(void)
{
    return trial.sig;
}

const Matrix &TimeVaryingECMaterial::getTangent(void)
{
    if (pristine && !trialSinceCommit)
        return this->getInitialTangent();
    return trial.Calg;
}

const Matrix &TimeVaryingECMaterial::getInitialTangent(void)
{
    double K = trial.K, G = trial.G, A = trial.A;
    if (pristine && !trialSinceCommit) {
        // before any loading (or after revertToStart): evaluate the law at the
        // current driver without changing the state
        double Kd, Gd, Ad;
        if (evaluateLaw(currentDriver(), Kd, Gd, Ad) == 0) {
            K = Kd; G = Gd; A = Ad;
        }
    }
    physicalStiffness(K, G, workMatrix);
    return workMatrix;
}

int TimeVaryingECMaterial::commitState(void)
{
    if (theReference == 0)
        return -1;
    if (!trialSinceCommit) {
        // nothing new: the reference is not committed again, so a stale
        // reference trial (after a failed trial or a revert) cannot leak in
        committedEpsStar = inputEpsStar;
        committedDriver = inputDriver;
        committedUseLocalDriver = useLocalDriver;
        return 0;
    }
    if (!referenceInSync) {
        int res = theReference->setTrialStrain(trial.epsP);
        if (res != 0) {
            opserr << "nDMaterial TimeVaryingEC (tag " << this->getTag()
                   << ") Error: the reference trial could not be restored; not committed\n";
            return res;
        }
        referenceInSync = true;
    }
    int res = theReference->commitState();
    if (res != 0)
        return res;
    committed = trial;
    recoverableStrain(committed);
    trial = committed;
    committedEpsStar = inputEpsStar;
    committedDriver = inputDriver;
    committedUseLocalDriver = useLocalDriver;
    trialSinceCommit = false;
    pristine = false;
    return 0;
}

int TimeVaryingECMaterial::revertToLastCommit(void)
{
    if (theReference == 0)
        return -1;
    // the reference is told, but its own revert is not relied upon: every
    // query is answered from `committed`, and its next trial integrates from
    // its committed state
    theReference->revertToLastCommit();
    trial = committed;
    inputEpsStar = committedEpsStar;
    inputDriver = committedDriver;
    useLocalDriver = committedUseLocalDriver;
    trialSinceCommit = false;
    referenceInSync = false;
    return 0;
}

int TimeVaryingECMaterial::revertToStart(void)
{
    if (theReference == 0)
        return -1;
    int res = theReference->revertToStart();
    resetToStart();
    return res;
}

NDMaterial *TimeVaryingECMaterial::getCopy(void)
{
    TimeVaryingECMaterial *c = new TimeVaryingECMaterial();
    c->setTag(this->getTag());
    c->theReference = theReference != 0 ? theReference->getCopy() : 0;
    c->law = law;
    c->C0 = C0;
    c->C0inv = C0inv;
    c->trial = trial;
    c->committed = committed;
    c->trialSinceCommit = trialSinceCommit;
    c->referenceInSync = referenceInSync;
    c->pristine = pristine;
    c->inputEpsStar = inputEpsStar;
    c->committedEpsStar = committedEpsStar;
    c->inputDriver = inputDriver;
    c->committedDriver = committedDriver;
    c->useLocalDriver = useLocalDriver;
    c->committedUseLocalDriver = committedUseLocalDriver;
    return c;
}

NDMaterial *TimeVaryingECMaterial::getCopy(const char *code)
{
    if (strcmp(code, "ThreeDimensional") == 0 || strcmp(code, "3D") == 0)
        return this->getCopy();
    opserr << "nDMaterial TimeVaryingEC: only ThreeDimensional is available, not " << code << "\n";
    return 0;
}

const char *TimeVaryingECMaterial::getType(void) const
{
    return "ThreeDimensional";
}

int TimeVaryingECMaterial::getOrder(void) const
{
    return theReference != 0 ? 6 : 0;
}

void TimeVaryingECMaterial::Print(OPS_Stream &s, int flag)
{
    s << "TimeVaryingEC (EXPERIMENTAL) tag: " << this->getTag() << endln;
    s << "  reference material tag: " << (theReference ? theReference->getTag() : -1) << endln;
    s << "  law points: " << (law ? law->t.Size() : 0) << endln;
    s << "  driver: " << trial.driver << "  K " << trial.K << "  G " << trial.G
      << "  A " << trial.A << endln;
    s << "  stress: " << trial.sig;
}

int TimeVaryingECMaterial::packState(const State &s, Vector &d, int pos)
{
    const Vector *v[7] = {&s.eps, &s.epsStar, &s.epsM, &s.epsP, &s.sigP, &s.sig, &s.eta};
    for (int k = 0; k < 7; k++)
        for (int i = 0; i < 6; i++)
            d(pos++) = (*v[k])(i);
    d(pos++) = s.driver; d(pos++) = s.K; d(pos++) = s.G; d(pos++) = s.A;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            d(pos++) = s.Calg(i, j);
    return pos;
}

int TimeVaryingECMaterial::unpackState(State &s, const Vector &d, int pos)
{
    Vector *v[7] = {&s.eps, &s.epsStar, &s.epsM, &s.epsP, &s.sigP, &s.sig, &s.eta};
    for (int k = 0; k < 7; k++)
        for (int i = 0; i < 6; i++)
            (*v[k])(i) = d(pos++);
    s.driver = d(pos++); s.K = d(pos++); s.G = d(pos++); s.A = d(pos++);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            s.Calg(i, j) = d(pos++);
    return pos;
}

int TimeVaryingECMaterial::sendSelf(int commitTag, Channel &theChannel)
{
    if (theReference == 0 || !law)
        return -1;
    int dataTag = this->getDbTag();
    if (referenceDbTag == 0)
        referenceDbTag = theChannel.getDbTag();
    int n = law->t.Size();

    static ID idData(11);
    idData(0) = TVEC_VERSION;
    idData(1) = this->getTag();
    idData(2) = theReference->getClassTag();
    idData(3) = referenceDbTag;
    idData(4) = n;
    idData(5) = trialSinceCommit ? 1 : 0;
    idData(6) = referenceInSync ? 1 : 0;
    idData(7) = useLocalDriver ? 1 : 0;
    idData(8) = committedUseLocalDriver ? 1 : 0;
    idData(9) = TVEC_STATE_SIZE;
    idData(10) = pristine ? 1 : 0;
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
        opserr << "TimeVaryingECMaterial::sendSelf - failed to send ID\n";
        return -1;
    }

    Vector data(4 * n + 72 + 2 * TVEC_STATE_SIZE + 14);
    int pos = 0;
    for (int i = 0; i < n; i++) {
        data(pos++) = law->t(i); data(pos++) = law->E(i);
        data(pos++) = law->K(i); data(pos++) = law->A(i);
    }
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) {
            data(pos++) = C0(i, j);
            data(pos++) = C0inv(i, j);
        }
    pos = packState(committed, data, pos);
    pos = packState(trial, data, pos);
    for (int i = 0; i < 6; i++) {
        data(pos++) = inputEpsStar(i);
        data(pos++) = committedEpsStar(i);
    }
    data(pos++) = inputDriver;
    data(pos++) = committedDriver;
    if (theChannel.sendVector(dataTag, commitTag, data) < 0) {
        opserr << "TimeVaryingECMaterial::sendSelf - failed to send data\n";
        return -1;
    }

    theReference->setDbTag(referenceDbTag);
    if (theReference->sendSelf(commitTag, theChannel) < 0) {
        opserr << "TimeVaryingECMaterial::sendSelf - the reference material failed to send itself\n";
        return -1;
    }
    return 0;
}

int TimeVaryingECMaterial::recvSelf(int commitTag, Channel &theChannel,
                                    FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    static ID idData(11);
    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
        opserr << "TimeVaryingECMaterial::recvSelf - failed to receive ID\n";
        return -1;
    }
    if (idData(0) != TVEC_VERSION || idData(9) != TVEC_STATE_SIZE) {
        opserr << "TimeVaryingECMaterial::recvSelf - incompatible data version "
               << idData(0) << " (expected " << TVEC_VERSION << ")\n";
        return -1;
    }
    this->setTag(idData(1));
    int refClassTag = idData(2);
    referenceDbTag = idData(3);
    int n = idData(4);
    if (n < 1) {
        opserr << "TimeVaryingECMaterial::recvSelf - invalid law size\n";
        return -1;
    }

    Vector data(4 * n + 72 + 2 * TVEC_STATE_SIZE + 14);
    if (theChannel.recvVector(dataTag, commitTag, data) < 0) {
        opserr << "TimeVaryingECMaterial::recvSelf - failed to receive data\n";
        return -1;
    }
    std::shared_ptr<Law> l = std::make_shared<Law>();
    l->t = Vector(n); l->E = Vector(n); l->K = Vector(n); l->A = Vector(n);
    int pos = 0;
    for (int i = 0; i < n; i++) {
        l->t(i) = data(pos++); l->E(i) = data(pos++);
        l->K(i) = data(pos++); l->A(i) = data(pos++);
    }
    law = l;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) {
            C0(i, j) = data(pos++);
            C0inv(i, j) = data(pos++);
        }
    pos = unpackState(committed, data, pos);
    pos = unpackState(trial, data, pos);
    for (int i = 0; i < 6; i++) {
        inputEpsStar(i) = data(pos++);
        committedEpsStar(i) = data(pos++);
    }
    inputDriver = data(pos++);
    committedDriver = data(pos++);
    trialSinceCommit = idData(5) != 0;
    useLocalDriver = idData(7) != 0;
    committedUseLocalDriver = idData(8) != 0;
    pristine = idData(10) != 0;

    if (theReference == 0 || theReference->getClassTag() != refClassTag) {
        if (theReference != 0)
            delete theReference;
        theReference = theBroker.getNewNDMaterial(refClassTag);
        if (theReference == 0) {
            opserr << "TimeVaryingECMaterial::recvSelf - the broker has no material with class tag "
                   << refClassTag << "\n";
            return -1;
        }
    }
    theReference->setDbTag(referenceDbTag);
    if (theReference->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "TimeVaryingECMaterial::recvSelf - the reference material failed to receive itself\n";
        return -1;
    }
    // the reference is transferred in its committed state only; a pending
    // trial is redone from it before the next commit
    referenceInSync = false;
    return 0;
}

int TimeVaryingECMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;
    // 4001: isotropic eigenstrain (normal components), as in TimeVarying
    if (strcmp(argv[0], "initNormalStrain") == 0) {
        param.setValue(inputEpsStar(0));
        return param.addObject(4001, this);
    }
    // 4002: local driver; once set it replaces domain time
    if (strcmp(argv[0], "evolutionVariable") == 0 || strcmp(argv[0], "localDriver") == 0) {
        param.setValue(inputDriver);
        return param.addObject(4002, this);
    }
    // parameters of the reference material are not forwarded: changing them
    // would invalidate C0
    return -1;
}

int TimeVaryingECMaterial::updateParameter(int parameterID, Information &info)
{
    switch (parameterID) {
    case 4001:
        if (!std::isfinite(info.theDouble))
            return -1;
        inputEpsStar.Zero();
        inputEpsStar(0) = info.theDouble;
        inputEpsStar(1) = info.theDouble;
        inputEpsStar(2) = info.theDouble;
        break;
    case 4002:
        if (!std::isfinite(info.theDouble))
            return -1;
        inputDriver = info.theDouble;
        useLocalDriver = true;
        break;
    default:
        return -1;
    }
    // Re-evaluate the trial at the current strain with the new input. An
    // analysis forms its first unbalance from the stress before it updates
    // the elements, so without this a step that only changes the driver or
    // the eigenstrain starts from the old stress, solves a zero increment
    // and passes a displacement-increment test with the domain out of
    // equilibrium.
    if (theReference != 0)
        return this->setTrialStrain(trial.eps) == 0 ? 0 : -1;
    return 0;
}

Response *TimeVaryingECMaterial::setResponse(const char **argv, int argc, OPS_Stream &s)
{
    if (argc < 1)
        return NDMaterial::setResponse(argv, argc, s);
    const char *r = argv[0];
    if (strcmp(r, "mechanicalStrain") == 0)
        return new MaterialResponse(this, 101, trial.epsM);
    if (strcmp(r, "eigenStrain") == 0)
        return new MaterialResponse(this, 102, trial.epsStar);
    if (strcmp(r, "projectedStrain") == 0)
        return new MaterialResponse(this, 103, trial.epsP);
    if (strcmp(r, "projectedStress") == 0)
        return new MaterialResponse(this, 104, trial.sigP);
    if (strcmp(r, "properties") == 0)
        return new MaterialResponse(this, 105, Vector(4));
    if (strcmp(r, "stress") == 0 || strcmp(r, "stresses") == 0 ||
        strcmp(r, "strain") == 0 || strcmp(r, "strains") == 0 ||
        strcmp(r, "Tangent") == 0 || strcmp(r, "tangent") == 0)
        return NDMaterial::setResponse(argv, argc, s);
    // other responses (e.g. plastic strain) are those of the reference
    // material, in the reference space
    return theReference != 0 ? theReference->setResponse(argv, argc, s) : 0;
}

int TimeVaryingECMaterial::getResponse(int responseID, Information &matInfo)
{
    switch (responseID) {
    case 101: return matInfo.setVector(trial.epsM);
    case 102: return matInfo.setVector(trial.epsStar);
    case 103: return matInfo.setVector(trial.epsP);
    case 104: return matInfo.setVector(trial.sigP);
    case 105: {
        static Vector p(4);
        p(0) = trial.K; p(1) = trial.G; p(2) = trial.A; p(3) = trial.driver;
        return matInfo.setVector(p);
    }
    default:
        return NDMaterial::getResponse(responseID, matInfo);
    }
}
