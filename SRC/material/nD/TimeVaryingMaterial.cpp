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

// José L. Larenas & José A. Abell (UANDES)
// Massimo Petracca - ASDEA Software, Italy
//
// A Wapper material that allow's time varying stiffness and
// resistance properties for the wrapped material.
//

#include <TimeVaryingMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <Parameter.h>
#include <Domain.h>
#include <ID.h>
#include <Information.h>
#include <classTags.h>
#include <cmath>


void *OPS_TimeVaryingMaterial(void)
{
    const char *usage =
        "nDMaterial TimeVarying $tag $projMatTag $N t1 .. tN E1 .. EN K1 .. KN A1 .. AN\n";

    // check arguments
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "nDMaterial TimeVarying Error: too few arguments.\n" << usage;
        return nullptr;
    }

    // get integer data
    int iData[2];
    int numData = 2;
    if (OPS_GetInt(&numData, iData) != 0)  {
        opserr << "nDMaterial TimeVarying Error: invalid nDMaterial tags.\n";
        return nullptr;
    }

    int Ndatapoints = 0;
    numData = 1;
    if (OPS_GetInt(&numData, &Ndatapoints) != 0 || Ndatapoints < 1)  {
        opserr << "nDMaterial TimeVarying Error: invalid number of data points N\n";
        return nullptr;
    }
    if (OPS_GetNumRemainingInputArgs() < 4 * Ndatapoints) {
        opserr << "nDMaterial TimeVarying Error: expected " << 4 * Ndatapoints
               << " values (t, E, K, A histories) for tag " << iData[0] << ".\n" << usage;
        return nullptr;
    }

    Vector t(Ndatapoints);
    Vector E(Ndatapoints);
    Vector K(Ndatapoints);
    Vector A(Ndatapoints);
    Vector *histories[4] = {&t, &E, &K, &A};
    const char *names[4] = {"t", "E", "K", "A"};
    for (int h = 0; h < 4; ++h) {
        numData = Ndatapoints;
        if (OPS_GetDouble(&numData, &(*histories[h])(0)) != 0) {
            opserr << "nDMaterial TimeVarying Error: reading " << names[h]
                   << " data for nDMaterial TimeVarying with tag " << iData[0] << ".\n";
            return nullptr;
        }
    }

    // validate histories
    for (int k = 0; k < Ndatapoints; ++k) {
        if (k > 0 && t(k) <= t(k - 1)) {
            opserr << "nDMaterial TimeVarying Error: times must be strictly increasing (tag " << iData[0] << ").\n";
            return nullptr;
        }
        if (A(k) <= 0.0 || K(k) <= 0.0 || E(k) <= 0.0 || E(k) >= 9.0 * K(k)) {
            opserr << "nDMaterial TimeVarying Error: need A > 0, K > 0 and 0 < E < 9K at every point (tag "
                   << iData[0] << ", point " << k + 1 << ").\n";
            return nullptr;
        }
    }

    // get the projected material to map
    NDMaterial *theProjMaterial = OPS_getNDMaterial(iData[1]);
    if (theProjMaterial == 0) {
        opserr << "WARNING: nDMaterial does not exist.\n";
        opserr << "nDMaterial: " << iData[1] << "\n";
        opserr << "nDMaterial TimeVarying: " << iData[0] << "\n";
        return nullptr;
    }

    {
        NDMaterial *test = theProjMaterial->getCopy("ThreeDimensional");
        if (test == 0) {
            opserr << "nDMaterial TimeVarying Error: material " << iData[1]
                   << " has no ThreeDimensional version.\n";
            return nullptr;
        }
        delete test;
    }

    // create the TimeVarying wrapper
    NDMaterial* theTimeVaryingMaterial = new TimeVaryingMaterial(
        iData[0],
        *theProjMaterial,
        t, E, K , A);
    if (theTimeVaryingMaterial == 0) {
        opserr << "nDMaterial TimeVarying Error: failed to allocate a new material.\n";
        return nullptr;
    }

    // done
    return theTimeVaryingMaterial;
}


TimeVaryingMaterial::TimeVaryingMaterial()
    : NDMaterial(0, ND_TAG_TimeVaryingMaterial)
{

}

TimeVaryingMaterial::TimeVaryingMaterial(
    int tag,
    NDMaterial &theProjMat,
    const Vector& t_, const Vector& E_, const Vector& K_, const Vector& A_)
    : NDMaterial(tag, ND_TAG_TimeVaryingMaterial)
{
    // copy the isotropic material
    theProjectedMaterial = theProjMat.getCopy("ThreeDimensional");
    if (theProjectedMaterial == 0) {
        opserr << "nDMaterial TimeVarying Error: failed to get a (3D) copy of the isotropic material\n";
    }

    std::shared_ptr<History> h = std::make_shared<History>();
    h->t = t_; h->E = E_; h->K = K_; h->A = A_;
    hist = h;

    new_time_step = true;
}

TimeVaryingMaterial::~TimeVaryingMaterial()
{
    if (theProjectedMaterial)
    {
        delete theProjectedMaterial;
        theProjectedMaterial = 0;
    }
}

double TimeVaryingMaterial::getRho(void)
{
    return theProjectedMaterial->getRho();
}

double TimeVaryingMaterial::currentDriver(void) const
{
    if (use_local_evolution_variable)
        return local_evolution_variable;
    Domain *dom = ops_TheActiveDomain != 0 ? ops_TheActiveDomain : OPS_GetDomain();
    return dom != 0 ? dom->getCurrentTime() : 0.0;
}

int TimeVaryingMaterial::computeMap(double current_driver)
{
    // Get the current parameters at current time
    getParameters(current_driver);

    int tag = this->getTag();
    double Ex  = current_E;  double Ey  = current_E;  double Ez  = current_E;
    double Gxy = current_G;  double Gyz = current_G;  double Gzx = current_G;
    double vxy = current_nu; double vyz = current_nu; double vzx = current_nu;

    // compute the initial orthotropic constitutive tensor
    static Matrix C0(6, 6);
    C0.Zero();
    double vyx = vxy * Ey / Ex;
    double vzy = vyz * Ez / Ey;
    double vxz = vzx * Ex / Ez;
    double d = (1.0 - vxy * vyx - vyz * vzy - vzx * vxz - 2.0 * vxy * vyz * vzx) / (Ex * Ey * Ez);
    C0(0, 0) = (1.0 - vyz * vzy) / (Ey * Ez * d);
    C0(1, 1) = (1.0 - vzx * vxz) / (Ez * Ex * d);
    C0(2, 2) = (1.0 - vxy * vyx) / (Ex * Ey * d);
    C0(1, 0) = (vxy + vxz * vzy) / (Ez * Ex * d);
    C0(0, 1) = C0(1, 0);
    C0(2, 0) = (vxz + vxy * vyz) / (Ex * Ey * d);
    C0(0, 2) = C0(2, 0);
    C0(2, 1) = (vyz + vxz * vyx) / (Ex * Ey * d);
    C0(1, 2) = C0(2, 1);
    C0(3, 3) = Gxy;
    C0(4, 4) = Gyz;
    C0(5, 5) = Gzx;

    // compute the Asigma and its inverse
    if (current_A <= 0 ) {
        opserr << "nDMaterial TimeVarying Error: A must be greater than 0 for tag = " << tag << "\n";
        opserr << "A = " << current_A << endln;
        return -1;
    }

    // compute the initial projected constitutive tensor and its inverse
    static Matrix C0proj(6, 6);
    static Matrix C0proj_inv(6, 6);
    C0proj = theProjectedMaterial->getInitialTangent();
    int res = C0proj.Invert(C0proj_inv);
    if (res < 0) {
        opserr << "nDMaterial TimeVarying Error: the isotropic material gave a singular initial tangent.\n";
        return -1;
    }

    // compute the strain tensor map inv(C0_iso) * Asigma * C0_ortho
    Aepsilon.addMatrixProduct(0.0, C0proj_inv, C0, current_A);
    mapReady = true;
    return 0;
}

int TimeVaryingMaterial::setTrialStrain(const Vector & strain)
{
    if (theProjectedMaterial == 0)
        return -1;
    for (int i = 0; i < strain.Size(); i++) {
        if (!std::isfinite(strain(i))) {
            opserr << "nDMaterial TimeVarying Error: non-finite trial strain (tag "
                   << this->getTag() << "); previous trial kept\n";
            return -1;
        }
    }

    // Compute the real strain increment. The members are only assigned once
    // the trial has succeeded, so a rejected trial leaves the previous one.
    static Vector eps_real_trial(6);
    static Vector depsilon_real(6);
    eps_real_trial = strain - epsilon_internal ; // epsilon_internal comes from thermal analysis
    depsilon_real = eps_real_trial - epsilon_real_n ; // depsilon_real = epsilon_real_new - epsilon_real_old

    int res = computeMap(currentDriver());
    if (res != 0)
        return res;

    //Compute the projected strain increment
    static Vector depsilon_proj(6);
    depsilon_proj.addMatrixVector(0.0, Aepsilon, depsilon_real, 1.0); // depsilon_proj = Aepsilon(t) * depsilon_real

    //Compute the projected total strain
    static Vector epsilon_proj(6);
    epsilon_proj = epsilon_proj_n;
    epsilon_proj.addVector(1.0, depsilon_proj, 1.0); // epsilon_proj = epsilon_proj_old + depsilon_proj

    // call projected material with the projected total strain
    res = theProjectedMaterial->setTrialStrain(epsilon_proj);
    if (res != 0) {
        opserr << "TimeVaryingMaterial::setTrialStrain\n";
        opserr << "nDMaterial Projected Material Error: the isotropic material failed in setTrialStrain.\n";
        return res;
    }

    epsilon_real = eps_real_trial;
    trialSinceCommit = true;
    reverted = false;

    //cool! we're done!
    return 0;
}

const Vector &TimeVaryingMaterial::getStrain(void)
{
    // the strain to report back is the real strain plus the internal,
    // the internal gets removed before converting it to "real"
    static Vector epsilon_return(6);
    epsilon_return = epsilon_real + epsilon_internal;
    return epsilon_return;
}

const Vector &TimeVaryingMaterial::getStress(void)
{
    // after a revert, the committed stress (the projected material's own
    // stress may still be that of the reverted trial)
    if (reverted) {
        sigma_real = sigma_real_n;
        return sigma_real;
    }

    // stress in projected space
    const Vector& sigma_proj = theProjectedMaterial->getStress();

    // compute the projected stress increment
    static Vector dsigma_proj(6);
    dsigma_proj = sigma_proj - sigma_proj_n; // dsigma_proj = sigma_proj_new - sigma_proj_old

    // rescale the stress increment to back to real space
    static Vector dsigma_real(6);
    // dsigma_real = Asigma^-1 * dsigma_proj
    dsigma_real = dsigma_proj / current_A;  // dsigma_real = Asigma^-1 * dsigma_proj

    //add real stress increment to the previous real stress
    sigma_real = sigma_real_n + dsigma_real;

    return sigma_real;
}

const Matrix &TimeVaryingMaterial::getTangent(void)
{
    // after a revert, the committed tangent (see getStress)
    if (reverted)
        return committed_tangent;
    if (!mapReady)
        return this->getInitialTangent();

    // tensor in isotropic space
    const Matrix &C_proj = theProjectedMaterial->getTangent();

    // compute orthotripic tangent
    static Matrix C_real(6, 6);
    C_real.addMatrixProduct(0.0, C_proj, Aepsilon, 1 / current_A);

    return C_real;
}

const Matrix &TimeVaryingMaterial::getInitialTangent(void)
{
    // before the first trial Aepsilon has not been built: build it at the
    // current driver (this only refreshes the cached properties)
    if (!mapReady && theProjectedMaterial != 0)
        computeMap(currentDriver());

    // elasticity tensor in projected space
    const Matrix& C_proj = theProjectedMaterial->getInitialTangent();

    // compute the real tangent from the projected one
    static Matrix C_real(6, 6);
    C_real.addMatrixProduct(0.0, C_proj, Aepsilon, 1 / current_A);
    return C_real;
}

int TimeVaryingMaterial::commitState(void)
{
    // Nothing new since the last commit (e.g. a commit right after a revert):
    // do not commit the projected material, whose trial may be the reverted one.
    if (!trialSinceCommit)
        return 0;

    // Materialize and snapshot a complete trial state before advancing either material,
    // so the committed state does not depend on whether getStress() was called first.
    Vector sigma_proj_trial(theProjectedMaterial->getStress());
    Vector dsigma_proj = sigma_proj_trial - sigma_proj_n;
    Vector sigma_real_trial = sigma_real_n + dsigma_proj / current_A;
    Vector epsilon_real_trial(epsilon_real);
    Vector epsilon_proj_trial(theProjectedMaterial->getStrain());

    int result = theProjectedMaterial->commitState();
    if (result == 0) {
        sigma_real_n = sigma_real_trial;
        sigma_proj_n = sigma_proj_trial;
        epsilon_real_n = epsilon_real_trial;
        epsilon_proj_n = epsilon_proj_trial;
        epsilon_new_n = epsilon_new;
        new_time_step = true;
        // committed copies for an immediate revert
        committed_E = current_E;
        committed_G = current_G;
        committed_nu = current_nu;
        committed_A = current_A;
        committed_Aepsilon = Aepsilon;
        committed_tangent.addMatrixProduct(0.0, theProjectedMaterial->getTangent(), Aepsilon, 1 / current_A);
        trialSinceCommit = false;
        reverted = false;
    }
    return result;
}

int TimeVaryingMaterial::revertToLastCommit(void)
{
    int result = theProjectedMaterial->revertToLastCommit();
    if (result == 0) {
        sigma_real = sigma_real_n;
        epsilon_real = epsilon_real_n;
        epsilon_new = epsilon_new_n;
        // restore the committed map; queries are answered from the committed
        // copies until the next trial
        current_E = committed_E;
        current_G = committed_G;
        current_nu = committed_nu;
        current_A = committed_A;
        Aepsilon = committed_Aepsilon;
        // the cached properties no longer match last_parameter_driver
        last_parameter_driver = std::numeric_limits<double>::quiet_NaN();
        trialSinceCommit = false;
        reverted = mapReady;
    }
    return result;
}

int TimeVaryingMaterial::revertToStart(void)
{
    int result = theProjectedMaterial->revertToStart();
    if (result == 0) {
        sigma_real.Zero();
        epsilon_real.Zero();
        sigma_real_n.Zero();
        sigma_proj_n.Zero();
        epsilon_real_n.Zero();
        epsilon_proj_n.Zero();
        epsilon_new_n.Zero();
        current_E = current_G = current_nu = current_A = 0.0;
        committed_E = committed_G = committed_nu = committed_A = 0.0;
        Aepsilon.Zero();
        committed_Aepsilon.Zero();
        committed_tangent.Zero();
        last_parameter_driver = std::numeric_limits<double>::quiet_NaN();
        new_time_step = true;
        trialSinceCommit = false;
        reverted = false;
        mapReady = false;
    }
    return result;
}

NDMaterial * TimeVaryingMaterial::getCopy(void)
{
    TimeVaryingMaterial *theCopy = new TimeVaryingMaterial();
    theCopy->setTag(getTag());
    theCopy->theProjectedMaterial = theProjectedMaterial->getCopy("ThreeDimensional");
    theCopy->Aepsilon = Aepsilon;
    theCopy->epsilon_internal = epsilon_internal;
    theCopy->sigma_real = sigma_real;
    theCopy->epsilon_real = epsilon_real;
    theCopy->sigma_real_n = sigma_real_n;
    theCopy->sigma_proj_n = sigma_proj_n;
    theCopy->epsilon_real_n = epsilon_real_n;
    theCopy->epsilon_proj_n = epsilon_proj_n;
    theCopy->epsilon_new = epsilon_new;
    theCopy->epsilon_new_n = epsilon_new_n;
    theCopy->current_E = current_E;
    theCopy->current_G = current_G;
    theCopy->current_nu = current_nu;
    theCopy->current_A = current_A;
    theCopy->new_time_step = new_time_step;
    theCopy->local_evolution_variable = local_evolution_variable;
    theCopy->use_local_evolution_variable = use_local_evolution_variable;
    theCopy->hist = hist;
    theCopy->last_parameter_driver = last_parameter_driver;
    theCopy->committed_E = committed_E;
    theCopy->committed_G = committed_G;
    theCopy->committed_nu = committed_nu;
    theCopy->committed_A = committed_A;
    theCopy->committed_Aepsilon = committed_Aepsilon;
    theCopy->committed_tangent = committed_tangent;
    theCopy->trialSinceCommit = trialSinceCommit;
    theCopy->reverted = reverted;
    theCopy->mapReady = mapReady;
    return theCopy;
}

NDMaterial* TimeVaryingMaterial::getCopy(const char* code)
{
    if (strcmp(code, "ThreeDimensional") == 0)
        return getCopy();
    return NDMaterial::getCopy(code);
}

const char* TimeVaryingMaterial::getType(void) const
{
    return "ThreeDimensional";
}

int TimeVaryingMaterial::getOrder(void) const
{
    return 6;
}

void TimeVaryingMaterial::Print(OPS_Stream &s, int flag)
{
    s << "Time Varying Material, tag: " << this->getTag() << "\n";
}

// sendSelf layout version (0 was the empty implementation)
static const int TV_SEND_VERSION = 1;

int TimeVaryingMaterial::sendSelf(int commitTag, Channel &theChannel)
{
    if (theProjectedMaterial == 0 || !hist)
        return -1;
    int dataTag = this->getDbTag();
    if (projectedDbTag == 0)
        projectedDbTag = theChannel.getDbTag();
    int n = hist->t.Size();

    static ID idData(10);
    idData(0) = TV_SEND_VERSION;
    idData(1) = this->getTag();
    idData(2) = theProjectedMaterial->getClassTag();
    idData(3) = projectedDbTag;
    idData(4) = n;
    idData(5) = use_local_evolution_variable ? 1 : 0;
    idData(6) = new_time_step ? 1 : 0;
    idData(7) = trialSinceCommit ? 1 : 0;
    idData(8) = reverted ? 1 : 0;
    idData(9) = mapReady ? 1 : 0;
    if (theChannel.sendID(dataTag, commitTag, idData) < 0) {
        opserr << "TimeVaryingMaterial::sendSelf - failed to send ID\n";
        return -1;
    }

    Vector data(4 * n + 9 * 6 + 3 * 36 + 10);
    int pos = 0;
    for (int i = 0; i < n; i++) {
        data(pos++) = hist->t(i); data(pos++) = hist->E(i);
        data(pos++) = hist->K(i); data(pos++) = hist->A(i);
    }
    const Vector *v[9] = {&epsilon_internal, &sigma_real, &epsilon_real, &epsilon_new,
                          &sigma_real_n, &sigma_proj_n, &epsilon_real_n, &epsilon_proj_n, &epsilon_new_n};
    for (int k = 0; k < 9; k++)
        for (int i = 0; i < 6; i++)
            data(pos++) = (*v[k])(i);
    const Matrix *m[3] = {&Aepsilon, &committed_Aepsilon, &committed_tangent};
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                data(pos++) = (*m[k])(i, j);
    data(pos++) = current_E; data(pos++) = current_G; data(pos++) = current_nu; data(pos++) = current_A;
    data(pos++) = committed_E; data(pos++) = committed_G; data(pos++) = committed_nu; data(pos++) = committed_A;
    data(pos++) = last_parameter_driver;
    data(pos++) = local_evolution_variable;
    if (theChannel.sendVector(dataTag, commitTag, data) < 0) {
        opserr << "TimeVaryingMaterial::sendSelf - failed to send data\n";
        return -1;
    }

    theProjectedMaterial->setDbTag(projectedDbTag);
    if (theProjectedMaterial->sendSelf(commitTag, theChannel) < 0) {
        opserr << "TimeVaryingMaterial::sendSelf - the projected material failed to send itself\n";
        return -1;
    }
    return 0;
}

int TimeVaryingMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
    int dataTag = this->getDbTag();
    static ID idData(10);
    if (theChannel.recvID(dataTag, commitTag, idData) < 0) {
        opserr << "TimeVaryingMaterial::recvSelf - failed to receive ID\n";
        return -1;
    }
    if (idData(0) != TV_SEND_VERSION) {
        opserr << "TimeVaryingMaterial::recvSelf - incompatible data version " << idData(0) << "\n";
        return -1;
    }
    this->setTag(idData(1));
    int projClassTag = idData(2);
    projectedDbTag = idData(3);
    int n = idData(4);
    if (n < 1) {
        opserr << "TimeVaryingMaterial::recvSelf - invalid history size\n";
        return -1;
    }
    use_local_evolution_variable = idData(5) != 0;
    new_time_step = idData(6) != 0;
    trialSinceCommit = idData(7) != 0;
    reverted = idData(8) != 0;
    mapReady = idData(9) != 0;

    Vector data(4 * n + 9 * 6 + 3 * 36 + 10);
    if (theChannel.recvVector(dataTag, commitTag, data) < 0) {
        opserr << "TimeVaryingMaterial::recvSelf - failed to receive data\n";
        return -1;
    }
    std::shared_ptr<History> h = std::make_shared<History>();
    h->t = Vector(n); h->E = Vector(n); h->K = Vector(n); h->A = Vector(n);
    int pos = 0;
    for (int i = 0; i < n; i++) {
        h->t(i) = data(pos++); h->E(i) = data(pos++);
        h->K(i) = data(pos++); h->A(i) = data(pos++);
    }
    hist = h;
    Vector *v[9] = {&epsilon_internal, &sigma_real, &epsilon_real, &epsilon_new,
                    &sigma_real_n, &sigma_proj_n, &epsilon_real_n, &epsilon_proj_n, &epsilon_new_n};
    for (int k = 0; k < 9; k++)
        for (int i = 0; i < 6; i++)
            (*v[k])(i) = data(pos++);
    Matrix *m[3] = {&Aepsilon, &committed_Aepsilon, &committed_tangent};
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                (*m[k])(i, j) = data(pos++);
    current_E = data(pos++); current_G = data(pos++); current_nu = data(pos++); current_A = data(pos++);
    committed_E = data(pos++); committed_G = data(pos++); committed_nu = data(pos++); committed_A = data(pos++);
    last_parameter_driver = data(pos++);
    local_evolution_variable = data(pos++);

    if (theProjectedMaterial == 0 || theProjectedMaterial->getClassTag() != projClassTag) {
        if (theProjectedMaterial != 0)
            delete theProjectedMaterial;
        theProjectedMaterial = theBroker.getNewNDMaterial(projClassTag);
        if (theProjectedMaterial == 0) {
            opserr << "TimeVaryingMaterial::recvSelf - the broker has no material with class tag "
                   << projClassTag << "\n";
            return -1;
        }
    }
    theProjectedMaterial->setDbTag(projectedDbTag);
    if (theProjectedMaterial->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "TimeVaryingMaterial::recvSelf - the projected material failed to receive itself\n";
        return -1;
    }
    // the projected material arrives in its committed state only (e.g.
    // J2Plasticity restores its committed internal variables, not its current
    // stress), so queries are answered from the committed copies until the
    // next trial; a pending trial is dropped
    trialSinceCommit = false;
    reverted = mapReady;
    sigma_real = sigma_real_n;
    epsilon_real = epsilon_real_n;
    return 0;
}

int TimeVaryingMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    // 4000 - init strain
    if (strcmp(argv[0], "initNormalStrain") == 0) {
        double initNormalStrain = epsilon_internal(0);
        param.setValue(initNormalStrain);
        return param.addObject(4001, this);
    }

    // 4002 - optional material-point evolution variable. Once updated, this
    // value replaces global domain time as the interpolation driver.
    if (strcmp(argv[0], "evolutionVariable") == 0 ||
            strcmp(argv[0], "localDriver") == 0) {
        param.setValue(local_evolution_variable);
        return param.addObject(4002, this);
    }

    // forward to the adapted (isotropic) material
    return theProjectedMaterial->setParameter(argv, argc, param);
}


int TimeVaryingMaterial::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {

    case 4001:
    {
        double initNormalStrain = info.theDouble; // alpha * (Temp(t)  - Temp())
        if (!std::isfinite(initNormalStrain))
            return -1;
        epsilon_internal.Zero();
        epsilon_internal(0) = initNormalStrain;
        epsilon_internal(1) = initNormalStrain;
        epsilon_internal(2) = initNormalStrain;
        break;
    }

    case 4002:
    {
        if (!std::isfinite(info.theDouble))
            return -1;
        local_evolution_variable = info.theDouble;
        use_local_evolution_variable = true;
        new_time_step = true;
        break;
    }

    // default
    default:
        return -1;
    }

    // KNOWN LIMITATION: the input is only stored; the stress changes at the
    // next setTrialStrain. With the Plain or Penalty handlers a step that only
    // changes this input starts from the old stress and can pass its
    // convergence test at iteration 1 with the domain out of equilibrium. The
    // Transformation handler updates the elements in applyLoad and is not
    // affected. Re-evaluating here was tried and rejected: the extra trials
    // change the path of trial-history-dependent wrapped materials
    // (ASDConcrete3D) and so the results of existing analyses. TimeVaryingEC
    // re-evaluates (its scope is elasticity and perfect J2).
    return 0;
}


Response* TimeVaryingMaterial::setResponse(const char** argv, int argc, OPS_Stream& s)
{
    if (argc > 0) {
        if (strcmp(argv[0], "stress") == 0 ||
                strcmp(argv[0], "stresses") == 0 ||
                strcmp(argv[0], "strain") == 0 ||
                strcmp(argv[0], "strains") == 0 ||
                strcmp(argv[0], "Tangent") == 0 ||
                strcmp(argv[0], "tangent") == 0) {
            // stresses, strain and tangent should be those of this adapter (orthotropic)
            return NDMaterial::setResponse(argv, argc, s);
        }

        else {
            // any other response should be obtained from the adapted (isotropic) material
            return theProjectedMaterial->setResponse(argv, argc, s);
        }
    }
    return NDMaterial::setResponse(argv, argc, s);
}

void TimeVaryingMaterial::getParameters(double time)
{
    int tag = this->getTag();
    // F14: the element's construction-time update() consumes new_time_step at the initial driver
    // value, so without the second test the first analysis step reuses properties evaluated there.
    if (new_time_step || !(time == last_parameter_driver)) {
        double K = 0;
        new_time_step = false;
        last_parameter_driver = time;
        // Find the interval in which 'time' falls within time_history
        int index = 0;
        for (; index < hist->t.Size(); ++index)
        {
            if (hist->t(index) >= time)
            {
                break;
            }
        }

        if (index == 0)
        {
            current_E = hist->E(0) ;
            current_A = hist->A(0) ;
            K = hist->K(0) ;

            current_G  = (3.0 * K * current_E) / (9.0 * K - current_E);
            current_nu = (3.0 * K - current_E) / (6.0 * K);
        }

        else if (index > 0 && index < hist->t.Size())
        {
            // Perform linear interpolation for the values of E, A, and  K
            double t1    = hist->t(index - 1);
            double t2    = hist->t(index);
            double alpha = (time - t1) / (t2 - t1);

            current_E = (1.0 - alpha) * hist->E(index - 1) + alpha * hist->E(index);
            current_A = (1.0 - alpha) * hist->A(index - 1) + alpha * hist->A(index);
            K = (1.0 - alpha) * hist->K(index - 1) + alpha * hist->K(index);

            current_G  = (3.0 * K * current_E) / (9.0 * K - current_E);
            current_nu = (3.0 * K - current_E) / (6.0 * K);
        }

        else {
            current_E = hist->E(hist->t.Size() - 1) ;
            current_A = hist->A(hist->t.Size() - 1) ;
            K = hist->K(hist->t.Size() - 1) ;

            current_G  = (3.0 * K * current_E) / (9.0 * K - current_E);
            current_nu = (3.0 * K - current_E) / (6.0 * K);
        }
    }
}
