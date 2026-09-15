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
#include <Element.h>


void *OPS_TimeVaryingMaterial(void)
{
    opserr << "Using TimeVaryingMaterial" << endln ;
    // check arguments
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 7) {
        opserr <<
               "nDMaterial TimeVarying Error: Few arguments (< 17).\n"
               "nDMaterial TimeVarying $tag $theProjectedMat $Nt t1 t2 t3 .. t_Nt E1 E2 E3 ... E_Nt K1 K2 K3...K_Nt  A1 A2 A3...A_Nt\n";
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
    if (OPS_GetInt(&numData, &Ndatapoints) != 0)  {
        opserr << "nDMaterial TimeVarying Error: error reading Ndatapoints\n";
        return nullptr;
    }

    Vector t(Ndatapoints);
    Vector E(Ndatapoints);
    Vector K(Ndatapoints);
    Vector A(Ndatapoints);

    // get double data
    double * dData = new double[Ndatapoints];

    numData = Ndatapoints;

    //Read time steps
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "nDMaterial TimeVarying Error: reading t data for nDMaterial TimeVarying with tag " << iData[0] << ".\n";
        return nullptr;
    }
    t.setData(dData, Ndatapoints);

    //Read E
    dData = new double[Ndatapoints];
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "nDMaterial TimeVarying Error: reading E data for nDMaterial TimeVarying with tag " << iData[0] << ".\n";
        return nullptr;
    }
    E.setData(dData, Ndatapoints);

    //Read K
    dData = new double[Ndatapoints];
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "nDMaterial TimeVarying Error: reading K data for nDMaterial TimeVarying with tag " << iData[0] << ".\n";
        return nullptr;
    }
    K.setData(dData, Ndatapoints);

    //Read A
    dData = new double[Ndatapoints];
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "nDMaterial TimeVarying Error: reading A data for nDMaterial TimeVarying with tag " << iData[0] << ".\n";
        return nullptr;
    }
    A.setData(dData, Ndatapoints);

    // get the projected material to map
    NDMaterial *theProjMaterial = OPS_getNDMaterial(iData[1]);
    if (theProjMaterial == 0) {
        opserr << "WARNING: nDMaterial does not exist.\n";
        opserr << "nDMaterial: " << iData[1] << "\n";
        opserr << "nDMaterial TimeVarying: " << iData[0] << "\n";
        return nullptr;
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

// int TimeVaryingMaterial::number_of_evolution_laws = 0;
// Vector* TimeVaryingMaterial::time_history = 0;
// Vector* TimeVaryingMaterial::E_history = 0;
// Vector* TimeVaryingMaterial::K_history = 0;
// Vector* TimeVaryingMaterial::A_history = 0;
std::map<int, Vector> TimeVaryingMaterial::time_histories;
std::map<int, Vector> TimeVaryingMaterial::E_histories;
std::map<int, Vector> TimeVaryingMaterial::K_histories;
std::map<int, Vector> TimeVaryingMaterial::A_histories;
std::map<int, bool> TimeVaryingMaterial::new_time_step ;
bool TimeVaryingMaterial::print_strain_once = true;
bool TimeVaryingMaterial::print_stress_once = true;
bool TimeVaryingMaterial::print_commit_once = true;
bool TimeVaryingMaterial::print_tang_once = true;
// double TimeVaryingMaterial::E = 0;
// double TimeVaryingMaterial::G = 0;
// double TimeVaryingMaterial::A = 0;
// double TimeVaryingMaterial::nu = 0;
std::map<int, double> TimeVaryingMaterial::E;
std::map<int, double> TimeVaryingMaterial::G;
std::map<int, double> TimeVaryingMaterial::A;
std::map<int, double> TimeVaryingMaterial::nu;

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
        opserr << "nDMaterial Orthotropic Error: failed to get a (3D) copy of the isotropic material\n";
        exit(-1);
    }

    int Ndatapoints = E_.Size();
    // time_history = new Vector(Ndatapoints);
    // E_history = new Vector(Ndatapoints);
    // K_history = new Vector(Ndatapoints);
    // A_history = new Vector(Ndatapoints);
    time_histories[tag] = Vector(Ndatapoints);
    E_histories[tag] = Vector(Ndatapoints);
    K_histories[tag] = Vector(Ndatapoints);
    A_histories[tag] = Vector(Ndatapoints);
    // *time_history = t_;
    // *E_history = E_;
    // *K_history = K_;
    // *A_history = A_;
    time_histories[tag] = t_;
    E_histories[tag] = E_;
    K_histories[tag] = K_;
    A_histories[tag] = A_;

    new_time_step[tag] = true;

    //This is a new material constructor, hence we save the evolution law
    // evolution_law_id = number_of_evolution_laws;
    //and advance the number of evolution laws.
    // number_of_evolution_laws++;

    opserr << "Created new TimeVaryingMaterial \n";
    // opserr << "  evolution_law_id          = " << evolution_law_id                << endln;
    opserr << "  tag          = " << tag                << endln;
    opserr << "  proj_tag     = " << theProjectedMaterial->getTag() << endln;
    opserr << "  Nt           = " << E_histories[tag].Size()  << endln;
    opserr << "  time_history = " << time_histories[tag]      << endln;
    opserr << "  E_history    = " << E_histories[tag]         << endln;
    opserr << "  K_history    = " << K_histories[tag]         << endln;
    opserr << "  A_history    = " << A_histories[tag]         << endln;
}

TimeVaryingMaterial::~TimeVaryingMaterial()
{
    if (theProjectedMaterial)
    {
        // delete time_history;
        // time_history = 0;
        // delete E_history;
        // E_history = 0;
        // delete K_history;
        // K_history = 0;
        // delete A_history;
        // A_history = 0;
        delete theProjectedMaterial;
        theProjectedMaterial = 0;
    }
}

double TimeVaryingMaterial::getRho(void)
{
    return theProjectedMaterial->getRho();
}

int TimeVaryingMaterial::setTrialStrain(const Vector & strain)
{

    //Compute the real strain increment
    static Vector depsilon_real(6);
    epsilon_real = strain - epsilon_internal ; // epsilon_internal comes from thermal analysis
    depsilon_real = epsilon_real - epsilon_real_n ; // depsilon_real = epsilon_real_new - epsilon_real_old

    // Get the current parameters at current time
    double current_time = OPS_GetDomain()->getCurrentTime();  //should not be used in displacement control user could set the current time as a parameter (check ASDConcrete3D)
    getParameters(current_time);

    int tag = this->getTag();
    double Ex  = E[tag];  double Ey  = E[tag];  double Ez  = E[tag];
    double Gxy = G[tag];  double Gyz = G[tag];  double Gzx = G[tag];
    double vxy = nu[tag]; double vyz = nu[tag]; double vzx = nu[tag];
    // double Asigmaxx   = A; double Asigmayy   = A; double Asigmazz   = A;
    // double Asigmaxyxy = A; double Asigmayzyz = A; double Asigmaxzxz = A;

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
    if (A[tag] <= 0 ) {
        opserr << "nDMaterial TimeVarying Error: A must be greater than 0 for tag = " << tag << "\n";
        opserr << "A = " << A[tag] << endln;
        exit(-1);
    }
    // static Matrix Asigma(6, 6);
    // Asigma.Zero();
    // Asigma(0, 0) = Asigmaxx;
    // Asigma(1, 1) = Asigmayy;
    // Asigma(2, 2) = Asigmazz;
    // Asigma(3, 3) = Asigmaxyxy;
    // Asigma(4, 4) = Asigmayzyz;
    // Asigma(5, 5) = Asigmaxzxz;
    // for (int i = 0; i < 6; ++i)
    //     Asigma_inv(i) = 1.0 / Asigma(i, i);

    // compute the initial projected constitutive tensor and its inverse
    static Matrix C0proj(6, 6);
    static Matrix C0proj_inv(6, 6);
    C0proj = theProjectedMaterial->getInitialTangent();
    int res = C0proj.Invert(C0proj_inv);
    if (res < 0) {
        opserr << "nDMaterial Orthotropic Error: the isotropic material gave a singular initial tangent.\n";
        exit(-1);
    }

    // compute the strain tensor map inv(C0_iso) * Asigma * C0_ortho
    // static Matrix Asigma_C0(6, 6);
    // Asigma_C0.addMatrixProduct(0.0, Asigma, C0, 1.0);
    // Asigma_C0 = A*C0;
    Aepsilon.addMatrixProduct(0.0, C0proj_inv, C0, A[tag]);

    //Compute the projected strain increment
    static Vector depsilon_proj(6);
    depsilon_proj.addMatrixVector(0.0, Aepsilon, depsilon_real, 1.0); // depsilon_proj = Aepsilon(t) * depsilon_real

    //Compute the projected total strain
    static Vector epsilon_proj(6);
    epsilon_proj = epsilon_proj_n;
    epsilon_proj.addVector(1.0, depsilon_proj, 1.0); // epsilon_proj = epsilon_proj_old + depsilon_proj

    my_element_tag = ops_TheActiveElement->getTag();
    // opserr << "my_element_tag = " << my_element_tag << endln;
    // if (my_element_tag == 83 )
    // {
    //     opserr << "@ TimeVaryingMaterial::setTrialStrain" << endln;
    //     opserr << "strain = " << strain << endln;
    //     opserr << "epsilon_internal = " << epsilon_internal << endln;
    //     opserr << "epsilon_real = " << epsilon_real << endln;
    //     opserr << "depsilon_real = " << depsilon_real << endln;
    //     opserr << "depsilon_proj = " << depsilon_proj << endln;
    //     opserr << "epsilon_proj = " << epsilon_proj << endln;
    //     print_strain_once = false;
    //     print_commit_once = true;
    // }

    // call projected material with the projected total strain
    res = theProjectedMaterial->setTrialStrain(epsilon_proj);
    if (res != 0) {
        opserr << "TimeVaryingMaterial::setTrialStrain\n";
        opserr << "nDMaterial Projected Material Error: the isotropic material failed in setTrialStrain.\n";
        return res;
    }

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
    // stress in projected space
    const Vector& sigma_proj = theProjectedMaterial->getStress();

    // compute the projected stress increment
    static Vector dsigma_proj(6);
    dsigma_proj = sigma_proj - sigma_proj_n; // dsigma_proj = sigma_proj_new - sigma_proj_old

    // rescale the stress increment to back to real space
    static Vector dsigma_real(6);
    // for (int i = 0; i < 6; ++i)
    //     dsigma_real(i) = dsigma_proj(i) / A;
    // dsigma_real = Asigma^-1 * dsigma_proj
    // dsigma_real(i) = Asigma_inv(i) * dsigma_proj(i); // dsigma_real = Asigma^-1 * dsigma_proj
    dsigma_real = dsigma_proj / A[this->getTag()];  // dsigma_real = Asigma^-1 * dsigma_proj

    //add real stress increment to the previous real stress
    sigma_real = sigma_real_n + dsigma_real;

    if (my_element_tag == 83)
    {
        opserr << " @TimeVaryingMaterial::getStress mattag = " << this->getTag() << endln;
        opserr << "sigma_real_n = " << sigma_real_n << endln;
        opserr << "dsigma_real = " << dsigma_real << endln;
        opserr << "sigma_real = " << sigma_real << endln;
        // opserr << "Aepsilon = " << Aepsilon << endln;
        print_stress_once = false;
    }

    return sigma_real;
}

const Matrix &TimeVaryingMaterial::getTangent(void)
{
    // tensor in isotropic space
    const Matrix &C_proj = theProjectedMaterial->getTangent();

    // compute orthotripic tangent
    static Matrix C_real(6, 6);
    // static Matrix temp(6, 6);
    // static Matrix cdiff(6, 6);
    // static Matrix invAsigma(6, 6);
    // invAsigma.Zero();
    // for (int i = 0; i < 6; ++i)
    //     invAsigma(i, i) = Asigma_inv(i);
    // temp.addMatrixProduct(0.0, C_proj, Aepsilon, 1.0);
    // C_real.addMatrixProduct(0.0, invAsigma, temp, 1.0);
    // C_real = temp / A;
    C_real.addMatrixProduct(0.0, C_proj, Aepsilon, 1 / A[this->getTag()]);


    // cdiff = C_proj - C_real;
    bool printnow = false;
    // for (int i = 0; i < 6; ++i)
    // {
    //     for (int j = 0; j < 6; ++j)
    //     {
    //         if (abs(cdiff(i, j)) > 1e-4 || cdiff(i, j) != cdiff(i, j) )
    //         {
    //             printnow = true;
    //         }
    //     }
    // }

    // if (my_element_tag == 83  )
    // {
    //     opserr << "@ getTangent" << endln;
    //     opserr << "C_real = " << C_real << endln;
    //     opserr << "C_proj = " << C_proj << endln;
    //     // opserr << "diff = " << cdiff  << endln;
    //     opserr << "Asigma_inv = " << 1 / A[this->getTag()]  << endln;
    //     opserr << "Aepsilon = " << Aepsilon  << endln;
    //     // print_tang_once = false;
    // }


    return C_real;
}

const Matrix &TimeVaryingMaterial::getInitialTangent(void)
{
    opserr << "@TimeVaryingMaterial::getInitialTangent" << endln;

    // elasticity tensor in projected space
    const Matrix& C_proj = theProjectedMaterial->getInitialTangent();

    // compute the real tangent from the projected one
    static Matrix C_real(6, 6);
    static Matrix temp(6, 6);
    // static Matrix invAsigma(6, 6);
    // invAsigma.Zero();
    // for (int i = 0; i < 6; ++i)
    //     invAsigma(i, i) = Asigma_inv(i);
    // temp.addMatrixProduct(0.0, C_proj, Aepsilon, 1.0);
    // C.addMatrixProduct(0.0, invAsigma, temp, 1.0);
    C_real.addMatrixProduct(0.0, C_proj, Aepsilon, 1 / A[this->getTag()]);
    return C_real;


}

int TimeVaryingMaterial::commitState(void)
{
    new_time_step[this->getTag()] = true;

    const Vector& sigma_proj = theProjectedMaterial->getStress();
    const Vector& epsilon_proj = theProjectedMaterial->getStrain();

    print_stress_once = true;
    print_strain_once = true;
    print_tang_once = true;
    sigma_real_n = sigma_real;
    sigma_proj_n = sigma_proj;
    epsilon_real_n = epsilon_real;
    epsilon_proj_n = epsilon_proj;
    epsilon_new_n = epsilon_new;

    if (my_element_tag == 83)
    {
        opserr << "@ TimeVaryingMaterial::commitState(void) " << endln;
        opserr << "sigma_real_n = " << sigma_real_n;
        opserr << "sigma_proj_n = " << sigma_proj_n;
        opserr << "epsilon_real_n = " << epsilon_real_n;
        opserr << "epsilon_proj_n = " << epsilon_proj_n;
        opserr << "epsilon_new_n = " << epsilon_new_n;
        // print_commit_once = false;
    }

    return theProjectedMaterial->commitState();
}

int TimeVaryingMaterial::revertToLastCommit(void)
{
    sigma_real = sigma_real_n ;
    // sigma_proj = sigma_proj_n ;
    epsilon_real = epsilon_real_n ;
    // epsilon_proj = epsilon_proj_n ;
    epsilon_new = epsilon_new_n ;
    return theProjectedMaterial->revertToLastCommit();
}

int TimeVaryingMaterial::revertToStart(void)
{
    sigma_real_n.Zero();
    sigma_proj_n.Zero();
    epsilon_real_n.Zero();
    epsilon_proj_n.Zero();
    epsilon_new_n .Zero();
    return theProjectedMaterial->revertToStart();
}

NDMaterial * TimeVaryingMaterial::getCopy(void)
{
    // opserr << "TimeVaryingMaterial::getCopy" << endln;
    TimeVaryingMaterial *theCopy = new TimeVaryingMaterial();
    theCopy->setTag(getTag());
    theCopy->theProjectedMaterial = theProjectedMaterial->getCopy("ThreeDimensional");
    // theCopy->evolution_law_id = evolution_lsaw_id;
    theCopy->Aepsilon = Aepsilon;
    theCopy->epsilon_internal = epsilon_internal;
    theCopy->sigma_real = sigma_real;
    // theCopy->sigma_proj = sigma_proj;
    theCopy->epsilon_real = epsilon_real;
    // theCopy->epsilon_proj = epsilon_proj;
    theCopy->sigma_real_n = sigma_real_n;
    theCopy->sigma_proj_n = sigma_proj_n;
    theCopy->epsilon_real_n = epsilon_real_n;
    theCopy->epsilon_proj_n = epsilon_proj_n;
    theCopy->epsilon_new = epsilon_new;
    theCopy->epsilon_new_n = epsilon_new_n;
    // theCopy->E = E;
    // theCopy->G = G;
    // theCopy->nu = nu;
    // theCopy->A = A;
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

int TimeVaryingMaterial::sendSelf(int commitTag, Channel &theChannel)
{

    return 0;
}

int TimeVaryingMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{

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

    // forward to the adapted (isotropic) material
    return theProjectedMaterial->setParameter(argv, argc, param);
}


int TimeVaryingMaterial::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {

    case 4001:
    {
        double initNormalStrain = info.theDouble; // alpha * (Temp(t)  - Temp())
        epsilon_internal.Zero();
        epsilon_internal(0) = initNormalStrain;
        epsilon_internal(1) = initNormalStrain;
        epsilon_internal(2) = initNormalStrain;
        return 0;
    }

    // default
    default:
        return -1;
    }
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

void TimeVaryingMaterial::getParameters(double time)//, double& E, double& G, double& nu, double& A)
{
    int tag = this->getTag();
    if (new_time_step[tag]) {
        double K = 0;
        // opserr << "Getting Parameters (E, G, nu, A)" << endln;
        new_time_step[tag] = false;

        // Find the interval in which 'time' falls within time_history
        int index = 0;
        for (; index < time_histories[tag].Size(); ++index)
        {
            if (time_histories[tag](index) >= time)
            {
                break;
            }
        }

        opserr << "index for time (" << time << ") = " << index << endln;

        if (index == 0)
        {
            E[tag] = E_histories[tag](0) ;
            A[tag] = A_histories[tag](0) ;
            K = K_histories[tag](0) ;

            G[tag]  = (3.0 * K * E[tag]) / (9.0 * K - E[tag]);
            nu[tag] = (3.0 * K - E[tag]) / (6.0 * K);
        }

        else if (index > 0 && index < time_histories[tag].Size())
        {
            // Perform linear interpolation for the values of E, A, and  K
            double t1    = time_histories[tag](index - 1);
            double t2    = time_histories[tag](index);
            double alpha = (time - t1) / (t2 - t1);
            // opserr << "t1    = " << t1    << endln ;
            // opserr << "t2    = " << t2    << endln ;
            // opserr << "alpha = " << alpha << endln ;

            E[tag] = (1.0 - alpha) * E_histories[tag](index - 1) + alpha * E_histories[tag](index);
            A[tag] = (1.0 - alpha) * A_histories[tag](index - 1) + alpha * A_histories[tag](index);
            K = (1.0 - alpha) * K_histories[tag](index - 1) + alpha * K_histories[tag](index);

            G[tag]  = (3.0 * K * E[tag]) / (9.0 * K - E[tag]);
            nu[tag] = (3.0 * K - E[tag]) / (6.0 * K);
        }

        else {
            E[tag] = E_histories[tag](time_histories[tag].Size() - 1) ;
            A[tag] = A_histories[tag](time_histories[tag].Size() - 1) ;
            K = K_histories[tag](time_histories[tag].Size() - 1) ;

            G[tag]  = (3.0 * K * E[tag]) / (9.0 * K - E[tag]);
            nu[tag] = (3.0 * K - E[tag]) / (6.0 * K);
        }

        opserr << " UPDATING PARAMETERS of mat = " << this->getTag() <<   " TO "  << endln;
        opserr << "  E  = " << E[tag]  << endln;
        opserr << "  G  = " << G[tag]  << endln;
        opserr << "  nu = " << nu[tag] << endln;
        opserr << "  A  = " << A[tag]  << endln;
    }
//     opserr << "current_time = " << time << " "
//            << "new_time_step = " << (int)new_time_step << " "
//            << "E = " << E << " "
//            << "G = " << G << " "
//            << "nu = " << nu << " "
//            << "A = " << A << " " << endln;
}
