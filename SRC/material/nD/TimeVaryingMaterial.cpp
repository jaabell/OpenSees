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

// $Revision: 1.4 $
// $Date: 2020-04-19 23:01:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/TimeVaryingMaterial.h,v $

// Davide Raino, Massimo Petracca - ASDEA Software, Italy
//
// A Generic Orthotropic Material Wrapper that can convert any
// nonlinear isotropic material into an orthotropic one by means of tensor
// mapping
//

#include <TimeVaryingMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <Parameter.h>


void *OPS_TimeVaryingMaterial(void)
{
    opserr << "Using TimeVaryingMaterial" << endln ;
    // check arguments
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 17) {
        opserr <<
               "nDMaterial TimeVarying Error: Few arguments (< 17).\n"
               "nDMaterial TimeVarying $tag $theIsoMat $Ex $Ey $Ez $Gxy $Gyz $Gzx $vxy $vyz $vzx $Asigmaxx $Asigmayy $Asigmazz $Asigmaxyxy $Asigmayzyz $Asigmaxzxz.\n";
        return nullptr;
    }

    // get integer data
    int iData[2];
    int numData = 2;
    if (OPS_GetInt(&numData, iData) != 0)  {
        opserr << "nDMaterial TimeVarying Error: invalid nDMaterial tags.\n";
        return nullptr;
    }

    // get double data
    double dData[15];
    numData = 15;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "nDMaterial TimeVarying Error: invalid data for nDMaterial TimeVarying with tah " << iData[0] << ".\n";
        return nullptr;
    }

    // get the isotropic material to map
    NDMaterial *theIsoMaterial = OPS_getNDMaterial(iData[1]);
    if (theIsoMaterial == 0) {
        opserr << "WARNING: nDMaterial does not exist.\n";
        opserr << "nDMaterial: " << iData[1] << "\n";
        opserr << "nDMaterial TimeVarying: " << iData[0] << "\n";
        return nullptr;
    }

    // Agregar lectura de curvas de tiempo para modulo de elasticidad y de resistencia


    // create the TimeVarying wrapper
    NDMaterial* theTimeVaryingMaterial = new TimeVaryingMaterial(
        iData[0],
        *theIsoMaterial,
        dData[0], dData[1], dData[2], dData[3],
        dData[4], dData[5], dData[6], dData[7],
        dData[8], dData[9], dData[10], dData[11],
        dData[12], dData[13], dData[14]);
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
    NDMaterial &theIsoMat,
    double Ex, double Ey, double Ez, double Gxy, double Gyz, double Gzx,
    double vxy, double vyz, double vzx,
    double Asigmaxx, double Asigmayy, double Asigmazz, double Asigmaxyxy, double Asigmayzyz, double Asigmaxzxz)
    : NDMaterial(tag, ND_TAG_TimeVaryingMaterial)
{
    // copy the isotropic material
    theIsotropicMaterial = theIsoMat.getCopy("ThreeDimensional");
    if (theIsotropicMaterial == 0) {
        opserr << "nDMaterial Orthotropic Error: failed to get a (3D) copy of the isotropic material\n";
        exit(-1);
    }

    opserr << "TimeVaryingMaterial paremeters:"<< endln ;
    opserr << "   Ex = "         << Ex         << endln ;
    opserr << "   Ey = "         << Ey         << endln ;
    opserr << "   Ez = "         << Ez         << endln ;
    opserr << "   Gxy = "        << Gxy        << endln ;
    opserr << "   Gyz = "        << Gyz        << endln ;
    opserr << "   Gzx = "        << Gzx        << endln ;
    opserr << "   vxy = "        << vxy        << endln ;
    opserr << "   vyz = "        << vyz        << endln ;
    opserr << "   vzx = "        << vzx        << endln ;
    opserr << "   Asigmaxx = "   << Asigmaxx   << endln ;
    opserr << "   Asigmayy = "   << Asigmayy   << endln ;
    opserr << "   Asigmazz = "   << Asigmazz   << endln ;
    opserr << "   Asigmaxyxy = " << Asigmaxyxy << endln ;
    opserr << "   Asigmayzyz = " << Asigmayzyz << endln ;
    opserr << "   Asigmaxzxz = " << Asigmaxzxz << endln ;

    parameters(0) = Ex;
    parameters(1) = Ey;
    parameters(2) = Ez;
    parameters(3) = Gxy;
    parameters(4) = Gyz;
    parameters(5) = Gzx;
    parameters(6) = vxy;
    parameters(7) = vyz;
    parameters(8) = vzx;
    parameters(9) = Asigmaxx;
    parameters(10) = Asigmayy;
    parameters(11) = Asigmazz;
    parameters(12) = Asigmaxyxy;
    parameters(13) = Asigmayzyz;
    parameters(14) = Asigmaxzxz;
}

TimeVaryingMaterial::~TimeVaryingMaterial()
{
    if (theIsotropicMaterial)
        delete theIsotropicMaterial;
}

double TimeVaryingMaterial::getRho(void)
{
    return theIsotropicMaterial->getRho();
}

int TimeVaryingMaterial::setTrialStrain(const Vector & strain)
{
    /* START --- THESE LINES ARE TO TEST FOR Δε

    static Vector depsilon(6);
    static Vector depsilon_internal(6);
    depsilon.Zero();
    depsilon_internal.Zero();
  
    depsilon = epsilon - epsilon_n;
    depsilon_internal = epsilon_internal - epsilon_internal_n;
   
    // commitState():
    //     epsilon_n = epsilon
    //     epsilon_internal_n = epsilon_internal
    // epsilon = (epsilon - epsilon_n) - (epsilon_internal - epsilon_internal_n)

    epsilon = depsilon - depsilon_internal ;
    END --- THESE LINES ARE TO TEST FOR Δε */ 

    // strain in orthotropic space (TOTAL STRAIN)
    epsilon = strain - epsilon_internal ;

    // Actualizar matrices de resistencia
    // Aepsilon
    // Asigma_inv
    // opserr << "ops_Dt = " << ops_Dt << endln; // ITS Δt NOT TOTAL CURRENT TIME

    // Lugar donde hay que variar los parametros elasticos
    // y de resistencia usando ops_Dt

    // Los de rigidez son rigideces "finales"
    // Los de resistencia son resistencia_inicial(t) / resistencia_actual(t) = A(t)
    // Asigmazz = 0.5 significa el doble de sigma_zz de resistencia

    double Ex         = parameters(0);  // E(t)
    double Ey         = parameters(1);  // E(t)
    double Ez         = parameters(2);  // E(t)
    double Gxy        = parameters(3);  // G(t) en funcion de Bulk(t) 
    double Gyz        = parameters(4);  // G(t) en funcion de Bulk(t) 
    double Gzx        = parameters(5);  // G(t) en funcion de Bulk(t) 
    double vxy        = parameters(6);  // nu(t) en funcion de Bulk(t) 
    double vyz        = parameters(7);  // nu(t) en funcion de Bulk(t) 
    double vzx        = parameters(8);  // nu(t) en funcion de Bulk(t) 
    double Asigmaxx   = parameters(9);  // A(t)
    double Asigmayy   = parameters(10); // A(t)
    double Asigmazz   = parameters(11); // A(t)
    double Asigmaxyxy = parameters(12); // A(t)
    double Asigmayzyz = parameters(13); // A(t)
    double Asigmaxzxz = parameters(14); // A(t)

    ///  Old code at constructor BEGIN ========================================================
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
    if (Asigmaxx <= 0 || Asigmayy <= 0 || Asigmazz <= 0 || Asigmaxyxy <= 0 || Asigmayzyz <= 0 || Asigmaxzxz <= 0) {
        opserr << "nDMaterial Orthotropic Error: Asigma11, Asigma22, Asigma33, Asigma12, Asigma23, Asigma13 must be greater than 0.\n";
        exit(-1);
    }
    static Matrix Asigma(6, 6);
    Asigma.Zero();
    Asigma(0, 0) = Asigmaxx;
    Asigma(1, 1) = Asigmayy;
    Asigma(2, 2) = Asigmazz;
    Asigma(3, 3) = Asigmaxyxy;
    Asigma(4, 4) = Asigmayzyz;
    Asigma(5, 5) = Asigmaxzxz;
    for (int i = 0; i < 6; ++i)
        Asigma_inv(i) = 1.0 / Asigma(i, i);

    // coompute the initial isotropic constitutive tensor and its inverse
    static Matrix C0iso(6, 6);
    static Matrix C0iso_inv(6, 6);
    C0iso = theIsotropicMaterial->getInitialTangent();
    int res = C0iso.Invert(C0iso_inv);
    if (res < 0) {
        opserr << "nDMaterial Orthotropic Error: the isotropic material gave a singular initial tangent.\n";
        exit(-1);
    }

    // compute the strain tensor map inv(C0_iso) * Asigma * C0_ortho
    static Matrix Asigma_C0(6, 6);
    Asigma_C0.addMatrixProduct(0.0, Asigma, C0, 1.0);
    Aepsilon.addMatrixProduct(0.0, C0iso_inv, Asigma_C0, 1.0);

    ///  Old code at constructor END ========================================================

    // move to isotropic space
    static Vector eps_iso(6);
    eps_iso.addMatrixVector(0.0, Aepsilon, epsilon, 1.0);

    // call isotropic material
    res = theIsotropicMaterial->setTrialStrain(eps_iso);
    if (res != 0) {
        opserr << "nDMaterial Orthotropic Error: the isotropic material failed in setTrialStrain.\n";
        return res;
    }
    return 0;
}

const Vector &TimeVaryingMaterial::getStrain(void)
{
    return epsilon;
}

const Vector &TimeVaryingMaterial::getStress(void)
{
    // stress in isotropic space
    const Vector& sigma_iso = theIsotropicMaterial->getStress();

    // move to orthotropic space
    static Vector sigma(6);
    for (int i = 0; i < 6; ++i)
        sigma(i) = Asigma_inv(i) * sigma_iso(i);
    return sigma;
}

const Matrix &TimeVaryingMaterial::getTangent(void)
{
    // tensor in isotropic space
    const Matrix &C_iso = theIsotropicMaterial->getTangent();

    // compute orthotripic tangent
    static Matrix C(6, 6);
    static Matrix temp(6, 6);
    static Matrix invAsigma(6, 6);
    invAsigma.Zero();
    for (int i = 0; i < 6; ++i)
        invAsigma(i, i) = Asigma_inv(i);
    temp.addMatrixProduct(0.0, C_iso, Aepsilon, 1.0);
    C.addMatrixProduct(0.0, invAsigma, temp, 1.0);
    return C;
}

const Matrix &TimeVaryingMaterial::getInitialTangent(void)
{
    // tensor in isotropic space
    const Matrix& C_iso = theIsotropicMaterial->getInitialTangent();

    // compute orthotripic tangent
    static Matrix C(6, 6);
    static Matrix temp(6, 6);
    static Matrix invAsigma(6, 6);
    invAsigma.Zero();
    for (int i = 0; i < 6; ++i)
        invAsigma(i, i) = Asigma_inv(i);
    temp.addMatrixProduct(0.0, C_iso, Aepsilon, 1.0);
    C.addMatrixProduct(0.0, invAsigma, temp, 1.0);
    return C;;
}

int TimeVaryingMaterial::commitState(void)
{
    epsilon_n=epsilon;
    epsilon_internal_n=epsilon_internal;
    return theIsotropicMaterial->commitState();
}

int TimeVaryingMaterial::revertToLastCommit(void)
{
    return theIsotropicMaterial->revertToLastCommit();
}

int TimeVaryingMaterial::revertToStart(void)
{
    return theIsotropicMaterial->revertToStart();
}

NDMaterial * TimeVaryingMaterial::getCopy(void)
{
    TimeVaryingMaterial *theCopy = new TimeVaryingMaterial();
    theCopy->setTag(getTag());
    theCopy->theIsotropicMaterial = theIsotropicMaterial->getCopy("ThreeDimensional");
    theCopy->epsilon = epsilon;
    theCopy->epsilon_n = epsilon_n;
    theCopy->Aepsilon = Aepsilon;
    theCopy->Asigma_inv = Asigma_inv;
    theCopy->parameters = parameters;
    theCopy->epsilon_internal = epsilon_internal;
    theCopy->epsilon_internal_n = epsilon_internal_n;
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
    // result
    int res = 0;

    // data
    static Vector data(48);
    int counter = 0;
    // store int values
    data(counter++) = static_cast<double>(getTag());
    data(counter++) = static_cast<double>(theIsotropicMaterial->getClassTag());
    int matDbTag = theIsotropicMaterial->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        theIsotropicMaterial->setDbTag(matDbTag);
    }
    data(counter++) = static_cast<double>(matDbTag);
    // store internal variables
    for (int i = 0; i < 6; ++i)
        data(counter++) = epsilon(i);
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            data(counter++) = Aepsilon(i, j);
    for (int i = 0; i < 6; ++i)
        data(counter++) = Asigma_inv(i);
    // send data
    res = theChannel.sendVector(getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "nDMaterial Orthotropic Error: failed to send vector data\n";
        return res;
    }

    // now send the materials data
    res = theIsotropicMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
        opserr << "nDMaterial Orthotropic Error: failed to send the isotropic material\n";
        return res;
    }

    // done
    return res;
}

int TimeVaryingMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
    // result
    int res = 0;

    // data
    static Vector data(48);
    int counter = 0;

    // receive data
    res = theChannel.recvVector(getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "nDMaterial Orthotropic Error: failed to send vector data\n";
        return res;
    }
    // get int values
    setTag(static_cast<int>(data(counter++)));
    int matClassTag = static_cast<int>(data(counter++));
    // if the associated material has not yet been created or is of the wrong type
    // create a new material for recvSelf later
    if ((theIsotropicMaterial == nullptr) || (theIsotropicMaterial->getClassTag() != matClassTag)) {
        if (theIsotropicMaterial)
            delete theIsotropicMaterial;
        theIsotropicMaterial = theBroker.getNewNDMaterial(matClassTag);
        if (theIsotropicMaterial == nullptr) {
            opserr << "nDMaterial Orthotropic Error: failed to get a material of type: " << matClassTag << endln;
            return -1;
        }
    }
    theIsotropicMaterial->setDbTag(static_cast<int>(data(counter++)));
    // store internal variables
    for (int i = 0; i < 6; ++i)
        epsilon(i) = data(counter++);
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            Aepsilon(i, j) = data(counter++);
    for (int i = 0; i < 6; ++i)
        Asigma_inv(i) = data(counter++);

    // now receive the associated materials data
    res = theIsotropicMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
        opserr << "nDMaterial Orthotropic Error: failed to receive the isotropic material\n";
        return res;
    }

    // done
    return res;
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
    return theIsotropicMaterial->setParameter(argv, argc, param);
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
            return theIsotropicMaterial->setResponse(argv, argc, s);
        }
    }
    return NDMaterial::setResponse(argv, argc, s);
}
