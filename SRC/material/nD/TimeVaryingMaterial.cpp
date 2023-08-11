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
#include <AnalysisModel.h>
#include <Parameter.h>


void *OPS_TimeVaryingMaterial(void)
{
    opserr << "Using TimeVaryingMaterial" << endln ;
    // check arguments
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 7) {
        opserr <<
               "nDMaterial TimeVarying Error: Few arguments (< 17).\n"
               "nDMaterial TimeVarying $tag $theIsoMat $Nt t1 t2 t3 .. t_Nt E1 E2 E3 ... E_Nt K1 K2 K3...K_Nt  A1 A2 A3...A_Nt\n";
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

// model basic -ndf 3 -ndm 3
// nDMaterial ElasticIsotropic 1 1 0.1
// nDMaterial TimeVarying 2 1 3 1 1 1 2 2 2 3 3 3 4 4 4

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
        t, E, K , A);
    if (theTimeVaryingMaterial == 0) {
        opserr << "nDMaterial TimeVarying Error: failed to allocate a new material.\n";
        return nullptr;
    }

    // done
    return theTimeVaryingMaterial;
}

Vector* TimeVaryingMaterial::time_history = 0;
Vector* TimeVaryingMaterial::E_history = 0;
Vector* TimeVaryingMaterial::K_history = 0;
Vector* TimeVaryingMaterial::A_history = 0;
bool TimeVaryingMaterial::new_time_step = false;

TimeVaryingMaterial::TimeVaryingMaterial()
    : NDMaterial(0, ND_TAG_TimeVaryingMaterial)
{
}

TimeVaryingMaterial::TimeVaryingMaterial(
    int tag,
    NDMaterial &theIsoMat,
    const Vector& t_, const Vector& E_, const Vector& K_, const Vector& A_)
    : NDMaterial(tag, ND_TAG_TimeVaryingMaterial)
{
    // copy the isotropic material
    theIsotropicMaterial = theIsoMat.getCopy("ThreeDimensional");
    if (theIsotropicMaterial == 0) {
        opserr << "nDMaterial Orthotropic Error: failed to get a (3D) copy of the isotropic material\n";
        exit(-1);
    }

    int Ndatapoints = E_.Size();
    time_history = new Vector(Ndatapoints);
    E_history = new Vector(Ndatapoints);
    K_history = new Vector(Ndatapoints);
    A_history = new Vector(Ndatapoints);

    *time_history = t_;
    *E_history = E_;
    *K_history = K_;
    *A_history = A_;

    opserr << "Created new TimeVaryingMaterial \n";
    opserr << "    tag = " << tag << endln;
    opserr << "    proj_tag = " << theIsoMat.getTag() << endln;
    opserr << "    Nt = " << E_history->Size() << endln;
    opserr << "    time_history = " << *time_history << endln;
    opserr << "    E_history = " << *E_history << endln;
    opserr << "    K_history = " << *K_history << endln;
    opserr << "    A_history = " << *A_history << endln;
}

TimeVaryingMaterial::~TimeVaryingMaterial()
{
    if (theIsotropicMaterial)
    {
    	delete time_history;
    	delete E_history;
    	delete K_history;
    	delete A_history;
        delete theIsotropicMaterial;
    }
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

    epsilon_real = strain;

    static Vector epsilon_use(6);
    epsilon_use.Zero();

    epsilon_use = strain - epsilon_internal ;

    // Actualizar matrices de resistencia
    // Aepsilon
    // Asigma_inv
    // opserr << "ops_Dt = " << ops_Dt << endln; // ITS Δt NOT TOTAL CURRENT TIME

    // Lugar donde hay que variar los parametros elasticos
    // y de resistencia usando ops_Dt

    // Los de rigidez son rigideces "finales"
    // Los de resistencia son resistencia_inicial(t) / resistencia_actual(t) = A(t)
    // Asigmazz = 0.5 significa el doble de sigma_zz de resistencia

    double current_time = (*OPS_GetAnalysisModel())->getCurrentDomainTime();

    opserr << "current_time = " << current_time << endln;

    double E, G, nu, A;
    getParameters(current_time, E, G, nu, A);

    double Ex         = E; //parameters(0);  // E(t)
    double Ey         = E; //parameters(1);  // E(t)
    double Ez         = E; //parameters(2);  // E(t)
    double Gxy        = G; //parameters(3);  // G(t) en funcion de Bulk(t) 
    double Gyz        = G; //parameters(4);  // G(t) en funcion de Bulk(t) 
    double Gzx        = G; //parameters(5);  // G(t) en funcion de Bulk(t) 
    double vxy        = nu; //parameters(6);  // nu(t) en funcion de Bulk(t) 
    double vyz        = nu; //parameters(7);  // nu(t) en funcion de Bulk(t) 
    double vzx        = nu; //parameters(8);  // nu(t) en funcion de Bulk(t) 
    double Asigmaxx   = A; //parameters(9);  // A(t)
    double Asigmayy   = A; //parameters(10); // A(t)
    double Asigmazz   = A; //parameters(11); // A(t)
    double Asigmaxyxy = A; //parameters(12); // A(t)
    double Asigmayzyz = A; //parameters(13); // A(t)
    double Asigmaxzxz = A; //parameters(14); // A(t)

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
    eps_iso.addMatrixVector(0.0, Aepsilon, epsilon_use, 1.0);

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
    return epsilon_real;
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
	sigma_real_n = sigma_real;
	sigma_proj_n = sigma_proj;
	epsilon_real_n = epsilon_real;
	epsilon_proj_n = epsilon_proj;
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
    // theCopy->epsilon = epsilon;
    // theCopy->epsilon_n = epsilon_n;
    theCopy->Aepsilon = Aepsilon;
    theCopy->Asigma_inv = Asigma_inv;
    // theCopy->parameters = parameters;
    theCopy->epsilon_internal = epsilon_internal;
    theCopy->sigma_real = sigma_real;
	theCopy->sigma_proj = sigma_proj;
	theCopy->epsilon_real = epsilon_real;
	theCopy->epsilon_proj = epsilon_proj;
	theCopy->sigma_real_n = sigma_real_n;
	theCopy->sigma_proj_n = sigma_proj_n;
	theCopy->epsilon_real_n = epsilon_real_n;
	theCopy->epsilon_proj_n = epsilon_proj_n;
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
    // int res = 0;

    // // data
    // static Vector data(48);
    // int counter = 0;
    // // store int values
    // data(counter++) = static_cast<double>(getTag());
    // data(counter++) = static_cast<double>(theIsotropicMaterial->getClassTag());
    // int matDbTag = theIsotropicMaterial->getDbTag();
    // if (matDbTag == 0) {
    //     matDbTag = theChannel.getDbTag();
    //     theIsotropicMaterial->setDbTag(matDbTag);
    // }
    // data(counter++) = static_cast<double>(matDbTag);
    // // store internal variables
    // for (int i = 0; i < 6; ++i)
    //     data(counter++) = epsilon(i);
    // for (int i = 0; i < 6; ++i)
    //     for (int j = 0; j < 6; ++j)
    //         data(counter++) = Aepsilon(i, j);
    // for (int i = 0; i < 6; ++i)
    //     data(counter++) = Asigma_inv(i);
    // // send data
    // res = theChannel.sendVector(getDbTag(), commitTag, data);
    // if (res < 0) {
    //     opserr << "nDMaterial Orthotropic Error: failed to send vector data\n";
    //     return res;
    // }

    // // now send the materials data
    // res = theIsotropicMaterial->sendSelf(commitTag, theChannel);
    // if (res < 0) {
    //     opserr << "nDMaterial Orthotropic Error: failed to send the isotropic material\n";
    //     return res;
    // }

    // // done
    // return res;
    return 0;
}

int TimeVaryingMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
    // // result
    // int res = 0;

    // // data
    // static Vector data(48);
    // int counter = 0;

    // // receive data
    // res = theChannel.recvVector(getDbTag(), commitTag, data);
    // if (res < 0) {
    //     opserr << "nDMaterial Orthotropic Error: failed to send vector data\n";
    //     return res;
    // }
    // // get int values
    // setTag(static_cast<int>(data(counter++)));
    // int matClassTag = static_cast<int>(data(counter++));
    // // if the associated material has not yet been created or is of the wrong type
    // // create a new material for recvSelf later
    // if ((theIsotropicMaterial == nullptr) || (theIsotropicMaterial->getClassTag() != matClassTag)) {
    //     if (theIsotropicMaterial)
    //         delete theIsotropicMaterial;
    //     theIsotropicMaterial = theBroker.getNewNDMaterial(matClassTag);
    //     if (theIsotropicMaterial == nullptr) {
    //         opserr << "nDMaterial Orthotropic Error: failed to get a material of type: " << matClassTag << endln;
    //         return -1;
    //     }
    // }
    // theIsotropicMaterial->setDbTag(static_cast<int>(data(counter++)));
    // // store internal variables
    // for (int i = 0; i < 6; ++i)
    //     epsilon(i) = data(counter++);
    // for (int i = 0; i < 6; ++i)
    //     for (int j = 0; j < 6; ++j)
    //         Aepsilon(i, j) = data(counter++);
    // for (int i = 0; i < 6; ++i)
    //     Asigma_inv(i) = data(counter++);

    // // now receive the associated materials data
    // res = theIsotropicMaterial->recvSelf(commitTag, theChannel, theBroker);
    // if (res < 0) {
    //     opserr << "nDMaterial Orthotropic Error: failed to receive the isotropic material\n";
    //     return res;
    // }

    // // done
    // return res;
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
    return theIsotropicMaterial->setParameter(argv, argc, param);
}


int TimeVaryingMaterial::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {

    case 4001:
    {
        double initNormalStrain = info.theDouble; // alpha * (Temp(t)  - Temp())
        opserr << "initNormalStrain = " << initNormalStrain << endln;
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


void TimeVaryingMaterial::getParameters(double time, double& E, double& G, double& nu, double& A)
{
	if(new_time_step)
	{
		//hacer calculos
		new_time_step = false;
		E = 0;   //interpolate
		G = 0;
		nu = 0;
		A = 0;
	}

	return;
}

