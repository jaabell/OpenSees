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

// Davide Raino, Massimo Petracca - ASDEA Software, Italy
//
// A Generic Orthotropic Material Wrapper that can convert any
// nonlinear isotropic material into an orthotropic one by means of tensor
// mapping
//

#ifndef TimeVaryingMaterial_h
#define TimeVaryingMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

class TimeVaryingMaterial : public NDMaterial 
{
public:
    // life-cycle
    TimeVaryingMaterial(
        int tag, 
        NDMaterial &theIsoMat,
        const Vector& time_history_in,
        const Vector& E_history_in,
        const Vector& K_history_in,
        const Vector& A_history_in);
    TimeVaryingMaterial();
    ~TimeVaryingMaterial();

    // info
    const char* getClassType(void) const { return "TimeVaryingMaterial"; };

    // density
    double getRho(void);

    // set state
    int setTrialStrain(const Vector &strain);

    // get state
    const Vector &getStrain(void);
    const Vector &getStress(void);
    const Matrix &getTangent(void);
    const Matrix &getInitialTangent(void);

    // handle state
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // copy and others...
    NDMaterial *getCopy(void);
    NDMaterial* getCopy(const char* code);
    const char *getType(void) const;
    int getOrder(void) const;
    void Print(OPS_Stream &s, int flag=0);

    // send/recv self
    virtual int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // parameters and responses
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

    Response* setResponse(const char** argv, int argc, OPS_Stream& s);



private:

    // double E(double time) const ;
    // double G(double time) const ;
    // double nu(double time) const ;
    // double A(double time) const ;
    void getParameters(double, double& E, double&G,double& nu,double& A);





    // the mapped isotropic material
    NDMaterial *theIsotropicMaterial = nullptr;
    // the strain in the real orthotropic space
    // Vector epsilon = Vector(6);
    // Vector epsilon_n = Vector(6);
    // strain tensor map
    Matrix Aepsilon = Matrix(6, 6);
    // inverse of stress tensor map (saved as vector, it's diagonal)
    Vector Asigma_inv = Vector(6);

    // double parameters[15];
    // Vector parameters = Vector(15);

    static Vector* time_history;
    static Vector* E_history;
    static Vector* K_history;
    static Vector* A_history;

    int evolution_law_id;

    Vector epsilon_internal = Vector(6);
    
    Vector sigma_real = Vector(6);    
    Vector sigma_proj = Vector(6);
    Vector epsilon_real = Vector(6);  
    Vector epsilon_proj = Vector(6);  

    Vector sigma_real_n = Vector(6);    // T
    Vector sigma_proj_n = Vector(6);
    Vector epsilon_real_n = Vector(6);  // Equivalent to epsilon_n
    Vector epsilon_proj_n = Vector(6);  //

    static bool new_time_step;
};
#endif