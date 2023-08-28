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

#ifndef TimeVaryingMaterial_h
#define TimeVaryingMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include <map>

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

    void getParameters(double);// Determines the current properties for E, G, nu and A for this material

    // the projected material pointer
    NDMaterial *theProjectedMaterial = nullptr;
    
    //Strain transformation map
    Matrix Aepsilon = Matrix(6, 6);

    //History of evolution
    // static Vector* time_history;
    // static Vector* E_history;
    // static Vector* K_history;
    // static Vector* A_history;
    static std::map<int, Vector> time_histories;
    static std::map<int, Vector> E_histories;
    static std::map<int, Vector> K_histories;
    static std::map<int, Vector> A_histories;

    //The "tag" of this material to the evolution laws (unused)
    // int evolution_law_id;
    // static int number_of_evolution_laws;

    //Strain offset due to thermal action
    Vector epsilon_internal = Vector(6);
    
    // State for the incremental model. 
    // _proj are the tensors for the projected material
    // _real are the tensors for the final (real) material 
    Vector sigma_real = Vector(6);    
    Vector sigma_proj = Vector(6);
    Vector epsilon_real = Vector(6);  
    Vector epsilon_proj = Vector(6);  
    Vector epsilon_new = Vector(6);  

    // State for the incremental model
    Vector sigma_real_n = Vector(6);    
    Vector sigma_proj_n = Vector(6);
    Vector epsilon_real_n = Vector(6);  
    Vector epsilon_proj_n = Vector(6);  
    Vector epsilon_new_n = Vector(6);  

    //global variables for all materials... should not be
    // static double E, G, nu, A;
    static std::map<int, double> E, G, nu, A;
    static std::map<int, bool> new_time_step;

    static bool print_strain_once;
    static bool print_stress_once;
    static bool print_commit_once;
    static bool print_tang_once;
};
#endif
