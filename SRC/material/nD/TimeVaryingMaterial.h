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
#include <limits>
#include <memory>

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
    // driver of this material point: the local evolution variable once set,
    // else the time of the domain being updated (ops_TheActiveDomain, so that
    // OpenSeesSP partitions use their subdomain's time), else OPS_GetDomain()
    double currentDriver(void) const;
    // evaluate the properties at the driver and build Aepsilon; same
    // arithmetic as before, shared by setTrialStrain and getInitialTangent
    int computeMap(double driver);

    // the projected material pointer
    NDMaterial *theProjectedMaterial = nullptr;

    //Strain transformation map
    Matrix Aepsilon = Matrix(6, 6);

    //History of evolution: immutable, owned per instance and shared by copies
    // (previously static maps keyed by tag)
    struct History {
        Vector t, E, K, A;
    };
    std::shared_ptr<const History> hist;

    //Strain offset due to thermal action
    Vector epsilon_internal = Vector(6);

    // State for the incremental model.
    // _proj are the tensors for the projected material
    // _real are the tensors for the final (real) material
    Vector sigma_real = Vector(6);
    Vector epsilon_real = Vector(6);
    Vector epsilon_new = Vector(6);

    // State for the incremental model
    Vector sigma_real_n = Vector(6);
    Vector sigma_proj_n = Vector(6);
    Vector epsilon_real_n = Vector(6);
    Vector epsilon_proj_n = Vector(6);
    Vector epsilon_new_n = Vector(6);

    //global variables for all materials... should not be
    // Current interpolated properties and refresh state belong to this
    // material-point copy. Evolution histories remain shared and immutable.
    double current_E = 0.0;
    double current_G = 0.0;
    double current_nu = 0.0;
    double current_A = 0.0;
    bool new_time_step = true;

    // Driver value the cached properties were evaluated at. NaN until the first
    // evaluation, so the first comparison always refreshes (audit F14).
    double last_parameter_driver = std::numeric_limits<double>::quiet_NaN();

    // Optional material-point evolution variable. Domain time remains the
    // driver until parameter 4002 (evolutionVariable) is updated.
    double local_evolution_variable = 0.0;
    bool use_local_evolution_variable = false;

    // committed copies, so that revertToLastCommit restores the mapped
    // response immediately (the wrapped material's own revert may be a no-op,
    // e.g. J2Plasticity, and its stress/tangent stay at the last trial)
    double committed_E = 0.0;
    double committed_G = 0.0;
    double committed_nu = 0.0;
    double committed_A = 0.0;
    Matrix committed_Aepsilon = Matrix(6, 6);
    Matrix committed_tangent = Matrix(6, 6);
    bool trialSinceCommit = false;  // a successful trial since the last commit
    bool reverted = false;          // queries answered from the committed copies
    bool mapReady = false;          // Aepsilon has been built at least once
    int projectedDbTag = 0;
};
#endif
