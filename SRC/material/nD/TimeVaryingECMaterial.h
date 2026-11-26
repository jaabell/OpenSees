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
//
// TimeVaryingEC: EXPERIMENTAL, opt-in wrapper that maps a reference
// small-strain 3D material with constant elastic tensor C0 onto a physical
// isotropic material whose stiffness C(xi) and strength map A(xi) are
// prescribed as functions of a driver xi (domain time, or a local value set
// through the parameter "evolutionVariable"). It is a separate class from
// TimeVarying: TimeVarying maps increments only, so a property change at fixed
// mechanical strain leaves the accepted stress unchanged; TimeVaryingEC adds
// the state-evolution term ("EC" correction):
//
//   eps_m  = eps - eps*                      (eps*: eigenstrain, "initNormalStrain")
//   S      = A C0^-1 C ,    eta_n = C_n^-1 sigma_n
//   epsP_{n+1} = epsP_n + S_{n+1}(eps_m,{n+1} - eps_m,n) + (S_{n+1} - S_n) eta_n
//   sig_{n+1}  = sig_n + [sigP_{n+1} - sigP_n - C0 (S_{n+1} - S_n) eta_n] / A_{n+1}
//                      + (C_{n+1} - C_n) eta_n
//   C_alg  = A_{n+1}^-1 Ctilde_alg S_{n+1}
//
// with eta_n, C_n, A_n, S_n always from the COMMITTED state n. The wrapper
// reads and writes no internal variable of the wrapped material.
//
// VALIDITY (first delivery): eta_n = C_n^-1 sigma_n is the recoverable strain
// only if the current unloading stiffness of the wrapped material, mapped,
// equals C_n. This holds for isotropic linear elasticity and for associated
// PERFECT J2 plasticity (configure J2Plasticity with sigInf = sig0 and
// delta = H = eta = 0). It is FALSE for stiffness-degrading materials: for a
// scalar damage model with unloading (1-d)C the correction is off by 1/(1-d)
// and under-corrects (measured 0.811 / 0.658 / 0.472 of the current-state
// target for d = 0.25 / 0.50 / 0.90). General hardening is also out of scope:
// stress identity does not imply equal hardening memory. Damage,
// ASDConcrete3D, hardening, creep, large strain and monolithic
// thermo-mechanical coupling are NOT covered. The wrapper cannot detect these
// features of the wrapped material; the user must not combine them.
//
// Further conventions and limits:
//  - Voigt order xx yy zz xy yz zx; engineering shear strains, tensor shear
//    stresses (the OpenSees convention of the 3D materials and elements).
//  - Initial state: eps_m,0 = 0 and sigma_0 = 0 at the first driver value; the
//    first eigenstrain change counts as an increment (as in TimeVarying). A
//    model that sets a nonzero initial stress (InitialStateAnalysis,
//    InitStressNDMaterial, prestress) needs eta_0 seeded consistently and is
//    out of the scope of this class.
//  - Driver fixed during the Newton iterations of a step: C_alg is the
//    mechanical derivative only; the derivative with respect to the driver is
//    not included.
//  - getStrain() returns the total input strain eps. Responses
//    "mechanicalStrain", "eigenStrain", "projectedStrain", "projectedStress"
//    and "properties" (K G A driver) are separate.
//  - State contract required from the wrapped material: setTrialStrain
//    integrates from its committed state, commitState advances it, getCopy
//    is deep. Its revertToLastCommit is called but not relied upon: every
//    query after a revert is answered from the copies kept by this wrapper.
//  - updateParameter (driver, eigenstrain) re-evaluates the trial at the
//    current strain, so the first unbalance of the next step already sees the
//    new input (otherwise a displacement-increment test can pass a step that
//    only changes the input, with the domain out of equilibrium). With the
//    domain-time driver the change arrives through setTrialStrain as usual.
//  - Domain-time driver: the time of the domain being updated
//    (ops_TheActiveDomain, else OPS_GetDomain()), so OpenSeesSP partitions
//    use their subdomain's time.
//  - getInitialTangent() returns the physical elastic stiffness C of the
//    current state; before the first trial (or after revertToStart) the law
//    is evaluated at the current driver value without changing the state.

#ifndef TimeVaryingECMaterial_h
#define TimeVaryingECMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <memory>

class TimeVaryingECMaterial : public NDMaterial
{
  public:
    // property law: E, K, A at strictly increasing driver values, linear
    // interpolation, constant extrapolation (same as TimeVarying)
    struct Law {
        Vector t, E, K, A;
    };

    TimeVaryingECMaterial(int tag, NDMaterial &theReferenceMaterial,
                          const Vector &t, const Vector &E,
                          const Vector &K, const Vector &A);
    TimeVaryingECMaterial();
    ~TimeVaryingECMaterial();

    const char *getClassType(void) const { return "TimeVaryingECMaterial"; }

    double getRho(void);

    int setTrialStrain(const Vector &strain);
    int setTrialStrain(const Vector &strain, const Vector &rate);
    const Vector &getStrain(void);
    const Vector &getStress(void);
    const Matrix &getTangent(void);
    const Matrix &getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *code);
    const char *getType(void) const;
    int getOrder(void) const;
    void Print(OPS_Stream &s, int flag = 0);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &matInfo);

  private:
    // complete state of one material point (trial or committed)
    struct State {
        State();
        Vector eps;      // total input strain
        Vector epsStar;  // eigenstrain
        Vector epsM;     // mechanical strain eps - eps*
        double driver;
        double K, G, A;  // physical properties at driver
        Vector epsP;     // projected (reference) strain
        Vector sigP;     // reference stress
        Vector sig;      // physical stress
        Matrix Calg;     // physical algorithmic tangent
        Vector eta;      // C^-1 sig (meaningful in the committed state)
    };

    double currentDriver(void) const;
    int evaluateLaw(double driver, double &K, double &G, double &A) const;
    void physicalStiffness(double K, double G, Matrix &C) const;
    void mapS(double K, double G, double A, Matrix &S) const;
    void recoverableStrain(State &s) const;
    void resetToStart(void);
    static int packState(const State &s, Vector &data, int pos);
    static int unpackState(State &s, const Vector &data, int pos);

    NDMaterial *theReference;           // wrapped reference material (owned)
    std::shared_ptr<const Law> law;     // immutable, shared between copies
    Matrix C0;                          // reference elastic tensor
    Matrix C0inv;

    State trial;
    State committed;
    bool trialSinceCommit;   // trial differs from committed
    bool referenceInSync;    // reference trial corresponds to `trial`
    bool pristine;           // no trial committed since construction / revertToStart

    // external inputs (pending, and as committed)
    Vector inputEpsStar, committedEpsStar;
    double inputDriver, committedDriver;
    bool useLocalDriver, committedUseLocalDriver;

    int referenceDbTag;

    static Vector workVector;
    static Matrix workMatrix;
};

#endif
