#include <ExplicitBathe.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <math.h>
#define OPS_Export

#include <NodeIter.h>
#include <LoadPatternIter.h>

void *OPS_ExplicitBathe(void) {
    // Pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = nullptr;

    // Check if enough arguments are provided
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "WARNING: Insufficient arguments for ExplicitBathe integrator. Expected 1 argument (p).\n";
        return nullptr;
    }

    // Read input parameters
    double p;
    if (OPS_GetDoubleInput(&numArgs, &p) < 0) {
        opserr << "WARNING: Invalid input for ExplicitBathe integrator.\n";
        return nullptr;
    }

    // Create the ExplicitBathe integrator with the provided parameters
    theIntegrator = new ExplicitBathe(p);

    if (theIntegrator == nullptr) {
        opserr << "WARNING - out of memory creating ExplicitBathe integrator\n";
    }
    return theIntegrator;
}

ExplicitBathe::ExplicitBathe()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(0.0), q0(0.0), q1(0.0), q2(0.0), s(0.0),
      U_t(0),V_t(0),A_t(0),R_t(0),
      U_tpdt(0),V_tpdt(0),A_tpdt(0),Rhat_tpdt(0),Rfunnyhat_tpdt(0),
      U_tdt(0),V_tdt(0),A_tdt(0),R_tdt(0), updateCount(0),
      a0(0.),a1(0.),a2(0.),a3(0.),a4(0.),a5(0.),a6(0.),a7(0.)
{}

ExplicitBathe::ExplicitBathe(double _p)//);, double _q0, double _q1, double _q2, double _s)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(_p),
      U_t(0),V_t(0),A_t(0),R_t(0),
      U_tpdt(0),V_tpdt(0),A_tpdt(0),Rhat_tpdt(0),Rfunnyhat_tpdt(0),
      U_tdt(0),V_tdt(0),A_tdt(0),R_tdt(0), updateCount(0),
      a0(0.),a1(0.),a2(0.),a3(0.),a4(0.),a5(0.),a6(0.),a7(0.)
{
    // Calculate the integration constants based on p
    q1 = (1 - 2*p)/(2*p*(1-p));
    q2 = 0.5 - p * q1;
    q0 = -q1 -q2 + 0.5;
    s = -1;
}

ExplicitBathe::~ExplicitBathe() {
    if (U_t) delete U_t;
    if (V_t) delete V_t;
    if (A_t) delete A_t;
    if (R_t) delete R_t;
    if (U_tpdt) delete U_tpdt;
    if (V_tpdt) delete V_tpdt;
    if (A_tpdt) delete A_tpdt;
    if (Rhat_tpdt) delete Rhat_tpdt;
    if (Rfunnyhat_tpdt) delete Rfunnyhat_tpdt;
    if (U_tdt) delete U_tdt;
    if (V_tdt) delete V_tdt;
    if (A_tdt) delete A_tdt;
    if (R_tdt) delete R_tdt;
}

int ExplicitBathe::newStep(double _deltaT) {
    deltaT = _deltaT;

    // A. Initial Calculations
    a0 = p * deltaT;
    a1 = 0.5 * pow(p * deltaT,2);
    a2 = 0.5 * a0;
    a3 = (1 - p) * deltaT;
    a4 = 0.5 * pow((1 - p) * deltaT,2);
    a5 = q0 * a3;
    a6 = (0.5 + q1) * a3;
    a7 = q2 * a3;

    // Ensure state variables are initialized
    // if (!Ut || !Utdot || !Utdotdot) {
    //     opserr << "ExplicitBathe::newStep() - state variables not initialized\n";
    //     return -1;
    // }

    // // Sub-step 1: Calculate intermediate displacement and velocity
    // // 1. Calculate displacements and effective loads at time t + pΔt
    // (*U_pdt) = *(U_t) + a0 * (*V_t) + a1 * (*A_t);
    // // Vector Rhat_pdt = *Ut * a5 + *Udotdot * a6; // Effective loads based on the given constants


    // LinearSOE *theLinSOE = this->getLinearSOE();
    // if (theLinSOE->setB(R_pdt) < 0) {
    //     opserr << "ExplicitBathe::newStep() - failed to set loads\n";
    //     return -1;
    // }
    // if (theLinSOE->solve() < 0) {
    //     opserr << "ExplicitBathe::newStep() - failed to solve for accelerations\n";
    //     return -1;
    // }
    // *Udotdot = theLinSOE->getX();

    // // 3. Calculate velocities at time t + pΔt
    // Vector V_pdt = *Utdot + a2 * (*Udotdot);

    // // Update intermediate state variables
    // *Utm1 = *Ut;
    // *Ut = U_pdt;
    // *Utdot = V_pdt;

    // // C. Second Sub-Step
    // // 1. Calculate displacements and effective loads at time t + Δt
    // Vector U_dt = *Ut + a3 * (*Utdot) + a4 * (*Udotdot);
    // Vector R_dt = *Ut * a5 + *Udotdot * a6 + *Utdot * a7;

    // // Solve for accelerations at t + Δt
    // if (theLinSOE->setB(R_dt) < 0) {
    //     opserr << "ExplicitBathe::newStep() - failed to set loads for second sub-step\n";
    //     return -1;
    // }
    // if (theLinSOE->solve() < 0) {
    //     opserr << "ExplicitBathe::newStep() - failed to solve for accelerations in second sub-step\n";
    //     return -1;
    // }
    // *Udotdot = theLinSOE->getX();

    // // 3. Calculate velocities at time t + Δt
    // Vector V_dt = *Utdot + a5 * (*Udotdot) + a6 * (*Udotdot);

    // // Update state variables for the final step
    // *Ut = U_dt;
    // *Utdot = V_dt;

    return 0;
}

int ExplicitBathe::domainChanged() {
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (!theModel || !theLinSOE) {
        opserr << "ExplicitBathe::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();

    if (size == 0) {
        opserr << "ExplicitBathe::domainChanged - invalid size\n";
        return -1;
    }

    // Allocate memory for state variables
    // if (!Ut || Ut->Size() != size) {
    //     if (U_t) delete U_t;
    //     if (V_t) delete V_t;
    //     if (A_t) delete A_t;
    //     if (U_tpdt) delete U_tpdt;
    //     if (V_tpdt) delete V_tpdt;
    //     if (A_tpdt) delete A_tpdt;
    //     if (Rhat_tpdt) delete Rhat_tpdt;
    //     if (Rfunnyhat_tpdt) delete Rfunnyhat_tpdt;
    //     if (U_tdt) delete U_tdt;
    //     if (V_tdt) delete V_tdt;
    //     if (A_tdt) delete A_tdt;

    //     U_t = new Vector(size);
    //     V_t = new Vector(size);
    //     A_t = new Vector(size);
    //     U_tpdt = new Vector(size);
    //     V_tpdt = new Vector(size);
    //     A_tpdt = new Vector(size);
    //     Rhat_tpdt = new Vector(size);
    //     Rfunnyhat_tpdt = new Vector(size);
    //     U_tdt = new Vector(size);
    //     V_tdt = new Vector(size);
    //     A_tdt = new Vector(size);

    //     if (!U_t || 
    //         !V_t || 
    //         !A_t || 
    //         !U_tpdt || 
    //         !V_tpdt || 
    //         !A_tpdt || 
    //         !Rhat_tpdt || 
    //         !Rfunnyhat_tpdt || 
    //         !U_tdt || 
    //         !V_tdt || 
    //         !A_tdt || 
    //         ) {
    //         opserr << "ExplicitBathe::domainChanged - out of memory\n";
    //         return -1;
    //     }
    // }

    // // Initialize state variables
    // DOF_GrpIter &theDOFs = theModel->getDOFs();
    // DOF_Group *dofPtr;
    // while ((dofPtr = theDOFs()) != nullptr) {
    //     const ID &id = dofPtr->getID();
    //     int idSize = id.Size();

    //     const Vector &disp = dofPtr->getCommittedDisp();
    //     for (int i = 0; i < idSize; ++i) {
    //         int loc = id(i);
    //         if (loc >= 0) {
    //             (*U_t)(loc) = disp(i);
    //             (*U_tpdt)(loc) = disp(i);  // Assume Ut + pdt = Ut initially
    //             (*U_tdt)(loc) = disp(i);  // Assume  Ut +  dt = Ut initially
    //         }
    //     }

    //     const Vector &vel = dofPtr->getCommittedVel();
    //     for (int i = 0; i < idSize; ++i) {
    //         int loc = id(i);
    //         if (loc >= 0) {
    //             (*V_t)(loc) = vel(i);
    //         }
    //     }

    //     const Vector &accel = dofPtr->getCommittedAccel();
    //     for (int i = 0; i < idSize; ++i) {
    //         int loc = id(i);
    //         if (loc >= 0) {
    //             (*A_t)(loc) = accel(i);
    //         }
    //     }
    // }

    return 0;
}




int 
ExplicitBathe::formRHS_tpdt()
{

    // Get pointer to the SOE
    LinearSOE *theSOE = this->getLinearSOE();

    theSOE->zeroB();


    // Code copied from Newmark...


    // // Get the analysis model
    // AnalysisModel *theModel = this->getAnalysisModel();

    // // Randomness in external load (including randomness in time series)
    // // Get domain
    // Domain *theDomain = theModel->getDomainPtr();

    // // Loop through nodes to zero the unbalaced load
    // Node *nodePtr;
    // NodeIter &theNodeIter = theDomain->getNodes();
    // while ((nodePtr = theNodeIter()) != 0)
    // nodePtr->zeroUnbalancedLoad();

    // // Loop through load patterns to add external load sensitivity
    // LoadPattern *loadPatternPtr;
    // LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    // double time;
    // while((loadPatternPtr = thePatterns()) != 0) {
    // time = theDomain->getCurrentTime();
    // loadPatternPtr->applyLoadSensitivity(time);
    // }

    // // Randomness in element/material contributions
    // // Loop through FE elements
    // FE_Element *elePtr;
    // FE_EleIter &theEles = theModel->getFEs();    
    // while((elePtr = theEles()) != 0) {
    // theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
    // }

    // // Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
    // DOF_Group *dofPtr;
    // DOF_GrpIter &theDOFs = theModel->getDOFs();
    // while((dofPtr = theDOFs()) != 0) {
    // theSOE->addB(  dofPtr->getUnbalance(this),  dofPtr->getID()  );
    // }

    // // Reset the sensitivity flag
    // sensitivityFlag = 0;

    return 0;
}






int ExplicitBathe::update(const Vector &U) {
  
   updateCount++;
    if (updateCount > 2)  {
        opserr << "WARNING ExplicitBathe::update() - called more than once -";
        opserr << " ExplicitBathe integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING ExplicitBathe::update() - no souAnalysisModel set\n";
        return -2;
    }

    // check domainChanged() has been called, i.e. Ut will not be zero
    if (U_t == 0)  {
        opserr << "WARNING ExplicitBathe::update() - domainChange() failed or not called\n";
        return -3;
    }

    // check Udotdot is of correct size
    if (A_t->Size() != A_tdt->Size()) {
        opserr << "WARNING ExplicitBathe::update() - Vectors of incompatible size ";
        opserr << " expecting " << A_t->Size() << " obtained " << A_tdt->Size() << endln;
        return -4;
    }

    int size = A_t->Size();


    // determine the response at t  + p *deltaT
    double pDt = p * deltaT ;

    // Get displacement at t  + p *deltaT
    *U_tpdt = *U_t;
    U_tpdt -> addVector(1.0, *V_t, a0);
    U_tpdt -> addVector(1.0, *A_t, a1);


    this->formRHS_tpdt();   // Agregar aqui las ecuacione sdel RHS
    
    // Solve for displacement sensitivity
    
    // theSOE->solve(); // solve is private here

    // A_tpdt = theSOE->getX();   // getX is private here

    // Get velocity  at t  + p *deltaT
    *V_tpdt = *V_t;
    V_tpdt -> addVector(1.0, *A_t, a2);
    V_tpdt -> addVector(1.0, *A_tpdt, a2);

    // Get displacement at t  + deltaT
    *U_tpdt = *U_tpdt;
    U_tpdt -> addVector(1.0, *V_tpdt, a3);
    U_tpdt -> addVector(1.0, *A_tpdt, a4);

    theModel->setResponse(*U_tpdt, *V_t, *A_t);

    if (theModel->updateDomain() < 0)  {
        opserr << "ExplicitBathe::update() - failed to update the domain\n";
        return -5;
    }
    // set response at t to be that at t+deltaT of previous step

    // (*Utdotdot) = Udotdot;
    // (*Utdotdot1) = Udotdot;

    return 0;
}

int ExplicitBathe::formEleTangent(FE_Element *theEle) {

    theEle->zeroTangent();
    theEle->addMtoTang();

    return 0;
}

int ExplicitBathe::formNodTangent(DOF_Group *theDof) {
    theDof->zeroTangent();
    theDof->addMtoTang();

    return 0;
}

int ExplicitBathe::commit() {
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == nullptr) {
        opserr << "ExplicitBathe::commit() - no AnalysisModel set\n";
        return -1;
    }

    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    theModel->setCurrentDomainTime(time);

    return theModel->commitDomain();
}

const Vector &ExplicitBathe::getVel() {
    // return *Utdot;
}

int ExplicitBathe::sendSelf(int cTag, Channel &theChannel) {
    Vector data(1);
    data(0) = p;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}

int ExplicitBathe::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    Vector data(1);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::recvSelf() - could not receive data\n";
        return -1;
    }

    p = data(0);

    // Recalculate integration constants based on received p
    q1 = (1 - 2*p)/(2*p*(1-p));
    q2 = 0.5 - p * q1;
    q0 = -q1 -q2 + 0.5;
    s = -1;

    return 0;
}

void ExplicitBathe::Print(OPS_Stream &stream, int flag) {
    stream << "Explicit Bathe Method:\n";
    stream << "  Time Step: " << deltaT << "\n";
    stream << "  p: " << p << ", q0: " << q0 << ", q1: " << q1 << ", q2: " << q2 << ", s: " << s << "\n";
}
