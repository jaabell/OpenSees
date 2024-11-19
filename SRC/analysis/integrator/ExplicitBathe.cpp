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
#define OPS_Export

void *OPS_ExplicitBathe(void) {
    // Pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = nullptr;

    // Check if enough arguments are provided
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "WARNING: Insufficient arguments for ExplicitBathe integrator. Expected 1 arguments (p, q0, q1, q2, s).\n";
        return nullptr;
    }

    // Read input parameters
    double data[1];
    if (OPS_GetDoubleInput(&numArgs, data) < 0) {
        opserr << "WARNING: Invalid input for ExplicitBathe integrator.\n";
        return nullptr;
    }

    // Create the ExplicitBathe integrator with the provided parameters
    theIntegrator = new ExplicitBathe(data[0]);//, data[1], data[2], data[3], data[4]);

    if (theIntegrator == nullptr) {
        opserr << "WARNING - out of memory creating ExplicitBathe integrator\n";
    }
    return theIntegrator;
}

ExplicitBathe::ExplicitBathe()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(0.0), q0(0.0), q1(0.0), q2(0.0), s(0.0),
      Utm1(0), Ut(0), Utdot(0), Utdotdot(0),
      Udot(0), Udotdot(0)
{}

ExplicitBathe::ExplicitBathe(double _p)//);, double _q0, double _q1, double _q2, double _s)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(_p),// q0(_q0), q1(_q1), q2(_q2), s(_s),
      Utm1(0), Ut(0), Utdot(0), Utdotdot(0),
      Udot(0), Udotdot(0)
{
    q1 = (1 - 2*p)/(2*p*(1-p));
    q2 = 0.5 - p*q1;
    q0 = -q1 -q2 + 0.5;
    s = -1;
}

ExplicitBathe::~ExplicitBathe() {
    if (Utm1) delete Utm1;
    if (Ut) delete Ut;
    if (Utdot) delete Utdot;
    if (Utdotdot) delete Utdotdot;
    if (Udot) delete Udot;
    if (Udotdot) delete Udotdot;
}

int ExplicitBathe::newStep(double _deltaT) {
    deltaT = _deltaT;
    double subStep1 = p * deltaT;
    double subStep2 = (1 - p) * deltaT;

    // Ensure state variables are initialized
    if (!Ut || !Utdot || !Utdotdot) {
        opserr << "ExplicitBathe::newStep() - state variables not initialized\n";
        return -1;
    }

    // Sub-step 1: Calculate intermediate displacement and velocity
    Vector U_half = *Ut + subStep1 * *Utdot + 0.5 * (subStep1 * subStep1) * *Utdotdot;
    Vector V_half = *Utdot + subStep1 * *Utdotdot;

    // Adjust for damping (Eq. 10)
    Vector U_hat = (1 - s) * *Ut + s * U_half;
    Vector V_hat = (1 - s) * *Utdot + s * V_half;

    // Sub-step 2: Calculate final displacement, velocity, and acceleration (Eq. 6-9)
    Vector U_final = U_hat + subStep2 * V_hat +
                     0.5 * (subStep2 * subStep2) * (*Utdotdot);

    Vector V_final = V_hat + subStep2 * *Utdotdot;

    // Update response variables
    *Utm1 = *Ut;  // Previous displacement
    *Ut = U_final;
    *Utdot = V_final;

    // Update acceleration (Eq. 9)
    *Utdotdot = q0 * (*Udot) + q1 * (*Utdotdot) + q2 * (*Udotdot);

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
    if (!Ut || Ut->Size() != size) {
        if (Utm1) delete Utm1;
        if (Ut) delete Ut;
        if (Utdot) delete Utdot;
        if (Utdotdot) delete Utdotdot;

        Utm1 = new Vector(size);
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);

        if (!Utm1 || !Ut || !Utdot || !Utdotdot) {
            opserr << "ExplicitBathe::domainChanged - out of memory\n";
            return -1;
        }
    }

    // Initialize state variables
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != nullptr) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*Ut)(loc) = disp(i);
                (*Utm1)(loc) = disp(i);  // Assume Ut-1 = Ut initially
            }
        }

        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*Utdot)(loc) = vel(i);
            }
        }

        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*Utdotdot)(loc) = accel(i);
            }
        }
    }

    return 0;
}

int ExplicitBathe::update(const Vector &U) {
    if (!Ut) {
        opserr << "ExplicitBathe::update() - state variables not initialized\n";
        return -1;
    }

    // Update displacement, velocity, and acceleration for this step
    *Ut = U;

    // Calculate velocity update
    double subStep1 = p * deltaT;
    double subStep2 = (1 - p) * deltaT;

    Vector V_updated = *Utdot + subStep1 * (*Utdotdot);
    *Utdot = V_updated;

    // Update the acceleration based on the new displacement
    Vector A_updated = q0 * (*Udot) + q1 * (*Utdotdot) + q2 * (*Udotdot);
    *Utdotdot = A_updated;

    return 0;
}

int ExplicitBathe::formEleTangent(FE_Element *theEle) {
    theEle->zeroTangent();

    double subStep1 = p * deltaT;
    double subStep2 = (1 - p) * deltaT;

    // Form the tangent matrix contributions based on sub-step values
    theEle->addMtoTang(subStep1);
    theEle->addCtoTang(subStep2);

    return 0;
}

int ExplicitBathe::formNodTangent(DOF_Group *theDof) {
    theDof->zeroTangent();

    double subStep1 = p * deltaT;
    double subStep2 = (1 - p) * deltaT;

    // Form the tangent contributions for nodal degrees of freedom
    theDof->addMtoTang(subStep1);
    theDof->addCtoTang(subStep2);

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
    return *Utdot;
}

int ExplicitBathe::sendSelf(int cTag, Channel &theChannel) {
    Vector data(5);
    data(0) = p;
    data(1) = q0;
    data(2) = q1;
    data(3) = q2;
    data(4) = s;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}

int ExplicitBathe::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    Vector data(5);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::recvSelf() - could not receive data\n";
        return -1;
    }

    p = data(0);
    q0 = data(1);
    q1 = data(2);
    q2 = data(3);
    s = data(4);

    return 0;
}

void ExplicitBathe::Print(OPS_Stream &stream, int flag) {
    stream << "Explicit Bathe Method:\n";
    stream << "  Time Step: " << deltaT << "\n";
    stream << "  p: " << p << ", q0: " << q0 << ", q1: " << q1 << ", q2: " << q2 << ", s: " << s << "\n";
}
