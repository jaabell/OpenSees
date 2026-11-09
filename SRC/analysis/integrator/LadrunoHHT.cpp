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

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell (Ladruño project)
//
// References:
//   H.M. Hilber, T.J.R. Hughes, R.L. Taylor (1977). Improved numerical
//     dissipation for time integration algorithms in structural dynamics.
//     Earthquake Engineering and Structural Dynamics 5(3):283-292.
//     doi:10.1002/eqe.4290050306
//   J.P. Conte, P.K. Vijalapura, M. Meghella (2003). Consistent finite-element
//     response sensitivity analysis. Journal of Engineering Mechanics
//     129(12):1380-1393. doi:10.1061/(ASCE)0733-9399(2003)129:12(1380)
//
// Description: HHT-alpha integrator with DDM response sensitivities. See
// LadrunoHHT.h for the derivation of the alpha-weighted sensitivity residual.

#include <LadrunoHHT.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <elementAPI.h>


void *
OPS_LadrunoHHT(void)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;

    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 1 && argc != 3) {
        opserr << "WARNING - incorrect number of args want LadrunoHHT $alpha <$gamma $beta>\n";
        return 0;
    }

    double dData[3];
    if (OPS_GetDouble(&argc, dData) != 0) {
        opserr << "WARNING - invalid args want LadrunoHHT $alpha <$gamma $beta>\n";
        return 0;
    }

    if (argc == 1)
        theIntegrator = new LadrunoHHT(dData[0]);
    else
        theIntegrator = new LadrunoHHT(dData[0], dData[2], dData[1]); // alpha, beta, gamma for the constructor

    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating LadrunoHHT integrator\n";

    return theIntegrator;
}


LadrunoHHT::LadrunoHHT()
    : HHT(INTEGRATOR_TAGS_LadrunoHHT, 1.0, 0.25, 0.5),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS()
{

}


LadrunoHHT::LadrunoHHT(double _alpha)
    // same beta and gamma that HHT derives from a single alpha
    : HHT(INTEGRATOR_TAGS_LadrunoHHT, _alpha, (2-_alpha)*(2-_alpha)*0.25, 1.5-_alpha),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS()
{

}


LadrunoHHT::LadrunoHHT(double _alpha, double _beta, double _gamma)
    : HHT(INTEGRATOR_TAGS_LadrunoHHT, _alpha, _beta, _gamma),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS()
{

}


LadrunoHHT::~LadrunoHHT()
{
    if (massMatrixMultiplicator != 0)
        delete massMatrixMultiplicator;
    if (dampingMatrixMultiplicator != 0)
        delete dampingMatrixMultiplicator;
}


int
LadrunoHHT::formEleResidual(FE_Element *theEle)
{
    if (sensitivityFlag == 0) {  // no sensitivity: same path as HHT
        this->TransientIntegrator::formEleResidual(theEle);
    } else {

        theEle->zeroResidual();

        // HHT residual at the intermediate state:
        //   R = F(t+alpha*dt) - P(Ualpha) - C*Ualphadot - M*Udotdot
        // U, Udot and Udotdot follow the Newmark relations, so a2..a8 are the
        // Newmark constants (see Newmark::formEleResidual):
        //   dUdotdot/dh = c3*dU/dh + a2*dUn + a3*dVn + a4*dAn
        //   dUdot/dh    = c2*dU/dh + a6*dUn + a7*dVn + a8*dAn
        double a2 = -c3;
        double a3 = -c2/gamma;
        double a4 = 1.0 - 1.0/(2.0*beta);
        double a6 = -c2;
        double a7 = 1.0 - gamma/beta;
        double dt = gamma/(beta*c2);
        double a8 = dt*(1.0 - gamma/(2.0*beta));

        // response sensitivities of the previous step dUn, dVn, dAn
        int vectorSize = U->Size();
        Vector dUn(vectorSize);
        Vector dVn(vectorSize);
        Vector dAn(vectorSize);
        int i, loc;

        AnalysisModel *myModel = this->getAnalysisModel();
        DOF_GrpIter &theDOFs = myModel->getDOFs();
        DOF_Group *dofPtr;
        while ((dofPtr = theDOFs()) != 0) {
            const ID &id = dofPtr->getID();
            int idSize = id.Size();

            const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
            for (i = 0; i < idSize; i++) {
                loc = id(i);
                if (loc >= 0)
                    dUn(loc) = dispSens(i);
            }
            const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
            for (i = 0; i < idSize; i++) {
                loc = id(i);
                if (loc >= 0)
                    dVn(loc) = velSens(i);
            }
            const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
            for (i = 0; i < idSize; i++) {
                loc = id(i);
                if (loc >= 0)
                    dAn(loc) = accelSens(i);
            }
        }

        // mass multiplier: M acts at t+dt, Newmark form (no alpha)
        //   tmp1 = a2*dUn + a3*dVn + a4*dAn
        Vector tmp1(vectorSize);
        tmp1.addVector(0.0, dUn, a2);
        tmp1.addVector(1.0, dVn, a3);
        tmp1.addVector(1.0, dAn, a4);

        // damping multiplier: C acts at Ualphadot, whose part independent of
        // dU/dh is (1-alpha)*dVn + alpha*(a6*dUn + a7*dVn + a8*dAn)
        Vector tmp2(vectorSize);
        tmp2.addVector(0.0, dUn, alpha*a6);
        tmp2.addVector(1.0, dVn, (1.0-alpha) + alpha*a7);
        tmp2.addVector(1.0, dAn, alpha*a8);

        if (massMatrixMultiplicator == 0)
            massMatrixMultiplicator = new Vector(tmp1.Size());
        if (dampingMatrixMultiplicator == 0)
            dampingMatrixMultiplicator = new Vector(tmp2.Size());

        (*massMatrixMultiplicator) = tmp1;
        (*dampingMatrixMultiplicator) = tmp2;

        // -dP/dh at fixed displacement (element state is Ualpha)
        theEle->addResistingForceSensitivity(gradNumber);

        // -dM/dh * acceleration at t+dt
        theEle->addM_ForceSensitivity(gradNumber, *Udotdot, -1.0);

        // -(1-alpha)*K*dUn, the (1-alpha) part of dUalpha/dh; addK_Force uses
        // the current tangent, consistent with the left-hand side
        theEle->addK_Force(dUn, -(1.0-alpha));

        // -M*(a2*dUn + a3*dVn + a4*dAn)
        theEle->addM_Force(*massMatrixMultiplicator, -1.0);

        // -C*(alpha-weighted damping multiplier)
        theEle->addD_Force(*dampingMatrixMultiplicator, -1.0);

        // -dC/dh * velocity at Ualphadot
        theEle->addD_ForceSensitivity(gradNumber, *Ualphadot, -1.0);
    }

    return 0;
}


int
LadrunoHHT::formNodUnbalance(DOF_Group *theDof)
{
    if (sensitivityFlag == 0) {  // no sensitivity: same path as HHT
        this->TransientIntegrator::formNodUnbalance(theDof);
    } else {

        theDof->zeroUnbalance();

        // The multipliers are built in formEleResidual. A model with nodal
        // mass or damping and no FE_Element never calls it, so allocate them
        // here (a new Vector is zero, the correct contribution in that case).
        if (massMatrixMultiplicator == 0)
            massMatrixMultiplicator = new Vector(U->Size());
        if (dampingMatrixMultiplicator == 0)
            dampingMatrixMultiplicator = new Vector(U->Size());

        // nodal contributions carry M and C only (no element stiffness)

        // -M*(a2*dUn + a3*dVn + a4*dAn)
        theDof->addM_Force(*massMatrixMultiplicator, -1.0);

        // -dM/dh * acceleration
        theDof->addM_ForceSensitivity(*Udotdot, -1.0);

        // -C*(alpha-weighted damping multiplier)
        theDof->addD_Force(*dampingMatrixMultiplicator, -1.0);

        // -dC/dh * velocity at Ualphadot
        theDof->addD_ForceSensitivity(*Ualphadot, -1.0);

        // nodal load sensitivities (formed by applyLoadSensitivity)
        theDof->addPtoUnbalance();
    }

    return 0;
}


int
LadrunoHHT::revertToStart()
{
    if (Ut != 0)        Ut->Zero();
    if (Utdot != 0)     Utdot->Zero();
    if (Utdotdot != 0)  Utdotdot->Zero();
    if (U != 0)         U->Zero();
    if (Udot != 0)      Udot->Zero();
    if (Udotdot != 0)   Udotdot->Zero();
    if (Ualpha != 0)    Ualpha->Zero();
    if (Ualphadot != 0) Ualphadot->Zero();

    return 0;
}


int
LadrunoHHT::formSensitivityRHS(int passedGradNumber)
{
    sensitivityFlag = 1;
    gradNumber = passedGradNumber;

    LinearSOE *theSOE = this->getLinearSOE();

    // optional independent part of the right-hand side
    if (assemblyFlag != 0)
        theSOE->setB(independentRHS);

    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = theModel->getDomainPtr();

    // zero the nodal unbalanced loads
    Node *nodePtr;
    NodeIter &theNodeIter = theDomain->getNodes();
    while ((nodePtr = theNodeIter()) != 0)
        nodePtr->zeroUnbalancedLoad();

    // external load sensitivity (load values and time series).
    // The HHT residual is formed at t_n + alpha*dt. With -computeAtEachStep
    // the sensitivities are computed before commit and getCurrentTime()
    // returns that time, so the load gradient matches the residual. After
    // commit (-computeByCommand) the domain time is t_n + dt; the offset of
    // (1-alpha)*dt only affects parameters with dF/dh != 0.
    LoadPattern *loadPatternPtr;
    LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    double time;
    while ((loadPatternPtr = thePatterns()) != 0) {
        time = theDomain->getCurrentTime();
        loadPatternPtr->applyLoadSensitivity(time);
    }

    // element and material contributions
    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();
    while ((elePtr = theEles()) != 0)
        theSOE->addB(elePtr->getResidual(this), elePtr->getID());

    // DOF groups (must be last)
    DOF_Group *dofPtr;
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    while ((dofPtr = theDOFs()) != 0)
        theSOE->addB(dofPtr->getUnbalance(this), dofPtr->getID());

    sensitivityFlag = 0;

    return 0;
}


int
LadrunoHHT::formIndependentSensitivityRHS()
{
    // as in Newmark: not used
    return 0;
}


int
LadrunoHHT::saveSensitivity(const Vector &vNew, int gradNum, int numGrads)
{
    // U, Udot and Udotdot follow the Newmark relations, so the update of the
    // response sensitivities is the Newmark one (no alpha)
    double a1 = c3;
    double a2 = -c3;
    double a3 = -c2/gamma;
    double a4 = 1.0 - 1.0/(2.0*beta);
    double a5 = c2;
    double a6 = -c2;
    double a7 = 1.0 - gamma/beta;
    double dt = gamma/(beta*c2);
    double a8 = dt*(1.0 - gamma/(2.0*beta));

    int vectorSize = U->Size();
    Vector dUn(vectorSize);
    Vector dVn(vectorSize);
    Vector dAn(vectorSize);
    int i, loc;

    AnalysisModel *myModel = this->getAnalysisModel();
    DOF_GrpIter &theDOFs = myModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();
        const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
        for (i = 0; i < idSize; i++) {
            loc = id(i);
            if (loc >= 0)
                dUn(loc) = dispSens(i);
        }
        const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
        for (i = 0; i < idSize; i++) {
            loc = id(i);
            if (loc >= 0)
                dVn(loc) = velSens(i);
        }
        const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
        for (i = 0; i < idSize; i++) {
            loc = id(i);
            if (loc >= 0)
                dAn(loc) = accelSens(i);
        }
    }

    Vector vdotNew(vectorSize);
    Vector vdotdotNew(vectorSize);
    // vdotdotNew = vNew*a1 + dUn*a2 + dVn*a3 + dAn*a4
    vdotdotNew.addVector(0.0, vNew, a1);
    vdotdotNew.addVector(1.0, dUn, a2);
    vdotdotNew.addVector(1.0, dVn, a3);
    vdotdotNew.addVector(1.0, dAn, a4);
    // vdotNew = vNew*a5 + dUn*a6 + dVn*a7 + dAn*a8
    vdotNew.addVector(0.0, vNew, a5);
    vdotNew.addVector(1.0, dUn, a6);
    vdotNew.addVector(1.0, dVn, a7);
    vdotNew.addVector(1.0, dAn, a8);

    DOF_GrpIter &theDOFGrps = myModel->getDOFs();
    DOF_Group *dofPtr1;
    while ((dofPtr1 = theDOFGrps()) != 0)
        dofPtr1->saveSensitivity(vNew, vdotNew, vdotdotNew, gradNum, numGrads);

    return 0;
}


int
LadrunoHHT::commitSensitivity(int gradNum, int numGrads)
{
    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();
    while ((elePtr = theEles()) != 0)
        elePtr->commitSensitivity(gradNum, numGrads);

    return 0;
}


int
LadrunoHHT::computeSensitivities()
{
    LinearSOE *theSOE = this->getLinearSOE();

    theSOE->zeroB();
    this->formIndependentSensitivityRHS();

    // Re-form the current tangent at the converged state once; the solves
    // for every parameter below reuse this factorization. The factorization
    // left by the primal algorithm can be stale (ModifiedNewton) or built from
    // the initial stiffness (-initial); with ModifiedNewton on a yielding
    // Steel01 truss reusing it gives a relative gradient error of about 1e-3.
    this->formTangent(CURRENT_TANGENT);

    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = theModel->getDomainPtr();
    ParameterIter &paramIter = theDomain->getParameters();
    Parameter *theParam;

    // de-activate all parameters
    while ((theParam = paramIter()) != 0)
        theParam->activate(false);

    int numGrads = theDomain->getNumParameters();
    paramIter = theDomain->getParameters();
    while ((theParam = paramIter()) != 0) {
        theParam->activate(true);

        theSOE->zeroB();
        int gradIndex = theParam->getGradIndex();
        this->formSensitivityRHS(gradIndex);
        theSOE->solve();
        this->saveSensitivity(theSOE->getX(), gradIndex, numGrads);
        this->commitSensitivity(gradIndex, numGrads);

        theParam->activate(false);
    }

    return 0;
}


const char *
LadrunoHHT::getClassType(void) const
{
    return "LadrunoHHT";
}


void
LadrunoHHT::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "LadrunoHHT (DDM sensitivity) - currentTime: " << currentTime << endln;
        s << "  alpha: " << alpha;
        s << "  beta: " << beta << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else
        s << "LadrunoHHT - no associated AnalysisModel\n";
}
