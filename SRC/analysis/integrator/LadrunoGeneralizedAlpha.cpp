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
//   J. Chung, G.M. Hulbert (1993). A time integration algorithm for structural
//     dynamics with improved numerical dissipation: the generalized-alpha
//     method. Journal of Applied Mechanics 60(2):371-375.
//     doi:10.1115/1.2900803
//   J.P. Conte, P.K. Vijalapura, M. Meghella (2003). Consistent finite-element
//     response sensitivity analysis. Journal of Engineering Mechanics
//     129(12):1380-1393. doi:10.1061/(ASCE)0733-9399(2003)129:12(1380)
//
// Description: generalized-alpha integrator with DDM response sensitivities.
// See LadrunoGeneralizedAlpha.h for the derivation. The sensitivity residual
// requires GeneralizedAlpha::update() to evaluate the inertia at
// Udotdot_aM = (1-alphaM)*Udotdot_n + alphaM*Udotdot_{n+1}; the mass terms
// below are the derivative of that residual.

#include <LadrunoGeneralizedAlpha.h>
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
OPS_LadrunoGeneralizedAlpha(void)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;

    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 2 && argc != 4) {
        opserr << "WARNING - incorrect number of args want LadrunoGeneralizedAlpha $alphaM $alphaF <$beta $gamma>\n";
        return 0;
    }

    double dData[4];
    if (OPS_GetDouble(&argc, dData) != 0) {
        opserr << "WARNING - invalid args want LadrunoGeneralizedAlpha $alphaM $alphaF <$beta $gamma>\n";
        return 0;
    }

    // argument order as in GeneralizedAlpha: alphaM, alphaF, beta, gamma
    if (argc == 2)
        theIntegrator = new LadrunoGeneralizedAlpha(dData[0], dData[1]);
    else
        theIntegrator = new LadrunoGeneralizedAlpha(dData[0], dData[1], dData[2], dData[3]);

    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating LadrunoGeneralizedAlpha integrator\n";

    return theIntegrator;
}


LadrunoGeneralizedAlpha::LadrunoGeneralizedAlpha()
    : GeneralizedAlpha(INTEGRATOR_TAGS_LadrunoGeneralizedAlpha, 0.0, 0.0, 0.25, 0.5),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS(), sensTangentFlag(0)
{

}


LadrunoGeneralizedAlpha::LadrunoGeneralizedAlpha(double _alphaM, double _alphaF)
    // same beta and gamma that GeneralizedAlpha derives from (alphaM, alphaF)
    : GeneralizedAlpha(INTEGRATOR_TAGS_LadrunoGeneralizedAlpha, _alphaM, _alphaF,
                       (1+_alphaM-_alphaF)*(1+_alphaM-_alphaF)*0.25, 0.5+_alphaM-_alphaF),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS(), sensTangentFlag(0)
{

}


LadrunoGeneralizedAlpha::LadrunoGeneralizedAlpha(double _alphaM, double _alphaF,
                                                 double _beta, double _gamma)
    : GeneralizedAlpha(INTEGRATOR_TAGS_LadrunoGeneralizedAlpha, _alphaM, _alphaF, _beta, _gamma),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS(), sensTangentFlag(0)
{

}


LadrunoGeneralizedAlpha::~LadrunoGeneralizedAlpha()
{
    if (massMatrixMultiplicator != 0)
        delete massMatrixMultiplicator;
    if (dampingMatrixMultiplicator != 0)
        delete dampingMatrixMultiplicator;
}


int
LadrunoGeneralizedAlpha::formEleResidual(FE_Element *theEle)
{
    if (sensitivityFlag == 0) {  // no sensitivity: same path as GeneralizedAlpha
        this->TransientIntegrator::formEleResidual(theEle);
    } else {

        theEle->zeroResidual();

        // K and C act at the alphaF state (Ualpha, Ualphadot), M at the
        // alphaM state (Ualphadotdot). U, Udot and Udotdot follow the Newmark
        // relations, so a2..a8 are the Newmark constants:
        //   dUdotdot/dh = c3*dU/dh + a2*dUn + a3*dVn + a4*dAn
        //   dUdot/dh    = c2*dU/dh + a6*dUn + a7*dVn + a8*dAn
        // and the intermediate-state sensitivities are
        //   dUalpha/dh       = (1-alphaF)*dUn + alphaF*dU/dh
        //   dUalphadot/dh    = (1-alphaF)*dVn + alphaF*dUdot/dh
        //   dUalphadotdot/dh = (1-alphaM)*dAn + alphaM*dUdotdot/dh
        // The dU/dh terms (alphaF*K, alphaF*c2*C, alphaM*c3*M) form the
        // left-hand side, identical to formEleTangent.
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

        // mass multiplier: part of dUalphadotdot/dh independent of dU/dh
        //   tmp1 = alphaM*a2*dUn + alphaM*a3*dVn + ((1-alphaM) + alphaM*a4)*dAn
        Vector tmp1(vectorSize);
        tmp1.addVector(0.0, dUn, alphaM*a2);
        tmp1.addVector(1.0, dVn, alphaM*a3);
        tmp1.addVector(1.0, dAn, (1.0-alphaM) + alphaM*a4);

        // damping multiplier: part of dUalphadot/dh independent of dU/dh
        //   tmp2 = alphaF*a6*dUn + ((1-alphaF) + alphaF*a7)*dVn + alphaF*a8*dAn
        Vector tmp2(vectorSize);
        tmp2.addVector(0.0, dUn, alphaF*a6);
        tmp2.addVector(1.0, dVn, (1.0-alphaF) + alphaF*a7);
        tmp2.addVector(1.0, dAn, alphaF*a8);

        if (massMatrixMultiplicator == 0)
            massMatrixMultiplicator = new Vector(tmp1.Size());
        if (dampingMatrixMultiplicator == 0)
            dampingMatrixMultiplicator = new Vector(tmp2.Size());

        (*massMatrixMultiplicator) = tmp1;
        (*dampingMatrixMultiplicator) = tmp2;

        // -dP/dh at fixed displacement (element state is Ualpha)
        theEle->addResistingForceSensitivity(gradNumber);

        // -dM/dh * acceleration at Ualphadotdot
        theEle->addM_ForceSensitivity(gradNumber, *Ualphadotdot, -1.0);

        // -(1-alphaF)*K*dUn, the (1-alphaF) part of dUalpha/dh; addK_Force
        // uses the current tangent, consistent with the left-hand side
        theEle->addK_Force(dUn, -(1.0-alphaF));

        // -M*(alphaM-weighted mass multiplier)
        theEle->addM_Force(*massMatrixMultiplicator, -1.0);

        // -C*(alphaF-weighted damping multiplier)
        theEle->addD_Force(*dampingMatrixMultiplicator, -1.0);

        // -dC/dh * velocity at Ualphadot
        theEle->addD_ForceSensitivity(gradNumber, *Ualphadot, -1.0);
    }

    return 0;
}


int
LadrunoGeneralizedAlpha::formNodUnbalance(DOF_Group *theDof)
{
    if (sensitivityFlag == 0) {  // no sensitivity: same path as GeneralizedAlpha
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

        // -M*(alphaM-weighted mass multiplier)
        theDof->addM_Force(*massMatrixMultiplicator, -1.0);

        // -dM/dh * acceleration at Ualphadotdot
        theDof->addM_ForceSensitivity(*Ualphadotdot, -1.0);

        // -C*(alphaF-weighted damping multiplier)
        theDof->addD_Force(*dampingMatrixMultiplicator, -1.0);

        // -dC/dh * velocity at Ualphadot
        theDof->addD_ForceSensitivity(*Ualphadot, -1.0);

        // nodal load sensitivities (formed by applyLoadSensitivity)
        theDof->addPtoUnbalance();
    }

    return 0;
}


int
LadrunoGeneralizedAlpha::formEleTangent(FE_Element *theEle)
{
    if (sensTangentFlag == 0)   // primal path: GeneralizedAlpha tangent
        return this->GeneralizedAlpha::formEleTangent(theEle);

    // sensitivity solve: Jacobian of the residual with respect to U_{n+1},
    // alphaF*K + alphaF*c2*C + alphaM*c3*M (the GeneralizedAlpha coefficients)
    theEle->zeroTangent();
    if (statusFlag == INITIAL_TANGENT)
        theEle->addKiToTang(alphaF*c1);
    else
        theEle->addKtToTang(alphaF*c1);
    theEle->addCtoTang(alphaF*c2);
    theEle->addMtoTang(alphaM*c3);

    return 0;
}


int
LadrunoGeneralizedAlpha::formNodTangent(DOF_Group *theDof)
{
    if (sensTangentFlag == 0)   // primal path: GeneralizedAlpha tangent
        return this->GeneralizedAlpha::formNodTangent(theDof);

    // sensitivity solve: nodal M and C with the GeneralizedAlpha coefficients
    theDof->zeroTangent();
    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaM*c3);

    return 0;
}


int
LadrunoGeneralizedAlpha::revertToStart()
{
    if (Ut != 0)           Ut->Zero();
    if (Utdot != 0)        Utdot->Zero();
    if (Utdotdot != 0)     Utdotdot->Zero();
    if (U != 0)            U->Zero();
    if (Udot != 0)         Udot->Zero();
    if (Udotdot != 0)      Udotdot->Zero();
    if (Ualpha != 0)       Ualpha->Zero();
    if (Ualphadot != 0)    Ualphadot->Zero();
    if (Ualphadotdot != 0) Ualphadotdot->Zero();

    return 0;
}


int
LadrunoGeneralizedAlpha::formSensitivityRHS(int passedGradNumber)
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
    // The residual is formed at t_n + alphaF*dt. With -computeAtEachStep the
    // sensitivities are computed before commit and getCurrentTime() returns
    // that time, so the load gradient matches the residual. After commit
    // (-computeByCommand) the domain time is t_n + dt; the offset of
    // (1-alphaF)*dt only affects parameters with dF/dh != 0.
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
LadrunoGeneralizedAlpha::formIndependentSensitivityRHS()
{
    // as in Newmark: not used
    return 0;
}


int
LadrunoGeneralizedAlpha::saveSensitivity(const Vector &vNew, int gradNum, int numGrads)
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
LadrunoGeneralizedAlpha::commitSensitivity(int gradNum, int numGrads)
{
    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();
    while ((elePtr = theEles()) != 0)
        elePtr->commitSensitivity(gradNum, numGrads);

    return 0;
}


int
LadrunoGeneralizedAlpha::computeSensitivities()
{
    LinearSOE *theSOE = this->getLinearSOE();

    theSOE->zeroB();
    this->formIndependentSensitivityRHS();

    // Re-form the current tangent at the converged state once; the solves
    // for every parameter below reuse this factorization. The factorization
    // left by the primal algorithm can be stale (ModifiedNewton) or built from
    // the initial stiffness (-initial).
    sensTangentFlag = 1;
    this->formTangent(CURRENT_TANGENT);
    sensTangentFlag = 0;

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
LadrunoGeneralizedAlpha::getClassType(void) const
{
    return "LadrunoGeneralizedAlpha";
}


void
LadrunoGeneralizedAlpha::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "LadrunoGeneralizedAlpha (DDM sensitivity) - currentTime: " << currentTime << endln;
        s << "  alphaM: " << alphaM << "  alphaF: " << alphaF;
        s << "  beta: " << beta << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else
        s << "LadrunoGeneralizedAlpha - no associated AnalysisModel\n";
}
