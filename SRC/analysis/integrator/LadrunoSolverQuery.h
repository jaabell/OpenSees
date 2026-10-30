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
// Description: subcommand dispatch of the runtime solver commands that
// query and steer an active integrator (`ladrunoArcLength`, `ladrunoDR`).
// The dispatch
// is shared by the classic Tcl interpreter (SRC/tcl/commands.cpp) and the
// Python and Tcl wrappers (SRC/interpreter). It is header-only so that each
// including translation unit binds the elementAPI calls (OPS_GetString,
// OPS_SetDoubleOutput, ...) to the backend of its own interpreter, with no
// link dependency between the interpreter libraries. The integrator is a
// parameter because each interpreter keeps its own integrator pointer.

#ifndef LadrunoSolverQuery_h
#define LadrunoSolverQuery_h

#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <LadrunoArcLength.h>
#include <LadrunoDynamicRelaxation.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <string.h>

namespace Ladruno {

// ladrunoArcLength <sub <value>>: runtime control of the active
// LadrunoArcLength, for script-driven cut-and-retry:
//   ok = analyze(1)
//   while ok != 0:
//       ladrunoArcLength('reduceStep', 0.5); ok = analyze(1)
// Mutators (one value, echo the result): reduceStep f | increaseStep f |
//   setArcLength v (echo the radius) | scaleCVisc f (echo the dissipation
//   ratio). Queries: arcLength | deltaLambdaStep | currentLambda | sign |
//   deltaUstepNorm | dissipationRatio | dissipatedEnergy | referenceEnergy.
// Without a subcommand the current radius is returned.
inline int arcLengthCommand(StaticIntegrator *si)
{
    if (si == 0 || si->getClassTag() != INTEGRATOR_TAGS_LadrunoArcLength) {
        opserr << "WARNING ladrunoArcLength - the active static integrator is not "
                  "a LadrunoArcLength (set `integrator LadrunoArcLength ...` first)\n";
        return -1;
    }
    LadrunoArcLength *la = (LadrunoArcLength *)si;

    if (OPS_GetNumRemainingInputArgs() < 1) {
        double al = la->getArcLength();
        int nd = 1;
        OPS_SetDoubleOutput(&nd, &al, true);
        return 0;
    }

    const char *sub = OPS_GetString();

    if (strcmp(sub, "reduceStep") == 0 || strcmp(sub, "increaseStep") == 0 ||
        strcmp(sub, "setArcLength") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING ladrunoArcLength " << sub << " - expects one value\n";
            return -1;
        }
        int nd = 1;
        double v = 0.0;
        if (OPS_GetDoubleInput(&nd, &v) < 0) {
            opserr << "WARNING ladrunoArcLength " << sub << " - failed to read value\n";
            return -1;
        }
        int res = 0;
        if (strcmp(sub, "reduceStep") == 0)        res = la->reduceStep(v);
        else if (strcmp(sub, "increaseStep") == 0) res = la->increaseStep(v);
        else                                       res = la->setArcLength(v);
        if (res < 0) return -1;
        double al = la->getArcLength();
        int n = 1;
        OPS_SetDoubleOutput(&n, &al, true);
        return 0;
    }

    if (strcmp(sub, "scaleCVisc") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING ladrunoArcLength scaleCVisc - expects one value\n";
            return -1;
        }
        int nd = 1;
        double v = 0.0;
        if (OPS_GetDoubleInput(&nd, &v) < 0) {
            opserr << "WARNING ladrunoArcLength scaleCVisc - failed to read value\n";
            return -1;
        }
        if (la->scaleCVisc(v) < 0) return -1;
        double r = la->getStabilizationDissipationRatio();
        int n = 1;
        OPS_SetDoubleOutput(&n, &r, true);
        return 0;
    }

    double out = 0.0;
    if      (strcmp(sub, "arcLength") == 0)        out = la->getArcLength();
    else if (strcmp(sub, "deltaLambdaStep") == 0)  out = la->getDeltaLambdaStep();
    else if (strcmp(sub, "currentLambda") == 0)    out = la->getCurrentLambda();
    else if (strcmp(sub, "sign") == 0)             out = (double)la->getSignLastDeltaLambdaStep();
    else if (strcmp(sub, "deltaUstepNorm") == 0)   out = la->getDeltaUstepNorm();
    else if (strcmp(sub, "dissipationRatio") == 0) out = la->getStabilizationDissipationRatio();
    else if (strcmp(sub, "dissipatedEnergy") == 0) out = la->getStabilizationDissipatedEnergy();
    else if (strcmp(sub, "referenceEnergy") == 0)  out = la->getReferenceStrainEnergy();
    else {
        opserr << "WARNING ladrunoArcLength - unknown subcommand '" << sub
               << "' (use reduceStep|increaseStep|setArcLength|scaleCVisc|arcLength|"
                  "deltaLambdaStep|currentLambda|sign|deltaUstepNorm|"
                  "dissipationRatio|dissipatedEnergy|referenceEnergy)\n";
        return -1;
    }
    int n = 1;
    OPS_SetDoubleOutput(&n, &out, true);
    return 0;
}

// ladrunoDR <sub>: read-only state of the active LadrunoDynamicRelaxation.
//   residualNorm    ||f_ext - f_int||_inf of the last step, the force
//                   criterion for rest
//   kineticEnergy   1/2 v^T M* v
//   stabilityMargin largest (omega_max dt / 2)^2 measured since the last
//                   domainChanged(); <= 1 inside the explicit stability
//                   limit; -1 for a lumped or unity mass, -2 with the probe
//                   off (not measured)
inline int drCommand(TransientIntegrator *ti)
{
    if (ti == 0 || ti->getClassTag() != INTEGRATOR_TAGS_LadrunoDynamicRelaxation) {
        opserr << "WARNING ladrunoDR - the active transient integrator is not a "
                  "LadrunoDynamicRelaxation (set `integrator LadrunoDynamicRelaxation "
                  "...` first)\n";
        return -1;
    }
    LadrunoDynamicRelaxation *dr = (LadrunoDynamicRelaxation *)ti;

    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING ladrunoDR - expects a subcommand "
                  "(residualNorm|kineticEnergy|stabilityMargin)\n";
        return -1;
    }
    const char *sub = OPS_GetString();
    double out = 0.0;
    if      (strcmp(sub, "residualNorm") == 0)    out = dr->getResidualNorm();
    else if (strcmp(sub, "kineticEnergy") == 0)   out = dr->getKineticEnergy();
    else if (strcmp(sub, "stabilityMargin") == 0) out = dr->getStabilityMargin();
    else {
        opserr << "WARNING ladrunoDR - unknown subcommand '" << sub
               << "' (use residualNorm|kineticEnergy|stabilityMargin)\n";
        return -1;
    }
    int n = 1;
    OPS_SetDoubleOutput(&n, &out, true);
    return 0;
}

} // namespace Ladruno

#endif
