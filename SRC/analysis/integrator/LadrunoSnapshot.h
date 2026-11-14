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

// Developed by: J.A. Abell (Ladruño project)
//
// Description: `ladrunoSnapshot save|restore|clear`, an in-memory copy of the
// committed state of the domain, so that a script can return to an earlier
// equilibrium (e.g. the last converged strength reduction factor after a
// runaway, or the state after gravity) without database save/restore:
//   - nodes: committed displacement, velocity and acceleration;
//   - materials: their committed state, through the parameter names
//     ladrunoSaveState / ladrunoRestoreState that elements pass to their
//     material points (implemented by ASDPlasticMaterial3D).
// Elements whose committed state is more than their nodes and materials, and
// materials that do not implement the two names, are not restored; the
// command reports how many material points took part. restore leaves the
// domain committed at the copy and marks it changed, so the next analyze
// re-seeds the integrator (for LadrunoDynamicRelaxation: velocities and the
// fictitious mass) from it. Material parameters (e.g. strength values) are not
// part of the copy. In a parallel run each process keeps its own subdomain.
// Header-only, shared by the interpreters as LadrunoSolverQuery.h.

#ifndef LadrunoSnapshot_h
#define LadrunoSnapshot_h

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <Parameter.h>
#include <Vector.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <map>
#include <string.h>

namespace Ladruno {

struct SnapshotNode {
    Vector disp, vel, accel;
};

inline std::map<int, SnapshotNode> &snapshotNodes()
{
    static std::map<int, SnapshotNode> nodes;
    return nodes;
}

// a Parameter that is never registered with the domain, used to reach the
// material points of every element by name
class SnapshotParameter : public Parameter
{
  public:
    SnapshotParameter() : Parameter(0) {}
    int count(void) const { return numObjects; }
};

// pass `name` to the material points of every element and invoke it;
// returns the number of material points that accepted it, < 0 on failure
inline int snapshotMaterials(Domain *dom, const char *name)
{
    SnapshotParameter p;
    const char *argv[1] = {name};
    ElementIter &theEles = dom->getElements();
    Element *ele;
    while ((ele = theEles()) != 0)
        p.addComponent(ele, argv, 1);
    if (p.update(1.0) != 0)
        return -1;
    return p.count();
}

inline int snapshotCommand(Domain *dom)
{
    if (dom == 0)
        return -1;
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING ladrunoSnapshot - expects save|restore|clear\n";
        return -1;
    }
    const char *sub = OPS_GetString();
    std::map<int, SnapshotNode> &nodes = snapshotNodes();

    if (strcmp(sub, "clear") == 0) {
        nodes.clear();
        int nm = 0;
        int n = 1;
        OPS_SetIntOutput(&n, &nm, true);
        return 0;
    }

    if (strcmp(sub, "save") == 0) {
        nodes.clear();
        NodeIter &theNodes = dom->getNodes();
        Node *nd;
        while ((nd = theNodes()) != 0) {
            SnapshotNode &s = nodes[nd->getTag()];
            s.disp = nd->getDisp();
            s.vel = nd->getVel();
            s.accel = nd->getAccel();
        }
        int nm = snapshotMaterials(dom, "ladrunoSaveState");
        if (nm < 0) {
            opserr << "WARNING ladrunoSnapshot save - a material point failed to save its state\n";
            return -1;
        }
        if (nm == 0)
            opserr << "WARNING ladrunoSnapshot save - no material point saved its state "
                      "(only ASDPlasticMaterial3D implements it); only the nodes are kept\n";
        int n = 1;
        OPS_SetIntOutput(&n, &nm, true);
        return 0;
    }

    if (strcmp(sub, "restore") == 0) {
        if (nodes.empty()) {
            opserr << "WARNING ladrunoSnapshot restore - no snapshot (use `ladrunoSnapshot save`)\n";
            return -1;
        }
        NodeIter &theNodes = dom->getNodes();
        Node *nd;
        int missing = 0;
        while ((nd = theNodes()) != 0) {
            std::map<int, SnapshotNode>::iterator it = nodes.find(nd->getTag());
            if (it == nodes.end() || it->second.disp.Size() != nd->getNumberDOF()) {
                missing++;
                continue;
            }
            nd->setTrialDisp(it->second.disp);
            nd->setTrialVel(it->second.vel);
            nd->setTrialAccel(it->second.accel);
            nd->commitState();
        }
        if (missing > 0) {
            opserr << "WARNING ladrunoSnapshot restore - " << missing
                   << " nodes were not in the snapshot; the domain changed since save\n";
            return -1;
        }
        int nm = snapshotMaterials(dom, "ladrunoRestoreState");
        if (nm < 0) {
            opserr << "WARNING ladrunoSnapshot restore - a material point had no saved state\n";
            return -1;
        }
        // The elements are not updated here: an update would integrate the
        // materials with their current parameters (e.g. the strengths of the
        // trial that failed), moving the restored stresses. The next step
        // updates them from the restored, committed state.
        // re-seed the analysis (integrator state) from the restored domain
        dom->domainChange();
        int n = 1;
        OPS_SetIntOutput(&n, &nm, true);
        return 0;
    }

    opserr << "WARNING ladrunoSnapshot - unknown subcommand '" << sub
           << "' (use save|restore|clear)\n";
    return -1;
}

} // namespace Ladruno

#endif
