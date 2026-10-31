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


#ifndef EnergyBalanceRecorder_h
#define EnergyBalanceRecorder_h

// Written: jaabell
// Created: 12/2024
// Revision: B (regions, closure residual, nodal mass and damping,
//              trapezoidal integration, -v2 channel layout;
//              N. Mora Bowen, P. Palacios, J.A. Abell, Ladruño project)
//
// References:
//   T. Belytschko, W.K. Liu, B. Moran (2000). Nonlinear Finite Elements for
//     Continua and Structures, Wiley, Section 6.2.3 (energy balance check).
//   C.-M. Uang, V.V. Bertero (1990). Evaluation of seismic energy in
//     structures. Earthquake Engineering and Structural Dynamics 19(1),
//     77-90. doi:10.1002/eqe.4290190108
//
// Description: records a structural-dynamics energy balance for the whole
// model and, optionally, for one or more MeshRegions. Each step writes
// EBR_NUM_ENERGY_COMPONENTS columns per block:
//
//   KE   kinetic energy            1/2 v^T M v          (instantaneous)
//   IE   internal energy           int F^T v dt         (accumulated)
//   DW   damping work              int v^T C v dt       (accumulated)
//   ULW  unbalanced-load work      int v^T Pext dt      (accumulated)
//   RES  residual                  ULW - (KE + IE + DW) (closure error)
//   ERR  relative residual [%]     100 * RES / Eref     (Eref = running
//                                                        max total energy)
//
// Column layout written each step:
//   [time]  (KE IE DW ULW RES ERR)_model  (KE IE DW ULW RES ERR)_region0  ...
//
// M and C include both element matrices (Element::getMass/getDamp, which
// carry the element Rayleigh aM*M + bK*K) and nodal contributions
// (Node::getMass, Node::getDamp = aM*nodalMass). Pext is the external applied
// nodal load (Node::getUnbalancedLoad after Domain::applyLoad). For an element
// whose damping force is not C v (a viscous material whose force is already
// inside the resisting force F), DW uses the power of the force the element
// adds beyond F and inertia, so the dashpot is counted once (in IE).
//
// Known coverage gaps (RES shows their effect):
//   * Modal damping is applied inside the integrator, not through any
//     element or node getDamp(), so DW does not contain it (the -v2 layout
//     reports it as E_modal when the integrator publishes it).
//   * Element loads (eleLoad) assemble into the element unbalance, not into
//     Node::unbalLoad, so they are not captured by ULW.
//   * Under _PARALLEL_PROCESSING each rank records only its own partition
//     (nodal terms counted on every rank holding a shared node). Reduce
//     offline for a global balance. A warning is printed at construction.
//
// The three work integrals use the trapezoidal rule
// (1/2 (rate_{n-1}+rate_n) dt); the first record() uses the rectangle rule.
//
// -v2 (opt-in): channel-aware layout. The model block becomes
//   KE_ele KE_nod IE DW_ele DW_nod ULW [E_inject] [E_lnvd] [E_modal] [E_hg]
//   RES ERR%
// The element/nodal split keeps the nodal terms, which are counted on every
// rank holding a shared boundary node, out of the element terms, which
// reduce exactly across ranks. E_inject is the incident-wave injection of
// absorbing boundaries, derived from the published resisting-sign leak L
// (IE shown = IE_raw - L - E_hg, E_inject = -L). E_lnvd is the local
// non-viscous damping dissipation. E_inject is published by LysmerTriangle
// and ASDAbsorbingBoundary2D/3D (bottom and lateral boundaries). E_modal is
// the modal-damping dissipation published by IncrementalIntegrator::commit()
// (integrators that override commit() without calling it do not publish, and
// their modal dissipation stays in RES). E_hg is the
// hourglass-stabilization energy,
// the sum of the element "hourglassEnergy" responses; the column is present
// when the model contains an element answering that response at the first
// record. Channel columns appear only when a producer declared them
// (LadrunoEnergyChannels.h). Region blocks keep the six-column layout in
// both modes. The channel set is fixed at the first record and reset by
// wipe(); parse -v2 output by the echoed or XML column names.
//
// -ownedNodes <regionTag> (opt-in): count the nodal terms (nodal-mass KE,
// nodal Rayleigh DW, nodal-load ULW) only for the nodes of the given
// MeshRegion. Under flat per-rank MPI (openseesmp) shared boundary nodes are
// built on every rank that touches them, so their nodal terms would be
// counted several times in a cross-rank sum; passing each rank the set of
// nodes it owns makes the per-rank nodal blocks sum exactly. Applies to the
// nodal terms of both layouts. Without the flag every node is counted. A
// missing region prints a warning and counts no nodes.

#include <Recorder.h>
#include <ID.h>
#include <Vector.h>
#include <TimeSeries.h>
#include <EnergyBalanceKernel.h>

#include <vector>
#include <unordered_map>
#include <unordered_set>

class Domain;
class FE_Datastore;
class Node;
class Element;
class MeshRegion;
class Response;

#define EBR_NUM_ENERGY_COMPONENTS 6

class EnergyBalanceRecorder: public Recorder
{
public:
    EnergyBalanceRecorder();
    EnergyBalanceRecorder(
        Domain &theDomain,
        OPS_Stream &theOutputHandler,
        bool echoTimeFlag = true,
        const ID &regions = ID(),
        bool v2Layout = false,
        int ownedNodesRegion = -1);

    ~EnergyBalanceRecorder();

    int record(int commitTag, double timeStamp);
    int flush();

    int domainChanged(void);
    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);


protected:

private:
    int initialize(void);
    void buildCache(void);

    Domain *theDomain;
    OPS_Stream *theOutputHandler;

    Vector response;

    bool echoTimeFlag;       // include the time column in the output
    bool initializationDone; // response and accumulators sized and zeroed
    bool firstRecord;        // first record() uses the rectangle rule

    double time_last;

    // whole-model accumulator (work integrals + closure/ERR), six columns
    ebkernel::EnergyAccumulator model_acc;

    // -v2 layout. Which channel columns exist is fixed at initialize() from
    // the registry's declared set (producers declare at construction or
    // setDomain, so define them before the first analyze).
    bool v2Layout;
    bool chInject;               // ABSORB_LEAK declared -> E_inject column
    bool chLnvd;                 // LNVD_WORK declared   -> E_lnvd column
    bool chModal;                // MODAL_WORK declared  -> E_modal column
    bool chHourglass;            // hourglassEnergy response found -> E_hg
    int  numModelCols;           // model-block width (6, or the -v2 width)
    ebkernel::EnergyAccumulatorV2 model_acc2;

    // Hourglass energy: Response objects found by probing every element for
    // a "hourglassEnergy" response, rebuilt (old ones deleted) with the
    // cache. Each entry keeps the element tag so record() can check that
    // getElement(tag) still returns the cached pointer before calling
    // getResponse(): `remove element` deletes the Element the Response
    // wraps, and recorders are not notified of domain changes.
    struct HgEntry {
        int tag;
        Element *ele;
        Response *res;
    };
    std::vector<HgEntry> hgCache;

    // -ownedNodes: region tag (-1 = inactive) and the resolved node tags;
    // nodal terms are counted only for members.
    int ownedRegionTag;
    std::unordered_set<int> ownedNodes;

    // per-region accumulators
    ID regionTags;
    int numRegions;
    std::vector<ebkernel::EnergyAccumulator> region_acc;

    // Cache of region membership, keyed by tag (a freed pointer could be
    // reused by a new entity; tags cannot alias). cachedNumElements/Nodes
    // are checked every record() and force a rebuild on any change.
    bool cacheValid;
    std::unordered_map<int, std::vector<int> > elemRegions;
    std::unordered_map<int, std::vector<int> > nodeRegions;
    int cachedNumElements;
    int cachedNumNodes;
};

#endif
