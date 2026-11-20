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

// Written: jaabell (12/2024)
// Revision: B (regions, closure residual, nodal mass and damping,
//              trapezoidal integration, -v2 channel layout;
//              N. Mora Bowen, P. Palacios, J.A. Abell, Ladruño project)
//
// Description: see EnergyBalanceRecorder.h for the energy-balance definition,
// the column layouts and the coverage gaps. The per-entity terms and the
// per-scope integration live in EnergyBalanceKernel.h.
//
// What: "@(#) EnergyBalanceRecorder.C, revB"

#include <EnergyBalanceRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <MeshRegion.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <TCP_Stream.h>
#include <DummyStream.h>

#include <Response.h>
#include <Information.h>

#include <elementAPI.h>
#include <classTags.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>

void*
OPS_EnergyBalanceRecorder()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: recorder EnergyBalance -file <fileName> <-time> "
                  "<-region tag ...> <-v2> <-ownedNodes regionTag>\n";
        return 0;
    }

    OPS_Stream *theOutputStream = 0;
    const char* filename = 0;

    const int STANDARD_STREAM = 0;
    const int DATA_STREAM = 1;
    const int XML_STREAM = 2;
    const int BINARY_STREAM = 4;
    const int DATA_STREAM_CSV = 5;
    const int TCP_STREAM = 6;
    const int DATA_STREAM_ADD = 7;

    int eMode = STANDARD_STREAM;

    bool echoTimeFlag = false;
    bool doScientific = false;
    int precision = 6;
    bool closeOnWrite = false;
    bool v2Layout = false;
    int ownedNodesRegion = -1;

    const char *inetAddr = 0;
    int inetPort = 0;

    ID regions(0, 6);

    while (OPS_GetNumRemainingInputArgs() > 0) {

        const char* option = OPS_GetString();

        if (strcmp(option, "-time") == 0) {
            echoTimeFlag = true;
        }
        else if (strcmp(option, "-load") == 0) {
            // alias of -time
            echoTimeFlag = true;
        }
        else if (strcmp(option, "-v2") == 0) {
            // channel-aware layout with the element/nodal split
            v2Layout = true;
        }
        else if (strcmp(option, "-scientific") == 0) {
            doScientific = true;
        }
        else if (strcmp(option, "-file") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = DATA_STREAM;
        }
        else if (strcmp(option, "-fileAdd") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = DATA_STREAM_ADD;
        }
        else if (strcmp(option, "-closeOnWrite") == 0) {
            closeOnWrite = true;
        }
        else if (strcmp(option, "-csv") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = DATA_STREAM_CSV;
        }
        else if (strcmp(option, "-tcp") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                inetAddr = OPS_GetString();
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &inetPort) < 0) {
                    opserr << "WARNING: EnergyBalance -tcp failed to read inetPort\n";
                    return 0;
                }
            }
            eMode = TCP_STREAM;
        }
        else if (strcmp(option, "-xml") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = XML_STREAM;
        }
        else if (strcmp(option, "-binary") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = BINARY_STREAM;
        }
        else if (strcmp(option, "-precision") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &precision) < 0) {
                    opserr << "WARNING: EnergyBalance -precision failed to read value\n";
                    return 0;
                }
            }
        }
        else if (strcmp(option, "-region") == 0) {
            // -region <tag> ; repeatable, e.g. ... -region 1 -region 2
            // Each named region adds one (KE IE DW ULW RES ERR) block after
            // the whole-model block. Region energy is summed over the
            // region's elements (KE, IE, DW) and nodes (KE, DW, ULW).
            int tag;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &tag) < 0) {
                    opserr << "WARNING: EnergyBalance -region requires an integer region tag\n";
                    return 0;
                }
                regions[regions.Size()] = tag;
            }
        }
        else if (strcmp(option, "-ownedNodes") == 0) {
            // count the nodal terms (KE_nod, DW_nod, ULW) only for nodes in
            // this MeshRegion (per-rank ownership under flat MPI, see header)
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &ownedNodesRegion) < 0) {
                    opserr << "WARNING: EnergyBalance -ownedNodes requires an "
                              "integer region tag\n";
                    return 0;
                }
            }
        }
    }

    // data handler
    if (eMode == DATA_STREAM && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else if (eMode == DATA_STREAM_ADD && filename != 0)
        theOutputStream = new DataFileStreamAdd(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else if (eMode == DATA_STREAM_CSV && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 1, closeOnWrite, precision, doScientific);
    else if (eMode == XML_STREAM && filename != 0)
        theOutputStream = new XmlFileStream(filename);
    else if (eMode == BINARY_STREAM && filename != 0)
        theOutputStream = new BinaryFileStream(filename);
    else if (eMode == TCP_STREAM && inetAddr != 0)
        theOutputStream = new TCP_Stream(inetPort, inetAddr);
    else
        theOutputStream = new StandardStream();

    theOutputStream->setPrecision(precision);

    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        delete theOutputStream;
        return 0;
    }

#ifdef _PARALLEL_PROCESSING
    opserr << "EnergyBalanceRecorder: WARNING under parallel processing the energy "
              "balance is computed per partition (local domain only); totals are not "
              "reduced across ranks and nodal terms are counted on every rank holding "
              "a shared node. Reduce the per-rank output offline for a global balance.\n";
#endif

    EnergyBalanceRecorder* recorder = new EnergyBalanceRecorder(
        *domain, *theOutputStream, echoTimeFlag, regions, v2Layout,
        ownedNodesRegion);

    return recorder;
}


EnergyBalanceRecorder::EnergyBalanceRecorder()
    : Recorder(RECORDER_TAGS_EnergyBalanceRecorder),
      theDomain(0), theOutputHandler(0),
      response(EBR_NUM_ENERGY_COMPONENTS),
      echoTimeFlag(true), initializationDone(false), firstRecord(true),
      time_last(0.),
      model_acc(),
      v2Layout(false), chInject(false), chLnvd(false),
      chModal(false), chHourglass(false),
      numModelCols(EBR_NUM_ENERGY_COMPONENTS), model_acc2(),
      hgCache(), ownedRegionTag(-1), ownedNodes(),
      regionTags(), numRegions(0),
      region_acc(),
      cacheValid(false), cachedNumElements(-1), cachedNumNodes(-1)
{

}


EnergyBalanceRecorder::EnergyBalanceRecorder(
    Domain &theDomain_,
    OPS_Stream &theOutputHandler_,
    bool echoTimeFlag_,
    const ID &regions,
    bool v2Layout_,
    int ownedNodesRegion)
    : Recorder(RECORDER_TAGS_EnergyBalanceRecorder),
      theDomain(&theDomain_), theOutputHandler(&theOutputHandler_),
      response((echoTimeFlag_ ? 1 : 0) + EBR_NUM_ENERGY_COMPONENTS * (1 + regions.Size())),
      echoTimeFlag(echoTimeFlag_), initializationDone(false), firstRecord(true),
      time_last(0.),
      model_acc(),
      v2Layout(v2Layout_), chInject(false), chLnvd(false),
      chModal(false), chHourglass(false),
      numModelCols(EBR_NUM_ENERGY_COMPONENTS), model_acc2(),
      hgCache(), ownedRegionTag(ownedNodesRegion), ownedNodes(),
      regionTags(regions), numRegions(regions.Size()),
      region_acc((size_t)regions.Size()),
      cacheValid(false), cachedNumElements(-1), cachedNumNodes(-1)
{
    response.Zero();
}


EnergyBalanceRecorder::~EnergyBalanceRecorder()
{
    // the hourglass cache owns its Response objects
    for (size_t i = 0; i < hgCache.size(); ++i)
        delete hgCache[i].res;
    hgCache.clear();

    if (theOutputHandler != 0) {
        if (initializationDone)
            theOutputHandler->endTag(); // closes the "Data" tag opened in initialize()
        delete theOutputHandler;
    }
}


void
EnergyBalanceRecorder::buildCache(void)
{
    elemRegions.clear();
    nodeRegions.clear();

    // region membership: element/node tag -> list of region indices
    for (int r = 0; r < numRegions; ++r) {
        MeshRegion *region = theDomain->getRegion(regionTags(r));
        if (region == 0) {
            opserr << "EnergyBalanceRecorder - WARNING region " << regionTags(r)
                   << " not found in the domain; its energy columns will be zero\n";
            continue;
        }
        const ID &elems = region->getElements();
        for (int i = 0; i < elems.Size(); ++i) {
            if (theDomain->getElement(elems(i)) != 0)
                elemRegions[elems(i)].push_back(r);
        }
        const ID &rnodes = region->getNodes();
        for (int i = 0; i < rnodes.Size(); ++i) {
            if (theDomain->getNode(rnodes(i)) != 0)
                nodeRegions[rnodes(i)].push_back(r);
        }
    }

    // -v2: probe every element for a "hourglassEnergy" response; the hits
    // form the E_hg cache. The probe stream is a DummyStream so the
    // setResponse tag output never reaches the real output.
    for (size_t i = 0; i < hgCache.size(); ++i)
        delete hgCache[i].res;
    hgCache.clear();
    if (v2Layout) {
        static DummyStream probeStream;
        static const char *hgArgv[1] = { "hourglassEnergy" };
        Element *ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0) {
            // skip the ShadowSubdomains of a PartitionedDomain (no node
            // pointers; a setResponse would go to the remote process)
            if (ele->isSubdomain())
                continue;
            Response *hgRes = ele->setResponse(hgArgv, 1, probeStream);
            if (hgRes != 0) {
                HgEntry e = { ele->getTag(), ele, hgRes };
                hgCache.push_back(e);
            }
        }
    }

    // structural staleness sentinels, compared with the live domain in every
    // record() (recorders are not notified of domain changes)
    cachedNumElements = theDomain->getNumElements();
    cachedNumNodes = theDomain->getNumNodes();

    // -ownedNodes: resolve the region into a node set. A missing region
    // warns and owns no nodes, so the nodal columns go to zero instead of
    // silently counting shared nodes several times.
    ownedNodes.clear();
    if (ownedRegionTag >= 0) {
        MeshRegion *owned = theDomain->getRegion(ownedRegionTag);
        if (owned == 0) {
            opserr << "EnergyBalanceRecorder - WARNING -ownedNodes region "
                   << ownedRegionTag << " not found; no nodal terms will be "
                      "counted\n";
        }
        else {
            const ID &onodes = owned->getNodes();
            for (int i = 0; i < onodes.Size(); ++i) {
                if (theDomain->getNode(onodes(i)) != 0)
                    ownedNodes.insert(onodes(i));
            }
        }
    }

    cacheValid = true;
}


int
EnergyBalanceRecorder::record(int commitTag, double timeStamp)
{
    if (theDomain == 0)
        return 0;

    if (theOutputHandler == 0) {
        opserr << "EnergyBalanceRecorder::record() - no output handler has been set\n";
        return -1;
    }

    if (initializationDone == false) {
        if (this->initialize() != 0) {
            opserr << "EnergyBalanceRecorder::record() - failed in initialize()\n";
            return -1;
        }
    }

    // Removal and addition safety. Recorders are not notified of domain
    // changes (Domain::hasDomainChanged() is consumed by the Analysis), so
    // re-validate structurally on every record: a change in the number of
    // elements or nodes, or a hourglass entry whose tag no longer resolves
    // to the same Element*, forces a rebuild before getResponse() is called
    // on a cached Response.
    if (cacheValid == true) {
        if (theDomain->getNumElements() != cachedNumElements ||
            theDomain->getNumNodes() != cachedNumNodes)
            cacheValid = false;
        else
            for (size_t i = 0; i < hgCache.size(); ++i)
                if (theDomain->getElement(hgCache[i].tag) != hgCache[i].ele) {
                    cacheValid = false;
                    break;
                }
    }

    if (cacheValid == false)
        this->buildCache();

    const double dT = timeStamp - time_last;

    int timeOffset = 0;
    if (echoTimeFlag == true) {
        timeOffset = 1;
        response(0) = timeStamp;
    }

    //
    // single sweep: compute each element/node contribution once, accumulate
    // into the whole-model total and bin into the regions it belongs to.
    //
    double g_ke = 0., g_internal_rate = 0., g_damping_rate = 0., g_unbalanced_rate = 0.;
    // -v2 element/nodal split sums
    double v2_ke_ele = 0., v2_ke_nod = 0., v2_dw_ele = 0., v2_dw_nod = 0.;
    std::vector<double> step_region_ke((size_t)numRegions, 0.0);
    std::vector<double> step_region_internal_rate((size_t)numRegions, 0.0);
    std::vector<double> step_region_damping_rate((size_t)numRegions, 0.0);
    std::vector<double> step_region_unbalanced_rate((size_t)numRegions, 0.0);

    {
        Element *ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0) {
            ebkernel::ElementEnergyTerms t;
            ebkernel::elementEnergyTerms(ele, t);
            g_ke += t.ke; g_internal_rate += t.ie; g_damping_rate += t.dw;
            if (v2Layout) { v2_ke_ele += t.ke; v2_dw_ele += t.dw; }

            if (numRegions > 0) {
                std::unordered_map<int, std::vector<int> >::const_iterator it =
                    elemRegions.find(ele->getTag());
                if (it != elemRegions.end())
                    for (size_t k = 0; k < it->second.size(); ++k) {
                        const int r = it->second[k];
                        step_region_ke[(size_t)r] += t.ke;
                        step_region_internal_rate[(size_t)r] += t.ie;
                        step_region_damping_rate[(size_t)r] += t.dw;
                    }
            }
        }
    }

    {
        Node *node;
        NodeIter &nodes = theDomain->getNodes();
        while ((node = nodes()) != 0) {
            // -ownedNodes gate (both layouts)
            if (ownedRegionTag >= 0 &&
                ownedNodes.find(node->getTag()) == ownedNodes.end())
                continue;
            double ke = 0., dw = 0., ulw = 0.;
            ebkernel::addNodeEnergy(node, ke, dw, ulw);
            g_ke += ke; g_damping_rate += dw; g_unbalanced_rate += ulw;
            if (v2Layout) { v2_ke_nod += ke; v2_dw_nod += dw; }

            if (numRegions > 0) {
                std::unordered_map<int, std::vector<int> >::const_iterator it =
                    nodeRegions.find(node->getTag());
                if (it != nodeRegions.end())
                    for (size_t k = 0; k < it->second.size(); ++k) {
                        const int r = it->second[k];
                        step_region_ke[(size_t)r] += ke;
                        step_region_damping_rate[(size_t)r] += dw;
                        step_region_unbalanced_rate[(size_t)r] += ulw;
                    }
            }
        }
    }

    //
    // integrate the work rates and close the balance; fill the response:
    // whole-model block, then one block per region.
    //
    if (v2Layout) {
        // channel totals are cumulative registry reads; the accumulator takes
        // its baseline on the first record so a new recorder starts at zero
        const Ladruno::EnergyChannelRegistry &reg =
            Ladruno::EnergyChannelRegistry::instance();

        // hourglass-stabilization energy: the elements accumulate at commit
        // cadence and the recorder samples the running scalar
        double hg_total = 0.0;
        for (size_t i = 0; i < hgCache.size(); ++i) {
            Response *hgRes = hgCache[i].res;
            if (hgRes->getResponse() >= 0) {
                const Information &info = hgRes->getInformation();
                if (info.theVector != 0 && info.theVector->Size() > 0)
                    hg_total += (*info.theVector)(0);
            }
        }

        double out[ebkernel::NUM_V2_SLOTS];
        model_acc2.step(dT, firstRecord,
                        v2_ke_ele, v2_ke_nod,
                        g_internal_rate, v2_dw_ele, v2_dw_nod,
                        g_unbalanced_rate,
                        reg.energy(Ladruno::EnergyChannelRegistry::ABSORB_LEAK),
                        reg.energy(Ladruno::EnergyChannelRegistry::LNVD_WORK),
                        reg.energy(Ladruno::EnergyChannelRegistry::MODAL_WORK),
                        hg_total,
                        out);
        int col = timeOffset;
        response(col++) = out[ebkernel::V2_KE_ELE];
        response(col++) = out[ebkernel::V2_KE_NOD];
        response(col++) = out[ebkernel::V2_IE];
        response(col++) = out[ebkernel::V2_DW_ELE];
        response(col++) = out[ebkernel::V2_DW_NOD];
        response(col++) = out[ebkernel::V2_ULW];
        if (chInject)
            response(col++) = out[ebkernel::V2_E_INJECT];
        if (chLnvd)
            response(col++) = out[ebkernel::V2_E_LNVD];
        if (chModal)
            response(col++) = out[ebkernel::V2_E_MODAL];
        if (chHourglass)
            response(col++) = out[ebkernel::V2_E_HG];
        response(col++) = out[ebkernel::V2_RES];
        response(col++) = out[ebkernel::V2_ERR];
    }
    else {
        double out[ebkernel::NUM_ENERGY_COMPONENTS];
        model_acc.step(dT, firstRecord,
                       g_ke, g_internal_rate, g_damping_rate, g_unbalanced_rate, out);
        int col = timeOffset;
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            response(col + c) = out[c];
    }

    for (int r = 0; r < numRegions; ++r) {
        double out[ebkernel::NUM_ENERGY_COMPONENTS];
        region_acc[(size_t)r].step(dT, firstRecord,
                                   step_region_ke[(size_t)r],
                                   step_region_internal_rate[(size_t)r],
                                   step_region_damping_rate[(size_t)r],
                                   step_region_unbalanced_rate[(size_t)r], out);
        // region blocks keep the six columns in both layouts; in -v2 the
        // model block may be wider, so offset by its actual width
        const int col = timeOffset + numModelCols + EBR_NUM_ENERGY_COMPONENTS * r;
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            response(col + c) = out[c];
    }

    firstRecord = false;

    theOutputHandler->write(response);

    time_last = timeStamp;

    return 0;
}


int
EnergyBalanceRecorder::setDomain(Domain &theDom)
{
    theDomain = &theDom;
    time_last = theDomain->getCurrentTime();
    firstRecord = true;
    cacheValid = false;
    return 0;
}


int
EnergyBalanceRecorder::domainChanged(void)
{
    // region membership may have changed: rebuild lazily
    cacheValid = false;
    return 0;
}


int
EnergyBalanceRecorder::sendSelf(int commitTag, Channel &theChannel)
{
    // each process builds its own recorder through OPS_EnergyBalanceRecorder
    return 0;
}


int
EnergyBalanceRecorder::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker)
{
    return 0;
}


int
EnergyBalanceRecorder::initialize(void)
{
    if (theDomain == 0) {
        opserr << "EnergyBalanceRecorder::initialize() - the domain has not been set\n";
        return -1;
    }

    // -v2: fix the channel columns from the registry's declared set.
    // Producers declare at construction or setDomain, so anything built
    // before analyze is visible here. E_hg is decided from the response
    // probe, which must run before the widths are fixed.
    if (v2Layout) {
        if (cacheValid == false)
            this->buildCache();
        const Ladruno::EnergyChannelRegistry &reg =
            Ladruno::EnergyChannelRegistry::instance();
        chInject = reg.declared(Ladruno::EnergyChannelRegistry::ABSORB_LEAK);
        chLnvd   = reg.declared(Ladruno::EnergyChannelRegistry::LNVD_WORK);
        chModal  = reg.declared(Ladruno::EnergyChannelRegistry::MODAL_WORK);
        chHourglass = !hgCache.empty();
        numModelCols = 6 + (chInject ? 1 : 0) + (chLnvd ? 1 : 0)
                     + (chModal ? 1 : 0) + (chHourglass ? 1 : 0) + 2;
    }
    else {
        numModelCols = EBR_NUM_ENERGY_COMPONENTS;
    }

    // size and zero the response and all per-scope accumulators
    const int timeOffset = (echoTimeFlag == true) ? 1 : 0;
    response.resize(timeOffset + numModelCols
                    + EBR_NUM_ENERGY_COMPONENTS * numRegions);
    response.Zero();

    model_acc.reset();
    model_acc2.reset();
    region_acc.assign((size_t)numRegions, ebkernel::EnergyAccumulator());

    // column metadata (used by the XML and binary streams; text streams
    // ignore tags, the layout is echoed once below)
    static const char *comp[EBR_NUM_ENERGY_COMPONENTS] =
        { "KE", "IE", "DW", "ULW", "RES", "ERR%" };

    if (echoTimeFlag == true) {
        theOutputHandler->tag("TimeOutput");
        theOutputHandler->tag("ResponseType", "time");
        theOutputHandler->endTag();
    }

    theOutputHandler->tag("EnergyOutput");
    theOutputHandler->attr("region", "model");
    if (v2Layout) {
        theOutputHandler->tag("ResponseType", "KE_ele");
        theOutputHandler->tag("ResponseType", "KE_nod");
        theOutputHandler->tag("ResponseType", "IE");
        theOutputHandler->tag("ResponseType", "DW_ele");
        theOutputHandler->tag("ResponseType", "DW_nod");
        theOutputHandler->tag("ResponseType", "ULW");
        if (chInject) theOutputHandler->tag("ResponseType", "E_inject");
        if (chLnvd)   theOutputHandler->tag("ResponseType", "E_lnvd");
        if (chModal)  theOutputHandler->tag("ResponseType", "E_modal");
        if (chHourglass) theOutputHandler->tag("ResponseType", "E_hg");
        theOutputHandler->tag("ResponseType", "RES");
        theOutputHandler->tag("ResponseType", "ERR%");
    } else {
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            theOutputHandler->tag("ResponseType", comp[c]);
    }
    theOutputHandler->endTag();

    for (int r = 0; r < numRegions; ++r) {
        theOutputHandler->tag("EnergyOutput");
        theOutputHandler->attr("region", regionTags(r));
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            theOutputHandler->tag("ResponseType", comp[c]);
        theOutputHandler->endTag();
    }

    theOutputHandler->tag("Data");

    // echo the column layout once for users of the plain-text output
    opserr << "EnergyBalanceRecorder: columns =";
    if (echoTimeFlag == true) opserr << " time";
    if (v2Layout) {
        opserr << " [model: KE_ele KE_nod IE DW_ele DW_nod ULW";
        if (chInject) opserr << " E_inject";
        if (chLnvd)   opserr << " E_lnvd";
        if (chModal)  opserr << " E_modal";
        if (chHourglass) opserr << " E_hg";
        opserr << " RES ERR%]";
    } else {
        opserr << " [model: KE IE DW ULW RES ERR%]";
    }
    for (int r = 0; r < numRegions; ++r)
        opserr << " [region " << regionTags(r) << ": KE IE DW ULW RES ERR%]";
    opserr << endln;

    initializationDone = true;

    return 0;
}


int EnergyBalanceRecorder::flush(void) {
    if (theOutputHandler != 0)
        return theOutputHandler->flush();
    return 0;
}
