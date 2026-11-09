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
//   Chopra, A.K. (2017). Dynamics of Structures, 5th ed. Pearson.
//   ASCE/SEI 7-22, Minimum Design Loads and Associated Criteria for
//     Buildings and Other Structures, Sec. 12.9.1 (modal response spectrum
//     analysis: each quantity is combined from its own per-mode values).
//
// Description: responseSpectrumPerMode, the per-mode response-spectrum
// resultants of every requested mode and direction in one call.
//
//   responseSpectrumPerMode (-ts $tag | -Tn $T.. -Sa $Sa.. | -fn $f.. -Sa $Sa..)
//       -dir $d1 <$d2 ..> <-scale $s> <-modes $m1 $m2 ..>
//       <-reactions <-about $x $y <$z>>>
//       <-reactionSet $name $node.. <-about ..>> ...
//       <-eleGroup $name -ele $ele.. <-cut $node..> <-about ..>> ...
//       <-inertia $name $node.. <-about ..>> ...
//
// Per mode k and direction d the amplitude is the one responseSpectrumAnalysis
// -mode k imposes:
//     u_k,d = phi_k s_k Gamma_k,d Sa(T_k) / lambda_k   about the committed state
// (s_k the eigenvector scale factor of modalProperties).
// Instead of updating the whole domain and sweeping every element for the
// reactions, the command updates ONLY the elements a requested quantity needs
// (the ones attached to a supported node, and the named element groups), on
// only their nodes, and returns
//     reactions        sum of the reactions of every SP-constrained node,
//                      on its SP-constrained DOFs only
//     reaction:<name>  the same over a named set of supported nodes
//     ele:<name>       sum of the resisting forces (eleForce) of the group's
//                      elements at its cut nodes (all nodes without -cut)
//     inertia:<name>   sum of the modal inertia forces M phi s Gamma Sa
// each by component (Fx Fy Fz Mx My Mz in 3D, Fx Fy Mz in 2D), as the
// INCREMENT over the committed state, flat [dir][mode][component]. -about adds
// r x F about the point to the moments of the group it follows.
// The domain is left at its committed state. Each per-mode quantity can then
// be combined by the rule of choice (SRSS, CQC, ...), as ASCE 7-22 12.9.1
// requires: combining displacements first and recovering forces from them is
// not the combination of the forces.

#include <ResponseSpectrumAnalysis.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <DomainModalProperties.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <TimeSeries.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <elementAPI.h>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

const char* const CMD = "responseSpectrumPerMode";

// read a run of numbers; stops (without consuming) at the first non-number
// (the responseSpectrumAnalysis list idiom: a failed read may consume)
int readOne(int* v) { int one = 1; return OPS_GetIntInput(&one, v); }
int readOne(double* v) { int one = 1; return OPS_GetDoubleInput(&one, v); }

void parseOne(const char* s, char** end, int* v) { *v = (int)strtol(s, end, 10); }
void parseOne(const char* s, char** end, double* v) { *v = strtod(s, end); }

// a whole list as ONE argument (classic Tcl `-Tn $list` unexpanded, as
// responseSpectrumAnalysis accepts): whitespace-separated numbers, all or none
template <class T>
bool parseListString(const std::string& text, std::vector<T>& out)
{
    std::vector<T> got;
    const char* p = text.c_str();
    while (true) {
        while (std::isspace((unsigned char)*p)) ++p;
        if (*p == 0) break;
        char* end = 0;
        T v;
        parseOne(p, &end, &v);
        if (end == p || (*end != 0 && !std::isspace((unsigned char)*end)))
            return false;
        got.push_back(v);
        p = end;
    }
    if (got.empty()) return false;
    out.insert(out.end(), got.begin(), got.end());
    return true;
}

template <class T>
void readList(std::vector<T>& out)
{
    const size_t start = out.size();
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const int before = OPS_GetNumRemainingInputArgs();
        T v;
        if (readOne(&v) < 0) {
            if (OPS_GetNumRemainingInputArgs() < before)
                OPS_ResetCurrentInputArg(-1);
            break;
        }
        out.push_back(v);
    }
    // nothing read: try the next argument as one list string
    if (out.size() == start && OPS_GetNumRemainingInputArgs() > 0) {
        // OPS_GetString, not OPS_GetStringFromAll: no buffer to truncate a
        // long node list (a number never gets here, readOne took it)
        const char* s = OPS_GetString();
        if (s == 0 || !parseListString(std::string(s), out))
            OPS_ResetCurrentInputArg(-1);
    }
}

struct Group {
    std::string key;                 // output key, e.g. "ele:pier1"
    std::vector<int> nodes;          // reactionSet / inertia / eleGroup -cut
    std::vector<int> eles;           // eleGroup
    bool hasCut = false;
    bool about = false;
    double p[3] = {0.0, 0.0, 0.0};
};

// add the nodal force block f (size ndfNode) of node nd to the resultant r
// (ncomp components), with the moment of f about p when `about`
void addBlock(double* r, int ncomp, const double* f, int ndfNode, Node* nd,
              const Group& g)
{
    const int nd3 = (ncomp == 6 && ndfNode == 4) ? 3 : ndfNode; // 3D U-P: skip p
    const int n = std::min(nd3, ncomp);
    for (int i = 0; i < n; ++i)
        r[i] += f[i];
    if (!g.about) return;
    const Vector& x = nd->getCrds();
    if (ncomp == 6) {
        const double rx = x(0) - g.p[0];
        const double ry = (x.Size() > 1 ? x(1) : 0.0) - g.p[1];
        const double rz = (x.Size() > 2 ? x(2) : 0.0) - g.p[2];
        const double fx = n > 0 ? f[0] : 0.0, fy = n > 1 ? f[1] : 0.0, fz = n > 2 ? f[2] : 0.0;
        r[3] += ry * fz - rz * fy;
        r[4] += rz * fx - rx * fz;
        r[5] += rx * fy - ry * fx;
    }
    else if (ncomp == 3) {
        const double rx = x(0) - g.p[0];
        const double ry = (x.Size() > 1 ? x(1) : 0.0) - g.p[1];
        const double fx = n > 0 ? f[0] : 0.0, fy = n > 1 ? f[1] : 0.0;
        r[2] += rx * fy - ry * fx;
    }
}

// one node block of a touched element: element force offset + node
struct Contrib {
    int ele;      // index into the touched-element list
    int off;      // offset of the node block in the element force vector
    int ndf;      // node DOFs
    Node* node;
};

} // namespace

int
OPS_LadrunoResponseSpectrumPerMode(void)
{
    AnalysisModel** modelPtr = OPS_GetAnalysisModel();
    AnalysisModel* theModel = modelPtr ? *modelPtr : 0;
    if (theModel == 0 || theModel->getDomainPtr() == 0) {
        opserr << CMD << ": no AnalysisModel - define an analysis and run eigen and modalProperties first\n";
        return -1;
    }
    Domain* domain = theModel->getDomainPtr();
    {
        // a partitioned domain holds Subdomain elements: not supported
        ElementIter& it = domain->getElements();
        Element* e;
        while ((e = it()) != 0)
            if (e->isSubdomain()) {
                opserr << CMD << ": a partitioned (subdomain) model is not supported\n";
                return -1;
            }
    }
    DomainModalProperties mp;
    if (domain->getModalProperties(mp) < 0) {
        opserr << CMD << ": eigen and modalProperties have not been called\n";
        return -1;
    }
    const int ncomp = mp.totalMass().Size();       // 3 (2D) or 6 (3D)

    // ---- parse ---------------------------------------------------------
    TimeSeries* ts = 0;
    std::vector<double> Tn, Sa;
    bool useFn = false;
    double scale = 1.0;
    std::vector<int> dirs, modes;
    bool haveModes = false;
    std::deque<Group> reactionSets, eleGroups, inertias;
    bool total = false;
    Group totalGroup;
    totalGroup.key = "reactions";
    Group* last = 0;
    std::set<std::string> keys;
    char tok[256];

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetStringFromAll(tok, sizeof(tok));
        if (opt == 0) {
            opserr << CMD << ": unreadable argument\n";
            return -1;
        }
        const std::string o(opt);
        if (o == "-ts") {
            int tag;
            int one = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetIntInput(&one, &tag) < 0) {
                opserr << CMD << ": -ts requires a timeSeries tag\n";
                return -1;
            }
            ts = OPS_getTimeSeries(tag);
            if (ts == 0) {
                opserr << CMD << ": timeSeries " << tag << " not found\n";
                return -1;
            }
        }
        else if (o == "-Tn" || o == "-fn") {
            Tn.clear();
            readList(Tn);
            useFn = (o == "-fn");
        }
        else if (o == "-Sa") {
            Sa.clear();
            readList(Sa);
        }
        else if (o == "-scale") {
            int one = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&one, &scale) < 0) {
                opserr << CMD << ": -scale requires a value\n";
                return -1;
            }
        }
        else if (o == "-dir") {
            readList(dirs);
        }
        else if (o == "-modes") {
            readList(modes);
            haveModes = true;
        }
        else if (o == "-reactions") {
            total = true;
            last = &totalGroup;
        }
        else if (o == "-reactionSet" || o == "-inertia" || o == "-eleGroup") {
            const char* name = OPS_GetNumRemainingInputArgs() > 0 ? OPS_GetStringFromAll(tok, sizeof(tok)) : 0;
            if (name == 0 || name[0] == '\0' || name[0] == '-') {
                opserr << CMD << ": " << o.c_str() << " requires a name\n";
                return -1;
            }
            Group g;
            std::deque<Group>* dst;
            if (o == "-reactionSet") { g.key = std::string("reaction:") + name; dst = &reactionSets; }
            else if (o == "-inertia") { g.key = std::string("inertia:") + name; dst = &inertias; }
            else { g.key = std::string("ele:") + name; dst = &eleGroups; }
            if (!keys.insert(g.key).second) {
                opserr << CMD << ": duplicate group " << g.key.c_str() << "\n";
                return -1;
            }
            if (o == "-eleGroup") {
                const char* e = OPS_GetNumRemainingInputArgs() > 0 ? OPS_GetStringFromAll(tok, sizeof(tok)) : 0;
                if (e == 0 || strcmp(e, "-ele") != 0) {
                    opserr << CMD << ": -eleGroup " << name << " requires -ele $tag ..\n";
                    return -1;
                }
                readList(g.eles);
                if (OPS_GetNumRemainingInputArgs() > 0) {
                    const char* c = OPS_GetStringFromAll(tok, sizeof(tok));
                    if (c != 0 && strcmp(c, "-cut") == 0) {
                        g.hasCut = true;
                        readList(g.nodes);
                        if (g.nodes.empty()) {
                            opserr << CMD << ": -eleGroup " << name << " -cut lists no node\n";
                            return -1;
                        }
                    }
                    else
                        OPS_ResetCurrentInputArg(-1);
                }
                if (g.eles.empty()) {
                    opserr << CMD << ": -eleGroup " << name << " lists no element\n";
                    return -1;
                }
            }
            else {
                readList(g.nodes);
                if (g.nodes.empty()) {
                    opserr << CMD << ": " << o.c_str() << " " << name << " lists no node\n";
                    return -1;
                }
            }
            dst->push_back(g);
            last = &dst->back();
        }
        else if (o == "-about") {
            if (last == 0) {
                opserr << CMD << ": -about must follow -reactions, -reactionSet, -eleGroup or -inertia\n";
                return -1;
            }
            std::vector<double> p;
            readList(p);
            const int ndm = OPS_GetNDM();
            if ((int)p.size() != ndm) {
                opserr << CMD << ": -about takes " << ndm << " coordinates (got " << (int)p.size() << ")\n";
                return -1;
            }
            last->about = true;
            for (int i = 0; i < ndm && i < 3; ++i) last->p[i] = p[i];
        }
        else {
            opserr << CMD << ": unknown option '" << opt << "'\n";
            return -1;
        }
    }

    // ---- validate --------------------------------------------------------
    if (ts == 0 && (Tn.empty() || Sa.empty())) {
        opserr << CMD << ": give the spectrum as -ts $tag or as -Tn (or -fn) and -Sa lists\n";
        return -1;
    }
    if (ts != 0 && (!Tn.empty() || !Sa.empty())) {
        opserr << CMD << ": -ts and -Tn/-Sa are mutually exclusive\n";
        return -1;
    }
    if (ts == 0) {
        if (useFn)
            for (double& v : Tn) if (v != 0.0) v = 1.0 / v;
        if (Tn.size() != Sa.size()) {
            opserr << CMD << ": -Tn and -Sa have different lengths (" << (int)Tn.size()
                   << " != " << (int)Sa.size() << ")\n";
            return -1;
        }
        for (size_t i = 0; i < Tn.size(); ++i) {
            if (Tn[i] < 0.0 || (i > 0 && Tn[i] <= Tn[i - 1])) {
                opserr << CMD << ": -Tn values must be non-negative and increasing\n";
                return -1;
            }
            if (Sa[i] < 0.0) {
                opserr << CMD << ": -Sa values must be non-negative\n";
                return -1;
            }
        }
    }
    if (dirs.empty()) {
        opserr << CMD << ": -dir $d .. is required\n";
        return -1;
    }
    for (size_t i = 0; i < dirs.size(); ++i) {
        if (dirs[i] < 1 || dirs[i] > ncomp) {
            opserr << CMD << ": direction " << dirs[i] << " out of range 1-" << ncomp << "\n";
            return -1;
        }
        for (size_t j = 0; j < i; ++j)
            if (dirs[j] == dirs[i]) {
                opserr << CMD << ": direction " << dirs[i] << " given twice\n";
                return -1;
            }
    }
    if (!total && reactionSets.empty() && eleGroups.empty() && inertias.empty()) {
        opserr << CMD << ": nothing requested (-reactions, -reactionSet, -eleGroup or -inertia)\n";
        return -1;
    }

    // the responseSpectrumAnalysis consistency check: modalProperties of THIS
    // mode set
    ResponseSpectrumAnalysis rsa0(theModel, ts, Tn, Sa, dirs[0], scale);
    if (rsa0.check() < 0)
        return -1;
    const int numEigen = domain->getEigenvalues().Size();
    if (!haveModes)
        for (int k = 1; k <= numEigen; ++k) modes.push_back(k);
    if (modes.empty()) {
        opserr << CMD << ": -modes lists no mode\n";
        return -1;
    }
    {
        std::set<int> seen;
        for (int k : modes) {
            if (k < 1 || k > numEigen) {
                opserr << CMD << ": mode " << k << " out of range 1-" << numEigen << "\n";
                return -1;
            }
            if (!seen.insert(k).second) {
                opserr << CMD << ": mode " << k << " given twice\n";
                return -1;
            }
            if (!(mp.eigenvalues()(k - 1) > 0.0)) {
                opserr << CMD << ": mode " << k << " has a non-positive eigenvalue ("
                       << mp.eigenvalues()(k - 1) << "); select modes with -modes\n";
                return -1;
            }
        }
    }

    // ---- supported nodes and the elements a quantity needs ------------------
    std::map<int, int> supIdx;            // node tag -> index into sup
    std::vector<Node*> sup;
    std::vector<char> supMask;            // [sup][6]: DOF carries an SP
    {
        SP_ConstraintIter& sps = domain->getDomainAndLoadPatternSPs();
        SP_Constraint* sp;
        while ((sp = sps()) != 0) {
            const int nt = sp->getNodeTag();
            Node* nd = domain->getNode(nt);
            if (nd == 0) continue;
            auto it = supIdx.find(nt);
            int s;
            if (it == supIdx.end()) {
                s = (int)sup.size();
                supIdx[nt] = s;
                sup.push_back(nd);
                supMask.resize(supMask.size() + 6, 0);
            }
            else
                s = it->second;
            const int dof = sp->getDOF_Number();
            if (dof >= 0 && dof < 6) supMask[s * 6 + dof] = 1;
        }
    }
    for (const Group& g : reactionSets)
        for (int nt : g.nodes) {
            if (domain->getNode(nt) == 0) {
                opserr << CMD << ": " << g.key.c_str() << ": node " << nt << " does not exist\n";
                return -1;
            }
            if (!supIdx.count(nt)) {
                opserr << CMD << ": " << g.key.c_str() << ": node " << nt << " has no SP constraint (not a supported node)\n";
                return -1;
            }
        }
    for (const Group& g : inertias)
        for (int nt : g.nodes)
            if (domain->getNode(nt) == 0) {
                opserr << CMD << ": " << g.key.c_str() << ": node " << nt << " does not exist\n";
                return -1;
            }
    const bool needReactions = total || !reactionSets.empty();

    std::vector<Element*> touched;
    std::map<int, int> touchedIdx;        // element tag -> index into touched
    std::vector<Contrib> supContrib;      // element blocks at supported nodes
    std::vector<int> supOf;               // supContrib -> sup index
    auto touch = [&](Element* e) -> int {
        auto it = touchedIdx.find(e->getTag());
        if (it != touchedIdx.end()) return it->second;
        const int idx = (int)touched.size();
        touched.push_back(e);
        touchedIdx[e->getTag()] = idx;
        return idx;
    };
    if (needReactions && !sup.empty()) {
        // one pass over the element connectivity (no state determination)
        ElementIter& it = domain->getElements();
        Element* e;
        while ((e = it()) != 0) {
            const int nn = e->getNumExternalNodes();
            Node** np = e->getNodePtrs();
            int off = 0, idx = -1;
            for (int a = 0; a < nn; ++a) {
                const int nd = np[a]->getNumberDOF();
                auto s = supIdx.find(np[a]->getTag());
                if (s != supIdx.end()) {
                    if (idx < 0) idx = touch(e);
                    supContrib.push_back(Contrib{idx, off, nd, np[a]});
                    supOf.push_back(s->second);
                }
                off += nd;
            }
        }
    }
    std::vector<std::vector<Contrib>> eleContrib(eleGroups.size());
    for (size_t gi = 0; gi < eleGroups.size(); ++gi) {
        const Group& g = eleGroups[gi];
        std::set<int> cut(g.nodes.begin(), g.nodes.end()), cutHit;
        std::set<int> seen;
        for (int et : g.eles) {
            Element* e = domain->getElement(et);
            if (e == 0) {
                opserr << CMD << ": " << g.key.c_str() << ": element " << et << " does not exist\n";
                return -1;
            }
            if (!seen.insert(et).second) {
                opserr << CMD << ": " << g.key.c_str() << ": element " << et << " given twice\n";
                return -1;
            }
            const int idx = touch(e);
            const int nn = e->getNumExternalNodes();
            Node** np = e->getNodePtrs();
            int off = 0;
            for (int a = 0; a < nn; ++a) {
                const int nd = np[a]->getNumberDOF();
                if (!g.hasCut || cut.count(np[a]->getTag())) {
                    eleContrib[gi].push_back(Contrib{idx, off, nd, np[a]});
                    cutHit.insert(np[a]->getTag());
                }
                off += nd;
            }
        }
        if (g.hasCut)
            for (int nt : g.nodes)
                if (!cutHit.count(nt)) {
                    opserr << CMD << ": " << g.key.c_str() << ": cut node " << nt
                           << " is on no element of the group\n";
                    return -1;
                }
    }

    // nodes whose trial displacement the touched elements read
    std::vector<Node*> tNodes;
    {
        std::set<int> seen;
        for (Element* e : touched) {
            const int nn = e->getNumExternalNodes();
            Node** np = e->getNodePtrs();
            for (int a = 0; a < nn; ++a)
                if (seen.insert(np[a]->getTag()).second) tNodes.push_back(np[a]);
        }
    }
    std::vector<std::vector<Node*>> inNodes(inertias.size());
    for (size_t gi = 0; gi < inertias.size(); ++gi)
        for (int nt : inertias[gi].nodes) inNodes[gi].push_back(domain->getNode(nt));

    // ---- per-mode scalars (the solveMode expressions) ----------------------
    const int nm = (int)modes.size(), ndir = (int)dirs.size();
    std::vector<double> period(nm), omega(nm);
    std::vector<double> SaV(ndir * nm), gamma(ndir * nm);
    for (int a = 0; a < nm; ++a) {
        const int k = modes[a] - 1;
        const double lambda = mp.eigenvalues()(k);
        const double w = std::sqrt(lambda);
        const double freq = w / 2.0 / M_PI;
        period[a] = 1.0 / freq;
        omega[a] = w;
        for (int d = 0; d < ndir; ++d) {
            SaV[d * nm + a] = rsa0.getSa(period[a]);
            gamma[d * nm + a] = mp.modalParticipationFactors()(k, dirs[d] - 1);
        }
    }

    // ---- start from the committed state ---------------------------------------
    if (theModel->revertDomainToLastCommit() < 0) {
        opserr << CMD << ": revertDomainToLastCommit failed\n";
        return -1;
    }
    // (the revert updated the whole domain once: F0 is the committed force)
    std::vector<Vector> F0(touched.size()), Ft(touched.size());
    for (size_t i = 0; i < touched.size(); ++i)
        F0[i] = touched[i]->getResistingForce();
    std::vector<double> R0(sup.size() * 6, 0.0), Rt(sup.size() * 6, 0.0);
    for (size_t c = 0; c < supContrib.size(); ++c) {
        const Contrib& b = supContrib[c];
        for (int i = 0; i < b.ndf && i < 6; ++i)
            R0[supOf[c] * 6 + i] += F0[b.ele](b.off + i);
    }
    std::vector<std::vector<double>> G0(eleGroups.size(), std::vector<double>(ncomp, 0.0));
    for (size_t gi = 0; gi < eleGroups.size(); ++gi)
        for (const Contrib& b : eleContrib[gi])
            addBlock(G0[gi].data(), ncomp, &F0[b.ele](b.off), b.ndf, b.node, eleGroups[gi]);

    // ---- outputs, flat [dir][mode][comp] -----------------------------------
    const size_t blk = (size_t)ndir * nm * ncomp;
    std::vector<double> outTotal(total ? blk : 0);
    std::vector<std::vector<double>> outSet(reactionSets.size(), std::vector<double>(blk, 0.0));
    std::vector<std::vector<double>> outEle(eleGroups.size(), std::vector<double>(blk, 0.0));
    std::vector<std::vector<double>> outIn(inertias.size(), std::vector<double>(blk, 0.0));
    std::vector<double> fblk(6), ablk(6), mvec(6);
    int rc = 0;

    for (int d = 0; d < ndir && rc == 0; ++d) {
        for (int a = 0; a < nm && rc == 0; ++a) {
            const int k = modes[a] - 1;
            const double lambda = mp.eigenvalues()(k);
            const double mga = SaV[d * nm + a];
            const double Vscale = mp.eigenVectorScaleFactors()(k);
            const double MPF = gamma[d * nm + a];
            const size_t o = ((size_t)d * nm + a) * ncomp;

            // on the touched objects only: revert and update at the committed
            // state (as Domain::revertToLastCommit does), then
            // trial = committed + u_k (solveMode) and update (updateDomain).
            // Elements that condense internal DOFs incrementally (ASDShellQ4)
            // depend on the previous trial, so the first update is not skipped.
            for (Node* nd : tNodes) nd->revertToLastCommit();
            for (size_t i = 0; i < touched.size() && rc == 0; ++i) {
                touched[i]->revertToLastCommit();
                if (touched[i]->update() < 0) {
                    opserr << CMD << ": element " << touched[i]->getTag() << " failed in update at mode "
                           << modes[a] << "\n";
                    rc = -1;
                }
            }
            if (rc) break;
            for (Node* nd : tNodes) {
                if (nd->getNumEigenvectors() <= k) continue;
                const Matrix& ev = nd->getEigenvectors();
                const int nndf = ev.noRows();
                const Vector& uc = nd->getDisp();
                for (int i = 0; i < std::min(nndf, ncomp); ++i) {
                    if (ncomp == 6 && nndf == 4 && i == 3) continue;
                    const double V = ev(i, k) * Vscale;
                    double u = V * MPF * mga / lambda;
                    u = uc(i) + 1.0 * u;
                    nd->setTrialDisp(u, i);
                }
            }
            for (size_t i = 0; i < touched.size(); ++i) {
                if (touched[i]->update() < 0) {
                    opserr << CMD << ": element " << touched[i]->getTag() << " failed in update at mode "
                           << modes[a] << "\n";
                    rc = -1;
                    break;
                }
                Ft[i] = touched[i]->getResistingForce();
            }
            if (rc) break;

            if (needReactions) {
                std::fill(Rt.begin(), Rt.end(), 0.0);
                for (size_t c = 0; c < supContrib.size(); ++c) {
                    const Contrib& b = supContrib[c];
                    for (int i = 0; i < b.ndf && i < 6; ++i)
                        Rt[supOf[c] * 6 + i] += Ft[b.ele](b.off + i);
                }
                auto addNode = [&](double* r, int s, const Group& g) {
                    const int nd = std::min(sup[s]->getNumberDOF(), 6);
                    // a reaction exists only on an SP-constrained DOF; at a free
                    // DOF the element sum is the mode's out-of-balance force
                    for (int i = 0; i < nd; ++i)
                        fblk[i] = supMask[s * 6 + i] ? Rt[s * 6 + i] - R0[s * 6 + i] : 0.0;
                    addBlock(r, ncomp, fblk.data(), nd, sup[s], g);
                };
                if (total)
                    for (size_t s = 0; s < sup.size(); ++s) addNode(&outTotal[o], (int)s, totalGroup);
                for (size_t gi = 0; gi < reactionSets.size(); ++gi)
                    for (int nt : reactionSets[gi].nodes)
                        addNode(&outSet[gi][o], supIdx[nt], reactionSets[gi]);
            }
            for (size_t gi = 0; gi < eleGroups.size(); ++gi) {
                double* r = &outEle[gi][o];
                for (const Contrib& b : eleContrib[gi])
                    addBlock(r, ncomp, &Ft[b.ele](b.off), b.ndf, b.node, eleGroups[gi]);
                for (int i = 0; i < ncomp; ++i) r[i] -= G0[gi][i];
            }
            // modal inertia M (phi s Gamma Sa) = lambda M u_k (no element work)
            for (size_t gi = 0; gi < inertias.size(); ++gi) {
                double* r = &outIn[gi][o];
                for (Node* nd : inNodes[gi]) {
                    if (nd->getNumEigenvectors() <= k) continue;
                    const Matrix& ev = nd->getEigenvectors();
                    const int nndf = std::min(ev.noRows(), 6);
                    const Matrix& M = nd->getMass();
                    for (int i = 0; i < nndf; ++i) {
                        const bool p = (ncomp == 6 && ev.noRows() == 4 && i == 3) || i >= ncomp;
                        ablk[i] = p ? 0.0 : ev(i, k) * Vscale * MPF * mga;
                    }
                    for (int i = 0; i < nndf; ++i) {
                        double s = 0.0;
                        for (int j = 0; j < nndf; ++j) s += M(i, j) * ablk[j];
                        mvec[i] = s;
                    }
                    addBlock(r, ncomp, mvec.data(), nndf, nd, inertias[gi]);
                }
            }
        }
    }

    // ---- leave the touched objects at the committed state -------------------
    for (Node* nd : tNodes) nd->revertToLastCommit();
    for (Element* e : touched) {
        e->revertToLastCommit();
        e->update();
    }
    if (rc) return rc;

    // ---- publish --------------------------------------------------------------
    std::map<const char*, std::vector<double>> data;
    std::deque<std::string> names;    // stable storage for the dynamic keys
    data["shape"] = {(double)ndir, (double)nm, (double)ncomp};
    for (int dd : dirs) data["dirs"].push_back(dd);
    for (int k : modes) data["modes"].push_back(k);
    data["period"] = period;
    data["omega"] = omega;
    data["Sa"] = SaV;
    data["gamma"] = gamma;
    data["numSupports"] = {(double)sup.size()};
    data["numUpdated"] = {(double)touched.size()};
    if (total) data["reactions"] = outTotal;
    auto put = [&](const Group& g, std::vector<double>& v) {
        names.push_back(g.key);
        data[names.back().c_str()] = v;
    };
    for (size_t gi = 0; gi < reactionSets.size(); ++gi) put(reactionSets[gi], outSet[gi]);
    for (size_t gi = 0; gi < eleGroups.size(); ++gi) put(eleGroups[gi], outEle[gi]);
    for (size_t gi = 0; gi < inertias.size(); ++gi) put(inertias[gi], outIn[gi]);
    if (OPS_SetDoubleDictListOutput(data) < 0) {
        opserr << CMD << ": failed to set the output\n";
        return -1;
    }
    return 0;
}
