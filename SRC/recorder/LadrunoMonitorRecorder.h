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
// Description: recorder Monitor (LadrunoMonitorRecorder), a live analysis
// monitor. It is an ordinary Recorder, so record() runs at every commit; it
// streams selected nodal scalars (disp, vel, accel, reaction) to the SWMR HDF5
// sink of LadrunoMonitorSink.h, which another process can read while the
// analysis runs. Two gates bound the data rate: -every K keeps one commit in
// K, -hz caps the wall-clock rate. Nodal scalar time histories only; in an
// interpreter-per-rank run each rank writes its own "<stem>.part-<rank><ext>".
//
// References:
//   The HDF Group, HDF5 Single-Writer/Multiple-Reader (SWMR) User's Guide,
//     https://docs.hdfgroup.org/hdf5/develop/_s_w_m_r.html

#ifndef LadrunoMonitorRecorder_h
#define LadrunoMonitorRecorder_h

#include <Recorder.h>
#include <ID.h>

#include <chrono>
#include <string>
#include <vector>

class Domain;
class Node;

namespace ladruno { class MonitorSink; }

class LadrunoMonitorRecorder : public Recorder {
public:
    LadrunoMonitorRecorder();
    LadrunoMonitorRecorder(const ID &theDofs,
                           const ID &theNodalTags,
                           const char *respName,
                           Domain &theDomain,
                           const char *filename,
                           int everyK,
                           double hz);
    ~LadrunoMonitorRecorder();

    int record(int commitTag, double timeStamp);
    int flush();

    int restart(void);
    int domainChanged(void);
    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private:
    int initialize(void);

    ID theDofs;            // 0-based dof indices
    ID theNodalTags;
    std::vector<int> theNodeTags;   // tags of resolved nodes that became columns;
                                    // re-looked-up each emit so a mid-run node
                                    // removal cannot dangle (see record()).
    int dataFlag;          // 0 disp, 1 vel, 2 accel; 7/8/9 reaction
                           // (static / inclInertia / +Rayleigh)
    std::string respName;
    std::string filename;
    Domain *theDomain;

    // rate gates
    int everyK;                 // >=1
    double minIntervalSec;      // from -hz; <=0 => no wall throttle
    long long recordCallCount;  // counts every record() call (for -every)
    std::chrono::steady_clock::time_point lastEmit;
    bool haveEmitted;

    ladruno::MonitorSink *sink; // owned
    bool initialized;
    bool initFailed;
    std::vector<double> frameBuf;
};

#endif
