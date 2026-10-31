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
// Description: SWMR (single-writer / multiple-reader) HDF5 sink of recorder
// Monitor. One file, written so that a separate process can read it while
// the analysis is still running. Layout:
//
//   /              attrs: FORMAT="ladruno-monitor", FORMAT_VERSION, GENERATOR, units (optional)
//     COLUMNS      [nCols]   variable-length strings (self-describing channels)
//     STEP         [T]       int32, chunked/unlimited  (commitTag per frame)
//     TIME         [T]       float64, chunked/unlimited (pseudo-time per frame)
//     FRAMES       [T,nCols] float64, chunked/unlimited (one row appended per frame)
//
// SWMR rules: the file is created with H5F_LIBVER_LATEST bounds, every
// dataset is chunked, the full object layout (COLUMNS included) exists before
// H5Fstart_swmr_write (no object may be created afterwards, so the channel set
// is fixed at the first record), and each append() calls H5Dflush so a reader
// sees the new frame. After the run the same file is an ordinary results file.
//
// References:
//   The HDF Group, HDF5 Single-Writer/Multiple-Reader (SWMR) User's Guide,
//     https://docs.hdfgroup.org/hdf5/develop/_s_w_m_r.html

#ifndef LadrunoMonitorSink_h
#define LadrunoMonitorSink_h

#include <string>
#include <vector>

#include "hdf5.h"

namespace ladruno {

class MonitorSink {
public:
    // columns: one self-describing label per scalar channel, e.g.
    //   "node5.disp.dof1". units: optional free-form string ("" = omit).
    MonitorSink(const std::string &filename,
                const std::vector<std::string> &columns,
                const std::string &units = "");
    ~MonitorSink();

    // Create the file + datasets and arm SWMR. Returns 0 on success, <0 on
    // failure (after which isOpen() is false and append() is a no-op).
    int open();

    // Append one frame. values must have columns().size() entries. Never
    // blocks on a reader. Returns 0 on success, <0 on a write error.
    int append(int step, double time, const std::vector<double> &values);

    // Flush + close. Safe to call more than once. Called from the dtor, never
    // from a static/atexit path (HDF5 may already be shut down there).
    void close();

    bool isOpen() const { return fileId >= 0; }
    size_t numColumns() const { return columns.size(); }
    long long numFramesWritten() const { return nFrames; }

private:
    MonitorSink(const MonitorSink &) = delete;
    MonitorSink &operator=(const MonitorSink &) = delete;

    int extendAndWrite(hid_t dset, hid_t memType, hsize_t ncols,
                       const void *rowData);

    std::string filename;
    std::vector<std::string> columns;
    std::string units;

    hid_t fileId;     // H5I_INVALID_HID until open()
    hid_t dsStep;
    hid_t dsTime;
    hid_t dsFrames;
    long long nFrames;
    bool writeFailed; // latched on first write error -> further appends no-op
};

} // namespace ladruno

#endif
