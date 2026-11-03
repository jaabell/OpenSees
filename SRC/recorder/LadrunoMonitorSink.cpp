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
// Description: SWMR HDF5 sink of recorder Monitor (see LadrunoMonitorSink.h).

#include "LadrunoMonitorSink.h"

#include <OPS_Globals.h>   // opserr

// SWMR writing (H5Fstart_swmr_write) exists from HDF5 1.10.0 on.
#if !defined(H5_VERSION_GE)
#error "recorder Monitor needs HDF5 >= 1.10 (SWMR)"
#elif !H5_VERSION_GE(1, 10, 0)
#error "recorder Monitor needs HDF5 >= 1.10 (SWMR)"
#endif

namespace ladruno {

namespace {

// FRAMES/STEP/TIME grow one row per frame; chunk along the time axis so a
// long explicit run does not thrash. 256 rows/chunk is a reasonable balance
// between flush latency and metadata overhead.
const hsize_t kTimeChunk = 256;

herr_t writeStrAttr(hid_t obj, const char *name, const std::string &val)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, val.size() + 1);
    H5Tset_strpad(atype, H5T_STR_NULLTERM);
    hid_t attr = H5Acreate(obj, name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
    herr_t st = H5Awrite(attr, atype, val.c_str());
    H5Aclose(attr);
    H5Tclose(atype);
    H5Sclose(space);
    return st;
}

herr_t writeIntAttr(hid_t obj, const char *name, int val)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(obj, name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT);
    herr_t st = H5Awrite(attr, H5T_NATIVE_INT, &val);
    H5Aclose(attr);
    H5Sclose(space);
    return st;
}

// Create a 1-D chunked, unlimited dataset. Returns the (open) dataset id, or
// H5I_INVALID_HID on failure.
hid_t makeTimeSeries1D(hid_t file, const char *name, hid_t fileType)
{
    hsize_t dims[1] = {0};
    hsize_t maxd[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {kTimeChunk};
    hid_t space = H5Screate_simple(1, dims, maxd);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 1, chunk);
    hid_t dset = H5Dcreate(file, name, fileType, space,
                           H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Pclose(dcpl);
    H5Sclose(space);
    return dset;
}

} // namespace

MonitorSink::MonitorSink(const std::string &fname,
                         const std::vector<std::string> &cols,
                         const std::string &u)
    : filename(fname), columns(cols), units(u),
      fileId(H5I_INVALID_HID), dsStep(H5I_INVALID_HID),
      dsTime(H5I_INVALID_HID), dsFrames(H5I_INVALID_HID),
      nFrames(0), writeFailed(false)
{
}

MonitorSink::~MonitorSink()
{
    this->close();
}

int MonitorSink::open()
{
    const hsize_t nCols = static_cast<hsize_t>(columns.size());

    // SWMR requires the latest format bounds.
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    fileId = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    if (fileId < 0) {
        opserr << "LadrunoMonitorSink: could not create file " << filename.c_str() << "\n";
        fileId = H5I_INVALID_HID;
        return -1;
    }

    // Root attributes (self-describing header).
    writeStrAttr(fileId, "FORMAT", "ladruno-monitor");
    writeIntAttr(fileId, "FORMAT_VERSION", 1);
    writeStrAttr(fileId, "GENERATOR", "Ladruno");
    if (!units.empty())
        writeStrAttr(fileId, "units", units);

    // COLUMNS: variable-length strings, one self-describing label per channel.
    {
        hsize_t cdim[1] = {nCols > 0 ? nCols : 1};
        hid_t cspace = H5Screate_simple(1, cdim, NULL);
        hid_t cstr = H5Tcopy(H5T_C_S1);
        H5Tset_size(cstr, H5T_VARIABLE);
        hid_t cset = H5Dcreate(fileId, "COLUMNS", cstr, cspace,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (cset >= 0 && nCols > 0) {
            std::vector<const char *> ptrs(nCols);
            for (hsize_t i = 0; i < nCols; ++i)
                ptrs[i] = columns[i].c_str();
            H5Dwrite(cset, cstr, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ptrs[0]);
            H5Dclose(cset);
        }
        H5Tclose(cstr);
        H5Sclose(cspace);
    }

    // Extensible per-frame datasets.
    dsStep = makeTimeSeries1D(fileId, "STEP", H5T_STD_I32LE);
    dsTime = makeTimeSeries1D(fileId, "TIME", H5T_IEEE_F64LE);
    {
        hsize_t dims[2] = {0, nCols};
        hsize_t maxd[2] = {H5S_UNLIMITED, nCols};
        hsize_t chunk[2] = {kTimeChunk, nCols > 0 ? nCols : 1};
        hid_t space = H5Screate_simple(2, dims, maxd);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(dcpl, 2, chunk);
        dsFrames = H5Dcreate(fileId, "FRAMES", H5T_IEEE_F64LE, space,
                             H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Pclose(dcpl);
        H5Sclose(space);
    }

    if (dsStep < 0 || dsTime < 0 || dsFrames < 0) {
        opserr << "LadrunoMonitorSink: failed to create datasets in " << filename.c_str() << "\n";
        this->close();
        return -1;
    }

    // Arm SWMR -- NO new objects may be created after this point.
    if (H5Fstart_swmr_write(fileId) < 0) {
        opserr << "LadrunoMonitorSink: H5Fstart_swmr_write failed for " << filename.c_str() << "\n";
        this->close();
        return -1;
    }

    return 0;
}

int MonitorSink::extendAndWrite(hid_t dset, hid_t memType, hsize_t ncols,
                                const void *rowData)
{
    const bool twoD = (ncols > 1) || (dset == dsFrames);
    const int rank = twoD ? 2 : 1;
    const hsize_t newLen = static_cast<hsize_t>(nFrames) + 1;

    hsize_t dims[2];
    hsize_t start[2];
    hsize_t count[2];
    if (twoD) {
        dims[0] = newLen;            dims[1] = ncols;
        start[0] = nFrames;          start[1] = 0;
        count[0] = 1;                count[1] = ncols;
    } else {
        dims[0] = newLen;
        start[0] = nFrames;
        count[0] = 1;
    }

    if (H5Dset_extent(dset, dims) < 0)
        return -1;

    hid_t fspace = H5Dget_space(dset);
    if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        H5Sclose(fspace);
        return -1;
    }
    hsize_t mdims[2] = {1, ncols};
    hid_t mspace = H5Screate_simple(rank, mdims, NULL);
    herr_t st = H5Dwrite(dset, memType, mspace, fspace, H5P_DEFAULT, rowData);
    H5Sclose(mspace);
    H5Sclose(fspace);
    return st < 0 ? -1 : 0;
}

int MonitorSink::append(int step, double time, const std::vector<double> &values)
{
    if (!isOpen() || writeFailed)
        return 0; // sink unavailable -> silently no-op, never stalls the solver

    if (values.size() != columns.size()) {
        opserr << "LadrunoMonitorSink: append got " << (int)values.size()
               << " values, expected " << (int)columns.size() << "\n";
        return -1;
    }

    int rc = 0;
    rc |= extendAndWrite(dsStep, H5T_NATIVE_INT, 1, &step);
    rc |= extendAndWrite(dsTime, H5T_NATIVE_DOUBLE, 1, &time);
    rc |= extendAndWrite(dsFrames, H5T_NATIVE_DOUBLE,
                         static_cast<hsize_t>(columns.size()),
                         values.empty() ? NULL : &values[0]);

    if (rc != 0) {
        opserr << "LadrunoMonitorSink: frame write failed; monitor disabled\n";
        writeFailed = true;
        return -1;
    }

    // SWMR: flush so a tailing reader sees this frame.
    H5Dflush(dsStep);
    H5Dflush(dsTime);
    H5Dflush(dsFrames);

    ++nFrames;
    return 0;
}

void MonitorSink::close()
{
    if (dsFrames >= 0) { H5Dclose(dsFrames); dsFrames = H5I_INVALID_HID; }
    if (dsTime >= 0)   { H5Dclose(dsTime);   dsTime = H5I_INVALID_HID; }
    if (dsStep >= 0)   { H5Dclose(dsStep);   dsStep = H5I_INVALID_HID; }
    if (fileId >= 0)   { H5Fclose(fileId);   fileId = H5I_INVALID_HID; }
}

} // namespace ladruno
