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
// Description: thin inline wrappers over the HDF5 C API used by the ladruno
// recorder (attributes, fixed datasets, chunked extensible time-series
// datasets, groups). Header-only, namespace ladruno::h5.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#ifndef Ladruno_Hdf5_h
#define Ladruno_Hdf5_h


#include "Ladruno_Types.h"

namespace ladruno {

/*ladruno - hdf5 interface utilities*/
namespace h5 {

	namespace attribute {

		// low level functions for c interface

		inline herr_t writei(hid_t obj, const char *attr_name, const int *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t attr = H5Acreate(obj, attr_name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, H5T_NATIVE_INT, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			return status;
		}
		inline herr_t writed(hid_t obj, const char *attr_name, const double *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t attr = H5Acreate(obj, attr_name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, H5T_NATIVE_DOUBLE, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			return status;
		}
		inline herr_t writes(hid_t obj, const char *attr_name, const char *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { 1 };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t atype = H5Tcopy(H5T_C_S1);
			status = H5Tset_size(atype, data_size);
			status = H5Tset_strpad(atype, H5T_STR_NULLTERM);
			hid_t attr = H5Acreate(obj, attr_name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, atype, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			H5Tclose(atype);
			return status;
		}

		// higher level c++ utils

		inline herr_t write(hid_t obj, const char *attr_name, int attr_data)
		{
			return writei(obj, attr_name, &attr_data, 1);
		}
		inline herr_t write(hid_t obj, const char *attr_name, const std::vector<int> &attr_data)
		{
			if (attr_data.size() > 0) {
				return writei(obj, attr_name, &attr_data[0], attr_data.size());
			}
			return 0;
		}
		inline herr_t write(hid_t obj, const char *attr_name, double attr_data)
		{
			return writed(obj, attr_name, &attr_data, 1);
		}
		inline herr_t write(hid_t obj, const char *attr_name, std::vector<double> &attr_data)
		{
			if (attr_data.size() > 0) {
				return writed(obj, attr_name, &attr_data[0], attr_data.size());
			}
			return 0;
		}
		inline herr_t write(hid_t obj, const char *attr_name, const std::string &attr_data)
		{
			if (!attr_data.empty()) {
				return writes(obj, attr_name, attr_data.c_str(), attr_data.size());
			}
			return 0;
		}

	}

	namespace group {

		// low level functions for c interface

		inline hid_t create(hid_t loc_id, const char *name, hid_t lcpl_id, hid_t gcpl_id, hid_t gapl_id) {
			return H5Gcreate(loc_id, name, lcpl_id, gcpl_id, gapl_id);
		}
		inline herr_t close(hid_t group_id) {
			if (group_id == HID_INVALID) return -1;
			return H5Gclose(group_id);
		}

		// higher level c++ utils

		inline hid_t createResultGroup(hid_t obj, hid_t gplist, const std::string &result_name,
			const std::string &disp_name, const std::string &components, int num_components,
			const std::string &dimension, const std::string &description,
			int result_type, int data_type)
		{
			herr_t status;
			hid_t h_gp_result = h5::group::create(obj, result_name.c_str(), H5P_DEFAULT, gplist, H5P_DEFAULT);
			status = h5::attribute::write(h_gp_result, "DISPLAY_NAME", disp_name);
			status = h5::attribute::write(h_gp_result, "COMPONENTS", components);
			status = h5::attribute::write(h_gp_result, "DIMENSION", dimension);
			status = h5::attribute::write(h_gp_result, "DESCRIPTION", description);
			status = h5::attribute::write(h_gp_result, "TYPE", result_type);
			status = h5::attribute::write(h_gp_result, "DATA_TYPE", data_type);
			return h_gp_result;
		}

	}

	namespace dataset {

		// low level functions for c interface

		inline herr_t close(hid_t dataset_id) {
			if (dataset_id == HID_INVALID) return -1;
			return H5Dclose(dataset_id);
		}
		inline hid_t createAndWrited1(hid_t obj, const char *name, const double *data, hsize_t data_size) {
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		inline hid_t createAndWrited2(hid_t obj, const char *name, const double *data, hsize_t rows, hsize_t cols)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[2] = { rows, cols };
			hid_t space = H5Screate_simple(2, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		inline hid_t createAndWritei1(hid_t obj, const char *name, const int *data, hsize_t data_size)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		inline hid_t createAndWritei2(hid_t obj, const char *name, const int *data, hsize_t rows, hsize_t cols)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[2] = { rows, cols };
			hid_t space = H5Screate_simple(2, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		inline hid_t createAndWrites(hid_t obj, const char *name, const char *data, hsize_t data_size)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { 1 };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t atype = H5Tcopy(H5T_C_S1);
			status = H5Tset_size(atype, data_size);
			status = H5Tset_strpad(atype, H5T_STR_NULLTERM);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			H5Tclose(atype);
			return dset;
		}

		// higher level c++ utils

		inline hid_t createAndWrite(hid_t obj, const char *name, double data)
		{
			return createAndWrited1(obj, name, &data, 1);
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<double> &data)
		{
			if (data.size() > 0) {
				return createAndWrited1(obj, name, &data[0], data.size());
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<double> &data, size_t rows, size_t cols)
		{
			if (data.size() > 0 && data.size() == rows*cols) {
				return createAndWrited2(obj, name, &data[0], rows, cols);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, int data)
		{
			return createAndWritei1(obj, name, &data, 1);
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<int> &data)
		{
			if (data.size() > 0) {
				return createAndWritei1(obj, name, &data[0], data.size());
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<int> &data, size_t rows, size_t cols)
		{
			if (data.size() > 0 && data.size() == rows*cols) {
				return createAndWritei2(obj, name, &data[0], rows, cols);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::string &data)
		{
			if (!data.empty()) {
				return createAndWrites(obj, name, data.c_str(), data.size());
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<utils::locax::quaternion> &data)
		{
			if (data.size() > 0) {
				return createAndWrited2(obj, name, data[0].data(), data.size(), 4);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<detail::element::FiberData> &data)
		{
			if (data.size() > 0) {
				return createAndWrited2(obj, name, data[0].data(), data.size(), 3);
			}
			return HID_INVALID;
		}
		inline hid_t createAndWrite(hid_t obj, const char *name, const std::vector<detail::element::ElemGaussPair> &data)
		{
			if (data.size() > 0) {
				return createAndWritei2(obj, name, data[0].data(), data.size(), 2);
			}
			return HID_INVALID;
		}

		// ---- extensible (chunked) time-series datasets ------------------------
		// One [T x nIds x nComp] dataset per result grows along the unlimited T
		// (time/mode) axis; one slab is appended per committed step. Chunked +
		// shuffle + deflate so the (compressible) time axis is compressed and a
		// single time-history is an O(1) dataset open instead of O(T) per-step
		// dataset opens. TIME[T] (double) / STEP[T] (int) are the matching
		// extensible 1-D axes.

		// choose how many T-slabs per chunk to land near ~256 KiB/chunk.
		inline hsize_t chunkSlabsPerBlock(hsize_t slab_elems, hsize_t elem_bytes) {
			const hsize_t target = 262144; // 256 KiB
			hsize_t per = slab_elems * elem_bytes;
			if (per == 0) return 1;
			hsize_t n = target / per;
			if (n < 1)    n = 1;
			if (n > 1024) n = 1024;
			return n;
		}

		// create a [0 x nIds x nComp] dataset, maxdims [UNLIM, nIds, nComp].
		// disk_type selects the ON-DISK element type (default f64 = lossless parity
		// path). Pass H5T_IEEE_F32LE for the opt-in `-precision f32` lossy mode:
		// appendSlab3d still writes from an in-memory double array and HDF5 narrows
		// f64->f32 on write, halving the payload at ~7 significant digits. The chunk
		// byte budget uses the on-disk element size so chunks stay ~256 KiB either way.
		inline hid_t createTimeSeries3d(hid_t obj, const char *name, hsize_t n_ids, hsize_t n_comp,
		                                hid_t disk_type = H5T_IEEE_F64LE) {
			if (n_ids < 1 || n_comp < 1) return HID_INVALID;  // unrepresentable
			hsize_t dims[3]    = { 0, n_ids, n_comp };
			hsize_t maxdims[3] = { H5S_UNLIMITED, n_ids, n_comp };
			hid_t space = H5Screate_simple(3, dims, maxdims);
			if (space < 0) return HID_INVALID;
			hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
			if (dcpl < 0) { H5Sclose(space); return HID_INVALID; }
			hsize_t elem_bytes = (H5Tequal(disk_type, H5T_IEEE_F32LE) > 0) ? 4 : 8;
			hsize_t ct = chunkSlabsPerBlock(n_ids * n_comp, elem_bytes);
			// HDF5 refuses a chunk of >= 4 GiB and the chunk spans every id, so a
			// slab above 2^29 f64 values cannot be created with one id tile. Above
			// a 1 GiB slab the id axis is tiled; below it one chunk holds the whole
			// slab for ct steps.
			hsize_t chunk_ids = n_ids;
			const hsize_t max_chunk_bytes = (hsize_t)1 << 30;  // 1 GiB
			const hsize_t row_bytes = n_comp * elem_bytes;
			if (n_ids * row_bytes > max_chunk_bytes) {
				chunk_ids = max_chunk_bytes / row_bytes;
				if (chunk_ids < 1) chunk_ids = 1;
			}
			hsize_t chunk[3] = { ct, chunk_ids, n_comp };
			herr_t st = H5Pset_chunk(dcpl, 3, chunk);
			if (st >= 0) st = H5Pset_shuffle(dcpl);
			if (st >= 0) st = H5Pset_deflate(dcpl, 4);
			hid_t dset = (st >= 0)
				? H5Dcreate(obj, name, disk_type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT)
				: HID_INVALID;
			H5Pclose(dcpl);
			H5Sclose(space);
			return (dset < 0) ? HID_INVALID : dset;
		}

		// ---- write-path layout for datasets that stay open for a stage -------
		//
		// A chunk of {ct, nIds, nComp} spans every id, so reading one entity's
		// history inflates the whole dataset; reopening DATA every step also
		// evicts the per-dataset chunk cache, so a partially filled chunk is
		// re-read, inflated and re-deflated on every step.
		//
		// Chunks target ~1 MiB. A small slab (<= 1 MiB) keeps every id in one
		// chunk and stacks up to 1024 steps in it; a large slab tiles the id axis
		// and stacks ct (<= 16) steps, so one entity's history touches T/ct chunks
		// of ~1 MiB. The dataset access list sizes the chunk cache to hold a whole
		// row of chunks (ct steps x every id tile), so with the handle kept open a
		// chunk is filled in memory and deflated ONCE, when it is complete.
		struct TimeSeriesPlan {
			hsize_t ct;          // steps per chunk
			hsize_t chunk_ids;   // ids per chunk
			size_t cache_bytes;  // chunk-cache size for the dataset access list
			size_t cache_slots;  // hash slots for the chunk cache
		};

		inline TimeSeriesPlan planTimeSeries(hsize_t n_ids, hsize_t n_comp, hsize_t elem_bytes) {
			const hsize_t target = (hsize_t)1 << 20;          // ~1 MiB per chunk
			const hsize_t row_budget = (hsize_t)32 << 20;     // ct x slab held in cache
			const hsize_t row_bytes = n_comp * elem_bytes;
			const hsize_t slab = n_ids * row_bytes;
			TimeSeriesPlan p;
			if (slab <= target) {
				p.chunk_ids = n_ids;
				p.ct = slab ? target / slab : 1;
				if (p.ct < 1) p.ct = 1;
				if (p.ct > 1024) p.ct = 1024;
			}
			else {
				p.ct = row_budget / slab;
				if (p.ct < 1) p.ct = 1;
				if (p.ct > 16) p.ct = 16;
				p.chunk_ids = target / (p.ct * row_bytes);
				if (p.chunk_ids < 1) p.chunk_ids = 1;
				if (p.chunk_ids > n_ids) p.chunk_ids = n_ids;
			}
			const hsize_t tiles = (n_ids + p.chunk_ids - 1) / p.chunk_ids;
			const hsize_t chunk_bytes = p.ct * p.chunk_ids * row_bytes;
			// A row of chunks (every id tile x ct steps) plus one spare chunk. With
			// ct == 1 every chunk is completed by a single write, so a cache that
			// cannot hold the row costs nothing; cap it at 256 MiB.
			hsize_t cache = tiles * chunk_bytes + chunk_bytes;
			const hsize_t cap = (hsize_t)256 << 20;
			if (cache > cap) cache = cap;
			if (cache < ((hsize_t)1 << 20)) cache = (hsize_t)1 << 20;
			p.cache_bytes = (size_t)cache;
			p.cache_slots = (size_t)(tiles * 10 + 521) | 1;
			return p;
		}

		// Create a [0 x nIds x nComp] extensible dataset with planTimeSeries() and a
		// sized chunk cache; deflate 0 disables the shuffle+deflate filters.
		inline hid_t createTimeSeries3dOpen(hid_t obj, const char *name, hsize_t n_ids,
		                                    hsize_t n_comp, hid_t disk_type, int deflate) {
			if (n_ids < 1 || n_comp < 1) return HID_INVALID;
			hsize_t elem_bytes = (H5Tequal(disk_type, H5T_IEEE_F32LE) > 0) ? 4 : 8;
			TimeSeriesPlan p = planTimeSeries(n_ids, n_comp, elem_bytes);
			hsize_t dims[3]    = { 0, n_ids, n_comp };
			hsize_t maxdims[3] = { H5S_UNLIMITED, n_ids, n_comp };
			hsize_t chunk[3]   = { p.ct, p.chunk_ids, n_comp };
			hid_t space = H5Screate_simple(3, dims, maxdims);
			if (space < 0) return HID_INVALID;
			hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
			hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
			herr_t st = (dcpl >= 0 && dapl >= 0) ? H5Pset_chunk(dcpl, 3, chunk) : -1;
			if (st >= 0 && deflate > 0) {
				st = H5Pset_shuffle(dcpl);
				if (st >= 0) st = H5Pset_deflate(dcpl, (unsigned)(deflate > 9 ? 9 : deflate));
			}
			// w0 = 1.0: evict fully written chunks first.
			if (st >= 0) st = H5Pset_chunk_cache(dapl, p.cache_slots, p.cache_bytes, 1.0);
			hid_t dset = (st >= 0)
				? H5Dcreate(obj, name, disk_type, space, H5P_DEFAULT, dcpl, dapl)
				: HID_INVALID;
			if (dapl >= 0) H5Pclose(dapl);
			if (dcpl >= 0) H5Pclose(dcpl);
			H5Sclose(space);
			return (dset < 0) ? HID_INVALID : dset;
		}

		// Write slab t of an open [T x nIds x nComp] dataset (extend to t+1 first).
		// `mspace` is a caller-owned [1 x nIds x nComp] memory space reused per step.
		inline herr_t writeSlab3dAt(hid_t dset, hsize_t t, const double *data,
		                            hsize_t n_ids, hsize_t n_comp, hid_t mspace) {
			hsize_t newdims[3] = { t + 1, n_ids, n_comp };
			if (H5Dset_extent(dset, newdims) < 0) return -1;
			hid_t fspace = H5Dget_space(dset);
			if (fspace < 0) return -1;
			hsize_t start[3] = { t, 0, 0 };
			hsize_t count[3] = { 1, n_ids, n_comp };
			herr_t status = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			if (status >= 0)
				status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, data);
			H5Sclose(fspace);
			return status;
		}

		// Write element t of an open [T] axis (extend to t+1 first).
		inline herr_t writeScalar1dAt(hid_t dset, hsize_t t, hid_t mem_type, const void *value,
		                              hid_t mspace1) {
			hsize_t newdims[1] = { t + 1 };
			if (H5Dset_extent(dset, newdims) < 0) return -1;
			hid_t fspace = H5Dget_space(dset);
			if (fspace < 0) return -1;
			hsize_t start[1] = { t }, count[1] = { 1 };
			herr_t status = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			if (status >= 0)
				status = H5Dwrite(dset, mem_type, mspace1, fspace, H5P_DEFAULT, value);
			H5Sclose(fspace);
			return status;
		}

		// current length of the unlimited (first) axis, 0 on error.
		inline hsize_t extent0(hid_t dset) {
			hid_t fspace = H5Dget_space(dset);
			if (fspace < 0) return 0;
			hsize_t cur[3] = { 0, 0, 0 };
			int nd = H5Sget_simple_extent_ndims(fspace);
			if (nd >= 1 && nd <= 3)
				H5Sget_simple_extent_dims(fspace, cur, NULL);
			H5Sclose(fspace);
			return cur[0];
		}

		// append one [1 x nIds x nComp] slab to a [T x nIds x nComp] dataset.
		inline herr_t appendSlab3d(hid_t dset, const double *data, hsize_t n_ids, hsize_t n_comp) {
			hid_t fspace = H5Dget_space(dset);
			hsize_t cur[3];
			H5Sget_simple_extent_dims(fspace, cur, NULL);
			H5Sclose(fspace);
			hsize_t t = cur[0];
			hsize_t newdims[3] = { t + 1, n_ids, n_comp };
			if (H5Dset_extent(dset, newdims) < 0) return -1;
			fspace = H5Dget_space(dset);
			hsize_t start[3] = { t, 0, 0 };
			hsize_t count[3] = { 1, n_ids, n_comp };
			H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			hsize_t mdims[3] = { 1, n_ids, n_comp };
			hid_t mspace = H5Screate_simple(3, mdims, NULL);
			herr_t status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, data);
			H5Sclose(mspace);
			H5Sclose(fspace);
			return status;
		}

		// create a [0] extensible 1-D axis of the given HDF5 file type.
		inline hid_t createTimeAxis1d(hid_t obj, const char *name, hid_t h5type) {
			hsize_t dims[1] = { 0 };
			hsize_t maxdims[1] = { H5S_UNLIMITED };
			hid_t space = H5Screate_simple(1, dims, maxdims);
			hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
			hsize_t chunk[1] = { 1024 };
			H5Pset_chunk(dcpl, 1, chunk);
			hid_t dset = H5Dcreate(obj, name, h5type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
			H5Pclose(dcpl);
			H5Sclose(space);
			return dset;
		}

		// append one scalar to a [T] extensible axis (mem type from the caller).
		inline herr_t appendScalar1d(hid_t dset, hid_t mem_type, const void *value) {
			hid_t fspace = H5Dget_space(dset);
			hsize_t cur[1];
			H5Sget_simple_extent_dims(fspace, cur, NULL);
			H5Sclose(fspace);
			hsize_t t = cur[0];
			hsize_t newdims[1] = { t + 1 };
			if (H5Dset_extent(dset, newdims) < 0) return -1;
			fspace = H5Dget_space(dset);
			hsize_t start[1] = { t }, count[1] = { 1 };
			H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
			hsize_t mdims[1] = { 1 };
			hid_t mspace = H5Screate_simple(1, mdims, NULL);
			herr_t status = H5Dwrite(dset, mem_type, mspace, fspace, H5P_DEFAULT, value);
			H5Sclose(mspace);
			H5Sclose(fspace);
			return status;
		}
		inline herr_t appendDouble1d(hid_t dset, double value) {
			return appendScalar1d(dset, H5T_NATIVE_DOUBLE, &value);
		}
		inline herr_t appendInt1d(hid_t dset, int value) {
			return appendScalar1d(dset, H5T_NATIVE_INT, &value);
		}

	}

	namespace file {

		// low level functions for c interface

		inline herr_t close(hid_t file_id) {
			return H5Fclose(file_id);
		}
		inline herr_t flush(hid_t file_id) {
			return H5Fflush(file_id, H5F_SCOPE_LOCAL);
		}
		inline hid_t create(const char *filename, hid_t create_plist, hid_t acc_plist) {
			return H5Fcreate(filename, H5F_ACC_TRUNC, create_plist, acc_plist);
		}
	}

	namespace plist {

		enum CreateOptions {
			FileCreate,
			FileAccess,
			GroupCreate
		};

		// low level functions for c interface

		inline herr_t close(hid_t plist_id) {
			if (plist_id == HID_INVALID) return -1;
			return H5Pclose(plist_id);
		}
		inline hid_t crate(int opt) {
			switch (opt) {
			case FileCreate:
				return H5Pcreate(H5P_FILE_CREATE);
			case FileAccess:
				return H5Pcreate(H5P_FILE_ACCESS);
			case GroupCreate:
				return H5Pcreate(H5P_GROUP_CREATE);
			default:
				return HID_INVALID;
			}
		}
		inline herr_t setLinkCreationOrder(hid_t plist_id, unsigned int crt_order_flags) {
			return H5Pset_link_creation_order(plist_id, crt_order_flags);
		}
		inline herr_t setLibVerBounds(hid_t plist_id, unsigned int minor, unsigned int major) {
			return H5Pset_libver_bounds(plist_id, (H5F_libver_t)minor, (H5F_libver_t)major);
		}
	}

}

} // namespace ladruno

#endif // Ladruno_Hdf5_h
