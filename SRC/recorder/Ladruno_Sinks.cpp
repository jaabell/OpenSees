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
// Description: persistence strategies of the ladruno recorder.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#include "Ladruno_Sinks.h"
#include "Ladruno_Hdf5.h" // pulls Ladruno_Types.h; provides ladruno::h5::*

#include "OPS_Globals.h"   // opserr (failure reports)

#include <sstream>
#include <cmath>

namespace ladruno {

	/* ===================================================================== */
	/* ResultFamily                                                          */
	/* ===================================================================== */

	const char* ResultFamily::groupName(ResultFamily::Enum f)
	{
		switch (f) {
		case ResultFamily::OnNodes:    return "ON_NODES";
		case ResultFamily::OnElements: return "ON_ELEMENTS";
		case ResultFamily::OnDomain:   return "ON_DOMAIN";
		case ResultFamily::OnRegions:  return "ON_REGIONS";
		default:                       return "ON_NODES";
		}
	}

	/* ===================================================================== */
	/* tree navigation helpers (local to this TU)                            */
	/*                                                                       */
	/* The recorder (LadrunoRecorder::writeModel) created the stage as   */
	/* "MODEL_STAGE[<current_model_stage_id>]" off the file root and an      */
	/* (empty) RESULTS child. The sinks resolve that same path here, opening */
	/* or creating each link in the chain so a sink works whether the        */
	/* recorder pre-created RESULTS or not, and whether this is the first    */
	/* result of the stage or a sibling.                                     */
	/* ===================================================================== */

	namespace {

		// MODEL_STAGE group name, matching LadrunoRecorder::writeModel exactly.
		std::string stageGroupName(const detail::ProcessInfo& info)
		{
			std::stringstream ss;
			ss << "MODEL_STAGE[" << info.current_model_stage_id << "]";
			return ss.str();
		}

		// Open child `name` under `loc` if it exists, else create it (creation
		// goes through the wrapper; opening an existing group is plain navigation).
		hid_t openOrCreateGroup(hid_t loc, const char* name, hid_t gplist)
		{
			htri_t exists = H5Lexists(loc, name, H5P_DEFAULT);
			if (exists > 0) {
				return H5Gopen2(loc, name, H5P_DEFAULT);
			}
			return h5::group::create(loc, name, H5P_DEFAULT, gplist, H5P_DEFAULT);
		}

		// Resolve MODEL_STAGE[..]/RESULTS/<family> and return its open handle.
		// Caller closes it. Returns HID_INVALID on failure.
		hid_t openFamilyGroup(detail::ProcessInfo& info, ResultFamily::Enum family)
		{
			hid_t h_stage = openOrCreateGroup(info.h_file_id,
				stageGroupName(info).c_str(), info.h_group_proplist);
			if (h_stage == HID_INVALID || h_stage < 0)
				return HID_INVALID;

			hid_t h_results = openOrCreateGroup(h_stage, "RESULTS", info.h_group_proplist);
			h5::group::close(h_stage);
			if (h_results == HID_INVALID || h_results < 0)
				return HID_INVALID;

			hid_t h_family = openOrCreateGroup(h_results,
				ResultFamily::groupName(family), info.h_group_proplist);
			h5::group::close(h_results);
			return h_family;
		}

	} // anonymous namespace

	/* ===================================================================== */
	/* StreamingSink                                                         */
	/* ===================================================================== */

	void StreamingSink::fail(const ResultSchema& schema, const char* why)
	{
		opserr << "LadrunoRecorder error: result \"" << schema.name.c_str()
		       << "\" " << why << "; it will NOT be recorded from here on\n";
		m_dead = true;
		m_initialized = true;   // never re-run begin() (no per-step retry/spam)
	}

	void StreamingSink::begin(detail::ProcessInfo& info, const ResultSource& src)
	{
		const ResultSchema& schema = src.schema();
		if (schema.num_components < 1)
			return;
		// An empty channel (no ids on this process) has nothing to write. Its
		// [T x 0 x C] DATA could never be created (chunk dim 1 > max dim 0) and
		// used to leave a DATA-less result group behind; write nothing instead.
		if (src.ids().empty()) {
			m_dead = true;
			m_initialized = true;
			return;
		}

		hid_t h_family = openFamilyGroup(info, m_family);
		if (h_family == HID_INVALID || h_family < 0) {
			fail(schema, "could not open its RESULTS family group");
			return;
		}

		// result group names are unique per channel within a MODEL_STAGE
		// (element buckets embed classTag + rule + header index), so a
		// pre-existing group means this result was requested twice (e.g.
		// `-N displacement displacement`) or two recorders target one file.
		// Appending into the existing group would double its rows silently;
		// refuse instead.
		if (H5Lexists(h_family, schema.name.c_str(), H5P_DEFAULT) > 0) {
			h5::group::close(h_family);
			fail(schema, "already exists in this MODEL_STAGE (duplicate request, "
			             "or two recorders writing the same file)");
			return;
		}

		// Result group with the self-describing attrs.
		hid_t h_gp_result = h5::group::createResultGroup(
			h_family, info.h_group_proplist,
			schema.name, schema.display_name,
			schema.components_csv, schema.num_components,
			schema.dimension, schema.description,
			(int)schema.result_type, (int)schema.data_type);
		if (h_gp_result < 0) {
			h5::group::close(h_family);
			fail(schema, "could not create its result group");
			return;
		}
		// how a reader combines this result across partition files.
		h5::attribute::write(h_gp_result, "PARTITION_REDUCTION",
			std::string(src.partitionReduction()));

		// ID dataset [nIds x 1] (same shape as in MPCORecorder).
		const std::vector<int>& ids = src.ids();
		const size_t n_ids = ids.size();
		const size_t n_comp = (size_t)schema.num_components;
		hid_t h_dset_id = h5::dataset::createAndWrite(
			h_gp_result, "ID", ids, n_ids, 1);

		// Chunked time-series layout: DATA is a single [T x nIds x nComp]
		// extensible dataset that grows one slab per step, with matching
		// TIME[T] (double) / STEP[T] (int) axes. accept() appends.
		// On-disk DATA type is f64 (default, lossless) or f32 (opt-in
		// `-precision f32` lossy mode); TIME/STEP stay f64/int either way.
		hid_t data_disk_type = info.store_data_f32 ? H5T_IEEE_F32LE : H5T_IEEE_F64LE;
		// chunk plan of planTimeSeries() + sized chunk cache, deflate level from
		// `-compress`; the handles are KEPT OPEN until the sink is destroyed.
		m_data = h5::dataset::createTimeSeries3dOpen(
			h_gp_result, "DATA", (hsize_t)n_ids, (hsize_t)n_comp, data_disk_type,
			info.deflate_level);
		m_time = h5::dataset::createTimeAxis1d(h_gp_result, "TIME", H5T_IEEE_F64LE);
		m_step = h5::dataset::createTimeAxis1d(h_gp_result, "STEP", H5T_STD_I32LE);
		hsize_t mdims[3] = { 1, (hsize_t)n_ids, (hsize_t)n_comp };
		m_mspace = H5Screate_simple(3, mdims, NULL);
		hsize_t mdims1[1] = { 1 };
		m_mspace1 = H5Screate_simple(1, mdims1, NULL);
		m_n_ids = n_ids;
		m_n_comp = n_comp;
		m_t = 0;
		const bool ok = (h_dset_id >= 0 && m_data >= 0 && m_time >= 0 && m_step >= 0 &&
		                 m_mspace >= 0 && m_mspace1 >= 0);

		if (h_dset_id >= 0) h5::dataset::close(h_dset_id);
		h5::group::close(h_gp_result);
		h5::group::close(h_family);

		if (!ok) {
			closeHandles();
			// reported: otherwise every later accept() would fail and the
			// whole result would be dropped silently.
			fail(schema, "could not create its ID/DATA/TIME/STEP datasets "
			             "(disk full or quota, or HDF5 refused the dataset)");
			return;
		}
		m_initialized = true;
	}

	void StreamingSink::accept(detail::ProcessInfo& info, const ResultSource& src,
	                           const std::vector<double>& buffer)
	{
		const ResultSchema& schema = src.schema();
		if (schema.num_components < 1)
			return;

		// Defensive: if accept() is reached without a begin() for this stage
		// (e.g. a re-begin path missed it), create the group lazily so a step
		// is never silently dropped.
		if (!m_initialized)
			begin(info, src);
		if (m_dead)
			return;

		const std::vector<int>& ids = src.ids();
		const size_t n_ids = ids.size();
		const size_t n_comp = (size_t)schema.num_components;

		// a short buffer skips the WHOLE step (DATA and TIME/STEP together)
		// so the three axes stay aligned.
		if (buffer.size() < n_ids * n_comp) {
			if (!m_warned_short) {
				opserr << "LadrunoRecorder warning: result \"" << schema.name.c_str()
				       << "\" got " << (int)buffer.size() << " values, expected "
				       << (int)(n_ids * n_comp) << "; step(s) skipped\n";
				m_warned_short = true;
			}
			return;
		}

		// write into the handles held open since begin(). The slab lands in
		// the dataset's chunk cache; a chunk is deflated once, when it is
		// complete or at a flush, instead of being re-read and re-deflated on
		// every step.
		herr_t st = h5::dataset::writeSlab3dAt(m_data, m_t, &buffer[0],
			(hsize_t)n_ids, (hsize_t)n_comp, m_mspace);
		if (st < 0) {
			// if DATA was extended before the write failed, that row
			// exists (fill value); give it its TIME/STEP so the axes stay aligned,
			// then stop the channel.
			const bool grew = h5::dataset::extent0(m_data) > m_t;
			fail(schema, "failed to append a step to DATA (disk full or quota?)");
			if (!grew)
				return;
		}
		const double t_val = info.current_time_step;
		const int s_val = info.current_time_step_id;
		h5::dataset::writeScalar1dAt(m_time, m_t, H5T_NATIVE_DOUBLE, &t_val, m_mspace1);
		h5::dataset::writeScalar1dAt(m_step, m_t, H5T_NATIVE_INT, &s_val, m_mspace1);
		++m_t;
	}

	void StreamingSink::closeHandles()
	{
		if (m_mspace1 >= 0) H5Sclose(m_mspace1);
		if (m_mspace >= 0) H5Sclose(m_mspace);
		if (m_step >= 0) H5Dclose(m_step);
		if (m_time >= 0) H5Dclose(m_time);
		if (m_data >= 0) H5Dclose(m_data);
		m_data = m_time = m_step = m_mspace = m_mspace1 = HID_INVALID;
	}

	StreamingSink::~StreamingSink()
	{
		closeHandles();
	}

	void StreamingSink::finalize(detail::ProcessInfo& /*info*/)
	{
		// Every step is handed to HDF5 in accept(); partially filled chunks sit in
		// the open dataset's chunk cache until the recorder's H5Fflush (its flush
		// cadence) or the handles close in the destructor. Nothing to do here.
	}

	/* ===================================================================== */
	/* EnvelopeSink                                                          */
	/* ===================================================================== */

	void EnvelopeSink::begin(detail::ProcessInfo& /*info*/, const ResultSource& src)
	{
		// Cache identity from the source (stable for its lifetime). Accumulators
		// are not allocated here; they are seeded on the first accept() so MIN/MAX
		// begin at the real first sample rather than +/-inf sentinels.
		const ResultSchema& schema = src.schema();
		m_name = schema.name;
		m_display_name = schema.display_name;
		m_components_csv = schema.components_csv;
		m_dimension = schema.dimension;
		m_description = schema.description;
		m_result_type = (int)schema.result_type;
		m_data_type = (int)schema.data_type;
		m_partition_reduction = src.partitionReduction();
		m_n_comp = (size_t)(schema.num_components < 0 ? 0 : schema.num_components);
		m_ids = src.ids();
		m_n_ids = m_ids.size();
	}

	void EnvelopeSink::accept(detail::ProcessInfo& info, const ResultSource& src,
	                          const std::vector<double>& buffer)
	{
		const ResultSchema& schema = src.schema();
		if (schema.num_components < 1)
			return;

		// Ensure identity + self-describing metadata are cached even if begin() was
		// skipped (the recorder drives accept() directly, no begin() call).
		if (m_name.empty()) {
			m_name = schema.name;
			m_display_name = schema.display_name;
			m_components_csv = schema.components_csv;
			m_dimension = schema.dimension;
			m_description = schema.description;
			m_result_type = (int)schema.result_type;
			m_data_type = (int)schema.data_type;
			m_partition_reduction = src.partitionReduction();   // begin() is skipped here
		}
		m_n_comp = (size_t)schema.num_components;
		if (m_ids.empty()) {
			m_ids = src.ids();
			m_n_ids = m_ids.size();
		}

		const size_t n = m_n_ids * m_n_comp;
		if (n == 0 || buffer.size() < n)
			return;

		const int step = info.current_time_step_id;

		if (!m_seeded) {
			// First accept: seed every accumulator from this step's values.
			m_min.assign(buffer.begin(), buffer.begin() + n);
			m_max.assign(buffer.begin(), buffer.begin() + n);
			m_absmax.resize(n);
			m_arg_step.assign(n, step);
			for (size_t i = 0; i < n; ++i)
				m_absmax[i] = std::abs(buffer[i]);
			m_seeded = true;
			return;
		}

		// Subsequent accepts: element-wise componentwise update (each component
		// tracked independently; ARG_STEP records the step at which the absolute
		// extreme for THAT component was last attained).
		for (size_t i = 0; i < n; ++i) {
			const double v = buffer[i];
			// NaN is STICKY. Ordered comparisons with NaN are false, so a run
			// that went NaN at step k would keep its finite pre-divergence
			// extremes and the envelope would look healthy. The first NaN
			// poisons MIN/MAX/ABSMAX of that component and ARG_STEP records the
			// step it appeared (a NaN first sample already sticks).
			if (v != v) {
				if (m_absmax[i] == m_absmax[i]) {
					m_min[i] = v; m_max[i] = v; m_absmax[i] = v;
					m_arg_step[i] = step;
				}
				continue;
			}
			if (v < m_min[i]) m_min[i] = v;
			if (v > m_max[i]) m_max[i] = v;
			const double a = std::abs(v);
			if (a > m_absmax[i]) {
				m_absmax[i] = a;
				m_arg_step[i] = step;
			}
		}
	}

	void EnvelopeSink::writeEnvelope(detail::ProcessInfo& info)
	{
		if (!m_seeded || m_name.empty() || m_n_ids == 0 || m_n_comp == 0)
			return;

		// after the first write, overwrite the four accumulator datasets in
		// place: no group walk, no delete/recreate, no attributes.
		if (m_created) {
			if (m_dmin >= 0) H5Dwrite(m_dmin, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_min[0]);
			if (m_dmax >= 0) H5Dwrite(m_dmax, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_max[0]);
			if (m_dabs >= 0) H5Dwrite(m_dabs, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_absmax[0]);
			if (m_darg >= 0) H5Dwrite(m_darg, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_arg_step[0]);
			return;
		}

		// ENVELOPES/<family>/<name> lives under the stage's RESULTS tree.
		hid_t h_stage = openOrCreateGroup(info.h_file_id,
			stageGroupName(info).c_str(), info.h_group_proplist);
		if (h_stage == HID_INVALID || h_stage < 0)
			return;
		hid_t h_results = openOrCreateGroup(h_stage, "RESULTS", info.h_group_proplist);
		h5::group::close(h_stage);
		if (h_results == HID_INVALID || h_results < 0)
			return;
		hid_t h_env = openOrCreateGroup(h_results, "ENVELOPES", info.h_group_proplist);
		h5::group::close(h_results);
		if (h_env == HID_INVALID || h_env < 0)
			return;
		hid_t h_family = openOrCreateGroup(h_env,
			ResultFamily::groupName(m_family), info.h_group_proplist);
		h5::group::close(h_env);
		if (h_family == HID_INVALID || h_family < 0)
			return;

		// Element result names are nested ("<display>/<bucket>", e.g.
		// "stress/204-FourNodeQuad[201:0:0]"); node/domain names are flat. The group
		// proplist carries no create-intermediate-group link property, so H5Gcreate
		// (in createResultGroup below) will NOT auto-create the "<display>" parent:
		// pre-create each ancestor prefix once.
		for (size_t slash = m_name.find('/'); slash != std::string::npos;
		     slash = m_name.find('/', slash + 1)) {
			std::string prefix = m_name.substr(0, slash);
			if (H5Lexists(h_family, prefix.c_str(), H5P_DEFAULT) <= 0) {
				hid_t h_mid = h5::group::create(h_family, prefix.c_str(),
					H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
				if (h_mid >= 0)
					h5::group::close(h_mid);
			}
		}

		// first write of this stage: create the group once. A pre-existing
		// group can only be a duplicate request, which the parser drops.
		if (H5Lexists(h_family, m_name.c_str(), H5P_DEFAULT) > 0) {
			opserr << "LadrunoRecorder error: envelope \"" << m_name.c_str()
			       << "\" already exists in this MODEL_STAGE; it will NOT be recorded\n";
			m_created = true;   // never retry (handles stay invalid)
			h5::group::close(h_family);
			return;
		}

		// Self-describing result group (same COMPONENTS/DISPLAY_NAME/DIMENSION attrs
		// as the time-series StreamingSink writes) so envelope output carries
		// authoritative component names in-file, not just bare arrays.
		hid_t h_name = h5::group::createResultGroup(
			h_family, info.h_group_proplist, m_name, m_display_name,
			m_components_csv, (int)m_n_comp, m_dimension, m_description,
			m_result_type, m_data_type);
		// same attribute as the streaming group. In a partitioned run the
		// recorder refuses to envelope anything but "NONE", so a written envelope is
		// always "NONE" there; serial files keep the source's own value.
		h5::attribute::write(h_name, "PARTITION_REDUCTION", m_partition_reduction);

		// ID [nIds x 1], and the four [nIds x nComp] accumulators; their
		// handles are kept for in-place rewrites.
		hid_t d_id = h5::dataset::createAndWrite(h_name, "ID", m_ids, m_n_ids, 1);
		m_dmin = h5::dataset::createAndWrite(h_name, "MIN", m_min, m_n_ids, m_n_comp);
		m_dmax = h5::dataset::createAndWrite(h_name, "MAX", m_max, m_n_ids, m_n_comp);
		m_dabs = h5::dataset::createAndWrite(h_name, "ABSMAX", m_absmax, m_n_ids, m_n_comp);
		m_darg = h5::dataset::createAndWrite(h_name, "ARG_STEP", m_arg_step, m_n_ids, m_n_comp);
		m_created = true;

		h5::dataset::close(d_id);
		h5::group::close(h_name);
		h5::group::close(h_family);
	}

	EnvelopeSink::~EnvelopeSink()
	{
		if (m_darg >= 0) H5Dclose(m_darg);
		if (m_dabs >= 0) H5Dclose(m_dabs);
		if (m_dmax >= 0) H5Dclose(m_dmax);
		if (m_dmin >= 0) H5Dclose(m_dmin);
	}

	void EnvelopeSink::flush(detail::ProcessInfo& info)
	{
		writeEnvelope(info);
	}

	void EnvelopeSink::finalize(detail::ProcessInfo& info)
	{
		writeEnvelope(info);
	}

	void EnvelopeSink::reset()
	{
		m_seeded = false;
		m_n_ids = 0;
		m_n_comp = 0;
		m_name.clear();
		m_ids.clear();
		m_min.clear();
		m_max.clear();
		m_absmax.clear();
		m_arg_step.clear();
		if (m_darg >= 0) H5Dclose(m_darg);
		if (m_dabs >= 0) H5Dclose(m_dabs);
		if (m_dmax >= 0) H5Dclose(m_dmax);
		if (m_dmin >= 0) H5Dclose(m_dmin);
		m_dmin = m_dmax = m_dabs = m_darg = HID_INVALID;
		m_created = false;
	}

} // namespace ladruno
