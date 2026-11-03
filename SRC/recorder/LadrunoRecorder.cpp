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
// Description: recorder ladruno (see LadrunoRecorder.h).
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#include "LadrunoRecorder.h"

// ladruno modules (all symbols namespaced to avoid an ODR clash with MPCORecorder)
#include "Ladruno_Hdf5.h"     // pulls Ladruno_Types.h
#include "Ladruno_ResultIO.h"
#include "Ladruno_NodeResults.h"
#include "Ladruno_ElementResults.h"
#include "Ladruno_DomainResults.h"
#include "Ladruno_Sinks.h"
#include "Ladruno_LaunchEnv.h"   // launcher rank probe

// OpenSees
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>

// NOTE on parallel output: this TU is part of the shared (sequentially compiled)
// OPS_Recorder library; it links into OpenSeesPy/OpenSeesMP/openseesmp alike and
// is NOT compiled with _PARALLEL_PROCESSING, so it cannot call MPI directly (as
// MPCORecorder, whose only _PARALLEL_PROCESSING bit, an Allreduce, is likewise
// compiled out here). The per-partition output contract is therefore
// driven WITHOUT MPI, by two always-available signals (see initialize()):
//   * PartitionedDomain / OpenSeesMP: the recorder is broadcast P0->workers via
//     sendSelf/recvSelf, so send_self_count != 0 and p_id is the assigned index;
//   * interpreter-per-rank / openseesmp: each rank builds its own recorder, so the
//     partition index is read from the launcher's PMI_RANK/PMI_SIZE env (Intel
//     MPI / MS-MPI set these per rank). Single process -> filename used verbatim.
#include <elementAPI.h>
#include <classTags.h>
#include <Response.h>
#include <CompositeResponse.h>
#include <MeshRegion.h>
#include <Pressure_ConstraintIter.h>
#include <Pressure_Constraint.h>
#include "section/SectionForceDeformation.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

namespace ladrunons = ladruno; // local alias

// A run identity shared by every rank of one launch, when one
// exists. Order: LADRUNO_RUN_ID (user/launcher-provided), SLURM_JOB_ID[.STEP],
// OpenMPI's job id; else a process-unique id (scope "process": readers must not
// cross-check it between part files).
static void ladrunoResolveRunId(std::string& id, std::string& scope)
{
	const char* v = std::getenv("LADRUNO_RUN_ID");
	if (v && *v) { id = v; scope = "user"; return; }
	v = std::getenv("SLURM_JOB_ID");
	if (v && *v) {
		id = std::string("slurm-") + v;
		const char* st = std::getenv("SLURM_STEP_ID");
		if (st && *st) id += std::string(".") + st;
		scope = "launcher";
		return;
	}
	v = std::getenv("OMPI_MCA_ess_base_jobid");
	if (v && *v) { id = std::string("ompi-") + v; scope = "launcher"; return; }
	std::stringstream ss;
	ss << "proc-" << std::hex
	   << (unsigned long long)std::chrono::system_clock::now().time_since_epoch().count()
	   << "-" << (unsigned long long)std::random_device{}();
	id = ss.str();
	scope = "process";
}

// The user stage barrier. The interpreters' user `domainChange` command
// (tcl/commands.cpp, interpreter/OpenSeesMiscCommands.cpp OPS_domainChange)
// calls ladrunoNoteUserStageBarrier(); the recorder compares the count at each
// recorded step. Internal Domain::domainChange() callers (an SP pattern,
// eleLoad, a contact re-sort, ...) do not touch it, so those re-stamps are
// folded into the open MODEL_STAGE. The count lives in this TU
// (the OPS_Recorder library every interpreter executable links), not in the
// Domain, and is process-local: on openseesmp every rank runs the same command.
static int s_ladruno_user_stage_barrier = 0;

void ladrunoNoteUserStageBarrier(void)
{
	++s_ladruno_user_stage_barrier;
}

int ladrunoUserStageBarrier(void)
{
	return s_ladruno_user_stage_barrier;
}

// The transient analysis the classic Tcl `eigen` command builds
// when none is defined (tcl/commands.cpp eigenAnalysis). openseespy deletes
// that analysis after the solve; classic Tcl keeps it, so an `eigen ...;
// record` snapshot would read "transient" in Tcl and nothing in openseespy.
// While the transient pointer is still this object, it does not count as a
// running analysis. `analysis` and `analyze` clear the mark (an analyze step
// on it is a real transient run), as do `wipe` and `wipeAnalysis`.
static const void* s_ladruno_eigen_only_analysis = 0;

void ladrunoNoteEigenOnlyAnalysis(const void* analysis)
{
	s_ladruno_eigen_only_analysis = analysis;
}

// The kind of the analysis running on this process right now:
// "static", "transient", or "" when it cannot be told (no analysis defined, or
// both pointers set). The interpreter-side analyses are read through the
// elementAPI accessors (classic Tcl: elementAPI_TCL.cpp; openseespy and the
// new Tcl wrapper: OpenSeesCommands.cpp). VariableTransient and PFEM set the
// transient pointer too. A PartitionedDomain (OpenSeesSP) worker runs no
// interpreter, so it sees neither pointer: its Subdomain holds the
// DomainDecompositionAnalysis P0 sent, but Subdomain::getDDAnalysis() is
// protected. There the kind stays unresolved (provisional "static", warned at
// stage close).
static std::string ladrunoRunningAnalysisKind()
{
	StaticAnalysis** sa = OPS_GetStaticAnalysis();
	DirectIntegrationAnalysis** ta = OPS_GetTransientAnalysis();
	const bool is_static = (sa != 0 && *sa != 0);
	const bool is_transient = (ta != 0 && *ta != 0 &&
	                           (const void*)(*ta) != s_ladruno_eigen_only_analysis);
	if (is_static && !is_transient)
		return "static";
	if (is_transient && !is_static)
		return "transient";
	return "";   // none, or both (ambiguous): unresolved, never a flip
}

// true when a modal (ModesOfVibration[Rotational]) channel is requested: such a
// recorder's stages are KIND "eigen" whatever analysis runs.
static bool ladrunoHasModalRequest(
	const std::vector<ladruno::detail::NodalResultType::Enum>& reqs)
{
	for (size_t i = 0; i < reqs.size(); ++i) {
		if (reqs[i] == ladruno::detail::NodalResultType::ModesOfVibration ||
		    reqs[i] == ladruno::detail::NodalResultType::ModesOfVibrationRotational)
			return true;
	}
	return false;
}

// Replace MODEL_STAGE[<id>]/@KIND (a provisional value written
// before any analysis existed). Delete + create: the string sizes differ.
static int ladrunoRewriteStageKind(hid_t h_file, int stage_id, const std::string& kind)
{
	std::stringstream ss;
	ss << "MODEL_STAGE[" << stage_id << "]";
	hid_t h_stage = H5Gopen2(h_file, ss.str().c_str(), H5P_DEFAULT);
	if (h_stage < 0)
		return -1;
	if (H5Aexists(h_stage, "KIND") > 0)
		H5Adelete(h_stage, "KIND");
	herr_t st = ladruno::h5::attribute::write(h_stage, "KIND", kind);
	H5Gclose(h_stage);
	return st < 0 ? -1 : 0;
}

// max number of iterations to guess the number of fibers (MPCO_MAX_TRIAL_NFIB in MPCORecorder)
#ifndef LADRUNO_MAX_TRIAL_NFIB
#define LADRUNO_MAX_TRIAL_NFIB 100000
#endif

/* ===================================================================== */
/* private_data                                                          */
/* ===================================================================== */

class LadrunoRecorder::private_data
{
public:
	private_data()
		: filename()
		, send_self_count(0)
		, p_id(0)
		, is_partitioned(false)
		, initialized(false)
		, first_domain_changed_done(false)
		, stage_failed(false)
		, init_failed(false)
		, flush_seconds(10.0)
		, last_flush()
		, last_seen_stamp(0)
		, topology_fp(0)
		, last_barrier(0)
		, stage_run_kind()
		, kind_pending(false)
		, worker_fold_warned(false)
		, energy_model_state()
		, energy_region_state()
		, run_id()
		, run_id_scope()
		, info()
		, output_freq()
		, has_region(false)
		, region_tag(0)
		, node_set()
		, elem_set()
		, energy_requested(false)
		, energy_region_tags()
		, stage_kind()
		, envelope_mode(false)
		, nodal_results_requests()
		, sens_grad_indices()
		, elemental_results_requests()
		, nodes()
		, elements()
		, elem_ngauss_nfiber_info()
		, elemental_responses()
	{}

	std::string filename;
	// Parallel partition bookkeeping (MPCORecorder convention):
	//   send_self_count == 0  -> sequential (one process, no partitioning)
	//   send_self_count  > 0  -> primary P0 (incremented once per worker channel)
	//   send_self_count  < 0  -> a secondary/worker process (set in recvSelf)
	// p_id is this process's partition index: 0 on P0, 1..N-1 on the workers
	// (the value P0's counter held when it serialized to that worker). The output
	// file becomes "<stem>.part-<p_id>.ladruno" whenever send_self_count != 0.
	int send_self_count;
	int p_id;
	// resolved in initialize() (either signal above); read at channel
	// creation to refuse -envelope of a source whose partitions are partials.
	bool is_partitioned;
	bool initialized;
	// first_domain_changed_done: false until the first record() resolves the
	// model-stage stamp; thereafter a stamp change triggers a model rebuild.
	bool first_domain_changed_done;
	// true after a MODEL_STAGE rebuild failed. writeModel() released
	// every source first, so nothing dangles; recording stays suspended (no
	// results, no per-step error spam) until the next stamp change rebuilds.
	bool stage_failed;
	// true after initialize() failed once: the recorder is then dead for the
	// run instead of retrying (and leaking) every step.
	bool init_failed;
	// wall-clock flush cadence (-flush <seconds>, 0 = every recorded step). A
	// per-step H5Fflush forces every dirty chunk to be deflated and written on
	// every step; between flushes HDF5's chunk cache keeps them.
	double flush_seconds;
	std::chrono::steady_clock::time_point last_flush;
	// the domain-change stamp moves for far more than topology (a contact
	// re-sort inside commit(), `eleLoad`, a pattern with SPs, an
	// `sp`/`imposedMotion`), and a new MODEL_STAGE is a full model copy plus an
	// envelope and energy reset. A stamp move rebuilds only when this
	// fingerprint of the node/element/pressure-constraint SET changed.
	int last_seen_stamp;
	unsigned long long topology_fp;
	// the user stage-barrier count (ladrunoUserStageBarrier()) when the
	// current MODEL_STAGE opened. A user `domainChange`
	// bumps it, and a moved count opens a new MODEL_STAGE even when the
	// topology fingerprint is unchanged. Process-local, never serialized.
	int last_barrier;
	// the running analysis kind ("static" /
	// "transient") the current MODEL_STAGE was opened under; "" while it is not
	// resolved yet (no analysis defined at the stage's first record). A
	// static<->transient flip opens a new MODEL_STAGE whatever the KIND label
	// (-kind, eigen), so no stage mixes static and transient steps.
	std::string stage_run_kind;
	// true while the open stage's KIND attr is a provisional "static" written
	// before any analysis existed; the first record() that sees the running
	// analysis rewrites it.
	bool kind_pending;
	// one-time warning on a PartitionedDomain worker, which cannot see the user
	// barrier (no interpreter runs there).
	bool worker_fold_warned;
	// energy integrals owned here, not by the per-stage source.
	ladrunons::EnergyState energy_model_state;
	ladrunons::EnergyState energy_region_state;
	// run identity written to INFO, so a reader can tell the
	// part files of ONE run from stale parts of an earlier one. Broadcast by P0
	// on the PartitionedDomain path (set in sendSelf, received in recvSelf).
	std::string run_id;
	std::string run_id_scope;
	ladrunons::detail::ProcessInfo info;

	// -T output frequency
	ladrunons::detail::OutputFrequency output_freq;

	// -R region restriction
	bool has_region;
	int region_tag;              // the -R region tag (0 if none); for MODEL/SETS
	std::vector<int> node_set;
	std::vector<int> elem_set;

	// -G energy balance: whole-model and/or per-region energy.
	bool energy_requested;               // -G energy   -> whole-model ON_DOMAIN
	std::vector<int> energy_region_tags; // -G <tag...> -> per-region ON_REGIONS

	// -kind <transient|static|eigen>: an explicit override of the MODEL_STAGE/KIND
	// LABEL for every stage; it does not stop a static<->transient flip from
	// opening a stage. Empty (the default) = automatic: each
	// stage's KIND is the analysis running at its first record() (Static ->
	// "static"; Transient / VariableTransient / PFEM -> "transient"). Either way,
	// a modal (ModesOfVibration) channel promotes KIND to "eigen" at writeModel.
	std::string stage_kind;              // "" = automatic, else the -kind value

	// -envelope: record componentwise running extremes instead of the
	// per-step time series. Every non-modal node/element/domain channel uses an
	// EnvelopeSink (MIN/MAX/ABSMAX + ARG_STEP under RESULTS/ENVELOPES) rather than a
	// StreamingSink. Off by default (time-series output unchanged).
	bool envelope_mode;

	// parsed requests
	std::vector<ladrunons::detail::NodalResultType::Enum> nodal_results_requests;
	std::vector<int> sens_grad_indices;
	std::vector<std::vector<std::string> > elemental_results_requests;

	// the node set (Node* in write order) used by all node sources
	std::vector<Node*> nodes;

	// element engine
	ladrunons::detail::element::ElementCollection elements;
	// for each element tag: per gauss point (fiber_base_index, num_fibers)
	std::map<int, std::vector<std::pair<int, int> > > elem_ngauss_nfiber_info;
	// all element responses, owned here (released in clearSources)
	std::vector<Response*> elemental_responses;

	/* --- node sources/sinks --- */
	struct NodeChannel {
		ladrunons::ResultSource* source;   // owned
		ladrunons::ResultSink* sink;       // owned (StreamingSink or, in -envelope mode, EnvelopeSink)
		int reaction_flag;               // -1 if not a reaction source
		bool is_modes;                   // ModesOfVibration(Rotational) special path
		NodeChannel() : source(0), sink(0), reaction_flag(-1), is_modes(false) {}
	};
	std::vector<NodeChannel> node_channels;

	/* --- element sources/sinks --- */
	struct ElemChannel {
		ladrunons::ElementResultSource* source;             // owned
		ladrunons::ResultSink* sink;                        // owned (Streaming or Envelope)
		const ladrunons::detail::element::OutputDescriptorHeader* header; // not owned
		bool column_map_written;
		ElemChannel() : source(0), sink(0), header(0), column_map_written(false) {}
	};
	std::vector<ElemChannel> elem_channels;
	// the per-request response-collection trees (own the OutputResponseCollection
	// buckets referenced by ElementResultSource); kept alive for source lifetime.
	ladrunons::detail::element::ResultRecorderCollection elemental_recorders;

	/* --- domain/region sources/sinks (energy balance) --- */
	struct DomainChannel {
		ladrunons::ResultSource* source;   // owned
		ladrunons::ResultSink* sink;       // owned (Streaming or Envelope)
		DomainChannel() : source(0), sink(0) {}
	};
	std::vector<DomainChannel> domain_channels;
};

/* ===================================================================== */
/* ctor / dtor                                                           */
/* ===================================================================== */

LadrunoRecorder::LadrunoRecorder()
	: Recorder(RECORDER_TAGS_LadrunoRecorder)
	, m_data(new private_data())
{
}

LadrunoRecorder::~LadrunoRecorder()
{
	if (m_data) {
		// the last stage never saw a running analysis.
		if (m_data->initialized && m_data->kind_pending)
			opserr << "LadrunoRecorder warning: the last MODEL_STAGE closed with no "
			          "analysis visible at any recorded step; its KIND stays the "
			          "provisional 'static' (pass -kind to state it)\n";
		// In envelope mode, write the final ENVELOPES datasets while the sinks +
		// file are still alive (EnvelopeSink defers all output to finalize()).
		// Inside H5E_BEGIN_TRY like the closes below: at process exit HDF5's
		// atexit handler may already have closed every ID, and this write would
		// print HDF5-DIAG stacks (no data is lost: the envelope was rewritten on
		// the last recorded step).
		if (m_data->envelope_mode && m_data->initialized &&
		    m_data->info.h_file_id != ladrunons::HID_INVALID) {
			H5E_BEGIN_TRY {
				finalizeAllSinks();
			} H5E_END_TRY;
		}
		// sinks hold DATA/TIME/STEP and envelope handles; closing them at
		// process exit (IDs already closed by HDF5's atexit) must be silent too.
		H5E_BEGIN_TRY {
			clearSources();
		} H5E_END_TRY;
		if (m_data->initialized && m_data->info.h_file_id != ladrunons::HID_INVALID) {
			// When the recorder dies at process exit (deck never wiped), HDF5's
			// own atexit handler has usually run first: it flushes + closes every
			// open ID (that is why the data survives) and the closes below then
			// operate on dead IDs, printing 4 HDF5-DIAG stacks per rank to stderr.
			// H5E_BEGIN_TRY is HDF5's canonical
			// attempt-silently wrapper for exactly this shutdown-path situation;
			// on a live library (normal `wipe`) these closes behave as before.
			H5E_BEGIN_TRY {
				ladrunons::h5::file::flush(m_data->info.h_file_id);
				ladrunons::h5::file::close(m_data->info.h_file_id);
				ladrunons::h5::plist::close(m_data->info.h_group_proplist);
				ladrunons::h5::plist::close(m_data->info.h_file_proplist);
			} H5E_END_TRY;
		}
		delete m_data;
	}
}

int LadrunoRecorder::clearSources()
{
	for (size_t i = 0; i < m_data->node_channels.size(); ++i) {
		delete m_data->node_channels[i].source;
		delete m_data->node_channels[i].sink;
	}
	m_data->node_channels.clear();
	for (size_t i = 0; i < m_data->elem_channels.size(); ++i) {
		delete m_data->elem_channels[i].source;
		delete m_data->elem_channels[i].sink;
	}
	m_data->elem_channels.clear();
	for (size_t i = 0; i < m_data->domain_channels.size(); ++i) {
		delete m_data->domain_channels[i].source;
		delete m_data->domain_channels[i].sink;
	}
	m_data->domain_channels.clear();
	// element responses are owned here
	for (size_t i = 0; i < m_data->elemental_responses.size(); ++i) {
		if (m_data->elemental_responses[i])
			delete m_data->elemental_responses[i];
	}
	m_data->elemental_responses.clear();
	m_data->elemental_recorders.clear();
	return 0;
}

void LadrunoRecorder::finalizeAllSinks()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	for (size_t i = 0; i < m_data->node_channels.size(); ++i)
		if (m_data->node_channels[i].sink)
			m_data->node_channels[i].sink->finalize(info);
	for (size_t i = 0; i < m_data->elem_channels.size(); ++i)
		if (m_data->elem_channels[i].sink)
			m_data->elem_channels[i].sink->finalize(info);
	for (size_t i = 0; i < m_data->domain_channels.size(); ++i)
		if (m_data->domain_channels[i].sink)
			m_data->domain_channels[i].sink->finalize(info);

	// Element ENVELOPES carry the same structured per-column COLUMN_MAP as the
	// time-series ON_ELEMENTS path (gauss/section/fiber column structure, beyond
	// the result-level COMPONENTS, which is empty for element results). The
	// EnvelopeSink creates its leaf group once and overwrites it in place, so
	// the COLUMN_MAP is written ONCE, right after the group first exists.
	// finalizeAllSinks is only ever called in envelope mode; the envelope_mode
	// guard keeps it correct if that changes.
	if (m_data->envelope_mode) {
		for (size_t i = 0; i < m_data->elem_channels.size(); ++i) {
			private_data::ElemChannel& ech = m_data->elem_channels[i];
			if (ech.column_map_written)
				continue;
			ladrunons::EnvelopeSink* es = dynamic_cast<ladrunons::EnvelopeSink*>(ech.sink);
			if (es != 0 && es->written()) {
				writeElementColumnMap((int)i);
				ech.column_map_written = true;
			}
		}
	}
}

/* ===================================================================== */
/* lifecycle                                                             */
/* ===================================================================== */

int LadrunoRecorder::setDomain(Domain& theDomain)
{
	m_data->info.domain = &theDomain;
	return 0;
}

int LadrunoRecorder::restart(void)
{
	return 0;
}

int LadrunoRecorder::domainChanged(void)
{
	return 0;
}

int LadrunoRecorder::record(int commitTag, double timeStamp)
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	// a recorder whose file could not be created is dead for
	// the run (the error was printed once by initialize()).
	if (m_data->init_failed)
		return 0;
	info.current_time_step_id = commitTag;
	info.current_time_step = timeStamp;

	// integrate the energy balance on EVERY commit, before the -T
	// gate, so IE/DW/ULW are trapezoids over the real steps, not over the -T
	// samples. advance() is idempotent per commit tag (the recorded-step
	// evaluate() calls it again) and sweeps the live domain (no cached
	// Element*), so it is safe even on a commit where the topology changed and
	// the stage has not been rebuilt yet.
	if (m_data->initialized && !m_data->stage_failed) {
		for (size_t i = 0; i < m_data->domain_channels.size(); ++i) {
			ladrunons::EnergyBalanceSource* es =
				dynamic_cast<ladrunons::EnergyBalanceSource*>(m_data->domain_channels[i].source);
			if (es)
				es->advance(info);
		}
	}

	// -T output-frequency gate.
	bool do_record = false;
	if (!m_data->initialized) {
		do_record = true; // always record the first time
	}
	else {
		if (m_data->output_freq.type == ladrunons::detail::OutputFrequency::NumberOfSteps) {
			if ((info.current_time_step_id - m_data->output_freq.last_step) >= m_data->output_freq.nsteps)
				do_record = true;
		}
		else if (m_data->output_freq.type == ladrunons::detail::OutputFrequency::DeltaTime) {
			// relative tolerance (as NodeRecorder's
			// relDeltaTTol), so a time accumulated by summing steps
			// (0.19999999999999996 vs 0.2) is not skipped; a time that moved
			// BACKWARD (setTime between stages) records and restarts the grid.
			const double dt = m_data->output_freq.dt;
			const double elapsed = info.current_time_step - m_data->output_freq.last_time;
			if (elapsed >= dt - 1.0e-5 * dt || elapsed < -1.0e-5 * dt)
				do_record = true;
		}
	}
	if (do_record) {
		m_data->output_freq.last_step = info.current_time_step_id;
		// snap `last_time` to the -T dt grid (last + k*dt)
		// instead of the actual recorded time, so the record times do not drift
		// by one analysis step per sample (dt_analysis 0.03, -T dt 0.1 recorded
		// at 0.12, 0.24, ... instead of on the nearest steps after 0.1, 0.2, ...).
		const double dt = m_data->output_freq.dt;
		const double elapsed = info.current_time_step - m_data->output_freq.last_time;
		if (m_data->initialized &&
		    m_data->output_freq.type == ladrunons::detail::OutputFrequency::DeltaTime &&
		    dt > 0.0 && elapsed >= 0.0) {
			const double k = std::floor((elapsed + 1.0e-5 * dt) / dt);
			m_data->output_freq.last_time += (k >= 1.0 ? k : 1.0) * dt;
		}
		else {
			m_data->output_freq.last_time = info.current_time_step;
		}
	}
	else {
		return 0;
	}

	if (!m_data->initialized) {
		if (initialize() != 0) {
			m_data->init_failed = true;
			opserr << "LadrunoRecorder: initialization failed; this recorder is "
			          "disabled for the rest of the run\n";
			return -1;
		}
		m_data->initialized = true;
	}

	// Multi-stage detection. On the first record, and again whenever the
	// domain-change stamp moves (a new MODEL_STAGE: wipeAnalysis + new
	// elements/nodes, staged analysis),
	// (re)write the model and rebuild every node/element source+sink. This is
	// what re-acquires fresh Element*/Response* pointers; without it the cached
	// pointers from the prior stage dangle and the new stage is never written.
	// Multi-stage stamp (per-process). MPCORecorder max-reduces this across ranks
	// via an MPI_Allreduce, but that is its only _PARALLEL_PROCESSING-gated code
	// and is compiled out of this shared sequential TU, so the local domain stamp
	// is used. Each partition writes its own part-N file, so cross-rank stamp
	// sync only matters for aligning MODEL_STAGE ids across partitions in a staged
	// analysis; that is a documented limitation here (single-stage is unaffected).
	bool rebuild_model = false;
	int new_stamp = info.domain->hasDomainChanged();
	// The two stage boundaries that leave topology unchanged.
	// (B) a user `domainChange` since the stage opened (the barrier count);
	// (A) a static<->transient flip of the running analysis, so that no
	// MODEL_STAGE mixes static and transient steps. The flip splits whatever
	// the KIND label is (-kind and eigen promotion only set the label). An
	// unresolved kind ("": no analysis yet) is never a flip; a repeated kind is
	// never a flip.
	const int barrier_now = ladrunoUserStageBarrier();
	const std::string run_kind = ladrunoRunningAnalysisKind();
	bool user_barrier = false;
	bool kind_flip = false;
	if (m_data->first_domain_changed_done && !m_data->stage_failed) {
		user_barrier = (barrier_now != m_data->last_barrier);
		kind_flip = !run_kind.empty() &&
		            !m_data->stage_run_kind.empty() &&
		            run_kind != m_data->stage_run_kind;
	}
	if (!m_data->first_domain_changed_done) {
		info.current_model_stage_id = new_stamp;
		m_data->last_seen_stamp = new_stamp;
		m_data->first_domain_changed_done = true;
		rebuild_model = true;
	}
	else if (!user_barrier && !kind_flip &&
	         new_stamp != m_data->last_seen_stamp &&
	         !m_data->stage_failed && topologyFingerprint() == m_data->topology_fp) {
		// the stamp moved but the node/element set did not (contact re-sort,
		// eleLoad, an SP pattern, ...): keep the MODEL_STAGE, its sources,
		// envelopes and energy integrals. Every cached Element*/Node* is still the
		// live object (the fingerprint hashes the pointers, not just the tags).
		m_data->last_seen_stamp = new_stamp;
		// a PartitionedDomain worker runs no interpreter, so a
		// user `domainChange` on P0 reaches it only as a stamp move; it cannot
		// be told from an internal re-stamp and is folded here. Say so once.
		if (m_data->send_self_count < 0 && !m_data->worker_fold_warned) {
			m_data->worker_fold_warned = true;
			opserr << "LadrunoRecorder warning: partition " << m_data->p_id
			       << " folded a domain-change stamp into MODEL_STAGE["
			       << info.current_model_stage_id << "]: a PartitionedDomain "
			          "worker cannot see a user `domainChange` stage barrier, so "
			          "this part file may hold fewer MODEL_STAGEs than part 0\n";
		}
	}
	else if (user_barrier || kind_flip || new_stamp != m_data->last_seen_stamp) {
		m_data->last_seen_stamp = new_stamp;
		// close out the ending stage BEFORE the stamp moves: envelopes are
		// rewritten on the flush cadence, not every step, so their last
		// extremes (and the element COLUMN_MAP, resolved from the stage stamp)
		// must land under the OLD MODEL_STAGE; then flush it to disk.
		if (!m_data->stage_failed) {
			if (m_data->envelope_mode)
				finalizeAllSinks();
			flushFile();
		}
		if (m_data->kind_pending) {
			opserr << "LadrunoRecorder warning: MODEL_STAGE["
			       << info.current_model_stage_id << "] closed with no analysis "
			          "visible at any recorded step; its KIND stays the "
			          "provisional 'static' (pass -kind to state it)\n";
		}
		// a barrier or a kind flip can arrive without a stamp
		// move, so the id is the stamp only when it is past the current id;
		// otherwise current+1. Every later stamp move takes max(id+1, stamp)
		// again, so the ids stay strictly increasing and never collide.
		info.current_model_stage_id =
			std::max(info.current_model_stage_id + 1, new_stamp);
		rebuild_model = true;
	}
	if (rebuild_model) {
		m_data->last_barrier = barrier_now;
		// writeModel() releases the prior stage's sources (and their
		// cached Element*/Response*) BEFORE writing anything, so a failed rebuild
		// leaves no dangling pointers behind. The stamp is already committed, so
		// the next record() would not retry: latch stage_failed and record nothing
		// until the next stamp change rebuilds (one error, not one per step).
		if (writeModel() != 0) {
			m_data->stage_failed = true;
			opserr << "LadrunoRecorder: MODEL_STAGE[" << info.current_model_stage_id
			       << "] could not be written; recording is suspended until the "
			          "next domain change\n";
			return -1;
		}
		m_data->stage_failed = false;
		m_data->topology_fp = topologyFingerprint();
	}
	if (m_data->stage_failed)
		return 0;

	// the stage opened before the running analysis could be told
	// (an explicit `record` first, an eigen snapshot): adopt the first kind a
	// recorded step sees. An automatic KIND label was written as a provisional
	// "static"; restate it.
	if (m_data->stage_run_kind.empty() && !run_kind.empty()) {
		m_data->stage_run_kind = run_kind;
		if (m_data->kind_pending) {
			if (ladrunoRewriteStageKind(info.h_file_id, info.current_model_stage_id,
			                            run_kind) != 0) {
				opserr << "LadrunoRecorder: could not rewrite MODEL_STAGE["
				       << info.current_model_stage_id << "]/KIND\n";
				return -1;
			}
			m_data->kind_pending = false;
		}
	}

	if (recordResultsOnNodes() != 0)
		return -1;
	if (recordResultsOnElements() != 0)
		return -1;
	if (recordResultsOnDomain() != 0)
		return -1;

	// on the -flush cadence (default 10 s of wall clock; 0 = every recorded
	// step) rewrite the envelopes in place and flush the file. Crash safety: at
	// most the last `-flush` seconds of output are lost; a per-step flush
	// forces every dirty chunk to be deflated and written on every step, which
	// dominates the cost of small or explicit runs.
	const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
	const double since = std::chrono::duration<double>(now - m_data->last_flush).count();
	const bool due = (m_data->flush_seconds <= 0.0 || since >= m_data->flush_seconds);
	if (m_data->envelope_mode) {
		// Small envelopes (<= 8 MiB of accumulators in total) are still
		// overwritten every recorded step (a plain in-place H5Dwrite, no flush),
		// so a deck that exits without `wipe` keeps its latest extremes (HDF5's
		// atexit close writes them out). Larger ones follow the -flush cadence.
		size_t env_bytes = 0;
		for (size_t i = 0; i < m_data->node_channels.size(); ++i) {
			ladrunons::EnvelopeSink* es =
				dynamic_cast<ladrunons::EnvelopeSink*>(m_data->node_channels[i].sink);
			if (es) env_bytes += es->payloadBytes();
		}
		for (size_t i = 0; i < m_data->elem_channels.size(); ++i) {
			ladrunons::EnvelopeSink* es =
				dynamic_cast<ladrunons::EnvelopeSink*>(m_data->elem_channels[i].sink);
			if (es) env_bytes += es->payloadBytes();
		}
		if (due || env_bytes <= ((size_t)8 << 20))
			finalizeAllSinks();
	}
	if (due)
		flushFile();

	return 0;
}

unsigned long long LadrunoRecorder::topologyFingerprint() const
{
	// FNV-1a over the node, element and pressure-constraint sets:
	// (tag, object address) per entry, plus the counts. Address, not only tag:
	// `remove element 5; element ... 5` re-uses the tag with a NEW object, and
	// the cached Response* of the old one must not survive. O(N) per stamp move
	// (not per step), cheap next to the full rebuild it avoids.
	const Domain* dom_c = m_data->info.domain;
	if (dom_c == 0)
		return 0ULL;
	Domain* dom = const_cast<Domain*>(dom_c);
	unsigned long long h = 1469598103934665603ULL;
	auto mix = [&h](unsigned long long v) {
		for (int b = 0; b < 8; ++b) {
			h ^= (v >> (8 * b)) & 0xFFULL;
			h *= 1099511628211ULL;
		}
	};
	mix((unsigned long long)dom->getNumNodes());
	mix((unsigned long long)dom->getNumElements());
	mix((unsigned long long)dom->getNumPCs());
	{
		NodeIter& it = dom->getNodes();
		Node* n = 0;
		while ((n = it()) != 0) {
			mix((unsigned long long)(long long)n->getTag());
			mix((unsigned long long)(uintptr_t)n);
		}
	}
	{
		ElementIter& it = dom->getElements();
		Element* e = 0;
		while ((e = it()) != 0) {
			mix((unsigned long long)(long long)e->getTag());
			mix((unsigned long long)(uintptr_t)e);
		}
	}
	{
		Pressure_ConstraintIter& it = dom->getPCs();
		Pressure_Constraint* pc = 0;
		while ((pc = it()) != 0)
			mix((unsigned long long)(uintptr_t)pc);
	}
	return h;
}

void LadrunoRecorder::flushFile()
{
	// the one place the recorder flushes (cadence, stage end, close).
	if (m_data->initialized && m_data->info.h_file_id != ladrunons::HID_INVALID)
		ladrunons::h5::file::flush(m_data->info.h_file_id);
	m_data->last_flush = std::chrono::steady_clock::now();
}

/* ===================================================================== */
/* file creation                                                         */
/* ===================================================================== */

int LadrunoRecorder::initialize()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;

	if (info.domain == 0) {
		opserr << "LadrunoRecorder error: domain is not defined\n";
		return -1;
	}

	// spatial dimension from the first node (enforce >0 and <=3)
	if (info.num_dimensions < 1) {
		NodeIter& nit = info.domain->getNodes();
		Node* nptr = 0;
		while ((nptr = nit()) != 0) {
			info.num_dimensions = nptr->getCrds().Size();
			break;
		}
		if (info.num_dimensions < 1 || info.num_dimensions > 3) {
			opserr << "LadrunoRecorder error: invalid spatial dimension "
			       << info.num_dimensions << "\n";
			return -1;
		}
	}

	// property lists (link creation order tracked/indexed, as in MPCORecorder)
	info.h_file_proplist = ladrunons::h5::plist::crate(ladrunons::h5::plist::FileCreate);
	ladrunons::h5::plist::setLinkCreationOrder(
		info.h_file_proplist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
	info.h_group_proplist = ladrunons::h5::plist::crate(ladrunons::h5::plist::GroupCreate);
	ladrunons::h5::plist::setLinkCreationOrder(
		info.h_group_proplist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

	// Resolve this process's partition identity for the parallel output contract
	// (no MPI: this TU is not _PARALLEL_PROCESSING; see the note near the includes).
	//   * PartitionedDomain/OpenSeesMP: the recorder was broadcast from P0
	//     (send_self_count != 0) and p_id is the assigned 0-based index;
	//   * interpreter-per-rank/openseesmp: each rank built the recorder itself
	//     (send_self_count == 0); index by the launcher's per-rank environment.
	//     Launchers disagree on the names: Intel MPI / MS-MPI export
	//     PMI_RANK+PMI_SIZE, OpenMPI exports OMPI_COMM_WORLD_RANK+_SIZE, and
	//     SLURM's srun exports SLURM_PROCID+SLURM_NTASKS under every MPI
	//     plugin, including pmix, where no PMI_* exist at all (without it every
	//     rank of a `srun --mpi=pmix_v3` job opens the same file). First pair with
	//     SIZE > 1 wins; SLURM is probed last because sbatch also exports
	//     SLURM_NTASKS into the batch shell itself, so a real MPI launcher's
	//     own variables should take precedence when both are present.
	//     (PMIX_RANK exists under pmix but has no reliable SIZE companion.)
	// Single process (no launcher env, or SIZE <= 1) -> not partitioned -> verbatim.
	int part_id = m_data->p_id;
	bool is_partitioned = (m_data->send_self_count != 0);
	int num_parts = 1;
	{
		// one probe (Ladruno_LaunchEnv.h): SLURM only inside an srun step, and
		// a SIZE > 1 without a valid RANK is refused
		// instead of every rank silently writing (truncating) part-0.
		int env_rank = 0, env_size = 1;
		std::string env_source, env_error;
		const ladrunons::launch::Status st =
			ladrunons::launch::detectRank(env_rank, env_size, env_source, env_error);
		if (st == ladrunons::launch::Inconsistent) {
			opserr << "LadrunoRecorder error: " << env_error.c_str()
			       << ": cannot tell which partition file this rank owns\n";
			ladrunons::h5::plist::close(info.h_group_proplist);
			ladrunons::h5::plist::close(info.h_file_proplist);
			info.h_group_proplist = ladrunons::HID_INVALID;
			info.h_file_proplist = ladrunons::HID_INVALID;
			return -1;
		}
		if (st == ladrunons::launch::Launched) {
			num_parts = env_size;
			if (!is_partitioned) {   // interpreter-per-rank: index by the launcher rank
				part_id = env_rank;
				is_partitioned = true;
			}
		}
	}
	m_data->p_id = part_id;
	m_data->is_partitioned = is_partitioned;

	// Per-partition output filename: "<stem>.part-<part_id>.ladruno" when this is
	// one rank of a partitioned run (so N processes never race to write one
	// file), else the filename verbatim. A reader stitches the set via
	// ^(?P<stem>.+?)\.part-(?P<idx>\d+)\.ladruno$ (0-based, contiguous).
	std::string the_filename = m_data->filename;
	if (is_partitioned) {
		std::string stem = m_data->filename;
		const std::string ext = ".ladruno";
		if (stem.size() >= ext.size() &&
		    stem.compare(stem.size() - ext.size(), ext.size(), ext) == 0)
			stem.erase(stem.size() - ext.size());
		std::stringstream ss;
		ss << stem << ".part-" << part_id << ".ladruno";
		the_filename = ss.str();
	}

	// create the file (no SWMR; H5P_DEFAULT file access)
	info.h_file_id = ladrunons::h5::file::create(
		the_filename.c_str(), info.h_file_proplist, H5P_DEFAULT);
	if (info.h_file_id == ladrunons::HID_INVALID) {
		opserr << "LadrunoRecorder error: cannot create file \""
		       << the_filename.c_str() << "\"\n";
		// release the property lists created above.
		ladrunons::h5::plist::close(info.h_group_proplist);
		ladrunons::h5::plist::close(info.h_file_proplist);
		info.h_group_proplist = ladrunons::HID_INVALID;
		info.h_file_proplist = ladrunons::HID_INVALID;
		return -1;
	}

	// INFO group
	hid_t h_info = ladrunons::h5::group::create(
		info.h_file_id, "INFO", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::attribute::write(h_info, "GENERATOR", std::string("Ladruno"));
	ladrunons::h5::attribute::write(h_info, "FORMAT_VERSION", (int)1);
	// minor 1 = per-stage KIND from the running analysis, and a
	// MODEL_STAGE also opens on a user `domainChange` or a static<->transient
	// flip. Absent (minor 0): KIND is the recorder-wide -kind (default
	// "static") and a stage opens only on a topology change.
	ladrunons::h5::attribute::write(h_info, "FORMAT_MINOR_VERSION", (int)1);
	ladrunons::h5::attribute::write(h_info, "SOLVER_NAME", std::string("OpenSees"));
	std::vector<int> solver_version;
	solver_version.push_back(3);
	solver_version.push_back(5);
	solver_version.push_back(1);
	ladrunons::h5::attribute::write(h_info, "SOLVER_VERSION", solver_version);
	ladrunons::h5::attribute::write(h_info, "SPATIAL_DIM", info.num_dimensions);
	// On-disk precision of per-step result DATA: "f64" (default, lossless) or
	// "f32" (opt-in `-precision f32` lossy mode, about 7 significant digits).
	// -precision f32 only affects the streaming per-step DATA datasets.
	// Envelope-mode output (MIN/MAX/ABSMAX) is always f64, so STORED_PRECISION
	// reports f64 there: the attribute describes what is actually on disk, and
	// a consumer relies on it to choose its comparison tolerance.
	bool effective_f32 = info.store_data_f32 && !m_data->envelope_mode;
	if (info.store_data_f32 && m_data->envelope_mode)
		opserr << "LadrunoRecorder warning: -precision f32 is ignored in -envelope mode "
		          "(envelope datasets are always f64); recording STORED_PRECISION=f64\n";
	ladrunons::h5::attribute::write(h_info, "STORED_PRECISION",
		std::string(effective_f32 ? "f32" : "f64"));

	// Partition manifest (parallel-output contract). Each ".part-N.ladruno" file
	// self-declares its 0-based partition index + the total count so a reader can
	// validate the contiguous-from-0 set it globs. PARTITIONED=0 marks an ordinary
	// single-file (sequential) output; NUM_PARTITIONS then stays 1. Values resolved
	// above (broadcast p_id for PartitionedDomain, MPI rank for interpreter-per-rank).
	ladrunons::h5::attribute::write(h_info, "PARTITION_ID", part_id);
	ladrunons::h5::attribute::write(h_info, "PARTITIONED", (int)(is_partitioned ? 1 : 0));
	ladrunons::h5::attribute::write(h_info, "NUM_PARTITIONS", num_parts);
	// RUN_ID + RUN_ID_SCOPE ("broadcast" | "user" | "launcher" |
	// "process"). Part files of one run carry the same RUN_ID unless the scope
	// is "process" (no shared identity was available on this launch path).
	if (m_data->run_id.empty())
		ladrunoResolveRunId(m_data->run_id, m_data->run_id_scope);
	ladrunons::h5::attribute::write(h_info, "RUN_ID", m_data->run_id);
	ladrunons::h5::attribute::write(h_info, "RUN_ID_SCOPE", m_data->run_id_scope);

	hid_t h_prov = ladrunons::h5::group::create(
		h_info, "PROVENANCE", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::group::close(h_prov);
	ladrunons::h5::group::close(h_info);

	// NOTE: model writing + source/sink building are NOT done here. They live in
	// writeModel(), driven by record()'s multi-stage rebuild_model block, so they
	// re-run on every MODEL_STAGE change (initialize() only sets up the file and
	// record() owns the model rebuild).

	return 0;
}

/* ===================================================================== */
/* model writing                                                         */
/* ===================================================================== */

int LadrunoRecorder::writeModel()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;

	// info.current_model_stage_id is set by record() before calling writeModel()
	// (the multi-stage rebuild_model block), so the stamp is consistent across
	// the MODEL_STAGE group and every result/source path resolved from it.

	// release the prior stage's sources, sinks and cached element
	// Response* FIRST. They wrap Element*/Node* of the previous domain state; if
	// any writer below fails and returns early, they must already be gone, or
	// the next record() (same stamp -> no rebuild) would call getResponse() on
	// deleted elements.
	if (clearSources() != 0)
		return -1;
	m_data->nodes.clear();

	// MODEL_STAGE[<stamp>] + MODEL + RESULTS skeleton
	std::stringstream ss_stage;
	ss_stage << "MODEL_STAGE[" << info.current_model_stage_id << "]";
	hid_t h_stage = ladrunons::h5::group::create(
		info.h_file_id, ss_stage.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::attribute::write(h_stage, "STEP", info.current_time_step_id);
	ladrunons::h5::attribute::write(h_stage, "TIME", info.current_time_step);
	// KIND: a modal result (ModesOfVibration[Rotational])
	// promotes it to "eigen"; else the -kind override; else the analysis running
	// at this, the stage's first record(). With no analysis visible yet it is a
	// provisional "static" that record() rewrites once one runs (kind_pending).
	// stage_run_kind (what the flip test compares) is the running analysis
	// whatever the label; "" until a recorded step sees one.
	const std::string run_kind = ladrunoRunningAnalysisKind();
	std::string stage_kind;
	m_data->kind_pending = false;
	m_data->stage_run_kind = run_kind;
	if (ladrunoHasModalRequest(m_data->nodal_results_requests)) {
		stage_kind = "eigen";
	}
	else if (!m_data->stage_kind.empty()) {
		stage_kind = m_data->stage_kind;
	}
	else if (!run_kind.empty()) {
		stage_kind = run_kind;
	}
	else {
		stage_kind = "static";
		m_data->kind_pending = true;
	}
	ladrunons::h5::attribute::write(h_stage, "KIND", stage_kind);

	hid_t h_model = ladrunons::h5::group::create(
		h_stage, "MODEL", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::group::close(h_model);

	hid_t h_results = ladrunons::h5::group::create(
		h_stage, "RESULTS", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	hid_t h_on_nodes = ladrunons::h5::group::create(
		h_results, "ON_NODES", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::group::close(h_on_nodes);
	hid_t h_on_elems = ladrunons::h5::group::create(
		h_results, "ON_ELEMENTS", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::group::close(h_on_elems);
	// ON_DOMAIN / ON_REGIONS skeleton (energy balance)
	hid_t h_on_domain = ladrunons::h5::group::create(
		h_results, "ON_DOMAIN", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::group::close(h_on_domain);
	hid_t h_on_regions = ladrunons::h5::group::create(
		h_results, "ON_REGIONS", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::group::close(h_on_regions);
	ladrunons::h5::group::close(h_results);

	ladrunons::h5::group::close(h_stage);

	if (writeModelNodes() != 0)
		return -1;
	if (writeModelElements() != 0)
		return -1;
	if (writeModelSets() != 0)
		return -1;
	if (writeModelLocalAxes() != 0)
		return -1;
	if (writeSections() != 0)
		return -1;

	// (Re)build the source/sink channels for this stage. The prior stage's
	// sources were released at the top of writeModel(); the fresh sinks are not
	// yet initialized, so they re-create their result groups under the new
	// MODEL_STAGE.
	if (initNodeSources() != 0)
		return -1;
	if (initElementSources() != 0)
		return -1;
	if (initDomainSources() != 0)
		return -1;

	return 0;
}

int LadrunoRecorder::writeModelNodes()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	const int ndim = info.num_dimensions;

	// gather the node set (skip pressure nodes)
	m_data->nodes.clear();
	ID pressure_node_tags(0, info.domain->getNumPCs());
	Pressure_ConstraintIter& pc_iter = info.domain->getPCs();
	Pressure_Constraint* pc = 0;
	while ((pc = pc_iter()) != 0) {
		Node* pn = pc->getPressureNode();
		if (pn)
			pressure_node_tags.insert(pn->getTag());
	}
	if (m_data->has_region) {
		for (size_t i = 0; i < m_data->node_set.size(); ++i) {
			int nd = m_data->node_set[i];
			if (pressure_node_tags.getLocationOrdered(nd) < 0) {
				Node* cn = info.domain->getNode(nd);
				if (cn)
					m_data->nodes.push_back(cn);
			}
		}
	}
	else {
		NodeIter& nit = info.domain->getNodes();
		Node* nptr = 0;
		while ((nptr = nit()) != 0) {
			int nd = nptr->getTag();
			if (pressure_node_tags.getLocationOrdered(nd) < 0)
				m_data->nodes.push_back(nptr);
		}
	}

	const size_t nnodes = m_data->nodes.size();
	if (nnodes == 0) {
		// no nodes on THIS process is a valid, empty partition: openseesmp
		// with a `-R` region that lives entirely on other ranks, or a staged
		// run that removed the whole region. A part file with no MODEL/NODES
		// would break the stitch, so write zero-length NODES/ID and
		// NODES/COORDINATES, mark the stage EMPTY_PARTITION=1, and carry on (the
		// node/element channels have nothing to record; ON_DOMAIN still works).
		std::stringstream ss_empty;
		ss_empty << "MODEL_STAGE[" << info.current_model_stage_id << "]";
		hid_t h_stage = H5Gopen2(info.h_file_id, ss_empty.str().c_str(), H5P_DEFAULT);
		if (h_stage >= 0) {
			ladrunons::h5::attribute::write(h_stage, "EMPTY_PARTITION", (int)1);
			ladrunons::h5::group::close(h_stage);
		}
		ss_empty << "/MODEL/NODES";
		hid_t h_nodes = ladrunons::h5::group::create(
			info.h_file_id, ss_empty.str().c_str(), H5P_DEFAULT,
			info.h_group_proplist, H5P_DEFAULT);
		if (h_nodes < 0)
			return -1;
		hsize_t d_ids[2] = { 0, 1 };
		hsize_t d_crd[2] = { 0, (hsize_t)ndim };
		hid_t s_ids = H5Screate_simple(2, d_ids, NULL);
		hid_t s_crd = H5Screate_simple(2, d_crd, NULL);
		hid_t ds_ids = H5Dcreate(h_nodes, "ID", H5T_STD_I32LE, s_ids, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t ds_crd = H5Dcreate(h_nodes, "COORDINATES", H5T_IEEE_F64LE, s_crd, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		const bool ok = (ds_ids >= 0 && ds_crd >= 0);
		if (ds_crd >= 0) H5Dclose(ds_crd);
		if (ds_ids >= 0) H5Dclose(ds_ids);
		H5Sclose(s_crd);
		H5Sclose(s_ids);
		ladrunons::h5::group::close(h_nodes);
		if (!ok) {
			opserr << "LadrunoRecorder Error: could not write an empty MODEL/NODES\n";
			return -1;
		}
		opserr << "LadrunoRecorder: MODEL_STAGE[" << info.current_model_stage_id
		       << "] has no nodes on this process; written as an empty partition\n";
		return 0;
	}

	std::vector<int> node_ids(nnodes);
	std::vector<double> node_coords(nnodes * (size_t)ndim);
	for (size_t i = 0; i < nnodes; ++i) {
		Node* inode = m_data->nodes[i];
		node_ids[i] = inode->getTag();
		const Vector& crds = inode->getCrds();
		size_t j = i * (size_t)ndim;
		const int nc = crds.Size();  // a node may carry fewer coords than ndim -> clamp
		node_coords[j] = (nc > 0) ? crds[0] : 0.0;
		if (ndim > 1) {
			node_coords[j + 1] = (nc > 1) ? crds[1] : 0.0;
			if (ndim > 2)
				node_coords[j + 2] = (nc > 2) ? crds[2] : 0.0;
		}
	}

	std::stringstream ss_nodes;
	ss_nodes << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/NODES";
	hid_t h_nodes = ladrunons::h5::group::create(
		info.h_file_id, ss_nodes.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);
	hid_t d_id = ladrunons::h5::dataset::createAndWrite(h_nodes, "ID", node_ids);
	ladrunons::h5::dataset::close(d_id);
	hid_t d_coord = ladrunons::h5::dataset::createAndWrite(
		h_nodes, "COORDINATES", node_coords, nnodes, (size_t)ndim);
	ladrunons::h5::dataset::close(d_coord);
	ladrunons::h5::group::close(h_nodes);

	return 0;
}

int LadrunoRecorder::writeModelElements()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;

	std::stringstream ss_elems;
	ss_elems << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/ELEMENTS";
	hid_t h_gp_elements = ladrunons::h5::group::create(
		info.h_file_id, ss_elems.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);

	// map all elements (also populates registered custom rules)
	m_data->elements.mapElements(info.domain, m_data->has_region, m_data->elem_set);

	namespace me = ladrunons::detail::element;
	for (me::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
		it1 != m_data->elements.items.end(); ++it1) {
		me::ElementWithSameClassTagCollection& elem_by_tag = it1->second;
		for (me::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
			it2 != elem_by_tag.items.end(); ++it2) {
			me::ElementWithSameIntRuleCollection& elem_by_rule = it2->second;
			for (me::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
				it3 != elem_by_rule.items.end(); ++it3) {
				me::ElementWithSameCustomIntRuleCollection& elem_by_custom_rule = it3->second;
				if (elem_by_custom_rule.items.empty())
					continue;
				// name <classTag>-<className>[<int_rule>:<custom_rule_index>]
				std::stringstream ss_dset;
				ss_dset << elem_by_tag.class_tag << "-" << elem_by_tag.class_name
					<< "[" << elem_by_rule.int_rule_type << ":"
					<< elem_by_custom_rule.custom_int_rule_index << "]";
				elem_by_custom_rule.name = ss_dset.str();

				// CONNECTIVITY [nElem x (1 + num_nodes)]: (elemTag, c1..cK)
				std::vector<int> buffer(elem_by_custom_rule.items.size() * (size_t)(1 + elem_by_tag.num_nodes));
				size_t offset = 0;
				for (me::ElementWithSameCustomIntRuleCollection::collection_type::iterator it4 = elem_by_custom_rule.items.begin();
					it4 != elem_by_custom_rule.items.end(); ++it4) {
					Element* elem = *it4;
					buffer[offset] = elem->getTag();
					for (int j = 0; j < elem_by_tag.num_nodes; ++j)
						buffer[(size_t)j + 1 + offset] = elem->getExternalNodes()(j);
					offset += (size_t)(1 + elem_by_tag.num_nodes);
				}
				// Element group <name> is a GROUP holding CONNECTIVITY + BASIS attrs
				// + a QUADRATURE child (a dataset cannot parent the QUADRATURE group).
				hid_t h_eg = ladrunons::h5::group::create(
					h_gp_elements, elem_by_custom_rule.name.c_str(), H5P_DEFAULT,
					info.h_group_proplist, H5P_DEFAULT);
				hid_t d_conn = ladrunons::h5::dataset::createAndWrite(
					h_eg, "CONNECTIVITY", buffer,
					elem_by_custom_rule.items.size(), (size_t)(1 + elem_by_tag.num_nodes));
				ladrunons::h5::dataset::close(d_conn);

				// BASIS / QUADRATURE descriptors DERIVED from the legacy
				// (geometry, integration rule) pair. The legacy
				// GEOMETRY/INTEGRATION_RULE attrs are kept (for tooling that reads them)
				// AND add the derived self-describing descriptors.
				ladrunons::detail::ElementGeometryType::Enum geom = elem_by_tag.geom_type;
				ladrunons::detail::ElementIntegrationRuleType::Enum irule = elem_by_rule.int_rule_type;
				ladrunons::h5::attribute::write(h_eg, "GEOMETRY", (int)geom);
				ladrunons::h5::attribute::write(h_eg, "INTEGRATION_RULE", (int)irule);

				// derive TOPOLOGY / FAMILY / ORDER / PARAM_DOMAIN / NUM_CTRL
				std::string topology = "custom";
				std::string param_domain = "[-1,1]";
				std::vector<int> order;
				switch (geom) {
				case ladrunons::detail::ElementGeometryType::Line_2N:
				case ladrunons::detail::ElementGeometryType::Line_3N:
					topology = "line"; param_domain = "[-1,1]"; order.push_back(1); break;
				case ladrunons::detail::ElementGeometryType::Triangle_3N:
				case ladrunons::detail::ElementGeometryType::Triangle_6N:
					topology = "tri"; param_domain = "bary"; order.push_back(1); break;
				case ladrunons::detail::ElementGeometryType::Quadrilateral_4N:
				case ladrunons::detail::ElementGeometryType::Quadrilateral_8N:
				case ladrunons::detail::ElementGeometryType::Quadrilateral_9N:
				case ladrunons::detail::ElementGeometryType::Quadrilateral_CohesiveBand_4N:
					topology = "quad"; param_domain = "[-1,1]"; order.push_back(1); order.push_back(1); break;
				case ladrunons::detail::ElementGeometryType::Tetrahedron_4N:
				case ladrunons::detail::ElementGeometryType::Tetrahedron_10N:
					topology = "tet"; param_domain = "bary"; order.push_back(1); break;
				case ladrunons::detail::ElementGeometryType::Hexahedron_8N:
				case ladrunons::detail::ElementGeometryType::Hexahedron_20N:
				case ladrunons::detail::ElementGeometryType::Hexahedron_27N:
					topology = "hex"; param_domain = "[-1,1]"; order.push_back(1); order.push_back(1); order.push_back(1); break;
				default:
					topology = "custom"; break;
				}
				// A self-describing element (one that answers the "basisInfo"
				// response) overrides the geometry-derived guesses above. Lagrange
				// elements don't answer (basis_info.valid == false) so the derived
				// defaults are kept; a quadratic Bernstein element instead records
				// FAMILY=bernstein, ORDER=2.
				const me::BasisInfo& binfo = elem_by_tag.basis_info;
				std::string family = "lagrange";
				int rational = 0;
				int num_ctrl = elem_by_tag.num_nodes;
				if (binfo.valid) {
					if (!binfo.family.empty())       family = binfo.family;
					if (!binfo.topology.empty())     topology = binfo.topology;
					if (!binfo.param_domain.empty()) param_domain = binfo.param_domain;
					if (binfo.rational >= 0)         rational = binfo.rational;
					if (binfo.num_ctrl > 0)          num_ctrl = binfo.num_ctrl;
					if (binfo.order > 0) {
						order.clear();
						// ORDER carries one entry per GP_PARAM column (the validator equates
						// len(ORDER) with the parametric dimension). For a multi-dim custom
						// rule (barycentric triangle: 2 free area coords) replicate the
						// declared total degree across the parametric directions.
						int ndir = 1;
						if (irule == ladrunons::detail::ElementIntegrationRuleType::CustomIntegrationRule
							&& elem_by_custom_rule.custom_int_rule_index != 0) {
							int d = m_data->elements.registered_custom_rules[
								elem_by_custom_rule.custom_int_rule_index].custom_rule_dimension;
							if (d > 1) ndir = d;
						}
						for (int k = 0; k < ndir; ++k) order.push_back(binfo.order);
					}
				}
				ladrunons::h5::attribute::write(h_eg, "TOPOLOGY", topology);
				ladrunons::h5::attribute::write(h_eg, "FAMILY", family);
				if (order.size() > 0)
					ladrunons::h5::attribute::write(h_eg, "ORDER", order);
				ladrunons::h5::attribute::write(h_eg, "PARAM_DOMAIN", param_domain);
				ladrunons::h5::attribute::write(h_eg, "RATIONAL", (int)rational);
				ladrunons::h5::attribute::write(h_eg, "NUM_CTRL", (int)num_ctrl);

				// ---- Resolve the quadrature rule: custom -> standard table -> none ----
				// QUADRATURE/{GP_PARAM[NUM_GP x NDIR], GP_WEIGHT[NUM_GP]} + NDIR + NUM_GP
				// are written for BOTH custom-rule elements AND standard Gauss-Legendre
				// elements (quad/hex/tri/tet/line). For standard rules the GP coords +
				// weights come from the built-in getStandardQuadrature table, written
				// file-side so a reader needs no per-element quadrature tables. The
				// bucket <name> is a GROUP (see above), so the QUADRATURE child group is
				// valid.
				std::vector<double> q_gp_param, q_gp_weight;
				int q_num_gp = 0, q_ndir = 0;
				bool have_rule = false;
				if (irule == ladrunons::detail::ElementIntegrationRuleType::CustomIntegrationRule) {
					ladrunons::h5::attribute::write(h_eg, "CUSTOM_INTEGRATION_RULE", elem_by_custom_rule.custom_int_rule_index);
					if (elem_by_custom_rule.custom_int_rule_index != 0) {
						me::ElementIntegrationRule& custom_rule =
							m_data->elements.registered_custom_rules[elem_by_custom_rule.custom_int_rule_index];
						// legacy GP_X attribute (flat coords) kept for tooling that reads it.
						ladrunons::h5::attribute::write(h_eg, "GP_X", custom_rule.x);
						ladrunons::h5::attribute::write(h_eg, "CUSTOM_INTEGRATION_RULE_DIMENSION", custom_rule.custom_rule_dimension);
						// reshape: legacy line rules store 1 coord/point (num_gp==0 => cols=1);
						// multi-dim parametric rules (e.g. barycentric triangles) store num_gp
						// rows of (x.size()/num_gp) coords each.
						q_num_gp = custom_rule.num_gp > 0 ? custom_rule.num_gp : (int)custom_rule.x.size();
						q_ndir = (custom_rule.num_gp > 0 && custom_rule.num_gp <= (int)custom_rule.x.size())
							? (int)custom_rule.x.size() / custom_rule.num_gp : 1;
						q_gp_param = custom_rule.x;
						q_gp_weight = custom_rule.w;
						have_rule = (q_num_gp > 0 && q_ndir > 0);
					}
				}
				else {
					// standard Gauss rule (returns false for untabulated rules ->
					// no QUADRATURE written; the reader tolerates its absence).
					have_rule = me::getStandardQuadrature(irule, q_gp_param, q_gp_weight, q_num_gp, q_ndir);
				}

				if (have_rule) {
					ladrunons::h5::attribute::write(h_eg, "NDIR", q_ndir);
					ladrunons::h5::attribute::write(h_eg, "NUM_GP", q_num_gp);
					// QUADRATURE child group: GP_PARAM (NUM_GP x NDIR) + optional GP_WEIGHT.
					hid_t h_q = ladrunons::h5::group::create(
						h_eg, "QUADRATURE", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
					hid_t d_gpp = ladrunons::h5::dataset::createAndWrite(
						h_q, "GP_PARAM", q_gp_param, (size_t)q_num_gp, (size_t)q_ndir);
					ladrunons::h5::dataset::close(d_gpp);
					if (!q_gp_weight.empty()) {
						hid_t d_gpw = ladrunons::h5::dataset::createAndWrite(
							h_q, "GP_WEIGHT", q_gp_weight);  // 1-D [nGP]
						ladrunons::h5::dataset::close(d_gpw);
					}
					ladrunons::h5::group::close(h_q);

					// ---- GLOBAL_GP_COORDS (write-side cross-check of GP_PARAM) -----
					// Static, reference-config, C++-computed x(ξ_g)=ΣN_i(ξ_g)X_i per
					// element via the write-side basis evaluator. Stored 2-D
					// [nElem x (NUM_GP*ndim)] because Ladruno_Hdf5 has no rank-3 writer;
					// reader reshapes to [nElem, NUM_GP, ndim]. Skipped (graceful) when
					// computeGlobalGP doesn't yet have the topology's basis. Reuses the
					// CONNECTIVITY node access pattern (domain->getNode->getCrds).
					int ndim = 0;
					{
						const ID& nd0 = elem_by_custom_rule.items.front()->getExternalNodes();
						Node* n0 = (nd0.Size() > 0) ? info.domain->getNode(nd0(0)) : 0;
						if (n0) ndim = n0->getCrds().Size();
					}
					if (ndim > 0) {
						size_t nElem = elem_by_custom_rule.items.size();
						std::vector<double> ggp;
						ggp.reserve(nElem * (size_t)q_num_gp * (size_t)ndim);
						bool ok = true;
						for (me::ElementWithSameCustomIntRuleCollection::collection_type::iterator it5 = elem_by_custom_rule.items.begin();
							it5 != elem_by_custom_rule.items.end() && ok; ++it5) {
							Element* elem = *it5;
							const ID& nd = elem->getExternalNodes();
							int nn = elem_by_tag.num_nodes;
							std::vector<std::vector<double> > node_xyz(
								(size_t)nn, std::vector<double>((size_t)ndim, 0.0));
							for (int j = 0; j < nn && ok; ++j) {
								Node* nj = (j < nd.Size()) ? info.domain->getNode(nd(j)) : 0;
								if (nj == 0) { ok = false; break; }
								const Vector& c = nj->getCrds();
								for (int d = 0; d < ndim && d < c.Size(); ++d)
									node_xyz[(size_t)j][(size_t)d] = c(d);
							}
							if (!ok) break;
							std::vector<double> ex;
							if (!me::computeGlobalGP(geom, nn, node_xyz, ndim,
								q_gp_param, q_num_gp, q_ndir, ex)) { ok = false; break; }
							ggp.insert(ggp.end(), ex.begin(), ex.end());
						}
						if (ok && !ggp.empty()) {
							hid_t d_g = ladrunons::h5::dataset::createAndWrite(
								h_eg, "GLOBAL_GP_COORDS", ggp,
								nElem, (size_t)q_num_gp * (size_t)ndim);
							ladrunons::h5::dataset::close(d_g);
						}
					}
				}
				ladrunons::h5::group::close(h_eg);
			}
		}
	}

	ladrunons::h5::group::close(h_gp_elements);
	return 0;
}

int LadrunoRecorder::writeModelSets()
{
	// MODEL/SETS: one SET_<tag> per region whose identity the file must be
	// self-describing about: the per-region energy tags plus the
	// -R restriction region (if any). Each SET_<tag> carries NODES + ELEMENTS
	// tag datasets and a TAG attr. If no region is referenced, the (empty) SETS
	// group is still created for a consistent layout.
	ladrunons::detail::ProcessInfo& info = m_data->info;

	// union of region tags to describe (preserve order, drop duplicates)
	std::vector<int> set_tags;
	{
		std::set<int> seen;
		if (m_data->has_region && m_data->region_tag != 0) {
			if (seen.insert(m_data->region_tag).second)
				set_tags.push_back(m_data->region_tag);
		}
		for (size_t i = 0; i < m_data->energy_region_tags.size(); ++i) {
			int t = m_data->energy_region_tags[i];
			if (seen.insert(t).second)
				set_tags.push_back(t);
		}
	}

	// Nothing to describe -> don't create an empty MODEL/SETS group (it adds
	// benign HDF5 open-probe noise on region-free models for no benefit).
	if (set_tags.empty())
		return 0;

	std::stringstream ss_sets;
	ss_sets << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/SETS";
	hid_t h_sets = ladrunons::h5::group::create(
		info.h_file_id, ss_sets.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);

	for (size_t i = 0; i < set_tags.size(); ++i) {
		int tag = set_tags[i];
		MeshRegion* reg = info.domain->getRegion(tag);
		if (reg == 0) {
			opserr << "LadrunoRecorder - WARNING region " << tag
			       << " referenced by the recorder was not found; its MODEL/SETS "
			          "entry is skipped\n";
			continue;
		}
		std::stringstream ss;
		ss << "MODEL_STAGE[" << info.current_model_stage_id
		   << "]/MODEL/SETS/SET_" << tag;
		hid_t h = ladrunons::h5::group::create(
			info.h_file_id, ss.str().c_str(), H5P_DEFAULT,
			info.h_group_proplist, H5P_DEFAULT);
		ladrunons::h5::attribute::write(h, "TAG", tag);

		const ID& rnodes = reg->getNodes();
		std::vector<int> node_ids((size_t)rnodes.Size());
		for (int k = 0; k < rnodes.Size(); ++k)
			node_ids[(size_t)k] = rnodes(k);
		const ID& relems = reg->getElements();
		std::vector<int> elem_ids((size_t)relems.Size());
		for (int k = 0; k < relems.Size(); ++k)
			elem_ids[(size_t)k] = relems(k);

		hid_t d_n = ladrunons::h5::dataset::createAndWrite(h, "NODES", node_ids);
		if (d_n != ladrunons::HID_INVALID)
			ladrunons::h5::dataset::close(d_n);
		hid_t d_e = ladrunons::h5::dataset::createAndWrite(h, "ELEMENTS", elem_ids);
		if (d_e != ladrunons::HID_INVALID)
			ladrunons::h5::dataset::close(d_e);

		ladrunons::h5::group::close(h);
	}

	ladrunons::h5::group::close(h_sets);
	return 0;
}

int LadrunoRecorder::writeModelLocalAxes()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	namespace me = ladrunons::detail::element;

	// Per-class MODEL/LOCAL_AXES/<classTag>-<ClassName>/{ID[E], FRAME[E×4]}.
	// The frame is sourced from the element's "localAxes" response (9 packed direction
	// cosines vx|vy|vz) and converted to a quaternion. Elements that don't answer are
	// skipped: a silent identity frame is NEVER written.
	const char* request[1] = { "localAxes" };
	hid_t h_la_root = ladrunons::HID_INVALID; // created lazily on the first answering class group

	for (me::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
		it1 != m_data->elements.items.end(); ++it1) {
		me::ElementWithSameClassTagCollection& elem_by_tag = it1->second;

		std::vector<int> ids;
		std::vector<ladrunons::utils::locax::quaternion> frames;
		Vector vx(3), vy(3), vz(3);

		for (me::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
			it2 != elem_by_tag.items.end(); ++it2) {
			me::ElementWithSameIntRuleCollection& elem_by_rule = it2->second;
			for (me::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
				it3 != elem_by_rule.items.end(); ++it3) {
				me::ElementWithSameCustomIntRuleCollection& ebcr = it3->second;
				for (std::vector<Element*>::iterator it4 = ebcr.items.begin();
					it4 != ebcr.items.end(); ++it4) {
					Element* elem = *it4;
					me::OutputDescriptor desc;
					me::OutputDescriptorStream stream(&desc);
					Response* resp = elem->setResponse(request, 1, stream);
					if (resp == 0)
						continue;
					resp->getResponse();
					const Vector& packed = resp->getInformation().getData();
					if (packed.Size() == 9) {
						for (int i = 0; i < 3; ++i) {
							vx(i) = packed(i); vy(i) = packed(i + 3); vz(i) = packed(i + 6);
						}
						ids.push_back(elem->getTag());
						frames.push_back(ladrunons::utils::locax::quatFromMat(vx, vy, vz));
					}
					delete resp;
				}
			}
		}

		if (ids.empty())
			continue;

		if (h_la_root == ladrunons::HID_INVALID) {
			std::stringstream ss_root;
			ss_root << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/LOCAL_AXES";
			h_la_root = ladrunons::h5::group::create(info.h_file_id, ss_root.str().c_str(),
				H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
		}
		std::stringstream ss_name;
		ss_name << elem_by_tag.class_tag << "-" << elem_by_tag.class_name;
		hid_t h_grp = ladrunons::h5::group::create(h_la_root, ss_name.str().c_str(),
			H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
		hid_t d_id = ladrunons::h5::dataset::createAndWrite(h_grp, "ID", ids);
		ladrunons::h5::dataset::close(d_id);
		hid_t d_fr = ladrunons::h5::dataset::createAndWrite(h_grp, "FRAME", frames); // [E × 4]
		ladrunons::h5::dataset::close(d_fr);
		ladrunons::h5::group::close(h_grp);
	}

	if (h_la_root != ladrunons::HID_INVALID)
		ladrunons::h5::group::close(h_la_root);
	return 0;
}

int LadrunoRecorder::writeSections()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	namespace me = ladrunons::detail::element;

	std::map<int, me::SectionAssignment> sec_assignments;
	std::map<me::FiberSectionData, me::SectionAssignment> aux_sec_assignments;

	m_data->elem_ngauss_nfiber_info.clear();

	std::stringstream ss_sec_dir;
	ss_sec_dir << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/SECTION_ASSIGNMENTS";
	hid_t h_gp_sec = ladrunons::h5::group::create(
		info.h_file_id, ss_sec_dir.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);

	// loop all element groups (section/fiber discovery, as in MPCORecorder)
	for (me::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
		it1 != m_data->elements.items.end(); ++it1) {
		me::ElementWithSameClassTagCollection& elem_by_tag = it1->second;
		for (me::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
			it2 != elem_by_tag.items.end(); ++it2) {
			me::ElementWithSameIntRuleCollection& elem_by_rule = it2->second;
			for (me::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
				it3 != elem_by_rule.items.end(); ++it3) {
				me::ElementWithSameCustomIntRuleCollection& elem_by_custom_rule = it3->second;
				for (std::vector<Element*>::iterator it4 = elem_by_custom_rule.items.begin();
					it4 != elem_by_custom_rule.items.end(); ++it4) {
					Element* elem = *it4;

					std::string request1 = ladrunons::utils::shell::isShellElementTag(elem->getClassTag()) ? "material" : "section";
					std::string request2 = "fiber";
					std::string request3 = "stress";
					int argc = 5;
					const char** argv = new const char*[argc];
					argv[0] = request1.c_str();
					argv[2] = request2.c_str();
					argv[4] = request3.c_str();
					int woagg_agrc = 6;
					const char** woagg_argv = new const char*[woagg_agrc];
					woagg_argv[0] = request1.c_str();
					woagg_argv[2] = request1.c_str();
					woagg_argv[3] = request2.c_str();
					woagg_argv[5] = request3.c_str();

					std::vector<int> elem_gauss_id;
					std::vector<int> elem_sec_id;
					std::vector<bool> elem_dummy_sec_flags;
					std::vector<me::FiberSectionData> elem_sections;
					std::vector<int> elem_fiber_base_index;

					int trial_num_sec = 0;
					while (true) {
						trial_num_sec++;
						if (trial_num_sec > LADRUNO_MAX_TRIAL_NSEC)
							break;
						std::stringstream ss_ts; ss_ts << trial_num_sec;
						std::string s_ts = ss_ts.str();
						argv[1] = s_ts.c_str();
						woagg_argv[1] = s_ts.c_str();

						int trial_num_fib = -1;
						bool break_sec_loop = false;
						bool first_fiber_done = false;
						bool do_workaround_for_aggregator = false;
						while (true) {
							trial_num_fib++;
							if (trial_num_fib > LADRUNO_MAX_TRIAL_NFIB)
								break;
							std::stringstream ss_tf; ss_tf << trial_num_fib;
							std::string s_tf = ss_tf.str();
							argv[3] = s_tf.c_str();
							woagg_argv[4] = s_tf.c_str();

							me::OutputDescriptor eo_descriptor;
							me::OutputDescriptorStream eo_stream(&eo_descriptor);
							Response* eo_response = do_workaround_for_aggregator ?
								elem->setResponse(woagg_argv, woagg_agrc, eo_stream) :
								elem->setResponse(argv, argc, eo_stream);
							eo_stream.finalizeSetResponse();
							if (eo_response)
								delete eo_response;

							std::vector<me::FiberData> trial_fiberdata;
							std::vector<int> trial_fiber_material_id;
							std::vector<int> trial_sec_id;
							std::vector<int> trial_gp_id;
							std::vector<bool> trial_dummy_flag;
							eo_descriptor.getFiberData(trial_fiberdata, trial_fiber_material_id,
								trial_sec_id, trial_gp_id, trial_dummy_flag);

							if (trial_gp_id.size() == 0) {
								break_sec_loop = true;
								break;
							}
							if (trial_sec_id.size() == 0)
								break;

							if (trial_sec_id.size() != trial_fiberdata.size()) trial_sec_id.resize(trial_fiberdata.size());
							if (trial_gp_id.size() != trial_fiberdata.size()) trial_gp_id.resize(trial_fiberdata.size());
							if (trial_dummy_flag.size() != trial_fiberdata.size()) trial_dummy_flag.resize(trial_fiberdata.size());

							if (trial_fiberdata.size() == 0 && trial_num_fib == 0)
								continue;

							if (!first_fiber_done) {
								elem_gauss_id.push_back(trial_gp_id[0]);
								elem_sec_id.push_back(trial_sec_id[0]);
								elem_dummy_sec_flags.push_back(trial_dummy_flag[0]);
								elem_sections.push_back(me::FiberSectionData());
								elem_fiber_base_index.push_back(trial_num_fib);
								first_fiber_done = true;
							}

							if (trial_fiberdata.size() == 0) {
								if (!do_workaround_for_aggregator && (trial_num_fib == 1)) {
									do_workaround_for_aggregator = true;
									if (first_fiber_done) {
										elem_gauss_id.pop_back();
										elem_sec_id.pop_back();
										elem_dummy_sec_flags.pop_back();
										elem_sections.pop_back();
										elem_fiber_base_index.pop_back();
									}
									trial_num_fib = -1;
									break_sec_loop = false;
									first_fiber_done = false;
									continue;
								}
								break;
							}

							if (trial_fiberdata.size() == 1) {
								me::FiberSectionData& curr = elem_sections.back();
								curr.fibers.push_back(trial_fiberdata[0]);
								curr.materials.push_back(trial_fiber_material_id[0]);
							}
						} // fiber loop
						if (break_sec_loop)
							break;
					} // section loop

					delete[] argv;
					delete[] woagg_argv;

					// store ngauss/nfiber info for the element-result section path
					std::vector<std::pair<int, int> >& nfpg = m_data->elem_ngauss_nfiber_info[elem->getTag()];
					nfpg.resize(elem_gauss_id.size());
					for (size_t igp = 0; igp < elem_gauss_id.size(); ++igp) {
						nfpg[igp].first = elem_fiber_base_index[igp];
						nfpg[igp].second = (int)elem_sections[igp].fibers.size();
					}

					size_t ngauss = elem_gauss_id.size();
					for (size_t igp = 0; igp < ngauss; ++igp) {
						int igp_sec_id = elem_sec_id[igp];
						bool igp_dummy = elem_dummy_sec_flags[igp];
						me::FiberSectionData& igp_data = elem_sections[igp];
						if (igp_dummy) {
							me::SectionAssignment& ca = aux_sec_assignments[igp_data];
							if (ca.is_new) { ca.fiber_section_data = igp_data; ca.is_new = false; }
							ca.assignments.push_back(me::ElemGaussPair(elem->getTag(), (int)igp));
						}
						else {
							me::SectionAssignment& ca = sec_assignments[igp_sec_id];
							if (ca.is_new) { ca.fiber_section_data = igp_data; ca.is_new = false; }
							ca.assignments.push_back(me::ElemGaussPair(elem->getTag(), (int)igp));
						}
					}
				} // each element
			}
		}
	}

	// resolve section class names
	for (std::map<int, me::SectionAssignment>::iterator it = sec_assignments.begin();
		it != sec_assignments.end(); ++it) {
		SectionForceDeformation* sfd = OPS_getSectionForceDeformation(it->first);
		if (sfd)
			it->second.name = sfd->getClassType();
	}

	// move dummy (floating-fiber) sections to tags above max
	int max_sec_tag = sec_assignments.size() > 0 ? (--sec_assignments.end())->first : 0;
	for (std::map<me::FiberSectionData, me::SectionAssignment>::iterator it = aux_sec_assignments.begin();
		it != aux_sec_assignments.end(); ++it) {
		sec_assignments[++max_sec_tag] = it->second;
	}

	// write each assignment (KIND=fiber)
	for (std::map<int, me::SectionAssignment>::iterator it = sec_assignments.begin();
		it != sec_assignments.end(); ++it) {
		me::SectionAssignment& sec_asn = it->second;
		if (sec_asn.assignments.size() == 0)
			continue;
		std::stringstream ss_name;
		ss_name << "SECTION_" << it->first << "[" << sec_asn.name << "]";
		hid_t h_isec = ladrunons::h5::group::create(
			h_gp_sec, ss_name.str().c_str(), H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
		ladrunons::h5::attribute::write(h_isec, "ID", it->first);
		ladrunons::h5::attribute::write(h_isec, "NAME", sec_asn.name);
		// a section that exposes no fibers (elastic, membrane-plate, ...) is a
		// stress-resultant section: KIND=resultant, no FIBER_DATA
		const bool has_fibers = !sec_asn.fiber_section_data.fibers.empty();
		ladrunons::h5::attribute::write(h_isec, "KIND",
			std::string(has_fibers ? "fiber" : "resultant"));
		{
			hid_t dset = ladrunons::h5::dataset::createAndWrite(h_isec, "ASSIGNMENT", sec_asn.assignments);
			ladrunons::h5::dataset::close(dset);
		}
		if (has_fibers) {
			hid_t dset = ladrunons::h5::dataset::createAndWrite(h_isec, "FIBER_DATA", sec_asn.fiber_section_data.fibers);
			ladrunons::h5::dataset::close(dset);
			dset = ladrunons::h5::dataset::createAndWrite(h_isec, "FIBER_MATERIALS", sec_asn.fiber_section_data.materials);
			ladrunons::h5::dataset::close(dset);
		}
		ladrunons::h5::group::close(h_isec);
	}

	ladrunons::h5::group::close(h_gp_sec);
	return 0;
}

/* ===================================================================== */
/* node sources                                                          */
/* ===================================================================== */

namespace {
	// in a partitioned run each partition file holds only its own share of
	// an additive result (e.g. the reaction at a support on a partition interface),
	// so the reader must combine the files. A per-partition EXTREME of a partial can
	// never be recombined (max(a+b) != max a + max b), so -envelope is refused for any
	// source whose PARTITION_REDUCTION is not "NONE". Returns true = skip the channel.
	bool refusePartitionedEnvelope(bool envelope_mode, bool is_partitioned,
	                               const ladrunons::ResultSource* src)
	{
		if (!envelope_mode || !is_partitioned ||
		    std::strcmp(src->partitionReduction(), "NONE") == 0)
			return false;
		opserr << "LadrunoRecorder WARNING: -envelope of '" << src->schema().name.c_str()
		       << "' is NOT recorded in this partitioned run: each partition holds only a "
		          "partial (PARTITION_REDUCTION=" << src->partitionReduction()
		       << "), and a per-partition extreme cannot be recombined. Record it without "
		          "-envelope; the reader combines the partition files.\n";
		return true;
	}
}

int LadrunoRecorder::initNodeSources()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	namespace NR = ladrunons::detail; // for NodalResultType

	// node-tag set passed to each source (matches writeModelNodes order)
	std::vector<int> node_ids(m_data->nodes.size());
	for (size_t i = 0; i < m_data->nodes.size(); ++i)
		node_ids[i] = m_data->nodes[i]->getTag();

	for (size_t i = 0; i < m_data->nodal_results_requests.size(); ++i) {
		NR::NodalResultType::Enum rtype = m_data->nodal_results_requests[i];
		int grad = m_data->sens_grad_indices[i];

		private_data::NodeChannel ch;
		ladrunons::ResultSource* src = 0;
		int reac = -1;
		bool is_modes = false;

		switch (rtype) {
		case NR::NodalResultType::Displacement:
			src = new ladrunons::DisplacementSource(info, node_ids); break;
		case NR::NodalResultType::Rotation:
			src = new ladrunons::RotationSource(info, node_ids); break;
		case NR::NodalResultType::Velocity:
			src = new ladrunons::VelocitySource(info, node_ids); break;
		case NR::NodalResultType::AngularVelocity:
			src = new ladrunons::AngularVelocitySource(info, node_ids); break;
		case NR::NodalResultType::Acceleration:
			src = new ladrunons::AccelerationSource(info, node_ids); break;
		case NR::NodalResultType::AngularAcceleration:
			src = new ladrunons::AngularAccelerationSource(info, node_ids); break;
		case NR::NodalResultType::Pressure:
			src = new ladrunons::PressureSource(info, node_ids); break;
		case NR::NodalResultType::ReactionForce: {
			ladrunons::ReactionForceSource* s = new ladrunons::ReactionForceSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::ReactionMoment: {
			ladrunons::ReactionMomentSource* s = new ladrunons::ReactionMomentSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::ReactionForceIncludingInertia: {
			ladrunons::ReactionForceIncInertiaSource* s = new ladrunons::ReactionForceIncInertiaSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::ReactionMomentIncludingInertia: {
			ladrunons::ReactionMomentIncInertiaSource* s = new ladrunons::ReactionMomentIncInertiaSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::RayleighForce: {
			ladrunons::RayleighForceSource* s = new ladrunons::RayleighForceSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::RayleighMoment: {
			ladrunons::RayleighMomentSource* s = new ladrunons::RayleighMomentSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::UnbalancedForce:
			src = new ladrunons::UnbalancedForceSource(info, node_ids); break;
		case NR::NodalResultType::UnbalancedMoment:
			src = new ladrunons::UnbalancedMomentSource(info, node_ids); break;
		case NR::NodalResultType::UnbalancedForceIncludingInertia:
			src = new ladrunons::UnbalancedForceIncInertiaSource(info, node_ids); break;
		case NR::NodalResultType::UnbalancedMomentIncludingInertia:
			src = new ladrunons::UnbalancedMomentIncInertiaSource(info, node_ids); break;
		case NR::NodalResultType::ModesOfVibration:
			src = new ladrunons::ModesOfVibrationSource(info, node_ids); is_modes = true; break;
		case NR::NodalResultType::ModesOfVibrationRotational:
			src = new ladrunons::ModesOfVibrationRotationalSource(info, node_ids); is_modes = true; break;
		case NR::NodalResultType::DisplacementSensitivity:
			src = new ladrunons::DisplacementSensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::RotationSensitivity:
			src = new ladrunons::RotationSensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::VelocitySensitivity:
			src = new ladrunons::VelocitySensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::AngularVelocitySensitivity:
			src = new ladrunons::AngularVelocitySensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::AccelerationSensitivity:
			src = new ladrunons::AccelerationSensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::AngularAccelerationSensitivity:
			src = new ladrunons::AngularAccelerationSensitivitySource(info, node_ids, grad); break;
		default:
			break;
		}

		if (src == 0)
			continue;
		if (!is_modes &&
		    refusePartitionedEnvelope(m_data->envelope_mode, m_data->is_partitioned, src)) {
			delete src;   // a partitioned envelope of a partial is unrecoverable
			continue;
		}

		ch.source = src;
		// Envelope mode swaps the sink to an EnvelopeSink for ordinary results;
		// ModesOfVibration (eigen) is never enveloped; it keeps the StreamingSink
		// mode-loop path.
		if (m_data->envelope_mode && !is_modes)
			ch.sink = new ladrunons::EnvelopeSink(ladrunons::ResultFamily::OnNodes);
		else
			ch.sink = new ladrunons::StreamingSink(ladrunons::ResultFamily::OnNodes);
		ch.reaction_flag = reac;
		ch.is_modes = is_modes;
		m_data->node_channels.push_back(ch);
	}

	return 0;
}

/* ===================================================================== */
/* element sources                                                       */
/* ===================================================================== */

int LadrunoRecorder::initElementSources()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	namespace me = ladrunons::detail::element;

	if (m_data->elemental_results_requests.empty())
		return 0;

	// Build the per-request response-collection trees (port of initElementRecorders).
	m_data->elemental_recorders.resize(m_data->elemental_results_requests.size());

	for (size_t cnt_request = 0; cnt_request < m_data->elemental_results_requests.size(); ++cnt_request) {
		const std::vector<std::string>& request = m_data->elemental_results_requests[cnt_request];

		// the printable form of the -E request,
		// i.e. the dot-joined token list exactly as the user typed it. Used only by
		// the warnings below.
		std::string request_str;
		for (size_t i = 0; i < request.size(); ++i) {
			if (i > 0) request_str += ".";
			request_str += request[i];
		}
		// Number of elements, anywhere in the domain, that answered this request.
		// Zero => the request produced no bucket at all (see the warning after the
		// class loop). Only counts elements whose response AND descriptor stream
		// were both accepted, which is exactly the bucket-fill condition.
		size_t num_request_answers = 0;

		std::vector<std::string> request_mod;
		bool do_all_materials = false;
		bool do_all_sections = false;
		bool do_all_fibers = false;
		size_t material_id_placeholder_index = 0;
		size_t section_id_placeholder_index = 0;
		size_t fiber_id_placeholder_index = 0;
		for (size_t i = 0; i < request.size(); ++i) {
			const std::string& cr = request[i];
			request_mod.push_back(cr);
			if (i == 0) {
				if (cr == "section") {
					request_mod.push_back("");
					section_id_placeholder_index = request_mod.size() - 1;
					do_all_sections = true;
				}
				else if (cr == "material") {
					request_mod.push_back("");
					material_id_placeholder_index = request_mod.size() - 1;
					do_all_materials = true;
				}
			}
			else {
				// a "fiber" sub-verb expands per-layer under BOTH a
				// "section.fiber.*" and a "material.fiber.*" request. The latter is the
				// natural per-layer verb for shells (whose layers are addressed through
				// the "material" keyword); without it do_all_fibers is never set, the
				// fiber index is never substituted and the request emits no bucket.
				if (cr == "fiber" && (do_all_sections || do_all_materials)) {
					request_mod.push_back("");
					fiber_id_placeholder_index = request_mod.size() - 1;
					do_all_fibers = true;
				}
			}
		}
		int argc = (int)request_mod.size();
		const char** argv = new const char*[argc];
		for (size_t i = 0; i < request_mod.size(); ++i)
			argv[i] = request_mod[i].c_str();

		std::string aux_section_keyword_standard = "section";
		std::string aux_section_keyword_for_shells = "material";

		me::ResultRecorder& recorder = m_data->elemental_recorders[cnt_request];
		recorder.result_request = request;

		for (me::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
			it1 != m_data->elements.items.end(); ++it1) {
			me::ElementWithSameClassTagCollection& elem_by_tag = it1->second;
			me::OutputWithSameClassTagCollection& eo_by_tag = recorder.response_map[it1->first];

			// aggregate, per element class, how many elements refused this
			// request. NEVER report per element: a
			// domain-wide unknown token would otherwise emit one line per element.
			size_t num_class_elements = 0;
			size_t num_class_no_response = 0;     // setResponse returned null
			size_t num_class_stream_error = 0;    // response built, descriptor stream rejected it

			bool standard_section_keyword_modified = false;
			if (do_all_sections && ladrunons::utils::shell::isShellElementTag(it1->first)) {
				argv[0] = aux_section_keyword_for_shells.c_str();
				standard_section_keyword_modified = true;
			}

			for (me::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
				it2 != elem_by_tag.items.end(); ++it2) {
				me::ElementWithSameIntRuleCollection& elem_by_rule = it2->second;
				me::OutputWithSameIntRuleCollection& eo_by_rule = eo_by_tag.items[it2->first];
				for (me::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
					it3 != elem_by_rule.items.end(); ++it3) {
					me::ElementWithSameCustomIntRuleCollection& elem_by_custom_rule = it3->second;
					me::OutputWithSameCustomIntRuleCollection& eo_by_custom_rule = eo_by_rule.items[it3->first];
					int header_local_index = 0;
					for (me::ElementWithSameCustomIntRuleCollection::collection_type::iterator it4 = elem_by_custom_rule.items.begin();
						it4 != elem_by_custom_rule.items.end(); ++it4) {
						Element* elem = *it4;
						++num_class_elements;
						me::OutputDescriptor eo_descriptor;
						Response* eo_response = 0;
						me::OutputDescriptorStream eo_stream(&eo_descriptor);

						if (do_all_sections) {
							std::map<int, std::vector<std::pair<int, int> > >::const_iterator it_info =
								m_data->elem_ngauss_nfiber_info.find(elem->getTag());
							if (it_info != m_data->elem_ngauss_nfiber_info.end()) {
								CompositeResponse* sec_comp_response = new CompositeResponse();
								int num_sec_responses = 0;
								for (size_t section_id = 0; section_id < it_info->second.size(); ++section_id) {
									std::stringstream ss_sid; ss_sid << section_id + 1;
									std::string s_sid = ss_sid.str();
									argv[section_id_placeholder_index] = s_sid.c_str();
									Response* sec_response = 0;
									if (do_all_fibers) {
										CompositeResponse* fib_comp_response = new CompositeResponse();
										int num_fib_responses = 0;
										int fiber_base_id = it_info->second[section_id].first;
										int num_fibers = it_info->second[section_id].second;
										for (int fiber_id = 0; fiber_id < num_fibers; ++fiber_id) {
											std::string s_fid = ladrunons::utils::strings::to_string(fiber_id + fiber_base_id);
											argv[fiber_id_placeholder_index] = s_fid.c_str();
											bool was_valid_before = (eo_stream.error_code == me::OutputDescriptorStream::ERROR_CODE_OK);
											Response* fib_response = elem->setResponse(argv, argc, eo_stream);
											if (fib_response) {
												num_fib_responses = fib_comp_response->addResponse(fib_response);
											}
											else {
												if (was_valid_before && (eo_stream.error_code == me::OutputDescriptorStream::ERROR_CODE_SECTION_AFTER_FIBER)) {
													eo_descriptor.fixSectionAfterFiberDueToFiberOutputFail();
													eo_stream.error_code = me::OutputDescriptorStream::ERROR_CODE_OK;
												}
											}
										}
										if (num_fib_responses == 0)
											delete fib_comp_response;
										else
											sec_response = fib_comp_response;
									}
									// a bare "section" request (no sub-verb) would forward a zero-length
									// arg list to the section's setResponse, which derefs argv[0] -> segfault.
									// Skip it (emit no bucket) instead of crashing.
									else if (argc > (int)section_id_placeholder_index + 1) {
										sec_response = elem->setResponse(argv, argc, eo_stream);
									}
									if (sec_response)
										num_sec_responses = sec_comp_response->addResponse(sec_response);
								}
								if (num_sec_responses == 0)
									delete sec_comp_response;
								else
									eo_response = sec_comp_response;
							}
						}
						else if (do_all_materials) {
							CompositeResponse* mat_comp_response = new CompositeResponse();
							int num_mat_responses = 0;
							if (do_all_fibers) {
								// "material.fiber.<resp>" per-layer expansion, mirroring the
								// do_all_sections path. For shells the "material" id addresses the Gauss
								// point and the "fiber" id the through-thickness layer; the per-(gp,layer)
								// counts come from the same elem_ngauss_nfiber_info discovery table.
								std::map<int, std::vector<std::pair<int, int> > >::const_iterator it_info =
									m_data->elem_ngauss_nfiber_info.find(elem->getTag());
								if (it_info != m_data->elem_ngauss_nfiber_info.end()) {
									for (size_t mat_id = 0; mat_id < it_info->second.size(); ++mat_id) {
										std::stringstream ss_mid; ss_mid << mat_id + 1;
										std::string s_mid = ss_mid.str();
										argv[material_id_placeholder_index] = s_mid.c_str();
										CompositeResponse* fib_comp_response = new CompositeResponse();
										int num_fib_responses = 0;
										int fiber_base_id = it_info->second[mat_id].first;
										int num_fibers = it_info->second[mat_id].second;
										for (int fiber_id = 0; fiber_id < num_fibers; ++fiber_id) {
											std::string s_fid = ladrunons::utils::strings::to_string(fiber_id + fiber_base_id);
											argv[fiber_id_placeholder_index] = s_fid.c_str();
											bool was_valid_before = (eo_stream.error_code == me::OutputDescriptorStream::ERROR_CODE_OK);
											Response* fib_response = elem->setResponse(argv, argc, eo_stream);
											if (fib_response) {
												num_fib_responses = fib_comp_response->addResponse(fib_response);
											}
											else {
												if (was_valid_before && (eo_stream.error_code == me::OutputDescriptorStream::ERROR_CODE_SECTION_AFTER_FIBER)) {
													eo_descriptor.fixSectionAfterFiberDueToFiberOutputFail();
													eo_stream.error_code = me::OutputDescriptorStream::ERROR_CODE_OK;
												}
											}
										}
										if (num_fib_responses == 0)
											delete fib_comp_response;
										else
											num_mat_responses = mat_comp_response->addResponse(fib_comp_response);
									}
								}
							}
							// skip a bare "material" request (no sub-verb) -> no crash.
							else if (argc > (int)material_id_placeholder_index + 1) {
							int material_id = 0;
							while (true) {
								material_id++;
								if (material_id > LADRUNO_MAX_TRIAL_NSEC)
									break;
								std::stringstream ss_mid; ss_mid << material_id;
								std::string s_mid = ss_mid.str();
								argv[material_id_placeholder_index] = s_mid.c_str();
								Response* mat_response = elem->setResponse(argv, argc, eo_stream);
								if (mat_response)
									num_mat_responses = mat_comp_response->addResponse(mat_response);
								else
									break;
							}
							}
							if (num_mat_responses == 0)
								delete mat_comp_response;
							else
								eo_response = mat_comp_response;
						}
						else {
							eo_response = elem->setResponse(argv, argc, eo_stream);
						}
						eo_stream.finalizeSetResponse();
						if (do_all_fibers)
							eo_descriptor.purge();

						if (eo_response && (eo_stream.error_code == me::OutputDescriptorStream::ERROR_CODE_OK)) {
							me::OutputDescriptorHeader header(eo_descriptor.makeHeader());
							header.workaroundForSizeInconsistency(eo_response->getInformation().getData().Size());
							header.workaroundForDuplicatedComponents();
							me::OutputResponseCollection& eo_by_header = eo_by_custom_rule.items[header];
							if (eo_by_header.is_new) {
								std::stringstream ss_buf;
								ss_buf << elem_by_tag.class_tag << "-" << elem_by_tag.class_name
									<< "[" << elem_by_rule.int_rule_type << ":"
									<< elem_by_custom_rule.custom_int_rule_index << ":"
									<< header_local_index++ << "]";
								eo_by_header.dir_name = ss_buf.str();
								eo_by_header.is_new = false;
							}
							eo_by_header.items.push_back(me::OutputResponse(elem, eo_response));
							m_data->elemental_responses.push_back(eo_response);
							++num_request_answers;
						}
						else if (eo_response) {
							delete eo_response;  // composite built but descriptor stream errored (e.g. GENERIC) -> owned by nobody
							++num_class_stream_error;
						}
						else {
							// setResponse returned null: this element does not know the
							// token, or does not support it in this configuration, and it
							// never appears in the file. Only COUNT here;
							// the aggregated message is emitted once per element class below.
							++num_class_no_response;
						}
					}
				}
			}
			if (standard_section_keyword_modified) {
				argv[0] = aux_section_keyword_standard.c_str();
				standard_section_keyword_modified = false;
			}

			// ONE line per (request, element class) when this class dropped any
			// element. Deliberately NOT a token vocabulary check: the element
			// response token space is per-element and
			// open-ended, so "did anyone answer" is the only sound signal.
			size_t num_class_dropped = num_class_no_response + num_class_stream_error;
			if (num_class_dropped > 0) {
				opserr << "LadrunoRecorder warning: -E " << request_str.c_str()
				       << " : " << (int)num_class_dropped << " of " << (int)num_class_elements
				       << " " << elem_by_tag.class_name.c_str()
				       << " element(s) returned no response (unknown or unsupported token";
				if (num_class_stream_error > 0)
					opserr << "; " << (int)num_class_stream_error
					       << " of them built a response the output descriptor rejected";
				opserr << "); nothing will be recorded for them\n";
			}
		}

		// a request nobody answered produces no HDF5 group at all; without this
		// warning a typo in -E looks exactly like a working recorder with an
		// empty result.
		if (num_request_answers == 0) {
			opserr << "LadrunoRecorder warning: -E " << request_str.c_str()
			       << " matched no element in the model; nothing will be recorded for it\n";
		}
		delete[] argv;
	}

	// Now build one ElementResultSource + StreamingSink per (request, class, rule,
	// custom-rule, header) bucket. The bucket references live in
	// m_data->elemental_recorders, kept alive for the source lifetime.
	for (size_t cnt_request = 0; cnt_request < m_data->elemental_recorders.size(); ++cnt_request) {
		me::ResultRecorder& recorder = m_data->elemental_recorders[cnt_request];
		const std::vector<std::string>& request = recorder.result_request;
		for (me::ResultRecorder::collection_type::iterator it1 = recorder.response_map.begin();
			it1 != recorder.response_map.end(); ++it1) {
			me::OutputWithSameClassTagCollection& eo_by_tag = it1->second;
			for (me::OutputWithSameClassTagCollection::collection_type::iterator it2 = eo_by_tag.items.begin();
				it2 != eo_by_tag.items.end(); ++it2) {
				me::OutputWithSameIntRuleCollection& eo_by_rule = it2->second;
				for (me::OutputWithSameIntRuleCollection::collection_type::iterator it3 = eo_by_rule.items.begin();
					it3 != eo_by_rule.items.end(); ++it3) {
					me::OutputWithSameCustomIntRuleCollection& eo_by_custom_rule = it3->second;
					for (me::OutputWithSameCustomIntRuleCollection::collection_type::iterator it4 = eo_by_custom_rule.items.begin();
						it4 != eo_by_custom_rule.items.end(); ++it4) {
						const me::OutputDescriptorHeader& header = it4->first;
						me::OutputResponseCollection& bucket = it4->second;
						if (bucket.items.empty())
							continue;
						private_data::ElemChannel ech;
						ech.source = new ladrunons::ElementResultSource(request, header, bucket);
						ech.sink = m_data->envelope_mode
							? (ladrunons::ResultSink*)new ladrunons::EnvelopeSink(ladrunons::ResultFamily::OnElements)
							: (ladrunons::ResultSink*)new ladrunons::StreamingSink(ladrunons::ResultFamily::OnElements);
						ech.header = &header;
						ech.column_map_written = false;
						m_data->elem_channels.push_back(ech);
						// The element source name is "<display>/<bucket>" (nested).
						// The StreamingSink creates that as a 2-level path under
						// ON_ELEMENTS, so the intermediate <display> group ("force",
						// "stress", ...) must already exist. Pre-create it here once.
						{
							std::stringstream ss_disp_grp;
							ss_disp_grp << "MODEL_STAGE[" << info.current_model_stage_id
								<< "]/RESULTS/ON_ELEMENTS/" << m_data->elem_channels.back().source->schema().display_name;
							if (H5Lexists(info.h_file_id, ss_disp_grp.str().c_str(), H5P_DEFAULT) <= 0) {
								hid_t h_disp = ladrunons::h5::group::create(
									info.h_file_id, ss_disp_grp.str().c_str(), H5P_DEFAULT,
									info.h_group_proplist, H5P_DEFAULT);
								ladrunons::h5::group::close(h_disp);
							}
						}
					}
				}
			}
		}
	}

	return 0;
}

/* ===================================================================== */
/* per-step recording                                                    */
/* ===================================================================== */

int LadrunoRecorder::recordResultsOnNodes()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;

	if (m_data->nodes.empty() || m_data->node_channels.empty())
		return 0;

	// EIGEN/MODES detection.
	info.record_eigen_on_this_step = false;
	// gate on the DOMAIN's spectrum, not the interpreter's numEigen. `wipe`
	// clears the domain eigenvalues but never resets numEigen (Tcl global /
	// OpenSeesCommands), so after `eigen; wipe; <new model>` the old count
	// survives and Domain::getEigenvalues() would exit the process.
	// getNumEigenvalues() is the non-exiting presence probe.
	int num_eigen = info.domain->getNumEigenvalues();
	if (num_eigen > 0) {
		bool eigen_requested = false;
		for (size_t i = 0; i < m_data->node_channels.size(); ++i) {
			if (m_data->node_channels[i].is_modes) { eigen_requested = true; break; }
		}
		if (eigen_requested) {
			double eigen_set_time = info.domain->getTimeEigenvaluesSet();
			const Vector& new_eigenvalues = info.domain->getEigenvalues();
			if (info.eigen_first_initialization_done) {
				if (std::abs(eigen_set_time - info.eigen_last_time_set) > 1.0e-10) {
					info.record_eigen_on_this_step = true;
				}
				else if (num_eigen != info.eigen_last_values.Size()) {
					info.record_eigen_on_this_step = true;
				}
				else {
					for (int i = 0; i < num_eigen; ++i) {
						if (std::abs(new_eigenvalues[i] - info.eigen_last_values[i]) > 1.0e-10) {
							info.record_eigen_on_this_step = true;
							break;
						}
					}
				}
			}
			else {
				info.record_eigen_on_this_step = true;
				info.eigen_first_initialization_done = true;
			}
			if (info.record_eigen_on_this_step) {
				info.eigen_last_time_set = eigen_set_time;
				info.eigen_last_values = new_eigenvalues;
			}
		}
	}

	// Drive each channel. Group by reaction flag,
	// priming Domain::calculateNodalReactions(flag) only when the flag changes.
	int previous_reac_type = -1;
	std::vector<double> buffer;
	for (size_t i = 0; i < m_data->node_channels.size(); ++i) {
		private_data::NodeChannel& ch = m_data->node_channels[i];

		// remember the last flag actually COMPUTED, so a non-reaction channel
		// (-1) between two reaction channels (`reactionForce displacement
		// reactionMoment`) does not sweep every element's resisting force twice.
		int curr_reac_type = ch.reaction_flag;
		if (curr_reac_type > -1 && curr_reac_type < 3 && curr_reac_type != previous_reac_type) {
			info.domain->calculateNodalReactions(curr_reac_type);
			previous_reac_type = curr_reac_type;
		}

		if (ch.is_modes) {
			// Modes special path: gated on record_eigen_on_this_step; drive the
			// mode loop, writing one MODE_<k> per mode with LAMBDA/OMEGA/...
			if (!info.record_eigen_on_this_step)
				continue;
			recordModeChannel((int)i);
		}
		else {
			ch.source->evaluate(info, buffer);
			ch.sink->accept(info, *ch.source, buffer);
		}
	}

	return 0;
}

/* Drive the mode loop for a ModesOfVibration(Rotational) source and write one
   MODE_<k> dataset per mode under STEP_<commitTag>, as MPCORecorder does. The
   StreamingSink writes one DATA dataset per source; modes need a STEP group
   containing MODE_<k> datasets, so that structure is written directly through
   the h5 wrapper here. */
void LadrunoRecorder::recordModeChannel(int node_channel_index)
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	private_data::NodeChannel& ch = m_data->node_channels[(size_t)node_channel_index];
	ladrunons::ModesOfVibrationSource* msrc =
		dynamic_cast<ladrunons::ModesOfVibrationSource*>(ch.source);
	if (msrc == 0)
		return;

	const ladrunons::ResultSchema& schema = msrc->schema();
	if (schema.num_components < 1)
		return;

	// Modal layout (as MPCORecorder):
	//   <stage>/RESULTS/ON_NODES/<name>/{ID, DATA/STEP_<step>/MODE_<k>}
	// DATA here is a GROUP (one MODE_<k> dataset per mode per step), NOT the chunked
	// time-series DATA *dataset* that StreamingSink::begin would create; the two
	// collide on the name "DATA" (a DATA/STEP_<step> group cannot be created under
	// an existing dataset). So modal channels
	// do NOT use the streaming sink: create the result group + ID + DATA group here,
	// once per stage (idempotent via H5Lexists, re-created when the stage path moves).
	std::stringstream ss_result;
	ss_result << "MODEL_STAGE[" << info.current_model_stage_id
		<< "]/RESULTS/ON_NODES/" << schema.name;
	if (H5Lexists(info.h_file_id, ss_result.str().c_str(), H5P_DEFAULT) <= 0) {
		std::stringstream ss_family;
		ss_family << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES";
		hid_t h_family = H5Gopen2(info.h_file_id, ss_family.str().c_str(), H5P_DEFAULT);
		if (h_family >= 0) {
			hid_t h_gp_result = ladrunons::h5::group::createResultGroup(
				h_family, info.h_group_proplist,
				schema.name, schema.display_name,
				schema.components_csv, schema.num_components,
				schema.dimension, schema.description,
				(int)schema.result_type, (int)schema.data_type);
			const std::vector<int>& mode_ids = msrc->ids();
			hid_t h_dset_id = ladrunons::h5::dataset::createAndWrite(
				h_gp_result, "ID", mode_ids, mode_ids.size(), (size_t)1);
			ladrunons::h5::dataset::close(h_dset_id);
			hid_t h_data_grp = ladrunons::h5::group::create(
				h_gp_result, "DATA", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
			ladrunons::h5::group::close(h_data_grp);
			ladrunons::h5::group::close(h_gp_result);
			ladrunons::h5::group::close(h_family);
		}
	}

	// create the per-step group under DATA: <name>/DATA/STEP_<step>
	std::stringstream ss_step;
	ss_step << ss_result.str() << "/DATA/STEP_" << info.current_time_step_id;
	hid_t h_gp_step = ladrunons::h5::group::create(
		info.h_file_id, ss_step.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);
	ladrunons::h5::attribute::write(h_gp_step, "STEP", info.current_time_step_id);
	ladrunons::h5::attribute::write(h_gp_step, "TIME", info.current_time_step);

	const std::vector<int>& ids = msrc->ids();
	const size_t n_ids = ids.size();
	const size_t n_comp = (size_t)schema.num_components;

	std::vector<double> buffer;
	int nmodes = msrc->numModes(info);
	for (int k = 0; k < nmodes; ++k) {
		msrc->setCurrentMode(k);
		msrc->evaluate(info, buffer);
		double lambda, omega, freq, period;
		msrc->modeInfo(info, k, lambda, omega, freq, period);
		std::stringstream ss_mode; ss_mode << "MODE_" << k;
		hid_t h_dset = ladrunons::h5::dataset::createAndWrite(
			h_gp_step, ss_mode.str().c_str(), buffer, n_ids, n_comp);
		ladrunons::h5::attribute::write(h_dset, "MODE", k);
		ladrunons::h5::attribute::write(h_dset, "LAMBDA", lambda);
		ladrunons::h5::attribute::write(h_dset, "OMEGA", omega);
		ladrunons::h5::attribute::write(h_dset, "FREQUENCY", freq);
		ladrunons::h5::attribute::write(h_dset, "PERIOD", period);
		ladrunons::h5::dataset::close(h_dset);
	}

	ladrunons::h5::group::close(h_gp_step);
}

int LadrunoRecorder::recordResultsOnElements()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	if (m_data->elem_channels.empty())
		return 0;

	std::vector<double> buffer;
	for (size_t i = 0; i < m_data->elem_channels.size(); ++i) {
		private_data::ElemChannel& ech = m_data->elem_channels[i];
		ech.source->evaluate(info, buffer);
		ech.sink->accept(info, *ech.source, buffer);
		// Time-series mode: write the COLUMN_MAP metadata once, after the
		// StreamingSink created the ON_ELEMENTS result group on first accept().
		// In envelope mode the EnvelopeSink defers all output to finalize() and
		// delete-recreates its leaf group each flush, so the COLUMN_MAP is written
		// (and re-stamped) from finalizeAllSinks() under ENVELOPES/ON_ELEMENTS
		// instead; see writeElementColumnMap().
		if (!m_data->envelope_mode && !ech.column_map_written) {
			writeElementColumnMap((int)i);
			ech.column_map_written = true;
		}
	}
	return 0;
}

/* ===================================================================== */
/* domain / region sources (energy balance)                              */
/* ===================================================================== */

int LadrunoRecorder::initDomainSources()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;

	// energy is PARTITION_REDUCTION=UNSUPPORTED (shared-node KE counted in
	// every partition that holds the node; RES/ERR derived). Streaming still writes each
	// partition's partial, flagged, so nothing is lost -- but say so once.
	if (m_data->is_partitioned && !m_data->envelope_mode &&
	    (m_data->energy_requested || !m_data->energy_region_tags.empty()))
		opserr << "LadrunoRecorder WARNING: energyBalance in a partitioned run is written "
		          "per partition and is NOT mergeable (PARTITION_REDUCTION=UNSUPPORTED): "
		          "the per-step reduction is not implemented.\n";

	// whole-model energy -> ON_DOMAIN/energyBalance
	if (m_data->energy_requested) {
		// the recorder-owned state outlives this stage's source.
		ladrunons::ResultSource* src =
			new ladrunons::EnergyBalanceSource(info, &m_data->energy_model_state);
		if (refusePartitionedEnvelope(m_data->envelope_mode, m_data->is_partitioned, src)) {
			delete src;
		} else {
			private_data::DomainChannel ch;
			ch.source = src;
			ch.sink = m_data->envelope_mode
				? (ladrunons::ResultSink*)new ladrunons::EnvelopeSink(ladrunons::ResultFamily::OnDomain)
				: (ladrunons::ResultSink*)new ladrunons::StreamingSink(ladrunons::ResultFamily::OnDomain);
			m_data->domain_channels.push_back(ch);
		}
	}

	// per-region energy -> ON_REGIONS/energyBalance
	if (!m_data->energy_region_tags.empty()) {
		ladrunons::ResultSource* src =
			new ladrunons::EnergyBalanceSource(info, m_data->energy_region_tags,
			                                   &m_data->energy_region_state);
		if (refusePartitionedEnvelope(m_data->envelope_mode, m_data->is_partitioned, src)) {
			delete src;
		} else {
			private_data::DomainChannel ch;
			ch.source = src;
			ch.sink = m_data->envelope_mode
				? (ladrunons::ResultSink*)new ladrunons::EnvelopeSink(ladrunons::ResultFamily::OnRegions)
				: (ladrunons::ResultSink*)new ladrunons::StreamingSink(ladrunons::ResultFamily::OnRegions);
			m_data->domain_channels.push_back(ch);
		}
	}

	return 0;
}

int LadrunoRecorder::recordResultsOnDomain()
{
	ladrunons::detail::ProcessInfo& info = m_data->info;

	std::vector<double> buffer;
	for (size_t i = 0; i < m_data->domain_channels.size(); ++i) {
		private_data::DomainChannel& ch = m_data->domain_channels[i];
		ch.source->evaluate(info, buffer);
		ch.sink->accept(info, *ch.source, buffer);
	}
	return 0;
}

/* Write the COLUMN_MAP child group + SECTION_MAP under the element result group
   Called once, after the StreamingSink has created
   ON_ELEMENTS/<display>/<bucket>. */
void LadrunoRecorder::writeElementColumnMap(int elem_channel_index)
{
	ladrunons::detail::ProcessInfo& info = m_data->info;
	private_data::ElemChannel& ech = m_data->elem_channels[(size_t)elem_channel_index];
	const ladrunons::detail::element::OutputDescriptorHeader& header = *ech.header;
	const ladrunons::ResultSchema& schema = ech.source->schema();

	// Open the bucket result group. In time-series mode that is
	// ON_ELEMENTS/<display>/<bucket> (created by the StreamingSink); in envelope
	// mode the EnvelopeSink writes the bucket under ENVELOPES/ON_ELEMENTS/<...>
	// instead, so the per-column structure must land there too (the result-level
	// COMPONENTS attr is empty for elements; COLUMN_MAP is the authoritative
	// gauss/section/fiber column map).
	std::stringstream ss_grp;
	ss_grp << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/";
	if (m_data->envelope_mode)
		ss_grp << "ENVELOPES/ON_ELEMENTS/" << schema.name;
	else
		ss_grp << "ON_ELEMENTS/" << schema.name;
	hid_t h_grp = H5Gopen2(info.h_file_id, ss_grp.str().c_str(), H5P_DEFAULT);
	if (h_grp < 0)
		return;

	ladrunons::h5::attribute::write(h_grp, "NUM_COLUMNS", header.num_columns);

	// COLUMN_MAP structured metadata. One row per descriptor block.
	const size_t k = header.num_components.size();
	std::vector<int> levels(k, 0);   // block path depth code = last element of path
	std::vector<int> gauss_id(k, -1);
	std::vector<int> section_tag(k, -1); // section tag not tracked at header level: -1
	std::vector<int> fiber_id(k, -1);    // fiber index not tracked at header level: -1
	std::vector<int> num_comp(k, 0);
	std::vector<int> multiplicity(k, 0);
	std::stringstream ss_comp_names;
	for (size_t i = 0; i < k; ++i) {
		// LEVELS: deepest path code in this block's components_path
		int depth_code = 0;
		if (i < header.components_path.size() && !header.components_path[i].empty())
			depth_code = header.components_path[i].back();
		levels[i] = depth_code;
		gauss_id[i] = (i < header.gauss_id.size()) ? header.gauss_id[i] : -1;
		num_comp[i] = header.num_components[i];
		multiplicity[i] = (i < header.multiplicity.size()) ? header.multiplicity[i] : 1;
		// COMP_NAMES: one line per block, CSV of component names
		if (i > 0) ss_comp_names << "\n";
		if (i < header.components.size()) {
			const std::vector<std::string>& cn = header.components[i];
			for (size_t j = 0; j < cn.size(); ++j) {
				if (j > 0) ss_comp_names << ",";
				ss_comp_names << cn[j];
			}
		}
	}

	hid_t h_cm = ladrunons::h5::group::create(
		h_grp, "COLUMN_MAP", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	if (k > 0) {
		hid_t d;
		d = ladrunons::h5::dataset::createAndWrite(h_cm, "LEVELS", levels, k, 1); ladrunons::h5::dataset::close(d);
		d = ladrunons::h5::dataset::createAndWrite(h_cm, "GAUSS_ID", gauss_id, k, 1); ladrunons::h5::dataset::close(d);
		d = ladrunons::h5::dataset::createAndWrite(h_cm, "SECTION_TAG", section_tag, k, 1); ladrunons::h5::dataset::close(d);
		d = ladrunons::h5::dataset::createAndWrite(h_cm, "FIBER_ID", fiber_id, k, 1); ladrunons::h5::dataset::close(d);
		d = ladrunons::h5::dataset::createAndWrite(h_cm, "NUM_COMP", num_comp, k, 1); ladrunons::h5::dataset::close(d);
		d = ladrunons::h5::dataset::createAndWrite(h_cm, "MULTIPLICITY", multiplicity, k, 1); ladrunons::h5::dataset::close(d);
	}
	ladrunons::h5::attribute::write(h_cm, "COMP_NAMES", ss_comp_names.str());
	ladrunons::h5::group::close(h_cm);

	ladrunons::h5::group::close(h_grp);
}

/* ===================================================================== */
/* parallel: broadcast the recorder CONFIG (not results) to worker ranks */
/* ===================================================================== */
//
// Under OpenSeesMP the recorder is created once (on P0) and shipped to each
// worker via sendSelf/recvSelf. Only the lightweight config is serialized: the
// requests, filename, frequency, region/sets, kind. Each worker then builds its
// own sources from ITS partition of the domain and writes its own
// "<stem>.part-<p_id>.ladruno" (result data is never communicated, as in
// MPCORecorder, whose per-step comm is just the one stage-stamp Allreduce in
// record()). p_id is carried as P0's send_self_count at send time.

namespace {
// Minimal self-contained binary (de)serializer. Same-architecture transport
// (MPCORecorder likewise ships raw Message bytes), so native layout is fine.
struct LadSer {
	std::string b;
	void put_i(int v) { b.append(reinterpret_cast<const char*>(&v), sizeof(int)); }
	void put_d(double v) { b.append(reinterpret_cast<const char*>(&v), sizeof(double)); }
	void put_s(const std::string& v) { put_i(static_cast<int>(v.size())); b.append(v); }
	void put_vi(const std::vector<int>& v) {
		put_i(static_cast<int>(v.size()));
		for (size_t k = 0; k < v.size(); ++k) put_i(v[k]);
	}
};
struct LadDeser {
	const char* p; const char* end; bool ok;
	LadDeser(const char* data, int n) : p(data), end(data + n), ok(true) {}
	int get_i() {
		if (!ok || p + sizeof(int) > end) { ok = false; return 0; }
		int v; std::memcpy(&v, p, sizeof(int)); p += sizeof(int); return v;
	}
	double get_d() {
		if (!ok || p + sizeof(double) > end) { ok = false; return 0.0; }
		double v; std::memcpy(&v, p, sizeof(double)); p += sizeof(double); return v;
	}
	std::string get_s() {
		int n = get_i();
		if (!ok || n < 0 || p + n > end) { ok = false; return std::string(); }
		std::string v(p, n); p += n; return v;
	}
	std::vector<int> get_vi() {
		int n = get_i(); std::vector<int> v;
		if (!ok || n < 0) { ok = false; return v; }
		v.reserve(n);
		for (int k = 0; k < n && ok; ++k) v.push_back(get_i());
		return v;
	}
};
} // anonymous namespace

int LadrunoRecorder::sendSelf(int commitTag, Channel& theChannel)
{
	if (theChannel.isDatastore() == 1) {
		opserr << "LadrunoRecorder::sendSelf() - does not send data to a datastore\n";
		return -1;
	}
	m_data->send_self_count++;
	if (m_data->run_id.empty()) {   // P0 owns the run identity
		ladrunoResolveRunId(m_data->run_id, m_data->run_id_scope);
		if (m_data->run_id_scope == "process")
			m_data->run_id_scope = "broadcast";
	}

	LadSer ser;
	ser.put_i(this->getTag());
	ser.put_i(m_data->send_self_count);   // becomes p_id in the receiver
	ser.put_s(m_data->filename);
	ser.put_i((int)m_data->output_freq.type);
	ser.put_d(m_data->output_freq.dt);
	ser.put_i(m_data->output_freq.nsteps);
	ser.put_i(m_data->has_region ? 1 : 0);
	ser.put_i(m_data->region_tag);
	ser.put_vi(m_data->node_set);
	ser.put_vi(m_data->elem_set);
	ser.put_i(m_data->energy_requested ? 1 : 0);
	ser.put_vi(m_data->energy_region_tags);
	ser.put_s(m_data->stage_kind);
	{
		std::vector<int> nr;
		nr.reserve(m_data->nodal_results_requests.size());
		for (size_t k = 0; k < m_data->nodal_results_requests.size(); ++k)
			nr.push_back((int)m_data->nodal_results_requests[k]);
		ser.put_vi(nr);
	}
	ser.put_vi(m_data->sens_grad_indices);
	ser.put_i((int)m_data->elemental_results_requests.size());
	for (size_t k = 0; k < m_data->elemental_results_requests.size(); ++k) {
		const std::vector<std::string>& req = m_data->elemental_results_requests[k];
		ser.put_i((int)req.size());
		for (size_t j = 0; j < req.size(); ++j) ser.put_s(req[j]);
	}

	ser.put_i(m_data->envelope_mode ? 1 : 0);    // lockstep with recvSelf
	ser.put_i(m_data->info.store_data_f32 ? 1 : 0);
	ser.put_i(m_data->info.deflate_level);
	ser.put_d(m_data->flush_seconds);
	ser.put_s(m_data->run_id);
	ser.put_s(m_data->run_id_scope);
	int msg_data_size = (int)ser.b.size();
	ID idata(1);
	idata(0) = msg_data_size;
	if (theChannel.sendID(0, commitTag, idata) < 0) {
		opserr << "LadrunoRecorder::sendSelf() - failed to send message size\n";
		return -1;
	}
	std::vector<char> msg_data(static_cast<size_t>(msg_data_size) + 1, '\0');
	std::copy(ser.b.begin(), ser.b.end(), msg_data.begin());
	Message msg(msg_data.data(), msg_data_size);
	if (theChannel.sendMsg(0, commitTag, msg) < 0) {
		opserr << "LadrunoRecorder::sendSelf() - failed to send message\n";
		return -1;
	}
	return 0;
}

int LadrunoRecorder::recvSelf(int commitTag, Channel& theChannel,
                                  FEM_ObjectBroker& /*theBroker*/)
{
	if (theChannel.isDatastore() == 1) {
		opserr << "LadrunoRecorder::recvSelf() - does not recv data from a datastore\n";
		return -1;
	}

	// fresh state: this process is a secondary/worker (send_self_count < 0)
	if (m_data) {
		clearSources();
		delete m_data;
	}
	m_data = new private_data();
	m_data->send_self_count = -1;

	ID idata(1);
	if (theChannel.recvID(0, commitTag, idata) < 0) {
		opserr << "LadrunoRecorder::recvSelf() - failed to recv message size\n";
		return -1;
	}
	int msg_data_size = idata(0);
	std::vector<char> msg_data(static_cast<size_t>(msg_data_size) + 1, '\0');
	Message msg(msg_data.data(), msg_data_size);
	if (theChannel.recvMsg(0, commitTag, msg) < 0) {
		opserr << "LadrunoRecorder::recvSelf() - failed to recv message\n";
		return -1;
	}

	LadDeser de(msg_data.data(), msg_data_size);
	int my_tag = de.get_i();
	m_data->p_id = de.get_i();             // P0's counter at send time -> partition idx
	m_data->filename = de.get_s();
	m_data->output_freq.type =
		(ladrunons::detail::OutputFrequency::IncrementType)de.get_i();
	m_data->output_freq.dt = de.get_d();
	m_data->output_freq.nsteps = de.get_i();
	m_data->has_region = de.get_i() != 0;
	m_data->region_tag = de.get_i();
	m_data->node_set = de.get_vi();
	m_data->elem_set = de.get_vi();
	m_data->energy_requested = de.get_i() != 0;
	m_data->energy_region_tags = de.get_vi();
	m_data->stage_kind = de.get_s();
	{
		std::vector<int> nr = de.get_vi();
		m_data->nodal_results_requests.clear();
		for (size_t k = 0; k < nr.size(); ++k)
			m_data->nodal_results_requests.push_back(
				(ladrunons::detail::NodalResultType::Enum)nr[k]);
	}
	m_data->sens_grad_indices = de.get_vi();
	{
		int n = de.get_i();
		m_data->elemental_results_requests.clear();
		for (int k = 0; k < n && de.ok; ++k) {
			int m = de.get_i();
			std::vector<std::string> req;
			for (int j = 0; j < m && de.ok; ++j) req.push_back(de.get_s());
			m_data->elemental_results_requests.push_back(req);
		}
	}
	m_data->envelope_mode = de.get_i() != 0;    // lockstep with sendSelf
	m_data->info.store_data_f32 = de.get_i() != 0;
	m_data->info.deflate_level = de.get_i();
	m_data->flush_seconds = de.get_d();
	m_data->run_id = de.get_s();
	m_data->run_id_scope = de.get_s();
	if (!de.ok) {
		opserr << "LadrunoRecorder::recvSelf() - failed to de-serialize config\n";
		return -1;
	}
	this->setTag(my_tag);
	return 0;
}

/* ===================================================================== */
/* command parser  (recorder ladruno <file> ...)                     */
/* ===================================================================== */

void* OPS_LadrunoRecorder()
{
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 1) {
		opserr << "LadrunoRecorder error: insufficient args; expected a filename\n"
		       << "  usage: recorder ladruno <file> [-N <res...>] [-NS <res...> <grad>] "
		          "[-E <res...>] [-R <regionTag>] [-G energy <regionTag...>] "
		          "[-T dt|nsteps <v>] [-envelope] [-kind <k>] [-precision f32|f64] "
		          "[-compress <0..9>] [-flush <s>]\n";
		return 0;
	}

	const char* filename = OPS_GetString();
	numdata--;

	Domain* domain = OPS_GetDomain();
	if (domain == 0) {
		opserr << "LadrunoRecorder error: domain is not defined\n";
		return 0;
	}

	// parse optional arguments (same syntax as OPS_MPCORecorder, plus -G and the options below)
	ladrunons::utils::parsing::option_type curr_opt = ladrunons::utils::parsing::opt_none;
	std::vector<ladrunons::detail::NodalResultType::Enum> nodal_results_requests;
	std::vector<int> sens_grad_indices;
	std::vector<std::vector<std::string> > elemental_results_requests;
	std::vector<std::string> tokens;
	ladrunons::detail::OutputFrequency output_freq;
	bool has_region = false;
	int region_tag = 0;                  // -R region tag (for MODEL/SETS)
	std::set<int> node_set;
	std::set<int> elem_set;
	bool energy_requested = false;       // -G energy
	std::vector<int> energy_region_tags; // -G <regionTag...>
	std::string stage_kind_opt;            // -kind <transient|static|eigen>; "" = automatic
	bool envelope_opt = false;             // -envelope flag
	int compress_opt = 1;                  // -compress <0..9> (default 1)
	double flush_opt = 10.0;               // -flush <seconds>
	bool store_data_f32 = false;           // -precision f32 (lossy) | f64 (default)
	int one_item = 1;

	while (numdata > 0) {
		const char* data = OPS_GetString();
		numdata--;
		if (strcmp(data, "-N") == 0) {
			curr_opt = ladrunons::utils::parsing::opt_result_on_nodes;
		}
		else if (strcmp(data, "-NS") == 0) {
			curr_opt = ladrunons::utils::parsing::opt_result_on_nodes_sens;
		}
		else if (strcmp(data, "-E") == 0) {
			curr_opt = ladrunons::utils::parsing::opt_result_on_elements;
		}
		else if (strcmp(data, "-T") == 0) {
			curr_opt = ladrunons::utils::parsing::opt_time;
			output_freq.reset();
		}
		else if (strcmp(data, "-R") == 0) {
			curr_opt = ladrunons::utils::parsing::opt_region;
			if (numdata > 0) {
				int rtag = 0;
				if (OPS_GetInt(&one_item, &rtag) != 0) {
					opserr << "LadrunoRecorder error: option -R (region) requires an int region tag\n";
					return 0;
				}
				MeshRegion* region = domain->getRegion(rtag);
				if (region == 0) {
					opserr << "LadrunoRecorder error: region " << rtag << " is null\n";
					return 0;
				}
				const ID& node_ids = region->getNodes();
				const ID& elem_ids = region->getElements();
				for (int i = 0; i < node_ids.Size(); ++i)
					node_set.insert(node_ids(i));
				for (int i = 0; i < elem_ids.Size(); ++i)
					elem_set.insert(elem_ids(i));
				has_region = true;
				region_tag = rtag; // retained for MODEL/SETS self-description
				numdata--;
			}
			else {
				opserr << "LadrunoRecorder error: option -R (region) requires an int region tag\n";
				return 0;
			}
			has_region = true;
		}
		else if (strcmp(data, "-G") == 0) {
			// -G energy <regionTag...> : whole-model energy balance (ON_DOMAIN)
			// + optional per-region balance (ON_REGIONS).
			// Consume eagerly: the "energy" keyword is a real string arg, but the
			// region tags are integers -> they MUST be read with OPS_GetIntInput,
			// not OPS_GetString (which returns "" for a numeric OpenSeesPy arg,
			// so atoi() would silently yield 0). Same idiom as -node.
			curr_opt = ladrunons::utils::parsing::opt_global;
			if (numdata > 0) {
				const char* gkind = OPS_GetString();
				numdata--;
				// `-G -T nsteps 10` must not eat "-T", read "nsteps"/"10" as
				// region tags and drop the output frequency: an option token is
				// handed back.
				if (gkind[0] == '-' && gkind[1] != '\0' &&
				    !(gkind[1] >= '0' && gkind[1] <= '9')) {
					OPS_ResetCurrentInputArg(-1);
					numdata++;
					opserr << "LadrunoRecorder warning: -G expects 'energy' "
					          "[regionTag...]; got option " << gkind << "\n";
					continue;
				}
				if (strcmp(gkind, "energy") == 0) {
					energy_requested = true;
				}
				else {
					// tolerate a tag passed as a string (Tcl path)
					int t = atoi(gkind);
					if (t != 0 && domain->getRegion(t) != 0)
						energy_region_tags.push_back(t);
					else
						opserr << "LadrunoRecorder warning: -G expects 'energy' "
						          "[regionTag...]; ignoring token (" << gkind << ")\n";
				}
			}
			// 0+ integer region tags follow; stop at the first non-int (an option
			// like -T) and hand it back to the outer parse loop.
			while (numdata > 0) {
				int tag = 0; int n = 1;
				int oldn = OPS_GetNumRemainingInputArgs();
				if (OPS_GetIntInput(&n, &tag) < 0) {
					// un-consume the flag only if the failed read consumed it:
					// openseespy advances the cursor on a failed int parse, the
					// classic Tcl elementAPI does not; a blind -1 rewind there
					// overshoots onto the previous token and desyncs numdata.
					if (OPS_GetNumRemainingInputArgs() < oldn)
						OPS_ResetCurrentInputArg(-1);
					break;
				}
				numdata--;
				MeshRegion* gregion = domain->getRegion(tag);
				if (gregion == 0)
					opserr << "LadrunoRecorder warning: -G region " << tag
					       << " not found in the domain; skipping its energy block\n";
				else
					energy_region_tags.push_back(tag);
			}
		}
		else if (strcmp(data, "-kind") == 0) {
			// -kind <transient|static|eigen> : the MODEL_STAGE/KIND attr.
			if (numdata > 0) {
				const char* k = OPS_GetString();
				numdata--;
				if (strcmp(k, "transient") == 0 || strcmp(k, "static") == 0 ||
				    strcmp(k, "eigen") == 0)
					stage_kind_opt = k;
				else
					opserr << "LadrunoRecorder warning: -kind expects "
					          "transient|static|eigen; got (" << k << "), keeping "
					       << (stage_kind_opt.empty() ? "the automatic per-stage KIND"
					                                  : stage_kind_opt.c_str()) << "\n";
			}
		}
		else if (strcmp(data, "-envelope") == 0) {
			// flag (no argument): record componentwise running extremes (MIN/MAX/
			// ABSMAX + ARG_STEP under RESULTS/ENVELOPES) instead of the
			// per-step time series. Does not change the active request mode.
			envelope_opt = true;
		}
		else if (strcmp(data, "-precision") == 0) {
			// -precision f32|f64 : on-disk element type for the per-step result
			// DATA. f64 (default) is the lossless parity path; f32 is the opt-in
			// lossy mode (~half the payload, ~7 significant digits). Anything else
			// warns and keeps f64.
			if (numdata > 0) {
				const char* p = OPS_GetString();
				numdata--;
				std::string ps(p);
				if (ps == "f32" || ps == "float32" || ps == "single")
					store_data_f32 = true;
				else if (ps == "f64" || ps == "float64" || ps == "double")
					store_data_f32 = false;
				else
					opserr << "LadrunoRecorder warning: -precision expects f32|f64; got ("
					       << p << "), keeping f64 (lossless)\n";
			}
		}
		else if (strcmp(data, "-compress") == 0) {
			// -compress <0..9> : deflate level of the time-series DATA (0 = no
			// shuffle/deflate filter). Default 1: on large slabs it costs 19 % less
			// recorder CPU than level 4 for a 0.6 % larger file.
			if (numdata > 0) {
				int lvl = 1;
				if (OPS_GetInt(&one_item, &lvl) != 0) {
					opserr << "LadrunoRecorder error: -compress requires an int level 0..9\n";
					return 0;
				}
				numdata--;
				if (lvl < 0 || lvl > 9) {
					opserr << "LadrunoRecorder warning: -compress level " << lvl
					       << " clamped to [0, 9]\n";
					lvl = (lvl < 0) ? 0 : 9;
				}
				compress_opt = lvl;
			}
			else {
				opserr << "LadrunoRecorder error: -compress requires a level 0..9\n";
				return 0;
			}
		}
		else if (strcmp(data, "-flush") == 0) {
			// -flush <seconds> : wall-clock interval between file flushes (and
			// envelope rewrites). 0 = every recorded step. The file is always
			// flushed at a stage change and on close.
			if (numdata > 0) {
				double sec = 10.0;
				if (OPS_GetDouble(&one_item, &sec) != 0) {
					opserr << "LadrunoRecorder error: -flush requires a number of seconds\n";
					return 0;
				}
				numdata--;
				flush_opt = (sec < 0.0) ? 0.0 : sec;
			}
			else {
				opserr << "LadrunoRecorder error: -flush requires a number of seconds\n";
				return 0;
			}
		}
		else {
			switch (curr_opt) {
			case ladrunons::utils::parsing::opt_result_on_nodes: {
				if (strcmp(data, "displacement") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::Displacement);
				else if (strcmp(data, "rotation") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::Rotation);
				else if (strcmp(data, "velocity") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::Velocity);
				else if (strcmp(data, "angularVelocity") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::AngularVelocity);
				else if (strcmp(data, "acceleration") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::Acceleration);
				else if (strcmp(data, "angularAcceleration") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::AngularAcceleration);
				else if (strcmp(data, "reactionForce") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::ReactionForce);
				else if (strcmp(data, "reactionMoment") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::ReactionMoment);
				else if (strcmp(data, "reactionForceIncludingInertia") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::ReactionForceIncludingInertia);
				else if (strcmp(data, "reactionMomentIncludingInertia") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::ReactionMomentIncludingInertia);
				else if (strcmp(data, "rayleighForce") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::RayleighForce);
				else if (strcmp(data, "rayleighMoment") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::RayleighMoment);
				else if (strcmp(data, "unbalancedForce") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::UnbalancedForce);
				else if (strcmp(data, "unbalancedMoment") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::UnbalancedMoment);
				else if (strcmp(data, "unbalancedForceIncludingInertia") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::UnbalancedForceIncludingInertia);
				else if (strcmp(data, "unbalancedMomentIncludingInertia") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::UnbalancedMomentIncludingInertia);
				else if (strcmp(data, "pressure") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::Pressure);
				else if (strcmp(data, "modesOfVibration") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::ModesOfVibration);
				else if (strcmp(data, "modesOfVibrationRotational") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::ModesOfVibrationRotational);
				else {
					opserr << "LadrunoRecorder error: option -N with unknown result type (" << data << ")\n";
					return 0;
				}
				sens_grad_indices.push_back(0); // placeholder
				break;
			}
			case ladrunons::utils::parsing::opt_result_on_nodes_sens: {
				if (strcmp(data, "displacementSensitivity") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::DisplacementSensitivity);
				else if (strcmp(data, "rotationSensitivity") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::RotationSensitivity);
				else if (strcmp(data, "velocitySensitivity") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::VelocitySensitivity);
				else if (strcmp(data, "angularVelocitySensitivity") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::AngularVelocitySensitivity);
				else if (strcmp(data, "accelerationSensitivity") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::AccelerationSensitivity);
				else if (strcmp(data, "angularAccelerationSensitivity") == 0)
					nodal_results_requests.push_back(ladrunons::detail::NodalResultType::AngularAccelerationSensitivity);
				else {
					opserr << "LadrunoRecorder error: option -NS with unknown result type (" << data << ")\n";
					return 0;
				}
				if (numdata > 0) {
					int grad_index;
					if (OPS_GetInt(&one_item, &grad_index) != 0) {
						opserr << "LadrunoRecorder error: option -NS requires an int sensitivity parameter index\n";
						return 0;
					}
					numdata--;
					sens_grad_indices.push_back(grad_index);
				}
				else {
					opserr << "LadrunoRecorder error: option -NS requires a sensitivity parameter index\n";
					return 0;
				}
				break;
			}
			case ladrunons::utils::parsing::opt_result_on_elements: {
				std::string temp(data);
				ladrunons::utils::strings::split(temp, '.', tokens, true);
				if (tokens.size() > 0)
					elemental_results_requests.push_back(tokens);
				break;
			}
			case ladrunons::utils::parsing::opt_time: {
				if (strcmp(data, "dt") == 0) {
					output_freq.type = ladrunons::detail::OutputFrequency::DeltaTime;
					output_freq.nsteps = 1;
					if (numdata > 0) {
						if (OPS_GetDouble(&one_item, &output_freq.dt) != 0) {
							opserr << "LadrunoRecorder error: invalid double argument for the delta time\n";
							return 0;
						}
						if (output_freq.dt < 0.0) output_freq.dt = 0.0;
						numdata--;
					}
					else {
						opserr << "LadrunoRecorder error: option -T dt requires a delta time argument\n";
						return 0;
					}
				}
				else if (strcmp(data, "nsteps") == 0) {
					output_freq.type = ladrunons::detail::OutputFrequency::NumberOfSteps;
					output_freq.dt = 0.0;
					if (numdata > 0) {
						if (OPS_GetInt(&one_item, &output_freq.nsteps) != 0) {
							opserr << "LadrunoRecorder error: invalid int argument for the number of steps\n";
							return 0;
						}
						if (output_freq.nsteps < 1) output_freq.nsteps = 1;
						numdata--;
					}
					else {
						opserr << "LadrunoRecorder error: option -T nsteps requires a number-of-steps argument\n";
						return 0;
					}
				}
				else {
					opserr << "LadrunoRecorder error: option -T with unknown frequency type (" << data << ")\n";
					return 0;
				}
				break;
			}
			case ladrunons::utils::parsing::opt_global: {
				// -G energy           -> whole-model energy balance (ON_DOMAIN)
				// -G energy <tag...>  -> additionally per-region (ON_REGIONS)
				if (strcmp(data, "energy") == 0) {
					energy_requested = true;
				}
				else {
					// interpret as a region tag (same int-parse style as -R)
					int tag = atoi(data);
					MeshRegion* region = domain->getRegion(tag);
					if (region == 0) {
						opserr << "LadrunoRecorder warning: option -G region " << tag
						       << " not found in the domain; skipping its energy block\n";
					}
					else {
						energy_region_tags.push_back(tag);
					}
				}
				break;
			}
			default: {
				opserr << "LadrunoRecorder error: unknown arg with option none " << data << "\n";
				return 0;
			}
			}
		}
	}

	// drop repeated requests at parse time (with a notice). Two channels for
	// one result share one HDF5 group name; the sink refuses the second, but
	// the request is a typo or an alias pair; say so here, once.
	{
		std::vector<ladrunons::detail::NodalResultType::Enum> n_keep;
		std::vector<int> g_keep;
		for (size_t i = 0; i < nodal_results_requests.size(); ++i) {
			const int g = (i < sens_grad_indices.size()) ? sens_grad_indices[i] : 0;
			bool dup = false;
			for (size_t j = 0; j < n_keep.size(); ++j)
				if (n_keep[j] == nodal_results_requests[i] && g_keep[j] == g) { dup = true; break; }
			if (dup) {
				opserr << "LadrunoRecorder warning: nodal result requested twice (an alias "
				          "or a repeated -N token); recording it once\n";
				continue;
			}
			n_keep.push_back(nodal_results_requests[i]);
			g_keep.push_back(g);
		}
		nodal_results_requests.swap(n_keep);
		sens_grad_indices.swap(g_keep);

		std::vector<std::vector<std::string> > e_keep;
		for (size_t i = 0; i < elemental_results_requests.size(); ++i) {
			bool dup = false;
			for (size_t j = 0; j < e_keep.size(); ++j)
				if (e_keep[j] == elemental_results_requests[i]) { dup = true; break; }
			if (dup) {
				opserr << "LadrunoRecorder warning: element result requested twice (-E";
				for (size_t k = 0; k < elemental_results_requests[i].size(); ++k)
					opserr << (k ? "." : " ") << elemental_results_requests[i][k].c_str();
				opserr << "); recording it once\n";
				continue;
			}
			e_keep.push_back(elemental_results_requests[i]);
		}
		elemental_results_requests.swap(e_keep);

		std::vector<int> r_keep;
		for (size_t i = 0; i < energy_region_tags.size(); ++i)
			if (std::find(r_keep.begin(), r_keep.end(), energy_region_tags[i]) == r_keep.end())
				r_keep.push_back(energy_region_tags[i]);
		energy_region_tags.swap(r_keep);
	}

	LadrunoRecorder* recorder = new LadrunoRecorder();
	recorder->m_data->filename = filename;
	recorder->m_data->output_freq = output_freq;
	recorder->m_data->nodal_results_requests.swap(nodal_results_requests);
	recorder->m_data->sens_grad_indices.swap(sens_grad_indices);
	recorder->m_data->elemental_results_requests.swap(elemental_results_requests);
	recorder->m_data->has_region = has_region;
	recorder->m_data->region_tag = region_tag;
	for (std::set<int>::const_iterator it = node_set.begin(); it != node_set.end(); ++it)
		recorder->m_data->node_set.push_back(*it);
	for (std::set<int>::const_iterator it = elem_set.begin(); it != elem_set.end(); ++it)
		recorder->m_data->elem_set.push_back(*it);
	recorder->m_data->energy_requested = energy_requested;
	recorder->m_data->energy_region_tags.swap(energy_region_tags);
	recorder->m_data->stage_kind = stage_kind_opt;
	recorder->m_data->envelope_mode = envelope_opt;
	recorder->m_data->info.deflate_level = compress_opt;
	recorder->m_data->flush_seconds = flush_opt;
	// -precision: carried on ProcessInfo (read by StreamingSink::begin and stamped
	// into INFO/STORED_PRECISION at initialize()). setDomain() does not overwrite it.
	recorder->m_data->info.store_data_f32 = store_data_f32;
	return recorder;
}
