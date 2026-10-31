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
// Description: persistence strategies of the ladruno recorder. StreamingSink
// appends one slab per recorded step to a chunked DATA[T x nIds x nComp]
// dataset with TIME[T] and STEP[T]; EnvelopeSink keeps the componentwise
// MIN, MAX, ABSMAX and the step of each extreme.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#ifndef Ladruno_Sinks_h
#define Ladruno_Sinks_h

#include "Ladruno_ResultIO.h"
#include "Ladruno_Types.h"

#include <string>
#include <vector>

namespace ladruno {

	/*
	Which RESULTS sub-tree a sink writes into (ON_NODES, ON_ELEMENTS,
	ON_DOMAIN, ON_REGIONS), so one sink class serves every source family.
	Also selects the ENVELOPES sub-tree.
	*/
	struct ResultFamily {
		enum Enum {
			OnNodes = 0,
			OnElements,
			OnDomain,
			OnRegions
		};
		// schema group name: "ON_NODES" / "ON_ELEMENTS" / "ON_DOMAIN" /
		// "ON_REGIONS"
		static const char* groupName(Enum f);
	};

	/*
	StreamingSink: one slab per recorded step.

	  begin():  create the RESULTS/<family>/<name> result group (DISPLAY_NAME /
	            COMPONENTS / DIMENSION / DESCRIPTION / TYPE / DATA_TYPE /
	            PARTITION_REDUCTION attrs), the ID dataset [nIds x 1] from
	            src.ids(), and the empty extensible DATA [0 x nIds x nComp],
	            TIME [0] and STEP [0].
	  accept(): append the buffer as slab t of DATA, the time to TIME and the
	            commit tag to STEP.
	  finalize(): nothing deferred (the handles stay open until destruction).
	*/
	class StreamingSink : public ResultSink {
	public:
		explicit StreamingSink(ResultFamily::Enum family)
			: m_family(family), m_initialized(false), m_dead(false),
			  m_warned_short(false), m_data(HID_INVALID), m_time(HID_INVALID),
			  m_step(HID_INVALID), m_mspace(HID_INVALID), m_mspace1(HID_INVALID),
			  m_t(0), m_n_ids(0), m_n_comp(0) {}
		// closes the DATA/TIME/STEP handles held open for the stage
		// (the recorder deletes sinks on a stage change and before H5Fclose).
		~StreamingSink() override;

		void begin(detail::ProcessInfo& info, const ResultSource& src) override;
		void accept(detail::ProcessInfo& info, const ResultSource& src,
		            const std::vector<double>& buffer) override;
		void finalize(detail::ProcessInfo& info) override;

		// New MODEL_STAGE / re-begin: forget that the group was created.
		void reset() { m_initialized = false; m_dead = false; m_warned_short = false; }

	private:
		// report a write failure ONCE and stop this channel, instead of
		// silently returning on every later step (the result would vanish from
		// the file with only HDF5-DIAG output on stderr).
		void fail(const ResultSchema& schema, const char* why);
		void closeHandles();

		ResultFamily::Enum m_family;
		bool m_initialized;   // group + ID written for the current stage
		bool m_dead;          // a create/append failed; channel stopped
		bool m_warned_short;  // short-buffer warning printed once
		// DATA/TIME/STEP stay open from begin() to destruction, so the
		// per-dataset chunk cache keeps the partially filled chunk in memory and
		// it is deflated once when complete (reopening would evict and
		// re-inflate it every step). m_t counts the slabs written.
		hid_t m_data, m_time, m_step;
		hid_t m_mspace, m_mspace1;   // reused [1 x nIds x nComp] and [1] memory spaces
		hsize_t m_t;
		size_t m_n_ids, m_n_comp;
	};

	/*
	EnvelopeSink: componentwise running extremes, per MODEL_STAGE.

	Accumulators are all shape [nIds x nComp], row-major:
	  m_min     : running minimum per component
	  m_max     : running maximum per component
	  m_absmax  : running max(|value|) per component
	  m_arg_step: commitTag at which m_absmax for that component was last set
	Init happens on the FIRST accept(): the accumulators are seeded from the
	first buffer, so MIN/MAX start at the real first sample, not +/-inf.

	  accept():   element-wise update of all four accumulators from the step
	              buffer; first call seeds them.
	  flush():    overwrite the ENVELOPES datasets in place for crash safety;
	              callable any time after the first accept.
	  finalize(): same as flush(); the final authoritative write for the stage.
	  reset():    drop accumulators (cheap) so the sink can be reused on the next
	              MODEL_STAGE.
	*/
	class EnvelopeSink : public ResultSink {
	public:
		explicit EnvelopeSink(ResultFamily::Enum family)
			: m_family(family), m_seeded(false), m_n_ids(0), m_n_comp(0),
			  m_result_type(0), m_data_type(0), m_created(false),
			  m_dmin(HID_INVALID), m_dmax(HID_INVALID), m_dabs(HID_INVALID),
			  m_darg(HID_INVALID) {}
		// closes the MIN/MAX/ABSMAX/ARG_STEP handles held for the stage.
		~EnvelopeSink() override;
		// true once the ENVELOPES/<name> group exists (the recorder writes the
		// element COLUMN_MAP once, after the first write).
		bool written() const { return m_created; }
		// bytes rewritten per envelope write (MIN/MAX/ABSMAX f64 + ARG_STEP i32).
		size_t payloadBytes() const { return m_n_ids * m_n_comp * 28; }

		void begin(detail::ProcessInfo& info, const ResultSource& src) override;
		void accept(detail::ProcessInfo& info, const ResultSource& src,
		            const std::vector<double>& buffer) override;
		void finalize(detail::ProcessInfo& info) override;

		// Periodic crash-safety rewrite of the envelope datasets. Identical
		// payload to finalize(); separated so the recorder can call it on its
		// flush cadence without the meaning of "stage done".
		void flush(detail::ProcessInfo& info);

		// Drop all accumulators for a new MODEL_STAGE. Cheap (vector clears).
		void reset();

	private:
		// shared writer for flush()/finalize()
		void writeEnvelope(detail::ProcessInfo& info);

	private:
		ResultFamily::Enum m_family;
		bool m_seeded;            // accumulators initialized from first accept
		size_t m_n_ids;
		size_t m_n_comp;
		std::string m_name;       // result name (schema().name); cached at begin/first accept
		// Self-describing result-group metadata, cached so the ENVELOPES/<name>
		// group carries the same COMPONENTS/DISPLAY_NAME/DIMENSION as the
		// time-series result group.
		std::string m_display_name, m_components_csv, m_dimension, m_description;
		int m_result_type, m_data_type;
		std::string m_partition_reduction;   // PARTITION_REDUCTION attr ("NONE"/"SUM"/"UNSUPPORTED")
		std::vector<int> m_ids;   // cached ids() snapshot for the ID dataset
		std::vector<double> m_min;
		std::vector<double> m_max;
		std::vector<double> m_absmax;
		std::vector<int> m_arg_step; // per-component step index of the abs-extreme
		// the group + datasets are created on the first write and then
		// OVERWRITTEN in place (H5Dwrite on held handles).
		bool m_created;
		hid_t m_dmin, m_dmax, m_dabs, m_darg;
	};

} // namespace ladruno

#endif // Ladruno_Sinks_h
