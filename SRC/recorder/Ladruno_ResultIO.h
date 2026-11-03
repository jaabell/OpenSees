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
// Description: the result-source / result-sink contract of the ladruno
// recorder. A ResultSource computes one result for the current step (nodal,
// element, whole-model or per-region); a ResultSink persists it (streaming
// time series or running envelope). Any source composes with any sink.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#ifndef Ladruno_ResultIO_h
#define Ladruno_ResultIO_h

#include "Ladruno_Types.h"

#include <string>
#include <vector>

namespace ladruno {

	// Self-describing column metadata for one result.
	struct ResultSchema {
		std::string name;            // on-disk group name (e.g. "DISPLACEMENT", "stress")
		std::string display_name;    // human label
		std::string components_csv;  // "Ux,Uy,Uz"
		std::string dimension;       // "L", "F", "F*L", ...
		std::string description;
		int num_components = 0;
		detail::ResultType::Enum result_type = detail::ResultType::Generic;
		detail::ResultDataType::Enum data_type = detail::ResultDataType::Scalar;
	};

	/*
	ResultSource: what to compute. A node / element / domain result implements
	this. The recorder pairs each source with a sink. A source yields, for the
	current step, a flat row-major buffer of shape [ ids().size() * schema.num_components ]
	(domain sources use a single synthetic id {0}). Node and element sources are
	stateless across steps; cross-step behavior (streaming vs envelope) is the
	sink's job. Exception: EnergyBalanceSource integrates work over time; its
	state is owned by the recorder and advanced on every commit.
	*/
	class ResultSource {
	public:
		virtual ~ResultSource() = default;

		// Constant for the lifetime of the source.
		virtual const ResultSchema& schema() const = 0;

		// Row identity: node tags / element-gp ids / {0} for domain scalars.
		virtual const std::vector<int>& ids() const = 0;

		// Fill `buffer` (resized by the source) with this step's values,
		// row-major over ids() x num_components.
		virtual void evaluate(const detail::ProcessInfo& info,
		                      std::vector<double>& buffer) = 0;

		// Parallel routing: true for additive/global quantities that need a
		// per-step partition reduction before a sink may accumulate them (e.g.
		// boundary-node reactions, base shear). Consistent quantities (kinematics,
		// element results) return false and stream/envelope locally with zero comm.
		virtual bool requiresPartitionReduction() const { return false; }

		// How a READER must combine this result's rows for an id present in several
		// partition files. Written verbatim to each result group's
		// PARTITION_REDUCTION attribute:
		//   "NONE"        consistent across partitions: any one copy is the value;
		//   "SUM"         additive: each partition holds the partial from its own
		//                 elements (e.g. boundary-node reactions); sum the copies;
		//   "UNSUPPORTED" partial per partition and NOT recoverable by a componentwise
		//                 sum (e.g. energy: shared-node KE counted twice, RES/ERR
		//                 derived); a reader must not present a merged value.
		// Also gates -envelope in partitioned runs: anything but "NONE" is refused there,
		// because a per-partition extreme of a partial cannot be recombined.
		virtual const char* partitionReduction() const
		{
			return requiresPartitionReduction() ? "SUM" : "NONE";
		}
	};

	/*
	ResultSink: how to persist. Wraps any source.
	  - StreamingSink appends one slab per recorded step to the chunked
	    RESULTS/.../DATA[T x nIds x nComp] plus TIME[T] / STEP[T],
	    holding the dataset handles open for the stage.
	  - EnvelopeSink keeps running min/max/absmax + arg-step and rewrites them in
	    place on the recorder's cadence and at stage end, per MODEL_STAGE.
	The sink owns the HDF5 group lifetime for its result. One sink lives for one
	MODEL_STAGE: the recorder deletes it at the stage change (after finalizing
	envelopes) and builds a fresh one.
	*/
	class ResultSink {
	public:
		virtual ~ResultSink() = default;

		// Create groups + the ID dataset for this source's result. The recorder
		// does not call this directly: both sinks begin lazily on their first
		// accept() (the sink is new for every MODEL_STAGE).
		virtual void begin(detail::ProcessInfo& info, const ResultSource& src) = 0;

		// Persist one step's buffer. NOT partition-reduced: each rank writes its
		// own part file and the reader combines partitions per the result's
		// PARTITION_REDUCTION attribute.
		virtual void accept(detail::ProcessInfo& info, const ResultSource& src,
		                    const std::vector<double>& buffer) = 0;

		// Flush any deferred datasets (envelope writes its accumulators here).
		// Called on stage change, periodically, and at recorder teardown.
		virtual void finalize(detail::ProcessInfo& info) = 0;
	};

} // namespace ladruno

#endif // Ladruno_ResultIO_h
