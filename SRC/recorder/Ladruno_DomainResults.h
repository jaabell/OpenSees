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
// Description: whole-model and per-region result sources of the ladruno
// recorder. EnergyBalanceSource evaluates the structural-dynamics energy
// balance (KE, IE, DW, ULW, RES, ERR) through the shared kernel of
// EnergyBalanceKernel.h, so its values equal those of EnergyBalanceRecorder.
//
// References:
//   T. Belytschko, W.K. Liu, B. Moran (2000). Nonlinear Finite Elements for
//     Continua and Structures, Wiley, Section 6.2.3 (energy balance check).
//   C.-M. Uang, V.V. Bertero (1990). Evaluation of seismic energy in
//     structures. Earthquake Engineering and Structural Dynamics 19(1),
//     77-90. doi:10.1002/eqe.4290190108
//

#ifndef Ladruno_DomainResults_h
#define Ladruno_DomainResults_h

#include "Ladruno_ResultIO.h"      // ladruno::ResultSource, ladruno::ResultSchema
#include "Ladruno_Types.h"         // ladruno::detail::ProcessInfo, enums

#include "EnergyBalanceKernel.h"   // ebkernel::EnergyAccumulator + per-entity terms

#include "Vector.h"

#include <string>
#include <vector>

namespace ladruno {

	/*
	EnergyBalanceSource: whole-model OR per-region energy balance.

	Two construction modes:
	  - Whole-model: ids() == {0}; evaluate() runs one accumulator over the
	    whole domain -> 6 doubles.
	  - Per-region : ids() == [regionTags]; evaluate() sweeps each region (in
	    the given tag order), stepping its own accumulator -> 6 doubles per
	    region, row-major [nRegions x 6].

	Unlike the stateless node/element sources, this source carries cross-step
	state (the trapezoidal work integrals and the closure live in the
	accumulators, plus the previous time stamp).

	That state lives in an EnergyState owned by the recorder, so it survives
	the source being rebuilt at a MODEL_STAGE change (a stage-local state would
	restart at zero and take a rate x t_total jump on the first step of every
	stage), and advance(), idempotent per commit tag, is called by the
	recorder on every commit, so the work integrals are not sampled only at
	the -T recorded steps. evaluate() = advance() + copy the latest values.
	*/
	struct EnergyState {
		ebkernel::EnergyAccumulator model;                // whole-model
		std::vector<ebkernel::EnergyAccumulator> regions; // per-region, tag order
		double prev_time = 0.0;      // time of the last integrated commit
		bool first = true;           // the first integrated commit seeds the rates
		bool has_last = false;
		int last_tag = 0;            // commit tag last integrated (idempotency key)
		std::vector<double> last_out;   // [nIds x 6] of the last integrated commit
	};

	class EnergyBalanceSource : public ResultSource {
	public:
		// Whole-model (ON_DOMAIN). ids() == {0}.
		// `state` (optional) is recorder-owned and outlives the source; without
		// it the source owns a private, stage-local one.
		explicit EnergyBalanceSource(const detail::ProcessInfo& info, EnergyState* state = 0);
		// Per-region (ON_REGIONS). ids() == region_tags.
		EnergyBalanceSource(const detail::ProcessInfo& info,
		                    const std::vector<int>& region_tags, EnergyState* state = 0);

		~EnergyBalanceSource() override = default;

		const ResultSchema& schema() const override { return m_schema; }
		const std::vector<int>& ids() const override { return m_ids; }

		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		// Integrate this commit (no-op if this commit tag was already
		// integrated). The recorder calls it on every commit.
		void advance(const detail::ProcessInfo& info);

		// Energy is an additive/global quantity that needs a per-step partition
		// reduction before a sink may accumulate it. The values are correct for
		// a serial run; partitioned runs write per-partition partials.
		bool requiresPartitionReduction() const override { return true; }
		// Not "SUM": KE at a node shared by two partitions is counted in both,
		// and RES/ERR are derived from the other components, so a componentwise
		// sum of the partition files is wrong.
		const char* partitionReduction() const override { return "UNSUPPORTED"; }

	private:
		void buildSchema();

	private:
		bool m_per_region;                 // false => whole-model {0}
		std::vector<int> m_ids;            // {0} or the region tags
		ResultSchema m_schema;

		// cross-step state (recorder-owned when given)
		EnergyState m_own;                 // used only when no state was passed
		EnergyState* m_state;              // -> recorder's state, or &m_own
	};

} // namespace ladruno

#endif // Ladruno_DomainResults_h
