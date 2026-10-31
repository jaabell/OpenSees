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
// Description: nodal result sources of the ladruno recorder (kinematics,
// reactions, unbalanced and Rayleigh forces, pressure, modes of vibration,
// sensitivities).
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#ifndef Ladruno_NodeResults_h
#define Ladruno_NodeResults_h

#include "Ladruno_ResultIO.h"  // ladruno::ResultSource, ladruno::ResultSchema
#include "Ladruno_Types.h"     // ladruno::detail::ProcessInfo, enums, utils

// OpenSees (only what the node sources need; NO h5 layer here)
#include "Node.h"
#include "NodeIter.h"
#include "Domain.h"
#include "Vector.h"

#include <string>
#include <vector>

namespace ladruno {

	/*
	NodeResultSource: common base for all nodal sources.

	Holds the schema, the spatial dimension (m_ndim) and the row-identity vector
	of node tags. The node set is gathered with a NodeIter over the ProcessInfo
	domain at construction (convenience ctor) or supplied explicitly.

	Subclasses fill in the schema in their ctor and override evaluate().
	*/
	class NodeResultSource : public ResultSource {
	public:
		// Gather all nodes from the domain (default node set).
		explicit NodeResultSource(const detail::ProcessInfo& info);
		// Or accept an explicit node-tag set (mirrors a region-restricted recorder).
		NodeResultSource(const detail::ProcessInfo& info, const std::vector<int>& node_ids);

		virtual ~NodeResultSource() = default;

		const ResultSchema& schema() const override { return m_schema; }
		const std::vector<int>& ids() const override { return m_ids; }

	protected:
		// Resolve a node tag -> Node* via the domain held in `info`.
		Node* getNode(const detail::ProcessInfo& info, int tag) const;
		// Node* of row i, resolved once per source. `getNode(tag)` is a std::map
		// walk (Domain's MapOfTaggedObjects) per node per channel per step. The
		// cache is valid for the source's lifetime: sources are rebuilt on every
		// domain-change stamp (any node add/remove bumps it) before the next
		// evaluate(), and a failed rebuild releases them.
		Node* nodeAt(const detail::ProcessInfo& info, size_t i) const;

	protected:
		int m_ndim;                 // info.num_dimensions
		std::vector<int> m_ids;     // node tags this source covers
		mutable std::vector<Node*> m_nodes;   // lazily resolved, m_ids order
		mutable bool m_nodes_resolved = false;
		ResultSchema m_schema;      // constant for the lifetime of the source
	};

	/* ===================================================================== */
	/* Kinematics: consistent quantities, requiresPartitionReduction()==false */
	/* ===================================================================== */

	class DisplacementSource : public NodeResultSource {
	public:
		explicit DisplacementSource(const detail::ProcessInfo& info);
		DisplacementSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class RotationSource : public NodeResultSource {
	public:
		explicit RotationSource(const detail::ProcessInfo& info);
		RotationSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class VelocitySource : public NodeResultSource {
	public:
		explicit VelocitySource(const detail::ProcessInfo& info);
		VelocitySource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class AngularVelocitySource : public NodeResultSource {
	public:
		explicit AngularVelocitySource(const detail::ProcessInfo& info);
		AngularVelocitySource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class AccelerationSource : public NodeResultSource {
	public:
		explicit AccelerationSource(const detail::ProcessInfo& info);
		AccelerationSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class AngularAccelerationSource : public NodeResultSource {
	public:
		explicit AngularAccelerationSource(const detail::ProcessInfo& info);
		AngularAccelerationSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	// PRESSURE: the u-p elements (FourNodeQuadUP, brickUP, SSP*UP, ...) carry
	// DOF ndm = integral of p dt, so the physical pressure is the velocity slot
	// ndm, as in MPCORecorder.
	class PressureSource : public NodeResultSource {
	public:
		explicit PressureSource(const detail::ProcessInfo& info);
		PressureSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	/* ===================================================================== */
	/* Reactions: additive at partition-boundary nodes.                     */
	/* requiresPartitionReduction()==true: the reaction at a node shared      */
	/* between partitions is the SUM of each partition's contribution, so a   */
	/* per-step partition reduction must precede sink accumulation. The flag  */
	/* returned by reactionFlag() selects which reaction the Node exposes     */
	/* (0=reaction, 1=incl. inertia, 2=Rayleigh); the recorder driver primes  */
	/* the Domain reactions with it before evaluate().                        */
	/* ===================================================================== */

	class ReactionForceSource : public NodeResultSource {
	public:
		explicit ReactionForceSource(const detail::ProcessInfo& info);
		ReactionForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		// reactions are additive across partition boundaries.
		bool requiresPartitionReduction() const override { return true; }
		// 0=reaction, 1=incl. inertia, 2=Rayleigh.
		virtual int reactionFlag() const { return 0; }
	protected:
		void build(const char* group, const char* display, const char* description);
	};

	class ReactionMomentSource : public NodeResultSource {
	public:
		explicit ReactionMomentSource(const detail::ProcessInfo& info);
		ReactionMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		bool requiresPartitionReduction() const override { return true; }
		virtual int reactionFlag() const { return 0; }
	protected:
		void build(const char* group, const char* display, const char* description);
	};

	// Inertia / Rayleigh variants only change schema name + reaction flag; the
	// numerical getReaction() path is identical (the Node returns the variant the
	// recorder driver primed via setReactionFlag-equivalent before record()).
	class ReactionForceIncInertiaSource : public ReactionForceSource {
	public:
		explicit ReactionForceIncInertiaSource(const detail::ProcessInfo& info);
		ReactionForceIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		int reactionFlag() const override { return 1; }
	};

	class ReactionMomentIncInertiaSource : public ReactionMomentSource {
	public:
		explicit ReactionMomentIncInertiaSource(const detail::ProcessInfo& info);
		ReactionMomentIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		int reactionFlag() const override { return 1; }
	};

	class RayleighForceSource : public ReactionForceSource {
	public:
		explicit RayleighForceSource(const detail::ProcessInfo& info);
		RayleighForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		int reactionFlag() const override { return 2; }
	};

	class RayleighMomentSource : public ReactionMomentSource {
	public:
		explicit RayleighMomentSource(const detail::ProcessInfo& info);
		RayleighMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		int reactionFlag() const override { return 2; }
	};

	/* ===================================================================== */
	/* Unbalanced forces/moments: also additive at partition boundaries.    */
	/* requiresPartitionReduction()==true.                                   */
	/* ===================================================================== */

	class UnbalancedForceSource : public NodeResultSource {
	public:
		explicit UnbalancedForceSource(const detail::ProcessInfo& info);
		UnbalancedForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		// unbalanced load is additive across partition boundaries.
		bool requiresPartitionReduction() const override { return true; }
	protected:
		void build(const char* group, const char* display, const char* description);
	};

	class UnbalancedMomentSource : public NodeResultSource {
	public:
		explicit UnbalancedMomentSource(const detail::ProcessInfo& info);
		UnbalancedMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		bool requiresPartitionReduction() const override { return true; }
	protected:
		void build(const char* group, const char* display, const char* description);
	};

	// IncInertia variants override evaluate() to use getUnbalancedLoadIncInertia().
	class UnbalancedForceIncInertiaSource : public UnbalancedForceSource {
	public:
		explicit UnbalancedForceIncInertiaSource(const detail::ProcessInfo& info);
		UnbalancedForceIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	};

	class UnbalancedMomentIncInertiaSource : public UnbalancedMomentSource {
	public:
		explicit UnbalancedMomentIncInertiaSource(const detail::ProcessInfo& info);
		UnbalancedMomentIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	};

	/* ===================================================================== */
	/* Modes of vibration: several modes per recorded step.                  */
	/*                                                                       */
	/* The ResultSource contract yields ONE flat buffer per evaluate(), so   */
	/* the source exposes:                                                   */
	/*   - evaluate(info, buf): the current mode (set via setCurrentMode())   */
	/*   - numModes(info), setCurrentMode(k), modeInfo(info, k, ...)          */
	/* The recorder drives the mode loop: for k in [0,numModes):             */
	/* setCurrentMode(k); evaluate(...); append one slab + the mode scalars. */
	/* requiresPartitionReduction()==false (eigenvectors are consistent).     */
	/* Recording is also gated on info.record_eigen_on_this_step by the sink. */
	/* ===================================================================== */

	class ModesOfVibrationSource : public NodeResultSource {
	public:
		explicit ModesOfVibrationSource(const detail::ProcessInfo& info);
		ModesOfVibrationSource(const detail::ProcessInfo& info, const std::vector<int>& ids);

		// Fill buffer for the CURRENT mode (translational components).
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;

		// Mode-aware API for the mode-loop driver (StreamingSink).
		int numModes(const detail::ProcessInfo& info) const;  // captured spectrum size
		void setCurrentMode(int k) { m_current_mode = k; }
		int currentMode() const { return m_current_mode; }
		// Per-mode scalars (lambda, omega = sqrt(lambda), f, T).
		void modeInfo(const detail::ProcessInfo& info, int k,
		              double& lambda, double& omega,
		              double& freq, double& period) const;

	protected:
		void build();
		int m_current_mode;
	};

	class ModesOfVibrationRotationalSource : public ModesOfVibrationSource {
	public:
		explicit ModesOfVibrationRotationalSource(const detail::ProcessInfo& info);
		ModesOfVibrationRotationalSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	/* ===================================================================== */
	/* Sensitivity sources: consistent, requiresPartitionReduction()==false. */
	/* Each carries a gradient index (m_grad); evaluate() guards m_grad        */
	/* against domain->getNumParameters().                                     */
	/* ===================================================================== */

	class DisplacementSensitivitySource : public NodeResultSource {
	public:
		DisplacementSensitivitySource(const detail::ProcessInfo& info, int grad);
		DisplacementSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class RotationSensitivitySource : public NodeResultSource {
	public:
		RotationSensitivitySource(const detail::ProcessInfo& info, int grad);
		RotationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class VelocitySensitivitySource : public NodeResultSource {
	public:
		VelocitySensitivitySource(const detail::ProcessInfo& info, int grad);
		VelocitySensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class AngularVelocitySensitivitySource : public NodeResultSource {
	public:
		AngularVelocitySensitivitySource(const detail::ProcessInfo& info, int grad);
		AngularVelocitySensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class AccelerationSensitivitySource : public NodeResultSource {
	public:
		AccelerationSensitivitySource(const detail::ProcessInfo& info, int grad);
		AccelerationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class AngularAccelerationSensitivitySource : public NodeResultSource {
	public:
		AngularAccelerationSensitivitySource(const detail::ProcessInfo& info, int grad);
		AngularAccelerationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

} // namespace ladruno

#endif // Ladruno_NodeResults_h
