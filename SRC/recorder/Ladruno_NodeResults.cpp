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
// Description: nodal result sources of the ladruno recorder.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#include "Ladruno_NodeResults.h"

#include "elementAPI.h"
#include "OPS_Globals.h"   // opserr
#include "Matrix.h"        // Node::getEigenvectors() returns const Matrix&

#include "ID.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ladruno {

	// Result-group name: MODEL_STAGE[<id>]/RESULTS/ON_NODES/<group>
	static std::string makeResultName(const detail::ProcessInfo& info, const std::string& group)
	{
		std::stringstream ss;
		ss << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/" << group;
		return ss.str();
	}

	/* ===================================================================== */
	/* NodeResultSource base                                                 */
	/* ===================================================================== */

	NodeResultSource::NodeResultSource(const detail::ProcessInfo& info)
		: m_ndim(info.num_dimensions)
	{
		// Gather all nodes from the domain (default node set, matching the
		// NodeIter walk of LadrunoRecorder::writeModel()).
		if (info.domain != 0) {
			NodeIter& nit = info.domain->getNodes();
			Node* nptr = 0;
			while ((nptr = nit()) != 0)
				m_ids.push_back(nptr->getTag());
		}
	}

	NodeResultSource::NodeResultSource(const detail::ProcessInfo& info, const std::vector<int>& node_ids)
		: m_ndim(info.num_dimensions)
		, m_ids(node_ids)
	{
	}

	Node* NodeResultSource::getNode(const detail::ProcessInfo& info, int tag) const
	{
		if (info.domain == 0)
			return 0;
		return info.domain->getNode(tag);
	}

	Node* NodeResultSource::nodeAt(const detail::ProcessInfo& info, size_t i) const
	{
		if (!m_nodes_resolved) {
			m_nodes.resize(m_ids.size());
			for (size_t k = 0; k < m_ids.size(); ++k)
				m_nodes[k] = getNode(info, m_ids[k]);
			m_nodes_resolved = true;
		}
		return (i < m_nodes.size()) ? m_nodes[i] : 0;
	}

	/* ===================================================================== */
	/* DisplacementSource                                                    */
	/* ===================================================================== */

	DisplacementSource::DisplacementSource(const detail::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	DisplacementSource::DisplacementSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void DisplacementSource::build()
	{
		m_schema.name = "DISPLACEMENT";
		m_schema.display_name = "Displacement";
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Ux";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Ux,Uy";
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Ux,Uy,Uz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L";
		m_schema.description = "Nodal displacement field";
		m_schema.result_type = detail::ResultType::Generic;
	}

	void DisplacementSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getTrialDisp(), buffer);
		}
	}

	/* ===================================================================== */
	/* RotationSource                                                        */
	/* ===================================================================== */

	RotationSource::RotationSource(const detail::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	RotationSource::RotationSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void RotationSource::build()
	{
		m_schema.name = "ROTATION";
		m_schema.display_name = "Rotation";
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "Rz";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "Rx,Ry,Rz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "A";
		m_schema.description = "Nodal rotation field";
		m_schema.result_type = detail::ResultType::Generic;
	}

	void RotationSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getTrialDisp(), buffer);
		}
	}

	/* ===================================================================== */
	/* VelocitySource                                                        */
	/* ===================================================================== */

	VelocitySource::VelocitySource(const detail::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	VelocitySource::VelocitySource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void VelocitySource::build()
	{
		m_schema.name = "VELOCITY";
		m_schema.display_name = "Velocity";
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Vx";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Vx,Vy";
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Vx,Vy,Vz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L/t";
		m_schema.description = "Nodal velocity field";
		m_schema.result_type = detail::ResultType::Generic;
	}

	void VelocitySource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getTrialVel(), buffer);
		}
	}

	/* ===================================================================== */
	/* AngularVelocitySource                                                 */
	/* ===================================================================== */

	AngularVelocitySource::AngularVelocitySource(const detail::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	AngularVelocitySource::AngularVelocitySource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void AngularVelocitySource::build()
	{
		m_schema.name = "ANGULAR_VELOCITY";
		m_schema.display_name = "Angular Velocity";
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "RVz";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "RVx,RVy,RVz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L/t*A";
		m_schema.description = "Nodal angular velocity field";
		m_schema.result_type = detail::ResultType::Generic;
	}

	void AngularVelocitySource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getTrialVel(), buffer);
		}
	}

	/* ===================================================================== */
	/* AccelerationSource                                                    */
	/* ===================================================================== */

	AccelerationSource::AccelerationSource(const detail::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	AccelerationSource::AccelerationSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void AccelerationSource::build()
	{
		m_schema.name = "ACCELERATION";
		m_schema.display_name = "Acceleration";
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Ax";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Ax,Ay";
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Ax,Ay,Az";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L/t^2";
		m_schema.description = "Nodal acceleration field";
		m_schema.result_type = detail::ResultType::Generic;
	}

	void AccelerationSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getTrialAccel(), buffer);
		}
	}

	/* ===================================================================== */
	/* AngularAccelerationSource                                             */
	/* ===================================================================== */

	AngularAccelerationSource::AngularAccelerationSource(const detail::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	AngularAccelerationSource::AngularAccelerationSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void AngularAccelerationSource::build()
	{
		m_schema.name = "ANGULAR_ACCELERATION";
		m_schema.display_name = "Angular Acceleration";
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "RAz";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "RAx,RAy,RAz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L/t^2*A";
		m_schema.description = "Nodal angular acceleration field";
		m_schema.result_type = detail::ResultType::Generic;
	}

	void AngularAccelerationSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getTrialAccel(), buffer);
		}
	}

	/* ===================================================================== */
	/* PressureSource                                                        */
	/* ===================================================================== */

	PressureSource::PressureSource(const detail::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	PressureSource::PressureSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void PressureSource::build()
	{
		m_schema.name = "PRESSURE";
		m_schema.display_name = "Pressure";
		m_schema.components_csv = "p";
		m_schema.num_components = 1;
		m_schema.data_type = detail::ResultDataType::Scalar;
		m_schema.dimension = "F/L^2";
		m_schema.description = "Nodal pressure field";
		m_schema.result_type = detail::ResultType::Generic;
	}

	void PressureSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		// Pressure is the extra DOF: idx 2 (2D u-p) or idx 3 (3D u-p), read from
		// the velocity vector (the u-p elements carry the integral of p dt).
		for (size_t i = 0; i < m_ids.size(); i++) {
			double pressure = 0.0;
			Node* inode = nodeAt(info, i);
			if (inode) {
				const Vector& resp = inode->getTrialVel();
				if (inode->getCrds().Size() == 2) {
					if (resp.Size() == 3) {
						pressure = resp[2];
					}
				}
				else if (inode->getCrds().Size() == 3) {
					if (resp.Size() == 4) {
						pressure = resp[3];
					}
				}
			}
			buffer[i] = pressure;
		}
	}

	/* ===================================================================== */
	/* ReactionForceSource (+ inertia/Rayleigh variants)                     */
	/* ===================================================================== */

	ReactionForceSource::ReactionForceSource(const detail::ProcessInfo& info)
		: NodeResultSource(info)
	{ build("REACTION_FORCE", "Reaction Force", "Nodal reaction force field"); }
	ReactionForceSource::ReactionForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids)
	{ build("REACTION_FORCE", "Reaction Force", "Nodal reaction force field"); }

	void ReactionForceSource::build(const char* group, const char* display, const char* description)
	{
		m_schema.name = group;
		m_schema.display_name = display;
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "RFx";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "RFx,RFy";
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "RFx,RFy,RFz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "F";
		m_schema.description = description;
		m_schema.result_type = detail::ResultType::Generic;
	}

	void ReactionForceSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getReaction(), buffer);
		}
	}

	ReactionForceIncInertiaSource::ReactionForceIncInertiaSource(const detail::ProcessInfo& info)
		: ReactionForceSource(info)
	{ build("REACTION_FORCE_INCLUDING_INERTIA", "Reaction Force Including Inertia",
	        "Nodal reaction force field including inertia"); }
	ReactionForceIncInertiaSource::ReactionForceIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: ReactionForceSource(info, ids)
	{ build("REACTION_FORCE_INCLUDING_INERTIA", "Reaction Force Including Inertia",
	        "Nodal reaction force field including inertia"); }

	RayleighForceSource::RayleighForceSource(const detail::ProcessInfo& info)
		: ReactionForceSource(info)
	{ build("RAYLEIGH_FORCE", "Rayleigh Force", "Nodal Rayleigh force field"); }
	RayleighForceSource::RayleighForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: ReactionForceSource(info, ids)
	{ build("RAYLEIGH_FORCE", "Rayleigh Force", "Nodal Rayleigh force field"); }

	/* ===================================================================== */
	/* ReactionMomentSource (+ inertia/Rayleigh variants)                    */
	/* ===================================================================== */

	ReactionMomentSource::ReactionMomentSource(const detail::ProcessInfo& info)
		: NodeResultSource(info)
	{ build("REACTION_MOMENT", "Reaction Moment", "Nodal reaction moment field"); }
	ReactionMomentSource::ReactionMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids)
	{ build("REACTION_MOMENT", "Reaction Moment", "Nodal reaction moment field"); }

	void ReactionMomentSource::build(const char* group, const char* display, const char* description)
	{
		m_schema.name = group;
		m_schema.display_name = display;
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "RMz";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "RMx,RMy,RMz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "F*L";
		m_schema.description = description;
		m_schema.result_type = detail::ResultType::Generic;
	}

	void ReactionMomentSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getReaction(), buffer);
		}
	}

	ReactionMomentIncInertiaSource::ReactionMomentIncInertiaSource(const detail::ProcessInfo& info)
		: ReactionMomentSource(info)
	{ build("REACTION_MOMENT_INCLUDING_INERTIA", "Reaction Moment Including Inertia",
	        "Nodal reaction moment field including inertia"); }
	ReactionMomentIncInertiaSource::ReactionMomentIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: ReactionMomentSource(info, ids)
	{ build("REACTION_MOMENT_INCLUDING_INERTIA", "Reaction Moment Including Inertia",
	        "Nodal reaction moment field including inertia"); }

	RayleighMomentSource::RayleighMomentSource(const detail::ProcessInfo& info)
		: ReactionMomentSource(info)
	{ build("RAYLEIGH_MOMENT", "Rayleigh Moment", "Nodal Rayleigh moment field"); }
	RayleighMomentSource::RayleighMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: ReactionMomentSource(info, ids)
	{ build("RAYLEIGH_MOMENT", "Rayleigh Moment", "Nodal Rayleigh moment field"); }

	/* ===================================================================== */
	/* UnbalancedForceSource (+ IncInertia variant)                          */
	/* ===================================================================== */

	UnbalancedForceSource::UnbalancedForceSource(const detail::ProcessInfo& info)
		: NodeResultSource(info)
	{ build("UNBALANCED_FORCE", "Unbalanced Force", "Nodal unbalanced force field"); }
	UnbalancedForceSource::UnbalancedForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids)
	{ build("UNBALANCED_FORCE", "Unbalanced Force", "Nodal unbalanced force field"); }

	void UnbalancedForceSource::build(const char* group, const char* display, const char* description)
	{
		m_schema.name = group;
		m_schema.display_name = display;
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Fx";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Fx,Fy";
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Fx,Fy,Fz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "F";
		m_schema.description = description;
		m_schema.result_type = detail::ResultType::Generic;
	}

	void UnbalancedForceSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getUnbalancedLoad(), buffer);
		}
	}

	UnbalancedForceIncInertiaSource::UnbalancedForceIncInertiaSource(const detail::ProcessInfo& info)
		: UnbalancedForceSource(info)
	{ build("UNBALANCED_FORCE_INCLUDING_INERTIA", "Unbalanced Force Including Inertia",
	        "Nodal unbalanced force field including inertia"); }
	UnbalancedForceIncInertiaSource::UnbalancedForceIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: UnbalancedForceSource(info, ids)
	{ build("UNBALANCED_FORCE_INCLUDING_INERTIA", "Unbalanced Force Including Inertia",
	        "Nodal unbalanced force field including inertia"); }

	void UnbalancedForceIncInertiaSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getUnbalancedLoadIncInertia(), buffer);
		}
	}

	/* ===================================================================== */
	/* UnbalancedMomentSource (+ IncInertia variant)                         */
	/* ===================================================================== */

	UnbalancedMomentSource::UnbalancedMomentSource(const detail::ProcessInfo& info)
		: NodeResultSource(info)
	{ build("UNBALANCED_MOMENT", "Unbalanced Moment", "Nodal unbalanced moment field"); }
	UnbalancedMomentSource::UnbalancedMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids)
	{ build("UNBALANCED_MOMENT", "Unbalanced Moment", "Nodal unbalanced moment field"); }

	void UnbalancedMomentSource::build(const char* group, const char* display, const char* description)
	{
		m_schema.name = group;
		m_schema.display_name = display;
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "Mz";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "Mx,My,Mz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "F*L";
		m_schema.description = description;
		m_schema.result_type = detail::ResultType::Generic;
	}

	void UnbalancedMomentSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getUnbalancedLoad(), buffer);
		}
	}

	UnbalancedMomentIncInertiaSource::UnbalancedMomentIncInertiaSource(const detail::ProcessInfo& info)
		: UnbalancedMomentSource(info)
	{ build("UNBALANCED_MOMENT_INCLUDING_INERTIA", "Unbalanced Moment Including Inertia",
	        "Nodal unbalanced moment field including inertia"); }
	UnbalancedMomentIncInertiaSource::UnbalancedMomentIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: UnbalancedMomentSource(info, ids)
	{ build("UNBALANCED_MOMENT_INCLUDING_INERTIA", "Unbalanced Moment Including Inertia",
	        "Nodal unbalanced moment field including inertia"); }

	void UnbalancedMomentIncInertiaSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getUnbalancedLoadIncInertia(), buffer);
		}
	}

	/* ===================================================================== */
	/* ModesOfVibrationSource (translational): SPECIAL multi-mode case.     */
	/* ===================================================================== */

	ModesOfVibrationSource::ModesOfVibrationSource(const detail::ProcessInfo& info)
		: NodeResultSource(info), m_current_mode(0) { build(); }
	ModesOfVibrationSource::ModesOfVibrationSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids), m_current_mode(0) { build(); }

	void ModesOfVibrationSource::build()
	{
		m_schema.name = "MODES_OF_VIBRATION(U)";
		m_schema.display_name = "Modes of Vibration (U)";
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Ux";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Ux,Uy";
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Ux,Uy,Uz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = "Eigenvector (translational components)";
		m_schema.result_type = detail::ResultType::Modal;
	}

	int ModesOfVibrationSource::numModes(const detail::ProcessInfo& info) const
	{
		// The spectrum the recorder captured from the domain
		// (recordResultsOnNodes sets eigen_last_values from getEigenvalues()), not
		// the interpreter's numEigen, which `wipe` never resets.
		return info.eigen_last_values.Size();
	}

	void ModesOfVibrationSource::modeInfo(const detail::ProcessInfo& info, int k,
	                                      double& lambda, double& omega,
	                                      double& freq, double& period) const
	{
		// Per-mode scalars.
		lambda = k < info.eigen_last_values.Size() ? info.eigen_last_values[k] : 0.0;
		omega = std::sqrt(lambda);
		freq = omega / (2.0 * M_PI);
		period = 1.0 / freq;
	}

	void ModesOfVibrationSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		// Translational eigenvector components for the CURRENT mode.
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (!n) continue;
			// Node::getEigenvectors() exits when unset (a node
			// added after `eigen`); probe first and leave such nodes zero.
			if (n->getNumEigenvectors() <= m_current_mode) continue;
			size_t j = i * m_ndim;
			const Matrix& inode_eigenvec = n->getEigenvectors();
			// A node with fewer DOFs than ndm (scalar node in a
			// 2D/3D model) has fewer rows; Matrix() is unchecked in release.
			const int nr = std::min(m_ndim, inode_eigenvec.noRows());
			for (int r = 0; r < nr; ++r)
				buffer[j + (size_t)r] = inode_eigenvec(r, m_current_mode);
		}
	}

	/* ===================================================================== */
	/* ModesOfVibrationRotationalSource                                      */
	/* ===================================================================== */

	ModesOfVibrationRotationalSource::ModesOfVibrationRotationalSource(const detail::ProcessInfo& info)
		: ModesOfVibrationSource(info) { build(); }
	ModesOfVibrationRotationalSource::ModesOfVibrationRotationalSource(const detail::ProcessInfo& info, const std::vector<int>& ids)
		: ModesOfVibrationSource(info, ids) { build(); }

	void ModesOfVibrationRotationalSource::build()
	{
		m_schema.name = "MODES_OF_VIBRATION(R)";
		m_schema.display_name = "Modes of Vibration (R)";
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "Rz";
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "Rx,Ry,Rz";
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = "Eigenvector (rotational components)";
		m_schema.result_type = detail::ResultType::Modal;
	}

	void ModesOfVibrationRotationalSource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		// Rotational eigenvector slice depends on the node's DOF count (rows>2 in
		// 2D, rows>5 in 3D); nodes without rotational DOFs stay zero.
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = nodeAt(info, i);
			if (!n) continue;
			if (n->getNumEigenvectors() <= m_current_mode) continue;  // node outside eigen
			const Matrix& inode_eigenvec = n->getEigenvectors();
			if (m_ndim == 2 && inode_eigenvec.noRows() > 2) {
				size_t j = i;
				buffer[j] = inode_eigenvec(2, m_current_mode);
			}
			else if (m_ndim == 3 && inode_eigenvec.noRows() > 5) {
				size_t j = i * 3;
				buffer[j] = inode_eigenvec(3, m_current_mode);
				buffer[j + 1] = inode_eigenvec(4, m_current_mode);
				buffer[j + 2] = inode_eigenvec(5, m_current_mode);
			}
		}
	}

	/* ===================================================================== */
	/* Sensitivity sources                                                   */
	/* ===================================================================== */

	DisplacementSensitivitySource::DisplacementSensitivitySource(const detail::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	DisplacementSensitivitySource::DisplacementSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void DisplacementSensitivitySource::build()
	{
		m_schema.name = LADRUNO_MAKE_STRING("DISPLACEMENT_SENSITIVITY_" << m_grad);
		m_schema.display_name = LADRUNO_MAKE_STRING("Displacement Sensitivity d_U/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Ux/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Ux/d_q" << m_grad << ",d_Uy/d_q" << m_grad);
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Ux/d_q" << m_grad << ",d_Uy/d_q" << m_grad << ",d_Uz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = LADRUNO_MAKE_STRING("Nodal displacement sensitivity field d_U/d_q" << m_grad);
		m_schema.result_type = detail::ResultType::Generic;
	}

	void DisplacementSensitivitySource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "LadrunoRecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = nodeAt(info, i);
			if (!inode) continue;
			size_t j = i * m_ndim;
			buffer[j] = inode->getDispSensitivity(1, m_grad - 1);
			if (m_ndim > 1) {
				buffer[j + 1] = inode->getDispSensitivity(2, m_grad - 1);
				if (m_ndim > 2) {
					buffer[j + 2] = inode->getDispSensitivity(3, m_grad - 1);
				}
			}
		}
	}

	RotationSensitivitySource::RotationSensitivitySource(const detail::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	RotationSensitivitySource::RotationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void RotationSensitivitySource::build()
	{
		m_schema.name = LADRUNO_MAKE_STRING("ROTATION_SENSITIVITY_" << m_grad);
		m_schema.display_name = LADRUNO_MAKE_STRING("Rotation Sensitivity d_R/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Rz/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Rx/d_q" << m_grad << ",d_Ry/d_q" << m_grad << ",d_Rz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = LADRUNO_MAKE_STRING("Nodal rotation sensitivity field d_R/d_q" << m_grad);
		m_schema.result_type = detail::ResultType::Generic;
	}

	void RotationSensitivitySource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "LadrunoRecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = nodeAt(info, i);
			if (!inode) continue;
			const Vector& iresponse = inode->getTrialDisp(); // aux, just to check for rotational components
			if (m_ndim == 2) {
				size_t j = i; // node_counter * 1 (1 rotation component in 2D)
				if (iresponse.Size() > 2)
					buffer[j] = inode->getDispSensitivity(3, m_grad - 1);
				else
					buffer[j] = 0.0;
			}
			else if (m_ndim == 3) {
				size_t j = i * 3; // node_counter * 3 (3 rotation components in 3D)
				if (iresponse.Size() > 5) {
					buffer[j] = inode->getDispSensitivity(4, m_grad - 1);
					buffer[j + 1] = inode->getDispSensitivity(5, m_grad - 1);
					buffer[j + 2] = inode->getDispSensitivity(6, m_grad - 1);
				}
				else {
					buffer[j] = 0.0;
					buffer[j + 1] = 0.0;
					buffer[j + 2] = 0.0;
				}
			}
		}
	}

	VelocitySensitivitySource::VelocitySensitivitySource(const detail::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	VelocitySensitivitySource::VelocitySensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void VelocitySensitivitySource::build()
	{
		m_schema.name = LADRUNO_MAKE_STRING("VELOCITY_SENSITIVITY_" << m_grad);
		m_schema.display_name = LADRUNO_MAKE_STRING("Velocity Sensitivity d_V/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Vx/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Vx/d_q" << m_grad << ",d_Vy/d_q" << m_grad);
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Vx/d_q" << m_grad << ",d_Vy/d_q" << m_grad << ",d_Vz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = LADRUNO_MAKE_STRING("Nodal velocity sensitivity field d_V/d_q" << m_grad);
		m_schema.result_type = detail::ResultType::Generic;
	}

	void VelocitySensitivitySource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "LadrunoRecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = nodeAt(info, i);
			if (!inode) continue;
			size_t j = i * m_ndim;
			buffer[j] = inode->getVelSensitivity(1, m_grad - 1);
			if (m_ndim > 1) {
				buffer[j + 1] = inode->getVelSensitivity(2, m_grad - 1);
				if (m_ndim > 2) {
					buffer[j + 2] = inode->getVelSensitivity(3, m_grad - 1);
				}
			}
		}
	}

	AngularVelocitySensitivitySource::AngularVelocitySensitivitySource(const detail::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	AngularVelocitySensitivitySource::AngularVelocitySensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void AngularVelocitySensitivitySource::build()
	{
		m_schema.name = LADRUNO_MAKE_STRING("ANGULAR_VELOCITY_SENSITIVITY_" << m_grad);
		m_schema.display_name = LADRUNO_MAKE_STRING("Angular Velocity Sensitivity d_RV/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_RVz/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_RVx/d_q" << m_grad << ",d_RVy/d_q" << m_grad << ",d_RVz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = LADRUNO_MAKE_STRING("Nodal angular velocity sensitivity field d_RV/d_q" << m_grad);
		m_schema.result_type = detail::ResultType::Generic;
	}

	void AngularVelocitySensitivitySource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "LadrunoRecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = nodeAt(info, i);
			if (!inode) continue;
			const Vector& iresponse = inode->getTrialVel(); // aux, just to check for angular velocityal components
			if (m_ndim == 2) {
				size_t j = i; // node_counter * 1 (1 angular velocity component in 2D)
				if (iresponse.Size() > 2)
					buffer[j] = inode->getVelSensitivity(3, m_grad - 1);
				else
					buffer[j] = 0.0;
			}
			else if (m_ndim == 3) {
				size_t j = i * 3; // node_counter * 3 (3 angular velocity components in 3D)
				if (iresponse.Size() > 5) {
					buffer[j] = inode->getVelSensitivity(4, m_grad - 1);
					buffer[j + 1] = inode->getVelSensitivity(5, m_grad - 1);
					buffer[j + 2] = inode->getVelSensitivity(6, m_grad - 1);
				}
				else {
					buffer[j] = 0.0;
					buffer[j + 1] = 0.0;
					buffer[j + 2] = 0.0;
				}
			}
		}
	}

	AccelerationSensitivitySource::AccelerationSensitivitySource(const detail::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	AccelerationSensitivitySource::AccelerationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void AccelerationSensitivitySource::build()
	{
		m_schema.name = LADRUNO_MAKE_STRING("ACCELERATION_SENSITIVITY_" << m_grad);
		m_schema.display_name = LADRUNO_MAKE_STRING("Acceleration Sensitivity d_A/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Ax/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Ax/d_q" << m_grad << ",d_Ay/d_q" << m_grad);
			m_schema.num_components = 2;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_Ax/d_q" << m_grad << ",d_Ay/d_q" << m_grad << ",d_Az/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = LADRUNO_MAKE_STRING("Nodal acceleration sensitivity field d_A/d_q" << m_grad);
		m_schema.result_type = detail::ResultType::Generic;
	}

	void AccelerationSensitivitySource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "LadrunoRecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		// Same as MPCORecorder: getVelSensitivity() is read here, not
		// getAccSensitivity(), so the two recorders return the same values.
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = nodeAt(info, i);
			if (!inode) continue;
			size_t j = i * m_ndim;
			buffer[j] = inode->getVelSensitivity(1, m_grad - 1);
			if (m_ndim > 1) {
				buffer[j + 1] = inode->getVelSensitivity(2, m_grad - 1);
				if (m_ndim > 2) {
					buffer[j + 2] = inode->getVelSensitivity(3, m_grad - 1);
				}
			}
		}
	}

	AngularAccelerationSensitivitySource::AngularAccelerationSensitivitySource(const detail::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	AngularAccelerationSensitivitySource::AngularAccelerationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void AngularAccelerationSensitivitySource::build()
	{
		m_schema.name = LADRUNO_MAKE_STRING("ANGULAR_ACCELERATION_SENSITIVITY_" << m_grad);
		m_schema.display_name = LADRUNO_MAKE_STRING("Angular Acceleration Sensitivity d_RA/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_RAz/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = detail::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = LADRUNO_MAKE_STRING("d_RAx/d_q" << m_grad << ",d_RAy/d_q" << m_grad << ",d_RAz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = detail::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = LADRUNO_MAKE_STRING("Nodal angular acceleration sensitivity field dRA/d_q" << m_grad);
		m_schema.result_type = detail::ResultType::Generic;
	}

	void AngularAccelerationSensitivitySource::evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "LadrunoRecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		// Same as MPCORecorder: getVelSensitivity() is read here, not the
		// acceleration sensitivity.
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = nodeAt(info, i);
			if (!inode) continue;
			const Vector& iresponse = inode->getTrialVel(); // aux, just to check for angular accelerational components
			if (m_ndim == 2) {
				size_t j = i; // node_counter * 1 (1 angular acceleration component in 2D)
				if (iresponse.Size() > 2)
					buffer[j] = inode->getVelSensitivity(3, m_grad - 1);
				else
					buffer[j] = 0.0;
			}
			else if (m_ndim == 3) {
				size_t j = i * 3; // node_counter * 3 (3 angular acceleration components in 3D)
				if (iresponse.Size() > 5) {
					buffer[j] = inode->getVelSensitivity(4, m_grad - 1);
					buffer[j + 1] = inode->getVelSensitivity(5, m_grad - 1);
					buffer[j + 2] = inode->getVelSensitivity(6, m_grad - 1);
				}
				else {
					buffer[j] = 0.0;
					buffer[j + 1] = 0.0;
					buffer[j + 2] = 0.0;
				}
			}
		}
	}

} // namespace ladruno
