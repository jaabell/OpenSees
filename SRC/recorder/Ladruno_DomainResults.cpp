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
// recorder (energy balance through EnergyBalanceKernel.h).

#include "Ladruno_DomainResults.h"

#include "EnergyBalanceKernel.h"

#include "Domain.h"
#include "Node.h"
#include "NodeIter.h"
#include "Element.h"
#include "ElementIter.h"
#include "MeshRegion.h"
#include "ID.h"

#include <cstddef>

namespace ladruno {

	// Instantaneous energy terms of a set of elements and nodes, summed with
	// the per-entity functions of EnergyBalanceKernel.h (the same functions
	// EnergyBalanceRecorder sums), so the two recorders cannot disagree.
	//   ke              = 1/2 v^T M v over element and nodal mass
	//   internal_rate   = F^T v of the element resisting forces
	//   damping_rate    = element damping power + nodal v^T C v
	//   unbalanced_rate = v^T P_ext of the nodal loads
	static void addElementTerms(Element* ele, double& ke, double& ir, double& dr)
	{
		ebkernel::ElementEnergyTerms t;
		ebkernel::elementEnergyTerms(ele, t);
		ke += t.ke; ir += t.ie; dr += t.dw;
	}

	static void addNodeTerms(Node* node, double& ke, double& dr, double& ur)
	{
		double nke = 0., ndw = 0., nulw = 0.;
		ebkernel::addNodeEnergy(node, nke, ndw, nulw);
		ke += nke; dr += ndw; ur += nulw;
	}

	static void sweepDomain(Domain* dom, double& ke, double& ir, double& dr, double& ur)
	{
		ke = 0.0; ir = 0.0; dr = 0.0; ur = 0.0;
		if (dom == 0)
			return;
		Element* ele;
		ElementIter& elements = dom->getElements();
		while ((ele = elements()) != 0)
			addElementTerms(ele, ke, ir, dr);
		Node* node;
		NodeIter& nodes = dom->getNodes();
		while ((node = nodes()) != 0)
			addNodeTerms(node, ke, dr, ur);
	}

	static void sweepRegion(Domain* dom, MeshRegion* reg, double& ke, double& ir,
	                        double& dr, double& ur)
	{
		ke = 0.0; ir = 0.0; dr = 0.0; ur = 0.0;
		if (dom == 0 || reg == 0)
			return;
		const ID& elems = reg->getElements();
		for (int i = 0; i < elems.Size(); ++i) {
			Element* ele = dom->getElement(elems(i));
			if (ele != 0)
				addElementTerms(ele, ke, ir, dr);
		}
		const ID& rnodes = reg->getNodes();
		for (int i = 0; i < rnodes.Size(); ++i) {
			Node* node = dom->getNode(rnodes(i));
			if (node != 0)
				addNodeTerms(node, ke, dr, ur);
		}
	}

	/* ===================================================================== */
	/* EnergyBalanceSource                                                   */
	/* ===================================================================== */

	EnergyBalanceSource::EnergyBalanceSource(const detail::ProcessInfo& /*info*/,
	                                         EnergyState* state)
		: m_per_region(false)
		, m_ids()
		, m_schema()
		, m_own()
		, m_state(state ? state : &m_own)
	{
		m_ids.push_back(0); // synthetic domain id {0}
		buildSchema();
	}

	EnergyBalanceSource::EnergyBalanceSource(const detail::ProcessInfo& /*info*/,
	                                         const std::vector<int>& region_tags,
	                                         EnergyState* state)
		: m_per_region(true)
		, m_ids(region_tags)
		, m_schema()
		, m_own()
		, m_state(state ? state : &m_own)
	{
		if (m_state->regions.size() != region_tags.size())
			m_state->regions.resize(region_tags.size());
		buildSchema();
	}

	void EnergyBalanceSource::buildSchema()
	{
		m_schema.name = "energyBalance";
		m_schema.display_name = "Energy Balance";
		m_schema.components_csv = "KE,IE,DW,ULW,RES,ERR";
		m_schema.num_components = ebkernel::NUM_ENERGY_COMPONENTS; // 6
		m_schema.dimension = "F*L";
		m_schema.description = "Structural-dynamics energy balance "
			"(kinetic/internal/damping/load-work + closure residual)";
		m_schema.data_type = detail::ResultDataType::Vectorial;
		m_schema.result_type = detail::ResultType::Generic;
	}

	void EnergyBalanceSource::advance(const detail::ProcessInfo& info)
	{
		EnergyState& st = *m_state;
		// once per commit: the recorder calls advance() on every commit and
		// evaluate() calls it again on a recorded one.
		if (st.has_last && st.last_tag == info.current_time_step_id)
			return;

		const size_t ncomp = (size_t)m_schema.num_components;
		st.last_out.assign(m_ids.size() * ncomp, 0.0);

		Domain* dom = info.domain;
		if (dom == 0)
			return;

		// dT from the previously integrated commit. The state outlives a stage
		// rebuild, so a new MODEL_STAGE continues the integrals instead of
		// restarting them with a rate x (t - 0) jump.
		const double dT = info.current_time_step - st.prev_time;

		if (!m_per_region) {
			double ke = 0., ir = 0., dr = 0., ur = 0.;
			sweepDomain(dom, ke, ir, dr, ur);
			double out[ebkernel::NUM_ENERGY_COMPONENTS];
			st.model.step(dT, st.first, ke, ir, dr, ur, out);
			for (size_t c = 0; c < ncomp; ++c)
				st.last_out[c] = out[c];
		}
		else {
			for (size_t i = 0; i < m_ids.size(); ++i) {
				MeshRegion* reg = dom->getRegion(m_ids[i]);
				double ke = 0., ir = 0., dr = 0., ur = 0.;
				sweepRegion(dom, reg, ke, ir, dr, ur);
				double out[ebkernel::NUM_ENERGY_COMPONENTS];
				st.regions[i].step(dT, st.first, ke, ir, dr, ur, out);
				const size_t base = i * ncomp;
				for (size_t c = 0; c < ncomp; ++c)
					st.last_out[base + c] = out[c];
			}
		}

		st.first = false;
		st.prev_time = info.current_time_step;
		st.last_tag = info.current_time_step_id;
		st.has_last = true;
	}

	void EnergyBalanceSource::evaluate(const detail::ProcessInfo& info,
	                                   std::vector<double>& buffer)
	{
		advance(info);
		const size_t n = m_ids.size() * (size_t)m_schema.num_components;
		buffer.assign(n, 0.0);
		for (size_t k = 0; k < n && k < m_state->last_out.size(); ++k)
			buffer[k] = m_state->last_out[k];
	}

} // namespace ladruno
