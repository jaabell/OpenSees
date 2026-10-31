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
// Description: element result source of the ladruno recorder.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#include "Ladruno_ElementResults.h"

#include <limits>
#include <sstream>

namespace ladruno {

	ElementResultSource::ElementResultSource(const std::vector<std::string>& request,
		const detail::element::OutputDescriptorHeader& header,
		detail::element::OutputResponseCollection& bucket)
		: m_header(header)
		, m_bucket(bucket)
		, m_size_warned(false)   // one-shot latch for the size-mismatch diagnostic
	{
		// Row identity = the element tag of each response in the bucket.
		m_ids.resize(m_bucket.items.size());
		for (size_t i = 0; i < m_bucket.items.size(); ++i)
			m_ids[i] = m_bucket.items[i].element->getTag();
		buildSchema(request);
	}

	void ElementResultSource::buildSchema(const std::vector<std::string>& request)
	{
		// The on-disk group name is the per-bucket dir_name produced by the
		// orchestrator (e.g. "204-FourNodeQuad[201:0:0]"), nested under the
		// per-result display name (e.g. "force" / "stress"). The recorder pairs
		// this source with a StreamingSink(OnElements) whose family root is
		// ON_ELEMENTS; the source name therefore embeds "<display>/<bucket>".
		std::stringstream ss_disp;
		for (size_t i = 0; i < request.size(); ++i) {
			if (i > 0) ss_disp << ".";
			ss_disp << request[i];
		}
		m_schema.display_name = ss_disp.str();
		// name = "<result_display>/<bucket_dir>" so the sink lands it at
		// ON_ELEMENTS/<result_display>/<bucket_dir>.
		std::stringstream ss_name;
		ss_name << m_schema.display_name << "/" << m_bucket.dir_name;
		m_schema.name = ss_name.str();

		// Element results: components/dimension are not a single CSV; the column
		// structure lives in the COLUMN_MAP (written separately by the recorder).
		// num_components is the flat NUM_COLUMNS so the sink writes [nElem x NCOL].
		m_schema.components_csv = "";
		m_schema.dimension = "";
		m_schema.description = "";
		m_schema.num_components = m_header.num_columns;
		m_schema.result_type = detail::ResultType::Generic;
		m_schema.data_type = detail::ResultDataType::Scalar;
	}

	void ElementResultSource::evaluate(const detail::ProcessInfo& /*info*/, std::vector<double>& buffer)
	{
		// For each element response in the bucket: getResponse(); copy its
		//   getInformation().getData() row into [i * NUM_COLUMNS .. +NUM_COLUMNS).
		const int num_columns = m_header.num_columns;
		const size_t num_rows = m_bucket.items.size();
		buffer.assign(num_rows * (size_t)num_columns, 0.0);
		for (size_t i = 0; i < num_rows; ++i) {
			detail::element::OutputResponse& current_response = m_bucket.items[i];
			current_response.response->getResponse();
			const Vector& current_data = current_response.response->getInformation().getData();
			if (current_data.Size() != num_columns) {
				// Reported once per bucket: the mismatch is a property of the
				// (element, bucket) pairing, so it recurs on every element on
				// every step. The row stays the declared width so the dataset
				// never desyncs, and it is NaN-filled, since a zero row is
				// indistinguishable from a real zero.
				if (!m_size_warned) {
					m_size_warned = true;
					opserr << "LadrunoRecorder warning: element "
					       << current_response.element->getTag()
					       << " in bucket " << m_bucket.dir_name.c_str()
					       << " returned " << current_data.Size()
					       << " values but the bucket declares " << num_columns
					       << " columns; that row is NaN-filled"
					          " (further size mismatches in this bucket are not reported)\n";
				}
				const size_t row = i * (size_t)num_columns;
				for (int j = 0; j < num_columns; ++j)
					buffer[row + (size_t)j] = std::numeric_limits<double>::quiet_NaN();
				continue;
			}
			size_t offset = i * (size_t)num_columns;
			for (int j = 0; j < num_columns; ++j)
				buffer[offset + (size_t)j] = current_data[j];
		}
	}

} // namespace ladruno
