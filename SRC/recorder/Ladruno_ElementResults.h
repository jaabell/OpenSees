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
// Description: element result source of the ladruno recorder. Walks the
// OutputDescriptor tree an element publishes through setResponse() and
// flattens one result (element, Gauss point, section, fiber or material
// level) into a row per element.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#ifndef Ladruno_ElementResults_h
#define Ladruno_ElementResults_h

#include "Ladruno_Types.h"
#include "Ladruno_ResultIO.h"

#include "Element.h"
#include "ElementIter.h"
#include "Domain.h"
#include "Response.h"
#include "Information.h"
#include "CompositeResponse.h"
#include "section/SectionForceDeformation.h"
#include "Vector.h"
#include "Matrix.h"
#include "ID.h"
#include "OPS_Stream.h"

// Distinct OPS_Stream tag so this descriptor stream cannot collide with
// MPCORecorder's (1001) when both link into the same binary.
#define OPS_STREAM_TAGS_Ladruno_ElementOutputDescriptorStream 1002

// Class tags of elements compiled only in some builds, as in MPCORecorder.cpp.
// All #ifndef-guarded so they never clash with identical definitions.
#ifndef LADRUNO_MAX_TRIAL_NSEC
#define LADRUNO_MAX_TRIAL_NSEC 100
#endif
#ifndef ELE_TAG_FourNodeQuadWithSensitivity
#define ELE_TAG_FourNodeQuadWithSensitivity 100000011
#endif
#ifndef ELE_TAG_DispBeamColumn2dWithSensitivity
#define ELE_TAG_DispBeamColumn2dWithSensitivity 102030
#endif
#ifndef ELE_TAG_DispBeamColumn3dWithSensitivity
#define ELE_TAG_DispBeamColumn3dWithSensitivity 1110000
#endif
#ifndef ELE_TAG_BbarBrickWithSensitivity
#define ELE_TAG_BbarBrickWithSensitivity 1984587234
#endif
#ifndef ELE_TAG_AC3D8HexWithSensitivity
#define ELE_TAG_AC3D8HexWithSensitivity 100001
#endif
#ifndef Ele_TAG_Elastic2dGNL
#define Ele_TAG_Elastic2dGNL -1
#endif
#ifndef TAG_InelasticYS2DGNL
#define TAG_InelasticYS2DGNL -1
#endif

namespace ladruno {

namespace detail {

	namespace element {

		struct OutputDescriptorHeader
		{
			OutputDescriptorHeader()
				: num_columns(0)
				, num_components()
				, gauss_id()
				, components_path()
				, components()
				, multiplicity()
			{}

			inline bool operator == (const OutputDescriptorHeader &other) const {
				if (this != &other) {
					if (*this < other) return false;
					if (*this > other) return false;
					return true;
				}
				return true;
			}

			inline bool operator < (const OutputDescriptorHeader &other) const {
				if (num_columns < other.num_columns) return true;
				if (num_columns > other.num_columns) return false;
				if (multiplicity.size() < other.multiplicity.size()) return true;
				if (multiplicity.size() > other.multiplicity.size()) return false;
				for (size_t i = 0; i < multiplicity.size(); i++) {
					if (multiplicity[i] < other.multiplicity[i]) return true;
					if (multiplicity[i] > other.multiplicity[i]) return false;
				}
				if (num_components.size() < other.num_components.size()) return true;
				if (num_components.size() > other.num_components.size()) return false;
				for (size_t i = 0; i < num_components.size(); i++) {
					if (num_components[i] < other.num_components[i]) return true;
					if (num_components[i] > other.num_components[i]) return false;
				}
				if (gauss_id.size() < other.gauss_id.size()) return true;
				if (gauss_id.size() > other.gauss_id.size()) return false;
				for (size_t i = 0; i < gauss_id.size(); i++) {
					if (gauss_id[i] < other.gauss_id[i]) return true;
					if (gauss_id[i] > other.gauss_id[i]) return false;
				}
				if (components_path.size() < other.components_path.size()) return true;
				if (components_path.size() > other.components_path.size()) return false;
				for (size_t i = 0; i < components_path.size(); i++) {
					if (components_path[i] < other.components_path[i]) return true;
					if (components_path[i] > other.components_path[i]) return false;
				}
				if (components.size() < other.components.size()) return true;
				if (components.size() > other.components.size()) return false;
				for (size_t i = 0; i < components.size(); i++) {
					if (components[i] < other.components[i]) return true;
					if (components[i] > other.components[i]) return false;
				}
				return false;
			}

			inline bool operator > (const OutputDescriptorHeader &other) const {
				if (num_columns > other.num_columns) return true;
				if (num_columns < other.num_columns) return false;
				if (multiplicity.size() > other.multiplicity.size()) return true;
				if (multiplicity.size() < other.multiplicity.size()) return false;
				for (size_t i = 0; i < multiplicity.size(); i++) {
					if (multiplicity[i] > other.multiplicity[i]) return true;
					if (multiplicity[i] < other.multiplicity[i]) return false;
				}
				if (num_components.size() > other.num_components.size()) return true;
				if (num_components.size() < other.num_components.size()) return false;
				for (size_t i = 0; i < num_components.size(); i++) {
					if (num_components[i] > other.num_components[i]) return true;
					if (num_components[i] < other.num_components[i]) return false;
				}
				if (gauss_id.size() > other.gauss_id.size()) return true;
				if (gauss_id.size() < other.gauss_id.size()) return false;
				for (size_t i = 0; i < gauss_id.size(); i++) {
					if (gauss_id[i] > other.gauss_id[i]) return true;
					if (gauss_id[i] < other.gauss_id[i]) return false;
				}
				if (components_path.size() > other.components_path.size()) return true;
				if (components_path.size() < other.components_path.size()) return false;
				for (size_t i = 0; i < components_path.size(); i++) {
					if (components_path[i] > other.components_path[i]) return true;
					if (components_path[i] < other.components_path[i]) return false;
				}
				if (components.size() > other.components.size()) return true;
				if (components.size() < other.components.size()) return false;
				for (size_t i = 0; i < components.size(); i++) {
					if (components[i] > other.components[i]) return true;
					if (components[i] < other.components[i]) return false;
				}
				return false;
			}

			inline std::string toString()const {
				std::stringstream ss;
				ss << "LadrunoRecorder Element Output Descriptor Header:\n";
				ss << "Columns: " << num_columns << "\n";
				ss << "N. Components:\n";
				for (size_t i = 0; i < num_components.size(); i++) {
					ss << num_components[i] << "   ";
				}
				ss << "\n\nPaths:\n";
				for (std::vector<std::vector<int> >::const_iterator
					it1 = components_path.begin(); it1 != components_path.end(); ++it1) {
					for (std::vector<int>::const_iterator
						it2 = it1->begin(); it2 != it1->end(); ++it2) {
						ss << detail::ElementOutputDescriptorType::toString((detail::ElementOutputDescriptorType::Enum)*it2) << ".";
					}
					ss << "\n";
				}
				ss << "\nComponents:\n";
				for (std::vector<std::vector<std::string> >::const_iterator
					it1 = components.begin(); it1 != components.end(); ++it1) {
					for (std::vector<std::string>::const_iterator
						it2 = it1->begin(); it2 != it1->end(); ++it2) {
						ss << *it2 << "   ";
					}
					ss << "\n";
				}
				return ss.str();
			}

			void workaroundForSizeInconsistency(int data_size) {
				/** WORKAROUND
				ouch! in some cases (see some material models), the material in the setResponse method, gives a vector of size N
				but without specifying the response type, which is what gives the components. In those cases there is
				a mismatch between response->getInformation().getData().size() and header.num_components.
				note: this workaround applies if and only if all items have zero components
				*/
				bool all_zero = true;
				for (size_t i = 0; i < num_components.size(); i++) {
					if (num_components[i] != 0) {
						all_zero = false;
						break;
					}
				}
				if (all_zero && data_size > 0) {
					if (num_components.size() == 0) { // empty
						num_components.resize(1);
						num_components[0] = data_size;
						gauss_id.resize(1);
						gauss_id[0] = -1;
						components_path.resize(1);
						components_path[0].push_back(detail::ElementOutputDescriptorType::Element);
						components.resize(1);
						components[0].resize((size_t)data_size);
						for (size_t i = 0; i < (size_t)data_size; i++) {
							std::stringstream aux; aux << "C" << i + 1; components[0][i] = aux.str();
						}
						multiplicity.resize(1);
						multiplicity[0] = 1;
						num_columns = data_size;
					}
					else {
						int num_items = (int)num_components.size();
						/* only if the data size divides equally by the number of items,
						i.e. each gauss point is assumed to have the same components
						*/
						if (data_size % num_items == 0) {
							int n_comp = data_size / num_items;
							for (size_t i = 0; i < num_components.size(); i++) {
								num_components[i] = n_comp;
								multiplicity[i] = 1;
								components[i].resize((size_t)n_comp);
								for (int j = 0; j < n_comp; j++) {
									std::stringstream aux; aux << "C" << j + 1; components[i][j] = aux.str();
								}
							}
							// gauss and comp paths are left as they are
							num_columns = data_size;
						}
					}
				}
			}

			void workaroundForDuplicatedComponents()
			{
				for (size_t i = 0; i < components.size(); i++) {
					std::vector<std::string> &i_components = components[i];
					if (i_components.size() > 0) {
						std::map<std::string, int> aux;
						for (size_t j = 0; j < i_components.size(); j++) {
							int &dupl_counter = aux[i_components[j]];
							if (dupl_counter > 0) {
								std::stringstream ss;
								ss << i_components[j] << "(" << dupl_counter << ")";
								i_components[j] = ss.str();
							}
							dupl_counter++;
						}
					}
				}
			}

			int num_columns;
			std::vector<int> num_components;
			std::vector<int> gauss_id; // note: 0-based
			std::vector<std::vector<int> > components_path;
			std::vector<std::vector<std::string> > components;
			std::vector<int> multiplicity;
		};

		// ─── Self-described element basis (Element-contract Part A) ──────
		// Filled from the OPTIONAL "basisInfo" probe. Lets a non-Lagrange
		// element (e.g. a 6-node Bezier triangle: family=bernstein, degree 2)
		// override the family/order/topology/... that would otherwise be
		// *guessed* from the legacy (geometry, integration-rule) pair.
		// Lagrange elements simply don't answer the probe ⇒ valid stays
		// false ⇒ the geometry-derived defaults are kept.
		struct BasisInfo
		{
			bool valid = false;
			std::string family;        // "bernstein", "lagrange", ...
			std::string topology;      // "tri", "quad", "hex", ...
			std::string param_domain;  // "bary", "[-1,1]", ...
			int order = -1;            // total polynomial degree (-1 = unset)
			int num_ctrl = -1;         // number of control points (-1 = unset)
			int num_gp = -1;           // number of result Gauss points (-1 = unset)
			int rational = -1;         // 0/1 (-1 = unset)
		};

		class OutputDescriptor
		{
		public:
			OutputDescriptor()
				: type(detail::ElementOutputDescriptorType::Element)
				, tag(0)
				, dummy_section_flag(false)
				, gp_number(0)
				, gp_eta(0.0)
				, gp_weight(0.0)
				, fib_y(0.0)
				, fib_z(0.0)
				, fib_a(0.0)
				, components()
				, items()
			{}

			OutputDescriptor(const OutputDescriptor &other)
				: type(other.type)
				, tag(other.tag)
				, dummy_section_flag(other.dummy_section_flag)
				, gp_number(other.gp_number)
				, gp_eta(other.gp_eta)
				, gp_weight(other.gp_weight)
				, fib_y(other.fib_y)
				, fib_z(other.fib_z)
				, fib_a(other.fib_a)
				, components(other.components)
				, items()
			{
				// make a deep copy of subitems
				items.resize(other.items.size());
				for (size_t i = 0; i < other.items.size(); i++)
					items[i] = new OutputDescriptor(*other.items[i]);
			}

			~OutputDescriptor() {
				for (size_t i = 0; i < items.size(); i++)
					if (items[i])
						delete items[i];
			}

			OutputDescriptor &operator = (const OutputDescriptor &other) {
				if (this != &other) {
					type = other.type;
					tag = other.tag;
					dummy_section_flag = other.dummy_section_flag;
					gp_number = other.gp_number;
					gp_eta = other.gp_eta;
					gp_weight = other.gp_weight;
					fib_y = other.fib_y;
					fib_z = other.fib_z;
					fib_a = other.fib_a;
					components = other.components;
					// make a deep copy of subitems
					items.resize(other.items.size());
					for (size_t i = 0; i < other.items.size(); i++)
						items[i] = new OutputDescriptor(*other.items[i]);
				}
				return *this;
			}

		public:
			// generic
			detail::ElementOutputDescriptorType::Enum type;
			// for material or section
			int tag;
			bool dummy_section_flag;
			// for gauss point
			int gp_number;
			double gp_eta; // use only eta (for 1D custom integration)
			double gp_weight;
			// for fibers
			double fib_y;
			double fib_z;
			double fib_a;
			// components (if items.size() == 0)
			std::vector<std::string> components;
			// subitems
			std::vector<OutputDescriptor*> items;

		public:
			std::string toString()const {
				std::stringstream ss;
				/*
				this prints a xml-like structure, just for debugging purposes
				*/
				ss << "----------------------------------------------------\n";
				ss << "OutputDescriptor info\n";
				ss << "----------------------------------------------------\n";
				ss << "XML-like structure\n";
				printInfo(0, ss);
				ss << "----------------------------------------------------\n";
				return ss.str();
			}

			void getGaussLocations(std::vector<double> &x) const {
				x.clear();
				appendGaussLocation(x);
				if (x.size() == 1) {
					x[0] = 0.0;
				}
				else {
					double xmin = std::numeric_limits<double>::max();
					double xmax = -xmin;
					for (size_t i = 0; i < x.size(); i++) {
						double ieta = x[i];
						if (ieta < xmin)
							xmin = ieta;
						else if (ieta > xmax)
							xmax = ieta;
					}
					double span = xmax - xmin;
					if (span == 0.0) {
						for (size_t i = 0; i < x.size(); i++)
							x[i] = 0.0;
					}
					else {
						for (size_t i = 0; i < x.size(); i++)
							x[i] = 2.0*(x[i] - xmin) / span - 1.0;
					}
				}
			}

			void appendGaussLocation(std::vector<double>& x) const {
				if (type == detail::ElementOutputDescriptorType::Gauss)
					x.push_back(gp_eta);
				for (size_t i = 0; i < items.size(); i++)
					items[i]->appendGaussLocation(x);
			}

			void appendGaussWeight(std::vector<double>& x) const {
				if (type == detail::ElementOutputDescriptorType::Gauss)
					x.push_back(gp_weight);
				for (size_t i = 0; i < items.size(); i++)
					items[i]->appendGaussWeight(x);
			}

			void getFiberData(std::vector<detail::element::FiberData> &data,
				std::vector<int> &data_mat_id,
				std::vector<int> &sec_id,
				std::vector<int> &gp_id,
				std::vector<bool> &dummy_sec_flags)const {
				data.clear();
				data_mat_id.clear();
				sec_id.clear();
				gp_id.clear();
				int *temp_gp = 0;
				int *temp_sec = 0;
				bool *temp_dummy = 0;
				appendFiberData(data, data_mat_id,
					sec_id, gp_id, dummy_sec_flags, temp_sec, temp_gp, temp_dummy);
				if (temp_gp)
					delete temp_gp;
				if (temp_sec)
					delete temp_sec;
				if (temp_dummy)
					delete temp_dummy;
			}

			detail::element::OutputDescriptorHeader makeHeader()const {
				detail::element::OutputDescriptorHeader header;
				std::list<int> temp_path;
				int temp_gp_id(-1);
				makeHeaderInternal(header, temp_path, temp_gp_id);
				return header;
			}

			void fixFloatingFiberOutput() {
				/*
				for some reason, some fiber-based cross section do not write the SectionOutput-tag before the FiberOutput-tag.
				this make things complicated and not robust. this is a workaround to fix this problem.
				*/
				fixFloatingFiberOutputInternal();
			}

			void fixSectionAfterFiberDueToFiberOutputFail() {
				/*
				due to a recent commit, the SectionForceDeformation first checks for fibers.
				If it fails when the fiber does not have the requested output, the SectionForceDeformation
				falls back to its setResponse, thus opening a section tag after a fiber tag without giving any result.
				This messes up everything!
				*/
				fixSectionAfterFiberDueToFiberOutputFailInternal();
			}

			int getNextGpTag() {
				int next_gp_tag = -1;
				getNextGpTagInternal(next_gp_tag);
				return next_gp_tag + 1;
			}

			void purge() {
				mergeGaussInternal();
				mergeSecInternal();
			}

		private:
			void printInfo(int level, std::stringstream &ss) const {
				std::stringstream ss_indent;
				for (int i = 0; i < level; i++) ss_indent << "\t";
				std::string indent = ss_indent.str();
				ss << indent << "<" << detail::ElementOutputDescriptorType::toString(this->type);
				if (this->type == detail::ElementOutputDescriptorType::Gauss) {
					ss << " number=\"" << this->gp_number << "\" eta=\"" << this->gp_eta << "\" weight=\"" << this->gp_weight << "\"";
				}
				else if (this->type == detail::ElementOutputDescriptorType::Section) {
					ss << " tag=\"" << this->tag << "\"";
				}
				else if (this->type == detail::ElementOutputDescriptorType::Material) {
					ss << " tag=\"" << this->tag << "\"";
				}
				ss << ">\n";
				for (int i = 0; i < this->components.size(); i++) {
					ss << indent << "\t" << components[i] << "\n";
				}
				for (int i = 0; i < items.size(); i++) {
					items[i]->printInfo(level + 1, ss);
				}
				ss << indent << "</" << detail::ElementOutputDescriptorType::toString(this->type) << ">\n";
			}

			void makeHeaderInternal(detail::element::OutputDescriptorHeader &header, std::list<int> &temp_path, int &temp_gp_id) const {
				if (type == detail::ElementOutputDescriptorType::Gauss)
					temp_gp_id = gp_number;
				temp_path.push_back((int)type);
				if (components.size() > 0 || items.size() == 0) {
					header.num_columns += (int)components.size();
					int next_num_components = (int)components.size();
					int next_gauss_id = temp_gp_id;
					std::vector<int> next_components_path(temp_path.begin(), temp_path.end());
					if (header.multiplicity.size() == 0) {
						header.num_components.push_back(next_num_components);
						header.gauss_id.push_back(next_gauss_id);
						header.components_path.push_back(next_components_path);
						header.components.push_back(components);
						header.multiplicity.push_back(1);
					}
					else {
						bool equal_to_previous = false;
						size_t last_index = header.multiplicity.size() - 1;
						if (header.num_components[last_index] == next_num_components) {
							if (header.gauss_id[last_index] == next_gauss_id) {
								if (utils::misc::areVectorsEqual(header.components_path[last_index], next_components_path)) {
									if (utils::misc::areVectorsEqual(header.components[last_index], components)) {
										equal_to_previous = true;
									}
								}
							}
						}
						if (equal_to_previous) {
							header.multiplicity[last_index]++;
						}
						else {
							header.num_components.push_back(next_num_components);
							header.gauss_id.push_back(next_gauss_id);
							header.components_path.push_back(next_components_path);
							header.components.push_back(components);
							header.multiplicity.push_back(1);
						}
					}
				}
				for (size_t i = 0; i < items.size(); i++)
					items[i]->makeHeaderInternal(header, temp_path, temp_gp_id);
				temp_path.pop_back();
			}

			void appendFiberData(std::vector<detail::element::FiberData> &data, std::vector<int> &data_mat_id,
				std::vector<int> &sec_id, std::vector<int> &gp_id, std::vector<bool> &dummy_sec_flags,
				int* &temp_sec, int* &temp_gp, bool* &temp_dummy) const {
				if (type == detail::ElementOutputDescriptorType::Gauss) {
					if (temp_gp == 0)
						temp_gp = new int();
					*temp_gp = gp_number;
				}
				else if (type == detail::ElementOutputDescriptorType::Section) {
					if (temp_sec == 0)
						temp_sec = new int();
					if (temp_dummy == 0)
						temp_dummy = new bool();
					*temp_sec = tag;
					*temp_dummy = dummy_section_flag;
				}
				else if (type == detail::ElementOutputDescriptorType::Fiber) {
					data.push_back(detail::element::FiberData(fib_y, fib_z, fib_a));
					int fiber_mat_tag = -1;
					if (items.size() == 1) {
						if (items[0]->type == detail::ElementOutputDescriptorType::Material) {
							fiber_mat_tag = items[0]->tag;
						}
					}
					data_mat_id.push_back(fiber_mat_tag);
				}
				if (type == detail::ElementOutputDescriptorType::Fiber || items.size() == 0) {
					if (temp_sec)
						sec_id.push_back(*temp_sec);
					if (temp_gp)
						gp_id.push_back(*temp_gp);
					if (temp_dummy)
						dummy_sec_flags.push_back(*temp_dummy);
					/** note: return here, without asking the sub-items.
					for example if this item is a fiber, it may have sub-times like material!
					*/
					return;
				}
				// ask subitems
				for (size_t i = 0; i < items.size(); i++) {
					items[i]->appendFiberData(data, data_mat_id,
						sec_id, gp_id, dummy_sec_flags, temp_sec, temp_gp, temp_dummy);
				}
			}

			void fixFloatingFiberOutputInternal() {
				if (items.size() > 0) {
					if (type != detail::ElementOutputDescriptorType::Section) {
						if (items[0]->type == detail::ElementOutputDescriptorType::Fiber) {
							/* check only the first one. items are of the same type... */
							OutputDescriptor *dummy_section_level = new OutputDescriptor();
							dummy_section_level->type = detail::ElementOutputDescriptorType::Section;
							dummy_section_level->tag = -123456;
							dummy_section_level->dummy_section_flag = true;
							dummy_section_level->items = items;
							items.clear();
							items.push_back(dummy_section_level);
						}
					}
					for (size_t i = 0; i < items.size(); i++)
						items[i]->fixFloatingFiberOutputInternal();
				}
			}

			void fixSectionAfterFiberDueToFiberOutputFailInternal() {
				if (items.size() > 0) {
					if (items[0]->type == detail::ElementOutputDescriptorType::Fiber) {
						/* check only the first one. items are of the same type...*/
						if (items.size() > 1) {
							if (items.back()->type != detail::ElementOutputDescriptorType::Fiber) {
								items.pop_back();
							}
						}
					}
					for (size_t i = 0; i < items.size(); i++)
						items[i]->fixSectionAfterFiberDueToFiberOutputFailInternal();
				}
			}

			void getNextGpTagInternal(int &next_gp_tag) {
				if (type == detail::ElementOutputDescriptorType::Gauss) {
					if (next_gp_tag < gp_number)
						next_gp_tag = gp_number;
				}
				else {
					for (size_t i = 0; i < items.size(); i++)
						items[i]->getNextGpTagInternal(next_gp_tag);
				}
			}

			void mergeGaussInternal() {
				/* if multiple gauss items with same id exist, merge their contents
				\todo: check if the order of fibers is preserved..
				note mandatory now, because when asking for multiple sections, duplicate gauss points are produced
				but only the first one is filled
				*/
				if (items.size() > 0) {
					if (items[0]->type == detail::ElementOutputDescriptorType::Gauss) {
						// all sub items are gauss descriptors, let's merge them
						std::map<int, OutputDescriptor*> aux;
						for (size_t i = 0; i < items.size(); i++) {
							OutputDescriptor* curr_item = items[i];
							std::map<int, OutputDescriptor*>::iterator it = aux.find(curr_item->gp_number);
							if (it == aux.end()) { // non existing
								aux[curr_item->gp_number] = curr_item;
							}
							else { // existing, move curr_item->items into existing_item->items
								OutputDescriptor* existing_item = it->second;
								for (size_t j = 0; j < curr_item->items.size(); j++)
									existing_item->items.push_back(curr_item->items[j]);
								curr_item->items.clear();
							}
						}
						items.clear();
						for (std::map<int, OutputDescriptor*>::iterator it = aux.begin(); it != aux.end(); ++it) {
							items.push_back(it->second);
						}
					}
					else {
						for (size_t i = 0; i < items.size(); i++)
							items[i]->mergeGaussInternal();
					}
				}
			}

			void mergeSecInternal() {
				if (items.size() > 0) {
					if (items[0]->type == detail::ElementOutputDescriptorType::Section) {
						// all sub items are section descriptors, let's merge them
						std::map<int, OutputDescriptor*> aux;
						for (size_t i = 0; i < items.size(); i++) {
							OutputDescriptor* curr_item = items[i];
							std::map<int, OutputDescriptor*>::iterator it = aux.find(curr_item->tag);
							if (it == aux.end()) { // non existing
								aux[curr_item->tag] = curr_item;
							}
							else { // existing, move curr_item->items into existing_item->items
								OutputDescriptor* existing_item = it->second;
								for (size_t j = 0; j < curr_item->items.size(); j++)
									existing_item->items.push_back(curr_item->items[j]);
								curr_item->items.clear();
							}
						}
						items.clear();
						for (std::map<int, OutputDescriptor*>::iterator it = aux.begin(); it != aux.end(); ++it) {
							items.push_back(it->second);
						}
					}
					else {
						for (size_t i = 0; i < items.size(); i++)
							items[i]->mergeSecInternal();
					}
				}
			}
		};

		class OutputDescriptorStream : public OPS_Stream
		{
		public:
			enum StreamErrorCode {
				ERROR_CODE_OK = 0,
				ERROR_CODE_SECTION_AFTER_FIBER,
				ERROR_CODE_GENERIC
			};
		public:
			OutputDescriptorStream(detail::element::OutputDescriptor * _d)
				: OPS_Stream(OPS_STREAM_TAGS_Ladruno_ElementOutputDescriptorStream)
				, descr(_d)
				, current_level(0)
				, pending_close_tag(false)
				, error_code(ERROR_CODE_OK)
				, in_basis(false)
			{}
			~OutputDescriptorStream() {}

			// Self-description captured from an optional "basisInfo" probe.
			// Stays default/invalid for the normal result-discovery path.
			detail::element::BasisInfo basis_info;
			bool in_basis;

			int tag(const char *name) {
				// get the element output descriptor at current level and id
				detail::element::OutputDescriptor *eo_curr_lev = descr;
				for (int i = 1; i <= current_level; i++) {
					if (eo_curr_lev->items.size() == 0) {
						opserr << "LadrunoRecorder Error: cannot set attribute(name, int), empty item list.\n";
						error_code = ERROR_CODE_GENERIC; return -1;  // drop this bucket; the run is not aborted
					}
					eo_curr_lev = eo_curr_lev->items[eo_curr_lev->items.size() - 1];
				}
				if (current_level == 0) {
					if (strcmp(name, "ElementOutput") == 0) {
						/** nothing to do. this is the root of the result tree*/
					}
					/* gauss output is the first entry */
					else if (strcmp(name, "GaussPoint") == 0 || strcmp(name, "GaussPointOutput") == 0) {
						detail::element::OutputDescriptor *eo_new_curr_lev = new detail::element::OutputDescriptor();
						eo_new_curr_lev->type = detail::ElementOutputDescriptorType::Gauss;
						ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
						eo_curr_lev->items.push_back(eo_new_curr_lev);
						current_level++;
					}
					/* self-describing basis probe (Element-contract Part A):
					   capture the <ElementBasis> attributes WITHOUT opening a
					   gauss level. Only emitted in response to a "basisInfo"
					   request, never during normal result discovery. */
					else if (strcmp(name, "ElementBasis") == 0) {
						in_basis = true;
						basis_info.valid = true;
					}
					else {
						/*opserr <<
						"LadrunoRecorder Error: invalid tag at level 0:\n"
						"expected \"GaussPoint\" or \"GaussPointOutput\", given \"" << name << "\"\n";
						exit(-1);*/
						/*
						let's try a last workaround for an inconsistency problem found in SSPbrick:
						in that element, in setResponse, the material's setResponse method is called without opening a GaussPoint tag...
						*/
						//simulate: tag("GaussPoint");
						{
							detail::element::OutputDescriptor *eo_new_curr_lev = new detail::element::OutputDescriptor();
							eo_new_curr_lev->type = detail::ElementOutputDescriptorType::Gauss;
							ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
							eo_curr_lev->items.push_back(eo_new_curr_lev);
							current_level++;
						}
						// create the attribute 
						attr("number", descr->getNextGpTag());
						/*
						close the gauss tag. It cannot be done here: pending_close_tag is set to true, so that
						when the 'real' tag gets closed, the 'manual' gauss tag is closed automatically
						*/
						pending_close_tag = true;
						// recursion: recall this request now inside a gauss point tag
						tag(name);
					}
				}
				else if (current_level > 0) {
					if (strcmp(name, "NdMaterialOutput") == 0 || strcmp(name, "UniaxialMaterialOutput") == 0) {
						// its parent can be anything but ElementOutput, ok if current_level > 1
						detail::element::OutputDescriptor *eo_new_curr_lev = new detail::element::OutputDescriptor();
						eo_new_curr_lev->type = detail::ElementOutputDescriptorType::Material;
						ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
						if (eo_curr_lev->items.size() > 0) {
							// multiple materials cannot be children of same gauss/fiber point. this happens when
							// an objects opens the tag, fails in getting response, and falls back to base class implementation,
							// which opens again the same tag
							for (detail::element::OutputDescriptor* sub_item : eo_curr_lev->items)
								delete sub_item;
							eo_curr_lev->items.clear();
						}
						eo_curr_lev->items.push_back(eo_new_curr_lev);
						current_level++;
					}
					else if (strcmp(name, "SectionOutput") == 0 || strcmp(name, "SectionForceDeformation") == 0) {
						// its parent can be GaussOutput or another SectionOutput
						if (!(eo_curr_lev->type == detail::ElementOutputDescriptorType::Gauss || eo_curr_lev->type == detail::ElementOutputDescriptorType::Section)) {
							opserr <<
								"LadrunoRecorder Error: invalid parent for \"" << name << "\" tag:\n"
								"expected \"GaussOutput\" or \"GaussPointOutput\""
								" or \"SectionOutput\" or \"SectionForceDeformation\", parent tag = \""
								<< detail::ElementOutputDescriptorType::toString(eo_curr_lev->type) << "\"\n";
							error_code = ERROR_CODE_GENERIC; return -1;  // drop this bucket; the run is not aborted
						}
						detail::element::OutputDescriptor *eo_new_curr_lev = new detail::element::OutputDescriptor();
						eo_new_curr_lev->type = detail::ElementOutputDescriptorType::Section;
						ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
						if (error_code == ERROR_CODE_OK) {
							if (eo_curr_lev->items.size() > 0) {
								// multiple sections cannot be children of same gauss point. this happens when
								// an objects opens the tag, fails in getting response, and falls back to base class implementation,
								// which opens again the same tag
								for (detail::element::OutputDescriptor* sub_item : eo_curr_lev->items)
									delete sub_item;
								eo_curr_lev->items.clear();
							}
							// do the above check only if there is no inconsistency with previous items!
						}
						eo_curr_lev->items.push_back(eo_new_curr_lev);
						current_level++;
					}
					else if (strcmp(name, "FiberOutput") == 0) {
						// its parent can be only a SectionOutput
						if (!(eo_curr_lev->type == detail::ElementOutputDescriptorType::Section || eo_curr_lev->type == detail::ElementOutputDescriptorType::Gauss)) {
							opserr <<
								"LadrunoRecorder Error: invalid parent for \"" << name << "\" tag:\n"
								"expected \"GaussOutput\" or \"GaussPointOutput\""
								" or \"SectionOutput\" or \"SectionForceDeformation\", parent tag = \""
								<< detail::ElementOutputDescriptorType::toString(eo_curr_lev->type) << "\"\n";
							error_code = ERROR_CODE_GENERIC; return -1;  // drop this bucket; the run is not aborted
						}
						detail::element::OutputDescriptor *eo_new_curr_lev = new detail::element::OutputDescriptor();
						eo_new_curr_lev->type = detail::ElementOutputDescriptorType::Fiber;
						ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
						eo_curr_lev->items.push_back(eo_new_curr_lev);
						current_level++;
					}
					else {
						/*
						let's try a last workaround for an inconsistency problem found in ForceBeamColumn3d:
						in that element, in setResponse, if multiple sections are requested the GaussOutput tag is not
						closed...
						*/
						if ((strcmp(name, "GaussPoint") == 0 || strcmp(name, "GaussPointOutput") == 0) && current_level == 1) {
							endTag();
							tag(name);
						}
						else {
							opserr <<
								"LadrunoRecorder Error: invalid tag at level " << current_level << ":\n"
								"expected \"NdMaterialOutput\" or \"UniaxialMaterialOutput\""
								" or \"SectionOutput\" or \"SectionForceDeformation\""
								" or \"FiberOutput\""
								", given \"" << name << "\"\n";
							error_code = ERROR_CODE_GENERIC; return -1;  // drop this bucket; the run is not aborted
						}
					}
				}
				return 0;
			};

			int tag(const char *name, const char *value) {
				detail::element::OutputDescriptor *eo_curr_lev = descr;
				for (int i = 1; i <= current_level; i++) {
					if (eo_curr_lev->items.size() == 0) {
						opserr << "LadrunoRecorder Error: cannot set attribute(name, int), empty item list.\n";
						error_code = ERROR_CODE_GENERIC; return -1;  // drop this bucket; the run is not aborted
					}
					eo_curr_lev = eo_curr_lev->items[eo_curr_lev->items.size() - 1];
				}
				if (strcmp(name, "ResponseType") == 0)
					eo_curr_lev->components.push_back(value);
				return 0;
			};

			int endTag() {
				if (in_basis) {
					// closing <ElementBasis>; it never opened a level
					in_basis = false;
					return 0;
				}
				if (current_level > 0) {
					// decrement the current level
					current_level--;
				}
				if (pending_close_tag) {
					// once more...
					if (current_level > 0) {
						// decrement the current level
						current_level--;
					}
					pending_close_tag = false;
				}
				return 0;
			};

			int attr(const char *name, int value) {
				if (in_basis) {
					if (strcmp(name, "rational") == 0) basis_info.rational = value;
					else if (strcmp(name, "numCtrl") == 0) basis_info.num_ctrl = value;
					else if (strcmp(name, "numGP") == 0) basis_info.num_gp = value;
					else if (strcmp(name, "orderU") == 0) basis_info.order = value;
					return 0;
				}
				if (current_level > 0) {
					detail::element::OutputDescriptor *eo_curr_lev = descr;
					for (int i = 1; i <= current_level; i++) {
						if (eo_curr_lev->items.size() == 0) {
							opserr << "LadrunoRecorder Error: cannot set attribute(name, int), empty item list.\n";
							error_code = ERROR_CODE_GENERIC; return -1;  // drop this bucket; the run is not aborted
						}
						eo_curr_lev = eo_curr_lev->items[eo_curr_lev->items.size() - 1];
					}
					if (eo_curr_lev->type == detail::ElementOutputDescriptorType::Gauss) {
						if (strcmp(name, "number") == 0)
							eo_curr_lev->gp_number = value - 1; // note: make it 0-based
					}
					else if (eo_curr_lev->type == detail::ElementOutputDescriptorType::Material) {
						if (strcmp(name, "tag") == 0 || strcmp(name, "matTag") == 0)
							eo_curr_lev->tag = value;
					}
					else if (eo_curr_lev->type == detail::ElementOutputDescriptorType::Section) {
						if (strcmp(name, "tag") == 0 || strcmp(name, "secTag") == 0)
							eo_curr_lev->tag = value;
					}
				}
				return 0;
			};

			int attr(const char *name, double value) {
				if (current_level > 0) {
					detail::element::OutputDescriptor *eo_curr_lev = descr;
					for (int i = 1; i <= current_level; i++) {
						if (eo_curr_lev->items.size() == 0) {
							opserr << "LadrunoRecorder Error: cannot set attribute(name, int), empty item list.\n";
							error_code = ERROR_CODE_GENERIC; return -1;  // drop this bucket; the run is not aborted
						}
						eo_curr_lev = eo_curr_lev->items[eo_curr_lev->items.size() - 1];
					}
					if (eo_curr_lev->type == detail::ElementOutputDescriptorType::Gauss) {
						if (strcmp(name, "eta") == 0)
							eo_curr_lev->gp_eta = value;
						if (strcmp(name, "weight") == 0)
							eo_curr_lev->gp_weight = value;
					}
					else if (eo_curr_lev->type == detail::ElementOutputDescriptorType::Fiber) {
						if (strcmp(name, "yLoc") == 0)
							eo_curr_lev->fib_y = value;
						else if (strcmp(name, "zLoc") == 0)
							eo_curr_lev->fib_z = value;
						else if ((strcmp(name, "area") == 0) || (strcmp(name, "thickness") == 0))
							eo_curr_lev->fib_a = value;
					}
				}
				return 0;
			};

			int attr(const char *name, const char *value) {
				if (in_basis) {
					if (strcmp(name, "family") == 0) basis_info.family = value ? value : "";
					else if (strcmp(name, "topology") == 0) basis_info.topology = value ? value : "";
					else if (strcmp(name, "paramDomain") == 0) basis_info.param_domain = value ? value : "";
				}
				return 0;
			};
			int write(Vector &data) { return 0; };

			OPS_Stream& write(const char *s, int n) { return *this; };
			OPS_Stream& write(const unsigned char *s, int n) { return *this; };
			OPS_Stream& write(const signed char *s, int n) { return *this; };
			OPS_Stream& write(const void *s, int n) { return *this; };
			OPS_Stream& operator<<(char c) { return *this; };
			OPS_Stream& operator<<(unsigned char c) { return *this; };
			OPS_Stream& operator<<(signed char c) { return *this; };
			OPS_Stream& operator<<(const char *s) { return *this; };
			OPS_Stream& operator<<(const unsigned char *s) { return *this; };
			OPS_Stream& operator<<(const signed char *s) { return *this; };
			OPS_Stream& operator<<(const void *p) { return *this; };
			OPS_Stream& operator<<(int n) { return *this; };
			OPS_Stream& operator<<(unsigned int n) { return *this; };
			OPS_Stream& operator<<(long n) { return *this; };
			OPS_Stream& operator<<(unsigned long n) { return *this; };
			OPS_Stream& operator<<(short n) { return *this; };
			OPS_Stream& operator<<(unsigned short n) { return *this; };
			OPS_Stream& operator<<(bool b) { return *this; };
			OPS_Stream& operator<<(double n) { return *this; };
			OPS_Stream& operator<<(float n) { return *this; };

			int sendSelf(int commitTag, Channel &theChannel) { return 0; };
			int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) { return 0; };

		public:
			void finalizeSetResponse() {
				while (current_level > 0) {
					endTag();
				}
				/*
				this method was originally automatically called when a call to endTag determines
				the end of an ElementOutput. Some elements (ForceBeamColumn3d) do not call endTag after a
				tag("GaussPointOutput") in case of output from multiple sections, so it is called manually.
				*/
				/*
				do some workaround...
				*/
				descr->fixFloatingFiberOutput();
			}

		private:
			void ensureItemsOfUniformType(detail::element::OutputDescriptor *parent, detail::element::OutputDescriptor *child) {
				if (parent->items.size() > 0) {
					if (child->type != parent->items.back()->type) {
						/*opserr << "LadrunoRecorder Error: (detail::element::OutputDescriptor) "
							"Responses at the same level of the response tree must be of the same type.\n"
							"Expected: " << detail::ElementOutputDescriptorType::toString(parent->items.back()->type)
							<< ", given: " << detail::ElementOutputDescriptorType::toString(child->type) << "\n";
						exit(-1);*/
						// M.Petracca - due to a recent commit (08/10/2021)
						// this one can be converted from a fatal error to a silent-skip...
						error_code = ERROR_CODE_GENERIC;
						if ((child->type == detail::ElementOutputDescriptorType::Section) &&
							(parent->items.back()->type == detail::ElementOutputDescriptorType::Fiber)) {
							error_code = ERROR_CODE_SECTION_AFTER_FIBER;
						}
					}
				}
			}

		public:
			detail::element::OutputDescriptor *descr;
			int current_level;
			bool pending_close_tag;
			StreamErrorCode error_code;
		};

		class OutputResponse
		{
		public:
			OutputResponse()
				: response(0)
				, element(0)
			{}
			OutputResponse(Element *_elem, Response *_resp)
				: response(_resp)
				, element(_elem)
			{}
			Response *response;
			Element *element;
		};

		struct OutputResponseCollection
		{
			OutputResponseCollection()
				: is_new(true)
				, dir_name("")
				, initialized(false)
				, items() {}
			bool is_new;
			std::string dir_name;
			bool initialized;
			std::vector<OutputResponse> items;
		};

		struct OutputWithSameCustomIntRuleCollection
		{
			typedef std::map<OutputDescriptorHeader, OutputResponseCollection> collection_type;
			collection_type items;
		};

		struct OutputWithSameIntRuleCollection
		{
			typedef std::map<int, OutputWithSameCustomIntRuleCollection> collection_type;
			collection_type items;
		};

		struct OutputWithSameClassTagCollection
		{
			typedef std::map<ElementIntegrationRuleType::Enum, OutputWithSameIntRuleCollection> collection_type;
			collection_type items;
		};

		struct ResultRecorder
		{
			typedef std::map<int, OutputWithSameClassTagCollection> collection_type;
			ResultRecorder()
				: initialized(false)
				, result_request()
				, response_map()
			{}
			bool initialized;
			std::vector<std::string> result_request;
			collection_type response_map;
		};

		typedef std::vector<ResultRecorder> ResultRecorderCollection;

		/*************************************************************************************

		utilities for element mapping based on:
		class tag
		integration rule
		default or custom integration rule
		<element group>

		In this way all elements in <element group> share the same:
		1) number of nodes, 2) number and location of integration points

		**************************************************************************************/

		struct ElementIntegrationRule
		{
			ElementIntegrationRule()
				: int_rule_type(ElementIntegrationRuleType::CustomIntegrationRule)
			{}
			ElementIntegrationRule(ElementIntegrationRuleType::Enum _int_rule_type)
				: int_rule_type(_int_rule_type)
			{}
			inline bool operator < (const ElementIntegrationRule &other) const {
				const double rel_tol = 1.0e-5;
				if (int_rule_type < other.int_rule_type) return true;
				if (int_rule_type > other.int_rule_type) return false;
				if (custom_rule_dimension < other.custom_rule_dimension) return true;
				if (custom_rule_dimension > other.custom_rule_dimension) return false;
				if (x.size() < other.x.size()) return true;
				if (x.size() > other.x.size()) return false;
				for (size_t i = 0; i < x.size(); i++) {
					double tol = std::max(std::abs(x[i]), std::abs(other.x[i]))*rel_tol;
					if (utils::misc::lessThanWithTol(x[i], other.x[i], rel_tol)) return true;
					if (utils::misc::greaterThanWithTol(x[i], other.x[i], rel_tol)) return false;
					// continue with the loop
				}
				return false; // everything is equal
			}
			inline bool operator > (const ElementIntegrationRule &other) const {
				const double rel_tol = 1.0e-5;
				if (int_rule_type > other.int_rule_type) return true;
				if (int_rule_type < other.int_rule_type) return false;
				if (custom_rule_dimension > other.custom_rule_dimension) return true;
				if (custom_rule_dimension < other.custom_rule_dimension) return false;
				if (x.size() > other.x.size()) return true;
				if (x.size() < other.x.size()) return false;
				for (size_t i = 0; i < x.size(); i++) {
					double tol = std::max(std::abs(x[i]), std::abs(other.x[i]))*rel_tol;
					if (utils::misc::greaterThanWithTol(x[i], other.x[i], rel_tol)) return true;
					if (utils::misc::lessThanWithTol(x[i], other.x[i], rel_tol)) return false;
					// continue with the loop
				}
				return false; // everything is equal
			}
			ElementIntegrationRuleType::Enum int_rule_type;
			std::vector<double> x;
			int custom_rule_dimension = 1;
			// num_gp == 0  ⇒ legacy line rule: x holds 1 coord per point
			//                (GP_PARAM written as [x.size() x 1]).
			// num_gp >  0  ⇒ multi-dim parametric rule: x is row-major
			//                [num_gp x (x.size()/num_gp)] (e.g. barycentric
			//                triangle coords). Excluded from operator</> below
			//                because it is derivable from (x, dimension).
			int num_gp = 0;
			std::vector<double> w; // optional quadrature weights, one per gp
		};

		/* ----------------------------------------------------------------- *
		 * Standard (non-custom) Gauss rule -> parametric GP coords + weights *
		 *                                                                     *
		 * GP_PARAM is row-major [num_gp x ndir]. The row order MUST match the *
		 * element's own integration-point loop (the GP order is per element,  *
		 * NOT a standard tensor sweep) so that GP_PARAM[k]                     *
		 * pairs with result gauss_id k. ndir = parametric dimension (decoupled *
		 * from polynomial ORDER): line 1, quad/tri 2, hex/tet 3. tri/tet are   *
		 * barycentric free-coords (first ndir; last = 1 - sum). Returns false  *
		 * for rules not tabulated -> caller omits QUADRATURE for that bucket.  *
		 * All abscissae/weights verified against the canonical element source. *
		 * ----------------------------------------------------------------- */
		inline bool getStandardQuadrature(
			ElementIntegrationRuleType::Enum rule,
			std::vector<double>& gp_param, std::vector<double>& gp_weight,
			int& num_gp, int& ndir)
		{
			gp_param.clear(); gp_weight.clear(); num_gp = 0; ndir = 0;
			const double a = 0.5773502691896258;   // 1/sqrt(3)      (2-pt GL)
			const double b = 0.7745966692414834;   // sqrt(0.6)      (3-pt GL)
			const double w5 = 5.0 / 9.0, w8 = 8.0 / 9.0;
			switch (rule) {
			// ---- line, domain [-1,1] -------------------------------------
			case ElementIntegrationRuleType::Line_GaussLegendre_1:
				ndir = 1; num_gp = 1; gp_param = { 0.0 }; gp_weight = { 2.0 }; return true;
			case ElementIntegrationRuleType::Line_GaussLegendre_2:
				ndir = 1; num_gp = 2; gp_param = { -a, a }; gp_weight = { 1.0, 1.0 }; return true;
			case ElementIntegrationRuleType::Line_GaussLegendre_3:
				ndir = 1; num_gp = 3; gp_param = { -b, 0.0, b }; gp_weight = { w5, w8, w5 }; return true;
			// ---- quad, domain [-1,1]^2 -----------------------------------
			case ElementIntegrationRuleType::Quadrilateral_GaussLegendre_1:
				ndir = 2; num_gp = 1; gp_param = { 0.0, 0.0 }; gp_weight = { 4.0 }; return true;
			case ElementIntegrationRuleType::Quadrilateral_GaussLegendre_2: // FourNodeQuad: CCW
				ndir = 2; num_gp = 4;
				gp_param = { -a,-a,  a,-a,  a, a,  -a, a };
				gp_weight = { 1.0, 1.0, 1.0, 1.0 }; return true;
			case ElementIntegrationRuleType::Quadrilateral_GaussLegendre_3: { // NineNodeQuad: corners,edges,center
				ndir = 2; num_gp = 9;
				gp_param = { -b,-b,  b,-b,  b, b,  -b, b,    // 4 CCW corners
				              0,-b,  b, 0,  0, b,  -b, 0,    // 4 CCW edge midpoints
				              0, 0 };                        // center
				const double wc = w5 * w5, wd = w8 * w5, w0 = w8 * w8;
				gp_weight = { wc, wc, wc, wc,  wd, wd, wd, wd,  w0 }; return true; }
			// ---- tri, barycentric (ndir=2) -------------------------------
			case ElementIntegrationRuleType::Triangle_GaussLegendre_1: // Tri31
				ndir = 2; num_gp = 1; gp_param = { 1.0 / 3.0, 1.0 / 3.0 }; gp_weight = { 0.5 }; return true;
			// ---- tet, barycentric (ndir=3) -------------------------------
			case ElementIntegrationRuleType::Tetrahedron_GaussLegendre_1: // FourNodeTetrahedron: 1 GP (not 4)
				ndir = 3; num_gp = 1; gp_param = { 0.25, 0.25, 0.25 }; gp_weight = { 1.0 / 6.0 }; return true;
			case ElementIntegrationRuleType::Tetrahedron_GaussLegendre_2: { // TenNodeTetrahedron: 4-pt rule
				// alpha=(5+3√5)/20, beta=(5−√5)/20; element GP order (ζ0,ζ1,ζ2) is
				// the (α,β,β),(β,α,β),(β,β,α),(β,β,β) cycle (TenNodeTetrahedron.cpp
				// gaussPoint[d]=sg[abs(d-k)] loop). w=1/24 each.
				const double al = 0.5854101966249685, be = 0.1381966011250105;
				ndir = 3; num_gp = 4;
				gp_param = { al, be, be,  be, al, be,  be, be, al,  be, be, be };
				const double tw = 1.0 / 24.0;
				gp_weight = { tw, tw, tw, tw }; return true; }
			// ---- hex, domain [-1,1]^3 ------------------------------------
			case ElementIntegrationRuleType::Hexahedron_GaussLegendre_1:
				ndir = 3; num_gp = 1; gp_param = { 0.0, 0.0, 0.0 }; gp_weight = { 8.0 }; return true;
			case ElementIntegrationRuleType::Hexahedron_GaussLegendre_2: { // Brick: i,j,k lexicographic (z fastest)
				ndir = 3; num_gp = 8;
				gp_param = { -a,-a,-a,  -a,-a, a,  -a, a,-a,  -a, a, a,
				              a,-a,-a,   a,-a, a,   a, a,-a,   a, a, a };
				gp_weight = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }; return true; }
			case ElementIntegrationRuleType::Hexahedron_GaussLegendre_3: { // Twenty_Node_Brick: 27-pt (brcshl)
				// GP order follows the serendipity node pattern (NOT i,j,k tensor):
				// GP[L] = b·(2·RA[L], 2·SA[L], 2·TA[L]) with the RA/SA/TA arrays from
				// shp3dv.cpp brcshl; 8 corners, 12 edge-mids, 6 face-centers, centroid.
				// This is the element's materialPointers[L] order, so GP_PARAM[L] pairs
				// with result gauss_id L. Weights = tensor {5/9,8/9} products (Σ=8).
				ndir = 3; num_gp = 27;
				gp_param = {
					-b,-b,-b,  b,-b,-b,  b, b,-b, -b, b,-b,    // 0-3  lower corners
					-b,-b, b,  b,-b, b,  b, b, b, -b, b, b,    // 4-7  upper corners
					 0,-b,-b,  b, 0,-b,  0, b,-b, -b, 0,-b,    // 8-11 lower edge mids
					 0,-b, b,  b, 0, b,  0, b, b, -b, 0, b,    // 12-15 upper edge mids
					-b,-b, 0,  b,-b, 0,  b, b, 0, -b, b, 0,    // 16-19 vertical edge mids
					 b, 0, 0,  0, b, 0,  0, 0, b, -b, 0, 0,    // 20-23 faces +r,+s,+t,-r
					 0,-b, 0,  0, 0,-b,  0, 0, 0 };            // 24-26 faces -s,-t, centroid
				const double q5 = 5.0 / 9.0, q8 = 8.0 / 9.0;
				const double wc = q5 * q5 * q5, wd = q5 * q5 * q8;
				const double wf = q5 * q8 * q8, w0 = q8 * q8 * q8;
				gp_weight = {
					wc, wc, wc, wc, wc, wc, wc, wc,            // 8 corners
					wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, // 12 edges
					wf, wf, wf, wf, wf, wf,                    // 6 faces
					w0 };                                      // centroid
				return true; }
			default:
				return false; // Tri_GL_2/2B/2C not yet tabulated
			}
		}

		/* ----------------------------------------------------------------- *
		 * Write-side reference-config global Gauss-point coordinates for ONE *
		 * element:   x(ξ_g) = Σ_i N_i(ξ_g) · X_i                              *
		 *                                                                     *
		 * The shape functions N_i are evaluated in the element's OWN node     *
		 * order (getExternalNodes()); every basis below was verified against *
		 * its canonical element source (FourNodeQuad / NineNodeQuad / Brick + *
		 * shp3d / Tri31 / FourNodeTetrahedron). GLOBAL_GP_COORDS both lets     *
		 * legacy elements be consumed with zero reader-side basis math AND     *
		 * catches GP-ordering/basis bugs at write time when paired with the    *
		 * round-trip check x(GP_PARAM[k]) vs GLOBAL_GP_COORDS[k].              *
		 *                                                                     *
		 * Inputs:                                                             *
		 *   geom      : element geometry enum (selects the basis)             *
		 *   num_nodes : number of element nodes                               *
		 *   node_xyz  : [num_nodes][ndim] reference node coords, in           *
		 *               getExternalNodes() order                              *
		 *   gp_param  : row-major [num_gp × ndir] natural coords (the rule,    *
		 *               in the element's GP order)                            *
		 *   ndim      : spatial dimension (1/2/3)                             *
		 * Output (on success): out is row-major [num_gp × ndim].             *
		 * Returns false for topologies whose write-side basis is not yet      *
		 *   implemented (caller then omits GLOBAL_GP_COORDS for that bucket;   *
		 *   the reader treats the dataset as optional). Coverage: line 2N,    *
		 *   quad 4N/9N, tri 3N, tet 4N, hex 8N.                               *
		 * ----------------------------------------------------------------- */
		inline bool computeGlobalGP(
			ElementGeometryType::Enum geom, int num_nodes,
			const std::vector<std::vector<double> >& node_xyz, int ndim,
			const std::vector<double>& gp_param, int num_gp, int ndir,
			std::vector<double>& out)
		{
			out.clear();
			if (ndim <= 0 || num_gp <= 0 || ndir <= 0) return false;
			if ((int)node_xyz.size() < num_nodes) return false;

			// Evaluate the shape functions for ONE Gauss point into N[num_nodes].
			// natc points at this point's ndir natural coords. Returns false if the
			// (geom, num_nodes) pair is not a supported write-side basis.
			struct ShapeEval {
				static bool eval(ElementGeometryType::Enum g, int nn,
					const double* natc, int nd, std::vector<double>& N)
				{
					N.assign(nn, 0.0);
					switch (g) {
					// ---- line, [-1,1] : 2-node linear -----------------------
					case ElementGeometryType::Line_2N: {
						if (nn < 2 || nd < 1) return false;
						double xi = natc[0];
						N[0] = 0.5 * (1.0 - xi);
						N[1] = 0.5 * (1.0 + xi);
						return true; }
					// ---- tri, barycentric (Tri31: N0=ξ,N1=η,N2=1-ξ-η) -------
					case ElementGeometryType::Triangle_3N: {
						if (nn < 3 || nd < 2) return false;
						double xi = natc[0], eta = natc[1];
						N[0] = xi; N[1] = eta; N[2] = 1.0 - xi - eta;
						return true; }
					// ---- quad, [-1,1]^2 (FourNodeQuad CCW bilinear) ---------
					case ElementGeometryType::Quadrilateral_4N: {
						if (nn < 4 || nd < 2) return false;
						double s = natc[0], t = natc[1];
						N[0] = 0.25 * (1.0 - s) * (1.0 - t);
						N[1] = 0.25 * (1.0 + s) * (1.0 - t);
						N[2] = 0.25 * (1.0 + s) * (1.0 + t);
						N[3] = 0.25 * (1.0 - s) * (1.0 + t);
						return true; }
					// ---- quad, [-1,1]^2 (NineNodeQuad biquadratic Lagrange) -
					// node order: 4 CCW corners, 4 CCW edge-mids, center.
					case ElementGeometryType::Quadrilateral_9N: {
						if (nn < 9 || nd < 2) return false;
						double s = natc[0], t = natc[1];
						N[0] =  (1.0 - s) * (1.0 - t) * s * t / 4.0;
						N[1] = -(1.0 + s) * (1.0 - t) * s * t / 4.0;
						N[2] =  (1.0 + s) * (1.0 + t) * s * t / 4.0;
						N[3] = -(1.0 - s) * (1.0 + t) * s * t / 4.0;
						N[4] = -(1.0 - s * s) * (1.0 - t) * t / 2.0;
						N[5] =  (1.0 + s) * (1.0 - t * t) * s / 2.0;
						N[6] =  (1.0 - s * s) * (1.0 + t) * t / 2.0;
						N[7] = -(1.0 - s) * (1.0 - t * t) * s / 2.0;
						N[8] =  (1.0 - s * s) * (1.0 - t * t);
						return true; }
					// ---- tet, barycentric (N0=r,N1=s,N2=t,N3=1-r-s-t) -------
					case ElementGeometryType::Tetrahedron_4N: {
						if (nn < 4 || nd < 3) return false;
						double r = natc[0], s = natc[1], u = natc[2];
						N[0] = r; N[1] = s; N[2] = u; N[3] = 1.0 - r - s - u;
						return true; }
					// ---- tet, quadratic (TenNodeTetrahedron, node 8/9 swapped) ---
					// ζ=(r,s,u); L4=1-r-s-u. Edge nodes: 4:0-1, 5:1-2, 6:2-0,
					// 7:0-3, 8:2-3, 9:1-3 (TenNodeTetrahedron.cpp shp3d).
					case ElementGeometryType::Tetrahedron_10N: {
						if (nn < 10 || nd < 3) return false;
						double r = natc[0], s = natc[1], u = natc[2];
						double L4 = 1.0 - r - s - u;
						N[0] = r * (2.0 * r - 1.0);
						N[1] = s * (2.0 * s - 1.0);
						N[2] = u * (2.0 * u - 1.0);
						N[3] = L4 * (2.0 * L4 - 1.0);
						N[4] = 4.0 * r * s;
						N[5] = 4.0 * s * u;
						N[6] = 4.0 * u * r;
						N[7] = 4.0 * r * L4;
						N[8] = 4.0 * u * L4;
						N[9] = 4.0 * s * L4;
						return true; }
					// ---- hex, [-1,1]^3 (Brick/shp3d trilinear) --------------
					// node order: bottom CCW (0:---,1:+--,2:++-,3:-+-),
					//             top    CCW (4:--+,5:+-+,6:+++,7:-++).
					case ElementGeometryType::Hexahedron_8N: {
						if (nn < 8 || nd < 3) return false;
						double r = natc[0], s = natc[1], u = natc[2];
						static const double sgn[8][3] = {
							{-1,-1,-1}, { 1,-1,-1}, { 1, 1,-1}, {-1, 1,-1},
							{-1,-1, 1}, { 1,-1, 1}, { 1, 1, 1}, {-1, 1, 1} };
						for (int i = 0; i < 8; ++i)
							N[i] = 0.125 * (1.0 + sgn[i][0] * r)
							             * (1.0 + sgn[i][1] * s)
							             * (1.0 + sgn[i][2] * u);
						return true; }
					// ---- hex, 20-node serendipity (Twenty_Node_Brick) -------
					// node order: 8 corners, then edge-mids 8-11 (lower edges
					// 1-2,2-3,3-4,4-1), 12-15 (upper), 16-19 (vertical 1-5..4-8)
					// per shp3dv.cpp brcshl. Standard serendipity basis (exact for
					// straight-sided bricks; reference-config GLOBAL_GP_COORDS).
					case ElementGeometryType::Hexahedron_20N: {
						if (nn < 20 || nd < 3) return false;
						double xyz[3] = { natc[0], natc[1], natc[2] };
						static const double nc[20][3] = {
							{-1,-1,-1}, { 1,-1,-1}, { 1, 1,-1}, {-1, 1,-1},
							{-1,-1, 1}, { 1,-1, 1}, { 1, 1, 1}, {-1, 1, 1},
							{ 0,-1,-1}, { 1, 0,-1}, { 0, 1,-1}, {-1, 0,-1},
							{ 0,-1, 1}, { 1, 0, 1}, { 0, 1, 1}, {-1, 0, 1},
							{-1,-1, 0}, { 1,-1, 0}, { 1, 1, 0}, {-1, 1, 0} };
						for (int i = 0; i < 20; ++i) {
							const double* c = nc[i];
							int zeros = 0, zdir = -1;
							for (int d = 0; d < 3; ++d)
								if (c[d] == 0.0) { ++zeros; zdir = d; }
							if (zeros == 0) { // corner
								double lin = xyz[0]*c[0] + xyz[1]*c[1] + xyz[2]*c[2] - 2.0;
								N[i] = 0.125 * (1.0 + xyz[0]*c[0]) * (1.0 + xyz[1]*c[1])
								             * (1.0 + xyz[2]*c[2]) * lin;
							} else { // edge mid: (1-u^2) along the zero direction
								double prod = 0.25 * (1.0 - xyz[zdir]*xyz[zdir]);
								for (int d = 0; d < 3; ++d)
									if (d != zdir) prod *= (1.0 + xyz[d]*c[d]);
								N[i] = prod;
							}
						}
						return true; }
					default:
						return false;
					}
				}
			};

			out.assign((size_t)num_gp * (size_t)ndim, 0.0);
			std::vector<double> N;
			for (int g = 0; g < num_gp; ++g) {
				if (!ShapeEval::eval(geom, num_nodes, &gp_param[(size_t)g * (size_t)ndir], ndir, N))
					return false;
				for (int i = 0; i < num_nodes; ++i) {
					double Ni = N[(size_t)i];
					const std::vector<double>& Xi = node_xyz[(size_t)i];
					for (int d = 0; d < ndim; ++d)
						out[(size_t)g * (size_t)ndim + (size_t)d] += Ni * Xi[(size_t)d];
				}
			}
			return true;
		}

		struct ElementWithSameCustomIntRuleCollection
		{
			typedef std::vector<Element*> collection_type;

			ElementWithSameCustomIntRuleCollection()
				: is_new(true), custom_int_rule_index(0), name(""), items()
			{}

			bool is_new;
			int custom_int_rule_index;
			std::string name;
			std::vector<Element*> items;
		};

		struct ElementWithSameIntRuleCollection
		{
			typedef std::map<int, ElementWithSameCustomIntRuleCollection> submap_type;

			ElementWithSameIntRuleCollection()
				: is_new(true), int_rule_type(ElementIntegrationRuleType::CustomIntegrationRule), items()
			{}

			bool is_new;
			ElementIntegrationRuleType::Enum int_rule_type;
			std::map<int, ElementWithSameCustomIntRuleCollection> items;
		};

		struct ElementWithSameClassTagCollection
		{
			typedef std::map<ElementIntegrationRuleType::Enum, ElementWithSameIntRuleCollection> submap_type;

			ElementWithSameClassTagCollection()
				: is_new(true), class_tag(0), class_name("unknown"), num_nodes(0), geom_type(ElementGeometryType::Custom), items() {}

			bool is_new;
			int class_tag;
			std::string class_name;
			int num_nodes;
			ElementGeometryType::Enum geom_type;
			BasisInfo basis_info; // from the optional "basisInfo" probe (one per class tag)
			std::map<detail::ElementIntegrationRuleType::Enum, ElementWithSameIntRuleCollection> items;
		};

		struct ElementCollection
		{
			typedef std::map<int, ElementWithSameClassTagCollection> submap_type;

			ElementCollection()
				: registered_custom_rules(), items()
			{}

			void mapElements(Domain *d, bool has_region, const std::vector<int> &subset) {
				/*
				utilties
				*/
				auto lam_get_num_ext_nodes = [](Element* elem) {
					switch (elem->getClassTag()) {
					case ELE_TAG_SFI_MVLEM_3D: return 4;
					case ELE_TAG_SFI_MVLEM: return 2;
					default: return elem->getNumExternalNodes();
					}
				};
				/*
				clear previous mappings
				*/
				registered_custom_rules.clear();
				items.clear();
				/*
				quick return
				*/
				if (d == 0) return;
				/*
				auxiliary map to register custom rules
				*/
				std::map<ElementIntegrationRule, int> aux_map_custom_rules;
				/*
				loop over all elements in the domain
				*/
				size_t subset_elem_counter(0);
				ElementIter* element_iter = &(d->getElements());
				Element* current_element = 0;
				//while ((current_element = (*element_iter)()) != 0) {
				while (true) {
					/*
					get next element
					*/
					if (has_region) {
						if (subset_elem_counter == subset.size())
							break;
						current_element = d->getElement(subset[subset_elem_counter++]);
						if (current_element == 0)
							continue; // skip null and go to next iteration
					}
					else {
						current_element = (*element_iter)();
						if (current_element == 0) 
							break; 
					}
					/*
					get class tag, geometry and integration rule type
					*/
					int elem_type = current_element->getClassTag();
					/*
					skip element classes that are not recorded
					*/
					if (elem_type == ELE_TAG_Subdomain ||
						elem_type == ELE_TAG_ASDEmbeddedNodeElement)
					{
						continue;
					}
					ElementGeometryType::Enum geom_type;
					ElementIntegrationRuleType::Enum int_rule_type;
					int custom_rule_dimension;
					getGeometryAndIntRuleByClassTag(elem_type, geom_type, int_rule_type, custom_rule_dimension);
					/*
					map by class tag
					*/
					ElementWithSameClassTagCollection &elem_coll_by_tag = items[elem_type];
					if (elem_coll_by_tag.is_new) {
						elem_coll_by_tag.class_tag = elem_type;
						elem_coll_by_tag.class_name = current_element->getClassType();
						elem_coll_by_tag.num_nodes = lam_get_num_ext_nodes(current_element);
						elem_coll_by_tag.geom_type = geom_type;
						elem_coll_by_tag.basis_info = getElementBasisInfo(current_element);
						elem_coll_by_tag.is_new = false;
					}
					/*
					make sure that every element with the same tag have the same number of nodes
					*/
					if (lam_get_num_ext_nodes(current_element) != elem_coll_by_tag.num_nodes) {
						opserr << "LadrunoRecorder Error while mapping elements: elements with different number of nodes "
							"exist within the same class tag. This is not supported\n";
						continue;  // skip the inconsistent element; the run is not aborted
					}
					/*
					create the integration rule
					*/
					ElementIntegrationRule int_rule(int_rule_type);
					if (int_rule_type == ElementIntegrationRuleType::CustomIntegrationRule) {
						// set the declared dimension FIRST so getCustomGaussPointLocations
						// can pick the multi-dim (e.g. barycentric) extraction path.
						int_rule.custom_rule_dimension = custom_rule_dimension;
						getCustomGaussPointLocations(current_element, int_rule);
					}
					/*
					if this is a custom rule, register it
					*/
					int custom_int_rule_index = 0;
					if (int_rule_type == ElementIntegrationRuleType::CustomIntegrationRule) {
						std::map<ElementIntegrationRule, int>::iterator reg_int_rule_iter = aux_map_custom_rules.find(int_rule);
						if (reg_int_rule_iter == aux_map_custom_rules.end()) {
							// new rule, define a new index starting from 1
							custom_int_rule_index = (int)aux_map_custom_rules.size() + 1;
							aux_map_custom_rules[int_rule] = custom_int_rule_index;
						}
						else {
							custom_int_rule_index = reg_int_rule_iter->second;
						}
					}
					/*
					map by integration rule
					*/
					ElementWithSameIntRuleCollection &elem_coll_by_rule = elem_coll_by_tag.items[int_rule_type];
					if (elem_coll_by_rule.is_new) {
						elem_coll_by_rule.int_rule_type = int_rule_type;
						elem_coll_by_rule.is_new = false;
					}
					/*
					map by custom integration rule index
					*/
					ElementWithSameCustomIntRuleCollection &elem_coll_by_custom_rule = elem_coll_by_rule.items[custom_int_rule_index];
					if (elem_coll_by_custom_rule.is_new) {
						elem_coll_by_custom_rule.custom_int_rule_index = custom_int_rule_index;
						elem_coll_by_custom_rule.is_new = false;
					}
					/*
					finally add this element
					*/
					elem_coll_by_custom_rule.items.push_back(current_element);
				}
				/*
				now fill the custom integration rule map
				*/
				for (std::map<ElementIntegrationRule, int>::iterator it = aux_map_custom_rules.begin();
					it != aux_map_custom_rules.end(); ++it) {
					registered_custom_rules[it->second] = it->first;
				}
			}

			void getGeometryAndIntRuleByClassTag(
				int elem_class_tag,
				ElementGeometryType::Enum &geom_type,
				ElementIntegrationRuleType::Enum &int_type,
				int &custom_rule_dimension) {
				/*
				set default values. custom geometry (i.e. point cloud)
				and no integration rule
				*/
				geom_type = ElementGeometryType::Custom;
				int_type = ElementIntegrationRuleType::NoIntegrationRule;
				custom_rule_dimension = 1;
				/*
				2-node line with 1 gp
				*/
				if (
					// ./adapter actuators
					elem_class_tag == ELE_TAG_Actuator ||
					elem_class_tag == ELE_TAG_ActuatorCorot ||
					// ./absorbentBoundaries
					elem_class_tag == ELE_TAG_FSIInterfaceElement2D ||
					elem_class_tag == ELE_TAG_FSIFluidBoundaryElement2D ||
					// ./truss
					elem_class_tag == ELE_TAG_Truss ||
					elem_class_tag == ELE_TAG_Truss2 ||
					elem_class_tag == ELE_TAG_TrussSection ||
					elem_class_tag == ELE_TAG_CorotTruss ||
					elem_class_tag == ELE_TAG_CorotTruss2 ||
					elem_class_tag == ELE_TAG_CorotTrussSection ||
					elem_class_tag == ELE_TAG_InertiaTruss ||
					// ./zeroLength
					elem_class_tag == ELE_TAG_ZeroLength ||
					elem_class_tag == ELE_TAG_ZeroLengthSection ||
					elem_class_tag == ELE_TAG_ZeroLengthND ||
					elem_class_tag == ELE_TAG_CoupledZeroLength ||
					elem_class_tag == ELE_TAG_ZeroLengthRocking ||
					elem_class_tag == ELE_TAG_ZeroLengthContact2D ||
					elem_class_tag == ELE_TAG_ZeroLengthContact3D ||
					elem_class_tag == ELE_TAG_ZeroLengthContactASDimplex ||
					elem_class_tag == ELE_Tag_ZeroLengthImpact3D ||
					// ./twoNodeLink
					elem_class_tag == ELE_TAG_TwoNodeLinkSection ||
					// ./elasticBeamColumn
					elem_class_tag == ELE_TAG_ElasticBeam2d ||
					elem_class_tag == ELE_TAG_ElasticBeam3d ||
					elem_class_tag == ELE_TAG_ElasticTimoshenkoBeam2d ||
					elem_class_tag == ELE_TAG_ElasticTimoshenkoBeam3d ||
					elem_class_tag == ELE_TAG_ModElasticBeam2d ||
					// .elastomericBearing
					elem_class_tag == ELE_TAG_ElastomericBearingBoucWen2d ||
					elem_class_tag == ELE_TAG_ElastomericBearingBoucWen3d ||
					elem_class_tag == ELE_TAG_ElastomericBearingBoucWenMod3d ||
					elem_class_tag == ELE_TAG_ElastomericBearingPlasticity2d ||
					elem_class_tag == ELE_TAG_ElastomericBearingPlasticity3d ||
					elem_class_tag == ELE_TAG_ElastomericBearingUFRP2d ||
					elem_class_tag == ELE_TAG_ElastomericBearingUFRP3d ||
					elem_class_tag == ELE_TAG_ElastomericX ||
					elem_class_tag == ELE_TAG_HDR ||
					elem_class_tag == ELE_TAG_LeadRubberX ||
					// ./ulBeamColumn
					/*warning: these two could go to the beam with custom integration rule, but they do not define everything properly! check in future versions */
					elem_class_tag == Ele_TAG_Elastic2dGNL || /*warning: no integrationPoints, no resp for all sections, no section response*/
					elem_class_tag == TAG_InelasticYS2DGNL /*warning: no integrationPoints, no resp for all sections, no section response*/
														   // todo... others
					) {
					geom_type = ElementGeometryType::Line_2N;
					int_type = ElementIntegrationRuleType::Line_GaussLegendre_1;
				}
				/*
				2-node beams with custom number of gp
				*/
				else if (
					// ./dispBeamColumn
					elem_class_tag == ELE_TAG_DispBeamColumn2d || /*warning: no integrationPoints*/
					elem_class_tag == ELE_TAG_DispBeamColumn3d || /*warning: no integrationPoints*/
					elem_class_tag == ELE_TAG_DispBeamColumn2dWithSensitivity || /*warning: no integrationPoints, no resp for all sections*/
					elem_class_tag == ELE_TAG_DispBeamColumn3dWithSensitivity || /*warning: no integrationPoints, no resp for all sections*/
																				 // ./dispBeamColumnInt
					elem_class_tag == ELE_TAG_DispBeamColumn2dInt || /*warning: no integrationPoints, no resp for all sections*/
					elem_class_tag == ELE_TAG_DispBeamColumn2dThermal || /*warning: no integrationPoints, no resp for all sections*/
																		 // ./forceBeamColumn
					elem_class_tag == ELE_TAG_ElasticForceBeamColumn2d || /*warning: no resp for all sections*/
					elem_class_tag == ELE_TAG_ElasticForceBeamColumn3d || /*warning: no resp for all sections*/
					elem_class_tag == ELE_TAG_ElasticForceBeamColumnWarping2d || /*warning: no resp for all sections*/
					elem_class_tag == ELE_TAG_ForceBeamColumn2d || /* <- OK! this one defines everything ! good job*/
					elem_class_tag == ELE_TAG_ForceBeamColumn3d || /* <- OK! this one defines everything ! good job*/
					elem_class_tag == ELE_TAG_ForceBeamColumnCBDI2d || /* <- OK! this one defines everything ! good job*/
					elem_class_tag == ELE_TAG_ForceBeamColumnWarping2d || /* <- OK! this one defines everything ! good job*/
					// ./mixedBeamColumn
					elem_class_tag == ELE_TAG_MixedBeamColumn2d ||
					elem_class_tag == ELE_TAG_MixedBeamColumn3d ||
					// ./mvlem (2D)
					elem_class_tag == ELE_TAG_MVLEM ||
					elem_class_tag == ELE_TAG_SFI_MVLEM
					) {
					geom_type = ElementGeometryType::Line_2N;
					int_type = ElementIntegrationRuleType::CustomIntegrationRule;
				}
				/*
				3-node triangle with 1 gp
				*/
				else if (
					// ./triangle
					elem_class_tag == ELE_TAG_Tri31
					) {
					geom_type = ElementGeometryType::Triangle_3N;
					int_type = ElementIntegrationRuleType::Triangle_GaussLegendre_1;
				}
				/*
				3-node triangle with 3 gp
				*/
				else if (
					// ./shell
					elem_class_tag == ELE_TAG_ASDShellT3
					) {
					geom_type = ElementGeometryType::Triangle_3N;
					int_type = ElementIntegrationRuleType::Triangle_GaussLegendre_2B;
				}
				/*
				3-node triangle with 4 gp
				*/
				else if (
					// ./shell
					elem_class_tag == ELE_TAG_ShellDKGT ||
					elem_class_tag == ELE_TAG_ShellNLDKGT
					) {
					geom_type = ElementGeometryType::Triangle_3N;
					int_type = ElementIntegrationRuleType::Triangle_GaussLegendre_2C;
				}
				/*
				4-node quadrilateral with 1 gp
				*/
				else if (
					// ./UWelements
					elem_class_tag == ELE_TAG_SSPquad ||
					elem_class_tag == ELE_TAG_SSPquadUP ||
					// ./absorbentBoundaries
					elem_class_tag == ELE_TAG_ASDAbsorbingBoundary2D ||
					elem_class_tag == ELE_TAG_FSIFluidElement2D
					)
				{
					geom_type = ElementGeometryType::Quadrilateral_4N;
					int_type = ElementIntegrationRuleType::Quadrilateral_GaussLegendre_1;
				}
				/*
				4-node quadrilateral with 2x2 gp
				*/
				else if (
					// ./quad
					elem_class_tag == ELE_TAG_ConstantPressureVolumeQuad ||
					elem_class_tag == ELE_TAG_EnhancedQuad ||
					elem_class_tag == ELE_TAG_FourNodeQuad ||
					elem_class_tag == ELE_TAG_FourNodeQuad3d ||
					elem_class_tag == ELE_TAG_FourNodeQuadWithSensitivity ||
					// ./shell
					elem_class_tag == ELE_TAG_ShellDKGQ ||
					elem_class_tag == ELE_TAG_ShellNLDKGQ ||
					elem_class_tag == ELE_TAG_ShellMITC4 ||
					elem_class_tag == ELE_TAG_ShellMITC4Thermal ||
					elem_class_tag == ELE_TAG_ASDShellQ4 ||
					// ./up
					elem_class_tag == ELE_TAG_BBarFourNodeQuadUP ||
					elem_class_tag == ELE_TAG_FourNodeQuadUP
					) {
					geom_type = ElementGeometryType::Quadrilateral_4N;
					int_type = ElementIntegrationRuleType::Quadrilateral_GaussLegendre_2;
				}
				/*
				4-node quadrilateral cohesive with custom rule
				*/
				else if (
					// ./mvlem
					elem_class_tag == ELE_TAG_MVLEM_3D ||
					elem_class_tag == ELE_TAG_SFI_MVLEM_3D ||
					elem_class_tag == ELE_TAG_E_SFI_MVLEM_3D
					) {
					geom_type = ElementGeometryType::Quadrilateral_CohesiveBand_4N;
					int_type = ElementIntegrationRuleType::CustomIntegrationRule;
					custom_rule_dimension = 2;
				}
				/*
				9-node quadrilateral with 3x3 gp
				*/
				else if (
					// ./quad
					elem_class_tag == ELE_TAG_NineNodeMixedQuad ||
					elem_class_tag == ELE_TAG_NineNodeQuad ||
					// ./shell
					elem_class_tag == ELE_TAG_ShellMITC9 ||
					// ./up
					elem_class_tag == ELE_TAG_Nine_Four_Node_QuadUP
					) {
					geom_type = ElementGeometryType::Quadrilateral_9N;
					int_type = ElementIntegrationRuleType::Quadrilateral_GaussLegendre_3;
				}
				/*
				4-node tetrahedron with 4 gp
				*/
				else if (
					// ./tetrahedron
					elem_class_tag == ELE_TAG_FourNodeTetrahedron
					)
				{
					geom_type = ElementGeometryType::Tetrahedron_4N;
					int_type = ElementIntegrationRuleType::Tetrahedron_GaussLegendre_1;
				}
				/*
				10-node tetrahedron with 1x1x1 gp
				*/
				else if (
					// ./tetrahedron
					elem_class_tag == ELE_TAG_TenNodeTetrahedron
					)
				{
					geom_type = ElementGeometryType::Tetrahedron_10N;
					int_type = ElementIntegrationRuleType::Tetrahedron_GaussLegendre_2;
				}
				/*
				8-node hexahedron with 1x1x1 gp
				*/
				else if (
					// ./UWelements
					elem_class_tag == ELE_TAG_SSPbrick ||
					elem_class_tag == ELE_TAG_SSPbrickUP ||
					// ./absorbentBoundaries
					elem_class_tag == ELE_TAG_ASDAbsorbingBoundary3D
					)
				{
					geom_type = ElementGeometryType::Hexahedron_8N;
					int_type = ElementIntegrationRuleType::Hexahedron_GaussLegendre_1;
				}
				/*
				8-node hexahedron with 2x2x2 gp
				*/
				else if (
					// ./brick
					elem_class_tag == ELE_TAG_BbarBrick ||
					elem_class_tag == ELE_TAG_BbarBrickWithSensitivity ||
					elem_class_tag == ELE_TAG_Brick ||
					// ./up
					elem_class_tag == ELE_TAG_BBarBrickUP ||
					elem_class_tag == ELE_TAG_BrickUP ||
					// ./XMUelements
					elem_class_tag == ELE_TAG_AC3D8HexWithSensitivity
					)
				{
					geom_type = ElementGeometryType::Hexahedron_8N;
					int_type = ElementIntegrationRuleType::Hexahedron_GaussLegendre_2;
				}
				/*
				20-node hexahedron with 3x3x3 gp
				*/
				else if (
					// ./brick
					elem_class_tag == ELE_TAG_Twenty_Node_Brick ||
					// ./up
					elem_class_tag == ELE_TAG_Twenty_Eight_Node_BrickUP
					)
				{
					geom_type = ElementGeometryType::Hexahedron_20N;
					int_type = ElementIntegrationRuleType::Hexahedron_GaussLegendre_3;
				}
			}

			// Probe an element's OPTIONAL "basisInfo" self-description.
			// Lagrange elements don't answer, in
			// which case the returned BasisInfo has valid == false and callers
			// fall back to geometry-derived guesses.
			BasisInfo getElementBasisInfo(Element *elem) {
				BasisInfo info;
				std::string request = "basisInfo";
				int argc = 1;
				const char **argv = new const char*[argc];
				argv[0] = request.c_str();
				OutputDescriptor eo_descriptor;
				OutputDescriptorStream eo_stream(&eo_descriptor);
				Response *eo_response = elem->setResponse(argv, argc, eo_stream);
				eo_stream.finalizeSetResponse();
				info = eo_stream.basis_info;
				if (eo_response)
					delete eo_response;
				delete[] argv;
				return info;
			}

			void getCustomGaussPointLocations(Element *elem, ElementIntegrationRule &rule) {
				/*
				clear any existing locations
				*/
				rule.x.clear();
				rule.w.clear();
				rule.num_gp = 0;
				/*
				Multi-dimensional parametric rules (custom_rule_dimension >= 2):
				the element returns "integrationPoints" as a Matrix [nGP x dim]
				of natural coordinates (e.g. a Bezier triangle: 3 GP x 2
				barycentric area coords). Store them verbatim, row-major, WITHOUT the 1D
				[-1,1] line normalization applied to line elements below; also
				grab matching "integrationWeights" if available. Elements that
				don't answer with a Matrix (e.g. MVLEM) fall through to the
				legacy 1D strategies, preserving their previous behavior.
				*/
				if (rule.custom_rule_dimension >= 2) {
					bool done = false;
					std::string request = "integrationPoints";
					int argc = 1;
					const char **argv = new const char*[argc];
					argv[0] = request.c_str();
					OutputDescriptor eo_descriptor;
					OutputDescriptorStream eo_stream(&eo_descriptor);
					Response *eo_response = elem->setResponse(argv, argc, eo_stream);
					eo_stream.finalizeSetResponse();
					if (eo_response) {
						eo_response->getResponse();
						Information &eo_info = eo_response->getInformation();
						if (eo_info.theType == MatrixType && eo_info.theMatrix != 0) {
							const Matrix &gp = *eo_info.theMatrix;
							int ngp = gp.noRows();
							int ncomp = gp.noCols();
							if (ngp > 0 && ncomp > 0) {
								rule.x.resize((size_t)ngp * (size_t)ncomp);
								for (int i = 0; i < ngp; i++)
									for (int j = 0; j < ncomp; j++)
										rule.x[(size_t)i * (size_t)ncomp + (size_t)j] = gp(i, j);
								rule.num_gp = ngp;
								rule.custom_rule_dimension = ncomp;
								done = true;
							}
						}
						delete eo_response;
					}
					delete[] argv;
					if (done) {
						// optional matching quadrature weights (one per gp)
						std::string wrequest = "integrationWeights";
						int wargc = 1;
						const char **wargv = new const char*[wargc];
						wargv[0] = wrequest.c_str();
						OutputDescriptor w_descriptor;
						OutputDescriptorStream w_stream(&w_descriptor);
						Response *w_response = elem->setResponse(wargv, wargc, w_stream);
						w_stream.finalizeSetResponse();
						if (w_response) {
							w_response->getResponse();
							const Vector &wv = w_response->getInformation().getData();
							if (wv.Size() == rule.num_gp) {
								rule.w.resize((size_t)wv.Size());
								for (int i = 0; i < wv.Size(); i++)
									rule.w[(size_t)i] = wv(i);
							}
							delete w_response;
						}
						delete[] wargv;
						return;
					}
					/* else: fall through to the legacy 1D strategies below */
				}
				/*
				ask for integrationPoints ...
				*/
				{
					//std::cout << "get custom gp: trying with \"integrationPoints\"...\n";
					bool done = false;
					std::string request = "integrationPoints";
					int argc = 1;
					const char **argv = new const char*[argc];
					argv[0] = request.c_str();
					OutputDescriptor eo_descriptor;
					OutputDescriptorStream eo_stream(&eo_descriptor);
					Response *eo_response = elem->setResponse(argv, argc, eo_stream);
					eo_stream.finalizeSetResponse();
					if (eo_response) {
						eo_response->getResponse();
						const Vector &data = eo_response->getInformation().getData();
						if (data.Size() > 0) {
							if (data.Size() == 1) {
								rule.x.resize(1);
								rule.x[0] = 0.0;
							}
							else {
								double x_min = std::numeric_limits<double>::max();
								double x_max = -x_min;
								for (int i = 0; i < data.Size(); i++) {
									double ieta = data[i];
									if (ieta < x_min)
										x_min = ieta;
									else if (ieta > x_max)
										x_max = ieta;
								}
								double span = x_max - x_min;
								if (span == 0.0) {
									rule.x.resize((size_t)data.Size(), 0.0);
								}
								else {
									rule.x.resize((size_t)data.Size());
									for (int i = 0; i < data.Size(); i++)
										rule.x[(size_t)i] = 2.0*(data[i] - x_min) / span - 1.0;
								}
							}
							done = true;
						}
						delete eo_response;
					}
					delete[] argv;
					if (done)
						return;
				}
				/*
				..., otherwise, ask for a dummy response on all sections ...
				*/
				{
					//std::cout << "get custom gp: trying with \"section(all)\"...\n";
					bool done = false;
					std::string request1 = "section";
					std::string request2 = "dummy";
					int argc = 2;
					const char **argv = new const char*[argc];
					argv[0] = request1.c_str();
					argv[1] = request2.c_str();
					OutputDescriptor eo_descriptor;
					OutputDescriptorStream eo_stream(&eo_descriptor);
					Response *eo_response = elem->setResponse(argv, argc, eo_stream);
					eo_stream.finalizeSetResponse();
					if (eo_response)
						delete eo_response; // not needed now
					eo_descriptor.getGaussLocations(rule.x);
					if (rule.x.size() > 0)
						done = true;
					delete[] argv;
					if (done)
						return;
				}
				/*
				..., otherwise, ask for a dummy response on all sections, finding out what is the number of gauss points
				*/
				{
					//std::cout << "get custom gp: trying with \"section(1,2,..,N)\"...\n";
					bool done = false;
					std::string request1 = "section";
					if (elem->getClassTag() == ELE_TAG_MVLEM ||
						elem->getClassTag() == ELE_TAG_SFI_MVLEM ||
						elem->getClassTag() == ELE_TAG_MVLEM_3D || 
						elem->getClassTag() == ELE_TAG_SFI_MVLEM_3D ||
						elem->getClassTag() == ELE_TAG_E_SFI_MVLEM_3D) {
						request1 = "material";
					}
					std::string request3 = "dummy";
					int argc = 3;
					const char **argv = new const char*[argc];
					argv[0] = request1.c_str();
					argv[2] = request3.c_str();

					int trial_num = 0;
					double rule_weight_sum = 0.0;
					while (true) {
						trial_num++;
						if (trial_num > LADRUNO_MAX_TRIAL_NSEC) {
							// not expected; the loop still needs a limit
							//opserr << "LadrunoRecorder warning: iterative guess of ngp: reached maximum number of iteration, giving up...\n";
							break;
						}
						std::stringstream ss_trial_num; ss_trial_num << trial_num;
						std::string s_trial_num = ss_trial_num.str();
						argv[1] = s_trial_num.c_str();
						OutputDescriptor eo_descriptor;
						OutputDescriptorStream eo_stream(&eo_descriptor);
						Response *eo_response = elem->setResponse(argv, argc, eo_stream);
						eo_stream.finalizeSetResponse();
						if (eo_response)
							delete eo_response; // not needed now
						std::vector<double> trial_x;
						std::vector<double> trial_w;
						eo_descriptor.appendGaussLocation(trial_x);
						eo_descriptor.appendGaussWeight(trial_w);
						if (trial_x.size() > 0) {
							if (trial_x.size() > 1) {
								// not expected
								opserr << "LadrunoRecorder warning: iterative guess of ngp: expected 1 trial gauss location, given = "
									<< (int)trial_x.size()
									<< "\nonly the first one will be considered\n";
							}
							rule.x.push_back(trial_x[0]);
							rule_weight_sum += trial_w[0];
						}
						else {
							// the maximum number of gauss points for this element is reached
							break;
						}
					}
					if (rule.x.size() > 0) {
						if (std::abs(rule_weight_sum - 2.0) > 1.0e-8) {
							// don't do auto-normalization if the integration weight is explicitly given
							// (only considered valid if the integration span is 2.0)
							if (rule.x.size() == 1) {
								rule.x[0] = 0.0;
							}
							else {
								double x_min = std::numeric_limits<double>::max();
								double x_max = -x_min;
								for (size_t i = 0; i < rule.x.size(); i++) {
									double ieta = rule.x[i];
									if (ieta < x_min)
										x_min = ieta;
									else if (ieta > x_max)
										x_max = ieta;
								}
								double span = x_max - x_min;
								if (span == 0.0) {
									for (size_t i = 0; i < rule.x.size(); i++)
										rule.x[i] = 0.0;
								}
								else {
									for (size_t i = 0; i < rule.x.size(); i++)
										rule.x[i] = 2.0 * (rule.x[i] - x_min) / span - 1.0;
								}
							}
						}
						done = true;
					}

					delete[] argv;
					if (done)
						return;
				}
				/*
				..., otherwise none of the strategies above worked: the element must implement
				"integrationPoints" or a per-section "integrationPoint" response
				*/
				opserr << "LadrunoRecorder warning: cannot get custom integration rule from element "
					<< elem->getTag() << "[class = " << elem->getClassType() << "]\n";
			}

			std::map<int, ElementIntegrationRule> registered_custom_rules;
			std::map<int, ElementWithSameClassTagCollection> items;
		};
	}

}

	/* ======================================================================
	ElementResultSource: wraps an already-built (header, response-collection)
	bucket as a ladruno::ResultSource. The bucket building and the per-step
	driving live in LadrunoRecorder (initElementSources /
	recordResultsOnElements).
	====================================================================== */
	class ElementResultSource : public ResultSource
	{
	public:
		ElementResultSource(const std::vector<std::string>& request,
			const detail::element::OutputDescriptorHeader& header,
			detail::element::OutputResponseCollection& bucket);

		const ResultSchema& schema() const override { return m_schema; }
		const std::vector<int>& ids() const override { return m_ids; }
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		bool requiresPartitionReduction() const override { return false; }

	private:
		void buildSchema(const std::vector<std::string>& request);

	private:
		const detail::element::OutputDescriptorHeader& m_header;
		detail::element::OutputResponseCollection& m_bucket;
		ResultSchema m_schema;
		std::vector<int> m_ids;
		// latch for the response-size mismatch diagnostic in evaluate(). The
		// mismatch is a property of the (element, bucket) pairing, so it repeats
		// identically on every element on every step; unlatched it floods the
		// log (and dominates the runtime) for the rest of the analysis.
		bool m_size_warned;
	};

	/* ======================================================================
	BASIS/QUADRATURE capture. The pieces:
	  - detail::element::BasisInfo                  (the captured descriptor)
	  - OutputDescriptorStream <ElementBasis> tag (captureBasisInfo)
	  - ElementCollection::getElementBasisInfo()  (probes "basisInfo")
	  - getCustomGaussPointLocations() multi-dim  (captureIntegrationPoints +
	                                                captureIntegrationWeights)
	  - getGeometryAndIntRuleByClassTag()         (deriveLegacyBasis fallback)
	LadrunoRecorder composes these: the probe overrides the legacy-derived
	FAMILY/ORDER/TOPOLOGY/PARAM_DOMAIN/RATIONAL/NUM_CTRL, and GP_PARAM /
	GP_WEIGHT come from the captured rule. Control-point weights of rational
	bases are not captured.
	====================================================================== */

} // namespace ladruno
#endif // Ladruno_ElementResults_h
