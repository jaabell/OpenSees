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
// Description: shared vocabulary of the ladruno HDF5 recorder (recorder
// ladruno): command-line option kinds, string and local-axis utilities,
// fiber-section and section-assignment metadata, the result enums, the output
// frequency and the per-process state. Every symbol lives in namespace
// ladruno, so this module links into the same binary as MPCORecorder, whose
// equivalent utilities carry the same unqualified names.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#ifndef Ladruno_Types_h
#define Ladruno_Types_h


// hdf5
#include "hdf5.h"

// std
#include <string>
#include <vector>
#include <list>
#include <map>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdint.h>
// needed by ported code below (Timer uses std::cout and clock_t)
#include <iostream>
#include <ctime>

// opensees
#include "OPS_Globals.h" // for opserr
#include "Vector.h"      // used by utils::locax, utils::misc, detail::ProcessInfo
#include "Domain.h"      // used by detail::ProcessInfo
#include "classTags.h"   // ELE_TAG_* used by utils::shell::isShellElementTag

namespace ladruno {

	// replaces the original `#define HID_INVALID -1`
	inline constexpr hid_t HID_INVALID = -1;

/*namespace for utilities*/
namespace utils {

	/*
	utilities for parsing
	*/
	namespace parsing {

		enum option_type {
			opt_none,
			opt_result_on_nodes,
			opt_result_on_nodes_sens,
			opt_result_on_elements,
			opt_time,
			opt_region,
			opt_global       // -G energy <regionTag...> : energy balance
		};

	}

	/*
	utilities for string
	*/
	namespace strings {

		template<class T>
		inline std::string to_string(T x) {
			std::stringstream ss;
			ss << x;
			return ss.str();
		}

		inline void split(const std::string &text, char sep, std::vector<std::string> &tokens, bool skip_empty = false)
		{
			if (tokens.size() > 0) tokens.clear();
			std::size_t start = 0, end = 0;
			while (true)
			{
				end = text.find(sep, start);
				if (end == std::string::npos)
				{
					if (start < text.size())
						tokens.push_back(text.substr(start));
					else if (!skip_empty)
						tokens.push_back("");
					break;
				}
				std::string subs = text.substr(start, end - start);
				if (skip_empty)
				{
					if (subs.size() > 0)
						tokens.push_back(subs);
				}
				else
				{
					tokens.push_back(subs);
				}
				start = end + 1;
			}
		}

		struct concatenator {
			template<class T>
			inline concatenator &operator << (const T &x) {
				ss << x;
				return *this;
			}
			std::stringstream ss;
		};

	}

	/*
	utilities for shells
	*/
	namespace shell {

		inline bool isShellElementTag(int ele_tag) {
			/**
			used to check if an element is a shell.
			used in the workaround for shell sections.
			note that the keyword "section" doesn't work with shells (only in beams)
			shells use "material" keyword, and then the "fiber" keyword
			*/
			return (
				ele_tag == ELE_TAG_ShellMITC4 ||
				ele_tag == ELE_TAG_ShellMITC4Thermal ||
				ele_tag == ELE_TAG_ShellMITC9 ||
				ele_tag == ELE_TAG_ShellDKGQ ||
				ele_tag == ELE_TAG_ShellNLDKGQ ||
				ele_tag == ELE_TAG_ShellNLDKGQThermal ||
				ele_tag == ELE_TAG_ShellDKGT ||
				ele_tag == ELE_TAG_ShellNLDKGT ||
				ele_tag == ELE_TAG_ASDShellQ4 ||
				ele_tag == ELE_TAG_ASDShellT3
				);
		}

	}

	/*
	utilities for local axes
	*/
	namespace locax {

		struct quaternion {
			double w;
			double x;
			double y;
			double z;
			inline const double *data()const { return &w; }
		};

		inline quaternion quatFromMat(const Vector &vx, const Vector &vy, const Vector &vz) {
			// vx = first column of rot matrix
			// vy = second column
			// vz = third column
			double xx = vx(0);
			double yy = vy(1);
			double zz = vz(2);
			double tr = xx + yy + zz;
			quaternion Q;
			if ((tr > xx) && (tr > yy) && (tr > zz))
			{
				double S = std::sqrt(tr + 1.0) * 2.0;
				Q.w = 0.25 * S;
				Q.x = (vz(1) - vy(2)) / S;
				Q.y = (vx(2) - vz(0)) / S;
				Q.z = (vy(0) - vx(1)) / S;
			}
			else if ((xx > yy) && (xx > zz))
			{
				double S = std::sqrt(1.0 + xx - yy - zz) * 2.0;
				Q.w = (vz(1) - vy(2)) / S;
				Q.x = 0.25 * S;
				Q.y = (vx(1) + vy(0)) / S;
				Q.z = (vx(2) + vz(0)) / S;
			}
			else if (yy > zz)
			{
				double S = std::sqrt(1.0 + yy - xx - zz) * 2.0;
				Q.w = (vx(2) - vz(0)) / S;
				Q.x = (vx(1) + vy(0)) / S;
				Q.y = 0.25 * S;
				Q.z = (vy(2) + vz(1)) / S;
			}
			else
			{
				double S = std::sqrt(1.0 + zz - xx - yy) * 2.0;
				Q.w = (vy(0) - vx(1)) / S;
				Q.x = (vx(2) + vz(0)) / S;
				Q.y = (vy(2) + vz(1)) / S;
				Q.z = 0.25 * S;
			}
			// safe normalization
			double squared_norm = Q.x*Q.x + Q.y*Q.y + Q.z*Q.z + Q.w*Q.w;
			if (squared_norm > 0.0 && squared_norm != 1.0) {
				double n = std::sqrt(squared_norm);
				Q.x /= n;
				Q.y /= n;
				Q.z /= n;
				Q.w /= n;
			}
			return Q;
		}
	}

	/*
	misc utilities
	*/
	namespace misc {

		template<class T>
		inline bool areVectorsEqual(const std::vector<T> &a, const std::vector<T> &b) {
			if (a.size() != b.size())
				return false;
			for (size_t i = 0; i < a.size(); i++)
				if (a[i] != b[i])
					return false;
			return true;
		}

		template<>
		inline bool areVectorsEqual<std::string>(const std::vector<std::string> &a, const std::vector<std::string> &b) {
			if (a.size() != b.size())
				return false;
			for (size_t i = 0; i < a.size(); i++) {
				const std::string &sa = a[i];
				const std::string &sb = b[i];
				if (sa.size() != sb.size())
					return false;
				if (sa != sb)
					return false;
			}
			return true;
		}

		inline bool lessThanWithTol(double a, double b, double tol) {
			return (a - b) < -tol;
		}

		inline bool greaterThanWithTol(double a, double b, double tol) {
			return (a - b) > tol;
		}

		inline void bufferNodeResponseVec3u(size_t node_counter, int ndim, const Vector &iresponse, std::vector<double> &buffer)
		{
			size_t j = node_counter * ndim;
			if (iresponse.Size() > 0)
				buffer[j] = iresponse[0];
			else
				buffer[j] = 0.0;
			if (ndim > 1) {
				if (iresponse.Size() > 1)
					buffer[j + 1] = iresponse[1];
				else
					buffer[j + 1] = 0.0;
				if (ndim > 2) {
					if (iresponse.Size() > 2)
						buffer[j + 2] = iresponse[2];
					else
						buffer[j + 2] = 0.0;
				}
			}
		}

		inline void bufferNodeResponseVec3r(size_t node_counter, int ndim, const Vector &iresponse, std::vector<double> &buffer)
		{
			if (ndim == 2) {
				size_t j = node_counter; // node_counter * 1 (1 rotation component in 2D)
				if (iresponse.Size() > 2)
					buffer[j] = iresponse[2];
				else
					buffer[j] = 0.0;
			}
			else if (ndim == 3) {
				size_t j = node_counter * 3; // node_counter * 3 (3 rotation components in 3D)
				if (iresponse.Size() > 5) {
					buffer[j] = iresponse[3];
					buffer[j + 1] = iresponse[4];
					buffer[j + 2] = iresponse[5];
				}
				else {
					buffer[j] = 0.0;
					buffer[j + 1] = 0.0;
					buffer[j + 2] = 0.0;
				}

			}
		}
	}

}

#define LADRUNO_MAKE_STRING(X) (ladruno::utils::strings::concatenator() << X).ss.str()

/*utilities for element results. need to stay here before the h5 namespace*/
namespace detail {

	namespace element {

		struct FiberData
		{
			FiberData()
				: y(0.0), z(0.0), a(0.0) {}
			FiberData(double _y, double _z, double _a)
				: y(_y), z(_z), a(_a) {}
			inline bool operator < (const FiberData &other) const {
				const double rel_tol = 1.0e-5;
				double tol = std::max(std::abs(y), std::abs(other.y))*rel_tol;
				if (utils::misc::lessThanWithTol(y, other.y, tol)) return true;
				if (utils::misc::greaterThanWithTol(y, other.y, tol)) return false;
				tol = std::max(std::abs(z), std::abs(other.z))*rel_tol;
				if (utils::misc::lessThanWithTol(z, other.z, tol)) return true;
				if (utils::misc::greaterThanWithTol(z, other.z, tol)) return false;
				tol = std::max(std::abs(a), std::abs(other.a))*rel_tol;
				if (utils::misc::lessThanWithTol(a, other.a, tol)) return true;
				if (utils::misc::greaterThanWithTol(a, other.a, tol)) return false;
				return false; // equal
			}
			inline bool operator > (const FiberData &other) const {
				const double rel_tol = 1.0e-5;
				double tol = std::max(std::abs(y), std::abs(other.y))*rel_tol;
				if (utils::misc::greaterThanWithTol(y, other.y, tol)) return true;
				if (utils::misc::lessThanWithTol(y, other.y, tol)) return false;
				tol = std::max(std::abs(z), std::abs(other.z))*rel_tol;
				if (utils::misc::greaterThanWithTol(z, other.z, tol)) return true;
				if (utils::misc::lessThanWithTol(z, other.z, tol)) return false;
				tol = std::max(std::abs(a), std::abs(other.a))*rel_tol;
				if (utils::misc::greaterThanWithTol(a, other.a, tol)) return true;
				if (utils::misc::lessThanWithTol(a, other.a, tol)) return false;
				return false; // equal
			}
			inline const double *data()const { return &y; }
			double y;
			double z;
			double a;
		};

		typedef std::vector<FiberData> FiberDataCollection;

		typedef std::vector<int> FiberDataMaterialCollection;

		struct FiberSectionData
		{
			FiberSectionData() : fibers() {}
			inline bool operator < (const FiberSectionData &other) const {
				if (fibers.size() < other.fibers.size()) return true;
				if (fibers.size() > other.fibers.size()) return false;
				for (size_t i = 0; i < fibers.size(); i++) {
					const detail::element::FiberData &a = fibers[i];
					const detail::element::FiberData &b = other.fibers[i];
					if (a < b) return true;
					if (a > b) return false;
				}
				if (materials.size() < other.materials.size()) return true;
				if (materials.size() > other.materials.size()) return false;
				for (size_t i = 0; i < materials.size(); i++) {
					int a = materials[i];
					int b = other.materials[i];
					if (a < b) return true;
					if (a > b) return false;
				}
				return false; // equal
			}
			inline bool operator > (const FiberSectionData &other) const {
				if (fibers.size() > other.fibers.size()) return true;
				if (fibers.size() < other.fibers.size()) return false;
				for (size_t i = 0; i < fibers.size(); i++) {
					const detail::element::FiberData &a = fibers[i];
					const detail::element::FiberData &b = other.fibers[i];
					if (a > b) return true;
					if (a < b) return false;
				}
				if (materials.size() > other.materials.size()) return true;
				if (materials.size() < other.materials.size()) return false;
				for (size_t i = 0; i < materials.size(); i++) {
					int a = materials[i];
					int b = other.materials[i];
					if (a > b) return true;
					if (a < b) return false;
				}
				return false; // equal
			}
			FiberDataCollection fibers;
			FiberDataMaterialCollection materials;
		};

		struct ElemGaussPair
		{
			ElemGaussPair()
				: elem_id(0), gauss_id(0) {}
			ElemGaussPair(int _elem_id, int _gauss_id)
				: elem_id(_elem_id), gauss_id(_gauss_id) {}
			inline const int *data() const { return &elem_id; }
			int elem_id;
			int gauss_id;
		};

		struct SectionAssignment
		{
			SectionAssignment()
				: is_new(true), name("UnknownClassType"), fiber_section_data(), assignments() {}
			bool is_new;
			std::string name;
			FiberSectionData fiber_section_data;
			std::vector<ElemGaussPair> assignments;
		};

	}

}

/*ladruno enums*/
namespace detail {

	struct ResultType {
		enum Enum {
			Generic = 0,
			Modal
		};
	};

	struct ResultDataType {
		enum Enum {
			Scalar = 0,
			Vectorial,
			Tensorial,
			BeamForceDeformation,
			ShellForceDeformation
		};
	};

	struct NodalResultType {
		enum Enum {
			Displacement,
			Rotation,
			Velocity,
			AngularVelocity,
			Acceleration,
			AngularAcceleration,
			ReactionForce,
			ReactionMoment,
			ReactionForceIncludingInertia,
			ReactionMomentIncludingInertia,
			RayleighForce,
			RayleighMoment,
			UnbalancedForce,
			UnbalancedForceIncludingInertia,
			UnbalancedMoment,
			UnbalancedMomentIncludingInertia,
			Pressure,
			ModesOfVibration,
			ModesOfVibrationRotational,
			// sensitivity
			DisplacementSensitivity,
			RotationSensitivity,
			VelocitySensitivity,
			AngularVelocitySensitivity,
			AccelerationSensitivity,
			AngularAccelerationSensitivity
		};
	};

	struct ElementGeometryType {
		enum Enum {
			Custom = 0,
			Line_2N = 1,
			Line_3N,
			Triangle_3N = 100,
			Triangle_6N,
			Quadrilateral_4N = 200,
			Quadrilateral_8N,
			Quadrilateral_9N,
			Tetrahedron_4N = 300,
			Tetrahedron_10N,
			Hexahedron_8N = 400,
			Hexahedron_20N,
			Hexahedron_27N,
			Quadrilateral_CohesiveBand_4N = 500
		};
	};

	struct ElementIntegrationRuleType {
		enum Enum {
			NoIntegrationRule = 0,
			Line_GaussLegendre_1 = 1,
			Line_GaussLegendre_2,
			Line_GaussLegendre_3,
			Triangle_GaussLegendre_1 = 100,
			Triangle_GaussLegendre_2,
			Triangle_GaussLegendre_2B,
			Triangle_GaussLegendre_2C,
			Quadrilateral_GaussLegendre_1 = 200,
			Quadrilateral_GaussLegendre_2,
			Quadrilateral_GaussLegendre_3,
			Tetrahedron_GaussLegendre_1 = 300,
			Tetrahedron_GaussLegendre_2,
			Hexahedron_GaussLegendre_1 = 400,
			Hexahedron_GaussLegendre_2,
			Hexahedron_GaussLegendre_3,
			CustomIntegrationRule = 1000
		};
	};

	struct ElementOutputDescriptorType {
		/**
		these paths are interpreted as
		results on elements (either constant over element (# metadata = 1) or per element node (# metadata = # element nodes)
		[Element]
		these paths are interpreted as
		results on gauss points (# metadata = # element gauss points)
		[Element.Gauss]
		[Element.Gauss.Material]
		[Element.Gauss.Section]
		[Element.Gauss.Section.Material]
		[Element.Material]
		[Element.Section]
		[Element.Section.Material]
		these paths are interpreted as
		results on fibers (sub-integration points) (# metadata = # element gauss points AND (for each metadata) multiplicity = # fibers)
		[Element.Gauss.Section.Fiber]
		[Element.Gauss.Section.Fiber.Material]
		[Element.Section.Fiber]
		[Element.Section.Fiber.Material]
		*/
		enum Enum {
			Element = 0,
			Gauss,
			Section,
			Fiber,
			Material
		};

		static const char* toString(detail::ElementOutputDescriptorType::Enum _type) {
			switch (_type) {
			case detail::ElementOutputDescriptorType::Element: return "ELEMENT";
			case detail::ElementOutputDescriptorType::Gauss: return "GAUSS";
			case detail::ElementOutputDescriptorType::Section: return "SECTION";
			case detail::ElementOutputDescriptorType::Fiber: return "FIBER";
			case detail::ElementOutputDescriptorType::Material: return "MATERIAL";
			default:
				return "UnknownOutput";
			}
		}
	};
}

/*ladruno utilities*/
namespace detail {

	/*
	holds information for output frequency
	*/
	struct OutputFrequency {
		enum IncrementType {
			DeltaTime,
			NumberOfSteps
		};

		IncrementType type;
		double dt;
		int nsteps;
		double last_time;
		int last_step;

		OutputFrequency() : type(DeltaTime), dt(0.0), nsteps(1), last_time(0.0), last_step(0) {}
		void reset() {
			type = DeltaTime;
			dt = 0.0;
			nsteps = 1;
		}
	};

	/*
	a simple timer
	*/
	class Timer
	{
	public:
		Timer(const std::string &task_name)
			: m_task_name(task_name), m_t0(0), m_t1(0) {}
		inline void start() {
			m_t0 = clock();
			m_t1 = m_t0;
		}
		inline void stop() {
			m_t1 = clock();
			std::cout << "          LadrunoRecorder. Task = \"" << m_task_name << "\". Elapsed time: " << double(m_t1 - m_t0) / double(CLOCKS_PER_SEC) << "\n";
		}
	private:
		std::string m_task_name;
		clock_t m_t0;
		clock_t m_t1;
	};

	/*
	holds current information
	*/
	class ProcessInfo
	{
	public:
		ProcessInfo()
			// domain and model information
			: domain(0)
			, num_dimensions(0)
			, current_model_stage_id(-1)
			// some handles in the hdf5 file
			, h_file_id(HID_INVALID)
			, h_file_proplist(HID_INVALID)
			, h_group_proplist(HID_INVALID)
			// time step info
			, current_time_step_id(0)
			, current_time_step(0.0)
			// misc
			, eigen_first_initialization_done(false)
			, record_eigen_on_this_step(false)
			, eigen_last_time_set(0.0)
			, eigen_last_values()
			, store_data_f32(false)
			, deflate_level(1)
		{}
	public:
		// domain and model information
		Domain *domain;
		int num_dimensions;
		int current_model_stage_id;
		// some handles in the hdf5 file
		hid_t h_file_id;
		hid_t h_file_proplist;
		hid_t h_group_proplist;
		// time step info
		int current_time_step_id;
		double current_time_step;
		// misc
		bool eigen_first_initialization_done;
		bool record_eigen_on_this_step;
		double eigen_last_time_set;
		Vector eigen_last_values;
		// Opt-in `-precision f32` lossy mode: StreamingSink writes per-step result
		// DATA as on-disk float32 (HDF5 narrows from the in-memory double slab),
		// ~half the payload at ~7 significant digits. Default false = lossless f64
		// parity path. Coords / TIME / STEP / envelope datasets stay f64/int.
		bool store_data_f32;
		// deflate level of the time-series DATA datasets, set by
		// `-compress <0..9>` (0 = no shuffle/deflate filter; default 1).
		// Compression runs once per chunk.
		int deflate_level;
	};

}

} // namespace ladruno

#endif // Ladruno_Types_h
