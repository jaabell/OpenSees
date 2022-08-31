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

// Written: Massimo Petracca 
// Created: 02/2020
// Revision: A

#include <PluginUniaxialMaterial.h>
#include <PluginFramework.h>
#include <elementAPI.h>
#include <Element.h>
#include <ID.h>
#include <Channel.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <actor/message/Message.h>
#include <string>
#include <exception>
#include <limits>
#include <algorithm>
#include <cmath>

namespace details_plugin_uniaxial_material {

	inline PluginMaterialDescriptor* getDescriptor(const std::string& plugin_library, const std::string& plugin_function)
	{
		PluginMaterialDescriptor* descriptor = 0;
		try {
			descriptor = PluginFramework::instance().getMaterialDescriptor(plugin_library, plugin_function);
		}
		catch (const std::exception & ex) {
			opserr << "PluginUniaxialMaterial internal error while getting the requested procedure:\n"
				<< ex.what() << "\n";
			return 0;
		}
		catch (...) {
			opserr << "PluginUniaxialMaterial unknown exception while getting the requested procedure:\n";
			return 0;
		}
		return descriptor;
	}

	inline double callGet(PluginMaterialData *data, PluginMaterialJobType ijob)
	{
		// set
		double res = 0.0;
		data->response_size = 1;
		data->response = &res;
		// call
		int rcode = 0;
		data->proc(data, ijob, &rcode);
		// unset
		data->response_size = 0;
		data->response = 0;
		// check
		if (rcode != 0) {
			opserr << "PluginUniaxialMaterial Error: Failed in job = " << PluginFramework::materialJobTypeToString(ijob) << "\n";
			exit(-1);
		}
		// done
		return res;
	}

	inline int callSimple(PluginMaterialData* data, PluginMaterialJobType ijob)
	{
		// call
		int rcode = 0;
		data->proc(data, ijob, &rcode);
		return rcode;
	}

	inline int callGetResponse(PluginMaterialData* data, const PluginResponseDescriptor& d, Information &info, PluginMaterialJobType ijob, int responseID, int gradIndex = 0)
	{
		// set
		data->response_id = responseID;
		data->sens_grad_index = gradIndex;
		double value_scalar = 0.0;
		Vector value_vector;
		if (d.components.size() == 1) {
			data->response_size = 1;
			data->response = &value_scalar;
		}
		else {
			value_vector.resize(static_cast<int>(d.components.size()));
			value_vector.Zero();
			data->response_size = d.components.size();
			data->response = &value_vector(0);
		}

		// call
		int rcode = 0;
		data->proc(data, ijob, &rcode);

		// unset
		data->response_id = 0;
		data->response_size = 0;
		data->sens_grad_index = 0;
		data->response = 0;

		// check
		if (rcode != 0) {
			return -1;
		}

		// done
		if (d.components.size() == 1)
			return info.setDouble(value_scalar);
		else
			return info.setVector(value_vector);
	}

}

namespace details = details_plugin_uniaxial_material;

void*
OPS_PluginUniaxialMaterial(void)
{
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 3) { // check for mandatory arguments (at least 3)
		opserr << "Want: uniaxialMaterial plugin $tag   $pluginLibrary $pluginFunction    <+ plugin dependent arguments>\n";
		return 0;
	}

	int idata[1];
	int num_data;

	// get tag
	num_data = 1;
	if (OPS_GetInt(&num_data, idata) != 0) {
		opserr << "WARNING invalid integer tag: uniaxialMaterial plugin \n";
		return 0;
	}
	int tag = idata[0];

	// get library and function names
	std::string plugin_library = OPS_GetString();
	std::string plugin_function = OPS_GetString();
	
	// get plugin material descriptor
	PluginMaterialDescriptor* descriptor = details::getDescriptor(plugin_library, plugin_function);
	if (descriptor == 0)
		return 0;

	// allocate the PluginMaterialData. now material data has the tag, the pointer to
	// the procedure, and all data properly inizialized to zero
	PluginMaterialData* data = PluginFramework::instance().makeMaterialData(descriptor->procedure, tag);
	if (data == 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to initialize material data\n";
		return 0;
	}
	data->mat_type = PF_MAT_TYPE_UNIAXIAL;

	// call to get initialization info. Here the function must set initialization info for proper allocation
	// (message)
	int rcode = 0;
	data->proc(data, PF_MAT_GET_INIT_INFO, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_INIT_INFO\n";
		delete data;
		return 0;
	}
	if (data->mat_type != PF_MAT_TYPE_UNIAXIAL) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_INIT_INFO\n"
			"The plugin function provided is not a uniaxial material type.\n";
		delete data;
		return 0;
	}

	// now we can parse the message and get information about this material plugin.
	// information such as the input arguments, the available results/parameters/variables.
	// note: do it only once for each plugin material type
	if (data->message == 0) {
		opserr << 
			"PluginUniaxialMaterial Error: message field is null.\n"
			"Please you need to fill the message field to parse input arguments\n";
		delete data;
		return 0;
	}
	if (descriptor->parseMessage(data->message) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to parse message.\n";
		delete data;
		return 0;
	}
	data->message = 0; // don't need it anymore, the plugin is free to do whatever it likes

	// get number of arguments for next allocation
	data->n_param = static_cast<int>(descriptor->arguments.size()); 

	// now we can allocate data for input arguments
	if (PluginFramework::instance().allocateData(data) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to allocate material data\n";
		delete data;
		return 0;
	}

	// we correctly parsed the message. now we can parse the tcl command
	if (PluginFramework::instance().parseTclCommand(descriptor, data) != 0) {
		opserr << "PluginUniaxialMaterial Error: cannot parse the TCL command.\n";
		PluginFramework::instance().releaseData(data);
		delete data;
		return 0;
	}

	// call for initialization. Here the function must initialize the sate of the plugin material
	rcode = 0;
	data->proc(data, PF_MAT_INITIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_INITIALIZE\n";
		PluginFramework::instance().releaseData(data);
		delete data;
		return 0;
	}

	// create the new material with the current data, it will complete 
	// the setup of pointers in data
	PluginUniaxialMaterial* material = new PluginUniaxialMaterial(descriptor, data);

	// done
	return material;
}

PluginUniaxialMaterial::PluginUniaxialMaterial()
	: UniaxialMaterial(0, MAT_TAG_PluginUniaxialMaterial)
	, m_descriptor(0)
	, m_data(0)
	, m_lch(1.0)
	, m_lch_calculated(false)
{
}

PluginUniaxialMaterial::PluginUniaxialMaterial(PluginMaterialDescriptor* descr, PluginMaterialData* d)
	: UniaxialMaterial(d->tag, MAT_TAG_PluginUniaxialMaterial)
	, m_descriptor(descr)
	, m_data(d)
	, m_lch(1.0)
	, m_lch_calculated(false)
{
}

PluginUniaxialMaterial::~PluginUniaxialMaterial()
{
	if (m_data) {
		int rcode = 0;
		m_data->proc(m_data, PF_MAT_FINALIZE, &rcode);
		if (rcode != 0) {
			opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_FINALIZE\n";
		}
		PluginFramework::instance().releaseData(m_data);
		delete m_data;
	}
}

void PluginUniaxialMaterial::Print(OPS_Stream& s, int flag)
{
	s << "PluginUniaxialMaterial - Tag: " << getTag() << "\n";
}

const char* PluginUniaxialMaterial::getClassType() const
{
	return "PluginUniaxialMaterial";
}

int PluginUniaxialMaterial::sendSelf(int commitTag, Channel& theChannel)
{
	// serialize everything in a string using the following format:
	//
	// plugin_library
	// plugin_function
	// tag
	// nparams
	// param1 param2 param3 ... paramN
	// custom_string_from_plugin_material
	//

	// the serializer
	PluginSerializer ser;

	// serialize all data we have here
	ser << m_descriptor->library->library_name
		<< m_descriptor->procedure_name
		<< getTag()
		<< m_data->n_param
		<< PluginSerializer::darray_wrapper(m_data->param, static_cast<std::size_t>(m_data->n_param));

	// serialize extra data from the plugin if any
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_SERIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_SERIALIZE\n";
		return -1;
	}
	if (m_data->message == 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_SERIALIZE\n";
		return -1;
	}
	
	// serialize message and reset the original pointer to 0
	ser << m_data->message;
	m_data->message = 0;

	// get message string and size
	std::string msg_string = ser.str();
	std::vector<char> msg_data(msg_string.size() + 1);
	std::copy(msg_string.begin(), msg_string.end(), msg_data.begin());
	msg_data.back() = '\0';
	int msg_data_size = static_cast<int>(msg_string.size());

	// send message size
	ID idata(1);
	idata(0) = msg_data_size;
	if (theChannel.sendID(0, commitTag, idata) < 0) {
		opserr << "PluginUniaxialMaterial::sendSelf() - failed to send message size\n";
		return -1;
	}

	// send message
	Message msg(msg_data.data(), msg_data_size);
	if (theChannel.sendMsg(0, commitTag, msg) < 0) {
		opserr << "PluginUniaxialMaterial::sendSelf() - failed to send message\n";
		return -1;
	}

	// done
	return 0;
}

int PluginUniaxialMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	// initial check
	if ((m_data != 0) || (m_descriptor != 0)) {
		opserr << "PluginUniaxialMaterial::recvSelf(): m_data and m_descriptor must be null when calling this method!\n";
		return -1;
	}

	// de-serialize everything from a string using the following format:
	//
	// plugin_library
	// plugin_function
	// tag
	// nparams
	// param1 param2 param3 ... paramN
	// custom_string_from_plugin_material
	//

	// recv message size
	ID idata(1);
	if (theChannel.recvID(0, commitTag, idata) < 0) {
		opserr << "PluginUniaxialMaterial::recvSelf() - failed to recv message size\n";
		return -1;
	}
	int msg_data_size = idata(0);

	// recv message
	std::vector<char> msg_data(static_cast<size_t>(msg_data_size) + 1);
	Message msg(msg_data.data(), msg_data_size);
	if (theChannel.recvMsg(0, commitTag, msg) < 0) {
		opserr << "PluginUniaxialMaterial::recvSelf() - failed to recv message\n";
		return -1;
	}
	msg_data.back() = '\0';

	// the serializer
	PluginSerializer ser;

	// initialize it with all the received message
	ser << msg_data.data();

	// get library and function names
	std::string plugin_library;
	std::string plugin_function;
	ser >> plugin_library >> plugin_function;

	// get plugin material descriptor
	m_descriptor = details::getDescriptor(plugin_library, plugin_function);
	if (m_descriptor == 0)
		return -1;

	// get tag
	int tag;
	ser >> tag;
	setTag(tag);

	// allocate the PluginMaterialData. now material data has the tag, the pointer to
	// the procedure, and all data properly inizialized to zero
	m_data = PluginFramework::instance().makeMaterialData(m_descriptor->procedure, tag);
	if (m_data == 0) {
		opserr << "PluginUniaxialMaterial::recvSelf() - Error: Failed to initialize material data\n";
		return -1;
	}
	m_data->mat_type = PF_MAT_TYPE_UNIAXIAL;

	// NOTE: we need to do the next 2 code blocks because this may be the first 
	// time this material is sent to another process where the plugin library
	// has not been loaded yet!

	// call to get initialization info. Here the function must set initialization info for proper allocation
	// (message)
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_INIT_INFO, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial::recvSelf() - Error: Failed in job = PF_MAT_GET_INIT_INFO\n";
		return -1;
	}

	// now we can parse the message and get information about this material plugin.
	// information such as the input arguments, the available results/parameters/variables.
	// note: do it only once for each plugin material type
	if (m_data->message == 0) {
		opserr <<
			"PluginUniaxialMaterial Error: message field is null.\n"
			"Please you need to fill the message field to parse input arguments\n";
		return -1;
	}
	if (m_descriptor->parseMessage(m_data->message) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to parse message.\n";
		return -1;
	}
	m_data->message = 0; // don't need it anymore, the plugin is free to do whatever it likes

	// get number of parameters
	ser >> m_data->n_param;

	// make sure they are the same as the number of arguments in message
	if (m_data->n_param != static_cast<int>(m_descriptor->arguments.size())) {
		opserr << "PluginUniaxialMaterial Error: n_param obtained from the de-serialization is different "
			"from the one in the info message\n";
		return -1;
	}

	// now we can allocate data for input arguments
	if (PluginFramework::instance().allocateData(m_data) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to allocate material data\n";
		return -1;
	}

	// now we do not call the parseTclCommand, instead we obtain the values of input
	// arguments from de-serialization
	{
		PluginSerializer::darray_wrapper dummy(m_data->param, static_cast<std::size_t>(m_data->n_param));
		ser >> dummy;
	}

	// call for initialization. Here the function must initialize the sate of the plugin material
	rcode = 0;
	m_data->proc(m_data, PF_MAT_INITIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_INITIALIZE\n";
		return -1;
	}

	// now we can get the remaining string from the stream and let the plugin material
	// deserialize its custom data (instead of a standard initialization)
	std::string plugin_message = ser.remaining();
	
	// set the message pointer
	m_data->message = plugin_message.data();

	// deserialize from string
	rcode = 0;
	m_data->proc(m_data, PF_MAT_DESERIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_DESERIALIZE\n";
		return -1;
	}
	
	// reset the original pointer to 0
	m_data->message = 0;

	// done
	return 0;
}

int PluginUniaxialMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	if (argc < 1)
		return 0;

	std::string iarg(argv[0]);
	for (std::size_t i = 0; i < m_descriptor->parameters.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->parameters[i];
		if (d.name == iarg) {
			return param.addObject(d.id, this);
		}
	}

	// call the base class implementation
	return UniaxialMaterial::setParameter(argv, argc, param);
}

int PluginUniaxialMaterial::updateParameter(int parameterID, Information& info)
{
	for (std::size_t i = 0; i < m_descriptor->responses.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->responses[i];
		if (d.id == parameterID) {
			m_data->response_id = parameterID;
			m_data->response_size = 1;
			m_data->response = &info.theDouble;
			int res = details::callSimple(m_data, PF_MAT_UPDATE_PARAMETER);
			m_data->response_id = 0;
			m_data->response_size = 0;
			m_data->response = 0;
			return res;
		}
	}

	// call the base class implementation
	return UniaxialMaterial::updateParameter(parameterID, info);
}

int PluginUniaxialMaterial::activateParameter(int parameterID)
{
	m_data->response_id = parameterID;
	int res = details::callSimple(m_data, PF_MAT_ACTIVATE_PARAMETER);
	m_data->response_id = 0;
	return res;
}

int PluginUniaxialMaterial::setVariable(const char* variable, Information& info)
{
	if (variable == 0)
		return -1;

	std::string iarg(variable);
	for (std::size_t i = 0; i < m_descriptor->variables.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->variables[i];
		if (d.name == iarg) {
			if (d.components.size() > 1) {
				if (info.theVector == 0) {
					opserr << "PluginUniaxialMaterial::setVariable - Error: theVector is null\n";
					exit(-1);
				}
				if (info.theVector->Size() != static_cast<int>(d.components.size())) {
					opserr << "PluginUniaxialMaterial::setVariable - Error: theVector.Size() is different from the number of components\n";
					exit(-1);
				}
				m_data->response = &((*info.theVector)(0));
			}
			else if (d.components.size() == 1) {
				m_data->response = &(info.theDouble);
			}
			else {
				return -1;
			}
			m_data->response_size = d.components.size();
			int res = details::callSimple(m_data, PF_MAT_SET_VARIABLE);
			m_data->response = 0;
			m_data->response_size = 0;
			return res;
		}
	}

	// call the base class implementation
	return UniaxialMaterial::setVariable(variable, info);
}

int PluginUniaxialMaterial::getVariable(const char* variable, Information& info)
{
	if (variable == 0)
		return -1;

	std::string iarg(variable);
	for (std::size_t i = 0; i < m_descriptor->variables.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->variables[i];
		if (d.name == iarg) {
			int res = -1;
			if (d.components.size() > 1) {
				Vector vector_data(static_cast<int>(d.components.size()));
				m_data->response = &(vector_data(0));
				m_data->response_size = d.components.size();
				res = details::callSimple(m_data, PF_MAT_GET_VARIABLE);
				m_data->response = 0;
				m_data->response_size = 0;
				info.setVector(vector_data);
			}
			else if (d.components.size() == 1) {
				double scalar_data(0.0);
				m_data->response = &scalar_data;
				m_data->response_size = d.components.size();
				res = details::callSimple(m_data, PF_MAT_GET_VARIABLE);
				m_data->response = 0;
				m_data->response_size = 0;
				info.setDouble(scalar_data);
			}
			return res;
		}
	}

	// call the base class implementation
	return UniaxialMaterial::getVariable(variable, info);
}

int PluginUniaxialMaterial::setTrialStrain(double strain, double strainRate)
{
	// compute lch only once
	if (!m_lch_calculated && (ops_TheActiveElement != 0)) {
		m_lch = ops_TheActiveElement->getCharacteristicLength();
		m_lch_calculated = true;
	}
	// set
	bool has_temp = m_data->temperature != 0;
	double dummy = 0.0;
	if (!has_temp)
		m_data->temperature = &dummy;
	m_data->strain = &strain;
	m_data->strain_rate = &strainRate;
	m_data->lch = m_lch;
	m_data->dT = ops_Dt;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_COMPUTE, &rcode);
	// unset
	if (!has_temp)
		m_data->temperature = 0;
	m_data->strain = 0;
	m_data->strain_rate = 0;
	m_data->lch = 0.0;
	m_data->dT = 0.0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_COMPUTE\n";
	}
	// done
	return rcode;
}

int PluginUniaxialMaterial::setTrialStrain(double strain, double temperature, double strainRate)
{
	// set
	m_data->temperature = &temperature;
	// call
	int res = setTrialStrain(strain, strainRate);
	// unset
	m_data->temperature = 0;
	// done
	return res;
}

double PluginUniaxialMaterial::getStrain()
{
	return details::callGet(m_data, PF_MAT_GET_STRAIN);
}

double PluginUniaxialMaterial::getStrainRate()
{
	return details::callGet(m_data, PF_MAT_GET_STRAIN_RATE);
}

double PluginUniaxialMaterial::getStress()
{
	return details::callGet(m_data, PF_MAT_GET_STRESS);
}

double PluginUniaxialMaterial::getTangent()
{
	return details::callGet(m_data, PF_MAT_GET_TANGENT);
}

double PluginUniaxialMaterial::getInitialTangent()
{
	return details::callGet(m_data, PF_MAT_GET_INITIAL_TANGENT);
}

double PluginUniaxialMaterial::getDampTangent()
{
	return details::callGet(m_data, PF_MAT_GET_DAMP_TANGENT);
}

double PluginUniaxialMaterial::getRho()
{
	return details::callGet(m_data, PF_MAT_GET_RHO);
}

int PluginUniaxialMaterial::commitState()
{
	return details::callSimple(m_data, PF_MAT_COMMIT);
}

int PluginUniaxialMaterial::revertToLastCommit()
{
	return details::callSimple(m_data, PF_MAT_REVERT);
}

int PluginUniaxialMaterial::revertToStart()
{
	return details::callSimple(m_data, PF_MAT_REVERT_TO_START);
}

UniaxialMaterial* PluginUniaxialMaterial::getCopy()
{
	// allocate the PluginMaterialData. now material data has the tag, the pointer to
	// the procedure, and all data properly inizialized to zero
	PluginMaterialData* data = PluginFramework::instance().makeMaterialData(m_descriptor->procedure, getTag());
	if (data == 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to initialize material data (in getCopy)\n";
		return 0;
	}

	// call for initialization info (PF_MAT_GET_INIT_INFO) avoided here. Can copy those values
	data->mat_type = m_data->mat_type;
	data->n_param = m_data->n_param;

	// now we can allocate data for input arguments and state variables
	if (PluginFramework::instance().allocateData(data) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to allocate material data (in getCopy)\n";
		delete data;
		return 0;
	}

	// copy data param and state. parseTclCommand avoided here. Can copy those values
	for (int i = 0; i < m_data->n_param; i++)
		data->param[i] = m_data->param[i];

	// call for initialization. Here the function must initialize the sate of the plugin material
	int rcode = 0;
	data->proc(data, PF_MAT_INITIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_INITIALIZE\n";
		delete data;
		return 0;
	}

	// create the new material with the current data, it will complete 
	// the setup of pointers in data
	PluginUniaxialMaterial* material = new PluginUniaxialMaterial(m_descriptor, data);

	// done
	return material;
}

Response* PluginUniaxialMaterial::setResponse(const char** argv, int argc, OPS_Stream& theOutputStream)
{
	if (argc < 1)
		return 0;

	std::string iarg(argv[0]);
	for (std::size_t i = 0; i < m_descriptor->responses.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->responses[i];
		if ((d.name == iarg) && (d.components.size() > 0)) {

			theOutputStream.tag("UniaxialMaterialOutput");
			theOutputStream.attr("matType", this->getClassType());
			theOutputStream.attr("matTag", this->getTag());

			for(std::size_t j = 0; j < d.components.size(); j++)
				theOutputStream.tag("ResponseType", d.components[j].c_str());

			if (d.components.size() == 1) {
				return new MaterialResponse(this, d.id, 0.0);
			}
			else {
				return new MaterialResponse(this, d.id, Vector(static_cast<int>(d.components.size())));
			}

			theOutputStream.endTag();
		}
	}
	
	// call the base class implementation
	return UniaxialMaterial::setResponse(argv, argc, theOutputStream);
}

int PluginUniaxialMaterial::getResponse(int responseID, Information& info)
{
	for (std::size_t i = 0; i < m_descriptor->responses.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->responses[i];
		if ((d.id == responseID) && (d.components.size() > 0)) {
			return details::callGetResponse(m_data, d, info, PF_MAT_GET_RESPONSE, responseID);
		}
	}

	// call the base class implementation
	return UniaxialMaterial::getResponse(responseID, info);
}

int PluginUniaxialMaterial::getResponseSensitivity(int responseID, int gradIndex, Information& info)
{
	for (std::size_t i = 0; i < m_descriptor->responses.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->responses[i];
		if ((d.id == responseID) && (d.components.size() > 0)) {
			return details::callGetResponse(m_data, d, info, PF_MAT_GET_RESPONSE_SENSITIVITY, responseID, gradIndex);
		}
	}
	
	// call the base class implementation
	return UniaxialMaterial::getResponseSensitivity(responseID, gradIndex, info);
}

bool PluginUniaxialMaterial::hasFailed()
{
	double res = details::callGet(m_data, PF_MAT_GET_IS_FAILED);
	return std::abs(res) > std::numeric_limits<double>::epsilon();
}

double PluginUniaxialMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
	m_data->sens_grad_index = gradIndex;
	m_data->sens_conditonal = conditional ? 1 : 0;
	double res = details::callGet(m_data, PF_MAT_GET_STRESS_SENSITIVITY);
	m_data->sens_grad_index = 0;
	m_data->sens_conditonal = 0;
	return res;
}

double PluginUniaxialMaterial::getStrainSensitivity(int gradIndex)
{
	m_data->sens_grad_index = gradIndex;
	double res = details::callGet(m_data, PF_MAT_GET_STRAIN_SENSITIVITY);
	m_data->sens_grad_index = 0;
	return res;
}

double PluginUniaxialMaterial::getTangentSensitivity(int gradIndex)
{
	m_data->sens_grad_index = gradIndex;
	double res = details::callGet(m_data, PF_MAT_GET_TANGENT_SENSITIVITY);
	m_data->sens_grad_index = 0;
	return res;
}

double PluginUniaxialMaterial::getInitialTangentSensitivity(int gradIndex)
{
	m_data->sens_grad_index = gradIndex;
	double res = details::callGet(m_data, PF_MAT_GET_INITIAL_TANGENT_SENSITIVITY);
	m_data->sens_grad_index = 0;
	return res;
}

double PluginUniaxialMaterial::getDampTangentSensitivity(int gradIndex)
{
	m_data->sens_grad_index = gradIndex;
	double res = details::callGet(m_data, PF_MAT_GET_DAMP_TANGENT_SENSITIVITY);
	m_data->sens_grad_index = 0;
	return res;
}

double PluginUniaxialMaterial::getRhoSensitivity(int gradIndex)
{
	m_data->sens_grad_index = gradIndex;
	double res = details::callGet(m_data, PF_MAT_GET_RHO_SENSITIVITY);
	m_data->sens_grad_index = 0;
	return res;
}

int PluginUniaxialMaterial::commitSensitivity(double strainGradient, int gradIndex, int numGrads)
{
	m_data->sens_strain_gradient = &strainGradient;
	m_data->sens_grad_index = gradIndex;
	m_data->sens_num_grad = numGrads;
	int res = details::callSimple(m_data, PF_MAT_COMMIT_SENSITIVITY);
	m_data->sens_strain_gradient = 0;
	m_data->sens_grad_index = 0;
	m_data->sens_num_grad = 0;
	return res;
}

double PluginUniaxialMaterial::getEnergy()
{
	return details::callGet(m_data, PF_MAT_GET_ENERGY);
}
