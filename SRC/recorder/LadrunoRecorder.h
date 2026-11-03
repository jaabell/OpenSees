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
// Description: recorder ladruno. Writes an HDF5 file (.ladruno) holding the
// model (nodes, elements, local axes, sections, sets) and the requested
// nodal, element, whole-model and per-region results, one MODEL_STAGE group
// per change of the node/element set, each result as a chunked
// [T x nIds x nComp] series or as an envelope. Under a partitioned run each
// process writes its own part file.
//
// References:
//   M. Petracca, MPCORecorder (SRC/recorder/MPCORecorder.cpp), the HDF5
//     recorder this module is derived from.
//   The HDF Group, HDF5 Library and File Format, https://www.hdfgroup.org/solutions/hdf5/
//

#ifndef LadrunoRecorder_h
#define LadrunoRecorder_h

#include <Recorder.h>

class Domain;
class Channel;
class FEM_ObjectBroker;

class LadrunoRecorder : public Recorder
{
	friend void* OPS_LadrunoRecorder();

public:
	LadrunoRecorder();
	~LadrunoRecorder();

	int record(int commitTag, double timeStamp);
	virtual int restart(void);
	virtual int domainChanged(void);
	virtual int setDomain(Domain& theDomain);
	virtual int sendSelf(int commitTag, Channel& theChannel);
	virtual int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

private:
	int initialize();
	int writeModel();
	int writeModelNodes();
	int writeModelElements();
	int writeModelSets();
	int writeModelLocalAxes();
	int writeSections();

	int initNodeSources();
	int initElementSources();
	int initDomainSources();
	int recordResultsOnNodes();
	int recordResultsOnElements();
	int recordResultsOnDomain();

	// helpers (operate on the channel at the given index in m_data, so the
	// header never has to name private_data's nested channel types)
	void recordModeChannel(int node_channel_index);
	void writeElementColumnMap(int elem_channel_index);

	int clearSources();
	// Write the deferred ENVELOPES datasets for every channel (EnvelopeSink::finalize;
	// a no-op for StreamingSink). Called on the -flush cadence in envelope mode (an
	// in-place rewrite), at a stage change and once more at teardown.
	void finalizeAllSinks();
	// H5Fflush + restart the -flush cadence clock.
	void flushFile();
	// hash of the node/element/pressure-constraint set (tags and
	// object addresses). A domain-change stamp move rebuilds the MODEL_STAGE only
	// when this changed.
	unsigned long long topologyFingerprint() const;

	class private_data;
	private_data* m_data;
};

#endif // LadrunoRecorder_h
