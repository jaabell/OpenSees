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

// Written: José A. Abell.
//
// Description: Elemental load that specifies a heat source.

#ifndef ThermalHeatSource_h
#define ThermalHeatSource_h

#include <ElementalLoad.h>

class ThermalHeatSource : public ElementalLoad
{
public:
	ThermalHeatSource(int tag, int eleTag, double q_);
	ThermalHeatSource();
	~ThermalHeatSource();

	const Vector &getData(int &type, double loadFactor);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s, int flag = 0);

protected:

private:
	int dir;
	static Vector data;
	double q;
};

#endif
