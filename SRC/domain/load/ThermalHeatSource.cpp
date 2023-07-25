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

#include <ThermalHeatSource.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector ThermalHeatSource::data(1);

ThermalHeatSource::ThermalHeatSource(int tag, int theElementTag, double q_)
    : ElementalLoad(tag, LOAD_TAG_ThermalHeatSource, theElementTag),
      q(q_)
{
    opserr << "Creating ThermalHeatSource for element " << theElementTag << " and rate q = " << q << endln;
}

ThermalHeatSource::ThermalHeatSource()
    : ElementalLoad(LOAD_TAG_ThermalHeatSource),
      q(0)
{
}

ThermalHeatSource::~ThermalHeatSource()
{
}

const Vector &
ThermalHeatSource::getData(int &type, double loadFactor)
{
    type = LOAD_TAG_ThermalHeatSource;

    data(0) = loadFactor * q;

    return data;
}

int
ThermalHeatSource::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(4);
    int dataTag = this->getDbTag();
    data(0) = this->getTag();
    data(1) = dataTag;
    data(2) = eleTag;
    data(3) = q;

    res = theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING ThermalHeatSource::sendSelf() - " << this->getTag() << " failed to send iddata\n";
        return res;
    }

    return -1;
}

int
ThermalHeatSource::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(4);
    int dataTag = this->getDbTag();

    res = theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING ThermalHeatSource::recvSelf() - " << this->getTag() << " failed to receive data\n";
        return res;
    }
    this->setTag(data(0));
    eleTag = data(2);

    dir = data(3);

    return res;
}

void
ThermalHeatSource::Print(OPS_Stream &s, int flag)
{
    s << "ThermalHeatSource...";
    s << "  element acted on: " << eleTag << " dir = " << dir << endln;
}

