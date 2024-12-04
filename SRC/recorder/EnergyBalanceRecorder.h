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


#ifndef EnergyBalanceRecorder_h
#define EnergyBalanceRecorder_h

// Written: jaabell
// Created: 12/2024
// Revision: A
//
// Description: This file contains the class definition for
// EnergyBalanceRecorder. A EnergyBalanceRecorder is used store
// the energy balance in the domain at each time step



#include <Recorder.h>
#include <ID.h>
#include <Vector.h>
#include <TimeSeries.h>

class Domain;
class FE_Datastore;
class Node;

#define EBR_NUM_ENERGY_COMPONENTS 4

class EnergyBalanceRecorder: public Recorder
{
public:
    EnergyBalanceRecorder();
    EnergyBalanceRecorder(
        Domain &theDomain,
        OPS_Stream &theOutputHandler,
        bool echoTimeFlag = true);

    ~EnergyBalanceRecorder();

    int record(int commitTag, double timeStamp);
    int flush();

    int domainChanged(void);
    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);
  

protected:

private:
    int initialize(void);

    Domain *theDomain;
    OPS_Stream *theOutputHandler;

    Vector response;

    bool echoTimeFlag;   // flag indicating whether time to be included in o/p
    bool initializationDone;

    double time_last;
    double internal_energy;
    double damping_work;
    double unbalanced_load_work;
};

#endif
