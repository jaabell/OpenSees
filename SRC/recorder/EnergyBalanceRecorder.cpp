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

// Written: fmk
//
// Description: This file contains the class definition for EnergyBalanceRecorder.
// A EnergyBalanceRecorder is used to record the specified dof responses
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) EnergyBalanceRecorder.C, revA"

#include <EnergyBalanceRecorder.h>
#include <Domain.h>
#include <Parameter.h>
#include <Node.h>
#include <NodeIter.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <FE_Datastore.h>
#include <FEM_ObjectBroker.h>
#include <Pressure_Constraint.h>
#include <MeshRegion.h>
#include <TimeSeries.h>

#include <Element.h>
#include <ElementIter.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <DatabaseStream.h>
#include <TCP_Stream.h>

#include <elementAPI.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>

void*
OPS_EnergyBalanceRecorder()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING: recorder EnergyBalance ";
        opserr << "-file <fileName> <-time>\n";
        opserr << "OPS_GetNumRemainingInputArgs() = " << OPS_GetNumRemainingInputArgs() << endln;
        return 0;
    }

    const char* responseID = 0;
    OPS_Stream *theOutputStream = 0;
    const char* filename = 0;

    const int STANDARD_STREAM = 0;
    const int DATA_STREAM = 1;
    const int XML_STREAM = 2;
    const int DATABASE_STREAM = 3;
    const int BINARY_STREAM = 4;
    const int DATA_STREAM_CSV = 5;
    const int TCP_STREAM = 6;
    const int DATA_STREAM_ADD = 7;

    int eMode = STANDARD_STREAM;

    bool echoTimeFlag = false;
    double dT = 0.0;
    double rTolDt = 0.00001;
    bool doScientific = false;

    int precision = 6;

    bool closeOnWrite = false;

    const char *inetAddr = 0;
    int inetPort;

    int gradIndex = -1;

    ID nodes(0, 6);
    ID dofs(0, 6);
    ID timeseries(0, 6);

    while (OPS_GetNumRemainingInputArgs() > 0) {

        const char* option = OPS_GetString();
        responseID = option;

        if (strcmp(option, "-time") == 0) {
            echoTimeFlag = true;
        }
        else if (strcmp(option, "-load") == 0) {
            echoTimeFlag = true;
        }
        else if (strcmp(option, "-scientific") == 0) {
            doScientific = true;
        }
        else if (strcmp(option, "-file") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = DATA_STREAM;
        }
        else if (strcmp(option, "-fileAdd") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = DATA_STREAM_ADD;
        }
        else if (strcmp(option, "-closeOnWrite") == 0) {
            closeOnWrite = true;
        }
        else if (strcmp(option, "-csv") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = DATA_STREAM_CSV;
        }
        else if (strcmp(option, "-tcp") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                inetAddr = OPS_GetString();
            }
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &inetPort) < 0) {
                    opserr << "WARNING: failed to read inetPort\n";
                    return 0;
                }
            }
            eMode = TCP_STREAM;
        }
        else if (strcmp(option, "-xml") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = XML_STREAM;
        }
        else if (strcmp(option, "-binary") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = BINARY_STREAM;
        }
        else if (strcmp(option, "-precision") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &precision) < 0) {
                    opserr << "WARNING: failed to read precision\n";
                    return 0;
                }
            }
        }
        else if (strcmp(option, "-region") == 0) {
            int tag;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &tag) < 0) {
                    opserr << "WARNING: failed to read region tag\n";
                    return 0;
                }
            }
            // Domain *domain = OPS_GetDomain();
            // MeshRegion *theRegion = domain->getRegion(tag);
            // if (theRegion == 0) {
            //     opserr << "WARNING: region does not exist\n";
            //     return 0;
            // }
            // const ID &nodeRegion = theRegion->getNodes();
            // int numNodes = 0;
            // for (int i = 0; i < nodeRegion.Size(); i++)
            //     nodes[numNodes++] = nodeRegion(i);
        }
    }

    // data handler
    if (eMode == DATA_STREAM && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else if (eMode == DATA_STREAM_ADD && filename != 0)
        theOutputStream = new DataFileStreamAdd(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else if (eMode == DATA_STREAM_CSV && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 1, closeOnWrite, precision, doScientific);
    else if (eMode == XML_STREAM && filename != 0)
        theOutputStream = new XmlFileStream(filename);
    //else if (eMode == DATABASE_STREAM && tableName != 0)
    //    theOutputStream = new DatabaseStream(theDatabase, tableName);
    else if (eMode == BINARY_STREAM && filename != 0)
        theOutputStream = new BinaryFileStream(filename);
    else if (eMode == TCP_STREAM && inetAddr != 0)
        theOutputStream = new TCP_Stream(inetPort, inetAddr);
    else
        theOutputStream = new StandardStream();

    theOutputStream->setPrecision(precision);

    Domain* domain = OPS_GetDomain();
    if (domain == 0)
        return 0;
    EnergyBalanceRecorder* recorder = new EnergyBalanceRecorder(
        *domain, *theOutputStream,
        echoTimeFlag);

    return recorder;
}


EnergyBalanceRecorder::EnergyBalanceRecorder()
    : Recorder(RECORDER_TAGS_EnergyBalanceRecorder),
      theDomain(0), theOutputHandler(0),
      response(EBR_NUM_ENERGY_COMPONENTS),
      echoTimeFlag(true),
      initializationDone(false),
      time_last(0.),
      internal_energy(0.),
      damping_work(0.),
      unbalanced_load_work(0.)
{

}


EnergyBalanceRecorder::EnergyBalanceRecorder(
    Domain &theDomain,
    OPS_Stream &theOutputHandler_,
    bool echoTimeFlag_)
    : Recorder(RECORDER_TAGS_EnergyBalanceRecorder),
      theDomain(&theDomain), theOutputHandler(&theOutputHandler_),
      response(EBR_NUM_ENERGY_COMPONENTS),
      echoTimeFlag(echoTimeFlag_),
      initializationDone(false),
      time_last(0.),
      internal_energy(0.),
      damping_work(0.),
      unbalanced_load_work(0.)
{

}


EnergyBalanceRecorder::~EnergyBalanceRecorder()
{
    if (theOutputHandler != 0) {
        theOutputHandler->endTag(); // Data
        delete theOutputHandler;
    }
}


int
EnergyBalanceRecorder::record(int commitTag, double timeStamp)
{
    if (theDomain == 0 ) {
        return 0;
    }

    if (theOutputHandler == 0) {
        opserr << "EnergyBalanceRecorder::record() - no DataOutputHandler has been set\n";
        return -1;
    }

    if (initializationDone == false) {
        if (this->initialize() != 0) {
            opserr << "EnergyBalanceRecorder::record() - failed in initialize()\n";
            return -1;
        }
    }

    // where relDeltaTTol is the maximum reliable ratio between analysis time step and deltaT
    // to provide adequate tolerance for floating point precision (default=1.0e-5)

    //
    // if need nodal reactions get the domain to calculate them
    // before we iterate over the nodes
    //

    // NEED TO GET REACTIONS?
    // if (dataFlag == 7)
    //     theDomain->calculateNodalReactions(0);
    // else if (dataFlag == 8)
    //     theDomain->calculateNodalReactions(1);
    // if (dataFlag == 9)
    //     theDomain->calculateNodalReactions(2);

    //
    // add time information if requested
    //

    int timeOffset = 0;
    if (echoTimeFlag == true) {
        timeOffset = 1;
        response(0) = timeStamp;
    }


    // Compute Kinetic Energy 
    double kinetic_energy = 0.;
    {   
        Element * ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0)
        {
            const Matrix &M = ele->getMass();
            const int numExternalNodes = ele->getNumExternalNodes();
            const int numDOF = ele->getNumDOF();
            Node **elenodes = ele->getNodePtrs();
            static Vector vel(numDOF);
            vel.resize(numDOF);
            int cnt = 0;
            for (int i = 0; i < numExternalNodes; ++i)
            {
                Node * node = elenodes[i];
                const Vector& node_vel = node->getVel();
                for (int j = 0; j < node_vel.Size(); ++j)
                {
                    vel(cnt) = node_vel(j);
                    cnt++;
                }
            }
            for (int i = 0; i < numDOF; ++i)
            {
                for (int j = 0; j < numDOF; ++j)
                {
                    kinetic_energy += M(i, j) * vel(i) * vel(j) / 2;
                }
            }
        }
    }
    response(timeOffset + 0) = kinetic_energy;


    // Compute Internal Energy 
    double internal_energy_increment = 0.;
    double damping_work_increment = 0.;
    double dT = timeStamp - time_last;
    {   
        Element * ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0)
        {
            const Matrix &C = ele->getDamp();
            const Vector &F = ele->getResistingForce();
            const int numExternalNodes = ele->getNumExternalNodes();
            const int numDOF = ele->getNumDOF();
            Node **elenodes = ele->getNodePtrs();
            static Vector vel(numDOF);
            vel.resize(numDOF);
            int cnt = 0;
            for (int i = 0; i < numExternalNodes; ++i)
            {
                Node * node = elenodes[i];
                const Vector& node_vel = node->getVel();
                for (int j = 0; j < node_vel.Size(); ++j)
                {
                    vel(cnt) = node_vel(j);
                    cnt++;
                }
            }
            for (int i = 0; i < numDOF; ++i)
            {
                for (int j = 0; j < numDOF; ++j)
                {
                    damping_work_increment += C(i, j) * vel(i) * vel(j) ;
                }
            }

            internal_energy_increment += F ^ vel;
        }
        internal_energy_increment *= dT;
        internal_energy += internal_energy_increment;
        
        damping_work_increment *= dT;
        damping_work += damping_work_increment;
    }
    response(timeOffset + 1) = internal_energy;
    response(timeOffset + 2) = damping_work;


    //Compute the work done by unbalanced forces (external loads)
    double unbalanced_load_work_increment = 0;
    {
        Node* node;
        NodeIter& nodes = theDomain->getNodes();
        while((node=nodes()) != 0)
        {
            const Vector& vel = node->getVel();
            const Vector& unbalanced_load = node->getUnbalancedLoad();

            unbalanced_load_work_increment += vel ^ unbalanced_load;
        }
        unbalanced_load_work_increment *= dT;
        unbalanced_load_work += unbalanced_load_work_increment;
    }
    response(timeOffset + 3) = unbalanced_load_work;


    //Tell the handler to write the output... now
    theOutputHandler->write(response);
    
    //Advance the step
    time_last = timeStamp;

    return 0;
}


int
EnergyBalanceRecorder::setDomain(Domain &theDom)
{
    theDomain = &theDom;
    time_last = theDomain->getCurrentTime();
    return 0;
}


int
EnergyBalanceRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  

    return 0;
}


int
EnergyBalanceRecorder::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker)
{
   

    return 0;
}


int
EnergyBalanceRecorder::domainChanged(void)
{
    return 0;
}


int
EnergyBalanceRecorder::initialize(void)
{
    if (theDomain == 0) {
        opserr << "EnergyBalanceRecorder::initialize() - the domain has not been set\n";
        return -1;
    }


    //
    // resize the response vector
    //

    int timeOffset = 0;
    if (echoTimeFlag == true)
    {
        timeOffset = 1;
        response.resize(EBR_NUM_ENERGY_COMPONENTS + 1);
    }


    // ID orderResponse(numValidResponse);

    // //
    // // need to create the data description, i.e. what each column of data is
    // //

    // char outputData[32];
    // char dataType[10];

    // if (dataFlag == 0) {
    //     strcpy(dataType, "D");
    // } else if (dataFlag == 1) {
    //     strcpy(dataType, "V");
    // } else if (dataFlag == 2) {
    //     strcpy(dataType, "A");
    // } else if (dataFlag == 3) {
    //     strcpy(dataType, "dD");
    // } else if (dataFlag == 4) {
    //     strcpy(dataType, "ddD");
    // } else if (dataFlag == 5) {
    //     strcpy(dataType, "U");
    // } else if (dataFlag == 6) {
    //     strcpy(dataType, "U");
    // } else if (dataFlag == 7) {
    //     strcpy(dataType, "R");
    // } else if (dataFlag == 8) {
    //     strcpy(dataType, "R");
    // } else if (dataFlag == 10000) {
    //     strcpy(dataType, "|D|");
    // } else if (dataFlag > 10) {
    //     sprintf(dataType, "E%d", dataFlag - 10);
    // } else
    //     strcpy(dataType, "Unknown");

    // /************************************************************
    // } else if ((strncmp(dataToStore, "sensitivity",11) == 0)) {
    //   int grad = atoi(&(dataToStore[11]));
    //   if (grad > 0)
    //     dataFlag = 1000 + grad;
    //   else
    //     dataFlag = 6;
    // } else if ((strncmp(dataToStore, "velSensitivity",14) == 0)) {
    //   int grad = atoi(&(dataToStore[14]));
    //   if (grad > 0)
    //     dataFlag = 2000 + grad;
    //   else
    //     dataFlag = 6;
    // } else if ((strncmp(dataToStore, "accSensitivity",14) == 0)) {
    //   int grad = atoi(&(dataToStore[14]));
    //   if (grad > 0)
    //     dataFlag = 3000 + grad;
    //   else
    //     dataFlag = 6;

    // ***********************************************************/
    // int numDOF = theDofs->Size();

    // // write out info to handler if parallel execution
    // //

    // ID xmlOrder(numValidNodes);

    // if (echoTimeFlag == true)
    //     xmlOrder.resize(numValidNodes + 1);

    // if (theNodalTags != 0 && addColumnInfo == 1) {

    //     int numNode = theNodalTags->Size();
    //     int count = 0;
    //     int nodeCount = 0;

    //     if (echoTimeFlag == true)  {
    //         orderResponse(count++) = 0;
    //         xmlOrder(nodeCount++) = 0;
    //     }

    //     for (int i = 0; i < numNode; i++) {
    //         int nodeTag = (*theNodalTags)(i);
    //         Node *theNode = theDomain->getNode(nodeTag);
    //         if (theNode != 0) {
    //             xmlOrder(nodeCount++) = i + 1;
    //             for (int j = 0; j < numDOF; j++)
    //                 orderResponse(count++) = i + 1;
    //         }
    //     }

    //     theOutputHandler->setOrder(xmlOrder);
    // }

    // char nodeCrdData[20];
    // sprintf(nodeCrdData, "coord");

    // if (echoTimeFlag == true) {
    //     if (theNodalTags != 0 && addColumnInfo == 1) {
    //         theOutputHandler->tag("TimeOutput");
    //         theOutputHandler->tag("ResponseType", "time");
    //         theOutputHandler->endTag();
    //     }
    // }

    // for (int i = 0; i < numValidNodes; i++) {
    //     int nodeTag = theNodes[i]->getTag();
    //     const Vector &nodeCrd = theNodes[i]->getCrds();
    //     int numCoord = nodeCrd.Size();


    //     theOutputHandler->tag("NodeOutput");
    //     theOutputHandler->attr("nodeTag", nodeTag);

    //     for (int j = 0; j < 3; j++) {
    //         sprintf(nodeCrdData, "coord%d", j + 1);
    //         if (j < numCoord)
    //             theOutputHandler->attr(nodeCrdData, nodeCrd(j));
    //         else
    //             theOutputHandler->attr(nodeCrdData, 0.0);
    //     }

    //     for (int k = 0; k < theDofs->Size(); k++) {
    //         sprintf(outputData, "%s%d", dataType, k + 1);
    //         theOutputHandler->tag("ResponseType", outputData);
    //     }

    //     theOutputHandler->endTag();
    // }

    // if (theNodalTags != 0 && addColumnInfo == 1) {
    //     theOutputHandler->setOrder(orderResponse);
    // }

    // theOutputHandler->tag("Data");
    // initializationDone = true;

    return 0;
}


int EnergyBalanceRecorder::flush(void) {
    if (theOutputHandler != 0) {
        return theOutputHandler->flush();
    }
    return 0;
}