/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

// ============================================================================
// 2021 By Jose Abell @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// Implements simple routines to test nDMaterials using direct strains.
//
// ============================================================================

#include "TimeSeries.h"
#include "map"
#include "string.h"
#include "elementAPI.h"
#include "vector"
#include "PathTimeSeries.h"
#include "PathSeries.h"
#include "Vector.h"
#include "Matrix.h"
#include "NDMaterial.h"
#include "Information.h"
#include "Response.h"
#include "MaterialResponse.h"

// Materials created by NDTest Copy / SendRecv live in this table (they keep
// the tag of their source, which cannot be changed from here); every NDTest
// command looks a tag up here first, then in the global material list.
static std::map<int, NDMaterial *> &ndTestCopies()
{
    static std::map<int, NDMaterial *> copies;
    return copies;
}

static NDMaterial *ndTestLookup(int tag)
{
    std::map<int, NDMaterial *>::iterator it = ndTestCopies().find(tag);
    if (it != ndTestCopies().end())
        return it->second;
    return OPS_getNDMaterial(tag);
}

static void ndTestStore(int tag, NDMaterial *m)
{
    std::map<int, NDMaterial *>::iterator it = ndTestCopies().find(tag);
    if (it != ndTestCopies().end())
        delete it->second;
    ndTestCopies()[tag] = m;
}

//Forward declaration for indexing purposes
void * OPS_NDSetStrain();               // To set a trial strain on the material
void * OPS_NDPrintStrain();             // Printing strains to the screen
void * OPS_NDPrintStress();             // Printing strains to the screen
void * OPS_NDGetStrain();               // Return current strains to a list variable
void * OPS_NDGetStress();               // Return current stresses to a list variable
void * OPS_NDGetTangentStiffness();     // Return current tangent stiffness to a list variable
void * OPS_NDCommitState();             // To commit the state of the material after last trial strain
void * OPS_NDUpdateIntegerParameter();  // To change an integer parameter in the material
void * OPS_NDUpdateDoubleParameter();   // To change a double parameter in the material
void * OPS_NDGetResponse();             // To get other responses

void* OPS_NDSetStrain()
{
    int tag = 0;
    double strains[6];

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0)
    {
        opserr << "OPS_NDSetStrain - got incorrect integer tag for material" << endln;
        return 0;
    }


    numdata = 6;
    if (OPS_GetDoubleInput(&numdata, strains) < 0) {
        opserr << "OPS_NDSetStrain - need 6 components of floating-point strains" << endln;
        return 0;
    }

    NDMaterial* mat = ndTestLookup(tag);

    if (mat == NULL)
    {
        opserr << "OPS_NDSetStrain - material with tag " << tag << " does not exist" << endln;
        return 0;
    }


    Vector new_strain(6);

    new_strain(0) = strains[0];
    new_strain(1) = strains[1];
    new_strain(2) = strains[2];
    new_strain(3) = strains[3];
    new_strain(4) = strains[4];
    new_strain(5) = strains[5];

    mat->setTrialStrain(new_strain);

    return 0;
}

void* OPS_NDPrintStrain()
{
    int tag = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = ndTestLookup(tag);

    if (mat == NULL)
    {
        opserr << "OPS_NDPrintStrain - material with tag " << tag << " does not exist" << endln;
        return 0;
    }

    const Vector &strain = mat->getStrain();

    return 0;
}

void* OPS_NDPrintStress()
{
    int tag = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = ndTestLookup(tag);

    if (mat == NULL)
    {
        opserr << "OPS_NDPrintStress - material with tag " << tag << " does not exist" << endln;
        return 0;
    }

    const Vector &stress = mat->getStress();

    return 0;
}


void* OPS_NDGetStrain()
{
    int tag = 0;
    int size = 6;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = ndTestLookup(tag);

    if (mat == NULL)
    {
        opserr << "OPS_NDGetStrain - material with tag " << tag << " does not exist" << endln;
        return 0;
    }

    const Vector &strain = mat->getStrain();

    std::vector<double> values(size);
    for (int i = 0; i < 6; i++) {
        values[i] = strain(i);
    }
    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
        opserr << "WARNING OPS_NDGetStress - failed to set double inputs\n";
        return 0;
    }


    return 0;
}



void* OPS_NDGetStress()
{
    int tag = 0;
    int size = 6;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0)
    {
        return 0;
    }


    NDMaterial* mat = ndTestLookup(tag);
    if (mat == NULL)
    {
        opserr << "OPS_NDGetStress() - Material tag " << tag << " not declared" << endln;
    }

    const Vector &stress = mat->getStress();

    std::vector<double> values(size);
    for (int i = 0; i < 6; i++) {
        values[i] = stress(i);
    }
    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
        opserr << "WARNING OPS_NDGetStress - failed to set double inputs\n";
        return 0;
    }


    return 0;
}

void* OPS_NDGetResponse()
{
    int tag = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0)
    {
        return 0;
    }


    NDMaterial* mat = ndTestLookup(tag);
    if (mat == NULL)
    {
        opserr << "OPS_NDGetStress() - Material tag " << tag << " not declared" << endln;
    }


    const char* response_type = OPS_GetString();

    Response* theResponse = mat->setResponse(&response_type, 1, opserr);
    theResponse->getResponse();
    Information& theInformation = theResponse->getInformation();

    const Vector &respose = *theInformation.theVector;

    int size = respose.Size();
    std::vector<double> values(size);
    for (int i = 0; i < size; i++) {
        values[i] = respose(i);
    }
    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
        opserr << "WARNING OPS_NDGetStress - failed to set double inputs\n";
        return 0;
    }


    return 0;
}



void* OPS_NDGetTangentStiffness()
{
    int tag = 0;
    int size = 36;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = ndTestLookup(tag);

    if (mat == NULL)
    {
        opserr << "OPS_NDGetTangentStiffness - material with tag " << tag << " does not exist" << endln;
        return 0;
    }

    const Matrix &stiffness = mat->getTangent();

    std::vector<double> values(size);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++)
        {
            double one_value = stiffness(i, j);
            values[6 * i + j] = one_value;
        }
    }
    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
        opserr << "WARNING OPS_NDGetStress - failed to set double inputs\n";
        return 0;
    }


    return 0;
}



void* OPS_NDCommitState()
{
    int tag = 0;
    double strains[6];
    int size = 6;
    double stressdata[6] = {1 , 2 , 3, 4, 5, 6};

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;


    NDMaterial* mat = ndTestLookup(tag);

    if (mat == NULL)
    {
        opserr << "OPS_NDCommitState - material with tag " << tag << " does not exist" << endln;
        return 0;
    }


    mat->commitState();

    return 0;
}


void* OPS_NDUpdateIntegerParameter()
{

    opserr << "OPS_NDUpdateIntegerParameter" << endln;

    int tag = 0;
    int responseID = 0;
    int theNewIntegerParameterValue = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;
    opserr << "tag = " << tag  << endln;
    if (OPS_GetIntInput(&numdata, &responseID) < 0) return 0;
    opserr << "responseID = " << responseID  << endln;
    if (OPS_GetIntInput(&numdata, &theNewIntegerParameterValue) < 0) return 0;
    opserr << "theNewIntegerParameterValue = " << theNewIntegerParameterValue  << endln;


    NDMaterial* mat = ndTestLookup(tag);

    if (mat == NULL)
    {
        opserr << "OPS_getNDMaterial - material with tag " << tag << " does not exist" << endln;
        return 0;
    }

    Information info;

    info.theInt = theNewIntegerParameterValue;

    mat->updateParameter(responseID, info);

    return 0;
}


void* OPS_NDUpdateDoubleParameter()
{
    int tag = 0;
    int responseID = 0;
    double theNewDoubleParameterValue = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;
    if (OPS_GetIntInput(&numdata, &responseID) < 0) return 0;
    if (OPS_GetDoubleInput(&numdata, &theNewDoubleParameterValue) < 0) return 0;


    NDMaterial* mat = ndTestLookup(tag);

    if (mat == NULL)
    {
        opserr << "OPS_NDUpdateDoubleParameter - material with tag " << tag << " does not exist" << endln;
        return 0;
    }

    Information info;

    info.theDouble = theNewDoubleParameterValue;

    mat->updateParameter(responseID, info);

    return 0;
}


// State-operation commands (added for TimeVaryingEC verification). Each
// returns an integer code (0 = success) so tests can check failures.
#include <FileDatastore.h>
#include <FEM_ObjectBrokerAllClasses.h>

static NDMaterial *ndTestMaterial(const char *who)
{
    int tag = 0;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << who << " - expects a material tag\n";
        return 0;
    }
    NDMaterial *mat = ndTestLookup(tag);
    if (mat == 0)
        opserr << who << " - material with tag " << tag << " does not exist\n";
    return mat;
}

static void ndTestReturn(int code)
{
    int one = 1;
    OPS_SetIntOutput(&one, &code, true);
}

// NDTest TrySetStrain tag e1..e6 -> setTrialStrain return code
void *OPS_NDTrySetStrain()
{
    NDMaterial *mat = ndTestMaterial("NDTest TrySetStrain");
    if (mat == 0) { ndTestReturn(-99); return 0; }
    double e[6];
    int numdata = 6;
    if (OPS_GetDoubleInput(&numdata, e) < 0) { ndTestReturn(-98); return 0; }
    Vector v(e, 6);
    ndTestReturn(mat->setTrialStrain(v));
    return 0;
}

// NDTest TryCommit tag -> commitState return code
void *OPS_NDTryCommit()
{
    NDMaterial *mat = ndTestMaterial("NDTest TryCommit");
    ndTestReturn(mat ? mat->commitState() : -99);
    return 0;
}

// NDTest RevertToLastCommit tag
void *OPS_NDRevertToLastCommit()
{
    NDMaterial *mat = ndTestMaterial("NDTest RevertToLastCommit");
    ndTestReturn(mat ? mat->revertToLastCommit() : -99);
    return 0;
}

// NDTest RevertToStart tag
void *OPS_NDRevertToStart()
{
    NDMaterial *mat = ndTestMaterial("NDTest RevertToStart");
    ndTestReturn(mat ? mat->revertToStart() : -99);
    return 0;
}

// NDTest Copy srcTag newTag ?type?: getCopy() (or getCopy(type)) of srcTag,
// reachable by NDTest as newTag
void *OPS_NDCopy()
{
    NDMaterial *mat = ndTestMaterial("NDTest Copy");
    int newTag = 0, numdata = 1;
    if (mat == 0 || OPS_GetIntInput(&numdata, &newTag) < 0) { ndTestReturn(-99); return 0; }
    // optional type: getCopy(type), e.g. ThreeDimensional for a base-class material
    NDMaterial *c = 0;
    if (OPS_GetNumRemainingInputArgs() > 0)
        c = mat->getCopy(OPS_GetString());
    else
        c = mat->getCopy();
    if (c == 0) { ndTestReturn(-1); return 0; }
    ndTestStore(newTag, c);
    ndTestReturn(0);
    return 0;
}

// NDTest SendRecv srcTag newTag fileBase: sendSelf of srcTag into a
// FileDatastore, then a new object from the broker by class tag and
// recvSelf; reachable by NDTest as newTag
void *OPS_NDSendRecv()
{
    NDMaterial *mat = ndTestMaterial("NDTest SendRecv");
    int newTag = 0, numdata = 1;
    if (mat == 0 || OPS_GetIntInput(&numdata, &newTag) < 0 || OPS_GetNumRemainingInputArgs() < 1) {
        ndTestReturn(-99);
        return 0;
    }
    const char *base = OPS_GetString();
    static FEM_ObjectBrokerAllClasses theBroker;
    Domain *dom = OPS_GetDomain();
    FileDatastore *ds = new FileDatastore(base, *dom, theBroker);
    int dbTag = ds->getDbTag();
    int commitTag = 1;
    mat->setDbTag(dbTag);
    if (mat->sendSelf(commitTag, *ds) < 0) { delete ds; ndTestReturn(-3); return 0; }
    NDMaterial *r = theBroker.getNewNDMaterial(mat->getClassTag());
    if (r == 0) { delete ds; ndTestReturn(-4); return 0; }
    r->setDbTag(dbTag);
    if (r->recvSelf(commitTag, *ds, theBroker) < 0) { delete r; delete ds; ndTestReturn(-5); return 0; }
    delete ds;
    ndTestStore(newTag, r);
    ndTestReturn(0);
    return 0;
}

// NDTest GetInitialTangent tag -> 36 values, row major
void *OPS_NDGetInitialTangent()
{
    NDMaterial *mat = ndTestMaterial("NDTest GetInitialTangent");
    if (mat == 0) return 0;
    const Matrix &C = mat->getInitialTangent();
    int size = C.noRows() * C.noCols();
    std::vector<double> values(size);
    for (int i = 0; i < C.noRows(); i++)
        for (int j = 0; j < C.noCols(); j++)
            values[i * C.noCols() + j] = C(i, j);
    OPS_SetDoubleOutput(&size, &values[0], false);
    return 0;
}

namespace {

struct char_cmp {
    bool operator () (const char *a, const char *b) const
    {
        return strcmp(a, b) < 0;
    }
};

typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;


static OPS_ParsingFunctionMap functionMap;

static int setUpFunctions(void)
{
    functionMap.insert(std::make_pair("SetStrain", &OPS_NDSetStrain));
    functionMap.insert(std::make_pair("CommitState", &OPS_NDCommitState));
    functionMap.insert(std::make_pair("PrintStress", &OPS_NDPrintStress));
    functionMap.insert(std::make_pair("PrintStrain", &OPS_NDPrintStrain));
    functionMap.insert(std::make_pair("GetStrain", &OPS_NDGetStrain));
    functionMap.insert(std::make_pair("GetStress", &OPS_NDGetStress));
    functionMap.insert(std::make_pair("GetTangentStiffness", &OPS_NDGetTangentStiffness));
    functionMap.insert(std::make_pair("UpdateIntegerParameter", &OPS_NDUpdateIntegerParameter));
    functionMap.insert(std::make_pair("UpdateDoubleParameter", &OPS_NDUpdateDoubleParameter));
    functionMap.insert(std::make_pair("GetResponse", &OPS_NDGetResponse));
    functionMap.insert(std::make_pair("TrySetStrain", &OPS_NDTrySetStrain));
    functionMap.insert(std::make_pair("TryCommit", &OPS_NDTryCommit));
    functionMap.insert(std::make_pair("RevertToLastCommit", &OPS_NDRevertToLastCommit));
    functionMap.insert(std::make_pair("RevertToStart", &OPS_NDRevertToStart));
    functionMap.insert(std::make_pair("Copy", &OPS_NDCopy));
    functionMap.insert(std::make_pair("SendRecv", &OPS_NDSendRecv));
    functionMap.insert(std::make_pair("GetInitialTangent", &OPS_NDGetInitialTangent));

    return 0;
}
}

int
OPS_NDTest()
{
    static bool initDone = false;
    if (initDone == false) {
        setUpFunctions();
        initDone = true;
    }

    // Identify what specific command of Patch we're calling
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING too few arguments: NDTest cmd? \n";
        opserr << " available commands: SetStrain|CommitState|GetStrain|GetStress \n";
        return -1;
    }

    const char* type = OPS_GetString();

    OPS_ParsingFunctionMap::const_iterator iter = functionMap.find(type);
    if (iter == functionMap.end()) {
        opserr << "WARNING NDTest type " << type << " is unknown\n";
        return -1;
    }

    // Call the function
    (*iter->second)();

    return 0;

}
