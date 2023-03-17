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

// ============================================================================
// 2023 By Jose Abell and José Larenas @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// Please read detailed description in SixNodeBoundryCondition.h.
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <SixNodeBoundryCondition.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <map>

void* OPS_SixNodeBoundryCondition()
{
    if (OPS_GetNumRemainingInputArgs() < 7)
    {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: element SixNodeBoundryCondition eleTag? Node1? Node2? Node3? Node4? Node5? Node6? beta, k, tamb, th\n";
        return 0;
    }

    int idata[7];
    int num = 7;
    if (OPS_GetIntInput(&num, idata) < 0)
    {
        opserr << "WARNING: invalid integer data\n";
        return 0;
    }

    double data[4] = {0, 0, 0, 0};
    num = OPS_GetNumRemainingInputArgs();

    if (num > 4)
    {
        num = 4;
    }
    if (num > 0)
    {
        if (OPS_GetDoubleInput(&num, data) < 0)
        {
            opserr << "WARNING: invalid double data\n";
            return 0;
        }
    }

    return new SixNodeBoundryCondition(idata[0], idata[1], idata[2], idata[3], idata[4], idata[5], idata[6], data[0], data[1], data[2], data[3]);
}

//static data
double  SixNodeBoundryCondition::xl[3][NumNodes]                     ;
Matrix  SixNodeBoundryCondition::stiff(NumDOFsTotal, NumDOFsTotal)   ;
Vector  SixNodeBoundryCondition::resid(NumDOFsTotal)                 ;
Matrix  SixNodeBoundryCondition::mass(NumDOFsTotal, NumDOFsTotal)    ;

//quadrature data
const double  SixNodeBoundryCondition::alpha = 2.0 / 3.0 ;
const double  SixNodeBoundryCondition::beta  = 1.0 / 6.0 ;
const double  SixNodeBoundryCondition::sg[]  = { alpha, beta, beta } ;
const double  SixNodeBoundryCondition::wg[]  = { 1.0 / 3.0 } ;

// static Matrix B(NumStressComponents, NumDOFsPerNode) ;
Matrix SixNodeBoundryCondition::B(NumStressComponents, NumDOFsPerNode) ;

//null constructor
SixNodeBoundryCondition::SixNodeBoundryCondition( )
    : Element( 0, ELE_TAG_SixNodeBoundryCondition ),
      connectedExternalNodes(NumNodes),
      applyLoad(0), load(0), Ki(0)
{
    B.Zero();

    inp_info[0] = 0.0;
    inp_info[1] = 0.0;
    inp_info[2] = 0.0;
    inp_info[3] = 1.0;

    for (int i = 0; i < NumNodes; i++ ) {
        nodePointers[i] = 0;
    }
}

//*********************************************************************

//full constructor
SixNodeBoundryCondition::SixNodeBoundryCondition(int tag,
        int node1,
        int node2,
        int node3,
        int node4,
        int node5,
        int node6,
        double beta,
        double k,
        double tamb,
        double th)
    : Element(tag, ELE_TAG_SixNodeBoundryCondition),
      connectedExternalNodes(NumNodes), applyLoad(0), load(0), Ki(0)
{
    B.Zero();
    connectedExternalNodes(0) = node1 ;
    connectedExternalNodes(1) = node2 ;
    connectedExternalNodes(2) = node3 ;
    connectedExternalNodes(3) = node4 ;
    connectedExternalNodes(4) = node5 ;
    connectedExternalNodes(5) = node6 ;

    for (int i = 0; i < NumNodes; i++ ) {
        nodePointers[i] = 0;
    }

    inp_info[0] = beta ;
    inp_info[1] = k ;
    inp_info[2] = tamb ;
    inp_info[3] = th ;
}

//******************************************************************

//destructor
SixNodeBoundryCondition::~SixNodeBoundryCondition( )
{
    if (load != 0)
        delete load;

    if (Ki != 0)
        delete Ki;
}

//set domain
void  SixNodeBoundryCondition::setDomain( Domain *theDomain )
{
    int i ;

    //node pointers
    for ( i = 0; i < NumNodes; i++ )
    {
        nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
    }

    this->DomainComponent::setDomain(theDomain);
}

//get the number of external nodes
int  SixNodeBoundryCondition::getNumExternalNodes( ) const
{
    return NumNodes ;
}

//return connected external nodes
const ID&  SixNodeBoundryCondition::getExternalNodes( )
{
    return connectedExternalNodes ;
}

Node **
SixNodeBoundryCondition::getNodePtrs(void)
{
    return nodePointers ;
}

//return number of dofs
int  SixNodeBoundryCondition::getNumDOF( )
{
    return NumDOFsTotal ;
}

//commit state
int  SixNodeBoundryCondition::commitState( )
{
    int success = 0 ;

    // call element commitState to do any base class stuff
    if ((success = this->Element::commitState()) != 0) {
        opserr << "SixNodeBoundryCondition::commitState () - failed in base class";
    }

    return success ;
}

//revert to last commit
int  SixNodeBoundryCondition::revertToLastCommit( )
{
    int success = 0 ;

    return success ;
}

//revert to start
int  SixNodeBoundryCondition::revertToStart( )
{
    int success = 0 ;

    return success ;
}

//print out element data
void  SixNodeBoundryCondition::Print(OPS_Stream &s, int flag)
{
    if (flag == 2) {
        s << "#SixNodeBoundryCondition\n";

        int i;
        const int numNodes = NumNodes;
        const int nstress = NumStressComponents;

        for (i = 0; i < numNodes; i++) {
            const Vector &nodeCrd = nodePointers[i]->getCrds();
            const Vector &nodeDisp = nodePointers[i]->getDisp();
            s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)
              << " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;
        }
    }

    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "Standard SixNodeBoundryCondition \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Nodes: " << connectedExternalNodes;
        s << endln;
        s << "Resisting Force (no inertia): " << this->getResistingForce();
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"SixNodeBoundryCondition\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        for (int i = 1; i < 2; i++)
            s << connectedExternalNodes(i) << ", ";
        s << connectedExternalNodes(3) << "], ";
    }
}

//return stiffness matrix
const Matrix&  SixNodeBoundryCondition::getTangentStiff( )
{
    int tang_flag = 1 ; //get the tangent

    //do tangent and residual here
    formResidAndTangent( tang_flag ) ;

    return stiff ;
}

//return initial matrix
const Matrix&  SixNodeBoundryCondition::getInitialStiff( )
{
    if (Ki != 0)
        return *Ki;


    static const int ndm = 3 ;

    static const int ndf = NumDOFsPerNode ;

    static const int numberNodes = NumNodes ;

    static const int numberGauss = NumGaussPoints ;

    static const int nShape = 3 ;

    static const int massIndex = nShape - 1 ;

    static const int stiffIndex = nShape - 1 ;

    static double volume ;

    static double xsj ;  // determinant jacaobian matrix

    static double dvol[numberGauss] ; //volume element

    static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

    static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

    static double gaussPoint[2] ;

    static Vector momentum(ndf) ;

    int i, j, k, p, q ;
    int jj, kk ;

    static double temp, stiffJK ;

    stiff.Zero( ) ;
    resid.Zero( ) ;

    //compute basis vectors and local nodal coordinates
    computeBasis( ) ;

    //gauss loop to compute and save shape functions
    int count = 0 ;
    volume = 0.0 ;

    for ( k = 0; k < numberGauss; k++ )
    {
        gaussPoint[0] = sg[k] ;
        gaussPoint[1] = sg[abs(1 - k)] ;

        //get shape functions
        shp3d( gaussPoint, xsj, shp, xl ) ;

        //save shape functions
        for ( p = 0; p < nShape; p++ )
        {
            for ( q = 0; q < numberNodes; q++ )
            {
                Shape[p][q][count] = shp[p][q] ;
            }
        } // end for p

        //volume element to also be saved
        dvol[count] = wg[0] * xsj ;
        volume += dvol[count] ;
        count++ ;
    } //end for k

    //gauss loop
    for ( i = 0; i < numberGauss; i++ )
    {
        //extract shape functions from saved array
        for ( p = 0; p < nShape; p++ )
        {
            for ( q = 0; q < numberNodes; q++ )
            {
                shp[p][q]  = Shape[p][q][i] ;
            }
        } // end for p

        //residual and tangent calculations node loops
        jj = 0 ;
        for ( j = 0; j < numberNodes; j++ )
        {
            temp = shp[stiffIndex][j] * dvol[i] * (inp_info[0] / inp_info[1]) * inp_info[3] ;

            //node-node mass
            kk = 0 ;
            for ( k = 0; k < numberNodes; k++ )
            {
                stiffJK = temp * shp[stiffIndex][k] ;
                for ( p = 0; p < ndf; p++ )
                {
                    stiff( jj + p, kk + p ) += stiffJK ;
                }
                kk += ndf ;
            } // end for k loop

            jj += ndf ;
        } // end for j loop
    } //end for i gauss loop

    Ki = new Matrix(stiff);

    return stiff ;
}

//return mass matrix
const Matrix&  SixNodeBoundryCondition::getMass( )
{
    int tangFlag = 1 ;

    formInertiaTerms( tangFlag ) ;

    return mass ;
}

void  SixNodeBoundryCondition::zeroLoad( )
{
    if (load != 0)
        load->Zero();

    applyLoad = 0 ;
    appliedQ = 0 ; 

    return ;
}

int
SixNodeBoundryCondition::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_BrickSelfWeight) {
        applyLoad = 1;
        appliedQ += loadFactor * inp_info[2] * (inp_info[0] / inp_info[1]);
        return 0;
    } else if (type == LOAD_TAG_SelfWeight) {
        // added compatibility with selfWeight class implemented for all continuum elements, C.McGann, U.W.
        applyLoad = 1;
        appliedQ += loadFactor * data(2) * inp_info[2] * (inp_info[0] / inp_info[1]);
        return 0;
    } else {
        opserr << "SixNodeBoundryCondition::addLoad() - ele with tag: " << this->getTag() << " does not deal with load type: " << type << "\n";
        return -1;
    }

    return -1;
}

int
SixNodeBoundryCondition::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

//get residual
const Vector&  SixNodeBoundryCondition::getResistingForce( )
{
    int tang_flag = 1 ; // get the tangent

    resid.Zero();
    formResidAndTangent( tang_flag ) ;

    if (load != 0)
        resid -= *load;

    return resid ;
}

//get residual with inertia terms
const Vector&  SixNodeBoundryCondition::getResistingForceIncInertia( )
{
    static Vector res(6); res.Zero();

    int tang_flag = 1 ; // get the tangent
    resid.Zero();

    //do tangent and residual here
    formResidAndTangent( tang_flag ) ;

    formInertiaTerms( tang_flag ) ;

    res = resid;

    if (load != 0)
        res -= *load;

    return res;
}

//*********************************************************************

//form inertia terms
void   SixNodeBoundryCondition::formInertiaTerms( int tangFlag )
{
    mass.Zero( ) ;
}

//form residual and tangent
void   SixNodeBoundryCondition::formResidAndTangent( int tangFlag )
{
    static const int ndm = 3 ;

    static const int ndf = NumDOFsPerNode ;

    static const int numberNodes = NumNodes ;

    static const int numberGauss = NumGaussPoints ;

    static const int nShape = 3 ;

    static const int massIndex = nShape - 1 ;

    static const int stiffIndex = nShape - 1 ;

    double xsj ;  // determinant jacaobian matrix

    double dvol[numberGauss] ; //volume element

    static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

    static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

    static double gaussPoint[2] ;

    static Vector momentum(ndf) ;

    int i, j, k, p, q ;
    int jj, kk ;

    double temp, stiffJK ;

    stiff.Zero( ) ;

    //compute basis vectors and local nodal coordinates
    computeBasis( ) ;

    //gauss loop to compute and save shape functions
    int count = 0 ;

    for ( k = 0; k < numberGauss; k++ )
    {
        gaussPoint[0] = sg[k] ;
        gaussPoint[1] = sg[abs(1 - k)] ;

        //get shape functions
        shp3d( gaussPoint, xsj, shp, xl ) ;

        //save shape functions
        for ( p = 0; p < nShape; p++ )
        {
            for ( q = 0; q < numberNodes; q++ )
            {
                Shape[p][q][count] = shp[p][q] ;
            }
        } // end for p

        //volume element to also be saved
        dvol[count] = wg[0] * xsj ;
        count++ ;

    } //end for k

    //gauss loop
    for ( i = 0; i < numberGauss; i++ )
    {
        //extract shape functions from saved array
        for ( p = 0; p < nShape; p++ )
        {
            for ( q = 0; q < numberNodes; q++ )
            {
                shp[p][q]  = Shape[p][q][i] ;
            }
        } // end for p

        //residual and tangent calculations node loops
        jj = 0 ;
        for ( j = 0; j < numberNodes; j++ )
        {
            resid( jj  ) -= dvol[i] * ( inp_info[2] * (inp_info[0] / inp_info[1]) ) * shp[2][j] ;
            temp = shp[stiffIndex][j] * dvol[i] * (inp_info[0] / inp_info[1]) * inp_info[3] ;

            if ( tangFlag == 1 )
            {
                //node-node mass
                kk = 0 ;
                for ( k = 0; k < numberNodes; k++ )
                {
                    stiffJK = temp * shp[stiffIndex][k] ;
                    for ( p = 0; p < ndf; p++ )
                    {
                        stiff( jj + p, kk + p ) += stiffJK ;
                    }
                    kk += ndf ;
                } // end for k loop
            } // end if tang_flag

            jj += ndf ;
        } // end for j loop
    } //end for i gauss loop

    Vector nodeTemp_dot(NumNodes);

    for (int i = 0; i < NumNodes; ++i)
    {
        const Vector &temp_dot_node_i = nodePointers[i]->getTrialDisp( ) ;
        nodeTemp_dot(i) = temp_dot_node_i(0);
    }

    resid.addMatrixVector(0.0, stiff, nodeTemp_dot,  1.0);
}

//*********************************************************************

//form residual and tangent
int
SixNodeBoundryCondition::update(void)
{
    return 0;
}


//************************************************************************

//compute local coordinates and basis
void   SixNodeBoundryCondition::computeBasis( )
{
    //nodal coordinates
    int i ;
    for ( i = 0; i < NumNodes; i++ ) {
        const Vector &coorI = nodePointers[i]->getCrds( ) ;
        xl[0][i] = coorI(0) ;
        xl[1][i] = coorI(1) ;
        xl[2][i] = coorI(2) ;
    }  //end for i
}

//*************************************************************************

//compute B
const Matrix&
SixNodeBoundryCondition::computeB( int node, const double shp[3][NumNodes] )
{
    //--------------------------------------------------------------------
    //
    //                -    -
    //               | N,x  |
    //   B       =   | N,y  | (2x1)
    //                -    -
    //
    //-------------------------------------------------------------------

    B(0, 0) = shp[0][node] ;
    B(1, 0) = shp[1][node] ;

    return B ;
}

//**********************************************************************

int  SixNodeBoundryCondition::sendSelf (int commitTag, Channel &theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // Quad packs its data into a Vector and sends this to theChannel
    // along with its dbTag and the commitTag passed in the arguments

    // Now quad sends the ids of its materials
    static ID idData(26);

    idData(20) = connectedExternalNodes(0);
    idData(21) = connectedExternalNodes(1);
    idData(22) = connectedExternalNodes(2);
    idData(23) = connectedExternalNodes(3);
    idData(24) = connectedExternalNodes(4);
    idData(25) = connectedExternalNodes(5);

    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING SixNodeBoundryCondition::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    static Vector dData(7);
    dData(0) = alphaM;
    dData(1) = betaK;
    dData(2) = betaK0;
    dData(3) = betaKc;
    dData(4) = inp_info[0];
    dData(5) = inp_info[1];
    dData(6) = inp_info[2];

    if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
        opserr << "SixNodeBoundryCondition::sendSelf() - failed to send double data\n";
        return -1;
    }

    return res;
}

int  SixNodeBoundryCondition::recvSelf (int commitTag,
                                   Channel &theChannel,
                                   FEM_ObjectBroker &theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    static ID idData(26);
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING SixNodeBoundryCondition::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    this->setTag(idData(26));

    static Vector dData(7);
    if (theChannel.recvVector(dataTag, commitTag, dData) < 0) {
        opserr << "DispBeamColumn2d::sendSelf() - failed to recv double data\n";
        return -1;
    }

    alphaM = dData(0);
    betaK  = dData(1);
    betaK0 = dData(2);
    betaKc = dData(3);
    inp_info[0] = dData(4);
    inp_info[1] = dData(5);
    inp_info[2] = dData(6);

    connectedExternalNodes(0) = idData(20);
    connectedExternalNodes(1) = idData(21);
    connectedExternalNodes(2) = idData(22);
    connectedExternalNodes(3) = idData(23);
    connectedExternalNodes(4) = idData(24);
    connectedExternalNodes(5) = idData(25);

    return res;
}

//**************************************************************************

int
SixNodeBoundryCondition::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    static Vector v5(3);
    static Vector v6(3);

    nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
    nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
    nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
    nodePointers[3]->getDisplayCrds(v4, fact, displayMode);
    nodePointers[4]->getDisplayCrds(v5, fact, displayMode);
    nodePointers[5]->getDisplayCrds(v6, fact, displayMode);

    // color vector
    static Vector values(6);
    values(0) = 0;
    values(1) = 0;
    values(2) = 0;
    values(3) = 0;
    values(4) = 0;
    values(5) = 0;

    // draw polygons for each tetrahedron face -ambaker1
    int res = 0;
    static Matrix coords(6, 3); // rows are face nodes

    // face 1 (1 3 2)
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v4(i);
        coords(2, i) = v2(i);
        coords(3, i) = v5(i);
        coords(4, i) = v3(i);
        coords(5, i) = v6(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    return res;
}

Response*
SixNodeBoundryCondition::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    char outputData[NumDOFsTotal];

    output.tag("ElementOutput");
    output.attr("eleType", "SixNodeBoundryCondition");
    output.attr("eleTag", this->getTag());
    for (int i = 1; i <= 6; i++)
    {
        sprintf(outputData, "node%d", i);
        output.attr(outputData, nodePointers[i - 1]->getTag());
    }

    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0)
    {
        for (int i = 1; i <= 10; i++)
        {
            sprintf(outputData, "P1_%d", i);
            output.tag("ResponseType", outputData);
            sprintf(outputData, "P2_%d", i);
            output.tag("ResponseType", outputData);
            sprintf(outputData, "P3_%d", i);
            output.tag("ResponseType", outputData);
        }
        theResponse = new ElementResponse(this, 1, resid);
    }

    else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
    {
        theResponse = new ElementResponse(this, 2, stiff);
    }

    else if (strcmp(argv[0],"mass") == 0)
    {
        theResponse = new ElementResponse(this, 3, mass);
    }

    else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0)
    {
        int pointNum = atoi(argv[1]);

        if (pointNum > 0 && pointNum <= 10)
        {
            output.tag("GaussPoint");
            output.attr("number", pointNum);
            output.endTag(); // GaussPoint
        }
    }

    else if (strcmp(argv[0], "stresses") == 0)
    {
        for (int i = 0; i < 4; i++)
        {
            output.tag("GaussPoint");
            output.attr("number", i + 1);

            output.tag("ResponseType", "sigma11");
            output.tag("ResponseType", "sigma22");
            output.tag("ResponseType", "sigma33");
            output.tag("ResponseType", "sigma12");
            output.tag("ResponseType", "sigma23");
            output.tag("ResponseType", "sigma13");

            output.endTag(); // GaussPoint
        }
        theResponse =  new ElementResponse(this, 4, Vector(6*10));
    }

    else if (strcmp(argv[0], "strains") == 0)
    {
        for (int i = 0; i < 4; i++)
        {
            output.tag("GaussPoint");
            output.attr("number", i + 1);

            output.tag("ResponseType", "eps11");
            output.tag("ResponseType", "eps22");
            output.tag("ResponseType", "eps33");
            output.tag("ResponseType", "eps12");
            output.tag("ResponseType", "eps23");
            output.tag("ResponseType", "eps13");

            output.endTag(); // GaussPoint
        }
        theResponse =  new ElementResponse(this, 5, Vector(6*10));
    }

    output.endTag(); // ElementOutput

    return theResponse;
}

int
SixNodeBoundryCondition::getResponse(int responseID, Information &eleInfo)
{
    static Vector stresses(3*6);

    if (responseID == 1)
        return eleInfo.setVector(this->getResistingForce());

    else if (responseID == 2)
        return eleInfo.setMatrix(this->getTangentStiff());

    else if (responseID == 3)
        return eleInfo.setMatrix(this->getMass());

    else
        return -1;
}

int
SixNodeBoundryCondition::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    int res = -1;

    return res;
}

int
SixNodeBoundryCondition::updateParameter(int parameterID, Information &info)
{
    int res = -1;
    int matRes = res;

    if (parameterID == res)
    {
        return -1;
    }

    else if (parameterID == 1313)
    {
        int doit = int(info.theDouble);
        if (doit == 1)
        {
            Domain * mydomain = this->getDomain();
            for ( int i = 0; i < NumNodes; i++ )
            {
                nodePointers[i] = mydomain->getNode( connectedExternalNodes(i) ) ;
            }
        }
        return 0;
    }

    else if (parameterID == 1414)
    {
        int new_do_update = int(info.theDouble);
        if (new_do_update == 0)
        {
            opserr << "6NBC::updateParameter - ele tag = " << this->getTag()  << " - will not update\n";
        }
        return 0;
    }

    else
    {
        return res;
    }
}

void
SixNodeBoundryCondition::shp3d( const double zeta[3], double &xsj, double shp[3][NumNodes], const double xl[3][NumNodes]   )
{
    // Mathematica formulation by Carlos Felippa.
    double zeta1 = zeta[0] ; double zeta2 = zeta[1] ; double zeta3 = 1.0 - zeta1 - zeta2 ;

    double x1 = xl[0][0] ; double y1 = xl[1][0] ; double z1 = xl[2][0] ;
    double x2 = xl[0][1] ; double y2 = xl[1][1] ; double z2 = xl[2][1] ;
    double x3 = xl[0][2] ; double y3 = xl[1][2] ; double z3 = xl[2][2] ;
    double x4 = xl[0][3] ; double y4 = xl[1][3] ; double z4 = xl[2][3] ;
    double x5 = xl[0][4] ; double y5 = xl[1][4] ; double z5 = xl[2][4] ;
    double x6 = xl[0][5] ; double y6 = xl[1][5] ; double z6 = xl[2][5] ;

    // dNi/dzeta0
    double dN0_dzeta0  = 4.0 * zeta1 - 1.0 ; double dN1_dzeta0 = 0.0               ;
    double dN2_dzeta0  = 0.0               ; double dN3_dzeta0 = 4.0 * zeta2       ;
    double dN4_dzeta0  = 0.0               ; double dN5_dzeta0 = 4.0 * zeta3       ;

    // dNi/dzeta1
    double dN0_dzeta1  = 0.0               ; double dN1_dzeta1 = 4.0 * zeta2 - 1.0 ;
    double dN2_dzeta1  = 0.0               ; double dN3_dzeta1 = 4.0 * zeta1       ;
    double dN4_dzeta1  = 4.0 * zeta3       ; double dN5_dzeta1 = 0.0               ;

    // dNi/dzeta3
    double dN0_dzeta3  = 0.0               ; double dN1_dzeta3 = 0.0               ;
    double dN2_dzeta3  = 4.0 * zeta3 - 1.0 ; double dN3_dzeta3 = 0.0               ;
    double dN4_dzeta3  = 4.0 * zeta2       ; double dN5_dzeta3 = 4.0 * zeta1       ;

    // IFEM.Ch24 (Eq. 24.7)
    // 
    //                       |  1      1      1  |
    // Jdet = (1 / 2) * det( | Jx1    Jx2    Jx3 | )
    //                       | Jy1    Jy2    Jy3 |
    //

    double dx4 = x4 - (x1 + x2) / 2.0 ;
    double dx5 = x5 - (x2 + x3) / 2.0 ;
    double dx6 = x6 - (x3 + x1) / 2.0 ;

    double dy4 = y4 - (y1 + y2) / 2.0 ;
    double dy5 = y5 - (y2 + y3) / 2.0 ;
    double dy6 = y6 - (y3 + y1) / 2.0 ;

    double Jx21 = (x2 - x1) + 4.0 * (dx4 * (zeta1 - zeta2) + (dx5 - dx6) * zeta3) ;
    double Jx32 = (x3 - x2) + 4.0 * (dx5 * (zeta2 - zeta3) + (dx6 - dx4) * zeta1) ;
    double Jx13 = (x1 - x3) + 4.0 * (dx6 * (zeta3 - zeta1) + (dx4 - dx5) * zeta2) ;

    double Jy12 = (y1 - y2) + 4.0 * (dy4 * (zeta2 - zeta1) + (dy6 - dy5) * zeta3) ;
    double Jy23 = (y2 - y3) + 4.0 * (dy5 * (zeta3 - zeta2) + (dy4 - dy6) * zeta1) ;
    double Jy31 = (y3 - y1) + 4.0 * (dy6 * (zeta1 - zeta3) + (dy5 - dy4) * zeta2) ;

    double Jdet = ( Jx21 * Jy31 - Jy12 * Jx13 ) ;

    // Saving the Jacobians Determinant
    xsj = Jdet / 2.0 ;

    shp[0][0] = (4.0 * zeta1 - 1.0) * Jy23 / Jdet ;
    shp[0][1] = (4.0 * zeta2 - 1.0) * Jy31 / Jdet ;
    shp[0][2] = (4.0 * zeta3 - 1.0) * Jy12 / Jdet ;
    shp[0][3] = 4.0 * (zeta2 * Jy23 + zeta1 * Jy31) / Jdet ;
    shp[0][4] = 4.0 * (zeta3 * Jy31 + zeta2 * Jy12) / Jdet ;
    shp[0][5] = 4.0 * (zeta1 * Jy12 + zeta3 * Jy23) / Jdet ;

    shp[1][0] = (4.0 * zeta1 - 1.0) * Jx32 / Jdet ;
    shp[1][1] = (4.0 * zeta2 - 1.0) * Jx13 / Jdet ;
    shp[1][2] = (4.0 * zeta3 - 1.0) * Jx21 / Jdet ;
    shp[1][3] = 4.0 * (zeta2 * Jx32 + zeta1 * Jx13) / Jdet ;
    shp[1][4] = 4.0 * (zeta3 * Jx13 + zeta2 * Jx21) / Jdet ;
    shp[1][5] = 4.0 * (zeta1 * Jx21 + zeta3 * Jx32) / Jdet ;

    // N1 - N6
    shp[2][0] = zeta1 * (2.0 * zeta1 - 1.0) ;
    shp[2][1] = zeta2 * (2.0 * zeta2 - 1.0) ;
    shp[2][2] = zeta3 * (2.0 * zeta3 - 1.0) ;
    shp[2][3] = 4.0 * zeta1 * zeta2 ;
    shp[2][4] = 4.0 * zeta2 * zeta3 ;
    shp[2][5] = 4.0 * zeta3 * zeta1 ;

    return ;
}

void
SixNodeBoundryCondition::onActivate()
{
    Domain* theDomain = this->getDomain();
    this->setDomain(theDomain);
    this->update();
}

void
SixNodeBoundryCondition::onDeactivate()
{

}