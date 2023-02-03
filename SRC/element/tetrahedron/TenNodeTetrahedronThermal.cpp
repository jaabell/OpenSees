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
// 2022 By Jose Abell and José Larenas @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// Please read detailed description in TenNodeTetrahedronThermal.h.
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
#include <TenNodeTetrahedronThermal.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <map>

void* OPS_TenNodeTetrahedronThermal()
{
    if (OPS_GetNumRemainingInputArgs() < 11)
    {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: element TenNodeTetrahedronThermal eleTag? Node1? Node2? Node3? Node4? Node5? Node6? Node7? Node8? Node9? Node10? kxx, kyy, kzz, rho, cp\n";
        return 0;
    }

    int idata[11];
    int num = 11;
    if (OPS_GetIntInput(&num, idata) < 0)
    {
        opserr << "WARNING: invalid integer data\n";
        return 0;
    }

    double data[5] = {0, 0, 0, 0, 0};
    num = OPS_GetNumRemainingInputArgs();

    if (num > 5)
    {
        num = 5;
    }
    if (num > 0)
    {
        if (OPS_GetDoubleInput(&num, data) < 0)
        {
            opserr << "WARNING: invalid double data\n";
            return 0;
        }
    }

    return new TenNodeTetrahedronThermal(idata[0], idata[1], idata[2], idata[3], idata[4], idata[5], idata[6], idata[7], idata[8], idata[9], idata[10], data[0], data[1], data[2], data[3], data[4]);
}

//static data
double  TenNodeTetrahedronThermal::xl[3][NumNodes]                     ;
Matrix  TenNodeTetrahedronThermal::stiff(NumDOFsTotal, NumDOFsTotal)   ;
Vector  TenNodeTetrahedronThermal::resid(NumDOFsTotal)                 ;
Matrix  TenNodeTetrahedronThermal::mass(NumDOFsTotal, NumDOFsTotal)    ;
Matrix  TenNodeTetrahedronThermal::damping(NumDOFsTotal, NumDOFsTotal) ;

//quadrature data
const double  TenNodeTetrahedronThermal::root3 = sqrt( 3.0 ) ;
const double  TenNodeTetrahedronThermal::one_over_root3 = 1.0 / root3 ;

const double  TenNodeTetrahedronThermal::alpha = ( 5.0 + 3.0 * sqrt( 5.0 ) ) / 20. ;
const double  TenNodeTetrahedronThermal::beta  = ( 5.0 - sqrt( 5.0 ) ) / 20.       ;

const double  TenNodeTetrahedronThermal::sg[]  = { alpha, beta, beta, beta } ;
const double  TenNodeTetrahedronThermal::wg[]  = { 1.0 / 4.0 }              ;
// const double  TenNodeTetrahedronThermal::wg[]  = { 1.0 / 24.0 }              ;

static Matrix B(NumStressComponents, NumDOFsPerNode) ;

//null constructor
TenNodeTetrahedronThermal::TenNodeTetrahedronThermal( )
    : Element( 0, ELE_TAG_TenNodeTetrahedronThermal ),
      connectedExternalNodes(NumNodes),
      load(0), Ki(0)
{
    B.Zero();

    krc[0] = 0.0;
    krc[1] = 0.0;
    krc[2] = 0.0;
    krc[3] = 0.0;
    krc[4] = 0.0;

    for (int i = 0; i < NumNodes; i++ ) {
        nodePointers[i] = 0;
    }
}

//*********************************************************************

//full constructor
TenNodeTetrahedronThermal::TenNodeTetrahedronThermal(int tag,
        int node1,
        int node2,
        int node3,
        int node4,
        int node5,
        int node6,
        int node7,
        int node8,
        int node9,
        int node10,
        double kxx,
        double kyy,
        double kzz,
        double rho,
        double cp)
    : Element(tag, ELE_TAG_TenNodeTetrahedronThermal),
      connectedExternalNodes(NumNodes), load(0), Ki(0)
{
    B.Zero();
    connectedExternalNodes(0) = node1 ;
    connectedExternalNodes(1) = node2 ;
    connectedExternalNodes(2) = node3 ;
    connectedExternalNodes(3) = node4 ;
    connectedExternalNodes(4) = node5 ;
    connectedExternalNodes(5) = node6 ;
    connectedExternalNodes(6) = node7 ;
    connectedExternalNodes(7) = node8 ;
    connectedExternalNodes(8) = node9 ;
    connectedExternalNodes(9) = node10 ;

    for (int i = 0; i < NumNodes; i++ ) {
        nodePointers[i] = 0;
    }

    krc[0] = kxx;
    krc[1] = kyy;
    krc[2] = kzz;
    krc[3] = rho;
    krc[4] = cp;
}

//******************************************************************

//destructor
TenNodeTetrahedronThermal::~TenNodeTetrahedronThermal( )
{
    if (load != 0)
        delete load;

    if (Ki != 0)
        delete Ki;
}

//set domain
void  TenNodeTetrahedronThermal::setDomain( Domain *theDomain )
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
int  TenNodeTetrahedronThermal::getNumExternalNodes( ) const
{
    return NumNodes ;
}

//return connected external nodes
const ID&  TenNodeTetrahedronThermal::getExternalNodes( )
{
    return connectedExternalNodes ;
}

Node **
TenNodeTetrahedronThermal::getNodePtrs(void)
{
    return nodePointers ;
}

//return number of dofs
int  TenNodeTetrahedronThermal::getNumDOF( )
{
    return NumDOFsTotal ;
}

//commit state
int  TenNodeTetrahedronThermal::commitState( )
{
    int success = 0 ;

    // call element commitState to do any base class stuff
    if ((success = this->Element::commitState()) != 0) {
        opserr << "TenNodeTetrahedronThermal::commitState () - failed in base class";
    }

    return success ;
}

//revert to last commit
int  TenNodeTetrahedronThermal::revertToLastCommit( )
{
    int success = 0 ;

    return success ;
}

//revert to start
int  TenNodeTetrahedronThermal::revertToStart( )
{
    int success = 0 ;

    return success ;
}

//print out element data
void  TenNodeTetrahedronThermal::Print(OPS_Stream &s, int flag)
{
    if (flag == 2) {
        s << "#TenNodeTetrahedronThermal\n";

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
        s << "Standard TenNodeTetrahedronThermal \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Nodes: " << connectedExternalNodes;
        s << endln;
        s << "Resisting Force (no inertia): " << this->getResistingForce();
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"TenNodeTetrahedronThermal\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        for (int i = 1; i < 2; i++)
            s << connectedExternalNodes(i) << ", ";
        s << connectedExternalNodes(3) << "], ";
    }
}

//return stiffness matrix
const Matrix&  TenNodeTetrahedronThermal::getTangentStiff( )
{
    int tang_flag = 1 ; //get the tangent

    //do tangent and residual here
    formResidAndTangent( tang_flag ) ;

    // opserr << "STIFF = " << stiff << endln ;

    // exit( -1 ) ;

    return stiff ;
}

//return initial matrix
const Matrix&  TenNodeTetrahedronThermal::getInitialStiff( )
{
    if (Ki != 0)
        return *Ki;

    //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31
    static const int ndm = 3 ;
    static const int ndf = NumDOFsPerNode ;
    static const int nstress = NumStressComponents ;
    static const int numberNodes = NumNodes ;
    static const int numberGauss = NumGaussPoints ;
    static const int nShape = 4 ;

    int i, j, k, p, q ;
    int jj, kk ;

    static double volume ;
    static double xsj ;  // determinant jacaobian matrix
    static double dvol[numberGauss] ; //volume element
    static double gaussPoint[ndm] ;
    static Vector strain(nstress) ;  //strain
    static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
    static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions
    static Matrix stiffJK(ndf, ndf) ; //nodeJK stiffness
    static Matrix dd(3, 3) ; //material tangent

    dd.Zero();

    //---------B-matrices------------------------------------

    static Matrix BJ(nstress, ndf) ;     // B matrix node J

    static Matrix BJtran(ndf, nstress) ;

    static Matrix BK(nstress, ndf) ;     // B matrix node k

    static Matrix BJtranD(ndf, nstress) ;

    //-------------------------------------------------------

    //zero stiffness and residual
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
        gaussPoint[2] = sg[abs(2 - k)] ;

        // sg = [alpha, beta, beta, beta]
        // k = 0: alpha beta beta
        // k = 1: beta alpha beta
        // k = 2: beta beta alpha
        // k = 3: beta beta beta

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

        dd(0, 0) = krc[0] * dvol[i];
        dd(1, 1) = krc[1] * dvol[i];
        dd(2, 2) = krc[2] * dvol[i];

        jj = 0;
        for ( j = 0; j < numberNodes; j++ )
        {
            BJ = computeB( j, shp ) ;
            //transpose
            for (p = 0; p < ndf; p++)
            {
                for (q = 0; q < nstress; q++)
                {
                    BJtran(p, q) = BJ(q, p) ;
                }
            }//end for p

            BJtranD.addMatrixProduct(0.0,  BJtran, dd, 1.0) ;

            kk = 0 ;
            for ( k = 0; k < numberNodes; k++ )
            {
                BK = computeB( k, shp ) ;

                stiffJK.addMatrixProduct(0.0,  BJtranD, BK, 1.0) ;

                for ( p = 0; p < ndf; p++ )
                {
                    for ( q = 0; q < ndf; q++ )
                    {
                        stiff( jj + p, kk + q ) += stiffJK( p, q ) ;
                    }
                } //end for p
                kk += ndf ;
            } // end for k loop

            jj += ndf ;
        } // end for j loop
    } //end for i gauss loop

    Ki = new Matrix(stiff);

    return stiff ;
}

//return mass matrix
const Matrix&  TenNodeTetrahedronThermal::getMass( )
{
    int tangFlag = 1 ;

    formInertiaTerms( tangFlag ) ;

    return mass ;
}

//return damping matrix
const Matrix&  TenNodeTetrahedronThermal::getDamp( )
{
    int tangFlag = 1 ;

    formDampingTerms( tangFlag ) ;

    // opserr << "DAMP = " << damping << endln ;

    // exit( -1 ) ;

    return damping ;
}

void  TenNodeTetrahedronThermal::zeroLoad( )
{
    if (load != 0)
        load->Zero();

    return ;
}

int
TenNodeTetrahedronThermal::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    return -1;
}

int
TenNodeTetrahedronThermal::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

//get residual
const Vector&  TenNodeTetrahedronThermal::getResistingForce( )
{
    int tang_flag = 1 ; // get the tangent

    resid.Zero();
    formResidAndTangent( tang_flag ) ;

    if (load != 0)
        resid -= *load;

    return resid ;
}

//get residual with inertia terms
const Vector&  TenNodeTetrahedronThermal::getResistingForceIncInertia( )
{
    static Vector res(10); res.Zero();

    int tang_flag = 1 ; // get the tangent
    resid.Zero();

    //do tangent and residual here
    formResidAndTangent( tang_flag ) ;

    formInertiaTerms( tang_flag ) ;

    formDampingTerms( tang_flag ) ;

    res = resid;

    return res;
}

//*********************************************************************

//form inertia terms
void   TenNodeTetrahedronThermal::formInertiaTerms( int tangFlag )
{
    mass.Zero( ) ;
}

void   TenNodeTetrahedronThermal::formDampingTerms( int tangFlag )
{
    static const int ndm = 3 ;

    static const int ndf = NumDOFsPerNode ;

    static const int numberNodes = NumNodes ;

    static const int numberGauss = NumGaussPoints ;

    static const int nShape = 4 ;

    static const int massIndex = nShape - 1 ;

    static const int dampingIndex = nShape - 1 ;

    double xsj ;  // determinant jacaobian matrix

    double dvol[numberGauss] ; //volume element

    static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

    static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

    static double gaussPoint[ndm] ;

    static Vector momentum(ndf) ;

    int i, j, k, p, q ;
    int jj, kk ;

    double temp, dampingJK ;

    damping.Zero( ) ;

    //compute basis vectors and local nodal coordinates
    computeBasis( ) ;

    //gauss loop to compute and save shape functions
    int count = 0 ;

    for ( k = 0; k < numberGauss; k++ )
    {
        gaussPoint[0] = sg[k] ;
        gaussPoint[1] = sg[abs(1 - k)] ;
        gaussPoint[2] = sg[abs(2 - k)] ;

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
            temp = shp[dampingIndex][j] * dvol[i] ;

            if ( tangFlag == 1 )
            {
                temp *= krc[3] * krc[4] ; // rho * cp

                //node-node mass
                kk = 0 ;
                for ( k = 0; k < numberNodes; k++ )
                {
                    dampingJK = temp * shp[dampingIndex][k] ;
                    for ( p = 0; p < ndf; p++ )
                    {
                        damping( jj + p, kk + p ) += dampingJK ;
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
        const Vector &temp_dot_node_i = nodePointers[i]->getTrialVel( ) ;
        nodeTemp_dot(i) = temp_dot_node_i(0);
    }

    resid.addMatrixVector(1.0, damping, nodeTemp_dot,  1.0);
}

//*********************************************************************

//form residual and tangent
int
TenNodeTetrahedronThermal::update(void)
{
    return 0;
}

//*********************************************************************

//form residual and tangent
void  TenNodeTetrahedronThermal::formResidAndTangent( int tang_flag )
{
    static const int ndm = 3 ;

    static const int ndf = NumDOFsPerNode ;

    static const int nstress = NumStressComponents ;

    static const int numberNodes = NumNodes ;

    static const int numberGauss = NumGaussPoints ;

    static const int nShape = 4 ;

    int i, j, k, p, q ;

    static double volume ;

    static double xsj ;  // determinant jacaobian matrix

    static double dvol[numberGauss] ; //volume element

    static double gaussPoint[ndm] ;

    static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

    static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

    static Vector residJ(ndf) ; //nodeJ residual

    static Matrix stiffJK(ndf, ndf) ; //nodeJK stiffness

    static Vector stress(nstress) ;  //stress

    static Matrix dd(3, 3) ; //material tangent

    //---------B-matrices------------------------------------

    static Matrix BJ(nstress, ndf) ;     // B matrix node J

    static Matrix BJtran(ndf, nstress) ;

    static Matrix BK(nstress, ndf) ;     // B matrix node k

    static Matrix BJtranD(ndf, nstress) ;

    //-------------------------------------------------------

    //zero stiffness and residual
    stiff.Zero( ) ;
    resid.Zero( ) ;

    //compute basis vectors and local nodal coordinates
    computeBasis( ) ;

    //gauss loop to compute and save shape functions
    volume = 0.0 ;
    int count = 0; 

    for ( k = 0; k < numberGauss; k++ )
    {
        gaussPoint[0] = sg[k] ;
        gaussPoint[1] = sg[abs(1 - k)] ;
        gaussPoint[2] = sg[abs(2 - k)] ;

        //get shape functions
        shp3d( gaussPoint, xsj, shp, xl ) ;

        //save shape functions
        for ( p = 0; p < nShape; p++ ) {
            for ( q = 0; q < numberNodes; q++ ) {
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

        if ( tang_flag == 1 ) {
            dd(0, 0) = krc[0] * dvol[i];
            dd(1, 1) = krc[1] * dvol[i];
            dd(2, 2) = krc[2] * dvol[i]; 
        } //end if tang_flag

        //residual and tangent calculations node loops
        int jj = 0 ;
        for ( j = 0; j < numberNodes; j++ )
        {
            BJ = computeB( j, shp ) ;

            //transpose
            for (p = 0; p < ndf; p++)
            {
                for (q = 0; q < nstress; q++)
                {
                    BJtran(p, q) = BJ(q, p) ;
                }
            }//end for p

            if ( tang_flag == 1 )
            {
                BJtranD.addMatrixProduct(0.0,  BJtran, dd, 1.0) ;

                int kk = 0 ;
                for ( k = 0; k < numberNodes; k++ )
                {
                    BK = computeB( k, shp ) ;
                    stiffJK.addMatrixProduct(0.0,  BJtranD, BK, 1.0) ;

                    for ( p = 0; p < ndf; p++ )
                    {
                        for ( q = 0; q < ndf; q++ )
                        {
                            stiff( jj + p, kk + q ) += stiffJK( p, q ) ;
                        }
                    } //end for p
                    kk += ndf ;
                } // end for k loop
            } // end if tang_flag
            jj += ndf ;
        } // end for j loop
    } //end for i gauss loop

    Vector nodeTemp(NumNodes);

    for (int i = 0; i < NumNodes; ++i)
    {
        const Vector &temperature_node_i = nodePointers[i]->getTrialDisp( ) ;
        nodeTemp(i) = temperature_node_i(0);
    }

    resid.addMatrixVector(0.0, stiff, nodeTemp, 1.0);

    return ;
}

//************************************************************************

//compute local coordinates and basis
void   TenNodeTetrahedronThermal::computeBasis( )
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
TenNodeTetrahedronThermal::computeB( int node, const double shp[4][NumNodes] )
{
    //--------------------------------------------------------------------
    //
    //                -    -
    //               | N,x  |
    //   B       =   | N,y  |
    //               | N,z  |   (3x1)
    //                -    -
    //
    //-------------------------------------------------------------------

    B(0, 0) = shp[0][node] ;
    B(1, 0) = shp[1][node] ;
    B(2, 0) = shp[2][node] ;

    return B ;
}

//**********************************************************************

int  TenNodeTetrahedronThermal::sendSelf (int commitTag, Channel &theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // Quad packs its data into a Vector and sends this to theChannel
    // along with its dbTag and the commitTag passed in the arguments

    // Now quad sends the ids of its materials
    static ID idData(30);

    idData(20) = connectedExternalNodes(0);
    idData(21) = connectedExternalNodes(1);
    idData(22) = connectedExternalNodes(2);
    idData(23) = connectedExternalNodes(3);
    idData(24) = connectedExternalNodes(4);
    idData(25) = connectedExternalNodes(5);
    idData(26) = connectedExternalNodes(6);
    idData(27) = connectedExternalNodes(7);
    idData(28) = connectedExternalNodes(8);
    idData(29) = connectedExternalNodes(9);

    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING TenNodeTetrahedronThermal::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    static Vector dData(9);
    dData(0) = alphaM;
    dData(1) = betaK;
    dData(2) = betaK0;
    dData(3) = betaKc;
    dData(4) = krc[0];
    dData(5) = krc[1];
    dData(6) = krc[2];
    dData(7) = krc[3];
    dData(8) = krc[4];

    if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
        opserr << "TenNodeTetrahedronThermal::sendSelf() - failed to send double data\n";
        return -1;
    }

    return res;
}

int  TenNodeTetrahedronThermal::recvSelf (int commitTag,
                                   Channel &theChannel,
                                   FEM_ObjectBroker &theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    static ID idData(30);
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING TenNodeTetrahedronThermal::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    this->setTag(idData(30));

    static Vector dData(9);
    if (theChannel.recvVector(dataTag, commitTag, dData) < 0) {
        opserr << "DispBeamColumn2d::sendSelf() - failed to recv double data\n";
        return -1;
    }

    alphaM = dData(0);
    betaK  = dData(1);
    betaK0 = dData(2);
    betaKc = dData(3);
    krc[0] = dData(4);
    krc[1] = dData(5);
    krc[2] = dData(6);
    krc[3] = dData(7);
    krc[4] = dData(8);

    connectedExternalNodes(0) = idData(20);
    connectedExternalNodes(1) = idData(21);
    connectedExternalNodes(2) = idData(22);
    connectedExternalNodes(3) = idData(23);
    connectedExternalNodes(4) = idData(24);
    connectedExternalNodes(5) = idData(25);
    connectedExternalNodes(6) = idData(26);
    connectedExternalNodes(7) = idData(27);
    connectedExternalNodes(8) = idData(28);
    connectedExternalNodes(9) = idData(29);

    return res;
}

//**************************************************************************

int
TenNodeTetrahedronThermal::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    static Vector v5(3);
    static Vector v6(3);
    static Vector v7(3);
    static Vector v8(3);
    static Vector v9(3);
    static Vector v10(3);

    nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
    nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
    nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
    nodePointers[3]->getDisplayCrds(v4, fact, displayMode);
    nodePointers[4]->getDisplayCrds(v5, fact, displayMode);
    nodePointers[5]->getDisplayCrds(v6, fact, displayMode);
    nodePointers[6]->getDisplayCrds(v7, fact, displayMode);
    nodePointers[7]->getDisplayCrds(v8, fact, displayMode);
    nodePointers[8]->getDisplayCrds(v9, fact, displayMode);
    nodePointers[9]->getDisplayCrds(v10, fact, displayMode);

    // color vector
    static Vector values(10);
    values(0) = 0;
    values(1) = 0;
    values(2) = 0;
    values(3) = 0;
    values(4) = 0;
    values(5) = 0;
    values(6) = 0;
    values(7) = 0;
    values(8) = 0;
    values(9) = 0;

    // draw polygons for each tetrahedron face -ambaker1
    int res = 0;
    static Matrix coords(10, 3); // rows are face nodes

    // face 1 (1 3 2)
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v5(i);
        coords(2, i) = v2(i);
        coords(3, i) = v6(i);
        coords(4, i) = v3(i);
        coords(5, i) = v7(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    // face 2 (1 2 4)
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v5(i);
        coords(2, i) = v2(i);
        coords(3, i) = v10(i);
        coords(4, i) = v4(i);
        coords(5, i) = v8(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    // face 3 (1 4 3)
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v7(i);
        coords(2, i) = v3(i);
        coords(3, i) = v9(i);
        coords(4, i) = v4(i);
        coords(5, i) = v8(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    // face 4 (2 3 4)
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v2(i);
        coords(1, i) = v6(i);
        coords(2, i) = v3(i);
        coords(3, i) = v9(i);
        coords(4, i) = v4(i);
        coords(5, i) = v10(i);
    }
    res += theViewer.drawPolygon(coords, values, this->getTag());

    return res;
}

Response*
TenNodeTetrahedronThermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    char outputData[NumDOFsTotal];

    output.tag("ElementOutput");
    output.attr("eleType", "TenNodeTetrahedronThermal");
    output.attr("eleTag", this->getTag());
    for (int i = 1; i <= 10; i++)
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

    else if (strcmp(argv[0],"damp") == 0)
    {
        theResponse = new ElementResponse(this, 4, damping);
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
        theResponse =  new ElementResponse(this, 5, Vector(6*10));
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
TenNodeTetrahedronThermal::getResponse(int responseID, Information &eleInfo)
{
    static Vector stresses(3*10);

    if (responseID == 1)
        return eleInfo.setVector(this->getResistingForce());

    else if (responseID == 2)
        return eleInfo.setMatrix(this->getTangentStiff());

    else if (responseID == 3)
        return eleInfo.setMatrix(this->getMass());

    else if (responseID == 4)
        return eleInfo.setMatrix(this->getDamp());

    else
        return -1;
}

int
TenNodeTetrahedronThermal::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    int res = -1;

    return res;
}

int
TenNodeTetrahedronThermal::updateParameter(int parameterID, Information &info)
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
            opserr << "10Ntet::updateParameter - ele tag = " << this->getTag()  << " - will not update\n";
        }
        return 0;
    }

    else
    {
        return res;
    }
}

void
TenNodeTetrahedronThermal::shp3d( const double zeta[4], double &xsj, double shp[4][NumNodes], const double xl[3][NumNodes]   )
{
    // Mathematica formulation by Carlos Felippa.
    // Modified by José Larenas so that it works when
    // switching N9 with N10.

    double zeta4 = 1 - zeta[0] - zeta[1] - zeta[2] ;

    double x1  = xl[0][0]; double y1  = xl[1][0] ; double z1  = xl[2][0] ;
    double x2  = xl[0][1]; double y2  = xl[1][1] ; double z2  = xl[2][1] ;
    double x3  = xl[0][2]; double y3  = xl[1][2] ; double z3  = xl[2][2] ;
    double x4  = xl[0][3]; double y4  = xl[1][3] ; double z4  = xl[2][3] ;
    double x5  = xl[0][4]; double y5  = xl[1][4] ; double z5  = xl[2][4] ;
    double x6  = xl[0][5]; double y6  = xl[1][5] ; double z6  = xl[2][5] ;
    double x7  = xl[0][6]; double y7  = xl[1][6] ; double z7  = xl[2][6] ;
    double x8  = xl[0][7]; double y8  = xl[1][7] ; double z8  = xl[2][7] ;
    double x9  = xl[0][8]; double y9  = xl[1][8] ; double z9  = xl[2][8] ;
    double x10 = xl[0][9]; double y10 = xl[1][9] ; double z10 = xl[2][9] ;

    double a1 = y2 * (z4 - z3) - y3 * (z4 - z2) + y4 * (z3 - z2)  ; double b1 = -x2 * (z4 - z3) + x3 * (z4 - z2) - x4 * (z3 - z2) ; double c1 = x2 * (y4 - y3) - x3 * (y4 - y2) + x4 * (y3 - y2)  ;
    double a2 = -y1 * (z4 - z3) + y3 * (z4 - z1) - y4 * (z3 - z1) ; double b2 = x1 * (z4 - z3) - x3 * (z4 - z1) + x4 * (z3 - z1)  ; double c2 = -x1 * (y4 - y3) + x3 * (y4 - y1) - x4 * (y3 - y1) ;
    double a3 = y1 * (z4 - z2) - y2 * (z4 - z1) + y4 * (z2 - z1)  ; double b3 = -x1 * (z4 - z2) + x2 * (z4 - z1) - x4 * (z2 - z1) ; double c3 = x1 * (y4 - y2) - x2 * (y4 - y1) + x4 * (y2 - y1)  ;
    double a4 = -y1 * (z3 - z2) + y2 * (z3 - z1) - y3 * (z2 - z1) ; double b4 = x1 * (z3 - z2) - x2 * (z3 - z1) + x3 * (z2 - z1)  ; double c4 = -x1 * (y3 - y2) + x2 * (y3 - y1) - x3 * (y2 - y1) ;

    // dNi/dzeta1
    double dN1_dzeta1  = 4 * zeta[0] - 1 ; double dN2_dzeta1 = 0         ;
    double dN3_dzeta1  = 0               ; double dN4_dzeta1 = 0         ;
    double dN5_dzeta1  = 4 * zeta[1]     ; double dN6_dzeta1 = 0         ;
    double dN7_dzeta1  = 4 * zeta[2]     ; double dN8_dzeta1 = 4 * zeta4 ;
    double dN10_dzeta1 = 0               ; double dN9_dzeta1 = 0         ;

    // dNi/dzeta2
    double dN1_dzeta2  = 0           ; double dN2_dzeta2 = 4 * zeta[1] - 1 ;
    double dN3_dzeta2  = 0           ; double dN4_dzeta2 = 0               ;
    double dN5_dzeta2  = 4 * zeta[0] ; double dN6_dzeta2 = 4 * zeta[2]     ;
    double dN7_dzeta2  = 0           ; double dN8_dzeta2 = 0               ;
    double dN10_dzeta2 = 4 * zeta4   ; double dN9_dzeta2 = 0               ;

    // dNi/dzeta3
    double dN1_dzeta3  = 0               ; double dN2_dzeta3 = 0           ;
    double dN3_dzeta3  = 4 * zeta[2] - 1 ; double dN4_dzeta3 = 0           ;
    double dN5_dzeta3  = 0               ; double dN6_dzeta3 = 4 * zeta[1] ;
    double dN7_dzeta3  = 4 * zeta[0]     ; double dN8_dzeta3 = 0           ;
    double dN10_dzeta3 = 0               ; double dN9_dzeta3 = 4 * zeta4   ;

    // dNi/dzeta3
    double dN1_dzeta4  = 0           ; double dN2_dzeta4 = 0             ;
    double dN3_dzeta4  = 0           ; double dN4_dzeta4 = 4 * zeta4 - 1 ;
    double dN5_dzeta4  = 0           ; double dN6_dzeta4 = 0             ;
    double dN7_dzeta4  = 0           ; double dN8_dzeta4 = 4 * zeta[0]   ;
    double dN10_dzeta4 = 4 * zeta[1] ; double dN9_dzeta4 = 4 * zeta[2]   ;

    double Jx1 = x1 * dN1_dzeta1 + x2 * dN2_dzeta1 + x3 * dN3_dzeta1 + x4 * dN4_dzeta1 + x5 * dN5_dzeta1 + x6 * dN6_dzeta1 + x7 * dN7_dzeta1 + x8 * dN8_dzeta1 + x9 * dN9_dzeta1 + x10 * dN10_dzeta1 ;
    double Jy1 = y1 * dN1_dzeta1 + y2 * dN2_dzeta1 + y3 * dN3_dzeta1 + y4 * dN4_dzeta1 + y5 * dN5_dzeta1 + y6 * dN6_dzeta1 + y7 * dN7_dzeta1 + y8 * dN8_dzeta1 + y9 * dN9_dzeta1 + y10 * dN10_dzeta1 ;
    double Jz1 = z1 * dN1_dzeta1 + z2 * dN2_dzeta1 + z3 * dN3_dzeta1 + z4 * dN4_dzeta1 + z5 * dN5_dzeta1 + z6 * dN6_dzeta1 + z7 * dN7_dzeta1 + z8 * dN8_dzeta1 + z9 * dN9_dzeta1 + z10 * dN10_dzeta1 ;

    double Jx2 = x1 * dN1_dzeta2 + x2 * dN2_dzeta2 + x3 * dN3_dzeta2 + x4 * dN4_dzeta2 + x5 * dN5_dzeta2 + x6 * dN6_dzeta2 + x7 * dN7_dzeta2 + x8 * dN8_dzeta2 + x9 * dN9_dzeta2 + x10 * dN10_dzeta2 ;
    double Jy2 = y1 * dN1_dzeta2 + y2 * dN2_dzeta2 + y3 * dN3_dzeta2 + y4 * dN4_dzeta2 + y5 * dN5_dzeta2 + y6 * dN6_dzeta2 + y7 * dN7_dzeta2 + y8 * dN8_dzeta2 + y9 * dN9_dzeta2 + y10 * dN10_dzeta2 ;
    double Jz2 = z1 * dN1_dzeta2 + z2 * dN2_dzeta2 + z3 * dN3_dzeta2 + z4 * dN4_dzeta2 + z5 * dN5_dzeta2 + z6 * dN6_dzeta2 + z7 * dN7_dzeta2 + z8 * dN8_dzeta2 + z9 * dN9_dzeta2 + z10 * dN10_dzeta2 ;

    double Jx3 = x1 * dN1_dzeta3 + x2 * dN2_dzeta3 + x3 * dN3_dzeta3 + x4 * dN4_dzeta3 + x5 * dN5_dzeta3 + x6 * dN6_dzeta3 + x7 * dN7_dzeta3 + x8 * dN8_dzeta3 + x9 * dN9_dzeta3 + x10 * dN10_dzeta3 ;
    double Jy3 = y1 * dN1_dzeta3 + y2 * dN2_dzeta3 + y3 * dN3_dzeta3 + y4 * dN4_dzeta3 + y5 * dN5_dzeta3 + y6 * dN6_dzeta3 + y7 * dN7_dzeta3 + y8 * dN8_dzeta3 + y9 * dN9_dzeta3 + y10 * dN10_dzeta3 ;
    double Jz3 = z1 * dN1_dzeta3 + z2 * dN2_dzeta3 + z3 * dN3_dzeta3 + z4 * dN4_dzeta3 + z5 * dN5_dzeta3 + z6 * dN6_dzeta3 + z7 * dN7_dzeta3 + z8 * dN8_dzeta3 + z9 * dN9_dzeta3 + z10 * dN10_dzeta3 ;

    double Jx4 = x1 * dN1_dzeta4 + x2 * dN2_dzeta4 + x3 * dN3_dzeta4 + x4 * dN4_dzeta4 + x5 * dN5_dzeta4 + x6 * dN6_dzeta4 + x7 * dN7_dzeta4 + x8 * dN8_dzeta4 + x9 * dN9_dzeta4 + x10 * dN10_dzeta4 ;
    double Jy4 = y1 * dN1_dzeta4 + y2 * dN2_dzeta4 + y3 * dN3_dzeta4 + y4 * dN4_dzeta4 + y5 * dN5_dzeta4 + y6 * dN6_dzeta4 + y7 * dN7_dzeta4 + y8 * dN8_dzeta4 + y9 * dN9_dzeta4 + y10 * dN10_dzeta4 ;
    double Jz4 = z1 * dN1_dzeta4 + z2 * dN2_dzeta4 + z3 * dN3_dzeta4 + z4 * dN4_dzeta4 + z5 * dN5_dzeta4 + z6 * dN6_dzeta4 + z7 * dN7_dzeta4 + z8 * dN8_dzeta4 + z9 * dN9_dzeta4 + z10 * dN10_dzeta4 ;

    // Terms in simplified Jacobian Matrix (3x3)
    double t1 = Jx2 - Jx1 ; double t2 = Jx3 - Jx1 ; double t3 = Jx4 - Jx1 ;
    double t4 = Jy2 - Jy1 ; double t5 = Jy3 - Jy1 ; double t6 = Jy4 - Jy1 ;
    double t7 = Jz2 - Jz1 ; double t8 = Jz3 - Jz1 ; double t9 = Jz4 - Jz1 ;

    // Assembling the Jacobians Determinant
    double Jdet = t1 * (t5 * t9 - t6 * t8) - t2 * (t4 * t9 - t6 * t7) + t3 * (t4 * t8 - t5 * t7) ;

    // opserr << "DET J" << Jdet << endln ;

    // Saving the Jacobians Determinant
    xsj = Jdet ;
    // xsj = Jdet / 6.0 ;

    // qx1 - qx10 (17.24)
    shp[0][0] = ( 1 / Jdet ) * (dN1_dzeta1 * a1  + dN1_dzeta2 * a2  + dN1_dzeta3 * a3  + dN1_dzeta4 * a4)  ;
    shp[0][1] = ( 1 / Jdet ) * (dN2_dzeta1 * a1  + dN2_dzeta2 * a2  + dN2_dzeta3 * a3  + dN2_dzeta4 * a4)  ;
    shp[0][2] = ( 1 / Jdet ) * (dN3_dzeta1 * a1  + dN3_dzeta2 * a2  + dN3_dzeta3 * a3  + dN3_dzeta4 * a4)  ;
    shp[0][3] = ( 1 / Jdet ) * (dN4_dzeta1 * a1  + dN4_dzeta2 * a2  + dN4_dzeta3 * a3  + dN4_dzeta4 * a4)  ;
    shp[0][4] = ( 1 / Jdet ) * (dN5_dzeta1 * a1  + dN5_dzeta2 * a2  + dN5_dzeta3 * a3  + dN5_dzeta4 * a4)  ;
    shp[0][5] = ( 1 / Jdet ) * (dN6_dzeta1 * a1  + dN6_dzeta2 * a2  + dN6_dzeta3 * a3  + dN6_dzeta4 * a4)  ;
    shp[0][6] = ( 1 / Jdet ) * (dN7_dzeta1 * a1  + dN7_dzeta2 * a2  + dN7_dzeta3 * a3  + dN7_dzeta4 * a4)  ;
    shp[0][7] = ( 1 / Jdet ) * (dN8_dzeta1 * a1  + dN8_dzeta2 * a2  + dN8_dzeta3 * a3  + dN8_dzeta4 * a4)  ;
    shp[0][8] = ( 1 / Jdet ) * (dN9_dzeta1 * a1  + dN9_dzeta2 * a2  + dN9_dzeta3 * a3  + dN9_dzeta4 * a4)  ;
    shp[0][9] = ( 1 / Jdet ) * (dN10_dzeta1 * a1 + dN10_dzeta2 * a2 + dN10_dzeta3 * a3 + dN10_dzeta4 * a4) ;

    // qy1 - qy10 (17.24)
    shp[1][0] = ( 1 / Jdet ) * (dN1_dzeta1 * b1  + dN1_dzeta2 * b2  + dN1_dzeta3 * b3  + dN1_dzeta4 * b4)  ;
    shp[1][1] = ( 1 / Jdet ) * (dN2_dzeta1 * b1  + dN2_dzeta2 * b2  + dN2_dzeta3 * b3  + dN2_dzeta4 * b4)  ;
    shp[1][2] = ( 1 / Jdet ) * (dN3_dzeta1 * b1  + dN3_dzeta2 * b2  + dN3_dzeta3 * b3  + dN3_dzeta4 * b4)  ;
    shp[1][3] = ( 1 / Jdet ) * (dN4_dzeta1 * b1  + dN4_dzeta2 * b2  + dN4_dzeta3 * b3  + dN4_dzeta4 * b4)  ;
    shp[1][4] = ( 1 / Jdet ) * (dN5_dzeta1 * b1  + dN5_dzeta2 * b2  + dN5_dzeta3 * b3  + dN5_dzeta4 * b4)  ;
    shp[1][5] = ( 1 / Jdet ) * (dN6_dzeta1 * b1  + dN6_dzeta2 * b2  + dN6_dzeta3 * b3  + dN6_dzeta4 * b4)  ;
    shp[1][6] = ( 1 / Jdet ) * (dN7_dzeta1 * b1  + dN7_dzeta2 * b2  + dN7_dzeta3 * b3  + dN7_dzeta4 * b4)  ;
    shp[1][7] = ( 1 / Jdet ) * (dN8_dzeta1 * b1  + dN8_dzeta2 * b2  + dN8_dzeta3 * b3  + dN8_dzeta4 * b4)  ;
    shp[1][8] = ( 1 / Jdet ) * (dN9_dzeta1 * b1  + dN9_dzeta2 * b2  + dN9_dzeta3 * b3  + dN9_dzeta4 * b4)  ;
    shp[1][9] = ( 1 / Jdet ) * (dN10_dzeta1 * b1 + dN10_dzeta2 * b2 + dN10_dzeta3 * b3 + dN10_dzeta4 * b4) ;

    // qz1 - qz10 (17.24)
    shp[2][0] = ( 1 / Jdet ) * (dN1_dzeta1 * c1  + dN1_dzeta2 * c2  + dN1_dzeta3 * c3  + dN1_dzeta4 * c4)  ;
    shp[2][1] = ( 1 / Jdet ) * (dN2_dzeta1 * c1  + dN2_dzeta2 * c2  + dN2_dzeta3 * c3  + dN2_dzeta4 * c4)  ;
    shp[2][2] = ( 1 / Jdet ) * (dN3_dzeta1 * c1  + dN3_dzeta2 * c2  + dN3_dzeta3 * c3  + dN3_dzeta4 * c4)  ;
    shp[2][3] = ( 1 / Jdet ) * (dN4_dzeta1 * c1  + dN4_dzeta2 * c2  + dN4_dzeta3 * c3  + dN4_dzeta4 * c4)  ;
    shp[2][4] = ( 1 / Jdet ) * (dN5_dzeta1 * c1  + dN5_dzeta2 * c2  + dN5_dzeta3 * c3  + dN5_dzeta4 * c4)  ;
    shp[2][5] = ( 1 / Jdet ) * (dN6_dzeta1 * c1  + dN6_dzeta2 * c2  + dN6_dzeta3 * c3  + dN6_dzeta4 * c4)  ;
    shp[2][6] = ( 1 / Jdet ) * (dN7_dzeta1 * c1  + dN7_dzeta2 * c2  + dN7_dzeta3 * c3  + dN7_dzeta4 * c4)  ;
    shp[2][7] = ( 1 / Jdet ) * (dN8_dzeta1 * c1  + dN8_dzeta2 * c2  + dN8_dzeta3 * c3  + dN8_dzeta4 * c4)  ;
    shp[2][8] = ( 1 / Jdet ) * (dN9_dzeta1 * c1  + dN9_dzeta2 * c2  + dN9_dzeta3 * c3  + dN9_dzeta4 * c4)  ;
    shp[2][9] = ( 1 / Jdet ) * (dN10_dzeta1 * c1 + dN10_dzeta2 * c2 + dN10_dzeta3 * c3 + dN10_dzeta4 * c4) ;

    // N1 - N10 (notice N9 and N10 are switched up)
    shp[3][0] = zeta[0] * (2 * zeta[0] - 1) ;
    shp[3][1] = zeta[1] * (2 * zeta[1] - 1) ;
    shp[3][2] = zeta[2] * (2 * zeta[2] - 1) ;
    shp[3][3] = zeta4 * (2 * zeta4 - 1)     ;
    shp[3][4] = 4 * zeta[0] * zeta[1]       ;
    shp[3][5] = 4 * zeta[1] * zeta[2]       ;
    shp[3][6] = 4 * zeta[2] * zeta[0]       ;
    shp[3][7] = 4 * zeta[0] * zeta4         ;
    shp[3][9] = 4 * zeta[1] * zeta4         ; 
    shp[3][8] = 4 * zeta[2] * zeta4         ; 

    return ;
}

void
TenNodeTetrahedronThermal::onActivate()
{
    Domain* theDomain = this->getDomain();
    this->setDomain(theDomain);
    this->update();
}

void
TenNodeTetrahedronThermal::onDeactivate()
{

}