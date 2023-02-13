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
// 2023 By Jose Abell and Jose Larenas @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// Implements a 6-node triangle element that accounts for boundry conditions
// for a thermal analysis.
//
// This element has 3 Gauss points of integration.
// ============================================================================


#ifndef SixNodeBoundryCondition_H
#define SixNodeBoundryCondition_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>


//Number of Gauss-points
// #define NumGaussPoints 3
// #define NumNodes 6
// #define NumDOFsPerNode 1
// #define NumStressComponents 3
// #define NumDOFsTotal NumNodes*NumDOFsPerNode

class SixNodeBoundryCondition : public Element {

public :

    //null constructor
    SixNodeBoundryCondition();

    //full constructor
    SixNodeBoundryCondition(int tag,
                       int node1,
                       int node2,
                       int node3,
                       int node4,
                       int node5,
                       int node6,
                       double beta = 0.0,
                       double k    = 0.0,
                       double tamb = 0.0);

    //destructor
    virtual ~SixNodeBoundryCondition( ) ;

    const char *getClassType(void) const {return "SixNodeBoundryCondition";};

    //set domain
    void setDomain( Domain *theDomain ) ;

    //get the number of external nodes
    int getNumExternalNodes( ) const ;

    //return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs(void);

    //return number of dofs
    int getNumDOF( ) ;

    //commit state
    int commitState( ) ;

    //revert to last commit
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    // update
    int update(void);

    //print out element data
    void Print( OPS_Stream &s, int flag ) ;

    //return stiffness matrix
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();
    const Matrix &getDamp();

    void zeroLoad( ) ;
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    //get residual
    const Vector &getResistingForce( ) ;

    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker
                  &theBroker);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);


    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    //plotting
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);
    void onActivate();
    void onDeactivate();

private :

    //Number of Gauss-points
    enum {NumGaussPoints=3} ;
    enum {NumNodes=6} ;
    enum {NumDOFsPerNode=1} ;
    enum {NumStressComponents=3} ;
    enum {NumDOFsTotal=NumNodes*NumDOFsPerNode} ;

    // Routine to compute shape functions and their derivatives. These get stored as follows. 
    // for node n:
    //   shp[0][n] --> dN_n / d x, 
    //   shp[1][n] --> dN_n / d y
    //   shp[2][n] --> dN_n / d z
    //   shp[3][n] --> N_n, shape function n value at the z values 
    void shp3d( 
        const double zeta[3],  // Tetrahedral coordinates  (input)
        double &xsj,         // Jacobian determinant (output)
        double shp[4][NumNodes], // Shape function and derivatives values at the tetrahedral coordinates (output) 
        const double xl[3][NumNodes]   ); // Node coordinates (input)

    //
    // private attributes
    //

    ID connectedExternalNodes ;  //four node numbers
    Node *nodePointers[NumNodes] ;      //pointers to eight nodes

    double inp_info[3];

    Vector *load;
    Matrix *Ki;

    //
    // static attributes
    //

    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ; 
    static Matrix damping ;

    static Matrix B ;

    //quadrature data
    static const double root3 ;
    static const double one_over_root3 ;
    static const double alpha ;
    static const double beta ;
    static const double sg[4] ;
    static const double wg[1] ;

    //local nodal coordinates, three coordinates for each of four nodes
    static double xl[3][NumNodes] ;

    //
    // private methods
    //

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;
    void formDampingTerms( int tangFlag ) ;

    //form residual and tangent
    void formResidAndTangent( int tang_flag ) ;

    //compute coordinate system
    void computeBasis( ) ;

    //compute B matrix
    const Matrix& computeB( int node, const double shp[4][NumNodes] ) ;
} ;

#endif
