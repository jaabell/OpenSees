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
                                                                       
#ifndef UnsaturatedSSPquadUP_h
#define UnsaturatedSSPquadUP_h

//
// Adapted: Luis Miranda, LNEC, Andre Barbosa, OSU; July 2015 - added boundary tractions (see LM changes)
//
// Created: C.McGann, UW, 05.2011
//
// Description: This file contains the class definition for UnsaturatedSSPquadUP
//                Stabilized Single-Point Quad element with a u-p formulation 
//                for plane strain analysis of saturated porous media
//
// Reference:   Zienkiewicz, O.C. and Shiomi, T. (1984). "Dynamic behavior of 
//                saturated porous media; the generalized Biot formulation and 
//                its numerical solution." International Journal for Numerical 
//                Methods in Geomechanics, 8, 71-96.

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <SSPquadUP.h>

// number of nodes per element
#define SQUP_NUM_NODE 4
// number of dimensions
#define SQUP_NUM_DIM  2
// degrees of freedom per element
#define SQUP_NUM_DOF  12

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;
class Response;

class UnsaturatedSSPquadUP : public SSPquadUP
{
  public:
    // LM change
    UnsaturatedSSPquadUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterial &theMat,
                       double thick, double Kf, double Rf, double k1, double k2,
                       double eVoid, double alpha, double b1 = 0.0, double b2 = 0.0,
                       double Pup = 0.0, double Plow = 0.0, double Pleft = 0.0, double Pright = 0.0);
    UnsaturatedSSPquadUP();

    const char* getClassType()  const { return "UnsaturatedSSPquadUP"; };
    
    const Matrix &getMass(void);
    const Vector &getResistingForce(void); 
 

  protected:

    void GetPermeabilityMatrix(void);          // compute permeability matrix ****

  private:

   
};

#endif
 
