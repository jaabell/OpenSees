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
// Description: interface between a constraint handler that enforces
// multi-point constraints by projecting the explicit acceleration onto the
// constraint manifold and the explicit integrator that performs the march.
// The handler builds the projector and passes it to the integrator through
// this interface (dynamic_cast on the Integrator), so the handler does not
// depend on a concrete integrator class.
//
// Contract for the integrator: store the projector (not owned; the handler
// deletes it), read the diagonal mass before the first factorization after a
// domain change, project the starting acceleration a_0 before the half-step
// velocity is seeded, and project each solved acceleration before it is used
// to form velocities.

#ifndef LadrunoProjectionConsumer_h
#define LadrunoProjectionConsumer_h

class LadrunoConstraintProjector;

class LadrunoProjectionConsumer
{
  public:
    virtual ~LadrunoProjectionConsumer() {}

    // Called by the handler after the DOFs are numbered. A null pointer
    // detaches the projector.
    virtual void setConstraintProjector(LadrunoConstraintProjector *theProjector) = 0;
};

#endif
