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
// References:
//   Wriggers, P. (2006). Computational Contact Mechanics, 2nd ed. Springer,
//     Berlin. doi:10.1007/978-3-540-32609-0
//   Laursen, T.A. (2002). Computational Contact and Impact Mechanics.
//     Springer, Berlin. doi:10.1007/978-3-662-04864-1
//   Hallquist, J.O. (2006). LS-DYNA Theory Manual, Livermore Software
//     Technology Corporation, Section 26 (Contact-Impact Algorithm).
//
// Description: LadrunoContactHandler, a ConstraintHandler that adds the
// contact FE_Element adapters (LadrunoContactFE) to the AnalysisModel.
// handle() is where FE_Elements are created (inside domainChanged(), after the
// AnalysisModel is cleared), so contact enters the assembly through a
// constraint handler. handle() builds the DOF_Groups and FE_Elements the way
// PlainHandler does, then appends the contact adapters with tags above the
// structural FE_Element tags (AnalysisModel::addFE_Element drops a duplicate
// tag).
//
// Like PlainHandler, single-point constraints are enforced by removing the
// DOF, so only homogeneous constraints are enforced; a non-homogeneous one,
// an MP_Constraint or a non-trivial EQ_Constraint is reported and not
// enforced.
//
// Command: constraints LadrunoContact

#ifndef LadrunoContactHandler_h
#define LadrunoContactHandler_h

#include <ConstraintHandler.h>

class FE_Element;
class DOF_Group;

class LadrunoContactHandler : public ConstraintHandler
{
  public:
    LadrunoContactHandler();
    ~LadrunoContactHandler();

    int handle(const ID *nodesNumberedLast = 0);
    void clearAll(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
};

#endif
