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
// Description: reduced-operator projection for complex modal analysis.
// See LadrunoDampingAssembler.h for the formulation and its contract.

#include "LadrunoDampingAssembler.h"

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>
#include <Matrix.h>
#include <Vector.h>
#include <OPS_Globals.h>

namespace {

// Gather the p retained mode columns of one node's eigenvector storage into
// rows [rowStart, rowStart+ndf) of PhiE. A node with NO stored eigenvectors
// (fully fixed: it never entered the eigen analysis) contributes zero rows,
// which is exactly its physical mode shape. Returns false when the node HAS
// eigenvectors but fewer than p (a stale or short eigen run; the caller errors out).
bool
gatherNodeRows(Node *nodePtr, int p, Matrix &PhiE, int rowStart)
{
    const int ndf = nodePtr->getNumberDOF();
    const int stored = nodePtr->getNumEigenvectors();
    if (stored == 0) {
        for (int i = 0; i < ndf; i++)
            for (int a = 0; a < p; a++)
                PhiE(rowStart + i, a) = 0.0;
        return true;
    }
    if (stored < p)
        return false;
    const Matrix &vecs = nodePtr->getEigenvectors();  // ndf x numEigen
    for (int i = 0; i < ndf; i++)
        for (int a = 0; a < p; a++)
            PhiE(rowStart + i, a) = vecs(i, a);
    return true;
}

} // namespace

int
LadrunoDampingAssembler::projectReduced(Domain &theDomain, int p,
                                        Matrix &Mt, Matrix &Ct, Matrix &Kt)
{
    if (p < 1 || Mt.noRows() != p || Mt.noCols() != p ||
        Ct.noRows() != p || Ct.noCols() != p ||
        Kt.noRows() != p || Kt.noCols() != p) {
        opserr << "LadrunoDampingAssembler::projectReduced - bad dimensions\n";
        return -1;
    }
    const int pAvail = theDomain.getNumEigenvalues();
    if (pAvail < p) {
        opserr << "LadrunoDampingAssembler::projectReduced - " << p
               << " modes requested but the domain holds " << pAvail << "\n";
        return -2;
    }
    const Vector &eigenvalues = theDomain.getEigenvalues();

    Mt.Zero();  Ct.Zero();  Kt.Zero();

    // ---- element contributions: Phi_e^T (C_e | M_e) Phi_e ----
    Element *elePtr;
    ElementIter &theElemIter = theDomain.getElements();
    while ((elePtr = theElemIter()) != 0) {
        const ID &extNodes = elePtr->getExternalNodes();
        Node **nodePtrs = elePtr->getNodePtrs();
        const int numNodes = extNodes.Size();

        int nde = 0;
        for (int n = 0; n < numNodes; n++) {
            if (nodePtrs[n] == 0) {
                opserr << "LadrunoDampingAssembler - element "
                       << elePtr->getTag() << " has a null node pointer\n";
                return -3;
            }
            nde += nodePtrs[n]->getNumberDOF();
        }

        // getDamp()/getMass() habitually return the SAME per-class static
        // scratch matrix: never hold references to both. Check dims with
        // separate calls, deep-copy C, then consume M straight away.
        if (elePtr->getDamp().noRows() != nde ||
            elePtr->getMass().noRows() != nde) {
            // an element whose matrix dimension is not the sum of its nodes'
            // DOFs (internal/condensed DOFs) cannot be projected from nodal
            // eigenvectors: refuse rather than misproject
            opserr << "LadrunoDampingAssembler - element " << elePtr->getTag()
                   << ": damping/mass dimension != sum of nodal DOFs (" << nde
                   << "); cannot project this model from nodal eigenvectors\n";
            return -4;
        }

        Matrix PhiE(nde, p);
        int row = 0;
        for (int n = 0; n < numNodes; n++) {
            if (!gatherNodeRows(nodePtrs[n], p, PhiE, row)) {
                opserr << "LadrunoDampingAssembler - node "
                       << nodePtrs[n]->getTag() << " stores fewer than " << p
                       << " eigenvectors; re-run 'eigen'\n";
                return -5;
            }
            row += nodePtrs[n]->getNumberDOF();
        }

        const Matrix Ce(elePtr->getDamp());              // deep copy
        Mt.addMatrixTripleProduct(1.0, PhiE, elePtr->getMass(), 1.0);
        Ct.addMatrixTripleProduct(1.0, PhiE, Ce, 1.0);
    }

    // ---- nodal contributions: lumped nodal mass + its alphaM Rayleigh ----
    Node *nodePtr;
    NodeIter &theNodeIter = theDomain.getNodes();
    while ((nodePtr = theNodeIter()) != 0) {
        if (nodePtr->getNumEigenvectors() == 0)
            continue;  // fixed/outside the analysis: phi = 0 rows
        const int ndf = nodePtr->getNumberDOF();
        Matrix PhiN(ndf, p);
        if (!gatherNodeRows(nodePtr, p, PhiN, 0)) {
            opserr << "LadrunoDampingAssembler - node " << nodePtr->getTag()
                   << " stores fewer than " << p
                   << " eigenvectors; re-run 'eigen'\n";
            return -5;
        }
        // same static-scratch aliasing hazard as elements: copy C first
        const Matrix Cn(nodePtr->getDamp());       // alphaM_node * Mn (copy)
        Mt.addMatrixTripleProduct(1.0, PhiN, nodePtr->getMass(), 1.0);
        Ct.addMatrixTripleProduct(1.0, PhiN, Cn, 1.0);
    }

    // ---- Kt from the spectrum: Kt = sym(Mt * diag(w^2)), exact, see .h ----
    for (int a = 0; a < p; a++)
        for (int b = 0; b < p; b++)
            Kt(a, b) = 0.5 * Mt(a, b) * (eigenvalues(a) + eigenvalues(b));

    return 0;
}
