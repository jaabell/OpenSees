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
// Description: LadrunoContactSurface, see LadrunoContactSurface.h.

#include "LadrunoContactSurface.h"
#include <Domain.h>
#include <Node.h>

LadrunoContactSurface::LadrunoContactSurface(int tag, Kind kind,
                                             const ID &nodeTags, int npseg)
  : myTag(tag), myKind(kind), theNodeTags(nodeTags), nodesPerSeg(npseg)
{
}

LadrunoContactSurface::~LadrunoContactSurface()
{
}

int
LadrunoContactSurface::getNumSegments(void) const
{
    if (myKind == MASTER_SEGMENTS && nodesPerSeg > 0)
        return theNodeTags.Size() / nodesPerSeg;   // faceted: count facets
    return theNodeTags.Size();   // SLAVE_NODES: one entry per slave node
}

int
LadrunoContactSurface::pruneMissingNodes(Domain *dom)
{
    if (dom == 0)
        return 0;

    const int n = theNodeTags.Size();
    if (n == 0)
        return 0;

    const bool faceted = (myKind == MASTER_SEGMENTS) && nodesPerSeg > 0;
    const int stride = faceted ? nodesPerSeg : 1;

    // Whole strides only: for a faceted surface the stride is the facet, so one
    // missing corner removes the facet; for SLAVE_NODES the stride is the node.
    ID kept(0, n);
    int nkept = 0;
    for (int base = 0; base + stride <= n; base += stride) {
        bool alive = true;
        for (int k = 0; k < stride; k++) {
            if (dom->getNode(theNodeTags(base + k)) == 0) {
                alive = false;
                break;
            }
        }
        if (!alive)
            continue;
        for (int k = 0; k < stride; k++)
            kept[nkept++] = theNodeTags(base + k);
    }

    const int dropped = n - nkept;
    if (dropped > 0) {
        ID tight(nkept);
        for (int i = 0; i < nkept; i++)
            tight(i) = kept(i);
        theNodeTags = tight;
    }
    return dropped;
}
