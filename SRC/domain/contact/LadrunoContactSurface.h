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
// Description: LadrunoContactSurface, the topology of one contact surface: a
// slave node set, or a master surface of tri3/quad4 segments stored as flat
// connectivity. Owned by LadrunoContactDomain.

#ifndef LadrunoContactSurface_h
#define LadrunoContactSurface_h

#include <ID.h>

class Domain;

class LadrunoContactSurface
{
  public:
    enum Kind { SLAVE_NODES = 0, MASTER_SEGMENTS = 1, SLAVE_SEGMENTS = 2 };

    // SLAVE_NODES: nodeTags are the slave nodes.
    // MASTER_SEGMENTS: nodeTags is the flat connectivity and nodesPerSeg the
    //   stride (3 = tri, 4 = quad); segment s spans nodeTags[s*nodesPerSeg .. +stride).
    // SLAVE_SEGMENTS: a faceted slave surface with the MASTER_SEGMENTS layout. The
    //   mortar formulation integrates the slave matrix D over these facets; a node
    //   set carries no facet and cannot give D.
    LadrunoContactSurface(int tag, Kind kind, const ID &nodeTags, int nodesPerSeg = 0);
    ~LadrunoContactSurface();

    int getTag(void) const { return myTag; }
    Kind getKind(void) const { return myKind; }
    const ID &getNodeTags(void) const { return theNodeTags; }
    int getNumSegments(void) const;       // faceted: segment count; SLAVE_NODES: node count
    int getNodesPerSeg(void) const { return nodesPerSeg; }

    // Drop the nodes the analysis has removed (`remove node`, element removal
    // by a collapse recorder) and return the number of node entries dropped.
    // This is sound because the contactSurface command refuses a node tag that
    // does not exist at declaration, so a node missing later can only have been
    // removed at run time. A faceted surface loses the whole segment when any of
    // its nodes is gone: a partial facet is not a facet, and the segment count is
    // derived as size/nodesPerSeg.
    int pruneMissingNodes(Domain *dom);

  private:
    int myTag;
    Kind myKind;
    ID theNodeTags;
    int nodesPerSeg;
};

#endif
