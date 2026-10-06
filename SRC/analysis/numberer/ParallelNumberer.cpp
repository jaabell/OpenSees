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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-02-14 18:43:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/numberer/ParallelNumberer.cpp,v $                                                                        

// Written: fmk 
// Revision: A
//
// Description: This file contains the implementation for ParallelNumberer.
//
// What: "@(#) ParallelNumberer.C, revA"

#include <ParallelNumberer.h>
#include <AnalysisModel.h>

#include <Domain.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <GraphNumberer.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <Node.h>
#include <ArrayOfTaggedObjects.h>
#include <stdlib.h>

#include <cstddef>
#include <unordered_map>
#include <vector>


ParallelNumberer::ParallelNumberer(int dTag, int numSub, Channel **theC) 
  :DOF_Numberer(NUMBERER_TAG_ParallelNumberer), 
   processID(dTag), numChannels(numSub)

{
  theChannels = new Channel *[numSub];
  for (int i=0; i<numSub; i++)
    theChannels[i] = theC[i];
}

ParallelNumberer::ParallelNumberer() 
  :DOF_Numberer(NUMBERER_TAG_ParallelNumberer), theNumberer(0),
   processID(0), numChannels(0), theChannels(0)
{
  
}


ParallelNumberer::ParallelNumberer(GraphNumberer &theGraphNumberer) 
  :DOF_Numberer(NUMBERER_TAG_ParallelNumberer), theNumberer(&theGraphNumberer),
   processID(0), numChannels(0), theChannels(0)
{
  
}


ParallelNumberer::~ParallelNumberer() 
{
  if (theChannels != 0)
    delete [] theChannels;

  if (theNumberer != 0)
    delete theNumberer;
}


int
ParallelNumberer::setProcessID(int dTag) 
{
  processID = dTag;
  return 0;
}

int
ParallelNumberer::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[numChannels];
  for (int i=0; i<numChannels; i++)
    theChannels[i] = theC[i];

  return 0;
}



namespace {

// Maps a reference tag (the node tag of a DOF_Group) to its position in the
// vertexRefs list in O(1), replacing the linear ID::getLocation() search of
// mergeSubGraph(). Node tags are usually close to contiguous, so a direct
// address table is used; if the tags turn out to be much sparser than the
// number of vertices the table is replaced by a hash map. Negative reference
// tags (DOF_Groups without a node) are always kept in the hash map. The first
// position recorded for a tag is kept, as getLocation() returns the first
// match.
class RefTagIndex
{
 public:
  explicit RefTagIndex(std::size_t expectedVertices)
    :expected(expectedVertices > 0 ? expectedVertices : 1), count(0),
     useHash(false)
  {}

  int find(int ref) const {
    if (ref >= 0 && !useHash) {
      if (static_cast<std::size_t>(ref) < dense.size())
	return dense[ref];
      return -1;
    }
    std::unordered_map<int, int>::const_iterator it = sparse.find(ref);
    return (it == sparse.end()) ? -1 : it->second;
  }

  void insert(int ref, int loc) {
    count++;
    if (ref < 0 || useHash) {
      sparse.insert(std::make_pair(ref, loc));
      return;
    }
    std::size_t need = static_cast<std::size_t>(ref) + 1;
    if (need > dense.size()) {
      std::size_t budget = 4 * (expected > count ? expected : count);
      if (need > budget && need > 1024) {
	// sparse tag space: move the table into the hash map
	useHash = true;
	sparse.reserve(2 * expected);
	for (std::size_t t = 0; t < dense.size(); t++)
	  if (dense[t] >= 0)
	    sparse.insert(std::make_pair(static_cast<int>(t), dense[t]));
	std::vector<int>().swap(dense);
	sparse.insert(std::make_pair(ref, loc));
	return;
      }
      std::size_t newSize = dense.empty() ? 1024 : dense.size();
      while (newSize < need)
	newSize *= 2;
      dense.resize(newSize, -1);
    }
    if (dense[ref] < 0)
      dense[ref] = loc;
  }

 private:
  std::vector<int> dense;
  std::unordered_map<int, int> sparse;
  std::size_t expected;
  std::size_t count;
  bool useHash;
};

// Same merge as ParallelNumberer::mergeSubGraph(): the vertices of the
// subgraph are visited in the same order, new merged vertices get the same
// tags in the same order and the edges are added in the same order, so the
// merged graph is identical. Only the searches differ: the reference tag
// lookup uses the RefTagIndex and the subgraph-to-merged tag lookup of the
// edge pass uses a hash map filled during the vertex pass, instead of linear
// searches of vertexRefs and theSubdomainMap, which made the merge
// quadratic in the number of vertices.
void
mergeSubGraphIndexed(Graph &theGraph, Graph &theSubGraph, ID &vertexTags,
		     ID &vertexRefs, ID &theSubdomainMap, RefTagIndex &refIndex)
{
  Vertex *subVertexPtr;
  VertexIter &theSubGraphIter1 = theSubGraph.getVertices();
  int count = 0;
  int numVertex = theGraph.getNumVertex();
  int numVertexSub = theSubGraph.getNumVertex();

  std::unordered_map<int, int> subToMerged;
  subToMerged.reserve(2 * static_cast<std::size_t>(numVertexSub));

  while ((subVertexPtr = theSubGraphIter1()) != 0) {
    int vertexTagSub = subVertexPtr->getTag();
    int vertexTagRef = subVertexPtr->getRef();
    int loc = refIndex.find(vertexTagRef);

    // see ParallelNumberer::mergeSubGraph(): DOF_Groups without a node
    // cannot be matched by reference tag
    if (vertexTagRef < 0 && loc >= 0) {
      opserr << "FATAL ParallelNumberer::numberDOF() - DOF_Group " << vertexTagSub
	     << " has no node (reference tag " << vertexTagRef << "), as for the"
	     << " Lagrange multipliers of constraints Lagrange; such DOF_Groups"
	     << " cannot be matched across processes. Use constraints"
	     << " Transformation, Penalty or Plain with a parallel numberer.\n";
      exit(-1);
    }

    int vertexTagMerged;
    if (loc < 0) {
      vertexTagMerged = theGraph.getFreeTag();
      vertexTags[numVertex] = vertexTagMerged;
      vertexRefs[numVertex] = vertexTagRef;
      Vertex *newVertex = new Vertex(vertexTagMerged, vertexTagRef,
				     subVertexPtr->getWeight(),
				     subVertexPtr->getColor());
      theGraph.addVertex(newVertex);
      refIndex.insert(vertexTagRef, numVertex);
      numVertex++;
    } else
      vertexTagMerged = vertexTags[loc];

    theSubdomainMap[count] = vertexTagSub;
    theSubdomainMap[count+numVertexSub] = vertexTagMerged;
    subToMerged[vertexTagSub] = vertexTagMerged;
    count++;
  }

  VertexIter &theSubGraphIter2 = theSubGraph.getVertices();
  while ((subVertexPtr = theSubGraphIter2()) != 0) {
    int vertexTagMerged = subToMerged[subVertexPtr->getTag()];
    const ID &adjacency = subVertexPtr->getAdjacency();
    for (int i=0; i<adjacency.Size(); i++) {
      std::unordered_map<int, int>::const_iterator it =
	subToMerged.find(adjacency(i));
      if (it == subToMerged.end()) {
	opserr << "FATAL ParallelNumberer::numberDOF() - subdomain graph vertex "
	       << subVertexPtr->getTag() << " is adjacent to unknown vertex "
	       << adjacency(i) << endln;
	exit(-1);
      }
      theGraph.addEdge(vertexTagMerged, it->second);
    }
  }
}

}  // namespace

// int numberDOF(void)
// The ParallelNumberer sitting on P0, collects each partition graph from P1 through Pn-1, 
// merges them into 1 large graph, & then numbers this graph. The ParallelNumberers sitting 
// on P1 through Pn-1 then receive the mapping info for the dof tag and dof numbering from P0.

int
ParallelNumberer::numberDOF(int lastDOF)
{
  int result = 0;

  // get a pointer to the model & check its not null
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  Domain *theDomain = 0;
  if (theModel != 0) theDomain = theModel->getDomainPtr();
  
  if (theModel == 0 || theDomain == 0) {
    opserr << "WARNING ParallelNumberer::numberDOF(int) -";
    opserr << " - no AnalysisModel - has setLinks() been invoked?\n";
    return -1;
  }
  
  if (lastDOF != -1) {
    opserr << "WARNING ParallelNumberer::numberDOF(int lastDOF):";
    opserr << " does not use the lastDOF as requested\n";
  }

  Graph &theGraph = theModel->getDOFGroupGraph();

  // if subdomain, collect graph, send it off, get 
  // ID back containing dof tags & start id numbers.
  if (processID != 0) {

    Channel *theChannel = theChannels[0];
    int numVertex = theGraph.getNumVertex();

    theGraph.sendSelf(0, *theChannel);

    // recv iD
    ID theID(2*numVertex);
    theChannel->recvID(0, 0, theID);

    // set vertex numbering based on ID received
    for (int i=0; i<numVertex; i ++) {
      int vertexTag = theID(i);
      int startID = theID(i+numVertex);
      //      Vertex *vertexPtr = theGraph.getVertexPtr(vertexTag);
      int dofTag = vertexTag;
      DOF_Group *dofPtr;	
      dofPtr = theModel->getDOF_GroupPtr(dofTag);
      if (dofPtr == 0) {
	opserr << "WARNING ParallelNumberer::numberDOF - ";
	opserr << "DOF_Group " << dofTag << "not in AnalysisModel!\n";
	result = -4;
      } else {
	const ID &theDOFID = dofPtr->getID();
	int idSize = theDOFID.Size();
	for (int j=0; j<idSize; j++)
	  if (theDOFID(j) == -2 || theDOFID(j) == -3) dofPtr->setID(j, startID++);
      }
    }

    theChannel->sendID(0, 0, theID);
  } 
  
  // if main domain, collect graphs from all subdomains,
  // merge into 1, number this one, send to subdomains the
  // id containing dof tags & start id's.
  else {

    // for P0 domain determine original vertex and ref tags
    int numVertex = theGraph.getNumVertex(); 
    int numVertexP0 = numVertex;

    ID vertexTags(numVertex);
    ID vertexRefs(numVertex);
    RefTagIndex refIndex(static_cast<std::size_t>(numVertex) *
			 static_cast<std::size_t>(numChannels + 1));
    Vertex *vertexPtr;
    int loc = 0;
    VertexIter &theVertices = theGraph.getVertices();
    bool contiguousTags = true;
    while ((vertexPtr = theVertices()) != 0) {
      vertexTags[loc] = vertexPtr->getTag();
      vertexRefs[loc] = vertexPtr->getRef();
      refIndex.insert(vertexPtr->getRef(), loc);
      if (vertexPtr->getTag() != loc)
	contiguousTags = false;
      loc++;
    }

    // The merged graph is built in a copy of the P0 graph held in an
    // ArrayOfTaggedObjects, where a vertex is found by indexing with its tag,
    // instead of in the P0 graph itself, whose map storage costs a tree search
    // for every vertex lookup of the merge and of the GraphNumberer. This
    // requires the vertex tags to equal their positions: the P0 tags are
    // 0..numVertexP0-1 and getFreeTag() then hands out the following
    // integers. The copy has the same vertices, tags, adjacency and
    // getFreeTag() sequence and the vertices are visited in the same
    // ascending tag order, so the merged graph and its numbering are
    // unchanged. If the P0 tags are not contiguous the P0 graph is used as
    // before.
    Graph *theMergedGraph = &theGraph;
    Graph *theDenseGraph = 0;
    if (contiguousTags == true) {
      std::size_t sizeHint = static_cast<std::size_t>(numVertex) *
	static_cast<std::size_t>(numChannels + 1);
      if (sizeHint < 1024)
	sizeHint = 1024;
      if (sizeHint > 0x7fffffff)
	sizeHint = 0x7fffffff;
      theDenseGraph = new Graph(*(new ArrayOfTaggedObjects(static_cast<int>(sizeHint))));
      VertexIter &theP0Vertices = theGraph.getVertices();
      while ((vertexPtr = theP0Vertices()) != 0)
	theDenseGraph->addVertex(new Vertex(vertexPtr->getTag(), vertexPtr->getRef(),
					    vertexPtr->getWeight(), vertexPtr->getColor()));
      VertexIter &theP0Edges = theGraph.getVertices();
      while ((vertexPtr = theP0Edges()) != 0) {
	const ID &adjacency = vertexPtr->getAdjacency();
	for (int i=0; i<adjacency.Size(); i++)
	  theDenseGraph->addEdge(vertexPtr->getTag(), adjacency(i));
      }
      theMergedGraph = theDenseGraph;
    }
    
    ID **theSubdomainIDs = new ID *[numChannels];
    FEM_ObjectBroker theBroker;

    // for each subdomain we receive graph, create an ID (to store
    // subdomain graph to merged graph vertex mapping and the final
    // subdoain graph vertex to startDOF mapping) and finally merge the
    // subdomain graph

    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      Graph *theSubGraph = new Graph();

      theSubGraph->recvSelf(0, *theChannel, theBroker);

      theSubdomainIDs[j] = new ID(theSubGraph->getNumVertex()*2);

      mergeSubGraphIndexed(*theMergedGraph, *theSubGraph, vertexTags, vertexRefs,
			   *theSubdomainIDs[j], refIndex);

      delete theSubGraph;
    }
    
    // we use graph numberer if one was provided in constructor,
    // otherwise we number based on subdomains (all in subdomain 1 numbered first, 
    // then  those in 2 not in 1 and so on till done.
    //    GraphNumberer *theNumberer = this->getGraphNumbererPtr();

    ID *theOrderedRefs = new ID(theMergedGraph->getNumVertex());

    if (theNumberer != 0) {

      // use the supplied graph numberer to number the merged graph
      *theOrderedRefs = theNumberer->number(*theMergedGraph, lastDOF);     

    } else {

      // assign numbers based on the subdomains. A vertex is skipped if it
      // has already been ordered. The flags give the same result as the
      // former search of theOrderedRefs: as that ID is zero-filled, the
      // vertex with tag 0 always counted as ordered and was left in the
      // last, zero-filled position. That order is kept.
      std::vector<char> ordered(static_cast<std::size_t>(theMergedGraph->getFreeTag()) + 1, 0);
      ordered[0] = 1;

      int loc = 0;
      for (int l=0; l<numChannels; l++) {
	const ID &theSubdomain = *theSubdomainIDs[l];
	int numVertexSubdomain = theSubdomain.Size()/2;

	for (int i=0; i<numVertexSubdomain; i++) {
	  int vertexTagMerged = theSubdomain(i+numVertexSubdomain);
	  if (ordered[vertexTagMerged] == 0) {
	    ordered[vertexTagMerged] = 1;
	    (*theOrderedRefs)[loc++] = vertexTagMerged;
	  }
	}
      }

      // now order those not yet ordered in p0
      for (int j=0; j<numVertexP0; j++) {
	int refTagP0 = vertexTags[j];
	if (ordered[refTagP0] == 0) {
	  ordered[refTagP0] = 1;
	  (*theOrderedRefs)[loc++] = refTagP0;
	}
      }
    }

    int count = 0;
    for (int i=0; i<theOrderedRefs->Size(); i++) {
      int vertexTag = (*theOrderedRefs)(i);
      //      int vertexTag = vertexTags[vertexRefs.getLocation(tag)];
      Vertex *vertexPtr = theMergedGraph->getVertexPtr(vertexTag);
      int numDOF= vertexPtr->getColor();
      vertexPtr->setTmp(count);
      count += numDOF;
    }

    delete theOrderedRefs;

    // number own dof's
    for (int i=0; i<numVertexP0; i++  ) {
      int vertexTag = vertexTags(i);
      Vertex *vertexPtr = theMergedGraph->getVertexPtr(vertexTag);

      int startID = vertexPtr->getTmp();
      int dofTag = vertexTag;
      DOF_Group *dofPtr;	
      dofPtr = theModel->getDOF_GroupPtr(dofTag);
      if (dofPtr == 0) {
	opserr << "WARNING ParallelNumberer::numberDOF - ";
	opserr << "DOF_Group (P0) " << dofTag << "not in AnalysisModel!\n";
	result = -4;
      } else {
	const ID &theDOFID = dofPtr->getID();
	int idSize = theDOFID.Size();
	for (int j=0; j<idSize; j++)
	  if (theDOFID(j) == -2 || theDOFID(j) == -3) dofPtr->setID(j, startID++);
      }
    }

    // now given the ordered refs we determine the mapping for each subdomain
    // and send the id with the information back to the subdomain, which it uses to order
    // it's own graph
    for (int k=0; k<numChannels; k++) {
      Channel *theChannel = theChannels[k];
      ID &theSubdomain = *theSubdomainIDs[k];
      int numVertexSubdomain = theSubdomain.Size()/2;

      for (int i=0; i<numVertexSubdomain; i++) {
	int vertexTagMerged = theSubdomain[numVertexSubdomain+i];
	Vertex *vertexPtr = theMergedGraph->getVertexPtr(vertexTagMerged);
	int startDOF = vertexPtr->getTmp();
	theSubdomain[i+numVertexSubdomain] = startDOF;
      }

      theChannel->sendID(0, 0, theSubdomain);
      theChannel->recvID(0, 0, theSubdomain);
      delete theSubdomainIDs[k];
    }      
    delete [] theSubdomainIDs;

    if (theDenseGraph != 0)
      delete theDenseGraph;
  }



  // iterate through the DOFs one last time setting any -4 values
  // iterate through  the DOFs second time setting -3 values
  // the MP_Constraints are indexed by constrained node in one pass; each
  // list keeps the getMPs() order, so a node constrained by several
  // MP_Constraints receives its IDs in the same order as with a full scan
  std::unordered_map<int, std::vector<MP_Constraint *> > mpIndex;
  {
    MP_ConstraintIter &theMPsAll = theDomain->getMPs();
    MP_Constraint *mpAll;
    while ((mpAll = theMPsAll()) != 0)
      mpIndex[mpAll->getNodeConstrained()].push_back(mpAll);
  }

  AnalysisModel *theAModel = this->getAnalysisModelPtr();
  DOF_GrpIter &tDOFs = theAModel->getDOFs();

  DOF_Group *dofPtr;
  while ((dofPtr = tDOFs()) != 0) {
    const ID &theID = dofPtr->getID();
    int have4s = 0;
    for (int i=0; i<theID.Size(); i++)
      if (theID(i) == -4) have4s = 1;
    
    if (have4s == 1) {
      int nodeID = dofPtr->getNodeTag();
      // the MP_Constraints of this node; the constraint matrix must be
      // diagonal with 1's on the diagonal
      std::unordered_map<int, std::vector<MP_Constraint *> >::iterator mpIt =
	mpIndex.find(nodeID);
      if (mpIt == mpIndex.end())
	continue;
      for (std::size_t m = 0; m < mpIt->second.size(); m++) {
	MP_Constraint *mpPtr = mpIt->second[m];
	{
	  int nodeRetained = mpPtr->getNodeRetained();
	  Node *nodeRetainedPtr = theDomain->getNode(nodeRetained);
	  DOF_Group *retainedDOF = nodeRetainedPtr->getDOF_GroupPtr();
	  const ID&retainedDOFIDs = retainedDOF->getID();
	  const ID&constrainedDOFs = mpPtr->getConstrainedDOFs();
	  const ID&retainedDOFs = mpPtr->getRetainedDOFs();
	  for (int i=0; i<constrainedDOFs.Size(); i++) {
	    int dofC = constrainedDOFs(i);
	    int dofR = retainedDOFs(i);
	    int dofID = retainedDOFIDs(dofR);
	    dofPtr->setID(dofC, dofID);
	  }
	}
      }		
    }	
  }

  // iterate through the FE_Element getting them to set their IDs
  FE_EleIter &theEle = theModel->getFEs();
  FE_Element *elePtr;
  while ((elePtr = theEle()) != 0)
    elePtr->setID();

  theModel->clearDOFGroupGraph();
  
  return result;
}


int
ParallelNumberer::mergeSubGraph(Graph &theGraph, Graph &theSubGraph, ID &vertexTags, ID &vertexRefs, ID &theSubdomainMap)
{  
  // for each vertex in the SubGraph we see if a vertex exists in the Graph which has the same
  // reference tag (Reference tag in the AnalysisModel graph is the node tag) .. if so this will be 
  // the new vertex tag for SubGraph vertex in new graph, otherwise we assign it some new vertex tag,
  // create a vertex for this new vertex tag & add it to the graph

  Vertex *subVertexPtr;
  VertexIter &theSubGraphIter1 = theSubGraph.getVertices();
  int count =0;
  int numVertex = theGraph.getNumVertex();
  int numVertexSub = theSubGraph.getNumVertex();

  while ((subVertexPtr = theSubGraphIter1()) != 0) {
    int vertexTagSub = subVertexPtr->getTag();
    int vertexTagRef = subVertexPtr->getRef();
    int loc = vertexRefs.getLocation(vertexTagRef);

    // a DOF_Group without a node (e.g. the multiplier DOFs of the Lagrange
    // constraint handler) has reference tag -1. Merging by reference tag would
    // fuse all of them into one vertex, so that different constraints share
    // the same equation numbers. Stop instead of numbering them wrongly.
    if (vertexTagRef < 0 && loc >= 0) {
      opserr << "FATAL ParallelNumberer::mergeSubGraph() - DOF_Group " << vertexTagSub
	     << " has no node (reference tag " << vertexTagRef << "), as for the"
	     << " Lagrange multipliers of constraints Lagrange; such DOF_Groups"
	     << " cannot be matched across processes. Use constraints"
	     << " Transformation, Penalty or Plain with a parallel numberer.\n";
      exit(-1);
    }

    int vertexTagMerged;
    if (loc < 0) {
      // if not already in, we will be creating a new vertex
      vertexTagMerged = theGraph.getFreeTag();
      vertexTags[numVertex] = vertexTagMerged;
      vertexRefs[numVertex] = vertexTagRef;
      Vertex *newVertex = new Vertex(vertexTagMerged, vertexTagRef, subVertexPtr->getWeight(), subVertexPtr->getColor());

      theGraph.addVertex(newVertex);
      numVertex++;
    } else
      vertexTagMerged = vertexTags[loc];

    // use the subgraphs ID to hold the mapping of vertex numbers between merged and original
    theSubdomainMap[count] = vertexTagSub;
    theSubdomainMap[count+numVertexSub] = vertexTagMerged;
    count++;
  }

  // for each vertex in subgraph, we add it's adjacenecy into the merged graph
  VertexIter &theSubGraphIter2 = theSubGraph.getVertices();
  while ((subVertexPtr = theSubGraphIter2()) != 0) {
    int vertexTagSub = subVertexPtr->getTag();
    int loc = theSubdomainMap.getLocation(vertexTagSub);
    int vertexTagMerged = theSubdomainMap[loc+numVertexSub];

    const ID &adjacency = subVertexPtr->getAdjacency();

    for (int i=0; i<adjacency.Size(); i++) {
      int vertexTagSubAdjacent = adjacency(i);
      int loc = theSubdomainMap.getLocation(vertexTagSubAdjacent);
      int vertexTagMergedAdjacent = theSubdomainMap[loc+numVertexSub];      
      theGraph.addEdge(vertexTagMerged, vertexTagMergedAdjacent);
    }
  }


  return 0;
}


int
ParallelNumberer::sendSelf(int cTag, Channel &theChannel)
{
  int sendID =0;

  // if P0 check if already sent. If already sent use old processID; if not allocate a new process 
  // id for remote part of object, enlarge channel * to hold a channel * for this remote object.

  // if not P0, send current processID

  if (processID == 0) {

    // check if already using this object
    bool found = false;
    for (int i=0; i<numChannels; i++)
      if (theChannels[i] == &theChannel) {
	sendID = i+1;
	found = true;
      }

    // if new object, enlarge Channel pointers to hold new channel * & allocate new ID
    if (found == false) {
      int nextNumChannels = numChannels + 1;
      Channel **nextChannels = new Channel *[nextNumChannels];
      if (nextNumChannels == 0) {
	opserr << "ParallelNumberer::sendSelf() - failed to allocate channel array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	nextChannels[i] = theChannels[i];
      nextChannels[numChannels] = &theChannel;
      
      numChannels = nextNumChannels;
      
      if (theChannels != 0)
	delete [] theChannels;
      
      theChannels = nextChannels;
      
      // allocate new processID for remote object
      sendID = numChannels;
    }

  } else 
    sendID = processID;


  // send remotes processID
  ID idData(1);
  idData(0) = sendID;
  
  int res = theChannel.sendID(0, cTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedSparseGenColLinSOE::sendSelf() - failed to send data\n";
    return -1;
  }

  return 0;
}

int
ParallelNumberer::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  ID idData(1);
  int res = theChannel.recvID(0, cTag, idData);
  if (res < 0) {
    opserr <<"WARNING Parallel::recvSelf() - failed to send data\n";
    return -1;
  }	      
  processID = idData(0);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  return 0;
}


int
ParallelNumberer::numberDOF(ID &lastDOFs)
{
  int result = 0;

  // get a pointer to the model & check its not null
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  Domain *theDomain = 0;
  if (theModel != 0) theDomain = theModel->getDomainPtr();
  
  if (theModel == 0 || theDomain == 0) {
    opserr << "WARNING ParallelNumberer::numberDOF(int) -";
    opserr << " - no AnalysisModel - has setLinks() been invoked?\n";
    return -1;
  }
  
  Graph &theGraph = theModel->getDOFGroupGraph();
  
  // if subdomain, collect graph, send it off, get 
  // ID back containing dof tags & start id numbers.
  if (processID != 0) {
    Channel *theChannel = theChannels[0];
    int numVertex = theGraph.getNumVertex();
    theGraph.sendSelf(0, *theChannel);
    ID theID(2*numVertex);
    theChannel->recvID(0, 0, theID);
    for (int i=0; i<numVertex; i += 2) {
      int dofTag = theID(i);
      int startID = theID(i+1);
      DOF_Group *dofPtr;	
      dofPtr = theModel->getDOF_GroupPtr(dofTag);
      if (dofPtr == 0) {
	opserr << "WARNING ParallelNumberer::numberDOF - ";
	opserr << "DOF_Group " << dofTag << "not in AnalysisModel!\n";
	result = -4;
      } else {
	const ID &theID = dofPtr->getID();
	int idSize = theID.Size();
	for (int j=0; j<idSize; j++)
	  if (theID(j) == -2) dofPtr->setID(j, startID++);
      }
    }
  } 
  
  // if main domain, collect graphs from all subdomains,
  // merge into 1, number this one, send to subdomains the
  // id containing dof tags & start id's.
  else {
    
    // determine original vertex and ref tags
    int numVertex = theGraph.getNumVertex();
    ID vertexTags(numVertex);
    ID vertexRefs(numVertex);
    Vertex *vertexPtr;
    int loc = 0;
    VertexIter &theVertices = theGraph.getVertices();
    while ((vertexPtr = theVertices()) != 0) {
      vertexTags[loc] = vertexPtr->getTag();
      vertexRefs[loc] = vertexPtr->getRef();
      loc++;
    }
    
    ID **theSubdomainIDs = new ID *[numChannels];
    FEM_ObjectBroker theBroker;

    // merge all subdomain graphs
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      Graph theSubGraph;
      theSubGraph.recvSelf(0, *theChannel, theBroker);
      
      theSubdomainIDs[j] = new ID(theSubGraph.getNumVertex()*2);
      this->mergeSubGraph(theGraph, theSubGraph, vertexTags, vertexRefs, *theSubdomainIDs[j]);
    }

    // number the merged graph
    //    result =  this->DOF_Numberer::number(theGraph);

    // send results of numbered back to subdomains
    for (int k=0; k<numChannels; k++) {
      Channel *theChannel = theChannels[k];
      // this->determineSubIDs
      theChannel->sendID(0, 0, *theSubdomainIDs[k]);
      delete theSubdomainIDs[k];
    }      
    delete [] theSubdomainIDs;
    
    // number own dof's
  }

  theModel->clearDOFGroupGraph();

  return result;
}
