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

// Written: M. Salehi opensees.net@gmail.com
// website : http://opensees.net
// Created: 02/19
// Revision: A


#include <PARDISOGenLinSOE.h>
#include <PARDISOGenLinSolver.h>

#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <stdlib.h>
#include <new>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

// The half-storage symmetry check in addA() runs in bursts: it inspects the
// element matrices of ASYM_CHECK_BUDGET tangent assemblies, then sleeps until
// the next multiple of ASYM_RESAMPLE_PERIOD assemblies (see zeroA()). This
// keeps its cost a small fixed fraction of assembly while still catching a
// tangent that turns unsymmetric part way through a run.
static const int ASYM_CHECK_BUDGET    = 3;
static const int ASYM_RESAMPLE_PERIOD = 64;

PARDISOGenLinSOE::PARDISOGenLinSOE(PARDISOGenLinSolver &the_Solver)
	:LinearSOE(the_Solver, LinSOE_TAGS_PARDISOGenLinSOE),
	size(0), nnz(0), A(0), B(0), X(0), colA(0), rowStartA(0),
	vectX(0), vectB(0),
	Asize(0), Bsize(0),
	factored(false), matType(0), asymWarned(0), asymBudget(0), asymPass(0),
	missWarned(0)
{
	the_Solver.setLinearSOE(*this);
}


PARDISOGenLinSOE::PARDISOGenLinSOE(PARDISOGenLinSolver &the_Solver, int _matType)
	:LinearSOE(the_Solver, LinSOE_TAGS_PARDISOGenLinSOE),
	size(0), nnz(0), A(0), B(0), X(0), colA(0), rowStartA(0),
	vectX(0), vectB(0),
	Asize(0), Bsize(0),
	factored(false), matType(_matType), asymWarned(0), asymBudget(0), asymPass(0),
	missWarned(0)
{
	if (matType < 0 || matType > 2) {
		opserr << "WARNING PARDISOGenLinSOE - unknown matrixType (" << matType
		       << "); unsymmetric storage assumed\n";
		matType = 0;
	}
	the_Solver.setLinearSOE(*this);
}


PARDISOGenLinSOE::~PARDISOGenLinSOE()
{
	if (A != 0) delete[] A;
	if (B != 0) delete[] B;
	if (X != 0) delete[] X;
	if (rowStartA != 0) delete[] rowStartA;
	if (colA != 0) delete[]colA;
	if (vectX != 0) delete vectX;
	if (vectB != 0) delete vectB;
}


int
PARDISOGenLinSOE::getNumEqn(void) const
{
	return size;
}

int
PARDISOGenLinSOE::getMatType(void) const
{
	return matType;
}

int
PARDISOGenLinSOE::setSize(Graph &theGraph)
{

	int result = 0;
	int oldSize = size;
	size = theGraph.getNumVertex();

	// fist itearte through the vertices of the graph to get nnz
	Vertex *theVertex;
	int newNNZ = 0;
	VertexIter &theVertices = theGraph.getVertices();
	while ((theVertex = theVertices()) != 0) {
		const ID &theAdjacency = theVertex->getAdjacency();
		newNNZ += theAdjacency.Size() + 1; // the +1 is for the diag entry
	}

	// Half storage: the graph adjacency is symmetric, so the off-diagonal
	// count halves and the diagonal is kept whole (as in MumpsSOE::setSize).
	if (matType != 0) {
		newNNZ -= size;
		newNNZ /= 2;
		newNNZ += size;
	}

	nnz = newNNZ;

	// arm the symmetry check for the first assemblies of this pattern
	asymPass = 0;
	if (matType != 0 && asymWarned == 0)
		asymBudget = ASYM_CHECK_BUDGET;

	// Plain new throws instead of returning 0, so the out-of-memory branch
	// below could never run. Use nothrow and leave no dangling pointers, so a
	// model that does not fit reports itself instead of terminating.
	if (newNNZ > Asize) { // we have to get more space for A and colA
		if (A != 0)
			delete[] A;
		if (colA != 0)
			delete[] colA;
		A = 0; colA = 0;

		A = new (std::nothrow) double[newNNZ];
		colA = new (std::nothrow) int[newNNZ];

		if (A == 0 || colA == 0) {
			opserr << "WARNING PARDISOGenLinSOE::setSize :";
			opserr << " ran out of memory for A and colA with nnz = ";
			opserr << newNNZ << " ("
			       << (newNNZ * (sizeof(double) + sizeof(int))) / (1024.0 * 1024.0)
			       << " MB requested)\n";
			if (A != 0) { delete[] A; A = 0; }
			if (colA != 0) { delete[] colA; colA = 0; }
			size = 0; Asize = 0; nnz = 0;
			return -1;
		}

		Asize = newNNZ;
	}

	// zero the matrix
	for (int i = 0; i < Asize; i++)
		A[i] = 0;

	factored = false;

	if (size > Bsize) { // we have to get space for the vectors

	// delete the old
		if (B != 0) delete[] B;
		if (X != 0) delete[] X;
		if (rowStartA != 0) delete[] rowStartA;
		B = 0; X = 0; rowStartA = 0;

		// create the new
		B = new (std::nothrow) double[size];
		X = new (std::nothrow) double[size];
		rowStartA = new (std::nothrow) int[size + 1];

		if (B == 0 || X == 0 || rowStartA == 0) {
			opserr << "WARNING PARDISOGenLinSOE::setSize :";
			opserr << " ran out of memory for vectors (size) (";
			opserr << size << ") \n";
			size = 0; Bsize = 0;
			return -1;
		}
		else
			Bsize = size;
	}

	// zero the vectors
	for (int j = 0; j < size; j++) {
		B[j] = 0;
		X[j] = 0;
	}

	// create new Vectors objects
	if (size != oldSize) {
		if (vectX != 0)
			delete vectX;

		if (vectB != 0)
			delete vectB;

		vectX = new Vector(X, size);
		vectB = new Vector(B, size);
	}

	// fill in rowStartA and colA
	if (size != 0) {
		rowStartA[0] = 0 + 1;
		int startLoc = 0;
		int lastLoc = 0;
		for (int a = 0; a < size; a++) {

			theVertex = theGraph.getVertexPtr(a);
			if (theVertex == 0) {
				opserr << "WARNING:PARDISOGenLinSOE::setSize :";
				opserr << " vertex " << a << " not in graph! - size set to 0\n";
				size = 0;
				return -1;
			}

			int vertexTag = theVertex->getTag();

			// nnz is a prediction that is exact only for a symmetric,
			// self-loop-free graph; bound the fill so a violation is an error
			// rather than a heap overflow of colA.
			if (lastLoc >= nnz) {
				opserr << "WARNING:PARDISOGenLinSOE::setSize : CSR fill would "
				          "overrun nnz=" << nnz << " at row " << a
				       << " (matType " << matType << ")\n";
				size = 0;
				return -1;
			}
			colA[lastLoc++] = vertexTag + 1; // place diag in first fortran index start at 1
			const ID &theAdjacency = theVertex->getAdjacency();
			int idSize = theAdjacency.Size();

			// now we have to place the entries in the ID into order in colA
			// (PARDISO requires ascending column indices within each row)
			for (int i = 0; i < idSize; i++) {

				int row = theAdjacency(i);

				// upper triangle only when half-storing; the diagonal is
				// already first and smaller than every entry kept here
				if (matType != 0 && row <= vertexTag)
					continue;

				if (lastLoc >= nnz) {
					opserr << "WARNING:PARDISOGenLinSOE::setSize : CSR fill "
					          "would overrun nnz=" << nnz << " in row " << a
					       << " (matType " << matType << ")\n";
					size = 0;
					return -1;
				}

				bool foundPlace = false;
				// find a place in colA for current col
				for (int j = startLoc; j < lastLoc; j++)
					if (colA[j] > row + 1) {
						// move the entries already there one further on
						// and place col in current location
						for (int k = lastLoc; k > j; k--)

							colA[k] = colA[k - 1];
						colA[j] = row + 1;
						foundPlace = true;
						j = lastLoc;
					}
				if (foundPlace == false) // put in at the end
					colA[lastLoc] = row + 1;

				lastLoc++;
			}
			rowStartA[a + 1] = lastLoc + 1;
			startLoc = lastLoc;
		}

		// Check the CSR layout PARDISO requires instead of assuming it: every
		// row non-empty, column indices strictly ascending, diagonal present
		// (even when zero). A violation can otherwise return a plausible wrong
		// solution rather than an error. The check is O(nnz).
		if (lastLoc != nnz) {
			opserr << "WARNING:PARDISOGenLinSOE::setSize : filled " << lastLoc
			       << " entries but nnz=" << nnz << " (matType " << matType
			       << ")\n";
			size = 0;
			return -1;
		}
		for (int a = 0; a < size; a++) {
			int rs = rowStartA[a] - 1, re = rowStartA[a + 1] - 1;
			bool sawDiag = false;
			if (rs >= re) {
				opserr << "WARNING:PARDISOGenLinSOE::setSize : row " << a
				       << " is empty - PARDISO requires a diagonal entry in every row\n";
				size = 0;
				return -1;
			}
			for (int k = rs; k < re; k++) {
				if (k > rs && colA[k] <= colA[k - 1]) {
					opserr << "WARNING:PARDISOGenLinSOE::setSize : row " << a
					       << " column indices are not strictly ascending at "
					       << k << " - PARDISO requires ascending CSR\n";
					size = 0;
					return -1;
				}
				if (colA[k] == a + 1)
					sawDiag = true;
			}
			if (sawDiag == false) {
				opserr << "WARNING:PARDISOGenLinSOE::setSize : row " << a
				       << " has no diagonal entry - PARDISO requires one, even if zero\n";
				size = 0;
				return -1;
			}
		}
	}

	// invoke setSize() on the Solver
	LinearSOESolver *the_Solver = this->getSolver();
	int solverOK = the_Solver->setSize();
	if (solverOK < 0) {
		opserr << "WARNING:PARDISOGenLinSOE::setSize :";
		opserr << " solver failed setSize()\n";
		return solverOK;
	}
	return result;
}

// Locate the one-based column col1 in the ascending run colA[lo, hi) and
// return its offset into A, or -1 if the pattern has no such entry. Binary
// search is valid because setSize() rejects any row whose columns are not
// strictly ascending. It replaces a linear scan of the row, which made the
// scatter of a 3D solid element cost O(n_dof^2 * row length).
static inline int
pardiso_findCol(const int *colA, int lo, int hi, int col1)
{
	while (lo < hi) {
		const int mid = lo + ((hi - lo) >> 1);
		const int c = colA[mid];
		if (c < col1)
			lo = mid + 1;
		else if (c > col1)
			hi = mid;
		else
			return mid;
	}
	return -1;
}


int
PARDISOGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
	// check for a quick return
	if (fact == 0.0)
		return 0;

	int idSize = id.Size();

	// check that m and id are of similar size
	if (idSize != m.noRows() && idSize != m.noCols()) {
		opserr << "PARDISOGenLinSOE::addA() ";
		opserr << " - Matrix and ID not of similar sizes\n";
		return -1;
	}

	// when half-storing, only the col >= row entries have a home in A
	const bool halfStore = (matType != 0);

	// Half storage discards the lower triangle. Compare it with its mirror
	// while it is in cache and warn once if the element matrix is
	// unsymmetric, since the run would then solve the reflected system. The
	// check is sampled (see zeroA), so it is a diagnostic, not a guarantee.
	if (halfStore && asymWarned == 0 && asymBudget > 0 &&
	    idSize == m.noRows() && idSize == m.noCols()) {
		double worstDev = 0.0, scale = 0.0;
		for (int i = 0; i < idSize; i++)
			for (int j = 0; j < i; j++) {
				double u = m(i, j), v = m(j, i);
				double au = u < 0.0 ? -u : u;
				double av = v < 0.0 ? -v : v;
				if (au > scale) scale = au;
				if (av > scale) scale = av;
				double d = u - v;
				if (d < 0.0) d = -d;
				if (d > worstDev) worstDev = d;
			}
		if (scale > 0.0 && worstDev > 1.0e-8 * scale) {
			asymWarned = 1;
			opserr << "WARNING PARDISOGenLinSOE: an assembled element matrix is "
			          "unsymmetric (max deviation " << worstDev << " vs scale "
			       << scale << ") at tangent assembly " << asymPass
			       << " of this pattern, but `system Pardiso -matrixType "
			       << matType << "` stores only the upper triangle.\n"
			          "     The lower half is discarded, so the run solves the "
			          "reflected system and may converge to a wrong answer.\n"
			          "     Use the default `system Pardiso` (-matrixType 0) for "
			          "contact, non-associated flow, follower loads or "
			          "corotational transformations. Reported once.\n";
		}
	}

	// An entry with no slot in the pattern would be dropped silently;
	// record it and report once below.
	int missing = 0;

	if (fact == 1.0) { // do not need to multiply
		for (int i = 0; i < idSize; i++) {
			int row = id(i);
			if (row < size && row >= 0) {
				int startRowLoc = rowStartA[row] - 1;
				int endRowLoc = rowStartA[row + 1] - 1;
				for (int j = 0; j < idSize; j++) {
					int col = id(j);
					if (col < size && col >= 0 && (!halfStore || col >= row)) {
						// find place in A using colA
						const int k = pardiso_findCol(colA, startRowLoc,
						                              endRowLoc, col + 1);
						if (k >= 0)
							A[k] += m(i, j);
						else
							missing = 1;
					}
				}  // for j
			}
		}  // for i
	}
	else {
		for (int i = 0; i < idSize; i++) {
			int row = id(i);
			if (row < size && row >= 0) {
				int startRowLoc = rowStartA[row] - 1;
				int endRowLoc = rowStartA[row + 1] - 1;
				for (int j = 0; j < idSize; j++) {
					int col = id(j);
					if (col < size && col >= 0 && (!halfStore || col >= row)) {
						// find place in A using colA
						const int k = pardiso_findCol(colA, startRowLoc,
						                              endRowLoc, col + 1);
						if (k >= 0)
							A[k] += fact * m(i, j);
						else
							missing = 1;
					}
				}  // for j
			}
		}  // for i
	}

	if (missing == 1 && missWarned == 0) {
		missWarned = 1;
		opserr << "WARNING PARDISOGenLinSOE::addA() - an element matrix entry has "
		          "no slot in the CSR pattern and is being discarded; the "
		          "sparsity pattern and the assembly disagree. Reported once.\n";
	}

	return 0;
}


int
PARDISOGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
	// check for a quick return
	if (fact == 0.0)  return 0;

	int idSize = id.Size();
	// check that m and id are of similar size
	if (idSize != v.Size()) {
		opserr << "PARDISOGenLinSOE::addB() ";
		opserr << " - Vector and ID not of similar sizes\n";
		return -1;
	}

	if (fact == 1.0) { // do not need to multiply if fact == 1.0
		for (int i = 0; i < idSize; i++) {
			int pos = id(i);
			if (pos < size && pos >= 0)
				B[pos] += v(i);
		}
	}
	else if (fact == -1.0) { // do not need to multiply if fact == -1.0
		for (int i = 0; i < idSize; i++) {
			int pos = id(i);
			if (pos < size && pos >= 0)
				B[pos] -= v(i);
		}
	}
	else {
		for (int i = 0; i < idSize; i++) {
			int pos = id(i);
			if (pos < size && pos >= 0)
				B[pos] += v(i) * fact;
		}
	}

	return 0;
}


int
PARDISOGenLinSOE::setB(const Vector &v, double fact)
{
	// check for a quick return
	if (fact == 0.0)  return 0;


	if (v.Size() != size) {
		opserr << "WARNING BandGenLinSOE::setB() -";
		opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
		return -1;
	}

	if (fact == 1.0) { // do not need to multiply if fact == 1.0
		for (int i = 0; i < size; i++) {
			B[i] = v(i);
		}
	}
	else if (fact == -1.0) {
		for (int i = 0; i < size; i++) {
			B[i] = -v(i);
		}
	}
	else {
		for (int i = 0; i < size; i++) {
			B[i] = v(i) * fact;
		}
	}
	return 0;
}

void
PARDISOGenLinSOE::zeroA(void)
{
	double *Aptr = A;
	for (int i = 0; i < Asize; i++)
		*Aptr++ = 0;

	factored = false;

	// zeroA runs once per tangent assembly: spend one unit of the symmetry
	// check budget here and re-arm it every ASYM_RESAMPLE_PERIOD assemblies
	// until an asymmetry has been reported.
	asymPass++;

	if (matType != 0 && asymWarned == 0 &&
	    (asymPass % ASYM_RESAMPLE_PERIOD) == 0)
		asymBudget = ASYM_CHECK_BUDGET;
	else if (asymBudget > 0)
		asymBudget--;
}

void
PARDISOGenLinSOE::zeroB(void)
{
	double *Bptr = B;
	for (int i = 0; i < size; i++)
		*Bptr++ = 0;
}

void
PARDISOGenLinSOE::setX(int loc, double value)
{
	if (loc < size && loc >= 0)
		X[loc] = value;
}

void
PARDISOGenLinSOE::setX(const Vector &x)
{
	if (x.Size() == size && vectX != 0)
		*vectX = x;
}

// getX()/getB() before setSize() has run (for example after an analysis
// that failed in the constraint handler) used to call exit(-1). Return an
// empty vector with a warning instead.
static const Vector &
pardiso_emptyVector(void)
{
	static Vector theEmptyVector;
	return theEmptyVector;
}

const Vector &
PARDISOGenLinSOE::getX(void)
{
	if (vectX == 0) {
		opserr << "WARNING PARDISOGenLinSOE::getX - the SOE has not been sized "
		          "yet; returning an empty Vector\n";
		return pardiso_emptyVector();
	}
	return *vectX;
}

const Vector &
PARDISOGenLinSOE::getB(void)
{
	if (vectB == 0) {
		opserr << "WARNING PARDISOGenLinSOE::getB - the SOE has not been sized "
		          "yet; returning an empty Vector\n";
		return pardiso_emptyVector();
	}
	return *vectB;
}

double
PARDISOGenLinSOE::normRHS(void)
{
	double norm = 0.0;
	for (int i = 0; i < size; i++) {
		double Yi = B[i];
		norm += Yi * Yi;
	}
	return sqrt(norm);

}


int
PARDISOGenLinSOE::setPARDISOGenLinSolver(PARDISOGenLinSolver &newSolver)
{
	newSolver.setLinearSOE(*this);

	if (size != 0) {
		int solverOK = newSolver.setSize();
		if (solverOK < 0) {
			opserr << "WARNING:PARDISOGenLinSOE::setSolver :";
			opserr << "the new solver could not setSeize() - staying with old\n";
			return -1;
		}
	}

	return this->LinearSOE::setSolver(newSolver);
}


int
PARDISOGenLinSOE::sendSelf(int cTag, Channel &theChannel)
{
	return 0;
}

int
PARDISOGenLinSOE::recvSelf(int cTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	return 0;
}
