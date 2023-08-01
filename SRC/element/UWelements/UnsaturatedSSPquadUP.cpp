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

//
// Subclassed by: José A. Abell (UANDES), Francisco Pinto (UChile), Ricardo Gallardo (PUCV)
//

#include "UnsaturatedSSPquadUP.h"

#include <elementAPI.h>
#include <NDMaterial.h>
#include <ElementResponse.h>


static int num_UnsaturatedSSPquadUP = 0;

void* OPS_UnsaturatedSSPquadUP(void)
{
	if (num_UnsaturatedSSPquadUP == 0) {
		num_UnsaturatedSSPquadUP++;
		opserr << "UnsaturatedSSPquadUP element - Written: José A. Abell (UANDES), Francisco Pinto (UChile), Ricardo Gallardo (PUCV)\n";
	}

	// Pointer to an element that will be returned
	Element *theElement = 0;

	int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

	if (numRemainingInputArgs < 21) {
		opserr << "Invalid #args, want: element UnsaturatedSSPquadUP eleTag? iNode? jNode? kNode? lNode? matTag? t? fBulk? fDen? k1? k2? e? alpha? Sres? Ssat? ga? gn? gc? gl? krel_min?  <b1? b2? Pup? Plow? Pleft? Pright?>\n";
		return 0;
	}

	int iData[6];
	double dData[20];
	dData[7]  = 0.0; // Sres?
	dData[8]  = 0.0; // Ssat?
	dData[9]  = 0.0; // ga?
	dData[10] = 0.0; // gn?
	dData[11] = 0.0; // gc?
	dData[12] = 0.0; // gl?
	dData[13] = 1e-4; // krel_min?
	dData[14] = 0.0; // <b1?>
	dData[15] = 0.0; // <b2?>
	dData[16] = 0.0; // <Pup?>
	dData[17] = 0.0; // <Plow?>
	dData[18] = 0.0; // <Pleft?>
	dData[19] = 0.0; // <Pright?>

	int numData = 6; // number of data to read
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid integer data: element UnsaturatedSSPquadUP " << iData[0] << endln;
		return 0;
	}
	numRemainingInputArgs -= numData; // update remaining args

	numData = 14; // update number of data to read
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING invalid double data: element UnsaturatedSSPquadUP " << iData[0] << endln;
		return 0;
	}
	numRemainingInputArgs -= numData; // update remaining args

	int matID = iData[5];
	NDMaterial *theMaterial = OPS_GetNDMaterial(matID);
	if (theMaterial == 0) {
		opserr << "WARNING element UnsaturatedSSPquadUP " << iData[0] << endln;
		opserr << " Material: " << matID << " not found\n";
		return 0;
	}

	if (numRemainingInputArgs >= 2) { // check if there are at least 2 more args
		numData = 2;
		if (OPS_GetDoubleInput(&numData, &dData[14]) != 0) {
			opserr << "WARNING invalid optional data: element UnsaturatedSSPquadUP " << iData[0] << endln;
			return 0;
		}
		numRemainingInputArgs -= numData; // update remaining args
	}

	if (numRemainingInputArgs >= 4) { // check if there are at least 4 more args
		numData = 4;
		if (OPS_GetDoubleInput(&numData, &dData[16]) != 0) {
			opserr << "WARNING invalid optional data: element UnsaturatedSSPquadUP " << iData[0] << endln;
			return 0;
		}
		numRemainingInputArgs -= numData; // update remaining args
	}

	theElement = new UnsaturatedSSPquadUP(iData[0], iData[1], iData[2], iData[3], iData[4], *theMaterial,
	                                      dData[0], dData[1], dData[2], dData[3], dData[4],  dData[5], dData[6],
	                                      dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
	                                      dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19]);

	if (theElement == 0) {
		opserr << "WARNING could not create element of type UnsaturatedSSPquadUP\n";
		return 0;
	}

	return theElement;
}

// full constructor
UnsaturatedSSPquadUP::UnsaturatedSSPquadUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterial &theMat,
        double thick, double Kf, double Rf, double k1, double k2,
        double eVoid, double alpha,
        double Sres, double Ssat, double ga, double gn, double gc, double gl, double krel_min,
        double b1, double b2,
        double Pup, double Plow, double Pleft, double Pright)
	: SSPquadUP(tag, Nd1, Nd2, Nd3, Nd4, theMat,
	            thick, Kf, Rf, k1, k2,
	            eVoid, alpha, b1, b2,
	            Pup, Plow, Pleft, Pright,
	            ELE_TAG_UnsaturatedSSPquadUP),
	  mSres(Sres), mSsat(Ssat), mga(ga), mgn(gn), mgc(gc), mgl(gl), mkrel_min(krel_min)
{

}

// null constructor
UnsaturatedSSPquadUP::UnsaturatedSSPquadUP()
	: SSPquadUP(ELE_TAG_UnsaturatedSSPquadUP),
	  mSres(0.0), mSsat(0.0), mga(0.0), mgn(0.0), mgc(0.0), mgl(0.0), mkrel_min(0.0)
{
}


Response*
UnsaturatedSSPquadUP::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;

	char outputData[32];

	output.tag("ElementOutput");
	output.attr("eleType", "BrickUP");
	output.attr("eleTag", this->getTag());
	for (int i = 1; i <= 4; i++)
	{
		sprintf(outputData, "node%d", i);
		output.attr(outputData, theNodes[i - 1]->getTag());
	}

	if (strcmp(argv[0], "saturation") == 0) {
		sprintf(outputData, "Saturation");
		output.tag("ResponseType", outputData);

		theResponse = new ElementResponse(this, 1010, Vector(1));
	}
	else
	{
		//If nothing here forward to the base class
		theResponse = SSPquadUP::setResponse(argv, argc, output);
	}

	output.endTag(); // ElementOutput

	return theResponse;
}

int
UnsaturatedSSPquadUP::getResponse(int responseID, Information &eleInfo)
{

	static Vector saturation_vec(1);

	if (responseID == 1010){
		//Compute fluid pressure as average of nodal values
		const double p0 = theNodes[0]->getTrialVel()(2);
		const double p1 = theNodes[1]->getTrialVel()(2);
		const double p2 = theNodes[2]->getTrialVel()(2);
		const double p3 = theNodes[3]->getTrialVel()(2);

		double pw = (p0 + p1 + p2 + p3) / 4;

		pw = pw < 0 ? pw : 0;

		// compute compressibility matrix term
		double glocal = sqrt(b[0] * b[0] + b[1] * b[1]);
		double psi = -pw / (fDens * glocal);
		double S = getRelativeSaturation(psi);
		saturation_vec(0) = S;
		return eleInfo.setVector(saturation_vec);
	}
	else
	{
		return SSPquadUP::getResponse(responseID, eleInfo);
	}

	return -1;
}


const Matrix &
UnsaturatedSSPquadUP::getMass(void)
{
	mMass.Zero();

	//Compute fluid pressure as average of nodal values
	const double p0 = theNodes[0]->getTrialVel()(2);
	const double p1 = theNodes[1]->getTrialVel()(2);
	const double p2 = theNodes[2]->getTrialVel()(2);
	const double p3 = theNodes[3]->getTrialVel()(2);

	double pw = (p0 + p1 + p2 + p3) / 4;

	pw = pw < 0 ? pw : 0;

	// compute compressibility matrix term
	double glocal = sqrt(b[0] * b[0] + b[1] * b[1]);
	double psi = -pw / (fDens * glocal);
	double S = getRelativeSaturation(psi);
	double dSdpw = getRelativeSaturationPressureDerivative(pw);
	double oneOverQ = -0.25 * J0 * mThickness * (mPorosity / fBulk * S - mPorosity * dSdpw); // Eq 58 Plaxis
	// double oneOverQ = -0.25 * J0 * mThickness * mPorosity / fBulk;

	// get mass density from the material
	double density = theMaterial->getRho();

	// transpose the shape function derivative array
	Matrix dNp(2, 4);
	dNp(0, 0) = dN(0, 0); dNp(0, 1) = dN(1, 0); dNp(0, 2) = dN(2, 0); dNp(0, 3) = dN(3, 0);
	dNp(1, 0) = dN(0, 1); dNp(1, 1) = dN(1, 1); dNp(1, 2) = dN(2, 1); dNp(1, 3) = dN(3, 1);

	// compute stabilization matrix for incompressible problems
	Matrix Kp(4, 4);
	Kp = -4.0 * mAlpha * J0 * mThickness * dN * dNp;

	// return zero matrix if density is zero
	if (density == 0.0) {
		return mMass;
	}

	// full mass matrix for the element [ M  0 ]
	//  includes M and S submatrices    [ 0 -S ]
	for (int i = 0; i < 4; i++) {

		int I    = 2 * i;
		int Ip1  = 2 * i + 1;
		int II   = 3 * i;
		int IIp1 = 3 * i + 1;
		int IIp2 = 3 * i + 2;

		for (int j = 0; j < 4; j++) {

			int J    = 2 * j;
			int Jp1  = 2 * j + 1;
			int JJ   = 3 * j;
			int JJp1 = 3 * j + 1;
			int JJp2 = 3 * j + 2;

			mMass(II, JJ)     = mSolidM(I, J);
			mMass(IIp1, JJ)   = mSolidM(Ip1, J);
			mMass(IIp1, JJp1) = mSolidM(Ip1, Jp1);
			mMass(II, JJp1)   = mSolidM(I, Jp1);

			// contribution of compressibility matrix
			mMass(IIp2, JJp2) = Kp(i, j) + oneOverQ;
		}
	}

	return mMass;
}


const Vector &
UnsaturatedSSPquadUP::getResistingForce(void)
// this function computes the resisting force vector for the element
{
	Vector f1(8);
	Vector f2(4);
	Vector mStress(3);

	// get stress from the material
	mStress = theMaterial->getStress();

	// get trial displacement
	const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
	const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
	const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
	const Vector &mDisp_4 = theNodes[3]->getTrialDisp();

	Vector d(8);
	d(0) = mDisp_1(0);
	d(1) = mDisp_1(1);
	d(2) = mDisp_2(0);
	d(3) = mDisp_2(1);
	d(4) = mDisp_3(0);
	d(5) = mDisp_3(1);
	d(6) = mDisp_4(0);
	d(7) = mDisp_4(1);

	// add stabilization force to internal force vector
	f1 = Kstab * d;

	// add internal force from the stress
	f1.addMatrixTransposeVector(1.0, Mmem, mStress, 4.0 * mThickness * J0);

	// get mass density from the material
	double density = theMaterial->getRho();

	// subtract body forces from internal force vector
	double xi[4];
	double eta[4];
	xi[0]  = -1.0; xi[1]  =  1.0; xi[2]  = 1.0; xi[3]  = -1.0;
	eta[0] = -1.0; eta[1] = -1.0; eta[2] = 1.0; eta[3] =  1.0;

	if (applyLoad == 0) {
		for (int i = 0; i < 4; i++) {
			f1(2 * i)   -= density * b[0] * mThickness * (J0 + J1 * xi[i] + J2 * eta[i]);
			f1(2 * i + 1) -= density * b[1] * mThickness * (J0 + J1 * xi[i] + J2 * eta[i]);
		}
	} else {
		for (int i = 0; i < 4; i++) {
			f1(2 * i)   -= density * appliedB[0] * mThickness * (J0 + J1 * xi[i] + J2 * eta[i]);
			f1(2 * i + 1) -= density * appliedB[1] * mThickness * (J0 + J1 * xi[i] + J2 * eta[i]);
		}
	}

	// account for fluid body forces
	Matrix k(2, 2);
	Vector body(2);
	// permeability tensor
	k(0, 0) = perm[0];
	k(1, 1) = perm[1];
	// body force vector
	if (applyLoad == 0) {
		body(0) = b[0];
		body(1) = b[1];
	} else {
		body(0) = appliedB[0];
		body(1) = appliedB[1];
	}

	//Compute fluid pressure as average of nodal values
	const double p0 = theNodes[0]->getTrialVel()(2);
	const double p1 = theNodes[1]->getTrialVel()(2);
	const double p2 = theNodes[2]->getTrialVel()(2);
	const double p3 = theNodes[3]->getTrialVel()(2);

	double pw = (p0 + p1 + p2 + p3) / 4;

	pw = pw < 0 ? pw : 0;

	// compute compressibility matrix term
	double glocal = sqrt(b[0] * b[0] + b[1] * b[1]);
	double psi = -pw / (fDens * glocal);
	double S = getRelativeSaturation(psi);
	double krel = getRelativePermeability(S);
	f2 = 4.0 * J0 * mThickness * fDens * dN * krel / (fDens * glocal) * k * body; // Eq. 59 Plaxis
	// f2 = 4.0 * J0 * mThickness * fDens * dN * k * body;

	// assemble full internal force vector for the element
	mInternalForces(0)  = f1(0);
	mInternalForces(1)  = f1(1);
	mInternalForces(2)  = f2(0);
	mInternalForces(3)  = f1(2);
	mInternalForces(4)  = f1(3);
	mInternalForces(5)  = f2(1);
	mInternalForces(6)  = f1(4);
	mInternalForces(7)  = f1(5);
	mInternalForces(8)  = f2(2);
	mInternalForces(9)  = f1(6);
	mInternalForces(10) = f1(7);
	mInternalForces(11) = f2(3);

	//LM change
	// Subtract pressure loading from internal force vector
	if (pressureUpperSide != 0.0 || pressureLowerSide != 0.0 || pressureLeftSide != 0.0 || pressureRightSide != 0.0) {
		mInternalForces.addVector(1.0, pressureLoad, -1.0);
	}
	//LM change

	// inertial unbalance load
	mInternalForces.addVector(1.0, Q, -1.0);

	return mInternalForces;
}

void
UnsaturatedSSPquadUP::GetPermeabilityMatrix(void)
// this function computes the permeability matrix for the element
{
	mPerm.Zero();
	Matrix k(2, 2);
	Matrix dNp(2, 4);

	// permeability tensor
	k(0, 0) = perm[0];
	k(1, 1) = perm[1];

	// transpose the shape function derivative array
	dNp(0, 0) = dN(0, 0); dNp(0, 1) = dN(1, 0); dNp(0, 2) = dN(2, 0); dNp(0, 3) = dN(3, 0);
	dNp(1, 0) = dN(0, 1); dNp(1, 1) = dN(1, 1); dNp(1, 2) = dN(2, 1); dNp(1, 3) = dN(3, 1);

	// compute permeability matrix
	const double p0 = theNodes[0]->getTrialVel()(2);
	const double p1 = theNodes[1]->getTrialVel()(2);
	const double p2 = theNodes[2]->getTrialVel()(2);
	const double p3 = theNodes[3]->getTrialVel()(2);

	double pw = (p0 + p1 + p2 + p3) / 4;

	pw = pw < 0 ? pw : 0;

	// compute compressibility matrix term
	double glocal = sqrt(b[0] * b[0] + b[1] * b[1]);
	double psi = -pw / (fDens * glocal);
	double S = getRelativeSaturation(psi);
	double krel = getRelativePermeability(S);


	mPerm.addMatrixTripleProduct(1.0, dNp, krel / (fDens * glocal)*k, 4.0 * J0 * mThickness); // Eq 57 Plaxis
	// mPerm.addMatrixTripleProduct(1.0, dNp, k, 4.0 * J0 * mThickness);

	return;
}
