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

// Created: Chris McGann, UW, 10.2011
//
// Description: This file contains the implementation of the UnsaturatedSSPbrickUP class
//                Stabilized Single-Point Brick element with a u-p formulation
//                for 3D analysis of saturated porous media
//
// Reference:   Zienkiewicz, O.C. and Shiomi, T. (1984). "Dynamic behavior of
//                saturated porous media; the generalized Biot formulation and
//                its numerical solution." International Journal for Numerical
//

//
// Subclassed by: José A. Abell (UANDES), Francisco Pinto (UChile), Ricardo Gallardo (PUCV)
//

#include "UnsaturatedSSPbrickUP.h"

#include <elementAPI.h>
#include <NDMaterial.h>
#include <ElementResponse.h>


static int num_UnsaturatedSSPbrickUP = 0;
void* OPS_UnsaturatedSSPbrickUP(void)
{
    if (num_UnsaturatedSSPbrickUP == 0) {
        num_UnsaturatedSSPbrickUP++;
        opserr << "UnsaturatedSSPbrickUP element - Written: José A. Abell (UANDES), Francisco Pinto (UChile), Ricardo Gallardo (PUCV)\n";
    }

    // Pointer to an element that will be returned
    Element *theElement = 0;

    int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingInputArgs < 25) {
        opserr << "Invalid #args, want: element UnsaturatedSSPbrickUP eleTag? iNode? jNode? kNode? lNode? mNode? nNode? oNode? pNode? matTag? Kf? Rf? k1? k2? k3? eVoid? alpha? Sres? Ssat? ga? gn? gc? gl? krel_min? <b1? b2? b3?>\n";
        return 0;
    }

    int iData[10];
    double dData[17];
    dData[13] = 1e-4; // krel_min?
    dData[14] = 0.0; // <b1?>
    dData[15] = 0.0; // <b2?>
    dData[16] = 0.0; // <b3?>

    int numData = 10; // number of data to read
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element UnsaturatedSSPbrickUP " << iData[0] << endln;
        return 0;
    }
    numRemainingInputArgs -= numData; // update remaining args

    numData = 14; // update number of data to read
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid double data: element UnsaturatedSSPbrickUP " << iData[0] << endln;
        return 0;
    }
    numRemainingInputArgs -= numData; // update remaining args

    int matID = iData[9];
    NDMaterial *theMaterial = OPS_GetNDMaterial(matID);
    if (theMaterial == 0) {
        opserr << "WARNING element UnsaturatedSSPbrickUP " << iData[0] << endln;
        opserr << " Material: " << matID << " not found\n";
        return 0;
    }

    if (numRemainingInputArgs >= 3) { // check if there are at least 3 more args
        numData = 3;
        if (OPS_GetDoubleInput(&numData, &dData[14]) != 0) {
            opserr << "WARNING invalid optional data: element UnsaturatedSSPbrickUP " << iData[0] << endln;
            return 0;
        }
        numRemainingInputArgs -= numData; // update remaining args
    }

    theElement = new UnsaturatedSSPbrickUP(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], iData[6], iData[7], iData[8], *theMaterial,
                                           dData[0], dData[1], dData[2], dData[3], dData[4],  dData[5], dData[6],
                                           dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
                                           dData[13], dData[14], dData[15]);

    if (theElement == 0) {
        opserr << "WARNING could not create element of type UnsaturatedSSPbrickUP\n";
        return 0;
    }

    return theElement;
}


// full constructor
UnsaturatedSSPbrickUP::UnsaturatedSSPbrickUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
        NDMaterial &theMat, double Kf, double Rf, double k1, double k2, double k3,
        double eVoid, double alpha,
        double Sres, double Ssat, double ga, double gn, double gc, double gl, double krel_min,
        double b1, double b2, double b3)
    : SSPbrickUP(tag,  Nd1,  Nd2,  Nd3,  Nd4,  Nd5,  Nd6,  Nd7,  Nd8,
                 theMat, Kf, Rf, k1, k2, k3,
                 eVoid, alpha, b1, b2, b3,
                 ELE_TAG_UnsaturatedSSPbrickUP),
      mSres(Sres), mSsat(Ssat), mga(ga), mgn(gn), mgc(gc), mgl(gl), mkrel_min(krel_min)
{

}

// null constructor
UnsaturatedSSPbrickUP::UnsaturatedSSPbrickUP()
    : SSPbrickUP(ELE_TAG_UnsaturatedSSPbrickUP),
      mSres(0.0), mSsat(0.0), mga(0.0), mgn(0.0), mgc(0.0), mgl(0.0), mkrel_min(0.0)
{
}




const Matrix &
UnsaturatedSSPbrickUP::getMass(void)
{
    mMass.Zero();

    //Compute fluid pressure as average of nodal values
    const double p0 = theNodes[0]->getTrialVel()(2);
    const double p1 = theNodes[1]->getTrialVel()(2);
    const double p2 = theNodes[2]->getTrialVel()(2);
    const double p3 = theNodes[3]->getTrialVel()(2);
    const double p4 = theNodes[4]->getTrialVel()(2);
    const double p5 = theNodes[5]->getTrialVel()(2);
    const double p6 = theNodes[6]->getTrialVel()(2);
    const double p7 = theNodes[7]->getTrialVel()(2);

    double pw = (p0 + p1 + p2 + p3 + p5 + p6 + p7) / 8;

    pw = pw < 0 ? pw : 0;

    // compute compressibilty matrix term S
    // double oneOverQ = -0.015625 * mVol * mPorosity / fBulk;
    double glocal = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    double psi = -pw / (fDens * glocal);
    double S = getRelativeSaturation(psi);
    double dSdpw = getRelativeSaturationPressureDerivative(pw);
    
    double oneOverQ = -0.015625 * mVol * (mPorosity / fBulk * S - mPorosity * dSdpw); // Eq 58 Plaxis

    // full mass matrix for the element [ M    0   ]
    //  includes M and S submatrices    [ 0 -(S+H) ]
    for (int i = 0; i < 8; i++) {

        int I    = 3 * i;
        int Ip1  = 3 * i + 1;
        int Ip2  = 3 * i + 2;
        int II   = 4 * i;
        int IIp1 = 4 * i + 1;
        int IIp2 = 4 * i + 2;
        int IIp3 = 4 * i + 3;

        for (int j = 0; j < 8; j++) {

            int J    = 3 * j;
            int Jp1  = 3 * j + 1;
            int Jp2  = 3 * j + 2;
            int JJ   = 4 * j;
            int JJp1 = 4 * j + 1;
            int JJp2 = 4 * j + 2;
            int JJp3 = 4 * j + 3;

            // contribution of solid phase mass
            mMass(II, JJ)     = mSolidM(I, J);
            mMass(IIp1, JJ)   = mSolidM(Ip1, J);
            mMass(IIp2, JJ)   = mSolidM(Ip2, J);
            mMass(IIp1, JJp1) = mSolidM(Ip1, Jp1);
            mMass(IIp2, JJp1) = mSolidM(Ip2, Jp1);
            mMass(IIp2, JJp2) = mSolidM(Ip2, Jp2);
            mMass(II, JJp1)   = mSolidM(I, Jp1);
            mMass(II, JJp2)   = mSolidM(I, Jp2);
            mMass(IIp1, JJp2) = mSolidM(Ip1, Jp2);

            // contribution of compressibility and stabilization matrices
            mMass(IIp3, JJp3) = oneOverQ - mPressStab(i, j);
        }
    }

    return mMass;
}


const Vector &
UnsaturatedSSPbrickUP::getResistingForce(void)
// this function computes the resisting force vector for the element
{
    Vector f1(24);
    f1.Zero();
    Vector f2(8);
    f2.Zero();

    // get stress from the material
    Vector mStress = theMaterial->getStress();

    // get trial displacement
    const Vector &mDisp_1 = theNodes[0]->getTrialDisp();
    const Vector &mDisp_2 = theNodes[1]->getTrialDisp();
    const Vector &mDisp_3 = theNodes[2]->getTrialDisp();
    const Vector &mDisp_4 = theNodes[3]->getTrialDisp();
    const Vector &mDisp_5 = theNodes[4]->getTrialDisp();
    const Vector &mDisp_6 = theNodes[5]->getTrialDisp();
    const Vector &mDisp_7 = theNodes[6]->getTrialDisp();
    const Vector &mDisp_8 = theNodes[7]->getTrialDisp();

    // assemble displacement vector
    Vector d(24);
    d(0) =  mDisp_1(0);
    d(1) =  mDisp_1(1);
    d(2) =  mDisp_1(2);
    d(3) =  mDisp_2(0);
    d(4) =  mDisp_2(1);
    d(5) =  mDisp_2(2);
    d(6) =  mDisp_3(0);
    d(7) =  mDisp_3(1);
    d(8) =  mDisp_3(2);
    d(9) =  mDisp_4(0);
    d(10) = mDisp_4(1);
    d(11) = mDisp_4(2);
    d(12) = mDisp_5(0);
    d(13) = mDisp_5(1);
    d(14) = mDisp_5(2);
    d(15) = mDisp_6(0);
    d(16) = mDisp_6(1);
    d(17) = mDisp_6(2);
    d(18) = mDisp_7(0);
    d(19) = mDisp_7(1);
    d(20) = mDisp_7(2);
    d(21) = mDisp_8(0);
    d(22) = mDisp_8(1);
    d(23) = mDisp_8(2);

    // add stabilization force to internal force vector
    f1 = Kstab * d;

    // add internal force from the stress  ->  fint = Kstab*d + 8*Jo*Bnot'*stress
    f1.addMatrixTransposeVector(1.0, Bnot, mStress, mVol);

    // subtract body forces from internal force vector
    double density = theMaterial->getRho();

    if (applyLoad == 0) {
        double polyJac = 0.0;
        for (int i = 0; i < 8; i++) {
            /*polyJac = J[0] + (J[1]*xi(i) + J[2]*et(i) + J[3]*ze(i) + J[7] + J[8] + J[9])/3.0
                     + (J[4]*hut(i) + J[5]*hus(i) + J[6]*hst(i) + J[10]*ze(i) + J[11]*et(i) + J[12]*xi(i) + J[13]*ze(i) + J[14]*et(i) + J[15]*xi(i))/9.0
                     + (J[16]*hstu(i) + J[17]*hut(i) + J[18]*hus(i) + J[19]*hst(i))/27.0;*/
            polyJac = J[0] * (1.0 + (J[1] * xi(i) + J[2] * et(i) + J[3] * ze(i) + J[7] + J[8] + J[9]) / 3.0
                              + (J[4] * hut(i) + J[5] * hus(i) + J[6] * hst(i) + J[10] * ze(i) + J[11] * et(i) + J[12] * xi(i) + J[13] * ze(i) + J[14] * et(i) + J[15] * xi(i)) / 9.0
                              + (J[16] * hstu(i) + J[17] * hut(i) + J[18] * hus(i) + J[19] * hst(i)) / 27.0);
            f1(3 * i)   -= density * b[0] * polyJac;
            f1(3 * i + 1) -= density * b[1] * polyJac;
            f1(3 * i + 2) -= density * b[2] * polyJac;
        }
    } else {
        double polyJac = 0.0;
        for (int i = 0; i < 8; i++) {
            /*polyJac = J[0] + (J[1]*xi(i) + J[2]*et(i) + J[3]*ze(i) + J[7] + J[8] + J[9])/3.0
                     + (J[4]*hut(i) + J[5]*hus(i) + J[6]*hst(i) + J[10]*ze(i) + J[11]*et(i) + J[12]*xi(i) + J[13]*ze(i) + J[14]*et(i) + J[15]*xi(i))/9.0
                     + (J[16]*hstu(i) + J[17]*hut(i) + J[18]*hus(i) + J[19]*hst(i))/27.0;*/
            polyJac = J[0] * (1.0 + (J[1] * xi(i) + J[2] * et(i) + J[3] * ze(i) + J[7] + J[8] + J[9]) / 3.0
                              + (J[4] * hut(i) + J[5] * hus(i) + J[6] * hst(i) + J[10] * ze(i) + J[11] * et(i) + J[12] * xi(i) + J[13] * ze(i) + J[14] * et(i) + J[15] * xi(i)) / 9.0
                              + (J[16] * hstu(i) + J[17] * hut(i) + J[18] * hus(i) + J[19] * hst(i)) / 27.0);
            f1(3 * i)   -= density * appliedB[0] * polyJac;
            f1(3 * i + 1) -= density * appliedB[1] * polyJac;
            f1(3 * i + 2) -= density * appliedB[2] * polyJac;
        }
    }

    // account for fluid body forces
    Matrix k(3, 3);
    Vector body(3);
    // permeability tensor
    k(0, 0) = perm[0];
    k(1, 1) = perm[1];
    k(2, 2) = perm[2];
    // body force vector
    if (applyLoad == 0) {
        body(0) = b[0];
        body(1) = b[1];
        body(2) = b[2];
    } else {
        body(0) = appliedB[0];
        body(1) = appliedB[1];
        body(2) = appliedB[2];
    }

    //Compute fluid pressure as average of nodal values
    const double p0 = theNodes[0]->getTrialVel()(2);
    const double p1 = theNodes[1]->getTrialVel()(2);
    const double p2 = theNodes[2]->getTrialVel()(2);
    const double p3 = theNodes[3]->getTrialVel()(2);
    const double p4 = theNodes[4]->getTrialVel()(2);
    const double p5 = theNodes[5]->getTrialVel()(2);
    const double p6 = theNodes[6]->getTrialVel()(2);
    const double p7 = theNodes[7]->getTrialVel()(2);

    double pw = (p0 + p1 + p2 + p3 + p5 + p6 + p7) / 8;

    pw = pw < 0 ? pw : 0;

    // compute compressibilty matrix term S
    // double oneOverQ = -0.015625 * mVol * mPorosity / fBulk;
    double glocal = sqrt(b[0] * b[0] + b[1] * b[1]);
    double psi = -pw / (fDens * glocal);
    double S = getRelativeSaturation(psi);
    double dSdpw = getRelativeSaturationPressureDerivative(pw);
    double krel = getRelativePermeability(S);

    f2 = mVol * fDens * dNmod * krel / (fDens * glocal) * k * body; // Eq. 59 Plaxis
    // f2 = mVol * fDens * dNmod * k * body;

    // assemble full internal force vector for the element
    mInternalForces(0)  = f1(0);
    mInternalForces(1)  = f1(1);
    mInternalForces(2)  = f1(2);
    mInternalForces(3)  = f2(0);
    mInternalForces(4)  = f1(3);
    mInternalForces(5)  = f1(4);
    mInternalForces(6)  = f1(5);
    mInternalForces(7)  = f2(1);
    mInternalForces(8)  = f1(6);
    mInternalForces(9)  = f1(7);
    mInternalForces(10) = f1(8);
    mInternalForces(11) = f2(2);
    mInternalForces(12) = f1(9);
    mInternalForces(13) = f1(10);
    mInternalForces(14) = f1(11);
    mInternalForces(15) = f2(3);
    mInternalForces(16) = f1(12);
    mInternalForces(17) = f1(13);
    mInternalForces(18) = f1(14);
    mInternalForces(19) = f2(4);
    mInternalForces(20) = f1(15);
    mInternalForces(21) = f1(16);
    mInternalForces(22) = f1(17);
    mInternalForces(23) = f2(5);
    mInternalForces(24) = f1(18);
    mInternalForces(25) = f1(19);
    mInternalForces(26) = f1(20);
    mInternalForces(27) = f2(6);
    mInternalForces(28) = f1(21);
    mInternalForces(29) = f1(22);
    mInternalForces(30) = f1(23);
    mInternalForces(31) = f2(7);

    // inertial unbalance load
    mInternalForces.addVector(1.0, Q, -1.0);

    return mInternalForces;
}

Response*
UnsaturatedSSPbrickUP::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    char outputData[32];

    output.tag("ElementOutput");
    output.attr("eleType", "OPS_UnsaturatedSSPbrickUP");
    output.attr("eleTag", this->getTag());
    for (int i = 1; i <= 8; i++)
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
        theResponse = SSPbrickUP::setResponse(argv, argc, output);
    }

    output.endTag(); // ElementOutput

    return theResponse;
}

int
UnsaturatedSSPbrickUP::getResponse(int responseID, Information &eleInfo)
{
    static Vector saturation_vec(1);

    if (responseID == 1010) {
        //Compute fluid pressure as average of nodal values
        const double p0 = theNodes[0]->getTrialVel()(2);
        const double p1 = theNodes[1]->getTrialVel()(2);
        const double p2 = theNodes[2]->getTrialVel()(2);
        const double p3 = theNodes[3]->getTrialVel()(2);
        const double p4 = theNodes[4]->getTrialVel()(2);
        const double p5 = theNodes[5]->getTrialVel()(2);
        const double p6 = theNodes[6]->getTrialVel()(2);
        const double p7 = theNodes[7]->getTrialVel()(2);

        double pw = (p0 + p1 + p2 + p3 + p5 + p6 + p7) / 8;

        // compute compressibility matrix term
        double glocal = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
        double psi = -pw / (fDens * glocal);
        double S = getRelativeSaturation(psi);
        saturation_vec(0) = S;
        return eleInfo.setVector(saturation_vec);
    }
    else
    {
        return SSPbrickUP::getResponse(responseID, eleInfo);
    }

    return -1;
}


void
UnsaturatedSSPbrickUP::GetPermeabilityMatrix(void)
// this function computes the permeability matrix for the element
{
    Matrix k(3, 3);
    double root3 = 8.0 / (sqrt(3.0));
    Vector s = root3 * xi;
    Vector t = root3 * et;
    Vector u = root3 * ze;
    Matrix dNloc(8, 3);
    Matrix Jmat(3, 3);
    Matrix Jinv(3, 3);
    Matrix dN(8, 3);
    Matrix dNT(3, 8);
    mPerm.Zero();
    mPressStab.Zero();

    // permeability tensor
    k(0, 0) = perm[0];
    k(1, 1) = perm[1];
    k(2, 2) = perm[2];

    for (int i = 0; i < 8; i++) {

        double dsN1 = -0.125 * (1 - t(i)) * (1 - u(i));
        double dsN2 =  0.125 * (1 - t(i)) * (1 - u(i));
        double dsN3 =  0.125 * (1 + t(i)) * (1 - u(i));
        double dsN4 = -0.125 * (1 + t(i)) * (1 - u(i));
        double dsN5 = -0.125 * (1 - t(i)) * (1 + u(i));
        double dsN6 =  0.125 * (1 - t(i)) * (1 + u(i));
        double dsN7 =  0.125 * (1 + t(i)) * (1 + u(i));
        double dsN8 = -0.125 * (1 + t(i)) * (1 + u(i));

        double dtN1 = -0.125 * (1 - s(i)) * (1 - u(i));
        double dtN2 = -0.125 * (1 + s(i)) * (1 - u(i));
        double dtN3 =  0.125 * (1 + s(i)) * (1 - u(i));
        double dtN4 =  0.125 * (1 - s(i)) * (1 - u(i));
        double dtN5 = -0.125 * (1 - s(i)) * (1 + u(i));
        double dtN6 = -0.125 * (1 + s(i)) * (1 + u(i));
        double dtN7 =  0.125 * (1 + s(i)) * (1 + u(i));
        double dtN8 =  0.125 * (1 - s(i)) * (1 + u(i));

        double duN1 = -0.125 * (1 - s(i)) * (1 - t(i));
        double duN2 = -0.125 * (1 + s(i)) * (1 - t(i));
        double duN3 = -0.125 * (1 + s(i)) * (1 + t(i));
        double duN4 = -0.125 * (1 - s(i)) * (1 + t(i));
        double duN5 =  0.125 * (1 - s(i)) * (1 - t(i));
        double duN6 =  0.125 * (1 + s(i)) * (1 - t(i));
        double duN7 =  0.125 * (1 + s(i)) * (1 + t(i));
        double duN8 =  0.125 * (1 - s(i)) * (1 + t(i));

        dNloc(0, 0) = dsN1; dNloc(0, 1) = dtN1; dNloc(0, 2) = duN1;
        dNloc(1, 0) = dsN2; dNloc(1, 1) = dtN2; dNloc(1, 2) = duN2;
        dNloc(2, 0) = dsN3; dNloc(2, 1) = dtN3; dNloc(2, 2) = duN3;
        dNloc(3, 0) = dsN4; dNloc(3, 1) = dtN4; dNloc(3, 2) = duN4;
        dNloc(4, 0) = dsN5; dNloc(4, 1) = dtN5; dNloc(4, 2) = duN5;
        dNloc(5, 0) = dsN6; dNloc(5, 1) = dtN6; dNloc(5, 2) = duN6;
        dNloc(6, 0) = dsN7; dNloc(6, 1) = dtN7; dNloc(6, 2) = duN7;
        dNloc(7, 0) = dsN8; dNloc(7, 1) = dtN8; dNloc(7, 2) = duN8;

        Jmat = mNodeCrd * dNloc;
        Jmat.Invert(Jinv);

        dN = dNloc * Jinv;

        dNT(0, 0) = dN(0, 0); dNT(0, 1) = dN(1, 0); dNT(0, 2) = dN(2, 0); dNT(0, 3) = dN(3, 0); dNT(0, 4) = dN(4, 0); dNT(0, 5) = dN(5, 0); dNT(0, 6) = dN(6, 0); dNT(0, 7) = dN(7, 0);
        dNT(1, 0) = dN(0, 1); dNT(1, 1) = dN(1, 1); dNT(1, 2) = dN(2, 1); dNT(1, 3) = dN(3, 1); dNT(1, 4) = dN(4, 1); dNT(1, 5) = dN(5, 1); dNT(1, 6) = dN(6, 1); dNT(1, 7) = dN(7, 1);
        dNT(2, 0) = dN(0, 2); dNT(2, 1) = dN(1, 2); dNT(2, 2) = dN(2, 2); dNT(2, 3) = dN(3, 2); dNT(2, 4) = dN(4, 2); dNT(2, 5) = dN(5, 2); dNT(2, 6) = dN(6, 2); dNT(2, 7) = dN(7, 2);

        double detJ =  Jmat(0, 0) * Jmat(1, 1) * Jmat(2, 2) + Jmat(0, 1) * Jmat(1, 2) * Jmat(2, 0) + Jmat(0, 2) * Jmat(1, 0) * Jmat(2, 1)
                       - Jmat(0, 2) * Jmat(1, 1) * Jmat(2, 0) - Jmat(0, 1) * Jmat(1, 0) * Jmat(2, 2) - Jmat(0, 0) * Jmat(1, 2) * Jmat(2, 1);


        //Compute fluid pressure as average of nodal values
        const double p0 = theNodes[0]->getTrialVel()(2);
        const double p1 = theNodes[1]->getTrialVel()(2);
        const double p2 = theNodes[2]->getTrialVel()(2);
        const double p3 = theNodes[3]->getTrialVel()(2);
        const double p4 = theNodes[4]->getTrialVel()(2);
        const double p5 = theNodes[5]->getTrialVel()(2);
        const double p6 = theNodes[6]->getTrialVel()(2);
        const double p7 = theNodes[7]->getTrialVel()(2);

        double pw = (p0 + p1 + p2 + p3 + p5 + p6 + p7) / 8;

        pw = pw < 0 ? pw : 0;

        // compute compressibilty matrix term S
        // double oneOverQ = -0.015625 * mVol * mPorosity / fBulk;
        double glocal = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
        double psi = -pw / (fDens * glocal);
        double S = getRelativeSaturation(psi);
        double krel = getRelativePermeability(S);

        mPerm.addMatrixTripleProduct(1.0, dNT, krel / (fDens * glocal)*k, detJ);
        // mPerm.addMatrixTripleProduct(1.0, dNT, k, detJ);
        mPressStab.addMatrixTransposeProduct(1.0, dNT, dNT, detJ * mAlpha);
    }

    return;
}
