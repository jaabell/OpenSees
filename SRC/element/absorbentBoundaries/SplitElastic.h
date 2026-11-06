// SPDX-License-Identifier: BSD-3-Clause
/*
Copyright (c) 2026, Jose Antonio Abell. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef SplitElastic_h
#define SplitElastic_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <array>

// Two half-property elastic continua in one quad/brick. A node with ndm DOFs
// is shared; a node with 2*ndm DOFs carries independent branch translations.
class SplitElastic : public Element {
public:
    SplitElastic(); // object broker
    SplitElastic(int tag, int ndm, const ID &nodes, double E, double nu,
                 double rho, double thickness, bool planeStress, bool lumped,
                 double coupling, const Vector &factors);
    const char *getClassType() const override { return "SplitElastic"; }
    int getNumExternalNodes() const override { return nodes.Size(); }
    const ID &getExternalNodes() override { return nodes; }
    Node **getNodePtrs() override { return nodePtrs.data(); }
    int getNumDOF() override { return numDOF; }
    void setDomain(Domain *domain) override;
    bool isValid() const { return valid; }
    int update() override { return valid ? 0 : -1; }
    int commitState() override { return valid ? 0 : -1; }
    int revertToLastCommit() override { return 0; }
    int revertToStart() override { return 0; }
    const Matrix &getTangentStiff() override;
    const Matrix &getInitialStiff() override { return getTangentStiff(); }
    const Matrix &getMass() override;
    const Matrix &getDamp() override;
    int setRayleighDampingFactors(double a, double b, double b0, double bc) override;
    void zeroLoad() override { load.Zero(); }
    int addLoad(ElementalLoad *, double) override;
    int addInertiaLoadToUnbalance(const Vector &accel) override;
    const Vector &getResistingForce() override;
    const Vector &getResistingForceIncInertia() override;
    int sendSelf(int commitTag, Channel &channel) override;
    int recvSelf(int commitTag, Channel &channel, FEM_ObjectBroker &broker) override;
    void Print(OPS_Stream &stream, int flag = 0) override;
    Response *setResponse(const char **argv, int argc, OPS_Stream &output) override;
    int getResponse(int responseID, Information &info) override;

private:
    bool integrate();
    bool shape(int gp, double N[8], double grad[8][3], double &weight) const;
    void constitutive(double D[6][6]) const;
    void strainMatrix(const double grad[8][3], double B[6][24]) const;
    void assemble(Matrix &out, double kScale, double mScale, bool relative) const;
    void force(bool inertia);
    void gather(int branch, int quantity, double *out) const;
    void scatter(int branch, const double *in, double scale, Vector &out) const;
    void massAction(const double *in, double *out) const;
    double relativePower() const;
    bool propertiesValid() const;

    int ndm = 0, numDOF = 0;
    ID nodes;
    std::array<Node *, 8> nodePtrs{};
    std::array<int, 24> mapA{}, mapB{};
    double E = 0., nu = 0., rho = 0., thickness = 1., coupling = 0.;
    bool planeStress = false, lumped = true, valid = false, allShared = false;
    Vector factors, volumes;
    Matrix stiffness, scalarMass;
    Vector resisting, load;
};
#endif
