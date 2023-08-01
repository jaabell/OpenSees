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
                                                                       
#ifndef UnsaturatedSSPbrickUP_h
#define UnsaturatedSSPbrickUP_h

//
// Subclassed by: José A. Abell (UANDES), Francisco Pinto (UChile), Ricardo Gallardo (PUCV)
//

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <SSPbrickUP.h>
#include <Response.h>

#include <cmath>

using namespace std;

class UnsaturatedSSPbrickUP : public SSPbrickUP
{
  public:
    UnsaturatedSSPbrickUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
                      NDMaterial &theMat, double Kf, double Rf, double k1, double k2, double k3,
					  double eVoid, double alpha, double Sres, double Ssat, double ga, double gn, double gc, double gl, double krel_min, double b1 = 0.0, double b2 = 0.0, double b3 = 0.0);
    UnsaturatedSSPbrickUP();

	const char* getClassType()  const { return "UnsaturatedSSPbrickUP"; };

    Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
    int getResponse(int responseID, Information &eleInformation);
    const Matrix &getMass(void);
    const Vector &getResistingForce(void); 
 

  protected:

    void GetPermeabilityMatrix(void);


    // Van Genuchten (1980) equations

    inline double getRelativeSaturation(double psi) const 
    {
        return mSres + (mSsat - mSres)*pow(1 + pow(mga*abs(psi), mgn), mgc);
    }

    inline double getRelativeSaturationPressureDerivative(double pw) const
    {   

        return (mSsat - mSres) * 
            mgc * 
            ( mgn * pow( mga / fDens, mgn) * pow(pw, mgn-1) ) * 
            pow(1 + pow(mga*abs(pw / fDens), mgn), mgc-1);
    }

    inline double getRelativePermeability(double S) const 
    {
        double Seff = (S - mSres)/(mSsat - mSres);
        return max(mkrel_min, 
            pow(Seff, mgl) * 
            pow(1 - pow(1 - pow(Seff, -1/mgc), -mgc),2)
            );
    }


    double mSres;           //A residual saturation, part of the fluid that remains in the pores even at high suction heads
    double mSsat;           //Saturated condition (Ssat=1, fully saturated which in general is not true)
    double mga;             //Van Genuchten fitting parameter related to air entry value of the soil, specific for material, units: 1/L
    double mgn;             //Van Genuchten fitting parameter related to the rate of water extraction from soil once air entry value is exceeded, specific for material, units: 1/L
    double mgc;             //Van Genuchten fitting parameter, specific for material, in absence of data can use gc = (1-gn)/gn
    double mgl;             //Van Genuchten fitting parameter for the relative permeability
    double mkrel_min;       //Minimum relative permeability (default 1e-4)


  private:
};

#endif
