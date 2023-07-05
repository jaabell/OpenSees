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
                                                                        
// ============================================================================
// 2023 By Jose Abell and Jose Larenas @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// 
// ============================================================================
                                          
                                                                        
#ifndef ThermalVolumetricLoadingPattern_h
#define ThermalVolumetricLoadingPattern_h



#include <LoadPattern.h>
#include <string>
#include <vector>
#include <fstream>
#include <iterator>
#include <iostream>



class GroundMotion;
class Vector;

class ThermalVolumetricLoadingPattern : public LoadPattern
{
  public:
    ThermalVolumetricLoadingPattern(int tag, double alpha_, double c1_, double c2_, double K_, std::string elements_filename_, std::string gausstemps_filename_);
    virtual ~ThermalVolumetricLoadingPattern();

    void setDomain(Domain *theDomain);

    virtual void applyLoad(double time);
    virtual bool addSP_Constraint(SP_Constraint *);
    virtual bool addNodalLoad(NodalLoad *);
    virtual bool addElementalLoad(ElementalLoad *);

    // methods for o/p
    virtual int sendSelf(int commitTag, Channel &theChannel) {return -1;};
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker) {return -1;};
    virtual void Print(OPS_Stream &s, int flag =0){s << "ThermalVolumetricLoadingPattern" << endln;};        

    // method to obtain a blank copy of the LoadPattern
    virtual LoadPattern *getCopy(void) {return 0;};;
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual void applyLoadSensitivity(double pseudoTime = 0.0);
    virtual int setParameter(const char **argv, int argc, Parameter &param);
    virtual int  updateParameter(int parameterID, Information &info);
    virtual int  activateParameter(int parameterID);
    // AddingSensitivity:END ///////////////////////////////////////////
    
 protected:

  private:
    double alpha;
    double c1, c2, K;
    std::string elements_filename; 
    std::string gausstemps_filename;
    double currentTime;

    std::vector<int> elementTags ;

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
