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
         

#include <ThermalVolumetricLoadingPattern.h>
#include <GroundMotion.h>

#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <stdlib.h>
#include <Channel.h>
#include <ErrorHandler.h>

#include <string.h>
#include <stdlib.h>

ThermalVolumetricLoadingPattern::ThermalVolumetricLoadingPattern(int tag, double alpha_, std::string elements_filename_, std::string gausstemps_filename_)
  :LoadPattern(tag, PATTERN_TAG_ThermalVolumetricLoadingPattern),  
  alpha(alpha_),
  elements_filename(elements_filename_),
  gausstemps_filename(gausstemps_filename_),
  currentTime(0.0), parameterID(0)
{

  opserr << "Creating ThermalVolumetricLoadingPattern" << endln;
  opserr << " alpha               = " << alpha << endln;
  opserr << " elements_filename   = " << elements_filename.c_str() << endln;
  opserr << " gausstemps_filename = " << gausstemps_filename.c_str() << endln;


}


ThermalVolumetricLoadingPattern::~ThermalVolumetricLoadingPattern()
{

}



void 
ThermalVolumetricLoadingPattern::applyLoad(double time)
{
    // buscar la linea que cooresponde al time en gausstemps_filename

    // para cada elemento en elements_filename, aplicarle los cambios de temperatura
    // que vienen en gausstemps_filename

    // iterar los puntos de gauss de cada elemento y llamar setParameter con la opcion initNormalStrain
}
    
void 
ThermalVolumetricLoadingPattern::applyLoadSensitivity(double time)
{
 
}


bool
ThermalVolumetricLoadingPattern::addSP_Constraint(SP_Constraint *)
{
  opserr << "ThermalVolumetricLoadingPattern::addSP_Constraint() - cannot add SP_Constraint to EQ pattern\n";
  return false;
}

bool
ThermalVolumetricLoadingPattern::addNodalLoad(NodalLoad *)
{
  opserr << "ThermalVolumetricLoadingPattern::addNodalLoad() - cannot add NodalLoad to EQ pattern\n";  
  return false;
}

bool
ThermalVolumetricLoadingPattern::addElementalLoad(ElementalLoad *)
{
  opserr << "ThermalVolumetricLoadingPattern::addElementalLoad() - cannot add ElementalLoad to EQ pattern\n";    
  return false;
}




// AddingSensitivity:BEGIN ////////////////////////////////////
int
ThermalVolumetricLoadingPattern::setParameter(const char **argv, int argc, Parameter &param)
{

    return 0;

}

int
ThermalVolumetricLoadingPattern::updateParameter(int pparameterID, Information &info)
{
  return 0;
}

int
ThermalVolumetricLoadingPattern::activateParameter(int pparameterID)
{
  parameterID = pparameterID;

  return 0;
}
// AddingSensitivity:END ////////////////////////////////////
