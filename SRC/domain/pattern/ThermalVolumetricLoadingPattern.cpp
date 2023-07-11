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
#include <Parameter.h>

#include <math.h>

#include <string>
#include <stdlib.h>
#include <sstream>

ThermalVolumetricLoadingPattern::ThermalVolumetricLoadingPattern(int tag, double alpha_, double c1_, double c2_, double K_, std::string elements_filename_, std::string gausstemps_filename_)
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


  // DEJAR BONITO (NOMBRES, STD::) 
  // Pasarselo a ChatGPT para que lo deje bien bonito y comentado


  std::ifstream file(elements_filename);
  if (!file) {
      std::cerr << "Could not open: " << elements_filename.c_str() << endln;
      return ;
  }

  std::vector<int> elementTags_((std::istream_iterator<int>(file)),
               std::istream_iterator<int>());

  elementTags = elementTags_ ;

  c1 = c1_ ; 
  c2 = c2_ ; 
  K = K_ ; 

  // for (int eleTag : elementTags) {
  //     std::cout << eleTag << ' ';
  // }

}


ThermalVolumetricLoadingPattern::~ThermalVolumetricLoadingPattern()
{

}

void 
ThermalVolumetricLoadingPattern::setDomain(Domain *theDomain)
{
  this->LoadPattern::setDomain(theDomain);
}

void 
ThermalVolumetricLoadingPattern::applyLoad(double time)
{
    // buscar la linea que cooresponde al time en gausstemps_filename

    // para cada elemento en elements_filename, aplicarle los cambios de temperatura
    // que vienen en gausstemps_filename

    // iterar los puntos de gauss de cada elemento y llamar setParameter con la opcion initNormalStrain

    // ------------------------------------------------------

    // para completar la funcion:
    // hacerlo para c/ punto de gauss
    // interpolar para distintos steps
    // modificar E y poisson (casi listo)

    double E  = c1 * log(time) + c2 ;
    double nu = (3 * K - E) / (6 * K) ;

    std::ifstream gausstemps_file(gausstemps_filename);
    std::string line;

    bool found_time = false;

    // Buscar la línea que corresponde al time en gausstemps_filename
    while (getline(gausstemps_file, line)) {
        double file_time;
        std::stringstream ss(line);
        ss >> file_time;

        if (abs(file_time - time) < 1e-3) {
            opserr << file_time ;
            found_time = true ;
            break;
        }
    }

    if (found_time)
      opserr << "Found time" ;
    else
      opserr << "Not found time" ;

    gausstemps_file.close() ;

    double tempChange = 1. ;
    double deltaEpsilon = 0. ;

    for (int eleTag : elementTags) {
      // std::cout << eleTag << ' ';

      Element *theElement = this->getDomain()->getElement(eleTag) ;


      for (int gp = 1; gp <= 4; gp++)
      {
        deltaEpsilon = alpha * tempChange ; // agarrar el tempchange de cada punto de gauss de cada elemento (restar la temp inicial con este)
        // siempre c/r al inicial no al anterior.
        // interpolar para distintos incrementos de tiempo (siempre sera la misma cantidad de tiempo, alcantidad de pasos no necesariamente)
        {
          const char *argv[3] = {"material", std::to_string(gp).c_str(), "initNormalStrain"};
          int argc = 3;
  
          Parameter param(0, theElement, argv, argc);
          param.update(deltaEpsilon);
        }

        {
          const char *argv[3] = {"material", std::to_string(gp).c_str(), "E"};
          int argc = 3;
          Parameter param(0, theElement, argv, argc);
          param.update(E);
        }

        {
          const char *argv[3] = {"material", std::to_string(gp).c_str(), "v"};
          int argc = 3;
          Parameter param(0, theElement, argv, argc);
          param.update(nu);
        }
      }
    }


    // // Para cada elemento en elements_filename, aplicarle los cambios de temperatura
    // // que vienen en gausstemps_filename
    // std::ifstream elements_file(elements_filename);
    // while (getline(elements_file, line)) {
    //     // Aquí necesitas alguna forma de obtener el elemento correspondiente a esta línea
    //     // Por ejemplo, si la línea es un ID de elemento, necesitarías una función para obtener el elemento a partir de su ID
    //     Element* element = get_element_from_line(line);

    //     // Aquí necesitas alguna forma de obtener la temperatura correspondiente a este elemento a partir de la línea de gausstemps_filename
    //     // Por ejemplo, si la temperatura está en la línea después del tiempo, podrías utilizar una función para obtener esa temperatura
    //     double temperature = get_temperature_from_gausstemps_line(line);

    //     element->setTemperature(temperature);
    // }

    // elements_file.close();

    // // Iterar los puntos de Gauss de cada elemento y llamar setParameter con la opción initNormalStrain
    // // Aquí estás suponiendo que tienes una lista de todos los elementos y que cada elemento tiene una lista de puntos de Gauss
    // for (Element* element : elements) {
    //     for (GaussPoint* gauss_point : element->getGaussPoints()) {
    //         const char* argv[] = {"initNormalStrain"};
    //         int argc = 1;
    //         Parameter param;

    //         gauss_point->setParameter(argv, argc, param);
    //     }
    // }

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
