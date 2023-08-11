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


/**
 * ThermalVolumetricLoadingPattern constructor.
 * Initializes member variables and reads the element tags from a file.
 *
 * @param tag Integer that represents the object tag
 * @param alpha_ Double, thermal expansion coefficient
 * @param c1_ Double, first coefficient for logarithmic function
 * @param c2_ Double, second coefficient for logarithmic function
 * @param K_ Double, bulk modulus
 * @param elements_filename_ String, name of the file containing element tags
 * @param gausstemps_filename_ String, name of the file containing Gauss temperature data
 */
ThermalVolumetricLoadingPattern::ThermalVolumetricLoadingPattern(int tag, double alpha_, double c1_, double c2_, double K_, std::string elements_filename_, std::string gausstemps_filename_)
  :LoadPattern(tag, PATTERN_TAG_ThermalVolumetricLoadingPattern),  
  alpha(alpha_),
  c1(c1_),
  c2(c2_),
  K(K_),
  elements_filename(elements_filename_),
  gausstemps_filename(gausstemps_filename_),
  currentTime(0.0), parameterID(0)
{

  opserr << "Creating ThermalVolumetricLoadingPattern" << endln;
  opserr << " alpha               = " << alpha << endln;
  opserr << " c1                  = " << c1 << endln;
  opserr << " c2                  = " << c2 << endln;
  opserr << " K                   = " << K << endln;
  opserr << " elements_filename   = " << elements_filename.c_str() << endln;
  opserr << " gausstemps_filename = " << gausstemps_filename.c_str() << endln;

  std::ifstream file(elements_filename);
  if (!file) {
      std::cerr << "Could not open: " << elements_filename.c_str() << endln;
      return ;
  }

  std::vector<int> elementTags_((std::istream_iterator<int>(file)),
               std::istream_iterator<int>());

  elementTags = elementTags_ ;
}


ThermalVolumetricLoadingPattern::~ThermalVolumetricLoadingPattern()
{

}

void 
ThermalVolumetricLoadingPattern::setDomain(Domain *theDomain)
{
  this->LoadPattern::setDomain(theDomain);
}

/**
 * Applies the load at a given time by updating strain and material properties.
 *
 * @param time Double, the time at which the load should be applied
 */
void
ThermalVolumetricLoadingPattern::applyLoad(double time)
{
    // Calculate Young's modulus and Poisson's ratio at the given time
    // double E = c2;

    // if (time >= 1)
    //     E = c1 * log(time) + c2;

    // double nu = (3 * K - E) / (6 * K);
    std::vector<double> initTemp;

    // Open the file containing Gauss temperature data
    std::ifstream gausstemps_file(gausstemps_filename);
    std::string line;

    std::vector<double> gaussDataEarlier, gaussDataLater;
    double earlierTime = 0, laterTime = 0;
    bool found_interval = false;
    
    // If the file is open and there's a line to read
    if (gausstemps_file && getline(gausstemps_file, line)) {
        std::stringstream ss(line);
        double file_time;
        ss >> file_time;
        earlierTime = file_time;
        
        // Now we read the first temperature data into initTemp
        double value;
        while (ss >> value) {
            initTemp.push_back(value);
        }
    }
    
    // Read the file to find the temperature data for the time interval
    while (getline(gausstemps_file, line)) {
        std::stringstream ss(line);
        double file_time;
        ss >> file_time;

        // Find the appropriate time interval in the file
        if (file_time >= time) {
            laterTime = file_time;
            double value;
            while (ss >> value) {
                gaussDataLater.push_back(value);
            }
            found_interval = true;
            break;
        }

        earlierTime = file_time;
        gaussDataEarlier.clear();
        double value;
        while (ss >> value) {
            gaussDataEarlier.push_back(value);
        }
    }

    // Close the file after reading
    gausstemps_file.close();

    // If time interval is not found, log the error and return
    if (!found_interval) {
        opserr << "Time interval not found." << endln;
        return;
    }

    double deltaEpsilon = 0.;
    int elementIndex = 0;

    // Loop through the element tags
    for (int eleTag : elementTags) {
        Element* theElement = this->getDomain()->getElement(eleTag);

        // Loop through the Gauss points (assumed 4 in this case)
        for (int gp = 1; gp <= 4; gp++) {
            int dataIndex = 4 * elementIndex + gp;
            if (dataIndex < gaussDataEarlier.size() && dataIndex < gaussDataLater.size()) {
                // Interpolate the temperature change
                double t = (time - earlierTime) / (laterTime - earlierTime);
                double tempChange = (1 - t) * gaussDataEarlier[dataIndex] + t * gaussDataLater[dataIndex] - initTemp[dataIndex];

                // Calculate the strain
                deltaEpsilon = alpha * tempChange;

                // Update the initial normal strain
                const char* argv[3] = {"material", std::to_string(gp).c_str(), "initNormalStrain"};
                int argc = 3;
                Parameter param(0, theElement, argv, argc);
                param.update(deltaEpsilon);
            }

            // {
            //     // Update Young's modulus
            //     const char* argv[3] = {"material", std::to_string(gp).c_str(), "E"};
            //     int argc = 3;
            //     Parameter param(0, theElement, argv, argc);
            //     param.update(E);
            // }

            // {
            //     // Update Poisson's ratio
            //     const char* argv[3] = {"material", std::to_string(gp).c_str(), "v"};
            //     int argc = 3;
            //     Parameter param(0, theElement, argv, argc);
            //     param.update(nu);
            // }
        }
        elementIndex++;
    }
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
