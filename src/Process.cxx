#ifndef __Process_cxx__
#define __Process_cxx__

#include "JetHadronSkimmer.h"
#include "CentralityAnalysis.h"
#include "JetSubtractedEnergy.h"
#include "TrackingPerformance.h"
#include "TrackMomentumResolution.h"
#include "JetEnergyResolution.h"
#include "Params.h"
#include "CentralityDefs.h"
#include "Process.h"

#include <string>
#include <iostream>

using namespace std;
using namespace JetHadronCorrelations;


int main (int argc, char** argv) {

  int argn                = 1;
  const string alg        = std::string (argv[argn++]);
  const char* subdir      = argv[argn++];
  const int dataSet       = atoi (argv[argn++]);

  TString collSys = TString (argv[argn++]);
  if (!SetCollisionSystem (collSys)) {
    std::cout << "In Process.cxx: Invalid collision system, exiting." << std::endl;
    return 1;
  }

  TString dType = TString (argv[argn++]);
  if (!SetDataType (dType)) {
    std::cout << "In Process.cxx: Invalid data type, exiting." << std::endl;
    return 2;
  }

  TString tType = TString (argv[argn++]);
  if (!SetTriggerType (tType)) {
    std::cout << "In Process.cxx: Invalid trigger type, exiting." << std::endl;
    return 3;
  }

  TString sFlag = TString (argv[argn++]);
  if (!SetSystFlag (sFlag)) {
    std::cout << "In Process.cxx: Invalid systematic flag, exiting." << std::endl;
    return 4;
  }

  const char* inFileName            = (argc > argn && argv[argn] ? argv[argn++] : "");
  if (!GetMCWeights (inFileName))
    std::cout << "In Process.cxx: Failed to get event weights, please check." << std::endl;

  std::cout << "Info: In Process.cxx: Configuration set to";
  std::cout << "\n\talg = " << alg;
  std::cout << "\n\tsubdir = " << subdir;
  std::cout << "\n\tdataSet = " << dataSet;
  std::cout << "\n\tCollisionSystem = " << ToTString (collisionSystem);
  std::cout << "\n\tDataType = " << ToTString (dataType);
  std::cout << "\n\tTriggerType = " << ToTString (triggerType);
  std::cout << "\n\tSystFlag = " << ToTString (systFlag);
  std::cout << "\n\tinFileName = " << inFileName;
  if (crossSectionPicoBarns != 0.)
    std::cout << "\n\t  --> deduced crossSectionPicoBarns = " << crossSectionPicoBarns;
  if (mcFilterEfficiency != 0.)
    std::cout << "\n\t  --> deduced mcFilterEfficiency = " << mcFilterEfficiency;
  if (mcNumberEvents != 0.)
    std::cout << "\n\t  --> deduced mcNumberEvents = " << mcNumberEvents;
  std::cout << std::endl;


  bool success = false;
  if (alg == "JetHadronSkimmer") {
    std::cout << "Info: In Process.cxx: Running JetHadronSkimmer algorithm..." << std::endl;
    success = JetHadronCorrelations::JetHadronSkimmer (subdir, dataSet, inFileName);
  }
  else if (alg == "CentralityAnalysis") {
    std::cout << "Info: In Process.cxx: Running CentralityAnalysis algorithm..." << std::endl;
    success = JetHadronCorrelations::CentralityAnalysis (subdir, dataSet, inFileName);
  }
  else if (alg == "JetSubtractedEnergy") {
    std::cout << "Info: In Process.cxx: Running JetSubtractedEnergy algorithm..." << std::endl;
    success = JetHadronCorrelations::JetSubtractedEnergy (subdir, dataSet, inFileName);
  }
  else if (alg == "TrackingPerformance") {
    std::cout << "Info: In Process.cxx: Running TrackingPerformance algorithm..." << std::endl;
    success = JetHadronCorrelations::TrackingPerformance (subdir, dataSet, inFileName);
  }
  else if (alg == "TrackMomentumResolution") {
    std::cout << "Info: In Process.cxx: Running TrackMomentumResolution algorithm..." << std::endl;
    success = JetHadronCorrelations::TrackMomentumResolution (subdir, dataSet, inFileName);
  }
  else if (alg == "JetEnergyResolution") {
    std::cout << "Info: In Process.cxx: Running JetEnergyResolution algorithm..." << std::endl;
    success = JetHadronCorrelations::JetEnergyResolution (subdir, dataSet, inFileName);
  }
  else {
    std::cout << "Error: In Process.cxx: Failed to recognize algorithm! Quitting." << std::endl;
    return 1;
  }

  delete [] zdcCentBins;
  zdcCentBins = nullptr;
  delete [] fcalCentBins;
  fcalCentBins = nullptr;
  delete [] fcalMixBins;
  fcalMixBins = nullptr;
  delete [] fineFcalCentBins;
  fineFcalCentBins = nullptr;
  delete [] ppMixBins;
  ppMixBins = nullptr;

  if (success) {
    std::cout << "Info: In Process.cxx: Finished algorithm!" << std::endl;
    return 0;
  }
  else {
    std::cout << "Error: In Process.cxx: Algorithm failed!" << std::endl;
    return 1;
  }
}

#endif
