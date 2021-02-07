#include "JetHadronSkimmer.h"
#include "ZDCAnalysis.h"
#include "Params.h"
#include "run.h"

#include <string>
#include <iostream>

using namespace std;
using namespace HadronYieldsAnalysis;


int main (int argc, char** argv) {

  int argn                = 1;
  const string alg        = string (argv[argn++]);
  const char* subdir      = argv[argn++];
  const int dataSet       = atoi (argv[argn++]);

  TString collSys = TString (argv[argn++]);
  if      (collSys == ToTString (CollisionSystem::pp15))        { collisionSystem = CollisionSystem::pp15;       }
  else if (collSys == ToTString (CollisionSystem::PbPb15))      { collisionSystem = CollisionSystem::PbPb15;     }
  else if (collSys == ToTString (CollisionSystem::pPb16s5TeV))  { collisionSystem = CollisionSystem::pPb16s5TeV; }
  else if (collSys == ToTString (CollisionSystem::pPb16))       { collisionSystem = CollisionSystem::pPb16;      }
  else if (collSys == ToTString (CollisionSystem::Pbp16))       { collisionSystem = CollisionSystem::Pbp16;      }
  else if (collSys == ToTString (CollisionSystem::XeXe17))      { collisionSystem = CollisionSystem::XeXe17;     }
  else if (collSys == ToTString (CollisionSystem::pp17))        { collisionSystem = CollisionSystem::pp17;       }
  else if (collSys == ToTString (CollisionSystem::PbPb18))      { collisionSystem = CollisionSystem::PbPb18;     }
  else {
    std::cout << "In Run.cxx: Invalid collision system, exiting." << std::endl;
    return 1;
  }

  TString dType = TString (argv[argn++]);
  if      (dType == ToTString (DataType::Collisions))       { dataType = DataType::Collisions;      }
  else if (dType == ToTString (DataType::MCSignal))         { dataType = DataType::MCSignal;        }
  else if (dType == ToTString (DataType::MCDataOverlay))    { dataType = DataType::MCDataOverlay;   }
  else if (dType == ToTString (DataType::MCHijing))         { dataType = DataType::MCHijing;        }
  else if (dType == ToTString (DataType::MCHijingOverlay))  { dataType = DataType::MCHijingOverlay; }
  else {
    std::cout << "In Run.cxx: Invalid data type, exiting." << std::endl;
    return 2;
  }

  TString tType = TString (argv[argn++]);
  if      (tType == ToTString (TriggerType::Jet))     { triggerType = TriggerType::Jet;     }
  else if (tType == ToTString (TriggerType::MinBias)) { triggerType = TriggerType::MinBias; }
  else if (IsCollisions ()) {
    std::cout << "In Run.cxx: Invalid trigger type in data, exiting." << std::endl;
    return 3;
  }

  doHITightVar            = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doPionsOnlyVar          = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);
  doWithPileupVar         = (argc > argn && argv[argn] ? string (argv[argn++]) == "true" : false);

  crossSectionPicoBarns             = (argc > argn && argv[argn] ? atof (argv[argn++]) : 0);
  mcFilterEfficiency                = (argc > argn && argv[argn] ? atof (argv[argn++]) : 0);
  mcNumberEvents                    = (argc > argn && argv[argn] ? atoi (argv[argn++]) : 0);
  const char* inFileName            = (argc > argn && argv[argn] ? argv[argn++] : "");
  const char* eventWeightsFileName  = (argc > argn && argv[argn] ? argv[argn++] : "");

  std::cout << "Info: In run.cxx: Configuration set to";
  std::cout << "\n\talg = " << alg;
  std::cout << "\n\tsubdir = " << subdir;
  std::cout << "\n\tdataSet = " << dataSet;
  std::cout << "\n\tCollisionSystem = " << ToTString (collisionSystem);
  std::cout << "\n\tDataType = " << ToTString (dataType);
  std::cout << "\n\tTriggerType = " << ToTString (triggerType);
  std::cout << "\n\tinFileName = " << inFileName;
  std::cout << "\n\tdoHITightVar = " << doHITightVar;
  std::cout << "\n\tdoPionsOnlyVar = " << doPionsOnlyVar;
  std::cout << "\n\tdoWithPileupVar = " << doWithPileupVar;
  if (crossSectionPicoBarns != 0.)
    std::cout << "\n\tcrossSectionPicoBarns = " << crossSectionPicoBarns;
  if (mcFilterEfficiency != 0.)
    std::cout << "\n\tmcFilterEfficiency = " << mcFilterEfficiency;
  if (mcNumberEvents != 0.)
    std::cout << "\n\tmcNumberEvents = " << mcNumberEvents;
  if (string (eventWeightsFileName) != "")
    std::cout << "\n\teventWeightsFileName = " << eventWeightsFileName;
  std::cout << std::endl;


  bool success = false;
  if (alg == "JetHadronSkimmer") {
    std::cout << "Info: In run.cxx: Running JetHadronSkimmer algorithm..." << std::endl;
    success = HadronYieldsAnalysis::JetHadronSkimmer (subdir, dataSet, inFileName);
  }
  else if (alg == "ZDCAnalysis") {
    std::cout << "Info: In run.cxx: Running ZDCAnalysis algorithm..." << std::endl;
    success = HadronYieldsAnalysis::ZDCAnalysis (subdir, dataSet, inFileName);
  }
  else {
    std::cout << "Error: In run.cxx: Failed to recognize algorithm! Quitting." << std::endl;
    return 1;
  }


  if (success) {
    std::cout << "Info: In run.cxx: Finished algorithm!" << std::endl;
    return 0;
  }
  else {
    std::cout << "Error: In run.cxx: Algorithm failed!" << std::endl;
    return 1;
  }
}
