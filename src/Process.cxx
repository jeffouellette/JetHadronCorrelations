#include "JetHadronSkimmer.h"
#include "ZDCAnalysis.h"
#include "Params.h"
#include "Process.h"

#include <string>
#include <iostream>

using namespace std;
using namespace JetHadronCorrelations;


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
    std::cout << "In Process.cxx: Invalid collision system, exiting." << std::endl;
    return 1;
  }

  TString dType = TString (argv[argn++]);
  if      (dType == ToTString (DataType::Collisions))       { dataType = DataType::Collisions;      }
  else if (dType == ToTString (DataType::MCSignal))         { dataType = DataType::MCSignal;        }
  else if (dType == ToTString (DataType::MCDataOverlay))    { dataType = DataType::MCDataOverlay;   }
  else if (dType == ToTString (DataType::MCHijing))         { dataType = DataType::MCHijing;        }
  else if (dType == ToTString (DataType::MCHijingOverlay))  { dataType = DataType::MCHijingOverlay; }
  else {
    std::cout << "In Process.cxx: Invalid data type, exiting." << std::endl;
    return 2;
  }

  TString tType = TString (argv[argn++]);
  if      (tType == ToTString (TriggerType::Jet))     { triggerType = TriggerType::Jet;     }
  else if (tType == ToTString (TriggerType::MinBias)) { triggerType = TriggerType::MinBias; }
  else if (IsCollisions ()) {
    std::cout << "In Process.cxx: Invalid trigger type in data, exiting." << std::endl;
    return 3;
  }

  TString sFlag = TString (argv[argn++]);
  if      (sFlag == ToTString (SystFlag::HITightVar))         { ToggleSyst (SystFlag::HITightVar);        }
  else if (sFlag == ToTString (SystFlag::PionsOnlyVar))       { ToggleSyst (SystFlag::PionsOnlyVar);      }
  else if (sFlag == ToTString (SystFlag::WithPileupVar))      { ToggleSyst (SystFlag::WithPileupVar);     }
  else if (sFlag == ToTString (SystFlag::JetES5PercUpVar))    { ToggleSyst (SystFlag::JetES5PercUpVar);   }
  else if (sFlag == ToTString (SystFlag::JetES5PercDownVar))  { ToggleSyst (SystFlag::JetES5PercDownVar); }
  else if (sFlag == ToTString (SystFlag::JetES2PercUpVar))    { ToggleSyst (SystFlag::JetES2PercUpVar);   }
  else if (sFlag == ToTString (SystFlag::JetES2PercDownVar))  { ToggleSyst (SystFlag::JetES2PercDownVar); }
  else {
    std::cout << "In Process.cxx: No systematics specified." << std::endl;
  }

  crossSectionPicoBarns             = (argc > argn && argv[argn] ? atof (argv[argn++]) : 0);
  mcFilterEfficiency                = (argc > argn && argv[argn] ? atof (argv[argn++]) : 0);
  mcNumberEvents                    = (argc > argn && argv[argn] ? atoi (argv[argn++]) : 0);
  const char* inFileName            = (argc > argn && argv[argn] ? argv[argn++] : "");
  const char* eventWeightsFileName  = (argc > argn && argv[argn] ? argv[argn++] : "");

  std::cout << "Info: In Process.cxx: Configuration set to";
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
  std::cout << "\n\tdoJetES5PercUpVar = " << doJetES5PercUpVar;
  std::cout << "\n\tdoJetES5PercDownVar = " << doJetES5PercDownVar;
  std::cout << "\n\tdoJetES2PercUpVar = " << doJetES2PercUpVar;
  std::cout << "\n\tdoJetES2PercDownVar = " << doJetES2PercDownVar;
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
    std::cout << "Info: In Process.cxx: Running JetHadronSkimmer algorithm..." << std::endl;
    success = JetHadronCorrelations::JetHadronSkimmer (subdir, dataSet, inFileName);
  }
  else if (alg == "ZDCAnalysis" && !Ispp ()) {
    std::cout << "Info: In Process.cxx: Running ZDCAnalysis algorithm..." << std::endl;
    success = JetHadronCorrelations::ZDCAnalysis (subdir, dataSet, inFileName);
  }
  else {
    std::cout << "Error: In Process.cxx: Failed to recognize algorithm! Quitting." << std::endl;
    return 1;
  }


  if (success) {
    std::cout << "Info: In Process.cxx: Finished algorithm!" << std::endl;
    return 0;
  }
  else {
    std::cout << "Error: In Process.cxx: Algorithm failed!" << std::endl;
    return 1;
  }
}
