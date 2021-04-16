#ifndef __LocalUtilities_cxx__
#define __LocalUtilities_cxx__

#include "LocalUtilities.h"
#include "TreeVariables.h"
#include "Params.h"
#include "Process.h"
#include "RunCorrelator.h"

#include <AtlasUtils.h>

#include <Utilities.h>

#include <TSystemDirectory.h>

#include <math.h>
#include <iostream>
#include <sys/stat.h>


namespace JetHadronCorrelations {


TString ToTString (const CollisionSystem collSys) {
  switch (collSys) {
    case CollisionSystem::pp15:       return TString ("pp15_5TeV");
    case CollisionSystem::PbPb15:     return TString ("PbPb15_5TeV");
    case CollisionSystem::pPb16s5TeV: return TString ("pPb16_5TeV");
    case CollisionSystem::pPb16:      return TString ("pPb16_8TeV");
    case CollisionSystem::Pbp16:      return TString ("Pbp16_8TeV");
    case CollisionSystem::XeXe17:     return TString ("XeXe17_5TeV");
    case CollisionSystem::pp17:       return TString ("pp17_5TeV");
    case CollisionSystem::PbPb18:     return TString ("PbPb18_5TeV");
    default:                          return TString ("???");
  }
}
 
 
 
TString ToTString (const DataType dType) {
  switch (dType) {
    case DataType::Collisions:      return TString ("Collisions");
    case DataType::MCSignal:        return TString ("MCSignal");
    case DataType::MCDataOverlay:   return TString ("MCDataOverlay");
    case DataType::MCHijing:        return TString ("MCHijing");
    case DataType::MCHijingOverlay: return TString ("MCHijingOverlay");
    default:                        return TString ("???");
  }
}



TString ToTString (const TriggerType tType) {
  switch (tType) {
    case TriggerType::Jet:      return TString ("Jet");
    case TriggerType::MinBias:  return TString ("MinBias");
    default:                    return TString ("???");
  }
}



TString ToTString (const SystFlag sFlag) {
  switch (sFlag) {
    case SystFlag::HITightVar:              return TString ("HITightVar");
    case SystFlag::PionsOnlyVar:            return TString ("PionsOnlyVar");
    case SystFlag::WithPileupVar:           return TString ("WithPileupVar");
    case SystFlag::JetES5PercUpVar:         return TString ("JetES5PercUpVar");
    case SystFlag::JetES5PercDownVar:       return TString ("JetES5PercDownVar");
    case SystFlag::JetES5PercSmearVar:      return TString ("JetES5PercSmearVar");
    case SystFlag::JetES2PercUpVar:         return TString ("JetES2PercUpVar");
    case SystFlag::JetES2PercDownVar:       return TString ("JetES2PercDownVar");
    case SystFlag::JetES2PercSmearVar:      return TString ("JetES2PercSmearVar");
    default:                                return TString ("???");
  }
}



bool IsIons (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::PbPb15 || collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16 || collSys == CollisionSystem::PbPb18);
}



bool IsPbPb (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::PbPb15 || collSys == CollisionSystem::PbPb18);
}



bool IsPbPb18 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::PbPb18);
}



bool IsPbPb15 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::PbPb15);
}



bool IsXeXe (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::XeXe17);
}



bool IspPb (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16);
}



bool IspPb16 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16);
}



bool Ispp (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pp15 || collSys == CollisionSystem::pp17);
}



bool Ispp15 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pp15);
}



bool Ispp17 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pp17);
}



bool IsPeriodA (const CollisionSystem collSys) {
  if (collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16) {
    return true;
  }
  else if (collSys == CollisionSystem::Pbp16) {
    return false;
  }
  std::cout << "In LocalUtilities.cxx::IsPeriodA (const CollisionSystem): Warning: collision system is symmetric, returning true by default." << std::endl;
  return true;
}



bool Is5TeV (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pp15 || collSys == CollisionSystem::PbPb15 || collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::XeXe17 || collSys == CollisionSystem::pp17 || collSys == CollisionSystem::PbPb18);
}



bool Is8TeV (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16);
}



bool Is2018 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::PbPb18);
}



bool Is2017 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::XeXe17 || collSys == CollisionSystem::pp17);
}



bool Is2016 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16);
}



bool Is2015 (const CollisionSystem collSys) {
  return (collSys == CollisionSystem::pp15 || collSys == CollisionSystem::PbPb15);
}



bool IsCollisions (const DataType dType) {
  return (dType == DataType::Collisions);
}



bool IsDataOverlay (const DataType dType) {
  return (dType == DataType::MCDataOverlay);
}



bool IsOverlay (const DataType dType) {
  return (dType == DataType::MCDataOverlay || dType == DataType::MCHijingOverlay);
}



bool IsHijing (const DataType dType) {
  return (dType == DataType::MCHijing || dType == DataType::MCHijingOverlay);
}



bool UseJetTriggers (const TriggerType tType) {
  return tType == TriggerType::Jet;
}



bool UseMinBiasTriggers (const TriggerType tType) {
  return tType == TriggerType::MinBias;
}



bool ToggleSyst (const SystFlag sFlag) {
  switch (sFlag) {
    case SystFlag::HITightVar:          { doHITightVar = true;          return true; }
    case SystFlag::PionsOnlyVar:        { doPionsOnlyVar = true;        return true; }
    case SystFlag::WithPileupVar:       { doWithPileupVar = true;       return true; }
    case SystFlag::JetES5PercUpVar:     { doJetES5PercUpVar = true;     return true; }
    case SystFlag::JetES5PercDownVar:   { doJetES5PercDownVar = true;   return true; }
    case SystFlag::JetES5PercSmearVar:  { doJetES5PercSmearVar = true;  return true; }
    case SystFlag::JetES2PercUpVar:     { doJetES2PercUpVar = true;     return true; }
    case SystFlag::JetES2PercDownVar:   { doJetES2PercDownVar = true;   return true; }
    case SystFlag::JetES2PercSmearVar:  { doJetES2PercSmearVar = true;  return true; }
    default:                                return false;
  }
}




bool IsIons () {
  return IsIons (collisionSystem);
}



bool IsPbPb () {
  return IsPbPb (collisionSystem);
}



bool IsPbPb18 () {
  return IsPbPb18 (collisionSystem);
}



bool IsPbPb15 () {
  return IsPbPb15 (collisionSystem);
}



bool IsXeXe () {
  return IsXeXe (collisionSystem);
}



bool IspPb () {
  return IspPb (collisionSystem);
}



bool IspPb16 () {
  return IspPb16 (collisionSystem);
}



bool Ispp () {
  return Ispp (collisionSystem);
}



bool Ispp15 () {
  return Ispp15 (collisionSystem);
}



bool Ispp17 () {
  return Ispp17 (collisionSystem);
}



bool IsPeriodA () {
  return IsPeriodA (collisionSystem);
}



bool Is5TeV () {
  return Is5TeV (collisionSystem);
}



bool Is8TeV () {
  return Is8TeV (collisionSystem);
}



bool Is2018 () {
  return Is2018 (collisionSystem);
}



bool Is2017 () {
  return Is2017 (collisionSystem);
}



bool Is2016 () {
  return Is2016 (collisionSystem);
}



bool Is2015 () {
  return Is2015 (collisionSystem);
}



bool IsCollisions () {
  return IsCollisions (dataType);
}



bool IsOverlay () {
  return IsOverlay (dataType);
}



bool IsDataOverlay () {
  return IsDataOverlay (dataType);
}



bool IsHijing () {
  return IsHijing (dataType);
}



bool UseJetTriggers () {
  return UseJetTriggers (triggerType);
}



bool UseMinBiasTriggers () {
  return UseMinBiasTriggers (triggerType);
}



/**
 * Returns the CoM boost relevant for asymmetric collision systems (i.e. p+Pb). 0 for everything else.
 */
double GetBoost (int rn) {
  double boost = 0;
  if (IspPb ()) {
    if (Is5TeV ())
      boost = -0.465;
    if (Is8TeV ())
      boost = (rn < 313500 ? -0.465 : 0.465);
  }
  return boost;
}



/**
 * Establishes path variables appropriately.
 */
void SetupDirectories (const TString dataSubDir, const bool addSubDir) {
  rootPath = extWorkPath + "rootFiles/" + dataSubDir + "/";
  plotPath = workPath + "Plots/" + dataSubDir + "/";

  if (addSubDir) {
    if (doHITightVar)
      rootPath = rootPath + "HITightVar/";
    else if (doPionsOnlyVar)
      rootPath = rootPath + "PionsOnlyVar/";
    else if (doWithPileupVar)
      rootPath = rootPath + "WithPileupVar/";
    else if (doJetES5PercUpVar)
      rootPath = rootPath + "JetES5PercUpVar/";
    else if (doJetES5PercDownVar)
      rootPath = rootPath + "JetES5PercDownVar/";
    else if (doJetES5PercSmearVar)
      rootPath = rootPath + "JetES5PercSmearVar/";
    else if (doJetES2PercUpVar)
      rootPath = rootPath + "JetES2PercUpVar/";
    else if (doJetES2PercDownVar)
      rootPath = rootPath + "JetES2PercDownVar/";
    else if (doJetES2PercSmearVar)
      rootPath = rootPath + "JetES2PercSmearVar/";
    else
      rootPath = rootPath + "Nominal/";
  }
}



/**
 * Returns a copy of the histogram detailing the Zdc cuts.
 */
TH1D* GetZdcCuts () {
  if (!IsPbPb ())
    return nullptr; // no cuts needed outside of Pb+Pb

  TFile* zdcFile = new TFile (Form ("%s/ZdcAnalysis/HIRun2PileUp_PbPb5p02_%s.root", rootPath.Data (), Is2015 () ? "2015" : "2018"), "read");
  TH1D* h_zdcCuts = (TH1D*) ((TH1D*) zdcFile->Get (Is2015 () ? "hCut" : "h_cuts")->Clone ("h_zdcCuts"));
  h_zdcCuts->SetDirectory (0);
  zdcFile->Close ();
  SaferDelete (&zdcFile);
  return h_zdcCuts;
}



/**
 * Returns the appropriate file in the given directory.
 * For MC, inFileName MUST be specified.
 * If dataSet == 0, will assume this is a "PhysCont" sample, i.e. an entire collision system in 1 data set.
 */
TFile* GetFile (const char* directory, const int dataSet, const char* inFileName) {
  TFile* file = nullptr;

  // First figure out the file we are looking for
  TString fileIdentifier;
  if (TString (inFileName) == "") {
    if (IsCollisions ())
      fileIdentifier = to_string (dataSet);
    else {
      cout << "Error: In LocalUtilities.cxx::GetFile (const char*, const int, const char*): Cannot identify this MC file! Will return null!" << endl;
      return nullptr;
    }
  }
  else
    fileIdentifier = inFileName;

  // Now get the list of files
  const TString dataPathTemp = dataPath + "/" + directory + "/";
  TSystemDirectory dir (dataPathTemp.Data (), dataPathTemp.Data ());
  TList* sysfiles = dir.GetListOfFiles ();
  if (!sysfiles) {
    cout << "Error: In LocalUtilities.cxx::GetFile (const char*, const int, const char*): Cannot get list of files! Will return null!" << endl;
    return nullptr;
  }
  TSystemFile* sysfile;
  TString fname;
  TIter next (sysfiles);

  while ((sysfile = (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
      if (fname.Contains (fileIdentifier)) {
        break;
      }
    }
  }

  if (!fname.Contains (fileIdentifier)) {
    cout << "Error: In LocalUtilities.cxx::GetFile (const char*, const int, const char*): TFile cannot be found for given data set. Will return null!" << endl;
    return nullptr;
  }

  const TString sourceName = dataPathTemp + fname;

  cout << "Info: In LocalUtilities.cxx::GetFile (const char*, const int, const char*): Resolved file at " << sourceName << endl;
  file = new TFile (sourceName, "read");
  
  if (!file) {
    cout << "Error: In LocalUtilities.cxx::GetFile (const char*, const int, const char*): TFile not obtained for given data set. Will return null!" << endl;
    return nullptr;
  }
  else return file;
} 


/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* directory, const char* inFileName) {
  TString id;
  if (IsCollisions ()) {
    id = to_string (dataSet);
    if (TString (directory).Contains ("pc"))
      id += "_pc";
    else if (TString (directory).Contains ("cc"))
      id += "_cc";
    return id;
  }

  if (IsPbPb ())
    id = (Is2015 () ? "PbPb15" : "PbPb18");
  else if (IspPb ())
    id = (IsPeriodA () ? "pPb16" : "Pbp16");
  else if (Ispp ())
    id = (Is2015 () ? "pp15" : "pp17");

  if (TString (inFileName).Contains ("Sherpa"))
    id = id + "_Sherpa";

  if (IsPbPb () || IspPb ()) {
    if (TString (inFileName).Contains ("pp_Zee"))
      id = id + "_pp_Zee";
    else if (TString (inFileName).Contains ("pn_Zee"))
      id = id + "_pn_Zee";
    else if (TString (inFileName).Contains ("np_Zee"))
      id = id + "_np_Zee";
    else if (TString (inFileName).Contains ("nn_Zee"))
      id = id + "_nn_Zee";
    else if (TString (inFileName).Contains ("Zee"))
      id = id + "_pp_Zee";
    else if (TString (inFileName).Contains ("pp_Zmumu"))
      id = id + "_pp_Zmumu";
    else if (TString (inFileName).Contains ("pn_Zmumu"))
      id = id + "_pn_Zmumu";
    else if (TString (inFileName).Contains ("np_Zmumu"))
      id = id + "_np_Zmumu";
    else if (TString (inFileName).Contains ("nn_Zmumu"))
      id = id + "_nn_Zmumu";
    else if (TString (inFileName).Contains ("Zmumu"))
      id = id + "_pp_Zmumu";
    else if (TString (inFileName).Contains ("Ztautau"))
      id = id + "_pp_Ztautau";
    else if (TString (inFileName).Contains ("ttbar"))
      id = id + "_pp_ttbar";
  } 
  else {
    if (TString (inFileName).Contains ("Zee"))
      id = id + "_Zee";
    else if (TString (inFileName).Contains ("Zmumu"))
      id = id + "_Zmumu";
    else if (TString (inFileName).Contains ("Ztautau"))
      id = id + "_Ztautau";
    else if (TString (inFileName).Contains ("ttbar"))
      id = id + "_ttbar";
  }
  

  if (IsHijing ()) {
    id = id + "_Hijing";
    if (TString (inFileName).Contains ("e4858") || TString (inFileName).Contains ("r11899") || TString (inFileName).Contains ("r11900") || TString (inFileName).Contains ("r11901") || TString (inFileName).Contains ("r11902") || TString (inFileName).Contains ("r11903"))
      id += "_SC";
    else if (TString (inFileName).Contains ("e4962") || TString (inFileName).Contains ("r11892") || TString (inFileName).Contains ("r11893") || TString (inFileName).Contains ("r11894") || TString (inFileName).Contains ("r11895") || TString (inFileName).Contains ("r11898"))
      id += "_MB";

    if (TString (inFileName).Contains ("r11892") || TString (inFileName).Contains ("r11899") || TString (inFileName).Contains ("s3531"))
      id += "_vtxz_n64";
    else if (TString (inFileName).Contains ("r11893") || TString (inFileName).Contains ("r11900") || TString (inFileName).Contains ("s3530"))
      id += "_vtxz_n24";
    else if (TString (inFileName).Contains ("r11894") || TString (inFileName).Contains ("r11901") || TString (inFileName).Contains ("s3520"))
      id += "_vtxz_0";
    else if (TString (inFileName).Contains ("r11895") || TString (inFileName).Contains ("r11902") || TString (inFileName).Contains ("s3528"))
      id += "_vtxz_p24";
    else if (TString (inFileName).Contains ("r11898") || TString (inFileName).Contains ("r11903") || TString (inFileName).Contains ("s3529"))
      id += "_vtxz_p64";
  }

  if (TString (inFileName).Contains ("ptmin25"))
    id = id + "_ptmin25";
  else if (TString (inFileName).Contains ("ygt175"))
    id = id + "_ygt175";

  return id;
}


/**
 * Returns the proper jet trigger luminosity for this data set in nb^-1
 */
double GetJetLuminosity () {
  if (IspPb ())
    return 0.355545; // in nb^-1
  else if (Ispp ())
    return 3574.87; // in nb^-1
  return -1;
}


/**
 * Returns true if this jet passes selection criteria.
 */
bool MeetsJetAcceptanceCuts (int iJ) {
  if (fabs (akt4_hi_jet_eta_xcalib[iJ]) > 2.8)
    return false;
  if (IspPb () && InDisabledHEC (akt4_hi_jet_eta_xcalib[iJ], akt4_hi_jet_phi[iJ]))
    return false;
  return true;
}


/**
 * Returns true if this jet pT passes pT cuts.
 */
bool MeetsJetPtCut (double jpt) {
  if (jet_min_pt > 0) {
    if (jpt < jet_min_pt)
      return false;
    if (jet_max_pt > jet_min_pt && jpt > jet_max_pt)
      return false;
  }
  if (jet_max_pt > 0 && jpt > jet_max_pt)
    return false;
  return true;
}


/**
 * Returns true if this track passes selection criteria.
 */
bool MeetsTrackCuts (int iTrk) {
  if (trk_pt[iTrk] < min_trk_pt)
    return false; // track minimum pT
  if (fabs (trk_eta[iTrk]) > 2.5)
    return false; // track maximum eta

  if (doHITightVar && !trk_HItight[iTrk])
    return false;
  else if (!doHITightVar && !trk_HIloose[iTrk])
    return false;

  if (IsPbPb ()) {
    if (fabs (trk_d0sig[iTrk]) > 3.0)
      return false; // d0 significance cut in Pb+Pb
    if (fabs (trk_z0sig[iTrk]) > 3.0)
      return false; // z0 significance cut in Pb+Pb
  }
  return true;
}


/**
 * Returns the appropriate per-jet reweighting factor. Takes in coordinates for an anti-kT R=0.4 HI jet (pT, eta, & phi).
 * Returns 0 if the jet is outside the acceptance.
 */
double GetAkt4JetWeight (const float jpt, const float jeta, const float jphi, const float jetr) {
  const double accept = ((IspPb () & InDisabledHEC (jeta, jphi)) || fabs (jeta) > 2.8 ? 0. : 1.);
  const double hecwgt = (IspPb () && jeta > 1.1 && jeta < 3.6 ? (2.*M_PI / (3.*M_PI/2. - 2*jetr)) : 1.);
  return accept * hecwgt;
}


} // end namespace

#endif
