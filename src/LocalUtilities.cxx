#ifndef __LocalUtilities_cxx__
#define __LocalUtilities_cxx__

#include "LocalUtilities.h"
#include "TreeVariables.h"
#include "Params.h"
#include "CentralityDefs.h"
#include "Process.h"
#include "RunCorrelator.h"

#include <AtlasUtils.h>

#include <Utilities.h>

#include <TSystemDirectory.h>

#include <math.h>
#include <iostream>
#include <fstream>


namespace JetHadronCorrelations {


TString workPath = TString (std::getenv ("JETHADRONCORR_PATH"));
TString extWorkPath = TString (std::getenv ("JETHADRONCORR_DATA_PATH")) + "/";
TString rootPath = extWorkPath + "rootFiles/";
TString dataPath = extWorkPath + "data/";



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
    case TriggerType::None:       return TString ("None");
    case TriggerType::Jet50GeV:   return TString ("Jet50GeV");
    case TriggerType::Jet100GeV:  return TString ("Jet100GeV");
    case TriggerType::MinBias:    return TString ("MinBias");
    default:                      return TString ("???");
  }
}



TString ToTString (const SystFlag sFlag) {
  switch (sFlag) {
    case SystFlag::None:                    return TString ("None");
    case SystFlag::HITightVar:              return TString ("HITightVar");
    case SystFlag::PionsOnlyVar:            return TString ("PionsOnlyVar");
    case SystFlag::WithPileupVar:           return TString ("WithPileupVar");
    case SystFlag::FcalCentVar:             return TString ("FcalCentVar");
    case SystFlag::FineFcalCentVar:         return TString ("FineFcalCentVar");
    case SystFlag::JetES5PercUpVar:         return TString ("JetES5PercUpVar");
    case SystFlag::JetES5PercDownVar:       return TString ("JetES5PercDownVar");
    case SystFlag::JetES5PercSmearVar:      return TString ("JetES5PercSmearVar");
    case SystFlag::JetES2PercUpVar:         return TString ("JetES2PercUpVar");
    case SystFlag::JetES2PercDownVar:       return TString ("JetES2PercDownVar");
    case SystFlag::JetES2PercSmearVar:      return TString ("JetES2PercSmearVar");
    default:                                return TString ("???");
  }
}



bool SetCollisionSystem (TString ts) {
  for (CollisionSystem cs : AllCollisionSystem) {
    if (ts == ToTString (cs)) { collisionSystem = cs; return true; }
  }
  return false;
}



bool SetDataType (TString ts) {
  for (DataType dt : AllDataType) {
    if (ts == ToTString (dt)) { dataType = dt; return true; }
  }
  return false;
}



bool SetTriggerType (TString ts) {
  for (TriggerType tt : AllTriggerType) {
    if (ts == ToTString (tt)) { triggerType = tt; return true; }
  }
  return false;
}



bool SetSystFlag (TString ts) {
  for (SystFlag sf : AllSystFlag) {
    if (ts == ToTString (sf)) { systFlag = sf; return true; }
  }
  return false;
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
  return UseJet50GeVTriggers (tType) || UseJet100GeVTriggers (tType);
}



bool UseJet50GeVTriggers (const TriggerType tType) {
  return tType == TriggerType::Jet50GeV;
}



bool UseJet100GeVTriggers (const TriggerType tType) {
  return tType == TriggerType::Jet100GeV;
}



bool UseMinBiasTriggers (const TriggerType tType) {
  return tType == TriggerType::MinBias;
}



bool DoHITightVar (const SystFlag sFlag) {
  return sFlag == SystFlag::HITightVar;
}



bool DoPionsOnlyVar (const SystFlag sFlag) {
  return sFlag == SystFlag::PionsOnlyVar;
}



bool DoWithPileupVar (const SystFlag sFlag) {
  return sFlag == SystFlag::WithPileupVar;
}



bool DoFcalCentVar (const SystFlag sFlag) {
  return sFlag == SystFlag::FcalCentVar;
}



bool DoFineFcalCentVar (const SystFlag sFlag) {
  return sFlag == SystFlag::FineFcalCentVar;
}



bool DoJetES5PercUpVar (const SystFlag sFlag) {
  return sFlag == SystFlag::JetES5PercUpVar;
}



bool DoJetES5PercDownVar (const SystFlag sFlag) {
  return sFlag == SystFlag::JetES5PercDownVar;
}



bool DoJetES5PercSmearVar (const SystFlag sFlag) {
  return sFlag == SystFlag::JetES5PercSmearVar;
}



bool DoJetES2PercUpVar (const SystFlag sFlag) {
  return sFlag == SystFlag::JetES2PercUpVar;
}



bool DoJetES2PercDownVar (const SystFlag sFlag) {
  return sFlag == SystFlag::JetES2PercDownVar;
}



bool DoJetES2PercSmearVar (const SystFlag sFlag) {
  return sFlag == SystFlag::JetES2PercSmearVar;
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



bool UseJet50GeVTriggers () {
  return UseJet50GeVTriggers (triggerType);
}



bool UseJet100GeVTriggers () {
  return UseJet100GeVTriggers (triggerType);
}



bool UseMinBiasTriggers () {
  return UseMinBiasTriggers (triggerType);
}



bool DoHITightVar () {
  return DoHITightVar (systFlag);
}



bool DoPionsOnlyVar () {
  return DoPionsOnlyVar (systFlag);
}



bool DoWithPileupVar () {
  return DoWithPileupVar (systFlag);
}



bool DoFcalCentVar () {
  return DoFcalCentVar (systFlag);
}



bool DoFineFcalCentVar () {
  return DoFineFcalCentVar (systFlag);
}



bool DoJetES5PercUpVar () {
  return DoJetES5PercUpVar (systFlag);
}



bool DoJetES5PercDownVar () {
  return DoJetES5PercDownVar (systFlag);
}



bool DoJetES5PercSmearVar () {
  return DoJetES5PercSmearVar (systFlag);
}



bool DoJetES2PercUpVar () {
  return DoJetES2PercUpVar (systFlag);
}


bool DoJetES2PercDownVar () {
  return DoJetES2PercDownVar (systFlag);
}



bool DoJetES2PercSmearVar () {
  return DoJetES2PercSmearVar (systFlag);
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

  if (addSubDir) {
    if (DoHITightVar ())
      rootPath = rootPath + "HITightVar/";
    else if (DoPionsOnlyVar ())
      rootPath = rootPath + "PionsOnlyVar/";
    else if (DoWithPileupVar ())
      rootPath = rootPath + "WithPileupVar/";
    else if (DoFcalCentVar ())
      rootPath = rootPath + "FcalCentVar/";
    else if (DoFineFcalCentVar ())
      rootPath = rootPath + "FineFcalCentVar/";
    else if (DoJetES5PercUpVar ())
      rootPath = rootPath + "JetES5PercUpVar/";
    else if (DoJetES5PercDownVar ())
      rootPath = rootPath + "JetES5PercDownVar/";
    else if (DoJetES5PercSmearVar ())
      rootPath = rootPath + "JetES5PercSmearVar/";
    else if (DoJetES2PercUpVar ())
      rootPath = rootPath + "JetES2PercUpVar/";
    else if (DoJetES2PercDownVar ())
      rootPath = rootPath + "JetES2PercDownVar/";
    else if (DoJetES2PercSmearVar ())
      rootPath = rootPath + "JetES2PercSmearVar/";
    else
      rootPath = rootPath + "Nominal/";
  }
}



/**
 * Looks up MC sample cross section, filter efficiency, and number of events.
 */
bool GetMCWeights (TString fname) {
  ifstream f_wgts;
  f_wgts.open ("MC_Weights.dat");

  string line;
  while (getline (f_wgts, line)) {

    vector <string> words = {};

    string word = "";
    for (auto x : line) {
        if (x == ' ') {
            words.push_back (word);
            word = "";
        }
        else {
            word = word + x;
        }
    }
    words.push_back (word);

    assert (words.size () == 4);

    if (fname.Contains (TString (words[0]))) {
      crossSectionPicoBarns = atof (words[1].c_str ());
      mcFilterEfficiency = atof (words[2].c_str ());
      mcNumberEvents = atof (words[3].c_str ());
      return true;
    }
  }
  return false;
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
 * Returns the probability histogram of each FCal ET value in 0-20% ZDC events.
 */
TH1D* GetFCalZdcWeights () {
  TDirectory* gdir = gDirectory;
  
  TString fname = Form ("%s/aux/CentralityDistributions.root", workPath.Data ());
  std::cout << "Trying to resolve centrality distributions file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  TH2D* h2 = (TH2D*) infile->Get ("h2_mb_Pb_fcal_et_zdc_calibE");
  TH1D* h = (TH1D*) h2->ProjectionX ("h_fcal_et_zdcWeights", h2->GetYaxis ()->FindBin (zdcCentBins[numZdcCentBins-1]), h2->GetYaxis ()->GetNbins ());

  //TH1D* h_fcal = (TH1D*) h2->ProjectionX ("h_fcal_et");
  //h->Divide (h_fcal);
  //SaferDelete (&h_fcal);

  h->Scale (1./h->Integral ());
  std::cout << "Loaded centrality distribution, closing file" << std::endl;

  h->SetDirectory (gdir);

  infile->Close ();
  SaferDelete (&infile);

  return h;
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

  if (TString (inFileName).Contains ("420010"))
    id = id + "_JZ0";
  else if (TString (inFileName).Contains ("420011"))
    id = id + "_JZ1";
  else if (TString (inFileName).Contains ("420012"))
    id = id + "_JZ2";
  else if (TString (inFileName).Contains ("420013"))
    id = id + "_JZ3";
  else if (TString (inFileName).Contains ("420014"))
    id = id + "_JZ4";
  else if (TString (inFileName).Contains ("420015"))
    id = id + "_JZ5";

  if (TString (inFileName).Contains ("e4108"))
    id = id + "_SampleA";
  else if (TString (inFileName).Contains ("e6608"))
    id = id + "_SampleB";

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

  if (DoHITightVar () && !trk_HItight[iTrk])
    return false;
  else if (!DoHITightVar () && !trk_HIloose[iTrk])
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



/**
 * Returns the tracking efficiency histograms.
 */
TH2D* LoadTrackingEfficiency () {
  TDirectory* gdir = gDirectory;

  TString fname = Form ("%s/TrackingPerformance/Nominal/outFile.root", rootPath.Data ());
  std::cout << "Trying to resolve tracking performance file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  TH2D* h2 = (TH2D*) infile->Get (Form ("h2_truth_matched_reco_tracks_%s", "pp"))->Clone ("h2_tracking_efficiency");
  h2->Divide ((TH2D*) infile->Get (Form ("h2_truth_tracks_%s", "pp")));
  std::cout << "Loaded tracking efficiencies, closing file" << std::endl;

  h2->SetDirectory (gdir);

  infile->Close ();
  SaferDelete (&infile);

  return h2;
}



/**
 * Returns the tracking purity histograms.
 */
TH2D* LoadTrackingPurity () {
  TDirectory* gdir = gDirectory;

  TString fname = Form ("%s/TrackingPerformance/Nominal/outFile.root", rootPath.Data ());
  std::cout << "Trying to resolve tracking performance file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  TH2D* h2 = (TH2D*) infile->Get (Form ("h2_primary_reco_tracks_%s", "pp"))->Clone ("h2_tracking_purity");
  h2->Divide ((TH2D*) infile->Get (Form ("h2_reco_tracks_%s", "pp")));
  std::cout << "Loaded tracking purities, closing file" << std::endl;

  h2->SetDirectory (gdir);

  infile->Close ();
  SaferDelete (&infile);

  return h2;
}



/**
 * Converts a TProfile to a TGraph assuming the x-axis of the TProfile is the y-axis of the TGraph.
 */
TGraphErrors* TProfY2TGE (TProfile* py) {
  TGraphErrors* g = new TGraphErrors ();
  for (int iX = 1; iX <= py->GetNbinsX (); iX++) {
    g->SetPoint (g->GetN (), py->GetBinContent (iX), py->GetBinCenter (iX));
    g->SetPointError (g->GetN ()-1, py->GetBinError (iX), py->GetBinWidth (iX) / 2.);
  }
  return g;
}

/**
 * Converts a TProfile to a TGraph assuming the x-axis of the TProfile is the x-axis of the TGraph.
 */
TGraphErrors* TProfX2TGE (TProfile* px) {
  TGraphErrors* g = new TGraphErrors ();
  for (int iX = 1; iX <= px->GetNbinsX (); iX++) {
    g->SetPoint (g->GetN (), px->GetBinCenter (iX), px->GetBinContent (iX));
    g->SetPointError (g->GetN ()-1, px->GetBinWidth (iX) / 2., px->GetBinError (iX));
  }
  return g;
}


} // end namespace

#endif
