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
#include <ArrayTemplates.h>

#include <TSystemDirectory.h>
#include <TLorentzVector.h>

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include <stdlib.h>
#include <stdio.h>


namespace JetHadronCorrelations {


TString workPath = TString (std::getenv ("JETHADRONCORR_PATH"));
TString extWorkPath = TString (std::getenv ("JETHADRONCORR_DATA_PATH")) + "/";
TString rootPath = extWorkPath + "rootFiles/";
TString dataPath = extWorkPath + "data/";



TString ToTString (const CollisionSystem& collSys) {
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
 
 
 
TString ToTString (const DataType& dType) {
  switch (dType) {
    case DataType::Collisions:      return TString ("Collisions");
    case DataType::MCSignal:        return TString ("MCSignal");
    case DataType::MCDataOverlay:   return TString ("MCDataOverlay");
    case DataType::MCHijing:        return TString ("MCHijing");
    case DataType::MCHijingOverlay: return TString ("MCHijingOverlay");
    default:                        return TString ("???");
  }
}



TString ToTString (const TriggerType& tType) {
  switch (tType) {
    case TriggerType::None:       return TString ("None");
    case TriggerType::J50:   return TString ("J50");
    case TriggerType::J100:  return TString ("J100");
    case TriggerType::MinBias:    return TString ("MinBias");
    default:                      return TString ("???");
  }
}



TString ToTString (const SystFlag& sFlag) {
  switch (sFlag) {
    case SystFlag::Nominal:                     return TString ("Nominal");
    case SystFlag::HITightVar:                  return TString ("HITightVar");
    case SystFlag::HILooseVar:                  return TString ("HILooseVar");
    case SystFlag::TrkEffVar:                   return TString ("TrkEffVar");
    case SystFlag::FakeRateVar:                 return TString ("FakeRateVar");
    case SystFlag::PrimFitVar:                  return TString ("PrimFitVar");
    case SystFlag::JetPrimFracVar:              return TString ("JetPrimFracVar");
    case SystFlag::PartSpcVar:                  return TString ("PartSpcVar");
    case SystFlag::FcalCentVar:                 return TString ("FcalCentVar");
    case SystFlag::FineFcalCentVar:             return TString ("FineFcalCentVar");
    case SystFlag::MixCatVar1:                  return TString ("MixCatVar1");
    case SystFlag::MixCatVar2:                  return TString ("MixCatVar2");
    case SystFlag::MixCatVar3:                  return TString ("MixCatVar3");
    case SystFlag::MixCatVar4:                  return TString ("MixCatVar4");
    case SystFlag::MixCatVar5:                  return TString ("MixCatVar5");
    case SystFlag::MixCatVar6:                  return TString ("MixCatVar6");
    case SystFlag::JESVar0:                     return TString ("JESVar0");
    case SystFlag::JESVar1:                     return TString ("JESVar1");
    case SystFlag::JESVar2:                     return TString ("JESVar2");
    case SystFlag::JESVar3:                     return TString ("JESVar3");
    case SystFlag::JESVar4:                     return TString ("JESVar4");
    case SystFlag::JESVar5:                     return TString ("JESVar5");
    case SystFlag::JESVar6:                     return TString ("JESVar6");
    case SystFlag::JESVar7:                     return TString ("JESVar7");
    case SystFlag::JESVar8:                     return TString ("JESVar8");
    case SystFlag::JESVar9:                     return TString ("JESVar9");
    case SystFlag::JESVar10:                    return TString ("JESVar10");
    case SystFlag::JESVar11:                    return TString ("JESVar11");
    case SystFlag::JESVar12:                    return TString ("JESVar12");
    case SystFlag::JESVar13:                    return TString ("JESVar13");
    case SystFlag::JESVar14:                    return TString ("JESVar14");
    case SystFlag::JESVar15:                    return TString ("JESVar15");
    case SystFlag::JESVar16:                    return TString ("JESVar16");
    case SystFlag::JESVar17:                    return TString ("JESVar17");
    case SystFlag::JESVar18:                    return TString ("JESVar18");
    case SystFlag::JESVar19:                    return TString ("JESVar19");
    case SystFlag::JESVar20:                    return TString ("JESVar20");
    case SystFlag::JERVar0:                     return TString ("JERVar0");
    case SystFlag::JERVar1:                     return TString ("JERVar1");
    case SystFlag::JERVar2:                     return TString ("JERVar2");
    case SystFlag::JERVar3:                     return TString ("JERVar3");
    case SystFlag::JERVar4:                     return TString ("JERVar4");
    case SystFlag::JERVar5:                     return TString ("JERVar5");
    case SystFlag::JERVar6:                     return TString ("JERVar6");
    case SystFlag::JERVar7:                     return TString ("JERVar7");
    case SystFlag::JERVar8:                     return TString ("JERVar8");
    case SystFlag::JERVar9:                     return TString ("JERVar9");
    case SystFlag::JERVar10:                    return TString ("JERVar10");
    case SystFlag::MCTruthJetsTruthParts:       return TString ("MCTruthJetsTruthParts");
    case SystFlag::MCRecoJetsTruthParts:        return TString ("MCRecoJetsTruthParts");
    case SystFlag::MCRecoJetsTruthMatchedParts: return TString ("MCRecoJetsTruthMatchedParts");
    case SystFlag::MCFCalWeighted:              return TString ("MCFCalWeighted");
    default:                                    return TString ("???");
  }
}



float GetRadius (const JetRadius& r) {
  switch (r) {
    case JetRadius::R0p2:     return 0.2;
    case JetRadius::R0p3:     return 0.3;
    case JetRadius::R0p4:     return 0.4;
    case JetRadius::R0p6:     return 0.6;
    case JetRadius::R0p8:     return 0.8;
    case JetRadius::R1p0:     return 1.0;
    case JetRadius::Invalid:  return 0.;
    default:                  return -1.;
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



bool IsIons (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::PbPb15 || collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16 || collSys == CollisionSystem::PbPb18);
}



bool IsPbPb (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::PbPb15 || collSys == CollisionSystem::PbPb18);
}



bool IsPbPb18 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::PbPb18);
}



bool IsPbPb15 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::PbPb15);
}



bool IsXeXe (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::XeXe17);
}



bool IspPb (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16);
}



bool IspPb16 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16);
}



bool Ispp (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pp15 || collSys == CollisionSystem::pp17);
}



bool Ispp15 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pp15);
}



bool Ispp17 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pp17);
}



bool IsPeriodA (const CollisionSystem& collSys) {
  if (collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16) {
    return true;
  }
  else if (collSys == CollisionSystem::Pbp16) {
    return false;
  }
  std::cout << "In LocalUtilities.cxx::IsPeriodA (const CollisionSystem): Warning: collision system is symmetric, returning true by default." << std::endl;
  return true;
}



bool Is5TeV (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pp15 || collSys == CollisionSystem::PbPb15 || collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::XeXe17 || collSys == CollisionSystem::pp17 || collSys == CollisionSystem::PbPb18);
}



bool Is8TeV (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16);
}



bool Is2018 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::PbPb18);
}



bool Is2017 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::XeXe17 || collSys == CollisionSystem::pp17);
}



bool Is2016 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pPb16s5TeV || collSys == CollisionSystem::pPb16 || collSys == CollisionSystem::Pbp16);
}



bool Is2015 (const CollisionSystem& collSys) {
  return (collSys == CollisionSystem::pp15 || collSys == CollisionSystem::PbPb15);
}



bool IsCollisions (const DataType& dType) {
  return (dType == DataType::Collisions);
}



bool IsDataOverlay (const DataType& dType) {
  return (dType == DataType::MCDataOverlay);
}



bool IsOverlay (const DataType& dType) {
  return (dType == DataType::MCDataOverlay || dType == DataType::MCHijingOverlay);
}



bool IsHijing (const DataType& dType) {
  return (dType == DataType::MCHijing || dType == DataType::MCHijingOverlay);
}



bool UseJetTriggers (const TriggerType& tType) {
  return UseJ50Triggers (tType) || UseJ100Triggers (tType);
}



bool UseJ50Triggers (const TriggerType& tType) {
  return tType == TriggerType::J50;
}



bool UseJ100Triggers (const TriggerType& tType) {
  return tType == TriggerType::J100;
}



bool UseMinBiasTriggers (const TriggerType& tType) {
  return tType == TriggerType::MinBias;
}



bool DoHITightVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::HITightVar;
}



bool DoHILooseVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::HILooseVar;
}



bool DoTrkEffVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::TrkEffVar;
}



bool DoFakeRateVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::FakeRateVar;
}



bool DoPrimFitVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::PrimFitVar;
}



bool DoJetPrimFracVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::JetPrimFracVar;
}



bool DoPartSpcVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::PartSpcVar;
}



bool DoFcalCentVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::FcalCentVar;
}



bool DoFineFcalCentVar (const SystFlag& sFlag) {
  return sFlag == SystFlag::FineFcalCentVar;
}



bool DoMixCatVar1 (const SystFlag& sFlag) {
  return sFlag == SystFlag::MixCatVar1;
}



bool DoMixCatVar2 (const SystFlag& sFlag) {
  return sFlag == SystFlag::MixCatVar2;
}



bool DoMixCatVar3 (const SystFlag& sFlag) {
  return sFlag == SystFlag::MixCatVar3;
}



bool DoMixCatVar4 (const SystFlag& sFlag) {
  return sFlag == SystFlag::MixCatVar4;
}



bool DoMixCatVar5 (const SystFlag& sFlag) {
  return sFlag == SystFlag::MixCatVar5;
}



bool DoMixCatVar6 (const SystFlag& sFlag) {
  return sFlag == SystFlag::MixCatVar6;
}



int GetNJESVar (const SystFlag& sFlag) {
  switch (sFlag) {
    case SystFlag::JESVar0:   return 0;
    case SystFlag::JESVar1:   return 1;
    case SystFlag::JESVar2:   return 2;
    case SystFlag::JESVar3:   return 3;
    case SystFlag::JESVar4:   return 4;
    case SystFlag::JESVar5:   return 5;
    case SystFlag::JESVar6:   return 6;
    case SystFlag::JESVar7:   return 7;
    case SystFlag::JESVar8:   return 8;
    case SystFlag::JESVar9:   return 9;
    case SystFlag::JESVar10:  return 10;
    case SystFlag::JESVar11:  return 11;
    case SystFlag::JESVar12:  return 12;
    case SystFlag::JESVar13:  return 13;
    case SystFlag::JESVar14:  return 14;
    case SystFlag::JESVar15:  return 15;
    case SystFlag::JESVar16:  return 16;
    case SystFlag::JESVar17:  return 17;
    case SystFlag::JESVar18:  return 18;
    case SystFlag::JESVar19:  return 19;
    case SystFlag::JESVar20:  return 20;
    default:                  return -1;
  }
}



int GetNJERVar (const SystFlag& sFlag) {
  switch (sFlag) {
    case SystFlag::JERVar0:   return 0;
    case SystFlag::JERVar1:   return 1;
    case SystFlag::JERVar2:   return 2;
    case SystFlag::JERVar3:   return 3;
    case SystFlag::JERVar4:   return 4;
    case SystFlag::JERVar5:   return 5;
    case SystFlag::JERVar6:   return 6;
    case SystFlag::JERVar7:   return 7;
    case SystFlag::JERVar8:   return 8;
    case SystFlag::JERVar9:   return 9;
    case SystFlag::JERVar10:  return 10;
    default:                  return -1;
  }
}



bool DoMCTruthJetsTruthParts (const SystFlag& sFlag) {
  return sFlag == SystFlag::MCTruthJetsTruthParts;
}



bool DoMCRecoJetsTruthParts (const SystFlag& sFlag) {
  return sFlag == SystFlag::MCRecoJetsTruthParts;
}



bool DoMCRecoJetsTruthMatchedParts (const SystFlag& sFlag) {
  return sFlag == SystFlag::MCRecoJetsTruthMatchedParts;
}



bool DoMCFCalWeights (const SystFlag& sFlag) {
  return sFlag == SystFlag::MCFCalWeighted;
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



bool UseJ50Triggers () {
  return UseJ50Triggers (triggerType);
}



bool UseJ100Triggers () {
  return UseJ100Triggers (triggerType);
}



bool UseMinBiasTriggers () {
  return UseMinBiasTriggers (triggerType);
}



bool DoHITightVar () {
  return DoHITightVar (systFlag);
}



bool DoHILooseVar () {
  return DoHILooseVar (systFlag);
}



bool DoTrkEffVar () {
  return DoTrkEffVar (systFlag);
}



bool DoFakeRateVar () {
  return DoFakeRateVar (systFlag);
}



bool DoPrimFitVar () {
  return DoPrimFitVar (systFlag);
}



bool DoJetPrimFracVar () {
  return DoJetPrimFracVar (systFlag);
}



bool DoPartSpcVar () {
  return DoPartSpcVar (systFlag);
}



bool DoFcalCentVar () {
  return DoFcalCentVar (systFlag);
}



bool DoFineFcalCentVar () {
  return DoFineFcalCentVar (systFlag);
}



bool DoMixCatVar1 () {
  return DoMixCatVar1 (systFlag);
}



bool DoMixCatVar2 () {
  return DoMixCatVar2 (systFlag);
}



bool DoMixCatVar3 () {
  return DoMixCatVar3 (systFlag);
}



bool DoMixCatVar4 () {
  return DoMixCatVar4 (systFlag);
}



bool DoMixCatVar5 () {
  return DoMixCatVar5 (systFlag);
}



bool DoMixCatVar6 () {
  return DoMixCatVar6 (systFlag);
}



int GetNJESVar () {
  return GetNJESVar (systFlag);
}



int GetNJERVar () {
  return GetNJERVar (systFlag);
}



bool DoMCTruthJetsTruthParts () {
  return DoMCTruthJetsTruthParts (systFlag);
}



bool DoMCRecoJetsTruthParts () {
  return DoMCRecoJetsTruthParts (systFlag);
}



bool DoMCRecoJetsTruthMatchedParts () {
  return DoMCRecoJetsTruthMatchedParts (systFlag);
}



bool DoMCFCalWeights () {
  return DoMCFCalWeights (systFlag);
}



bool UseTruthJets () {
  return DoMCTruthJetsTruthParts ();
}



bool UseTruthParticles () {
  return DoMCTruthJetsTruthParts () || DoMCRecoJetsTruthParts ();
}



bool UseTruthMatchedParticles () {
  return DoMCRecoJetsTruthMatchedParts ();
}



bool UseMCFCalWeights () {
  return DoMCFCalWeights () && IsDataOverlay ();
}





/**
 * Returns the CoM boost relevant for asymmetric collision systems (i.e. p+Pb). 0 for everything else.
 */
double GetBoost (int rn) {
  double boost = 0;
  if (IspPb ()) {
    boost = (IsPeriodA () ? -0.465 : 0.465);
  }
  return boost;
}



/**
 * Establishes path variables appropriately.
 */
void SetupDirectories (const TString& dataSubDir, const bool addSubDir) {
  rootPath = extWorkPath + "rootFiles/" + dataSubDir + "/";

  if (addSubDir)
    rootPath = rootPath + ToTString (systFlag) + "/";
}



/**
 * Looks up MC sample cross section, filter efficiency, and number of events.
 */
bool GetMCWeights (const TString& fname) {
  ifstream f_wgts;
  f_wgts.open (Form ("%s/aux/MC_Weights.dat", workPath.Data ()));

  std::string line;
  while (getline (f_wgts, line)) {

    std::istringstream lineStream (line);

    std::vector <std::string> words = {};
    std::string word;

    while (lineStream >> word)
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
 * Returns the minimum anti-kT R=0.4 truth jet pT for this JZXR04 MC sample.
 */
float GetJZXR04MinPt (const TString& fname) {
  if (!IsCollisions ()) {
    if (TString (fname).Contains ("JZ0") || TString (fname).Contains ("Hijing"))      return 0;
    else if (TString (fname).Contains ("JZ1")) return 20;
    else if (TString (fname).Contains ("JZ2")) return 60;
    else if (TString (fname).Contains ("JZ3")) return 160;
    else if (TString (fname).Contains ("JZ4")) return 400;
    else if (TString (fname).Contains ("JZ5")) return 800;
  }
  return 0;
}



/**
 * Returns the maximum anti-kT R=0.4 truth jet pT for this JZXR04 MC sample.
 */
float GetJZXR04MaxPt (const TString& fname) {
  if (!IsCollisions ()) {
    if (TString (fname).Contains ("JZ0") || TString (fname).Contains ("Hijing"))      return 20;
    else if (TString (fname).Contains ("JZ1")) return 60;
    else if (TString (fname).Contains ("JZ2")) return 160;
    else if (TString (fname).Contains ("JZ3")) return 400;
    else if (TString (fname).Contains ("JZ4")) return 800;
    else if (TString (fname).Contains ("JZ5")) return 1300;
  }
  return FLT_MAX;
}



/**
 * Returns the residual JZ scaling factor to smooth the spectrum transition at 60 GeV from SoftQCD to HardQCD
 */
float GetJZScaleFactor (const TString& fname) {
  if (!IsCollisions () && !(fname.Contains ("JZ0") || fname.Contains ("JZ1") || fname.Contains ("Hijing")))
    return 0.859275;
  return 1.;
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
 * Returns a copy of the histogram detailing the probability of sampling a given MC event.
 */
TH1D* GetFCalResamplingProbs () {
  if (!IsDataOverlay () || !IspPb ())
    return nullptr;

  TFile* probsFile = new TFile (Form ("%s/aux/MCResampling.root", workPath.Data ()), "read");
  TH1D* h_probs = (TH1D*) ((TH1D*) probsFile->Get ("hjet_ratio")->Clone ("h_MCResamplingProbs"));
  h_probs->SetDirectory (0);
  probsFile->Close ();
  SaferDelete (&probsFile);
  return h_probs;
}



/**
 * Returns a copy of the histogram with data/MC ratios of the FCal Et distribution
 */
TH1D* GetMCFCalWeights () {
  if (!IsDataOverlay () || !IspPb ())
    return nullptr;

  TFile* wgtsFile = new TFile (Form ("%s/aux/MCFCalWeights.root", workPath.Data ()), "read");
  TH1D* h_wgts = (TH1D*) ((TH1D*) wgtsFile->Get ("h_rat")->Clone ("h_MCFCalWeights"));
  h_wgts->SetDirectory (0);
  wgtsFile->Close ();
  SaferDelete (&wgtsFile);
  return h_wgts;
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
  TH1D* h = (TH1D*) h2->ProjectionX ("h_fcal_et_zdcWeights", h2->GetYaxis ()->FindBin (zdcCentBins[nZdcCentBins-1]), h2->GetYaxis ()->GetNbins ());

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
 * Returns a map from event numbers to pure overlay A-side FCal ET values.
 */
std::map <const unsigned int, float>* GetOverlayFCalMap () {

  TString fname = Form ("%s/aux/OverlayFCalEtMap.dat", workPath.Data ());
  std::cout << "Loading FCal Et map for p+Pb data overlay from " << fname << std::endl;

  ifstream inFile;
  inFile.open (fname);

  char fcalA_et[20];
  char evtNum[20];

  char* endPtr;

  std::map <const unsigned int, float>* overlay_fcalA_et_map = new std::map <const unsigned int, float> ();

  while (inFile) {

    inFile >> evtNum >> fcalA_et;

    overlay_fcalA_et_map->insert (std::pair <const unsigned int, float> (std::strtoul (evtNum, &endPtr, 10), std::strtof (fcalA_et, &endPtr)));

  }

  std::cout << "Loaded map, closing file" << std::endl;

  inFile.close ();

  return overlay_fcalA_et_map;
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

  if ((IsPbPb () || IspPb ()) && !IsCollisions () && !IsOverlay () && !IsHijing ())
    id = id + "_so";

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

  // FROM Z+h analysis (now deprecated)
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

    // FROM Z+h analysis Hijing samples (now deprecated)
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

  // FROM Z+h analysis Hijing samples (now deprecated)
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
 * Returns true if this truth jet passes selection criteria.
 */
bool MeetsTruthJetAcceptanceCuts (int iTJ, const JetRadius& radius) {
  assert (GetAktTruthJetN (radius) > iTJ);
  if (std::fabs (GetAktTruthJetEta (iTJ, radius)) > 2.8)
    return false;
  if (IspPb () && InDisabledHEC (GetAktTruthJetEta (iTJ, radius), GetAktTruthJetPhi (iTJ, radius), GetRadius (radius)))
    return false; // cut out jets in the disabled HEC
  return true;
}



/**
 * Returns true if this jet passes selection criteria.
 */
bool MeetsJetAcceptanceCuts (int iJ, const JetRadius& radius, const int nJESVar) {
  if (UseTruthJets ())
    return MeetsTruthJetAcceptanceCuts (iJ, radius); // if working at truth level then we should actually be checking the truth jet
  assert (GetAktHIJetN (radius) > iJ);
  assert (nJESVar >= -1 && nJESVar < nJESSys);
  if (std::fabs (GetAktHIJetEta (iJ, radius, nJESVar)) > 2.8)
    return false; // cut out forward jets
  if (IspPb () && InDisabledHEC (GetAktHIJetEta (iJ, radius, nJESVar), GetAktHIJetPhi (iJ, radius, nJESVar), GetRadius (radius)))
    return false; // cut out jets in the disabled HEC
  if (IsCollisions () && Ispp () && UseMinBiasTriggers () && GetAktHIJetTiming (iJ, radius) > 10)
    return false; // cut out jets with bad timing in pp if using MinBias trigger
  if (!GetAktHIJetCleaning (iJ, radius))
    return false; // cut on jets not meeting jet cleaning cut
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
bool MeetsTrackCuts (int iTrk, const int nWPVar) {
  if (trk_pt[iTrk] < min_trk_pt)
    return false; // track minimum pT
  if (std::fabs (trk_eta[iTrk]) > 2.5)
    return false; // track maximum eta
  if (trk_charge[iTrk] == 0)
    return false; // cut on neutrals

  assert (nWPVar >= 0 && nWPVar < (int)trackWPs.size ());

  if (!UseTruthParticles ()) {
    if (!trackWPs[nWPVar][iTrk])
      return false;
    if (UseTruthMatchedParticles () && trk_prob_truth[iTrk] < 0.5) 
      return false;
  }
  else if (trk_truth_barcode[iTrk] <= 0 || 200000 <= trk_truth_barcode[iTrk] || !trk_truth_isHadron[iTrk])
    return false;

  // Heavy ion additional vertex matching cuts
  //if (IsPbPb ()) {
  //  if (std::fabs (trk_d0sig[iTrk]) > 3.0)
  //    return false; // d0 significance cut in Pb+Pb
  //  if (std::fabs (trk_z0sig[iTrk]) > 3.0)
  //    return false; // z0 significance cut in Pb+Pb
  //}
  return true;
}



///**
// * Returns the truth-particle corrected FCal ET values.
// */
//FCalEt GetTruthCorrectedFCal (FCalEt values) {
//  TLorentzVector tlv;
//  for (int iTTrk = 0; iTTrk < truth_trk_n; iTTrk++) {
//    const float eta = truth_trk_eta[iTTrk];
//    const float aeta = std::fabs (eta);
//
//    if (aeta > 4.9 || aeta < 3.2)
//      continue;
//
//    if (truth_trk_barcode[iTTrk] <= 0 || truth_trk_barcode[iTTrk] >= 200000)
//      continue; // cut on secondaries -- we only want particles from the Pythia event! (Otherwise we double count some energy)
//
//    //if (truth_trk_charge[iTTrk] != 0 && truth_trk_pt[iTTrk] < 0.0225) // cut on loopers (they never hit the FCal); formula is pT^max = e * B * R^max where R^max ~ 0.5 * 75 mm for the FCal and B ~ 2.0T. For systematics, could use B ~ 1.7T and R^max = 0.5*70mm or 80mm. Try 2 * pT^max?
//    //  continue;
//
//    if (truth_trk_pt[iTTrk] < 0.120) // cut on particles within 2 sigma of noise.
//      continue;
//
//    //float sf = 1; // area-based correction factor for jets only partially in the FCal.
//    //if (aeta > 4.9) {
//    //  const float delta = aeta - 4.9;
//    //  //sf = std::acos (1 - delta/0.4) / M_PI - (0.4 - delta) * std::sqrt (delta * (0.8 - delta)) / (M_PI * 0.16);
//    //}
//    //else if (aeta < 3.2) {
//    //  const float delta = 3.2 - aeta;
//    //  //sf = std::acos (1 - delta/0.4) / M_PI - (0.4 - delta) * std::sqrt (delta * (0.8 - delta)) / (M_PI * 0.16);
//    //}
//
//    tlv.SetPtEtaPhiM (truth_trk_pt[iTTrk], eta, truth_trk_phi[iTTrk], 0.137); // assume pion mass?
//    const float subet = tlv.Et ();
//
//    if (eta > 0)
//      values.first = values.first - subet;
//    else
//      values.second = values.second - subet;
//  }
//  return values;
//}



///**
// * NOW DEPRECATED
// * Returns the matched truth jet within DR < pi to this HI jet.
// * Returns -1 if no truth jet is matched within this DR range, or the radius is invalid.
// */
//int GetAktTruthJetMatch (const int iJ, const JetRadius& radius, const int nJESVar) {
//  assert (GetAktHIJetN (radius) > iJ);
//  float mindr = M_PI;
//  int match = -1;
//  const float jeta = GetAktHIJetEta (iJ, radius, nJESVar);
//  const float jphi = GetAktHIJetPhi (iJ, radius, nJESVar);
//  const int nTJ = GetAktTruthJetN (radius);
//  for (int iTJ = 0; iTJ < nTJ; iTJ++) {
//    const float dr = DeltaR (jeta, GetAktTruthJetEta (iTJ, radius), jphi, GetAktTruthJetPhi (iTJ, radius));
//    if (dr < mindr) {
//      match = iTJ;
//      mindr = dr;
//    }
//  }
//  if (mindr > GetAktTruthMatchMaxDR (radius))
//    return -1;
//  return match;
//}



/**
 * Returns a vector storing the index of the reco jet match for each truth jet.
 * Elements are -1 if the truth jet has no reco match.
 */
std::vector <short> GetAktRecoJetMatches (const JetRadius& radius, const short nJESVar) {
  if (IsCollisions ())
    return std::vector <short> (0);

  const short nRJ = GetAktHIJetN (radius);

  std::vector <short> recoMatches (GetAktTruthJetN (radius));

  for (short iTJ = 0; iTJ < (short) recoMatches.size (); iTJ++) {
    const float tjeta = GetAktTruthJetEta (iTJ, radius);
    const float tjphi = GetAktTruthJetPhi (iTJ, radius);

    float mindr = GetAktTruthMatchMaxDR (radius);
    short recoMatch = -1;

    for (short iRJ = 0; iRJ < nRJ; iRJ++) {
      if (GetAktHIJetPt (iRJ, radius, nJESVar) < GetAktTruthMatchMinRecoPt (radius))
        continue;
      const float dr = DeltaR (tjeta, GetAktHIJetEta (iRJ, radius, nJESVar), tjphi, GetAktHIJetPhi (iRJ, radius, nJESVar));
      if (dr < mindr) {
        recoMatch = iRJ;
        mindr = dr;
      }
    } // end loop over iRJ

    recoMatches[iTJ] = recoMatch;
  } // end loop over iTJ

  return recoMatches; 
}


/**
 * Inverts the vector storing the index of the reco jet match for each truth jet.
 * If a reco jet has a truth match according to the input array, the truth index is stored at the reco index.
 * Elements are -1 if the reco jet has no truth match.
 */
std::vector <short> GetAktTruthJetMatches (const std::vector <short>& recoMatches, const JetRadius& radius) {
  if (IsCollisions ())
    return std::vector <short> (0);

  std::vector <short> truthMatches (GetAktHIJetN (radius));

  for (short iRJ = 0; iRJ < (short) truthMatches.size (); iRJ++) {
    short iTJ = 0;
    while (iTJ < (short) recoMatches.size () && recoMatches[iTJ] != iRJ) iTJ++;
    truthMatches[iRJ] = (iTJ == (short) recoMatches.size () ? -1 : iTJ);
  }

  return truthMatches;
}



/**
 * Determines the correct truth jet count for this radius.
 * Returns -1 if radius was not recognized.
 */
int GetAktTruthJetN (const JetRadius& radius) {
  switch (radius) {
  case JetRadius::R0p4: {
    assert (akt4_truth_jet_n >= 0);
    return akt4_truth_jet_n;
  }
  case JetRadius::R0p2: {
    assert (akt2_truth_jet_n >= 0);
    return akt2_truth_jet_n;
  }
  default: return -1;
  } // end switch
}



/**
 * Returns the appropriate truth jet pT for the given jet radius.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthJetPt (const int iJ, const JetRadius& radius) {
  assert (GetAktTruthJetN (radius) > iJ);
  switch (radius) {
  case JetRadius::R0p4: return akt4_truth_jet_pt[iJ];
  case JetRadius::R0p2: return akt2_truth_jet_pt[iJ];
  default: return std::nan ("");
  } // end switch
}



/**
 * Returns the appropriate truth jet eta for the given jet radius.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthJetEta (const int iJ, const JetRadius& radius) {
  assert (GetAktTruthJetN (radius) > iJ);
  switch (radius) {
  case JetRadius::R0p4: return akt4_truth_jet_eta[iJ];
  case JetRadius::R0p2: return akt2_truth_jet_eta[iJ];
  default: return std::nan ("");
  } // end switch
}



/**
 * Returns the appropriate truth jet phi for the given jet radius.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthJetPhi (const int iJ, const JetRadius& radius) {
  assert (GetAktTruthJetN (radius) > iJ);
  switch (radius) {
  case JetRadius::R0p4: return akt4_truth_jet_phi[iJ];
  case JetRadius::R0p2: return akt2_truth_jet_phi[iJ];
  default: return std::nan ("");
  } // end switch
}



/**
 * Returns the appropriate truth jet energy for the given jet radius.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthJetEn (const int iJ, const JetRadius& radius) {
  assert (GetAktTruthJetN (radius) > iJ);
  switch (radius) {
  case JetRadius::R0p4: return akt4_truth_jet_e[iJ];
  case JetRadius::R0p2: return akt2_truth_jet_e[iJ];
  default: return std::nan ("");
  } // end switch
}



/**
 * Returns the minimum delta R to any other truth jet with a pT above some threshold.
 * Will cause an error due to NaN returns if the radius is not supported.
 */
float GetAktTruthJetIso (const int iJ, const JetRadius& radius) {
  assert (GetAktTruthJetN (radius) > iJ);

  int jn = GetAktTruthJetN (radius);
  float maxdr = FLT_MAX;
  float jeta = GetAktTruthJetEta (iJ, radius);
  float jphi = GetAktTruthJetPhi (iJ, radius);

  for (int iJp = 0; iJp < jn; iJp++) {
    if (iJp == iJ)
      continue; // skip the jet of interest
    if (GetAktTruthJetPt (iJp, radius) < GetAktTruthIsoMinPtCut (radius))
      continue; // minimum jet pT cut
    float dr = DeltaR (jeta, GetAktTruthJetEta (iJp, radius), jphi, GetAktTruthJetPhi (iJp, radius));
    maxdr = std::fmin (dr, maxdr); // set isolation
  }
  return maxdr;
}



/**
 * Returns the minimum truth jet pT for consideration in isolation calculation.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthIsoMinPtCut (const JetRadius& radius) {
  switch (radius) {
  case JetRadius::R0p4: return akt4_truth_IsoMinPt;
  case JetRadius::R0p2: return akt2_truth_IsoMinPt;
  default: return std::nan ("");
  } // end switch
}


/**
 * Returns the minimum DR for considering a truth jet isolated.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthIsoMinDR (const JetRadius& radius) {
  switch (radius) {
  case JetRadius::R0p4: return akt4_truth_IsoMinDR;
  case JetRadius::R0p2: return akt2_truth_IsoMinDR;
  default: return std::nan ("");
  } // end switch
}


/**
 * Returns the minimum pT for reco. jets in the truth-matching procedure.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthMatchMinRecoPt (const JetRadius& radius) {
  switch (radius) {
  case JetRadius::R0p4: return akt4_hi_TruthMatchMinRecoPt;
  case JetRadius::R0p2: return akt2_hi_TruthMatchMinRecoPt;
  default: return std::nan ("");
  } // end switch
}



///**
// * NOW DEPRECATED
// * Returns the matched HI jet within DR < 1 to this truth jet.
// * Returns -1 if no HI jet is matched within this DR range, or the radius is invalid.
// */
//int GetAktHIJetMatch (const int iTJ, const JetRadius& radius, const int nJESVar) {
//  assert (GetAktTruthJetN (radius) > iTJ);
//  assert (nJESVar >= -1 && nJESVar < nJESSys);
//  float mindr = M_PI;
//  int match = -1;
//  const float tjeta = GetAktTruthJetEta (iTJ, radius);
//  const float tjphi = GetAktTruthJetPhi (iTJ, radius);
//  //const float minpt = GetAktTruthMatchMinRecoPt (radius);
//  const int nJ = GetAktHIJetN (radius);
//  for (int iJ = 0; iJ < nJ; iJ++) {
//    //if (GetAktHIJetPt (iJ, radius, nJESVar) < minpt)
//    //  continue; 
//    const float dr = DeltaR (tjeta, GetAktHIJetEta (iJ, radius, nJESVar), tjphi, GetAktHIJetPhi (iJ, radius, nJESVar));
//    if (dr < mindr) {
//      match = iJ;
//      mindr = dr;
//    }
//  }
//  if (mindr > GetAktTruthMatchMaxDR (radius))
//    return -1;
//  return match;
//}



/**
 * Determines the correct HI jet count for this radius.
 * Returns -1 if radius was not recognized.
 */
int GetAktHIJetN (const JetRadius& radius) {
  switch (radius) {
  case JetRadius::R0p4: {
    assert (akt4_hi_jet_n >= 0);
    return akt4_hi_jet_n;
  }
  case JetRadius::R0p2: {
    assert (akt2_hi_jet_n >= 0);
    return akt2_hi_jet_n;
  }
  default: return -1;
  } // end switch
}



/**
 * Determines the optimal jet pT to return (EtaJES or Cross-calibrated).
 * Jets in data must be cross-calibrated and jets in MC must not be, but jets in MC + data overlay should be cross-calibrated if they are not truth-matched.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetPt (const int iJ, const JetRadius& radius, const int nJESVar, const short scale) {
  assert (GetAktHIJetN (radius) > iJ);
  assert (nJESVar >= -1 && nJESVar < nJESSys);
  assert (scale >= -1 && scale <= 2);
  switch (radius) {
  case JetRadius::R0p4: {
    const float jesVar = 1. + ((nJESVar == -1 || nJESVar == 18 || nJESVar == 19) ? 0 : akt4_hi_jet_pt_sys_JES_ALL[nJESVar][iJ]); // default to no variation, or if flavour uncertainty evaluate in-loop
    if (IsCollisions () || scale == 0)                  return akt4_hi_jet_pt_xcalib[iJ] * jesVar;
    if (!IsDataOverlay () || scale == 1)                return akt4_hi_jet_pt_etajes[iJ] * jesVar;
    if (scale == 2)                                     return akt4_hi_jet_pt_precalib[iJ] * jesVar;
    //if (GetAktTruthJetMatch (iJ, radius, nJESVar) >= 0) return akt4_hi_jet_pt_etajes[iJ] * jesVar;
    else                                                return akt4_hi_jet_pt_xcalib[iJ] * jesVar;
  }
  case JetRadius::R0p2: {
    const float jesVar = 1. + (nJESVar == -1 ? 0 : akt2_hi_jet_pt_sys_JES_ALL[nJESVar][iJ]);
    if (IsCollisions () || scale == 0)                  return akt2_hi_jet_pt_xcalib[iJ] * jesVar;
    if (!IsDataOverlay () || scale == 1)                return akt2_hi_jet_pt_etajes[iJ] * jesVar;
    if (scale == 2)                                     return akt2_hi_jet_pt_precalib[iJ] * jesVar;
    //if (GetAktTruthJetMatch (iJ, radius, nJESVar) >= 0) return akt2_hi_jet_pt_etajes[iJ] * jesVar;
    else                                                return akt2_hi_jet_pt_xcalib[iJ] * jesVar;
  }
  default: 
    return std::nan ("");
  } // end switch
}



/**
 * Determines the optimal jet eta to return (EtaJES or Cross-calibrated).
 * The cross-calibration does nothing to jet eta so this function is trivial.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetEta (const int iJ, const JetRadius& radius, const int nJESVar) {
  assert (GetAktHIJetN (radius) > iJ);
  assert (nJESVar >= -1 && nJESVar < nJESSys);
  switch (radius) {
  case JetRadius::R0p4: return akt4_hi_jet_eta_etajes[iJ];
  case JetRadius::R0p2: return akt2_hi_jet_eta_etajes[iJ];
  default:              return std::nan ("");
  } // end switch
}



/**
 * Determines the optimal jet phi to return.
 * The EtaJES and cross-calibration do nothing to jet phi so this function is trivial.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetPhi (const int iJ, const JetRadius& radius, const int nJESVar) {
  assert (GetAktHIJetN (radius) > iJ);
  assert (nJESVar >= -1 && nJESVar < nJESSys);
  switch (radius) {
  case JetRadius::R0p4: return akt4_hi_jet_phi[iJ];
  case JetRadius::R0p2: return akt2_hi_jet_phi[iJ];
  default:              return std::nan ("");
  } // end switch
}



/**
 * Determines the optimal jet energy to return (EtaJES or Cross-calibrated).
 * Jets in data must be cross-calibrated and jets in MC must not be, but jets in MC + data overlay should be cross-calibrated if they are not truth-matched.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetEn (const int iJ, const JetRadius& radius, const int nJESVar, const short scale) {
  assert (GetAktHIJetN (radius) > iJ);
  assert (nJESVar >= -1 && nJESVar < nJESSys);
  assert (scale >= -1 && scale <= 2);
  switch (radius) {
  case JetRadius::R0p4: {
    const float jesVar = 1. + (nJESVar == -1 ? 0 : akt4_hi_jet_pt_sys_JES_ALL[nJESVar][iJ]);
    if (IsCollisions () || scale == 0)                  return akt4_hi_jet_e_xcalib[iJ] * jesVar;
    if (!IsDataOverlay () || scale == 1)                return akt4_hi_jet_e_etajes[iJ] * jesVar;
    if (scale == 2)                                     return akt4_hi_jet_e_precalib[iJ] * jesVar;
    //if (GetAktTruthJetMatch (iJ, radius, nJESVar) >= 0) return akt4_hi_jet_e_etajes[iJ] * jesVar;
    else                                                return akt4_hi_jet_e_xcalib[iJ] * jesVar;
  }
  case JetRadius::R0p2: {
    const float jesVar = 1. + (nJESVar == -1 ? 0 : akt2_hi_jet_pt_sys_JES_ALL[nJESVar][iJ]);
    if (IsCollisions () || scale == 0)                  return akt2_hi_jet_e_xcalib[iJ] * jesVar;
    if (!IsDataOverlay () || scale == 1)                return akt2_hi_jet_e_etajes[iJ] * jesVar;
    if (scale == 2)                                     return akt2_hi_jet_e_precalib[iJ] * jesVar;
    //if (GetAktTruthJetMatch (iJ, radius, nJESVar) >= 0) return akt2_hi_jet_e_etajes[iJ] * jesVar;
    else                                                return akt2_hi_jet_e_xcalib[iJ] * jesVar;
  }

  default:
    return std::nan ("");
  } // end switch
}



/**
 * Returns the maximum delta R for a truth-reco. jet match.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthMatchMaxDR (const JetRadius& radius) {
  switch (radius) {
  case JetRadius::R0p4: return akt4_TruthMatchMaxDR;
  case JetRadius::R0p2: return akt2_TruthMatchMaxDR;
  default: return std::nan ("");
  } // end switch
}



/**
 * Returns the minimum reconstructed jet pT for consideration in isolation calculation.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIIsoMinPtCut (const JetRadius& radius) {
  switch (radius) {
  case JetRadius::R0p4: return akt4_hi_IsoMinPt;
  case JetRadius::R0p2: return akt2_hi_IsoMinPt;
  default: return std::nan ("");
  } // end switch
}



/**
 * Returns the minimum DR for considering a reco. jet isolated.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIIsoMinDR (const JetRadius& radius) {
  switch (radius) {
  case JetRadius::R0p4: return akt4_hi_IsoMinDR;
  case JetRadius::R0p2: return akt2_hi_IsoMinDR;
  default: return std::nan ("");
  } // end switch
}



/**
 * Returns the minimum delta R to any other reconstructed jet with a pT above some threshold.
 * Will cause an error due to NaN returns if the radius is not supported.
 */
float GetAktHIJetIso (const int iJ, const JetRadius& radius, const int nJESVar) {
  assert (GetAktHIJetN (radius) > iJ);
  assert (nJESVar >= -1 && nJESVar < nJESSys);

  int jn = GetAktHIJetN (radius);
  float maxdr = FLT_MAX;
  float jeta = GetAktHIJetEta (iJ, radius, nJESVar);
  float jphi = GetAktHIJetPhi (iJ, radius, nJESVar);

  for (int iJp = 0; iJp < jn; iJp++) {
    if (iJp == iJ)
      continue; // skip the jet of interest
    if (GetAktHIJetPt (iJp, radius, nJESVar, 2) < GetAktHIIsoMinPtCut (radius))
      continue; // minimum jet pT cut
    float dr = DeltaR (jeta, GetAktHIJetEta (iJp, radius, nJESVar), jphi, GetAktHIJetPhi (iJp, radius, nJESVar));
    maxdr = std::fmin (dr, maxdr); // set isolation
  }
  return maxdr;
}



/**
 * Determines the optimal jet timing to return (depends on jet radius).
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetTiming (const int iJ, const JetRadius& radius) {
  assert (GetAktHIJetN (radius) > iJ);
  switch (radius) {
  case JetRadius::R0p4: return akt4_hi_jet_timing[iJ];
  case JetRadius::R0p2: return akt2_hi_jet_timing[iJ];
  default:              return std::nan ("");
  } // end switch
}



/**
 * Checks the jet cleaning boolean for this jet.
 * Returns false if radius was not recognized.
 */
bool GetAktHIJetCleaning (const int iJ, const JetRadius& radius) {
  assert (GetAktHIJetN (radius) > iJ);
  switch (radius) {
  case JetRadius::R0p4: return akt4_hi_jet_LooseBad[iJ];
  case JetRadius::R0p2: return akt2_hi_jet_LooseBad[iJ];
  default:              return false;
  }
}



/**
 * Returns the appropriate per-jet reweighting factor. Takes in coordinates for an anti-kT HI jet (pT, eta, & phi).
 * Returns 0 if the jet is outside the acceptance.
 */
double GetAktJetWeight (const float jpt, const float jeta, const float jphi, const JetRadius& jetr) {
  //if (UseTruthJets ())
  //  return 1;
  const double accept = ((IspPb () & InDisabledHEC (jeta, jphi, GetRadius (jetr))) || std::fabs (jeta) > 2.8 ? 0. : 1.);
  const double hecwgt = (IspPb () && jeta > 1.1 && jeta < 3.6 ? (2.*M_PI / (3.*M_PI/2. - 2*GetRadius (jetr))) : 1.);
  return accept * hecwgt;
}



/**
 * Returns the Pb-going Q2 vector, or 0 vector if there is no Pb beam.
 * If both beams are Pb, the A and C side values are summed.
 * Optionally will return values from the matched event instead of the trigger event.
 */
QnVector GetPbQ2Vec (const bool getMatching) {
  float q2x = 0, q2y = 0, tot = 0;
  if (IsPbPb () || (IspPb () && IsPeriodA ())) {
    q2x += getMatching ? fcalA_et_Cos2_matching : fcalA_et_Cos2;
    q2y += getMatching ? fcalA_et_Sin2_matching : fcalA_et_Sin2;
    tot += getMatching ? fcalA_et_matching      : fcalA_et;
  }
  if (IsPbPb () || (IspPb () && !IsPeriodA ())) {
    q2x += getMatching ? fcalC_et_Cos2_matching : fcalC_et_Cos2;
    q2y += getMatching ? fcalC_et_Sin2_matching : fcalC_et_Sin2;
    tot += getMatching ? fcalC_et_matching      : fcalC_et;
  }
  if (tot > 0)
    return QnVector (q2x / tot, q2y / tot);
  return QnVector (0, 0);
}



/**
 * Returns the proton-going Q2 vector, or 0 vector if there is no proton beam.
 * If both beams are protons, the A and C side values are summed.
 * Optionally will return values from the matched event instead of the trigger event.
 */
QnVector GetProtonQ2Vec (const bool getMatching) {
  float q2x = 0, q2y = 0, tot = 0;
  if (Ispp () || (IspPb () && !IsPeriodA ())) {
    q2x += getMatching ? fcalA_et_Cos2_matching : fcalA_et_Cos2;
    q2y += getMatching ? fcalA_et_Sin2_matching : fcalA_et_Sin2;
    tot += getMatching ? fcalA_et_matching      : fcalA_et;
  }
  if (Ispp () || (IspPb () && IsPeriodA ())) {
    q2x += getMatching ? fcalC_et_Cos2_matching : fcalC_et_Cos2;
    q2y += getMatching ? fcalC_et_Sin2_matching : fcalC_et_Sin2;
    tot += getMatching ? fcalC_et_matching      : fcalC_et;
  }
  if (tot > 0)
    return QnVector (q2x / tot, q2y / tot);
  return QnVector (0, 0);
}



/**
 * Returns the jet pT weight functions for MC.
 */
TF1** LoadJetPtWgtFuncs () {
  //TDirectory* gdir = gDirectory;

  TString fname = Form ("%s/aux/JetPtWeights.root", workPath.Data ());
  std::cout << "Trying to resolve MC jet pT weights file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  const short nBins = (Ispp () ? 1 : nZdcCentBins+1);

  TF1** farr = new TF1*[nBins];
  for (int iBin = 0; iBin < nBins; iBin++) {

    TF1* f = (TF1*) infile->Get (Form ("f_jet_pt_datamcScaled_ratio_%s_JZ0123", Ispp () ? "ref" : (iBin < nZdcCentBins ? Form ("iCent%i", iBin) : "allCent")))->Clone (Form ("f_jet_pt_weights_%s_JZ0123", Ispp () ? "ref" : (iBin < nZdcCentBins ? Form ("iCent%i", iBin) : "allCent")));

    farr[iBin] = f;

    //f->SetDirectory (gdir);
  }
  std::cout << "Loaded MC jet pT weights, closing file" << std::endl;

  infile->Close ();
  SaferDelete (&infile);

  return farr;
}



/**
 * Returns the jet pT weight histograms for MC.
 */
TH1D** LoadJetPtWgtHists () {
  TDirectory* gdir = gDirectory;

  TString fname = Form ("%s/aux/JetPtWeights.root", workPath.Data ());
  std::cout << "Trying to resolve MC jet pT weights file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  const short nBins = (Ispp () ? 1 : nZdcCentBins+1);

  TH1D** harr = new TH1D*[nBins];
  for (int iBin = 0; iBin < nBins; iBin++) {

    TH1D* h = (TH1D*) infile->Get (Form ("h_jet_pt_datamcScaled_ratio_%s_JZ0123", Ispp () ? "ref" : (iBin < nZdcCentBins ? Form ("iCent%i", iBin) : "allCent")))->Clone (Form ("h_jet_pt_weights_%s_JZ0123", Ispp () ? "ref" : (iBin < nZdcCentBins ? Form ("iCent%i", iBin) : "allCent")));

    harr[iBin] = h;

    h->SetDirectory (gdir);
  }
  std::cout << "Loaded MC jet pT weights, closing file" << std::endl;

  infile->Close ();
  SaferDelete (&infile);

  return harr;
}



/**
 * Returns the tracking efficiency histograms.
 */
TH2D** LoadTrackingEfficiency () {
  TDirectory* gdir = gDirectory;

  //TString fname = Form ("%s/TrackingPerformance/Nominal/outFile.root", rootPath.Data ());
  //std::cout << "Trying to resolve tracking performance file in " << fname.Data () << std::endl;
  //TFile* infile = new TFile (fname, "read");
  TString fname = Form ("%s/aux/TrackingPerformance.root", workPath.Data ());
  std::cout << "Trying to resolve tracking performance file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  const std::string wp = (DoHITightVar () ? "trk_HItight" : (DoHILooseVar () ? "trk_HIloose" : "trk_TightPrimary"));
  //const int iMult = nMultBins-1; // TODO change me for mult. unc.
  const std::string sys = Ispp () ? "pp" : "pPb";
  const int PID = (DoPartSpcVar () ? 211 : 0);

  TH2D** h2arr = new TH2D*[nMultBins];
  for (int iMult = 0; iMult < nMultBins; iMult++) {

    TH2D* h2 = (TH2D*) infile->Get (Form ("h2_efficiency_%s_PID%i_%s_iMult%i", sys.c_str (), PID, wp.c_str (), iMult))->Clone (Form ("h2_tracking_efficiency_iMult%i", iMult));
    h2arr[iMult] = h2;

    h2->SetDirectory (gdir);

    if (DoTrkEffVar ()) {
      for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
        for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) * (1.005 + (0.013/2.35)*(std::fabs (h2->GetXaxis ()->GetBinCenter (iX)) - 0.05))); // uncertainty on tracking efficiency, linear in |eta|
        }
      }
    }
  }
  std::cout << "Loaded tracking efficiencies, closing file" << std::endl;

  infile->Close ();
  SaferDelete (&infile);

  return h2arr;
}



/**
 * Returns the tracking purity histograms (stored as TGAEs).
 */
TGAE** LoadTrackingPurity (const bool useHybridPrimFrac) {
  TString fname = Form ("%s/aux/TrackingPerformance.root", workPath.Data ());
  std::cout << "Trying to resolve tracking performance file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  const std::string wp = (DoHITightVar () ? "trk_HItight" : (DoHILooseVar () ? "trk_HIloose" : "trk_TightPrimary"));
  const int iDR = 3;
  const std::string sys = (Ispp () ? "pp" : (useHybridPrimFrac ? "hybrid" : "pPb"));

  TGAE** g_trk_pur = Get1DArray <TGAE*> (nEtaTrkBins);

  for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
    g_trk_pur[iEta] = (TGAE*) infile->Get (Form ("g_primary_rate_%s_%s_iDR%i_iEta%i", sys.c_str (), wp.c_str (), iDR, iEta))->Clone (Form ("g_primary_rate_iDR%i_iEta%i", iDR, iEta));
    g_trk_pur[iEta]->SetBit (TGraph::kIsSortedX, true); // specifies that x-axis is sorted from low to high; improves Eval function time by allowing binary search
  }

  std::cout << "Loaded tracking purities, closing file" << std::endl;

  infile->Close ();
  SaferDelete (&infile);

  return g_trk_pur;
}



/**
 * Returns array of TGAEs of fits to the tracking purity.
 */
TGAE** LoadTrackingPurityFuncs (const bool useHybridPrimFrac) {
  TString fname = Form ("%s/aux/TrackingPerformance.root", workPath.Data ());
  std::cout << "Trying to resolve tracking performance file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  const std::string wp = (DoHITightVar () ? "trk_HItight" : (DoHILooseVar () ? "trk_HIloose" : "trk_TightPrimary"));
  const int iDR = 3;
  const std::string sys = (Ispp () ? "pp" : (useHybridPrimFrac ? "hybrid" : "pPb"));

  TGAE** gf_trk_pur = Get1DArray <TGAE*> (nEtaTrkBins);

  for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
    TString name = Form ("gf_primary_rate%s_%s_%s_iDR%i_iEta%i", DoFakeRateVar () ? "_fakes_p100" : "", sys.c_str (), wp.c_str (), iDR, iEta);
    if (DoFakeRateVar () && sys == std::string ("hybrid"))
      name = Form ("gf_primary_rate_%s%s_%s_iDR%i_iEta%i", sys.c_str (), DoFakeRateVar () ? "_fakes_p100" : "", wp.c_str (), iDR, iEta); // accidentally switched order of names with hybrid and fake rate plus 100% variation... whoops.
    gf_trk_pur[iEta] = (TGAE*) ((TGAE*) infile->Get (name.Data ()))->Clone (Form ("gf_primary_rate_iDR%i_iEta%i", iDR, iEta));
    assert (gf_trk_pur[iEta] != nullptr);
    gf_trk_pur[iEta]->SetBit (TGraph::kIsSortedX, true); // specifies that x-axis is sorted from low to high; improves Eval function time by allowing binary search
  }

  std::cout << "Loaded tracking purity functions, closing file" << std::endl;

  infile->Close ();
  SaferDelete (&infile);

  return gf_trk_pur;
}



/**
 * Returns the jet energy resolution function for smearing truth jet pT values
 */
TF1* LoadJetEnergyResFunction () {
  TString fname = Form ("%s/aux/JetEnergyResolution.root", workPath.Data ());
  std::cout << "Trying to resolve jet energy resolutions file in " << fname.Data () << std::endl;
  TFile* infile = new TFile (fname, "read");

  const std::string sys = Ispp () ? "pp" : "pPb";
  const short iAllEta = 10;

  TF1* f_jer = (TF1*) ((TF1*) infile->Get (Form ("f_r4_avg_jer_%s_iEta%i", sys.c_str (), iAllEta)))->Clone ("f_r4_avg_jer");
  double xmin, xmax;
  f_jer->GetRange (xmin, xmax);
  f_jer->SetRange (0.001, xmax);

  std::cout << "Loaded JER function" << std::endl;

  infile->Close ();
  SaferDelete (&infile);

  return f_jer;
}



/**
 * Converts a TProfile to a TGraph assuming the x-axis of the TProfile is the y-axis of the TGraph.
 */
TGAE* TProfY2TGAE (TProfile* py) {
  TGAE* g = new TGAE ();
  for (int iX = 1; iX <= py->GetNbinsX (); iX++) {
    g->SetPoint (g->GetN (), py->GetBinContent (iX), py->GetBinCenter (iX));
    g->SetPointEXhigh (g->GetN ()-1, py->GetBinError (iX));
    g->SetPointEXlow (g->GetN ()-1, py->GetBinError (iX));
    g->SetPointEYhigh (g->GetN ()-1, py->GetBinWidth (iX) / 2.);
    g->SetPointEYlow (g->GetN ()-1, py->GetBinWidth (iX) / 2.);
  }
  return g;
}



/**
 * Converts a TProfile to a TGraph assuming the x-axis of the TProfile is the x-axis of the TGraph.
 */
TGAE* TProfX2TGAE (TProfile* px) {
  TGAE* g = new TGAE ();
  for (int iX = 1; iX <= px->GetNbinsX (); iX++) {
    g->SetPoint (g->GetN (), px->GetBinCenter (iX), px->GetBinContent (iX));
    g->SetPointEXhigh (g->GetN ()-1, px->GetBinWidth (iX) / 2.);
    g->SetPointEXlow (g->GetN ()-1, px->GetBinWidth (iX) / 2.);
    g->SetPointEYhigh (g->GetN ()-1, px->GetBinError (iX));
    g->SetPointEYlow (g->GetN ()-1, px->GetBinError (iX));
  }
  return g;
}



/**
 * Sets central values in g according to values in centralValues, but keeps relative uncertainties the same.
 * Values must be positive-definite! I.e. not negative or zero. Otherwise no uncertainty resetting is done for those bins.
 */
void SetCentralValuesKeepRelativeErrors (TGAE* g, TH1D* centralValues) {
  assert (g->GetN () == centralValues->GetNbinsX ());

  double /*xelo, xehi, */yelo, yehi, x, y, ynew;
  for (int i = 0; i < g->GetN (); i++) {
    //xelo = g->GetErrorXlow (i);
    //xehi = g->GetErrorXhigh (i);
    yelo = g->GetErrorYlow (i);
    yehi = g->GetErrorYhigh (i);
    g->GetPoint (i, x, y);

    ynew = centralValues->GetBinContent (i+1);

    if (y > 0) {
      g->SetPoint (i, x, ynew);
      //g->SetPointEXlow (i, xelo);
      //g->SetPointEXhigh (i, xehi);
      g->SetPointEYlow (i, yelo*ynew/y);
      g->SetPointEYhigh (i, yehi*ynew/y);
    }
  }
}



/**
 * Takes the TGAE and sets y --> -y.
 */
void FlipTGAE (TGAE* g) {
  double x, y;
  for (int i = 0; i < g->GetN (); i++) {
    g->GetPoint (i, x, y);
    g->SetPoint (i, x, -y);
  }
  return;
}



///**
// * Performs a bin-by-bin unfold on a TH1D y^meas(x) using a given TF1 with unfolding factors f(x) such that y^unfold(x) = y^meas(x) * f(x).
// */
//void BinByBinUnfold (TH1D* h, TF1* f, const float mult) {
//  for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
//    const float uf = f->Eval (h->GetBinCenter (iX));
//    h->SetBinContent (iX, h->GetBinContent (iX) * uf);
//    h->SetBinError (iX, h->GetBinError (iX) * uf);
//  }
//  return;
//}


/**
 * Multiplies a target histogram by a given TF1 with an optional multiplier on the function.
 */
void MultiplyByTF1 (TH1D* h, TF1* f, const float mult) {
  for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
    //float sf = f->Eval (h->GetBinCenter (iX));
    //if (sf > 1)
    //  sf = 1-mult + mult*sf;//mult * (sf - 1) + 1;
    //else
    //  sf = 1-mult + mult*sf;//1 - mult * (1 - sf);
    const float sf = 1-mult + mult * f->Eval (h->GetBinCenter (iX)); // check: if mult=0.5 and f->Eval()=1.1, then sf=1-0.5+0.5*1.1=0.5+0.55=1.05.
                                                                     // note that if f evals to 0 when it should be ~1, then sf=1-0.5=0.5 incorrectly!
    h->SetBinContent (iX, h->GetBinContent (iX) * sf);
    h->SetBinError (iX, h->GetBinError (iX) * sf);
  }
  return;
}


/**
 * Divides a target histogram by a given TF1 with an optional multiplier on the function.
 */
void DivideByTF1 (TH1D* h, TF1* f, const float mult) {
  for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
    //float sf = f->Eval (h->GetBinCenter (iX));
    //if (sf > 1)
    //  sf = 1-mult + mult*sf;//mult * (sf - 1) + 1;
    //else
    //  sf = 1-mult + mult*sf;//1 - mult * (1 - sf);
    const float sf = 1-mult + mult * f->Eval (h->GetBinCenter (iX)); // check: if mult=0.5 and f->Eval()=1.1, then sf=1-0.5+0.5*1.1=0.5+0.55=1.05.
                                                                     // note that if f evals to 0 when it should be ~1, then sf=1-0.5=0.5 incorrectly!
    h->SetBinContent (iX, h->GetBinContent (iX) / sf);
    h->SetBinError (iX, h->GetBinError (iX) / sf);
  }
  return;
}


/**
 * Divides a histogram by another without propagating uncertainties.
 */
void DivideNoErrors (TH1D* h, const TH1D* hd) {
  assert (hd->GetNbinsX () == h->GetNbinsX ());
  for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
    assert (hd->GetBinCenter (iX) == h->GetBinCenter (iX));
    if (hd->GetBinContent (iX) != 0) {
      h->SetBinContent (iX, h->GetBinContent (iX) / hd->GetBinContent (iX));
      h->SetBinError (iX, h->GetBinError (iX) / hd->GetBinContent (iX));
    }
  }
  return;
}


/**
 * Extension of CalcSystematics (TGAE* sys, TH1D* nom, TH1D* var) for smoothing uncertainties. 
 */
void SmoothSystematics (TGAE* sys, TF1* func, TH1D* nom, TH1D* var) {
  TGraphErrors* tg = new TGraphErrors ();
  double x, y;
  double xlo = DBL_MAX, xhi = DBL_MIN;
  for (int i = 0; i < sys->GetN (); i++) {
    sys->GetPoint (i, x, y);
    xlo = std::fmin (x, xlo);
    xhi = std::fmax (x, xhi);
    if (y != 0) {
      tg->SetPoint (i, x, (sys->GetErrorYhigh (i) - sys->GetErrorYlow (i)) / y);
      double yv = var->GetBinContent (i+1);
      double yve = var->GetBinError (i+1);
      double yn = nom->GetBinContent (i+1);
      double yne = nom->GetBinError (i+1);
      tg->SetPointError (i, 0.5*nom->GetBinWidth (i+1), std::sqrt (std::pow (yv*yne/(yn*yn), 2) + std::pow (yve/yn, 2)));
    } 
  } 
  
  tg->Fit (func, "RN0Q");
  
  for (int i = 0; i < sys->GetN (); i++) {
    sys->GetPoint (i, x, y);
    double newErr = func->Eval (x) * y;
    if (newErr > 0) {
      sys->SetPointEYhigh (i, std::fabs (newErr));
      sys->SetPointEYlow (i, 0);
    } 
    else {
      sys->SetPointEYlow (i, std::fabs (newErr));
      sys->SetPointEYhigh (i, 0);
    } 
  } 
  SaferDelete (&tg);
  return;
}



/**
 * Returns the covariance matrix contained in inFileName.
 * Has dimensions (nPtJBins*nPtChBins)^2
 */
TMatrixD* GetCovarianceMatrix (const TString inFileName) {

  //const int nBins = (nPtJBins+2)*(nPtChBins+2);
  const int nBins = (nPtJBins)*(nPtChBins);

  TMatrixD* cov = new TMatrixD (nBins, nBins);
  ifstream f;
  f.open (inFileName.Data ());
  
  double value;
  int ix = 0;
  int iy = 0;

  f >> value;
  while (f) {
    if (ix >= nBins*nBins || iy >= nBins*nBins) {
      std::cout << "problem at " << ix%nBins << ", " << iy/nBins << std::endl;
      continue;
    }
    (*cov)[ix%nBins][iy/nBins] = value;
    ix++;
    iy++;
    f >> value;
  }

  //cov.Print ();
    
  return cov;
}



/**
 * Returns 2D histogram with relative uncertainty on the flavour fraction.
 */
TH2D* GetFlavorFractionUnc (const JetRadius& r) {
  TFile* inFile = new TFile (Form ("%s/aux/FlavorJESUncertainty_%s%s.root", workPath.Data (), r == JetRadius::R0p4 ? "R0p4" : "R0p2", IspPb () ? "_pPb" : ""), "read");

  TH2D* h2 = (TH2D*) inFile->Get (IspPb () ? "term2" : "termAbs2")->Clone ("h2_flavourFracUnc");
  h2->SetDirectory (0);

  inFile->Close ();
  SaferDelete (&inFile);
  return h2;
}



/**
 * Returns 2D histogram with relative uncertainty on the flavour response.
 */
TH2D* GetFlavorResponseUnc (const JetRadius& r) {
  TFile* inFile = new TFile (Form ("%s/aux/FlavorJESUncertainty_%s%s.root", workPath.Data (), r == JetRadius::R0p4 ? "R0p4" : "R0p2", IspPb () ? "_pPb" : ""), "read");

  TH2D* h2 = (TH2D*) inFile->Get (IspPb () ? "term1" : "termAbs1")->Clone ("h2_flavourRespUnc");
  h2->SetDirectory (0);

  inFile->Close ();
  SaferDelete (&inFile);
  return h2;
}



} // end namespace

#endif
