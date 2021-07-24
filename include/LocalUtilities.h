#ifndef __LocalUtilities_h__
#define __LocalUtilities_h__

#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TProfile.h>

#include <vector>

namespace JetHadronCorrelations {

enum class CollisionSystem { pp15, PbPb15, pPb16s5TeV, pPb16, Pbp16, XeXe17, pp17, PbPb18 }; // run 2 HI data sets
static const std::vector <CollisionSystem> AllCollisionSystem = { CollisionSystem::pp15, CollisionSystem::PbPb15, CollisionSystem::pPb16s5TeV, CollisionSystem::pPb16, CollisionSystem::Pbp16, CollisionSystem::XeXe17, CollisionSystem::pp17, CollisionSystem::PbPb18 };

enum class DataType { Collisions, MCSignal, MCDataOverlay, MCHijing, MCHijingOverlay }; // data types used in HI
static const std::vector <DataType> AllDataType = { DataType::Collisions, DataType::MCSignal, DataType:: MCDataOverlay, DataType::MCHijing, DataType::MCHijingOverlay };

enum class TriggerType { None, Jet50GeV, Jet100GeV, MinBias }; // types of triggers in this analysis
static const std::vector <TriggerType> AllTriggerType = { TriggerType::None, TriggerType::Jet50GeV, TriggerType::Jet100GeV, TriggerType::MinBias };

enum class SystFlag {
  Nominal,
  HITightVar,
  HILooseVar,
  FcalCentVar,
  FineFcalCentVar,
  JESVar0,
  JESVar1,
  JESVar2,
  JESVar3,
  JESVar4,
  JESVar5,
  JESVar6,
  JESVar7,
  JESVar8,
  JESVar9,
  JESVar10,
  JESVar11,
  JESVar12,
  JESVar13,
  JESVar14,
  JESVar15,
  JESVar16,
  JESVar17,
  JESVar18,
  JESVar19,
  JESVar20,
  JERVar0,
  JERVar1,
  JERVar2,
  JERVar3,
  JERVar4,
  JERVar5,
  JERVar6,
  JERVar7,
  JERVar8,
  JERVar9,
  JERVar10
}; // types of systematic variations
static const std::vector <SystFlag> AllSystFlag = {
  SystFlag::Nominal,
  SystFlag::HITightVar, 
  SystFlag::HILooseVar,
  SystFlag::FcalCentVar,
  SystFlag::FineFcalCentVar,
  SystFlag::JESVar0,
  SystFlag::JESVar1,
  SystFlag::JESVar2,
  SystFlag::JESVar3,
  SystFlag::JESVar4,
  SystFlag::JESVar5,
  SystFlag::JESVar6,
  SystFlag::JESVar7,
  SystFlag::JESVar8,
  SystFlag::JESVar9,
  SystFlag::JESVar10,
  SystFlag::JESVar11,
  SystFlag::JESVar12,
  SystFlag::JESVar13,
  SystFlag::JESVar14,
  SystFlag::JESVar15,
  SystFlag::JESVar16,
  SystFlag::JESVar17,
  SystFlag::JESVar18,
  SystFlag::JESVar19,
  SystFlag::JESVar20,
  SystFlag::JERVar0,
  SystFlag::JERVar1,
  SystFlag::JERVar2,
  SystFlag::JERVar3,
  SystFlag::JERVar4,
  SystFlag::JERVar5,
  SystFlag::JERVar6,
  SystFlag::JERVar7,
  SystFlag::JERVar8,
  SystFlag::JERVar9,
  SystFlag::JERVar10
};

enum class JetRadius { R0p2, R0p3, R0p4, R0p6, R0p8, R1p0, Invalid };


TString ToTString (const CollisionSystem& collSys);
TString ToTString (const DataType& dType);
TString ToTString (const TriggerType& tType);
TString ToTString (const SystFlag& sFlag);
float   GetRadius (const JetRadius& r);


bool SetCollisionSystem (TString ts);
bool SetDataType (TString ts);
bool SetTriggerType (TString ts);
bool SetSystFlag (TString ts);


bool IsIons (const CollisionSystem& collSys);
bool IsPbPb (const CollisionSystem& collSys);
bool IsPbPb18 (const CollisionSystem& collSys);
bool IsPbPb15 (const CollisionSystem& collSys);
bool IsXeXe (const CollisionSystem& collSys);
bool IspPb (const CollisionSystem& collSys);
bool IspPb16 (const CollisionSystem& collSys);
bool Ispp (const CollisionSystem& collSys);
bool Ispp15 (const CollisionSystem& collSys);
bool Ispp17 (const CollisionSystem& collSys);
bool IsPeriodA (const CollisionSystem& collSys);
bool Is5TeV (const CollisionSystem& collSys);
bool Is8TeV (const CollisionSystem& collSys);
bool Is2018 (const CollisionSystem& collSys);
bool Is2017 (const CollisionSystem& collSys);
bool Is2016 (const CollisionSystem& collSys);
bool Is2015 (const CollisionSystem& collSys);

bool IsCollisions (const DataType& dType);
bool IsOverlay (const DataType& dType);
bool IsDataOverlay (const DataType& dType);
bool IsHijing (const DataType& dType);

bool UseJetTriggers (const TriggerType& tType);
bool UseJet50GeVTriggers (const TriggerType& tType);
bool UseJet100GeVTriggers (const TriggerType& tType);
bool UseMinBiasTriggers (const TriggerType& tType);

bool DoHITightVar (const SystFlag& sFlag);
bool DoHILooseVar (const SystFlag& sFlag);
bool DoFcalCentVar (const SystFlag& sFlag);
bool DoFineFcalCentVar (const SystFlag& sFlag);

int GetNJESVar (const SystFlag& sFlag);
int GetNJERVar (const SystFlag& sFlag);



bool IsIons ();
bool IsPbPb ();
bool IsPbPb18 ();
bool IsPbPb15 ();
bool IsXeXe ();
bool IspPb ();
bool IspPb16 ();
bool Ispp ();
bool Ispp15 ();
bool Ispp17 ();
bool IsPeriodA ();
bool Is5TeV ();
bool Is8TeV ();
bool Is2018 ();
bool Is2017 ();
bool Is2016 ();
bool Is2015 ();

bool IsCollisions ();
bool IsOverlay ();
bool IsDataOverlay ();
bool IsHijing ();

bool UseJetTriggers ();
bool UseJet50GeVTriggers ();
bool UseJet100GeVTriggers ();
bool UseMinBiasTriggers ();

bool DoHITightVar ();
bool DoHILooseVar ();
bool DoFcalCentVar ();
bool DoFineFcalCentVar ();

int GetNJESVar ();
int GetNJERVar ();


/**
 * Returns the CoM boost relevant for asymmetric collision systems (i.e. p+Pb). 0 for everything else.
 */
double GetBoost (int rn);


/**
 * Establishes path variables appropriately.
 */
void SetupDirectories (const TString& dataSubDir, const bool addSubDir = true);


/**
 * Looks up MC sample cross section, filter efficiency, and number of events.
 */
bool GetMCWeights (const TString& fname);


/**
 * Returns the minimum anti-kT R=0.4 truth jet pT for this JZXR04 MC sample.
 */
float GetJZXR04MinPt (const TString& fname);


/**
 * Returns the maximum anti-kT R=0.4 truth jet pT for this JZXR04 MC sample.
 */
float GetJZXR04MaxPt (const TString& fname);


/**
 * Returns a copy of the histogram detailing the Zdc cuts.
 */
TH1D* GetZdcCuts ();


/**
 * Returns the probability histogram of each FCal ET value in 0-20% ZDC events.
 */
TH1D* GetFCalZdcWeights ();


/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* directory, const char* inFileName);


/**
 * Returns the proper jet trigger luminosity for this data set in nb^-1
 */
double GetJetLuminosity ();


/**
 * Returns true if this jet passes selection criteria.
 */
bool MeetsJetAcceptanceCuts (int iJ, const JetRadius& radius, const int nJESVar = -1);


/**
 * Returns true if this jet pT passes pT cuts.
 */
bool MeetsJetPtCut (double jpt);


/**
 * Returns true if this track passes selection criteria.
 */
bool MeetsTrackCuts (int iTrk);


/**
 * Returns the matched truth jet within DR < 1 to this HI jet.
 * Returns -1 if no truth jet is matched within this DR range, or the radius is invalid.
 */
int GetAktTruthJetMatch (const int iJ, const JetRadius& radius, const int nJESVar = -1);


/**
 * Determines the correct truth jet count for this radius.
 * Returns -1 if radius was not recognized.
 */
int GetAktTruthJetN (const JetRadius& radius);


/**
 * Returns the appropriate truth jet pT for the given jet radius.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthJetPt (const int iJ, const JetRadius& radius);


/**
 * Returns the appropriate truth jet eta for the given jet radius.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthJetEta (const int iJ, const JetRadius& radius);


/**
 * Returns the appropriate truth jet phi for the given jet radius.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthJetPhi (const int iJ, const JetRadius& radius);


/**
 * Returns the appropriate truth jet energy for the given jet radius.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthJetEn (const int iJ, const JetRadius& radius);


/**
 * Returns the minimum delta R to any other truth jet with a pT above some threshold.
 * Will cause an error due to NaN returns if the radius is not supported.
 */
float GetAktTruthJetIso (const int iJ, const JetRadius& radius);


/**
 * Returns the minimum truth jet pT for consideration in isolation calculation.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthIsoMinPtCut (const JetRadius& radius);


/**
 * Returns the minimum DR for considering a truth jet isolated.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthIsoMinDR (const JetRadius& radius);


/**
 * Returns the minimum pT for reco. jets in the truth-matching procedure.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthMatchMinRecoPt (const JetRadius& radius);


/**
 * Returns the matched HI jet within DR < 1 to this truth jet.
 * Returns -1 if no HI jet is matched within this DR range, or the radius is invalid.
 */
int GetAktHIJetMatch (const int iTJ, const JetRadius& radius, const int nJESVar = -1);


/**
 * Determines the correct HI jet count for this radius.
 * Returns -1 if radius was not recognized.
 */
int GetAktHIJetN (const JetRadius& radius);


/**
 * Determines the optimal jet pT to return (EtaJES or Cross-calibrated).
 * Jets in data must be cross-calibrated and jets in MC must not be, but jets in MC + data overlay should be cross-calibrated if they are not truth-matched.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetPt (const int iJ, const JetRadius& radius, const int nJESVar = -1, const short scale = -1);


/**
 * Determines the optimal jet eta to return (EtaJES or Cross-calibrated).
 * The cross-calibration does nothing to jet eta so this function is trivial.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetEta (const int iJ, const JetRadius& radius, const int nJESVar = -1);


/**
 * Determines the optimal jet phi to return.
 * The EtaJES and cross-calibration do nothing to jet phi so this function is trivial.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetPhi (const int iJ, const JetRadius& radius, const int nJESVar = -1);


/**
 * Determines the optimal jet energy to return (EtaJES or Cross-calibrated).
 * Jets in data must be cross-calibrated and jets in MC must not be, but jets in MC + data overlay should be cross-calibrated if they are not truth-matched.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetEn (const int iJ, const JetRadius& radius, const int nJESVar = -1, const short scale = -1);


/**
 * Returns the maximum delta R for a truth-reco. jet match.
 * Returns NaN if radius was not recognized.
 */
float GetAktTruthMatchMaxDR (const JetRadius& radius);


/**
 * Returns the minimum reconstructed jet pT for consideration in isolation calculation.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIIsoMinPtCut (const JetRadius& radius);


/**
 * Returns the minimum DR for considering a reco. jet isolated.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIIsoMinDR (const JetRadius& radius);


/**
 * Returns the minimum delta R to any other reconstructed jet with a pT above some threshold.
 * Will cause an error due to NaN returns if the radius is not supported.
 */
float GetAktHIJetIso (const int iJ, const JetRadius& radius, const int nJESVar = -1);


/**
 * Determines the optimal jet timing to return (depends on jet radius).
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetTiming (const int iJ, const JetRadius& radius);


/**
 * Returns the appropriate per-jet reweighting factor. Takes in coordinates for an anti-kT HI jet (pT, eta, & phi).
 * Returns 0 if the jet is outside the acceptance.
 */
double GetAktJetWeight (const float jpt, const float jeta, const float jphi, const JetRadius& jetr);


/**
 * Returns the tracking efficiency histograms.
 */
TH2D* LoadTrackingEfficiency ();


/**
 * Returns the tracking purity histograms.
 */
TH2D* LoadTrackingPurity ();


/**
 * Returns array of functions that fit the tracking purity.
 */
TF1** LoadTrackingPurityFuncs ();


/**
 * Converts a TProfile to a TGraph assuming the x-axis of the TProfile is the y-axis of the TGraph.
 */
TGraphErrors* TProfY2TGE (TProfile* py);


/**
 * Converts a TProfile to a TGraph assuming the x-axis of the TProfile is the x-axis of the TGraph.
 */
TGraphErrors* TProfX2TGE (TProfile* px);


} // end namespace

#endif
