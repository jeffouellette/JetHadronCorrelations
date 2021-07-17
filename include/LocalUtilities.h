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

enum class SystFlag { Nominal, HITightVar, PionsOnlyVar, WithPileupVar, FcalCentVar, FineFcalCentVar, JetES5PercUpVar, JetES5PercDownVar, JetES5PercSmearVar, JetES2PercUpVar, JetES2PercDownVar, JetES2PercSmearVar }; // types of systematic variations
static const std::vector <SystFlag> AllSystFlag = { SystFlag::Nominal, SystFlag::HITightVar, SystFlag::PionsOnlyVar, SystFlag::WithPileupVar, SystFlag::FcalCentVar, SystFlag::FineFcalCentVar, SystFlag::JetES5PercUpVar, SystFlag::JetES5PercDownVar, SystFlag::JetES5PercSmearVar, SystFlag::JetES2PercUpVar, SystFlag::JetES2PercDownVar, SystFlag::JetES2PercSmearVar };

enum class JetRadius { R0p2, R0p3, R0p4, R0p6, R0p8, R1p0, Invalid };


TString ToTString (const CollisionSystem collSys);
TString ToTString (const DataType dType);
TString ToTString (const TriggerType tType);
TString ToTString (const SystFlag sFlag);
float   GetRadius (const JetRadius r);


bool SetCollisionSystem (TString ts);
bool SetDataType (TString ts);
bool SetTriggerType (TString ts);
bool SetSystFlag (TString ts);


bool IsIons (const CollisionSystem collSys);
bool IsPbPb (const CollisionSystem collSys);
bool IsPbPb18 (const CollisionSystem collSys);
bool IsPbPb15 (const CollisionSystem collSys);
bool IsXeXe (const CollisionSystem collSys);
bool IspPb (const CollisionSystem collSys);
bool IspPb16 (const CollisionSystem collSys);
bool Ispp (const CollisionSystem collSys);
bool Ispp15 (const CollisionSystem collSys);
bool Ispp17 (const CollisionSystem collSys);
bool IsPeriodA (const CollisionSystem collSys);
bool Is5TeV (const CollisionSystem collSys);
bool Is8TeV (const CollisionSystem collSys);
bool Is2018 (const CollisionSystem collSys);
bool Is2017 (const CollisionSystem collSys);
bool Is2016 (const CollisionSystem collSys);
bool Is2015 (const CollisionSystem collSys);

bool IsCollisions (const DataType dType);
bool IsOverlay (const DataType dType);
bool IsDataOverlay (const DataType dType);
bool IsHijing (const DataType dType);

bool UseJetTriggers (const TriggerType tType);
bool UseJet50GeVTriggers (const TriggerType tType);
bool UseJet100GeVTriggers (const TriggerType tType);
bool UseMinBiasTriggers (const TriggerType tType);

bool DoHITightVar (const SystFlag sFlag);
bool DoPionsOnlyVar (const SystFlag sFlag);
bool DoWithPileupVar (const SystFlag sFlag);
bool DoFcalCentVar (const SystFlag sFlag);
bool DoFineFcalCentVar (const SystFlag sFlag);
bool DoJetES5PercUpVar (const SystFlag sFlag);
bool DoJetES5PercDownVar (const SystFlag sFlag);
bool DoJetES5PercSmearVar (const SystFlag sFlag);
bool DoJetES2PercUpVar (const SystFlag sFlag);
bool DoJetES2PercDownVar (const SystFlag sFlag);
bool DoJetES2PercSmearVar (const SystFlag sFlag);



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
bool DoPionsOnlyVar ();
bool DoWithPileupVar ();
bool DoFcalCentVar ();
bool DoFineFcalCentVar ();
bool DoJetES5PercUpVar ();
bool DoJetES5PercDownVar ();
bool DoJetES5PercSmearVar ();
bool DoJetES2PercUpVar ();
bool DoJetES2PercDownVar ();
bool DoJetES2PercSmearVar ();


/**
 * Returns the CoM boost relevant for asymmetric collision systems (i.e. p+Pb). 0 for everything else.
 */
double GetBoost (int rn);


/**
 * Establishes path variables appropriately.
 */
void SetupDirectories (const TString dataSubDir, const bool addSubDir = true);


/**
 * Looks up MC sample cross section, filter efficiency, and number of events.
 */
bool GetMCWeights (TString fname);


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
bool MeetsJetAcceptanceCuts (int iJ, const JetRadius radius = JetRadius::R0p4);


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
int GetAktTruthJetMatch (const int iJ, const JetRadius radius = JetRadius::R0p4);


/**
 * Determines the optimal jet pT to return (EtaJES or Cross-calibrated).
 * Jets in data must be cross-calibrated and jets in MC must not be, but jets in MC + data overlay should be cross-calibrated if they are not truth-matched.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetPt (const int iJ, const JetRadius radius = JetRadius::R0p4);


/**
 * Determines the optimal jet eta to return (EtaJES or Cross-calibrated).
 * The cross-calibration does nothing to jet eta so this function is trivial.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetEta (const int iJ, const JetRadius radius = JetRadius::R0p4);


/**
 * Determines the optimal jet phi to return.
 * The EtaJES and cross-calibration do nothing to jet phi so this function is trivial.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetPhi (const int iJ, const JetRadius radius = JetRadius::R0p4);


/**
 * Determines the optimal jet energy to return (EtaJES or Cross-calibrated).
 * Jets in data must be cross-calibrated and jets in MC must not be, but jets in MC + data overlay should be cross-calibrated if they are not truth-matched.
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetEn (const int iJ, const JetRadius radius = JetRadius::R0p4);


/**
 * Determines the optimal jet timing to return (depends on jet radius).
 * Returns NaN if radius was not recognized.
 */
float GetAktHIJetTiming (const int iJ, const JetRadius radius);


/**
 * Returns the appropriate per-jet reweighting factor. Takes in coordinates for an anti-kT HI jet (pT, eta, & phi).
 * Returns 0 if the jet is outside the acceptance.
 */
double GetAktJetWeight (const float jpt, const float jeta, const float jphi, const JetRadius jetr = JetRadius::R0p4);


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
