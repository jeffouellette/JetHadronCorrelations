#ifndef __LocalUtilities_h__
#define __LocalUtilities_h__

#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>
#include <TMatrixD.h>

#include <vector>
#include <map>

typedef TGraphAsymmErrors TGAE;
typedef std::pair <float, float> QnVector;
//typedef std::pair <float, float> FCalEt;

namespace JetHadronCorrelations {

// run 2 HI data sets
enum class CollisionSystem { pp15, PbPb15, pPb16s5TeV, pPb16, Pbp16, XeXe17, pp17, PbPb18 };
static const std::vector <CollisionSystem> AllCollisionSystem = { CollisionSystem::pp15, CollisionSystem::PbPb15, CollisionSystem::pPb16s5TeV, CollisionSystem::pPb16, CollisionSystem::Pbp16, CollisionSystem::XeXe17, CollisionSystem::pp17, CollisionSystem::PbPb18 };

// data types used in HI
enum class DataType { Collisions, MCSignal, MCDataOverlay, MCHijing, MCHijingOverlay };
static const std::vector <DataType> AllDataType = { DataType::Collisions, DataType::MCSignal, DataType:: MCDataOverlay, DataType::MCHijing, DataType::MCHijingOverlay };

// types of triggers in this analysis
enum class TriggerType { None, J50, J100, MinBias };
static const std::vector <TriggerType> AllTriggerType = { TriggerType::None, TriggerType::J50, TriggerType::J100, TriggerType::MinBias };

// types of systematic variations
enum class SystFlag {
  Nominal,
  HITightVar,       // tracking quality variations
  HILooseVar,
  TrkEffVar,        // tracking efficiency uncertainty
  FakeRateVar,      // fake rate uncertainty
  PrimFitVar,       // primary fraction uncertainty from fit
  JetPrimFracVar,   // use pp primary fraction instead of Hijing in p+Pb signal component
  PartSpcVar,       // particle species uncertainty
  FcalCentVar,      // centrality variations
  FineFcalCentVar,
  MixCatVar1,       // mixing variations
  MixCatVar2,
  MixCatVar3,
  MixCatVar4,       // 2nd order EP matching in pPb
  MixCatVar5,       // 2nd order EP matching in pp
  MixCatVar6,       // less fine fcal binning in peripheral?
  JESVar0,          // jet energy scale variations
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
  JESVar16,         // 2018 Eta-intercalibration non-closure uncertainty  (not applicable)
  JESVar17,         // cross-calibration JES uncertainty
  JESVar18,         // flavor response JES uncertainty
  JESVar19,         // flavor fraction JES uncertainty
  JESVar20,         // R-scan uncertainty                                 (not applicable)
  JERVar0,          // jet energy resolution variations
  JERVar1,
  JERVar2,
  JERVar3,
  JERVar4,
  JERVar5,
  JERVar6,
  JERVar7,
  JERVar8,
  JERVar9,
  JERVar10,
  MCTruthJetsTruthParts,  // does the analysis completely at the MC-truth level
  MCRecoJetsTruthParts,    // selects truth charged particles with reco jets for bin-by-bin unfolding factors (numerator). Denominator is in MCTruthJetsTruthParts
  MCRecoJetsTruthMatchedParts, // does the nominal analysis but only with truth-matched particles (so no background component).
  MCFCalWeighted, // applies MC sum FCal Et weights
};
static const std::vector <SystFlag> AllSystFlag = {
  SystFlag::Nominal,
  SystFlag::HITightVar, 
  SystFlag::HILooseVar,
  SystFlag::TrkEffVar,
  SystFlag::FakeRateVar,
  SystFlag::PrimFitVar,
  SystFlag::JetPrimFracVar,
  SystFlag::PartSpcVar,
  SystFlag::FcalCentVar,
  SystFlag::FineFcalCentVar,
  SystFlag::MixCatVar1,
  SystFlag::MixCatVar2,
  SystFlag::MixCatVar3,
  SystFlag::MixCatVar4,
  SystFlag::MixCatVar5,
  SystFlag::MixCatVar6,
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
  SystFlag::JERVar10,
  SystFlag::MCTruthJetsTruthParts,
  SystFlag::MCRecoJetsTruthParts,
  SystFlag::MCRecoJetsTruthMatchedParts,
  SystFlag::MCFCalWeighted,
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
bool UseJ50Triggers (const TriggerType& tType);
bool UseJ100Triggers (const TriggerType& tType);
bool UseMinBiasTriggers (const TriggerType& tType);

bool DoHITightVar (const SystFlag& sFlag);
bool DoHILooseVar (const SystFlag& sFlag);
bool DoTrkEffVar (const SystFlag& sFlag);
bool DoFakeRateVar (const SystFlag& sFlag);
bool DoPrimFitVar (const SystFlag& sFlag);
bool DoJetPrimFracVar (const SystFlag& sFlag);
bool DoPartSpcVar (const SystFlag& sFlag);
bool DoFcalCentVar (const SystFlag& sFlag);
bool DoFineFcalCentVar (const SystFlag& sFlag);
bool DoMixCatVar1 (const SystFlag& sFlag);
bool DoMixCatVar2 (const SystFlag& sFlag);
bool DoMixCatVar3 (const SystFlag& sFlag);
bool DoMixCatVar4 (const SystFlag& sFlag);
bool DoMixCatVar5 (const SystFlag& sFlag);
bool DoMixCatVar6 (const SystFlag& sFlag);

int GetNJESVar (const SystFlag& sFlag);
int GetNJERVar (const SystFlag& sFlag);

bool DoMCTruthJetsTruthParts (const SystFlag& sFlag);
bool DoMCRecoJetsTruthParts (const SystFlag& sFlag);
bool DoMCRecoJetsTruthMatchedParts (const SystFlag& sFlag);
bool DoMCFCalWeights (const SystFlag& sFlag);



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
bool UseJ50Triggers ();
bool UseJ100Triggers ();
bool UseMinBiasTriggers ();

bool DoHITightVar ();
bool DoHILooseVar ();
bool DoTrkEffVar ();
bool DoFakeRateVar ();
bool DoPrimFitVar ();
bool DoJetPrimFracVar ();
bool DoPartSpcVar ();
bool DoFcalCentVar ();
bool DoFineFcalCentVar ();
bool DoMixCatVar1 ();
bool DoMixCatVar2 ();
bool DoMixCatVar3 ();
bool DoMixCatVar4 ();
bool DoMixCatVar5 ();
bool DoMixCatVar6 ();

int GetNJESVar ();
int GetNJERVar ();

bool DoMCTruthJetsTruthParts ();
bool DoMCRecoJetsTruthParts ();
bool DoMCRecoJetsTruthMatchedParts ();
bool DoMCFCalWeights ();

bool UseTruthJets ();
bool UseTruthParticles ();
bool UseTruthMatchedParticles ();
bool UseMCFCalWeights ();


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
 * Returns the residual JZ scaling factor to smooth the spectrum transition at 60 GeV from SoftQCD to HardQCD
 */
float GetJZScaleFactor (const TString& fname);


/**
 * Returns a copy of the histogram detailing the Zdc cuts.
 */
TH1D* GetZdcCuts ();


/**
 * Returns a copy of the histogram detailing the probability of sampling a given MC event.
 */
TH1D* GetFCalResamplingProbs ();


/**
 * Returns a copy of the histogram with data/MC ratios of the FCal Et distribution
 */
TH1D* GetMCFCalWeights ();


/**
 * Returns the probability histogram of each FCal ET value in 0-20% ZDC events.
 */
TH1D* GetFCalZdcWeights ();


/**
 * Returns a map from event numbers to pure overlay A-side FCal ET values.
 */
std::map <const unsigned int, float>* GetOverlayFCalMap ();


/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* directory, const char* inFileName);


/**
 * Returns the proper jet trigger luminosity for this data set in nb^-1
 */
double GetJetLuminosity ();


/**
 * Returns true if this truth jet passes selection criteria.
 */
bool MeetsTruthJetAcceptanceCuts (int iTJ, const JetRadius& radius);


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
bool MeetsTrackCuts (int iTrk, const int nWPVar = 0);


///**
// * Returns the truth-particle corrected FCal ET values.
// */
//FCalEt GetTruthCorrectedFCal (FCalEt values);


///**
// * NOW DEPRECATED
// * Returns the matched truth jet within DR < 1 to this HI jet.
// * Returns -1 if no truth jet is matched within this DR range, or the radius is invalid.
// */
//int GetAktTruthJetMatch (const int iJ, const JetRadius& radius, const int nJESVar = -1);


/**
 * Returns a vector storing the index of the reco jet match for each truth jet.
 * Elements are -1 if the truth jet has no reco match.
 */
std::vector <short> GetAktRecoJetMatches (const JetRadius& radius, const short nJESVar = -1);


/**
 * Inverts the vector storing the index of the reco jet match for each truth jet.
 * If a reco jet has a truth match according to the input array, the truth index is stored at the reco index.
 * Elements are -1 if the reco jet has no truth match.
 */
std::vector <short> GetAktTruthJetMatches (const std::vector <short>& recoMatches, const JetRadius& radius);


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


///**
// * NOW DEPRECATED
// * Returns the matched HI jet within DR < 1 to this truth jet.
// * Returns -1 if no HI jet is matched within this DR range, or the radius is invalid.
// */
//int GetAktHIJetMatch (const int iTJ, const JetRadius& radius, const int nJESVar = -1);


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
 * Checks the jet cleaning boolean for this jet.
 * Returns false if radius was not recognized.
 */
bool GetAktHIJetCleaning (const int iJ, const JetRadius& radius);


/**
 * Returns the appropriate per-jet reweighting factor. Takes in coordinates for an anti-kT HI jet (pT, eta, & phi).
 * Returns 0 if the jet is outside the acceptance.
 */
double GetAktJetWeight (const float jpt, const float jeta, const float jphi, const JetRadius& jetr);


/**
 * Returns the Pb-going Q2 vector, or 0 vector if there is no Pb beam.
 * If both beams are Pb, the A and C side values are summed.
 * Optionally will return values from the matched event instead of the trigger event.
 */
QnVector GetPbQ2Vec (const bool getMatching = false);


/**
 * Returns the proton-going Q2 vector, or 0 vector if there is no proton beam.
 * If both beams are protons, the A and C side values are summed.
 * Optionally will return values from the matched event instead of the trigger event.
 */
QnVector GetProtonQ2Vec (const bool getMatching = false);


/**
 * Returns the jet pT weight functions for MC.
 */
TF1** LoadJetPtWgtFuncs ();


/**
 * Returns the jet pT weight histograms for MC.
 */
TH1D** LoadJetPtWgtHists ();


/**
 * Returns the tracking efficiency histograms.
 */
TH2D** LoadTrackingEfficiency ();


/**
 * Returns the tracking purity histograms (stored as TGAEs).
 */
TGAE** LoadTrackingPurity (const bool useHybridPrimFrac);


/**
 * Returns array of TGAEs of fits to the tracking purity.
 */
TGAE** LoadTrackingPurityFuncs (const bool useHybridPrimFrac);


/**
 * Returns the jet energy resolution function for smearing truth jet pT values
 */
TF1* LoadJetEnergyResFunction ();


/**
 * Converts a TProfile to a TGraph assuming the x-axis of the TProfile is the y-axis of the TGraph.
 */
TGAE* TProfY2TGAE (TProfile* py);


/**
 * Converts a TProfile to a TGraph assuming the x-axis of the TProfile is the x-axis of the TGraph.
 */
TGAE* TProfX2TGAE (TProfile* px);


/**
 * Sets central values in g according to values in centralValues, but keeps relative uncertainties the same.
 * Values must be positive-definite! I.e. not negative or zero. Otherwise no uncertainty resetting is done for those bins.
 */
void SetCentralValuesKeepRelativeErrors (TGAE* g, TH1D* centralValues);


/**
 * Takes the TGAE and sets y --> -y.
 */
void FlipTGAE (TGAE* g);


///**
// * Performs a bin-by-bin unfold on a TH1D y^meas(x) using a given TF1 with unfolding factors f(x) such that y^unfold(x) = y^meas(x) * f(x).
// */
//void BinByBinUnfold (TH1D* h, TF1* f);


/**
 * Multiplies a target histogram by a given TF1 with an optional multiplier on the function.
 */
void MultiplyByTF1 (TH1D* h, TF1* f, const float mult = 1.);


/**
 * Divides a target histogram by a given TF1 with an optional multiplier on the function.
 */
void DivideByTF1 (TH1D* h, TF1* f, const float mult = 1.);


/**
 * Divides a histogram by another without propagating uncertainties.
 */
void DivideNoErrors (TH1D* h, const TH1D* hd);


/**
 * Extension of CalcSystematics (TGAE* sys, TH1D* nom, TH1D* var) for smoothing uncertainties. 
 */
void SmoothSystematics (TGAE* sys, TF1* func, TH1D* nom, TH1D* var);


/**
 * Returns the covariance matrix contained in inFileName.
 * Has dimensions (nPtJBins*nPtChBins)^2
 */
TMatrixD GetCovarianceMatrix (const TString inFileName);


/**
 * Returns 2D histogram with relative uncertainty on the flavour fraction.
 */
TH2D* GetFlavorFractionUnc (const JetRadius& r);


/**
 * Returns 2D histogram with relative uncertainty on the flavour response.
 */
TH2D* GetFlavorResponseUnc (const JetRadius& r);


} // end namespace

#endif
