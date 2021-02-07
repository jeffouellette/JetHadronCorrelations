#ifndef __LocalUtilities_h__
#define __LocalUtilities_h__


#include <TString.h>
#include <TH1D.h>
#include <TFile.h>

namespace HadronYieldsAnalysis {

enum class CollisionSystem { pp15, PbPb15, pPb16s5TeV, pPb16, Pbp16, XeXe17, pp17, PbPb18 }; // run 2 HI data sets
enum class DataType { Collisions, MCSignal, MCDataOverlay, MCHijing, MCHijingOverlay }; // data types used in HI
enum class TriggerType { Jet, MinBias }; // types of triggers in this analysis

TString ToTString (const CollisionSystem collSys);
TString ToTString (const DataType dType);
TString ToTString (const TriggerType tType);


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
bool UseMinBiasTriggers (const TriggerType tType);


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
bool UseMinBiasTriggers ();


/**
 * Returns the CoM boost relevant for asymmetric collision systems (i.e. p+Pb). 0 for everything else.
 */
double GetBoost (int rn);


/**
 * Establishes path variables appropriately.
 */
void SetupDirectories (const TString dataSubDir, const bool addSubDir = true);


/**
 * Returns a copy of the histogram detailing the Zdc cuts.
 */
TH1D* GetZdcCuts ();


/**
 * Returns the appropriate file in the given directory.
 * For MC, inFileName MUST be specified.
 */
TFile* GetFile (const char* directory, const int dataSet, const char* inFileName);


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
bool MeetsJetCuts (int iJ);


/**
 * Returns true if this track passes selection criteria.
 */
bool MeetsTrackCuts (int iTrk);

} // end namespace

#endif
