#ifndef __Params_h__
#define __Params_h__

#include "LocalUtilities.h"

#include <Utilities.h>

#include <TString.h>

#include <fstream>
#include <string>
#include <set>
#include <math.h>

typedef std::pair <float, float> QnVector;

namespace JetHadronCorrelations { 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global variable declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const double electron_mass = 0.000511;
const double muon_mass = 0.105658;


const float min_trk_pt = 0.4;

const float akt2_TruthMatchMaxDR = 0.15; // Truth-matching maximum dR for R=0.2 jets
const float akt4_TruthMatchMaxDR = 0.3; // Truth-matching maximum dR for R=0.4 jets
const float akt2_hi_TruthMatchMinRecoPt = 0; // Minimum HI jet pT for reco. matching of R=0.2 truth jets
const float akt4_hi_TruthMatchMinRecoPt = 0; // Minimum HI jet pT for reco. matching of R=0.4 truth jets
const float akt2_hi_IsoMinPt = 7; // Minimum pT for isolation calculation of R=0.2 HI jets
const float akt4_hi_IsoMinPt = 7; // Minimum pT for isolation calculation of R=0.4 HI jets
const float akt2_truth_IsoMinPt = 7; // Minimum pT for isolation calculation of R=0.2 truth jets
const float akt4_truth_IsoMinPt = 7; // Minimum pT for isolation calculation of R=0.4 truth jets
const float akt2_hi_IsoMinDR = 0.0; // Minimum isolation DR for R=0.2 HI jets
const float akt4_hi_IsoMinDR = 0.0; // Minimum isolation DR for R=0.4 HI jets
const float akt2_truth_IsoMinDR = 0.5; // Minimum isolation DR for R=0.2 truth jets
const float akt4_truth_IsoMinDR = 1.0; // Minimum isolation DR for R=0.4 truth jets

const float nJERSigma = 3; // Number of sigma on the JER to cut at in MC

extern TString workPath;
extern TString extWorkPath;
extern TString rootPath;
extern TString dataPath;

extern float crossSectionPicoBarns;
extern float mcFilterEfficiency;
extern int mcNumberEvents;

// systematics configuration variables
extern bool doHITightVar;
extern bool doPionsOnlyVar;
extern bool doWithPileupVar;
extern bool doFcalCentVar;
extern bool doJetES5PercUpVar;
extern bool doJetES5PercDownVar;
extern bool doJetES5PercSmearVar;
extern bool doJetES2PercUpVar;
extern bool doJetES2PercDownVar;
extern bool doJetES2PercSmearVar;

extern CollisionSystem collisionSystem;
extern DataType dataType; 
extern TriggerType triggerType;


// histogram level jet pT cuts
extern float jet_min_pt;
extern float jet_max_pt;

// list of track working point names
const std::vector <std::string> trackWPNames = {"trk_TightPrimary", "trk_HItight", "trk_HIloose"};

TString GetJetPtStr (const char* tag) {
  if (strcmp (tag, "15GeVJets") == 0)
    return "15 GeV";
  if (strcmp (tag, "30GeVJets") == 0)
    return "30 GeV";
  if (strcmp (tag, "60GeVJets") == 0)
    return "60 GeV";
  if (strcmp (tag, "120GeVJets") == 0)
    return "120 GeV";
  return "??? GeV";
}


const std::vector <TString> directions = {"ns", "perp", "as"};
const int nDir = (int) directions.size ();

const std::vector <TString> pTChSelections = {"gt0p5_lt1", "gt1_lt1p5", "gt1p5_lt2", "gt2_lt4", "gt4_lt6", "gt6_lt8", "gt8_lt10", "gt10_lt15", "gt15_lt20", "gt20_lt30"};
const int nPtChSelections = pTChSelections.size ();

std::map <TString, TString> pTChStrs = {
  {"gt0p5_lt1", "#it{p}_{T}^{ch} = 0.5-1 GeV"},
  {"gt1_lt1p5", "#it{p}_{T}^{ch} = 1-1.5 GeV"},
  {"gt1p5_lt2", "#it{p}_{T}^{ch} = 1.5-2 GeV"},
  {"gt2_lt4", "#it{p}_{T}^{ch} = 2-4 GeV"},
  {"gt4_lt6", "#it{p}_{T}^{ch} = 4-6 GeV"},
  {"gt6_lt8", "#it{p}_{T}^{ch} = 6-8 GeV"},
  {"gt8_lt10", "#it{p}_{T}^{ch} = 8-10 GeV"},
  {"gt10_lt15", "#it{p}_{T}^{ch} = 10-15 GeV"},
  {"gt15_lt20", "#it{p}_{T}^{ch} = 15-20 GeV"},
  {"gt20_lt30", "#it{p}_{T}^{ch} = 20-30 GeV"}
};

std::map <TString, std::pair <double, double>> pTChStrCuts = {
  {"gt0p5_lt1", std::pair <double, double> (0.5, 1)},
  {"gt1_lt1p5", std::pair <double, double> (1, 1.5)},
  {"gt1p5_lt2", std::pair <double, double> (1.5, 2)},
  {"gt2_lt4", std::pair <double, double> (2, 4)},
  {"gt4_lt6", std::pair <double, double> (4, 6)},
  {"gt6_lt8", std::pair <double, double> (6, 8)},
  {"gt8_lt10", std::pair <double, double> (8, 10)},
  {"gt10_lt15", std::pair <double, double> (10, 15)},
  {"gt15_lt20", std::pair <double, double> (15, 20)},
  {"gt20_lt30", std::pair <double, double> (20, 30)}
};


const int pPbRuns[] = {312796, 312837, 312937, 312945, 312968, 314199};
const short npPbRuns = sizeof (pPbRuns) / sizeof(pPbRuns[0]);

const int ppRuns[] = {340644, 340683, 340697, 340718, 340814, 340849, 340850, 340910, 340925, 340973, 341027, 341123, 341184};
const short nppRuns = sizeof (ppRuns) / sizeof (ppRuns[0]);

const double multBins[] = {-0.5, 60.5, 120.5, 999.5};
const short nMultBins = sizeof (multBins) / sizeof (multBins[0]); // last bin is inclusive in multiplicity

const double drBins[] = {0, 0.2, 0.4, M_PI};
const short nDRBins = sizeof (drBins) / sizeof (drBins[0]); // last bin is inclusive in dR

const double etaTrkBins[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const short nEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

// position in array       0                                                11                          18                                                                                    36                  40
const double pTJBins[] = {10, 11, 12, 13, 14, 15, 17.5, 20, 22.5, 25, 27.5, 30, 33, 36, 40, 45, 50, 55, 60, 65, 70, 75, 82.5, 90, 100, 110, 120, 130, 145, 160, 180, 200, 220, 240, 260, 280, 300, 325, 350, 375, 400};
//double pTJBins[] = {7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70, 72.5, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 240, 260, 280, 300, 325, 350, 375, 400};
//double pTJBins[] = {15, 20, 30, 45, 60, 90, 120, 160, 200, 240, 300, 350, 400};
const short nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;

short GetPtJBin (const float jpt) {
  if (jpt < pTJBins[0])
    return -1;
  short iPtJ = 0;
  while (iPtJ < nPtJBins) {
    if (jpt < pTJBins[iPtJ+1])
      break;
    else
      iPtJ++;
  }
  return iPtJ;
}


const double dPhiBins[] = {0, M_PI/24., M_PI/12., M_PI/8., M_PI/6., 5*M_PI/24., M_PI/4., 7*M_PI/24., M_PI/3., 9*M_PI/24., 5*M_PI/12., 11*M_PI/24., M_PI/2., 13*M_PI/24., 7*M_PI/12., 5*M_PI/8., 2*M_PI/3., 17*M_PI/24., 3*M_PI/4., 19*M_PI/24., 5*M_PI/6., 7*M_PI/8., 11*M_PI/12., 23*M_PI/24., M_PI};
const short nDPhiBins = sizeof (dPhiBins) / sizeof (dPhiBins[0]) - 1;
short GetDPhiBin (const float dphi) {
  short iDPhi = 0;
  while (iDPhi < nDPhiBins) {
    if (dphi < dPhiBins[iDPhi+1])
      break;
    else
      iDPhi++;
  }
  return iDPhi;
}


//const double pTChBins[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60};//, 70, 80, 90, 100, 120}; // old binning
const double pTChBins[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 50, 60, 75, 90, 120}; // new binning
const short nPtChBins = sizeof (pTChBins) / sizeof (pTChBins[0]) - 1;
short GetPtChBin (const float ptch) {
  if (ptch < pTChBins[0])
    return -1;
  short iPtCh = 0;
  while (iPtCh < nPtChBins) {
    if (ptch < pTChBins[iPtCh+1])
      break;
    else
      iPtCh++;
  }
  return iPtCh;
}


} // end namespace

#endif
