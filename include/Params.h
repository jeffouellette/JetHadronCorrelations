#ifndef __Params_h__
#define __Params_h__

#include "LocalUtilities.h"

#include <Utilities.h>

#include <TString.h>

#include <string>
#include <set>
#include <math.h>

using namespace std;

namespace JetHadronCorrelations { 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global variable declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const double pi = atan (1)*4;

const double electron_mass = 0.000511;
const double muon_mass = 0.105658;


const float min_trk_pt = 0.5;
const float min_akt4_hi_jet_pt = 30;

TString workPath = TString (std::getenv ("JETHADRONCORR_PATH")) + "/Analysis/";
TString extWorkPath = TString (std::getenv ("JETHADRONCORR_DATA_PATH")) + "/";
TString rootPath = extWorkPath + "rootFiles/";
TString dataPath = extWorkPath + "data/";
TString plotPath = workPath + "Plots/"; 

extern float crossSectionPicoBarns;
extern float mcFilterEfficiency;
extern int mcNumberEvents;

// systematics configuration variables
extern bool doHITightVar;
extern bool doPionsOnlyVar;
extern bool doWithPileupVar;
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
extern double jet_min_pt;
extern double jet_max_pt;


// Centrality classes for mixing events
// See centralities defined in Kurt's note: https://cds.cern.ch/record/2301540/files/ATL-COM-DAPR-2018-002.pdf
const float CentBins[12] = {
 -1000,     // 100%
     3.54,  // 90%
     7.60,  // 80%
    12.27,  // 70%
    17.51,  // 60%
    23.45,  // 50%
    30.29,  // 40%
    38.49,  // 30%
    49.03,  // 20%
    65.08,  // 10%
    79.45,  // 5%
   109.01   // 1%
};
const int numCentBins = sizeof (CentBins) / sizeof (CentBins[0]) - 1; // no bin needed for pp

/**
 * Returns the bin corresponding to this sum fcal et bin.
 * Returns -1 for >90% peripheral collisions (i.e. FCal Et is less than energy for 90% centrality).
 */
short GetCentBin (const float fcal_et) {
  if (fcal_et < CentBins[0])
    return -1;
  short i = 0;
  while (i < numCentBins) {
    i++;
    if (fcal_et < CentBins[i])
      break;
  }
  return i-1;
}


const double zdcCentBins[3] = {
/* UNCALIBRATED ADC VALUES
   0,         // 100%
   0.725967,  // 90%
   2.94087,   // 80%
   5.22591,   // 70%
   7.12421,   // 60%
   8.7794,    // 50%
  10.2273,    // 40%
  11.5395,    // 30%
  12.8396,    // 20%
  14.4122,    // 10%
  //15.6168,    // 5%
  //16.9208,    // 2%
  //17.7623,    // 1%
  //18.5288,    // 0.5%
  //19.4441,    // 0.2%
  //20.051,     // 0.1%
  26.95       // 0%
*/

/* CALIBRATED ENERGIES (over all runs) */
   0,         // "100%"
 //  3.62017,   // 90%
 // 14.9817,    // 80%
 // 26.5908,    // 70%
 // 36.3317,    // 60%
 // 44.8672,    // 50%
 // 52.3714,    // 40%
 // 59.1911,    // 30%
  65.9228,    // 20%
 // 73.9471,    // 10%
 // 79.9879,    // 5%
 // 86.4855,    // 2%
 // 90.7422,    // 1%
 // 94.6261,    // 0.5%
 // 99.3721,    // 0.2%
 //102.71,      // 0.1%
 140          // "0%"
};
const int numZdcCentBins = sizeof (zdcCentBins) / sizeof (zdcCentBins[0]) - 1; // no bin needed for pp


short GetZdcCentBin (const float zdc_E) {
  if (zdc_E < zdcCentBins[0])
    return -1;
  short i = 0;
  while (i < numZdcCentBins) {
    i++;
    if (zdc_E < zdcCentBins[i])
      break;
  }
  return i-1;
}


// Run vs approx. luminosity in nb^-1 for 2016 p+Pb data at 8.16 TeV
// p+Pb (period A) 313063 313067 313100 313107 313136 313187 313259 313285 313295 313333 313435
//                 0.03   1.24   9.66   11.92  10.40  3.67   5.12   4.74   10.69  4.13   0.39
// Pb+p (period B) 313572 313574 313575 313603 313629 313630 313688 313695 313833 313878 313929 313935 313984 314014 314077 314105 314112 314157 314170
//                 0.01   1.33   7.54   8.69   6.86   7.90   7.96   4.53   5.11   2.16   0.63   10.96  2.40   7.36   10.19  6.50   10.49  9.83   4.92

// 
const set<int> groupA = {313063, 313067, 313100, 313107};
const set<int> groupB = {313136, 313187, 313259, 313285};
const set<int> groupC = {313295, 313333, 313435};
const set<int> groupD = {313572, 313574, 313575, 313603};
const set<int> groupE = {313629, 313630, 313688};
const set<int> groupF = {313695, 313833, 313878};
const set<int> groupG = {313929, 313935, 313984, 314014};
const set<int> groupH = {314077, 314105};
const set<int> groupI = {314112, 314157, 314170};

const vector <pair <string, const set<int>*>> runGroups = {
  {"GroupA", &groupA},
  {"GroupB", &groupB},
  {"GroupC", &groupC},
  {"GroupD", &groupD},
  {"GroupE", &groupE},
  {"GroupF", &groupF},
  {"GroupG", &groupG},
  {"GroupH", &groupH},
  {"GroupI", &groupI},
};
const int numRunGroups = runGroups.size ();


TString GetRunGroupTString (int rn) {
  for (const auto& group : runGroups) {
    if (group.second->find (rn) != group.second->end ())
      return TString (group.first);
  }
  return "";
}


short GetRunGroup (int rn) {
  short rg = 0;
  for (const auto& group : runGroups) {
    if (group.second->find (rn) == group.second->end ())
      rg++;
    else
      return rg;
  }
  return -1;
}


const int nPtJBins = 60;
double* pTJBins = logspace (30, 450, nPtJBins);

short GetPtJBin (const float jpt) {
  short iPtJ = 0;
  while (iPtJ < nPtJBins) {
    if (jpt < pTJBins[iPtJ+1])
      break;
    else
      iPtJ++;
  }
  return iPtJ;
}


const int nDPhiBins = 24;
short GetDPhiBin (const float dphi) {
  short iDPhi = 0;
  while (iDPhi < nDPhiBins) {
    if (dphi < (iDPhi+1)*(pi/nDPhiBins))
      break;
    else
      iDPhi++;
  }
  return iDPhi;
}

const int nPtChBins = 50;
double pTChBins[nPtChBins+1] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60};
short GetPtChBin (const float ptch) {
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
