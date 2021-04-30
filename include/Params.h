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


const float min_trk_pt = 0.4;
const float min_akt4_hi_jet_pt = 30;

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
extern double jet_min_pt;
extern double jet_max_pt;


// Centrality classes for mixing events
// See centralities defined in Kurt's note: https://cds.cern.ch/record/2301540/files/ATL-COM-DAPR-2018-002.pdf
double fcalCentBins[] = {
 -1000,     // 100%
   //  2.21,  // 90%
   //  4.94,  // 80%
   //  8.24,  // 70%
   // 12.14,  // 60%
   // 16.71,  // 50%
   // 22.11,  // 40%
   // 28.64,  // 30%
    37.07,  // 20%
   // 39.13,  // 18%
   // 41.40,  // 16%
   // 43.90,  // 14%
   // 46.71,  // 12%
   // 49.94,  // 10%
   // 53.78,  // 8%
   // 58.54,  // 6%
   // 64.96,  // 4%
   // 75.34,  // 2%
   // 85.08,  // 1%
   220     // "0%" -- really just an upper bound for sanity
};
int numFcalCentBins = sizeof (fcalCentBins) / sizeof (fcalCentBins[0]) - 1;


/**
 * Returns the bin corresponding to this sum fcal et bin.
 * Returns -1 for >90% peripheral collisions (i.e. FCal Et is less than energy for 90% centrality).
 */
short GetFcalCentBin (const float fcal_et) {
  if (fcal_et < fcalCentBins[0])
    return -1;
  short i = 0;
  while (i < numFcalCentBins) {
    i++;
    if (fcal_et < fcalCentBins[i])
      break;
  }
  return i-1;
}

double fineFcalCentBins[] = {
 -1000,     // 100%
   //  2.21,  // 90%
   //  4.94,  // 80%
   //  8.24,  // 70%
   // 12.14,  // 60%
   // 16.71,  // 50%
   // 22.11,  // 40%
   // 28.64,  // 30%
    37.07,  // 20%
    39.13,  // 18%
    41.40,  // 16%
    43.90,  // 14%
    46.71,  // 12%
    49.94,  // 10%
    53.78,  // 8%
    58.54,  // 6%
    64.96,  // 4%
    75.34,  // 2%
   // 85.08,  // 1%
   220     // "0%" -- really just an upper bound for sanity
};
int numFineFcalCentBins = sizeof (fineFcalCentBins) / sizeof (fineFcalCentBins[0]) - 1;


/**
 * Returns the bin corresponding to this sum fcal et bin.
 * Returns -1 for >90% peripheral collisions (i.e. FCal Et is less than energy for 90% centrality).
 */
short GetFineFcalCentBin (const float fcal_et) {
  if (fcal_et < fineFcalCentBins[0])
    return -1;
  short i = 0;
  while (i < numFineFcalCentBins) {
    i++;
    if (fcal_et < fineFcalCentBins[i])
      break;
  }
  return i-1;
}


double zdcCentBins[] = {
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
int numZdcCentBins = sizeof (zdcCentBins) / sizeof (zdcCentBins[0]) - 1; // no bin needed for pp


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


double* centBins = nullptr;
int numCentBins = 0;

short GetCentBin (const float val) {
  if (numCentBins == 0 || val < centBins[0])
    return -1;
  short i = 0;
  while (i < numCentBins) {
    i++;
    if (val < centBins[i])
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

const double pTChBins[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60};
const int nPtChBins = sizeof (pTChBins) / sizeof (pTChBins[0]) - 1;
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
