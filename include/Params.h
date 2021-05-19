#ifndef __Params_h__
#define __Params_h__

#include "LocalUtilities.h"

#include <Utilities.h>

#include <TString.h>

#include <fstream>
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


// Centrality classes for mixing events
// See Run 2 centrality discussion in Kurt's note: https://cds.cern.ch/record/2301540/
// See Run 1 centrality overview in Dennis' note: https://cds.cern.ch/record/1545591/
int numFcalCentBins = 5;

double* InitFCalCentBins () {
  std::ifstream cutsfile;
  cutsfile.open (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()));

  double* cuts = new double[numFcalCentBins+1];
  int iCent = 0;
  std::string perc, dummy;
  double cut;
  while (cutsfile && iCent <= numFcalCentBins) {
    cutsfile >> perc >> dummy;
    cut = atof (dummy.c_str ());
    if (strcmp (perc.c_str (), "100%") == 0 || 
        strcmp (perc.c_str (), "80%") == 0 ||
        strcmp (perc.c_str (), "60%") == 0 ||
        strcmp (perc.c_str (), "40%") == 0 ||
        strcmp (perc.c_str (), "20%") == 0 ||
        strcmp (perc.c_str (), "0%") == 0)
      cuts[iCent++] = cut;
  }

  return cuts;
}

double* fcalCentBins = InitFCalCentBins ();
//extern double* fcalCentBins;
int fcalCentPercs[] = {100, 80, 60, 40, 20, 0};


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


int numFineFcalCentBins = 100;

double* InitFineFCalCentBins () {
  std::ifstream cutsfile;
  cutsfile.open (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()));

  double* cuts = new double[numFineFcalCentBins+1];
  int iCent = 0;
  std::string perc, dummy;
  double cut;
  while (cutsfile && iCent <= numFineFcalCentBins) {
    cutsfile >> perc >> dummy;
    cut = atof (dummy.c_str ());
    //if (strcmp (perc.c_str (), "100%") == 0 || 
    //    strcmp (perc.c_str (), "80%") == 0 ||
    //    strcmp (perc.c_str (), "60%") == 0 ||
    //    strcmp (perc.c_str (), "40%") == 0 ||
    //    strcmp (perc.c_str (), "20%") == 0 ||
    //    strcmp (perc.c_str (), "19%") == 0 ||
    //    strcmp (perc.c_str (), "18%") == 0 ||
    //    strcmp (perc.c_str (), "17%") == 0 ||
    //    strcmp (perc.c_str (), "16%") == 0 ||
    //    strcmp (perc.c_str (), "15%") == 0 ||
    //    strcmp (perc.c_str (), "14%") == 0 ||
    //    strcmp (perc.c_str (), "13%") == 0 ||
    //    strcmp (perc.c_str (), "12%") == 0 ||
    //    strcmp (perc.c_str (), "11%") == 0 ||
    //    strcmp (perc.c_str (), "10%") == 0 ||
    //    strcmp (perc.c_str (), "9%") == 0 ||
    //    strcmp (perc.c_str (), "8%") == 0 ||
    //    strcmp (perc.c_str (), "7%") == 0 ||
    //    strcmp (perc.c_str (), "6%") == 0 ||
    //    strcmp (perc.c_str (), "5%") == 0 ||
    //    strcmp (perc.c_str (), "4%") == 0 ||
    //    strcmp (perc.c_str (), "3%") == 0 ||
    //    strcmp (perc.c_str (), "2%") == 0 ||
    //    strcmp (perc.c_str (), "1%") == 0 ||
    //    strcmp (perc.c_str (), "0%") == 0)
    cuts[iCent++] = cut;
  }

  return cuts;
}

double* fineFcalCentBins = InitFineFCalCentBins ();
//extern double* fineFcalCentBins;
int fineFcalCentPercs[] = {100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};


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


int numZdcCentBins = 5;

double* InitZdcCentBins () {
  std::ifstream cutsfile;
  cutsfile.open (Form ("%s/aux/ZdcCentCuts.dat", workPath.Data ()));

  double* cuts = new double[numZdcCentBins+1];
  int iCent = 0;
  std::string perc, dummy;
  double cut;
  std::getline (cutsfile, dummy); // first line is a comment
  while (cutsfile && iCent <= numZdcCentBins) {
    cutsfile >> perc >> dummy;
    cut = atof (dummy.c_str ());
    if (strcmp (perc.c_str (), "100%") == 0 || 
        strcmp (perc.c_str (), "80%") == 0 ||
        strcmp (perc.c_str (), "60%") == 0 ||
        strcmp (perc.c_str (), "40%") == 0 ||
        strcmp (perc.c_str (), "20%") == 0 ||
        strcmp (perc.c_str (), "0%") == 0)
      cuts[iCent++] = cut;
  }

  return cuts;
}

double* zdcCentBins = InitZdcCentBins ();
//extern double* zdcCentBins;
int zdcCentPercs[] = {100, 80, 60, 40, 20, 0};


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
