#ifndef __CentralityDefs_h__
#define __CentralityDefs_h__

#include <TString.h>

#include <fstream>
#include <string>

namespace JetHadronCorrelations { 

extern TString workPath;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variables for defining centrality classes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FCal-derived centrality classes for categorizing events
// See Run 2 centrality discussion in Kurt's note: https://cds.cern.ch/record/2301540/
// See Run 1 centrality overview in Dennis' note: https://cds.cern.ch/record/1545591/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int fcalCentPercs[] = {100, 80, 60, 40, 20, 0};
int numFcalCentBins = sizeof (fcalCentPercs) / sizeof (fcalCentPercs[0]) - 1;

/**
 * Initializes fcal centrality cuts from latest definitions (stored in aux directory).
 */
double* InitFCalCentBins () {
  std::ifstream cutsfile;
  cutsfile.open (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()));

  double* cuts = new double[numFcalCentBins+1];
  int iCent = 0;
  std::string perc, dummy;
  double cut;
  while (cutsfile && iCent <= numFcalCentBins) {
    cutsfile >> perc >> dummy;
    cut = std::atof (dummy.c_str ());
    if (strcmp (perc.c_str (), Form ("%i%%", fcalCentPercs[iCent])) == 0)
      cuts[iCent++] = cut;
  }

  return cuts;
}

double* fcalCentBins = InitFCalCentBins ();


/**
 * Returns the bin corresponding to this sum fcal et.
 * Returns -1 for >100% peripheral collisions (i.e. FCal Et is less than energy for 100% centrality).
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FCal-derived centrality classes for mixing p+Pb events
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int fcalMixPercs[] = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0};
int numFcalMixBins = sizeof (fcalMixPercs) / sizeof (fcalMixPercs[0]) - 1;

/**
 * Initializes fcal centrality cuts from latest definitions (stored in aux directory).
 */
double* InitFCalMixBins () {
  std::ifstream cutsfile;
  cutsfile.open (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()));

  double* cuts = new double[numFcalMixBins+1];
  int iCent = 0;
  std::string perc, dummy;
  double cut;
  while (cutsfile && iCent <= numFcalMixBins) {
    cutsfile >> perc >> dummy;
    cut = std::atof (dummy.c_str ());
    if (strcmp (perc.c_str (), Form ("%i%%", fcalMixPercs[iCent])) == 0)
      cuts[iCent++] = cut;
  }

  return cuts;
}

double* fcalMixBins = InitFCalMixBins ();


/**
 * Returns the bin corresponding to this sum fcal et.
 * Returns -1 for >100% peripheral collisions (i.e. FCal Et is less than energy for 100% centrality).
 */
short GetFcalMixBin (const float fcal_et) {
  if (fcal_et < fcalMixBins[0])
    return -1;
  short i = 0;
  while (i < numFcalMixBins) {
    i++;
    if (fcal_et < fcalMixBins[i])
      break;
  }
  return i-1;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Very granular (1%) FCal-derived centrality classes for studying bias selection
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int fineFcalCentPercs[] = {100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
int numFineFcalCentBins = sizeof (fineFcalCentPercs) / sizeof (fineFcalCentPercs[0]) - 1;

/**
 * Initializes fcal centrality cuts from latest definitions (stored in aux directory).
 */
double* InitFineFCalCentBins () {
  std::ifstream cutsfile;
  cutsfile.open (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()));

  double* cuts = new double[numFineFcalCentBins+1];
  int iCent = 0;
  std::string perc, dummy;
  double cut;
  while (cutsfile && iCent <= numFineFcalCentBins) {
    cutsfile >> perc >> dummy;
    cut = std::atof (dummy.c_str ());
    if (strcmp (perc.c_str (), Form ("%i%%", fineFcalCentPercs[iCent])) == 0)
      cuts[iCent++] = cut;
  }

  return cuts;
}

double* fineFcalCentBins = InitFineFCalCentBins ();


/**
 * Returns the bin corresponding to this sum fcal et.
 * Returns -1 for >100% peripheral collisions (i.e. FCal Et is less than energy for 100% centrality).
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Zdc-derived centrality classes for categorizing events
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int zdcCentPercs[] = {100, 80, 60, 40, 20, 0};
int numZdcCentBins = sizeof (zdcCentPercs) / sizeof (zdcCentPercs[0]) - 1;

/**
 * Initializes Zdc centrality cuts from latest definitions (stored in aux directory).
 */
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
    cut = std::atof (dummy.c_str ());
    if (strcmp (perc.c_str (), Form ("%i%%", zdcCentPercs[iCent])) == 0)
      cuts[iCent++] = cut;
  }

  return cuts;
}

double* zdcCentBins = InitZdcCentBins ();


/**
 * Returns the bin corresponding to this Pb-Zdc energy.
 * Returns -1 for >100% peripheral collisions (i.e. Pb-Zdc energy is less than energy for 100% centrality).
 */
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FCal-derived centrality classes for mixing pp events
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ppMixPercs[] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0};
int numppMixBins = sizeof (ppMixPercs) / sizeof (ppMixPercs[0]) - 1;

/**
 * Initializes fcal centrality cuts from latest definitions (stored in aux directory).
 */
double* InitppMixBins () {
  std::ifstream cutsfile;
  cutsfile.open (Form ("%s/aux/ppMixCuts.dat", workPath.Data ()));

  double* cuts = new double[numppMixBins+1];
  int iCent = 0;
  std::string perc, dummy;
  double cut;
  while (cutsfile && iCent <= numppMixBins) {
    cutsfile >> perc >> dummy;
    cut = std::atof (dummy.c_str ());
    if (strcmp (perc.c_str (), Form ("%i%%", ppMixPercs[iCent])) == 0)
      cuts[iCent++] = cut;
  }

  return cuts;
}

double* ppMixBins = InitFCalMixBins ();


/**
 * Returns the bin corresponding to this sum fcal et.
 * Returns -1 for >100% peripheral collisions (i.e. FCal Et is less than energy for 100% centrality).
 */
short GetppMixBin (const float fcal_et) {
  if (fcal_et < ppMixBins[0])
    return -1;
  short i = 0;
  while (i < numppMixBins) {
    i++;
    if (fcal_et < ppMixBins[i])
      break;
  }
  return i-1;
}


} // end namespace

#endif
