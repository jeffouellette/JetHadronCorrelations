#ifndef __CentralityDefs_h__
#define __CentralityDefs_h__

#include <TString.h>

#include <fstream>
#include <string>
#include <cassert>

namespace JetHadronCorrelations { 

extern TString workPath;


/**
 * Finds the bin corresponding to this value in the specified array.
 * Note: array must be in ascending order.
 * Returns -1 for values lower than the 0th element in the array.
 */
short GetBin (double* bins, const int nBins, const float val) {
  if (val < bins[0])
    return -1;
  short i = 0;
  while (i++ < nBins) {
    if (val < bins[i])
      break;
  }
  return i-1;
}


/**
 * Initializes centrality cuts from latest definitions (stored in aux directory).
 */
double* InitCentBins (const char* fName, const int* percs, const int nBins) {
  assert (percs != nullptr);
  assert (nBins > 0);

  std::ifstream cutsfile;
  cutsfile.open (fName);

  double* cuts = new double[nBins+1];
  int iCent = 0;
  std::string perc, dummy;
  double cut;

  while (cutsfile && iCent <= nBins) {
    cutsfile >> perc >> dummy;
    cut = std::atof (dummy.c_str ());
    if (strcmp (perc.c_str (), Form ("%i%%", percs[iCent])) == 0)
      cuts[iCent++] = cut;
  }

  //for (int i = 0; i < nBins; i++)
  //  std::cout << cuts[i] << ", ";
  //std::cout << cuts[nBins] << std::endl;

  return cuts;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variables for defining centrality classes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FCal-derived centrality classes for categorizing events
// See Run 2 centrality discussion in Kurt's note: https://cds.cern.ch/record/2301540/
// See Run 1 centrality overview in Dennis' note: https://cds.cern.ch/record/1545591/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int fcalCentPercs[] = {100, 80, 60, 40, 20, 0};
int nFcalCentBins = sizeof (fcalCentPercs) / sizeof (fcalCentPercs[0]) - 1;

double* fcalCentBins = InitCentBins (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()), fcalCentPercs, nFcalCentBins);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FCal-derived centrality classes for mixing p+Pb events
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int fcalMixPercs[] = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0};
int nFcalMixBins = sizeof (fcalMixPercs) / sizeof (fcalMixPercs[0]) - 1;

double* fcalMixBins = InitCentBins (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()), fcalMixPercs, nFcalMixBins);


int fcalMixVar2Percs[] = {100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
int nFcalMixVar2Bins = sizeof (fcalMixVar2Percs) / sizeof (fcalMixVar2Percs[0]) - 1;

double* fcalMixVar2Bins = InitCentBins (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()), fcalMixVar2Percs, nFcalMixVar2Bins);


int fcalMixVar6Percs[] = {100, 80, 60, 40, 20, 0};
int nFcalMixVar6Bins = sizeof (fcalMixVar6Percs) / sizeof (fcalMixVar6Percs[0]) - 1;

double* fcalMixVar6Bins = InitCentBins (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()), fcalMixVar6Percs, nFcalMixVar6Bins);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Very granular (1%) FCal-derived centrality classes for studying bias selection
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int fineFcalCentPercs[] = {100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
int nFineFcalCentBins = sizeof (fineFcalCentPercs) / sizeof (fineFcalCentPercs[0]) - 1;

double* fineFcalCentBins = InitCentBins (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()), fineFcalCentPercs, nFineFcalCentBins);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Zdc-derived centrality classes for categorizing events
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const double zdcNcollValues[] = {1.24, 3.04, 6.53, 9.99, 13.6};
const double zdcNcollErrors[] = {0.10, 0.18, 0.50, 0.90, 1.48};
int zdcCentPercs[] = {100, 80, 60, 40, 20, 0};
int nZdcCentBins = sizeof (zdcCentPercs) / sizeof (zdcCentPercs[0]) - 1;

double* zdcCentBins = InitCentBins (Form ("%s/aux/ZdcCentCuts.dat", workPath.Data ()), zdcCentPercs, nZdcCentBins);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FCal-derived centrality classes for mixing pp events
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ppMixPercs[] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0};
int nppMixBins = sizeof (ppMixPercs) / sizeof (ppMixPercs[0]) - 1;

double* ppMixBins = InitCentBins (Form ("%s/aux/ppMixCuts.dat", workPath.Data ()), ppMixPercs, nppMixBins);


int ppMixVar1Percs[] = {100, 98, 96, 94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 74, 72, 70, 68, 66, 64, 62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0};
int nppMixVar1Bins = sizeof (ppMixVar1Percs) / sizeof (ppMixVar1Percs[0]) - 1;

double* ppMixVar1Bins = InitCentBins (Form ("%s/aux/ppMixCuts.dat", workPath.Data ()), ppMixVar1Percs, nppMixVar1Bins);


int ppMixVar3Percs[] = {100, 95, 85, 75, 65, 55, 45, 35, 25, 15, 5, 0};
int nppMixVar3Bins = sizeof (ppMixVar3Percs) / sizeof (ppMixVar3Percs[0]) - 1;

double* ppMixVar3Bins = InitCentBins (Form ("%s/aux/ppMixCuts.dat", workPath.Data ()), ppMixVar3Percs, nppMixVar3Bins);



} // end namespace

#endif
