#ifndef __MakeLatexCentralityTable_C__
#define __MakeLatexCentralityTable_C__

#include <iostream>
#include <fstream>

#include "include/Params.h"
#include "include/LocalUtilities.h"
#include "include/CentralityDefs.h"

using namespace JetHadronCorrelations;

bool ArrContains (const int* arr, const int len, const double val) {
  int i = 0;
  while (i < len) {
    if (arr[i] == val) return true;
    i++;
  }
  return false;
}

void MakeLatexCentralityTable () {

  int percs[] = {100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  int nPercs = sizeof (percs) / sizeof (percs[0]) - 1;

  double* pPbBins = InitCentBins (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()), percs, nPercs);
  double* ppBins = InitCentBins (Form ("%s/aux/ppMixCuts.dat", workPath.Data ()), percs, nPercs);

  ofstream tfile ("LatexCentralityTable.tex");

  tfile << "\\documentclass[IntNote.tex]{subfiles}" << std::endl;
  tfile << std::endl;
  tfile << "\\begin{document}" << std::endl;
  tfile << std::endl;
  tfile << "\\begin{center}" << std::endl;
  tfile << "\\begin{table}[]" << std::endl;
  tfile << "  \\small" << std::endl;
  tfile << "  \\centering" << std::endl;
  tfile << "  \\begin{tabular}{l|c|c||l|c|c} \\hline" << std::endl;
  tfile << "      Percentile & \\pPb, $\\sumet$ [\\GeV] & \\pp, $\\sumet$ [\\GeV] &  & \\pPb, $\\sumet$ [\\GeV] & \\pp, $\\sumet$ [\\GeV]\\\\ \\hline" << std::endl;

  for (int iPerc = 0; iPerc < nPercs/2 + 1; iPerc++) {
    tfile << percs[iPerc] << "\\\% & ";
    if (ArrContains (fcalMixPercs, nFcalMixBins+1, percs[iPerc]))
      tfile << "\\textbf{" << pPbBins[iPerc] << "} & ";
    else tfile << pPbBins[iPerc] << " & ";
    if (ArrContains (ppMixPercs, nppMixBins+1, percs[iPerc]))
      tfile << "\\textbf{" << ppBins[iPerc] << "} & ";
    else tfile << ppBins[iPerc] << " & ";

    if (iPerc == 0)
      tfile << "& & \\\\" << std::endl;
    else {
      tfile << percs[iPerc + nPercs/2] << "\\\% & ";
      if (ArrContains (fcalMixPercs, nFcalMixBins+1, percs[iPerc + nPercs/2]))
        tfile << "\\textbf{" << pPbBins[iPerc + nPercs/2] << "} & ";
      else tfile << pPbBins[iPerc + nPercs/2] << " & ";

      if (ArrContains (ppMixPercs, nppMixBins+1, percs[iPerc + nPercs/2]))
        tfile << "\\textbf{" << ppBins[iPerc + nPercs/2] << "} \\\\" << std::endl;
      else tfile << ppBins[iPerc + nPercs/2] << " \\\\" << std::endl;
    }
  }

  tfile << "    \\end{tabular}" << std::endl;
  tfile << "  \\caption{Forward calorimeter (FCal) $\\sumet$ centiles used in this analysis, for \\pPb (A-side only) and \\pp (A+C-side sum). Bolded values correspond to the bin edges used in the nominal analysis results; other values are used in variations on the matching.}" << std::endl;
  tfile << "  \\label{tab:fcaletbins}" << std::endl;
  tfile << "\\end{table}" << std::endl;
  tfile << "\\end{center}" << std::endl;
  tfile << std::endl;
  tfile << "\\end{document}" << std::endl;

  tfile.close ();
  

}

#endif
