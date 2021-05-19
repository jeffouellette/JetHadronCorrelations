#ifndef __JetHadronCorrelatorCombineFineFCalPtCh_C__
#define __JetHadronCorrelatorCombineFineFCalPtCh_C__

#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <math.h>

using namespace JetHadronCorrelations;

void CombineFineFCalPtCh (const char* fileTag) {

  TFile* file = nullptr;

  TH1D** h_jet_trk_pt_ns_iaa      = Get1DArray <TH1D*> (numFineFcalCentBins);
  TH1D** h_jet_trk_pt_as_iaa      = Get1DArray <TH1D*> (numFineFcalCentBins);

  TH1D*  h_jet_trk_pt_ns_iaa_comb = nullptr;
  TH1D*  h_jet_trk_pt_as_iaa_comb = nullptr;


  TString fileName = fileTag;
  fileName.ReplaceAll (".root", "");
  fileName = Form ("%s/Results/PlotFineFCalPtCh_%s.root", rootPath.Data (), fileName.Data ());
  std::cout << "Writing " << fileName.Data () << std::endl;
  file = new TFile (fileName.Data (), "update");

  file->Delete ("h_jet_trk_pt_ns_iaa_FineFcalComb;*");
  file->Delete ("h_jet_trk_pt_as_iaa_FineFcalComb;*");


  for (int iCent = 0; iCent < numFineFcalCentBins; iCent++) {
    h_jet_trk_pt_ns_iaa[iCent] = (TH1D*) file->Get (Form ("h_jet_trk_pt_ns_iaa_FineFcalCent%i", iCent));
    h_jet_trk_pt_as_iaa[iCent] = (TH1D*) file->Get (Form ("h_jet_trk_pt_as_iaa_FineFcalCent%i", iCent));
  }


  h_jet_trk_pt_ns_iaa_comb = (TH1D*) h_jet_trk_pt_ns_iaa[0]->Clone ("h_jet_trk_pt_ns_iaa_FineFcalComb"); 
  h_jet_trk_pt_ns_iaa_comb->Reset ();

  h_jet_trk_pt_as_iaa_comb = (TH1D*) h_jet_trk_pt_as_iaa[0]->Clone ("h_jet_trk_pt_as_iaa_FineFcalComb"); 
  h_jet_trk_pt_as_iaa_comb->Reset ();

  TH1D* h_wgts = GetFCalZdcWeights ();

  TH1D* h_wgts_rebinned = new TH1D ("h_wgts_rebinned", ";#Sigma#it{E}_{T}^{FCal} percentile [%];Weight", 100, 0, 100);

  double sum_wgts = 0.;
  for (int iCent = 0; iCent < numFineFcalCentBins; iCent++) {

    const double plo = fineFcalCentBins[iCent];
    const double phi = fineFcalCentBins[iCent+1];

    double wgt = 0.;
    double n = 0.;
    for (int iX = 1; iX <= h_wgts->GetNbinsX (); iX++) {
      double ci = h_wgts->GetBinLowEdge (iX);
      double cip1 = h_wgts->GetBinLowEdge (iX) + h_wgts->GetBinWidth (iX);

      if (!(plo >= cip1 || phi <= ci)) {
        const double width = std::fmin (phi, cip1) - std::fmax (plo, ci);
        wgt += h_wgts->GetBinContent (iX) * width;
        n += width;
      }
    }
    assert (n > 0.);
    wgt = wgt / n;

    assert (!std::isnan (wgt));

    h_wgts_rebinned->SetBinContent (h_wgts_rebinned->FindBin (0.5 * (fineFcalCentPercs[iCent] + fineFcalCentPercs[iCent+1])), wgt);

    h_jet_trk_pt_ns_iaa_comb->Add (h_jet_trk_pt_ns_iaa[iCent], wgt);
    h_jet_trk_pt_as_iaa_comb->Add (h_jet_trk_pt_as_iaa[iCent], wgt);

    //for (int iX = 1; iX <= h_jet_trk_pt_ns_iaa_comb->GetNbinsX (); iX++) {
    //  if (std::isnan (h_jet_trk_pt_ns_iaa[iCent]->GetBinContent (iX)) || std::isnan (h_jet_trk_pt_ns_iaa[iCent]->GetBinError (iX))) {
    //    std::cout << "NaN in iaa? iCent = " << iCent << ", pT = " << h_jet_trk_pt_ns_iaa[iCent]->GetBinCenter (iX) << std::endl;
    //    continue; // indicates no yield in pp
    //  }
    //  h_jet_trk_pt_ns_iaa_comb->SetBinContent (iX, h_jet_trk_pt_ns_iaa_comb->GetBinContent (iX) * wgt);
    //  h_jet_trk_pt_ns_iaa_comb->SetBinError (iX, std::sqrt (std::pow (h_jet_trk_pt_ns_iaa_comb->GetBinError (iX), 2) + std::pow (h_jet_trk_pt_ns_iaa[iCent]->GetBinError (iX) * wgt, 2)));
    //}

    //for (int iX = 1; iX <= h_jet_trk_pt_as_iaa_comb->GetNbinsX (); iX++) {
    //  if (std::isnan (h_jet_trk_pt_as_iaa[iCent]->GetBinContent (iX)) || std::isnan (h_jet_trk_pt_as_iaa[iCent]->GetBinError (iX))) {
    //    std::cout << "NaN in iaa? iCent = " << iCent << ", pT = " << h_jet_trk_pt_as_iaa[iCent]->GetBinCenter (iX) << std::endl;
    //    continue; // indicates no yield in pp
    //  }
    //  h_jet_trk_pt_as_iaa_comb->SetBinContent (iX, h_jet_trk_pt_as_iaa_comb->GetBinContent (iX) * wgt);
    //  h_jet_trk_pt_as_iaa_comb->SetBinError (iX, std::sqrt (std::pow (h_jet_trk_pt_as_iaa_comb->GetBinError (iX), 2) + std::pow (h_jet_trk_pt_as_iaa[iCent]->GetBinError (iX) * wgt, 2)));
    //}

    sum_wgts += wgt;
  }

  h_jet_trk_pt_ns_iaa_comb->Scale (1./sum_wgts);
  h_jet_trk_pt_as_iaa_comb->Scale (1./sum_wgts);


  file->cd ();

  h_wgts_rebinned->Write ();

  h_jet_trk_pt_ns_iaa_comb->Write ();
  h_jet_trk_pt_as_iaa_comb->Write ();

  file->Close ();
  
}


#endif
