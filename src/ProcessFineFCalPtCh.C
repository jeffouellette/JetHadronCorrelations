#ifndef __JetHadronCorrelatorProcessFineFCalPtCh_C__
#define __JetHadronCorrelatorProcessFineFCalPtCh_C__

#include "CentralityDefs.h"
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

void ProcessFineFCalPtCh (const char* tag, const char* outFileTag) {

  TFile* inFile = nullptr;

  TH1D*  h_evt_counts_ref         = nullptr;
  TH1D*  h_jet_counts_ref         = nullptr;
  TH1D*  h_evt_counts_ref_bkg     = nullptr;
  TH1D*  h_jet_counts_ref_bkg     = nullptr;
  TH1D** h_evt_counts             = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH1D** h_jet_counts             = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH1D** h_evt_counts_bkg         = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH1D** h_jet_counts_bkg         = Get1DArray <TH1D*> (nFineFcalCentBins);

  TH1D*  h_jet_trk_pt_ns_ref            = nullptr;
  TH2D*  h2_jet_trk_pt_ns_cov_ref       = nullptr;
  TH1D*  h_jet_trk_pt_perp_ref          = nullptr;
  TH2D*  h2_jet_trk_pt_perp_cov_ref     = nullptr;
  TH1D*  h_jet_trk_pt_as_ref            = nullptr;
  TH2D*  h2_jet_trk_pt_as_cov_ref       = nullptr;
  TH1D*  h_jet_trk_pt_ns_ref_bkg        = nullptr;
  TH2D*  h2_jet_trk_pt_ns_cov_ref_bkg   = nullptr;
  TH1D*  h_jet_trk_pt_perp_ref_bkg      = nullptr;
  TH2D*  h2_jet_trk_pt_perp_cov_ref_bkg = nullptr;
  TH1D*  h_jet_trk_pt_as_ref_bkg        = nullptr;
  TH2D*  h2_jet_trk_pt_as_cov_ref_bkg   = nullptr;
  TH1D** h_jet_trk_pt_ns                = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH2D** h2_jet_trk_pt_ns_cov           = Get1DArray <TH2D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_perp              = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH2D** h2_jet_trk_pt_perp_cov         = Get1DArray <TH2D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_as                = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH2D** h2_jet_trk_pt_as_cov           = Get1DArray <TH2D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_ns_bkg            = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH2D** h2_jet_trk_pt_ns_cov_bkg       = Get1DArray <TH2D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_perp_bkg          = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH2D** h2_jet_trk_pt_perp_cov_bkg     = Get1DArray <TH2D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_as_bkg            = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH2D** h2_jet_trk_pt_as_cov_bkg       = Get1DArray <TH2D*> (nFineFcalCentBins);

  TH1D*  h_jet_trk_pt_ns_ref_sig    = nullptr;
  TH1D*  h_jet_trk_pt_perp_ref_sig  = nullptr;
  TH1D*  h_jet_trk_pt_as_ref_sig    = nullptr;
  TH1D** h_jet_trk_pt_ns_sig        = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_perp_sig      = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_as_sig        = Get1DArray <TH1D*> (nFineFcalCentBins);

  TH1D** h_jet_trk_pt_ns_iaa        = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_perp_iaa      = Get1DArray <TH1D*> (nFineFcalCentBins);
  TH1D** h_jet_trk_pt_as_iaa        = Get1DArray <TH1D*> (nFineFcalCentBins);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/PlotFineFCalPtCh_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  {
    TString inFileName = Form ("%s/Histograms/%s/JetsHists/FineFcalCentVar/data17_5TeV_hists.root", rootPath.Data (), tag);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");
    outFile->cd ();

    h_evt_counts_ref = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", tag))->Clone ("h_evt_counts_ref");
    h_jet_counts_ref = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", tag))->Clone ("h_jet_counts_ref");

    h_jet_trk_pt_ns_ref = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data17", tag))->Clone ("h_jet_trk_pt_ns_ref");
    h2_jet_trk_pt_ns_cov_ref = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data17", tag))->Clone ("h2_jet_trk_pt_ns_cov_ref");
    h_jet_trk_pt_perp_ref = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_data17", tag))->Clone ("h_jet_trk_pt_perp_ref");
    h2_jet_trk_pt_perp_cov_ref = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_data17", tag))->Clone ("h2_jet_trk_pt_perp_cov_ref");
    h_jet_trk_pt_as_ref = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data17", tag))->Clone ("h_jet_trk_pt_as_ref");
    h2_jet_trk_pt_as_cov_ref = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data17", tag))->Clone ("h2_jet_trk_pt_as_cov_ref");

    inFile->Close ();

    h_jet_trk_pt_ns_ref->Rebin (5);
    h2_jet_trk_pt_ns_cov_ref->RebinX (5);
    h2_jet_trk_pt_ns_cov_ref->RebinY (5);
    h_jet_trk_pt_perp_ref->Rebin (5);
    h2_jet_trk_pt_perp_cov_ref->RebinX (5);
    h2_jet_trk_pt_perp_cov_ref->RebinY (5);
    h_jet_trk_pt_as_ref->Rebin (5);
    h2_jet_trk_pt_as_cov_ref->RebinX (5);
    h2_jet_trk_pt_as_cov_ref->RebinY (5);

    CalcUncertainties (h_jet_trk_pt_ns_ref, h2_jet_trk_pt_ns_cov_ref, h_jet_counts_ref);
    CalcUncertainties (h_jet_trk_pt_perp_ref, h2_jet_trk_pt_perp_cov_ref, h_jet_counts_ref);
    CalcUncertainties (h_jet_trk_pt_as_ref, h2_jet_trk_pt_as_cov_ref, h_jet_counts_ref);
  }


  {
    TString inFileName = Form ("%s/Histograms/%s/MixedHists/FineFcalCentVar/data17_5TeV_hists.root", rootPath.Data (), tag);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");
    outFile->cd ();

    h_evt_counts_ref_bkg = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_mixed_data17", tag))->Clone ("h_evt_counts_ref_bkg");
    h_jet_counts_ref_bkg = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_data17", tag))->Clone ("h_jet_counts_ref_bkg");

    h_jet_trk_pt_ns_ref_bkg = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_mixed_data17", tag))->Clone ("h_jet_trk_pt_ns_ref_bkg");
    h2_jet_trk_pt_ns_cov_ref_bkg = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_mixed_data17", tag))->Clone ("h2_jet_trk_pt_ns_cov_ref_bkg");
    h_jet_trk_pt_perp_ref_bkg = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_mixed_data17", tag))->Clone ("h_jet_trk_pt_perp_ref_bkg");
    h2_jet_trk_pt_perp_cov_ref_bkg = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_mixed_data17", tag))->Clone ("h2_jet_trk_pt_perp_cov_ref_bkg");
    h_jet_trk_pt_as_ref_bkg = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_mixed_data17", tag))->Clone ("h_jet_trk_pt_as_ref_bkg");
    h2_jet_trk_pt_as_cov_ref_bkg = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_mixed_data17", tag))->Clone ("h2_jet_trk_pt_as_cov_ref_bkg");

    inFile->Close ();

    h_jet_trk_pt_ns_ref_bkg->Rebin (5);
    h2_jet_trk_pt_ns_cov_ref_bkg->RebinX (5);
    h2_jet_trk_pt_ns_cov_ref_bkg->RebinY (5);
    h_jet_trk_pt_perp_ref_bkg->Rebin (5);
    h2_jet_trk_pt_perp_cov_ref_bkg->RebinX (5);
    h2_jet_trk_pt_perp_cov_ref_bkg->RebinY (5);
    h_jet_trk_pt_as_ref_bkg->Rebin (5);
    h2_jet_trk_pt_as_cov_ref_bkg->RebinX (5);
    h2_jet_trk_pt_as_cov_ref_bkg->RebinY (5);

    CalcUncertainties (h_jet_trk_pt_ns_ref_bkg, h2_jet_trk_pt_ns_cov_ref_bkg, h_jet_counts_ref_bkg);
    CalcUncertainties (h_jet_trk_pt_perp_ref_bkg, h2_jet_trk_pt_perp_cov_ref_bkg, h_jet_counts_ref_bkg);
    CalcUncertainties (h_jet_trk_pt_as_ref_bkg, h2_jet_trk_pt_as_cov_ref_bkg, h_jet_counts_ref_bkg);
  }


  for (int iCent = 0; iCent < nFineFcalCentBins; iCent++) {
    TString inFileName = Form ("%s/Histograms/%s/JetsHists/FineFcalCentVar/data16_5TeV_iCent%i_hists.root", rootPath.Data (), tag, iCent);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");
    outFile->cd ();

    h_evt_counts[iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_pPb_FineFcalCent%i", iCent));
    h_jet_counts[iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", tag))->Clone (Form ("h_jet_counts_pPb_FineFcalCent%i", iCent));

    h_jet_trk_pt_ns[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data16", tag))->Clone (Form ("h_jet_trk_pt_ns_pPb_FineFcalCent%i", iCent));
    h2_jet_trk_pt_ns_cov[iCent] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data16", tag))->Clone (Form ("h2_jet_trk_pt_ns_cov_pPb_FineFcalCent%i", iCent));
    h_jet_trk_pt_perp[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_data16", tag))->Clone (Form ("h_jet_trk_pt_perp_pPb_FineFcalCent%i", iCent));
    h2_jet_trk_pt_perp_cov[iCent] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_data16", tag))->Clone (Form ("h2_jet_trk_pt_perp_cov_pPb_FineFcalCent%i", iCent));
    h_jet_trk_pt_as[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data16", tag))->Clone (Form ("h_jet_trk_pt_as_pPb_FineFcalCent%i", iCent));
    h2_jet_trk_pt_as_cov[iCent] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data16", tag))->Clone (Form ("h2_jet_trk_pt_as_cov_pPb_FineFcalCent%i", iCent));

    inFile->Close ();

    h_jet_trk_pt_ns[iCent]->Rebin (5);
    h2_jet_trk_pt_ns_cov[iCent]->RebinX (5);
    h2_jet_trk_pt_ns_cov[iCent]->RebinY (5);
    h_jet_trk_pt_perp[iCent]->Rebin (5);
    h2_jet_trk_pt_perp_cov[iCent]->RebinX (5);
    h2_jet_trk_pt_perp_cov[iCent]->RebinY (5);
    h_jet_trk_pt_as[iCent]->Rebin (5);
    h2_jet_trk_pt_as_cov[iCent]->RebinX (5);
    h2_jet_trk_pt_as_cov[iCent]->RebinY (5);

    CalcUncertainties (h_jet_trk_pt_ns[iCent], h2_jet_trk_pt_ns_cov[iCent], h_jet_counts[iCent]);
    CalcUncertainties (h_jet_trk_pt_perp[iCent], h2_jet_trk_pt_perp_cov[iCent], h_jet_counts[iCent]);
    CalcUncertainties (h_jet_trk_pt_as[iCent], h2_jet_trk_pt_as_cov[iCent], h_jet_counts[iCent]);
  }


  for (int iCent = 0; iCent < nFineFcalCentBins; iCent++) {
    TString inFileName = Form ("%s/Histograms/%s/MixedHists/FineFcalCentVar/data16_5TeV_iCent%i_hists.root", rootPath.Data (), tag, iCent);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");
    outFile->cd ();

    h_evt_counts_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_mixed_data16", tag))->Clone (Form ("h_evt_counts_pPb_bkg_FineFcalCent%i", iCent));
    h_jet_counts_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_data16", tag))->Clone (Form ("h_jet_counts_pPb_bkg_FineFcalCent%i", iCent));

    h_jet_trk_pt_ns_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_mixed_data16", tag))->Clone (Form ("h_jet_trk_pt_ns_pPb_bkg_FineFcalCent%i", iCent));
    h2_jet_trk_pt_ns_cov_bkg[iCent] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_mixed_data16", tag))->Clone (Form ("h2_jet_trk_pt_ns_cov_pPb_bkg_FineFcalCent%i", iCent));
    h_jet_trk_pt_perp_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_mixed_data16", tag))->Clone (Form ("h_jet_trk_pt_perp_pPb_bkg_FineFcalCent%i", iCent));
    h2_jet_trk_pt_perp_cov_bkg[iCent] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_mixed_data16", tag))->Clone (Form ("h2_jet_trk_pt_perp_cov_pPb_bkg_FineFcalCent%i", iCent));
    h_jet_trk_pt_as_bkg[iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_mixed_data16", tag))->Clone (Form ("h_jet_trk_pt_as_pPb_bkg_FineFcalCent%i", iCent));
    h2_jet_trk_pt_as_cov_bkg[iCent] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_mixed_data16", tag))->Clone (Form ("h2_jet_trk_pt_as_cov_pPb_bkg_FineFcalCent%i", iCent));

    inFile->Close ();

    h_jet_trk_pt_ns_bkg[iCent]->Rebin (5);
    h2_jet_trk_pt_ns_cov_bkg[iCent]->RebinX (5);
    h2_jet_trk_pt_ns_cov_bkg[iCent]->RebinY (5);
    h_jet_trk_pt_perp_bkg[iCent]->Rebin (5);
    h2_jet_trk_pt_perp_cov_bkg[iCent]->RebinX (5);
    h2_jet_trk_pt_perp_cov_bkg[iCent]->RebinY (5);
    h_jet_trk_pt_as_bkg[iCent]->Rebin (5);
    h2_jet_trk_pt_as_cov_bkg[iCent]->RebinX (5);
    h2_jet_trk_pt_as_cov_bkg[iCent]->RebinY (5);

    CalcUncertainties (h_jet_trk_pt_ns_bkg[iCent], h2_jet_trk_pt_ns_cov_bkg[iCent], h_jet_counts_bkg[iCent]);
    CalcUncertainties (h_jet_trk_pt_perp_bkg[iCent], h2_jet_trk_pt_perp_cov_bkg[iCent], h_jet_counts_bkg[iCent]);
    CalcUncertainties (h_jet_trk_pt_as_bkg[iCent], h2_jet_trk_pt_as_cov_bkg[iCent], h_jet_counts_bkg[iCent]);
  }


  {
    h_jet_trk_pt_ns_ref_sig = (TH1D*) h_jet_trk_pt_ns_ref->Clone (Form ("h_jet_trk_pt_ns_ref_sig"));
    h_jet_trk_pt_ns_ref_sig->Add (h_jet_trk_pt_ns_ref_bkg, -1);

    h_jet_trk_pt_perp_ref_sig = (TH1D*) h_jet_trk_pt_perp_ref->Clone (Form ("h_jet_trk_pt_perp_ref_sig"));
    h_jet_trk_pt_perp_ref_sig->Add (h_jet_trk_pt_perp_ref_bkg, -1);

    h_jet_trk_pt_as_ref_sig = (TH1D*) h_jet_trk_pt_as_ref->Clone (Form ("h_jet_trk_pt_as_ref_sig"));
    h_jet_trk_pt_as_ref_sig->Add (h_jet_trk_pt_as_ref_bkg, -1);

    for (int iCent = 0; iCent < nFineFcalCentBins; iCent++) {
      h_jet_trk_pt_ns_sig[iCent] = (TH1D*) h_jet_trk_pt_ns[iCent]->Clone (Form ("h_jet_trk_pt_ns_pPb_sig_FineFcalCent%i", iCent));
      h_jet_trk_pt_ns_sig[iCent]->Add (h_jet_trk_pt_ns_bkg[iCent], -1);

      h_jet_trk_pt_perp_sig[iCent] = (TH1D*) h_jet_trk_pt_perp[iCent]->Clone (Form ("h_jet_trk_pt_perp_pPb_sig_FineFcalCent%i", iCent));
      h_jet_trk_pt_perp_sig[iCent]->Add (h_jet_trk_pt_perp_bkg[iCent], -1);

      h_jet_trk_pt_as_sig[iCent] = (TH1D*) h_jet_trk_pt_as[iCent]->Clone (Form ("h_jet_trk_pt_as_pPb_sig_FineFcalCent%i", iCent));
      h_jet_trk_pt_as_sig[iCent]->Add (h_jet_trk_pt_as_bkg[iCent], -1);

      h_jet_trk_pt_ns_iaa[iCent] = (TH1D*) h_jet_trk_pt_ns_sig[iCent]->Clone (Form ("h_jet_trk_pt_ns_iaa_FineFcalCent%i", iCent));
      h_jet_trk_pt_ns_iaa[iCent]->Divide (h_jet_trk_pt_ns_ref_sig);

      h_jet_trk_pt_perp_iaa[iCent] = (TH1D*) h_jet_trk_pt_perp_sig[iCent]->Clone (Form ("h_jet_trk_pt_perp_iaa_FineFcalCent%i", iCent));
      h_jet_trk_pt_perp_iaa[iCent]->Divide (h_jet_trk_pt_perp_ref_sig);

      h_jet_trk_pt_as_iaa[iCent] = (TH1D*) h_jet_trk_pt_as_sig[iCent]->Clone (Form ("h_jet_trk_pt_as_iaa_FineFcalCent%i", iCent));
      h_jet_trk_pt_as_iaa[iCent]->Divide (h_jet_trk_pt_as_ref_sig);
    }
  }


  {
    outFile->cd ();

    h_evt_counts_ref->Write ();
    h_jet_counts_ref->Write ();

    h_evt_counts_ref_bkg->Write ();
    h_jet_counts_ref_bkg->Write ();

    h_jet_trk_pt_ns_ref->Write ();
    h_jet_trk_pt_perp_ref->Write ();
    h_jet_trk_pt_as_ref->Write ();

    h_jet_trk_pt_ns_ref_bkg->Write ();
    h_jet_trk_pt_perp_ref_bkg->Write ();
    h_jet_trk_pt_as_ref_bkg->Write ();

    for (int iCent = 0; iCent < nFineFcalCentBins; iCent++) {

      h_evt_counts[iCent]->Write ();
      h_jet_counts[iCent]->Write ();

      h_evt_counts_bkg[iCent]->Write ();
      h_jet_counts_bkg[iCent]->Write ();

      h_jet_trk_pt_ns[iCent]->Write ();
      h_jet_trk_pt_perp[iCent]->Write ();
      h_jet_trk_pt_as[iCent]->Write ();

      h_jet_trk_pt_ns_bkg[iCent]->Write ();
      h_jet_trk_pt_perp_bkg[iCent]->Write ();
      h_jet_trk_pt_as_bkg[iCent]->Write ();

      h_jet_trk_pt_ns_sig[iCent]->Write ();
      h_jet_trk_pt_perp_sig[iCent]->Write ();
      h_jet_trk_pt_as_sig[iCent]->Write ();

      h_jet_trk_pt_ns_iaa[iCent]->Write ();
      h_jet_trk_pt_perp_iaa[iCent]->Write ();
      h_jet_trk_pt_as_iaa[iCent]->Write ();
    }

    outFile->Close ();
  }
}


#endif
