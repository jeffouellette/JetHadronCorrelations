#ifndef __JetHadronCorrelatorProcessPtCh_C__
#define __JetHadronCorrelatorProcessPtCh_C__

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
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;

void ProcessPtCh (const char* outFileTag, const char* tag1, const char* tag2 = nullptr) {

  bool doMix = (tag2 == nullptr);

  TFile* inFile = nullptr;

  TH1D**  h_evt_counts_ref = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_counts_ref = Get1DArray <TH1D*> (nVar);
  TH1D**  h_evt_counts_ref_bkg = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_counts_ref_bkg = Get1DArray <TH1D*> (nVar);
  TH1D*** h_evt_counts = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_counts = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_evt_counts_bkg = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_counts_bkg = Get2DArray <TH1D*> (nZdcCentBins, nVar);

  TH1D**  h_jet_trk_pt_ns_ref = Get1DArray <TH1D*> (nVar);
  TH2D**  h2_jet_trk_pt_ns_cov_ref = Get1DArray <TH2D*> (nVar);
  TH1D**  h_jet_trk_pt_perp_ref = Get1DArray <TH1D*> (nVar);
  TH2D**  h2_jet_trk_pt_perp_cov_ref = Get1DArray <TH2D*> (nVar);
  TH1D**  h_jet_trk_pt_as_ref = Get1DArray <TH1D*> (nVar);
  TH2D**  h2_jet_trk_pt_as_cov_ref = Get1DArray <TH2D*> (nVar);

  TH1D**  h_jet_trk_pt_ns_ref_bkg = Get1DArray <TH1D*> (nVar);
  TH2D**  h2_jet_trk_pt_ns_cov_ref_bkg = Get1DArray <TH2D*> (nVar);
  TH1D**  h_jet_trk_pt_perp_ref_bkg = Get1DArray <TH1D*> (nVar);
  TH2D**  h2_jet_trk_pt_perp_cov_ref_bkg = Get1DArray <TH2D*> (nVar);
  TH1D**  h_jet_trk_pt_as_ref_bkg = Get1DArray <TH1D*> (nVar);
  TH2D**  h2_jet_trk_pt_as_cov_ref_bkg = Get1DArray <TH2D*> (nVar);

  TH1D*** h_jet_trk_pt_ns = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH2D*** h2_jet_trk_pt_ns_cov = Get2DArray <TH2D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_perp = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH2D*** h2_jet_trk_pt_perp_cov = Get2DArray <TH2D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_as = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH2D*** h2_jet_trk_pt_as_cov = Get2DArray <TH2D*> (nZdcCentBins, nVar);

  TH1D*** h_jet_trk_pt_ns_bkg = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH2D*** h2_jet_trk_pt_ns_cov_bkg = Get2DArray <TH2D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_perp_bkg = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH2D*** h2_jet_trk_pt_perp_cov_bkg = Get2DArray <TH2D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_as_bkg = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH2D*** h2_jet_trk_pt_as_cov_bkg = Get2DArray <TH2D*> (nZdcCentBins, nVar);

  TH1D**  h_jet_trk_pt_ns_ref_sig = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_trk_pt_perp_ref_sig = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_trk_pt_as_ref_sig = Get1DArray <TH1D*> (nVar);

  TH1D*** h_jet_trk_pt_ns_sig = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_perp_sig = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_as_sig = Get2DArray <TH1D*> (nZdcCentBins, nVar);

  TH1D*** h_jet_trk_pt_ns_iaa = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_perp_iaa = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_trk_pt_as_iaa = Get2DArray <TH1D*> (nZdcCentBins, nVar);


  TGAE**  g_jet_trk_pt_ns_ref_syst = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_perp_ref_syst = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_as_ref_syst = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_ns_ref_bkg_syst = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_perp_ref_bkg_syst = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_as_ref_bkg_syst = Get1DArray <TGAE*> (nVar);
  TGAE*** g_jet_trk_pt_ns_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_perp_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_as_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_ns_bkg_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_perp_bkg_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_as_bkg_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);

  TGAE**  g_jet_trk_pt_ns_ref_sig_syst = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_perp_ref_sig_syst = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_as_ref_sig_syst = Get1DArray <TGAE*> (nVar);
  TGAE*** g_jet_trk_pt_ns_sig_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_perp_sig_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_as_sig_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);

  TGAE*** g_jet_trk_pt_ns_iaa_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_perp_iaa_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_as_iaa_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/PlotPtCh_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  for (int iVar = 0; iVar < nVar; iVar++) {

    {
      TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data17_5TeV_hists.root", rootPath.Data (), tag1, variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts_ref[iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", tag1))->Clone (Form ("h_evt_counts_ref_%s", variations[iVar].Data ()));
      h_jet_counts_ref[iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", tag1))->Clone (Form ("h_jet_counts_ref_%s", variations[iVar].Data ()));

      h_jet_trk_pt_ns_ref[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data17", tag1))->Clone (Form ("h_jet_trk_pt_ns_ref_%s", variations[iVar].Data ()));
      h2_jet_trk_pt_ns_cov_ref[iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data17", tag1))->Clone (Form ("h2_jet_trk_pt_ns_cov_ref_%s", variations[iVar].Data ()));
      h_jet_trk_pt_perp_ref[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_data17", tag1))->Clone (Form ("h_jet_trk_pt_perp_ref_%s", variations[iVar].Data ()));
      h2_jet_trk_pt_perp_cov_ref[iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_data17", tag1))->Clone (Form ("h2_jet_trk_pt_perp_cov_ref_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_ref[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data17", tag1))->Clone (Form ("h_jet_trk_pt_as_ref_%s", variations[iVar].Data ()));
      h2_jet_trk_pt_as_cov_ref[iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data17", tag1))->Clone (Form ("h2_jet_trk_pt_as_cov_ref_%s", variations[iVar].Data ()));

      inFile->Close ();

      CalcUncertainties (h_jet_trk_pt_ns_ref[iVar], h2_jet_trk_pt_ns_cov_ref[iVar], h_jet_counts_ref[iVar]);
      CalcUncertainties (h_jet_trk_pt_perp_ref[iVar], h2_jet_trk_pt_perp_cov_ref[iVar], h_jet_counts_ref[iVar]);
      CalcUncertainties (h_jet_trk_pt_as_ref[iVar], h2_jet_trk_pt_as_cov_ref[iVar], h_jet_counts_ref[iVar]);
    }



    {
      TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/data17_5TeV_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts_ref_bkg[iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_evt_counts_ref_bkg_%s", variations[iVar].Data ()));
      h_jet_counts_ref_bkg[iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_counts_ref_bkg_%s", variations[iVar].Data ()));

      h_jet_trk_pt_ns_ref_bkg[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_ns_ref_bkg_%s", variations[iVar].Data ()));
      h2_jet_trk_pt_ns_cov_ref_bkg[iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_ns_cov_ref_bkg_%s", variations[iVar].Data ()));
      h_jet_trk_pt_perp_ref_bkg[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_perp_ref_bkg_%s", variations[iVar].Data ()));
      h2_jet_trk_pt_perp_cov_ref_bkg[iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_perp_cov_ref_bkg_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_ref_bkg[iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_as_ref_bkg_%s", variations[iVar].Data ()));
      h2_jet_trk_pt_as_cov_ref_bkg[iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_as_cov_ref_bkg_%s", variations[iVar].Data ()));

      inFile->Close ();

      CalcUncertainties (h_jet_trk_pt_ns_ref_bkg[iVar], h2_jet_trk_pt_ns_cov_ref_bkg[iVar], h_jet_counts_ref_bkg[iVar]);
      CalcUncertainties (h_jet_trk_pt_perp_ref_bkg[iVar], h2_jet_trk_pt_perp_cov_ref_bkg[iVar], h_jet_counts_ref_bkg[iVar]);
      CalcUncertainties (h_jet_trk_pt_as_ref_bkg[iVar], h2_jet_trk_pt_as_cov_ref_bkg[iVar], h_jet_counts_ref_bkg[iVar]);
    }



    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), tag1, variations[iVar].Data (), iCent);
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", tag1))->Clone (Form ("h_evt_counts_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_counts[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", tag1))->Clone (Form ("h_jet_counts_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));

      h_jet_trk_pt_ns[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data16", tag1))->Clone (Form ("h_jet_trk_pt_ns_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h2_jet_trk_pt_ns_cov[iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data16", tag1))->Clone (Form ("h2_jet_trk_pt_ns_cov_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_trk_pt_perp[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_data16", tag1))->Clone (Form ("h_jet_trk_pt_perp_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h2_jet_trk_pt_perp_cov[iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_data16", tag1))->Clone (Form ("h2_jet_trk_pt_perp_cov_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_trk_pt_as[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data16", tag1))->Clone (Form ("h_jet_trk_pt_as_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h2_jet_trk_pt_as_cov[iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data16", tag1))->Clone (Form ("h2_jet_trk_pt_as_cov_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));

      inFile->Close ();

      CalcUncertainties (h_jet_trk_pt_ns[iCent][iVar], h2_jet_trk_pt_ns_cov[iCent][iVar], h_jet_counts[iCent][iVar]);
      CalcUncertainties (h_jet_trk_pt_perp[iCent][iVar], h2_jet_trk_pt_perp_cov[iCent][iVar], h_jet_counts[iCent][iVar]);
      CalcUncertainties (h_jet_trk_pt_as[iCent][iVar], h2_jet_trk_pt_as_cov[iCent][iVar], h_jet_counts[iCent][iVar]);
    }



    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", variations[iVar].Data (), iCent);
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts_bkg[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_evt_counts_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_counts_bkg[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_counts_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));

      h_jet_trk_pt_ns_bkg[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_ns_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));
      h2_jet_trk_pt_ns_cov_bkg[iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_ns_cov_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_trk_pt_perp_bkg[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_perp_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));
      h2_jet_trk_pt_perp_cov_bkg[iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_perp_cov_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_trk_pt_as_bkg[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_as_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));
      h2_jet_trk_pt_as_cov_bkg[iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_as_cov_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));

      inFile->Close ();

      CalcUncertainties (h_jet_trk_pt_ns_bkg[iCent][iVar], h2_jet_trk_pt_ns_cov_bkg[iCent][iVar], h_jet_counts_bkg[iCent][iVar]);
      CalcUncertainties (h_jet_trk_pt_perp_bkg[iCent][iVar], h2_jet_trk_pt_perp_cov_bkg[iCent][iVar], h_jet_counts_bkg[iCent][iVar]);
      CalcUncertainties (h_jet_trk_pt_as_bkg[iCent][iVar], h2_jet_trk_pt_as_cov_bkg[iCent][iVar], h_jet_counts_bkg[iCent][iVar]);
    }



    {
      h_jet_trk_pt_ns_ref_sig[iVar] = (TH1D*) h_jet_trk_pt_ns_ref[iVar]->Clone (Form ("h_jet_trk_pt_ns_ref_sig_%s", variations[iVar].Data ()));
      h_jet_trk_pt_ns_ref_sig[iVar]->Add (h_jet_trk_pt_ns_ref_bkg[iVar], -1);

      h_jet_trk_pt_perp_ref_sig[iVar] = (TH1D*) h_jet_trk_pt_perp_ref[iVar]->Clone (Form ("h_jet_trk_pt_perp_ref_sig_%s", variations[iVar].Data ()));
      h_jet_trk_pt_perp_ref_sig[iVar]->Add (h_jet_trk_pt_perp_ref_bkg[iVar], -1);

      h_jet_trk_pt_as_ref_sig[iVar] = (TH1D*) h_jet_trk_pt_as_ref[iVar]->Clone (Form ("h_jet_trk_pt_as_ref_sig_%s", variations[iVar].Data ()));
      h_jet_trk_pt_as_ref_sig[iVar]->Add (h_jet_trk_pt_as_ref_bkg[iVar], -1);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_jet_trk_pt_ns_sig[iCent][iVar] = (TH1D*) h_jet_trk_pt_ns[iCent][iVar]->Clone (Form ("h_jet_trk_pt_ns_pPb_sig_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_ns_sig[iCent][iVar]->Add (h_jet_trk_pt_ns_bkg[iCent][iVar], -1);

        h_jet_trk_pt_perp_sig[iCent][iVar] = (TH1D*) h_jet_trk_pt_perp[iCent][iVar]->Clone (Form ("h_jet_trk_pt_perp_pPb_sig_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_perp_sig[iCent][iVar]->Add (h_jet_trk_pt_perp_bkg[iCent][iVar], -1);

        h_jet_trk_pt_as_sig[iCent][iVar] = (TH1D*) h_jet_trk_pt_as[iCent][iVar]->Clone (Form ("h_jet_trk_pt_as_pPb_sig_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_as_sig[iCent][iVar]->Add (h_jet_trk_pt_as_bkg[iCent][iVar], -1);

        h_jet_trk_pt_ns_iaa[iCent][iVar] = (TH1D*) h_jet_trk_pt_ns_sig[iCent][iVar]->Clone (Form ("h_jet_trk_pt_ns_iaa_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_ns_iaa[iCent][iVar]->Divide (h_jet_trk_pt_ns_ref_sig[iVar]);

        h_jet_trk_pt_perp_iaa[iCent][iVar] = (TH1D*) h_jet_trk_pt_perp_sig[iCent][iVar]->Clone (Form ("h_jet_trk_pt_perp_iaa_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_perp_iaa[iCent][iVar]->Divide (h_jet_trk_pt_perp_ref_sig[iVar]);

        h_jet_trk_pt_as_iaa[iCent][iVar] = (TH1D*) h_jet_trk_pt_as_sig[iCent][iVar]->Clone (Form ("h_jet_trk_pt_as_iaa_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_trk_pt_as_iaa[iCent][iVar]->Divide (h_jet_trk_pt_as_ref_sig[iVar]);

      }
    }
  }



  {
    g_jet_trk_pt_ns_ref_syst[0] = make_graph (h_jet_trk_pt_ns_ref[0]);
    g_jet_trk_pt_perp_ref_syst[0] = make_graph (h_jet_trk_pt_perp_ref[0]);
    g_jet_trk_pt_as_ref_syst[0] = make_graph (h_jet_trk_pt_as_ref[0]);

    g_jet_trk_pt_ns_ref_bkg_syst[0] = make_graph (h_jet_trk_pt_ns_ref_bkg[0]);
    g_jet_trk_pt_perp_ref_bkg_syst[0] = make_graph (h_jet_trk_pt_perp_ref_bkg[0]);
    g_jet_trk_pt_as_ref_bkg_syst[0] = make_graph (h_jet_trk_pt_as_ref_bkg[0]);

    g_jet_trk_pt_ns_ref_sig_syst[0] = make_graph (h_jet_trk_pt_ns_ref_sig[0]);
    g_jet_trk_pt_perp_ref_sig_syst[0] = make_graph (h_jet_trk_pt_perp_ref_sig[0]);
    g_jet_trk_pt_as_ref_sig_syst[0] = make_graph (h_jet_trk_pt_as_ref_sig[0]);

    ResetTGAEErrors (g_jet_trk_pt_ns_ref_syst[0]);
    ResetTGAEErrors (g_jet_trk_pt_perp_ref_syst[0]);
    ResetTGAEErrors (g_jet_trk_pt_as_ref_syst[0]);

    ResetTGAEErrors (g_jet_trk_pt_ns_ref_bkg_syst[0]);
    ResetTGAEErrors (g_jet_trk_pt_perp_ref_bkg_syst[0]);
    ResetTGAEErrors (g_jet_trk_pt_as_ref_bkg_syst[0]);

    ResetTGAEErrors (g_jet_trk_pt_ns_ref_sig_syst[0]);
    ResetTGAEErrors (g_jet_trk_pt_perp_ref_sig_syst[0]);
    ResetTGAEErrors (g_jet_trk_pt_as_ref_sig_syst[0]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      g_jet_trk_pt_ns_syst[iCent][0] = make_graph (h_jet_trk_pt_ns[iCent][0]);
      g_jet_trk_pt_perp_syst[iCent][0] = make_graph (h_jet_trk_pt_perp[iCent][0]);
      g_jet_trk_pt_as_syst[iCent][0] = make_graph (h_jet_trk_pt_as[iCent][0]);

      g_jet_trk_pt_ns_bkg_syst[iCent][0] = make_graph (h_jet_trk_pt_ns_bkg[iCent][0]);
      g_jet_trk_pt_perp_bkg_syst[iCent][0] = make_graph (h_jet_trk_pt_perp_bkg[iCent][0]);
      g_jet_trk_pt_as_bkg_syst[iCent][0] = make_graph (h_jet_trk_pt_as_bkg[iCent][0]);

      g_jet_trk_pt_ns_sig_syst[iCent][0] = make_graph (h_jet_trk_pt_ns_sig[iCent][0]);
      g_jet_trk_pt_perp_sig_syst[iCent][0] = make_graph (h_jet_trk_pt_perp_sig[iCent][0]);
      g_jet_trk_pt_as_sig_syst[iCent][0] = make_graph (h_jet_trk_pt_as_sig[iCent][0]);

      g_jet_trk_pt_ns_iaa_syst[iCent][0] = make_graph (h_jet_trk_pt_ns_iaa[iCent][0]);
      g_jet_trk_pt_perp_iaa_syst[iCent][0] = make_graph (h_jet_trk_pt_perp_iaa[iCent][0]);
      g_jet_trk_pt_as_iaa_syst[iCent][0] = make_graph (h_jet_trk_pt_as_iaa[iCent][0]);

      ResetTGAEErrors (g_jet_trk_pt_ns_syst[iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_perp_syst[iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_as_syst[iCent][0]);

      ResetTGAEErrors (g_jet_trk_pt_ns_bkg_syst[iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_perp_bkg_syst[iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_as_bkg_syst[iCent][0]);

      ResetTGAEErrors (g_jet_trk_pt_ns_sig_syst[iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_perp_sig_syst[iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_as_sig_syst[iCent][0]);

      ResetTGAEErrors (g_jet_trk_pt_ns_iaa_syst[iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_perp_iaa_syst[iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_as_iaa_syst[iCent][0]);

    }
  }



  for (int iVar = 1; iVar < nVar; iVar++) {

    g_jet_trk_pt_ns_ref_syst[iVar] = new TGAE ();
    g_jet_trk_pt_perp_ref_syst[iVar] = new TGAE ();
    g_jet_trk_pt_as_ref_syst[iVar] = new TGAE ();

    g_jet_trk_pt_ns_ref_bkg_syst[iVar] = new TGAE ();
    g_jet_trk_pt_perp_ref_bkg_syst[iVar] = new TGAE ();
    g_jet_trk_pt_as_ref_bkg_syst[iVar] = new TGAE ();

    g_jet_trk_pt_ns_ref_sig_syst[iVar] = new TGAE ();
    g_jet_trk_pt_perp_ref_sig_syst[iVar] = new TGAE ();
    g_jet_trk_pt_as_ref_sig_syst[iVar] = new TGAE ();

    CalcSystematics (g_jet_trk_pt_ns_ref_syst[iVar], h_jet_trk_pt_ns_ref[0], h_jet_trk_pt_ns_ref[iVar]);
    CalcSystematics (g_jet_trk_pt_perp_ref_syst[iVar], h_jet_trk_pt_perp_ref[0], h_jet_trk_pt_perp_ref[iVar]);
    CalcSystematics (g_jet_trk_pt_as_ref_syst[iVar], h_jet_trk_pt_as_ref[0], h_jet_trk_pt_as_ref[iVar]);

    CalcSystematics (g_jet_trk_pt_ns_ref_bkg_syst[iVar], h_jet_trk_pt_ns_ref_bkg[0], h_jet_trk_pt_ns_ref_bkg[iVar]);
    CalcSystematics (g_jet_trk_pt_perp_ref_bkg_syst[iVar], h_jet_trk_pt_perp_ref_bkg[0], h_jet_trk_pt_perp_ref_bkg[iVar]);
    CalcSystematics (g_jet_trk_pt_as_ref_bkg_syst[iVar], h_jet_trk_pt_as_ref_bkg[0], h_jet_trk_pt_as_ref_bkg[iVar]);

    CalcSystematics (g_jet_trk_pt_ns_ref_sig_syst[iVar], h_jet_trk_pt_ns_ref_sig[0], h_jet_trk_pt_ns_ref_sig[iVar]);
    CalcSystematics (g_jet_trk_pt_perp_ref_sig_syst[iVar], h_jet_trk_pt_perp_ref_sig[0], h_jet_trk_pt_perp_ref_sig[iVar]);
    CalcSystematics (g_jet_trk_pt_as_ref_sig_syst[iVar], h_jet_trk_pt_as_ref_sig[0], h_jet_trk_pt_as_ref_sig[iVar]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      g_jet_trk_pt_ns_syst[iCent][iVar] = new TGAE ();
      g_jet_trk_pt_perp_syst[iCent][iVar] = new TGAE ();
      g_jet_trk_pt_as_syst[iCent][iVar] = new TGAE ();

      g_jet_trk_pt_ns_bkg_syst[iCent][iVar] = new TGAE ();
      g_jet_trk_pt_perp_bkg_syst[iCent][iVar] = new TGAE ();
      g_jet_trk_pt_as_bkg_syst[iCent][iVar] = new TGAE ();

      g_jet_trk_pt_ns_sig_syst[iCent][iVar] = new TGAE ();
      g_jet_trk_pt_perp_sig_syst[iCent][iVar] = new TGAE ();
      g_jet_trk_pt_as_sig_syst[iCent][iVar] = new TGAE ();

      g_jet_trk_pt_ns_iaa_syst[iCent][iVar] = new TGAE ();
      g_jet_trk_pt_perp_iaa_syst[iCent][iVar] = new TGAE ();
      g_jet_trk_pt_as_iaa_syst[iCent][iVar] = new TGAE ();

      CalcSystematics (g_jet_trk_pt_ns_syst[iCent][iVar], h_jet_trk_pt_ns[iCent][0], h_jet_trk_pt_ns[iCent][iVar]);
      CalcSystematics (g_jet_trk_pt_perp_syst[iCent][iVar], h_jet_trk_pt_perp[iCent][0], h_jet_trk_pt_perp[iCent][iVar]);
      CalcSystematics (g_jet_trk_pt_as_syst[iCent][iVar], h_jet_trk_pt_as[iCent][0], h_jet_trk_pt_as[iCent][iVar]);

      CalcSystematics (g_jet_trk_pt_ns_bkg_syst[iCent][iVar], h_jet_trk_pt_ns_bkg[iCent][0], h_jet_trk_pt_ns_bkg[iCent][iVar]);
      CalcSystematics (g_jet_trk_pt_perp_bkg_syst[iCent][iVar], h_jet_trk_pt_perp_bkg[iCent][0], h_jet_trk_pt_perp_bkg[iCent][iVar]);
      CalcSystematics (g_jet_trk_pt_as_bkg_syst[iCent][iVar], h_jet_trk_pt_as_bkg[iCent][0], h_jet_trk_pt_as_bkg[iCent][iVar]);

      CalcSystematics (g_jet_trk_pt_ns_sig_syst[iCent][iVar], h_jet_trk_pt_ns_sig[iCent][0], h_jet_trk_pt_ns_sig[iCent][iVar]);
      CalcSystematics (g_jet_trk_pt_perp_sig_syst[iCent][iVar], h_jet_trk_pt_perp_sig[iCent][0], h_jet_trk_pt_perp_sig[iCent][iVar]);
      CalcSystematics (g_jet_trk_pt_as_sig_syst[iCent][iVar], h_jet_trk_pt_as_sig[iCent][0], h_jet_trk_pt_as_sig[iCent][iVar]);

      CalcSystematics (g_jet_trk_pt_ns_iaa_syst[iCent][iVar], h_jet_trk_pt_ns_iaa[iCent][0], h_jet_trk_pt_ns_iaa[iCent][iVar]);
      CalcSystematics (g_jet_trk_pt_perp_iaa_syst[iCent][iVar], h_jet_trk_pt_perp_iaa[iCent][0], h_jet_trk_pt_perp_iaa[iCent][iVar]);
      CalcSystematics (g_jet_trk_pt_as_iaa_syst[iCent][iVar], h_jet_trk_pt_as_iaa[iCent][0], h_jet_trk_pt_as_iaa[iCent][iVar]);
    }
  }



  // takes the maximum variation of the jet pT ES up/down variations
  {
    const int syst = 7; // MixCatVar1

    AddMaxSystematic (g_jet_trk_pt_ns_ref_syst[0], g_jet_trk_pt_ns_ref_syst[syst], g_jet_trk_pt_ns_ref_syst[syst]);
    AddMaxSystematic (g_jet_trk_pt_perp_ref_syst[0], g_jet_trk_pt_perp_ref_syst[syst], g_jet_trk_pt_perp_ref_syst[syst]);
    AddMaxSystematic (g_jet_trk_pt_as_ref_syst[0], g_jet_trk_pt_as_ref_syst[syst], g_jet_trk_pt_as_ref_syst[syst]);

    AddMaxSystematic (g_jet_trk_pt_ns_ref_bkg_syst[0], g_jet_trk_pt_ns_ref_bkg_syst[syst], g_jet_trk_pt_ns_ref_bkg_syst[syst]);
    AddMaxSystematic (g_jet_trk_pt_perp_ref_bkg_syst[0], g_jet_trk_pt_perp_ref_bkg_syst[syst], g_jet_trk_pt_perp_ref_bkg_syst[syst]);
    AddMaxSystematic (g_jet_trk_pt_as_ref_bkg_syst[0], g_jet_trk_pt_as_ref_bkg_syst[syst], g_jet_trk_pt_as_ref_bkg_syst[syst]);

    AddMaxSystematic (g_jet_trk_pt_ns_ref_sig_syst[0], g_jet_trk_pt_ns_ref_sig_syst[syst], g_jet_trk_pt_ns_ref_sig_syst[syst]);
    AddMaxSystematic (g_jet_trk_pt_perp_ref_sig_syst[0], g_jet_trk_pt_perp_ref_sig_syst[syst], g_jet_trk_pt_perp_ref_sig_syst[syst]);
    AddMaxSystematic (g_jet_trk_pt_as_ref_sig_syst[0], g_jet_trk_pt_as_ref_sig_syst[syst], g_jet_trk_pt_as_ref_sig_syst[syst]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      AddMaxSystematic (g_jet_trk_pt_ns_syst[iCent][0], g_jet_trk_pt_ns_syst[iCent][syst], g_jet_trk_pt_ns_syst[iCent][syst]);
      AddMaxSystematic (g_jet_trk_pt_perp_syst[iCent][0], g_jet_trk_pt_perp_syst[iCent][syst], g_jet_trk_pt_perp_syst[iCent][syst]);
      AddMaxSystematic (g_jet_trk_pt_as_syst[iCent][0], g_jet_trk_pt_as_syst[iCent][syst], g_jet_trk_pt_as_syst[iCent][syst]);

      AddMaxSystematic (g_jet_trk_pt_ns_bkg_syst[iCent][0], g_jet_trk_pt_ns_bkg_syst[iCent][syst], g_jet_trk_pt_ns_bkg_syst[iCent][syst]);
      AddMaxSystematic (g_jet_trk_pt_perp_bkg_syst[iCent][0], g_jet_trk_pt_perp_bkg_syst[iCent][syst], g_jet_trk_pt_perp_bkg_syst[iCent][syst]);
      AddMaxSystematic (g_jet_trk_pt_as_bkg_syst[iCent][0], g_jet_trk_pt_as_bkg_syst[iCent][syst], g_jet_trk_pt_as_bkg_syst[iCent][syst]);

      AddMaxSystematic (g_jet_trk_pt_ns_sig_syst[iCent][0], g_jet_trk_pt_ns_sig_syst[iCent][syst], g_jet_trk_pt_ns_sig_syst[iCent][syst]);
      AddMaxSystematic (g_jet_trk_pt_perp_sig_syst[iCent][0], g_jet_trk_pt_perp_sig_syst[iCent][syst], g_jet_trk_pt_perp_sig_syst[iCent][syst]);
      AddMaxSystematic (g_jet_trk_pt_as_sig_syst[iCent][0], g_jet_trk_pt_as_sig_syst[iCent][syst], g_jet_trk_pt_as_sig_syst[iCent][syst]);

      AddMaxSystematic (g_jet_trk_pt_ns_iaa_syst[iCent][0], g_jet_trk_pt_ns_iaa_syst[iCent][syst], g_jet_trk_pt_ns_iaa_syst[iCent][syst]);
      AddMaxSystematic (g_jet_trk_pt_perp_iaa_syst[iCent][0], g_jet_trk_pt_perp_iaa_syst[iCent][syst], g_jet_trk_pt_perp_iaa_syst[iCent][syst]);
      AddMaxSystematic (g_jet_trk_pt_as_iaa_syst[iCent][0], g_jet_trk_pt_as_iaa_syst[iCent][syst], g_jet_trk_pt_as_iaa_syst[iCent][syst]);

    }

    //const int syst = 1; // 5% JES systematic

    //AddMaxSystematic (g_jet_trk_pt_ns_ref_syst[0], g_jet_trk_pt_ns_ref_syst[syst], g_jet_trk_pt_ns_ref_syst[syst+1]);
    //AddMaxSystematic (g_jet_trk_pt_perp_ref_syst[0], g_jet_trk_pt_perp_ref_syst[syst], g_jet_trk_pt_perp_ref_syst[syst+1]);
    //AddMaxSystematic (g_jet_trk_pt_as_ref_syst[0], g_jet_trk_pt_as_ref_syst[syst], g_jet_trk_pt_as_ref_syst[syst+1]);

    //AddMaxSystematic (g_jet_trk_pt_ns_ref_bkg_syst[0], g_jet_trk_pt_ns_ref_bkg_syst[syst], g_jet_trk_pt_ns_ref_bkg_syst[syst+1]);
    //AddMaxSystematic (g_jet_trk_pt_perp_ref_bkg_syst[0], g_jet_trk_pt_perp_ref_bkg_syst[syst], g_jet_trk_pt_perp_ref_bkg_syst[syst+1]);
    //AddMaxSystematic (g_jet_trk_pt_as_ref_bkg_syst[0], g_jet_trk_pt_as_ref_bkg_syst[syst], g_jet_trk_pt_as_ref_bkg_syst[syst+1]);

    //AddMaxSystematic (g_jet_trk_pt_ns_ref_sig_syst[0], g_jet_trk_pt_ns_ref_sig_syst[syst], g_jet_trk_pt_ns_ref_sig_syst[syst+1]);
    //AddMaxSystematic (g_jet_trk_pt_perp_ref_sig_syst[0], g_jet_trk_pt_perp_ref_sig_syst[syst], g_jet_trk_pt_perp_ref_sig_syst[syst+1]);
    //AddMaxSystematic (g_jet_trk_pt_as_ref_sig_syst[0], g_jet_trk_pt_as_ref_sig_syst[syst], g_jet_trk_pt_as_ref_sig_syst[syst+1]);

    //for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

    //  AddMaxSystematic (g_jet_trk_pt_ns_syst[iCent][0], g_jet_trk_pt_ns_syst[iCent][syst], g_jet_trk_pt_ns_syst[iCent][syst+1]);
    //  AddMaxSystematic (g_jet_trk_pt_perp_syst[iCent][0], g_jet_trk_pt_perp_syst[iCent][syst], g_jet_trk_pt_perp_syst[iCent][syst+1]);
    //  AddMaxSystematic (g_jet_trk_pt_as_syst[iCent][0], g_jet_trk_pt_as_syst[iCent][syst], g_jet_trk_pt_as_syst[iCent][syst+1]);

    //  AddMaxSystematic (g_jet_trk_pt_ns_bkg_syst[iCent][0], g_jet_trk_pt_ns_bkg_syst[iCent][syst], g_jet_trk_pt_ns_bkg_syst[iCent][syst+1]);
    //  AddMaxSystematic (g_jet_trk_pt_perp_bkg_syst[iCent][0], g_jet_trk_pt_perp_bkg_syst[iCent][syst], g_jet_trk_pt_perp_bkg_syst[iCent][syst+1]);
    //  AddMaxSystematic (g_jet_trk_pt_as_bkg_syst[iCent][0], g_jet_trk_pt_as_bkg_syst[iCent][syst], g_jet_trk_pt_as_bkg_syst[iCent][syst+1]);

    //  AddMaxSystematic (g_jet_trk_pt_ns_sig_syst[iCent][0], g_jet_trk_pt_ns_sig_syst[iCent][syst], g_jet_trk_pt_ns_sig_syst[iCent][syst+1]);
    //  AddMaxSystematic (g_jet_trk_pt_perp_sig_syst[iCent][0], g_jet_trk_pt_perp_sig_syst[iCent][syst], g_jet_trk_pt_perp_sig_syst[iCent][syst+1]);
    //  AddMaxSystematic (g_jet_trk_pt_as_sig_syst[iCent][0], g_jet_trk_pt_as_sig_syst[iCent][syst], g_jet_trk_pt_as_sig_syst[iCent][syst+1]);

    //  AddMaxSystematic (g_jet_trk_pt_ns_iaa_syst[iCent][0], g_jet_trk_pt_ns_iaa_syst[iCent][syst], g_jet_trk_pt_ns_iaa_syst[iCent][syst+1]);
    //  AddMaxSystematic (g_jet_trk_pt_perp_iaa_syst[iCent][0], g_jet_trk_pt_perp_iaa_syst[iCent][syst], g_jet_trk_pt_perp_iaa_syst[iCent][syst+1]);
    //  AddMaxSystematic (g_jet_trk_pt_as_iaa_syst[iCent][0], g_jet_trk_pt_as_iaa_syst[iCent][syst], g_jet_trk_pt_as_iaa_syst[iCent][syst+1]);

    //}
  }



  {
    outFile->cd ();

    g_jet_trk_pt_ns_ref_syst[0]->Write ("g_jet_trk_pt_ns_ref_syst");
    g_jet_trk_pt_perp_ref_syst[0]->Write ("g_jet_trk_pt_perp_ref_syst");
    g_jet_trk_pt_as_ref_syst[0]->Write ("g_jet_trk_pt_as_ref_syst");

    g_jet_trk_pt_ns_ref_bkg_syst[0]->Write ("g_jet_trk_pt_ns_ref_bkg_syst");
    g_jet_trk_pt_perp_ref_bkg_syst[0]->Write ("g_jet_trk_pt_perp_ref_bkg_syst");
    g_jet_trk_pt_as_ref_bkg_syst[0]->Write ("g_jet_trk_pt_as_ref_bkg_syst");

    g_jet_trk_pt_ns_ref_sig_syst[0]->Write ("g_jet_trk_pt_ns_ref_sig_syst");
    g_jet_trk_pt_perp_ref_sig_syst[0]->Write ("g_jet_trk_pt_perp_ref_sig_syst");
    g_jet_trk_pt_as_ref_sig_syst[0]->Write ("g_jet_trk_pt_as_ref_sig_syst");

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      g_jet_trk_pt_ns_syst[iCent][0]->Write (Form ("g_jet_trk_pt_ns_syst_iCent%i", iCent));
      g_jet_trk_pt_perp_syst[iCent][0]->Write (Form ("g_jet_trk_pt_perp_syst_iCent%i", iCent));
      g_jet_trk_pt_as_syst[iCent][0]->Write (Form ("g_jet_trk_pt_as_syst_iCent%i", iCent));

      g_jet_trk_pt_ns_bkg_syst[iCent][0]->Write (Form ("g_jet_trk_pt_ns_bkg_syst_iCent%i", iCent));
      g_jet_trk_pt_perp_bkg_syst[iCent][0]->Write (Form ("g_jet_trk_pt_perp_bkg_syst_iCent%i", iCent));
      g_jet_trk_pt_as_bkg_syst[iCent][0]->Write (Form ("g_jet_trk_pt_as_bkg_syst_iCent%i", iCent));

      g_jet_trk_pt_ns_sig_syst[iCent][0]->Write (Form ("g_jet_trk_pt_ns_sig_syst_iCent%i", iCent));
      g_jet_trk_pt_perp_sig_syst[iCent][0]->Write (Form ("g_jet_trk_pt_perp_sig_syst_iCent%i", iCent));
      g_jet_trk_pt_as_sig_syst[iCent][0]->Write (Form ("g_jet_trk_pt_as_sig_syst_iCent%i", iCent));

      g_jet_trk_pt_ns_iaa_syst[iCent][0]->Write (Form ("g_jet_trk_pt_ns_iaa_syst_iCent%i", iCent));
      g_jet_trk_pt_perp_iaa_syst[iCent][0]->Write (Form ("g_jet_trk_pt_perp_iaa_syst_iCent%i", iCent));
      g_jet_trk_pt_as_iaa_syst[iCent][0]->Write (Form ("g_jet_trk_pt_as_iaa_syst_iCent%i", iCent));

    }

    for (int iVar = 0; iVar < nVar; iVar++) {

      h_evt_counts_ref[iVar]->Write ();
      h_jet_counts_ref[iVar]->Write ();

      h_evt_counts_ref_bkg[iVar]->Write ();
      h_jet_counts_ref_bkg[iVar]->Write ();

      h_jet_trk_pt_ns_ref[iVar]->Write ();
      h_jet_trk_pt_perp_ref[iVar]->Write ();
      h_jet_trk_pt_as_ref[iVar]->Write ();

      h_jet_trk_pt_ns_ref_bkg[iVar]->Write ();
      h_jet_trk_pt_perp_ref_bkg[iVar]->Write ();
      h_jet_trk_pt_as_ref_bkg[iVar]->Write ();

      h_jet_trk_pt_ns_ref_sig[iVar]->Write ();
      h_jet_trk_pt_perp_ref_sig[iVar]->Write ();
      h_jet_trk_pt_as_ref_sig[iVar]->Write ();

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_evt_counts[iCent][iVar]->Write ();
        h_jet_counts[iCent][iVar]->Write ();

        h_jet_trk_pt_ns[iCent][iVar]->Write ();
        h_jet_trk_pt_perp[iCent][iVar]->Write ();
        h_jet_trk_pt_as[iCent][iVar]->Write ();

        h_jet_trk_pt_ns_bkg[iCent][iVar]->Write ();
        h_jet_trk_pt_perp_bkg[iCent][iVar]->Write ();
        h_jet_trk_pt_as_bkg[iCent][iVar]->Write ();

        h_jet_trk_pt_ns_sig[iCent][iVar]->Write ();
        h_jet_trk_pt_perp_sig[iCent][iVar]->Write ();
        h_jet_trk_pt_as_sig[iCent][iVar]->Write ();

        h_jet_trk_pt_ns_iaa[iCent][iVar]->Write ();
        h_jet_trk_pt_perp_iaa[iCent][iVar]->Write ();
        h_jet_trk_pt_as_iaa[iCent][iVar]->Write ();

      }
    }

    outFile->Close ();
  }
}


#endif
