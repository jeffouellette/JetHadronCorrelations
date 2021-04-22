#ifndef __JetHadronCorrelatorProcessPtCh_C__
#define __JetHadronCorrelatorProcessPtCh_C__

#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <math.h>

using namespace JetHadronCorrelations;

const int nVar = 8;
std::vector <TString> variations = {"Nominal", "JetES5PercUpVar", "JetES5PercDownVar", "JetES5PercSmearVar", "JetES2PercUpVar", "JetES2PercDownVar", "JetES2PercSmearVar", "FcalCentVar"};

void ProcessDPhi (const char* tag, const char* outFileTag) {

  SetupDirectories ("Data");

  TFile* inFile = nullptr;

  TH1D*** h_evt_counts = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_jet_counts = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_counts = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_counts = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_evt_counts_bkg = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_jet_counts_bkg = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_counts_bkg = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_counts_bkg = Get2DArray <TH1D*> (2, nVar);

  TH1D*** h_jet_trk_dphi = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_jet_trk_dphi_cov = Get2DArray <TH2D*> (2, nVar);
  TH1D*** h_ljet_trk_dphi = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_ljet_trk_dphi_cov = Get2DArray <TH2D*> (2, nVar);
  TH1D*** h_sljet_trk_dphi = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_sljet_trk_dphi_cov = Get2DArray <TH2D*> (2, nVar);
  TH1D*** h_jet_trk_dphi_bkg = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_jet_trk_dphi_cov_bkg = Get2DArray <TH2D*> (2, nVar);
  TH1D*** h_ljet_trk_dphi_bkg = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_ljet_trk_dphi_cov_bkg = Get2DArray <TH2D*> (2, nVar);
  TH1D*** h_sljet_trk_dphi_bkg = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_sljet_trk_dphi_cov_bkg = Get2DArray <TH2D*> (2, nVar);

  TH1D*** h_jet_trk_dphi_sig = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_ljet_trk_dphi_sig = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_sljet_trk_dphi_sig = Get2DArray <TH1D*> (2, nVar);

  TH1D** h_jet_trk_dphi_iaa = Get1DArray <TH1D*> (nVar);
  TH1D** h_ljet_trk_dphi_iaa = Get1DArray <TH1D*> (nVar);
  TH1D** h_sljet_trk_dphi_iaa = Get1DArray <TH1D*> (nVar);

  TGAE*** g_jet_trk_dphi_syst = Get2DArray <TGAE*> (2, nVar);
  TGAE*** g_ljet_trk_dphi_syst = Get2DArray <TGAE*> (2, nVar);
  TGAE*** g_sljet_trk_dphi_syst = Get2DArray <TGAE*> (2, nVar);
  TGAE*** g_jet_trk_dphi_bkg_syst = Get2DArray <TGAE*> (2, nVar);
  TGAE*** g_ljet_trk_dphi_bkg_syst = Get2DArray <TGAE*> (2, nVar);
  TGAE*** g_sljet_trk_dphi_bkg_syst = Get2DArray <TGAE*> (2, nVar);

  TGAE*** g_jet_trk_dphi_sig_syst = Get2DArray <TGAE*> (2, nVar);
  TGAE*** g_ljet_trk_dphi_sig_syst = Get2DArray <TGAE*> (2, nVar);
  TGAE*** g_sljet_trk_dphi_sig_syst = Get2DArray <TGAE*> (2, nVar);

  TGAE** g_jet_trk_dphi_iaa_syst = Get1DArray <TGAE*> (nVar);
  TGAE** g_ljet_trk_dphi_iaa_syst = Get1DArray <TGAE*> (nVar);
  TGAE** g_sljet_trk_dphi_iaa_syst = Get1DArray <TGAE*> (nVar);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("./rootFiles/Results/PlotDPhi_%s.root", outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  for (int iVar = 0; iVar < nVar; iVar++) {
    {
      TString inFileName = Form ("./rootFiles/%s/JetsHists/%s/data16_5TeV_hists.root", tag, variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts[1][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_%s_data16_%s", tag, variations[iVar].Data ()));
      h_jet_counts[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_%s_data16_%s", tag, variations[iVar].Data ()));
      h_ljet_counts[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_%s_data16_%s", tag, variations[iVar].Data ()));
      h_sljet_counts[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_data16", tag))->Clone (Form ("h_sljet_counts_%s_data16_%s", tag, variations[iVar].Data ()));
      h_jet_trk_dphi[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_data16", tag))->Clone (Form ("h_jet_trk_dphi_%s_data16_%s", tag, variations[iVar].Data ()));
      h2_jet_trk_dphi_cov[1][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_data16", tag))->Clone (Form ("h2_jet_trk_dphi_cov_%s_data16_%s", tag, variations[iVar].Data ()));
      h_ljet_trk_dphi[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_dphi_%s_data16", tag))->Clone (Form ("h_ljet_trk_dphi_%s_data16_%s", tag, variations[iVar].Data ()));
      h2_ljet_trk_dphi_cov[1][iVar] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_dphi_cov_%s_data16", tag))->Clone (Form ("h2_ljet_trk_dphi_cov_%s_data16_%s", tag, variations[iVar].Data ()));
      h_sljet_trk_dphi[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_dphi_%s_data16", tag))->Clone (Form ("h_sljet_trk_dphi_%s_data16_%s", tag, variations[iVar].Data ()));
      h2_sljet_trk_dphi_cov[1][iVar] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_dphi_cov_%s_data16", tag))->Clone (Form ("h2_sljet_trk_dphi_cov_%s_data16_%s", tag, variations[iVar].Data ()));

      inFile->Close ();

      CalcUncertainties (h_jet_trk_dphi[1][iVar], h2_jet_trk_dphi_cov[1][iVar], h_jet_counts[1][iVar]);
      CalcUncertainties (h_ljet_trk_dphi[1][iVar], h2_ljet_trk_dphi_cov[1][iVar], h_ljet_counts[1][iVar]);
      CalcUncertainties (h_sljet_trk_dphi[1][iVar], h2_sljet_trk_dphi_cov[1][iVar], h_sljet_counts[1][iVar]);
    }



    {
      TString inFileName = Form ("./rootFiles/%s/JetsHists/%s/data17_5TeV_hists.root", tag, variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_%s_data17_%s", tag, variations[iVar].Data ()));
      h_jet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_%s_data17_%s", tag, variations[iVar].Data ()));
      h_ljet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_%s_data17_%s", tag, variations[iVar].Data ()));
      h_sljet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_data17", tag))->Clone (Form ("h_sljet_counts_%s_data17_%s", tag, variations[iVar].Data ()));
      h_jet_trk_dphi[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_data17", tag))->Clone (Form ("h_jet_trk_dphi_%s_data17_%s", tag, variations[iVar].Data ()));
      h2_jet_trk_dphi_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_data17", tag))->Clone (Form ("h2_jet_trk_dphi_cov_%s_data17_%s", tag, variations[iVar].Data ()));
      h_ljet_trk_dphi[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_dphi_%s_data17", tag))->Clone (Form ("h_ljet_trk_dphi_%s_data17_%s", tag, variations[iVar].Data ()));
      h2_ljet_trk_dphi_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_dphi_cov_%s_data17", tag))->Clone (Form ("h2_ljet_trk_dphi_cov_%s_data17_%s", tag, variations[iVar].Data ()));
      h_sljet_trk_dphi[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_dphi_%s_data17", tag))->Clone (Form ("h_sljet_trk_dphi_%s_data17_%s", tag, variations[iVar].Data ()));
      h2_sljet_trk_dphi_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_dphi_cov_%s_data17", tag))->Clone (Form ("h2_sljet_trk_dphi_cov_%s_data17_%s", tag, variations[iVar].Data ()));

      inFile->Close ();

      CalcUncertainties (h_jet_trk_dphi[0][iVar], h2_jet_trk_dphi_cov[0][iVar], h_jet_counts[0][iVar]);
      CalcUncertainties (h_ljet_trk_dphi[0][iVar], h2_ljet_trk_dphi_cov[0][iVar], h_ljet_counts[0][iVar]);
      CalcUncertainties (h_sljet_trk_dphi[0][iVar], h2_sljet_trk_dphi_cov[0][iVar], h_sljet_counts[0][iVar]);
    }



    {
      TString inFileName = Form ("./rootFiles/%s/MixedHists/%s/data16_5TeV_hists.root", tag, variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_mixed_data16", tag))->Clone (Form ("h_evt_counts_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h_jet_counts_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_data16", tag))->Clone (Form ("h_evt_counts_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h_ljet_counts_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_mixed_data16", tag))->Clone (Form ("h_evt_counts_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h_sljet_counts_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_mixed_data16", tag))->Clone (Form ("h_sljet_counts_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h_jet_trk_dphi_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_mixed_data16", tag))->Clone (Form ("h_jet_trk_dphi_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h2_jet_trk_dphi_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_mixed_data16", tag))->Clone (Form ("h2_jet_trk_dphi_cov_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h_ljet_trk_dphi_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_dphi_%s_mixed_data16", tag))->Clone (Form ("h_ljet_trk_dphi_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h2_ljet_trk_dphi_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_dphi_cov_%s_mixed_data16", tag))->Clone (Form ("h2_ljet_trk_dphi_cov_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h_sljet_trk_dphi_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_dphi_%s_mixed_data16", tag))->Clone (Form ("h_sljet_trk_dphi_%s_mixed_data16_%s", tag, variations[iVar].Data ()));
      h2_sljet_trk_dphi_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_dphi_cov_%s_mixed_data16", tag))->Clone (Form ("h2_sljet_trk_dphi_cov_%s_mixed_data16_%s", tag, variations[iVar].Data ()));

      inFile->Close ();

      CalcUncertainties (h_jet_trk_dphi_bkg[1][iVar], h2_jet_trk_dphi_cov_bkg[1][iVar], h_jet_counts_bkg[1][iVar]);
      CalcUncertainties (h_ljet_trk_dphi_bkg[1][iVar], h2_ljet_trk_dphi_cov_bkg[1][iVar], h_ljet_counts_bkg[1][iVar]);
      CalcUncertainties (h_sljet_trk_dphi_bkg[1][iVar], h2_sljet_trk_dphi_cov_bkg[1][iVar], h_sljet_counts_bkg[1][iVar]);
    }



    {
      h_jet_trk_dphi_sig[0][iVar] = (TH1D*) h_jet_trk_dphi[0][iVar]->Clone (Form ("h_jet_trk_dphi_data17_sig_%s", variations[iVar].Data ()));
      h_jet_trk_dphi_sig[1][iVar] = (TH1D*) h_jet_trk_dphi[1][iVar]->Clone (Form ("h_jet_trk_dphi_data16_sig_%s", variations[iVar].Data ()));
      h_jet_trk_dphi_sig[1][iVar]->Add (h_jet_trk_dphi_bkg[1][iVar], -1);

      h_ljet_trk_dphi_sig[0][iVar] = (TH1D*) h_ljet_trk_dphi[0][iVar]->Clone (Form ("h_ljet_trk_dphi_data17_sig_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi_sig[1][iVar] = (TH1D*) h_ljet_trk_dphi[1][iVar]->Clone (Form ("h_ljet_trk_dphi_data16_sig_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi_sig[1][iVar]->Add (h_ljet_trk_dphi_bkg[1][iVar], -1);

      h_sljet_trk_dphi_sig[0][iVar] = (TH1D*) h_sljet_trk_dphi[0][iVar]->Clone (Form ("h_sljet_trk_dphi_data17_sig_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi_sig[1][iVar] = (TH1D*) h_sljet_trk_dphi[1][iVar]->Clone (Form ("h_sljet_trk_dphi_data16_sig_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi_sig[1][iVar]->Add (h_sljet_trk_dphi_bkg[1][iVar], -1);

      h_jet_trk_dphi_iaa[iVar] = (TH1D*) h_jet_trk_dphi_sig[1][iVar]->Clone (Form ("h_jet_trk_dphi_data16_iaa_%s", variations[iVar].Data ()));
      h_jet_trk_dphi_iaa[iVar]->Divide (h_jet_trk_dphi_sig[0][iVar]);

      h_ljet_trk_dphi_iaa[iVar] = (TH1D*) h_ljet_trk_dphi_sig[1][iVar]->Clone (Form ("h_ljet_trk_dphi_data16_iaa_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi_iaa[iVar]->Divide (h_ljet_trk_dphi_sig[0][iVar]);

      h_sljet_trk_dphi_iaa[iVar] = (TH1D*) h_sljet_trk_dphi_sig[1][iVar]->Clone (Form ("h_sljet_trk_dphi_data16_iaa_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi_iaa[iVar]->Divide (h_sljet_trk_dphi_sig[0][iVar]);
    }
  }



  {
    g_jet_trk_dphi_syst[0][0] = make_graph (h_jet_trk_dphi[0][0]);
    g_ljet_trk_dphi_syst[0][0] = make_graph (h_ljet_trk_dphi[0][0]);
    g_sljet_trk_dphi_syst[0][0] = make_graph (h_sljet_trk_dphi[0][0]);
    g_jet_trk_dphi_syst[1][0] = make_graph (h_jet_trk_dphi[1][0]);
    g_ljet_trk_dphi_syst[1][0] = make_graph (h_ljet_trk_dphi[1][0]);
    g_sljet_trk_dphi_syst[1][0] = make_graph (h_sljet_trk_dphi[1][0]);
    g_jet_trk_dphi_bkg_syst[1][0] = make_graph (h_jet_trk_dphi_bkg[1][0]);
    g_ljet_trk_dphi_bkg_syst[1][0] = make_graph (h_ljet_trk_dphi_bkg[1][0]);
    g_sljet_trk_dphi_bkg_syst[1][0] = make_graph (h_sljet_trk_dphi_bkg[1][0]);

    g_jet_trk_dphi_sig_syst[0][0] = make_graph (h_jet_trk_dphi_sig[0][0]);
    g_ljet_trk_dphi_sig_syst[0][0] = make_graph (h_ljet_trk_dphi_sig[0][0]);
    g_sljet_trk_dphi_sig_syst[0][0] = make_graph (h_sljet_trk_dphi_sig[0][0]);
    g_jet_trk_dphi_sig_syst[1][0] = make_graph (h_jet_trk_dphi_sig[1][0]);
    g_ljet_trk_dphi_sig_syst[1][0] = make_graph (h_ljet_trk_dphi_sig[1][0]);
    g_sljet_trk_dphi_sig_syst[1][0] = make_graph (h_sljet_trk_dphi_sig[1][0]);

    g_jet_trk_dphi_iaa_syst[0] = make_graph (h_jet_trk_dphi_iaa[0]);
    g_ljet_trk_dphi_iaa_syst[0] = make_graph (h_ljet_trk_dphi_iaa[0]);
    g_sljet_trk_dphi_iaa_syst[0] = make_graph (h_sljet_trk_dphi_iaa[0]);

    ResetTGAEErrors (g_jet_trk_dphi_syst[0][0]);
    ResetTGAEErrors (g_ljet_trk_dphi_syst[0][0]);
    ResetTGAEErrors (g_sljet_trk_dphi_syst[0][0]);
    ResetTGAEErrors (g_jet_trk_dphi_syst[1][0]);
    ResetTGAEErrors (g_ljet_trk_dphi_syst[1][0]);
    ResetTGAEErrors (g_sljet_trk_dphi_syst[1][0]);
    ResetTGAEErrors (g_jet_trk_dphi_bkg_syst[1][0]);
    ResetTGAEErrors (g_ljet_trk_dphi_bkg_syst[1][0]);
    ResetTGAEErrors (g_sljet_trk_dphi_bkg_syst[1][0]);

    ResetTGAEErrors (g_jet_trk_dphi_sig_syst[0][0]);
    ResetTGAEErrors (g_ljet_trk_dphi_sig_syst[0][0]);
    ResetTGAEErrors (g_sljet_trk_dphi_sig_syst[0][0]);
    ResetTGAEErrors (g_jet_trk_dphi_sig_syst[1][0]);
    ResetTGAEErrors (g_ljet_trk_dphi_sig_syst[1][0]);
    ResetTGAEErrors (g_sljet_trk_dphi_sig_syst[1][0]);

    ResetTGAEErrors (g_jet_trk_dphi_iaa_syst[0]);
    ResetTGAEErrors (g_ljet_trk_dphi_iaa_syst[0]);
    ResetTGAEErrors (g_sljet_trk_dphi_iaa_syst[0]);
  }



  for (int iVar = 1; iVar < nVar; iVar++) {
    g_jet_trk_dphi_syst[0][iVar] = new TGAE ();
    g_ljet_trk_dphi_syst[0][iVar] = new TGAE ();
    g_sljet_trk_dphi_syst[0][iVar] = new TGAE ();
    g_jet_trk_dphi_syst[1][iVar] = new TGAE ();
    g_ljet_trk_dphi_syst[1][iVar] = new TGAE ();
    g_sljet_trk_dphi_syst[1][iVar] = new TGAE ();
    g_jet_trk_dphi_bkg_syst[1][iVar] = new TGAE ();
    g_ljet_trk_dphi_bkg_syst[1][iVar] = new TGAE ();
    g_sljet_trk_dphi_bkg_syst[1][iVar] = new TGAE ();

    g_jet_trk_dphi_sig_syst[0][iVar] = new TGAE ();
    g_ljet_trk_dphi_sig_syst[0][iVar] = new TGAE ();
    g_sljet_trk_dphi_sig_syst[0][iVar] = new TGAE ();
    g_jet_trk_dphi_sig_syst[1][iVar] = new TGAE ();
    g_ljet_trk_dphi_sig_syst[1][iVar] = new TGAE ();
    g_sljet_trk_dphi_sig_syst[1][iVar] = new TGAE ();

    g_jet_trk_dphi_iaa_syst[iVar] = new TGAE ();
    g_ljet_trk_dphi_iaa_syst[iVar] = new TGAE ();
    g_sljet_trk_dphi_iaa_syst[iVar] = new TGAE ();

    CalcSystematics (g_jet_trk_dphi_syst[0][iVar], h_jet_trk_dphi[0][0], h_jet_trk_dphi[0][iVar]);
    CalcSystematics (g_ljet_trk_dphi_syst[0][iVar], h_ljet_trk_dphi[0][0], h_ljet_trk_dphi[0][iVar]);
    CalcSystematics (g_sljet_trk_dphi_syst[0][iVar], h_sljet_trk_dphi[0][0], h_sljet_trk_dphi[0][iVar]);
    CalcSystematics (g_jet_trk_dphi_syst[1][iVar], h_jet_trk_dphi[1][0], h_jet_trk_dphi[1][iVar]);
    CalcSystematics (g_ljet_trk_dphi_syst[1][iVar], h_ljet_trk_dphi[1][0], h_ljet_trk_dphi[1][iVar]);
    CalcSystematics (g_sljet_trk_dphi_syst[1][iVar], h_sljet_trk_dphi[1][0], h_sljet_trk_dphi[1][iVar]);
    CalcSystematics (g_jet_trk_dphi_bkg_syst[1][iVar], h_jet_trk_dphi_bkg[1][0], h_jet_trk_dphi_bkg[1][iVar]);
    CalcSystematics (g_ljet_trk_dphi_bkg_syst[1][iVar], h_ljet_trk_dphi_bkg[1][0], h_ljet_trk_dphi_bkg[1][iVar]);
    CalcSystematics (g_sljet_trk_dphi_bkg_syst[1][iVar], h_sljet_trk_dphi_bkg[1][0], h_sljet_trk_dphi_bkg[1][iVar]);

    CalcSystematics (g_jet_trk_dphi_sig_syst[0][iVar], h_jet_trk_dphi_sig[0][0], h_jet_trk_dphi_sig[0][iVar]);
    CalcSystematics (g_ljet_trk_dphi_sig_syst[0][iVar], h_ljet_trk_dphi_sig[0][0], h_ljet_trk_dphi_sig[0][iVar]);
    CalcSystematics (g_sljet_trk_dphi_sig_syst[0][iVar], h_sljet_trk_dphi_sig[0][0], h_sljet_trk_dphi_sig[0][iVar]);
    CalcSystematics (g_jet_trk_dphi_sig_syst[1][iVar], h_jet_trk_dphi_sig[1][0], h_jet_trk_dphi_sig[1][iVar]);
    CalcSystematics (g_ljet_trk_dphi_sig_syst[1][iVar], h_ljet_trk_dphi_sig[1][0], h_ljet_trk_dphi_sig[1][iVar]);
    CalcSystematics (g_sljet_trk_dphi_sig_syst[1][iVar], h_sljet_trk_dphi_sig[1][0], h_sljet_trk_dphi_sig[1][iVar]);

    CalcSystematics (g_jet_trk_dphi_iaa_syst[iVar], h_jet_trk_dphi_iaa[0], h_jet_trk_dphi_iaa[iVar]);
    CalcSystematics (g_ljet_trk_dphi_iaa_syst[iVar], h_ljet_trk_dphi_iaa[0], h_ljet_trk_dphi_iaa[iVar]);
    CalcSystematics (g_sljet_trk_dphi_iaa_syst[iVar], h_sljet_trk_dphi_iaa[0], h_sljet_trk_dphi_iaa[iVar]);
  }



  // takes the maximum variation of the jet pT ES up/down variations
  // TODO
  {
    const int syst = 1;

    AddMaxSystematic (g_jet_trk_dphi_syst[0][0], g_jet_trk_dphi_syst[0][syst], g_jet_trk_dphi_syst[0][syst+1]);
    AddMaxSystematic (g_ljet_trk_dphi_syst[0][0], g_ljet_trk_dphi_syst[0][syst], g_ljet_trk_dphi_syst[0][syst+1]);
    AddMaxSystematic (g_sljet_trk_dphi_syst[0][0], g_sljet_trk_dphi_syst[0][syst], g_sljet_trk_dphi_syst[0][syst+1]);
    AddMaxSystematic (g_jet_trk_dphi_syst[1][0], g_jet_trk_dphi_syst[1][syst], g_jet_trk_dphi_syst[1][syst+1]);
    AddMaxSystematic (g_ljet_trk_dphi_syst[1][0], g_ljet_trk_dphi_syst[1][syst], g_ljet_trk_dphi_syst[1][syst+1]);
    AddMaxSystematic (g_sljet_trk_dphi_syst[1][0], g_sljet_trk_dphi_syst[1][syst], g_sljet_trk_dphi_syst[1][syst+1]);
    AddMaxSystematic (g_jet_trk_dphi_bkg_syst[1][0], g_jet_trk_dphi_bkg_syst[1][syst], g_jet_trk_dphi_bkg_syst[1][syst+1]);
    AddMaxSystematic (g_ljet_trk_dphi_bkg_syst[1][0], g_ljet_trk_dphi_bkg_syst[1][syst], g_ljet_trk_dphi_bkg_syst[1][syst+1]);
    AddMaxSystematic (g_sljet_trk_dphi_bkg_syst[1][0], g_sljet_trk_dphi_bkg_syst[1][syst], g_sljet_trk_dphi_bkg_syst[1][syst+1]);

    AddMaxSystematic (g_jet_trk_dphi_sig_syst[0][0], g_jet_trk_dphi_sig_syst[0][syst], g_jet_trk_dphi_sig_syst[0][syst+1]);
    AddMaxSystematic (g_ljet_trk_dphi_sig_syst[0][0], g_ljet_trk_dphi_sig_syst[0][syst], g_ljet_trk_dphi_sig_syst[0][syst+1]);
    AddMaxSystematic (g_sljet_trk_dphi_sig_syst[0][0], g_sljet_trk_dphi_sig_syst[0][syst], g_sljet_trk_dphi_sig_syst[0][syst+1]);
    AddMaxSystematic (g_jet_trk_dphi_sig_syst[1][0], g_jet_trk_dphi_sig_syst[1][syst], g_jet_trk_dphi_sig_syst[1][syst+1]);
    AddMaxSystematic (g_ljet_trk_dphi_sig_syst[1][0], g_ljet_trk_dphi_sig_syst[1][syst], g_ljet_trk_dphi_sig_syst[1][syst+1]);
    AddMaxSystematic (g_sljet_trk_dphi_sig_syst[1][0], g_sljet_trk_dphi_sig_syst[1][syst], g_sljet_trk_dphi_sig_syst[1][syst+1]);

    AddMaxSystematic (g_jet_trk_dphi_iaa_syst[0], g_jet_trk_dphi_iaa_syst[syst], g_jet_trk_dphi_iaa_syst[syst+1]);
    AddMaxSystematic (g_ljet_trk_dphi_iaa_syst[0], g_ljet_trk_dphi_iaa_syst[syst], g_ljet_trk_dphi_iaa_syst[syst+1]);
    AddMaxSystematic (g_sljet_trk_dphi_iaa_syst[0], g_sljet_trk_dphi_iaa_syst[syst], g_sljet_trk_dphi_iaa_syst[syst+1]);
  }



  {
    outFile->cd ();

    h_evt_counts[0][0]->Write ("h_evt_counts_ref");
    h_jet_counts[0][0]->Write ("h_jet_counts_ref");
    h_ljet_counts[0][0]->Write ("h_ljet_counts_ref");
    h_sljet_counts[0][0]->Write ("h_sljet_counts_ref");
    h_evt_counts[1][0]->Write ("h_evt_counts");
    h_jet_counts[1][0]->Write ("h_jet_counts");
    h_ljet_counts[1][0]->Write ("h_ljet_counts");
    h_sljet_counts[1][0]->Write ("h_sljet_counts");
    h_evt_counts_bkg[1][0]->Write ("h_evt_counts_bkg");
    h_jet_counts_bkg[1][0]->Write ("h_jet_counts_bkg");
    h_ljet_counts_bkg[1][0]->Write ("h_ljet_counts_bkg");
    h_sljet_counts_bkg[1][0]->Write ("h_sljet_counts_bkg");

    h_jet_trk_dphi[0][0]->Write ("h_jet_trk_dphi_ref");
    h_ljet_trk_dphi[0][0]->Write ("h_ljet_trk_dphi_ref");
    h_sljet_trk_dphi[0][0]->Write ("h_sljet_trk_dphi_ref");
    h_jet_trk_dphi[1][0]->Write ("h_jet_trk_dphi");
    h_ljet_trk_dphi[1][0]->Write ("h_ljet_trk_dphi");
    h_sljet_trk_dphi[1][0]->Write ("h_sljet_trk_dphi");
    h_jet_trk_dphi_bkg[1][0]->Write ("h_jet_trk_dphi_bkg");
    h_ljet_trk_dphi_bkg[1][0]->Write ("h_ljet_trk_dphi_bkg");
    h_sljet_trk_dphi_bkg[1][0]->Write ("h_sljet_trk_dphi_bkg");

    h_jet_trk_dphi_sig[0][0]->Write ("h_jet_trk_dphi_ref_sig");
    h_ljet_trk_dphi_sig[0][0]->Write ("h_ljet_trk_dphi_ref_sig");
    h_sljet_trk_dphi_sig[0][0]->Write ("h_sljet_trk_dphi_ref_sig");
    h_jet_trk_dphi_sig[1][0]->Write ("h_jet_trk_dphi_sig");
    h_ljet_trk_dphi_sig[1][0]->Write ("h_ljet_trk_dphi_sig");
    h_sljet_trk_dphi_sig[1][0]->Write ("h_sljet_trk_dphi_sig");

    h_jet_trk_dphi_iaa[0]->Write ("h_jet_trk_dphi_iaa");
    h_ljet_trk_dphi_iaa[0]->Write ("h_ljet_trk_dphi_iaa");
    h_sljet_trk_dphi_iaa[0]->Write ("h_sljet_trk_dphi_iaa");

    g_jet_trk_dphi_syst[0][0]->Write ("g_jet_trk_dphi_ref_syst");
    g_ljet_trk_dphi_syst[0][0]->Write ("g_ljet_trk_dphi_ref_syst");
    g_sljet_trk_dphi_syst[0][0]->Write ("g_sljet_trk_dphi_ref_syst");
    g_jet_trk_dphi_syst[1][0]->Write ("g_jet_trk_dphi_syst");
    g_ljet_trk_dphi_syst[1][0]->Write ("g_ljet_trk_dphi_syst");
    g_sljet_trk_dphi_syst[1][0]->Write ("g_sljet_trk_dphi_syst");
    g_jet_trk_dphi_bkg_syst[1][0]->Write ("g_jet_trk_dphi_bkg_syst");
    g_ljet_trk_dphi_bkg_syst[1][0]->Write ("g_ljet_trk_dphi_bkg_syst");
    g_sljet_trk_dphi_bkg_syst[1][0]->Write ("g_sljet_trk_dphi_bkg_syst");

    g_jet_trk_dphi_sig_syst[0][0]->Write ("g_jet_trk_dphi_ref_sig_syst");
    g_ljet_trk_dphi_sig_syst[0][0]->Write ("g_ljet_trk_dphi_ref_sig_syst");
    g_sljet_trk_dphi_sig_syst[0][0]->Write ("g_sljet_trk_dphi_ref_sig_syst");
    g_jet_trk_dphi_sig_syst[1][0]->Write ("g_jet_trk_dphi_sig_syst");
    g_ljet_trk_dphi_sig_syst[1][0]->Write ("g_ljet_trk_dphi_sig_syst");
    g_sljet_trk_dphi_sig_syst[1][0]->Write ("g_sljet_trk_dphi_sig_syst");

    g_jet_trk_dphi_iaa_syst[0]->Write ("g_jet_trk_dphi_iaa_syst");
    g_ljet_trk_dphi_iaa_syst[0]->Write ("g_ljet_trk_dphi_iaa_syst");
    g_sljet_trk_dphi_iaa_syst[0]->Write ("g_sljet_trk_dphi_iaa_syst");


    for (int iVar = 1; iVar < nVar; iVar++) {
      h_jet_trk_dphi[0][iVar]->Write (Form ("h_jet_trk_dphi_ref_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi[0][iVar]->Write (Form ("h_ljet_trk_dphi_ref_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi[0][iVar]->Write (Form ("h_sljet_trk_dphi_ref_%s", variations[iVar].Data ()));
      h_jet_trk_dphi[1][iVar]->Write (Form ("h_jet_trk_dphi_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi[1][iVar]->Write (Form ("h_ljet_trk_dphi_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi[1][iVar]->Write (Form ("h_sljet_trk_dphi_%s", variations[iVar].Data ()));
      h_jet_trk_dphi_bkg[1][iVar]->Write (Form ("h_jet_trk_dphi_bkg_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi_bkg[1][iVar]->Write (Form ("h_ljet_trk_dphi_bkg_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi_bkg[1][iVar]->Write (Form ("h_sljet_trk_dphi_bkg_%s", variations[iVar].Data ()));

      h_jet_trk_dphi_sig[0][iVar]->Write (Form ("h_jet_trk_dphi_ref_sig_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi_sig[0][iVar]->Write (Form ("h_ljet_trk_dphi_ref_sig_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi_sig[0][iVar]->Write (Form ("h_sljet_trk_dphi_ref_sig_%s", variations[iVar].Data ()));
      h_jet_trk_dphi_sig[1][iVar]->Write (Form ("h_jet_trk_dphi_sig_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi_sig[1][iVar]->Write (Form ("h_ljet_trk_dphi_sig_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi_sig[1][iVar]->Write (Form ("h_sljet_trk_dphi_sig_%s", variations[iVar].Data ()));

      h_jet_trk_dphi_iaa[iVar]->Write (Form ("h_jet_trk_dphi_iaa_%s", variations[iVar].Data ()));
      h_ljet_trk_dphi_iaa[iVar]->Write (Form ("h_ljet_trk_dphi_iaa_%s", variations[iVar].Data ()));
      h_sljet_trk_dphi_iaa[iVar]->Write (Form ("h_sljet_trk_dphi_iaa_%s", variations[iVar].Data ()));
    }

    outFile->Close ();
  }
}


#endif
