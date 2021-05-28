#ifndef __JetHadronCorrelatorProcessFCalvsZDCPtCh_C__
#define __JetHadronCorrelatorProcessFCalvsZDCPtCh_C__

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

void ProcessFCalvsZDCPtCh (const char* tag, const char* outFileTag, const int iCent) {

  TFile* inFile = nullptr;

  TH1D*** h_evt_counts = Get2DArray <TH1D*> (2, 2);
  TH1D*** h_jet_counts = Get2DArray <TH1D*> (2, 2);
  //TH1D*** h_ljet_counts = Get2DArray <TH1D*> (2, 2);
  //TH1D*** h_sljet_counts = Get2DArray <TH1D*> (2, 2);
  TH1D*** h_evt_counts_bkg = Get2DArray <TH1D*> (2, 2);
  TH1D*** h_jet_counts_bkg = Get2DArray <TH1D*> (2, 2);
  //TH1D*** h_ljet_counts_bkg = Get2DArray <TH1D*> (2, 2);
  //TH1D*** h_sljet_counts_bkg = Get2DArray <TH1D*> (2, 2);

  TH1D*** h_jet_trk_pt_ns = Get2DArray <TH1D*> (2, 2);
  TH2D*** h2_jet_trk_pt_ns_cov = Get2DArray <TH2D*> (2, 2);
  //TH1D*** h_ljet_trk_pt_ns = Get2DArray <TH1D*> (2, 2);
  //TH2D*** h2_ljet_trk_pt_ns_cov = Get2DArray <TH2D*> (2, 2);
  //TH1D*** h_sljet_trk_pt_ns = Get2DArray <TH1D*> (2, 2);
  //TH2D*** h2_sljet_trk_pt_ns_cov = Get2DArray <TH2D*> (2, 2);
  TH1D*** h_jet_trk_pt_as = Get2DArray <TH1D*> (2, 2);
  TH2D*** h2_jet_trk_pt_as_cov = Get2DArray <TH2D*> (2, 2);
  //TH1D*** h_ljet_trk_pt_as = Get2DArray <TH1D*> (2, 2);
  //TH2D*** h2_ljet_trk_pt_as_cov = Get2DArray <TH2D*> (2, 2);
  //TH1D*** h_sljet_trk_pt_as = Get2DArray <TH1D*> (2, 2);
  //TH2D*** h2_sljet_trk_pt_as_cov = Get2DArray <TH2D*> (2, 2);
  TH1D*** h_jet_trk_pt_ns_bkg = Get2DArray <TH1D*> (2, 2);
  TH2D*** h2_jet_trk_pt_ns_cov_bkg = Get2DArray <TH2D*> (2, 2);
  //TH1D*** h_ljet_trk_pt_ns_bkg = Get2DArray <TH1D*> (2, 2);
  //TH2D*** h2_ljet_trk_pt_ns_cov_bkg = Get2DArray <TH2D*> (2, 2);
  //TH1D*** h_sljet_trk_pt_ns_bkg = Get2DArray <TH1D*> (2, 2);
  //TH2D*** h2_sljet_trk_pt_ns_cov_bkg = Get2DArray <TH2D*> (2, 2);
  TH1D*** h_jet_trk_pt_as_bkg = Get2DArray <TH1D*> (2, 2);
  TH2D*** h2_jet_trk_pt_as_cov_bkg = Get2DArray <TH2D*> (2, 2);
  //TH1D*** h_ljet_trk_pt_as_bkg = Get2DArray <TH1D*> (2, 2);
  //TH2D*** h2_ljet_trk_pt_as_cov_bkg = Get2DArray <TH2D*> (2, 2);
  //TH1D*** h_sljet_trk_pt_as_bkg = Get2DArray <TH1D*> (2, 2);
  //TH2D*** h2_sljet_trk_pt_as_cov_bkg = Get2DArray <TH2D*> (2, 2);

  TH1D*** h_jet_trk_pt_ns_sig = Get2DArray <TH1D*> (2, 2);
  //TH1D*** h_ljet_trk_pt_ns_sig = Get2DArray <TH1D*> (2, 2);
  //TH1D*** h_sljet_trk_pt_ns_sig = Get2DArray <TH1D*> (2, 2);
  TH1D*** h_jet_trk_pt_as_sig = Get2DArray <TH1D*> (2, 2);
  //TH1D*** h_ljet_trk_pt_as_sig = Get2DArray <TH1D*> (2, 2);
  //TH1D*** h_sljet_trk_pt_as_sig = Get2DArray <TH1D*> (2, 2);

  TH1D** h_jet_trk_pt_ns_iaa = Get1DArray <TH1D*> (2);
  //TH1D** h_ljet_trk_pt_ns_iaa = Get1DArray <TH1D*> (2);
  //TH1D** h_sljet_trk_pt_ns_iaa = Get1DArray <TH1D*> (2);
  TH1D** h_jet_trk_pt_as_iaa = Get1DArray <TH1D*> (2);
  //TH1D** h_ljet_trk_pt_as_iaa = Get1DArray <TH1D*> (2);
  //TH1D** h_sljet_trk_pt_as_iaa = Get1DArray <TH1D*> (2);


  TGAE*** g_jet_trk_pt_ns_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_ljet_trk_pt_ns_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_sljet_trk_pt_ns_syst = Get2DArray <TGAE*> (2, 2);
  TGAE*** g_jet_trk_pt_as_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_ljet_trk_pt_as_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_sljet_trk_pt_as_syst = Get2DArray <TGAE*> (2, 2);
  TGAE*** g_jet_trk_pt_ns_bkg_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_ljet_trk_pt_ns_bkg_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_sljet_trk_pt_ns_bkg_syst = Get2DArray <TGAE*> (2, 2);
  TGAE*** g_jet_trk_pt_as_bkg_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_ljet_trk_pt_as_bkg_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_sljet_trk_pt_as_bkg_syst = Get2DArray <TGAE*> (2, 2);

  TGAE*** g_jet_trk_pt_ns_sig_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_ljet_trk_pt_ns_sig_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_sljet_trk_pt_ns_sig_syst = Get2DArray <TGAE*> (2, 2);
  TGAE*** g_jet_trk_pt_as_sig_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_ljet_trk_pt_as_sig_syst = Get2DArray <TGAE*> (2, 2);
  //TGAE*** g_sljet_trk_pt_as_sig_syst = Get2DArray <TGAE*> (2, 2);

  TGAE** g_jet_trk_pt_ns_iaa_syst = Get1DArray <TGAE*> (2);
  //TGAE** g_ljet_trk_pt_ns_iaa_syst = Get1DArray <TGAE*> (2);
  //TGAE** g_sljet_trk_pt_ns_iaa_syst = Get1DArray <TGAE*> (2);
  TGAE** g_jet_trk_pt_as_iaa_syst = Get1DArray <TGAE*> (2);
  //TGAE** g_ljet_trk_pt_as_iaa_syst = Get1DArray <TGAE*> (2);
  //TGAE** g_sljet_trk_pt_as_iaa_syst = Get1DArray <TGAE*> (2);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/PlotPtCh_%s_iCent%i.root", rootPath.Data (), outFileName.Data (), iCent);
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  {
    TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), tag, variations[0].Data (), iCent);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");
    outFile->cd ();

    h_evt_counts[1][0] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_%s_pPb_Nominal", tag));
    h_jet_counts[1][0] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", tag))->Clone (Form ("h_jet_counts_%s_pPb_Nominal", tag));
    //h_ljet_counts[1][0] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_%s_pPb_Nominal", tag));
    //h_sljet_counts[1][0] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_data16", tag))->Clone (Form ("h_sljet_counts_%s_pPb_%s", tag));
    h_jet_trk_pt_ns[1][0] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data16", tag))->Clone (Form ("h_jet_trk_pt_ns_%s_pPb_%s", tag));
    h2_jet_trk_pt_ns_cov[1][0] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data16", tag))->Clone (Form ("h2_jet_trk_pt_ns_cov_%s_pPb_%s", tag));
    //h_ljet_trk_pt_ns[1][0] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_%s_data16", tag))->Clone (Form ("h_ljet_trk_pt_ns_%s_pPb_%s", tag));
    //h2_ljet_trk_pt_ns_cov[1][0] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_pt_ns_cov_%s_data16", tag))->Clone (Form ("h2_ljet_trk_pt_ns_cov_%s_pPb_%s", tag));
    //h_sljet_trk_pt_ns[1][0] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_%s_data16", tag))->Clone (Form ("h_sljet_trk_pt_ns_%s_pPb_%s", tag));
    //h2_sljet_trk_pt_ns_cov[1][0] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_pt_ns_cov_%s_data16", tag))->Clone (Form ("h2_sljet_trk_pt_ns_cov_%s_pPb_%s", tag));
    h_jet_trk_pt_as[1][0] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data16", tag))->Clone (Form ("h_jet_trk_pt_as_%s_pPb_%s", tag));
    h2_jet_trk_pt_as_cov[1][0] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data16", tag))->Clone (Form ("h2_jet_trk_pt_as_cov_%s_pPb_%s", tag));
    //h_ljet_trk_pt_as[1][0] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_%s_data16", tag))->Clone (Form ("h_ljet_trk_pt_as_%s_pPb_%s", tag));
    //h2_ljet_trk_pt_as_cov[1][0] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_pt_as_cov_%s_data16", tag))->Clone (Form ("h2_ljet_trk_pt_as_cov_%s_pPb_%s", tag));
    //h_sljet_trk_pt_as[1][0] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_%s_data16", tag))->Clone (Form ("h_sljet_trk_pt_as_%s_pPb_%s", tag));
    //h2_sljet_trk_pt_as_cov[1][0] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_pt_as_cov_%s_data16", tag))->Clone (Form ("h2_sljet_trk_pt_as_cov_%s_pPb_%s", tag));

    inFile->Close ();

    CalcUncertainties (h_jet_trk_pt_ns[1][0], h2_jet_trk_pt_ns_cov[1][0], h_jet_counts[1][0]);
    //CalcUncertainties (h_ljet_trk_pt_ns[1][0], h2_ljet_trk_pt_ns_cov[1][0], h_ljet_counts[1][0]);
    //CalcUncertainties (h_sljet_trk_pt_ns[1][0], h2_sljet_trk_pt_ns_cov[1][0], h_sljet_counts[1][0]);
    CalcUncertainties (h_jet_trk_pt_as[1][0], h2_jet_trk_pt_as_cov[1][0], h_jet_counts[1][0]);
    //CalcUncertainties (h_ljet_trk_pt_as[1][0], h2_ljet_trk_pt_as_cov[1][0], h_ljet_counts[1][0]);
    //CalcUncertainties (h_sljet_trk_pt_as[1][0], h2_sljet_trk_pt_as_cov[1][0], h_sljet_counts[1][0]);
  }



  {
    TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data17_5TeV_hists.root", rootPath.Data (), tag);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");
    outFile->cd ();

    h_evt_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_%s_pp_%s", tag));
    h_jet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_%s_pp_%s", tag));
    //h_ljet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_%s_pp_%s", tag));
    //h_sljet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_data17", tag))->Clone (Form ("h_sljet_counts_%s_pp_%s", tag));
    h_jet_trk_pt_ns[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data17", tag))->Clone (Form ("h_jet_trk_pt_ns_%s_pp_%s", tag));
    h2_jet_trk_pt_ns_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data17", tag))->Clone (Form ("h2_jet_trk_pt_ns_cov_%s_pp_%s", tag));
    //h_ljet_trk_pt_ns[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_%s_data17", tag))->Clone (Form ("h_ljet_trk_pt_ns_%s_pp_%s", tag));
    //h2_ljet_trk_pt_ns_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_pt_ns_cov_%s_data17", tag))->Clone (Form ("h2_ljet_trk_pt_ns_cov_%s_pp_%s", tag));
    //h_sljet_trk_pt_ns[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_%s_data17", tag))->Clone (Form ("h_sljet_trk_pt_ns_%s_pp_%s", tag));
    //h2_sljet_trk_pt_ns_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_pt_ns_cov_%s_data17", tag))->Clone (Form ("h2_sljet_trk_pt_ns_cov_%s_pp_%s", tag));
    h_jet_trk_pt_as[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data17", tag))->Clone (Form ("h_jet_trk_pt_as_%s_pp_%s", tag));
    h2_jet_trk_pt_as_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data17", tag))->Clone (Form ("h2_jet_trk_pt_as_cov_%s_pp_%s", tag));
    //h_ljet_trk_pt_as[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_%s_data17", tag))->Clone (Form ("h_ljet_trk_pt_as_%s_pp_%s", tag));
    //h2_ljet_trk_pt_as_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_pt_as_cov_%s_data17", tag))->Clone (Form ("h2_ljet_trk_pt_as_cov_%s_pp_%s", tag));
    //h_sljet_trk_pt_as[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_%s_data17", tag))->Clone (Form ("h_sljet_trk_pt_as_%s_pp_%s", tag));
    //h2_sljet_trk_pt_as_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_pt_as_cov_%s_data17", tag))->Clone (Form ("h2_sljet_trk_pt_as_cov_%s_pp_%s", tag));

    inFile->Close ();

    CalcUncertainties (h_jet_trk_pt_ns[0][iVar], h2_jet_trk_pt_ns_cov[0][iVar], h_jet_counts[0][iVar]);
    //CalcUncertainties (h_ljet_trk_pt_ns[0][iVar], h2_ljet_trk_pt_ns_cov[0][iVar], h_ljet_counts[0][iVar]);
    //CalcUncertainties (h_sljet_trk_pt_ns[0][iVar], h2_sljet_trk_pt_ns_cov[0][iVar], h_sljet_counts[0][iVar]);
    CalcUncertainties (h_jet_trk_pt_as[0][iVar], h2_jet_trk_pt_as_cov[0][iVar], h_jet_counts[0][iVar]);
    //CalcUncertainties (h_ljet_trk_pt_as[0][iVar], h2_ljet_trk_pt_as_cov[0][iVar], h_ljet_counts[0][iVar]);
    //CalcUncertainties (h_sljet_trk_pt_as[0][iVar], h2_sljet_trk_pt_as_cov[0][iVar], h_sljet_counts[0][iVar]);
  }



  {
    TString inFileName = Form ("%s/Histograms/%s/MixedHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), tag, iCent);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");
    outFile->cd ();

    h_evt_counts_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_mixed_data16", tag))->Clone (Form ("h_evt_counts_%s_mixed_pPb_%s", tag));
    h_jet_counts_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_data16", tag))->Clone (Form ("h_evt_counts_%s_mixed_pPb_%s", tag));
    //h_ljet_counts_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_mixed_data16", tag))->Clone (Form ("h_evt_counts_%s_mixed_pPb_%s", tag));
    //h_sljet_counts_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_mixed_data16", tag))->Clone (Form ("h_sljet_counts_%s_mixed_pPb_%s", tag));
    h_jet_trk_pt_ns_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_mixed_data16", tag))->Clone (Form ("h_jet_trk_pt_ns_%s_mixed_pPb_%s", tag));
    h2_jet_trk_pt_ns_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_mixed_data16", tag))->Clone (Form ("h2_jet_trk_pt_ns_cov_%s_mixed_pPb_%s", tag));
    //h_ljet_trk_pt_ns_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_ns_%s_mixed_data16", tag))->Clone (Form ("h_ljet_trk_pt_ns_%s_mixed_pPb_%s", tag));
    //h2_ljet_trk_pt_ns_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_pt_ns_cov_%s_mixed_data16", tag))->Clone (Form ("h2_ljet_trk_pt_ns_cov_%s_mixed_pPb_%s", tag));
    //h_sljet_trk_pt_ns_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_ns_%s_mixed_data16", tag))->Clone (Form ("h_sljet_trk_pt_ns_%s_mixed_pPb_%s", tag));
    //h2_sljet_trk_pt_ns_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_pt_ns_cov_%s_mixed_data16", tag))->Clone (Form ("h2_sljet_trk_pt_ns_cov_%s_mixed_pPb_%s", tag));
    h_jet_trk_pt_as_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_mixed_data16", tag))->Clone (Form ("h_jet_trk_pt_as_%s_mixed_pPb_%s", tag));
    h2_jet_trk_pt_as_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_mixed_data16", tag))->Clone (Form ("h2_jet_trk_pt_as_cov_%s_mixed_pPb_%s", tag));
    //h_ljet_trk_pt_as_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_trk_pt_as_%s_mixed_data16", tag))->Clone (Form ("h_ljet_trk_pt_as_%s_mixed_pPb_%s", tag));
    //h2_ljet_trk_pt_as_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_ljet_trk_pt_as_cov_%s_mixed_data16", tag))->Clone (Form ("h2_ljet_trk_pt_as_cov_%s_mixed_pPb_%s", tag));
    //h_sljet_trk_pt_as_bkg[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_trk_pt_as_%s_mixed_data16", tag))->Clone (Form ("h_sljet_trk_pt_as_%s_mixed_pPb_%s", tag));
    //h2_sljet_trk_pt_as_cov_bkg[1][iVar] = (TH2D*) inFile->Get (Form ("h2_sljet_trk_pt_as_cov_%s_mixed_data16", tag))->Clone (Form ("h2_sljet_trk_pt_as_cov_%s_mixed_pPb_%s", tag));

    inFile->Close ();

    CalcUncertainties (h_jet_trk_pt_ns_bkg[1][iVar], h2_jet_trk_pt_ns_cov_bkg[1][iVar], h_jet_counts_bkg[1][iVar]);
    //CalcUncertainties (h_ljet_trk_pt_ns_bkg[1][iVar], h2_ljet_trk_pt_ns_cov_bkg[1][iVar], h_ljet_counts_bkg[1][iVar]);
    //CalcUncertainties (h_sljet_trk_pt_ns_bkg[1][iVar], h2_sljet_trk_pt_ns_cov_bkg[1][iVar], h_sljet_counts_bkg[1][iVar]);
    CalcUncertainties (h_jet_trk_pt_as_bkg[1][iVar], h2_jet_trk_pt_as_cov_bkg[1][iVar], h_jet_counts_bkg[1][iVar]);
    //CalcUncertainties (h_ljet_trk_pt_as_bkg[1][iVar], h2_ljet_trk_pt_as_cov_bkg[1][iVar], h_ljet_counts_bkg[1][iVar]);
    //CalcUncertainties (h_sljet_trk_pt_as_bkg[1][iVar], h2_sljet_trk_pt_as_cov_bkg[1][iVar], h_sljet_counts_bkg[1][iVar]);
  }



  {
    h_jet_trk_pt_ns_sig[0][iVar] = (TH1D*) h_jet_trk_pt_ns[0][iVar]->Clone (Form ("h_jet_trk_pt_ns_pp_sig_%s"));
    h_jet_trk_pt_ns_sig[1][iVar] = (TH1D*) h_jet_trk_pt_ns[1][iVar]->Clone (Form ("h_jet_trk_pt_ns_pPb_sig_%s"));
    h_jet_trk_pt_ns_sig[1][iVar]->Add (h_jet_trk_pt_ns_bkg[1][iVar], -1);

    //h_ljet_trk_pt_ns_sig[0][iVar] = (TH1D*) h_ljet_trk_pt_ns[0][iVar]->Clone (Form ("h_ljet_trk_pt_ns_pp_sig_%s"));
    //h_ljet_trk_pt_ns_sig[1][iVar] = (TH1D*) h_ljet_trk_pt_ns[1][iVar]->Clone (Form ("h_ljet_trk_pt_ns_pPb_sig_%s"));
    //h_ljet_trk_pt_ns_sig[1][iVar]->Add (h_ljet_trk_pt_ns_bkg[1][iVar], -1);

    //h_sljet_trk_pt_ns_sig[0][iVar] = (TH1D*) h_sljet_trk_pt_ns[0][iVar]->Clone (Form ("h_sljet_trk_pt_ns_pp_sig_%s"));
    //h_sljet_trk_pt_ns_sig[1][iVar] = (TH1D*) h_sljet_trk_pt_ns[1][iVar]->Clone (Form ("h_sljet_trk_pt_ns_pPb_sig_%s"));
    //h_sljet_trk_pt_ns_sig[1][iVar]->Add (h_sljet_trk_pt_ns_bkg[1][iVar], -1);

    h_jet_trk_pt_as_sig[0][iVar] = (TH1D*) h_jet_trk_pt_as[0][iVar]->Clone (Form ("h_jet_trk_pt_as_pp_sig_%s"));
    h_jet_trk_pt_as_sig[1][iVar] = (TH1D*) h_jet_trk_pt_as[1][iVar]->Clone (Form ("h_jet_trk_pt_as_pPb_sig_%s"));
    h_jet_trk_pt_as_sig[1][iVar]->Add (h_jet_trk_pt_as_bkg[1][iVar], -1);

    //h_ljet_trk_pt_as_sig[0][iVar] = (TH1D*) h_ljet_trk_pt_as[0][iVar]->Clone (Form ("h_ljet_trk_pt_as_pp_sig_%s"));
    //h_ljet_trk_pt_as_sig[1][iVar] = (TH1D*) h_ljet_trk_pt_as[1][iVar]->Clone (Form ("h_ljet_trk_pt_as_pPb_sig_%s"));
    //h_ljet_trk_pt_as_sig[1][iVar]->Add (h_ljet_trk_pt_as_bkg[1][iVar], -1);

    //h_sljet_trk_pt_as_sig[0][iVar] = (TH1D*) h_sljet_trk_pt_as[0][iVar]->Clone (Form ("h_sljet_trk_pt_as_pp_sig_%s"));
    //h_sljet_trk_pt_as_sig[1][iVar] = (TH1D*) h_sljet_trk_pt_as[1][iVar]->Clone (Form ("h_sljet_trk_pt_as_pPb_sig_%s"));
    //h_sljet_trk_pt_as_sig[1][iVar]->Add (h_sljet_trk_pt_as_bkg[1][iVar], -1);


    h_jet_trk_pt_ns_iaa[iVar] = (TH1D*) h_jet_trk_pt_ns_sig[1][iVar]->Clone (Form ("h_jet_trk_pt_ns_iaa_%s"));
    h_jet_trk_pt_ns_iaa[iVar]->Divide (h_jet_trk_pt_ns_sig[0][iVar]);

    //h_ljet_trk_pt_ns_iaa[iVar] = (TH1D*) h_ljet_trk_pt_ns_sig[1][iVar]->Clone (Form ("h_ljet_trk_pt_ns_iaa_%s"));
    //h_ljet_trk_pt_ns_iaa[iVar]->Divide (h_ljet_trk_pt_ns_sig[0][iVar]);

    //h_sljet_trk_pt_ns_iaa[iVar] = (TH1D*) h_sljet_trk_pt_ns_sig[1][iVar]->Clone (Form ("h_sljet_trk_pt_ns_iaa_%s"));
    //h_sljet_trk_pt_ns_iaa[iVar]->Divide (h_sljet_trk_pt_ns_sig[0][iVar]);

    h_jet_trk_pt_as_iaa[iVar] = (TH1D*) h_jet_trk_pt_as_sig[1][iVar]->Clone (Form ("h_jet_trk_pt_as_iaa_%s"));
    h_jet_trk_pt_as_iaa[iVar]->Divide (h_jet_trk_pt_as_sig[0][iVar]);

    //h_ljet_trk_pt_as_iaa[iVar] = (TH1D*) h_ljet_trk_pt_as_sig[1][iVar]->Clone (Form ("h_ljet_trk_pt_as_iaa_%s"));
    //h_ljet_trk_pt_as_iaa[iVar]->Divide (h_ljet_trk_pt_as_sig[0][iVar]);

    //h_sljet_trk_pt_as_iaa[iVar] = (TH1D*) h_sljet_trk_pt_as_sig[1][iVar]->Clone (Form ("h_sljet_trk_pt_as_iaa_%s"));
    //h_sljet_trk_pt_as_iaa[iVar]->Divide (h_sljet_trk_pt_as_sig[0][iVar]);
  }



  {
    g_jet_trk_pt_ns_syst[0][0] = make_graph (h_jet_trk_pt_ns[0][0]);
    //g_ljet_trk_pt_ns_syst[0][0] = make_graph (h_ljet_trk_pt_ns[0][0]);
    //g_sljet_trk_pt_ns_syst[0][0] = make_graph (h_sljet_trk_pt_ns[0][0]);
    g_jet_trk_pt_as_syst[0][0] = make_graph (h_jet_trk_pt_as[0][0]);
    //g_ljet_trk_pt_as_syst[0][0] = make_graph (h_ljet_trk_pt_as[0][0]);
    //g_sljet_trk_pt_as_syst[0][0] = make_graph (h_sljet_trk_pt_as[0][0]);
    g_jet_trk_pt_ns_syst[1][0] = make_graph (h_jet_trk_pt_ns[1][0]);
    //g_ljet_trk_pt_ns_syst[1][0] = make_graph (h_ljet_trk_pt_ns[1][0]);
    //g_sljet_trk_pt_ns_syst[1][0] = make_graph (h_sljet_trk_pt_ns[1][0]);
    g_jet_trk_pt_as_syst[1][0] = make_graph (h_jet_trk_pt_as[1][0]);
    //g_ljet_trk_pt_as_syst[1][0] = make_graph (h_ljet_trk_pt_as[1][0]);
    //g_sljet_trk_pt_as_syst[1][0] = make_graph (h_sljet_trk_pt_as[1][0]);
    g_jet_trk_pt_ns_bkg_syst[1][0] = make_graph (h_jet_trk_pt_ns_bkg[1][0]);
    //g_ljet_trk_pt_ns_bkg_syst[1][0] = make_graph (h_ljet_trk_pt_ns_bkg[1][0]);
    //g_sljet_trk_pt_ns_bkg_syst[1][0] = make_graph (h_sljet_trk_pt_ns_bkg[1][0]);
    g_jet_trk_pt_as_bkg_syst[1][0] = make_graph (h_jet_trk_pt_as_bkg[1][0]);
    //g_ljet_trk_pt_as_bkg_syst[1][0] = make_graph (h_ljet_trk_pt_as_bkg[1][0]);
    //g_sljet_trk_pt_as_bkg_syst[1][0] = make_graph (h_sljet_trk_pt_as_bkg[1][0]);

    g_jet_trk_pt_ns_sig_syst[0][0] = make_graph (h_jet_trk_pt_ns_sig[0][0]);
    //g_ljet_trk_pt_ns_sig_syst[0][0] = make_graph (h_ljet_trk_pt_ns_sig[0][0]);
    //g_sljet_trk_pt_ns_sig_syst[0][0] = make_graph (h_sljet_trk_pt_ns_sig[0][0]);
    g_jet_trk_pt_as_sig_syst[0][0] = make_graph (h_jet_trk_pt_as_sig[0][0]);
    //g_ljet_trk_pt_as_sig_syst[0][0] = make_graph (h_ljet_trk_pt_as_sig[0][0]);
    //g_sljet_trk_pt_as_sig_syst[0][0] = make_graph (h_sljet_trk_pt_as_sig[0][0]);
    g_jet_trk_pt_ns_sig_syst[1][0] = make_graph (h_jet_trk_pt_ns_sig[1][0]);
    //g_ljet_trk_pt_ns_sig_syst[1][0] = make_graph (h_ljet_trk_pt_ns_sig[1][0]);
    //g_sljet_trk_pt_ns_sig_syst[1][0] = make_graph (h_sljet_trk_pt_ns_sig[1][0]);
    g_jet_trk_pt_as_sig_syst[1][0] = make_graph (h_jet_trk_pt_as_sig[1][0]);
    //g_ljet_trk_pt_as_sig_syst[1][0] = make_graph (h_ljet_trk_pt_as_sig[1][0]);
    //g_sljet_trk_pt_as_sig_syst[1][0] = make_graph (h_sljet_trk_pt_as_sig[1][0]);

    g_jet_trk_pt_ns_iaa_syst[0] = make_graph (h_jet_trk_pt_ns_iaa[0]);
    //g_ljet_trk_pt_ns_iaa_syst[0] = make_graph (h_ljet_trk_pt_ns_iaa[0]);
    //g_sljet_trk_pt_ns_iaa_syst[0] = make_graph (h_sljet_trk_pt_ns_iaa[0]);
    g_jet_trk_pt_as_iaa_syst[0] = make_graph (h_jet_trk_pt_as_iaa[0]);
    //g_ljet_trk_pt_as_iaa_syst[0] = make_graph (h_ljet_trk_pt_as_iaa[0]);
    //g_sljet_trk_pt_as_iaa_syst[0] = make_graph (h_sljet_trk_pt_as_iaa[0]);

    ResetTGAEErrors (g_jet_trk_pt_ns_syst[0][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_ns_syst[0][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_ns_syst[0][0]);
    ResetTGAEErrors (g_jet_trk_pt_as_syst[0][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_as_syst[0][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_as_syst[0][0]);
    ResetTGAEErrors (g_jet_trk_pt_ns_syst[1][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_ns_syst[1][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_ns_syst[1][0]);
    ResetTGAEErrors (g_jet_trk_pt_as_syst[1][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_as_syst[1][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_as_syst[1][0]);
    ResetTGAEErrors (g_jet_trk_pt_ns_bkg_syst[1][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_ns_bkg_syst[1][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_ns_bkg_syst[1][0]);
    ResetTGAEErrors (g_jet_trk_pt_as_bkg_syst[1][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_as_bkg_syst[1][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_as_bkg_syst[1][0]);

    ResetTGAEErrors (g_jet_trk_pt_ns_sig_syst[0][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_ns_sig_syst[0][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_ns_sig_syst[0][0]);
    ResetTGAEErrors (g_jet_trk_pt_as_sig_syst[0][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_as_sig_syst[0][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_as_sig_syst[0][0]);
    ResetTGAEErrors (g_jet_trk_pt_ns_sig_syst[1][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_ns_sig_syst[1][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_ns_sig_syst[1][0]);
    ResetTGAEErrors (g_jet_trk_pt_as_sig_syst[1][0]);
    //ResetTGAEErrors (g_ljet_trk_pt_as_sig_syst[1][0]);
    //ResetTGAEErrors (g_sljet_trk_pt_as_sig_syst[1][0]);

    ResetTGAEErrors (g_jet_trk_pt_ns_iaa_syst[0]);
    //ResetTGAEErrors (g_ljet_trk_pt_ns_iaa_syst[0]);
    //ResetTGAEErrors (g_sljet_trk_pt_ns_iaa_syst[0]);
    ResetTGAEErrors (g_jet_trk_pt_as_iaa_syst[0]);
    //ResetTGAEErrors (g_ljet_trk_pt_as_iaa_syst[0]);
    //ResetTGAEErrors (g_sljet_trk_pt_as_iaa_syst[0]);
  }



  for (int iVar = 1; iVar < nVar; iVar++) {
    g_jet_trk_pt_ns_syst[0][iVar] = new TGAE ();
    //g_ljet_trk_pt_ns_syst[0][iVar] = new TGAE ();
    //g_sljet_trk_pt_ns_syst[0][iVar] = new TGAE ();
    g_jet_trk_pt_as_syst[0][iVar] = new TGAE ();
    //g_ljet_trk_pt_as_syst[0][iVar] = new TGAE ();
    //g_sljet_trk_pt_as_syst[0][iVar] = new TGAE ();
    g_jet_trk_pt_ns_syst[1][iVar] = new TGAE ();
    //g_ljet_trk_pt_ns_syst[1][iVar] = new TGAE ();
    //g_sljet_trk_pt_ns_syst[1][iVar] = new TGAE ();
    g_jet_trk_pt_as_syst[1][iVar] = new TGAE ();
    //g_ljet_trk_pt_as_syst[1][iVar] = new TGAE ();
    //g_sljet_trk_pt_as_syst[1][iVar] = new TGAE ();
    g_jet_trk_pt_ns_bkg_syst[1][iVar] = new TGAE ();
    //g_ljet_trk_pt_ns_bkg_syst[1][iVar] = new TGAE ();
    //g_sljet_trk_pt_ns_bkg_syst[1][iVar] = new TGAE ();
    g_jet_trk_pt_as_bkg_syst[1][iVar] = new TGAE ();
    //g_ljet_trk_pt_as_bkg_syst[1][iVar] = new TGAE ();
    //g_sljet_trk_pt_as_bkg_syst[1][iVar] = new TGAE ();

    g_jet_trk_pt_ns_sig_syst[0][iVar] = new TGAE ();
    //g_ljet_trk_pt_ns_sig_syst[0][iVar] = new TGAE ();
    //g_sljet_trk_pt_ns_sig_syst[0][iVar] = new TGAE ();
    g_jet_trk_pt_as_sig_syst[0][iVar] = new TGAE ();
    //g_ljet_trk_pt_as_sig_syst[0][iVar] = new TGAE ();
    //g_sljet_trk_pt_as_sig_syst[0][iVar] = new TGAE ();
    g_jet_trk_pt_ns_sig_syst[1][iVar] = new TGAE ();
    //g_ljet_trk_pt_ns_sig_syst[1][iVar] = new TGAE ();
    //g_sljet_trk_pt_ns_sig_syst[1][iVar] = new TGAE ();
    g_jet_trk_pt_as_sig_syst[1][iVar] = new TGAE ();
    //g_ljet_trk_pt_as_sig_syst[1][iVar] = new TGAE ();
    //g_sljet_trk_pt_as_sig_syst[1][iVar] = new TGAE ();

    g_jet_trk_pt_ns_iaa_syst[iVar] = new TGAE ();
    //g_ljet_trk_pt_ns_iaa_syst[iVar] = new TGAE ();
    //g_sljet_trk_pt_ns_iaa_syst[iVar] = new TGAE ();
    g_jet_trk_pt_as_iaa_syst[iVar] = new TGAE ();
    //g_ljet_trk_pt_as_iaa_syst[iVar] = new TGAE ();
    //g_sljet_trk_pt_as_iaa_syst[iVar] = new TGAE ();

    CalcSystematics (g_jet_trk_pt_ns_syst[0][iVar], h_jet_trk_pt_ns[0][0], h_jet_trk_pt_ns[0][iVar]);
    //CalcSystematics (g_ljet_trk_pt_ns_syst[0][iVar], h_ljet_trk_pt_ns[0][0], h_ljet_trk_pt_ns[0][iVar]);
    //CalcSystematics (g_sljet_trk_pt_ns_syst[0][iVar], h_sljet_trk_pt_ns[0][0], h_sljet_trk_pt_ns[0][iVar]);
    CalcSystematics (g_jet_trk_pt_as_syst[0][iVar], h_jet_trk_pt_as[0][0], h_jet_trk_pt_as[0][iVar]);
    //CalcSystematics (g_ljet_trk_pt_as_syst[0][iVar], h_ljet_trk_pt_as[0][0], h_ljet_trk_pt_as[0][iVar]);
    //CalcSystematics (g_sljet_trk_pt_as_syst[0][iVar], h_sljet_trk_pt_as[0][0], h_sljet_trk_pt_as[0][iVar]);
    CalcSystematics (g_jet_trk_pt_ns_syst[1][iVar], h_jet_trk_pt_ns[1][0], h_jet_trk_pt_ns[1][iVar]);
    //CalcSystematics (g_ljet_trk_pt_ns_syst[1][iVar], h_ljet_trk_pt_ns[1][0], h_ljet_trk_pt_ns[1][iVar]);
    //CalcSystematics (g_sljet_trk_pt_ns_syst[1][iVar], h_sljet_trk_pt_ns[1][0], h_sljet_trk_pt_ns[1][iVar]);
    CalcSystematics (g_jet_trk_pt_as_syst[1][iVar], h_jet_trk_pt_as[1][0], h_jet_trk_pt_as[1][iVar]);
    //CalcSystematics (g_ljet_trk_pt_as_syst[1][iVar], h_ljet_trk_pt_as[1][0], h_ljet_trk_pt_as[1][iVar]);
    //CalcSystematics (g_sljet_trk_pt_as_syst[1][iVar], h_sljet_trk_pt_as[1][0], h_sljet_trk_pt_as[1][iVar]);
    CalcSystematics (g_jet_trk_pt_ns_bkg_syst[1][iVar], h_jet_trk_pt_ns_bkg[1][0], h_jet_trk_pt_ns_bkg[1][iVar]);
    //CalcSystematics (g_ljet_trk_pt_ns_bkg_syst[1][iVar], h_ljet_trk_pt_ns_bkg[1][0], h_ljet_trk_pt_ns_bkg[1][iVar]);
    //CalcSystematics (g_sljet_trk_pt_ns_bkg_syst[1][iVar], h_sljet_trk_pt_ns_bkg[1][0], h_sljet_trk_pt_ns_bkg[1][iVar]);
    CalcSystematics (g_jet_trk_pt_as_bkg_syst[1][iVar], h_jet_trk_pt_as_bkg[1][0], h_jet_trk_pt_as_bkg[1][iVar]);
    //CalcSystematics (g_ljet_trk_pt_as_bkg_syst[1][iVar], h_ljet_trk_pt_as_bkg[1][0], h_ljet_trk_pt_as_bkg[1][iVar]);
    //CalcSystematics (g_sljet_trk_pt_as_bkg_syst[1][iVar], h_sljet_trk_pt_as_bkg[1][0], h_sljet_trk_pt_as_bkg[1][iVar]);

    CalcSystematics (g_jet_trk_pt_ns_sig_syst[0][iVar], h_jet_trk_pt_ns_sig[0][0], h_jet_trk_pt_ns_sig[0][iVar]);
    //CalcSystematics (g_ljet_trk_pt_ns_sig_syst[0][iVar], h_ljet_trk_pt_ns_sig[0][0], h_ljet_trk_pt_ns_sig[0][iVar]);
    //CalcSystematics (g_sljet_trk_pt_ns_sig_syst[0][iVar], h_sljet_trk_pt_ns_sig[0][0], h_sljet_trk_pt_ns_sig[0][iVar]);
    CalcSystematics (g_jet_trk_pt_as_sig_syst[0][iVar], h_jet_trk_pt_as_sig[0][0], h_jet_trk_pt_as_sig[0][iVar]);
    //CalcSystematics (g_ljet_trk_pt_as_sig_syst[0][iVar], h_ljet_trk_pt_as_sig[0][0], h_ljet_trk_pt_as_sig[0][iVar]);
    //CalcSystematics (g_sljet_trk_pt_as_sig_syst[0][iVar], h_sljet_trk_pt_as_sig[0][0], h_sljet_trk_pt_as_sig[0][iVar]);
    CalcSystematics (g_jet_trk_pt_ns_sig_syst[1][iVar], h_jet_trk_pt_ns_sig[1][0], h_jet_trk_pt_ns_sig[1][iVar]);
    //CalcSystematics (g_ljet_trk_pt_ns_sig_syst[1][iVar], h_ljet_trk_pt_ns_sig[1][0], h_ljet_trk_pt_ns_sig[1][iVar]);
    //CalcSystematics (g_sljet_trk_pt_ns_sig_syst[1][iVar], h_sljet_trk_pt_ns_sig[1][0], h_sljet_trk_pt_ns_sig[1][iVar]);
    CalcSystematics (g_jet_trk_pt_as_sig_syst[1][iVar], h_jet_trk_pt_as_sig[1][0], h_jet_trk_pt_as_sig[1][iVar]);
    //CalcSystematics (g_ljet_trk_pt_as_sig_syst[1][iVar], h_ljet_trk_pt_as_sig[1][0], h_ljet_trk_pt_as_sig[1][iVar]);
    //CalcSystematics (g_sljet_trk_pt_as_sig_syst[1][iVar], h_sljet_trk_pt_as_sig[1][0], h_sljet_trk_pt_as_sig[1][iVar]);

    CalcSystematics (g_jet_trk_pt_ns_iaa_syst[iVar], h_jet_trk_pt_ns_iaa[0], h_jet_trk_pt_ns_iaa[iVar]);
    //CalcSystematics (g_ljet_trk_pt_ns_iaa_syst[iVar], h_ljet_trk_pt_ns_iaa[0], h_ljet_trk_pt_ns_iaa[iVar]);
    //CalcSystematics (g_sljet_trk_pt_ns_iaa_syst[iVar], h_sljet_trk_pt_ns_iaa[0], h_sljet_trk_pt_ns_iaa[iVar]);
    CalcSystematics (g_jet_trk_pt_as_iaa_syst[iVar], h_jet_trk_pt_as_iaa[0], h_jet_trk_pt_as_iaa[iVar]);
    //CalcSystematics (g_ljet_trk_pt_as_iaa_syst[iVar], h_ljet_trk_pt_as_iaa[0], h_ljet_trk_pt_as_iaa[iVar]);
    //CalcSystematics (g_sljet_trk_pt_as_iaa_syst[iVar], h_sljet_trk_pt_as_iaa[0], h_sljet_trk_pt_as_iaa[iVar]);
  }



  // takes the maximum variation of the jet pT ES up/down variations
  {
    const int syst = 1;

    AddMaxSystematic (g_jet_trk_pt_ns_syst[0][0], g_jet_trk_pt_ns_syst[0][syst], g_jet_trk_pt_ns_syst[0][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_ns_syst[0][0], g_ljet_trk_pt_ns_syst[0][syst], g_ljet_trk_pt_ns_syst[0][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_ns_syst[0][0], g_sljet_trk_pt_ns_syst[0][syst], g_sljet_trk_pt_ns_syst[0][syst+1]);
    AddMaxSystematic (g_jet_trk_pt_as_syst[0][0], g_jet_trk_pt_as_syst[0][syst], g_jet_trk_pt_as_syst[0][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_as_syst[0][0], g_ljet_trk_pt_as_syst[0][syst], g_ljet_trk_pt_as_syst[0][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_as_syst[0][0], g_sljet_trk_pt_as_syst[0][syst], g_sljet_trk_pt_as_syst[0][syst+1]);
    AddMaxSystematic (g_jet_trk_pt_ns_syst[1][0], g_jet_trk_pt_ns_syst[1][syst], g_jet_trk_pt_ns_syst[1][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_ns_syst[1][0], g_ljet_trk_pt_ns_syst[1][syst], g_ljet_trk_pt_ns_syst[1][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_ns_syst[1][0], g_sljet_trk_pt_ns_syst[1][syst], g_sljet_trk_pt_ns_syst[1][syst+1]);
    AddMaxSystematic (g_jet_trk_pt_as_syst[1][0], g_jet_trk_pt_as_syst[1][syst], g_jet_trk_pt_as_syst[1][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_as_syst[1][0], g_ljet_trk_pt_as_syst[1][syst], g_ljet_trk_pt_as_syst[1][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_as_syst[1][0], g_sljet_trk_pt_as_syst[1][syst], g_sljet_trk_pt_as_syst[1][syst+1]);
    AddMaxSystematic (g_jet_trk_pt_ns_bkg_syst[1][0], g_jet_trk_pt_ns_bkg_syst[1][syst], g_jet_trk_pt_ns_bkg_syst[1][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_ns_bkg_syst[1][0], g_ljet_trk_pt_ns_bkg_syst[1][syst], g_ljet_trk_pt_ns_bkg_syst[1][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_ns_bkg_syst[1][0], g_sljet_trk_pt_ns_bkg_syst[1][syst], g_sljet_trk_pt_ns_bkg_syst[1][syst+1]);
    AddMaxSystematic (g_jet_trk_pt_as_bkg_syst[1][0], g_jet_trk_pt_as_bkg_syst[1][syst], g_jet_trk_pt_as_bkg_syst[1][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_as_bkg_syst[1][0], g_ljet_trk_pt_as_bkg_syst[1][syst], g_ljet_trk_pt_as_bkg_syst[1][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_as_bkg_syst[1][0], g_sljet_trk_pt_as_bkg_syst[1][syst], g_sljet_trk_pt_as_bkg_syst[1][syst+1]);

    AddMaxSystematic (g_jet_trk_pt_ns_sig_syst[0][0], g_jet_trk_pt_ns_sig_syst[0][syst], g_jet_trk_pt_ns_sig_syst[0][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_ns_sig_syst[0][0], g_ljet_trk_pt_ns_sig_syst[0][syst], g_ljet_trk_pt_ns_sig_syst[0][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_ns_sig_syst[0][0], g_sljet_trk_pt_ns_sig_syst[0][syst], g_sljet_trk_pt_ns_sig_syst[0][syst+1]);
    AddMaxSystematic (g_jet_trk_pt_as_sig_syst[0][0], g_jet_trk_pt_as_sig_syst[0][syst], g_jet_trk_pt_as_sig_syst[0][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_as_sig_syst[0][0], g_ljet_trk_pt_as_sig_syst[0][syst], g_ljet_trk_pt_as_sig_syst[0][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_as_sig_syst[0][0], g_sljet_trk_pt_as_sig_syst[0][syst], g_sljet_trk_pt_as_sig_syst[0][syst+1]);
    AddMaxSystematic (g_jet_trk_pt_ns_sig_syst[1][0], g_jet_trk_pt_ns_sig_syst[1][syst], g_jet_trk_pt_ns_sig_syst[1][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_ns_sig_syst[1][0], g_ljet_trk_pt_ns_sig_syst[1][syst], g_ljet_trk_pt_ns_sig_syst[1][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_ns_sig_syst[1][0], g_sljet_trk_pt_ns_sig_syst[1][syst], g_sljet_trk_pt_ns_sig_syst[1][syst+1]);
    AddMaxSystematic (g_jet_trk_pt_as_sig_syst[1][0], g_jet_trk_pt_as_sig_syst[1][syst], g_jet_trk_pt_as_sig_syst[1][syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_as_sig_syst[1][0], g_ljet_trk_pt_as_sig_syst[1][syst], g_ljet_trk_pt_as_sig_syst[1][syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_as_sig_syst[1][0], g_sljet_trk_pt_as_sig_syst[1][syst], g_sljet_trk_pt_as_sig_syst[1][syst+1]);

    AddMaxSystematic (g_jet_trk_pt_ns_iaa_syst[0], g_jet_trk_pt_ns_iaa_syst[syst], g_jet_trk_pt_ns_iaa_syst[syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_ns_iaa_syst[0], g_ljet_trk_pt_ns_iaa_syst[syst], g_ljet_trk_pt_ns_iaa_syst[syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_ns_iaa_syst[0], g_sljet_trk_pt_ns_iaa_syst[syst], g_sljet_trk_pt_ns_iaa_syst[syst+1]);
    AddMaxSystematic (g_jet_trk_pt_as_iaa_syst[0], g_jet_trk_pt_as_iaa_syst[syst], g_jet_trk_pt_as_iaa_syst[syst+1]);
    //AddMaxSystematic (g_ljet_trk_pt_as_iaa_syst[0], g_ljet_trk_pt_as_iaa_syst[syst], g_ljet_trk_pt_as_iaa_syst[syst+1]);
    //AddMaxSystematic (g_sljet_trk_pt_as_iaa_syst[0], g_sljet_trk_pt_as_iaa_syst[syst], g_sljet_trk_pt_as_iaa_syst[syst+1]);
  }



  {
    outFile->cd ();

    h_evt_counts[0][0]->Write ("h_evt_counts_ref");
    h_jet_counts[0][0]->Write ("h_jet_counts_ref");
    //h_ljet_counts[0][0]->Write ("h_ljet_counts_ref");
    //h_sljet_counts[0][0]->Write ("h_sljet_counts_ref");
    h_evt_counts[1][0]->Write ("h_evt_counts");
    h_jet_counts[1][0]->Write ("h_jet_counts");
    //h_ljet_counts[1][0]->Write ("h_ljet_counts");
    //h_sljet_counts[1][0]->Write ("h_sljet_counts");
    h_evt_counts_bkg[1][0]->Write ("h_evt_counts_bkg");
    h_jet_counts_bkg[1][0]->Write ("h_jet_counts_bkg");
    //h_ljet_counts_bkg[1][0]->Write ("h_ljet_counts_bkg");
    //h_sljet_counts_bkg[1][0]->Write ("h_sljet_counts_bkg");

    h_jet_trk_pt_ns[0][0]->Write ("h_jet_trk_pt_ns_ref");
    //h_ljet_trk_pt_ns[0][0]->Write ("h_ljet_trk_pt_ns_ref");
    //h_sljet_trk_pt_ns[0][0]->Write ("h_sljet_trk_pt_ns_ref");
    h_jet_trk_pt_as[0][0]->Write ("h_jet_trk_pt_as_ref");
    //h_ljet_trk_pt_as[0][0]->Write ("h_ljet_trk_pt_as_ref");
    //h_sljet_trk_pt_as[0][0]->Write ("h_sljet_trk_pt_as_ref");
    h_jet_trk_pt_ns[1][0]->Write ("h_jet_trk_pt_ns");
    //h_ljet_trk_pt_ns[1][0]->Write ("h_ljet_trk_pt_ns");
    //h_sljet_trk_pt_ns[1][0]->Write ("h_sljet_trk_pt_ns");
    h_jet_trk_pt_as[1][0]->Write ("h_jet_trk_pt_as");
    //h_ljet_trk_pt_as[1][0]->Write ("h_ljet_trk_pt_as");
    //h_sljet_trk_pt_as[1][0]->Write ("h_sljet_trk_pt_as");
    h_jet_trk_pt_ns_bkg[1][0]->Write ("h_jet_trk_pt_ns_bkg");
    //h_ljet_trk_pt_ns_bkg[1][0]->Write ("h_ljet_trk_pt_ns_bkg");
    //h_sljet_trk_pt_ns_bkg[1][0]->Write ("h_sljet_trk_pt_ns_bkg");
    h_jet_trk_pt_as_bkg[1][0]->Write ("h_jet_trk_pt_as_bkg");
    //h_ljet_trk_pt_as_bkg[1][0]->Write ("h_ljet_trk_pt_as_bkg");
    //h_sljet_trk_pt_as_bkg[1][0]->Write ("h_sljet_trk_pt_as_bkg");

    h_jet_trk_pt_ns_sig[0][0]->Write ("h_jet_trk_pt_ns_ref_sig");
    //h_ljet_trk_pt_ns_sig[0][0]->Write ("h_ljet_trk_pt_ns_ref_sig");
    //h_sljet_trk_pt_ns_sig[0][0]->Write ("h_sljet_trk_pt_ns_ref_sig");
    h_jet_trk_pt_as_sig[0][0]->Write ("h_jet_trk_pt_as_ref_sig");
    //h_ljet_trk_pt_as_sig[0][0]->Write ("h_ljet_trk_pt_as_ref_sig");
    //h_sljet_trk_pt_as_sig[0][0]->Write ("h_sljet_trk_pt_as_ref_sig");
    h_jet_trk_pt_ns_sig[1][0]->Write ("h_jet_trk_pt_ns_sig");
    //h_ljet_trk_pt_ns_sig[1][0]->Write ("h_ljet_trk_pt_ns_sig");
    //h_sljet_trk_pt_ns_sig[1][0]->Write ("h_sljet_trk_pt_ns_sig");
    h_jet_trk_pt_as_sig[1][0]->Write ("h_jet_trk_pt_as_sig");
    //h_ljet_trk_pt_as_sig[1][0]->Write ("h_ljet_trk_pt_as_sig");
    //h_sljet_trk_pt_as_sig[1][0]->Write ("h_sljet_trk_pt_as_sig");

    h_jet_trk_pt_ns_iaa[0]->Write ("h_jet_trk_pt_ns_iaa");
    //h_ljet_trk_pt_ns_iaa[0]->Write ("h_ljet_trk_pt_ns_iaa");
    //h_sljet_trk_pt_ns_iaa[0]->Write ("h_sljet_trk_pt_ns_iaa");
    h_jet_trk_pt_as_iaa[0]->Write ("h_jet_trk_pt_as_iaa");
    //h_ljet_trk_pt_as_iaa[0]->Write ("h_ljet_trk_pt_as_iaa");
    //h_sljet_trk_pt_as_iaa[0]->Write ("h_sljet_trk_pt_as_iaa");

    g_jet_trk_pt_ns_syst[0][0]->Write ("g_jet_trk_pt_ns_ref_syst");
    //g_ljet_trk_pt_ns_syst[0][0]->Write ("g_ljet_trk_pt_ns_ref_syst");
    //g_sljet_trk_pt_ns_syst[0][0]->Write ("g_sljet_trk_pt_ns_ref_syst");
    g_jet_trk_pt_as_syst[0][0]->Write ("g_jet_trk_pt_as_ref_syst");
    //g_ljet_trk_pt_as_syst[0][0]->Write ("g_ljet_trk_pt_as_ref_syst");
    //g_sljet_trk_pt_as_syst[0][0]->Write ("g_sljet_trk_pt_as_ref_syst");
    g_jet_trk_pt_ns_syst[1][0]->Write ("g_jet_trk_pt_ns_syst");
    //g_ljet_trk_pt_ns_syst[1][0]->Write ("g_ljet_trk_pt_ns_syst");
    //g_sljet_trk_pt_ns_syst[1][0]->Write ("g_sljet_trk_pt_ns_syst");
    g_jet_trk_pt_as_syst[1][0]->Write ("g_jet_trk_pt_as_syst");
    //g_ljet_trk_pt_as_syst[1][0]->Write ("g_ljet_trk_pt_as_syst");
    //g_sljet_trk_pt_as_syst[1][0]->Write ("g_sljet_trk_pt_as_syst");
    g_jet_trk_pt_ns_bkg_syst[1][0]->Write ("g_jet_trk_pt_ns_bkg_syst");
    //g_ljet_trk_pt_ns_bkg_syst[1][0]->Write ("g_ljet_trk_pt_ns_bkg_syst");
    //g_sljet_trk_pt_ns_bkg_syst[1][0]->Write ("g_sljet_trk_pt_ns_bkg_syst");
    g_jet_trk_pt_as_bkg_syst[1][0]->Write ("g_jet_trk_pt_as_bkg_syst");
    //g_ljet_trk_pt_as_bkg_syst[1][0]->Write ("g_ljet_trk_pt_as_bkg_syst");
    //g_sljet_trk_pt_as_bkg_syst[1][0]->Write ("g_sljet_trk_pt_as_bkg_syst");

    g_jet_trk_pt_ns_sig_syst[0][0]->Write ("g_jet_trk_pt_ns_ref_sig_syst");
    //g_ljet_trk_pt_ns_sig_syst[0][0]->Write ("g_ljet_trk_pt_ns_ref_sig_syst");
    //g_sljet_trk_pt_ns_sig_syst[0][0]->Write ("g_sljet_trk_pt_ns_ref_sig_syst");
    g_jet_trk_pt_as_sig_syst[0][0]->Write ("g_jet_trk_pt_as_ref_sig_syst");
    //g_ljet_trk_pt_as_sig_syst[0][0]->Write ("g_ljet_trk_pt_as_ref_sig_syst");
    //g_sljet_trk_pt_as_sig_syst[0][0]->Write ("g_sljet_trk_pt_as_ref_sig_syst");
    g_jet_trk_pt_ns_sig_syst[1][0]->Write ("g_jet_trk_pt_ns_sig_syst");
    //g_ljet_trk_pt_ns_sig_syst[1][0]->Write ("g_ljet_trk_pt_ns_sig_syst");
    //g_sljet_trk_pt_ns_sig_syst[1][0]->Write ("g_sljet_trk_pt_ns_sig_syst");
    g_jet_trk_pt_as_sig_syst[1][0]->Write ("g_jet_trk_pt_as_sig_syst");
    //g_ljet_trk_pt_as_sig_syst[1][0]->Write ("g_ljet_trk_pt_as_sig_syst");
    //g_sljet_trk_pt_as_sig_syst[1][0]->Write ("g_sljet_trk_pt_as_sig_syst");

    g_jet_trk_pt_ns_iaa_syst[0]->Write ("g_jet_trk_pt_ns_iaa_syst");
    //g_ljet_trk_pt_ns_iaa_syst[0]->Write ("g_ljet_trk_pt_ns_iaa_syst");
    //g_sljet_trk_pt_ns_iaa_syst[0]->Write ("g_sljet_trk_pt_ns_iaa_syst");
    g_jet_trk_pt_as_iaa_syst[0]->Write ("g_jet_trk_pt_as_iaa_syst");
    //g_ljet_trk_pt_as_iaa_syst[0]->Write ("g_ljet_trk_pt_as_iaa_syst");
    //g_sljet_trk_pt_as_iaa_syst[0]->Write ("g_sljet_trk_pt_as_iaa_syst");


    for (int iVar = 1; iVar < nVar; iVar++) {
      h_jet_trk_pt_ns[0][iVar]->Write (Form ("h_jet_trk_pt_ns_ref_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_ns[0][iVar]->Write (Form ("h_ljet_trk_pt_ns_ref_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_ns[0][iVar]->Write (Form ("h_sljet_trk_pt_ns_ref_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_as[0][iVar]->Write (Form ("h_jet_trk_pt_as_ref_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_as[0][iVar]->Write (Form ("h_ljet_trk_pt_as_ref_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_as[0][iVar]->Write (Form ("h_sljet_trk_pt_as_ref_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_ns[1][iVar]->Write (Form ("h_jet_trk_pt_ns_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_ns[1][iVar]->Write (Form ("h_ljet_trk_pt_ns_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_ns[1][iVar]->Write (Form ("h_sljet_trk_pt_ns_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_as[1][iVar]->Write (Form ("h_jet_trk_pt_as_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_as[1][iVar]->Write (Form ("h_ljet_trk_pt_as_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_as[1][iVar]->Write (Form ("h_sljet_trk_pt_as_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_ns_bkg[1][iVar]->Write (Form ("h_jet_trk_pt_ns_bkg_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_ns_bkg[1][iVar]->Write (Form ("h_ljet_trk_pt_ns_bkg_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_ns_bkg[1][iVar]->Write (Form ("h_sljet_trk_pt_ns_bkg_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_as_bkg[1][iVar]->Write (Form ("h_jet_trk_pt_as_bkg_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_as_bkg[1][iVar]->Write (Form ("h_ljet_trk_pt_as_bkg_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_as_bkg[1][iVar]->Write (Form ("h_sljet_trk_pt_as_bkg_FineFcalCentVar_%i", iCent));

      h_jet_trk_pt_ns_sig[0][iVar]->Write (Form ("h_jet_trk_pt_ns_ref_sig_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_ns_sig[0][iVar]->Write (Form ("h_ljet_trk_pt_ns_ref_sig_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_ns_sig[0][iVar]->Write (Form ("h_sljet_trk_pt_ns_ref_sig_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_as_sig[0][iVar]->Write (Form ("h_jet_trk_pt_as_ref_sig_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_as_sig[0][iVar]->Write (Form ("h_ljet_trk_pt_as_ref_sig_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_as_sig[0][iVar]->Write (Form ("h_sljet_trk_pt_as_ref_sig_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_ns_sig[1][iVar]->Write (Form ("h_jet_trk_pt_ns_sig_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_ns_sig[1][iVar]->Write (Form ("h_ljet_trk_pt_ns_sig_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_ns_sig[1][iVar]->Write (Form ("h_sljet_trk_pt_ns_sig_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_as_sig[1][iVar]->Write (Form ("h_jet_trk_pt_as_sig_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_as_sig[1][iVar]->Write (Form ("h_ljet_trk_pt_as_sig_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_as_sig[1][iVar]->Write (Form ("h_sljet_trk_pt_as_sig_FineFcalCentVar_%i", iCent));

      h_jet_trk_pt_ns_iaa[iVar]->Write (Form ("h_jet_trk_pt_ns_iaa_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_ns_iaa[iVar]->Write (Form ("h_ljet_trk_pt_ns_iaa_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_ns_iaa[iVar]->Write (Form ("h_sljet_trk_pt_ns_iaa_FineFcalCentVar_%i", iCent));
      h_jet_trk_pt_as_iaa[iVar]->Write (Form ("h_jet_trk_pt_as_iaa_FineFcalCentVar_%i", iCent));
      //h_ljet_trk_pt_as_iaa[iVar]->Write (Form ("h_ljet_trk_pt_as_iaa_FineFcalCentVar_%i", iCent));
      //h_sljet_trk_pt_as_iaa[iVar]->Write (Form ("h_sljet_trk_pt_as_iaa_FineFcalCentVar_%i", iCent));
    }

    outFile->Close ();
  }
}


#endif
