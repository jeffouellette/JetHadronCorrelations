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
#include <vector>
#include <math.h>

using namespace JetHadronCorrelations;

void ProcessDPhi (const char* outFileTag, const char* tag1, const char* tag2 = nullptr) {

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

  TH1D***  h_jet_trk_dphi_ref = Get2DArray <TH1D*> (nPtChSelections, nVar);
  TH2D***  h2_jet_trk_dphi_cov_ref = Get2DArray <TH2D*> (nPtChSelections, nVar);
  TH1D***  h_jet_trk_dphi_ref_bkg = Get2DArray <TH1D*> (nPtChSelections, nVar);
  TH2D***  h2_jet_trk_dphi_cov_ref_bkg = Get2DArray <TH2D*> (nPtChSelections, nVar);
  TH1D**** h_jet_trk_dphi= Get3DArray <TH1D*> (nZdcCentBins, nPtChSelections, nVar);
  TH2D**** h2_jet_trk_dphi_cov = Get3DArray <TH2D*> (nZdcCentBins, nPtChSelections, nVar);
  TH1D**** h_jet_trk_dphi_bkg = Get3DArray <TH1D*> (nZdcCentBins, nPtChSelections, nVar);
  TH2D**** h2_jet_trk_dphi_cov_bkg = Get3DArray <TH2D*> (nZdcCentBins, nPtChSelections, nVar);

  TH1D***  h_jet_trk_dphi_ref_sig = Get2DArray <TH1D*> (nPtChSelections, nVar);
  TH1D**** h_jet_trk_dphi_sig = Get3DArray <TH1D*> (nZdcCentBins, nPtChSelections, nVar);

  TH1D**** h_jet_trk_dphi_iaa = Get3DArray <TH1D*> (nZdcCentBins, nPtChSelections, nVar);

  TGAE***  g_jet_trk_dphi_ref_syst = Get2DArray <TGAE*> (nPtChSelections, nVar);
  TGAE***  g_jet_trk_dphi_ref_bkg_syst = Get2DArray <TGAE*> (nPtChSelections, nVar);
  TGAE**** g_jet_trk_dphi_syst = Get3DArray <TGAE*> (nZdcCentBins, nPtChSelections, nVar);
  TGAE**** g_jet_trk_dphi_bkg_syst = Get3DArray <TGAE*> (nZdcCentBins, nPtChSelections, nVar);

  TGAE***  g_jet_trk_dphi_ref_sig_syst = Get2DArray <TGAE*> (nPtChSelections, nVar);
  TGAE**** g_jet_trk_dphi_sig_syst = Get3DArray <TGAE*> (nZdcCentBins, nPtChSelections, nVar);

  TGAE**** g_jet_trk_dphi_iaa_syst = Get3DArray <TGAE*> (nZdcCentBins, nPtChSelections, nVar);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/PlotDPhi_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  // loop over all variations of analysis
  for (int iVar = 0; iVar < nVar; iVar++) {

    // read in pp histograms
    {
      TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data17_5TeV_hists.root", rootPath.Data (), tag1, variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts_ref[iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", tag1))->Clone (Form ("h_evt_counts_pp_%s", variations[iVar].Data ()));
      h_jet_counts_ref[iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", tag1))->Clone (Form ("h_jet_counts_pp_%s", variations[iVar].Data ()));

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
        h_jet_trk_dphi_ref[iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_data17", pTChSelections[iPtCh].Data (), tag1))->Clone (Form ("h_jet_trk_dphi_%s_pp_%s", pTChSelections[iPtCh].Data (), variations[iVar].Data ()));
        h2_jet_trk_dphi_cov_ref[iPtCh][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_%s_cov_%s_data17", pTChSelections[iPtCh].Data (), tag1))->Clone (Form ("h2_jet_trk_dphi_%s_cov_pp_%s", pTChSelections[iPtCh].Data (), variations[iVar].Data ()));
      }

      inFile->Close ();

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++)
        CalcUncertainties (h_jet_trk_dphi_ref[iPtCh][iVar], h2_jet_trk_dphi_cov_ref[iPtCh][iVar], h_jet_counts_ref[iVar]);
    }



    // read in pp histograms
    {
      TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/data17_5TeV_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts_ref_bkg[iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_evt_counts_ref_bkg_%s", variations[iVar].Data ()));
      h_jet_counts_ref_bkg[iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_counts_ref_bkg_%s", variations[iVar].Data ()));

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
        h_jet_trk_dphi_ref_bkg[iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_data17", pTChSelections[iPtCh].Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_dphi_%s_ref_bkg_%s", pTChSelections[iPtCh].Data (), variations[iVar].Data ()));
        h2_jet_trk_dphi_cov_ref_bkg[iPtCh][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_%s_cov_%s_data17", pTChSelections[iPtCh].Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_dphi_%s_cov_ref_bkg_%s", pTChSelections[iPtCh].Data (), variations[iVar].Data ()));
      }

      inFile->Close ();

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++)
        CalcUncertainties (h_jet_trk_dphi_ref_bkg[iPtCh][iVar], h2_jet_trk_dphi_cov_ref_bkg[iPtCh][iVar], h_jet_counts_ref_bkg[iVar]);
    }



    // read in all p+Pb centralities histograms
    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), tag1, variations[iVar].Data (), iCent);
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", tag1))->Clone (Form ("h_evt_counts_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_counts[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", tag1))->Clone (Form ("h_jet_counts_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
        h_jet_trk_dphi[iCent][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_data16", pTChSelections[iPtCh].Data (), tag1))->Clone (Form ("h_jet_trk_dphi_%s_pPb_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
        h2_jet_trk_dphi_cov[iCent][iPtCh][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_%s_cov_%s_data16", pTChSelections[iPtCh].Data (), tag1))->Clone (Form ("h2_jet_trk_dphi_%s_cov_pPb_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
      }

      inFile->Close ();

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++)
        CalcUncertainties (h_jet_trk_dphi[iCent][iPtCh][iVar], h2_jet_trk_dphi_cov[iCent][iPtCh][iVar], h_jet_counts[iCent][iVar]);
    }



    // read in all p+Pb centralities mixed event histograms
    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", variations[iVar].Data (), iCent);
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts_bkg[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_evt_counts_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_counts_bkg[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_counts_pPb_bkg_iCent%i_%s", iCent, variations[iVar].Data ()));

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
        h_jet_trk_dphi_bkg[iCent][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_data16", pTChSelections[iPtCh].Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_dphi_%s_pPb_bkg_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
        h2_jet_trk_dphi_cov_bkg[iCent][iPtCh][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_%s_cov_%s_data16", pTChSelections[iPtCh].Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_dphi_%s_cov_pPb_bkg_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
      }

      inFile->Close ();

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++)
        CalcUncertainties (h_jet_trk_dphi_bkg[iCent][iPtCh][iVar], h2_jet_trk_dphi_cov_bkg[iCent][iPtCh][iVar], h_jet_counts_bkg[iCent][iVar]);
    }



    for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      h_jet_trk_dphi_ref_sig[iPtCh][iVar] = (TH1D*) h_jet_trk_dphi_ref[iPtCh][iVar]->Clone (Form ("h_jet_trk_dphi_%s_ref_sig_%s", pTChSelections[iPtCh].Data (), variations[iVar].Data ()));
      h_jet_trk_dphi_ref_sig[iPtCh][iVar]->Add (h_jet_trk_dphi_ref_bkg[iPtCh][iVar], -1);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_jet_trk_dphi_sig[iCent][iPtCh][iVar] = (TH1D*) h_jet_trk_dphi[iCent][iPtCh][iVar]->Clone (Form ("h_jet_trk_dphi_%s_pPb_sig_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
        h_jet_trk_dphi_sig[iCent][iPtCh][iVar]->Add (h_jet_trk_dphi_bkg[iCent][iPtCh][iVar], -1);

        h_jet_trk_dphi_iaa[iCent][iPtCh][iVar] = (TH1D*) h_jet_trk_dphi_sig[iCent][iPtCh][iVar]->Clone (Form ("h_jet_trk_dphi_%s_iaa_iCent%i_%s", pTChSelections[iPtCh].Data (), iCent, variations[iVar].Data ()));
        h_jet_trk_dphi_iaa[iCent][iPtCh][iVar]->Divide (h_jet_trk_dphi_ref_sig[iPtCh][iVar]);

      }
    }
  }



  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

    g_jet_trk_dphi_ref_syst[iPtCh][0] = make_graph (h_jet_trk_dphi_ref[iPtCh][0]);
    g_jet_trk_dphi_ref_bkg_syst[iPtCh][0] = make_graph (h_jet_trk_dphi_ref_bkg[iPtCh][0]);
    g_jet_trk_dphi_ref_sig_syst[iPtCh][0] = make_graph (h_jet_trk_dphi_ref_sig[iPtCh][0]);

    ResetTGAEErrors (g_jet_trk_dphi_ref_syst[iPtCh][0]);
    ResetTGAEErrors (g_jet_trk_dphi_ref_bkg_syst[iPtCh][0]);
    ResetTGAEErrors (g_jet_trk_dphi_ref_sig_syst[iPtCh][0]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      g_jet_trk_dphi_syst[iCent][iPtCh][0] = make_graph (h_jet_trk_dphi[iCent][iPtCh][0]);
      g_jet_trk_dphi_bkg_syst[iCent][iPtCh][0] = make_graph (h_jet_trk_dphi_bkg[iCent][iPtCh][0]);
      g_jet_trk_dphi_sig_syst[iCent][iPtCh][0] = make_graph (h_jet_trk_dphi_sig[iCent][iPtCh][0]);
      g_jet_trk_dphi_iaa_syst[iCent][iPtCh][0] = make_graph (h_jet_trk_dphi_iaa[iCent][iPtCh][0]);

      ResetTGAEErrors (g_jet_trk_dphi_syst[iCent][iPtCh][0]);
      ResetTGAEErrors (g_jet_trk_dphi_bkg_syst[iCent][iPtCh][0]);
      ResetTGAEErrors (g_jet_trk_dphi_sig_syst[iCent][iPtCh][0]);
      ResetTGAEErrors (g_jet_trk_dphi_iaa_syst[iCent][iPtCh][0]);

    }
  }



  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    for (int iVar = 1; iVar < nVar; iVar++) {

      g_jet_trk_dphi_ref_syst[iPtCh][iVar] = new TGAE ();
      g_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar] = new TGAE ();
      g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar] = new TGAE ();

      CalcSystematics (g_jet_trk_dphi_ref_syst[iPtCh][iVar], h_jet_trk_dphi_ref[iPtCh][0], h_jet_trk_dphi_ref[iPtCh][iVar]);
      CalcSystematics (g_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar], h_jet_trk_dphi_ref_bkg[iPtCh][0], h_jet_trk_dphi_ref_bkg[iPtCh][iVar]);
      CalcSystematics (g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar], h_jet_trk_dphi_ref_sig[iPtCh][0], h_jet_trk_dphi_ref_sig[iPtCh][iVar]);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        g_jet_trk_dphi_syst[iCent][iPtCh][iVar] = new TGAE ();
        g_jet_trk_dphi_bkg_syst[iCent][iPtCh][iVar] = new TGAE ();
        g_jet_trk_dphi_sig_syst[iCent][iPtCh][iVar] = new TGAE ();
        g_jet_trk_dphi_iaa_syst[iCent][iPtCh][iVar] = new TGAE ();

        CalcSystematics (g_jet_trk_dphi_syst[iCent][iPtCh][iVar], h_jet_trk_dphi[iCent][iPtCh][0], h_jet_trk_dphi[iCent][iPtCh][iVar]);
        CalcSystematics (g_jet_trk_dphi_bkg_syst[iCent][iPtCh][iVar], h_jet_trk_dphi_bkg[iCent][iPtCh][0], h_jet_trk_dphi_bkg[iCent][iPtCh][iVar]);
        CalcSystematics (g_jet_trk_dphi_sig_syst[iCent][iPtCh][iVar], h_jet_trk_dphi_sig[iCent][iPtCh][0], h_jet_trk_dphi_sig[iCent][iPtCh][iVar]);
        CalcSystematics (g_jet_trk_dphi_iaa_syst[iCent][iPtCh][iVar], h_jet_trk_dphi_iaa[iCent][iPtCh][0], h_jet_trk_dphi_iaa[iCent][iPtCh][iVar]);
      }
    }
  }



/*
  // takes the maximum variation of the jet pT ES up/down variations
  // TODO
  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

    const int syst = 7; // MixCatVar1

    AddMaxSystematic (g_jet_trk_dphi_ref_syst[iPtCh][0], g_jet_trk_dphi_ref_syst[iPtCh][syst], g_jet_trk_dphi_ref_syst[iPtCh][syst]);
    AddMaxSystematic (g_jet_trk_dphi_ref_bkg_syst[iPtCh][0], g_jet_trk_dphi_ref_bkg_syst[iPtCh][syst], g_jet_trk_dphi_ref_bkg_syst[iPtCh][syst]);
    AddMaxSystematic (g_jet_trk_dphi_ref_sig_syst[iPtCh][0], g_jet_trk_dphi_ref_sig_syst[iPtCh][syst], g_jet_trk_dphi_ref_sig_syst[iPtCh][syst]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      AddMaxSystematic (g_jet_trk_dphi_syst[iCent][iPtCh][0], g_jet_trk_dphi_syst[iCent][iPtCh][syst], g_jet_trk_dphi_syst[iCent][iPtCh][syst]);
      AddMaxSystematic (g_jet_trk_dphi_bkg_syst[iCent][iPtCh][0], g_jet_trk_dphi_bkg_syst[iCent][iPtCh][syst], g_jet_trk_dphi_bkg_syst[iCent][iPtCh][syst]);
      AddMaxSystematic (g_jet_trk_dphi_sig_syst[iCent][iPtCh][0], g_jet_trk_dphi_sig_syst[iCent][iPtCh][syst], g_jet_trk_dphi_sig_syst[iCent][iPtCh][syst]);
      AddMaxSystematic (g_jet_trk_dphi_iaa_syst[iCent][iPtCh][0], g_jet_trk_dphi_iaa_syst[iCent][iPtCh][syst], g_jet_trk_dphi_iaa_syst[iCent][iPtCh][syst]);

    }


    //const int syst = 1; // 5% JES systematic

    //AddMaxSystematic (g_jet_trk_dphi_ref_syst[iPtCh][0], g_jet_trk_dphi_ref_syst[iPtCh][syst], g_jet_trk_dphi_ref_syst[iPtCh][syst+1]);
    //AddMaxSystematic (g_jet_trk_dphi_ref_bkg_syst[iPtCh][0], g_jet_trk_dphi_ref_bkg_syst[iPtCh][syst], g_jet_trk_dphi_ref_bkg_syst[iPtCh][syst+1]);
    //AddMaxSystematic (g_jet_trk_dphi_ref_sig_syst[iPtCh][0], g_jet_trk_dphi_ref_sig_syst[iPtCh][syst], g_jet_trk_dphi_ref_sig_syst[iPtCh][syst+1]);

    //for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

    //  AddMaxSystematic (g_jet_trk_dphi_syst[iCent][iPtCh][0], g_jet_trk_dphi_syst[iCent][iPtCh][syst], g_jet_trk_dphi_syst[iCent][iPtCh][syst+1]);
    //  AddMaxSystematic (g_jet_trk_dphi_bkg_syst[iCent][iPtCh][0], g_jet_trk_dphi_bkg_syst[iCent][iPtCh][syst], g_jet_trk_dphi_bkg_syst[iCent][iPtCh][syst+1]);
    //  AddMaxSystematic (g_jet_trk_dphi_sig_syst[iCent][iPtCh][0], g_jet_trk_dphi_sig_syst[iCent][iPtCh][syst], g_jet_trk_dphi_sig_syst[iCent][iPtCh][syst+1]);
    //  AddMaxSystematic (g_jet_trk_dphi_iaa_syst[iCent][iPtCh][0], g_jet_trk_dphi_iaa_syst[iCent][iPtCh][syst], g_jet_trk_dphi_iaa_syst[iCent][iPtCh][syst+1]);

    //}
  }
*/



  {
    outFile->cd ();

    for (int iVar = 0; iVar < nVar; iVar++) {

      h_evt_counts_ref[iVar]->Write ();
      h_jet_counts_ref[iVar]->Write ();

      h_evt_counts_ref_bkg[iVar]->Write ();
      h_jet_counts_ref_bkg[iVar]->Write ();

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_evt_counts[iCent][iVar]->Write ();
        h_jet_counts[iCent][iVar]->Write ();

        h_evt_counts_bkg[iCent][iVar]->Write ();
        h_jet_counts_bkg[iCent][iVar]->Write ();

      }
    }

    for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
      for (int iVar = 0; iVar < nVar; iVar++) {

        h_jet_trk_dphi_ref[iPtCh][iVar]->Write ();
        h_jet_trk_dphi_ref_bkg[iPtCh][iVar]->Write ();
        h_jet_trk_dphi_ref_sig[iPtCh][iVar]->Write ();

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          h_jet_trk_dphi[iCent][iPtCh][iVar]->Write ();
          h_jet_trk_dphi_bkg[iCent][iPtCh][iVar]->Write ();
          h_jet_trk_dphi_sig[iCent][iPtCh][iVar]->Write ();
          h_jet_trk_dphi_iaa[iCent][iPtCh][iVar]->Write ();

        }
      }

      g_jet_trk_dphi_ref_syst[iPtCh][0]->Write (Form ("g_jet_trk_dphi_%s_ref_syst", pTChSelections[iPtCh].Data ()));
      g_jet_trk_dphi_ref_bkg_syst[iPtCh][0]->Write (Form ("g_jet_trk_dphi_%s_ref_bkg_syst", pTChSelections[iPtCh].Data ()));
      g_jet_trk_dphi_ref_sig_syst[iPtCh][0]->Write (Form ("g_jet_trk_dphi_%s_ref_sig_syst", pTChSelections[iPtCh].Data ()));

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        g_jet_trk_dphi_syst[iCent][iPtCh][0]->Write (Form ("g_jet_trk_dphi_%s_syst_iCent%i", pTChSelections[iPtCh].Data (), iCent));
        g_jet_trk_dphi_bkg_syst[iCent][iPtCh][0]->Write (Form ("g_jet_trk_dphi_%s_bkg_syst_iCent%i", pTChSelections[iPtCh].Data (), iCent));
        g_jet_trk_dphi_sig_syst[iCent][iPtCh][0]->Write (Form ("g_jet_trk_dphi_%s_sig_syst_iCent%i", pTChSelections[iPtCh].Data (), iCent));
        g_jet_trk_dphi_iaa_syst[iCent][iPtCh][0]->Write (Form ("g_jet_trk_dphi_%s_iaa_syst_iCent%i", pTChSelections[iPtCh].Data (), iCent));

      }
    }

    outFile->Close ();
  }
}


#endif
