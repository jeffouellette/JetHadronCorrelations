#ifndef __JetHadronCorrelator_ProcessCorrelations_C__
#define __JetHadronCorrelator_ProcessCorrelations_C__

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


void ProcessCorrelations (const char* outFileTag, const char* tag1, const char* tag2 = nullptr) {

  bool doMix = (tag2 == nullptr);

  TFile* inFile = nullptr;

  TH1D***  h_evt_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_jet_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_evt_counts_ref_bkg = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_jet_counts_ref_bkg = Get2DArray <TH1D*> (2, nVar);
  TH1D**** h_evt_counts         = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_counts         = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_evt_counts_bkg     = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_counts_bkg     = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);

  TH1D****  h_jet_trk_pt_ref        = Get3DArray <TH1D*> (2, nDir, nVar);
  TH2D****  h2_jet_trk_pt_cov_ref   = Get3DArray <TH2D*> (2, nDir, nVar);
  TH1D****  h_jet_trk_dphi_ref      = Get3DArray <TH1D*> (2, nPtChSelections, nVar);
  TH2D****  h2_jet_trk_dphi_cov_ref = Get3DArray <TH2D*> (2, nPtChSelections, nVar);

  TH1D****  h_jet_trk_pt_ref_bkg        = Get3DArray <TH1D*> (2, nDir, nVar);
  TH2D****  h2_jet_trk_pt_cov_ref_bkg   = Get3DArray <TH2D*> (2, nDir, nVar);
  TH1D****  h_jet_trk_dphi_ref_bkg      = Get3DArray <TH1D*> (2, nPtChSelections, nVar);
  TH2D****  h2_jet_trk_dphi_cov_ref_bkg = Get3DArray <TH2D*> (2, nPtChSelections, nVar);

  TH1D***** h_jet_trk_pt                = Get4DArray <TH1D*> (2, nDir, nZdcCentBins, nVar);
  TH2D***** h2_jet_trk_pt_cov           = Get4DArray <TH2D*> (2, nDir, nZdcCentBins, nVar);
  TH1D***** h_jet_trk_dphi              = Get4DArray <TH1D*> (2, nPtChSelections, nZdcCentBins, nVar);
  TH2D***** h2_jet_trk_dphi_cov         = Get4DArray <TH2D*> (2, nPtChSelections, nZdcCentBins, nVar);

  TH1D***** h_jet_trk_pt_bkg            = Get4DArray <TH1D*> (2, nDir, nZdcCentBins, nVar);
  TH2D***** h2_jet_trk_pt_cov_bkg       = Get4DArray <TH2D*> (2, nDir, nZdcCentBins, nVar);
  TH1D***** h_jet_trk_dphi_bkg          = Get4DArray <TH1D*> (2, nPtChSelections, nZdcCentBins, nVar);
  TH2D***** h2_jet_trk_dphi_cov_bkg     = Get4DArray <TH2D*> (2, nPtChSelections, nZdcCentBins, nVar);

  TH1D****  h_jet_trk_pt_ref_sig        = Get3DArray <TH1D*> (2, nDir, nVar);
  TH1D***** h_jet_trk_pt_sig            = Get4DArray <TH1D*> (2, nDir, nZdcCentBins, nVar);
  TH1D****  h_jet_trk_dphi_ref_sig      = Get3DArray <TH1D*> (2, nPtChSelections, nVar);
  TH1D***** h_jet_trk_dphi_sig          = Get4DArray <TH1D*> (2, nPtChSelections, nZdcCentBins, nVar);

  TH1D**    h_jet_trk_pt_ref_bbb        = Get1DArray <TH1D*> (nDir);
  TH1D***   h_jet_trk_pt_bbb            = Get2DArray <TH1D*> (nDir, nZdcCentBins);
  TH1D**    h_jet_trk_dphi_ref_bbb      = Get1DArray <TH1D*> (nPtChSelections);
  TH1D***   h_jet_trk_dphi_bbb          = Get2DArray <TH1D*> (nPtChSelections, nZdcCentBins);

  TF1**    f_jet_trk_pt_ref_bbb         = Get1DArray <TF1*> (nDir);
  TF1***   f_jet_trk_pt_bbb             = Get2DArray <TF1*> (nDir, nZdcCentBins);
  TF1**    f_jet_trk_dphi_ref_bbb       = Get1DArray <TF1*> (nPtChSelections);
  TF1***   f_jet_trk_dphi_bbb           = Get2DArray <TF1*> (nPtChSelections, nZdcCentBins);

  TH1D***** h_jet_trk_pt_iaa            = Get4DArray <TH1D*> (2, nDir, nZdcCentBins, nVar);
  TH1D***** h_jet_trk_dphi_iaa          = Get4DArray <TH1D*> (2, nPtChSelections, nZdcCentBins, nVar);


  TGAE***  g_jet_trk_pt_ref_syst        = Get2DArray <TGAE*> (nDir, nVar);
  TGAE***  g_jet_trk_dphi_ref_syst      = Get2DArray <TGAE*> (nPtChSelections, nVar);
  TGAE***  g_jet_trk_pt_ref_bkg_syst    = Get2DArray <TGAE*> (nDir, nVar);
  TGAE***  g_jet_trk_dphi_ref_bkg_syst  = Get2DArray <TGAE*> (nPtChSelections, nVar);
  TGAE**** g_jet_trk_pt_syst            = Get3DArray <TGAE*> (nDir, nZdcCentBins, nVar);
  TGAE**** g_jet_trk_dphi_syst          = Get3DArray <TGAE*> (nPtChSelections, nZdcCentBins, nVar);
  TGAE**** g_jet_trk_pt_bkg_syst        = Get3DArray <TGAE*> (nDir, nZdcCentBins, nVar);
  TGAE**** g_jet_trk_dphi_bkg_syst      = Get3DArray <TGAE*> (nPtChSelections, nZdcCentBins, nVar);

  TGAE***  g_jet_trk_pt_ref_sig_syst    = Get2DArray <TGAE*> (nDir, nVar);
  TGAE***  g_jet_trk_dphi_ref_sig_syst  = Get2DArray <TGAE*> (nPtChSelections, nVar);
  TGAE**** g_jet_trk_pt_sig_syst        = Get3DArray <TGAE*> (nDir, nZdcCentBins, nVar);
  TGAE**** g_jet_trk_dphi_sig_syst      = Get3DArray <TGAE*> (nPtChSelections, nZdcCentBins, nVar);

  TGAE**** g_jet_trk_pt_iaa_syst        = Get3DArray <TGAE*> (nDir, nZdcCentBins, nVar);
  TGAE**** g_jet_trk_dphi_iaa_syst      = Get3DArray <TGAE*> (nPtChSelections, nZdcCentBins, nVar);


  TGAE***  g_jet_trk_pt_ref_systTot       = Get2DArray <TGAE*> (nDir, 3);
  TGAE***  g_jet_trk_dphi_ref_systTot     = Get2DArray <TGAE*> (nPtChSelections, 3);
  TGAE***  g_jet_trk_pt_ref_bkg_systTot   = Get2DArray <TGAE*> (nDir, 3);
  TGAE***  g_jet_trk_dphi_ref_bkg_systTot = Get2DArray <TGAE*> (nPtChSelections, 3);
  TGAE**** g_jet_trk_pt_systTot           = Get3DArray <TGAE*> (nDir, nZdcCentBins, 3);
  TGAE**** g_jet_trk_dphi_systTot         = Get3DArray <TGAE*> (nPtChSelections, nZdcCentBins, 3);
  TGAE**** g_jet_trk_pt_bkg_systTot       = Get3DArray <TGAE*> (nDir, nZdcCentBins, 3);
  TGAE**** g_jet_trk_dphi_bkg_systTot     = Get3DArray <TGAE*> (nPtChSelections, nZdcCentBins, 3);

  TGAE***  g_jet_trk_pt_ref_sig_systTot   = Get2DArray <TGAE*> (nDir, 3);
  TGAE***  g_jet_trk_dphi_ref_sig_systTot = Get2DArray <TGAE*> (nPtChSelections, 3);
  TGAE**** g_jet_trk_pt_sig_systTot       = Get3DArray <TGAE*> (nDir, nZdcCentBins, 3);
  TGAE**** g_jet_trk_dphi_sig_systTot     = Get3DArray <TGAE*> (nPtChSelections, nZdcCentBins, 3);

  TGAE**** g_jet_trk_pt_iaa_systTot       = Get3DArray <TGAE*> (nDir, nZdcCentBins, 3);
  TGAE**** g_jet_trk_dphi_iaa_systTot     = Get3DArray <TGAE*> (nPtChSelections, nZdcCentBins, 3);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");

  for (int iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (int iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      {
        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s17_5TeV_hists.root", rootPath.Data (), tag1, var.Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts_ref[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h_evt_counts_ref_%s_%s", dType.Data (), var.Data ()));
        h_jet_counts_ref[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h_jet_counts_ref_%s_%s", dType.Data (), var.Data ()));

        for (int iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jet_trk_pt_ref[iDType][iDir][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_%s17",       dir.Data (), tag1, dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_ref_%s_%s",       dir.Data (), dType.Data (), var.Data ()));
          h2_jet_trk_pt_cov_ref[iDType][iDir][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_%s_cov_%s_%s17",  dir.Data (), tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_%s_cov_ref_%s_%s",  dir.Data (), dType.Data (), var.Data ()));

          CalcUncertainties (h_jet_trk_pt_ref[iDType][iDir][iVar],   h2_jet_trk_pt_cov_ref[iDType][iDir][iVar],   h_evt_counts_ref[iDType][iVar]);

        } // end loop over iDir

        for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          h_jet_trk_dphi_ref[iDType][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_%s17", ptch.Data (), tag1, dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_ref_%s_%s", ptch.Data (), dType.Data (), var.Data ()));
          h2_jet_trk_dphi_cov_ref[iDType][iPtCh][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_%s_cov_%s_%s17", ptch.Data (), tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_%s_cov_ref_%s_%s", ptch.Data (), dType.Data (), var.Data ()));

          CalcUncertainties (h_jet_trk_dphi_ref[iDType][iPtCh][iVar], h2_jet_trk_dphi_cov_ref[iDType][iPtCh][iVar], h_evt_counts_ref[iDType][iVar]);

        } // end loop over iPtCh

        inFile->Close ();

      }



      {
        TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/%s17_5TeV_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", var.Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts_ref_bkg[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_evt_counts_ref_bkg_%s_%s", dType.Data (), var.Data ()));
        h_jet_counts_ref_bkg[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_counts_ref_bkg_%s_%s", dType.Data (), var.Data ()));

        for (int iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jet_trk_pt_ref_bkg[iDType][iDir][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_%s17",       dir.Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_ref_bkg_%s_%s",       dir.Data (), dType.Data (), var.Data ()));
          h2_jet_trk_pt_cov_ref_bkg[iDType][iDir][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_%s_cov_%s_%s17",  dir.Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_%s_cov_ref_bkg_%s_%s",  dir.Data (), dType.Data (), var.Data ()));

          CalcUncertainties (h_jet_trk_pt_ref_bkg[iDType][iDir][iVar],   h2_jet_trk_pt_cov_ref_bkg[iDType][iDir][iVar],   h_evt_counts_ref_bkg[iDType][iVar]);

        } // end loop over iDir

        for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          h_jet_trk_dphi_ref_bkg[iDType][iPtCh][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_%s17",       ptch.Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_ref_bkg_%s_%s",       ptch.Data (), dType.Data (), var.Data ()));
          h2_jet_trk_dphi_cov_ref_bkg[iDType][iPtCh][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_%s_cov_%s_%s17",  ptch.Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_%s_cov_ref_bkg_%s_%s",  ptch.Data (), dType.Data (), var.Data ()));

          CalcUncertainties (h_jet_trk_dphi_ref_bkg[iDType][iPtCh][iVar], h2_jet_trk_dphi_cov_ref_bkg[iDType][iPtCh][iVar], h_evt_counts_ref_bkg[iDType][iVar]);

        } // end loop over iPtCh

        inFile->Close ();

      }



      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s16_5TeV_iCent%i_hists.root", rootPath.Data (), tag1, var.Data (), dType.Data (), iCent);
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s16", tag1, dType.Data ()))->Clone (Form ("h_evt_counts_pPb_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));
        h_jet_counts[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s16", tag1, dType.Data ()))->Clone (Form ("h_jet_counts_pPb_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));

        for (int iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jet_trk_pt[iDType][iDir][iCent][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_%s16",       dir.Data (), tag1, dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_pPb_iCent%i_%s_%s",       dir.Data (), iCent, dType.Data (), var.Data ()));
          h2_jet_trk_pt_cov[iDType][iDir][iCent][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_%s_cov_%s_%s16",  dir.Data (), tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_%s_cov_pPb_iCent%i_%s_%s",  dir.Data (), iCent, dType.Data (), var.Data ()));

          CalcUncertainties (h_jet_trk_pt[iDType][iDir][iCent][iVar],    h2_jet_trk_pt_cov[iDType][iDir][iCent][iVar],    h_evt_counts[iDType][iCent][iVar]);

        } // end loop over iDir

        for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          h_jet_trk_dphi[iDType][iPtCh][iCent][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_%s16",       ptch.Data (), tag1, dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_pPb_iCent%i_%s_%s",      ptch.Data (), iCent, dType.Data (), var.Data ()));
          h2_jet_trk_dphi_cov[iDType][iPtCh][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_%s_cov_%s_%s16",  ptch.Data (), tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_%s_cov_pPb_iCent%i_%s_%s", ptch.Data (), iCent, dType.Data (), var.Data ()));

          CalcUncertainties (h_jet_trk_dphi[iDType][iPtCh][iCent][iVar], h2_jet_trk_dphi_cov[iDType][iPtCh][iCent][iVar], h_evt_counts[iDType][iCent][iVar]);

        } // end loop over iPtCh

        inFile->Close ();

      } // end loop over iCent



      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/%s16_5TeV_iCent%i_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", var.Data (), dType.Data (), iCent);
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts_bkg[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_evt_counts_pPb_bkg_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));
        h_jet_counts_bkg[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_counts_pPb_bkg_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));

        for (int iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jet_trk_pt_bkg[iDType][iDir][iCent][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_%s16",       dir.Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_pPb_bkg_iCent%i_%s_%s",      dir.Data (), iCent, dType.Data (), var.Data ()));
          h2_jet_trk_pt_cov_bkg[iDType][iDir][iCent][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_%s_cov_%s_%s16",  dir.Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_%s_cov_pPb_bkg_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), var.Data ()));

          CalcUncertainties (h_jet_trk_pt_bkg[iDType][iDir][iCent][iVar],    h2_jet_trk_pt_cov_bkg[iDType][iDir][iCent][iVar],    h_evt_counts_bkg[iDType][iCent][iVar]);

        } // end loop over iDir

        for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          h_jet_trk_dphi_bkg[iDType][iPtCh][iCent][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_%s16",       ptch.Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_pPb_bkg_iCent%i_%s_%s",       ptch.Data (), iCent, dType.Data (), var.Data ()));
          h2_jet_trk_dphi_cov_bkg[iDType][iPtCh][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_%s_cov_%s_%s16",  ptch.Data (), doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_%s_cov_pPb_bkg_iCent%i_%s_%s",  ptch.Data (), iCent, dType.Data (), var.Data ()));

          CalcUncertainties (h_jet_trk_dphi_bkg[iDType][iPtCh][iCent][iVar], h2_jet_trk_dphi_cov_bkg[iDType][iPtCh][iCent][iVar], h_evt_counts_bkg[iDType][iCent][iVar]);

        } // end loop over iPtCh

        inFile->Close ();

      } // end loop over iCent

    } // end loop over iVar

  } // end loop over iDType



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // DERIVE BIN-BY-BIN UNFOLDING FACTORS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    int iVarNum = 0;
    while (iVarNum < nVar && strcmp (variations[iVarNum], "MCTruthLevel") != 0) iVarNum++;
    if (iVarNum == nVar) {
      std::cout << "Cannot find MC truth jet result? Please check!" << std::endl;
      return;
    }
    int iVarDen = 0;
    while (iVarDen < nVar && strcmp (variations[iVarDen], "MCBbyBReco") != 0) iVarDen++;
    if (iVarDen == nVar) {
      std::cout << "Cannot find MC reco. jet result? Please check!" << std::endl;
      return;
    }

    TF1* f = nullptr;

    const TString funcStr = "[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)";

    for (int iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      h_jet_trk_pt_ref_bbb[iDir] = (TH1D*) h_jet_trk_pt_ref[1][iDir][iVarNum]->Clone (Form ("h_jet_trk_pt_%s_ref_bbb", dir.Data ()));
      h_jet_trk_pt_ref_bbb[iDir]->Divide (h_jet_trk_pt_ref[1][iDir][iVarDen]);

      f = new TF1 (Form ("f_jet_trk_pt_%s_ref_bbb", dir.Data ()), funcStr, h_jet_trk_pt_ref_bbb[iDir]->GetXaxis ()->GetBinLowEdge (1), h_jet_trk_pt_ref_bbb[iDir]->GetXaxis ()->GetBinLowEdge (h_jet_trk_pt_ref_bbb[iDir]->GetNbinsX () + 1));
      h_jet_trk_pt_ref_bbb[iDir]->Fit (f, "RN0Q");
      f_jet_trk_pt_ref_bbb[iDir] = f;

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_jet_trk_pt_bbb[iDir][iCent] = (TH1D*) h_jet_trk_pt[1][iDir][iCent][iVarNum]->Clone (Form ("h_jet_trk_pt_%s_pPb_bbb_iCent%i", dir.Data (), iCent));
        h_jet_trk_pt_bbb[iDir][iCent]->Divide (h_jet_trk_pt[1][iDir][iCent][iVarDen]);

        f = new TF1 (Form ("f_jet_trk_pt_%s_pPb_bbb_iCent%i", dir.Data (), iCent), funcStr, h_jet_trk_pt_bbb[iDir][iCent]->GetXaxis ()->GetBinLowEdge (1), h_jet_trk_pt_bbb[iDir][iCent]->GetXaxis ()->GetBinLowEdge (h_jet_trk_pt_bbb[iDir][iCent]->GetNbinsX () + 1));
        h_jet_trk_pt_bbb[iDir][iCent]->Fit (f, "RN0Q");
        f_jet_trk_pt_bbb[iDir][iCent] = f;

      } // end loop over iCent

    } // end loop over iDir
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CALCULATE SIGNAL YIELDS BY SUBTRACTING THE BACKGROUND COMPONENT, UNFOLD, AND CALCULATE RATIOS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (int iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      const bool hasRefBkg = (variationsWithNoppBkgd.count (var) == 0);
      const bool hasBkg = (variationsWithNopPbBkgd.count (var) == 0);
      const bool doUnfold = (variationsWithNoUnfold.count (var) == 0);

      for (int iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        if (doUnfold) {
          BinByBinUnfold (h_jet_trk_pt_ref[iDType][iDir][iVar],     f_jet_trk_pt_ref_bbb[iDir]);
          BinByBinUnfold (h_jet_trk_pt_ref_bkg[iDType][iDir][iVar], f_jet_trk_pt_ref_bbb[iDir]);
        }

        h_jet_trk_pt_ref_sig[iDType][iDir][iVar] = (TH1D*) h_jet_trk_pt_ref[iDType][iDir][iVar]->Clone (Form ("h_jet_trk_pt_%s_ref_sig_%s_%s", dir.Data (), dType.Data (), var.Data ()));
        if (hasRefBkg)
          h_jet_trk_pt_ref_sig[iDType][iDir][iVar]->Add (h_jet_trk_pt_ref_bkg[iDType][iDir][iVar], -1);

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          if (doUnfold) {
            BinByBinUnfold (h_jet_trk_pt[iDType][iDir][iCent][iVar],      f_jet_trk_pt_bbb[iDir][iCent]);
            BinByBinUnfold (h_jet_trk_pt_bkg[iDType][iDir][iCent][iVar],  f_jet_trk_pt_bbb[iDir][iCent]);
          }

          h_jet_trk_pt_sig[iDType][iDir][iCent][iVar] = (TH1D*) h_jet_trk_pt[iDType][iDir][iCent][iVar]->Clone (Form ("h_jet_trk_pt_%s_pPb_sig_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), var.Data ()));
          if (hasBkg)
            h_jet_trk_pt_sig[iDType][iDir][iCent][iVar]->Add (h_jet_trk_pt_bkg[iDType][iDir][iCent][iVar], -1);

          h_jet_trk_pt_iaa[iDType][iDir][iCent][iVar] = (TH1D*) h_jet_trk_pt_sig[iDType][iDir][iCent][iVar]->Clone (Form ("h_jet_trk_pt_%s_iaa_iCent%i_%s_%s", dir.Data (), iCent, dType.Data (), var.Data ()));
          h_jet_trk_pt_iaa[iDType][iDir][iCent][iVar]->Divide (h_jet_trk_pt_ref_sig[iDType][iDir][iVar]);

        } // end loop over iCent

      } // end loop over iDir

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        const TString ptch = pTChSelections[iPtCh].Data ();

        h_jet_trk_dphi_ref_sig[iDType][iPtCh][iVar] = (TH1D*) h_jet_trk_dphi_ref[iDType][iPtCh][iVar]->Clone (Form ("h_jet_trk_dphi_%s_ref_sig_%s_%s", ptch.Data (), dType.Data (), var.Data ()));
        if (hasRefBkg)
          h_jet_trk_dphi_ref_sig[iDType][iPtCh][iVar]->Add (h_jet_trk_dphi_ref_bkg[iDType][iPtCh][iVar], -1);

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          h_jet_trk_dphi_sig[iDType][iPtCh][iCent][iVar] = (TH1D*) h_jet_trk_dphi[iDType][iPtCh][iCent][iVar]->Clone (Form ("h_jet_trk_dphi_%s_pPb_sig_iCent%i_%s_%s", ptch.Data (), iCent, dType.Data (), var.Data ()));
          if (hasBkg)
            h_jet_trk_dphi_sig[iDType][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi_bkg[iDType][iPtCh][iCent][iVar], -1);

          h_jet_trk_dphi_iaa[iDType][iPtCh][iCent][iVar] = (TH1D*) h_jet_trk_dphi_sig[iDType][iPtCh][iCent][iVar]->Clone (Form ("h_jet_trk_dphi_%s_iaa_iCent%i_%s_%s", ptch.Data (), iCent, dType.Data (), var.Data ()));
          h_jet_trk_dphi_iaa[iDType][iPtCh][iCent][iVar]->Divide (h_jet_trk_dphi_ref_sig[iDType][iPtCh][iVar]);

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iVar

  } // end loop over iDType



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE TOTAL SYSTEMATIC UNCERTAINTIES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iDir = 0; iDir < nDir; iDir++) {

    g_jet_trk_pt_ref_syst[iDir][0]     = make_graph (h_jet_trk_pt_ref[0][iDir][0]);
    g_jet_trk_pt_ref_bkg_syst[iDir][0] = make_graph (h_jet_trk_pt_ref_bkg[0][iDir][0]);
    g_jet_trk_pt_ref_sig_syst[iDir][0] = make_graph (h_jet_trk_pt_ref_sig[0][iDir][0]);

    ResetTGAEErrors (g_jet_trk_pt_ref_syst[iDir][0]);
    ResetTGAEErrors (g_jet_trk_pt_ref_bkg_syst[iDir][0]);
    ResetTGAEErrors (g_jet_trk_pt_ref_sig_syst[iDir][0]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      g_jet_trk_pt_syst[iDir][iCent][0]        = make_graph (h_jet_trk_pt[0][iDir][iCent][0]);
      g_jet_trk_pt_bkg_syst[iDir][iCent][0]    = make_graph (h_jet_trk_pt_bkg[0][iDir][iCent][0]);
      g_jet_trk_pt_sig_syst[iDir][iCent][0]    = make_graph (h_jet_trk_pt_sig[0][iDir][iCent][0]);
      g_jet_trk_pt_iaa_syst[iDir][iCent][0]    = make_graph (h_jet_trk_pt_iaa[0][iDir][iCent][0]);

      ResetTGAEErrors (g_jet_trk_pt_syst[iDir][iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_bkg_syst[iDir][iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_sig_syst[iDir][iCent][0]);
      ResetTGAEErrors (g_jet_trk_pt_iaa_syst[iDir][iCent][0]);

    } // end loop over iCent

  } // end loop over iDir

  for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

    g_jet_trk_dphi_ref_syst[iPtCh][0]     = make_graph (h_jet_trk_dphi_ref[0][iPtCh][0]);
    g_jet_trk_dphi_ref_bkg_syst[iPtCh][0] = make_graph (h_jet_trk_dphi_ref_bkg[0][iPtCh][0]);
    g_jet_trk_dphi_ref_sig_syst[iPtCh][0] = make_graph (h_jet_trk_dphi_ref_sig[0][iPtCh][0]);

    ResetTGAEErrors (g_jet_trk_dphi_ref_syst[iPtCh][0]);
    ResetTGAEErrors (g_jet_trk_dphi_ref_bkg_syst[iPtCh][0]);
    ResetTGAEErrors (g_jet_trk_dphi_ref_sig_syst[iPtCh][0]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      g_jet_trk_dphi_syst[iPtCh][iCent][0]      = make_graph (h_jet_trk_dphi[0][iPtCh][iCent][0]);
      g_jet_trk_dphi_bkg_syst[iPtCh][iCent][0]  = make_graph (h_jet_trk_dphi_bkg[0][iPtCh][iCent][0]);
      g_jet_trk_dphi_sig_syst[iPtCh][iCent][0]  = make_graph (h_jet_trk_dphi_sig[0][iPtCh][iCent][0]);
      g_jet_trk_dphi_iaa_syst[iPtCh][iCent][0]  = make_graph (h_jet_trk_dphi_iaa[0][iPtCh][iCent][0]);

      ResetTGAEErrors (g_jet_trk_dphi_syst[iPtCh][iCent][0]);
      ResetTGAEErrors (g_jet_trk_dphi_bkg_syst[iPtCh][iCent][0]);
      ResetTGAEErrors (g_jet_trk_dphi_sig_syst[iPtCh][iCent][0]);
      ResetTGAEErrors (g_jet_trk_dphi_iaa_syst[iPtCh][iCent][0]);

    } // end loop over iCent

  } // end loop over iPtCh



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE SEPARATELY.
  // THEN CALCULATES SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE BY TAKING DIFFERENCES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iDType = 0; iDType < 2; iDType++) {

    for (int iVar = 1; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
        continue;

      for (int iDir = 0; iDir < nDir; iDir++) {

        g_jet_trk_pt_ref_syst[iDir][iVar]     = new TGAE ();
        g_jet_trk_pt_ref_bkg_syst[iDir][iVar] = new TGAE ();
        g_jet_trk_pt_ref_sig_syst[iDir][iVar] = new TGAE ();

        CalcSystematics (g_jet_trk_pt_ref_syst[iDir][iVar],     h_jet_trk_pt_ref[iDType][iDir][0],      h_jet_trk_pt_ref[iDType][iDir][iVar]);
        CalcSystematics (g_jet_trk_pt_ref_bkg_syst[iDir][iVar], h_jet_trk_pt_ref_bkg[iDType][iDir][0],  h_jet_trk_pt_ref_bkg[iDType][iDir][iVar]);
        CalcSystematics (g_jet_trk_pt_ref_sig_syst[iDir][iVar], h_jet_trk_pt_ref_sig[iDType][iDir][0],  h_jet_trk_pt_ref_sig[iDType][iDir][iVar]);

        // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
        if (variationsThatDontCancelInSig.count (var) != 0) {
          ResetTGAEErrors (g_jet_trk_pt_ref_sig_syst[iDir][iVar]);
          AddErrorsInQuadrature (g_jet_trk_pt_ref_sig_syst[iDir][iVar], g_jet_trk_pt_ref_syst[iDir][iVar], false, false);
          AddErrorsInQuadrature (g_jet_trk_pt_ref_sig_syst[iDir][iVar], g_jet_trk_pt_ref_bkg_syst[iDir][iVar], false, false);
        }

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          g_jet_trk_pt_syst[iDir][iCent][iVar]     = new TGAE ();
          g_jet_trk_pt_bkg_syst[iDir][iCent][iVar] = new TGAE ();
          g_jet_trk_pt_sig_syst[iDir][iCent][iVar] = new TGAE ();
          g_jet_trk_pt_iaa_syst[iDir][iCent][iVar] = new TGAE ();

          CalcSystematics (g_jet_trk_pt_syst[iDir][iCent][iVar],     h_jet_trk_pt[iDType][iDir][iCent][0],     h_jet_trk_pt[iDType][iDir][iCent][iVar]);
          CalcSystematics (g_jet_trk_pt_bkg_syst[iDir][iCent][iVar], h_jet_trk_pt_bkg[iDType][iDir][iCent][0], h_jet_trk_pt_bkg[iDType][iDir][iCent][iVar]);
          CalcSystematics (g_jet_trk_pt_sig_syst[iDir][iCent][iVar], h_jet_trk_pt_sig[iDType][iDir][iCent][0], h_jet_trk_pt_sig[iDType][iDir][iCent][iVar]);
          CalcSystematics (g_jet_trk_pt_iaa_syst[iDir][iCent][iVar], h_jet_trk_pt_iaa[iDType][iDir][iCent][0], h_jet_trk_pt_iaa[iDType][iDir][iCent][iVar]);

          // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
          if (variationsThatDontCancelInSig.count (var) != 0) {
            ResetTGAEErrors (g_jet_trk_pt_sig_syst[iDir][iCent][iVar]);
            AddErrorsInQuadrature (g_jet_trk_pt_sig_syst[iDir][iCent][iVar], g_jet_trk_pt_syst[iDir][iCent][iVar], false, false);
            AddErrorsInQuadrature (g_jet_trk_pt_sig_syst[iDir][iCent][iVar], g_jet_trk_pt_bkg_syst[iDir][iCent][iVar], false, false);
          }
          // allow some uncertainties to not cancel in the p+Pb / pp ratio by overwriting the current uncertainties
          if (variationsThatDontCancelInRatio.count (var) != 0) {
            ResetTGAEErrors (g_jet_trk_pt_iaa_syst[iDir][iCent][iVar]);
            AddRelErrorsInQuadrature (g_jet_trk_pt_iaa_syst[iDir][iCent][iVar], g_jet_trk_pt_ref_sig_syst[iDir][iVar], false, false);
            AddRelErrorsInQuadrature (g_jet_trk_pt_iaa_syst[iDir][iCent][iVar], g_jet_trk_pt_sig_syst[iDir][iCent][iVar], false, false);
          }

        } // end loop over iCent

      } // end loop over iDir

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        g_jet_trk_dphi_ref_syst[iPtCh][iVar]      = new TGAE ();
        g_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar]  = new TGAE ();
        g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar]  = new TGAE ();

        CalcSystematics (g_jet_trk_dphi_ref_syst[iPtCh][iVar],      h_jet_trk_dphi_ref[iDType][iPtCh][0],     h_jet_trk_dphi_ref[iDType][iPtCh][iVar]);
        CalcSystematics (g_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar],  h_jet_trk_dphi_ref_bkg[iDType][iPtCh][0], h_jet_trk_dphi_ref_bkg[iDType][iPtCh][iVar]);
        CalcSystematics (g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar],  h_jet_trk_dphi_ref_sig[iDType][iPtCh][0], h_jet_trk_dphi_ref_sig[iDType][iPtCh][iVar]);

        // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
        if (variationsThatDontCancelInSig.count (var) != 0) {
          ResetTGAEErrors (g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar]);
          AddErrorsInQuadrature (g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar], g_jet_trk_dphi_ref_syst[iPtCh][iVar], false, false);
          AddErrorsInQuadrature (g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar], g_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar], false, false);
        }

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          g_jet_trk_dphi_syst[iPtCh][iCent][iVar]     = new TGAE ();
          g_jet_trk_dphi_bkg_syst[iPtCh][iCent][iVar] = new TGAE ();
          g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar] = new TGAE ();
          g_jet_trk_dphi_iaa_syst[iPtCh][iCent][iVar] = new TGAE ();

          CalcSystematics (g_jet_trk_dphi_syst[iPtCh][iCent][iVar],     h_jet_trk_dphi[iDType][iPtCh][iCent][0],      h_jet_trk_dphi[iDType][iPtCh][iCent][iVar]);
          CalcSystematics (g_jet_trk_dphi_bkg_syst[iPtCh][iCent][iVar], h_jet_trk_dphi_bkg[iDType][iPtCh][iCent][0],  h_jet_trk_dphi_bkg[iDType][iPtCh][iCent][iVar]);
          CalcSystematics (g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar], h_jet_trk_dphi_sig[iDType][iPtCh][iCent][0],  h_jet_trk_dphi_sig[iDType][iPtCh][iCent][iVar]);
          CalcSystematics (g_jet_trk_dphi_iaa_syst[iPtCh][iCent][iVar], h_jet_trk_dphi_iaa[iDType][iPtCh][iCent][0],  h_jet_trk_dphi_iaa[iDType][iPtCh][iCent][iVar]);

          // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
          if (variationsThatDontCancelInSig.count (var) != 0) {
            ResetTGAEErrors (g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar]);
            AddErrorsInQuadrature (g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar], g_jet_trk_dphi_syst[iPtCh][iCent][iVar], false, false);
            AddErrorsInQuadrature (g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar], g_jet_trk_dphi_bkg_syst[iPtCh][iCent][iVar], false, false);
          }
          // allow some uncertainties to not cancel in the p+Pb / pp ratio by overwriting the current uncertainties
          if (variationsThatDontCancelInRatio.count (var) != 0) {
            ResetTGAEErrors (g_jet_trk_dphi_iaa_syst[iPtCh][iCent][iVar]);
            AddRelErrorsInQuadrature (g_jet_trk_dphi_iaa_syst[iPtCh][iCent][iVar], g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar], false, false);
            AddRelErrorsInQuadrature (g_jet_trk_dphi_iaa_syst[iPtCh][iCent][iVar], g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar], false, false);
          }

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iVar

  } // end loop over iDType



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // SYSTEMATIC UNCERTAINTIES DERIVED IN MC MUST HAVE CENTRAL VALUES SET BY CENTRAL VALUES IN DATA
  // THE FINAL UNCERTAINTY IS ASSIGNED TO MATCH THE FRACTIONAL UNCERTAINTY IN MC
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iVar = 1; iVar < nVar; iVar++) {

    const TString var = variations[iVar];

    if (dataVariations.count (var) > 0 || mcVariations.count (var) == 0)
      continue; // skip variations already evaluated in data or that are not evaluated in MC

    for (int iDir = 0; iDir < nDir; iDir++) {

      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ref_syst[iDir][iVar],      h_jet_trk_pt_ref[0][iDir][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ref_bkg_syst[iDir][iVar],  h_jet_trk_pt_ref_bkg[0][iDir][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ref_sig_syst[iDir][iVar],  h_jet_trk_pt_ref_sig[0][iDir][0]);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_syst[iDir][iCent][iVar],     h_jet_trk_pt[0][iDir][iCent][0]);
        SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_bkg_syst[iDir][iCent][iVar], h_jet_trk_pt_bkg[0][iDir][iCent][0]);
        SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_sig_syst[iDir][iCent][iVar], h_jet_trk_pt_sig[0][iDir][iCent][0]);
        SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_iaa_syst[iDir][iCent][iVar], h_jet_trk_pt_iaa[0][iDir][iCent][0]);

      } // end loop over iCent

    } // end loop over iDir

    for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      SetCentralValuesKeepRelativeErrors (g_jet_trk_dphi_ref_syst[iPtCh][iVar],     h_jet_trk_dphi_ref[0][iPtCh][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar], h_jet_trk_dphi_ref_bkg[0][iPtCh][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar], h_jet_trk_dphi_ref_sig[0][iPtCh][0]);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        SetCentralValuesKeepRelativeErrors (g_jet_trk_dphi_syst[iPtCh][iCent][iVar],      h_jet_trk_dphi[0][iPtCh][iCent][0]);
        SetCentralValuesKeepRelativeErrors (g_jet_trk_dphi_bkg_syst[iPtCh][iCent][iVar],  h_jet_trk_dphi_bkg[0][iPtCh][iCent][0]);
        SetCentralValuesKeepRelativeErrors (g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar],  h_jet_trk_dphi_sig[0][iPtCh][iCent][0]);
        SetCentralValuesKeepRelativeErrors (g_jet_trk_dphi_iaa_syst[iPtCh][iCent][iVar],  h_jet_trk_dphi_iaa[0][iPtCh][iCent][0]);

      } // end loop over iCent

    } // end loop over iPtCh
    
  } // end loop over iVar



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // THESE GRAPHS STORE SUMMARY SYSTEMATIC UNCERTAINTIES FOR EACH CATEGORY: TRACKING, JETS, & MIXING.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iTotVar = 0; iTotVar < 3; iTotVar++) {

    for (int iDir = 0; iDir < nDir; iDir++) {

      g_jet_trk_pt_ref_systTot[iDir][iTotVar]     = (TGAE*) g_jet_trk_pt_ref_syst[iDir][0]->Clone ();
      g_jet_trk_pt_ref_bkg_systTot[iDir][iTotVar] = (TGAE*) g_jet_trk_pt_ref_bkg_syst[iDir][0]->Clone ();
      g_jet_trk_pt_ref_sig_systTot[iDir][iTotVar] = (TGAE*) g_jet_trk_pt_ref_sig_syst[iDir][0]->Clone ();

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        g_jet_trk_pt_systTot[iDir][iCent][iTotVar]     = (TGAE*) g_jet_trk_pt_syst[iDir][iCent][0]->Clone ();
        g_jet_trk_pt_bkg_systTot[iDir][iCent][iTotVar] = (TGAE*) g_jet_trk_pt_bkg_syst[iDir][iCent][0]->Clone ();
        g_jet_trk_pt_sig_systTot[iDir][iCent][iTotVar] = (TGAE*) g_jet_trk_pt_sig_syst[iDir][iCent][0]->Clone ();
        g_jet_trk_pt_iaa_systTot[iDir][iCent][iTotVar] = (TGAE*) g_jet_trk_pt_iaa_syst[iDir][iCent][0]->Clone ();

      } // end loop over iCent

    } // end loop over iDir

    for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      g_jet_trk_dphi_ref_systTot[iPtCh][iTotVar]     = (TGAE*) g_jet_trk_dphi_ref_syst[iPtCh][0]->Clone ();
      g_jet_trk_dphi_ref_bkg_systTot[iPtCh][iTotVar] = (TGAE*) g_jet_trk_dphi_ref_bkg_syst[iPtCh][0]->Clone ();
      g_jet_trk_dphi_ref_sig_systTot[iPtCh][iTotVar] = (TGAE*) g_jet_trk_dphi_ref_sig_syst[iPtCh][0]->Clone ();

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        g_jet_trk_dphi_systTot[iPtCh][iCent][iTotVar]     = (TGAE*) g_jet_trk_dphi_syst[iPtCh][iCent][0]->Clone ();
        g_jet_trk_dphi_bkg_systTot[iPtCh][iCent][iTotVar] = (TGAE*) g_jet_trk_dphi_bkg_syst[iPtCh][iCent][0]->Clone ();
        g_jet_trk_dphi_sig_systTot[iPtCh][iCent][iTotVar] = (TGAE*) g_jet_trk_dphi_sig_syst[iPtCh][iCent][0]->Clone ();
        g_jet_trk_dphi_iaa_systTot[iPtCh][iCent][iTotVar] = (TGAE*) g_jet_trk_dphi_iaa_syst[iPtCh][iCent][0]->Clone ();

      } // end loop over iCent

    } // end loop over iPtCh

  } // end loop over iTotVar



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM EACH CATEGORY INTO ONE OF THREE GRAPHS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (std::vector <TString> vgroup : variationGroups) {

    std::vector <int> iVars = {};
    for (TString s : vgroup) {
      const int iVar = GetVariationN (s);
      if (0 <= iVar && iVar < nVar)
        iVars.push_back (iVar);
    }

    if (iVars.size () == 0)
      continue;

    else if (iVars.size () == 1) {

      const int iVar = iVars[0];
      const TString var = variations[iVar];
      const int iTotVar = (IsJetsVariation (var) ? 2 : (IsTrackingVariation (var) ? 1 : 0));

      for (int iDir = 0; iDir < nDir; iDir++) {

        AddErrorsInQuadrature (g_jet_trk_pt_ref_systTot[iDir][iTotVar],      g_jet_trk_pt_ref_syst[iDir][iVar], false, true);
        AddErrorsInQuadrature (g_jet_trk_pt_ref_bkg_systTot[iDir][iTotVar],  g_jet_trk_pt_ref_bkg_syst[iDir][iVar], false, true);
        AddErrorsInQuadrature (g_jet_trk_pt_ref_sig_systTot[iDir][iTotVar],  g_jet_trk_pt_ref_sig_syst[iDir][iVar], false, true);

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          AddErrorsInQuadrature (g_jet_trk_pt_systTot[iDir][iCent][iTotVar],     g_jet_trk_pt_syst[iDir][iCent][iVar], false, true);
          AddErrorsInQuadrature (g_jet_trk_pt_bkg_systTot[iDir][iCent][iTotVar], g_jet_trk_pt_bkg_syst[iDir][iCent][iVar], false, true);
          AddErrorsInQuadrature (g_jet_trk_pt_sig_systTot[iDir][iCent][iTotVar], g_jet_trk_pt_sig_syst[iDir][iCent][iVar], false, true);
          AddErrorsInQuadrature (g_jet_trk_pt_iaa_systTot[iDir][iCent][iTotVar], g_jet_trk_pt_iaa_syst[iDir][iCent][iVar], false, true);

        } // end loop over iCent

      } // end loop over iDir

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        AddErrorsInQuadrature (g_jet_trk_dphi_ref_systTot[iPtCh][iTotVar],     g_jet_trk_dphi_ref_syst[iPtCh][iVar], false, true);
        AddErrorsInQuadrature (g_jet_trk_dphi_ref_bkg_systTot[iPtCh][iTotVar], g_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar], false, true);
        AddErrorsInQuadrature (g_jet_trk_dphi_ref_sig_systTot[iPtCh][iTotVar], g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar], false, true);

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          AddErrorsInQuadrature (g_jet_trk_dphi_systTot[iPtCh][iCent][iTotVar],      g_jet_trk_dphi_syst[iPtCh][iCent][iVar], false, true);
          AddErrorsInQuadrature (g_jet_trk_dphi_bkg_systTot[iPtCh][iCent][iTotVar],  g_jet_trk_dphi_bkg_syst[iPtCh][iCent][iVar], false, true);
          AddErrorsInQuadrature (g_jet_trk_dphi_sig_systTot[iPtCh][iCent][iTotVar],  g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar], false, true);
          AddErrorsInQuadrature (g_jet_trk_dphi_iaa_systTot[iPtCh][iCent][iTotVar],  g_jet_trk_dphi_iaa_syst[iPtCh][iCent][iVar], false, true);

        } // end loop over iCent

      } // end loop over iPtCh

    }

    else {

      const TString var = variations[iVars[0]];
      const int iTotVar = (IsJetsVariation (var) ? 2 : (IsTrackingVariation (var) ? 1 : 0));

      for (int iDir = 0; iDir < nDir; iDir++) {

        AddErrorsInQuadrature (g_jet_trk_pt_ref_systTot[iDir][iTotVar],      g_jet_trk_pt_ref_syst[iDir],     &iVars, true);
        AddErrorsInQuadrature (g_jet_trk_pt_ref_bkg_systTot[iDir][iTotVar],  g_jet_trk_pt_ref_bkg_syst[iDir], &iVars, true);
        AddErrorsInQuadrature (g_jet_trk_pt_ref_sig_systTot[iDir][iTotVar],  g_jet_trk_pt_ref_sig_syst[iDir], &iVars, true);

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          AddErrorsInQuadrature (g_jet_trk_pt_systTot[iDir][iCent][iTotVar],     g_jet_trk_pt_syst[iDir][iCent],      &iVars, true);
          AddErrorsInQuadrature (g_jet_trk_pt_bkg_systTot[iDir][iCent][iTotVar], g_jet_trk_pt_bkg_syst[iDir][iCent],  &iVars, true);
          AddErrorsInQuadrature (g_jet_trk_pt_sig_systTot[iDir][iCent][iTotVar], g_jet_trk_pt_sig_syst[iDir][iCent],  &iVars, true);
          AddErrorsInQuadrature (g_jet_trk_pt_iaa_systTot[iDir][iCent][iTotVar], g_jet_trk_pt_iaa_syst[iDir][iCent],  &iVars, true);

        } // end loop over iCent

      } // end loop over iDir

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        AddErrorsInQuadrature (g_jet_trk_dphi_ref_systTot[iPtCh][iTotVar],     g_jet_trk_dphi_ref_syst[iPtCh],      &iVars, true);
        AddErrorsInQuadrature (g_jet_trk_dphi_ref_bkg_systTot[iPtCh][iTotVar], g_jet_trk_dphi_ref_bkg_syst[iPtCh],  &iVars, true);
        AddErrorsInQuadrature (g_jet_trk_dphi_ref_sig_systTot[iPtCh][iTotVar], g_jet_trk_dphi_ref_sig_syst[iPtCh],  &iVars, true);

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          AddErrorsInQuadrature (g_jet_trk_dphi_systTot[iPtCh][iCent][iTotVar],      g_jet_trk_dphi_syst[iPtCh][iCent],     &iVars, true);
          AddErrorsInQuadrature (g_jet_trk_dphi_bkg_systTot[iPtCh][iCent][iTotVar],  g_jet_trk_dphi_bkg_syst[iPtCh][iCent], &iVars, true);
          AddErrorsInQuadrature (g_jet_trk_dphi_sig_systTot[iPtCh][iCent][iTotVar],  g_jet_trk_dphi_sig_syst[iPtCh][iCent], &iVars, true);
          AddErrorsInQuadrature (g_jet_trk_dphi_iaa_systTot[iPtCh][iCent][iTotVar],  g_jet_trk_dphi_iaa_syst[iPtCh][iCent], &iVars, true);

        } // end loop over iCent

      } // end loop over iPtCh

    }

  } // end loop over iVar



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM ALL SOURCES IN QUADRATURE, STORING RESULTS IN A SINGLE GRAPH
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iTotVar = 0; iTotVar < 3; iTotVar++) {

    for (int iDir = 0; iDir < nDir; iDir++) {
  
      AddErrorsInQuadrature (g_jet_trk_pt_ref_syst[iDir][0],      g_jet_trk_pt_ref_systTot[iDir][iTotVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_ref_bkg_syst[iDir][0],  g_jet_trk_pt_ref_bkg_systTot[iDir][iTotVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_ref_sig_syst[iDir][0],  g_jet_trk_pt_ref_sig_systTot[iDir][iTotVar]);
  
      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  
        AddErrorsInQuadrature (g_jet_trk_pt_syst[iDir][iCent][0],     g_jet_trk_pt_systTot[iDir][iCent][iTotVar]);
        AddErrorsInQuadrature (g_jet_trk_pt_bkg_syst[iDir][iCent][0], g_jet_trk_pt_bkg_systTot[iDir][iCent][iTotVar]);
        AddErrorsInQuadrature (g_jet_trk_pt_sig_syst[iDir][iCent][0], g_jet_trk_pt_sig_systTot[iDir][iCent][iTotVar]);
        AddErrorsInQuadrature (g_jet_trk_pt_iaa_syst[iDir][iCent][0], g_jet_trk_pt_iaa_systTot[iDir][iCent][iTotVar]);
  
      } // end loop over iCent
  
    } // end loop over iDir
  
    for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
  
      AddErrorsInQuadrature (g_jet_trk_dphi_ref_syst[iPtCh][0],     g_jet_trk_dphi_ref_systTot[iPtCh][iTotVar]);
      AddErrorsInQuadrature (g_jet_trk_dphi_ref_bkg_syst[iPtCh][0], g_jet_trk_dphi_ref_bkg_systTot[iPtCh][iTotVar]);
      AddErrorsInQuadrature (g_jet_trk_dphi_ref_sig_syst[iPtCh][0], g_jet_trk_dphi_ref_sig_systTot[iPtCh][iTotVar]);
  
      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  
        AddErrorsInQuadrature (g_jet_trk_dphi_syst[iPtCh][iCent][0],      g_jet_trk_dphi_systTot[iPtCh][iCent][iTotVar]);
        AddErrorsInQuadrature (g_jet_trk_dphi_bkg_syst[iPtCh][iCent][0],  g_jet_trk_dphi_bkg_systTot[iPtCh][iCent][iTotVar]);
        AddErrorsInQuadrature (g_jet_trk_dphi_sig_syst[iPtCh][iCent][0],  g_jet_trk_dphi_sig_systTot[iPtCh][iCent][iTotVar]);
        AddErrorsInQuadrature (g_jet_trk_dphi_iaa_syst[iPtCh][iCent][0],  g_jet_trk_dphi_iaa_systTot[iPtCh][iCent][iTotVar]);
  
      } // end loop over iCent
  
    } // end loop over iPtCh

  } // end loop over iTotVar



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // FINALLY WRITE OUT EVERYTHING TO A SINGLE ROOT FILE WITH ALL THE RESULTS.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    outFile->cd ();

    for (int iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if (dataVariations.count (var) == 0 && mcVariations.count (var) == 0)
        continue;

      for (int iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        g_jet_trk_pt_ref_syst[iDir][iVar]->Write     (Form ("g_jet_trk_pt_%s_ref_syst_%s",      dir.Data (), var.Data ()));
        g_jet_trk_pt_ref_bkg_syst[iDir][iVar]->Write (Form ("g_jet_trk_pt_%s_ref_bkg_syst_%s",  dir.Data (), var.Data ()));
        g_jet_trk_pt_ref_sig_syst[iDir][iVar]->Write (Form ("g_jet_trk_pt_%s_ref_sig_syst_%s",  dir.Data (), var.Data ()));

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          g_jet_trk_pt_syst[iDir][iCent][iVar]->Write      (Form ("g_jet_trk_pt_%s_syst_iCent%i_%s",      dir.Data (), iCent, var.Data ()));
          g_jet_trk_pt_bkg_syst[iDir][iCent][iVar]->Write  (Form ("g_jet_trk_pt_%s_bkg_syst_iCent%i_%s",  dir.Data (), iCent, var.Data ()));
          g_jet_trk_pt_sig_syst[iDir][iCent][iVar]->Write  (Form ("g_jet_trk_pt_%s_sig_syst_iCent%i_%s",  dir.Data (), iCent, var.Data ()));
          g_jet_trk_pt_iaa_syst[iDir][iCent][iVar]->Write  (Form ("g_jet_trk_pt_%s_iaa_syst_iCent%i_%s",  dir.Data (), iCent, var.Data ()));

        } // end loop over iCent

      } // end loop over iDir

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        const TString ptch = pTChSelections[iPtCh].Data ();

        g_jet_trk_dphi_ref_syst[iPtCh][iVar]->Write      (Form ("g_jet_trk_dphi_%s_ref_syst_%s",      ptch.Data (), var.Data ()));
        g_jet_trk_dphi_ref_bkg_syst[iPtCh][iVar]->Write  (Form ("g_jet_trk_dphi_%s_ref_bkg_syst_%s",  ptch.Data (), var.Data ()));
        g_jet_trk_dphi_ref_sig_syst[iPtCh][iVar]->Write  (Form ("g_jet_trk_dphi_%s_ref_sig_syst_%s",  ptch.Data (), var.Data ()));

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          g_jet_trk_dphi_syst[iPtCh][iCent][iVar]->Write     (Form ("g_jet_trk_dphi_%s_syst_iCent%i_%s",      ptch.Data (), iCent, var.Data ()));
          g_jet_trk_dphi_bkg_syst[iPtCh][iCent][iVar]->Write (Form ("g_jet_trk_dphi_%s_bkg_syst_iCent%i_%s",  ptch.Data (), iCent, var.Data ()));
          g_jet_trk_dphi_sig_syst[iPtCh][iCent][iVar]->Write (Form ("g_jet_trk_dphi_%s_sig_syst_iCent%i_%s",  ptch.Data (), iCent, var.Data ()));
          g_jet_trk_dphi_iaa_syst[iPtCh][iCent][iVar]->Write (Form ("g_jet_trk_dphi_%s_iaa_syst_iCent%i_%s",  ptch.Data (), iCent, var.Data ()));

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iVar


    for (int iTotVar = 0; iTotVar < 3; iTotVar++) {

      const TString totVar = totalVariations[iTotVar];

      for (int iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        g_jet_trk_pt_ref_systTot[iDir][iTotVar]->Write     (Form ("g_jet_trk_pt_%s_ref_%s_systTot",      dir.Data (), totVar.Data ()));
        g_jet_trk_pt_ref_bkg_systTot[iDir][iTotVar]->Write (Form ("g_jet_trk_pt_%s_ref_bkg_%s_systTot",  dir.Data (), totVar.Data ()));
        g_jet_trk_pt_ref_sig_systTot[iDir][iTotVar]->Write (Form ("g_jet_trk_pt_%s_ref_sig_%s_systTot",  dir.Data (), totVar.Data ()));

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          g_jet_trk_pt_systTot[iDir][iCent][iTotVar]->Write      (Form ("g_jet_trk_pt_%s_%s_systTot_iCent%i",      dir.Data (), totVar.Data (), iCent));
          g_jet_trk_pt_bkg_systTot[iDir][iCent][iTotVar]->Write  (Form ("g_jet_trk_pt_%s_bkg_%s_systTot_iCent%i",  dir.Data (), totVar.Data (), iCent));
          g_jet_trk_pt_sig_systTot[iDir][iCent][iTotVar]->Write  (Form ("g_jet_trk_pt_%s_sig_%s_systTot_iCent%i",  dir.Data (), totVar.Data (), iCent));
          g_jet_trk_pt_iaa_systTot[iDir][iCent][iTotVar]->Write  (Form ("g_jet_trk_pt_%s_iaa_%s_systTot_iCent%i",  dir.Data (), totVar.Data (), iCent));

        } // end loop over iCent

      } // end loop over iDir

      for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        const TString ptch = pTChSelections[iPtCh].Data ();

        g_jet_trk_dphi_ref_systTot[iPtCh][iTotVar]->Write      (Form ("g_jet_trk_dphi_%s_ref_%s_systTot",      ptch.Data (), totVar.Data ()));
        g_jet_trk_dphi_ref_bkg_systTot[iPtCh][iTotVar]->Write  (Form ("g_jet_trk_dphi_%s_ref_bkg_%s_systTot",  ptch.Data (), totVar.Data ()));
        g_jet_trk_dphi_ref_sig_systTot[iPtCh][iTotVar]->Write  (Form ("g_jet_trk_dphi_%s_ref_sig_%s_systTot",  ptch.Data (), totVar.Data ()));

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          g_jet_trk_dphi_systTot[iPtCh][iCent][iTotVar]->Write     (Form ("g_jet_trk_dphi_%s_%s_systTot_iCent%i",      ptch.Data (), totVar.Data (), iCent));
          g_jet_trk_dphi_bkg_systTot[iPtCh][iCent][iTotVar]->Write (Form ("g_jet_trk_dphi_%s_bkg_%s_systTot_iCent%i",  ptch.Data (), totVar.Data (), iCent));
          g_jet_trk_dphi_sig_systTot[iPtCh][iCent][iTotVar]->Write (Form ("g_jet_trk_dphi_%s_sig_%s_systTot_iCent%i",  ptch.Data (), totVar.Data (), iCent));
          g_jet_trk_dphi_iaa_systTot[iPtCh][iCent][iTotVar]->Write (Form ("g_jet_trk_dphi_%s_iaa_%s_systTot_iCent%i",  ptch.Data (), totVar.Data (), iCent));

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iTotVar


    for (int iDType = 0; iDType < 2; iDType++) {

      for (int iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;

        h_evt_counts_ref[iDType][iVar]->Write ();
        h_jet_counts_ref[iDType][iVar]->Write ();

        h_evt_counts_ref_bkg[iDType][iVar]->Write ();
        h_jet_counts_ref_bkg[iDType][iVar]->Write ();

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          h_evt_counts[iDType][iCent][iVar]->Write ();
          h_jet_counts[iDType][iCent][iVar]->Write ();

        } // end loop over iCent


        for (int iDir = 0; iDir < nDir; iDir++) {

          h_jet_trk_pt_ref[iDType][iDir][iVar]->Write ();
          h_jet_trk_pt_ref_bkg[iDType][iDir][iVar]->Write ();
          h_jet_trk_pt_ref_sig[iDType][iDir][iVar]->Write ();

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            h_jet_trk_pt[iDType][iDir][iCent][iVar]->Write ();
            h_jet_trk_pt_bkg[iDType][iDir][iCent][iVar]->Write ();
            h_jet_trk_pt_sig[iDType][iDir][iCent][iVar]->Write ();
            h_jet_trk_pt_iaa[iDType][iDir][iCent][iVar]->Write ();

          } // end loop over iCent

        } // end loop over iDir

        for (int iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          h_jet_trk_dphi_ref[iDType][iPtCh][iVar]->Write ();
          h_jet_trk_dphi_ref_bkg[iDType][iPtCh][iVar]->Write ();
          h_jet_trk_dphi_ref_sig[iDType][iPtCh][iVar]->Write ();

          for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

            h_jet_trk_dphi[iDType][iPtCh][iCent][iVar]->Write ();
            h_jet_trk_dphi_bkg[iDType][iPtCh][iCent][iVar]->Write ();
            h_jet_trk_dphi_sig[iDType][iPtCh][iCent][iVar]->Write ();
            h_jet_trk_dphi_iaa[iDType][iPtCh][iCent][iVar]->Write ();

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iVar

    } // end loop over iDType


    for (int iDir = 0; iDir < nDir; iDir++) {

      h_jet_trk_pt_ref_bbb[iDir]->Write ();
      f_jet_trk_pt_ref_bbb[iDir]->Write ();

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_jet_trk_pt_bbb[iDir][iCent]->Write ();
        f_jet_trk_pt_bbb[iDir][iCent]->Write ();

      } // end loop over iCent

    } // end loop over iDir


    outFile->Close ();
  }
}


#endif
