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

  TH1D***  h_evt_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_jet_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_evt_counts_ref_bkg = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_jet_counts_ref_bkg = Get2DArray <TH1D*> (2, nVar);
  TH1D**** h_evt_counts         = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_counts         = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_evt_counts_bkg     = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_counts_bkg     = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);

  TH1D***  h_jet_trk_pt_ns_ref        = Get2DArray <TH1D*> (2, nVar);
  TH2D***  h2_jet_trk_pt_ns_cov_ref   = Get2DArray <TH2D*> (2, nVar);
  TH1D***  h_jet_trk_pt_perp_ref      = Get2DArray <TH1D*> (2, nVar);
  TH2D***  h2_jet_trk_pt_perp_cov_ref = Get2DArray <TH2D*> (2, nVar);
  TH1D***  h_jet_trk_pt_as_ref        = Get2DArray <TH1D*> (2, nVar);
  TH2D***  h2_jet_trk_pt_as_cov_ref   = Get2DArray <TH2D*> (2, nVar);

  TH1D***  h_jet_trk_pt_ns_ref_bkg        = Get2DArray <TH1D*> (2, nVar);
  TH2D***  h2_jet_trk_pt_ns_cov_ref_bkg   = Get2DArray <TH2D*> (2, nVar);
  TH1D***  h_jet_trk_pt_perp_ref_bkg      = Get2DArray <TH1D*> (2, nVar);
  TH2D***  h2_jet_trk_pt_perp_cov_ref_bkg = Get2DArray <TH2D*> (2, nVar);
  TH1D***  h_jet_trk_pt_as_ref_bkg        = Get2DArray <TH1D*> (2, nVar);
  TH2D***  h2_jet_trk_pt_as_cov_ref_bkg   = Get2DArray <TH2D*> (2, nVar);

  TH1D**** h_jet_trk_pt_ns        = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH2D**** h2_jet_trk_pt_ns_cov   = Get3DArray <TH2D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_perp      = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH2D**** h2_jet_trk_pt_perp_cov = Get3DArray <TH2D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_as        = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH2D**** h2_jet_trk_pt_as_cov   = Get3DArray <TH2D*> (2, nZdcCentBins, nVar);

  TH1D**** h_jet_trk_pt_ns_bkg        = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH2D**** h2_jet_trk_pt_ns_cov_bkg   = Get3DArray <TH2D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_perp_bkg      = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH2D**** h2_jet_trk_pt_perp_cov_bkg = Get3DArray <TH2D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_as_bkg        = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH2D**** h2_jet_trk_pt_as_cov_bkg   = Get3DArray <TH2D*> (2, nZdcCentBins, nVar);

  TH1D***  h_jet_trk_pt_ns_ref_sig    = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_jet_trk_pt_perp_ref_sig  = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_jet_trk_pt_as_ref_sig    = Get2DArray <TH1D*> (2, nVar);

  TH1D**** h_jet_trk_pt_ns_sig    = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_perp_sig  = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_as_sig    = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);

  TH1D**** h_jet_trk_pt_ns_iaa    = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_perp_iaa  = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_trk_pt_as_iaa    = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);


  TGAE**  g_jet_trk_pt_ns_ref_syst        = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_perp_ref_syst      = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_as_ref_syst        = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_ns_ref_bkg_syst    = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_perp_ref_bkg_syst  = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_as_ref_bkg_syst    = Get1DArray <TGAE*> (nVar);
  TGAE*** g_jet_trk_pt_ns_syst            = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_perp_syst          = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_as_syst            = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_ns_bkg_syst        = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_perp_bkg_syst      = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_as_bkg_syst        = Get2DArray <TGAE*> (nZdcCentBins, nVar);

  TGAE**  g_jet_trk_pt_ns_ref_sig_syst    = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_perp_ref_sig_syst  = Get1DArray <TGAE*> (nVar);
  TGAE**  g_jet_trk_pt_as_ref_sig_syst    = Get1DArray <TGAE*> (nVar);
  TGAE*** g_jet_trk_pt_ns_sig_syst        = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_perp_sig_syst      = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_as_sig_syst        = Get2DArray <TGAE*> (nZdcCentBins, nVar);

  TGAE*** g_jet_trk_pt_ns_iaa_syst    = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_perp_iaa_syst  = Get2DArray <TGAE*> (nZdcCentBins, nVar);
  TGAE*** g_jet_trk_pt_as_iaa_syst    = Get2DArray <TGAE*> (nZdcCentBins, nVar);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/PlotPtCh_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");

  for (int iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (int iVar = 0; iVar < nVar; iVar++) {

      if ((iDType == 0 && dataVariations.count (variations[iVar]) == 0) || (iDType == 1 && mcVariations.count (variations[iVar]) == 0))
        continue;

      {
        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s17_5TeV_hists.root", rootPath.Data (), tag1, variations[iVar].Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts_ref[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h_evt_counts_ref_%s_%s", dType.Data (), variations[iVar].Data ()));
        h_jet_counts_ref[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h_jet_counts_ref_%s_%s", dType.Data (), variations[iVar].Data ()));

        h_jet_trk_pt_ns_ref[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h_jet_trk_pt_ns_ref_%s_%s", dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_ns_cov_ref[iDType][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_ns_cov_ref_%s_%s", dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_perp_ref[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h_jet_trk_pt_perp_ref_%s_%s", dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_perp_cov_ref[iDType][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_perp_cov_ref_%s_%s", dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_as_ref[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h_jet_trk_pt_as_ref_%s_%s", dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_as_cov_ref[iDType][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_%s17", tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_as_cov_ref_%s_%s", dType.Data (), variations[iVar].Data ()));

        inFile->Close ();

        CalcUncertainties (h_jet_trk_pt_ns_ref[iDType][iVar],   h2_jet_trk_pt_ns_cov_ref[iDType][iVar],   h_evt_counts_ref[iDType][iVar]);
        CalcUncertainties (h_jet_trk_pt_perp_ref[iDType][iVar], h2_jet_trk_pt_perp_cov_ref[iDType][iVar], h_evt_counts_ref[iDType][iVar]);
        CalcUncertainties (h_jet_trk_pt_as_ref[iDType][iVar],   h2_jet_trk_pt_as_cov_ref[iDType][iVar],   h_evt_counts_ref[iDType][iVar]);
      }



      {
        TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/data17_5TeV_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", variations[iVar].Data ());
        //TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/%s17_5TeV_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", variations[iVar].Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts_ref_bkg[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_evt_counts_ref_bkg_%s_%s", dType.Data (), variations[iVar].Data ()));
        h_jet_counts_ref_bkg[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_counts_ref_bkg_%s_%s", dType.Data (), variations[iVar].Data ()));
        //h_evt_counts_ref_bkg[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_evt_counts_ref_bkg_%s_%s", dType.Data (), variations[iVar].Data ()));
        //h_jet_counts_ref_bkg[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s17", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_counts_ref_bkg_%s_%s", dType.Data (), variations[iVar].Data ()));

        h_jet_trk_pt_ns_ref_bkg[iDType][iVar]         = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data17",         doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_ns_ref_bkg_%s_%s",         dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_ns_cov_ref_bkg[iDType][iVar]    = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data17",    doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_ns_cov_ref_bkg_%s_%s",    dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_perp_ref_bkg[iDType][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_data17",       doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_perp_ref_bkg_%s_%s",       dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_perp_cov_ref_bkg[iDType][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_data17",  doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_perp_cov_ref_bkg_%s_%s",  dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_as_ref_bkg[iDType][iVar]         = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data17",         doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_as_ref_bkg_%s_%s",         dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_as_cov_ref_bkg[iDType][iVar]    = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data17",    doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_as_cov_ref_bkg_%s_%s",    dType.Data (), variations[iVar].Data ()));
        //h_jet_trk_pt_ns_ref_bkg[iDType][iVar]         = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_%s17",         doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_pt_ns_ref_bkg_%s_%s",         dType.Data (), variations[iVar].Data ()));
        //h2_jet_trk_pt_ns_cov_ref_bkg[iDType][iVar]    = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_%s17",    doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_ns_cov_ref_bkg_%s_%s",    dType.Data (), variations[iVar].Data ()));
        //h_jet_trk_pt_perp_ref_bkg[iDType][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_%s17",       doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_pt_perp_ref_bkg_%s_%s",       dType.Data (), variations[iVar].Data ()));
        //h2_jet_trk_pt_perp_cov_ref_bkg[iDType][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_%s17",  doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_perp_cov_ref_bkg_%s_%s",  dType.Data (), variations[iVar].Data ()));
        //h_jet_trk_pt_as_ref_bkg[iDType][iVar]         = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_%s17",         doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_pt_as_ref_bkg_%s_%s",         dType.Data (), variations[iVar].Data ()));
        //h2_jet_trk_pt_as_cov_ref_bkg[iDType][iVar]    = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_%s17",    doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_as_cov_ref_bkg_%s_%s",    dType.Data (), variations[iVar].Data ()));

        inFile->Close ();

        CalcUncertainties (h_jet_trk_pt_ns_ref_bkg[iDType][iVar],   h2_jet_trk_pt_ns_cov_ref_bkg[iDType][iVar],   h_evt_counts_ref_bkg[iDType][iVar]);
        CalcUncertainties (h_jet_trk_pt_perp_ref_bkg[iDType][iVar], h2_jet_trk_pt_perp_cov_ref_bkg[iDType][iVar], h_evt_counts_ref_bkg[iDType][iVar]);
        CalcUncertainties (h_jet_trk_pt_as_ref_bkg[iDType][iVar],   h2_jet_trk_pt_as_cov_ref_bkg[iDType][iVar],   h_evt_counts_ref_bkg[iDType][iVar]);
      }



      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s16_5TeV_iCent%i_hists.root", rootPath.Data (), tag1, variations[iVar].Data (), dType.Data (), iCent);
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s16", tag1, dType.Data ()))->Clone (Form ("h_evt_counts_pPb_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
        h_jet_counts[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s16", tag1, dType.Data ()))->Clone (Form ("h_jet_counts_pPb_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));

        h_jet_trk_pt_ns[iDType][iCent][iVar]        = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_%s16",         tag1, dType.Data ()))->Clone (Form ("h_jet_trk_pt_ns_pPb_iCent%i_%s_%s",        iCent, dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_ns_cov[iDType][iCent][iVar]   = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_%s16",    tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_ns_cov_pPb_iCent%i_%s_%s",   iCent, dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_perp[iDType][iCent][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_%s16",       tag1, dType.Data ()))->Clone (Form ("h_jet_trk_pt_perp_pPb_iCent%i_%s_%s",      iCent, dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_perp_cov[iDType][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_%s16",  tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_perp_cov_pPb_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_as[iDType][iCent][iVar]        = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_%s16",         tag1, dType.Data ()))->Clone (Form ("h_jet_trk_pt_as_pPb_iCent%i_%s_%s",        iCent, dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_as_cov[iDType][iCent][iVar]   = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_%s16",    tag1, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_as_cov_pPb_iCent%i_%s_%s",   iCent, dType.Data (), variations[iVar].Data ()));

        inFile->Close ();

        CalcUncertainties (h_jet_trk_pt_ns[iDType][iCent][iVar],    h2_jet_trk_pt_ns_cov[iDType][iCent][iVar],    h_evt_counts[iDType][iCent][iVar]);
        CalcUncertainties (h_jet_trk_pt_perp[iDType][iCent][iVar],  h2_jet_trk_pt_perp_cov[iDType][iCent][iVar],  h_evt_counts[iDType][iCent][iVar]);
        CalcUncertainties (h_jet_trk_pt_as[iDType][iCent][iVar],    h2_jet_trk_pt_as_cov[iDType][iCent][iVar],    h_evt_counts[iDType][iCent][iVar]);
      } // end loop over iCent



      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
        //TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", variations[iVar].Data (), iCent); // for now just use the data background for MC
        TString inFileName = Form ("%s/Histograms/%s/%sHists/%s/%s16_5TeV_iCent%i_hists.root", rootPath.Data (), doMix ? tag1 : tag2, doMix ? "Mixed" : "Jets", variations[iVar].Data (), dType.Data (), iCent);
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        //h_evt_counts_bkg[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_evt_counts_pPb_bkg_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
        //h_jet_counts_bkg[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_counts_pPb_bkg_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
        h_evt_counts_bkg[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_evt_counts_pPb_bkg_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
        h_jet_counts_bkg[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s16", doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_counts_pPb_bkg_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));

        //h_jet_trk_pt_ns_bkg[iDType][iCent][iVar]        = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_data16",         doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_ns_pPb_bkg_iCent%i_%s_%s",        iCent, dType.Data (), variations[iVar].Data ()));
        //h2_jet_trk_pt_ns_cov_bkg[iDType][iCent][iVar]   = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_data16",    doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_ns_cov_pPb_bkg_iCent%i_%s_%s",   iCent, dType.Data (), variations[iVar].Data ()));
        //h_jet_trk_pt_perp_bkg[iDType][iCent][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_data16",       doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_perp_pPb_bkg_iCent%i_%s_%s",      iCent, dType.Data (), variations[iVar].Data ()));
        //h2_jet_trk_pt_perp_cov_bkg[iDType][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_data16",  doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_perp_cov_pPb_bkg_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
        //h_jet_trk_pt_as_bkg[iDType][iCent][iVar]        = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_data16",         doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h_jet_trk_pt_as_pPb_bkg_iCent%i_%s_%s",        iCent, dType.Data (), variations[iVar].Data ()));
        //h2_jet_trk_pt_as_cov_bkg[iDType][iCent][iVar]   = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_data16",    doMix ? (std::string (tag1) + "_mixed").c_str () : tag2))->Clone (Form ("h2_jet_trk_pt_as_cov_pPb_bkg_iCent%i_%s_%s",   iCent, dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_ns_bkg[iDType][iCent][iVar]        = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_ns_%s_%s16",         doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_pt_ns_pPb_bkg_iCent%i_%s_%s",         iCent, dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_ns_cov_bkg[iDType][iCent][iVar]   = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_ns_cov_%s_%s16",    doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_ns_cov_pPb_bkg_iCent%i_%s_%s",    iCent, dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_perp_bkg[iDType][iCent][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_perp_%s_%s16",       doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_pt_perp_pPb_bkg_iCent%i_%s_%s",       iCent, dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_perp_cov_bkg[iDType][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_perp_cov_%s_%s16",  doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_perp_cov_pPb_bkg_iCent%i_%s_%s",  iCent, dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_as_bkg[iDType][iCent][iVar]        = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_as_%s_%s16",         doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h_jet_trk_pt_as_pPb_bkg_iCent%i_%s_%s",         iCent, dType.Data (), variations[iVar].Data ()));
        h2_jet_trk_pt_as_cov_bkg[iDType][iCent][iVar]   = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_as_cov_%s_%s16",    doMix ? (std::string (tag1) + "_mixed").c_str () : tag2, dType.Data ()))->Clone (Form ("h2_jet_trk_pt_as_cov_pPb_bkg_iCent%i_%s_%s",    iCent, dType.Data (), variations[iVar].Data ()));

        inFile->Close ();

        CalcUncertainties (h_jet_trk_pt_ns_bkg[iDType][iCent][iVar],    h2_jet_trk_pt_ns_cov_bkg[iDType][iCent][iVar],    h_evt_counts_bkg[iDType][iCent][iVar]);
        CalcUncertainties (h_jet_trk_pt_perp_bkg[iDType][iCent][iVar],  h2_jet_trk_pt_perp_cov_bkg[iDType][iCent][iVar],  h_evt_counts_bkg[iDType][iCent][iVar]);
        CalcUncertainties (h_jet_trk_pt_as_bkg[iDType][iCent][iVar],    h2_jet_trk_pt_as_cov_bkg[iDType][iCent][iVar],    h_evt_counts_bkg[iDType][iCent][iVar]);
      } // end loop over iCent



      {
        h_jet_trk_pt_ns_ref_sig[iDType][iVar] = (TH1D*) h_jet_trk_pt_ns_ref[iDType][iVar]->Clone (Form ("h_jet_trk_pt_ns_ref_sig_%s_%s", dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_ns_ref_sig[iDType][iVar]->Add (h_jet_trk_pt_ns_ref_bkg[iDType][iVar], -1);

        h_jet_trk_pt_perp_ref_sig[iDType][iVar] = (TH1D*) h_jet_trk_pt_perp_ref[iDType][iVar]->Clone (Form ("h_jet_trk_pt_perp_ref_sig_%s_%s", dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_perp_ref_sig[iDType][iVar]->Add (h_jet_trk_pt_perp_ref_bkg[iDType][iVar], -1);

        h_jet_trk_pt_as_ref_sig[iDType][iVar] = (TH1D*) h_jet_trk_pt_as_ref[iDType][iVar]->Clone (Form ("h_jet_trk_pt_as_ref_sig_%s_%s", dType.Data (), variations[iVar].Data ()));
        h_jet_trk_pt_as_ref_sig[iDType][iVar]->Add (h_jet_trk_pt_as_ref_bkg[iDType][iVar], -1);

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          h_jet_trk_pt_ns_sig[iDType][iCent][iVar] = (TH1D*) h_jet_trk_pt_ns[iDType][iCent][iVar]->Clone (Form ("h_jet_trk_pt_ns_pPb_sig_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
          h_jet_trk_pt_ns_sig[iDType][iCent][iVar]->Add (h_jet_trk_pt_ns_bkg[iDType][iCent][iVar], -1);

          h_jet_trk_pt_perp_sig[iDType][iCent][iVar] = (TH1D*) h_jet_trk_pt_perp[iDType][iCent][iVar]->Clone (Form ("h_jet_trk_pt_perp_pPb_sig_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
          h_jet_trk_pt_perp_sig[iDType][iCent][iVar]->Add (h_jet_trk_pt_perp_bkg[iDType][iCent][iVar], -1);

          h_jet_trk_pt_as_sig[iDType][iCent][iVar] = (TH1D*) h_jet_trk_pt_as[iDType][iCent][iVar]->Clone (Form ("h_jet_trk_pt_as_pPb_sig_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
          h_jet_trk_pt_as_sig[iDType][iCent][iVar]->Add (h_jet_trk_pt_as_bkg[iDType][iCent][iVar], -1);


          h_jet_trk_pt_ns_iaa[iDType][iCent][iVar] = (TH1D*) h_jet_trk_pt_ns_sig[iDType][iCent][iVar]->Clone (Form ("h_jet_trk_pt_ns_iaa_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
          h_jet_trk_pt_ns_iaa[iDType][iCent][iVar]->Divide (h_jet_trk_pt_ns_ref_sig[iDType][iVar]);

          h_jet_trk_pt_perp_iaa[iDType][iCent][iVar] = (TH1D*) h_jet_trk_pt_perp_sig[iDType][iCent][iVar]->Clone (Form ("h_jet_trk_pt_perp_iaa_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
          h_jet_trk_pt_perp_iaa[iDType][iCent][iVar]->Divide (h_jet_trk_pt_perp_ref_sig[iDType][iVar]);

          h_jet_trk_pt_as_iaa[iDType][iCent][iVar] = (TH1D*) h_jet_trk_pt_as_sig[iDType][iCent][iVar]->Clone (Form ("h_jet_trk_pt_as_iaa_iCent%i_%s_%s", iCent, dType.Data (), variations[iVar].Data ()));
          h_jet_trk_pt_as_iaa[iDType][iCent][iVar]->Divide (h_jet_trk_pt_as_ref_sig[iDType][iVar]);

        } // end loop over iCent
      }
    } // end loop over iVar
  } // end loop over iDType



  {
    g_jet_trk_pt_ns_ref_syst[0]       = make_graph (h_jet_trk_pt_ns_ref[0][0]);
    g_jet_trk_pt_perp_ref_syst[0]     = make_graph (h_jet_trk_pt_perp_ref[0][0]);
    g_jet_trk_pt_as_ref_syst[0]       = make_graph (h_jet_trk_pt_as_ref[0][0]);

    g_jet_trk_pt_ns_ref_bkg_syst[0]   = make_graph (h_jet_trk_pt_ns_ref_bkg[0][0]);
    g_jet_trk_pt_perp_ref_bkg_syst[0] = make_graph (h_jet_trk_pt_perp_ref_bkg[0][0]);
    g_jet_trk_pt_as_ref_bkg_syst[0]   = make_graph (h_jet_trk_pt_as_ref_bkg[0][0]);

    g_jet_trk_pt_ns_ref_sig_syst[0]   = make_graph (h_jet_trk_pt_ns_ref_sig[0][0]);
    g_jet_trk_pt_perp_ref_sig_syst[0] = make_graph (h_jet_trk_pt_perp_ref_sig[0][0]);
    g_jet_trk_pt_as_ref_sig_syst[0]   = make_graph (h_jet_trk_pt_as_ref_sig[0][0]);

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

      g_jet_trk_pt_ns_syst[iCent][0]        = make_graph (h_jet_trk_pt_ns[0][iCent][0]);
      g_jet_trk_pt_perp_syst[iCent][0]      = make_graph (h_jet_trk_pt_perp[0][iCent][0]);
      g_jet_trk_pt_as_syst[iCent][0]        = make_graph (h_jet_trk_pt_as[0][iCent][0]);

      g_jet_trk_pt_ns_bkg_syst[iCent][0]    = make_graph (h_jet_trk_pt_ns_bkg[0][iCent][0]);
      g_jet_trk_pt_perp_bkg_syst[iCent][0]  = make_graph (h_jet_trk_pt_perp_bkg[0][iCent][0]);
      g_jet_trk_pt_as_bkg_syst[iCent][0]    = make_graph (h_jet_trk_pt_as_bkg[0][iCent][0]);

      g_jet_trk_pt_ns_sig_syst[iCent][0]    = make_graph (h_jet_trk_pt_ns_sig[0][iCent][0]);
      g_jet_trk_pt_perp_sig_syst[iCent][0]  = make_graph (h_jet_trk_pt_perp_sig[0][iCent][0]);
      g_jet_trk_pt_as_sig_syst[iCent][0]    = make_graph (h_jet_trk_pt_as_sig[0][iCent][0]);

      g_jet_trk_pt_ns_iaa_syst[iCent][0]    = make_graph (h_jet_trk_pt_ns_iaa[0][iCent][0]);
      g_jet_trk_pt_perp_iaa_syst[iCent][0]  = make_graph (h_jet_trk_pt_perp_iaa[0][iCent][0]);
      g_jet_trk_pt_as_iaa_syst[iCent][0]    = make_graph (h_jet_trk_pt_as_iaa[0][iCent][0]);

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



  for (int iDType = 0; iDType < 2; iDType++) {

    for (int iVar = 1; iVar < nVar; iVar++) {

      if ((iDType == 0 && dataVariations.count (variations[iVar]) == 0) || (iDType == 1 && mcVariations.count (variations[iVar]) == 0))
        continue;

      g_jet_trk_pt_ns_ref_syst[iVar]        = new TGAE ();
      g_jet_trk_pt_perp_ref_syst[iVar]      = new TGAE ();
      g_jet_trk_pt_as_ref_syst[iVar]        = new TGAE ();

      g_jet_trk_pt_ns_ref_bkg_syst[iVar]    = new TGAE ();
      g_jet_trk_pt_perp_ref_bkg_syst[iVar]  = new TGAE ();
      g_jet_trk_pt_as_ref_bkg_syst[iVar]    = new TGAE ();

      g_jet_trk_pt_ns_ref_sig_syst[iVar]    = new TGAE ();
      g_jet_trk_pt_perp_ref_sig_syst[iVar]  = new TGAE ();
      g_jet_trk_pt_as_ref_sig_syst[iVar]    = new TGAE ();

      CalcSystematics (g_jet_trk_pt_ns_ref_syst[iVar],    h_jet_trk_pt_ns_ref[iDType][0],    h_jet_trk_pt_ns_ref[iDType][iVar]);
      CalcSystematics (g_jet_trk_pt_perp_ref_syst[iVar],  h_jet_trk_pt_perp_ref[iDType][0],  h_jet_trk_pt_perp_ref[iDType][iVar]);
      CalcSystematics (g_jet_trk_pt_as_ref_syst[iVar],    h_jet_trk_pt_as_ref[iDType][0],    h_jet_trk_pt_as_ref[iDType][iVar]);

      CalcSystematics (g_jet_trk_pt_ns_ref_bkg_syst[iVar],    h_jet_trk_pt_ns_ref_bkg[iDType][0],    h_jet_trk_pt_ns_ref_bkg[iDType][iVar]);
      CalcSystematics (g_jet_trk_pt_perp_ref_bkg_syst[iVar],  h_jet_trk_pt_perp_ref_bkg[iDType][0],  h_jet_trk_pt_perp_ref_bkg[iDType][iVar]);
      CalcSystematics (g_jet_trk_pt_as_ref_bkg_syst[iVar],    h_jet_trk_pt_as_ref_bkg[iDType][0],    h_jet_trk_pt_as_ref_bkg[iDType][iVar]);

      CalcSystematics (g_jet_trk_pt_ns_ref_sig_syst[iVar],    h_jet_trk_pt_ns_ref_sig[iDType][0],    h_jet_trk_pt_ns_ref_sig[iDType][iVar]);
      CalcSystematics (g_jet_trk_pt_perp_ref_sig_syst[iVar],  h_jet_trk_pt_perp_ref_sig[iDType][0],  h_jet_trk_pt_perp_ref_sig[iDType][iVar]);
      CalcSystematics (g_jet_trk_pt_as_ref_sig_syst[iVar],    h_jet_trk_pt_as_ref_sig[iDType][0],    h_jet_trk_pt_as_ref_sig[iDType][iVar]);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        g_jet_trk_pt_ns_syst[iCent][iVar]       = new TGAE ();
        g_jet_trk_pt_perp_syst[iCent][iVar]     = new TGAE ();
        g_jet_trk_pt_as_syst[iCent][iVar]       = new TGAE ();

        g_jet_trk_pt_ns_bkg_syst[iCent][iVar]   = new TGAE ();
        g_jet_trk_pt_perp_bkg_syst[iCent][iVar] = new TGAE ();
        g_jet_trk_pt_as_bkg_syst[iCent][iVar]   = new TGAE ();

        g_jet_trk_pt_ns_sig_syst[iCent][iVar]   = new TGAE ();
        g_jet_trk_pt_perp_sig_syst[iCent][iVar] = new TGAE ();
        g_jet_trk_pt_as_sig_syst[iCent][iVar]   = new TGAE ();

        g_jet_trk_pt_ns_iaa_syst[iCent][iVar]   = new TGAE ();
        g_jet_trk_pt_perp_iaa_syst[iCent][iVar] = new TGAE ();
        g_jet_trk_pt_as_iaa_syst[iCent][iVar]   = new TGAE ();

        CalcSystematics (g_jet_trk_pt_ns_syst[iCent][iVar],       h_jet_trk_pt_ns[iDType][iCent][0],       h_jet_trk_pt_ns[iDType][iCent][iVar]);
        CalcSystematics (g_jet_trk_pt_perp_syst[iCent][iVar],     h_jet_trk_pt_perp[iDType][iCent][0],     h_jet_trk_pt_perp[iDType][iCent][iVar]);
        CalcSystematics (g_jet_trk_pt_as_syst[iCent][iVar],       h_jet_trk_pt_as[iDType][iCent][0],       h_jet_trk_pt_as[iDType][iCent][iVar]);

        CalcSystematics (g_jet_trk_pt_ns_bkg_syst[iCent][iVar],   h_jet_trk_pt_ns_bkg[iDType][iCent][0],   h_jet_trk_pt_ns_bkg[iDType][iCent][iVar]);
        CalcSystematics (g_jet_trk_pt_perp_bkg_syst[iCent][iVar], h_jet_trk_pt_perp_bkg[iDType][iCent][0], h_jet_trk_pt_perp_bkg[iDType][iCent][iVar]);
        CalcSystematics (g_jet_trk_pt_as_bkg_syst[iCent][iVar],   h_jet_trk_pt_as_bkg[iDType][iCent][0],   h_jet_trk_pt_as_bkg[iDType][iCent][iVar]);

        CalcSystematics (g_jet_trk_pt_ns_sig_syst[iCent][iVar],   h_jet_trk_pt_ns_sig[iDType][iCent][0],   h_jet_trk_pt_ns_sig[iDType][iCent][iVar]);
        CalcSystematics (g_jet_trk_pt_perp_sig_syst[iCent][iVar], h_jet_trk_pt_perp_sig[iDType][iCent][0], h_jet_trk_pt_perp_sig[iDType][iCent][iVar]);
        CalcSystematics (g_jet_trk_pt_as_sig_syst[iCent][iVar],   h_jet_trk_pt_as_sig[iDType][iCent][0],   h_jet_trk_pt_as_sig[iDType][iCent][iVar]);

        CalcSystematics (g_jet_trk_pt_ns_iaa_syst[iCent][iVar],   h_jet_trk_pt_ns_iaa[iDType][iCent][0],   h_jet_trk_pt_ns_iaa[iDType][iCent][iVar]);
        CalcSystematics (g_jet_trk_pt_perp_iaa_syst[iCent][iVar], h_jet_trk_pt_perp_iaa[iDType][iCent][0], h_jet_trk_pt_perp_iaa[iDType][iCent][iVar]);
        CalcSystematics (g_jet_trk_pt_as_iaa_syst[iCent][iVar],   h_jet_trk_pt_as_iaa[iDType][iCent][0],   h_jet_trk_pt_as_iaa[iDType][iCent][iVar]);

      } // end loop over iCent

    } // end loop over iVar

  } // end loop over iDType



  for (int iVar = 1; iVar < nVar; iVar++) {

    if (dataVariations.count (variations[iVar]) > 0 || mcVariations.count (variations[iVar]) == 0)
      continue;

    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ns_ref_syst[iVar], h_jet_trk_pt_ns_ref[0][0]);
    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_perp_ref_syst[iVar], h_jet_trk_pt_perp_ref[0][0]);
    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_as_ref_syst[iVar], h_jet_trk_pt_as_ref[0][0]);

    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ns_ref_bkg_syst[iVar], h_jet_trk_pt_ns_ref_bkg[0][0]);
    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_perp_ref_bkg_syst[iVar], h_jet_trk_pt_perp_ref_bkg[0][0]);
    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_as_ref_bkg_syst[iVar], h_jet_trk_pt_as_ref_bkg[0][0]);

    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ns_ref_sig_syst[iVar], h_jet_trk_pt_ns_ref_sig[0][0]);
    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_perp_ref_sig_syst[iVar], h_jet_trk_pt_perp_ref_sig[0][0]);
    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_as_ref_sig_syst[iVar], h_jet_trk_pt_as_ref_sig[0][0]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ns_syst[iCent][iVar], h_jet_trk_pt_ns[0][iCent][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_perp_syst[iCent][iVar], h_jet_trk_pt_perp[0][iCent][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_as_syst[iCent][iVar], h_jet_trk_pt_as[0][iCent][0]);

      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ns_bkg_syst[iCent][iVar], h_jet_trk_pt_ns_bkg[0][iCent][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_perp_bkg_syst[iCent][iVar], h_jet_trk_pt_perp_bkg[0][iCent][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_as_bkg_syst[iCent][iVar], h_jet_trk_pt_as_bkg[0][iCent][0]);

      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ns_sig_syst[iCent][iVar], h_jet_trk_pt_ns_sig[0][iCent][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_perp_sig_syst[iCent][iVar], h_jet_trk_pt_perp_sig[0][iCent][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_as_sig_syst[iCent][iVar], h_jet_trk_pt_as_sig[0][iCent][0]);

      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ns_iaa_syst[iCent][iVar], h_jet_trk_pt_ns_iaa[0][iCent][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_perp_iaa_syst[iCent][iVar], h_jet_trk_pt_perp_iaa[0][iCent][0]);
      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_as_iaa_syst[iCent][iVar], h_jet_trk_pt_as_iaa[0][iCent][0]);

    } // end loop over iCent
    
  } // end loop over iVar



  for (int iVar = 1; iVar < nVar; iVar++) {

    if (dataVariations.count (variations[iVar]) == 0)
      continue;

    AddErrorsInQuadrature (g_jet_trk_pt_ns_ref_syst[0], g_jet_trk_pt_ns_ref_syst[iVar]);
    AddErrorsInQuadrature (g_jet_trk_pt_perp_ref_syst[0], g_jet_trk_pt_perp_ref_syst[iVar]);
    AddErrorsInQuadrature (g_jet_trk_pt_as_ref_syst[0], g_jet_trk_pt_as_ref_syst[iVar]);

    AddErrorsInQuadrature (g_jet_trk_pt_ns_ref_bkg_syst[0], g_jet_trk_pt_ns_ref_bkg_syst[iVar]);
    AddErrorsInQuadrature (g_jet_trk_pt_perp_ref_bkg_syst[0], g_jet_trk_pt_perp_ref_bkg_syst[iVar]);
    AddErrorsInQuadrature (g_jet_trk_pt_as_ref_bkg_syst[0], g_jet_trk_pt_as_ref_bkg_syst[iVar]);

    AddErrorsInQuadrature (g_jet_trk_pt_ns_ref_sig_syst[0], g_jet_trk_pt_ns_ref_sig_syst[iVar]);
    AddErrorsInQuadrature (g_jet_trk_pt_perp_ref_sig_syst[0], g_jet_trk_pt_perp_ref_sig_syst[iVar]);
    AddErrorsInQuadrature (g_jet_trk_pt_as_ref_sig_syst[0], g_jet_trk_pt_as_ref_sig_syst[iVar]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      AddErrorsInQuadrature (g_jet_trk_pt_ns_syst[iCent][0], g_jet_trk_pt_ns_syst[iCent][iVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_perp_syst[iCent][0], g_jet_trk_pt_perp_syst[iCent][iVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_as_syst[iCent][0], g_jet_trk_pt_as_syst[iCent][iVar]);
  
      AddErrorsInQuadrature (g_jet_trk_pt_ns_bkg_syst[iCent][0], g_jet_trk_pt_ns_bkg_syst[iCent][iVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_perp_bkg_syst[iCent][0], g_jet_trk_pt_perp_bkg_syst[iCent][iVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_as_bkg_syst[iCent][0], g_jet_trk_pt_as_bkg_syst[iCent][iVar]);
  
      AddErrorsInQuadrature (g_jet_trk_pt_ns_sig_syst[iCent][0], g_jet_trk_pt_ns_sig_syst[iCent][iVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_perp_sig_syst[iCent][0], g_jet_trk_pt_perp_sig_syst[iCent][iVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_as_sig_syst[iCent][0], g_jet_trk_pt_as_sig_syst[iCent][iVar]);

      AddErrorsInQuadrature (g_jet_trk_pt_ns_iaa_syst[iCent][0], g_jet_trk_pt_ns_iaa_syst[iCent][iVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_perp_iaa_syst[iCent][0], g_jet_trk_pt_perp_iaa_syst[iCent][iVar]);
      AddErrorsInQuadrature (g_jet_trk_pt_as_iaa_syst[iCent][0], g_jet_trk_pt_as_iaa_syst[iCent][iVar]);

    } // end loop over iCent

  } // end loop over iVar



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

    } // end loop over iCent

    for (int iDType = 0; iDType < 2; iDType++) {

      for (int iVar = 0; iVar < nVar; iVar++) {

        if ((iDType == 0 && dataVariations.count (variations[iVar]) == 0) || (iDType == 1 && mcVariations.count (variations[iVar]) == 0))
          continue;

        h_evt_counts_ref[iDType][iVar]->Write ();
        h_jet_counts_ref[iDType][iVar]->Write ();

        h_evt_counts_ref_bkg[iDType][iVar]->Write ();
        h_jet_counts_ref_bkg[iDType][iVar]->Write ();

        h_jet_trk_pt_ns_ref[iDType][iVar]->Write ();
        h_jet_trk_pt_perp_ref[iDType][iVar]->Write ();
        h_jet_trk_pt_as_ref[iDType][iVar]->Write ();

        h_jet_trk_pt_ns_ref_bkg[iDType][iVar]->Write ();
        h_jet_trk_pt_perp_ref_bkg[iDType][iVar]->Write ();
        h_jet_trk_pt_as_ref_bkg[iDType][iVar]->Write ();

        h_jet_trk_pt_ns_ref_sig[iDType][iVar]->Write ();
        h_jet_trk_pt_perp_ref_sig[iDType][iVar]->Write ();
        h_jet_trk_pt_as_ref_sig[iDType][iVar]->Write ();

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          h_evt_counts[iDType][iCent][iVar]->Write ();
          h_jet_counts[iDType][iCent][iVar]->Write ();

          h_jet_trk_pt_ns[iDType][iCent][iVar]->Write ();
          h_jet_trk_pt_perp[iDType][iCent][iVar]->Write ();
          h_jet_trk_pt_as[iDType][iCent][iVar]->Write ();

          h_jet_trk_pt_ns_bkg[iDType][iCent][iVar]->Write ();
          h_jet_trk_pt_perp_bkg[iDType][iCent][iVar]->Write ();
          h_jet_trk_pt_as_bkg[iDType][iCent][iVar]->Write ();

          h_jet_trk_pt_ns_sig[iDType][iCent][iVar]->Write ();
          h_jet_trk_pt_perp_sig[iDType][iCent][iVar]->Write ();
          h_jet_trk_pt_as_sig[iDType][iCent][iVar]->Write ();

          h_jet_trk_pt_ns_iaa[iDType][iCent][iVar]->Write ();
          h_jet_trk_pt_perp_iaa[iDType][iCent][iVar]->Write ();
          h_jet_trk_pt_as_iaa[iDType][iCent][iVar]->Write ();

        } // end loop over iCent

      } // end loop over iVar

    } // end loop over iDType

    outFile->Close ();
  }
}


#endif
