#ifndef __JetHadronCorrelator_ProcessUnfolding_C__
#define __JetHadronCorrelator_ProcessUnfolding_C__

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"
#include "PiecewisePolynomialConstantFunc.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TH1D.h>
#include <TH2D.h>

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

#include <iostream>
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;


short GetJetSpectraNIters (const short iPtJInt, const short iCent) {
  switch (iPtJInt) {
    case 0: // pTJ > 30 GeV
    switch (iCent) {
      case 5:  return 13; // all centralities
      case 4:  return 13; // 0-20%
      case 3:  return 11; // 20-40%...
      case 2:  return 11;
      case 1:  return 9;
      case 0:  return 9;  // ... 80-100%
      case -1: return 15; // pp
      default: return 0;
    }
    case 1: // pTJ > 60 GeV
    switch (iCent) {
      case 5:  return 17; // all centralities
      case 4:  return 15; // 0-20%
      case 3:  return 13; // 20-40% ...
      case 2:  return 13;
      case 1:  return 11;
      case 0:  return 9;  // ... 80-100%
      case -1: return 15; // pp
      default: return 0;
    }
    default: return 0;
  }
}


short GetTrkSpectraNIters (const short iPtJInt, const short iDir, const short iCent) {
  switch (iPtJInt) {
    case 0: // pTJ > 30 GeV
    switch (iDir) {
      case 0: // near-side
      switch (iCent) {
        case 5:  return 11; // all centralities
        case 4:  return 9; // 0-20%
        case 3:  return 9; // 20-40%...
        case 2:  return 9;
        case 1:  return 7;
        case 0:  return 7; // ... 80-100%
        case -1: return 7; // pp
        default: return 0;
      }
      case 1: // perpendicular
      switch (iCent) {
        case 5:  return 1;
        case 4:  return 1;
        case 3:  return 1;
        case 2:  return 1;
        case 1:  return 1;
        case 0:  return 1;
        case -1: return 1;
        default: return 0;
      }
      case 2: // away-side
      switch (iCent) {
        case 5:  return 9; // all centralities
        case 4:  return 9; // 0-20%
        case 3:  return 9; // 20-40%...
        case 2:  return 7;
        case 1:  return 7;
        case 0:  return 7; // ... 80-100%
        case -1: return 7; // pp
        default: return 0;
      }
      default: return 0;
    }
    case 1: // pTJ > 60 GeV
    switch (iDir) {
      case 0: // near-side
      switch (iCent) {
        case 5:  return 17; // all centralities
        case 4:  return 15; // 0-20%
        case 3:  return 15; // 20-40%...
        case 2:  return 13;
        case 1:  return 11;
        case 0:  return 9;  // ... 80-100%
        case -1: return 15; // pp
        default: return 0;
      }
      case 1: // perpendicular
      switch (iCent) {
        case 5:  return 1;
        case 4:  return 1;
        case 3:  return 1;
        case 2:  return 1;
        case 1:  return 1;
        case 0:  return 1;
        case -1: return 3;
        default: return 0;
      }
      case 2: // away-side
      switch (iCent) {
        case 5:  return 15; // all centralities
        case 4:  return 13; // 0-20%
        case 3:  return 13; // 20-40%...
        case 2:  return 11;
        case 1:  return 11; // "nominal" version
        //case 1:  return 9; //TODO DEBUGGING -- remove me
        case 0:  return 9;  // ... 80-100%
        case -1: return 15; // pp
        default: return 0;
      }
      default: return 0;
    }
    default: return 0;
  }
}



void ProcessUnfolding (const char* inFileTag, const char* outFileTag) {

  TFile* inFile = nullptr;

  TH1D***   h_jet_pt_ref                          = Get2DArray <TH1D*> (2, nVar);
  TH1D****  h_jet_pt                              = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);

  TH1D*****   h_jet_trk_pt_ref_sig                = Get4DArray <TH1D*> (2, nPtJBins, nDir, nVar);
  TH1D******  h_jet_trk_pt_sig                    = Get5DArray <TH1D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_pt_ref_sig             = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt_sig                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_dphi_ref_sig           = Get4DArray <TH1D*> (2, 2, nPtChSelections, nVar);
  TH1D******  h_jetInt_trk_dphi_sig               = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);

  TH1D****    h_jet_pt_ref_unf                    = Get3DArray <TH1D*> (2, 2, nVar);          // iDType, iPtJInt, nVar
  TH1D*****   h_jet_pt_unf                        = Get4DArray <TH1D*> (2, 2, nZdcCentBins+1, nVar);

  // now the pTJet-integrated histograms (e.g. > 30 GeV and > 60 GeV)
  TH1D*****   h_jetInt_trk_pt_ref_unf             = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt_unf                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);

  TH1D******  h_jetInt_trk_pt_iaa                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);
  TH1D******  h_jetInt_trk_pt_iaaNoUnf            = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);

  TH1D******  h_jetInt_trk_dphi_iaa               = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);

  TH1D***     h_jetInt_trk_pt_ref_unf_closure     = Get2DArray <TH1D*> (2, nDir);
  TGAE***     g_jetInt_trk_pt_ref_unf_closure     = Get2DArray <TGAE*> (2, nDir);
  TH1D****    h_jetInt_trk_pt_unf_closure         = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);
  TGAE****    g_jetInt_trk_pt_unf_closure         = Get3DArray <TGAE*> (2, nDir, nZdcCentBins+1);
  TH1D****    h_jetInt_trk_pt_iaa_closure         = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);
  TGAE****    g_jetInt_trk_pt_iaa_closure         = Get3DArray <TGAE*> (2, nDir, nZdcCentBins+1);


  RooUnfoldResponse**    rooUnfResp_jet_pt_ref                    = Get1DArray <RooUnfoldResponse*> (3); // with (1) vs. without (0) jet pT weights, or (2) for no weights ("full closure")
  RooUnfoldResponse***   rooUnfResp_jet_pt                        = Get2DArray <RooUnfoldResponse*> (3, nFcalCentBins+1);
  RooUnfoldResponse***   rooUnfResp_jet_trk_pt_ref_sig            = Get2DArray <RooUnfoldResponse*> (3, nDir);
  RooUnfoldResponse****  rooUnfResp_jet_trk_pt_sig                = Get3DArray <RooUnfoldResponse*> (3, nDir, nFcalCentBins+1);


  TGAE****  g_jet_trk_pt_ref_sig_syst         = Get3DArray <TGAE*> (nPtJBins, nDir, nVar);

  TGAE***** g_jet_trk_pt_sig_syst             = Get4DArray <TGAE*> (nPtJBins, nDir, nZdcCentBins+1, nVar);


  TGAE****  g_jetInt_trk_pt_ref_unf_syst      = Get3DArray <TGAE*> (2, nDir, nVar);

  TGAE***** g_jetInt_trk_pt_unf_syst          = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_pt_iaa_syst          = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);


  TGAE****  g_jetInt_trk_dphi_ref_sig_syst    = Get3DArray <TGAE*> (2, nPtChSelections, nVar);

  TGAE***** g_jetInt_trk_dphi_sig_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_dphi_iaa_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);


  TGAE****  g_jetInt_trk_pt_ref_unf_systTot   = Get3DArray <TGAE*> (2, nDir, 3);

  TGAE***** g_jetInt_trk_pt_unf_systTot       = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_pt_iaa_systTot       = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);

 
  TGAE***** g_jetInt_trk_dphi_iaa_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, 3);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/ProcessUnfolding_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // BEGIN READING IN HISTOGRAMS AND GRAPHS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {

    const TString inFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), inFileTag);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");


    std::cout << "Reading in lots of histograms & graphs..." << std::endl;
    for (short iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];
        h_jet_pt_ref[iDType][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_%s_%s", dType.Data (), var.Data ()));

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          h_jet_pt[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_%s_%s", cent.Data (), dType.Data (), var.Data ()));

        } // end loop over iCent


        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_%s",  dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
    
              h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

            } // end loop over iCent

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            const TString ptch = pTChSelections[iPtCh].Data ();

            h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_ref_sig_%s_%s_%s",  ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
    
              h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_dphi_%s_pPb_sig_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

            } // end loop over iCent

          } // end loop over iPtCh

        } // end loop over iPtJInt

      } // end loop over iVar


      {
        const int iVar = 0;
        const TString var = "Nominal";

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar]  = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_sig_%s_%s_%s",  dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
    
              h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_sig_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            } // end loop over iCent

          } // end loop over iDir

        } // end loop over iPtJ

      } // end iVar = 0 scope

    } // end loop over iDType



    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (short iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar] = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_ref_sig_syst_%s_%s",  dir.Data (), pTJ.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jet_trk_pt_%s_sig_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJ.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

      } // end loop over iDir

    } // end loop over iPtJ


    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        const TString ptch = pTChSelections[iPtCh].Data ();

        for (short iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar] = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_ref_sig_syst_%s_%s",  ptch.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jetInt_trk_dphi_%s_sig_syst_%s_%s_%s", ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iVar

      } // end loop over iPtCh

    } // end loop over iPtJInt
    std::cout << "Done reading in histograms & graphs!" << std::endl;

  }
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // DONE READING IN HISTOGRAMS AND GRAPHS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // READ IN RESPONSE MATRICES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Reading in response matrices..." << std::endl;
  {
    inFile = new TFile (Form ("%s/MakeResponseMatrix/Nominal/allSamples.root", rootPath.Data ()), "read");

    rooUnfResp_jet_pt_ref[2] = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_fullClosure");
    rooUnfResp_jet_pt_ref[1] = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_wgts");
    rooUnfResp_jet_pt_ref[0] = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_altwgts");

    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      rooUnfResp_jet_trk_pt_ref_sig[2][iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_fullClosure", dir.Data ()));
      rooUnfResp_jet_trk_pt_ref_sig[1][iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_wgts", dir.Data ()));
      rooUnfResp_jet_trk_pt_ref_sig[0][iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_altwgts", dir.Data ()));

    } // end loop over iDir

    for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

      const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

      rooUnfResp_jet_pt[2][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_fullClosure", cent));
      rooUnfResp_jet_pt[1][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_wgts", cent));
      rooUnfResp_jet_pt[0][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_altwgts", cent));

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        rooUnfResp_jet_trk_pt_sig[2][iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_fullClosure", dir.Data (), cent));
        rooUnfResp_jet_trk_pt_sig[1][iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_wgts", dir.Data (), cent));
        rooUnfResp_jet_trk_pt_sig[0][iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_altwgts", dir.Data (), cent));

      } // end loop over iFile

    } // end loop over iDir
  }




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // READ IN CLOSURE STUDY HISTOGRAMS AND CORRECTION FUNCTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Reading in closure study histograms & fit functions..." << std::endl;
  {
    inFile = new TFile (Form ("%s/aux/MCClosureHists.root", workPath.Data ()), "read");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        h_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_unf_%s_fullClosure", dir.Data (), pTJInt.Data ()));
        //f_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir] = (TF1*)  inFile->Get (Form ("f_jetInt_trk_pt_%s_ref_sig_unf_%s_fullClosure", dir.Data (), pTJInt.Data ()));
        g_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir] = make_graph (h_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir]);

        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

          const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

          h_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_%s_sig_unf_%s_fullClosure", dir.Data (), cent, pTJInt.Data ()));
          //f_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent] = (TF1*)  inFile->Get (Form ("f_jetInt_trk_pt_%s_%s_sig_unf_%s_fullClosure", dir.Data (), cent, pTJInt.Data ()));
          g_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent] = make_graph (h_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent]);

          h_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_%s_iaa_unf_%s_fullClosure", dir.Data (), cent, pTJInt.Data ()));
          //f_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent] = (TF1*)  inFile->Get (Form ("f_jetInt_trk_pt_%s_%s_iaa_unf_%s_fullClosure", dir.Data (), cent, pTJInt.Data ()));
          g_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent] = make_graph (h_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent]);

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over iPtJInt
  }




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE 2D RESULT HISTOGRAMS, DO NOMINAL 2D BAYES UNFOLD
  // THEN CONVERT BACK TO 1D HISTOGRAMS FOR INTEGRATED JET PT SELECTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Unfolding a lot of histograms..." << std::endl;
  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      const bool doUnfold = (variationsWithNoUnfold.count (var) == 0);
      const short jetWgtsType = (iDType == 1 ? 2 : ((variationsWithUnwgtdRespMatrix.count (var) == 0) ? 1 : 0));

      if (doUnfold) {
        std::cout << "Doing 1D unfold for var = " << var << std::endl;

        for (short iPtJInt : {0, 1}) { 
          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          {
            RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_pt_ref[jetWgtsType], h_jet_pt_ref[iDType][iVar], GetJetSpectraNIters (iPtJInt, -1));
            bayesUnf->SetVerbose (-1);
            h_jet_pt_ref_unf[iDType][iPtJInt][iVar] = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_ref_unf_%s_%s_%s", dType.Data (), pTJInt.Data (), var.Data ()));
            SaferDelete (&bayesUnf);
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  
            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
  
            RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_pt[jetWgtsType][iCent], h_jet_pt[iDType][iCent][iVar], GetJetSpectraNIters (iPtJInt, iCent));
            bayesUnf->SetVerbose (-1);
            h_jet_pt_unf[iDType][iPtJInt][iCent][iVar] = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_unf_%s_%s_%s_%s", dType.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            SaferDelete (&bayesUnf);

          } // end loop over iCent

        } // end loop over iPtJInt
      }

      else {
        for (short iPtJInt : {0, 1}) { 

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
          h_jet_pt_ref_unf[iDType][iPtJInt][iVar] = (TH1D*) h_jet_pt_ref[iDType][iVar]->Clone (Form ("h_jet_pt_ref_unf_%s_%s_%s", dType.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
            h_jet_pt_unf[iDType][iPtJInt][iCent][iVar] = (TH1D*) h_jet_pt[iDType][iCent][iVar]->Clone (Form ("h_jet_pt_unf_%s_%s_%s_%s", dType.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iPtJInt
      }


      if (doUnfold) {

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          TH2D* h2 = new TH2D ("h2_temp", "", nPtChBins, pTChBins, nPtJBins, pTJBins);
          h2->Sumw2 ();

          std::cout << "Doing 2D unfold for var = " << var << ", iCent = pp, iDir = " << iDir << std::endl;

          for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

            TGAE* g = g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar];
            TH1D* h = h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][0];
            const float nJet = h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1);

            for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
              if (g != nullptr && g->GetErrorYhigh (iPtCh) * g->GetErrorYlow (iPtCh) != 0)
                std::cout << "--> nonzero error in both x & y for iDType = " << iDType << ", var = " << var << ", iPtJ = " << iPtJ << ", iDir = " << iDir << ", iPtCh = " << iPtCh << std::endl;
              const float shift = (g == nullptr ? 1. : g->GetErrorYhigh (iPtCh) - g->GetErrorYlow (iPtCh)); // only one of these is nonzero, so this picks up the correct sign!
              h2->SetBinContent (iPtCh+1, iPtJ+1, nJet * (h->GetBinContent (iPtCh+1) + shift) * h->GetBinWidth (iPtCh+1));
              h2->SetBinError   (iPtCh+1, iPtJ+1, nJet * h->GetBinError   (iPtCh+1) * h->GetBinWidth (iPtCh+1)); // assumes no additional statistical uncertainty from variation... probably ok
            } // end loop over iPtCh

          } // end loop over iPtJ


          for (short iPtJInt : {0, 1}) {

            const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
            const float minJetPt = (iPtJInt == 0 ? 30. : 60.);
            const float maxJetPt = 300;

            //RooUnfoldBayes* bayesUnf1D = new RooUnfoldBayes (rooUnfResp_jet_pt_ref[jetWgtsType], h_jet_pt_ref[iDType][iVar], GetJetSpectraNIters (iPtJInt, -1));
            //bayesUnf1D->SetVerbose (-1);
            //TH1D* h_unf_njet = (TH1D*) bayesUnf1D->Hreco ()->Clone ("h_unf_njet");
            //SaferDelete (&bayesUnf1D);

            RooUnfoldBayes* bayesUnf2D = new RooUnfoldBayes (rooUnfResp_jet_trk_pt_ref_sig[jetWgtsType][iDir], h2, GetTrkSpectraNIters (iPtJInt, iDir, -1));
            bayesUnf2D->SetVerbose (-1);
            TH2D* h2_unf = (TH2D*) bayesUnf2D->Hreco ()->Clone ("h2_unf");
            SaferDelete (&bayesUnf2D);

            TH1D* h = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            h->Sumw2 ();
            h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar] = h;

            double totalJetsUF = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              totalJetsUF += h_jet_pt_ref_unf[iDType][iPtJInt][iVar]->GetBinContent (iPtJ+1);
              //totalJetsUF += h_unf_njet->GetBinContent (iPtJ+1);

              for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
                h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) + h2_unf->GetBinContent (iPtCh+1, iPtJ+1));
                h->SetBinError (iPtCh+1, h->GetBinError (iPtCh+1) + std::pow (h2_unf->GetBinError (iPtCh+1, iPtJ+1), 2));
              } // end loop over iPtCh

            } // end loop over iPtJ

            //SaferDelete (&h_unf_njet);
            SaferDelete (&h2_unf);

            for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
              h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
              h->SetBinError (iPtCh+1, std::sqrt (h->GetBinError (iPtCh+1)) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
            } // end loop over iPtCh

          } // end loop over iPtJInt

          SaferDelete (&h2);

        } // end loop over iDir


        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent)); 

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            TH2D* h2 = new TH2D ("h2_temp", "", nPtChBins, pTChBins, nPtJBins, pTJBins);
            h2->Sumw2 ();

            std::cout << "Doing 2D unfold for var = " << var << ", iCent = " << iCent << " , iDir = " << iDir << std::endl;

            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              TGAE* g = g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar];
              TH1D* h = h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][0];
              const float nJet = h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1);

              for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
                if (g != nullptr && g->GetErrorYhigh (iPtCh) * g->GetErrorYlow (iPtCh) != 0)
                  std::cout << "--> nonzero error in both x & y for iDType = " << iDType << ", var = " << var << ", iCent = " << iCent << ", iPtJ = " << iPtJ << ", iDir = " << iDir << ", iPtCh = " << iPtCh << std::endl;
                const float shift = (g == nullptr ? 1. : g->GetErrorYhigh (iPtCh) - g->GetErrorYlow (iPtCh)); // only one of these is nonzero, so this picks up the correct sign!
                h2->SetBinContent (iPtCh+1, iPtJ+1, nJet * (h->GetBinContent (iPtCh+1) + shift) * h->GetBinWidth (iPtCh+1));
                h2->SetBinError   (iPtCh+1, iPtJ+1, nJet * h->GetBinError   (iPtCh+1) * h->GetBinWidth (iPtCh+1)); // assumes no additional statistical uncertainty from variation... probably ok
              } // end loop over iPtCh

            } // end loop over iPtJ


            for (short iPtJInt : {0, 1}) {

              const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
              const float minJetPt = (iPtJInt == 0 ? 30. : 60.);
              const float maxJetPt = 300;
  
              //RooUnfoldBayes* bayesUnf1D = new RooUnfoldBayes (rooUnfResp_jet_pt[jetWgtsType][iCent], h_jet_pt[iDType][iCent][iVar], GetJetSpectraNIters (iPtJInt, iCent));
              //bayesUnf1D->SetVerbose (-1);
              //TH1D* h_unf_njet = (TH1D*) bayesUnf1D->Hreco ()->Clone ("h_unf_njet");
              //SaferDelete (&bayesUnf1D);

              RooUnfoldBayes* bayesUnf2D = new RooUnfoldBayes (rooUnfResp_jet_trk_pt_sig[jetWgtsType][iDir][iCent], h2, GetTrkSpectraNIters (iPtJInt, iDir, iCent));
              bayesUnf2D->SetVerbose (-1);
              TH2D* h2_unf = (TH2D*) bayesUnf2D->Hreco ()->Clone ("h2_unf");
              SaferDelete (&bayesUnf2D);

              TH1D* h = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
              h->Sumw2 ();
              h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar] = h;

              double totalJetsUF = 0;
              for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

                if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

                totalJetsUF += h_jet_pt_unf[iDType][iPtJInt][iCent][iVar]->GetBinContent (iPtJ+1);
                //totalJetsUF += h_unf_njet->GetBinContent (iPtJ+1);

                for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
                  h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) + h2_unf->GetBinContent (iPtCh+1, iPtJ+1));
                  h->SetBinError (iPtCh+1, h->GetBinError (iPtCh+1) + std::pow (h2_unf->GetBinError (iPtCh+1, iPtJ+1), 2));
                } // end loop over iPtCh

              } // end loop over iPtJ

              //SaferDelete (&h_unf_njet);
              SaferDelete (&h2_unf);

              for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
                h->SetBinContent (iPtCh+1, h->GetBinContent (iPtCh+1) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
                h->SetBinError (iPtCh+1, std::sqrt (h->GetBinError (iPtCh+1)) / (totalJetsUF * h->GetBinWidth (iPtCh+1)));
              } // end loop over iPtCh

            } // end loop over iPtJInt

            SaferDelete (&h2);

          } // end loop over iDir

        } // end loop over iCent

      } else {

        for (short iPtJInt : {0, 1}) {

          const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar] = (TH1D*) h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent)); 

              h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar] = (TH1D*) h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));

            } // end loop over iCent

          } // end loop over iDir

        } // end loop over iPtJInt
      }

    } // end loop over iVar

  } // end loop over iDType
  std::cout << "Unfolding finished!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // REBIN BY 2
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Rebinning histograms..." << std::endl;
//                       = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 50, 60, 75, 90, 120};
  double rebinPtChBins[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 20, 40, 60, 90, 120};
  int nRebinPtChBins = sizeof (rebinPtChBins) / sizeof (rebinPtChBins[0]) - 1;
  
  for (short iDType = 0; iDType < 2; iDType++) {

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      for (short iPtJInt : {0, 1}) {

        for (short iDir = 0; iDir < nDir; iDir++) {

          ////UnscaleWidth (h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]);
          ////h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]->Rebin (2);
          RebinSomeBins (&(h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]), nRebinPtChBins, rebinPtChBins, true);
          ////ScaleWidth (h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            if ((iPtJInt == 0 && iCent <= 5) || (iPtJInt == 1 && iCent <= 5)) {
              //UnscaleWidth (h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]);
              RebinSomeBins (&(h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]), nRebinPtChBins, rebinPtChBins, true);
              //ScaleWidth (h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]);
            }

          } // end loop over iCent

        } // end loop over iDir

      } // end loop over iPtJInt

    } // end loop over iVar

  } // end loop over iDType
  std::cout << "Rebinning finished!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CALCULATE IAA RATIOS FOR INTEGRATED JET PT SELECTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Calculating I_pPb ratios..." << std::endl;
  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar] = (TH1D*) h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_iaa_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            TH1D* href = h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar];
            if ((iPtJInt == 0 && iCent <= 5) || (iPtJInt == 1 && iCent <= 5)) {
              href = (TH1D*) href->Clone ("href_temp");
              RebinSomeBins (&href, nRebinPtChBins, rebinPtChBins, true);
              h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar]->Divide (href);
              SaferDelete (&href);
            }
            else {
              h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar]->Divide (href);
            }

            h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent][iVar] = (TH1D*) h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_iaaNoUnf_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            href = h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar];
            if ((iPtJInt == 0 && iCent <= 5) || (iPtJInt == 1 && iCent <= 5)) {
              RebinSomeBins (&h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent][iVar], nRebinPtChBins, rebinPtChBins, true);
              href = (TH1D*) href->Clone ("href_temp");
              RebinSomeBins (&href, nRebinPtChBins, rebinPtChBins, true);
              h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent][iVar]->Divide (href);
              SaferDelete (&href);
            }
            else {
              h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent][iVar]->Divide (href);
            }

          } // end loop over iCent

        } // end loop over iDir


        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][iVar] = (TH1D*) h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Clone (Form ("h_jetInt_trk_dphi_%s_iaa_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][iVar]->Divide (h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]);

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iPtJInt

    } // end loop over iVar

  } // end loop over iDType
  std::cout << "Finished calculating I_pPb ratios!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CORRECT HISTOGRAMS BY NON-CLOSURE FITS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Correcting by non-closure..." << std::endl;
  for (short iDType = 0; iDType < 2; iDType++) {

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((variationsWithNoUnfold.count (var) > 0) || (iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      //const float scaleFactor = (var == "NonClosureVar" ? 1.0 : 0.5); // for non-closure systematic correct by 100% of non-closure; nominal only corrects by 50% OLDER VERSION
      const float scaleFactor = (var == "NonClosureVar" ? 1.0 : 0.0); // for non-closure systematic correct by 100% of non-closure; nominal only corrects by 50%

      if (scaleFactor == 0)
        continue;

      for (short iPtJInt : {0, 1}) {

        for (short iDir = 0; iDir < nDir; iDir++) {

          //DivideByTF1 (h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar], f_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir], scaleFactor);
          TH1D* h = h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar];
          TGAE* g = g_jetInt_trk_pt_ref_unf_closure[iPtJInt][iDir];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            double den = g->Eval (h->GetBinCenter (iX)) * scaleFactor  - scaleFactor + 1;
            h->SetBinContent (iX, h->GetBinContent (iX) / den);
            h->SetBinError   (iX, h->GetBinError   (iX) / den);
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            //DivideByTF1 (h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar], f_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent], scaleFactor);
            h = h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar];
            g = g_jetInt_trk_pt_unf_closure[iPtJInt][iDir][iCent];
            for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
              double den = g->Eval (h->GetBinCenter (iX)) * scaleFactor  - scaleFactor + 1;
              h->SetBinContent (iX, h->GetBinContent (iX) / den);
              h->SetBinError   (iX, h->GetBinError   (iX) / den);
            }

            //DivideByTF1 (h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar], f_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent], scaleFactor);
            h = h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar];
            g = g_jetInt_trk_pt_iaa_closure[iPtJInt][iDir][iCent];
            for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
              double den = g->Eval (h->GetBinCenter (iX)) * scaleFactor  - scaleFactor + 1;
              h->SetBinContent (iX, h->GetBinContent (iX) / den);
              h->SetBinError   (iX, h->GetBinError   (iX) / den);
            }

          } // end loop over iCent

        } // end loop over iDir

      } // end loop over iPtJInt

    } // end loop over iVar

  } // end loop over iDType
  std::cout << "Finished correcting by non-closure!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE TOTAL SYSTEMATIC UNCERTAINTIES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Computing systematic uncertainties..." << std::endl;
  for (short iPtJInt : {0, 1}) {

    for (short iDir = 0; iDir < nDir; iDir++) {

      g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0] = make_graph (h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir][0]);

      ResetTGAEErrors (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0]);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]    = make_graph (h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent][0]);
        g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]    = make_graph (h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent][0]);

        ResetTGAEErrors (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]);

      } // end loop over iCent

    } // end loop over iDir

    for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0]  = make_graph (h_jetInt_trk_dphi_iaa[0][iPtJInt][iPtCh][iCent][0]);

        ResetTGAEErrors (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0]);

      } // end loop over iCent

    } // end loop over iPtCh

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE SEPARATELY.
  // THEN CALCULATES SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE BY TAKING DIFFERENCES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Taking differences between variations & nominal results..." << std::endl;
  for (short iDType = 0; iDType < 2; iDType++) {

    for (short iPtJInt : {0, 1}) {

      for (short iVar = 1; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
          continue;

        for (short iDir = 0; iDir < nDir; iDir++) {

          g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar] = new TGAE ();

          CalcSystematics (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar], h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][0],  h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]);

          if (variationsToSmooth.count (var) > 0) {
            SmoothSystematics (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar]);
          }
          
          //if (variationsToSmooth.count (var) > 0) {
          //  if (variationSmoothFuncs[var] == "ppcf") {
          //    PiecewisePolynomialConstantFunc ppcf;
          //    ppcf.SetNDeriv (2);
          //    ppcf.SetDegree (4);
          //    double params[] = {10, -0.1/(2*std::pow (std::log (10), 3)), 0};

          //    TF1* fitFunc = new TF1 ("fitFunc", &ppcf, 0.5, pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], ppcf.NDF ());
          //    fitFunc->SetParameters (params);
          //    SmoothSystematics (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar], fitFunc, h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][0],  h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]);
          //    SaferDelete (&fitFunc);
          //  }
          //  else {
          //    SmoothSystematics (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar], h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][0],  h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar], variationSmoothFuncs[var]);
          //  }
          //}

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar] = new TGAE ();
            g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar] = new TGAE ();

            CalcSystematics (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]);
            CalcSystematics (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar]);

            if (variationsToSmooth.count (var) > 0) {
              SmoothSystematics (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]);
              SmoothSystematics (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]);
            }

            //if (variationsToSmooth.count (var) > 0) {
            //  if (variationSmoothFuncs[var] == "ppcf") {
            //    PiecewisePolynomialConstantFunc ppcf;
            //    ppcf.SetNDeriv (2);
            //    ppcf.SetDegree (4);
            //    double params[] = {10, -0.1/(2*std::pow (std::log (10), 3)), 0};  // params are constant, polyCoeff[3], polyCoeff[4].

            //    TF1* fitFunc = new TF1 ("fitFunc", &ppcf, 0.5, pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], ppcf.NDF ());
            //    fitFunc->SetParameters (params);
            //    SmoothSystematics (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], fitFunc, h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]);
            //    SaferDelete (&fitFunc);

            //    fitFunc = new TF1 ("fitFunc", &ppcf, 0.5, pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], ppcf.NDF ());
            //    fitFunc->SetParameters (params);
            //    SmoothSystematics (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], fitFunc, h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar]);
            //    SaferDelete (&fitFunc);
            //  }
            //  else {
            //    SmoothSystematics (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar], variationSmoothFuncs[var]);
            //    SmoothSystematics (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar], variationSmoothFuncs[var]);
            //  }
            //}

            // allow some uncertainties to not cancel in the p+Pb / pp ratio by overwriting the current uncertainties
            if (variationsThatDontCancelInRatio.count (var) != 0) {
              ResetTGAEErrors (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]);
              AddRelErrorsInQuadrature (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar], false, false);
              AddRelErrorsInQuadrature (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], false, false);
            }

          } // end loop over iCent

        } // end loop over iDir

 
       for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar] = new TGAE ();

            CalcSystematics (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar], h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][0],  h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][iVar]);

            // allow some uncertainties to not cancel in the p+Pb / pp ratio by overwriting the current uncertainties
            if (variationsThatDontCancelInRatio.count (var) != 0) {
              ResetTGAEErrors (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar]);
              AddRelErrorsInQuadrature (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar], g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], false, false);
              AddRelErrorsInQuadrature (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar], g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], false, false);
            }

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iVar

    } // end loop over iPtJInt

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // SYSTEMATIC UNCERTAINTIES DERIVED IN MC MUST HAVE CENTRAL VALUES SET BY CENTRAL VALUES IN DATA
  // THE FINAL UNCERTAINTY IS ASSIGNED TO MATCH THE FRACTIONAL UNCERTAINTY IN MC
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Scaling uncertainties from MC to data..." << std::endl;
  for (short iVar = 1; iVar < nVar; iVar++) {

    const TString var = variations[iVar];

    if (dataVariations.count (var) > 0 || mcVariations.count (var) == 0)
      continue; // skip variations already evaluated in data or that are not evaluated in MC

    for (short iPtJInt : {0, 1}) {

      for (short iDir = 0; iDir < nDir; iDir++) {

        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar],  h_jetInt_trk_pt_ref_unf[0][iPtJInt][iDir][0]);

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_unf[0][iPtJInt][iDir][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent][0]);

        } // end loop over iCent

      } // end loop over iDir

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar],  h_jetInt_trk_dphi_iaa[0][iPtJInt][iPtCh][iCent][0]);

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iPtJInt
    
  } // end loop over iVar




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // THESE GRAPHS STORE SUMMARY SYSTEMATIC UNCERTAINTIES FOR EACH CATEGORY: TRACKING, JETS, & MIXING.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Summing uncertainties into categories..." << std::endl;
  for (short iPtJInt : {0, 1}) {

    for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

      for (short iDir = 0; iDir < nDir; iDir++) {

        g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0]->Clone ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar] = (TGAE*) g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0]->Clone ();
          g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar] = (TGAE*) g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0]->Clone ();

        } // end loop over iCent

      } // end loop over iDir

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar] = (TGAE*) g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0]->Clone ();

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iTotVar

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM EACH CATEGORY INTO ONE OF THREE GRAPHS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    for (std::vector <TString> vgroup : variationGroups) {

      std::vector <int> iVars = {};
      for (TString s : vgroup) {
        const short iVar = GetVarN (s);
        if (0 <= iVar && iVar < nVar)
          iVars.push_back (iVar);
      }

      if (iVars.size () == 0)
        continue;

      else if (iVars.size () == 1) {

        const short iVar = iVars[0];
        const short iTotVar = GetTotVarN (GetTotVar (variations[iVar]));
        std::cout << "For var " << variations[iVar] << " assigning to " << totalVariations[iTotVar] << std::endl;

        for (short iDir = 0; iDir < nDir; iDir++) {

          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar], false, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            AddErrorsInQuadrature (g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar], false, true);

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            AddErrorsInQuadrature (g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar], false, true);

          } // end loop over iCent

        } // end loop over iPtCh

      }

      else {

        const short iTotVar = GetTotVarN (GetTotVar (variations[iVars[0]]));
        for (int iVar : iVars)
          std::cout << "For var " << variations[iVar] << " assigning to " << GetTotVar (variations[iVar]) << std::endl;

        for (short iDir = 0; iDir < nDir; iDir++) {

          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir], &iVars, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            AddErrorsInQuadrature (g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent],  &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent],  &iVars, true);

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            AddErrorsInQuadrature (g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent], &iVars, true);

          } // end loop over iCent

        } // end loop over iPtCh

      }

    } // end loop over iVar

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM ALL SOURCES IN QUADRATURE, STORING RESULTS IN A SINGLE GRAPH
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Summing categories into total systematic uncertainty..." << std::endl;
  for (short iPtJInt : {0, 1}) {

    for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

      for (short iDir = 0; iDir < nDir; iDir++) {
    
        AddErrorsInQuadrature (g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][0],  g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar]);
    
        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
          AddErrorsInQuadrature (g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][0], g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][0], g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar]);
    
        } // end loop over iCent
    
      } // end loop over iDir
    
      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    
        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
          AddErrorsInQuadrature (g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][0],  g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar]);
    
        } // end loop over iCent
    
      } // end loop over iPtCh

    } // end loop over iTotVar

  } // end loop over iPtJInt
  std::cout << "Finished computing systematic uncertainties!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // FINALLY WRITE OUT EVERYTHING TO A SINGLE ROOT FILE WITH ALL THE RESULTS.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    outFile->cd ();

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if (dataVariations.count (var) == 0 && mcVariations.count (var) == 0)
        continue;

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          g_jetInt_trk_pt_ref_unf_syst[iPtJInt][iDir][iVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_unf_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_pt_unf_syst[iPtJInt][iDir][iCent][iVar]->Write  (Form ("g_jetInt_trk_pt_%s_unf_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_iaa_syst[iPtJInt][iDir][iCent][iVar]->Write  (Form ("g_jetInt_trk_pt_%s_iaa_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_dphi_iaa_syst[iPtJInt][iPtCh][iCent][iVar]->Write (Form ("g_jetInt_trk_dphi_%s_iaa_syst_%s_%s_%s", ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iPtJInt

    } // end loop over iVar


    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

        const TString totVar = totalVariations[iTotVar];

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          g_jetInt_trk_pt_ref_unf_systTot[iPtJInt][iDir][iTotVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_unf_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_pt_unf_systTot[iPtJInt][iDir][iCent][iTotVar]->Write  (Form ("g_jetInt_trk_pt_%s_unf_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_iaa_systTot[iPtJInt][iDir][iCent][iTotVar]->Write  (Form ("g_jetInt_trk_pt_%s_iaa_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jetInt_trk_dphi_iaa_systTot[iPtJInt][iPtCh][iCent][iTotVar]->Write (Form ("g_jetInt_trk_dphi_%s_iaa_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iTotVar

    } // end loop over iPtJInt



    for (short iDType = 0; iDType < 2; iDType++) {

      for (short iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;

        for (short iPtJInt : {0, 1})
          h_jet_pt_ref_unf[iDType][iPtJInt][iVar]->Write ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          for (short iPtJInt : {0, 1})
            h_jet_pt_unf[iDType][iPtJInt][iCent][iVar]->Write ();

        } // end loop over iCent

      } // end loop over iVar


      for (short iPtJInt : {0, 1}) {

        for (short iVar = 0; iVar < nVar; iVar++) {

          const TString var = variations[iVar];

          if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
            continue;

          for (short iDir = 0; iDir < nDir; iDir++) {

            h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir][iVar]->Write ();

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent][iVar]->Write ();

            } // end loop over iCent

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              h_jetInt_trk_dphi_iaa[iDType][iPtJInt][iPtCh][iCent][iVar]->Write ();

            } // end loop over iCent

          } // end loop over iPtCh

        } // end loop over iVar

      } // end loop over iPtJInt

    } // end loop over iDType


    outFile->Close ();
  }
}


int main  (int argn, char** argv) {
  assert (argn >= 3); 
  ProcessUnfolding (argv[1], argv[2]);
}


#endif
