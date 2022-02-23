#ifndef __JetHadronCorrelator_ProcessNIters_C__
#define __JetHadronCorrelator_ProcessNIters_C__

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

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

#include <iostream>
#include <thread>
#include <string.h>
#include <math.h>
#include <chrono>

using namespace JetHadronCorrelations;


void Unfold2DRef (
    RooUnfoldResponse* resp,
    TH2D* h2_raw,
    TMatrixD* h2_cov,
    TH1D****& h_unf_arr,
    TH1D****& h_rfld_arr,
    const short iDir,
    const TString baseName,
    const double* nItersVals,
    const short nItersLen
) {

  for (short iIter = 0; iIter < nItersLen; iIter++) {
    const short nIters = (short) nItersVals[iIter];

    TString bayesName = baseName;
    bayesName.ReplaceAll ("h_", "bayes_");
    bayesName.ReplaceAll ("nIters", Form ("nIters%i", nIters));
    bayesName.ReplaceAll ("REPLACE_", "");
    bayesName.ReplaceAll ("pTJInt_", "");

    RooUnfoldBayes* bayesUnf2D = new RooUnfoldBayes (bayesName, bayesName);
    bayesUnf2D->SetResponse (resp);
    bayesUnf2D->SetMeasured (h2_raw);
    bayesUnf2D->SetMeasuredCov (*h2_cov);
    bayesUnf2D->SetIterations (nIters);
    bayesUnf2D->SetSmoothing (false);
    bayesUnf2D->SetVerbose (-1);

    TH2D* h2_unf = (TH2D*) bayesUnf2D->Hreco ()->Clone (Form ("h2_unf_%iIters", nIters));

    delete bayesUnf2D;
    //SaferDelete (&bayesUnf2D);

    TH2D* h2_rfld = (TH2D*) resp->ApplyToTruth (h2_unf)->Clone (Form ("h2_rfld_%iIters", nIters));

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
      const float minJetPt = (iPtJInt == 0 ? 30. : 60.);
      const float maxJetPt = 300;

      TString unfName = baseName;
      unfName.ReplaceAll ("REPLACE", "unf");
      unfName.ReplaceAll ("nIters", Form ("nIters%i", nIters));
      unfName.ReplaceAll ("pTJInt", pTJInt);

      TH1D* h_unf = new TH1D (unfName, "", nPtChBins, pTChBins);
      h_unf->Sumw2 ();
      h_unf_arr[iPtJInt][iDir][iIter] = h_unf;

      TString rfldName = baseName;
      rfldName.ReplaceAll ("REPLACE", "rfld");
      rfldName.ReplaceAll ("nIters", Form ("nIters%i", nIters));
      rfldName.ReplaceAll ("pTJInt", pTJInt);

      TH1D* h_rfld = new TH1D (rfldName, "", nPtChBins, pTChBins);
      h_rfld->Sumw2 ();
      h_rfld_arr[iPtJInt][iDir][iIter] = h_rfld;
  
      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
    
        if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;
    
        for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
          h_unf->SetBinContent (iPtCh+1, h_unf->GetBinContent (iPtCh+1) + h2_unf->GetBinContent (iPtCh+1, iPtJ+1));
          h_unf->SetBinError (iPtCh+1, h_unf->GetBinError (iPtCh+1) + std::pow (h2_unf->GetBinError (iPtCh+1, iPtJ+1), 2));

          h_rfld->SetBinContent (iPtCh+1, h_rfld->GetBinContent (iPtCh+1) + h2_rfld->GetBinContent (iPtCh+1, iPtJ+1));
          h_rfld->SetBinError (iPtCh+1, h_rfld->GetBinError (iPtCh+1) + std::pow (h2_rfld->GetBinError (iPtCh+1, iPtJ+1), 2));
        } // end loop over iPtCh
    
      } // end loop over iPtJ

      for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
        h_unf->SetBinContent (iPtCh+1, h_unf->GetBinContent (iPtCh+1) / (h_unf->GetBinWidth (iPtCh+1)));
        h_unf->SetBinError (iPtCh+1, std::sqrt (h_unf->GetBinError (iPtCh+1)) / (h_unf->GetBinWidth (iPtCh+1)));

        h_rfld->SetBinContent (iPtCh+1, h_rfld->GetBinContent (iPtCh+1) / (h_rfld->GetBinWidth (iPtCh+1)));
        h_rfld->SetBinError (iPtCh+1, std::sqrt (h_rfld->GetBinError (iPtCh+1)) / (h_rfld->GetBinWidth (iPtCh+1)));
      } // end loop over iPtCh

    } // end loop over iPtJInt

    delete h2_unf;
    delete h2_rfld;

  } // end loop over iIter

  delete h2_raw;
  delete h2_cov;

  return;
}



void Unfold2D (
    RooUnfoldResponse* resp,
    TH2D* h2_raw,
    TMatrixD* h2_cov,
    TH1D*****& h_unf_arr,
    TH1D*****& h_rfld_arr,
    const short iDir,
    const short iCent,
    const TString baseName,
    const double* nItersVals,
    const short nItersLen
) {

  for (short iIter = 0; iIter < nItersLen; iIter++) {
    const short nIters = (short) nItersVals[iIter];

    TString bayesName = baseName;
    bayesName.ReplaceAll ("h_", "bayes_");
    bayesName.ReplaceAll ("nIters", Form ("nIters%i", nIters));
    bayesName.ReplaceAll ("REPLACE_", "");
    bayesName.ReplaceAll ("pTJInt_", "");

    RooUnfoldBayes* bayesUnf2D = new RooUnfoldBayes (bayesName, bayesName);
    bayesUnf2D->SetResponse (resp);
    bayesUnf2D->SetMeasured (h2_raw);
    bayesUnf2D->SetMeasuredCov (*h2_cov);
    bayesUnf2D->SetIterations (nIters);
    bayesUnf2D->SetSmoothing (false);
    bayesUnf2D->SetVerbose (-1);

    TH2D* h2_unf = (TH2D*) bayesUnf2D->Hreco ()->Clone (Form ("h2_unf_%iIters", nIters));

    delete bayesUnf2D;
    //SaferDelete (&bayesUnf2D);

    TH2D* h2_rfld = (TH2D*) resp->ApplyToTruth (h2_unf)->Clone (Form ("h2_rfld_%iIters", nIters));

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
      const float minJetPt = (iPtJInt == 0 ? 30. : 60.);
      const float maxJetPt = 300;

      TString unfName = baseName;
      unfName.ReplaceAll ("REPLACE", "unf");
      unfName.ReplaceAll ("nIters", Form ("nIters%i", nIters));
      unfName.ReplaceAll ("pTJInt", pTJInt);

      TH1D* h_unf = new TH1D (unfName, "", nPtChBins, pTChBins);
      h_unf->Sumw2 ();
      h_unf_arr[iPtJInt][iDir][iCent][iIter] = h_unf;

      TString rfldName = baseName;
      rfldName.ReplaceAll ("REPLACE", "rfld");
      rfldName.ReplaceAll ("nIters", Form ("nIters%i", nIters));
      rfldName.ReplaceAll ("pTJInt", pTJInt);

      TH1D* h_rfld = new TH1D (rfldName, "", nPtChBins, pTChBins);
      h_rfld->Sumw2 ();
      h_rfld_arr[iPtJInt][iDir][iCent][iIter] = h_rfld;
  
      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
    
        if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;
    
        for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
          h_unf->SetBinContent (iPtCh+1, h_unf->GetBinContent (iPtCh+1) + h2_unf->GetBinContent (iPtCh+1, iPtJ+1));
          h_unf->SetBinError (iPtCh+1, h_unf->GetBinError (iPtCh+1) + std::pow (h2_unf->GetBinError (iPtCh+1, iPtJ+1), 2));

          h_rfld->SetBinContent (iPtCh+1, h_rfld->GetBinContent (iPtCh+1) + h2_rfld->GetBinContent (iPtCh+1, iPtJ+1));
          h_rfld->SetBinError (iPtCh+1, h_rfld->GetBinError (iPtCh+1) + std::pow (h2_rfld->GetBinError (iPtCh+1, iPtJ+1), 2));
        } // end loop over iPtCh
    
      } // end loop over iPtJ

      for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
        h_unf->SetBinContent (iPtCh+1, h_unf->GetBinContent (iPtCh+1) / (h_unf->GetBinWidth (iPtCh+1)));
        h_unf->SetBinError (iPtCh+1, std::sqrt (h_unf->GetBinError (iPtCh+1)) / (h_unf->GetBinWidth (iPtCh+1)));

        h_rfld->SetBinContent (iPtCh+1, h_rfld->GetBinContent (iPtCh+1) / (h_rfld->GetBinWidth (iPtCh+1)));
        h_rfld->SetBinError (iPtCh+1, std::sqrt (h_rfld->GetBinError (iPtCh+1)) / (h_rfld->GetBinWidth (iPtCh+1)));
      } // end loop over iPtCh

    } // end loop over iPtJInt

    delete h2_unf;
    delete h2_rfld;

  } // end loop over iIter

  delete h2_raw;
  delete h2_cov;

}


void ProcessNIters (const char* rawTag, const char* outFileTag, const short nItersMax = 20, const short nIters1DMax = 1000) {

  if (nIters1DMax > 100) {
    std::cout << "About to process nIters = 1, ..., " << nIters1DMax << " for jet spectra (time complexity O(" << nIters1DMax << "^2)). Continue? [Y/N] ";
    std::string in;
    std::cin >> in;
    if (in.c_str ()[0] != 'Y' && in.c_str ()[0] != 'y') {
      std::cout << "Quitting gracefully." << std::endl;
      return;
    }
  }

  const short nItersMin = 1;
  const double* nItersVals = linspace (nItersMin, nItersMax, nItersMax-nItersMin);

  const short nIters1DMin = nItersMin;
  const double* nIters1DVals = linspace (nIters1DMin, nIters1DMax, nIters1DMax-nIters1DMin);

  const bool useJetWgts = true;
  const bool useCentDiffUnf = true;

  TH1D*     h_jet_pt_ref                  = nullptr;
  TH1D**    h_jet_pt                      = Get1DArray <TH1D*> (nZdcCentBins+1);
  TH2D*     h2_jet_pt_ref_cov             = nullptr;
  TH2D**    h2_jet_pt_cov                 = Get1DArray <TH2D*> (nZdcCentBins+1);

  //TH1D**     h_jet_pt_ref_unf             = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+1);
  //TH1D***    h_jet_pt_unf                 = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+1);

  //TH1D**     h_jet_pt_ref_rfld            = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+1);
  //TH1D***    h_jet_pt_rfld                = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+1);

  TH1D***     h_jet_trk_pt_ref_sig        = Get2DArray <TH1D*> (nPtJBins, nDir);
  TH1D****    h_jet_trk_pt_sig            = Get3DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1);
  TH1D***     h_jetInt_trk_pt_ref_sig     = Get2DArray <TH1D*> (2, nDir);
  TH1D****    h_jetInt_trk_pt_sig         = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);

  TH1D**      h_jet_pt_ref_unf_nIters     = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+1);
  TH1D***     h_jet_pt_unf_nIters         = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+1);

  TH1D**      h_jet_pt_ref_rfld_nIters     = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+1);
  TH1D***     h_jet_pt_rfld_nIters         = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+1);

  // now the pTJet-integrated histograms (e.g. > 30 GeV and > 60 GeV)
  TH1D****    h_jetInt_trk_pt_ref_unf_nIters      = Get3DArray <TH1D*> (nPtJBins, nDir, nItersMax-nItersMin+1);
  TH1D*****   h_jetInt_trk_pt_unf_nIters          = Get4DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1, nItersMax-nItersMin+1);

  TH1D****    h_jetInt_trk_pt_ref_rfld_nIters     = Get3DArray <TH1D*> (nPtJBins, nDir, nItersMax-nItersMin+1);
  TH1D*****   h_jetInt_trk_pt_rfld_nIters         = Get4DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1, nItersMax-nItersMin+1);


  RooUnfoldResponse*    rooUnfResp_jet_pt_ref                     = nullptr;
  RooUnfoldResponse**   rooUnfResp_jet_pt                         = Get1DArray <RooUnfoldResponse*> (nFcalCentBins+1);
  RooUnfoldResponse**   rooUnfResp_jet_trk_pt_ref_sig             = Get1DArray <RooUnfoldResponse*> (nDir);
  RooUnfoldResponse***  rooUnfResp_jet_trk_pt_sig                 = Get2DArray <RooUnfoldResponse*> (nDir, nFcalCentBins+1);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/ProcessNIters_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // BEGIN READING IN HISTOGRAMS AND GRAPHS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {

    const TString inFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), rawTag);
    std::cout << "Reading " << inFileName.Data () << std::endl;
    TFile* inFile = new TFile (inFileName.Data (), "read");


    h_jet_pt_ref = (TH1D*) inFile->Get ("h_jet_pt_ref_data_Nominal");
    h2_jet_pt_ref_cov = (TH2D*) inFile->Get ("h2_jet_pt_cov_ref_data_Nominal");

    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      h_jet_pt[iCent] = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_data_Nominal", cent.Data ()));
      h2_jet_pt_cov[iCent] = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_pPb_%s_data_Nominal", cent.Data ()));

    } // end loop over iCent

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        h_jetInt_trk_pt_ref_sig[iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_data_%s_Nominal",  dir.Data (), pTJInt.Data ()));

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
    
          h_jetInt_trk_pt_sig[iPtJInt][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_data_%s_Nominal", dir.Data (), cent.Data (), pTJInt.Data ()));

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over iPtJInt


    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        h_jet_trk_pt_ref_sig[iPtJ][iDir]  = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_ref_sig_data_%s_Nominal",  dir.Data (), pTJ.Data ()));

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
    
          h_jet_trk_pt_sig[iPtJ][iDir][iCent] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_pPb_sig_%s_data_%s_Nominal", dir.Data (), cent.Data (), pTJ.Data ()));

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over iPtJ
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // DONE READING IN HISTOGRAMS AND GRAPHS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // READ IN RESPONSE MATRICES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    TFile* inFile = new TFile (Form ("%s/MakeResponseMatrix/Nominal/allSamples.root", rootPath.Data ()), "read");

    if (useJetWgts) rooUnfResp_jet_pt_ref = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_wgts");
    else            rooUnfResp_jet_pt_ref = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_fullClosure");

    //for (short iDir = 0; iDir < nDir; iDir++) {
    for (short iDir : {0, 2}) {

      const TString dir = directions[iDir];

      if (useJetWgts) rooUnfResp_jet_trk_pt_ref_sig[iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_wgts", dir.Data ()));
      else            rooUnfResp_jet_trk_pt_ref_sig[iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_fullClosure", dir.Data ()));

    } // end loop over iDir

    for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

      const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

      if (useJetWgts) rooUnfResp_jet_pt[iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_wgts", cent));
      else            rooUnfResp_jet_pt[iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_fullClosure", cent));

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        if (useJetWgts) rooUnfResp_jet_trk_pt_sig[iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_wgts", dir.Data (), cent));
        else            rooUnfResp_jet_trk_pt_sig[iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_fullClosure", dir.Data (), cent));

      } // end loop over iFile

    } // end loop over iDir
  }




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE 2D RESULT HISTOGRAMS AND DO 2D BAYES UNFOLD STUDY VS # OF ITERATIONS
  // THEN CONVERT BACK TO 1D HISTOGRAMS FOR INTEGRATED JET PT SELECTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    {
      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {

        const short nIters = (short) nIters1DVals[iIter];

        RooUnfoldResponse* resp = rooUnfResp_jet_pt_ref;
        resp->UseOverflow (0);
        RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (resp, h_jet_pt_ref, nIters);
        bayesUnf->SetVerbose (-1);
        //bayesUnf->SetMeasuredCov (TMatrixD (nPtJBins, nPtJBins, h2_jet_pt_ref_cov->GetArray ()));
        TH1D* h_unf = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_ref_unf_data_Nominal_nIters%i", nIters));
        TH1D* h_rfld = (TH1D*) resp->ApplyToTruth (h_unf)->Clone (Form ("h_jet_pt_ref_rfld_data_Nominal_nIters%i", nIters));
        SaferDelete (&bayesUnf);

        h_jet_pt_ref_unf_nIters[iIter] = h_unf;
        h_jet_pt_ref_rfld_nIters[iIter] = h_rfld;

      } // end loop over iIter
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
      const short iUnfCent = (useCentDiffUnf ? iCent : nZdcCentBins);

      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
  
        const short nIters = (short) nIters1DVals[iIter];

        RooUnfoldResponse* resp = rooUnfResp_jet_pt[iUnfCent];
        resp->UseOverflow (0);
        RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (resp, h_jet_pt[iCent], nIters);
        bayesUnf->SetVerbose (-1);
        //bayesUnf->SetMeasuredCov (TMatrixD (nPtJBins, nPtJBins, h2_jet_pt_cov[iCent]->GetArray ()));
        TH1D* h_unf = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_unf_data_%s_Nominal_nIters%i", cent.Data (), nIters));
        TH1D* h_rfld = (TH1D*) resp->ApplyToTruth (h_unf)->Clone (Form ("h_jet_pt_rfld_data_%s_Nominal_nIters%i", cent.Data (), nIters));
        SaferDelete (&bayesUnf);

        h_jet_pt_unf_nIters[iCent][iIter] = h_unf;
        h_jet_pt_rfld_nIters[iCent][iIter] = h_rfld;

      } // end loop over iIter

    } // end loop over iCent


    //for (short iDir = 0; iDir < nDir; iDir++) {
    for (short iDir : {0, 2}) {

      std::vector <std::thread> threads;

      const TString dir = directions[iDir];


      // Designate response matrices to not use histogram overflow bins
      {
        rooUnfResp_jet_trk_pt_ref_sig[iDir]->UseOverflow (0);

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          rooUnfResp_jet_trk_pt_sig[iDir][iCent]->UseOverflow (0);

        } // end loop over iCent
      }


      auto start = std::chrono::steady_clock::now (); // time at beginning of 2D unfolds



      {
        std::cout << "Doing 2D unfold for iCent = pp, iDir = " << iDir << std::endl;

        TH2D* h2_raw = new TH2D ("h2_raw_ref", "", nPtChBins, pTChBins, nPtJBins, pTJBins);
        h2_raw->Sumw2 ();

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          TH1D* h = h_jet_trk_pt_ref_sig[iPtJ][iDir];

          for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
            h2_raw->SetBinContent (iPtCh+1, iPtJ+1, h->GetBinContent (iPtCh+1) * h->GetBinWidth (iPtCh+1));
            h2_raw->SetBinError   (iPtCh+1, iPtJ+1, h->GetBinError   (iPtCh+1) * h->GetBinWidth (iPtCh+1));
          } // end loop over iPtCh

        } // end loop over iPtJ

        TMatrixD* h2_cov = GetCovarianceMatrix (rootPath + "/Results/ProcessCovarianceMatrices_AllJets/" + TString ("h2_jetAll_trk_pt_cov_") + dir + TString ("_ref_sig_data.txt"));
        TString baseName = Form ("h_jetInt_trk_pt_%s_ref_REPLACE_data_pTJInt_Nominal_nIters", dir.Data ());
        //Unfold2DRef (rooUnfResp_jet_trk_pt_ref_sig[iDir], h2_raw, h2_cov, h_jetInt_trk_pt_ref_unf_nIters, h_jetInt_trk_pt_ref_rfld_nIters, iDir, baseName, nItersVals, nItersMax-nItersMin+1);

        std::thread t (Unfold2DRef, rooUnfResp_jet_trk_pt_ref_sig[iDir], h2_raw, h2_cov, std::ref (h_jetInt_trk_pt_ref_unf_nIters), std::ref (h_jetInt_trk_pt_ref_rfld_nIters), iDir, baseName, nItersVals, nItersMax-nItersMin+1);
        threads.push_back (std::move (t));
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent)); 
        const short iUnfCent = (useCentDiffUnf ? iCent : nZdcCentBins);

        const TString dir = directions[iDir];

        TH2D* h2_raw = new TH2D (Form ("h2_raw_pPb_%s", cent.Data ()), "", nPtChBins, pTChBins, nPtJBins, pTJBins);
        h2_raw->Sumw2 ();

        std::cout << "Doing 2D unfold for iCent = " << iCent << " , iDir = " << iDir << std::endl;

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          TH1D* h = h_jet_trk_pt_sig[iPtJ][iDir][iCent];

          for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
            h2_raw->SetBinContent (iPtCh+1, iPtJ+1, h->GetBinContent (iPtCh+1) * h->GetBinWidth (iPtCh+1));
            h2_raw->SetBinError   (iPtCh+1, iPtJ+1, h->GetBinError   (iPtCh+1) * h->GetBinWidth (iPtCh+1));
          } // end loop over iPtCh

        } // end loop over iPtJ

        TMatrixD* h2_cov = GetCovarianceMatrix (rootPath + "/Results/ProcessCovarianceMatrices_AllJets/" + TString ("h2_jetAll_trk_pt_cov_") + dir + TString ("_pPb_sig_") + cent + TString ("_data.txt"));
        TString baseName = Form ("h_jetInt_trk_pt_%s_pPb_REPLACE_%s_data_pTJInt_Nominal_nIters", dir.Data (), cent.Data ());
        //Unfold2D (rooUnfResp_jet_trk_pt_sig[iDir][iUnfCent], h2_raw, h2_cov, h_jetInt_trk_pt_unf_nIters, h_jetInt_trk_pt_rfld_nIters, iDir, iCent, baseName, nItersVals, nItersMax-nItersMin+1);

        std::thread t (Unfold2D, rooUnfResp_jet_trk_pt_sig[iDir][iUnfCent], h2_raw, h2_cov, std::ref (h_jetInt_trk_pt_unf_nIters), std::ref (h_jetInt_trk_pt_rfld_nIters), iDir, iCent, baseName, nItersVals, nItersMax-nItersMin+1);
        threads.push_back (std::move (t));

      } // end loop over iCent

      // wait for threads to finish...
      for (std::thread& t : threads) {
	      t.join ();
      }


      auto end = std::chrono::steady_clock::now (); // time at end of correlator function.
      auto nseconds = std::chrono::duration_cast <chrono::seconds> (end - start).count ();
      std::cout << "2D unfolds for iDir = " << iDir << " completed in " << nseconds << " seconds (avg. of " << (float)nseconds / ((float)(nItersMax-nItersMin+1)) << " seconds per number of iterations)" << std::endl;

    } // end loop over iDir



    //for (short iDir = 0; iDir < nDir; iDir++) {
    for (short iDir : {0, 2}) {

      for (short iIter = 0; iIter < nItersMax-nItersMin+1; iIter++) {
  
        for (short iPtJInt : {0, 1}) {

          const float minJetPt = (iPtJInt == 0 ? 30. : 60.);
          const float maxJetPt = 300;

          TH1D* h_unf = h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter];

          TH1D* h_rfld = h_jetInt_trk_pt_ref_rfld_nIters[iPtJInt][iDir][iIter];
  
          double totalJets = 0, totalJetsUF = 0;
          for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
   
            if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;
    
            totalJets   += h_jet_pt_ref->GetBinContent (iPtJ+1);
            totalJetsUF += h_jet_pt_ref_unf_nIters[iIter]->GetBinContent (iPtJ+1);
    
          } // end loop over iPtJ

          for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
            h_unf->SetBinContent (iPtCh+1, h_unf->GetBinContent (iPtCh+1) / totalJetsUF);
            h_unf->SetBinError (iPtCh+1, h_unf->GetBinError (iPtCh+1) / totalJetsUF);

            h_rfld->SetBinContent (iPtCh+1, h_rfld->GetBinContent (iPtCh+1) / totalJets);
            h_rfld->SetBinError (iPtCh+1, h_rfld->GetBinError (iPtCh+1) / totalJets);
          } // end loop over iPtCh

        } // end loop over iPtJInt

      } // end loop over iIter


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        for (short iIter = 0; iIter < nItersMax-nItersMin+1; iIter++) {

          for (short iPtJInt : {0, 1}) {

            const float minJetPt = (iPtJInt == 0 ? 30. : 60.);
            const float maxJetPt = 300;

            TH1D* h_unf = h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter];

            TH1D* h_rfld = h_jetInt_trk_pt_rfld_nIters[iPtJInt][iDir][iCent][iIter];

            double totalJets = 0, totalJetsUF = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
    
              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;
    
              totalJets   += h_jet_pt[iCent]->GetBinContent (iPtJ+1);
              totalJetsUF += h_jet_pt_unf_nIters[iCent][iIter]->GetBinContent (iPtJ+1);
    
            } // end loop over iPtJ

            for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
              h_unf->SetBinContent (iPtCh+1, h_unf->GetBinContent (iPtCh+1) / totalJetsUF);
              h_unf->SetBinError (iPtCh+1, h_unf->GetBinError (iPtCh+1) / totalJetsUF);

              h_rfld->SetBinContent (iPtCh+1, h_rfld->GetBinContent (iPtCh+1) / totalJets);
              h_rfld->SetBinError (iPtCh+1, h_rfld->GetBinError (iPtCh+1) / totalJets);
            } // end loop over iPtCh

          } // end loop over iPtJInt

        } // end loop over iIter

      } // end loop over iCent

    } // end loop over iDir

  }




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // FINALLY WRITE OUT EVERYTHING TO A SINGLE ROOT FILE WITH ALL THE RESULTS.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    outFile->cd ();

    for (short iPtJInt : {0, 1}) {

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        h_jetInt_trk_pt_ref_sig[iPtJInt][iDir]->Write ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          h_jetInt_trk_pt_sig[iPtJInt][iDir][iCent]->Write ();

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over iPtJInt


    for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {

      h_jet_pt_ref_unf_nIters[iIter]->Write ();
      h_jet_pt_ref_rfld_nIters[iIter]->Write ();

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  
        h_jet_pt_unf_nIters[iCent][iIter]->Write ();
        h_jet_pt_rfld_nIters[iCent][iIter]->Write ();
  
      } // end loop over iCent

    }


    for (short iIter = 0; iIter < nItersMax-nItersMin+1; iIter++) {

      for (short iPtJInt : {0, 1}) {

        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {
    
          h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter]->Write ();
          h_jetInt_trk_pt_ref_rfld_nIters[iPtJInt][iDir][iIter]->Write ();

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter]->Write ();
            h_jetInt_trk_pt_rfld_nIters[iPtJInt][iDir][iCent][iIter]->Write ();

          } // end loop over iCent
  
        } // end loop over iDir

      } // end loop over iPtJInt

    } // end loop over iIter



    outFile->Close ();
  }
}


int main  (int argn, char** argv) {
  assert (argn >= 3);
  int nItersMax = 20;
  int nIters1DMax = 100;
  if (argn >= 4)
    nItersMax = atoi (argv[3]);
  if (argn >= 5)
    nIters1DMax = atoi (argv[4]);
  ProcessNIters (argv[1], argv[2], nItersMax, nIters1DMax);
  return 0;
}


#endif
