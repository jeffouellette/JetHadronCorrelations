#ifndef __JetHadronCorrelations_PlotResponseMatrix_cxx__
#define __JetHadronCorrelations_PlotResponseMatrix_cxx__

#include "Params.h"
#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"
#include "Process.h"
#include "PiecewisePolynomialConstantFunc.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TString.h>
#include <TFitResult.h>
#include <TLine.h>
#include <TLatex.h>

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

#include <iostream>
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;

PiecewisePolynomialConstantFunc polyFunc;


TF1* DoClosureFit (const TString name, TH1D* h, const float xmin, const float xmax, short polyDeg = 4, const short maxPolyDeg = 8) {
  bool succeeded = false;
  TF1* f = nullptr;
  while (polyDeg < maxPolyDeg && !succeeded) {
    polyFunc.SetNDeriv (2);
    polyFunc.SetDegree (polyDeg);
    const int ndf = polyFunc.NDF ();

    if (f)
      delete f;
    f = new TF1 (name, &polyFunc, xmin, xmax, ndf);

    f->SetParameter (0, 10);
    for (short param = 1; param <= polyDeg-2; param++) {
      f->SetParameter (param, 0);
    }

    TFitResultPtr res = h->Fit (f, "RN0QS");
    succeeded = (res->IsValid ());
    polyDeg++;
  }

  if (succeeded)
    return f;
  else {
    delete f;
    return nullptr;
  }
}


void PlotResponseMatrix (const bool doPrimTracksOnly = false) {

  //const double pTJBins[] = {20, 30, 45, 60, 80, 100, 130, 160, 200, 240, 320};
  //const short nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;
  const Color_t cols[10] = {myCyan, myLiteBlue, myLitePurple, myPurple, myRed, myMaroon, myOrange, myLiteYellow, myLiteGreen, myGreen};

  const short maxIters = 100;
  const short max2DIters = std::min ((short)6, maxIters);

  // jet yield histograms
  TH1D***       h_jet_pt_ref                      = Get2DArray <TH1D*> (2, 2);
  TH2D***       h2_jet_pt_cov_ref                 = Get2DArray <TH2D*> (2, 2);
  TH1D****      h_jet_pt                          = Get3DArray <TH1D*> (2, nFcalCentBins+1, 2);
  TH2D****      h2_jet_pt_cov                     = Get3DArray <TH2D*> (2, nFcalCentBins+1, 2);

  // total particle yield histograms
  TH2D****      h2_jet_trk_pt_ref_tot             = Get3DArray <TH2D*> (2, nDir, 2);
  TH2D*****     h2_jet_trk_pt_tot                 = Get4DArray <TH2D*> (2, nDir, nFcalCentBins+1, 2);

  TH1D*****     h_jet_trk_pt_ref_tot              = Get4DArray <TH1D*> (2, nPtJBins, nDir, 2);
  TH1D******    h_jet_trk_pt_tot                  = Get5DArray <TH1D*> (2, nPtJBins, nDir, nFcalCentBins+1, 2);
  TH2D*****     h2_jet_trk_pt_cov_ref_tot         = Get4DArray <TH2D*> (2, nPtJBins, nDir, 2);
  TH2D******    h2_jet_trk_pt_cov_tot             = Get5DArray <TH2D*> (2, nPtJBins, nDir, nFcalCentBins+1, 2);

  // particle yield denominator histograms (number of jets)
  TH1D****      h_jet_counts_ref_tot              = Get3DArray <TH1D*> (2, nPtJBins, 2);
  TH1D*****     h_jet_counts_tot                  = Get4DArray <TH1D*> (2, nPtJBins, nFcalCentBins+1, 2);

  TH1D***       h_jet_counts_bkg                  = Get2DArray <TH1D*> (nPtJBins, nFcalCentBins+1);

  // UE estimate particle yield histograms
  TH1D****      h_jet_trk_pt_bkg                  = Get3DArray <TH1D*> (nPtJBins, nDir, nFcalCentBins+1);
  TH2D****      h2_jet_trk_pt_cov_bkg             = Get3DArray <TH2D*> (nPtJBins, nDir, nFcalCentBins+1);

  // UE subtracted particle yield histograms
  TH2D* ***     h2_jet_trk_pt_ref_sig             = Get3DArray <TH2D*> (2, nDir, 2);
  TH2D* ****    h2_jet_trk_pt_sig                 = Get4DArray <TH2D*> (2, nDir, nFcalCentBins+1, 2);
  TH1D* ****    h_jet_trk_pt_ref_sig              = Get4DArray <TH1D*> (2, nPtJBins, nDir, 2);
  TH1D* *****   h_jet_trk_pt_sig                  = Get5DArray <TH1D*> (2, nPtJBins, nDir, nFcalCentBins+1, 2);
  TH1D* ****    h_jetInt_trk_pt_ref_sig           = Get4DArray <TH1D*> (2, 2, nDir, 2);
  TH1D* *****   h_jetInt_trk_pt_sig               = Get5DArray <TH1D*> (2, 2, nDir, nFcalCentBins+1, 2);
  TH1D* *****   h_jetInt_trk_pt_iaa               = Get5DArray <TH1D*> (2, 2, nDir, nFcalCentBins+1, 2);

  // unfolded jet yield histograms
  TH1D***       h_jet_pt_ref_unf                  = Get2DArray <TH1D*> (2, maxIters);
  TH1D****      h_jet_pt_unf                      = Get3DArray <TH1D*> (2, nFcalCentBins+1, maxIters);

  TH1D***       h_njet_ref_unf  = Get2DArray <TH1D*> (2, 2); // iEvFrac, iPtJInt
  TH1D****      h_njet_unf      = Get3DArray <TH1D*> (2, 2, nFcalCentBins+1); // iEvFrac, iPtJInt, iCent
  TH1D***       h_ptjet_ref_unf  = Get2DArray <TH1D*> (2, 2); // iEvFrac, iPtJInt
  TH1D****      h_ptjet_unf      = Get3DArray <TH1D*> (2, 2, nFcalCentBins+1); // iEvFrac, iPtJInt, iCent

  // unfolded UE subtracted particle yield histograms
  TH2D****      h2_jet_trk_pt_ref_sig_unf         = Get3DArray <TH2D*> (2, nDir, max2DIters);
  TH2D****      h2_jet_trk_pt_ref_sig_rfld        = Get3DArray <TH2D*> (2, nDir, max2DIters);
  TH2D*****     h2_jet_trk_pt_sig_unf             = Get4DArray <TH2D*> (2, nDir, nFcalCentBins+1, max2DIters);
  TH2D*****     h2_jet_trk_pt_sig_rfld            = Get4DArray <TH2D*> (2, nDir, nFcalCentBins+1, max2DIters);

  TH1D*****     h_jet_trk_pt_ref_sig_unf          = Get4DArray <TH1D*> (2, nPtJBins, nDir, max2DIters);
  TH1D******    h_jet_trk_pt_sig_unf              = Get5DArray <TH1D*> (2, nPtJBins, nDir, nFcalCentBins+1, max2DIters);

  TH1D*****     h_jetInt_trk_pt_ref_sig_unf       = Get4DArray <TH1D*> (2, 2, nDir, max2DIters);
  TH1D*****     h_jetInt_trk_pt_ref_sig_rfld      = Get4DArray <TH1D*> (2, 2, nDir, max2DIters);
  TH1D******    h_jetInt_trk_pt_sig_unf           = Get5DArray <TH1D*> (2, 2, nDir, nFcalCentBins+1, max2DIters);
  TH1D******    h_jetInt_trk_pt_sig_rfld          = Get5DArray <TH1D*> (2, 2, nDir, nFcalCentBins+1, max2DIters);
  TH1D******    h_jetInt_trk_pt_iaa_unf           = Get5DArray <TH1D*> (2, 2, nDir, nFcalCentBins+1, max2DIters);
  TH1D******    h_jetInt_trk_pt_iaa_rfld          = Get5DArray <TH1D*> (2, 2, nDir, nFcalCentBins+1, max2DIters);


  RooUnfoldResponse**   rooUnfResp_jet_pt_ref         = Get1DArray <RooUnfoldResponse*> (3);
  RooUnfoldResponse***  rooUnfResp_jet_pt             = Get2DArray <RooUnfoldResponse*> (3, nFcalCentBins+1);
  RooUnfoldResponse***  rooUnfResp_jet_trk_pt_ref_sig = Get2DArray <RooUnfoldResponse*> (3, nDir);
  RooUnfoldResponse**** rooUnfResp_jet_trk_pt_sig     = Get3DArray <RooUnfoldResponse*> (3, nDir, nFcalCentBins+1);


  TH1D***   h_jetInt_trk_pt_ref_sig_unf_fullClosure = Get2DArray <TH1D*> (2, nDir);
  TF1***    f_jetInt_trk_pt_ref_sig_unf_fullClosure = Get2DArray <TF1*>  (2, nDir);
  TH1D****  h_jetInt_trk_pt_sig_unf_fullClosure     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);
  TF1****   f_jetInt_trk_pt_sig_unf_fullClosure     = Get3DArray <TF1*>  (2, nDir, nZdcCentBins+1);
  TH1D****  h_jetInt_trk_pt_iaa_unf_fullClosure     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);
  TF1****   f_jetInt_trk_pt_iaa_unf_fullClosure     = Get3DArray <TF1*>  (2, nDir, nZdcCentBins+1);

  double rebinPtChBins[] = {0.4, 0.5, 0.7, 0.9, 1.2, 1.6, 2, 3, 4, 5, 6, 8, 10, 12, 20, 40, 60, 90, 120};
  int nRebinPtChBins = sizeof (rebinPtChBins) / sizeof (rebinPtChBins[0]) - 1;


  {
    TFile* inFile = new TFile (Form ("%s/MakeResponseMatrix/%s/allSamples.root", rootPath.Data (), doPrimTracksOnly ? "MCRecoJetsTruthMatchedParts" : "Nominal"), "read");

    
    //////////////////////////////////////////////////////////////////////////////////////////////////// 
    // LOAD FULL AND HALF-CLOSURE TRAINING, PLUS TEST HISTOGRAMS
    //////////////////////////////////////////////////////////////////////////////////////////////////// 
    for (short iEvFrac : {0, 1}) {

      const TString evFrac = (iEvFrac == 0 ? "half" : "full");

      for (short iVar = 0; iVar < 2; iVar++) {

        const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");

        h_jet_pt_ref[iEvFrac][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%sClosure_ref_mc_%s", evFrac.Data (), var.Data ()));

        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {

          const TString dir = directions[iDir];
      
          h2_jet_trk_pt_ref_tot[iEvFrac][iDir][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_%sClosure_%s_ref_sig_mc_%s", evFrac.Data (), dir.Data (), var.Data ()))->Clone (Form ("h2_jet_trk_pt_%s_ref_tot_mc_%s", dir.Data (), var.Data ()));

        } // end loop over iDir


        for (short iCent = 0; iCent < nFcalCentBins; iCent++) {

          const char* cent = Form ("pPb_iCent%i", iCent);

          h_jet_pt[iEvFrac][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%sClosure_%s_mc_%s", evFrac.Data (), cent, var.Data ()));

          //for (short iDir = 0; iDir < nDir; iDir++) {
          for (short iDir : {0, 2}) {

            const TString dir = directions[iDir];
      
            h2_jet_trk_pt_tot[iEvFrac][iDir][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_%sClosure_%s_%s_sig_mc_%s", evFrac.Data (), dir.Data (), cent, var.Data ()))->Clone (Form ("h2_jet_trk_pt_%s_%s_tot_mc_%s", dir.Data (), cent, var.Data ()));

          } // end loop over iDir

        } // end loop over iCent

        {
          const char* cent = "pPb_allCent";

          h_jet_pt[iEvFrac][nFcalCentBins][iVar] = (TH1D*) h_jet_pt[iEvFrac][0][iVar]->Clone (Form ("h_jet_pt_%sClosure_%s_mc_%s", evFrac.Data (), cent, var.Data ()));

          for (short iCent = 1; iCent < nFcalCentBins; iCent++) {

            h_jet_pt[iEvFrac][nFcalCentBins][iVar]->Add (h_jet_pt[iEvFrac][iCent][iVar]);

          } // end loop over iCent

          //for (short iDir = 0; iDir < nDir; iDir++) {
          for (short iDir : {0, 2}) {

            const TString dir = directions[iDir];
      
            h2_jet_trk_pt_tot[iEvFrac][iDir][nFcalCentBins][iVar] = (TH2D*) h2_jet_trk_pt_tot[iEvFrac][iDir][0][iVar]->Clone (Form ("h2_jet_trk_pt_%sClosure_%s_%s_sig_mc_%s", evFrac.Data (), dir.Data (), cent, var.Data ()));

            for (short iCent = 1; iCent < nFcalCentBins; iCent++) {

              h2_jet_trk_pt_tot[iEvFrac][iDir][nFcalCentBins][iVar]->Add (h2_jet_trk_pt_tot[iEvFrac][iDir][iCent][iVar]);

            } // end loop over iCent

          } // end loop over iDir

        }

      } // end loop over iVar

    } // end loop over iEvFrac




    //////////////////////////////////////////////////////////////////////////////////////////////////// 
    // LOAD FULL AND HALF-CLOSURE RESPONSE MATRICES
    //////////////////////////////////////////////////////////////////////////////////////////////////// 
    for (short iEvFrac : {0, 1, 2}) {

      const TString evFrac = (iEvFrac == 0 ? "halfClosure" : (iEvFrac == 1 ? "fullClosure" : "wgts"));

      rooUnfResp_jet_pt_ref[iEvFrac] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_ref_mc_%s", evFrac.Data ()));

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        rooUnfResp_jet_trk_pt_ref_sig[iEvFrac][iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_%s", dir.Data (), evFrac.Data ()));

      } // end loop over iDir

      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

        rooUnfResp_jet_pt[iEvFrac][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_%s", cent, evFrac.Data ()));

        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {

          const TString dir = directions[iDir];

          rooUnfResp_jet_trk_pt_sig[iEvFrac][iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_%s", dir.Data (), cent, evFrac.Data ()));

        } // end loop over iFile

      } // end loop over iDir

    } // end loop over iEvFrac
  }




  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// LOAD HISTOGRAMS STORING TOTAL PER-JET CHARGED PARTICLE YIELD  ESTIMATE AND SCALE BY NJET
  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //for (short iVar = 0; iVar < 2; iVar++) {

  //  const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");

  //  {
  //    TFile* inFile = new TFile (Form ("%s/Histograms/All/JetsHists/%s/mc17_5TeV_hists.root", rootPath.Data (), var.Data ()), "read");

  //    h_evt_counts_ref_tot[iVar] = (TH1D*) inFile->Get ("h_evt_counts_mc17")->Clone (Form ("h_evt_counts_ref_tot_mc_%s", var.Data ()));

  //    h_jet_pt_ref[iVar]      = (TH1D*) inFile->Get ("h_jet_pt_mc17")->Clone      (Form ("h_jet_pt_ref_mc_%s",      var.Data ()));
  //    h2_jet_pt_cov_ref[iVar] = (TH2D*) inFile->Get ("h2_jet_pt_cov_mc17")->Clone (Form ("h2_jet_pt_cov_ref_mc_%s", var.Data ()));
  //    CalcUncertainties (h_jet_pt_ref[iVar], h2_jet_pt_cov_ref[iVar], h_evt_counts_ref_tot[iVar]); // calculates correct statistical uncertainties and provides (dN_jet/N_evt)

  //    h_jet_pt_ref[iVar]->Scale (h_evt_counts_ref_tot[iVar]->GetBinContent (2));
  //    //UnscaleWidth (h_jet_pt_ref[iVar], h_evt_counts_ref_tot[iVar]->GetBinContent (2));

  //    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

  //      const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

  //      h_jet_counts_ref_tot[iPtJ][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mc17", pTJ.Data ()))->Clone (Form ("h_jet_counts_ref_tot_%s_mc_%s", pTJ.Data (), var.Data ()));

  //      //for (short iDir = 0; iDir < nDir; iDir++) {
  //      for (short iDir : {0, 2}) {

  //        const TString dir = directions[iDir];

  //        h_jet_trk_pt_ref_tot[iPtJ][iDir][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_mc17",       dir.Data (), pTJ.Data ()))->Clone (Form ("h_jet_trk_pt_%s_ref_tot_%s_mc_%s",      dir.Data (), pTJ.Data (), var.Data ()));
  //        h2_jet_trk_pt_cov_ref_tot[iPtJ][iDir][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_mc17",  dir.Data (), pTJ.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_ref_tot_%s_mc_%s", dir.Data (), pTJ.Data (), var.Data ()));

  //        CalcUncertainties (h_jet_trk_pt_ref_tot[iPtJ][iDir][iVar], h2_jet_trk_pt_cov_ref_tot[iPtJ][iDir][iVar], h_jet_counts_ref_tot[iPtJ][iVar]); // calculates correct statistical uncertainties and provides (dN_ch/N_jet)

  //        h_jet_trk_pt_ref_tot[iPtJ][iDir][iVar]->Scale (h_jet_counts_ref_tot[iPtJ][iVar]->GetBinContent (2));
  //        //UnscaleWidth (h_jet_trk_pt_ref_tot[iPtJ][iDir][iVar], h_jet_counts_ref_tot[iPtJ][iVar]->GetBinContent (2));

  //      } // end loop over iDir

  //    } // end loop over iPtJ
  //  }


  //  for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

  //    const TString cent = (iCent == nFcalCentBins ? "allCent" : Form ("iCent%i", iCent));

  //    TFile* inFile = new TFile (Form ("%s/Histograms/All/JetsHists/%s/mc16_5TeV_%s_hists.root", rootPath.Data (), var.Data (), cent.Data ()), "read");

  //    h_evt_counts_tot[iCent][iVar] = (TH1D*) inFile->Get ("h_evt_counts_mc16")->Clone (Form ("h_evt_counts_pPb_tot_%s_mc_%s", cent.Data (), var.Data ()));

  //    h_jet_pt[iCent][iVar]       = (TH1D*) inFile->Get ("h_jet_pt_mc16")->Clone      (Form ("h_jet_pt_pPb_%s_mc_%s",       cent.Data (), var.Data ()));
  //    h2_jet_pt_cov[iCent][iVar]  = (TH2D*) inFile->Get ("h2_jet_pt_cov_mc16")->Clone (Form ("h2_jet_pt_cov_pPb_%s_mc_%s",  cent.Data (), var.Data ()));
  //    CalcUncertainties (h_jet_pt[iCent][iVar], h2_jet_pt_cov[iCent][iVar], h_evt_counts_tot[iCent][iVar]); // calculates correct statistical uncertainties and provides (dN_jet/N_evt)

  //    h_jet_pt[iCent][iVar]->Scale (h_evt_counts_tot[iCent][iVar]->GetBinContent (2));
  //    //UnscaleWidth (h_jet_pt[iCent][iVar], h_evt_counts_tot[iCent][iVar]->GetBinContent (2));

  //    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

  //      const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

  //      h_jet_counts_tot[iPtJ][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mc16", pTJ.Data ()))->Clone (Form ("h_jet_counts_pPb_tot_%s_%s_mc_%s", cent.Data (), pTJ.Data (), var.Data ()));

  //      //for (short iDir = 0; iDir < nDir; iDir++) {
  //      for (short iDir : {0, 2}) {

  //        const TString dir = directions[iDir];

  //        h_jet_trk_pt_tot[iPtJ][iDir][iCent][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_mc16",       dir.Data (), pTJ.Data ()))->Clone (Form ("h_jet_trk_pt_%s_pPb_tot_%s_%s_mc_%s",      dir.Data (), cent.Data (), pTJ.Data (), var.Data ()));
  //        h2_jet_trk_pt_cov_tot[iPtJ][iDir][iCent][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_mc16",  dir.Data (), pTJ.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_pPb_tot_%s_%s_mc_%s", dir.Data (), cent.Data (), pTJ.Data (), var.Data ()));

  //        CalcUncertainties (h_jet_trk_pt_tot[iPtJ][iDir][iCent][iVar], h2_jet_trk_pt_cov_tot[iPtJ][iDir][iCent][iVar], h_jet_counts_tot[iPtJ][iCent][iVar]); // calculates correct statistical uncertainties and provides (dN_ch/N_jet)

  //        h_jet_trk_pt_tot[iPtJ][iDir][iCent][iVar]->Scale (h_jet_counts_tot[iPtJ][iCent][iVar]->GetBinContent (2));
  //        //UnscaleWidth (h_jet_trk_pt_tot[iPtJ][iDir][iCent][iVar], h_jet_counts_tot[iPtJ][iCent][iVar]->GetBinContent (2));

  //      } // end loop over iDir

  //    } // end loop over iPtJ

  //  } // end loop over iCent

  //} // end loop over iVar




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // LOAD HISTOGRAMS STORING PER-JET UE ESTIMATE AND SCALE BY NJET IN TOTAL
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

    const TString cent = (iCent == nFcalCentBins ? "allCent" : Form ("iCent%i", iCent));

    TFile* inFile = new TFile (Form ("%s/Histograms/All/MixedHists/Nominal/mc16_5TeV_%s_hists.root", rootPath.Data (), cent.Data ()), "read");
    //TFile* inFile = new TFile (Form ("%s/Histograms/All/MixedHists/MixCatVar2/mc16_5TeV_%s_hists.root", rootPath.Data (), cent.Data ()), "read");
    //TFile* inFile = new TFile (Form ("%s/Histograms/All/MixedHists/MixCatVar6/mc16_5TeV_%s_hists.root", rootPath.Data (), cent.Data ()), "read");

    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

      h_jet_counts_bkg[iPtJ][iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_mc16", pTJ.Data ()))->Clone (Form ("h_jet_counts_pPb_bkg_%s_%s_mc_Nominal", cent.Data (), pTJ.Data ()));

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        h_jet_trk_pt_bkg[iPtJ][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_mixed_mc16",       dir.Data (), pTJ.Data ()))->Clone (Form ("h_jet_trk_pt_%s_pPb_bkg_%s_%s_mc_Nominal",      dir.Data (), cent.Data (), pTJ.Data ()));
        h2_jet_trk_pt_cov_bkg[iPtJ][iDir][iCent]  = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_mixed_mc16",  dir.Data (), pTJ.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_pPb_bkg_%s_%s_mc_Nominal", dir.Data (), cent.Data (), pTJ.Data ()));

        CalcUncertainties (h_jet_trk_pt_bkg[iPtJ][iDir][iCent], h2_jet_trk_pt_cov_bkg[iPtJ][iDir][iCent], h_jet_counts_bkg[iPtJ][iCent]); // calculates correct statistical uncertainties and provides (dN_ch/N_jet)

        h_jet_trk_pt_bkg[iPtJ][iDir][iCent]->Scale (1., "width");

      } // end loop over iDir

    } // end loop over iPtJ

  } // end loop over iCent




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE SIGNAL HISTOGRAMS AND SUBTRACT OFF UNDERLYING EVENT IN P+PB
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iEvFrac : {0, 1}) {

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    for (short iVar = 0; iVar < 2; iVar++) {

      const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        {
          h2_jet_trk_pt_ref_sig[iEvFrac][iDir][iVar] = (TH2D*) h2_jet_trk_pt_ref_tot[iEvFrac][iDir][iVar]->Clone (Form ("h2_jet_trk_pt_%s_ref_sig_%s_mc_%s", dir.Data (), evFrac.Data (), var.Data ()));
        }


        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

          const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

          h2_jet_trk_pt_sig[iEvFrac][iDir][iCent][iVar] = (TH2D*) h2_jet_trk_pt_tot[iEvFrac][iDir][iCent][iVar]->Clone (Form ("h2_jet_trk_pt_%s_%s_sig_%s_mc_%s", dir.Data (), cent, evFrac.Data (), var.Data ()));

          TH2D* h2tot = h2_jet_trk_pt_tot[iEvFrac][iDir][iCent][iVar];
          TH2D* h2sig = h2_jet_trk_pt_sig[iEvFrac][iDir][iCent][iVar];

          for (short iPtJY = 1; iPtJY <= nPtJBins; iPtJY++) {

            const double jpt = 0.5 * (pTJBins[iPtJY] + pTJBins[iPtJY-1]);

            TH1D* hbkg = h_jet_trk_pt_bkg[iPtJY-1][iDir][iCent];
    
            const double nJetSF = h_jet_pt[iEvFrac][iCent][iVar]->GetBinContent (h_jet_pt[iEvFrac][iCent][iVar]->FindBin (jpt));
    
            for (short iPtChX = 1; iPtChX <= nPtChBins; iPtChX++) {

              if (doPrimTracksOnly) {
                h2sig->SetBinContent (iPtChX, iPtJY, h2tot->GetBinContent (iPtChX, iPtJY));
                h2sig->SetBinError   (iPtChX, iPtJY, h2tot->GetBinError (iPtChX, iPtJY));
              } else {
                h2sig->SetBinContent (iPtChX, iPtJY, h2tot->GetBinContent (iPtChX, iPtJY) - (iVar == 0 ? nJetSF * hbkg->GetBinContent (iPtChX) * hbkg->GetBinWidth (iPtChX) : 0));
                h2sig->SetBinError   (iPtChX, iPtJY, std::hypot (h2tot->GetBinError (iPtChX, iPtJY), (iVar == 0 ? nJetSF * hbkg->GetBinError (iPtChX) * hbkg->GetBinWidth (iPtChX) : 0)));
              }
    
            } // end loop over iPtChY
    
          } // end loop over iPtJX

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over iVar

  } // end loop over iEvFrac




/*
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // UNFOLD FULL AND HALF-CLOSURE HISTOGRAMS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Unfolding a lot of histograms many times..." << std::endl;
  for (short iEvFrac : {0, 1}) {
    if (iEvFrac == 0)
      std::cout << "|--> working on half closure test." << std::endl;
    else if (iEvFrac == 1)
      std::cout << "|--> working on full closure test." << std::endl;

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    RooUnfoldBayes* bayesUnf = nullptr;
    TH1D* h = nullptr;

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      h_njet_ref_unf[iEvFrac][iPtJInt] = new TH1D (Form ("h_njet_ref_unf_%s_%s_mc", evFrac.Data (), pTJInt.Data ()), ";Iterations;N_{jet}^{unfolded} / N_{jet}^{truth}", maxIters, 0.5, maxIters + 0.5);
      h_njet_ref_unf[iEvFrac][iPtJInt]->Sumw2 ();
      h_ptjet_ref_unf[iEvFrac][iPtJInt] = new TH1D (Form ("h_ptjet_ref_unf_%s_%s_mc", evFrac.Data (), pTJInt.Data ()), ";Iterations;#LT#it{p}_{T}^{jet}#GT^{unfolded} / #LT#it{p}_{T}^{jet}#GT^{truth}", maxIters, 0.5, maxIters + 0.5);
      h_ptjet_ref_unf[iEvFrac][iPtJInt]->Sumw2 ();

      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));
        h_njet_unf[iEvFrac][iPtJInt][iCent] = new TH1D (Form ("h_njet_unf_%s_%s_%s_mc", evFrac.Data (), cent, pTJInt.Data ()), ";Iterations;N_{jet}^{unfolded} / N_{jet}^{truth}", maxIters, 0.5, maxIters + 0.5);
        h_njet_unf[iEvFrac][iPtJInt][iCent]->Sumw2 ();
        h_ptjet_unf[iEvFrac][iPtJInt][iCent] = new TH1D (Form ("h_ptjet_unf_%s_%s_%s_mc", evFrac.Data (), cent, pTJInt.Data ()), ";Iterations;#LT#it{p}_{T}^{jet}#GT^{unfolded} / #LT#it{p}_{T}^{jet}#GT^{truth}", maxIters, 0.5, maxIters + 0.5);
        h_ptjet_unf[iEvFrac][iPtJInt][iCent]->Sumw2 ();

      } // end loop over iCent

    } // end loop over iPtJInt

    for (short nIter = 0; nIter < maxIters; nIter++) {
      std::cout << "  |--> on 1D unfold, with niter = " << nIter+1 << "..." << std::endl;

      RooUnfoldResponse* resp = rooUnfResp_jet_pt_ref[iEvFrac];
      resp->UseOverflow (0);

      TH1D* h_raw = (TH1D*) h_jet_pt_ref[iEvFrac][0]->Clone ("h_raw");
      for (int iX = 1; iX <= h_raw->GetNbinsX (); iX++) {
        if (h_raw->GetBinCenter (iX) < 20) {
          h_raw->SetBinContent (iX, 0);
          h_raw->SetBinError (iX, 0);
        }
      }

      bayesUnf = new RooUnfoldBayes (resp, h_raw, nIter+1);
      bayesUnf->SetVerbose (-1);
      h = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_ref_unf_%s_mc_nIter%i", evFrac.Data (), nIter+1));
      h_jet_pt_ref_unf[iEvFrac][nIter] = h;

      SaferDelete (&h_raw);
      SaferDelete (&bayesUnf);

      for (short iPtJInt : {0, 1}) {
        const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
        const double maxJetPt = 300;
        TH1D* htemp = (TH1D*) h->Clone ("htemp"); // ensures SetRange doesn't mess anything up later
        double err = 0;
        double njet = htemp->IntegralAndError (htemp->FindBin (minJetPt+0.01), htemp->FindBin (maxJetPt-0.01), err);
        h_njet_ref_unf[iEvFrac][iPtJInt]->SetBinContent (nIter+1, njet);
        h_njet_ref_unf[iEvFrac][iPtJInt]->SetBinError   (nIter+1, err);
        htemp->GetXaxis ()->SetRange (iPtJInt == 0  ? 3 : 4, htemp->GetNbinsX () - 2);
        std::cout << "    |--> N_jet (pp)    = " << njet << std::endl;
        std::cout << "    |--> <pT^jet> (pp) = " << htemp->GetMean () << " GeV" << std::endl;
        h_ptjet_ref_unf[iEvFrac][iPtJInt]->SetBinContent (nIter+1, htemp->GetMean ());
        h_ptjet_ref_unf[iEvFrac][iPtJInt]->SetBinError   (nIter+1, htemp->GetMeanError ());
        SaferDelete (&htemp);
      }

      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

        const short iUnfCent = iCent;

        RooUnfoldResponse* resp = rooUnfResp_jet_pt[iEvFrac][iUnfCent];
        resp->UseOverflow (0);
  
        TH1D* h_raw = (TH1D*) h_jet_pt[iEvFrac][iCent][0]->Clone ("h_raw");
        for (int iX = 1; iX <= h_raw->GetNbinsX (); iX++) {
          if (h_raw->GetBinCenter (iX) < 20) {
            h_raw->SetBinContent (iX, 0);
            h_raw->SetBinError (iX, 0);
          }
        }
  
        bayesUnf = new RooUnfoldBayes (resp, h_raw, nIter+1);
        bayesUnf->SetVerbose (-1);
        h = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_unf_%s_%s_mc_nIter%i", evFrac.Data (), cent, nIter+1));
        h_jet_pt_unf[iEvFrac][iCent][nIter] = h;
  
        SaferDelete (&h_raw);
        SaferDelete (&bayesUnf);

        for (short iPtJInt : {0, 1}) {
          const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
          const double maxJetPt = 300;
          TH1D* htemp = (TH1D*) h->Clone ("htemp"); // ensures SetRange doesn't mess anything up later
          double err = 0;
          h_njet_unf[iEvFrac][iPtJInt][iCent]->SetBinContent (nIter+1, htemp->IntegralAndError (htemp->FindBin (minJetPt+0.01), htemp->FindBin (maxJetPt-0.01), err));
          h_njet_unf[iEvFrac][iPtJInt][iCent]->SetBinError   (nIter+1, err);
          htemp->GetXaxis ()->SetRange (iPtJInt == 0  ? 3 : 4, htemp->GetNbinsX () - 2);
          h_ptjet_unf[iEvFrac][iPtJInt][iCent]->SetBinContent (nIter+1, htemp->GetMean ());
          h_ptjet_unf[iEvFrac][iPtJInt][iCent]->SetBinError   (nIter+1, htemp->GetMeanError ());
          SaferDelete (&htemp);
        }

      } // end loop over iCent

    } // end loop over nIter


    for (short nIter = 0; nIter < max2DIters; nIter++) {
      std::cout << "  |--> on 2D unfold, with niter = " << nIter+1 << "..." << std::endl;

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        RooUnfoldResponse* resp = rooUnfResp_jet_trk_pt_ref_sig[iEvFrac][iDir];
        resp->UseOverflow (0);
  
        TH2D* h2_raw = (TH2D*) h2_jet_trk_pt_ref_sig[iEvFrac][iDir][0]->Clone ("h2_raw");
        for (int iY = 1; iY <= h2_raw->GetNbinsY (); iY++) {
          if (h2_raw->GetYaxis ()->GetBinCenter (iY) < 20) {
            for (int iX = 1; iX <= h2_raw->GetNbinsX (); iX++) {
              h2_raw->SetBinContent (iX, iY, 0);
              h2_raw->SetBinError (iX, iY, 0);
            }
          }
        }
  
        bayesUnf = new RooUnfoldBayes (resp, h2_raw, nIter+1);
        bayesUnf->SetVerbose (-1);
        h2_jet_trk_pt_ref_sig_unf[iEvFrac][iDir][nIter] = (TH2D*) bayesUnf->Hreco ()->Clone (Form ("h2_jet_trk_pt_%s_ref_sig_unf_%s_mc_nIter%i", dir.Data (), evFrac.Data (), nIter+1));

        SaferDelete (&h2_raw);
        SaferDelete (&bayesUnf);

        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
    
          const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

          const short iUnfCent = iCent;

          RooUnfoldResponse* resp = rooUnfResp_jet_trk_pt_sig[iEvFrac][iDir][iUnfCent];
          resp->UseOverflow (0);
    
          TH2D* h2_raw = (TH2D*) h2_jet_trk_pt_sig[iEvFrac][iDir][iCent][0]->Clone ("h2_raw");
          for (int iY = 1; iY <= h2_raw->GetNbinsY (); iY++) {
            if (h2_raw->GetYaxis ()->GetBinCenter (iY) < 20) {
              for (int iX = 1; iX <= h2_raw->GetNbinsX (); iX++) {
                h2_raw->SetBinContent (iX, iY, 0);
                h2_raw->SetBinError (iX, iY, 0);
              }
            }
          }
    
          bayesUnf = new RooUnfoldBayes (resp, h2_raw, nIter+1);
          bayesUnf->SetVerbose (-1);
          h2_jet_trk_pt_sig_unf[iEvFrac][iDir][iCent][nIter] = (TH2D*) bayesUnf->Hreco ()->Clone (Form ("h2_jet_trk_pt_%s_%s_sig_unf_%s_mc_nIter%i", dir.Data (), cent, evFrac.Data (), nIter+1));

          SaferDelete (&h2_raw);
          SaferDelete (&bayesUnf);

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over nIter

  } // end loop over iEvFrac
  std::cout << "Finished unfolding!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // SCALE TEST HISTOGRAMS BY NJET (AND BIN WIDTH) AFTER UNFOLDING
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Scaling histograms by N_jet after unfolding..." << std::endl;
  for (short iEvFrac : {0, 1}) {

    for (short iVar = 0; iVar < 2; iVar++) {

      TH1D* h = h_jet_pt_ref[iEvFrac][iVar]; 

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        TH2D* h2 = h2_jet_trk_pt_ref_sig[iEvFrac][iDir][iVar];
        for (short iX = 1; iX <= h2->GetNbinsX (); iX++) {
          for (short iY = 1; iY <= h2->GetNbinsY (); iY++) {
            h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
            h2->SetBinError   (iX, iY, h2->GetBinError   (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
          }
        }
        h2 = h2_jet_trk_pt_ref_tot[iEvFrac][iDir][iVar];
        for (short iX = 1; iX <= h2->GetNbinsX (); iX++) {
          for (short iY = 1; iY <= h2->GetNbinsY (); iY++) {
            h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
            h2->SetBinError   (iX, iY, h2->GetBinError   (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
          }
        }

      } // end loop over iDir


      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        h = h_jet_pt[iEvFrac][iCent][iVar]; 

        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {
    
          TH2D* h2 = h2_jet_trk_pt_sig[iEvFrac][iDir][iCent][iVar];
          for (short iX = 1; iX <= h2->GetNbinsX (); iX++) {
            for (short iY = 1; iY <= h2->GetNbinsY (); iY++) {
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
              h2->SetBinError   (iX, iY, h2->GetBinError   (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
            }
          }
          h2 = h2_jet_trk_pt_tot[iEvFrac][iDir][iCent][iVar];
          for (short iX = 1; iX <= h2->GetNbinsX (); iX++) {
            for (short iY = 1; iY <= h2->GetNbinsY (); iY++) {
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
              h2->SetBinError   (iX, iY, h2->GetBinError   (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
            }
          }

        } // end loop over iDir

      } // end loop over iCent

    } // end loop over iVar


    for (short nIter = 0; nIter < max2DIters; nIter++) {

      TH1D* h = h_jet_pt_ref_unf[iEvFrac][1];

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        TH2D* h2 = h2_jet_trk_pt_ref_sig_unf[iEvFrac][iDir][nIter];
        for (short iX = 1; iX <= h2->GetNbinsX (); iX++) {
          for (short iY = 1; iY <= h2->GetNbinsY (); iY++) {
            h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
            h2->SetBinError   (iX, iY, h2->GetBinError   (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
          }
        }

      } // end loop over iDir


      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        h = h_jet_pt_unf[iEvFrac][iCent][1];

        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {
    
          TH2D* h2 = h2_jet_trk_pt_sig_unf[iEvFrac][iDir][iCent][nIter];
          for (short iX = 1; iX <= h2->GetNbinsX (); iX++) {
            for (short iY = 1; iY <= h2->GetNbinsY (); iY++) {
              h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
              h2->SetBinError   (iX, iY, h2->GetBinError   (iX, iY) / (h->GetBinContent (iY) * h2->GetXaxis ()->GetBinWidth (iX)));
            }
          }

        } // end loop over iDir

      } // end loop over iCent

    } // end loop over nIter

  } // end loop over iEvFrac
  std::cout << "Finished scaling histograms!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // GET 1D PROJECTIONS OF JET-TAGGED HADRONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Taking 1D projections of unfolded histograms..." << std::endl;
  for (short iEvFrac : {0, 1}) {

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    //for (short iDir = 0; iDir < nDir; iDir++) {
    for (short iDir : {0, 2}) {

      const TString dir = directions[iDir];

      for (short iVar = 0; iVar < 2; iVar++) {

        const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
          h_jet_trk_pt_ref_tot[iEvFrac][iPtJ][iDir][iVar] = h2_jet_trk_pt_ref_tot[iEvFrac][iDir][iVar]->ProjectionX (Form ("h_jet_trk_pt_%s_ref_tot_%s_%g-%gGeVJets_mc_%s", dir.Data (), evFrac.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1], var.Data ()), iPtJ+1, iPtJ+1);
          h_jet_trk_pt_ref_sig[iEvFrac][iPtJ][iDir][iVar] = h2_jet_trk_pt_ref_sig[iEvFrac][iDir][iVar]->ProjectionX (Form ("h_jet_trk_pt_%s_ref_sig_%s_%g-%gGeVJets_mc_%s", dir.Data (), evFrac.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1], var.Data ()), iPtJ+1, iPtJ+1);
        }

        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

          const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

          for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
            h_jet_trk_pt_tot[iEvFrac][iPtJ][iDir][iCent][iVar] = h2_jet_trk_pt_tot[iEvFrac][iDir][iCent][iVar]->ProjectionX (Form ("h_jet_trk_pt_%s_%s_tot_%s_%g-%gGeVJets_mc_%s", dir.Data (), cent, evFrac.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1], var.Data ()), iPtJ+1, iPtJ+1);
            h_jet_trk_pt_sig[iEvFrac][iPtJ][iDir][iCent][iVar] = h2_jet_trk_pt_sig[iEvFrac][iDir][iCent][iVar]->ProjectionX (Form ("h_jet_trk_pt_%s_%s_sig_%s_%g-%gGeVJets_mc_%s", dir.Data (), cent, evFrac.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1], var.Data ()), iPtJ+1, iPtJ+1);
          }

        } // end loop over iCent

      } // end loop over iVar
    
      for (short nIter = 0; nIter < max2DIters; nIter++) {

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++)
          h_jet_trk_pt_ref_sig_unf[iEvFrac][iPtJ][iDir][nIter] = h2_jet_trk_pt_ref_sig_unf[iEvFrac][iDir][nIter]->ProjectionX (Form ("h_jet_trk_pt_%s_ref_sig_%s_%g-%gGeVJets_mc_nIter%i", dir.Data (), evFrac.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1], nIter+1), iPtJ+1, iPtJ+1);

        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

          const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

          for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++)
            h_jet_trk_pt_sig_unf[iEvFrac][iPtJ][iDir][iCent][nIter] = h2_jet_trk_pt_sig_unf[iEvFrac][iDir][iCent][nIter]->ProjectionX (Form ("h_jet_trk_pt_%s_%s_sig_%s_%g-%gGeVJets_mc_nIter%i", dir.Data (), cent, evFrac.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1], nIter+1), iPtJ+1, iPtJ+1);

        } // end loop over iCent

      } // end loop over nIter

    } // end loop over iDir

  } // end loop over iEvFrac
  std::cout << "Finished taking 1D projections!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE >30 GeV AND >60 GeV HISTOGRAMS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Integrating over jet pT bins..." << std::endl;
  for (short iEvFrac : {0, 1}) {

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
      const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
      const double maxJetPt = 300;

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {
    
        const TString dir = directions[iDir];
    
        for (short iVar = 0; iVar < 2; iVar++) {
    
          const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");
    
          {
            h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_mc_%s", dir.Data (), evFrac.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
    
            double norm = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
              const double jpt = 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]);
              if (jpt < minJetPt || maxJetPt < jpt) continue;
              const double nJets = h_jet_pt_ref[iEvFrac][iVar]->GetBinContent (h_jet_pt_ref[iEvFrac][iVar]->FindBin (jpt));
              h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][iVar]->Add (h_jet_trk_pt_ref_sig[iEvFrac][iPtJ][iDir][iVar], nJets);
              norm += nJets;
            }
            h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][iVar]->Scale (1./norm);
          }
    
    
          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
    
            const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));
    
            h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_%s_sig_%s_%s_mc_%s", dir.Data (), cent, evFrac.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
    
            double norm = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
              const double jpt = 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]);
              if (jpt < minJetPt || maxJetPt < jpt) continue;
              const double nJets = h_jet_pt[iEvFrac][iCent][iVar]->GetBinContent (h_jet_pt[iEvFrac][iCent][iVar]->FindBin (jpt));
              h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][iVar]->Add (h_jet_trk_pt_sig[iEvFrac][iPtJ][iDir][iCent][iVar], nJets);
              norm += nJets;
            }
            h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][iVar]->Scale (1./norm);
    
          } // end loop over iCent
    
        } // end loop over iVar
     
     
        for (short nIter = 0; nIter < max2DIters; nIter++) {
    
          {
            h_jetInt_trk_pt_ref_sig_unf[iEvFrac][iPtJInt][iDir][nIter] = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_sig_unf_%s_%s_mc_nIter%i", dir.Data (), evFrac.Data (), pTJInt.Data (), nIter+1), "", nPtChBins, pTChBins);
            TH1D* h = h_jet_pt_ref_unf[iEvFrac][1];
    
            double norm = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
              const double jpt = 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]);
              if (jpt < minJetPt || maxJetPt < jpt) continue;
              const double nJets = h->GetBinContent (h->FindBin (jpt));
              h_jetInt_trk_pt_ref_sig_unf[iEvFrac][iPtJInt][iDir][nIter]->Add (h_jet_trk_pt_ref_sig_unf[iEvFrac][iPtJ][iDir][nIter], nJets);
              norm += nJets;
            }
            h_jetInt_trk_pt_ref_sig_unf[iEvFrac][iPtJInt][iDir][nIter]->Scale (1./norm);
          }
    
    
          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
    
            const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));
    
            h_jetInt_trk_pt_sig_unf[iEvFrac][iPtJInt][iDir][iCent][nIter] = new TH1D (Form ("h_jetInt_trk_pt_%s_%s_sig_unf_%s_%s_mc_nIter%i", dir.Data (), cent, evFrac.Data (), pTJInt.Data (), nIter+1), "", nPtChBins, pTChBins);
            TH1D* h = h_jet_pt_unf[iEvFrac][iCent][1];
    
            double norm = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
              const double jpt = 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]);
              if (jpt < minJetPt || maxJetPt < jpt) continue;
              const double nJets = h->GetBinContent (h->FindBin (jpt));
              h_jetInt_trk_pt_sig_unf[iEvFrac][iPtJInt][iDir][iCent][nIter]->Add (h_jet_trk_pt_sig_unf[iEvFrac][iPtJ][iDir][iCent][nIter], nJets);
              norm += nJets;
            }
            h_jetInt_trk_pt_sig_unf[iEvFrac][iPtJInt][iDir][iCent][nIter]->Scale (1./norm);
    
          } // end loop over iCent
    
        } // end loop over nIter
    
      } // end loop over iDir

    } // end loop over iPtJInt

  } // end loop over iEvFrac
  std::cout << "Finished integrating over jet pT bins!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // COMPUTE IAA FOR >30 GeV AND >60 GeV HISTOGRAMS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Computing I_pPb ratios..." << std::endl;
  for (short iEvFrac : {0, 1}) {

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        for (short iVar = 0; iVar < 2; iVar++) {

          const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");

          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

            const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

            h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][iVar] = (TH1D*) h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_%s_iaa_%s_%s_mc_%s", dir.Data (), cent, evFrac.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][iVar]->Divide (h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][iVar]);

          } // end loop over iCent

        } // end loop over iVar
 
 
        for (short nIter = 0; nIter < max2DIters; nIter++) {

          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

            const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

            h_jetInt_trk_pt_iaa_unf[iEvFrac][iPtJInt][iDir][iCent][nIter] = (TH1D*) h_jetInt_trk_pt_sig_unf[iEvFrac][iPtJInt][iDir][iCent][nIter]->Clone (Form ("h_jetInt_trk_pt_%s_%s_iaa_unf_%s_%s_mc_nIter%i", dir.Data (), cent, evFrac.Data (), pTJInt.Data (), nIter+1));
            h_jetInt_trk_pt_iaa_unf[iEvFrac][iPtJInt][iDir][iCent][nIter]->Divide (h_jetInt_trk_pt_ref_sig_unf[iEvFrac][iPtJInt][iDir][3]);

          } // end loop over iCent

        } // end loop over nIter

      } // end loop over iDir

    } // end loop over iPtJInt

  } // end loop over iEvFrac
  std::cout << "Finished computing I_pPb ratios!" << std::endl;




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // COMPUTE CLOSURE HISTOGRAMS -- USE NITER = 2 RESULTS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::cout << "Computing closure of final results..." << std::endl;
  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
    //const TString funcStr = "[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)";
    std::cout << "|--> pTJInt = " << pTJInt << std::endl;

    //for (short iDir = 0; iDir < nDir; iDir++) {
    for (short iDir : {0, 2}) {

      const TString dir = directions[iDir];
      std::cout << "  |--> dir = " << dir << std::endl;

      h_jetInt_trk_pt_ref_sig_unf_fullClosure[iPtJInt][iDir] = (TH1D*) h_jetInt_trk_pt_ref_sig_unf[0][iPtJInt][iDir][1]->Clone (Form ("h_jetInt_trk_pt_%s_ref_sig_unf_%s_fullClosure", dir.Data (), pTJInt.Data ()));
      DivideNoErrors (h_jetInt_trk_pt_ref_sig_unf_fullClosure[iPtJInt][iDir], h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir][1]);

      h_jetInt_trk_pt_ref_sig_unf_fullClosure[iPtJInt][iDir]->Smooth (3);

      //TF1* f = DoClosureFit (Form ("f_jetInt_trk_pt_%s_ref_sig_unf_%s_fullClosure", dir.Data (), pTJInt.Data ()), h_jetInt_trk_pt_ref_sig_unf_fullClosure[iPtJInt][iDir], pTChBins[1], pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], 4, 8);
      //if (!f) std::cout << "Fit failed, please investigate! iPtJint = " << iPtJInt << ", iDir = " << iDir << std::endl;
      //f_jetInt_trk_pt_ref_sig_unf_fullClosure[iPtJInt][iDir] = f;


      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));
        std::cout << "  |--> cent = " << cent << std::endl;

        h_jetInt_trk_pt_sig_unf_fullClosure[iPtJInt][iDir][iCent] = (TH1D*) h_jetInt_trk_pt_sig_unf[0][iPtJInt][iDir][iCent][1]->Clone (Form ("h_jetInt_trk_pt_%s_%s_sig_unf_%s_fullClosure", dir.Data (), cent, pTJInt.Data ()));
        DivideNoErrors (h_jetInt_trk_pt_sig_unf_fullClosure[iPtJInt][iDir][iCent], h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent][1]);

        h_jetInt_trk_pt_sig_unf_fullClosure[iPtJInt][iDir][iCent]->Smooth (3);

        //TF1* f = DoClosureFit (Form ("f_jetInt_trk_pt_%s_%s_sig_unf_%s_fullClosure", dir.Data (), cent, pTJInt.Data ()), h_jetInt_trk_pt_sig_unf_fullClosure[iPtJInt][iDir][iCent], pTChBins[1], pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], 4, 8);
        //if (!f) std::cout << "Fit failed, please investigate! iPtJint = " << iPtJInt << ", iDir = " << iDir << ", iCent = " << iCent << std::endl;
        //f_jetInt_trk_pt_sig_unf_fullClosure[iPtJInt][iDir][iCent] = f;


        h_jetInt_trk_pt_iaa_unf_fullClosure[iPtJInt][iDir][iCent] = (TH1D*) h_jetInt_trk_pt_iaa_unf[0][iPtJInt][iDir][iCent][1]->Clone (Form ("h_jetInt_trk_pt_%s_%s_iaa_unf_%s_fullClosure", dir.Data (), cent, pTJInt.Data ()));
        DivideNoErrors (h_jetInt_trk_pt_iaa_unf_fullClosure[iPtJInt][iDir][iCent], h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent][1]);

        h_jetInt_trk_pt_iaa_unf_fullClosure[iPtJInt][iDir][iCent]->Smooth (3);

        //f = DoClosureFit (Form ("f_jetInt_trk_pt_%s_%s_iaa_unf_%s_fullClosure", dir.Data (), cent, pTJInt.Data ()), h_jetInt_trk_pt_iaa_unf_fullClosure[iPtJInt][iDir][iCent], pTChBins[1], pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], 4, 8);
        //if (!f) std::cout << "Fit failed, please investigate! iPtJint = " << iPtJInt << ", iDir = " << iDir << ", iCent = " << iCent << std::endl;
        //f_jetInt_trk_pt_iaa_unf_fullClosure[iPtJInt][iDir][iCent] = f;

      } // end loop over iCent
 
    } // end loop over iDir

  } // end loop over iPtJInt 
  std::cout << "Finished computing closure of final results!" << std::endl;



  {
    TFile* outFile = new TFile (Form ("%s/aux/MCClosureHists.root", workPath.Data ()), "recreate");

    for (short iPtJInt : {0, 1}) {

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        h_jetInt_trk_pt_ref_sig_unf_fullClosure[iPtJInt][iDir]->Write ();
        //f_jetInt_trk_pt_ref_sig_unf_fullClosure[iPtJInt][iDir]->Write ();

        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

          h_jetInt_trk_pt_sig_unf_fullClosure[iPtJInt][iDir][iCent]->Write ();
          //f_jetInt_trk_pt_sig_unf_fullClosure[iPtJInt][iDir][iCent]->Write ();
      
          h_jetInt_trk_pt_iaa_unf_fullClosure[iPtJInt][iDir][iCent]->Write ();
          //f_jetInt_trk_pt_iaa_unf_fullClosure[iPtJInt][iDir][iCent]->Write ();

        } // end loop over iCent
 
      } // end loop over iDir

    } // end loop over iPtJInt
  }
*/



  std::cout << "Finished main computational code, proceeded to plotting routines." << std::endl;
  {
    TLine* l = new TLine ();
    TLatex* tl = new TLatex ();


/*
    for (short iEvFrac : {1, 0}) {
  
      const TString evFrac = (iEvFrac == 0 ? "half" : "full");
  
      {
        const char* canvasName = Form ("c_jet_pt_%sClosure", evFrac.Data ());
        TCanvas* c = new TCanvas (canvasName, "", 1600, 800);
        c->Divide (4, 2);
  
        TH1D* h = nullptr;
        TGAE* g = nullptr;
  
        {
          c->cd (7);
  
          gPad->SetLogx ();
  
          h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];Ratio", 1, pTJBins[0], pTJBins[nPtJBins]);
          h->GetXaxis ()->SetMoreLogLabels ();
          //h->GetYaxis ()->SetRangeUser (0.0, 2.0);
          h->GetYaxis ()->SetRangeUser (0.6, 1.8);
          h->GetYaxis ()->CenterTitle ();
          h->SetBinContent (1, 1);
          h->SetLineStyle (2);
          h->SetLineWidth (2);
          h->SetLineColor (kBlack);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);
  
          l->SetLineWidth (2);
          l->SetLineColor (kGray+1);
          l->SetLineStyle (2);
          l->DrawLine (pTJBins[0], 1.05, pTJBins[nPtJBins], 1.05);
          l->DrawLine (pTJBins[0], 0.95, pTJBins[nPtJBins], 0.95);

          l->SetLineWidth (2);
          l->SetLineColor (kBlack);
          l->SetLineStyle (2);
          l->DrawLine (30, 0.6, 30, 1.6);
          l->DrawLine (300, 0.6, 300, 1.8);
 
          short iCol = 9;
          for (short nIter : {17, 15, 13, 11, 9, 7, 5, 3, 1}) {
          //short iCol = 3;
          //for (short nIter : {5, 3, 1}) {
            const Color_t col = cols[iCol--];
            g = make_graph (h_jet_pt_ref_unf[iEvFrac][nIter]);
            ScaleGraph (g, h_jet_pt_ref[iEvFrac][1]);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
          }
  
          g = make_graph (h_jet_pt_ref[iEvFrac][0]);
          ScaleGraph (g, h_jet_pt_ref[iEvFrac][1]);
          myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
  
          myText (0.2, 0.84, kBlack, "#bf{#it{pp}}", 0.06);
        }
  
        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
          c->cd (nFcalCentBins+1-iCent);
  
          gPad->SetLogx ();
  
          h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];Ratio", 1, pTJBins[0], pTJBins[nPtJBins]);
          h->GetXaxis ()->SetMoreLogLabels ();
          //h->GetYaxis ()->SetRangeUser (0.0, 2.0);
          h->GetYaxis ()->SetRangeUser (0.6, 1.8);
          h->GetYaxis ()->CenterTitle ();
          h->SetBinContent (1, 1);
          h->SetLineStyle (2);
          h->SetLineWidth (2);
          h->SetLineColor (kBlack);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);
  
          l->SetLineWidth (2);
          l->SetLineColor (kGray+1);
          l->SetLineStyle (2);
          l->DrawLine (pTJBins[0], 1.05, pTJBins[nPtJBins], 1.05);
          l->DrawLine (pTJBins[0], 0.95, pTJBins[nPtJBins], 0.95);

          l->SetLineWidth (2);
          l->SetLineColor (kBlack);
          l->SetLineStyle (2);
          l->DrawLine (30, 0.6, 30, 1.6);
          l->DrawLine (300, 0.6, 300, 1.8);
 
          short iCol = 9;
          for (short nIter : {17, 15, 13, 11, 9, 7, 5, 3, 1}) {
          //short iCol = 3;
          //for (short nIter : {5, 3, 1}) {
            const Color_t col = cols[iCol--];
            g = make_graph (h_jet_pt_unf[iEvFrac][iCent][nIter]);
            ScaleGraph (g, h_jet_pt[iEvFrac][iCent][1]);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
          }
  
          g = make_graph (h_jet_pt[iEvFrac][iCent][0]);
          ScaleGraph (g, h_jet_pt[iEvFrac][iCent][1]);
          myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
  
          if (iCent < nFcalCentBins)
            myText (0.2, 0.84, kBlack, Form ("#bf{#it{p}+Pb, FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
          else
            myText (0.2, 0.84, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);
  
        } // end loop over iCent
  
        c->cd (8);
        myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
        myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
        myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
        myLineText2 (0.15, 0.56, myCyan,        kFullCircle, "Raw / truth", 1.2, 0.06);
        myLineText2 (0.15, 0.48, myLiteBlue,    kOpenCircle, "2 iters.", 1.2, 0.06);
        myLineText2 (0.15, 0.40, myLitePurple,  kOpenCircle, "4 iters.", 1.2, 0.06);
        myLineText2 (0.15, 0.32, myPurple,      kOpenCircle, "6 iters.", 1.2, 0.06);
        myLineText2 (0.15, 0.24, myRed,         kOpenCircle, "8 iters.", 1.2, 0.06);
        myLineText2 (0.15, 0.16, myMaroon,      kOpenCircle, "10 iters.", 1.2, 0.06);
        myLineText2 (0.55, 0.48, myOrange,      kOpenCircle, "12 iters.", 1.2, 0.06);
        myLineText2 (0.55, 0.40, myLiteYellow,  kOpenCircle, "14 iters.", 1.2, 0.06);
        myLineText2 (0.55, 0.32, myLiteGreen,   kOpenCircle, "16 iters.", 1.2, 0.06);
        myLineText2 (0.55, 0.24, myGreen,       kOpenCircle, "18 iters.", 1.2, 0.06);
  
        c->SaveAs (Form ("%s/Plots/Unfolding/MC_JetPt_%sClosure.pdf", workPath.Data (), evFrac.Data ()));
  
      }
  
  
  

      for (short iPtJInt : {0, 1}) {
  
        const int minJetPt = (iPtJInt == 0 ? 30 : 60);
  
        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {
  
          const TString dir = directions[iDir];
  
          const char* canvasName = Form ("c_jet_trk_pt_%s_%iGeVJets_%sClosure_6Iters", dir.Data (), minJetPt, evFrac.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);
  
          TH1D* h = nullptr;
          TGAE* g = nullptr;
  
          {
            gPad->SetLogx ();
  
            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Unfolded yield / truth", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (0.9, 1.2);
            h->GetYaxis ()->CenterTitle ();
            h->SetBinContent (1, 1);
            h->SetLineStyle (2);
            h->SetLineWidth (2);
            h->SetLineColor (kBlack);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);
  
            l->SetLineWidth (2);
            l->SetLineColor (kGray+1);
            l->SetLineStyle (2);
            l->DrawLine (pTChBins[1], 1.05, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)], 1.05);
            l->DrawLine (pTChBins[1], 0.95, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)], 0.95);
  
            const Color_t col = colorfulColors[0];
            TH1D* h = (TH1D*) h_jetInt_trk_pt_ref_sig_unf[iEvFrac][iPtJInt][iDir][5]->Clone ("htemp");
            RebinSomeBins (&h, nRebinPtChBins, rebinPtChBins, true);
            g = make_graph (h);
            SaferDelete (&h);
            h = (TH1D*) h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][1]->Clone ("htemp");
            RebinSomeBins (&h, nRebinPtChBins, rebinPtChBins, true);
            ScaleGraph (g, h);
            SaferDelete (&h);
            TrimGraph (g, 0.5, iPtJInt == 0 ? 60 : 90);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "LP", false);
            SaferDelete (&g);
          }
  
          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
            const Color_t col = colorfulColors[iCent+1];
            TH1D* h = (TH1D*) h_jetInt_trk_pt_sig_unf[iEvFrac][iPtJInt][iDir][iCent][5]->Clone ("htemp");
            RebinSomeBins (&h, nRebinPtChBins, rebinPtChBins, true);
            g = make_graph (h);
            SaferDelete (&h);
            h = (TH1D*) h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][1]->Clone ("htemp");
            RebinSomeBins (&h, nRebinPtChBins, rebinPtChBins, true);
            ScaleGraph (g, h);
            SaferDelete (&h);
            TrimGraph (g, 0.5, iPtJInt == 0 ? 60 : 90);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "LP", false);
            SaferDelete (&g);
          } // end loop over iCent
  
          myText (0.58, 0.89, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.032);
          myText (0.58, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          myText (0.58, 0.81, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.032);
          myText (0.58, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, %s", minJetPt, (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side"))), 0.032);
          myText (0.58, 0.73, kBlack, "6 iterations", 0.032);

          myLineText2 (0.27, 0.25, colorfulColors[0], kOpenCircle, "#bf{#it{pp}}", 1.4, 0.032);
          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
            if (iCent < nFcalCentBins)
              myLineText2 (0.27+0.18*((iCent+1)/2), 0.25-0.04*((iCent+1)%2), colorfulColors[iCent+1], kOpenCircle, Form ("#bf{%i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 1.4, 0.032);
            else
              myLineText2 (0.27+0.18*((iCent+1)/2), 0.25-0.04*((iCent+1)%2), colorfulColors[iCent+1], kOpenCircle, "#bf{0-100%}", 1.4, 0.032);
          } // end loop over iCent
  
          c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_%s_%iGeVJets_%sClosure_6Iters.pdf", workPath.Data (), dir.Data (), minJetPt, evFrac.Data ()));
  
        } // end loop over iDir
  
      } // end loop over iPtJInt




      for (short iPtJInt : {0, 1}) {
  
        const int minJetPt = (iPtJInt == 0 ? 30 : 60);
  
        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {
  
          const TString dir = directions[iDir];
  
          const char* canvasName = Form ("c_jet_trk_pt_iaa_%s_%iGeVJets_%sClosure_6Iters", dir.Data (), minJetPt, evFrac.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);
  
          TH1D* h = nullptr;
          TGAE* g = nullptr;
  
          {
            gPad->SetLogx ();
  
            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Unfolded #it{I}_{#it{p}Pb} / truth", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (0.9, 1.2);
            h->GetYaxis ()->CenterTitle ();
            h->SetBinContent (1, 1);
            h->SetLineStyle (2);
            h->SetLineWidth (2);
            h->SetLineColor (kBlack);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);
  
            l->SetLineWidth (2);
            l->SetLineColor (kGray+1);
            l->SetLineStyle (2);
            l->DrawLine (pTChBins[1], 1.05, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)], 1.05);
            l->DrawLine (pTChBins[1], 0.95, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)], 0.95);
          }
  
          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
            const Color_t col = colorfulColors[iCent+1];
            TH1D* h = (TH1D*) h_jetInt_trk_pt_iaa_unf[iEvFrac][iPtJInt][iDir][iCent][5]->Clone ("htemp");
            RebinSomeBins (&h, nRebinPtChBins, rebinPtChBins, true);
            g = make_graph (h);
            SaferDelete (&h);
            h = (TH1D*) h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][1]->Clone ("htemp");
            RebinSomeBins (&h, nRebinPtChBins, rebinPtChBins, true);
            ScaleGraph (g, h);
            SaferDelete (&h);
            TrimGraph (g, 0.5, iPtJInt == 0 ? 60 : 90);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "LP", false);
            SaferDelete (&g);
          } // end loop over iCent
  
          myText (0.58, 0.89, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.032);
          myText (0.58, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
          myText (0.58, 0.81, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.032);
          myText (0.58, 0.77, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, %s", minJetPt, (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side"))), 0.032);
          myText (0.58, 0.73, kBlack, "6 iterations (#it{p}+Pb)", 0.032);
          myText (0.58, 0.69, kBlack, "4 iterations (#it{pp})", 0.032);

          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
            if (iCent < nFcalCentBins)
              myLineText2 (0.27+0.20*((iCent)/2), 0.25-0.04*((iCent)%2), colorfulColors[iCent+1], kOpenCircle, Form ("#bf{%i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 1.4, 0.032);
            else
              myLineText2 (0.27+0.20*((iCent)/2), 0.25-0.04*((iCent)%2), colorfulColors[iCent+1], kOpenCircle, "#bf{0-100%}", 1.4, 0.032);
          } // end loop over iCent
  
          c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_IAA_%s_%iGeVJets_%sClosure_6Iters.pdf", workPath.Data (), dir.Data (), minJetPt, evFrac.Data ()));
  
        } // end loop over iDir
  
      } // end loop over iPtJInt




      for (short iPtJInt : {0, 1}) {
  
        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
        const int minJetPt = (iPtJInt == 0 ? 30 : 60);
        const int maxJetPt = 300;

        const char* canvasName = Form ("c_njet_%iGeV_%sClosure", minJetPt, evFrac.Data ());
        TCanvas* c = new TCanvas (canvasName, "", 1600, 800);
        c->Divide (4, 2);
  
        TH1D* h = nullptr;
        TGAE* g = nullptr;
  
        {
          c->cd (7);

          h = h_jet_pt_ref[iEvFrac][1];
          const double njet_truth = h->Integral (h->FindBin (minJetPt+0.01), h->FindBin (maxJetPt-0.01));
          std::cout << "Truth nJet (pp): " << njet_truth << std::endl;
  
          h = new TH1D ("h", ";Iterations;N_{jet}^{unfolded} / N_{jet}^{truth}", 1, 0.5, maxIters+0.5);
          h->GetYaxis ()->SetRangeUser (0.90, 1.25);
          h->GetYaxis ()->CenterTitle ();
          h->SetBinContent (1, 1);
          h->SetLineStyle (2);
          h->SetLineWidth (2);
          h->SetLineColor (kBlack);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);
  
          l->SetLineWidth (2);
          l->SetLineColor (kGray+1);
          l->SetLineStyle (2);
          l->DrawLine (0.5, 1.05, maxIters+0.5, 1.05);
          l->DrawLine (0.5, 0.95, maxIters+0.5, 0.95);

          g = make_graph (h_njet_ref_unf[iEvFrac][iPtJInt]);
          ScaleGraph (g, nullptr, 1./njet_truth);
          myDraw (g, colorfulColors[0], kFullCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
  
          myText (0.2, 0.84, kBlack, "#bf{#it{pp}}", 0.06);
        }
  
        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

          c->cd (nFcalCentBins+1-iCent);

          h = h_jet_pt[iEvFrac][iCent][1];
          const double njet_truth = h->Integral (h->FindBin (minJetPt+0.01), h->FindBin (maxJetPt-0.01));
          std::cout << "Truth nJet (iCent = " << iCent << "): " << njet_truth << std::endl;
  
          h = new TH1D ("h", ";Iterations;N_{jet}^{unfolded} / N_{jet}^{truth}", 1, 0.5, maxIters+0.5);
          h->GetYaxis ()->SetRangeUser (0.90, 1.25);
          h->GetYaxis ()->CenterTitle ();
          h->SetBinContent (1, 1);
          h->SetLineStyle (2);
          h->SetLineWidth (2);
          h->SetLineColor (kBlack);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);
  
          l->SetLineWidth (2);
          l->SetLineColor (kGray+1);
          l->SetLineStyle (2);
          l->DrawLine (0.5, 1.05, maxIters+0.5, 1.05);
          l->DrawLine (0.5, 0.95, maxIters+0.5, 0.95);

          g = make_graph (h_njet_unf[iEvFrac][iPtJInt][iCent]);
          ScaleGraph (g, nullptr, 1./njet_truth);
          myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
  
          if (iCent < nFcalCentBins)
            myText (0.2, 0.84, kBlack, Form ("#bf{#it{p}+Pb, FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
          else
            myText (0.2, 0.84, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);
  
        } // end loop over iCent
  
        c->cd (8);
        myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
        myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
        myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
        myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, %i)", minJetPt, maxJetPt), 0.07);
  
        c->SaveAs (Form ("%s/Plots/Unfolding/MC_NJet_%s_%sClosure.pdf", workPath.Data (), pTJInt.Data (), evFrac.Data ()));
  
      } // end loop over iPtJInt




      for (short iPtJInt : {0, 1}) {
  
        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
        const int minJetPt = (iPtJInt == 0 ? 30 : 60);
        const int maxJetPt = 300;

        const char* canvasName = Form ("c_ptjet_%iGeV_%sClosure", minJetPt, evFrac.Data ());
        TCanvas* c = new TCanvas (canvasName, "", 1600, 800);
        c->Divide (4, 2);
  
        TH1D* h = nullptr;
        TGAE* g = nullptr;
  
        {
          c->cd (7);

          h = (TH1D*) h_jet_pt_ref[iEvFrac][1]->Clone ("htemp");
          h->GetXaxis ()->SetRange (iPtJInt == 0 ? 3 : 4, h->GetNbinsX () - 2);
          const double ptjet_truth = h->GetMean ();
          SaferDelete (&h);
          std::cout << "Truth <pT^jet> (pp): " << ptjet_truth << " GeV" << std::endl;
  
          h = new TH1D ("h", ";Iterations;#LT#it{p}_{T}^{jet}#GT^{unfolded} / #LT#it{p}_{T}^{jet}#GT^{truth}", 1, 0.5, maxIters+0.5);
          h->GetYaxis ()->SetRangeUser (0.90, 1.25);
          h->GetYaxis ()->CenterTitle ();
          h->SetBinContent (1, 1);
          h->SetLineStyle (2);
          h->SetLineWidth (2);
          h->SetLineColor (kBlack);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);
  
          l->SetLineWidth (2);
          l->SetLineColor (kGray+1);
          l->SetLineStyle (2);
          l->DrawLine (0.5, 1.05, maxIters+0.5, 1.05);
          l->DrawLine (0.5, 0.95, maxIters+0.5, 0.95);

          g = make_graph (h_ptjet_ref_unf[iEvFrac][iPtJInt]);
          ScaleGraph (g, nullptr, 1./ptjet_truth);
          myDraw (g, colorfulColors[0], kFullCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
  
          myText (0.2, 0.84, kBlack, "#bf{#it{pp}}", 0.06);
        }
  
        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

          c->cd (nFcalCentBins+1-iCent);

          h = (TH1D*) h_jet_pt[iEvFrac][iCent][1]->Clone ("htemp");
          h->GetXaxis ()->SetRange (iPtJInt == 0 ? 3 : 4, h->GetNbinsX () - 2);
          const double ptjet_truth = h->GetMean ();
          SaferDelete (&h);
          std::cout << "Truth <pT^jet> (iCent = " << iCent << "): " << ptjet_truth << " GeV" << std::endl;
  
          h = new TH1D ("h", ";Iterations;#LT#it{p}_{T}^{jet}#GT^{unfolded} / #LT#it{p}_{T}^{jet}#GT^{truth}", 1, 0.5, maxIters+0.5);
          h->GetYaxis ()->SetRangeUser (0.90, 1.25);
          h->GetYaxis ()->CenterTitle ();
          h->SetBinContent (1, 1);
          h->SetLineStyle (2);
          h->SetLineWidth (2);
          h->SetLineColor (kBlack);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);
  
          l->SetLineWidth (2);
          l->SetLineColor (kGray+1);
          l->SetLineStyle (2);
          l->DrawLine (0.5, 1.05, maxIters+0.5, 1.05);
          l->DrawLine (0.5, 0.95, maxIters+0.5, 0.95);

          g = make_graph (h_ptjet_unf[iEvFrac][iPtJInt][iCent]);
          ScaleGraph (g, nullptr, 1./ptjet_truth);
          myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
  
          if (iCent < nFcalCentBins)
            myText (0.2, 0.84, kBlack, Form ("#bf{#it{p}+Pb, FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
          else
            myText (0.2, 0.84, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);
  
        } // end loop over iCent
  
        c->cd (8);
        myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
        myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
        myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
        myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, %i)", minJetPt, maxJetPt), 0.07);
  
        c->SaveAs (Form ("%s/Plots/Unfolding/MC_AvgPtJet_%s_%sClosure.pdf", workPath.Data (), pTJInt.Data (), evFrac.Data ()));
  
      } // end loop over iPtJInt




      for (short iPtJInt : {0, 1}) {
  
        const int minJetPt = (iPtJInt == 0 ? 30 : 60);
  
        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {
  
          const TString dir = directions[iDir];
  
          const char* canvasName = Form ("c_jet_trk_pt_%s_%iGeVJets_%sClosure", dir.Data (), minJetPt, evFrac.Data ());
          //TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
          TCanvas* c = new TCanvas (canvasName, "", 1600, 800);
          c->Divide (4, 2);
  
          TH1D* h = nullptr;
          TGAE* g = nullptr;
  
          {
            c->cd (7);
  
            gPad->SetLogx ();
  
            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Ratio", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (0.7, 1.3);
            h->GetYaxis ()->CenterTitle ();
            h->SetBinContent (1, 1);
            h->SetLineStyle (2);
            h->SetLineWidth (2);
            h->SetLineColor (kBlack);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);
  
            l->SetLineWidth (2);
            l->SetLineColor (kGray+1);
            l->SetLineStyle (2);
            l->DrawLine (pTChBins[1], 1.05, pTChBins[nPtChBins-4], 1.05);
            l->DrawLine (pTChBins[1], 0.95, pTChBins[nPtChBins-4], 0.95);
  
            for (short nIter : {1, 3, 5}) {
            //for (short nIter : {1}) {
              const Color_t col = (nIter == 1 ? myLiteBlue : (nIter == 3 ? myLitePurple : myPurple));
              g = make_graph (h_jetInt_trk_pt_ref_sig_unf[iEvFrac][iPtJInt][iDir][nIter]);
              ScaleGraph (g, h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][1]);
              myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
              SaferDelete (&g);
            }
  
            //if (iEvFrac == 1) 
            //  myDraw (f_jetInt_trk_pt_ref_sig_unf_fullClosure[iPtJInt][iDir], kBlack, 2, 2);
  
            g = make_graph (h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][0]);
            ScaleGraph (g, h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][1]);
            myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
  
            myText (0.2, 0.84, kBlack, "#bf{#it{pp}}", 0.06);
          }
  
          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
            c->cd (nFcalCentBins+1-iCent);
  
            gPad->SetLogx ();
  
            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Ratio", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (0.7, 1.3);
            h->GetYaxis ()->CenterTitle ();
            h->SetBinContent (1, 1);
            h->SetLineStyle (2);
            h->SetLineWidth (2);
            h->SetLineColor (kBlack);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);
  
            l->SetLineWidth (2);
            l->SetLineColor (kGray+1);
            l->SetLineStyle (2);
            l->DrawLine (pTChBins[1], 1.05, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)], 1.05);
            l->DrawLine (pTChBins[1], 0.95, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)], 0.95);
  
            for (short nIter : {1, 3, 5}) {
            //for (short nIter : {1}) {
              const Color_t col = (nIter == 1 ? myLiteBlue : (nIter == 3 ? myLitePurple : myPurple));
              g = make_graph (h_jetInt_trk_pt_sig_unf[iEvFrac][iPtJInt][iDir][iCent][nIter]);
              ScaleGraph (g, h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][1]);
              myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
              SaferDelete (&g);
            }
  
            //if (iEvFrac == 1) 
            //  myDraw (f_jetInt_trk_pt_sig_unf_fullClosure[iPtJInt][iDir][iCent], kBlack, 2, 2);
  
            g = make_graph (h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][0]);
            ScaleGraph (g, h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][1]);
            myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
  
            if (iCent < nFcalCentBins)
              myText (0.2, 0.84, kBlack, Form ("#bf{#it{p}+Pb, FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
            else
              myText (0.2, 0.84, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);
  
          } // end loop over iCent
  
          c->cd (8);
          myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
          myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
          myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
          myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, %s", minJetPt, (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side"))), 0.07);
          myLineText2 (0.15, 0.48, myCyan,        kFullCircle, "Raw / truth", 1.2, 0.06);
          myLineText2 (0.15, 0.40, myLiteBlue,    kOpenCircle, "Unfolded / truth (2 iterations)", 1.2, 0.06);
          myLineText2 (0.15, 0.32, myLitePurple,  kOpenCircle, "Unfolded / truth (4 iterations)", 1.2, 0.06);
          myLineText2 (0.15, 0.24, myPurple,      kOpenCircle, "Unfolded / truth (6 iterations)", 1.2, 0.06);
  
          c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_%s_%iGeVJets_%sClosure.pdf", workPath.Data (), dir.Data (), minJetPt, evFrac.Data ()));
  
        } // end loop over iDir
  
      } // end loop over iPtJInt
  
  
  
  
      for (short iPtJInt : {0, 1}) {
  
        const int minJetPt = (iPtJInt == 0 ? 30 : 60);
  
        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {
  
          const TString dir = directions[iDir];
  
          const char* canvasName = Form ("c_jet_trk_pt_iaa_%s_%iGeVJets_%sClosure", dir.Data (), minJetPt, evFrac.Data ());
          //TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
          TCanvas* c = new TCanvas (canvasName, "", 1600, 800);
          c->Divide (4, 2);
  
          TH1D* h = nullptr;
          TGAE* g = nullptr;
  
          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
            c->cd (nFcalCentBins+1-iCent);
  
            gPad->SetLogx ();
  
            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Ratio", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
            h->GetXaxis ()->SetMoreLogLabels ();
            h->GetYaxis ()->SetRangeUser (0.7, 1.3);
            h->GetYaxis ()->CenterTitle ();
            h->SetBinContent (1, 1);
            h->SetLineStyle (2);
            h->SetLineWidth (2);
            h->SetLineColor (kBlack);
            h->DrawCopy ("hist ][");
            SaferDelete (&h);
  
            l->SetLineWidth (2);
            l->SetLineColor (kGray+1);
            l->SetLineStyle (2);
            l->DrawLine (pTChBins[1], 1.05, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)], 1.05);
            l->DrawLine (pTChBins[1], 0.95, pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)], 0.95);
  
            for (short nIter : {1, 3, 5}) {
            //for (short nIter : {1}) {
              const Color_t col = (nIter == 1 ? myLiteBlue : (nIter == 3 ? myLitePurple : myPurple));
              g = make_graph (h_jetInt_trk_pt_iaa_unf[iEvFrac][iPtJInt][iDir][iCent][nIter]);
              ScaleGraph (g, h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][1]);
              myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
              SaferDelete (&g);
            }
 
            //if (iEvFrac == 1)
            //  myDraw (f_jetInt_trk_pt_iaa_unf_fullClosure[iPtJInt][iDir][iCent], kBlack, 2, 2);
  
            g = make_graph (h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][0]);
            ScaleGraph (g, h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][1]);
            myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
  
            if (iCent < nFcalCentBins)
              myText (0.2, 0.84, kBlack, Form ("#bf{#it{p}+Pb, FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
            else
              myText (0.2, 0.84, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);
  
          } // end loop over iCent
  
          c->cd (7);
          myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
          myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
          myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
          myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, %s", minJetPt, (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side"))), 0.07);
          myLineText2 (0.15, 0.48, myCyan,        kFullCircle, "Raw / truth", 1.2, 0.06);
          myLineText2 (0.15, 0.40, myLiteBlue,    kOpenCircle, "Unfolded / truth (2 iterations)", 1.2, 0.06);
          myLineText2 (0.15, 0.32, myLitePurple,  kOpenCircle, "Unfolded / truth (4 iterations)", 1.2, 0.06);
          myLineText2 (0.15, 0.24, myPurple,      kOpenCircle, "Unfolded / truth (6 iterations)", 1.2, 0.06);
  
          c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_IAA_%s_%iGeVJets_%sClosure.pdf", workPath.Data (), dir.Data (), minJetPt, evFrac.Data ()));
  
        } // end loop over iDir
  
      } // end loop over iPtJInt

    } // end loop over iEvFrac
*/




    //for (short iDir = 0; iDir < nDir; iDir++) {
    for (short iDir : {0, 2}) {

      const TString dir = directions[iDir];

      const char* cname = Form ("c_2d_respMatrix_ref_%s", dir.Data ());

      TCanvas* c = new TCanvas (cname, "", 880, 800);
      const double lMargin = 0.17;
      const double rMargin = 0.11;
      const double bMargin = 0.17;
      const double tMargin = 0.04;

      c->SetLeftMargin (lMargin);
      c->SetRightMargin (rMargin);
      c->SetBottomMargin (bMargin);
      c->SetTopMargin (tMargin);

      c->SetLogz ();

      TH2D* h2 = (TH2D*) rooUnfResp_jet_trk_pt_ref_sig[2][iDir]->HresponseNoOverflow ()->Clone ("htemp");

      TAxis* xax = h2->GetXaxis ();
      TAxis* yax = h2->GetYaxis ();
      TAxis* zax = h2->GetZaxis ();

      xax->SetTitle ("Reco.-level #it{p}_{T}^{jet} [GeV]");
      yax->SetTitle ("Truth-level #it{p}_{T}^{jet} [GeV]");
      zax->SetTitle ("");

      xax->SetTitleFont (43);
      xax->SetTitleSize (32);
      yax->SetTitleFont (43);
      yax->SetTitleSize (32);
      zax->SetTitleFont (43);
      zax->SetTitleSize (32);
      xax->SetLabelFont (43);
      xax->SetLabelSize (0);
      yax->SetLabelFont (43);
      yax->SetLabelSize (0);
      zax->SetLabelFont (43);
      zax->SetLabelSize (24);

      xax->SetTitleOffset (2.0);
      yax->SetTitleOffset (2.0);

      xax->SetNdivisions (-40);
      yax->SetNdivisions (-40);

      h2->DrawCopy ("colz");
      SaferDelete (&h2);

      tl->SetTextFont (43);
      tl->SetTextSize (14);
      tl->SetTextAlign (33);
      tl->SetTextAngle (60);

      tl->DrawLatex (std::floor (0.5*nPtChBins), -6, "10-11");
      tl->DrawLatex (std::floor (1.5*nPtChBins), -6, "11-12");
      tl->DrawLatex (std::floor (2.5*nPtChBins), -6, "12-13");
      tl->DrawLatex (std::floor (3.5*nPtChBins), -6, "13-14");
      tl->DrawLatex (std::floor (4.5*nPtChBins), -6, "14-15");
      tl->DrawLatex (std::floor (5.5*nPtChBins), -6, "15-17.5");
      tl->DrawLatex (std::floor (6.5*nPtChBins), -6, "17.5-20");
      tl->DrawLatex (std::floor (7.5*nPtChBins), -6, "20-22.5");
      tl->DrawLatex (std::floor (8.5*nPtChBins), -6, "22.5-25");
      tl->DrawLatex (std::floor (9.5*nPtChBins), -6, "25-27.5");
      tl->DrawLatex (std::floor (10.5*nPtChBins), -6, "27.5-30");
      tl->DrawLatex (std::floor (11.5*nPtChBins), -6, "30-33");
      tl->DrawLatex (std::floor (12.5*nPtChBins), -6, "33-36");
      tl->DrawLatex (std::floor (13.5*nPtChBins), -6, "36-40");
      tl->DrawLatex (std::floor (14.5*nPtChBins), -6, "40-45");
      tl->DrawLatex (std::floor (15.5*nPtChBins), -6, "45-50");
      tl->DrawLatex (std::floor (16.5*nPtChBins), -6, "50-55");
      tl->DrawLatex (std::floor (17.5*nPtChBins), -6, "55-60");
      tl->DrawLatex (std::floor (18.5*nPtChBins), -6, "60-65");
      tl->DrawLatex (std::floor (19.5*nPtChBins), -6, "65-70");
      tl->DrawLatex (std::floor (20.5*nPtChBins), -6, "70-75");
      tl->DrawLatex (std::floor (21.5*nPtChBins), -6, "75-82.5");
      tl->DrawLatex (std::floor (22.5*nPtChBins), -6, "82.5-90");
      tl->DrawLatex (std::floor (23.5*nPtChBins), -6, "90-100");
      tl->DrawLatex (std::floor (24.5*nPtChBins), -6, "100-110");
      tl->DrawLatex (std::floor (25.5*nPtChBins), -6, "110-120");
      tl->DrawLatex (std::floor (26.5*nPtChBins), -6, "120-130");
      tl->DrawLatex (std::floor (27.5*nPtChBins), -6, "130-145");
      tl->DrawLatex (std::floor (28.5*nPtChBins), -6, "145-160");
      tl->DrawLatex (std::floor (29.5*nPtChBins), -6, "160-180");
      tl->DrawLatex (std::floor (30.5*nPtChBins), -6, "180-200");
      tl->DrawLatex (std::floor (31.5*nPtChBins), -6, "200-220");
      tl->DrawLatex (std::floor (32.5*nPtChBins), -6, "220-240");
      tl->DrawLatex (std::floor (33.5*nPtChBins), -6, "240-260");
      tl->DrawLatex (std::floor (34.5*nPtChBins), -6, "260-280");
      tl->DrawLatex (std::floor (35.5*nPtChBins), -6, "280-300");
      tl->DrawLatex (std::floor (36.5*nPtChBins), -6, "300-325");
      tl->DrawLatex (std::floor (37.5*nPtChBins), -6, "325-350");
      tl->DrawLatex (std::floor (38.5*nPtChBins), -6, "350-375");
      tl->DrawLatex (std::floor (39.5*nPtChBins), -6, "375-400");

      tl->SetTextAlign (32);
      tl->SetTextAngle (330);

      tl->DrawLatex (-6, std::floor (0.5*nPtChBins), "10-11");
      tl->DrawLatex (-6, std::floor (1.5*nPtChBins), "11-12");
      tl->DrawLatex (-6, std::floor (2.5*nPtChBins), "12-13");
      tl->DrawLatex (-6, std::floor (3.5*nPtChBins), "13-14");
      tl->DrawLatex (-6, std::floor (4.5*nPtChBins), "14-15");
      tl->DrawLatex (-6, std::floor (5.5*nPtChBins), "15-17.5");
      tl->DrawLatex (-6, std::floor (6.5*nPtChBins), "17.5-20");
      tl->DrawLatex (-6, std::floor (7.5*nPtChBins), "20-22.5");
      tl->DrawLatex (-6, std::floor (8.5*nPtChBins), "22.5-25");
      tl->DrawLatex (-6, std::floor (9.5*nPtChBins), "25-27.5");
      tl->DrawLatex (-6, std::floor (10.5*nPtChBins), "27.5-30");
      tl->DrawLatex (-6, std::floor (11.5*nPtChBins), "30-33");
      tl->DrawLatex (-6, std::floor (12.5*nPtChBins), "33-36");
      tl->DrawLatex (-6, std::floor (13.5*nPtChBins), "36-40");
      tl->DrawLatex (-6, std::floor (14.5*nPtChBins), "40-45");
      tl->DrawLatex (-6, std::floor (15.5*nPtChBins), "45-50");
      tl->DrawLatex (-6, std::floor (16.5*nPtChBins), "50-55");
      tl->DrawLatex (-6, std::floor (17.5*nPtChBins), "55-60");
      tl->DrawLatex (-6, std::floor (18.5*nPtChBins), "60-65");
      tl->DrawLatex (-6, std::floor (19.5*nPtChBins), "65-70");
      tl->DrawLatex (-6, std::floor (20.5*nPtChBins), "70-75");
      tl->DrawLatex (-6, std::floor (21.5*nPtChBins), "75-82.5");
      tl->DrawLatex (-6, std::floor (22.5*nPtChBins), "82.5-90");
      tl->DrawLatex (-6, std::floor (23.5*nPtChBins), "90-100");
      tl->DrawLatex (-6, std::floor (24.5*nPtChBins), "100-110");
      tl->DrawLatex (-6, std::floor (25.5*nPtChBins), "110-120");
      tl->DrawLatex (-6, std::floor (26.5*nPtChBins), "120-130");
      tl->DrawLatex (-6, std::floor (27.5*nPtChBins), "130-145");
      tl->DrawLatex (-6, std::floor (28.5*nPtChBins), "145-160");
      tl->DrawLatex (-6, std::floor (29.5*nPtChBins), "160-180");
      tl->DrawLatex (-6, std::floor (30.5*nPtChBins), "180-200");
      tl->DrawLatex (-6, std::floor (31.5*nPtChBins), "200-220");
      tl->DrawLatex (-6, std::floor (32.5*nPtChBins), "220-240");
      tl->DrawLatex (-6, std::floor (33.5*nPtChBins), "240-260");
      tl->DrawLatex (-6, std::floor (34.5*nPtChBins), "260-280");
      tl->DrawLatex (-6, std::floor (35.5*nPtChBins), "280-300");
      tl->DrawLatex (-6, std::floor (36.5*nPtChBins), "300-325");
      tl->DrawLatex (-6, std::floor (37.5*nPtChBins), "325-350");
      tl->DrawLatex (-6, std::floor (38.5*nPtChBins), "350-375");
      tl->DrawLatex (-6, std::floor (39.5*nPtChBins), "375-400");

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      tl->SetTextAngle (0);
      
      tl->SetTextColor (kBlack);
      tl->SetTextAlign (12);

      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side")));

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kGreen+2);
      l->DrawLine (nPtChBins*0,  nPtChBins*11, nPtChBins*36, nPtChBins*11);
      l->DrawLine (nPtChBins*11, nPtChBins*0,  nPtChBins*11, nPtChBins*36);

      l->SetLineColor (kMagenta+2);
      l->DrawLine (nPtChBins*0,  nPtChBins*18, nPtChBins*36, nPtChBins*18);
      l->DrawLine (nPtChBins*18, nPtChBins*0,  nPtChBins*18, nPtChBins*36);

      l->SetLineColor (kBlack);
      l->DrawLine (nPtChBins*36, nPtChBins*11, nPtChBins*36, nPtChBins*36);
      l->DrawLine (nPtChBins*11, nPtChBins*36, nPtChBins*36, nPtChBins*36);

      c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_ResponseMatrix_ref_2D_%s.png", workPath.Data (), dir.Data ()));
      c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_ResponseMatrix_ref_2D_%s.pdf", workPath.Data (), dir.Data ()));
    } // end loop over iDir




    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        const char* cname = Form ("c_2d_respMatrix_%s_%s", dir.Data (), cent);

        TCanvas* c = new TCanvas (cname, "", 880, 800);
        const double lMargin = 0.17;
        const double rMargin = 0.11;
        const double bMargin = 0.17;
        const double tMargin = 0.04;

        c->SetLeftMargin (lMargin);
        c->SetRightMargin (rMargin);
        c->SetBottomMargin (bMargin);
        c->SetTopMargin (tMargin);

        c->SetLogz ();

        TH2D* h2 = (TH2D*) rooUnfResp_jet_trk_pt_sig[2][iDir][iCent]->HresponseNoOverflow ()->Clone ("htemp");

        TAxis* xax = h2->GetXaxis ();
        TAxis* yax = h2->GetYaxis ();
        TAxis* zax = h2->GetZaxis ();

        xax->SetTitle ("Reco.-level #it{p}_{T}^{jet} [GeV]");
        yax->SetTitle ("Truth-level #it{p}_{T}^{jet} [GeV]");
        zax->SetTitle ("");

        xax->SetTitleFont (43);
        xax->SetTitleSize (32);
        yax->SetTitleFont (43);
        yax->SetTitleSize (32);
        zax->SetTitleFont (43);
        zax->SetTitleSize (32);
        xax->SetLabelFont (43);
        xax->SetLabelSize (0);
        yax->SetLabelFont (43);
        yax->SetLabelSize (0);
        zax->SetLabelFont (43);
        zax->SetLabelSize (24);

        xax->SetTitleOffset (2.0);
        yax->SetTitleOffset (2.0);

        xax->SetNdivisions (-40);
        yax->SetNdivisions (-40);

        h2->DrawCopy ("colz");
        SaferDelete (&h2);

        tl->SetTextFont (43);
        tl->SetTextSize (14);
        tl->SetTextAlign (33);
        tl->SetTextAngle (60);

        tl->DrawLatex (std::floor (0.5*nPtChBins), -6, "10-11");
        tl->DrawLatex (std::floor (1.5*nPtChBins), -6, "11-12");
        tl->DrawLatex (std::floor (2.5*nPtChBins), -6, "12-13");
        tl->DrawLatex (std::floor (3.5*nPtChBins), -6, "13-14");
        tl->DrawLatex (std::floor (4.5*nPtChBins), -6, "14-15");
        tl->DrawLatex (std::floor (5.5*nPtChBins), -6, "15-17.5");
        tl->DrawLatex (std::floor (6.5*nPtChBins), -6, "17.5-20");
        tl->DrawLatex (std::floor (7.5*nPtChBins), -6, "20-22.5");
        tl->DrawLatex (std::floor (8.5*nPtChBins), -6, "22.5-25");
        tl->DrawLatex (std::floor (9.5*nPtChBins), -6, "25-27.5");
        tl->DrawLatex (std::floor (10.5*nPtChBins), -6, "27.5-30");
        tl->DrawLatex (std::floor (11.5*nPtChBins), -6, "30-33");
        tl->DrawLatex (std::floor (12.5*nPtChBins), -6, "33-36");
        tl->DrawLatex (std::floor (13.5*nPtChBins), -6, "36-40");
        tl->DrawLatex (std::floor (14.5*nPtChBins), -6, "40-45");
        tl->DrawLatex (std::floor (15.5*nPtChBins), -6, "45-50");
        tl->DrawLatex (std::floor (16.5*nPtChBins), -6, "50-55");
        tl->DrawLatex (std::floor (17.5*nPtChBins), -6, "55-60");
        tl->DrawLatex (std::floor (18.5*nPtChBins), -6, "60-65");
        tl->DrawLatex (std::floor (19.5*nPtChBins), -6, "65-70");
        tl->DrawLatex (std::floor (20.5*nPtChBins), -6, "70-75");
        tl->DrawLatex (std::floor (21.5*nPtChBins), -6, "75-82.5");
        tl->DrawLatex (std::floor (22.5*nPtChBins), -6, "82.5-90");
        tl->DrawLatex (std::floor (23.5*nPtChBins), -6, "90-100");
        tl->DrawLatex (std::floor (24.5*nPtChBins), -6, "100-110");
        tl->DrawLatex (std::floor (25.5*nPtChBins), -6, "110-120");
        tl->DrawLatex (std::floor (26.5*nPtChBins), -6, "120-130");
        tl->DrawLatex (std::floor (27.5*nPtChBins), -6, "130-145");
        tl->DrawLatex (std::floor (28.5*nPtChBins), -6, "145-160");
        tl->DrawLatex (std::floor (29.5*nPtChBins), -6, "160-180");
        tl->DrawLatex (std::floor (30.5*nPtChBins), -6, "180-200");
        tl->DrawLatex (std::floor (31.5*nPtChBins), -6, "200-220");
        tl->DrawLatex (std::floor (32.5*nPtChBins), -6, "220-240");
        tl->DrawLatex (std::floor (33.5*nPtChBins), -6, "240-260");
        tl->DrawLatex (std::floor (34.5*nPtChBins), -6, "260-280");
        tl->DrawLatex (std::floor (35.5*nPtChBins), -6, "280-300");
        tl->DrawLatex (std::floor (36.5*nPtChBins), -6, "300-325");
        tl->DrawLatex (std::floor (37.5*nPtChBins), -6, "325-350");
        tl->DrawLatex (std::floor (38.5*nPtChBins), -6, "350-375");
        tl->DrawLatex (std::floor (39.5*nPtChBins), -6, "375-400");

        tl->SetTextAlign (32);
        tl->SetTextAngle (330);

        tl->DrawLatex (-6, std::floor (0.5*nPtChBins), "10-11");
        tl->DrawLatex (-6, std::floor (1.5*nPtChBins), "11-12");
        tl->DrawLatex (-6, std::floor (2.5*nPtChBins), "12-13");
        tl->DrawLatex (-6, std::floor (3.5*nPtChBins), "13-14");
        tl->DrawLatex (-6, std::floor (4.5*nPtChBins), "14-15");
        tl->DrawLatex (-6, std::floor (5.5*nPtChBins), "15-17.5");
        tl->DrawLatex (-6, std::floor (6.5*nPtChBins), "17.5-20");
        tl->DrawLatex (-6, std::floor (7.5*nPtChBins), "20-22.5");
        tl->DrawLatex (-6, std::floor (8.5*nPtChBins), "22.5-25");
        tl->DrawLatex (-6, std::floor (9.5*nPtChBins), "25-27.5");
        tl->DrawLatex (-6, std::floor (10.5*nPtChBins), "27.5-30");
        tl->DrawLatex (-6, std::floor (11.5*nPtChBins), "30-33");
        tl->DrawLatex (-6, std::floor (12.5*nPtChBins), "33-36");
        tl->DrawLatex (-6, std::floor (13.5*nPtChBins), "36-40");
        tl->DrawLatex (-6, std::floor (14.5*nPtChBins), "40-45");
        tl->DrawLatex (-6, std::floor (15.5*nPtChBins), "45-50");
        tl->DrawLatex (-6, std::floor (16.5*nPtChBins), "50-55");
        tl->DrawLatex (-6, std::floor (17.5*nPtChBins), "55-60");
        tl->DrawLatex (-6, std::floor (18.5*nPtChBins), "60-65");
        tl->DrawLatex (-6, std::floor (19.5*nPtChBins), "65-70");
        tl->DrawLatex (-6, std::floor (20.5*nPtChBins), "70-75");
        tl->DrawLatex (-6, std::floor (21.5*nPtChBins), "75-82.5");
        tl->DrawLatex (-6, std::floor (22.5*nPtChBins), "82.5-90");
        tl->DrawLatex (-6, std::floor (23.5*nPtChBins), "90-100");
        tl->DrawLatex (-6, std::floor (24.5*nPtChBins), "100-110");
        tl->DrawLatex (-6, std::floor (25.5*nPtChBins), "110-120");
        tl->DrawLatex (-6, std::floor (26.5*nPtChBins), "120-130");
        tl->DrawLatex (-6, std::floor (27.5*nPtChBins), "130-145");
        tl->DrawLatex (-6, std::floor (28.5*nPtChBins), "145-160");
        tl->DrawLatex (-6, std::floor (29.5*nPtChBins), "160-180");
        tl->DrawLatex (-6, std::floor (30.5*nPtChBins), "180-200");
        tl->DrawLatex (-6, std::floor (31.5*nPtChBins), "200-220");
        tl->DrawLatex (-6, std::floor (32.5*nPtChBins), "220-240");
        tl->DrawLatex (-6, std::floor (33.5*nPtChBins), "240-260");
        tl->DrawLatex (-6, std::floor (34.5*nPtChBins), "260-280");
        tl->DrawLatex (-6, std::floor (35.5*nPtChBins), "280-300");
        tl->DrawLatex (-6, std::floor (36.5*nPtChBins), "300-325");
        tl->DrawLatex (-6, std::floor (37.5*nPtChBins), "325-350");
        tl->DrawLatex (-6, std::floor (38.5*nPtChBins), "350-375");
        tl->DrawLatex (-6, std::floor (39.5*nPtChBins), "375-400");

        tl->SetTextFont (43);
        tl->SetTextSize (32);
        tl->SetTextAlign (21);
        tl->SetTextAngle (0);
        
        tl->SetTextColor (kBlack);
        tl->SetTextAlign (12);

        tl->SetTextSize (26);
        tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
        tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp} + #it{p}+Pb overlay, #sqrt{s} = 5.02 TeV");
        tl->DrawLatexNDC (0.26, 0.810, iCent == nZdcCentBins ? "0-100%" : Form ("FCal %i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]));
        tl->DrawLatexNDC (0.26, 0.770, (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side")));

        l->SetLineStyle (2);
        l->SetLineWidth (2);
        l->SetLineColor (kGreen+2);
        l->DrawLine (nPtChBins*0,  nPtChBins*11, nPtChBins*36, nPtChBins*11);
        l->DrawLine (nPtChBins*11, nPtChBins*0,  nPtChBins*11, nPtChBins*36);

        l->SetLineColor (kMagenta+2);
        l->DrawLine (nPtChBins*0,  nPtChBins*18, nPtChBins*36, nPtChBins*18);
        l->DrawLine (nPtChBins*18, nPtChBins*0,  nPtChBins*18, nPtChBins*36);

        l->SetLineColor (kBlack);
        l->DrawLine (nPtChBins*36, nPtChBins*11, nPtChBins*36, nPtChBins*36);
        l->DrawLine (nPtChBins*11, nPtChBins*36, nPtChBins*36, nPtChBins*36);

        c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_ResponseMatrix_%s_2D_%s.png", workPath.Data (), cent, dir.Data ()));
        c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_ResponseMatrix_%s_2D_%s.pdf", workPath.Data (), cent, dir.Data ()));

      } // end loop over iDir

    } // end loop over iCent




    {
      const char* cname = Form ("c_2d_respMatrix_ref");

      TCanvas* c = new TCanvas (cname, "", 880, 800);
      const double lMargin = 0.17;
      const double rMargin = 0.11;
      const double bMargin = 0.17;
      const double tMargin = 0.04;

      c->SetLeftMargin (lMargin);
      c->SetRightMargin (rMargin);
      c->SetBottomMargin (bMargin);
      c->SetTopMargin (tMargin);

      c->SetLogx ();
      c->SetLogy ();
      c->SetLogz ();

      TH2D* h2 = (TH2D*) rooUnfResp_jet_pt_ref[2]->HresponseNoOverflow ()->Clone ("htemp");

      TAxis* xax = h2->GetXaxis ();
      TAxis* yax = h2->GetYaxis ();
      TAxis* zax = h2->GetZaxis ();

      xax->SetTitle ("Reco.-level #it{p}_{T}^{jet} [GeV]");
      yax->SetTitle ("Truth-level #it{p}_{T}^{jet} [GeV]");
      zax->SetTitle ("");

      xax->SetMoreLogLabels ();
      yax->SetMoreLogLabels ();

      xax->SetTitleFont (43);
      xax->SetTitleSize (32);
      yax->SetTitleFont (43);
      yax->SetTitleSize (32);
      zax->SetTitleFont (43);
      zax->SetTitleSize (32);
      xax->SetLabelFont (43);
      xax->SetLabelSize (32);
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      zax->SetLabelFont (43);
      zax->SetLabelSize (24);

      xax->SetTitleOffset (2.0);
      yax->SetTitleOffset (2.0);

      xax->SetNdivisions (-12);
      yax->SetNdivisions (-12);

      h2->DrawCopy ("colz");
      SaferDelete (&h2);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      tl->SetTextAngle (0);
      
      tl->SetTextColor (kBlack);
      tl->SetTextAlign (12);

      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kGreen+2);
      l->DrawLine (30, 10, 30, 300);
      l->DrawLine (10, 30, 300, 30);

      l->SetLineColor (kMagenta+2);
      l->DrawLine (60, 10, 60, 300);
      l->DrawLine (10, 60, 300, 60);

      l->SetLineColor (kBlack);
      l->DrawLine (300, 30, 300, 300);
      l->DrawLine (30, 300, 300, 300);

      c->SaveAs (Form ("%s/Plots/Unfolding/MC_JetPt_ResponseMatrix_ref_2D.png", workPath.Data ()));
      c->SaveAs (Form ("%s/Plots/Unfolding/MC_JetPt_ResponseMatrix_ref_2D.pdf", workPath.Data ()));
    }




    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

      const char* cname = Form ("c_2d_respMatrix_%s", cent);

      TCanvas* c = new TCanvas (cname, "", 880, 800);
      const double lMargin = 0.17;
      const double rMargin = 0.11;
      const double bMargin = 0.17;
      const double tMargin = 0.04;

      c->SetLeftMargin (lMargin);
      c->SetRightMargin (rMargin);
      c->SetBottomMargin (bMargin);
      c->SetTopMargin (tMargin);

      c->SetLogx ();
      c->SetLogy ();
      c->SetLogz ();

      TH2D* h2 = (TH2D*) rooUnfResp_jet_pt[2][iCent]->HresponseNoOverflow ()->Clone ("htemp");

      TAxis* xax = h2->GetXaxis ();
      TAxis* yax = h2->GetYaxis ();
      TAxis* zax = h2->GetZaxis ();

      xax->SetTitle ("Reco.-level #it{p}_{T}^{jet} [GeV]");
      yax->SetTitle ("Truth-level #it{p}_{T}^{jet} [GeV]");
      zax->SetTitle ("");

      xax->SetMoreLogLabels ();
      yax->SetMoreLogLabels ();

      xax->SetTitleFont (43);
      xax->SetTitleSize (32);
      yax->SetTitleFont (43);
      yax->SetTitleSize (32);
      zax->SetTitleFont (43);
      zax->SetTitleSize (32);
      xax->SetLabelFont (43);
      xax->SetLabelSize (32);
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      zax->SetLabelFont (43);
      zax->SetLabelSize (24);

      xax->SetTitleOffset (2.0);
      yax->SetTitleOffset (2.0);

      xax->SetNdivisions (-12);
      yax->SetNdivisions (-12);

      h2->DrawCopy ("colz");
      SaferDelete (&h2);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      tl->SetTextAngle (0);
      
      tl->SetTextColor (kBlack);
      tl->SetTextAlign (12);

      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp} + #it{p}+Pb overlay, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, iCent == nZdcCentBins ? "0-100%" : Form ("FCal %i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]));

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kGreen+2);
      l->DrawLine (30, 10, 30, 300);
      l->DrawLine (10, 30, 300, 30);

      l->SetLineColor (kMagenta+2);
      l->DrawLine (60, 10, 60, 300);
      l->DrawLine (10, 60, 300, 60);

      l->SetLineColor (kBlack);
      l->DrawLine (300, 30, 300, 300);
      l->DrawLine (30, 300, 300, 300);

      c->SaveAs (Form ("%s/Plots/Unfolding/MC_JetPt_ResponseMatrix_%s_2D.png", workPath.Data (), cent));
      c->SaveAs (Form ("%s/Plots/Unfolding/MC_JetPt_ResponseMatrix_%s_2D.pdf", workPath.Data (), cent));
    } // end loop over iCent


    /*
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  
      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {
  
        const TString dir = directions[iDir];
  
        const char* canvasName = Form ("c_jet_trk_pt_%s_%g-%gGeVJets_%sClosure", dir.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1]);
        TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
        c->Divide (4, 2);
  
        TH1D* h = nullptr;
        TGAE* g = nullptr;
  
        {
          c->cd (7);
  
          gPad->SetLogx ();
  
          h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Ratio", 1, pTChBins[0], pTChBins[nPtChBins]);
          h->GetXaxis ()->SetMoreLogLabels ();
          h->GetYaxis ()->SetRangeUser (0.7, 1.3);
          h->GetYaxis ()->CenterTitle ();
          h->SetBinContent (1, 1);
          h->SetLineStyle (2);
          h->SetLineWidth (2);
          h->SetLineColor (kBlack);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);
  
          l->SetLineWidth (2);
          l->SetLineColor (kGray+1);
          l->SetLineStyle (2);
          l->DrawLine (pTChBins[0], 1.05, pTChBins[nPtChBins], 1.05);
          l->DrawLine (pTChBins[0], 0.95, pTChBins[nPtChBins], 0.95);
  
          for (short nIter : {1, 3, 5}) {
            const Color_t col = (nIter == 1 ? myLiteBlue : (nIter == 3 ? myLitePurple : myPurple));
            g = make_graph (h_jet_trk_pt_ref_sig_unf[1][iPtJ][iDir][nIter]);
            ScaleGraph (g, h_jet_trk_pt_ref_sig[1][iPtJ][iDir][1]);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
          }
  
          g = make_graph (h_jet_trk_pt_ref_sig[1][iPtJ][iDir][0]);
          ScaleGraph (g, h_jet_trk_pt_ref_sig[1][iPtJ][iDir][1]);
          myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
  
          myText (0.2, 0.84, kBlack, "#bf{#it{pp}}", 0.06);
        }
  
        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
          c->cd (nFcalCentBins+1-iCent);
  
          gPad->SetLogx ();
  
          h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Ratio", 1, pTChBins[0], pTChBins[nPtChBins]);
          h->GetXaxis ()->SetMoreLogLabels ();
          h->GetYaxis ()->SetRangeUser (0.7, 1.3);
          h->GetYaxis ()->CenterTitle ();
          h->SetBinContent (1, 1);
          h->SetLineStyle (2);
          h->SetLineWidth (2);
          h->SetLineColor (kBlack);
          h->DrawCopy ("hist ][");
          SaferDelete (&h);
  
          l->SetLineWidth (2);
          l->SetLineColor (kGray+1);
          l->SetLineStyle (2);
          l->DrawLine (pTChBins[0], 1.05, pTChBins[nPtChBins], 1.05);
          l->DrawLine (pTChBins[0], 0.95, pTChBins[nPtChBins], 0.95);
  
          for (short nIter : {1, 3, 5}) {
            const Color_t col = (nIter == 1 ? myLiteBlue : (nIter == 3 ? myLitePurple : myPurple));
            g = make_graph (h_jet_trk_pt_sig_unf[1][iPtJ][iDir][iCent][nIter]);
            ScaleGraph (g, h_jet_trk_pt_sig[1][iPtJ][iDir][iCent][1]);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
          }
  
          g = make_graph (h_jet_trk_pt_sig[1][iPtJ][iDir][iCent][0]);
          ScaleGraph (g, h_jet_trk_pt_sig[1][iPtJ][iDir][iCent][1]);
          myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
  
          if (iCent < nFcalCentBins)
            myText (0.2, 0.84, kBlack, Form ("#bf{FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
          else {
            myText (0.2, 0.84, kBlack, "#bf{0-100%}", 0.06);
            myText (0.2, 0.76, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.06);
          }
  
        } // end loop over iCent
  
        c->cd (8);
        myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
        myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
        myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
        myText (0.1, 0.57, kBlack, Form ("%g < #it{p}_{T}^{jet} < %g GeV, %s", pTJBins[iPtJ], pTJBins[iPtJ+1], (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side"))), 0.07);
        myLineText2 (0.15, 0.48, myCyan,        kFullCircle, "Raw / truth", 1.2, 0.06);
        myLineText2 (0.15, 0.40, myLiteBlue,    kOpenCircle, "Unfolded / truth (2 iterations)", 1.2, 0.06);
        myLineText2 (0.15, 0.32, myLitePurple,  kOpenCircle, "Unfolded / truth (4 iterations)", 1.2, 0.06);
        myLineText2 (0.15, 0.24, myPurple,      kOpenCircle, "Unfolded / truth (6 iterations)", 1.2, 0.06);
  
        c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_%s_%g-%gGeVJets_FullClosure.pdf", workPath.Data (), dir.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1]));
  
      } // end loop over iDir
  
    } // end loop over iPtJ
    */

  }

}

#endif
