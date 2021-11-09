#ifndef __JetHadronCorrelations_PlotResponseMatrix_cxx__
#define __JetHadronCorrelations_PlotResponseMatrix_cxx__

#include "Params.h"
#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"
#include "Process.h"
#include "PrimaryFractionFit.h"

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

PrimaryFractionFit polyFunc;


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
  }

  if (succeeded)
    return f;
  else {
    delete f;
    return nullptr;
  }
}


void PlotResponseMatrix () {

  //const double pTJBins[] = {20, 30, 45, 60, 80, 100, 130, 160, 200, 240, 320};
  //const short nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;
  const Color_t cols[10] = {myCyan, myLiteBlue, myLitePurple, myPurple, myRed, myMaroon, myOrange, myLiteYellow, myLiteGreen, myGreen};

  const short maxIters = 20;

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

  // unfolded UE subtracted particle yield histograms
  TH2D****      h2_jet_trk_pt_ref_sig_unf         = Get3DArray <TH2D*> (2, nDir, maxIters);
  TH2D*****     h2_jet_trk_pt_sig_unf             = Get4DArray <TH2D*> (2, nDir, nFcalCentBins+1, maxIters);
  TH1D*****     h_jet_trk_pt_ref_sig_unf          = Get4DArray <TH1D*> (2, nPtJBins, nDir, maxIters);
  TH1D******    h_jet_trk_pt_sig_unf              = Get5DArray <TH1D*> (2, nPtJBins, nDir, nFcalCentBins+1, maxIters);
  TH1D*****     h_jetInt_trk_pt_ref_sig_unf       = Get4DArray <TH1D*> (2, 2, nDir, maxIters);
  TH1D******    h_jetInt_trk_pt_sig_unf           = Get5DArray <TH1D*> (2, 2, nDir, nFcalCentBins+1, maxIters);
  TH1D******    h_jetInt_trk_pt_iaa_unf           = Get5DArray <TH1D*> (2, 2, nDir, nFcalCentBins+1, maxIters);


  RooUnfoldResponse**   rooUnfResp_jet_pt_ref         = Get1DArray <RooUnfoldResponse*> (2);
  RooUnfoldResponse***  rooUnfResp_jet_pt             = Get2DArray <RooUnfoldResponse*> (2, nFcalCentBins+1);
  RooUnfoldResponse***  rooUnfResp_jet_trk_pt_ref_sig = Get2DArray <RooUnfoldResponse*> (2, nDir);
  RooUnfoldResponse**** rooUnfResp_jet_trk_pt_sig     = Get3DArray <RooUnfoldResponse*> (2, nDir, nFcalCentBins+1);


  TH1D***   h_jetInt_trk_pt_ref_sig_unf_halfClosure = Get2DArray <TH1D*> (2, nDir);
  TF1***    f_jetInt_trk_pt_ref_sig_unf_halfClosure = Get2DArray <TF1*>  (2, nDir);
  TH1D****  h_jetInt_trk_pt_sig_unf_halfClosure     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);
  TF1****   f_jetInt_trk_pt_sig_unf_halfClosure     = Get3DArray <TF1*>  (2, nDir, nZdcCentBins+1);
  TH1D****  h_jetInt_trk_pt_iaa_unf_halfClosure     = Get3DArray <TH1D*> (2, nDir, nZdcCentBins+1);
  TF1****   f_jetInt_trk_pt_iaa_unf_halfClosure     = Get3DArray <TF1*>  (2, nDir, nZdcCentBins+1);


  {
    TFile* inFile = new TFile (Form ("%s/MakeResponseMatrix/Nominal/allSamples.root", rootPath.Data ()), "read");
    //TFile* inFile = new TFile (Form ("%s/MakeResponseMatrix/Nominal/allSamples_primTracksOnly.root", rootPath.Data ()), "read");

    
    //////////////////////////////////////////////////////////////////////////////////////////////////// 
    // LOAD HALF-CLOSURE TRAINING AND TEST HISTOGRAMS
    //////////////////////////////////////////////////////////////////////////////////////////////////// 
    for (short iEvFrac : {0, 1}) {

      const TString evFrac = (iEvFrac == 0 ? "half" : "full");

      for (short iVar = 0; iVar < 2; iVar++) {

        const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");

        h_jet_pt_ref[iEvFrac][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%sClosure_ref_mc_%s", evFrac.Data (), var.Data ()));

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];
      
          h2_jet_trk_pt_ref_tot[iEvFrac][iDir][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_%sClosure_%s_ref_sig_mc_%s", evFrac.Data (), dir.Data (), var.Data ()))->Clone (Form ("h2_jet_trk_pt_%s_ref_tot_mc_%s", dir.Data (), var.Data ()));

        } // end loop over iDir


        for (short iCent = 0; iCent < nFcalCentBins; iCent++) {

          const char* cent = Form ("pPb_iCent%i", iCent);

          h_jet_pt[iEvFrac][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%sClosure_%s_mc_%s", evFrac.Data (), cent, var.Data ()));

          for (short iDir = 0; iDir < nDir; iDir++) {

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

          for (short iDir = 0; iDir < nDir; iDir++) {

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
    // LOAD HALF-CLOSURE RESPONSE MATRICES
    //////////////////////////////////////////////////////////////////////////////////////////////////// 
    for (short iEvFrac : {0, 1}) {

      const TString evFrac = (iEvFrac == 0 ? "half" : "full");

      rooUnfResp_jet_pt_ref[iEvFrac] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_ref_mc_%sClosure", evFrac.Data ()));

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        rooUnfResp_jet_trk_pt_ref_sig[iEvFrac][iDir] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_ref_sig_mc_%sClosure", dir.Data (), evFrac.Data ()));

      } // end loop over iDir

      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

        rooUnfResp_jet_pt[iEvFrac][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_%sClosure", cent, evFrac.Data ()));

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          rooUnfResp_jet_trk_pt_sig[iEvFrac][iDir][iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_trk_pt_%s_%s_sig_mc_%sClosure", dir.Data (), cent, evFrac.Data ()));

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

  //      for (short iDir = 0; iDir < nDir; iDir++) {

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

  //      for (short iDir = 0; iDir < nDir; iDir++) {

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

      for (short iDir = 0; iDir < nDir; iDir++) {

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

      for (short iDir = 0; iDir < nDir; iDir++) {

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

              //h2sig->SetBinContent (iPtChX, iPtJY, h2tot->GetBinContent (iPtChX, iPtJY));
              //h2sig->SetBinError   (iPtChX, iPtJY, h2tot->GetBinError (iPtChX, iPtJY));
              h2sig->SetBinContent (iPtChX, iPtJY, h2tot->GetBinContent (iPtChX, iPtJY) - (iVar == 0 ? nJetSF * hbkg->GetBinContent (iPtChX) * hbkg->GetBinWidth (iPtChX) : 0));
              h2sig->SetBinError   (iPtChX, iPtJY, std::hypot (h2tot->GetBinError (iPtChX, iPtJY), (iVar == 0 ? nJetSF * hbkg->GetBinError (iPtChX) * hbkg->GetBinWidth (iPtChX) : 0)));
    
            } // end loop over iPtChY
    
          } // end loop over iPtJX

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over iVar

  } // end loop over iEvFrac




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // UNFOLD HALF-CLOSURE HISTOGRAMS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  RooUnfoldBayes* bayesUnf;
  for (short iEvFrac : {0, 1}) {

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    for (short nIter = 0; nIter < maxIters; nIter++) {

      bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_pt_ref[iEvFrac], h_jet_pt_ref[iEvFrac][0], nIter+1);
      bayesUnf->SetVerbose (-1);
      h_jet_pt_ref_unf[iEvFrac][nIter] = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_ref_unf_%s_mc_nIter%i", evFrac.Data (), nIter+1));
      SaferDelete (&bayesUnf);

      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

        bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_pt[iEvFrac][iCent], h_jet_pt[iEvFrac][iCent][0], nIter+1);
        bayesUnf->SetVerbose (-1);
        h_jet_pt_unf[iEvFrac][iCent][nIter] = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_unf_%s_%s_mc_nIter%i", evFrac.Data (), cent, nIter+1));
        SaferDelete (&bayesUnf);

      } // end loop over iCent


      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_trk_pt_ref_sig[iEvFrac][iDir], h2_jet_trk_pt_ref_sig[iEvFrac][iDir][0], nIter+1);
        bayesUnf->SetVerbose (-1);
        h2_jet_trk_pt_ref_sig_unf[iEvFrac][iDir][nIter] = (TH2D*) bayesUnf->Hreco ()->Clone (Form ("h2_jet_trk_pt_%s_ref_sig_unf_%s_mc_nIter%i", dir.Data (), evFrac.Data (), nIter+1));
        SaferDelete (&bayesUnf);

        for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
    
          const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

          bayesUnf = new RooUnfoldBayes (rooUnfResp_jet_trk_pt_sig[iEvFrac][iDir][iCent], h2_jet_trk_pt_sig[iEvFrac][iDir][iCent][0], nIter+1);
          bayesUnf->SetVerbose (-1);
          h2_jet_trk_pt_sig_unf[iEvFrac][iDir][iCent][nIter] = (TH2D*) bayesUnf->Hreco ()->Clone (Form ("h2_jet_trk_pt_%s_%s_sig_unf_%s_mc_nIter%i", dir.Data (), cent, evFrac.Data (), nIter+1));
          SaferDelete (&bayesUnf);

        } // end loop over iCent

      } // end loop over iDir

    } // end loop over nIter

  } // end loop over iEvFrac




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // SCALE TEST HISTOGRAMS BY NJET AFTER UNFOLDING
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iEvFrac : {0, 1}) {

    for (short iVar = 0; iVar < 2; iVar++) {

      TH1D* h = h_jet_pt_ref[iEvFrac][iVar]; 

      for (short iDir = 0; iDir < nDir; iDir++) {

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

        for (short iDir = 0; iDir < nDir; iDir++) {
    
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


    for (short nIter = 0; nIter < maxIters; nIter++) {

      TH1D* h = h_jet_pt_ref_unf[iEvFrac][1];

      for (short iDir = 0; iDir < nDir; iDir++) {

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

        for (short iDir = 0; iDir < nDir; iDir++) {
    
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




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // GET 1D PROJECTIONS OF JET-TAGGED HADRONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iEvFrac : {0, 1}) {

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    for (short iDir = 0; iDir < nDir; iDir++) {

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
    
      for (short nIter = 0; nIter < maxIters; nIter++) {

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




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE >30 GeV AND >60 GeV HISTOGRAMS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iEvFrac : {0, 1}) {

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
      const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
      const double maxJetPt = 300;

      for (short iDir = 0; iDir < nDir; iDir++) {
    
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
     
     
        for (short nIter = 0; nIter < maxIters; nIter++) {
    
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




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // COMPUTE IAA FOR >30 GeV AND >60 GeV HISTOGRAMS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iEvFrac : {0, 1}) {

    const TString evFrac = (iEvFrac == 0 ? "half" : "full");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (short iVar = 0; iVar < 2; iVar++) {

          const TString var = (iVar == 0 ? "Nominal" : "MCTruthJetsTruthParts");

          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

            const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

            h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][iVar] = (TH1D*) h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][iVar]->Clone (Form ("h_jetInt_trk_pt_%s_%s_iaa_%s_%s_mc_%s", dir.Data (), cent, evFrac.Data (), pTJInt.Data (), var.Data ()));
            h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][iVar]->Divide (h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][iVar]);

          } // end loop over iCent

        } // end loop over iVar
 
 
        for (short nIter = 0; nIter < maxIters; nIter++) {

          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

            const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

            h_jetInt_trk_pt_iaa_unf[iEvFrac][iPtJInt][iDir][iCent][nIter] = (TH1D*) h_jetInt_trk_pt_sig_unf[iEvFrac][iPtJInt][iDir][iCent][nIter]->Clone (Form ("h_jetInt_trk_pt_%s_%s_iaa_unf_%s_%s_mc_nIter%i", dir.Data (), cent, evFrac.Data (), pTJInt.Data (), nIter+1));
            h_jetInt_trk_pt_iaa_unf[iEvFrac][iPtJInt][iDir][iCent][nIter]->Divide (h_jetInt_trk_pt_ref_sig_unf[iEvFrac][iPtJInt][iDir][3]);

          } // end loop over iCent

        } // end loop over nIter

      } // end loop over iDir

    } // end loop over iPtJInt

  } // end loop over iEvFrac




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // COMPUTE CLOSURE HISTOGRAMS -- USE NITER = 2 RESULTS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
    //const TString funcStr = "[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)";

    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      h_jetInt_trk_pt_ref_sig_unf_halfClosure[iPtJInt][iDir] = (TH1D*) h_jetInt_trk_pt_ref_sig_unf[0][iPtJInt][iDir][1]->Clone (Form ("h_jetInt_trk_pt_%s_ref_sig_unf_%s_halfClosure", dir.Data (), pTJInt.Data ()));
      DivideNoErrors (h_jetInt_trk_pt_ref_sig_unf_halfClosure[iPtJInt][iDir], h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir][1]);

      TF1* f = DoClosureFit (Form ("f_jetInt_trk_pt_%s_ref_sig_unf_%s_halfClosure", dir.Data (), pTJInt.Data ()), h_jetInt_trk_pt_ref_sig_unf_halfClosure[iPtJInt][iDir], pTChBins[1], pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], 4, 8);
      if (!f) std::cout << "Fit failed, please investigate! iPtJint = " << iPtJInt << ", iDir = " << iDir << std::endl;
      f_jetInt_trk_pt_ref_sig_unf_halfClosure[iPtJInt][iDir] = f;


      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

        const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

        h_jetInt_trk_pt_sig_unf_halfClosure[iPtJInt][iDir][iCent] = (TH1D*) h_jetInt_trk_pt_sig_unf[0][iPtJInt][iDir][iCent][1]->Clone (Form ("h_jetInt_trk_pt_%s_%s_sig_unf_%s_halfClosure", dir.Data (), cent, pTJInt.Data ()));
        DivideNoErrors (h_jetInt_trk_pt_sig_unf_halfClosure[iPtJInt][iDir][iCent], h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent][1]);

        TF1* f = DoClosureFit (Form ("f_jetInt_trk_pt_%s_%s_sig_unf_%s_halfClosure", dir.Data (), cent, pTJInt.Data ()), h_jetInt_trk_pt_sig_unf_halfClosure[iPtJInt][iDir][iCent], pTChBins[1], pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], 4, 8);
        if (!f) std::cout << "Fit failed, please investigate! iPtJint = " << iPtJInt << ", iDir = " << iDir << ", iCent = " << iCent << std::endl;
        f_jetInt_trk_pt_sig_unf_halfClosure[iPtJInt][iDir][iCent] = f;


        h_jetInt_trk_pt_iaa_unf_halfClosure[iPtJInt][iDir][iCent] = (TH1D*) h_jetInt_trk_pt_iaa_unf[0][iPtJInt][iDir][iCent][1]->Clone (Form ("h_jetInt_trk_pt_%s_%s_iaa_unf_%s_halfClosure", dir.Data (), cent, pTJInt.Data ()));
        DivideNoErrors (h_jetInt_trk_pt_iaa_unf_halfClosure[iPtJInt][iDir][iCent], h_jetInt_trk_pt_iaa[0][iPtJInt][iDir][iCent][1]);

        f = DoClosureFit (Form ("f_jetInt_trk_pt_%s_%s_iaa_unf_%s_halfClosure", dir.Data (), cent, pTJInt.Data ()), h_jetInt_trk_pt_iaa_unf_halfClosure[iPtJInt][iDir][iCent], pTChBins[1], pTChBins[nPtChBins - (iPtJInt == 0 ? 3 : 1)], 4, 8);
        if (!f) std::cout << "Fit failed, please investigate! iPtJint = " << iPtJInt << ", iDir = " << iDir << ", iCent = " << iCent << std::endl;
        f_jetInt_trk_pt_iaa_unf_halfClosure[iPtJInt][iDir][iCent] = f;

      } // end loop over iCent
 
    } // end loop over iDir

  } // end loop over iPtJInt 



  //{
  //  TFile* outFile = new TFile (Form ("%s/aux/MCClosureHists.root", workPath.Data ()), "recreate");

  //  for (short iPtJInt : {0, 1}) {

  //    for (short iDir = 0; iDir < nDir; iDir++) {

  //      h_jetInt_trk_pt_ref_sig_unf_halfClosure[iPtJInt][iDir]->Write ();
  //      f_jetInt_trk_pt_ref_sig_unf_halfClosure[iPtJInt][iDir]->Write ();

  //      for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

  //        h_jetInt_trk_pt_sig_unf_halfClosure[iPtJInt][iDir][iCent]->Write ();
  //        f_jetInt_trk_pt_sig_unf_halfClosure[iPtJInt][iDir][iCent]->Write ();
  //    
  //        h_jetInt_trk_pt_iaa_unf_halfClosure[iPtJInt][iDir][iCent]->Write ();
  //        f_jetInt_trk_pt_iaa_unf_halfClosure[iPtJInt][iDir][iCent]->Write ();

  //      } // end loop over iCent
 
  //    } // end loop over iDir

  //  } // end loop over iPtJInt
  //}




  {
    TLine* l = new TLine ();
    TLatex* tl = new TLatex ();

    for (short iEvFrac : {0, 1}) {
  
      const TString evFrac = (iEvFrac == 0 ? "half" : "full");
  
      {
        const char* canvasName = Form ("c_jet_pt_%sClosure", evFrac.Data ());
        TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
        c->Divide (4, 2);
  
        TH1D* h = nullptr;
        TGAE* g = nullptr;
  
        double x, y;
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
            myText (0.2, 0.84, kBlack, Form ("#bf{FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
          else
            myText (0.2, 0.84, kBlack, "#bf{All centralities}", 0.06);
  
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
  
        for (short iDir = 0; iDir < 3; iDir++) {
  
          const TString dir = directions[iDir];
  
          const char* canvasName = Form ("c_jet_trk_pt_%s_%iGeVJets_%sClosure_6Iters", dir.Data (), minJetPt, evFrac.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 800, 800);
  
          TH1D* h = nullptr;
          TGAE* g = nullptr;
  
          double x, y;
          {
            gPad->SetLogx ();
  
            h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Ratio", 1, pTChBins[1], pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);
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
            g = make_graph (h_jetInt_trk_pt_ref_sig_unf[iEvFrac][iPtJInt][iDir][5]);
            ScaleGraph (g, h_jetInt_trk_pt_ref_sig[iEvFrac][iPtJInt][iDir][1]);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
          }
  
          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
            const Color_t col = colorfulColors[iCent+1];
            g = make_graph (h_jetInt_trk_pt_sig_unf[iEvFrac][iPtJInt][iDir][iCent][5]);
            ScaleGraph (g, h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][1]);
            myDraw (g, col, kOpenCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
          } // end loop over iCent
  
          //myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
          //myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
          //myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
          //myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{jet} > %i GeV, %s", minJetPt, (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side"))), 0.07);
          myText (0.65, 0.88, colorfulColors[0], "#bf{#it{pp}}", 0.032);

          for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {
            if (iCent < nFcalCentBins)
              myText (0.65, 0.88-0.04*(iCent+1), colorfulColors[iCent+1], Form ("#bf{FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
            else
              myText (0.65, 0.88-0.04*(iCent+1), colorfulColors[iCent+1], "#bf{All centralities}", 0.032);
          } // end loop over iCent
  
          c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_%s_%iGeVJets_%sClosure_6Iters.pdf", workPath.Data (), dir.Data (), minJetPt, evFrac.Data ()));
  
        } // end loop over iDir
  
      } // end loop over iPtJInt




      for (short iPtJInt : {0, 1}) {
  
        const int minJetPt = (iPtJInt == 0 ? 30 : 60);
  
        for (short iDir = 0; iDir < 3; iDir++) {
  
          const TString dir = directions[iDir];
  
          const char* canvasName = Form ("c_jet_trk_pt_%s_%iGeVJets_%sClosure", dir.Data (), minJetPt, evFrac.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
          c->Divide (4, 2);
  
          TH1D* h = nullptr;
          TGAE* g = nullptr;
  
          double x, y;
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
  
            if (iEvFrac == 0) 
              myDraw (f_jetInt_trk_pt_ref_sig_unf_halfClosure[iPtJInt][iDir], kBlack, 2, 2);
  
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
  
            if (iEvFrac == 0) 
              myDraw (f_jetInt_trk_pt_sig_unf_halfClosure[iPtJInt][iDir][iCent], kBlack, 2, 2);
  
            g = make_graph (h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][0]);
            ScaleGraph (g, h_jetInt_trk_pt_sig[iEvFrac][iPtJInt][iDir][iCent][1]);
            myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
  
            if (iCent < nFcalCentBins)
              myText (0.2, 0.84, kBlack, Form ("#bf{FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
            else
              myText (0.2, 0.84, kBlack, "#bf{All centralities}", 0.06);
  
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
  
        for (short iDir = 0; iDir < 3; iDir++) {
  
          const TString dir = directions[iDir];
  
          const char* canvasName = Form ("c_jet_trk_pt_iaa_%s_%iGeVJets_%sClosure", dir.Data (), minJetPt, evFrac.Data ());
          TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
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
 
            if (iEvFrac == 0) 
              myDraw (f_jetInt_trk_pt_iaa_unf_halfClosure[iPtJInt][iDir][iCent], kBlack, 2, 2);
  
            g = make_graph (h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][0]);
            ScaleGraph (g, h_jetInt_trk_pt_iaa[iEvFrac][iPtJInt][iDir][iCent][1]);
            myDraw (g, myCyan, kFullCircle, 1.0, 1, 2, "P", false);
            SaferDelete (&g);
  
            if (iCent < nFcalCentBins)
              myText (0.2, 0.84, kBlack, Form ("#bf{FCal %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
            else
              myText (0.2, 0.84, kBlack, "#bf{All centralities}", 0.06);
  
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




    for (short iDir = 0; iDir < 3; iDir++) {

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

      TH2D* h2 = (TH2D*) rooUnfResp_jet_trk_pt_ref_sig[1][iDir]->HresponseNoOverflow ()->Clone ("htemp");

      TAxis* xax = h2->GetXaxis ();
      TAxis* yax = h2->GetYaxis ();
      TAxis* zax = h2->GetZaxis ();

      xax->SetTitle ("Truth-level #it{p}_{T}^{jet}");
      yax->SetTitle ("Reco.-level #it{p}_{T}^{jet}");
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

      xax->SetNdivisions (-12);
      yax->SetNdivisions (-12);

      h2->DrawCopy ("colz");
      SaferDelete (&h2);

      tl->SetTextFont (43);
      tl->SetTextSize (16);
      tl->SetTextAlign (33);
      tl->SetTextAngle (30);

      tl->DrawLatex (std::floor (0.5*nPtChBins), -3, "15-20 GeV");
      tl->DrawLatex (std::floor (1.5*nPtChBins), -3, "20-30 GeV");
      tl->DrawLatex (std::floor (2.5*nPtChBins), -3, "30-45 GeV");
      tl->DrawLatex (std::floor (3.5*nPtChBins), -3, "45-60 GeV");
      tl->DrawLatex (std::floor (4.5*nPtChBins), -3, "60-90 GeV");
      tl->DrawLatex (std::floor (5.5*nPtChBins), -3, "90-120 GeV");
      tl->DrawLatex (std::floor (6.5*nPtChBins), -3, "120-160 GeV");
      tl->DrawLatex (std::floor (7.5*nPtChBins), -3, "160-200 GeV");
      tl->DrawLatex (std::floor (8.5*nPtChBins), -3, "200-240 GeV");
      tl->DrawLatex (std::floor (9.5*nPtChBins), -3, "240-300 GeV");
      tl->DrawLatex (std::floor (10.5*nPtChBins), -3, "300-350 GeV");
      tl->DrawLatex (std::floor (11.5*nPtChBins), -3, "350-400 GeV");

      tl->SetTextAlign (32);
      tl->SetTextAngle (330);

      tl->DrawLatex (-3, 1+std::floor (0.5*nPtChBins), "15-20 GeV");
      tl->DrawLatex (-3, 1+std::floor (1.5*nPtChBins), "20-30 GeV");
      tl->DrawLatex (-3, 1+std::floor (2.5*nPtChBins), "30-45 GeV");
      tl->DrawLatex (-3, 1+std::floor (3.5*nPtChBins), "45-60 GeV");
      tl->DrawLatex (-3, 1+std::floor (4.5*nPtChBins), "60-90 GeV");
      tl->DrawLatex (-3, 1+std::floor (5.5*nPtChBins), "90-120 GeV");
      tl->DrawLatex (-3, 1+std::floor (6.5*nPtChBins), "120-160 GeV");
      tl->DrawLatex (-3, 1+std::floor (7.5*nPtChBins), "160-200 GeV");
      tl->DrawLatex (-3, 1+std::floor (8.5*nPtChBins), "200-240 GeV");
      tl->DrawLatex (-3, 1+std::floor (9.5*nPtChBins), "240-300 GeV");
      tl->DrawLatex (-3, 1+std::floor (10.5*nPtChBins), "300-350 GeV");
      tl->DrawLatex (-3, 1+std::floor (11.5*nPtChBins), "350-400 GeV");

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

      l->SetLineColor (kGreen+2);
      l->DrawLine (0, nPtChBins*2, nPtChBins*10, nPtChBins*2);
      l->DrawLine (nPtChBins*2, 0, nPtChBins*2, nPtChBins*10);

      l->SetLineColor (kMagenta+2);
      l->DrawLine (nPtChBins*4, nPtChBins*4, nPtChBins*10, nPtChBins*4);
      l->DrawLine (nPtChBins*4, nPtChBins*4, nPtChBins*4, nPtChBins*10);

      l->SetLineColor (kBlack);
      l->DrawLine (nPtChBins*10, nPtChBins*2, nPtChBins*10, nPtChBins*10);
      l->DrawLine (nPtChBins*2, nPtChBins*10, nPtChBins*10, nPtChBins*10);

      c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_ResponseMatrix_ref_2D_%s.png", workPath.Data (), dir.Data ()));
      c->SaveAs (Form ("%s/Plots/Unfolding/MC_Jet_TrkPt_ResponseMatrix_ref_2D_%s.pdf", workPath.Data (), dir.Data ()));
    } // end loop over iDir




    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

      for (short iDir = 0; iDir < 3; iDir++) {

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

        TH2D* h2 = (TH2D*) rooUnfResp_jet_trk_pt_sig[1][iDir][iCent]->HresponseNoOverflow ()->Clone ("htemp");

        TAxis* xax = h2->GetXaxis ();
        TAxis* yax = h2->GetYaxis ();
        TAxis* zax = h2->GetZaxis ();

        xax->SetTitle ("Truth-level #it{p}_{T}^{jet}");
        yax->SetTitle ("Reco.-level #it{p}_{T}^{jet}");
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

        xax->SetNdivisions (-12);
        yax->SetNdivisions (-12);

        h2->DrawCopy ("colz");
        SaferDelete (&h2);

        tl->SetTextFont (43);
        tl->SetTextSize (16);
        tl->SetTextAlign (33);
        tl->SetTextAngle (30);

        tl->DrawLatex (std::floor (0.5*nPtChBins), -3, "15-20 GeV");
        tl->DrawLatex (std::floor (1.5*nPtChBins), -3, "20-30 GeV");
        tl->DrawLatex (std::floor (2.5*nPtChBins), -3, "30-45 GeV");
        tl->DrawLatex (std::floor (3.5*nPtChBins), -3, "45-60 GeV");
        tl->DrawLatex (std::floor (4.5*nPtChBins), -3, "60-90 GeV");
        tl->DrawLatex (std::floor (5.5*nPtChBins), -3, "90-120 GeV");
        tl->DrawLatex (std::floor (6.5*nPtChBins), -3, "120-160 GeV");
        tl->DrawLatex (std::floor (7.5*nPtChBins), -3, "160-200 GeV");
        tl->DrawLatex (std::floor (8.5*nPtChBins), -3, "200-240 GeV");
        tl->DrawLatex (std::floor (9.5*nPtChBins), -3, "240-300 GeV");
        tl->DrawLatex (std::floor (10.5*nPtChBins), -3, "300-350 GeV");
        tl->DrawLatex (std::floor (11.5*nPtChBins), -3, "350-400 GeV");

        tl->SetTextAlign (32);
        tl->SetTextAngle (330);

        tl->DrawLatex (-3, 1+std::floor (0.5*nPtChBins), "15-20 GeV");
        tl->DrawLatex (-3, 1+std::floor (1.5*nPtChBins), "20-30 GeV");
        tl->DrawLatex (-3, 1+std::floor (2.5*nPtChBins), "30-45 GeV");
        tl->DrawLatex (-3, 1+std::floor (3.5*nPtChBins), "45-60 GeV");
        tl->DrawLatex (-3, 1+std::floor (4.5*nPtChBins), "60-90 GeV");
        tl->DrawLatex (-3, 1+std::floor (5.5*nPtChBins), "90-120 GeV");
        tl->DrawLatex (-3, 1+std::floor (6.5*nPtChBins), "120-160 GeV");
        tl->DrawLatex (-3, 1+std::floor (7.5*nPtChBins), "160-200 GeV");
        tl->DrawLatex (-3, 1+std::floor (8.5*nPtChBins), "200-240 GeV");
        tl->DrawLatex (-3, 1+std::floor (9.5*nPtChBins), "240-300 GeV");
        tl->DrawLatex (-3, 1+std::floor (10.5*nPtChBins), "300-350 GeV");
        tl->DrawLatex (-3, 1+std::floor (11.5*nPtChBins), "350-400 GeV");

        tl->SetTextFont (43);
        tl->SetTextSize (32);
        tl->SetTextAlign (21);
        tl->SetTextAngle (0);
        
        tl->SetTextColor (kBlack);
        tl->SetTextAlign (12);

        tl->SetTextSize (26);
        tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
        tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp} + #it{p}+Pb overlay, #sqrt{s} = 5.02 TeV");
        tl->DrawLatexNDC (0.26, 0.810, iCent == nZdcCentBins ? "All centralities" : Form ("FCal %i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]));
        tl->DrawLatexNDC (0.26, 0.770, (dir == "ns" ? "Near-side" : (dir == "perp" ? "Perpendicular" : "Away-side")));

        l->SetLineColor (kGreen+2);
        l->DrawLine (nPtChBins*2, nPtChBins*2, nPtChBins*10, nPtChBins*2);
        l->DrawLine (nPtChBins*2, nPtChBins*2, nPtChBins*2, nPtChBins*10);

        l->SetLineColor (kMagenta+2);
        l->DrawLine (nPtChBins*4, nPtChBins*4, nPtChBins*10, nPtChBins*4);
        l->DrawLine (nPtChBins*4, nPtChBins*4, nPtChBins*4, nPtChBins*10);

        l->SetLineColor (kBlack);
        l->DrawLine (nPtChBins*10, nPtChBins*2, nPtChBins*10, nPtChBins*10);
        l->DrawLine (nPtChBins*2, nPtChBins*10, nPtChBins*10, nPtChBins*10);

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

      TH2D* h2 = (TH2D*) rooUnfResp_jet_pt_ref[1]->HresponseNoOverflow ()->Clone ("htemp");

      TAxis* xax = h2->GetXaxis ();
      TAxis* yax = h2->GetYaxis ();
      TAxis* zax = h2->GetZaxis ();

      xax->SetTitle ("Truth-level #it{p}_{T}^{jet}");
      yax->SetTitle ("Reco.-level #it{p}_{T}^{jet}");
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

      l->SetLineColor (kGreen+2);
      l->DrawLine (30, 30, 30, 300);
      l->DrawLine (30, 30, 300, 30);

      l->SetLineColor (kMagenta+2);
      l->DrawLine (60, 60, 60, 300);
      l->DrawLine (60, 60, 300, 60);

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

      TH2D* h2 = (TH2D*) rooUnfResp_jet_pt[1][iCent]->HresponseNoOverflow ()->Clone ("htemp");

      TAxis* xax = h2->GetXaxis ();
      TAxis* yax = h2->GetYaxis ();
      TAxis* zax = h2->GetZaxis ();

      xax->SetTitle ("Truth-level #it{p}_{T}^{jet}");
      yax->SetTitle ("Reco.-level #it{p}_{T}^{jet}");
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
      tl->DrawLatexNDC (0.26, 0.810, iCent == nZdcCentBins ? "All centralities" : Form ("FCal %i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]));

      l->SetLineColor (kGreen+2);
      l->DrawLine (30, 30, 30, 300);
      l->DrawLine (30, 30, 300, 30);

      l->SetLineColor (kMagenta+2);
      l->DrawLine (60, 60, 60, 300);
      l->DrawLine (60, 60, 300, 60);

      l->SetLineColor (kBlack);
      l->DrawLine (300, 30, 300, 300);
      l->DrawLine (30, 300, 300, 300);

      c->SaveAs (Form ("%s/Plots/Unfolding/MC_JetPt_ResponseMatrix_%s_2D.png", workPath.Data (), cent));
      c->SaveAs (Form ("%s/Plots/Unfolding/MC_JetPt_ResponseMatrix_%s_2D.pdf", workPath.Data (), cent));
    } // end loop over iCent


    /*
    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  
      for (short iDir = 0; iDir < 3; iDir++) {
  
        const TString dir = directions[iDir];
  
        const char* canvasName = Form ("c_jet_trk_pt_%s_%g-%gGeVJets_%sClosure", dir.Data (), pTJBins[iPtJ], pTJBins[iPtJ+1]);
        TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
        c->Divide (4, 2);
  
        TH1D* h = nullptr;
        TGAE* g = nullptr;
  
        double x, y;
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
            myText (0.2, 0.84, kBlack, "#bf{All centralities}", 0.06);
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
