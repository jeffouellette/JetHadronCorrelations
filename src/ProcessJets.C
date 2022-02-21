#ifndef __JetHadronCorrelatorPlotJets_C__
#define __JetHadronCorrelatorPlotJets_C__

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <math.h>

using namespace JetHadronCorrelations;

const bool doJZ123 = false;


TString GetSamp (const short iDType, const short iSamp) {
  if (iDType == 0)
    return "AllTrigs";
  switch (iSamp) {
    case 0: return "JZ0";
    case 1: return "JZ1";
    case 2: return "JZ2";
    case 3: return "JZ3";
    case 4: return "JZ123";
    case 5: return "JZ0123";
  }
  return "???";
}



void ProcessJets (const char* tag, const char* outFileTag) {

  const TString var = variations[0];
  const short nSamps = 6; // JZ0, 1, 2, 3, JZ1-3, JZ0-3


  const short nIters1DMax = 100;
  const short nIters1DMin = 1;
  const double* nIters1DVals = linspace (nIters1DMin, nIters1DMax, nIters1DMax-nIters1DMin);


  const bool doUnfOnMC = false;
  const bool useJetWgts = !doUnfOnMC;
  const bool useCentDiffUnf = true;

  const bool truncateLowPtReco = true;


  TFile* inFile = nullptr;

  TH1D***   h_evt_counts_ref                = Get2DArray <TH1D*> (2, nSamps);
  TH1D****  h_evt_counts                    = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH1D***   h_jet_pt_ref                    = Get2DArray <TH1D*> (2, nSamps);
  TH2D***   h2_jet_pt_cov_ref               = Get2DArray <TH2D*> (2, nSamps);

  TH1D****  h_jet_pt                        = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);
  TH2D****  h2_jet_pt_cov                   = Get3DArray <TH2D*> (2, nZdcCentBins+1, nSamps);

  TH1D****  h_jet_pt_ratio                  = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH1D**    h_jet_pt_datamc_ratio_ref       = Get1DArray <TH1D*> (nSamps);
  TH1D***   h_jet_pt_datamc_ratio           = Get2DArray <TH1D*> (nZdcCentBins+1, nSamps);

  TF1**     f_jet_pt_datamc_ratio_ref       = Get1DArray <TF1*> (nSamps);
  TF1***    f_jet_pt_datamc_ratio           = Get2DArray <TF1*> (nZdcCentBins+1, nSamps);


  // unfolded jet spectra
  TH1D**    h_jet_pt_ref_unf_nIters         = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+2);
  TH1D***   h_jet_pt_unf_nIters             = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+1);

  // refolded jet spectra
  TH1D**    h_jet_pt_ref_rfld_nIters        = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+2);
  TH1D***   h_jet_pt_rfld_nIters            = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+1);


  TH2D****  h2_jet_eta_phi_ref              = Get3DArray <TH2D*> (2, nPtJBins, nSamps);
  TH2D***** h2_jet_eta_phi                  = Get4DArray <TH2D*> (2, nPtJBins, nZdcCentBins+1, nSamps);


  TH2D***   h2_jet_pt_eta_jer_frac_num_ref  = Get2DArray <TH2D*> (2, nSamps);
  TH2D****  h2_jet_pt_eta_jer_frac_num      = Get3DArray <TH2D*> (2, nZdcCentBins+1, nSamps);

  TH1D***   h_jet_pt_jer_frac_num_ref       = Get2DArray <TH1D*> (2, nSamps);
  TH1D****  h_jet_pt_jer_frac_num           = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH2D***   h2_jet_pt_eta_jer_frac_den_ref  = Get2DArray <TH2D*> (2, nSamps);
  TH2D****  h2_jet_pt_eta_jer_frac_den      = Get3DArray <TH2D*> (2, nZdcCentBins+1, nSamps);

  TH1D***   h_jet_pt_jer_frac_den_ref       = Get2DArray <TH1D*> (2, nSamps);
  TH1D****  h_jet_pt_jer_frac_den           = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH2D***   h2_jet_pt_eta_jer_frac_ref      = Get2DArray <TH2D*> (2, nSamps);
  TH2D****  h2_jet_pt_eta_jer_frac          = Get3DArray <TH2D*> (2, nZdcCentBins+1, nSamps);

  TH1D***   h_jet_pt_jer_frac_ref           = Get2DArray <TH1D*> (2, nSamps);
  TH1D****  h_jet_pt_jer_frac               = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);
  TGAE***   g_jet_pt_jer_frac_ref           = Get2DArray <TGAE*> (2, nSamps);
  TGAE****  g_jet_pt_jer_frac               = Get3DArray <TGAE*> (2, nZdcCentBins+1, nSamps);


  RooUnfoldResponse*    rooUnfResp_jet_pt_ref                     = nullptr;
  RooUnfoldResponse**   rooUnfResp_jet_pt                         = Get1DArray <RooUnfoldResponse*> (nFcalCentBins+1);



  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/ProcessJets_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iSamp = 0; iSamp < nSamps; iSamp++) {

      if (iDType == 0 && iSamp > 0)
        continue;

      const TString samp = GetSamp (iDType, iSamp);

      {
        TString inFileName = Form ("%s/JetPtWeights/%s/%s17_5TeV_%s.root", rootPath.Data (), var.Data (), dType.Data (), samp.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_jet_pt_ref[iDType][iSamp]       = (TH1D*) inFile->Get ("h_jet_pt")->Clone (Form ("h_jet_pt_ref_%s_%s",      dType.Data (), samp.Data ()));
        h2_jet_pt_cov_ref[iDType][iSamp]  = (TH2D*) inFile->Get ("h2_jet_pt_cov")->Clone (Form ("h2_jet_pt_cov_ref_%s_%s", dType.Data (), samp.Data ()));

        h_evt_counts_ref[iDType][iSamp]   = (TH1D*) inFile->Get ("h_evt_counts")->Clone (Form ("h_evt_counts_ref_%s_%s",  dType.Data (), samp.Data ()));

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        //  h2_jet_eta_phi_ref[iDType][iPtJ][iSamp]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_%s17", pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_eta_phi_ref_%s_%s", dType.Data (), samp.Data ()));

        //} // end loop over iPtJ

        h2_jet_pt_eta_jer_frac_num_ref[iDType][iSamp] = (TH2D*) inFile->Get ("h2_jet_pt_eta_jer_frac_num")->Clone     (Form ("h2_jet_pt_eta_jer_frac_num_ref_%s_%s",  dType.Data (), samp.Data ()));
        h2_jet_pt_eta_jer_frac_den_ref[iDType][iSamp] = (TH2D*) inFile->Get ("h2_jet_pt_eta_jer_frac_den")->Clone     (Form ("h2_jet_pt_eta_jer_frac_den_ref_%s_%s",  dType.Data (), samp.Data ()));

        h_jet_pt_jer_frac_num_ref[iDType][iSamp]      = h2_jet_pt_eta_jer_frac_num_ref[iDType][iSamp]->ProjectionX    (Form ("h_jet_pt_jer_frac_num_ref_%s_%s",  dType.Data (), samp.Data ()));
        h_jet_pt_jer_frac_den_ref[iDType][iSamp]      = h2_jet_pt_eta_jer_frac_den_ref[iDType][iSamp]->ProjectionX    (Form ("h_jet_pt_jer_frac_den_ref_%s_%s",  dType.Data (), samp.Data ()));

        h2_jet_pt_eta_jer_frac_ref[iDType][iSamp]     = (TH2D*) h2_jet_pt_eta_jer_frac_num_ref[iDType][iSamp]->Clone  (Form ("h2_jet_pt_eta_jer_frac_ref_%s_%s",      dType.Data (), samp.Data ()));
        h_jet_pt_jer_frac_ref[iDType][iSamp]          = (TH1D*) h_jet_pt_jer_frac_num_ref[iDType][iSamp]->Clone       (Form ("h_jet_pt_jer_frac_ref_%s_%s",           dType.Data (), samp.Data ()));
        if (iDType == 1) { // for iDType == 0 (i.e. data) this quantity is not defined, and calling divide will result in divide by zero
          h2_jet_pt_eta_jer_frac_ref[iDType][iSamp]->Divide (h2_jet_pt_eta_jer_frac_den_ref[iDType][iSamp]);
          h_jet_pt_jer_frac_ref[iDType][iSamp]->Divide (h_jet_pt_jer_frac_den_ref[iDType][iSamp]);
        }
        TEfficiency* e = new TEfficiency (*(h_jet_pt_jer_frac_num_ref[iDType][iSamp]), *(h_jet_pt_jer_frac_den_ref[iDType][iSamp]));
        g_jet_pt_jer_frac_ref[iDType][iSamp] = e->CreateGraph ();
        delete e;

        inFile->Close ();

        //const double nEvts = h_evt_counts_ref[iDType][iSamp]->GetBinContent (2); // total number of accepted evts

        CalcUncertainties (h_jet_pt_ref[iDType][iSamp], h2_jet_pt_cov_ref[iDType][iSamp], h_evt_counts_ref[iDType][iSamp]);

        if (iDType == 1 && iSamp < 4)
          h_jet_pt_ref[iDType][iSamp]->Scale (h_evt_counts_ref[iDType][iSamp]->GetBinContent (2), "width");
        else
          h_jet_pt_ref[iDType][iSamp]->Scale (h_evt_counts_ref[iDType][iSamp]->GetBinContent (2));
        //else
        //  h_jet_pt_ref[iDType][iSamp]->Scale (1./h_jet_pt_ref[iDType][iSamp]->Integral (), "width");

        //if (iDType == 0) {
        //  TH1D* h = h_jet_pt_ref[iDType][iSamp];
        //  TH2D* h2 = h2_jet_pt_cov_ref[iDType][iSamp];

        //  //float normFactor = 0;
        //  //float minjpt = (strcmp (tag, "30GeVJets") == 0 ? 30 : 60);
        //  //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
        //  //  if (h->GetBinLowEdge (iX) > minjpt)
        //  //    normFactor += h->GetBinContent (iX) * h->GetBinWidth (iX);
        //  //}
        //  //h->Scale (1./normFactor);
        //  //h2->Scale (1./std::pow (normFactor, 2));
        //  //h->Scale (nEvts/nJets);
        //  //h2->Scale (std::pow (nEvts/nJets, 2));

        //  float norm = 0, avgptj = 0, avgptjerr = 0;
        //  norm = h->Integral ();
        //  std::cout << "norm = " << norm << std::endl;
        //  for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
        //    avgptj += h->GetBinCenter (iX) * h->GetBinContent (iX) / norm;
        //    for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
        //      avgptjerr += h->GetBinCenter (iX) * h2->GetBinContent (iX, iY) * h->GetBinCenter (iY) / (norm*norm);
        //    }
        //  }
        //  std::cout << "avg pt = " << avgptj << " +/-" << std::sqrt (avgptjerr) << std::endl;
        //}

        //RebinSomeBins (&(h_jet_pt_ref[iDType][iSamp]), nAltPtJBins, (double*) altPtJBins, true);

      }



      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        TString inFileName = Form ("%s/JetPtWeights/%s/%s16_5TeV_%s_%s.root", rootPath.Data (), var.Data (), dType.Data (), samp.Data (), cent.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts[iDType][iCent][iSamp]   = (TH1D*) inFile->Get ("h_evt_counts")->Clone (Form ("h_evt_counts_pPb_%s_%s_%s",    cent.Data(), dType.Data (), samp.Data ()));

        h_jet_pt[iDType][iCent][iSamp]       = (TH1D*) inFile->Get ("h_jet_pt")->Clone (Form ("h_jet_pt_pPb_%s_%s_%s",            cent.Data (), dType.Data (), samp.Data ()));
        h2_jet_pt_cov[iDType][iCent][iSamp]  = (TH2D*) inFile->Get ("h2_jet_pt_cov")->Clone (Form ("h2_jet_pt_cov_pPb_%s_%s_%s",  cent.Data (), dType.Data (), samp.Data ()));

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        //  h2_jet_eta_phi[iDType][iPtJ][iCent][iSamp]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_%s16", pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_eta_phi_pPb_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), samp.Data ()));

        //} // end loop over iPtJ

        h2_jet_pt_eta_jer_frac_num[iDType][iCent][iSamp]  = (TH2D*) inFile->Get ("h2_jet_pt_eta_jer_frac_num")->Clone       (Form ("h2_jet_pt_eta_jer_frac_num_%s_%s_%s", cent.Data (), dType.Data (), samp.Data ()));
        h2_jet_pt_eta_jer_frac_den[iDType][iCent][iSamp]  = (TH2D*) inFile->Get ("h2_jet_pt_eta_jer_frac_den")->Clone       (Form ("h2_jet_pt_eta_jer_frac_den_%s_%s_%s", cent.Data (), dType.Data (), samp.Data ()));

        h_jet_pt_jer_frac_num[iDType][iCent][iSamp]       = h2_jet_pt_eta_jer_frac_num[iDType][iCent][iSamp]->ProjectionX   (Form ("h_jet_pt_jer_frac_num_%s_%s_%s",      cent.Data (), dType.Data (), samp.Data ()));
        h_jet_pt_jer_frac_den[iDType][iCent][iSamp]       = h2_jet_pt_eta_jer_frac_den[iDType][iCent][iSamp]->ProjectionX   (Form ("h_jet_pt_jer_frac_den_%s_%s_%s",      cent.Data (), dType.Data (), samp.Data ()));

        h2_jet_pt_eta_jer_frac[iDType][iCent][iSamp]      = (TH2D*) h2_jet_pt_eta_jer_frac_num[iDType][iCent][iSamp]->Clone (Form ("h2_jet_pt_eta_jer_frac_%s_%s_%s",     cent.Data (), dType.Data (), samp.Data ()));
        h_jet_pt_jer_frac[iDType][iCent][iSamp]           = (TH1D*) h_jet_pt_jer_frac_num[iDType][iCent][iSamp]->Clone      (Form ("h_jet_pt_jer_frac_%s_%s_%s",          cent.Data (), dType.Data (), samp.Data ()));
        if (iDType == 1) { // for iDType == 0 (i.e. data) this quantity is not defined, and calling divide will result in divide by zero
          h2_jet_pt_eta_jer_frac[iDType][iCent][iSamp]->Divide (h2_jet_pt_eta_jer_frac_den[iDType][iCent][iSamp]);
          h_jet_pt_jer_frac[iDType][iCent][iSamp]->Divide (h_jet_pt_jer_frac_den[iDType][iCent][iSamp]);
        }
        TEfficiency* e = new TEfficiency (*(h_jet_pt_jer_frac_num[iDType][iCent][iSamp]), *(h_jet_pt_jer_frac_den[iDType][iCent][iSamp]));
        g_jet_pt_jer_frac[iDType][iCent][iSamp] = e->CreateGraph ();
        delete e;

        inFile->Close ();
      
        //const double nEvts = h_evt_counts[iDType][iCent][iSamp]->GetBinContent (2); // total number of accepted evts

        CalcUncertainties (h_jet_pt[iDType][iCent][iSamp], h2_jet_pt_cov[iDType][iCent][iSamp], h_evt_counts[iDType][iCent][iSamp]);

        if (iDType == 1 && iSamp < 4)
          h_jet_pt[iDType][iCent][iSamp]->Scale (h_evt_counts[iDType][iCent][iSamp]->GetBinContent (2), "width");
        //else
        //  h_jet_pt[iDType][iCent][iSamp]->Scale (1./h_jet_pt[iDType][iCent][iSamp]->Integral (), "width");
        else
          h_jet_pt[iDType][iCent][iSamp]->Scale (h_evt_counts[iDType][iCent][iSamp]->GetBinContent (2));

        //if (iDType == 0) {
        //  TH1D* h = h_jet_pt[iDType][iCent][iSamp];
        //  TH2D* h2 = h2_jet_pt_cov[iDType][iCent][iSamp];

        //  //float normFactor = 0;
        //  //float minjpt = (strcmp (tag, "30GeVJets") == 0 ? 30 : 60);
        //  //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
        //  //  if (h->GetBinLowEdge (iX) > minjpt)
        //  //    normFactor += h->GetBinContent (iX) * h->GetBinWidth (iX);
        //  //}
        //  //h->Scale (1./normFactor);
        //  //h2->Scale (1./std::pow (normFactor, 2));
        //  //h->Scale (nEvts/nJets);
        //  //h2->Scale (std::pow (nEvts/nJets, 2));

        //  float norm = 0, avgptj = 0, avgptjerr = 0;
        //  norm = h->Integral ();
        //  std::cout << "norm = " << norm << std::endl;
        //  for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
        //    avgptj += h->GetBinCenter (iX) * h->GetBinContent (iX) / norm;
        //    for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
        //      avgptjerr += h->GetBinCenter (iX) * h2->GetBinContent (iX, iY) * h->GetBinCenter (iY) / (norm*norm);
        //    }
        //  }
        //  std::cout << "avg pt = " << avgptj << " +/-" << std::sqrt (avgptjerr) << std::endl;
        //}

        //RebinSomeBins (&(h_jet_pt[iDType][iCent][iSamp]), nAltPtJBins, (double*) altPtJBins, true);

      } // end loop over iCent



      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        h_jet_pt_ratio[iDType][iCent][iSamp] = (TH1D*) h_jet_pt[iDType][iCent][iSamp]->Clone (Form ("h_jet_pt_ratio_%s_%s_%s", cent.Data (), dType.Data (), samp.Data ()));
        h_jet_pt_ratio[iDType][iCent][iSamp]->Divide (h_jet_pt_ref[iDType][iSamp]);

      } // end loop over iCent

    } // end loop over iSamp

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // READ IN RESPONSE MATRICES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    TFile* inFile = new TFile (Form ("%s/MakeResponseMatrix/Nominal/allSamples_finePtJBins%s.root", rootPath.Data (), truncateLowPtReco ? "_truncated" : ""), "read");

    if (useJetWgts) rooUnfResp_jet_pt_ref = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_altwgts");
    else            rooUnfResp_jet_pt_ref = (RooUnfoldResponse*) inFile->Get ("rooUnfResp_jet_pt_ref_mc_fullClosure");

    for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

      const char* cent = (iCent == nFcalCentBins ? "pPb_allCent" : Form ("pPb_iCent%i", iCent));

      if (useJetWgts) rooUnfResp_jet_pt[iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_altwgts", cent));
      else            rooUnfResp_jet_pt[iCent] = (RooUnfoldResponse*) inFile->Get (Form ("rooUnfResp_jet_pt_%s_mc_fullClosure", cent));

    } // end loop over iCent
  }




  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// CREATE DATA/MC HISTOGRAMS
  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iSamp = 0; iSamp < nSamps; iSamp++) {

    const TString samp = GetSamp (1, iSamp);

    const TString funcStr = "[0]+exp([1]+[2]*log(x))";

    {
      h_jet_pt_datamc_ratio_ref[iSamp] = (TH1D*) h_jet_pt_ref[0][0]->Clone (Form ("h_jet_pt_datamc_ratio_ref_%s", samp.Data ()));
      h_jet_pt_datamc_ratio_ref[iSamp]->Divide (h_jet_pt_ref[1][iSamp]);
      h_jet_pt_datamc_ratio_ref[iSamp]->Scale (h_jet_pt_ref[1][iSamp]->Integral () / h_jet_pt_ref[0][0]->Integral ());

      TF1* f = new TF1 (Form ("f_jet_pt_datamc_ratio_ref_%s", samp.Data ()), funcStr.Data (), 15/*pTJBins[0]*/, pTJBins[nPtJBins]);
      h_jet_pt_datamc_ratio_ref[iSamp]->Fit (f, "RN0QI");
      f_jet_pt_datamc_ratio_ref[iSamp] = f;
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      h_jet_pt_datamc_ratio[iCent][iSamp] = (TH1D*) h_jet_pt[0][iCent][0]->Clone (Form ("h_jet_pt_datamc_ratio_%s_%s", cent.Data (), samp.Data ()));
      h_jet_pt_datamc_ratio[iCent][iSamp]->Divide (h_jet_pt[1][iCent][iSamp]);
      h_jet_pt_datamc_ratio[iCent][iSamp]->Scale (h_jet_pt[1][iCent][iSamp]->Integral () / h_jet_pt[0][iCent][0]->Integral ());

      TF1* f = new TF1 (Form ("f_jet_pt_datamc_ratio_%s_%s", cent.Data (), samp.Data ()), funcStr.Data (), 15/*pTJBins[0]*/, pTJBins[nPtJBins]);
      h_jet_pt_datamc_ratio[iCent][iSamp]->Fit (f, "RN0QI");
      f_jet_pt_datamc_ratio[iCent][iSamp] = f;

    } // end loop over iCent

  }




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // DO UNFOLD STUDIES VS NUMBER OF ITERATIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    const short iDType = (doUnfOnMC ? 1 : 0);
    const short iSamp = (doUnfOnMC ? nSamps-1 : 0); 

    {
      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {

        const short nIters = (short) nIters1DVals[iIter];
        std::cout << "|--> Doing 1D unfold in pp at " << nIters << " iterations... " << std::endl;

        TH1D* h_raw = (TH1D*) h_jet_pt_ref[iDType][iSamp]->Clone ("h_raw");
        if (truncateLowPtReco) {
          for (int iX = 1; iX <= h_raw->GetNbinsX (); iX++) {
            if (h_raw->GetBinCenter (iX) < 20) {
              h_raw->SetBinContent (iX, 0);
              h_raw->SetBinError (iX, 0);
            }
          }
        }

        RooUnfoldResponse* resp = rooUnfResp_jet_pt_ref;
        resp->UseOverflow (0);
        RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (resp, h_raw, nIters);
        bayesUnf->SetVerbose (-1);
        //bayesUnf->SetMeasuredCov (TMatrixD (nPtJBins, nPtJBins, h2_jet_pt_ref_cov->GetArray ()));
        TH1D* h_unf = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_ref_unf_data_nIters%i", nIters));
        TH1D* h_rfld = (TH1D*) resp->ApplyToTruth (h_unf)->Clone (Form ("h_jet_pt_ref_rfld_data_nIters%i", nIters));

        SaferDelete (&h_raw);
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
        std::cout << "|--> Doing 1D unfold in p+Pb (" << cent.Data () << ") at " << nIters << " iterations... " << std::endl;

        TH1D* h_raw = (TH1D*) h_jet_pt[iDType][iCent][iSamp]->Clone ("h_raw");
        if (truncateLowPtReco) {
          for (int iX = 1; iX <= h_raw->GetNbinsX (); iX++) {
            if (h_raw->GetBinCenter (iX) < 20) {
              h_raw->SetBinContent (iX, 0);
              h_raw->SetBinError (iX, 0);
            }
          }
        }

        RooUnfoldResponse* resp = rooUnfResp_jet_pt[iUnfCent];
        resp->UseOverflow (0);
        RooUnfoldBayes* bayesUnf = new RooUnfoldBayes (resp, h_raw, nIters);
        bayesUnf->SetVerbose (-1);
        //bayesUnf->SetMeasuredCov (TMatrixD (nPtJBins, nPtJBins, h2_jet_pt_cov[iCent]->GetArray ()));
        TH1D* h_unf = (TH1D*) bayesUnf->Hreco ()->Clone (Form ("h_jet_pt_unf_data_%s_nIters%i", cent.Data (), nIters));
        TH1D* h_rfld = (TH1D*) resp->ApplyToTruth (h_unf)->Clone (Form ("h_jet_pt_rfld_data_%s_nIters%i", cent.Data (), nIters));
        SaferDelete (&h_raw);
        SaferDelete (&bayesUnf);

        h_jet_pt_unf_nIters[iCent][iIter] = h_unf;
        h_jet_pt_rfld_nIters[iCent][iIter] = h_rfld;

      } // end loop over iIter

    } // end loop over iCent
  }




  {
    outFile->cd ();

    for (short iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iSamp = 0; iSamp < nSamps; iSamp++) {

        const TString samp = GetSamp (iDType, iSamp);

        if (iDType == 0 && iSamp > 0)
          continue;

        h_evt_counts_ref[iDType][iSamp]->Write ();

        h_jet_pt_ref[iDType][iSamp]->Write ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          h_evt_counts[iDType][iCent][iSamp]->Write ();

          h_jet_pt[iDType][iCent][iSamp]->Write ();
          h_jet_pt_ratio[iDType][iCent][iSamp]->Write ();

        } // end loop over iCent

        h2_jet_pt_eta_jer_frac_num_ref[iDType][iSamp]->Write ();
        h2_jet_pt_eta_jer_frac_den_ref[iDType][iSamp]->Write ();
        h2_jet_pt_eta_jer_frac_ref[iDType][iSamp]->Write ();

        h_jet_pt_jer_frac_num_ref[iDType][iSamp]->Write ();
        h_jet_pt_jer_frac_den_ref[iDType][iSamp]->Write ();
        h_jet_pt_jer_frac_ref[iDType][iSamp]->Write ();
        g_jet_pt_jer_frac_ref[iDType][iSamp]->Write (Form ("g_jet_pt_jer_frac_ref_%s_%s", dType.Data (), samp.Data ()));

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          h2_jet_pt_eta_jer_frac_num[iDType][iCent][iSamp]->Write ();
          h2_jet_pt_eta_jer_frac_den[iDType][iCent][iSamp]->Write ();
          h2_jet_pt_eta_jer_frac[iDType][iCent][iSamp]->Write ();

          h_jet_pt_jer_frac_num[iDType][iCent][iSamp]->Write ();
          h_jet_pt_jer_frac_den[iDType][iCent][iSamp]->Write ();
          h_jet_pt_jer_frac[iDType][iCent][iSamp]->Write ();
          g_jet_pt_jer_frac[iDType][iCent][iSamp]->Write (Form ("g_jet_pt_jer_frac_%s_%s_%s", cent.Data (), dType.Data (), samp.Data ()));

        } // end loop over iCent

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  h2_jet_eta_phi_ref[iDType][iPtJ][iSamp]->Write ();

        //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        //    h2_jet_eta_phi[iDType][iPtJ][iCent][iSamp]->Write ();

        //  } // end loop over iCent

        //} // end loop over iPtJ

      } // end loop over iSamp

    } // end loop over iDType


    for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {

      h_jet_pt_ref_unf_nIters[iIter]->Write ();
      h_jet_pt_ref_rfld_nIters[iIter]->Write ();

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  
        h_jet_pt_unf_nIters[iCent][iIter]->Write ();
        h_jet_pt_rfld_nIters[iCent][iIter]->Write ();
  
      } // end loop over iCent

    }


    for (short iSamp = 0; iSamp < nSamps; iSamp++) {

      h_jet_pt_datamc_ratio_ref[iSamp]->Write ();
      f_jet_pt_datamc_ratio_ref[iSamp]->Write ();

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        h_jet_pt_datamc_ratio[iCent][iSamp]->Write ();
        f_jet_pt_datamc_ratio[iCent][iSamp]->Write ();

      } // end loop over iCent

    } // end loop over iSamp



    TFile* wgtsFile = new TFile (Form ("%s/aux/JetPtWeights.root", workPath.Data ()), "recreate");
    wgtsFile->cd ();

    for (short iSamp = 0; iSamp < nSamps; iSamp++) {

      h_jet_pt_datamc_ratio_ref[iSamp]->Write ();
      f_jet_pt_datamc_ratio_ref[iSamp]->Write ();

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        h_jet_pt_datamc_ratio[iCent][iSamp]->Write ();
        f_jet_pt_datamc_ratio[iCent][iSamp]->Write ();

      } // end loop over iCent

    } // end loop over iSamp


    outFile->Close ();
    wgtsFile->Close ();

  }
}


#endif
