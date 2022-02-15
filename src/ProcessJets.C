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

  TFile* inFile = nullptr;

  TH1D***   h_evt_counts_ref     = Get2DArray <TH1D*> (2, nSamps);
  TH1D****  h_evt_counts         = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH1D***   h_jet_pt_ref          = Get2DArray <TH1D*> (3, nSamps);
  TH2D***   h2_jet_pt_cov_ref     = Get2DArray <TH2D*> (2, nSamps);

  TH1D****  h_jet_pt             = Get3DArray <TH1D*> (3, nZdcCentBins+1, nSamps);
  TH2D****  h2_jet_pt_cov        = Get3DArray <TH2D*> (2, nZdcCentBins+1, nSamps);

  TH1D****  h_jet_pt_ratio       = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH1D***   h_jet_pt_datamc_ratio_ref  = Get2DArray <TH1D*> (2, nSamps);
  TH1D****  h_jet_pt_datamc_ratio      = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TF1***    f_jet_pt_datamc_ratio_ref  = Get2DArray <TF1*> (2, nSamps);
  TF1****   f_jet_pt_datamc_ratio      = Get3DArray <TF1*> (2, nZdcCentBins+1, nSamps);

  //TH2D****  h2_jet_eta_phi_ref   = Get3DArray <TH2D*> (2, nPtJBins, nSamps);
  //TH2D***** h2_jet_eta_phi       = Get4DArray <TH2D*> (2, nPtJBins, nZdcCentBins+1, nSamps);

  //const int nAltPtJBins = (strcmp (tag, "30GeVJets") == 0 ? 13 : 16);
  //double* altPtJBins = new double[nAltPtJBins+1];
  //if (strcmp (tag, "30GeVJets") == 0) {
  //  altPtJBins[0] = pTJBins[0];
  //  altPtJBins[1] = pTJBins[2];
  //  altPtJBins[2] = pTJBins[4];
  //  altPtJBins[3] = pTJBins[6];
  //  altPtJBins[4] = pTJBins[8];
  //  altPtJBins[5] = pTJBins[10];
  //  altPtJBins[6] = pTJBins[12];
  //  altPtJBins[7] = pTJBins[14];
  //  altPtJBins[8] = pTJBins[18];
  //  altPtJBins[9] = pTJBins[24];
  //  altPtJBins[10] = pTJBins[30];
  //  altPtJBins[11] = pTJBins[40];
  //  altPtJBins[12] = pTJBins[50];
  //  altPtJBins[13] = pTJBins[60];
  //}
  //else {
  //  altPtJBins[0] = pTJBins[0];
  //  altPtJBins[1] = pTJBins[3];
  //  altPtJBins[2] = pTJBins[6];
  //  altPtJBins[3] = pTJBins[9];
  //  altPtJBins[4] = pTJBins[12];
  //  altPtJBins[5] = pTJBins[15];
  //  altPtJBins[6] = pTJBins[18];
  //  altPtJBins[7] = pTJBins[21];
  //  altPtJBins[8] = pTJBins[24];
  //  altPtJBins[9] = pTJBins[27];
  //  altPtJBins[10] = pTJBins[30];
  //  altPtJBins[11] = pTJBins[33];
  //  altPtJBins[12] = pTJBins[36];
  //  altPtJBins[13] = pTJBins[42];
  //  altPtJBins[14] = pTJBins[48];
  //  altPtJBins[15] = pTJBins[54];
  //  altPtJBins[16] = pTJBins[60];
  //}
  


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

        inFile->Close ();

        //const double nEvts = h_evt_counts_ref[iDType][iSamp]->GetBinContent (2); // total number of accepted evts

        CalcUncertainties (h_jet_pt_ref[iDType][iSamp], h2_jet_pt_cov_ref[iDType][iSamp], h_evt_counts_ref[iDType][iSamp]);

        h_jet_pt_ref[iDType][iSamp]->Scale (1./h_jet_pt_ref[iDType][iSamp]->Integral (), "width");

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

        h_evt_counts[iDType][iCent][iSamp]   = (TH1D*) inFile->Get ("h_evt_counts")->Clone (Form ("h_evt_counts_pPb_%s_%s_%s", cent.Data(), dType.Data (), samp.Data ()));

        h_jet_pt[iDType][iCent][iSamp]       = (TH1D*) inFile->Get ("h_jet_pt")->Clone (Form ("h_jet_pt_pPb_%s_%s_%s",       cent.Data (), dType.Data (), samp.Data ()));
        h2_jet_pt_cov[iDType][iCent][iSamp]  = (TH2D*) inFile->Get ("h2_jet_pt_cov")->Clone (Form ("h2_jet_pt_cov_pPb_%s_%s_%s",  cent.Data (), dType.Data (), samp.Data ()));

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        //  h2_jet_eta_phi[iDType][iPtJ][iCent][iSamp]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_%s16", pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_eta_phi_pPb_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), samp.Data ()));

        //} // end loop over iPtJ

        inFile->Close ();
      
        //const double nEvts = h_evt_counts[iDType][iCent][iSamp]->GetBinContent (2); // total number of accepted evts

        CalcUncertainties (h_jet_pt[iDType][iCent][iSamp], h2_jet_pt_cov[iDType][iCent][iSamp], h_evt_counts[iDType][iCent][iSamp]);

        h_jet_pt[iDType][iCent][iSamp]->Scale (1./h_jet_pt[iDType][iCent][iSamp]->Integral (), "width");

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




  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// SMOOTH TRANSITION AT 60 GEV (IF USING JZ1 and 2) AND CREATE DATA/MC HISTOGRAMS
  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iSamp = 0; iSamp < nSamps; iSamp++) {

    const TString samp = GetSamp (1, iSamp);

    const TString funcStr = "[0]+exp([1]+[2]*log(x))";

    {

      h_jet_pt_ref[2][iSamp] = (TH1D*) h_jet_pt_ref[1][iSamp]->Clone (Form ("h_jet_pt_ref_mcScaled_%s", samp.Data ()));

      if (iSamp >= 4) {
        TF1* f_mb   = new TF1 ("f_mb",  "exp([0]+[1]*log(x))", 30, 50);
        TF1* f_jet  = new TF1 ("f_jet", "exp([0]+[1]*log(x))", 70, 200);
        short bin = h_jet_pt_ref[1][iSamp]->FindBin (30);
        f_mb->SetParameter (0, std::log (h_jet_pt_ref[1][iSamp]->GetBinContent (bin)));
        f_mb->SetParameter (1, -5.5);
        bin = h_jet_pt_ref[1][iSamp]->FindBin (60);
        f_jet->SetParameter (0, std::log (h_jet_pt_ref[1][iSamp]->GetBinContent (bin)));
        f_jet->SetParameter (1, -5.5);

        h_jet_pt_ref[1][iSamp]->Fit (f_mb, "RN0QI");
        h_jet_pt_ref[1][iSamp]->Fit (f_jet, "RN0QI");

        const double sf = f_mb->Eval (60) / f_jet->Eval (60);
        std::cout << "sf = " << f_mb->Eval (60) << " / " << f_jet->Eval (60) << " = " << sf << std::endl;

        for (int iX = 1; iX <= h_jet_pt_ref[2][iSamp]->GetNbinsX (); iX++) {
          if (h_jet_pt_ref[2][iSamp]->GetBinCenter (iX) < 60) continue;

          h_jet_pt_ref[2][iSamp]->SetBinContent (iX, h_jet_pt_ref[2][iSamp]->GetBinContent (iX) * sf);
        }

        SaferDelete (&f_mb);
        SaferDelete (&f_jet);
      }


      h_jet_pt_datamc_ratio_ref[0][iSamp] = (TH1D*) h_jet_pt_ref[0][0]->Clone (Form ("h_jet_pt_datamc_ratio_ref_%s", samp.Data ()));
      h_jet_pt_datamc_ratio_ref[0][iSamp]->Divide (h_jet_pt_ref[1][iSamp]);

      TF1* f = new TF1 (Form ("f_jet_pt_datamc_ratio_ref_%s", samp.Data ()), funcStr.Data (), pTJBins[0], pTJBins[nPtJBins]);
      h_jet_pt_datamc_ratio_ref[0][iSamp]->Fit (f, "RN0QI");
      f_jet_pt_datamc_ratio_ref[0][iSamp] = f;


      h_jet_pt_datamc_ratio_ref[1][iSamp] = (TH1D*) h_jet_pt_ref[0][0]->Clone (Form ("h_jet_pt_datamcScaled_ratio_ref_%s", samp.Data ()));
      h_jet_pt_datamc_ratio_ref[1][iSamp]->Divide (h_jet_pt_ref[2][iSamp]);

      f = new TF1 (Form ("f_jet_pt_datamcScaled_ratio_ref_%s", samp.Data ()), funcStr.Data (), pTJBins[0], pTJBins[nPtJBins]);
      h_jet_pt_datamc_ratio_ref[1][iSamp]->Fit (f, "RN0QI");
      f_jet_pt_datamc_ratio_ref[1][iSamp] = f;
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      h_jet_pt[2][iCent][iSamp] = (TH1D*) h_jet_pt[1][iCent][iSamp]->Clone (Form ("h_jet_pt_pPb_%s_mcScaled_%s", cent.Data (), samp.Data ()));

      if (iSamp >= 4) {
        TF1* f_mb   = new TF1 ("f_mb",  "exp([0]+[1]*log(x))", 30, 45);
        TF1* f_jet  = new TF1 ("f_jet", "exp([0]+[1]*log(x))", 60, 200);
        short bin = h_jet_pt_ref[1][iSamp]->FindBin (30);
        f_mb->SetParameter (0, std::log (h_jet_pt[1][iCent][iSamp]->GetBinContent (bin)));
        f_mb->SetParameter (1, -5.5);
        bin = h_jet_pt[1][iCent][iSamp]->FindBin (60);
        f_jet->SetParameter (0, std::log (h_jet_pt[1][iCent][iSamp]->GetBinContent (bin)));
        f_jet->SetParameter (1, -5.5);

        h_jet_pt[1][iCent][iSamp]->Fit (f_mb, "RN0QI");
        h_jet_pt[1][iCent][iSamp]->Fit (f_jet, "RN0QI");

        const double sf = f_mb->Eval (60) / f_jet->Eval (60);
        std::cout << "sf = " << f_mb->Eval (60) << " / " << f_jet->Eval (60) << " = " << sf << std::endl;
  
        for (int iX = 1; iX <= h_jet_pt[2][iCent][iSamp]->GetNbinsX (); iX++) {
          if (h_jet_pt[2][iCent][iSamp]->GetBinCenter (iX) < 60) continue;
  
          h_jet_pt[2][iCent][iSamp]->SetBinContent (iX, h_jet_pt[2][iCent][iSamp]->GetBinContent (iX) * sf);
        }

        SaferDelete (&f_mb);
        SaferDelete (&f_jet);
      }


      h_jet_pt_datamc_ratio[0][iCent][iSamp] = (TH1D*) h_jet_pt[0][iCent][0]->Clone (Form ("h_jet_pt_datamc_ratio_%s_%s", cent.Data (), samp.Data ()));
      h_jet_pt_datamc_ratio[0][iCent][iSamp]->Divide (h_jet_pt[1][iCent][iSamp]);

      TF1* f = new TF1 (Form ("f_jet_pt_datamc_ratio_%s_%s", cent.Data (), samp.Data ()), funcStr.Data (), pTJBins[0], pTJBins[nPtJBins]);
      h_jet_pt_datamc_ratio[0][iCent][iSamp]->Fit (f, "RN0QI");
      f_jet_pt_datamc_ratio[0][iCent][iSamp] = f;


      h_jet_pt_datamc_ratio[1][iCent][iSamp] = (TH1D*) h_jet_pt[0][iCent][0]->Clone (Form ("h_jet_pt_datamcScaled_ratio_%s_%s", cent.Data (), samp.Data ()));
      h_jet_pt_datamc_ratio[1][iCent][iSamp]->Divide (h_jet_pt[2][iCent][iSamp]);

      f = new TF1 (Form ("f_jet_pt_datamcScaled_ratio_%s_%s", cent.Data (), samp.Data ()), funcStr.Data (), pTJBins[0], pTJBins[nPtJBins]);
      h_jet_pt_datamc_ratio[1][iCent][iSamp]->Fit (f, "RN0QI");
      f_jet_pt_datamc_ratio[1][iCent][iSamp] = f;

    } // end loop over iCent

  }




  {
    outFile->cd ();

    for (short iDType = 0; iDType < 2; iDType++) {

      for (short iSamp = 0; iSamp < nSamps; iSamp++) {

        if (iDType == 0 && iSamp > 0)
          continue;

        h_evt_counts_ref[iDType][iSamp]->Write ();

        h_jet_pt_ref[iDType][iSamp]->Write ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          h_evt_counts[iDType][iCent][iSamp]->Write ();

          h_jet_pt[iDType][iCent][iSamp]->Write ();
          h_jet_pt_ratio[iDType][iCent][iSamp]->Write ();

        } // end loop over iCent

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  h2_jet_eta_phi_ref[iDType][iPtJ][iSamp]->Write ();

        //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        //    h2_jet_eta_phi[iDType][iPtJ][iCent][iSamp]->Write ();

        //  } // end loop over iCent

        //} // end loop over iPtJ

      } // end loop over iSamp

    } // end loop over iDType


    for (short iSamp = 0; iSamp < nSamps; iSamp++) {

      h_jet_pt_ref[2][iSamp]->Write ();

      h_jet_pt_datamc_ratio_ref[0][iSamp]->Write ();
      f_jet_pt_datamc_ratio_ref[0][iSamp]->Write ();

      h_jet_pt_datamc_ratio_ref[1][iSamp]->Write ();
      f_jet_pt_datamc_ratio_ref[1][iSamp]->Write ();

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        h_jet_pt[2][iCent][iSamp]->Write ();

        h_jet_pt_datamc_ratio[0][iCent][iSamp]->Write ();
        f_jet_pt_datamc_ratio[0][iCent][iSamp]->Write ();

        h_jet_pt_datamc_ratio[1][iCent][iSamp]->Write ();
        f_jet_pt_datamc_ratio[1][iCent][iSamp]->Write ();

      } // end loop over iCent

    } // end loop over iSamp



    TFile* wgtsFile = new TFile (Form ("%s/aux/JetPtWeights.root", workPath.Data ()), "recreate");
    wgtsFile->cd ();

    for (short iSamp = 0; iSamp < nSamps; iSamp++) {

      h_jet_pt_datamc_ratio_ref[0][iSamp]->Write ();
      f_jet_pt_datamc_ratio_ref[0][iSamp]->Write ();

      h_jet_pt_datamc_ratio_ref[1][iSamp]->Write ();
      f_jet_pt_datamc_ratio_ref[1][iSamp]->Write ();

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        h_jet_pt_datamc_ratio[0][iCent][iSamp]->Write ();
        f_jet_pt_datamc_ratio[0][iCent][iSamp]->Write ();

        h_jet_pt_datamc_ratio[1][iCent][iSamp]->Write ();
        f_jet_pt_datamc_ratio[1][iCent][iSamp]->Write ();

      } // end loop over iCent

    } // end loop over iSamp


    outFile->Close ();
    wgtsFile->Close ();

  }
}


#endif
