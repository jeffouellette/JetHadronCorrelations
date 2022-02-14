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

void ProcessJets (const char* tag, const char* outFileTag) {

  TFile* inFile = nullptr;

  TH1D***   h_evt_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  //TH1D****  h_jet_counts_ref     = Get3DArray <TH1D*> (2, nPtJBins, nVar);
  TH1D****  h_evt_counts         = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  //TH1D***** h_jet_counts         = Get4DArray <TH1D*> (2, nPtJBins, nZdcCentBins+1, nVar);

  TH1D***   h_jet_pt_ref          = Get2DArray <TH1D*> (2, nVar);
  TH2D***   h2_jet_pt_cov_ref     = Get2DArray <TH2D*> (2, nVar);

  TH1D****  h_jet_pt             = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH2D****  h2_jet_pt_cov        = Get3DArray <TH2D*> (2, nZdcCentBins+1, nVar);

  TH1D****  h_jet_pt_ratio       = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);

  TH1D**    h_jet_pt_datamc_ratio_ref  = Get1DArray <TH1D*> (nVar);
  TH1D***   h_jet_pt_datamc_ratio      = Get2DArray <TH1D*> (nZdcCentBins+1, nVar);

  TF1**     f_jet_pt_datamc_ratio_ref  = Get1DArray <TF1*> (nVar);
  TF1***    f_jet_pt_datamc_ratio      = Get2DArray <TF1*> (nZdcCentBins+1, nVar);

  //TH2D****  h2_jet_eta_phi_ref   = Get3DArray <TH2D*> (2, nPtJBins, nVar);
  //TH2D***** h2_jet_eta_phi       = Get4DArray <TH2D*> (2, nPtJBins, nZdcCentBins+1, nVar);

  TGAE**    g_jet_pt_ref_syst     = Get1DArray <TGAE*> (nVar);
  TGAE***   g_jet_pt_syst         = Get2DArray <TGAE*> (nZdcCentBins+1, nVar);

  TGAE***   g_jet_pt_ratio_syst   = Get2DArray <TGAE*> (nZdcCentBins+1, nVar);

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


  for (int iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    //for (int iVar = 0; iVar < nVar; iVar++) {
    for (int iVar : {0}) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      {
        TString inFileName = Form ("%s/JetPtWeights/%s/%s17_5TeV%s.root", rootPath.Data (), var.Data (), dType.Data (), iDType == 1 && doJZ123 ? "_JZ123" : "");
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_jet_pt_ref[iDType][iVar]        = (TH1D*) inFile->Get ("h_jet_pt")->Clone (Form ("h_jet_pt_ref_%s_%s",      dType.Data (), var.Data ()));
        h2_jet_pt_cov_ref[iDType][iVar]   = (TH2D*) inFile->Get ("h2_jet_pt_cov")->Clone (Form ("h2_jet_pt_cov_ref_%s_%s", dType.Data (), var.Data ()));

        h_evt_counts_ref[iDType][iVar]    = (TH1D*) inFile->Get ("h_evt_counts")->Clone (Form ("h_evt_counts_ref_%s_%s",  dType.Data (), var.Data ()));

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        //  h_jet_counts_ref[iDType][iPtJ][iVar]    = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s17",  pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_ref_%s_%s_%s",  dType.Data (), pTJ.Data (), var.Data ()));

        //  h2_jet_eta_phi_ref[iDType][iPtJ][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_%s17", pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_eta_phi_ref_%s_%s", dType.Data (), var.Data ()));

        //} // end loop over iPtJ

        inFile->Close ();

        //const double nEvts = h_evt_counts_ref[iDType][iVar]->GetBinContent (2); // total number of accepted evts
        //const double nJets = h_jet_counts_ref[iDType][iVar]->GetBinContent (2); // total number of triggered jets

        //CalcUncertainties (h_jet_pt_ref[iDType][iVar], h2_jet_pt_cov_ref[iDType][iVar], h_jet_counts_ref[iDType][iVar]);
        CalcUncertainties (h_jet_pt_ref[iDType][iVar], h2_jet_pt_cov_ref[iDType][iVar], h_evt_counts_ref[iDType][iVar]);

        h_jet_pt_ref[iDType][iVar]->Scale (1./h_jet_pt_ref[iDType][iVar]->Integral (), "width");

        if (iDType == 0 && (iVar == 0)) {
          TH1D* h = h_jet_pt_ref[iDType][iVar];
          TH2D* h2 = h2_jet_pt_cov_ref[iDType][iVar];

          //float normFactor = 0;
          //float minjpt = (strcmp (tag, "30GeVJets") == 0 ? 30 : 60);
          //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          //  if (h->GetBinLowEdge (iX) > minjpt)
          //    normFactor += h->GetBinContent (iX) * h->GetBinWidth (iX);
          //}
          //h->Scale (1./normFactor);
          //h2->Scale (1./std::pow (normFactor, 2));
          //h->Scale (nEvts/nJets);
          //h2->Scale (std::pow (nEvts/nJets, 2));

          float norm = 0, avgptj = 0, avgptjerr = 0;
          norm = h->Integral ();
          std::cout << "norm = " << norm << std::endl;
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            avgptj += h->GetBinCenter (iX) * h->GetBinContent (iX) / norm;
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
              avgptjerr += h->GetBinCenter (iX) * h2->GetBinContent (iX, iY) * h->GetBinCenter (iY) / (norm*norm);
            }
          }
          std::cout << "avg pt = " << avgptj << " +/-" << std::sqrt (avgptjerr) << std::endl;
        }

        //RebinSomeBins (&(h_jet_pt_ref[iDType][iVar]), nAltPtJBins, (double*) altPtJBins, true);

      }



      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        TString inFileName = Form ("%s/JetPtWeights/%s/%s16_5TeV_%s%s.root", rootPath.Data (), var.Data (), dType.Data (), cent.Data (), iDType == 1 && doJZ123 ? "_JZ123" : "");
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts[iDType][iCent][iVar]   = (TH1D*) inFile->Get ("h_evt_counts")->Clone (Form ("h_evt_counts_pPb_%s_%s_%s", cent.Data(), dType.Data (), var.Data ()));

        h_jet_pt[iDType][iCent][iVar]       = (TH1D*) inFile->Get ("h_jet_pt")->Clone (Form ("h_jet_pt_pPb_%s_%s_%s",       cent.Data (), dType.Data (), var.Data ()));
        h2_jet_pt_cov[iDType][iCent][iVar]  = (TH2D*) inFile->Get ("h2_jet_pt_cov")->Clone (Form ("h2_jet_pt_cov_pPb_%s_%s_%s",  cent.Data (), dType.Data (), var.Data ()));

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        //  h_jet_counts[iDType][iPtJ][iCent][iVar]   = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s16",  pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_pPb_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

        //  h2_jet_eta_phi[iDType][iPtJ][iCent][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_%s16", pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_eta_phi_pPb_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

        //} // end loop over iPtJ

        inFile->Close ();
    
        //const double nEvts = h_evt_counts[iDType][iCent][iVar]->GetBinContent (2); // total number of accepted evts
        //const double nJets = h_jet_counts[iDType][iCent][iVar]->GetBinContent (2); // total number of triggered jets

        //CalcUncertainties (h_jet_pt[iDType][iCent][iVar], h2_jet_pt_cov[iDType][iCent][iVar], h_jet_counts[iDType][iCent][iVar]);
        CalcUncertainties (h_jet_pt[iDType][iCent][iVar], h2_jet_pt_cov[iDType][iCent][iVar], h_evt_counts[iDType][iCent][iVar]);

        h_jet_pt[iDType][iCent][iVar]->Scale (1./h_jet_pt[iDType][iCent][iVar]->Integral (), "width");

        if (iDType == 0 && (iVar == 0)) {
          TH1D* h = h_jet_pt[iDType][iCent][iVar];
          TH2D* h2 = h2_jet_pt_cov[iDType][iCent][iVar];

          //float normFactor = 0;
          //float minjpt = (strcmp (tag, "30GeVJets") == 0 ? 30 : 60);
          //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          //  if (h->GetBinLowEdge (iX) > minjpt)
          //    normFactor += h->GetBinContent (iX) * h->GetBinWidth (iX);
          //}
          //h->Scale (1./normFactor);
          //h2->Scale (1./std::pow (normFactor, 2));
          //h->Scale (nEvts/nJets);
          //h2->Scale (std::pow (nEvts/nJets, 2));

          float norm = 0, avgptj = 0, avgptjerr = 0;
          norm = h->Integral ();
          std::cout << "norm = " << norm << std::endl;
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            avgptj += h->GetBinCenter (iX) * h->GetBinContent (iX) / norm;
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
              avgptjerr += h->GetBinCenter (iX) * h2->GetBinContent (iX, iY) * h->GetBinCenter (iY) / (norm*norm);
            }
          }
          std::cout << "avg pt = " << avgptj << " +/-" << std::sqrt (avgptjerr) << std::endl;
        }

        //RebinSomeBins (&(h_jet_pt[iDType][iCent][iVar]), nAltPtJBins, (double*) altPtJBins, true);

      } // end loop over iCent



      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        h_jet_pt_ratio[iDType][iCent][iVar] = (TH1D*) h_jet_pt[iDType][iCent][iVar]->Clone (Form ("h_jet_pt_ratio_%s_%s_%s", cent.Data (), dType.Data (), var.Data ()));
        h_jet_pt_ratio[iDType][iCent][iVar]->Divide (h_jet_pt_ref[iDType][iVar]);

      } // end loop over iCent

    } // end loop over iVar

  } // end loop over iDType




  {
    //const TString funcStr = "[0]+[1]*pow(log(x),-1)+[2]*pow(log(x),-2)+[3]*pow(log(x),-3)+[4]*pow(log(x),-4)";
    //const TString funcStr = "[0]+[1]*pow(log(x),-1)+[2]*pow(log(x),-2)+[3]*pow(log(x),-3)+[4]*pow(log(x),-4)+[5]*pow(log(x),-5)+[6]*pow(log(x),-6)";
    //const TString funcStr = "[0]+[1]*pow(log(x),-1)+[2]*pow(log(x),-2)";
    //const TString funcStr = "(exp([0]*x) + [1]*pow(x,[2])) / ([3]*exp([4]*x) + [5]*pow(x,[6]))";
    //const TString funcStr = "[0] * exp([1]*x) + [2] * pow (x, [3])";
    const TString funcStr = "[0]+exp([1]+[2]*log(x))";

    //for (int iVar = 0; iVar < nVar; iVar++) {
    for (int iVar : {0}) {

      const TString var = variations[iVar];

      {
        TF1* f_mb   = new TF1 ("f_mb",  "exp([0]+[1]*log(x))", 30, 60);
        TF1* f_jet  = new TF1 ("f_jet", "exp([0]+[1]*log(x))", 60, 120);
        f_mb->SetParameter (0, std::log (h_jet_pt_ref[1][iVar]->GetBinContent (1)));
        f_mb->SetParameter (1, (std::log (h_jet_pt_ref[1][iVar]->GetBinContent (2)) - std::log (h_jet_pt_ref[1][iVar]->GetBinContent (1))) / (0.5 * h_jet_pt_ref[1][iVar]->GetBinWidth (1) + 0.5 * h_jet_pt_ref[1][iVar]->GetBinWidth (2)));
        f_jet->SetParameter (0, std::log (h_jet_pt_ref[1][iVar]->GetBinContent (16)));
        f_jet->SetParameter (1, (std::log (h_jet_pt_ref[1][iVar]->GetBinContent (17)) - std::log (h_jet_pt_ref[1][iVar]->GetBinContent (16))) / (0.5 * h_jet_pt_ref[1][iVar]->GetBinWidth (16) + 0.5 * h_jet_pt_ref[1][iVar]->GetBinWidth (17)));

        h_jet_pt_ref[1][iVar]->Fit (f_mb, "RN0QI");
        h_jet_pt_ref[1][iVar]->Fit (f_jet, "RN0QI");

        const double sf = f_mb->Eval (60) / f_jet->Eval (60);
        std::cout << "sf = " << f_mb->Eval (60) << " / " << f_jet->Eval (60) << " = " << sf << std::endl;

        for (int iX = 1; iX <= h_jet_pt_ref[1][iVar]->GetNbinsX (); iX++) {
          if (h_jet_pt_ref[1][iVar]->GetBinCenter (iX) < 60) continue;

          h_jet_pt_ref[1][iVar]->SetBinContent (iX, h_jet_pt_ref[1][iVar]->GetBinContent (iX) * sf);
        }

        h_jet_pt_datamc_ratio_ref[iVar] = (TH1D*) h_jet_pt_ref[0][iVar]->Clone (Form ("h_jet_pt_datamc_ratio_ref_%s", var.Data ()));
        h_jet_pt_datamc_ratio_ref[iVar]->Divide (h_jet_pt_ref[1][iVar]);

        //TF1* f_mb   = new TF1 ("f_mb",  "[0]+[1]*log(x)", 30, 60);
        //TF1* f_jet  = new TF1 ("f_jet", "[0]+[1]*log(x)", 60, 120);

        //h_jet_pt_datamc_ratio_ref[iVar]->Fit (f_mb, "RN0QI");
        //h_jet_pt_datamc_ratio_ref[iVar]->Fit (f_jet, "RN0QI");

        //const double sf = f_mb->Eval (60) / f_jet->Eval (60);
        //std::cout << "sf = " << f_mb->Eval (60) << " / " << f_jet->Eval (60) << " = " << sf << std::endl;

        //for (int iX = 1; iX <= h_jet_pt_datamc_ratio_ref[iVar]->GetNbinsX (); iX++) {
        //  if (h_jet_pt_datamc_ratio_ref[iVar]->GetBinCenter (iX) < 60) continue;

        //  h_jet_pt_datamc_ratio_ref[iVar]->SetBinContent (iX, h_jet_pt_datamc_ratio_ref[iVar]->GetBinContent (iX) * sf);
        //}

        SaferDelete (&f_mb);
        SaferDelete (&f_jet);


        //h_jet_pt_datamc_ratio_ref[iVar]->Scale (h_jet_pt_ref[1][iVar]->Integral () / h_jet_pt_ref[0][iVar]->Integral ());

        TF1* f = new TF1 (Form ("f_jet_pt_datamc_ratio_ref_%s", var.Data ()), funcStr.Data (), pTJBins[0], pTJBins[nPtJBins]);
        h_jet_pt_datamc_ratio_ref[iVar]->Fit (f, "RN0QI");
        f_jet_pt_datamc_ratio_ref[iVar] = f;
      }

      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        TF1* f_mb   = new TF1 ("f_mb",  "exp([0]+[1]*log(x))", 30, 60);
        TF1* f_jet  = new TF1 ("f_jet", "exp([0]+[1]*log(x))", 60, 120);
        f_mb->SetParameter (0, std::log (h_jet_pt[1][iCent][iVar]->GetBinContent (1)));
        f_mb->SetParameter (1, (std::log (h_jet_pt[1][iCent][iVar]->GetBinContent (2)) - std::log (h_jet_pt[1][iCent][iVar]->GetBinContent (1))) / (0.5 * h_jet_pt[1][iCent][iVar]->GetBinWidth (1) + 0.5 * h_jet_pt[1][iCent][iVar]->GetBinWidth (2)));
        f_jet->SetParameter (0, std::log (h_jet_pt[1][iCent][iVar]->GetBinContent (16)));
        f_jet->SetParameter (1, (std::log (h_jet_pt[1][iCent][iVar]->GetBinContent (17)) - std::log (h_jet_pt[1][iCent][iVar]->GetBinContent (16))) / (0.5 * h_jet_pt[1][iCent][iVar]->GetBinWidth (16) + 0.5 * h_jet_pt[1][iCent][iVar]->GetBinWidth (17)));

        h_jet_pt[1][iCent][iVar]->Fit (f_mb, "RN0QI");
        h_jet_pt[1][iCent][iVar]->Fit (f_jet, "RN0QI");

        const double sf = f_mb->Eval (60) / f_jet->Eval (60);
        std::cout << "sf = " << f_mb->Eval (60) << " / " << f_jet->Eval (60) << " = " << sf << std::endl;
  
        for (int iX = 1; iX <= h_jet_pt[1][iCent][iVar]->GetNbinsX (); iX++) {
          if (h_jet_pt[1][iCent][iVar]->GetBinCenter (iX) < 60) continue;
  
          h_jet_pt[1][iCent][iVar]->SetBinContent (iX, h_jet_pt[1][iCent][iVar]->GetBinContent (iX) * sf);
        }

        h_jet_pt_datamc_ratio[iCent][iVar] = (TH1D*) h_jet_pt[0][iCent][iVar]->Clone (Form ("h_jet_pt_datamc_ratio_%s_%s", cent.Data (), var.Data ()));
        h_jet_pt_datamc_ratio[iCent][iVar]->Divide (h_jet_pt[1][iCent][iVar]);

        //TF1* f_mb   = new TF1 ("f_mb",  "[0]+[1]*log(x)", 30, 60);
        //TF1* f_jet  = new TF1 ("f_jet", "[0]+[1]*log(x)", 60, 120);
  
        //h_jet_pt_datamc_ratio[iCent][iVar]->Fit (f_mb, "RN0QI");
        //h_jet_pt_datamc_ratio[iCent][iVar]->Fit (f_jet, "RN0QI");
  
        //const double sf = f_mb->Eval (60) / f_jet->Eval (60);
        //std::cout << "sf = " << f_mb->Eval (60) << " / " << f_jet->Eval (60) << " = " << sf << std::endl;
  
        //for (int iX = 1; iX <= h_jet_pt_datamc_ratio[iCent][iVar]->GetNbinsX (); iX++) {
        //  if (h_jet_pt_datamc_ratio[iCent][iVar]->GetBinCenter (iX) < 60) continue;
  
        //  h_jet_pt_datamc_ratio[iCent][iVar]->SetBinContent (iX, h_jet_pt_datamc_ratio[iCent][iVar]->GetBinContent (iX) * sf);
        //}
  
        SaferDelete (&f_mb);
        SaferDelete (&f_jet);

        //h_jet_pt_datamc_ratio[iCent][iVar]->Scale (h_jet_pt_ref[1][iVar]->Integral () / h_jet_pt_ref[0][iVar]->Integral ());

        TF1* f = new TF1 (Form ("f_jet_pt_datamc_ratio_%s_%s", cent.Data (), var.Data ()), funcStr.Data (), pTJBins[0], pTJBins[nPtJBins]);
        h_jet_pt_datamc_ratio[iCent][iVar]->Fit (f, "RN0QI");
        f_jet_pt_datamc_ratio[iCent][iVar] = f;

      } // end loop over iCent

    } // end loop over iVar
  }



  {
    g_jet_pt_ref_syst[0] = make_graph (h_jet_pt_ref[0][0]);

    ResetTGAEErrors (g_jet_pt_ref_syst[0]);

    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      g_jet_pt_syst[iCent][0]       = make_graph (h_jet_pt[0][iCent][0]);
      g_jet_pt_ratio_syst[iCent][0] = make_graph (h_jet_pt_ratio[0][iCent][0]);

      ResetTGAEErrors (g_jet_pt_syst[iCent][0]);
      ResetTGAEErrors (g_jet_pt_ratio_syst[iCent][0]);

    } // end loop over iCent
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// CREATE GRAPHS THAT WILL STORE SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE SEPARATELY.
  //// THEN CALCULATES SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE BY TAKING DIFFERENCES
  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //for (int iDType = 0; iDType < 2; iDType++) {

  //  for (int iVar = 1; iVar < nVar; iVar++) {

  //    const TString var = variations[iVar];

  //    if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
  //      continue;

  //    g_jet_pt_ref_syst[iVar] = new TGAE ();

  //    CalcSystematics (g_jet_pt_ref_syst[iVar], h_jet_pt_ref[iDType][0], h_jet_pt_ref[iDType][iVar]);

  //    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

  //      g_jet_pt_syst[iCent][iVar] = new TGAE ();
  //      g_jet_pt_ratio_syst[iCent][iVar] = new TGAE ();

  //      CalcSystematics (g_jet_pt_syst[iCent][iVar],        h_jet_pt[iDType][iCent][0],       h_jet_pt[iDType][iCent][iVar]);
  //      CalcSystematics (g_jet_pt_ratio_syst[iCent][iVar],  h_jet_pt_ratio[iDType][iCent][0], h_jet_pt_ratio[iDType][iCent][iVar]);

  //      // allow some uncertainties to not cancel in the p+Pb / pp ratio by overwriting the current uncertainties
  //      if (variationsThatDontCancelInRatio.count (var) != 0) {
  //        ResetTGAEErrors (g_jet_pt_ratio_syst[iCent][iVar]);
  //        AddRelErrorsInQuadrature (g_jet_pt_ratio_syst[iCent][iVar], g_jet_pt_ref_syst[iVar], false, true);
  //        AddRelErrorsInQuadrature (g_jet_pt_ratio_syst[iCent][iVar], g_jet_pt_syst[iCent][iVar], false, true);
  //      }

  //    } // end loop over iCent

  //  } // end loop over iVar

  //} // end loop over iDType


  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// SYSTEMATIC UNCERTAINTIES DERIVED IN MC MUST HAVE CENTRAL VALUES SET BY CENTRAL VALUES IN DATA
  //// THE FINAL UNCERTAINTY IS ASSIGNED TO MATCH THE FRACTIONAL UNCERTAINTY IN MC
  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //for (int iVar = 1; iVar < nVar; iVar++) {

  //  const TString var = variations[iVar];

  //  if (dataVariations.count (var) > 0 || mcVariations.count (var) == 0)
  //    continue; // skip variations already evaluated in data or that are not evaluated in MC

  //  for (int iDir = 0; iDir < nDir; iDir++) {

  //    SetCentralValuesKeepRelativeErrors (g_jet_pt_ref_syst[iVar], h_jet_pt_ref[0][0]);

  //    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

  //      SetCentralValuesKeepRelativeErrors (g_jet_pt_syst[iCent][iVar],       h_jet_pt[0][iCent][0]);
  //      SetCentralValuesKeepRelativeErrors (g_jet_pt_ratio_syst[iCent][iVar], h_jet_pt_ratio[0][iCent][0]);

  //    } // end loop over iCent

  //  } // end loop over iDir

  //} // end loop over iVar



  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// ADD SYSTEMATIC UNCERTAINTIES FROM ALL SOURCES IN QUADRATURE, STORING RESULTS IN A SINGLE GRAPH
  ////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //for (int iVar = 1; iVar < nVar; iVar++) {

  //  const TString var = variations[iVar];

  //  //if (dataVariations.count (var) == 0 && mcVariations.count (var) == 0)
  //  if (mcVariations.count (var) == 0)
  //    continue;

  //  for (int iDir = 0; iDir < nDir; iDir++) {
  //
  //    AddErrorsInQuadrature (g_jet_pt_ref_syst[0], g_jet_pt_ref_syst[iVar], false, true);
  //
  //    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  //
  //      AddErrorsInQuadrature (g_jet_pt_syst[iCent][0],       g_jet_pt_syst[iCent][iVar], false, true);
  //      AddErrorsInQuadrature (g_jet_pt_ratio_syst[iCent][0], g_jet_pt_ratio_syst[iCent][iVar], false, true);
  //
  //    } // end loop over iCent
  //
  //  } // end loop over iDir

  //} // end loop over iVar
  


  {
    outFile->cd ();


    //for (int iVar = 0; iVar < nVar; iVar++) {
    for (int iVar : {0}) {

      const TString var = variations[iVar];

      //if (dataVariations.count (var) == 0 && mcVariations.count (var) == 0)
      if (mcVariations.count (var) == 0)
        continue;

      g_jet_pt_ref_syst[iVar]->Write (Form ("g_jet_pt_ref_syst_%s", var.Data ()));

      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        g_jet_pt_syst[iCent][iVar]->Write (Form ("g_jet_pt_syst_pPb_%s_%s", cent.Data (), var.Data ()));
        g_jet_pt_ratio_syst[iCent][iVar]->Write (Form ("g_jet_pt_ratio_syst_%s_%s", cent.Data (), var.Data ()));

      } // end loop over iCent

    } // end loop over iVar


    for (int iDType = 0; iDType < 2; iDType++) {

      //for (int iVar = 0; iVar < nVar; iVar++) {
      for (int iVar : {0}) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;

        h_evt_counts_ref[iDType][iVar]->Write ();

        h_jet_pt_ref[iDType][iVar]->Write ();

        for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          h_evt_counts[iDType][iCent][iVar]->Write ();

          h_jet_pt[iDType][iCent][iVar]->Write ();
          h_jet_pt_ratio[iDType][iCent][iVar]->Write ();

        } // end loop over iCent

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  h_jet_counts_ref[iDType][iPtJ][iVar]->Write ();

        //  h2_jet_eta_phi_ref[iDType][iPtJ][iVar]->Write ();

        //  for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        //    h_jet_counts[iDType][iPtJ][iCent][iVar]->Write ();

        //    h2_jet_eta_phi[iDType][iPtJ][iCent][iVar]->Write ();

        //  } // end loop over iCent

        //} // end loop over iPtJ

      } // end loop over iVar

    } // end loop over iDType


    //for (int iVar = 0; iVar < nVar; iVar++) {
    for (int iVar : {0}) {

      h_jet_pt_datamc_ratio_ref[iVar]->Write ();
      f_jet_pt_datamc_ratio_ref[iVar]->Write ();

      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        h_jet_pt_datamc_ratio[iCent][iVar]->Write ();
        f_jet_pt_datamc_ratio[iCent][iVar]->Write ();

      } // end loop over iCent

    } // end loop over iVar



    TFile* wgtsFile = new TFile (Form ("%s/aux/JetPtWeights.root", workPath.Data ()), "recreate");
    wgtsFile->cd ();

    //for (int iVar = 0; iVar < nVar; iVar++) {
    for (int iVar : {0}) {

      h_jet_pt_datamc_ratio_ref[iVar]->Write ();
      f_jet_pt_datamc_ratio_ref[iVar]->Write ();

      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        h_jet_pt_datamc_ratio[iCent][iVar]->Write ();
        f_jet_pt_datamc_ratio[iCent][iVar]->Write ();

      } // end loop over iCent

    } // end loop over iVar


    outFile->Close ();
    wgtsFile->Close ();

  }
}


#endif
