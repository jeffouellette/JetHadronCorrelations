#ifndef __JetHadronCorrelatorPlotNIters_C__
#define __JetHadronCorrelatorPlotNIters_C__

#include <iostream>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

#include <ArrayTemplates.h>
#include <Utilities.h>
#include <MyStyle.h>
#include <MyColors.h>

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"


using namespace JetHadronCorrelations;


void PlotNIters (const char* inFileTag) {

  TLine* l = new TLine ();
  TLatex* tl = new TLatex ();

  TFile* inFile = nullptr;

  const short nItersMin = 1;
  const short nItersMax = 20;
  const double* nItersVals = linspace (nItersMin, nItersMax, nItersMax-nItersMin);

  TH1D**    h_jet_pt_ref              = Get1DArray <TH1D*> (2);
  TH1D***   h_jet_pt                  = Get2DArray <TH1D*> (2, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_sig   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_sig       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_unf   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_unf       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_unf_nIters = Get3DArray <TH1D*> (nPtJBins, nDir, nItersMax-nItersMin+2);
  TH1D***** h_jetInt_trk_pt_unf_nIters     = Get4DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1, nItersMax-nItersMin+2);

  TH1D***** h_jetInt_trk_pt_iaa           = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);
  TH1D***** h_jetInt_trk_pt_iaaNoUnf      = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D**    h_jet_pt_ref_unf_nIters       = Get1DArray <TH1D*> (nItersMax-nItersMin+2);
  TH1D***   h_jet_pt_unf_nIters           = Get2DArray <TH1D*> (nZdcCentBins+1, nItersMax-nItersMin+2);

  TGraph**    g_jet_pt_ref_unfIterUnc     = Get1DArray <TGraph*> (2);                 // sums of iterations uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfIterUnc         = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of iterations uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfSumUnc      = Get1DArray <TGraph*> (2);                 // sums of statistical uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfSumUnc          = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of statistical uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfTotUnc      = Get1DArray <TGraph*> (2);                 // sums of total uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfTotUnc          = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of total uncertainties as a function of nIter -- data only

  TGraph**    g_jet_pt_ref_unfIterRelUnc  = Get1DArray <TGraph*> (2);                 // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfIterRelUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfSumRelUnc   = Get1DArray <TGraph*> (2);                 // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfSumRelUnc       = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfTotRelUnc   = Get1DArray <TGraph*> (2);                 // sums of total relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfTotRelUnc       = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of total relative uncertainties as a function of nIter -- data only


  TGraph***   g_jetInt_trk_pt_ref_unfIterUnc    = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterUnc        = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfSumUnc     = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfSumUnc         = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotUnc     = Get2DArray <TGraph*> (2, nDir);                           // sums of total uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotUnc         = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total uncertainties as a function of nIter -- data only

  TGraph***   g_jetInt_trk_pt_ref_unfIterRelUnc = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterRelUnc     = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfSumRelUnc  = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfSumRelUnc      = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotRelUnc  = Get2DArray <TGraph*> (2, nDir);                           // sums of total relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotRelUnc      = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total relative uncertainties as a function of nIter -- data only


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (short iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      h_jet_pt_ref[iDType] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_%s_Nominal", dType.Data ()));

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        h_jet_pt[iDType][iCent] = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_%s_Nominal", cent.Data (), dType.Data ()));

      } // end loop over iCent

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));
          h_jetInt_trk_pt_ref_unf[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));

        } // end loop over iDir


        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_unf[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_iaa[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaa_%s_%s_%s_Nominal",     dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));
            h_jetInt_trk_pt_iaaNoUnf[iDType][iPtJInt][iDir][iCent]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_iaaNoUnf_%s_%s_%s_Nominal",  dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

          } // end loop over iDir

        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType


    for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
      const short nIters = (short) nItersVals[iIter];

      h_jet_pt_ref_unf_nIters[iIter] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_unf_data_Nominal_nIters%i", nIters));

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        h_jet_pt_unf_nIters[iCent][iIter] = (TH1D*) inFile-> Get (Form ("h_jet_pt_unf_data_%s_Nominal_nIters%i", cent.Data (), nIters));
      } // end loop over iCent

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_data_%s_Nominal_nIters%i", dir.Data (), pTJInt.Data (), nIters));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
  
            h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_data_%s_Nominal_nIters%i", dir.Data (), cent.Data (), pTJInt.Data (), nIters));
  
          } // end loop over iCent

        } // end loop over iDir

      } // end loop over iPtJInt

    } // end loop over iIter

  }




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS SUMMARIZING UNCERTAINTIES VS. # OF ITERATIONS ON THE UNFOLDING
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    const short iDType = 0;
    const TString dType = "data";

    for (short iPtJInt : {0, 1}) {

      const double minJetPt = (iPtJInt == 0 ? 30. : 60.);

      TH1D* h_unf = h_jet_pt_ref[iDType];
      TH1D* h_unf_prev = h_jet_pt_ref[iDType];

      g_jet_pt_ref_unfSumUnc[iPtJInt]      = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfIterUnc[iPtJInt]     = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfTotUnc[iPtJInt]      = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfSumRelUnc[iPtJInt]   = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfIterRelUnc[iPtJInt]  = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfTotRelUnc[iPtJInt]   = new TGraph (nItersMax - nItersMin + 1);


      {
        double tot = 0, totErr = 0;
        for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
          if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
          tot += h_unf->GetBinContent (iX);
          totErr += std::pow (h_unf->GetBinError (iX), 2);
        } // end loop over iX
        double totRelErr = (tot > 0 ? totErr / (tot*tot) : 0);

        g_jet_pt_ref_unfSumUnc[iPtJInt]->SetPoint    (g_jet_pt_ref_unfSumUnc[iPtJInt]->GetN (),     0, std::sqrt (totErr));
        g_jet_pt_ref_unfTotUnc[iPtJInt]->SetPoint    (g_jet_pt_ref_unfTotUnc[iPtJInt]->GetN (),     0, std::sqrt (totErr));
        g_jet_pt_ref_unfSumRelUnc[iPtJInt]->SetPoint (g_jet_pt_ref_unfSumRelUnc[iPtJInt]->GetN (),  0, std::sqrt (totRelErr));
        g_jet_pt_ref_unfTotRelUnc[iPtJInt]->SetPoint (g_jet_pt_ref_unfTotRelUnc[iPtJInt]->GetN (),  0, std::sqrt (totRelErr));
      }


      for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {

        const short nIters = (short) nItersVals[iIter];

        h_unf = h_jet_pt_ref_unf_nIters[iIter];

        double tot = 0, totErr = 0;
        for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
          if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
          tot += h_unf->GetBinContent (iX);
          totErr += std::pow (h_unf->GetBinError (iX), 2);
        } // end loop over iX
        double totRelErr = (tot > 0 ? totErr / (tot*tot) : 0);

        double iterErr = 0;
        for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
          if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
          iterErr += std::pow (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX), 2);
        } // end loop over iX
        double iterRelErr = (tot > 0 ? iterErr / (tot*tot) : 0);

        g_jet_pt_ref_unfSumUnc[iPtJInt]->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totErr));
        g_jet_pt_ref_unfIterUnc[iPtJInt]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (iterErr));
        g_jet_pt_ref_unfTotUnc[iPtJInt]->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totErr + iterErr));
        g_jet_pt_ref_unfSumRelUnc[iPtJInt]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelErr));
        g_jet_pt_ref_unfIterRelUnc[iPtJInt]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (iterRelErr));
        g_jet_pt_ref_unfTotRelUnc[iPtJInt]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelErr + iterRelErr));

        h_unf_prev = h_unf;
      } // end loop over iIter

    } // end loop over iPtJInt


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      for (short iPtJInt : {0, 1}) {

        const double minJetPt = (iPtJInt == 0 ? 30. : 60.);

        TH1D* h_unf = h_jet_pt[iDType][iCent];
        TH1D* h_unf_prev = h_jet_pt[iDType][iCent];

        g_jet_pt_unfSumUnc[iPtJInt][iCent]     = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfIterUnc[iPtJInt][iCent]    = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfTotUnc[iPtJInt][iCent]     = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfSumRelUnc[iPtJInt][iCent]  = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfIterRelUnc[iPtJInt][iCent] = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfTotRelUnc[iPtJInt][iCent]  = new TGraph (nItersMax - nItersMin + 1);


        {
          double tot = 0, totErr = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
            tot += h_unf->GetBinContent (iX);
            totErr += std::pow (h_unf->GetBinError (iX), 2);
          } // end loop over iX
          double totRelErr = (tot > 0 ? totErr / (tot*tot) : 0);

          g_jet_pt_unfSumUnc[iPtJInt][iCent]->SetPoint     (g_jet_pt_unfSumUnc[iPtJInt][iCent]->GetN (),    0, std::sqrt (totErr));
          g_jet_pt_unfTotUnc[iPtJInt][iCent]->SetPoint     (g_jet_pt_unfTotUnc[iPtJInt][iCent]->GetN (),    0, std::sqrt (totErr));
          g_jet_pt_unfSumRelUnc[iPtJInt][iCent]->SetPoint  (g_jet_pt_unfSumRelUnc[iPtJInt][iCent]->GetN (), 0, std::sqrt (totRelErr));
          g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->SetPoint  (g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->GetN (), 0, std::sqrt (totRelErr));
        }

 
        for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
          const short nIters = (short) nItersVals[iIter];

          h_unf = h_jet_pt_unf_nIters[iCent][iIter];

          double tot = 0, totErr = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
            tot += h_unf->GetBinContent (iX);
            totErr += std::pow (h_unf->GetBinError (iX), 2);
          } // end loop over iX
          double totRelErr = (tot > 0 ? totErr / (tot*tot) : 0);
  
          double iterErr = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || 300 < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
            iterErr += std::pow (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX), 2);
          } // end loop over iX
          double iterRelErr = (tot > 0 ? iterErr / (tot*tot) : 0);
  
          g_jet_pt_unfSumUnc[iPtJInt][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totErr));
          g_jet_pt_unfIterUnc[iPtJInt][iCent]->SetPoint    (nIters - nItersMin, nIters, std::sqrt (iterErr));
          g_jet_pt_unfTotUnc[iPtJInt][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totErr + iterErr));
          g_jet_pt_unfSumRelUnc[iPtJInt][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelErr));
          g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->SetPoint (nIters - nItersMin, nIters, std::sqrt (iterRelErr));
          g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelErr + iterRelErr));

          h_unf_prev = h_unf;

        } // end loop over iIter

      } // end loop over iPtJInt

    } // end loop over iCent


    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      for (short iPtJInt : {0, 1}) {

        const double maxPtCh = (iPtJInt == 0 ? 40 : 75);

        TH1D* h_unf = h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir];
        TH1D* h_unf_prev = h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir];

        g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]      = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]     = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]      = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]   = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]  = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]   = new TGraph (nItersMax - nItersMin + 1);


        {
          double tot = 0, totErr = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
            if (ptch < 5 || maxPtCh < ptch) continue;
            tot += h_unf->GetBinContent (iX) * h_unf->GetBinWidth (iX);
            totErr += std::pow (h_unf->GetBinError (iX) * h_unf->GetBinWidth (iX), 2);
          } // end loop over iX
          double totRelErr = (tot > 0 ? totErr / (tot*tot) : 0);

          g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]->SetPoint    (g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]->GetN (),     0, std::sqrt (totErr));
          g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->SetPoint    (g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->GetN (),     0, std::sqrt (totErr));
          g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]->SetPoint (g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]->GetN (),  0, std::sqrt (totRelErr));
          g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->SetPoint (g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->GetN (),  0, std::sqrt (totRelErr));
        }


        for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
          const short nIters = (short) nItersVals[iIter];

          h_unf = h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter];

          double tot = 0, totErr = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
            if (ptch < 5 || maxPtCh < ptch) continue;
            tot += h_unf->GetBinContent (iX) * h_unf->GetBinWidth (iX);
            totErr += std::pow (h_unf->GetBinError (iX) * h_unf->GetBinWidth (iX), 2);
          } // end loop over iX
          double totRelErr = (tot > 0 ? totErr / (tot*tot) : 0);
  
          double iterErr = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
            if (ptch < 5 || maxPtCh < ptch) continue;
            iterErr += std::pow ((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) * h_unf->GetBinWidth (iX), 2);
          } // end loop over iX
          double iterRelErr = (tot > 0 ? iterErr / (tot*tot) : 0);
  
          g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir]->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totErr));
          g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (iterErr));
          g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totErr + iterErr));
          g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelErr));
          g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (iterRelErr));
          g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelErr + iterRelErr));

          h_unf_prev = h_unf;

        } // end loop over iIter

      } // end loop over iPtJInt

    } // end loop over iDir


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent)); 

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (short iPtJInt : {0, 1}) {

          const double maxPtCh = (iPtJInt == 0 ? 40 : 75);
  
          TH1D* h_unf = h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent];
          TH1D* h_unf_prev = h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent];

          g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]     = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]    = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]     = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]  = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent] = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]  = new TGraph (nItersMax - nItersMin + 1);


          {
            double tot = 0, totErr = 0;
            for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
              const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
              if (ptch < 5 || maxPtCh < ptch) continue;
              tot += h_unf->GetBinContent (iX) * h_unf->GetBinWidth (iX);
              totErr += std::pow (h_unf->GetBinError (iX) * h_unf->GetBinWidth (iX), 2);
            } // end loop over iX
            double totRelErr = (tot > 0 ? totErr / (tot*tot) : 0);

            g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]->SetPoint     (g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]->GetN (),    0, std::sqrt (totErr));
            g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->SetPoint     (g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->GetN (),    0, std::sqrt (totErr));
            g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]->SetPoint  (g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]->GetN (), 0, std::sqrt (totRelErr));
            g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->SetPoint  (g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->GetN (), 0, std::sqrt (totRelErr));
          }


          for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {

            const short nIters = (short) nItersVals[iIter];

            h_unf = h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter];

            double tot = 0, totErr = 0;
            for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
              const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
              if (ptch < 5 || maxPtCh < ptch) continue;
              tot += h_unf->GetBinContent (iX) * h_unf->GetBinWidth (iX);
              totErr += std::pow (h_unf->GetBinError (iX) * h_unf->GetBinWidth (iX), 2);
            } // end loop over iX
            double totRelErr = (tot > 0 ? totErr / (tot*tot) : 0);
    
            double iterErr = 0;
            for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
              const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
              if (ptch < 5 || maxPtCh < ptch) continue;
              iterErr += std::pow ((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) * h_unf->GetBinWidth (iX), 2);
            } // end loop over iX
            double iterRelErr = (tot > 0 ? iterErr / (tot*tot) : 0);

            g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totErr));
            g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->SetPoint    (nIters - nItersMin, nIters, std::sqrt (iterErr));
            g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totErr + iterErr));
            g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelErr));
            g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->SetPoint (nIters - nItersMin, nIters, std::sqrt (iterRelErr));
            g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelErr + iterRelErr));

            h_unf_prev = h_unf;

          } // end loop over iPtJInt

        } // end loop over iIter

      } // end loop over iDir

    } // end loop over iCent

  } // end iDType = 0



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    const char* canvasName = Form ("c_jet_pt_unfUncs_%s", pTJInt.Data ());
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);

    TGraph* g = nullptr;

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);

    double ymax = 0, ymin = 0;
    {
      c->cd (7);

      gPad->SetLogy ();

      double x, y, xmin;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_ref_unfTotUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) {
          xmin = x; ymin = y;
        }
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= 5;
      
      g = g_jet_pt_ref_unfSumUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_ref_unfIterUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      ymin *= 0.2;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{J}^{2}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_ref_unfTotUnc[iPtJInt];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfSumUnc[iPtJInt];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfIterUnc[iPtJInt];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xmin, ymin, xmin, ymax/5);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      gPad->SetLogy ();

      double x, y, xmin;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_unfTotUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) {
          xmin = x; ymin = y;
        }
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= 5;

      g = g_jet_pt_unfSumUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_unfIterUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      ymin *= 0.2;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{J}^{2}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_unfTotUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfSumUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfIterUnc[iPtJInt][iCent];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);
      l->DrawLine (xmin, ymin, xmin, ymax/5);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
    myLineText2 (0.15, 0.56, kBlack,        kOpenCircle, "Stat. unc. #oplus Iter. diff.", 1.2, 0.06);
    myLineText2 (0.15, 0.48, kBlue,         kOpenCircle, "Stat. unc. only", 1.2, 0.06);
    myLineText2 (0.15, 0.40, kRed,          kOpenCircle, "Iter. diff. only, i #rightarrow i+1", 1.2, 0.06);

    c->SaveAs (Form ("%s/Plots/PtCh/UnfUncs_Summary_JetSpectra_%s.pdf", workPath.Data (), pTJInt.Data ()));

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    const char* canvasName = Form ("c_jet_pt_unfRelUncs_%s", pTJInt.Data ());
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);

    TGraph* g = nullptr;

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);

    double ymax = 0, ymin = 0;
    {
      c->cd (7);

      //gPad->SetLogy ();

      double x, y, xmin;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_ref_unfTotRelUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) {
          xmin = x; ymin = y;
        }
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= 1.2;

      ymin = 0;  
      //g = g_jet_pt_ref_unfSumRelUnc[iPtJInt];
      //for (short i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  if (x != 0 && y < ymin && y > 0) ymin = y;
      //}
      //g = g_jet_pt_ref_unfIterRelUnc[iPtJInt];
      //for (short i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  if (x != 0 && y < ymin && y > 0) ymin = y;
      //}
      //ymin *= 0.2;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{J}^{2}} / #Sigma N_{J}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_ref_unfTotRelUnc[iPtJInt];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfSumRelUnc[iPtJInt];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfIterRelUnc[iPtJInt];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xmin, ymin, xmin, ymax/1.2);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      //gPad->SetLogy ();

      double x, y, xmin;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_unfTotRelUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) {
          xmin = x; ymin = y;
        }
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= 1.2;

      ymin = 0;  
      //g = g_jet_pt_unfSumRelUnc[iPtJInt][iCent];
      //for (short i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  if (x != 0 && y < ymin && y > 0) ymin = y;
      //}
      //g = g_jet_pt_unfIterRelUnc[iPtJInt][iCent];
      //for (short i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  if (x != 0 && y < ymin && y > 0) ymin = y;
      //}
      //ymin *= 0.2;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{J}^{2}} / #Sigma N_{J}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_unfTotRelUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfSumRelUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfIterRelUnc[iPtJInt][iCent];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);
      l->DrawLine (xmin, ymin, xmin, ymax/1.2);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
    myLineText2 (0.15, 0.56, kBlack,        kOpenCircle, "Stat. unc. #oplus Iter. diff.", 1.2, 0.06);
    myLineText2 (0.15, 0.48, kBlue,         kOpenCircle, "Stat. unc. only", 1.2, 0.06);
    myLineText2 (0.15, 0.40, kRed,          kOpenCircle, "Iter. diff. only, i #rightarrow i+1", 1.2, 0.06);

    c->SaveAs (Form ("%s/Plots/PtCh/UnfRelUncs_Summary_JetSpectra_%s.pdf", workPath.Data (), pTJInt.Data ()));

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_unfUncs_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGraph* g = nullptr;

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);

      double ymax = 0, ymin = 0;
      {
        c->cd (7);

        gPad->SetLogy ();

        double x, y, xmin;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) {
            xmin = x; ymin = y;
          }
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= 5;

        g = g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        ymin *= 0.2;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{Y}^{2}}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfSumUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xmin, ymin, xmin, ymax/5);
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogy ();

        double x, y, xmin;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) {
            xmin = x; ymin = y;
          }
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= 5;

        g = g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        ymin *= 0.2;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{Y}^{2}}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfSumUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);
        l->DrawLine (xmin, ymin, xmin, ymax/5);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (5, %i)", iPtJInt == 0 ? 40 : 75), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kOpenCircle, "Stat. unc. #oplus Iter. diff.", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "Stat. unc. only", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenCircle, "Iter. diff. only, i #rightarrow i+1", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_unfRelUncs_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGraph* g = nullptr;

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);

      double ymax = 0, ymin = 0;
      {
        c->cd (7);

        gPad->SetLogy ();

        double x, y, xmin;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) {
            xmin = x; ymin = y;
          }
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= 5;

        g = g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        ymin *= 0.2;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{Y}^{2}} / #Sigma Y", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfSumRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xmin, ymin, xmin, ymax/5);
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogy ();

        double x, y, xmin;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) {
            xmin = x; ymin = y;
          }
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= 5;

        g = g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        ymin *= 0.2;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{Y}^{2}} / #Sigma Y", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfSumRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);
        l->DrawLine (xmin, ymin, xmin, ymax/5);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (5, %i)", iPtJInt == 0 ? 40 : 75), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kOpenCircle, "Stat. unc. #oplus Iter. diff.", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "Stat. unc. only", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenCircle, "Iter. diff. only, i #rightarrow i+1", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfRelUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  {
    const char* canvasName = "c_jet_pt_sigVsUnf";
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);

    TGAE* g = nullptr;

    {
      c->cd (7);

      gPad->SetLogx ();

      //TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];N_{iter} iterations / No unfold", 1, pTJBins[0], pTJBins[nPtJBins]);
      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];N_{iter} iterations / N_{iter} - 1 iterations", 1, pTJBins[0], pTJBins[nPtJBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.00, 2.00);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);


      for (short iIter : {0, 1, 2, 3, 4, 5}) {
        h = h_jet_pt_ref_unf_nIters[iIter];
        g = make_graph (h);
        ScaleGraph (g, (iIter == 0 ? h_jet_pt_ref[0] : h_jet_pt_ref_unf_nIters[iIter-1]));
        //ScaleGraph (g, h_jet_pt_ref[0]);
        //ResetXErrors (g);
        myDraw (g, colorfulColors[iIter], kOpenCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);
      }

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      gPad->SetLogx ();

      //TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];N_{iter} iterations / No unfold", 1, pTJBins[0], pTJBins[nPtJBins]);
      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];N_{iter} iterations / N_{iter} - 1 iterations", 1, pTJBins[0], pTJBins[nPtJBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.00, 2.00);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      for (short iIter : {0, 1, 2, 3, 4, 5}) {
        h = h_jet_pt_unf_nIters[iCent][iIter];
        g = make_graph (h);
        ScaleGraph (g, (iIter == 0 ? h_jet_pt[0][iCent] : h_jet_pt_unf_nIters[iCent][iIter-1]));
        //ScaleGraph (g, h_jet_pt[0][iCent]);
        //ResetXErrors (g);
        myDraw (g, colorfulColors[iIter], kOpenCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);
      }

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    for (short iIter : {0, 1, 2, 3, 4, 5}) {
      myLineText2 (0.15, 0.66-0.07*iIter, colorfulColors[iIter], kOpenCircle, Form ("%i iterations", iIter+1), 1.2, 0.055);
    }

    c->SaveAs (Form ("%s/Plots/PtCh/UnfComp_Summary_JetSpectrum.pdf", workPath.Data ()));
  }




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_sigVsUnf_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGAE* g = nullptr;

      {
        c->cd (7);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];N_{iter} iterations / N_{iter}-1 iterations", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 4 : 2)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.88, 1.24);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);


        for (short iIter : {0, 1, 2, 3, 4, 5}) {
          h = h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter];
          g = make_graph (h);
          ScaleGraph (g, (iIter == 0 ? h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir] : h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter-1]));
          ResetXErrors (g);
          if (iDir == 1)
            TrimGraph (g, 0, 10);
          myDraw (g, colorfulColors[iIter], kOpenCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
        }

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];N_{iter} iterations / N_{iter}-1 iterations", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 4 : 2)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.88, 1.24);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        for (short iIter : {0, 1, 2, 3, 4, 5}) {
          h = h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter];
          g = make_graph (h);
          ScaleGraph (g, (iIter == 0 ? h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent] : h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter-1]));
          ResetXErrors (g);
          if (iDir == 1)
            TrimGraph (g, 0, 10);
          myDraw (g, colorfulColors[iIter], kOpenCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
        }

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{All centralities}", 0.05);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      for (short iIter : {0, 1, 2, 3, 4, 5}) {
        myLineText2 (0.15, 0.48-0.07*iIter, colorfulColors[iIter], kOpenCircle, Form ("%i iterations", iIter+1), 1.2, 0.055);
      }

      c->SaveAs (Form ("%s/Plots/PtCh/UnfComp_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt


}


#endif
