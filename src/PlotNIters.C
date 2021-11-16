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

#include "ProcessUnfolding.C" // gets NIter functions for plotting purposes


using namespace JetHadronCorrelations;


void PlotNIters (const char* inFileTag) {

  TLine* l = new TLine ();
  TLatex* tl = new TLatex ();

  TFile* inFile = nullptr;

  const short nItersMin = 1;
  const short nItersMax = 20;
  const double* nItersVals = linspace (nItersMin, nItersMax, nItersMax-nItersMin);

  const bool doLogY = true;
  const float ymaxSF = (doLogY ? 5 : 1.2);

  TH1D**    h_jet_pt_ref              = Get1DArray <TH1D*> (2);
  TH1D***   h_jet_pt                  = Get2DArray <TH1D*> (2, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_sig   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_sig       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_unf   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_unf       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_unf_nIters = Get3DArray <TH1D*> (nPtJBins, nDir, nItersMax-nItersMin+2);
  TH1D***** h_jetInt_trk_pt_unf_nIters     = Get4DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1, nItersMax-nItersMin+2);

  TH1D**    h_jet_pt_ref_unf_nIters       = Get1DArray <TH1D*> (nItersMax-nItersMin+2);
  TH1D***   h_jet_pt_unf_nIters           = Get2DArray <TH1D*> (nZdcCentBins+1, nItersMax-nItersMin+2);

  TGraph**    g_jet_pt_ref_unfIterUnc     = Get1DArray <TGraph*> (2);                 // sums of iterations uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfIterUnc         = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of iterations uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfStatUnc     = Get1DArray <TGraph*> (2);                 // sums of statistical uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfStatUnc         = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of statistical uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfTotUnc      = Get1DArray <TGraph*> (2);                 // sums of total uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfTotUnc          = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of total uncertainties as a function of nIter -- data only

  TGraph**    g_jet_pt_ref_unfIterRelUnc  = Get1DArray <TGraph*> (2);                 // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfIterRelUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfStatRelUnc  = Get1DArray <TGraph*> (2);                 // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfStatRelUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfTotRelUnc   = Get1DArray <TGraph*> (2);                 // sums of total relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfTotRelUnc       = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of total relative uncertainties as a function of nIter -- data only

  TGraph**    g_jet_pt_ref_unfIterHybUnc  = Get1DArray <TGraph*> (2);                 // sums of iterations hybrid uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfIterHybUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of iterations hybrid uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfStatHybUnc  = Get1DArray <TGraph*> (2);                 // sums of statistical hybrid uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfStatHybUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of statistical hybrid uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfTotHybUnc   = Get1DArray <TGraph*> (2);                 // sums of total hybrid uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfTotHybUnc       = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of total hybrid uncertainties as a function of nIter -- data only


  TGraph***   g_jetInt_trk_pt_ref_unfIterUnc    = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterUnc        = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfStatUnc    = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfStatUnc        = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotUnc     = Get2DArray <TGraph*> (2, nDir);                           // sums of total uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotUnc         = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total uncertainties as a function of nIter -- data only

  TGraph***   g_jetInt_trk_pt_ref_unfIterRelUnc = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterRelUnc     = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfStatRelUnc = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfStatRelUnc     = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotRelUnc  = Get2DArray <TGraph*> (2, nDir);                           // sums of total relative uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotRelUnc      = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total relative uncertainties as a function of nIter -- data only

  TGraph***   g_jetInt_trk_pt_ref_unfIterHybUnc = Get2DArray <TGraph*> (2, nDir);                           // sums of iterations hybrid uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfIterHybUnc     = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of iterations hybrid uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfStatHybUnc = Get2DArray <TGraph*> (2, nDir);                           // sums of statistical hybrid uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfStatHybUnc     = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of statistical hybrid uncertainties as a function of nIter -- data only
  TGraph***   g_jetInt_trk_pt_ref_unfTotHybUnc  = Get2DArray <TGraph*> (2, nDir);                           // sums of total hybrid uncertainties as a function of nIter -- data only
  TGraph****  g_jetInt_trk_pt_unfTotHybUnc      = Get3DArray <TGraph*> (2, nDir, nZdcCentBins+1);           // sums of total hybrid uncertainties as a function of nIter -- data only


  int***  jetInt_trk_pt_ref_niter_opt = Get3DArray <int> (2, nDir, 3);
  int**** jetInt_trk_pt_niter_opt     = Get4DArray <int> (2, nDir, nZdcCentBins+1, 3);


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

        const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
        const double maxJetPt = 300;

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_data_%s_Nominal_nIters%i", dir.Data (), pTJInt.Data (), nIters));
          h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter]->Scale (h_jet_pt_ref_unf_nIters[iIter]->Integral (h_jet_pt_ref_unf_nIters[iIter]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[iIter]->FindBin (maxJetPt-0.01)));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
  
            h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_data_%s_Nominal_nIters%i", dir.Data (), cent.Data (), pTJInt.Data (), nIters));
            h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter]->Scale (h_jet_pt_unf_nIters[iCent][iIter]->Integral (h_jet_pt_unf_nIters[iCent][iIter]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][iIter]->FindBin (maxJetPt-0.01)));
  
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
      const double maxJetPt = 300;

      TH1D* h_unf = h_jet_pt_ref[iDType];
      TH1D* h_unf_prev = h_jet_pt_ref[iDType];

      g_jet_pt_ref_unfStatUnc[iPtJInt]     = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfIterUnc[iPtJInt]     = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfTotUnc[iPtJInt]      = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfStatHybUnc[iPtJInt]  = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfIterHybUnc[iPtJInt]  = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfTotHybUnc[iPtJInt]   = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfStatRelUnc[iPtJInt]  = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfIterRelUnc[iPtJInt]  = new TGraph (nItersMax - nItersMin + 1);
      g_jet_pt_ref_unfTotRelUnc[iPtJInt]   = new TGraph (nItersMax - nItersMin + 1);


      {
        double tot = 0, totVar = 0, totHybVar = 0, totRelVar = 0;
        for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
          if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || maxJetPt < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
          tot += h_unf->GetBinContent (iX);
          totVar += std::pow (h_unf->GetBinError (iX), 2);
          totHybVar += std::pow (h_unf->GetBinError (iX), 2) / h_unf->GetBinContent (iX);
          totRelVar += std::pow (h_unf->GetBinError (iX) / h_unf->GetBinContent (iX), 2);
        } // end loop over iX

        g_jet_pt_ref_unfStatUnc[iPtJInt]->SetPoint    (g_jet_pt_ref_unfStatUnc[iPtJInt]->GetN (),     0, std::sqrt (totVar));
        g_jet_pt_ref_unfTotUnc[iPtJInt]->SetPoint     (g_jet_pt_ref_unfTotUnc[iPtJInt]->GetN (),      0, std::sqrt (totVar));
        g_jet_pt_ref_unfStatHybUnc[iPtJInt]->SetPoint (g_jet_pt_ref_unfStatHybUnc[iPtJInt]->GetN (),  0, std::sqrt (totHybVar));
        g_jet_pt_ref_unfTotHybUnc[iPtJInt]->SetPoint  (g_jet_pt_ref_unfTotHybUnc[iPtJInt]->GetN (),   0, std::sqrt (totHybVar));
        g_jet_pt_ref_unfStatRelUnc[iPtJInt]->SetPoint (g_jet_pt_ref_unfStatRelUnc[iPtJInt]->GetN (),  0, std::sqrt (totRelVar));
        g_jet_pt_ref_unfTotRelUnc[iPtJInt]->SetPoint  (g_jet_pt_ref_unfTotRelUnc[iPtJInt]->GetN (),   0, std::sqrt (totRelVar));
      }


      for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {

        const short nIters = (short) nItersVals[iIter];

        h_unf = h_jet_pt_ref_unf_nIters[iIter];

        double tot = 0, totVar = 0, totHybVar = 0, totRelVar = 0;
        for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
          if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || maxJetPt < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
          tot += h_unf->GetBinContent (iX);
          totVar += std::pow (h_unf->GetBinError (iX), 2);
          totHybVar += std::pow (h_unf->GetBinError (iX), 2) / h_unf->GetBinContent (iX);
          totRelVar += std::pow (h_unf->GetBinError (iX) / h_unf->GetBinContent (iX), 2);
        } // end loop over iX

        double iterVar = 0, iterHybVar = 0, iterRelVar = 0;
        for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
          if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || maxJetPt < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
          iterVar += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX));
          iterHybVar += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / std::sqrt (h_unf_prev->GetBinContent (iX));
          iterRelVar += std::fabs ((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / h_unf_prev->GetBinContent (iX));
        } // end loop over iX
        iterVar *= iterVar;
        iterHybVar *= iterHybVar;
        iterRelVar *= iterRelVar;

        g_jet_pt_ref_unfStatUnc[iPtJInt]->SetPoint    (nIters - nItersMin, nIters, std::sqrt (totVar));
        g_jet_pt_ref_unfIterUnc[iPtJInt]->SetPoint    (nIters - nItersMin, nIters, std::sqrt (iterVar));
        g_jet_pt_ref_unfTotUnc[iPtJInt]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totVar + iterVar));
        g_jet_pt_ref_unfStatHybUnc[iPtJInt]->SetPoint (nIters - nItersMin, nIters, std::sqrt (totHybVar));
        g_jet_pt_ref_unfIterHybUnc[iPtJInt]->SetPoint (nIters - nItersMin, nIters, std::sqrt (iterHybVar));
        g_jet_pt_ref_unfTotHybUnc[iPtJInt]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totHybVar + iterHybVar));
        g_jet_pt_ref_unfStatRelUnc[iPtJInt]->SetPoint (nIters - nItersMin, nIters, std::sqrt (totRelVar));
        g_jet_pt_ref_unfIterRelUnc[iPtJInt]->SetPoint (nIters - nItersMin, nIters, std::sqrt (iterRelVar));
        g_jet_pt_ref_unfTotRelUnc[iPtJInt]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelVar + iterRelVar));

        h_unf_prev = h_unf;
      } // end loop over iIter

    } // end loop over iPtJInt


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      for (short iPtJInt : {0, 1}) {

        const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
        const double maxJetPt = 300;

        TH1D* h_unf = h_jet_pt[iDType][iCent];
        TH1D* h_unf_prev = h_jet_pt[iDType][iCent];

        g_jet_pt_unfStatUnc[iPtJInt][iCent]     = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfIterUnc[iPtJInt][iCent]     = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfTotUnc[iPtJInt][iCent]      = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfStatHybUnc[iPtJInt][iCent]  = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfIterHybUnc[iPtJInt][iCent]  = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfTotHybUnc[iPtJInt][iCent]   = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfStatRelUnc[iPtJInt][iCent]  = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfIterRelUnc[iPtJInt][iCent]  = new TGraph (nItersMax - nItersMin + 1);
        g_jet_pt_unfTotRelUnc[iPtJInt][iCent]   = new TGraph (nItersMax - nItersMin + 1);


        {
          double tot = 0, totVar = 0, totHybVar = 0, totRelVar = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || maxJetPt < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
            tot += h_unf->GetBinContent (iX);
            totVar += std::pow (h_unf->GetBinError (iX), 2);
            totHybVar += std::pow (h_unf->GetBinError (iX), 2) / h_unf->GetBinContent (iX);
            totRelVar += std::pow (h_unf->GetBinError (iX) / h_unf->GetBinContent (iX), 2);
          } // end loop over iX

          g_jet_pt_unfStatUnc[iPtJInt][iCent]->SetPoint     (g_jet_pt_unfStatUnc[iPtJInt][iCent]->GetN (),    0, std::sqrt (totVar));
          g_jet_pt_unfTotUnc[iPtJInt][iCent]->SetPoint      (g_jet_pt_unfTotUnc[iPtJInt][iCent]->GetN (),     0, std::sqrt (totVar));
          g_jet_pt_unfStatHybUnc[iPtJInt][iCent]->SetPoint  (g_jet_pt_unfStatHybUnc[iPtJInt][iCent]->GetN (), 0, std::sqrt (totHybVar));
          g_jet_pt_unfTotHybUnc[iPtJInt][iCent]->SetPoint   (g_jet_pt_unfTotHybUnc[iPtJInt][iCent]->GetN (),  0, std::sqrt (totHybVar));
          g_jet_pt_unfStatRelUnc[iPtJInt][iCent]->SetPoint  (g_jet_pt_unfStatRelUnc[iPtJInt][iCent]->GetN (), 0, std::sqrt (totRelVar));
          g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->SetPoint   (g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->GetN (),  0, std::sqrt (totRelVar));
        }

 
        for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
          const short nIters = (short) nItersVals[iIter];

          h_unf = h_jet_pt_unf_nIters[iCent][iIter];

          double tot = 0, totVar = 0, totHybVar = 0, totRelVar = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || maxJetPt < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
            tot += h_unf->GetBinContent (iX);
            totVar += std::pow (h_unf->GetBinError (iX), 2);
            totHybVar += std::pow (h_unf->GetBinError (iX), 2) / h_unf->GetBinContent (iX);
            totRelVar += std::pow (h_unf->GetBinError (iX) / h_unf->GetBinContent (iX), 2);
          } // end loop over iX

          double iterVar = 0, iterHybVar = 0, iterRelVar = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            if (h_unf->GetXaxis ()->GetBinCenter (iX) < minJetPt || maxJetPt < h_unf->GetXaxis ()->GetBinCenter (iX)) continue;
            iterVar += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX));
            iterHybVar += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / std::sqrt (h_unf_prev->GetBinContent (iX));
            iterRelVar += std::fabs ((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / h_unf_prev->GetBinContent (iX));
          } // end loop over iX
          iterVar *= iterVar;
          iterHybVar *= iterHybVar;
          iterRelVar *= iterRelVar;

          g_jet_pt_unfStatUnc[iPtJInt][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totVar));
          g_jet_pt_unfIterUnc[iPtJInt][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (iterVar));
          g_jet_pt_unfTotUnc[iPtJInt][iCent]->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totVar + iterVar));
          g_jet_pt_unfStatHybUnc[iPtJInt][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totHybVar));
          g_jet_pt_unfIterHybUnc[iPtJInt][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (iterHybVar));
          g_jet_pt_unfTotHybUnc[iPtJInt][iCent]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totHybVar + iterHybVar));
          g_jet_pt_unfStatRelUnc[iPtJInt][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelVar));
          g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (iterRelVar));
          g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelVar + iterRelVar));

          h_unf_prev = h_unf;

        } // end loop over iIter

      } // end loop over iPtJInt

    } // end loop over iCent


    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      for (short iPtJInt : {0, 1}) {

        const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
        const double maxJetPt = 300;

        const double minPtCh = 4;
        const double maxPtCh = (iPtJInt == 0 ? 60 : 90);

        TH1D* h_unf = h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir];
        TH1D* h_unf_prev = h_unf;
        TH1D* h_j_unf = h_jet_pt_ref[iDType];
        TH1D* h_j_unf_prev = h_j_unf;

        g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]     = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]     = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]      = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]  = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]  = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir]   = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]  = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]  = new TGraph (nItersMax - nItersMin + 1);
        g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]   = new TGraph (nItersMax - nItersMin + 1);


        {
          const double nJets = 1;//h_j_unf->Integral (h_j_unf->FindBin (minJetPt+0.01), h_j_unf->FindBin (maxJetPt-0.01));

          double tot = 0, totVar = 0, totHybVar = 0, totRelVar = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
            if (ptch < minPtCh || maxPtCh < ptch) continue;
            const double err = h_unf->GetBinError (iX) * h_unf->GetBinWidth (iX) * nJets;
            const double wgt = h_unf->GetBinContent (iX) * h_unf->GetBinWidth (iX) * nJets;
            tot += wgt;
            totVar += std::pow (err, 2);
            totHybVar += std::pow (err, 2) / wgt;
            totRelVar += std::pow (err / wgt, 2);
          } // end loop over iX

          g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]->SetPoint     (g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]->GetN (),    0, std::sqrt (totVar));
          g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->SetPoint      (g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->GetN (),     0, std::sqrt (totVar));
          g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]->SetPoint  (g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]->GetN (), 0, std::sqrt (totHybVar));
          g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir]->SetPoint   (g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir]->GetN (),  0, std::sqrt (totHybVar));
          g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]->SetPoint  (g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]->GetN (), 0, std::sqrt (totRelVar));
          g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->SetPoint   (g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->GetN (),  0, std::sqrt (totRelVar));
        }


        for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {
  
          const short nIters = (short) nItersVals[iIter];

          h_unf = h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter];

          h_j_unf = h_jet_pt_ref_unf_nIters[iIter];
          const double nJets = 1;//h_j_unf->Integral (h_j_unf->FindBin (minJetPt+0.01), h_j_unf->FindBin (maxJetPt-0.01));
          const double nJets_prev = (iIter == 0 ? h_j_unf_prev->Integral (h_j_unf_prev->FindBin (minJetPt+0.01), h_j_unf_prev->FindBin (maxJetPt-0.01)) : 1);
          h_j_unf_prev = h_j_unf;

          double tot = 0, totVar = 0, totHybVar = 0, totRelVar = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
            if (ptch < minPtCh || maxPtCh < ptch) continue;
            const double err = h_unf->GetBinError (iX) * h_unf->GetBinWidth (iX) * nJets;
            const double wgt = h_unf->GetBinContent (iX) * h_unf->GetBinWidth (iX) * nJets;
            tot += wgt;
            totVar += std::pow (err, 2);
            totHybVar += std::pow (err, 2) / wgt;
            totRelVar += std::pow (err / wgt, 2);
          } // end loop over iX
  
          double iterVar = 0, iterHybVar = 0, iterRelVar = 0;
          for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
            const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
            if (ptch < minPtCh || maxPtCh < ptch) continue;
            const double diff = std::fabs ((h_unf->GetBinContent (iX) * nJets - h_unf_prev->GetBinContent (iX) * nJets_prev) * h_unf->GetBinWidth (iX));
            const double wgt = h_unf_prev->GetBinContent (iX) * h_unf_prev->GetBinWidth (iX) * nJets_prev;
            iterVar += diff;
            iterHybVar += diff / std::sqrt (wgt);
            iterRelVar += diff / wgt;
          } // end loop over iX
          iterVar *= iterVar;
          iterHybVar *= iterHybVar;
          iterRelVar *= iterRelVar;
  
          g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totVar));
          g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (iterVar));
          g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->SetPoint      (nIters - nItersMin, nIters, std::sqrt (totVar + iterVar));
          g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totHybVar));
          g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (iterHybVar));
          g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totHybVar + iterHybVar));
          g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelVar));
          g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (iterRelVar));
          g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->SetPoint   (nIters - nItersMin, nIters, std::sqrt (totRelVar + iterRelVar));

          h_unf_prev = h_unf;

        } // end loop over iIter

      } // end loop over iPtJInt

    } // end loop over iDir


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent)); 

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        for (short iPtJInt : {0, 1}) {

          const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
          const double maxJetPt = 300;

          const double minPtCh = 4;
          const double maxPtCh = (iPtJInt == 0 ? 60 : 90);
  
          TH1D* h_unf = h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent];
          TH1D* h_unf_prev = h_unf;
          TH1D* h_j_unf = h_jet_pt[iDType][iCent];
          TH1D* h_j_unf_prev = h_j_unf;

          g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]    = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]    = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]     = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent] = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent] = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent]  = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent] = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent] = new TGraph (nItersMax - nItersMin + 1);
          g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]  = new TGraph (nItersMax - nItersMin + 1);


          {
            const double nJets = 1;//h_j_unf->Integral (h_j_unf->FindBin (minJetPt+0.01), h_j_unf->FindBin (maxJetPt-0.01));

            double tot = 0, totVar = 0, totHybVar = 0, totRelVar = 0;
            for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
              const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
              if (ptch < minPtCh || maxPtCh < ptch) continue;
              const double err = h_unf->GetBinError (iX) * h_unf->GetBinWidth (iX) * nJets;
              const double wgt = h_unf->GetBinContent (iX) * h_unf->GetBinWidth (iX) * nJets;
              tot += wgt;
              totVar += std::pow (err, 2);
              totHybVar += std::pow (err, 2) / wgt;
              totRelVar += std::pow (err / wgt, 2);
            } // end loop over iX

            g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]->SetPoint    (g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]->GetN (),     0, std::sqrt (totVar));
            g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->SetPoint     (g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->GetN (),      0, std::sqrt (totVar));
            g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent]->SetPoint (g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent]->GetN (),  0, std::sqrt (totHybVar));
            g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent]->SetPoint  (g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent]->GetN (),   0, std::sqrt (totHybVar));
            g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent]->SetPoint (g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent]->GetN (),  0, std::sqrt (totRelVar));
            g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->SetPoint  (g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->GetN (),   0, std::sqrt (totRelVar));
          }


          for (short iIter = 0; iIter < nItersMax-nItersMin+2; iIter++) {

            const short nIters = (short) nItersVals[iIter];

            h_unf = h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter];

            h_j_unf = h_jet_pt_unf_nIters[iCent][iIter];
            const double nJets = 1;//h_j_unf->Integral (h_j_unf->FindBin (minJetPt+0.01), h_j_unf->FindBin (maxJetPt-0.01));
            const double nJets_prev = (iIter == 0 ? h_j_unf_prev->Integral (h_j_unf_prev->FindBin (minJetPt+0.01), h_j_unf_prev->FindBin (maxJetPt-0.01)) : 1);
            h_j_unf_prev = h_j_unf;

            double tot = 0, totVar = 0, totHybVar = 0, totRelVar = 0;
            for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
              const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
              if (ptch < minPtCh || maxPtCh < ptch) continue;
              const double err = h_unf->GetBinError (iX) * h_unf->GetBinWidth (iX) * nJets;
              const double wgt = h_unf->GetBinContent (iX) * h_unf->GetBinWidth (iX) * nJets;
              tot += wgt;
              totVar += std::pow (err, 2);
              totHybVar += std::pow (err, 2) / wgt;
              totRelVar += std::pow (err / wgt, 2);
            } // end loop over iX
    
            double iterVar = 0, iterHybVar = 0, iterRelVar = 0;
            for (short iX = 1; iX <= h_unf->GetNbinsX (); iX++) {
              const double ptch = h_unf->GetXaxis ()->GetBinCenter (iX);
              if (ptch < minPtCh || maxPtCh < ptch) continue;
              const double diff = std::fabs ((h_unf->GetBinContent (iX) * nJets - h_unf_prev->GetBinContent (iX) * nJets_prev) * h_unf->GetBinWidth (iX));
              const double wgt = h_unf_prev->GetBinContent (iX) * h_unf_prev->GetBinWidth (iX) * nJets_prev;
              iterVar += diff;
              iterHybVar += diff / std::sqrt (wgt);
              iterRelVar += diff / wgt;
            } // end loop over iX
            iterVar *= iterVar;
            iterHybVar *= iterHybVar;
            iterRelVar *= iterRelVar;

            g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]->SetPoint    (nIters - nItersMin, nIters, std::sqrt (totVar));
            g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->SetPoint    (nIters - nItersMin, nIters, std::sqrt (iterVar));
            g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->SetPoint     (nIters - nItersMin, nIters, std::sqrt (totVar + iterVar));
            g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent]->SetPoint (nIters - nItersMin, nIters, std::sqrt (totHybVar));
            g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent]->SetPoint (nIters - nItersMin, nIters, std::sqrt (iterHybVar));
            g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totHybVar + iterHybVar));
            g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent]->SetPoint (nIters - nItersMin, nIters, std::sqrt (totRelVar));
            g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->SetPoint (nIters - nItersMin, nIters, std::sqrt (iterRelVar));
            g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->SetPoint  (nIters - nItersMin, nIters, std::sqrt (totRelVar + iterRelVar));

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

        if (doLogY) gPad->SetLogy ();

        double x, y, xopt = -1;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= ymaxSF;

        g = g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        for (short i = 0; i < g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->GetN (); i++) {
          g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]->GetPoint (i, x, y);
          double ydum = y;
          g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->GetPoint (i, x, y);
          if ((xopt == -1 && ydum > y) || i == g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->GetN () - 1)
            xopt = x;
          else if (ydum < y)
            xopt = -1;
          i++;
        }

        if (doLogY) ymin *= 0.2;
        else        ymin = 0;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2}}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kFullCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenSquare);
        g->Draw ("P");

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (iPtJInt, iDir, -1);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][0] = xopt;
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        if (doLogY) gPad->SetLogy ();

        double x, y, xopt = -1;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= ymaxSF;

        g = g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        short i = 0;
        while (xopt == -1 && i < g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->GetN ()) {
          g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]->GetPoint (i, x, y);
          double ydum = y;
          g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->GetPoint (i, x, y);
          if (ydum > y)
            xopt = x;
          i++;
        }
        for (short i = 0; i < g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->GetN (); i++) {
          g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]->GetPoint (i, x, y);
          double ydum = y;
          g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->GetPoint (i, x, y);
          if ((xopt == -1 && ydum > y) || i == g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->GetN () - 1)
            xopt = x;
          else if (ydum < y)
            xopt = -1;
          i++;
        }

        if (doLogY) ymin *= 0.2;
        else        ymin = 0;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2}}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kFullCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenSquare);
        g->Draw ("P");

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (iPtJInt, iDir, iCent);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_niter_opt[iPtJInt][iDir][iCent][0] = xopt;

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (4, %i)", iPtJInt == 0 ? 60 : 90), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kFullCircle, "#sigma_{stat} #oplus #sigma_{iter}, #alpha = 0", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "#sigma_{stat} only, #alpha = 0", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenSquare, "#sigma_{iter} only, #alpha = 0", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_unfHybUncs_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGraph* g = nullptr;

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);

      double ymax = 0, ymin = 0;
      {
        c->cd (7);

        if (doLogY) gPad->SetLogy ();

        double x, y, xopt = -1;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= ymaxSF;

        g = g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        for (short i = 0; i < g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]->GetN (); i++) {
          g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]->GetPoint (i, x, y);
          double ydum = y;
          g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]->GetPoint (i, x, y);
          if ((xopt == -1 && ydum > y) || i == g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]->GetN () - 1)
            xopt = x;
          else if (ydum < y)
            xopt = -1;
          i++;
        }

        if (doLogY) ymin *= 0.2;
        else        ymin = 0;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2} / N}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kFullCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenSquare);
        g->Draw ("P");

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (iPtJInt, iDir, -1);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][1] = xopt;
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        if (doLogY) gPad->SetLogy ();

        double x, y, xopt = -1;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= ymaxSF;

        g = g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        for (short i = 0; i < g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent]->GetN (); i++) {
          g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent]->GetPoint (i, x, y);
          double ydum = y;
          g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent]->GetPoint (i, x, y);
          if ((xopt == -1 && ydum > y) || i == g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent]->GetN () - 1)
            xopt = x;
          else if (ydum < y)
            xopt = -1;
          i++;
        }

        if (doLogY) ymin *= 0.2;
        else        ymin = 0;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2} / N}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kFullCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenSquare);
        g->Draw ("P");

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (iPtJInt, iDir, iCent);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_niter_opt[iPtJInt][iDir][iCent][1] = xopt;

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (4, %i)", iPtJInt == 0 ? 60 : 90), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kFullCircle, "#sigma_{stat} #oplus #sigma_{iter}, #alpha = 0.5", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "#sigma_{stat} only, #alpha = 0.5", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenSquare, "#sigma_{iter} only, #alpha = 0.5", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfHybUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

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

        if (doLogY) gPad->SetLogy ();

        double x, y, xopt = -1;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= ymaxSF;

        g = g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        for (short i = 0; i < g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->GetN (); i++) {
          g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]->GetPoint (i, x, y);
          double ydum = y;
          g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->GetPoint (i, x, y);
          if ((xopt == -1 && ydum > y) || i == g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->GetN () - 1)
            xopt = x;
          else if (ydum < y)
            xopt = -1;
          i++;
        }

        if (doLogY) ymin *= 0.2;
        else        ymin = 0;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2} / N^{2}}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kFullCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenSquare);
        g->Draw ("P");

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (iPtJInt, iDir, -1);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][2] = xopt;
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        if (doLogY) gPad->SetLogy ();

        double x, y, xopt = -1;
        ymax = 0, ymin = DBL_MAX;
        g = g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (ymax > 0 && y > 2. * ymax) continue;
          ymax = std::fmax (y, ymax);
        }
        ymax *= ymaxSF;

        g = g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        g = g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent];
        for (short i = 0; i < g->GetN (); i++) {
          g->GetPoint (i, x, y);
          if (x != 0 && y < ymin && y > 0) ymin = y;
        }
        for (short i = 0; i < g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->GetN (); i++) {
          g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent]->GetPoint (i, x, y);
          double ydum = y;
          g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->GetPoint (i, x, y);
          if ((xopt == -1 && ydum > y) || i == g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->GetN () - 1)
            xopt = x;
          else if (ydum < y)
            xopt = -1;
          i++;
        }

        if (doLogY) ymin *= 0.2;
        else        ymin = 0;

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2} / N^{2}}", 1, 0, 20);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        g = g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlack);
        g->SetMarkerStyle (kFullCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kBlue);
        g->SetMarkerStyle (kOpenCircle);
        g->Draw ("P");

        g = g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent];
        g->SetMarkerColor (kRed);
        g->SetMarkerStyle (kOpenSquare);
        g->Draw ("P");

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (iPtJInt, iDir, iCent);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_niter_opt[iPtJInt][iDir][iCent][2] = xopt;

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (4, %i)", iPtJInt == 0 ? 60 : 90), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kFullCircle, "#sigma_{stat} #oplus #sigma_{iter}, #alpha = 1", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "#sigma_{stat} only, #alpha = 1", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenSquare, "#sigma_{iter} only, #alpha = 1", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfRelUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




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

      if (doLogY) gPad->SetLogy ();

      double x, y, xopt = -1;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_ref_unfTotUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= ymaxSF;

      g = g_jet_pt_ref_unfStatUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_ref_unfIterUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      for (short i = 0; i < g_jet_pt_ref_unfIterUnc[iPtJInt]->GetN (); i++) {
        g_jet_pt_ref_unfStatUnc[iPtJInt]->GetPoint (i, x, y);
        double ydum = y;
        g_jet_pt_ref_unfIterUnc[iPtJInt]->GetPoint (i, x, y);
        if ((xopt == -1 && ydum > y) || i == g_jet_pt_ref_unfIterUnc[iPtJInt]->GetN () - 1)
          xopt = x;
        else if (ydum < y)
          xopt = -1;
        i++;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_ref_unfTotUnc[iPtJInt];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kFullCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfStatUnc[iPtJInt];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfIterUnc[iPtJInt];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenSquare);
      g->Draw ("P");

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (iPtJInt, -1);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      if (doLogY) gPad->SetLogy ();

      double x, y, xopt = -1;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_unfTotUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= ymaxSF;

      g = g_jet_pt_unfStatUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_unfIterUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      for (short i = 0; i < g_jet_pt_unfIterUnc[iPtJInt][iCent]->GetN (); i++) {
        g_jet_pt_unfStatUnc[iPtJInt][iCent]->GetPoint (i, x, y);
        double ydum = y;
        g_jet_pt_unfIterUnc[iPtJInt][iCent]->GetPoint (i, x, y);
        if ((xopt == -1 && ydum > y) || i == g_jet_pt_unfIterUnc[iPtJInt][iCent]->GetN () - 1)
          xopt = x;
        else if (ydum < y)
          xopt = -1;
        i++;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_unfTotUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kFullCircle);
      g->Draw ("P");

      g = g_jet_pt_unfStatUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfIterUnc[iPtJInt][iCent];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenSquare);
      g->Draw ("P");

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (iPtJInt, iCent);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
    myLineText2 (0.15, 0.56, kBlack,        kFullCircle, "#sigma_{stat} #oplus #sigma_{iter}, #alpha = 0", 1.2, 0.06);
    myLineText2 (0.15, 0.48, kBlue,         kOpenCircle, "#sigma_{stat} only, #alpha = 0", 1.2, 0.06);
    myLineText2 (0.15, 0.40, kRed,          kOpenSquare, "#sigma_{iter} only, #alpha = 0", 1.2, 0.06);

    c->SaveAs (Form ("%s/Plots/PtCh/UnfUncs_Summary_JetSpectra_%s.pdf", workPath.Data (), pTJInt.Data ()));

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    const char* canvasName = Form ("c_jet_pt_unfHybUncs_%s", pTJInt.Data ());
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);

    TGraph* g = nullptr;

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);

    double ymax = 0, ymin = 0;
    {
      c->cd (7);

      if (doLogY) gPad->SetLogy ();

      double x, y, xopt = -1;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_ref_unfTotHybUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= ymaxSF;

      g = g_jet_pt_ref_unfStatHybUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_ref_unfIterHybUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      for (short i = 0; i < g_jet_pt_ref_unfIterHybUnc[iPtJInt]->GetN (); i++) {
        g_jet_pt_ref_unfStatHybUnc[iPtJInt]->GetPoint (i, x, y);
        double ydum = y;
        g_jet_pt_ref_unfIterHybUnc[iPtJInt]->GetPoint (i, x, y);
        if ((xopt == -1 && ydum > y) || i == g_jet_pt_ref_unfIterHybUnc[iPtJInt]->GetN () - 1)
          xopt = x;
        else if (ydum < y)
          xopt = -1;
        i++;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_ref_unfTotHybUnc[iPtJInt];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kFullCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfStatHybUnc[iPtJInt];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfIterHybUnc[iPtJInt];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenSquare);
      g->Draw ("P");

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (iPtJInt, -1);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      if (doLogY) gPad->SetLogy ();

      double x, y, xopt = -1;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_unfTotHybUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= ymaxSF;

      g = g_jet_pt_unfStatHybUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_unfIterHybUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      for (short i = 0; i < g_jet_pt_unfIterHybUnc[iPtJInt][iCent]->GetN (); i++) {
        g_jet_pt_unfStatHybUnc[iPtJInt][iCent]->GetPoint (i, x, y);
        double ydum = y;
        g_jet_pt_unfIterHybUnc[iPtJInt][iCent]->GetPoint (i, x, y);
        if ((xopt == -1 && ydum > y) || i == g_jet_pt_unfIterHybUnc[iPtJInt][iCent]->GetN () - 1)
          xopt = x;
        else if (ydum < y)
          xopt = -1;
        i++;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_unfTotHybUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kFullCircle);
      g->Draw ("P");

      g = g_jet_pt_unfStatHybUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfIterHybUnc[iPtJInt][iCent];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenSquare);
      g->Draw ("P");

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (iPtJInt, iCent);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
    myLineText2 (0.15, 0.56, kBlack,        kFullCircle, "#sigma_{stat} #oplus #sigma_{iter}, #alpha = 0.5", 1.2, 0.06);
    myLineText2 (0.15, 0.48, kBlue,         kOpenCircle, "#sigma_{stat} only, #alpha = 0.5", 1.2, 0.06);
    myLineText2 (0.15, 0.40, kRed,          kOpenSquare, "#sigma_{iter} only, #alpha = 0.5", 1.2, 0.06);

    c->SaveAs (Form ("%s/Plots/PtCh/UnfHybUncs_Summary_JetSpectra_%s.pdf", workPath.Data (), pTJInt.Data ()));

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

      if (doLogY) gPad->SetLogy ();

      double x, y, xopt = -1;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_ref_unfTotRelUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= ymaxSF;

      g = g_jet_pt_ref_unfStatRelUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_ref_unfIterRelUnc[iPtJInt];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      for (short i = 0; i < g_jet_pt_ref_unfIterRelUnc[iPtJInt]->GetN (); i++) {
        g_jet_pt_ref_unfStatRelUnc[iPtJInt]->GetPoint (i, x, y);
        double ydum = y;
        g_jet_pt_ref_unfIterRelUnc[iPtJInt]->GetPoint (i, x, y);
        if ((xopt == -1 && ydum > y) || i == g_jet_pt_ref_unfIterRelUnc[iPtJInt]->GetN () - 1)
          xopt = x;
        else if (ydum < y)
          xopt = -1;
        i++;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n^{2}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_ref_unfTotRelUnc[iPtJInt];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kFullCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfStatRelUnc[iPtJInt];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_ref_unfIterRelUnc[iPtJInt];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenSquare);
      g->Draw ("P");

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (iPtJInt, -1);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      if (doLogY) gPad->SetLogy ();

      double x, y, xopt = -1;
      ymax = 0, ymin = DBL_MAX;
      g = g_jet_pt_unfTotRelUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (ymax > 0 && y > 2. * ymax) continue;
        ymax = std::fmax (y, ymax);
      }
      ymax *= ymaxSF;

      g = g_jet_pt_unfStatRelUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      g = g_jet_pt_unfIterRelUnc[iPtJInt][iCent];
      for (short i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (x != 0 && y < ymin && y > 0) ymin = y;
      }
      for (short i = 0; i < g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->GetN (); i++) {
        g_jet_pt_unfStatRelUnc[iPtJInt][iCent]->GetPoint (i, x, y);
        double ydum = y;
        g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->GetPoint (i, x, y);
        if ((xopt == -1 && ydum > y) || i == g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->GetN () - 1)
          xopt = x;
        else if (ydum < y)
          xopt = -1;
        i++;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n^{2}}", 1, 0, 20);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      g = g_jet_pt_unfTotRelUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlack);
      g->SetMarkerStyle (kFullCircle);
      g->Draw ("P");

      g = g_jet_pt_unfStatRelUnc[iPtJInt][iCent];
      g->SetMarkerColor (kBlue);
      g->SetMarkerStyle (kOpenCircle);
      g->Draw ("P");

      g = g_jet_pt_unfIterRelUnc[iPtJInt][iCent];
      g->SetMarkerColor (kRed);
      g->SetMarkerStyle (kOpenSquare);
      g->Draw ("P");

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (iPtJInt, iCent);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
    myLineText2 (0.15, 0.56, kBlack,        kFullCircle, "#sigma_{stat} #oplus #sigma_{iter}, #alpha = 1", 1.2, 0.06);
    myLineText2 (0.15, 0.48, kBlue,         kOpenCircle, "#sigma_{stat} only, #alpha = 1", 1.2, 0.06);
    myLineText2 (0.15, 0.40, kRed,          kOpenSquare, "#sigma_{iter} only, #alpha = 1", 1.2, 0.06);

    c->SaveAs (Form ("%s/Plots/PtCh/UnfRelUncs_Summary_JetSpectra_%s.pdf", workPath.Data (), pTJInt.Data ()));

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


      short iCol = 7;
      for (short iIter : {11, 9, 7, 5, 3, 2, 1, 0}) {
        h = h_jet_pt_ref_unf_nIters[iIter];
        g = make_graph (h);
        ScaleGraph (g, (iIter == 0 ? h_jet_pt_ref[0] : h_jet_pt_ref_unf_nIters[iIter-1]));
        //ResetXErrors (g);
        myDraw (g, colorfulColors[iCol--], kOpenCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);
      }

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (30, 0.25, 30, 1.75);
      l->DrawLine (300, 0.25, 300, 1.75);

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

      short iCol = 7;
      for (short iIter : {11, 9, 7, 5, 3, 2, 1, 0}) {
        h = h_jet_pt_unf_nIters[iCent][iIter];
        g = make_graph (h);
        ScaleGraph (g, (iIter == 0 ? h_jet_pt[0][iCent] : h_jet_pt_unf_nIters[iCent][iIter-1]));
        //ResetXErrors (g);
        myDraw (g, colorfulColors[iCol--], kOpenCircle, 1.0, 1, 2, "P", false);
        SaferDelete (&g);
      }

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      l->SetLineColor (kBlack);
      l->DrawLine (30, 0.25, 30, 1.75);
      l->DrawLine (300, 0.25, 300, 1.75);

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
    short iCol = 7;
    for (short iIter : {11, 9, 7, 5}) {
      myLineText2 (0.55, 0.48-0.07*(iCol-4), colorfulColors[iCol], kOpenCircle, Form ("%i iterations", iIter+1), 1.2, 0.055); iCol--;
    }
    for (short iIter : {3, 2, 1, 0}) {
      myLineText2 (0.15, 0.48-0.07*(iCol), colorfulColors[iCol], kOpenCircle, Form ("%i iterations", iIter+1), 1.2, 0.055); iCol--;
    }

    c->SaveAs (Form ("%s/Plots/PtCh/UnfComp_Summary_JetSpectrum.pdf", workPath.Data ()));
  }




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
    const double maxJetPt = 300;

    for (short iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_sigVsUnf_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGAE* g = nullptr;

      {
        c->cd (7);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];N_{iter} iterations / N_{iter}-1 iterations", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.72, 1.44);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);


        short iCol = 7;
        for (short iIter : {11, 9, 7, 5, 3, 2, 1, 0}) {
          h = h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter];
          g = make_graph (h);
          ScaleGraph (g, (iIter == 0 ? h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir] : h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter-1]), iIter == 0 ? 1. / h_jet_pt_ref[0]->Integral (h_jet_pt_ref[0]->FindBin (minJetPt+0.01), h_jet_pt_ref[0]->FindBin (maxJetPt-0.01)) : 1);
          ResetXErrors (g);
          if (iDir == 1) TrimGraph (g, 0, 10);
          myDraw (g, colorfulColors[iCol--], kOpenCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
        }

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];N_{iter} iterations / N_{iter}-1 iterations", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.72, 1.44);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        short iCol = 7;
        for (short iIter : {11, 9, 7, 5, 3, 2, 1, 0}) {
          h = h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter];
          g = make_graph (h);
          ScaleGraph (g, (iIter == 0 ? h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent] : h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter-1]), iIter == 0 ? 1. / h_jet_pt[0][iCent]->Integral (h_jet_pt[0][iCent]->FindBin (minJetPt+0.01), h_jet_pt[0][iCent]->FindBin (maxJetPt-0.01)) : 1);
          ResetXErrors (g);
          if (iDir == 1) TrimGraph (g, 0, 10);
          myDraw (g, colorfulColors[iCol--], kOpenCircle, 1.0, 1, 2, "P", false);
          SaferDelete (&g);
        }

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      short iCol = 7;
      for (short iIter : {11, 9, 7, 5}) {
        myLineText2 (0.55, 0.48-0.07*(iCol-4), colorfulColors[iCol], kOpenCircle, Form ("%i iterations", iIter+1), 1.2, 0.055); iCol--;
      }
      for (short iIter : {3, 2, 1, 0}) {
        myLineText2 (0.15, 0.48-0.07*(iCol), colorfulColors[iCol], kOpenCircle, Form ("%i iterations", iIter+1), 1.2, 0.055); iCol--;
      }

      c->SaveAs (Form ("%s/Plots/PtCh/UnfComp_All_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    for (short iDir : {0, 1, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_sigVsUnf_Summary_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGAE* g = nullptr;

      {
        c->cd (7);

        gPad->SetLogx ();

        const short nIterOpt = GetTrkSpectraNIters (iPtJInt, iDir, -1);//jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][iPtJInt == 0 ? 1 : 0];
        const short nIterVar = nIterOpt + 1;//jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][0];

        TH1D* h = new TH1D ("h", Form (";#it{p}_{T}^{ch} [GeV];%i vs. %i iterations", nIterOpt, nIterVar), 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.97, 1.03);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);


        h = h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][nIterOpt];
        //h->Add (h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][nIterVar], -1);
        g = make_graph (h);
        ResetTGAEErrors (g);
        ResetXErrors (g);
        CalcSystematics (g, h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][nIterVar], false);
        ScaleGraph (g, h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][nIterOpt]);
        if (iDir == 1) TrimGraph (g, 0, 10);

        TGAE* gup = new TGAE ();
        TGAE* gdown = new TGAE ();
        MakeGupAndGdown (g, gup, gdown);
        myDraw (gup, colorfulColors[0], kDot, 0, 1, 2, "L");
        myDraw (gdown, colorfulColors[0], kDot, 0, 1, 2, "L");
        myDrawFill (gup, gdown, colorfulColors[0], 0.3);
        SaferDelete (&g);
        SaferDelete (&gup);
        SaferDelete (&gdown);

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        gPad->SetLogx ();

        const short nIterOpt = GetTrkSpectraNIters (iPtJInt, iDir, iCent);//jetInt_trk_pt_niter_opt[iPtJInt][iDir][iCent][iPtJInt == 0 ? 1 : 0];
        const short nIterVar = nIterOpt + 1;//jetInt_trk_pt_niter_opt[iPtJInt][iDir][iCent][0];

        TH1D* h = new TH1D ("h", Form (";#it{p}_{T}^{ch} [GeV];%i vs. %i iterations", nIterOpt, nIterVar), 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.97, 1.03);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        // for Niter / Niter-1 ratio plot
        h = h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][nIterOpt];
        g = make_graph (h);
        ResetTGAEErrors (g);
        ResetXErrors (g);
        CalcSystematics (g, h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][nIterVar], false);
        ScaleGraph (g, h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][nIterOpt]);
        if (iDir == 1) TrimGraph (g, 0, 10);

        TGAE* gup = new TGAE ();
        TGAE* gdown = new TGAE ();
        MakeGupAndGdown (g, gup, gdown);
        myDraw (gup, colorfulColors[iCent+1], kDot, 0, 1, 2, "L");
        myDraw (gdown, colorfulColors[iCent+1], kDot, 0, 1, 2, "L");
        myDrawFill (gup, gdown, colorfulColors[iCent+1], 0.3);
        SaferDelete (&g);
        SaferDelete (&gup);
        SaferDelete (&gdown);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);

      } // end loop over iCent

      c->cd (8);
      myText (0.1, 0.93, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
      myText (0.1, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.75, kBlack, "#it{p}+Pb, #sqrt{s} = 5.02 TeV", 0.07);
      myText (0.1, 0.66, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, 300)", iPtJInt == 0 ? 30 : 60), 0.065);
      myText (0.1, 0.57, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfComp_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt


}


#endif
