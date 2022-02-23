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


double GetMinPtCh (short iPtJInt) {
  switch (iPtJInt) {
  case 0: return 4;//8;
  case 1: return 4;//8;
  }
  return 0.4;
}
double GetMaxPtCh (short iPtJInt) {
  switch (iPtJInt) {
  case 0: return 60;//30;
  case 1: return 90;//60;
  }
  return 120;
}



double GetTGraphMin (const TGraph* g, double ymin = DBL_MAX) {
  double x = 0, y = 0;
  for (short i = 0; i < g->GetN (); i++) {
    g->GetPoint (i, x, y);
    if (x != 0 && y > 0) ymin = std::fmin (y, ymin);
  }
  return ymin;
}



double GetTGraphMax (const TGraph* g, double ymax = DBL_MIN) {
  double x = 0, y = 0;
  for (short i = 0; i < g->GetN (); i++) {
    g->GetPoint (i, x, y);
    if (x != 0 && y > 0) ymax = std::fmax (y, ymax);
  }
  return ymax;
}



int GetCrossOverPoint (const TGraph* gstat, const TGraph* giter, const int offset = 1) {
  double x = 0, y = 0, xopt = 0;
  for (short i = 0; i < giter->GetN (); i++) {
    gstat->GetPoint (i+1, x, y);
    double ydum = y;
    giter->GetPoint (i, x, y);
    if (xopt == -1 && (ydum > y || i == giter->GetN () - 1))
      xopt = x;
    else if (ydum < y)
      xopt = -1;
  }
  return xopt;
}



void PlotNIters (const char* rawTag, const char* nitersTag, const short nItersMax = 20, const short nIters1DMax = 100) {

  TLine* l = new TLine ();
  TLatex* tl = new TLatex ();

  TFile* inFile = nullptr;

  const short nItersMin = 1;
  const double* nItersVals = linspace (nItersMin, nItersMax, nItersMax-nItersMin);
  const short nIters1DMin = 1;
  const double* nIters1DVals = linspace (nIters1DMin, nIters1DMax, nIters1DMax-nIters1DMin);

  const bool doLogY = true;
  const float ymaxSF = (doLogY ? 5 : 1.2);

  TH1D**    h_jet_pt_ref              = Get1DArray <TH1D*> (2);
  TH1D***   h_jet_pt                  = Get2DArray <TH1D*> (2, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_sig   = Get3DArray <TH1D*> (2, 2, nDir);
  TH1D***** h_jetInt_trk_pt_sig       = Get4DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D****  h_jetInt_trk_pt_ref_unf_nIters  = Get3DArray <TH1D*> (nPtJBins, nDir, nItersMax-nItersMin+1);
  TH1D***** h_jetInt_trk_pt_unf_nIters      = Get4DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1, nItersMax-nItersMin+1);

  TH1D****  h_jetInt_trk_pt_ref_rfld_nIters = Get3DArray <TH1D*> (nPtJBins, nDir, nItersMax-nItersMin+1);
  TH1D***** h_jetInt_trk_pt_rfld_nIters     = Get4DArray <TH1D*> (nPtJBins, nDir, nZdcCentBins+1, nItersMax-nItersMin+1);

  TH1D**    h_jet_pt_ref_unf_nIters       = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+1);
  TH1D**    h_jet_pt_ref_rfld_nIters      = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+1);
  TH1D***   h_jet_pt_unf_nIters           = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+1);
  TH1D***   h_jet_pt_rfld_nIters          = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+1);

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
    TString inFileName = rawTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    //for (short iDType = 0; iDType < 2; iDType++) {
    for (short iDType = 0; iDType < 1; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      h_jet_pt_ref[iDType] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_%s_Nominal", dType.Data ()));

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        h_jet_pt[iDType][iCent] = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_%s_Nominal", cent.Data (), dType.Data ()));

      } // end loop over iCent

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));

        } // end loop over iDir


        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          //for (short iDir = 0; iDir < nDir; iDir++) {
          for (short iDir : {0, 2}) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

          } // end loop over iDir

        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType

  }



  {
    TString inFileName = nitersTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessNIters_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    //for (short iDType = 0; iDType < 2; iDType++) {
    for (short iDType = 0; iDType < 1; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir]  = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_Nominal",  dir.Data (), dType.Data (), pTJInt.Data ()));

        } // end loop over iDir


        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          //for (short iDir = 0; iDir < nDir; iDir++) {
          for (short iDir : {0, 2}) {

            const TString dir = directions[iDir];

            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent]       = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_Nominal", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()));

          } // end loop over iDir

        } // end loop over iCent

      } // end loop over iPtJInt

    } // end loop over iDType


    for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
  
      const short nIters = (short) nIters1DVals[iIter];

      h_jet_pt_ref_unf_nIters[iIter] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_unf_data_Nominal_nIters%i", nIters));
      h_jet_pt_ref_rfld_nIters[iIter] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_rfld_data_Nominal_nIters%i", nIters));

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        h_jet_pt_unf_nIters[iCent][iIter] = (TH1D*) inFile-> Get (Form ("h_jet_pt_unf_data_%s_Nominal_nIters%i", cent.Data (), nIters));
        h_jet_pt_rfld_nIters[iCent][iIter] = (TH1D*) inFile-> Get (Form ("h_jet_pt_rfld_data_%s_Nominal_nIters%i", cent.Data (), nIters));

      } // end loop over iCent

    } // end loop over iIter


    for (short iIter = 0; iIter < nItersMax-nItersMin+1; iIter++) {
  
      const short nIters = (short) nItersVals[iIter];

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
        const double maxJetPt = 300;

        //for (short iDir = 0; iDir < nDir; iDir++) {
        for (short iDir : {0, 2}) {

          const TString dir = directions[iDir];

          h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_unf_data_%s_Nominal_nIters%i", dir.Data (), pTJInt.Data (), nIters));
          h_jetInt_trk_pt_ref_rfld_nIters[iPtJInt][iDir][iIter] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_ref_rfld_data_%s_Nominal_nIters%i", dir.Data (), pTJInt.Data (), nIters));
          //h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter]->Scale (h_jet_pt_ref_unf_nIters[iIter]->Integral (h_jet_pt_ref_unf_nIters[iIter]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[iIter]->FindBin (maxJetPt-0.01)));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
  
            h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_unf_%s_data_%s_Nominal_nIters%i", cent.Data (), dir.Data (), pTJInt.Data (), nIters));
            h_jetInt_trk_pt_rfld_nIters[iPtJInt][iDir][iCent][iIter] = (TH1D*) inFile->Get (Form ("h_jetInt_trk_pt_%s_pPb_rfld_%s_data_%s_Nominal_nIters%i", cent.Data (), dir.Data (), pTJInt.Data (), nIters));
            //h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter]->Scale (h_jet_pt_unf_nIters[iCent][iIter]->Integral (h_jet_pt_unf_nIters[iCent][iIter]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][iIter]->FindBin (maxJetPt-0.01)));
  
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

      g_jet_pt_ref_unfStatUnc[iPtJInt]     = new TGraph ();
      g_jet_pt_ref_unfIterUnc[iPtJInt]     = new TGraph ();
      g_jet_pt_ref_unfTotUnc[iPtJInt]      = new TGraph ();
      g_jet_pt_ref_unfStatHybUnc[iPtJInt]  = new TGraph ();
      g_jet_pt_ref_unfIterHybUnc[iPtJInt]  = new TGraph ();
      g_jet_pt_ref_unfTotHybUnc[iPtJInt]   = new TGraph ();
      g_jet_pt_ref_unfStatRelUnc[iPtJInt]  = new TGraph ();
      g_jet_pt_ref_unfIterRelUnc[iPtJInt]  = new TGraph ();
      g_jet_pt_ref_unfTotRelUnc[iPtJInt]   = new TGraph ();


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


      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {

        const short nIters = (short) nIters1DVals[iIter];

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
          iterVar += std::pow (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX), 2);
          iterHybVar += std::pow ((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / std::sqrt (h_unf_prev->GetBinContent (iX)), 2);
          iterRelVar += std::pow (((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / h_unf_prev->GetBinContent (iX)), 2);
        } // end loop over iX

        g_jet_pt_ref_unfStatUnc[iPtJInt]->SetPoint    (g_jet_pt_ref_unfStatUnc[iPtJInt]->GetN (),     nIters, std::sqrt (totVar));
        g_jet_pt_ref_unfIterUnc[iPtJInt]->SetPoint    (g_jet_pt_ref_unfIterUnc[iPtJInt]->GetN (),     nIters, std::sqrt (iterVar));
        g_jet_pt_ref_unfTotUnc[iPtJInt]->SetPoint     (g_jet_pt_ref_unfTotUnc[iPtJInt]->GetN (),      nIters, std::sqrt (totVar + iterVar));
        g_jet_pt_ref_unfStatHybUnc[iPtJInt]->SetPoint (g_jet_pt_ref_unfStatHybUnc[iPtJInt]->GetN (),  nIters, std::sqrt (totHybVar));
        g_jet_pt_ref_unfIterHybUnc[iPtJInt]->SetPoint (g_jet_pt_ref_unfIterHybUnc[iPtJInt]->GetN (),  nIters, std::sqrt (iterHybVar));
        g_jet_pt_ref_unfTotHybUnc[iPtJInt]->SetPoint  (g_jet_pt_ref_unfTotHybUnc[iPtJInt]->GetN (),   nIters, std::sqrt (totHybVar + iterHybVar));
        g_jet_pt_ref_unfStatRelUnc[iPtJInt]->SetPoint (g_jet_pt_ref_unfStatRelUnc[iPtJInt]->GetN (),  nIters, std::sqrt (totRelVar));
        g_jet_pt_ref_unfIterRelUnc[iPtJInt]->SetPoint (g_jet_pt_ref_unfIterRelUnc[iPtJInt]->GetN (),  nIters, std::sqrt (iterRelVar));
        g_jet_pt_ref_unfTotRelUnc[iPtJInt]->SetPoint  (g_jet_pt_ref_unfTotRelUnc[iPtJInt]->GetN (),   nIters, std::sqrt (totRelVar + iterRelVar));

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

        g_jet_pt_unfStatUnc[iPtJInt][iCent]     = new TGraph ();
        g_jet_pt_unfIterUnc[iPtJInt][iCent]     = new TGraph ();
        g_jet_pt_unfTotUnc[iPtJInt][iCent]      = new TGraph ();
        g_jet_pt_unfStatHybUnc[iPtJInt][iCent]  = new TGraph ();
        g_jet_pt_unfIterHybUnc[iPtJInt][iCent]  = new TGraph ();
        g_jet_pt_unfTotHybUnc[iPtJInt][iCent]   = new TGraph ();
        g_jet_pt_unfStatRelUnc[iPtJInt][iCent]  = new TGraph ();
        g_jet_pt_unfIterRelUnc[iPtJInt][iCent]  = new TGraph ();
        g_jet_pt_unfTotRelUnc[iPtJInt][iCent]   = new TGraph ();


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

 
        for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
  
          const short nIters = (short) nIters1DVals[iIter];

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
            iterVar += std::pow (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX), 2);
            iterHybVar += std::pow ((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / std::sqrt (h_unf_prev->GetBinContent (iX)), 2);
            iterRelVar += std::pow (((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / h_unf_prev->GetBinContent (iX)), 2);
          } // end loop over iX

          g_jet_pt_unfStatUnc[iPtJInt][iCent]->SetPoint     (g_jet_pt_unfStatUnc[iPtJInt][iCent]->GetN (),    nIters, std::sqrt (totVar));
          g_jet_pt_unfIterUnc[iPtJInt][iCent]->SetPoint     (g_jet_pt_unfIterUnc[iPtJInt][iCent]->GetN (),    nIters, std::sqrt (iterVar));
          g_jet_pt_unfTotUnc[iPtJInt][iCent]->SetPoint      (g_jet_pt_unfTotUnc[iPtJInt][iCent]->GetN (),     nIters, std::sqrt (totVar + iterVar));
          g_jet_pt_unfStatHybUnc[iPtJInt][iCent]->SetPoint  (g_jet_pt_unfStatHybUnc[iPtJInt][iCent]->GetN (), nIters, std::sqrt (totHybVar));
          g_jet_pt_unfIterHybUnc[iPtJInt][iCent]->SetPoint  (g_jet_pt_unfIterHybUnc[iPtJInt][iCent]->GetN (), nIters, std::sqrt (iterHybVar));
          g_jet_pt_unfTotHybUnc[iPtJInt][iCent]->SetPoint   (g_jet_pt_unfTotHybUnc[iPtJInt][iCent]->GetN (),  nIters, std::sqrt (totHybVar + iterHybVar));
          g_jet_pt_unfStatRelUnc[iPtJInt][iCent]->SetPoint  (g_jet_pt_unfStatRelUnc[iPtJInt][iCent]->GetN (), nIters, std::sqrt (totRelVar));
          g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->SetPoint  (g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->GetN (), nIters, std::sqrt (iterRelVar));
          g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->SetPoint   (g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->GetN (),  nIters, std::sqrt (totRelVar + iterRelVar));

          h_unf_prev = h_unf;

        } // end loop over iIter

      } // end loop over iPtJInt

    } // end loop over iCent


    //for (short iDir = 0; iDir < nDir; iDir++) {
    for (short iDir : {0, 2}) {

      const TString dir = directions[iDir];

      for (short iPtJInt : {0, 1}) {

        const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
        const double maxJetPt = 300;

        const double minPtCh = GetMinPtCh (iPtJInt);
        const double maxPtCh = GetMaxPtCh (iPtJInt);

        TH1D* h_unf = h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir];
        TH1D* h_unf_prev = h_unf;
        TH1D* h_j_unf = h_jet_pt_ref[iDType];
        TH1D* h_j_unf_prev = h_j_unf;

        g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]     = new TGraph ();
        g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]     = new TGraph ();
        g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]      = new TGraph ();
        g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]  = new TGraph ();
        g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]  = new TGraph ();
        g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir]   = new TGraph ();
        g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]  = new TGraph ();
        g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]  = new TGraph ();
        g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]   = new TGraph ();


        {
          const double nJets = h_j_unf->Integral (h_j_unf->FindBin (minJetPt+0.01), h_j_unf->FindBin (maxJetPt-0.01));

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


        for (short iIter = 0; iIter < nItersMax-nItersMin+1; iIter++) {
  
          const short nIters = (short) nItersVals[iIter];

          h_unf = h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter];

          h_j_unf = h_jet_pt_ref_unf_nIters[iIter];
          const double nJets = h_j_unf->Integral (h_j_unf->FindBin (minJetPt+0.01), h_j_unf->FindBin (maxJetPt-0.01));
          //const double nJets_prev = (iIter == 0 ? h_j_unf_prev->Integral (h_j_unf_prev->FindBin (minJetPt+0.01), h_j_unf_prev->FindBin (maxJetPt-0.01)) : 1);
          const double nJets_prev = h_j_unf_prev->Integral (h_j_unf_prev->FindBin (minJetPt+0.01), h_j_unf_prev->FindBin (maxJetPt-0.01));
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
            iterVar += std::pow (diff, 2);
            iterHybVar += std::pow (diff, 2) / wgt;
            iterRelVar += std::pow (diff / wgt, 2);
          } // end loop over iX
  
          g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]->SetPoint     (g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]->GetN (),    nIters, std::sqrt (totVar));
          g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->SetPoint     (g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]->GetN (),    nIters, std::sqrt (iterVar));
          g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->SetPoint      (g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir]->GetN (),     nIters, std::sqrt (totVar + iterVar));
          g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]->SetPoint  (g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]->GetN (), nIters, std::sqrt (totHybVar));
          g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]->SetPoint  (g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]->GetN (), nIters, std::sqrt (iterHybVar));
          g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir]->SetPoint   (g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir]->GetN (),  nIters, std::sqrt (totHybVar + iterHybVar));
          g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]->SetPoint  (g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]->GetN (), nIters, std::sqrt (totRelVar));
          g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->SetPoint  (g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]->GetN (), nIters, std::sqrt (iterRelVar));
          g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->SetPoint   (g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir]->GetN (),  nIters, std::sqrt (totRelVar + iterRelVar));

          h_unf_prev = h_unf;

        } // end loop over iIter

      } // end loop over iPtJInt

    } // end loop over iDir


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent)); 

      //for (short iDir = 0; iDir < nDir; iDir++) {
      for (short iDir : {0, 2}) {

        const TString dir = directions[iDir];

        for (short iPtJInt : {0, 1}) {

          const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
          const double maxJetPt = 300;

          const double minPtCh = GetMinPtCh (iPtJInt);
          const double maxPtCh = GetMaxPtCh (iPtJInt);
  
          TH1D* h_unf = h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent];
          TH1D* h_unf_prev = h_unf;
          TH1D* h_j_unf = h_jet_pt[iDType][iCent];
          TH1D* h_j_unf_prev = h_j_unf;

          g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]    = new TGraph ();
          g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]    = new TGraph ();
          g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]     = new TGraph ();
          g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent] = new TGraph ();
          g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent] = new TGraph ();
          g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent]  = new TGraph ();
          g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent] = new TGraph ();
          g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent] = new TGraph ();
          g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]  = new TGraph ();


          {
            const double nJets = h_j_unf->Integral (h_j_unf->FindBin (minJetPt+0.01), h_j_unf->FindBin (maxJetPt-0.01));

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


          for (short iIter = 0; iIter < nItersMax-nItersMin+1; iIter++) {

            const short nIters = (short) nItersVals[iIter];

            h_unf = h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter];

            h_j_unf = h_jet_pt_unf_nIters[iCent][iIter];
            const double nJets = h_j_unf->Integral (h_j_unf->FindBin (minJetPt+0.01), h_j_unf->FindBin (maxJetPt-0.01));
            //const double nJets_prev = (iIter == 0 ? h_j_unf_prev->Integral (h_j_unf_prev->FindBin (minJetPt+0.01), h_j_unf_prev->FindBin (maxJetPt-0.01)) : 1);
            const double nJets_prev = h_j_unf_prev->Integral (h_j_unf_prev->FindBin (minJetPt+0.01), h_j_unf_prev->FindBin (maxJetPt-0.01));
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
              iterVar += std::pow (diff, 2);
              iterHybVar += std::pow (diff, 2) / wgt;
              iterRelVar += std::pow (diff / wgt, 2);
            } // end loop over iX

            g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]->SetPoint    (g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]->GetN (),     nIters, std::sqrt (totVar));
            g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->SetPoint    (g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]->GetN (),     nIters, std::sqrt (iterVar));
            g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->SetPoint     (g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent]->GetN (),      nIters, std::sqrt (totVar + iterVar));
            g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent]->SetPoint (g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent]->GetN (),  nIters, std::sqrt (totHybVar));
            g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent]->SetPoint (g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent]->GetN (),  nIters, std::sqrt (iterHybVar));
            g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent]->SetPoint  (g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent]->GetN (),   nIters, std::sqrt (totHybVar + iterHybVar));
            g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent]->SetPoint (g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent]->GetN (),  nIters, std::sqrt (totRelVar));
            g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->SetPoint (g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]->GetN (),  nIters, std::sqrt (iterRelVar));
            g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->SetPoint  (g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent]->GetN (),   nIters, std::sqrt (totRelVar + iterRelVar));

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

    //for (short iDir : {0, 1, 2}) {
    for (short iDir : {0, 2}) {

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

        // Find max y value
        ymax = ymaxSF * GetTGraphMax (g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir], 0);

        // Find min y value
        if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir]), GetTGraphMin (g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir]));
        else        ymin = 0;

        // Calculate cross-over point to optimize Niter
        const double xopt = GetCrossOverPoint (g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir], g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir], 1);

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2}}", 1, 0, nItersMax);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        myDraw (g_jetInt_trk_pt_ref_unfTotUnc[iPtJInt][iDir], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_ref_unfStatUnc[iPtJInt][iDir], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_ref_unfIterUnc[iPtJInt][iDir], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (false, iPtJInt, iDir, -1);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][0] = xopt;
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        if (doLogY) gPad->SetLogy ();

        // Find max y value
        ymax = ymaxSF * GetTGraphMax (g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent], 0);

        // Find min y value
        if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent]), GetTGraphMin (g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent]));
        else        ymin = 0;

        // Calculate cross-over point to optimize Niter
        const double xopt = GetCrossOverPoint (g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent], g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent], 1);

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2}}", 1, 0, nItersMax);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        myDraw (g_jetInt_trk_pt_unfTotUnc[iPtJInt][iDir][iCent], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_unfStatUnc[iPtJInt][iDir][iCent], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_unfIterUnc[iPtJInt][iDir][iCent], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (false, iPtJInt, iDir, iCent);
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
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (%g, %g)", GetMinPtCh (iPtJInt), GetMaxPtCh (iPtJInt)), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kFullCircle, "#sigma_{stat} #oplus #sigma_{iter}, #alpha = 0", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "#sigma_{stat} only, #alpha = 0", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenSquare, "#sigma_{iter} only, #alpha = 0", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    //for (short iDir : {0, 1, 2}) {
    for (short iDir : {0, 2}) {

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

        // Find max y value
        ymax = ymaxSF * GetTGraphMax (g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir], 0);

        // Find min y value
        if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir]), GetTGraphMin (g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir]));
        else        ymin = 0;

        // Calculate cross-over point to optimize Niter
        const double xopt = GetCrossOverPoint (g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir], g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir], 1);

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2} / N}", 1, 0, nItersMax);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        myDraw (g_jetInt_trk_pt_ref_unfTotHybUnc[iPtJInt][iDir], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_ref_unfStatHybUnc[iPtJInt][iDir], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_ref_unfIterHybUnc[iPtJInt][iDir], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (false, iPtJInt, iDir, -1);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][1] = xopt;
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        if (doLogY) gPad->SetLogy ();

        // Find max y value
        ymax = ymaxSF * GetTGraphMax (g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent], 0);

        // Find min y value
        if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent]), GetTGraphMin (g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent]));
        else        ymin = 0;

        // Calculate cross-over point to optimize Niter
        const double xopt = GetCrossOverPoint (g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent], g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent], 1);

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2} / N}", 1, 0, nItersMax);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        myDraw (g_jetInt_trk_pt_unfTotHybUnc[iPtJInt][iDir][iCent], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_unfStatHybUnc[iPtJInt][iDir][iCent], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_unfIterHybUnc[iPtJInt][iDir][iCent], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (false, iPtJInt, iDir, iCent);
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
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (%g, %g)", GetMinPtCh (iPtJInt), GetMaxPtCh (iPtJInt)), 0.065);
      myText (0.1, 0.48, kBlack, Form ("#Delta#phi_{ch,jet} %s", directions[iDir] == "ns" ? "< #pi/8" : (directions[iDir] == "as" ? "> 7#pi/8" : "#in (#pi/3, 2#pi/3)")), 0.065);
      myLineText2 (0.15, 0.40, kBlack,        kFullCircle, "#sigma_{stat} #oplus #sigma_{iter}, #alpha = 0.5", 1.2, 0.06);
      myLineText2 (0.15, 0.31, kBlue,         kOpenCircle, "#sigma_{stat} only, #alpha = 0.5", 1.2, 0.06);
      myLineText2 (0.15, 0.22, kRed,          kOpenSquare, "#sigma_{iter} only, #alpha = 0.5", 1.2, 0.06);

      c->SaveAs (Form ("%s/Plots/PtCh/UnfHybUncs_Summary_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    //for (short iDir : {0, 1, 2}) {
    for (short iDir : {0, 2}) {

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

        // Find max y value
        ymax = ymaxSF * GetTGraphMax (g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir], 0);

        // Find min y value
        if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir]), GetTGraphMin (g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir]));
        else        ymin = 0;

        // Calculate cross-over point to optimize Niter
        const double xopt = GetCrossOverPoint (g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir], g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir], 1);

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2} / N^{2}}", 1, 0, nItersMax);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        myDraw (g_jetInt_trk_pt_ref_unfTotRelUnc[iPtJInt][iDir], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_ref_unfStatRelUnc[iPtJInt][iDir], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_ref_unfIterRelUnc[iPtJInt][iDir], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

        myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (false, iPtJInt, iDir, -1);
        l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
        l->SetLineColor (kBlack);
        l->SetLineStyle (2);

        jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][2] = xopt;
      }


      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
        c->cd (nZdcCentBins+1-iCent);

        if (doLogY) gPad->SetLogy ();

        // Find max y value
        ymax = ymaxSF * GetTGraphMax (g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent], 0);

        // Find min y value
        if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent]), GetTGraphMin (g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent]));
        else        ymin = 0;

        // Calculate cross-over point to optimize Niter
        const double xopt = GetCrossOverPoint (g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent], g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent], 1);

        TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{N}^{2} / N^{2}}", 1, 0, nItersMax);
        h->GetYaxis ()->SetRangeUser (ymin, ymax);
        h->SetLineWidth (0);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        myDraw (g_jetInt_trk_pt_unfTotRelUnc[iPtJInt][iDir][iCent], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_unfStatRelUnc[iPtJInt][iDir][iCent], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
        myDraw (g_jetInt_trk_pt_unfIterRelUnc[iPtJInt][iDir][iCent], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

        if (iCent < nZdcCentBins)
          myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
        else
          myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
        l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

        l->SetLineColor (myGreen);
        l->SetLineStyle (3);
        const int xopt_2d = GetTrkSpectraNIters (false, iPtJInt, iDir, iCent);
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
      myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{ch} [GeV] #in (%g, %g)", GetMinPtCh (iPtJInt), GetMaxPtCh (iPtJInt)), 0.065);
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

      // Find max y value
      ymax = ymaxSF * GetTGraphMax (g_jet_pt_ref_unfTotUnc[iPtJInt], 0);

      // Find min y value
      if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jet_pt_ref_unfStatUnc[iPtJInt]), GetTGraphMin (g_jet_pt_ref_unfIterUnc[iPtJInt]));
      else        ymin = 0;

      // Calculate cross-over point to optimize Niter
      const double xopt = GetCrossOverPoint (g_jet_pt_ref_unfStatUnc[iPtJInt], g_jet_pt_ref_unfIterUnc[iPtJInt], 1);

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2}}", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (g_jet_pt_ref_unfTotUnc[iPtJInt], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_ref_unfStatUnc[iPtJInt], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_ref_unfIterUnc[iPtJInt], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, -1);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      if (doLogY) gPad->SetLogy ();

      // Find max y value
      ymax = ymaxSF * GetTGraphMax (g_jet_pt_unfTotUnc[iPtJInt][iCent], 0);

      // Find min y value
      if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jet_pt_unfStatUnc[iPtJInt][iCent]), GetTGraphMin (g_jet_pt_unfIterUnc[iPtJInt][iCent]));
      else        ymin = 0;

      // Calculate cross-over point to optimize Niter
      const double xopt = GetCrossOverPoint (g_jet_pt_unfStatUnc[iPtJInt][iCent], g_jet_pt_unfIterUnc[iPtJInt][iCent], 1);

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2}}", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (g_jet_pt_unfTotUnc[iPtJInt][iCent], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_unfStatUnc[iPtJInt][iCent], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_unfIterUnc[iPtJInt][iCent], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, iCent);
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

      // Find max y value
      ymax = ymaxSF * GetTGraphMax (g_jet_pt_ref_unfTotHybUnc[iPtJInt], 0);

      // Find min y value
      if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jet_pt_ref_unfStatHybUnc[iPtJInt]), GetTGraphMin (g_jet_pt_ref_unfIterHybUnc[iPtJInt]));
      else        ymin = 0;

      // Calculate cross-over point to optimize Niter
      const double xopt = GetCrossOverPoint (g_jet_pt_ref_unfStatHybUnc[iPtJInt], g_jet_pt_ref_unfIterHybUnc[iPtJInt], 1);

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n}", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (g_jet_pt_ref_unfTotHybUnc[iPtJInt], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_ref_unfStatHybUnc[iPtJInt], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_ref_unfIterHybUnc[iPtJInt], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, -1);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      if (doLogY) gPad->SetLogy ();

      // Find max y value
      ymax = ymaxSF * GetTGraphMax (g_jet_pt_unfTotHybUnc[iPtJInt][iCent], 0);

      // Find min y value
      if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jet_pt_unfStatHybUnc[iPtJInt][iCent]), GetTGraphMin (g_jet_pt_unfIterHybUnc[iPtJInt][iCent]));
      else        ymin = 0;

      // Calculate cross-over point to optimize Niter
      const double xopt = GetCrossOverPoint (g_jet_pt_unfStatHybUnc[iPtJInt][iCent], g_jet_pt_unfIterHybUnc[iPtJInt][iCent], 1);

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n}", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (g_jet_pt_unfTotHybUnc[iPtJInt][iCent], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_unfStatHybUnc[iPtJInt][iCent], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_unfIterHybUnc[iPtJInt][iCent], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, iCent);
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
        if (xopt == -1 && (ydum > y || i == g_jet_pt_ref_unfIterRelUnc[iPtJInt]->GetN () - 1))
          xopt = x;
        else if (ydum < y)
          xopt = -1;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n^{2}}", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (g_jet_pt_ref_unfTotRelUnc[iPtJInt], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_ref_unfStatRelUnc[iPtJInt], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_ref_unfIterRelUnc[iPtJInt], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

      myText (0.2, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, -1);
      l->DrawLine (xopt_2d, ymin, xopt_2d, ymax/ymaxSF);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
    }


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      if (doLogY) gPad->SetLogy ();

      // Find max y value
      ymax = ymaxSF * GetTGraphMax (g_jet_pt_unfTotRelUnc[iPtJInt][iCent], 0);

      // Find min y value
      if (doLogY) ymin = 0.5 * std::fmin (GetTGraphMin (g_jet_pt_unfStatRelUnc[iPtJInt][iCent]), GetTGraphMin (g_jet_pt_unfIterRelUnc[iPtJInt][iCent]));
      else        ymin = 0;

      // Calculate cross-over point to optimize Niter
      const double xopt = GetCrossOverPoint (g_jet_pt_unfStatRelUnc[iPtJInt][iCent], g_jet_pt_unfIterRelUnc[iPtJInt][iCent], 1);

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n^{2}}", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetLineWidth (0);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (g_jet_pt_unfTotRelUnc[iPtJInt][iCent], kBlack, kFullCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_unfStatRelUnc[iPtJInt][iCent], kBlue, kOpenCircle, 1.0, 1, 2, "PL", false);
      myDraw (g_jet_pt_unfIterRelUnc[iPtJInt][iCent], kRed, kOpenSquare, 1.0, 1, 2, "PL", false);

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);
      l->DrawLine (xopt, ymin, xopt, ymax/ymaxSF);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, iCent);
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
    const int minJetPt = (iPtJInt == 0 ? 30 : 60);
    const int maxJetPt = 300;

    const char* canvasName = Form ("c_njet_%iGeV_Converge", minJetPt);
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);
  
    TH1D* h = nullptr;
    TGAE* g = nullptr;

    const float ymin = 0.66;//0.95;
    const float ymax = 1.17;//1.075;

    {
      c->cd (7);

      g = new TGAE (nIters1DMax - nIters1DMin + 1);

      int nIter1p = -1;

      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
        double ratio = 0, ratio_err = 0;
        for (short iX = h_jet_pt_ref_unf_nIters[iIter]->FindBin (minJetPt+0.01); iX <= h_jet_pt_ref_unf_nIters[iIter]->FindBin (maxJetPt-0.01); iX++) {
          ratio += h_jet_pt_ref_unf_nIters[iIter]->GetBinContent (iX);
          ratio_err += std::pow (h_jet_pt_ref_unf_nIters[iIter]->GetBinError (iX), 2); 
        }
        ratio_err = std::sqrt (ratio_err);
        //double ratio = h_jet_pt_ref_unf_nIters[iIter]->IntegralAndError (h_jet_pt_ref_unf_nIters[iIter]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[iIter]->FindBin (maxJetPt-0.01), ratio_err);

        double den;
        if (iIter == 0)
          den = h_jet_pt_ref[0]->Integral (h_jet_pt_ref[0]->FindBin (minJetPt+0.01), h_jet_pt_ref[0]->FindBin (maxJetPt-0.01));
        else
          den = h_jet_pt_ref_unf_nIters[iIter-1]->Integral (h_jet_pt_ref_unf_nIters[iIter-1]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[iIter-1]->FindBin (maxJetPt-0.01));

        ratio = ratio / den;
        ratio_err = ratio_err / den;

        if (nIter1p == -1) {
          if (std::fabs (ratio - 1) < 0.01)
          //if (ratio - ratio_err < 1 && ratio + ratio_err > 1)
            nIter1p = nIters1DVals[iIter];
        }
        //else if (ratio - ratio_err > 1 || ratio + ratio_err < 1)
        else if (std::fabs (ratio - 1) > 0.01)
          nIter1p = -1;

        g->SetPoint       (iIter, nIters1DVals[iIter], ratio);
        g->SetPointEYhigh (iIter, ratio_err);
        g->SetPointEYlow  (iIter, ratio_err);
      }
      std::cout << "For minJetPt = " << minJetPt << ", pp, opt nIter = " << nIter1p << std::endl;

      h = new TH1D ("h", ";Iterations;N_{jet} at N_{iter.} / N_{iter.}-1", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);
  
      l->SetLineWidth (2);
      l->SetLineColor (kGray+1);
      l->SetLineStyle (2);
      l->DrawLine (0, 1.01, nIters1DMax, 1.01);
      l->DrawLine (0, 0.99, nIters1DMax, 0.99);

      TGAE* gup = new TGAE ();
      TGAE* gdown = new TGAE ();
      MakeGupAndGdown (g, gup, gdown);
      myDrawFill (gup, gdown, colorfulSystColors[0], 0.7);
      ResetTGAEErrors (g);
      ResetXErrors (g);
      myDraw (g, colorfulColors[0], kDot, 0, 1, 2, "L");
      SaferDelete (&g);
      SaferDelete (&gup);
      SaferDelete (&gdown);

      //myDraw (g, colorfulColors[0], kFullCircle, 1.0, 1, 2, "P", false);
      //SaferDelete (&g);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      l->DrawLine (GetJetSpectraNIters (false, iPtJInt, -1), ymin, GetJetSpectraNIters (false, iPtJInt, -1), ymax);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
      l->DrawLine (nIter1p, ymin, nIter1p, ymax);
  
      myText (0.2, 0.84, kBlack, "#bf{#it{pp}}", 0.06);
    }
  
    for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

      c->cd (nFcalCentBins+1-iCent);

      g = new TGAE (nIters1DMax - nIters1DMin + 1);

      int nIter1p = -1;

      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
        double ratio = 0, ratio_err = 0;
        for (short iX = h_jet_pt_unf_nIters[iCent][iIter]->FindBin (minJetPt+0.01); iX <= h_jet_pt_unf_nIters[iCent][iIter]->FindBin (maxJetPt-0.01); iX++) {
          ratio += h_jet_pt_unf_nIters[iCent][iIter]->GetBinContent (iX);
          ratio_err += std::pow (h_jet_pt_unf_nIters[iCent][iIter]->GetBinError (iX), 2); 
        }
        ratio_err = std::sqrt (ratio_err);
        //double ratio = h_jet_pt_unf_nIters[iCent][iIter]->IntegralAndError (h_jet_pt_unf_nIters[iCent][iIter]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][iIter]->FindBin (maxJetPt-0.01), ratio_err);

        double den;
        if (iIter == 0)
          den = h_jet_pt[0][iCent]->Integral (h_jet_pt[0][iCent]->FindBin (minJetPt+0.01), h_jet_pt[0][iCent]->FindBin (maxJetPt-0.01));
        else
          den = h_jet_pt_unf_nIters[iCent][iIter-1]->Integral (h_jet_pt_unf_nIters[iCent][iIter-1]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][iIter-1]->FindBin (maxJetPt-0.01));

        ratio = ratio / den;
        ratio_err = ratio_err / den;

        if (nIter1p == -1) {
          if (std::fabs (ratio - 1) < 0.01)
          //if (ratio - ratio_err < 1 && ratio + ratio_err > 1)
            nIter1p = nIters1DVals[iIter];
        }
        //else if (ratio - ratio_err > 1 || ratio + ratio_err < 1)
        else if (std::fabs (ratio - 1) > 0.01)
          nIter1p = -1;

        g->SetPoint       (iIter, nIters1DVals[iIter], ratio);
        g->SetPointEYhigh (iIter, ratio_err);
        g->SetPointEYlow  (iIter, ratio_err);
      }
      std::cout << "For minJetPt = " << minJetPt << ", " << zdcCentPercs[iCent+1] << "-" << zdcCentPercs[iCent] << "% central p+Pb, opt nIter = " << nIter1p << std::endl;
  
      h = new TH1D ("h", ";Iterations;N_{jet} at N_{iter.} / N_{iter.}-1", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);
  
      l->SetLineWidth (2);
      l->SetLineColor (kGray+1);
      l->SetLineStyle (2);
      l->DrawLine (0, 1.01, nIters1DMax, 1.01);
      l->DrawLine (0, 0.99, nIters1DMax, 0.99);

      TGAE* gup = new TGAE ();
      TGAE* gdown = new TGAE ();
      MakeGupAndGdown (g, gup, gdown);
      myDrawFill (gup, gdown, colorfulSystColors[iCent+1], 0.7);
      ResetTGAEErrors (g);
      ResetXErrors (g);
      myDraw (g, colorfulColors[iCent+1], kDot, 0, 1, 2, "L");
      SaferDelete (&g);
      SaferDelete (&gup);
      SaferDelete (&gdown);

      //myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.0, 1, 2, "P", false);
      //SaferDelete (&g);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      l->DrawLine (GetJetSpectraNIters (false, iPtJInt, iCent), ymin, GetJetSpectraNIters (false, iPtJInt, iCent), ymax);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
      l->DrawLine (nIter1p, ymin, nIter1p, ymax);
  
      if (iCent < nFcalCentBins)
        myText (0.2, 0.84, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
      else
        myText (0.2, 0.84, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);
  
    } // end loop over iCent
  
    c->cd (8);
    myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.07);
    myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, %i)", minJetPt, maxJetPt), 0.07);
  
    c->SaveAs (Form ("%s/Plots/PtCh/UnfRatioTest_NJet_%s.pdf", workPath.Data (), pTJInt.Data ()));
  
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {
  
    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
    const int minJetPt = (iPtJInt == 0 ? 30 : 60);
    const int maxJetPt = 300;

    const char* canvasName = Form ("c_njet_%iGeV_Converge", minJetPt);
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);
  
    TH1D* h = nullptr;
    TGAE* g = nullptr;

    const float ymin = (iPtJInt == 0 ? 0.30 : 0.66);//0.95;
    const float ymax = (iPtJInt == 0 ? 1.35 : 1.17);//1.075;

    {
      c->cd (7);

      g = new TGAE (nIters1DMax - nIters1DMin + 1);

      int nIter1p = -1;

      double den = h_jet_pt_ref_unf_nIters[nIters1DMax-1]->Integral (h_jet_pt_ref_unf_nIters[nIters1DMax-1]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[nIters1DMax-1]->FindBin (maxJetPt-0.01));
      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
        double ratio = 0, ratio_err = 0;
        for (short iX = h_jet_pt_ref_unf_nIters[iIter]->FindBin (minJetPt+0.01); iX <= h_jet_pt_ref_unf_nIters[iIter]->FindBin (maxJetPt-0.01); iX++) {
          ratio += h_jet_pt_ref_unf_nIters[iIter]->GetBinContent (iX);
          ratio_err += std::pow (h_jet_pt_ref_unf_nIters[iIter]->GetBinError (iX), 2); 
        }
        ratio_err = std::sqrt (ratio_err);
        //double ratio = h_jet_pt_ref_unf_nIters[iIter]->IntegralAndError (h_jet_pt_ref_unf_nIters[iIter]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[iIter]->FindBin (maxJetPt-0.01), ratio_err);

        //double den;
        //if (iIter == 0)
        //  den = h_jet_pt_ref[0]->Integral (h_jet_pt_ref[0]->FindBin (minJetPt+0.01), h_jet_pt_ref[0]->FindBin (maxJetPt-0.01));
        //else
        //  den = h_jet_pt_ref_unf_nIters[iIter-1]->Integral (h_jet_pt_ref_unf_nIters[iIter-1]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[iIter-1]->FindBin (maxJetPt-0.01));

        ratio = ratio / den;
        ratio_err = ratio_err / den;

        if (nIter1p == -1) {
          if (std::fabs (ratio - 1) < 0.01)
          //if (ratio - ratio_err < 1 && ratio + ratio_err > 1)
            nIter1p = nIters1DVals[iIter];
        }
        //else if (ratio - ratio_err > 1 || ratio + ratio_err < 1)
        else if (std::fabs (ratio - 1) > 0.01)
          nIter1p = -1;

        g->SetPoint       (iIter, nIters1DVals[iIter], ratio);
        g->SetPointEYhigh (iIter, ratio_err);
        g->SetPointEYlow  (iIter, ratio_err);
      }
      std::cout << "For minJetPt = " << minJetPt << ", pp, opt nIter = " << nIter1p << std::endl;

      //h = new TH1D ("h", ";Iterations;N_{jet} at N_{iter.} / N_{iter.}-1", 1, 0, nIters1DMax);
      h = new TH1D ("h", Form (";Iterations;N_{jet} at N_{iter.} / %i iters.", nIters1DMax), 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);
  
      l->SetLineWidth (2);
      l->SetLineColor (kGray+1);
      l->SetLineStyle (2);
      l->DrawLine (0, 1.01, nIters1DMax, 1.01);
      l->DrawLine (0, 0.99, nIters1DMax, 0.99);

      TGAE* gup = new TGAE ();
      TGAE* gdown = new TGAE ();
      MakeGupAndGdown (g, gup, gdown);
      myDrawFill (gup, gdown, colorfulSystColors[0], 0.7);
      ResetTGAEErrors (g);
      ResetXErrors (g);
      myDraw (g, colorfulColors[0], kDot, 0, 1, 2, "L");
      SaferDelete (&g);
      SaferDelete (&gup);
      SaferDelete (&gdown);

      //myDraw (g, colorfulColors[0], kFullCircle, 1.0, 1, 2, "P", false);
      //SaferDelete (&g);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      l->DrawLine (GetJetSpectraNIters (false, iPtJInt, -1), ymin, GetJetSpectraNIters (false, iPtJInt, -1), ymax);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
      l->DrawLine (nIter1p, ymin, nIter1p, ymax);
  
      myText (0.2, 0.84, kBlack, "#bf{#it{pp}}", 0.06);
    }
  
    for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

      c->cd (nFcalCentBins+1-iCent);

      g = new TGAE (nIters1DMax - nIters1DMin + 1);

      int nIter1p = -1;

      double den = h_jet_pt_unf_nIters[iCent][nIters1DMax-1]->Integral (h_jet_pt_unf_nIters[iCent][nIters1DMax-1]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][nIters1DMax-1]->FindBin (maxJetPt-0.01));
      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
        double ratio = 0, ratio_err = 0;
        for (short iX = h_jet_pt_unf_nIters[iCent][iIter]->FindBin (minJetPt+0.01); iX <= h_jet_pt_unf_nIters[iCent][iIter]->FindBin (maxJetPt-0.01); iX++) {
          ratio += h_jet_pt_unf_nIters[iCent][iIter]->GetBinContent (iX);
          ratio_err += std::pow (h_jet_pt_unf_nIters[iCent][iIter]->GetBinError (iX), 2); 
        }
        ratio_err = std::sqrt (ratio_err);
        //double ratio = h_jet_pt_unf_nIters[iCent][iIter]->IntegralAndError (h_jet_pt_unf_nIters[iCent][iIter]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][iIter]->FindBin (maxJetPt-0.01), ratio_err);

        //double den;
        //if (iIter == 0)
        //  den = h_jet_pt[0][iCent]->Integral (h_jet_pt[0][iCent]->FindBin (minJetPt+0.01), h_jet_pt[0][iCent]->FindBin (maxJetPt-0.01));
        //else
        //  den = h_jet_pt_unf_nIters[iCent][iIter-1]->Integral (h_jet_pt_unf_nIters[iCent][iIter-1]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][iIter-1]->FindBin (maxJetPt-0.01));

        ratio = ratio / den;
        ratio_err = ratio_err / den;

        if (nIter1p == -1) {
          if (std::fabs (ratio - 1) < 0.01)
          //if (ratio - ratio_err < 1 && ratio + ratio_err > 1)
            nIter1p = nIters1DVals[iIter];
        }
        //else if (ratio - ratio_err > 1 || ratio + ratio_err < 1)
        else if (std::fabs (ratio - 1) > 0.01)
          nIter1p = -1;

        g->SetPoint       (iIter, nIters1DVals[iIter], ratio);
        g->SetPointEYhigh (iIter, ratio_err);
        g->SetPointEYlow  (iIter, ratio_err);
      }
      std::cout << "For minJetPt = " << minJetPt << ", " << zdcCentPercs[iCent+1] << "-" << zdcCentPercs[iCent] << "% central p+Pb, opt nIter = " << nIter1p << std::endl;
  
      //h = new TH1D ("h", ";Iterations;N_{jet} at N_{iter.} / N_{iter.}-1", 1, 0, nIters1DMax);
      h = new TH1D ("h", Form (";Iterations;N_{jet} at N_{iter.} / %i iters.", nIters1DMax), 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);
  
      l->SetLineWidth (2);
      l->SetLineColor (kGray+1);
      l->SetLineStyle (2);
      l->DrawLine (0, 1.01, nIters1DMax, 1.01);
      l->DrawLine (0, 0.99, nIters1DMax, 0.99);

      TGAE* gup = new TGAE ();
      TGAE* gdown = new TGAE ();
      MakeGupAndGdown (g, gup, gdown);
      myDrawFill (gup, gdown, colorfulSystColors[iCent+1], 0.7);
      ResetTGAEErrors (g);
      ResetXErrors (g);
      myDraw (g, colorfulColors[iCent+1], kDot, 0, 1, 2, "L");
      SaferDelete (&g);
      SaferDelete (&gup);
      SaferDelete (&gdown);

      //myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.0, 1, 2, "P", false);
      //SaferDelete (&g);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      l->DrawLine (GetJetSpectraNIters (false, iPtJInt, iCent), ymin, GetJetSpectraNIters (false, iPtJInt, iCent), ymax);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
      l->DrawLine (nIter1p, ymin, nIter1p, ymax);
  
      if (iCent < nFcalCentBins)
        myText (0.2, 0.84, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
      else
        myText (0.2, 0.84, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);
  
    } // end loop over iCent
  
    c->cd (8);
    myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.07);
    myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, %i)", minJetPt, maxJetPt), 0.07);
  
    c->SaveAs (Form ("%s/Plots/PtCh/UnfConv_NJet_%s.pdf", workPath.Data (), pTJInt.Data ()));
  
  } // end loop over iPtJInt



  for (short iPtJInt : {0, 1}) {
  
    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
    const int minJetPt = (iPtJInt == 0 ? 30 : 60);
    const int maxJetPt = 300;

    const char* canvasName = Form ("c_njet_%iGeV_RefoldingComp", minJetPt);
    TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
    c->Divide (4, 2);
  
    TH1D* h = nullptr;
    TGAE* g = nullptr;

    const float ymin = 0.78;//0.95;
    const float ymax = 1.17;//1.075;

    {
      c->cd (7);

      g = new TGAE (nIters1DMax - nIters1DMin + 1);

      int nIter1p = -1;

      double den = h_jet_pt_ref[0]->Integral (h_jet_pt_ref[0]->FindBin (minJetPt+0.01), h_jet_pt_ref[0]->FindBin (maxJetPt-0.01));
      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
        double ratio = 0, ratio_err = 0;
        ratio = h_jet_pt_ref_rfld_nIters[iIter]->IntegralAndError (h_jet_pt_ref_rfld_nIters[iIter]->FindBin (minJetPt+0.01), h_jet_pt_ref_rfld_nIters[iIter]->FindBin (maxJetPt-0.01), ratio_err);

        ratio = ratio / den;
        ratio_err = ratio_err / den;

        if (nIter1p == -1) {
          if (std::fabs (ratio - 1) < 0.01)
            nIter1p = nIters1DVals[iIter];
        }
        else if (std::fabs (ratio - 1) > 0.01)
          nIter1p = -1;

        g->SetPoint       (iIter, nIters1DVals[iIter], ratio);
        g->SetPointEYhigh (iIter, ratio_err);
        g->SetPointEYlow  (iIter, ratio_err);
      }
      std::cout << "For minJetPt = " << minJetPt << ", pp, opt nIter = " << nIter1p << std::endl;

      h = new TH1D ("h", ";Iterations;Refolded N_{jet} / No unfold", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);
  
      l->SetLineWidth (2);
      l->SetLineColor (kGray+1);
      l->SetLineStyle (2);
      l->DrawLine (0, 1.01, nIters1DMax, 1.01);
      l->DrawLine (0, 0.99, nIters1DMax, 0.99);

      TGAE* gup = new TGAE ();
      TGAE* gdown = new TGAE ();
      MakeGupAndGdown (g, gup, gdown);
      myDrawFill (gup, gdown, colorfulSystColors[0], 0.7);
      ResetTGAEErrors (g);
      ResetXErrors (g);
      myDraw (g, colorfulColors[0], kDot, 0, 1, 2, "L");
      SaferDelete (&g);
      SaferDelete (&gup);
      SaferDelete (&gdown);

      //myDraw (g, colorfulColors[0], kFullCircle, 1.0, 1, 2, "P", false);
      //SaferDelete (&g);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      l->DrawLine (GetJetSpectraNIters (false, iPtJInt, -1), ymin, GetJetSpectraNIters (false, iPtJInt, -1), ymax);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
      l->DrawLine (nIter1p, ymin, nIter1p, ymax);
  
      myText (0.2, 0.84, kBlack, "#bf{#it{pp}}", 0.06);
    }
  
    for (short iCent = 0; iCent < nFcalCentBins+1; iCent++) {

      c->cd (nFcalCentBins+1-iCent);

      g = new TGAE (nIters1DMax - nIters1DMin + 1);

      int nIter1p = -1;

      double den = h_jet_pt[0][iCent]->Integral (h_jet_pt[0][iCent]->FindBin (minJetPt+0.01), h_jet_pt[0][iCent]->FindBin (maxJetPt-0.01));
      for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+1; iIter++) {
        double ratio = 0, ratio_err = 0;
        ratio = h_jet_pt_rfld_nIters[iCent][iIter]->IntegralAndError (h_jet_pt_rfld_nIters[iCent][iIter]->FindBin (minJetPt+0.01), h_jet_pt_rfld_nIters[iCent][iIter]->FindBin (maxJetPt-0.01), ratio_err);

        ratio = ratio / den;
        ratio_err = ratio_err / den;

        if (nIter1p == -1) {
          if (std::fabs (ratio - 1) < 0.01)
            nIter1p = nIters1DVals[iIter];
        }
        else if (std::fabs (ratio - 1) > 0.01)
          nIter1p = -1;

        g->SetPoint       (iIter, nIters1DVals[iIter], ratio);
        g->SetPointEYhigh (iIter, ratio_err);
        g->SetPointEYlow  (iIter, ratio_err);
      }
      std::cout << "For minJetPt = " << minJetPt << ", " << zdcCentPercs[iCent+1] << "-" << zdcCentPercs[iCent] << "% central p+Pb, opt nIter = " << nIter1p << std::endl;
  
      h = new TH1D ("h", ";Iterations;Refolded N_{jet} / No unfold", 1, 0, nIters1DMax);
      h->GetYaxis ()->SetRangeUser (ymin, ymax);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);
  
      l->SetLineWidth (2);
      l->SetLineColor (kGray+1);
      l->SetLineStyle (2);
      l->DrawLine (0, 1.01, nIters1DMax, 1.01);
      l->DrawLine (0, 0.99, nIters1DMax, 0.99);

      TGAE* gup = new TGAE ();
      TGAE* gdown = new TGAE ();
      MakeGupAndGdown (g, gup, gdown);
      myDrawFill (gup, gdown, colorfulSystColors[iCent+1], 0.7);
      ResetTGAEErrors (g);
      ResetXErrors (g);
      myDraw (g, colorfulColors[iCent+1], kDot, 0, 1, 2, "L");
      SaferDelete (&g);
      SaferDelete (&gup);
      SaferDelete (&gdown);

      //myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.0, 1, 2, "P", false);
      //SaferDelete (&g);

      l->SetLineColor (myGreen);
      l->SetLineStyle (3);
      l->DrawLine (GetJetSpectraNIters (false, iPtJInt, iCent), ymin, GetJetSpectraNIters (false, iPtJInt, iCent), ymax);
      l->SetLineColor (kBlack);
      l->SetLineStyle (2);
      l->DrawLine (nIter1p, ymin, nIter1p, ymax);
  
      if (iCent < nFcalCentBins)
        myText (0.2, 0.84, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.06);
      else
        myText (0.2, 0.84, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.06);
  
    } // end loop over iCent
  
    c->cd (8);
    myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Internal", 0.07);
    myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.07);
    myText (0.1, 0.57, kBlack, Form ("#it{p}_{T}^{jet} [GeV] #in (%i, %i)", minJetPt, maxJetPt), 0.07);
  
    c->SaveAs (Form ("%s/Plots/PtCh/Refolded_NJet_%s.pdf", workPath.Data (), pTJInt.Data ()));
  
  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
    const double maxJetPt = 300;

    //for (short iDir : {0, 1, 2}) {
    for (short iDir : {0, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_sigVsUnf_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGAE* g = nullptr;

      {
        c->cd (7);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];N_{iter} iterations / N_{iter}-1 iterations", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.86, 1.22);
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

          TH1D* hprev = (iIter == 0 ? h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir] : h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIter-1]);
          TH1D* hjetiter = h_jet_pt_ref_unf_nIters[iIter];
          TH1D* hjetprev = (iIter == 0 ? h_jet_pt_ref[0] : h_jet_pt_ref_unf_nIters[iIter-1]);

          const double nJetsIter = hjetiter->Integral (hjetiter->FindBin (minJetPt+0.01), hjetiter->FindBin (maxJetPt-0.01));
          const double nJetsPrev = hjetprev->Integral (hjetprev->FindBin (minJetPt+0.01), hjetprev->FindBin (maxJetPt-0.01));

          ScaleGraph (g, hprev, nJetsIter / nJetsPrev);
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
        h->GetYaxis ()->SetRangeUser (0.86, 1.22);
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

          TH1D* hprev = (iIter == 0 ? h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent] : h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIter-1]);
          TH1D* hjetiter = h_jet_pt_unf_nIters[iCent][iIter];
          TH1D* hjetprev = (iIter == 0 ? h_jet_pt[0][iCent] : h_jet_pt_unf_nIters[iCent][iIter-1]);

          const double nJetsIter = hjetiter->Integral (hjetiter->FindBin (minJetPt+0.01), hjetiter->FindBin (maxJetPt-0.01));
          const double nJetsPrev = hjetprev->Integral (hjetprev->FindBin (minJetPt+0.01), hjetprev->FindBin (maxJetPt-0.01));

          ScaleGraph (g, hprev, nJetsIter / nJetsPrev);
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

    const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
    const double maxJetPt = 300;

    //for (short iDir : {0, 1, 2}) {
    for (short iDir : {0, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_sigVsRfld_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGAE* g = nullptr;

      //std::vector <int> iIters = {
      //  iPtJInt == 0 ? 19 : 11,
      //  iPtJInt == 0 ? 16 : 9,
      //  iPtJInt == 0 ? 12 : 7,
      //  iPtJInt == 0 ? 9 : 5,
      //  iPtJInt == 0 ? 6 : 3,
      //  iPtJInt == 0 ? 4 : 2,
      //  iPtJInt == 0 ? 2 : 1,
      //  0
      //};
      std::vector <int> iIters = {11, 9, 7, 5, 3, 2, 1, 0};

      {
        c->cd (7);

        gPad->SetLogx ();

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Refolded N_{iter} iterations / No unfold", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.86, 1.22);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);


        short iCol = 7;
        for (short iIter : iIters) {
          //if (iIter > GetTrkSpectraNIters (false, iPtJInt, iDir, -1)) {
          //  iCol--;
          //  continue;
          //}

          h = (TH1D*) h_jetInt_trk_pt_ref_rfld_nIters[iPtJInt][iDir][iIter]->Clone ("h");
          TH1D* hprev = h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir];

          h->Divide (hprev);
          g = make_graph (h);
          SaferDelete (&h);

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

        TH1D* h = new TH1D ("h", ";#it{p}_{T}^{ch} [GeV];Refolded N_{iter} iterations / No unfold", 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.86, 1.22);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        short iCol = 7;
        for (short iIter : iIters) {
          //if (iIter > GetTrkSpectraNIters (false, iPtJInt, iDir, iCent)) {
          //  iCol--;
          //  continue;
          //}

          h = (TH1D*) h_jetInt_trk_pt_rfld_nIters[iPtJInt][iDir][iCent][iIter]->Clone ("h");
          TH1D* hprev = h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent];

          h->Divide (hprev);
          g = make_graph (h);
          SaferDelete (&h);

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
      for (short iIter : iIters) {
        myLineText2 (0.15+0.4*(iCol / 4), 0.48-0.07*(iCol%4), colorfulColors[iCol], kOpenCircle, Form ("%i iterations", iIter+1), 1.2, 0.055); iCol--;
      }

      c->SaveAs (Form ("%s/Plots/PtCh/Refolded_All_%iGeVJets_%s.pdf", workPath.Data (), iPtJInt == 0 ? 30 : 60, directions[iDir] == "ns" ? "nearside" : (directions[iDir] == "as" ? "awayside" : "perpendicular")));

    } // end loop over iDir

  } // end loop over iPtJInt




  for (short iPtJInt : {0, 1}) {

    const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

    const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
    const double maxJetPt = 300;

    //for (short iDir : {0, 1, 2}) {
    for (short iDir : {0, 2}) {

      const char* canvasName = Form ("c_jetInt_trk_pt_sigVsUnf_Summary_iDir%i_%s", iDir, pTJInt.Data ());
      TCanvas* c = new TCanvas (canvasName, "", 1400, 700);
      c->Divide (4, 2);

      TGAE* g = nullptr;

      {
        c->cd (7);

        gPad->SetLogx ();

        const short nIterOpt = GetTrkSpectraNIters (false, iPtJInt, iDir, -1);//jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][iPtJInt == 0 ? 1 : 0];
        const short nIterVar = nIterOpt + 1;//jetInt_trk_pt_ref_niter_opt[iPtJInt][iDir][0];

        short iIterOpt = 0;
        while (iIterOpt < nItersMax-nItersMin+1 && nItersVals[iIterOpt] != nIterOpt) iIterOpt++;
        short iIterVar = iIterOpt+1;

        TH1D* h = new TH1D ("h", Form (";#it{p}_{T}^{ch} [GeV];%i vs. %i iterations", nIterOpt, nIterVar), 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.97, 1.03);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        const double nJetsOpt = h_jet_pt_ref_unf_nIters[iIterOpt]->Integral (h_jet_pt_ref_unf_nIters[iIterOpt]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[iIterOpt]->FindBin (maxJetPt-0.01));
        const double nJetsVar = h_jet_pt_ref_unf_nIters[iIterVar]->Integral (h_jet_pt_ref_unf_nIters[iIterVar]->FindBin (minJetPt+0.01), h_jet_pt_ref_unf_nIters[iIterVar]->FindBin (maxJetPt-0.01));

        TH1D* hopt = (TH1D*) h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIterOpt]->Clone ("hden");
        hopt->Scale (nJetsOpt);
        TH1D* hvar = (TH1D*) h_jetInt_trk_pt_ref_unf_nIters[iPtJInt][iDir][iIterVar]->Clone ("hnum");
        hvar->Scale (nJetsVar);

        g = make_graph (hvar);
        ResetTGAEErrors (g);
        ResetXErrors (g);
        CalcSystematics (g, hopt, false);
        ScaleGraph (g, hopt);
        if (iDir == 1) TrimGraph (g, 0, 10);

        SaferDelete (&hopt);
        SaferDelete (&hvar);

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

        const short nIterOpt = GetTrkSpectraNIters (false, iPtJInt, iDir, iCent);//jetInt_trk_pt_niter_opt[iPtJInt][iDir][iCent][iPtJInt == 0 ? 1 : 0];
        const short nIterVar = nIterOpt + 1;//jetInt_trk_pt_niter_opt[iPtJInt][iDir][iCent][0];

        short iIterOpt = 0;
        while (iIterOpt < nItersMax-nItersMin+1 && nItersVals[iIterOpt] != nIterOpt) iIterOpt++;
        short iIterVar = iIterOpt+1;

        TH1D* h = new TH1D ("h", Form (";#it{p}_{T}^{ch} [GeV];%i vs. %i iterations", nIterOpt, nIterVar), 1, pTChBins[1], iDir == 1 ? 10 : pTChBins[nPtChBins-(iPtJInt == 0 ? 3 : 1)]);//pTChBins[nPtChBins]);
        h->GetXaxis ()->SetMoreLogLabels ();
        h->GetYaxis ()->SetRangeUser (0.97, 1.03);
        h->SetBinContent (1, 1);
        h->SetLineStyle (2);
        h->SetLineWidth (2);
        h->SetLineColor (kBlack);
        h->DrawCopy ("hist ][");
        SaferDelete (&h);

        const double nJetsOpt = h_jet_pt_unf_nIters[iCent][iIterOpt]->Integral (h_jet_pt_unf_nIters[iCent][iIterOpt]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][iIterOpt]->FindBin (maxJetPt-0.01));
        const double nJetsVar = h_jet_pt_unf_nIters[iCent][iIterVar]->Integral (h_jet_pt_unf_nIters[iCent][iIterVar]->FindBin (minJetPt+0.01), h_jet_pt_unf_nIters[iCent][iIterVar]->FindBin (maxJetPt-0.01));

        // for Niter / Niter-1 ratio plot
        TH1D* hopt = (TH1D*) h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIterOpt]->Clone ("hden");
        hopt->Scale (nJetsOpt);
        TH1D* hvar = (TH1D*) h_jetInt_trk_pt_unf_nIters[iPtJInt][iDir][iCent][iIterVar]->Clone ("hnum");
        hvar->Scale (nJetsVar);

        g = make_graph (hvar);
        ResetTGAEErrors (g);
        ResetXErrors (g);
        CalcSystematics (g, hopt, false);
        ScaleGraph (g, hopt);
        if (iDir == 1) TrimGraph (g, 0, 10);

        SaferDelete (&hopt);
        SaferDelete (&hvar);

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
