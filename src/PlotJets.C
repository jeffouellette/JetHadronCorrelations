#ifndef __JetHadronCorrelatorPlotJets_C__
#define __JetHadronCorrelatorPlotJets_C__

#include <iostream>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLorentzVector.h>

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
#include "ProcessUnfolding.C"

using namespace JetHadronCorrelations;


TLine* l = new TLine ();
TLatex* tl = new TLatex ();


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



void PlotJets (const char* tag, const char* inFileTag) {

  const TString var = variations[0];
  const short nSamps = 6; // JZ0, 1, 2, 3, JZ1-3, JZ0-3


  const short nIters1DMax = 20;
  const short nIters1DMin = 1;
  const double* nIters1DVals = linspace (nIters1DMin, nIters1DMax, nIters1DMax-nIters1DMin);


  const bool doLogY = true;
  const float ymaxSF = (doLogY ? 5 : 1.2);


  TFile* inFile = nullptr;

  TH1D***  h_evt_counts_ref     = Get2DArray <TH1D*> (2, nSamps);
  TH1D**** h_evt_counts         = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH1D*** h_jet_pt_ref          = Get2DArray <TH1D*> (3, nSamps);

  TH1D**** h_jet_pt             = Get3DArray <TH1D*> (3, nZdcCentBins+1, nSamps);

  TH1D**** h_jet_pt_ratio       = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TH1D***   h_jet_pt_datamc_ratio_ref  = Get2DArray <TH1D*> (2, nSamps);
  TH1D****  h_jet_pt_datamc_ratio      = Get3DArray <TH1D*> (2, nZdcCentBins+1, nSamps);

  TF1***    f_jet_pt_datamc_ratio_ref  = Get2DArray <TF1*> (2, nSamps);
  TF1****   f_jet_pt_datamc_ratio      = Get3DArray <TF1*> (2, nZdcCentBins+1, nSamps);


  // unfolded jet spectra
  TH1D**    h_jet_pt_ref_unf_nIters         = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+2);
  TH1D***   h_jet_pt_unf_nIters             = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+2);

  // refolded jet spectra
  TH1D**    h_jet_pt_ref_rfld_nIters        = Get1DArray <TH1D*> (nIters1DMax-nIters1DMin+2);
  TH1D***   h_jet_pt_rfld_nIters            = Get2DArray <TH1D*> (nZdcCentBins+1, nIters1DMax-nIters1DMin+2);


  //TH2D****  h2_jet_eta_phi_ref   = Get3DArray <TH2D*> (2, nPtJBins, nSamps);
  //TH2D***** h2_jet_eta_phi       = Get4DArray <TH2D*> (2, nPtJBins, nZdcCentBins+1, nSamps);

  TH2D***   h2_jet_pt_eta_jer_frac_num_ref  = Get2DArray <TH2D*> (2, nSamps);
  TH2D****  h2_jet_pt_eta_jer_frac_num      = Get3DArray <TH2D*> (2, nZdcCentBins+1, nSamps);

  TH2D***   h2_jet_pt_eta_jer_frac_den_ref  = Get2DArray <TH2D*> (2, nSamps);
  TH2D****  h2_jet_pt_eta_jer_frac_den      = Get3DArray <TH2D*> (2, nZdcCentBins+1, nSamps);

  TH2D***   h2_jet_pt_eta_jer_frac_ref      = Get2DArray <TH2D*> (2, nSamps);
  TH2D****  h2_jet_pt_eta_jer_frac          = Get3DArray <TH2D*> (2, nZdcCentBins+1, nSamps);


  // alpha=0 (absolute) uncertainty studies
  TGraph**    g_jet_pt_ref_unfIterUnc     = Get1DArray <TGraph*> (2);                 // sums of iterations uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfIterUnc         = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of iterations uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfStatUnc     = Get1DArray <TGraph*> (2);                 // sums of statistical uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfStatUnc         = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of statistical uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfTotUnc      = Get1DArray <TGraph*> (2);                 // sums of total uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfTotUnc          = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of total uncertainties as a function of nIter -- data only

  // alpha=0.5 (hybrid) uncertainty studies
  TGraph**    g_jet_pt_ref_unfIterRelUnc  = Get1DArray <TGraph*> (2);                 // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfIterRelUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of iterations relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfStatRelUnc  = Get1DArray <TGraph*> (2);                 // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfStatRelUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of statistical relative uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfTotRelUnc   = Get1DArray <TGraph*> (2);                 // sums of total relative uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfTotRelUnc       = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of total relative uncertainties as a function of nIter -- data only

  // alpha=1 (relative) uncertainty studies
  TGraph**    g_jet_pt_ref_unfIterHybUnc  = Get1DArray <TGraph*> (2);                 // sums of iterations hybrid uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfIterHybUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of iterations hybrid uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfStatHybUnc  = Get1DArray <TGraph*> (2);                 // sums of statistical hybrid uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfStatHybUnc      = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of statistical hybrid uncertainties as a function of nIter -- data only
  TGraph**    g_jet_pt_ref_unfTotHybUnc   = Get1DArray <TGraph*> (2);                 // sums of total hybrid uncertainties as a function of nIter -- data only
  TGraph***   g_jet_pt_unfTotHybUnc       = Get2DArray <TGraph*> (2, nZdcCentBins+1); // sums of total hybrid uncertainties as a function of nIter -- data only


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessJets_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (short iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (short iSamp = 0; iSamp < nSamps; iSamp++) {

        if (iDType == 0 && iSamp > 0)
          continue;

        const TString samp = GetSamp (iDType, iSamp);

        h_evt_counts_ref[iDType][iSamp] = (TH1D*) inFile->Get (Form ("h_evt_counts_ref_%s_%s",    dType.Data (), samp.Data ()));

        h_jet_pt_ref[iDType][iSamp]     = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_%s_%s",        dType.Data (), samp.Data ()));

        for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+2; iIter++) {

          const short nIters = (short) nIters1DVals[iIter];
          h_jet_pt_ref_unf_nIters[iIter] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_unf_data_nIters%i", nIters));
          h_jet_pt_ref_rfld_nIters[iIter] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_rfld_data_nIters%i", nIters));

        } // end loop over iIter

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        //  h2_jet_eta_phi_ref[iDType][iPtJ][iSamp]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_ref_%s_%s", dType.Data (), samp.Data ()));

        //} // end loop over iPtJ

        h2_jet_pt_eta_jer_frac_num_ref[iDType][iSamp] = (TH2D*) inFile->Get (Form ("h2_jet_pt_eta_jer_frac_num_ref_%s_%s",  dType.Data (), samp.Data ()));
        h2_jet_pt_eta_jer_frac_den_ref[iDType][iSamp] = (TH2D*) inFile->Get (Form ("h2_jet_pt_eta_jer_frac_den_ref_%s_%s",  dType.Data (), samp.Data ()));
        h2_jet_pt_eta_jer_frac_ref[iDType][iSamp]     = (TH2D*) inFile->Get (Form ("h2_jet_pt_eta_jer_frac_ref_%s_%s",      dType.Data (), samp.Data ()));
  

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
  
          h_evt_counts[iDType][iCent][iSamp]    = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_%s_%s_%s",   cent.Data (), dType.Data (), samp.Data ()));

          h_jet_pt[iDType][iCent][iSamp]        = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_%s_%s",       cent.Data (), dType.Data (), samp.Data ()));
          h_jet_pt_ratio[iDType][iCent][iSamp]  = (TH1D*) inFile->Get (Form ("h_jet_pt_ratio_%s_%s_%s",     cent.Data (), dType.Data (), samp.Data ()));

          for (short iIter = 0; iIter < nIters1DMax-nIters1DMin+2; iIter++) {
  
            const short nIters = (short) nIters1DVals[iIter];
            h_jet_pt_unf_nIters[iCent][iIter]   = (TH1D*) inFile->Get (Form ("h_jet_pt_unf_data_%s_nIters%i",       cent.Data (), nIters));
            h_jet_pt_rfld_nIters[iCent][iIter]  = (TH1D*) inFile->Get (Form ("h_jet_pt_rfld_data_%s_nIters%i",  cent.Data (), nIters));
  
          } // end loop over iIter

          //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          //  const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          //  h2_jet_eta_phi[iDType][iPtJ][iCent][iSamp]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_pPb_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), samp.Data ()));

          //} // end loop over iPtJ

          h2_jet_pt_eta_jer_frac_num[iDType][iCent][iSamp] = (TH2D*) inFile->Get  (Form ("h2_jet_pt_eta_jer_frac_num_%s_%s_%s", cent.Data (), dType.Data (), samp.Data ()));
          h2_jet_pt_eta_jer_frac_den[iDType][iCent][iSamp] = (TH2D*) inFile->Get  (Form ("h2_jet_pt_eta_jer_frac_den_%s_%s_%s", cent.Data (), dType.Data (), samp.Data ()));
          h2_jet_pt_eta_jer_frac[iDType][iCent][iSamp]     = (TH2D*) inFile->Get  (Form ("h2_jet_pt_eta_jer_frac_%s_%s_%s",     cent.Data (), dType.Data (), samp.Data ()));
  
        } // end loop over iCent

      } // end loop over iSamp

    } // end loop over iDType


    for (short iSamp = 0; iSamp < nSamps; iSamp++) {

      const TString samp = GetSamp (1, iSamp);

      h_jet_pt_ref[2][iSamp]              = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_mcScaled_%s",            samp.Data ()));

      h_jet_pt_datamc_ratio_ref[0][iSamp] = (TH1D*) inFile->Get (Form ("h_jet_pt_datamc_ratio_ref_%s",        samp.Data ()));
      f_jet_pt_datamc_ratio_ref[0][iSamp] = (TF1*)  inFile->Get (Form ("f_jet_pt_datamc_ratio_ref_%s",        samp.Data ()));

      h_jet_pt_datamc_ratio_ref[1][iSamp] = (TH1D*) inFile->Get (Form ("h_jet_pt_datamcScaled_ratio_ref_%s",  samp.Data ()));
      f_jet_pt_datamc_ratio_ref[1][iSamp] = (TF1*)  inFile->Get (Form ("f_jet_pt_datamcScaled_ratio_ref_%s",  samp.Data ()));

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        h_jet_pt[2][iCent][iSamp]               = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_mcScaled_%s",       cent.Data (), samp.Data ()));

        h_jet_pt_datamc_ratio[0][iCent][iSamp]  = (TH1D*) inFile->Get (Form ("h_jet_pt_datamc_ratio_%s_%s",       cent.Data (), samp.Data ()));
        f_jet_pt_datamc_ratio[0][iCent][iSamp]  = (TF1*)  inFile->Get (Form ("f_jet_pt_datamc_ratio_%s_%s",       cent.Data (), samp.Data ()));

        h_jet_pt_datamc_ratio[1][iCent][iSamp]  = (TH1D*) inFile->Get (Form ("h_jet_pt_datamcScaled_ratio_%s_%s", cent.Data (), samp.Data ()));
        f_jet_pt_datamc_ratio[1][iCent][iSamp]  = (TF1*)  inFile->Get (Form ("f_jet_pt_datamcScaled_ratio_%s_%s", cent.Data (), samp.Data ()));

      } // end loop over iCent

    } // end loop over iSamp

  }




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS SUMMARIZING UNCERTAINTIES VS. # OF ITERATIONS ON THE UNFOLDING
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    const short iDType = 0;
    const TString dType = "data";

    const short iSamp = 0;

    for (short iPtJInt : {0, 1}) {

      const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
      const double maxJetPt = 300;

      TH1D* h_unf = h_jet_pt_ref[iDType][iSamp];
      TH1D* h_unf_prev = h_jet_pt_ref[iDType][iSamp];

      g_jet_pt_ref_unfStatUnc[iPtJInt]     = new TGraph (nIters1DMax - nIters1DMin + 1);
      g_jet_pt_ref_unfIterUnc[iPtJInt]     = new TGraph (nIters1DMax - nIters1DMin + 1);
      g_jet_pt_ref_unfTotUnc[iPtJInt]      = new TGraph (nIters1DMax - nIters1DMin + 1);
      g_jet_pt_ref_unfStatHybUnc[iPtJInt]  = new TGraph (nIters1DMax - nIters1DMin + 1);
      g_jet_pt_ref_unfIterHybUnc[iPtJInt]  = new TGraph (nIters1DMax - nIters1DMin + 1);
      g_jet_pt_ref_unfTotHybUnc[iPtJInt]   = new TGraph (nIters1DMax - nIters1DMin + 1);
      g_jet_pt_ref_unfStatRelUnc[iPtJInt]  = new TGraph (nIters1DMax - nIters1DMin + 1);
      g_jet_pt_ref_unfIterRelUnc[iPtJInt]  = new TGraph (nIters1DMax - nIters1DMin + 1);
      g_jet_pt_ref_unfTotRelUnc[iPtJInt]   = new TGraph (nIters1DMax - nIters1DMin + 1);


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
          iterVar += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX));
          iterHybVar += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / std::sqrt (h_unf_prev->GetBinContent (iX));
          iterRelVar += std::fabs ((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / h_unf_prev->GetBinContent (iX));
        } // end loop over iX
        iterVar *= iterVar;
        iterHybVar *= iterHybVar;
        iterRelVar *= iterRelVar;

        g_jet_pt_ref_unfStatUnc[iPtJInt]->SetPoint    (nIters - nIters1DMin, nIters, std::sqrt (totVar));
        g_jet_pt_ref_unfIterUnc[iPtJInt]->SetPoint    (nIters - nIters1DMin, nIters, std::sqrt (iterVar));
        g_jet_pt_ref_unfTotUnc[iPtJInt]->SetPoint     (nIters - nIters1DMin, nIters, std::sqrt (totVar + iterVar));
        g_jet_pt_ref_unfStatHybUnc[iPtJInt]->SetPoint (nIters - nIters1DMin, nIters, std::sqrt (totHybVar));
        g_jet_pt_ref_unfIterHybUnc[iPtJInt]->SetPoint (nIters - nIters1DMin, nIters, std::sqrt (iterHybVar));
        g_jet_pt_ref_unfTotHybUnc[iPtJInt]->SetPoint  (nIters - nIters1DMin, nIters, std::sqrt (totHybVar + iterHybVar));
        g_jet_pt_ref_unfStatRelUnc[iPtJInt]->SetPoint (nIters - nIters1DMin, nIters, std::sqrt (totRelVar));
        g_jet_pt_ref_unfIterRelUnc[iPtJInt]->SetPoint (nIters - nIters1DMin, nIters, std::sqrt (iterRelVar));
        g_jet_pt_ref_unfTotRelUnc[iPtJInt]->SetPoint  (nIters - nIters1DMin, nIters, std::sqrt (totRelVar + iterRelVar));

        h_unf_prev = h_unf;
      } // end loop over iIter

    } // end loop over iPtJInt


    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

      for (short iPtJInt : {0, 1}) {

        const double minJetPt = (iPtJInt == 0 ? 30. : 60.);
        const double maxJetPt = 300;

        TH1D* h_unf = h_jet_pt[iDType][iCent][iSamp];
        TH1D* h_unf_prev = h_jet_pt[iDType][iCent][iSamp];

        g_jet_pt_unfStatUnc[iPtJInt][iCent]     = new TGraph (nIters1DMax - nIters1DMin + 1);
        g_jet_pt_unfIterUnc[iPtJInt][iCent]     = new TGraph (nIters1DMax - nIters1DMin + 1);
        g_jet_pt_unfTotUnc[iPtJInt][iCent]      = new TGraph (nIters1DMax - nIters1DMin + 1);
        g_jet_pt_unfStatHybUnc[iPtJInt][iCent]  = new TGraph (nIters1DMax - nIters1DMin + 1);
        g_jet_pt_unfIterHybUnc[iPtJInt][iCent]  = new TGraph (nIters1DMax - nIters1DMin + 1);
        g_jet_pt_unfTotHybUnc[iPtJInt][iCent]   = new TGraph (nIters1DMax - nIters1DMin + 1);
        g_jet_pt_unfStatRelUnc[iPtJInt][iCent]  = new TGraph (nIters1DMax - nIters1DMin + 1);
        g_jet_pt_unfIterRelUnc[iPtJInt][iCent]  = new TGraph (nIters1DMax - nIters1DMin + 1);
        g_jet_pt_unfTotRelUnc[iPtJInt][iCent]   = new TGraph (nIters1DMax - nIters1DMin + 1);


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
            iterVar += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX));
            iterHybVar += std::fabs (h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / std::sqrt (h_unf_prev->GetBinContent (iX));
            iterRelVar += std::fabs ((h_unf->GetBinContent (iX) - h_unf_prev->GetBinContent (iX)) / h_unf_prev->GetBinContent (iX));
          } // end loop over iX
          iterVar *= iterVar;
          iterHybVar *= iterHybVar;
          iterRelVar *= iterRelVar;

          g_jet_pt_unfStatUnc[iPtJInt][iCent]->SetPoint     (nIters - nIters1DMin, nIters, std::sqrt (totVar));
          g_jet_pt_unfIterUnc[iPtJInt][iCent]->SetPoint     (nIters - nIters1DMin, nIters, std::sqrt (iterVar));
          g_jet_pt_unfTotUnc[iPtJInt][iCent]->SetPoint      (nIters - nIters1DMin, nIters, std::sqrt (totVar + iterVar));
          g_jet_pt_unfStatHybUnc[iPtJInt][iCent]->SetPoint  (nIters - nIters1DMin, nIters, std::sqrt (totHybVar));
          g_jet_pt_unfIterHybUnc[iPtJInt][iCent]->SetPoint  (nIters - nIters1DMin, nIters, std::sqrt (iterHybVar));
          g_jet_pt_unfTotHybUnc[iPtJInt][iCent]->SetPoint   (nIters - nIters1DMin, nIters, std::sqrt (totHybVar + iterHybVar));
          g_jet_pt_unfStatRelUnc[iPtJInt][iCent]->SetPoint  (nIters - nIters1DMin, nIters, std::sqrt (totRelVar));
          g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->SetPoint  (nIters - nIters1DMin, nIters, std::sqrt (iterRelVar));
          g_jet_pt_unfTotRelUnc[iPtJInt][iCent]->SetPoint   (nIters - nIters1DMin, nIters, std::sqrt (totRelVar + iterRelVar));

          h_unf_prev = h_unf;

        } // end loop over iIter

      } // end loop over iPtJInt

    } // end loop over iCent
  }





  //for (float trigpt : {30., 60.}) {
//  float maxpt = (trigpt == 30. ? 60. : 300.);
  //  float* njet = new float[nZdcCentBins+2];

  //  std::cout << "---------------" << std::endl << "JETS IN DATA > " << trigpt << " GeV" << std::endl << "---------------" << std::endl;

  //  njet[0] = 0;
  //  for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  //    if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
  //      njet[0] += h_jet_counts_ref[0][iPtJ][0]->GetBinContent (1);
  //  } // end loop over iPtJ
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  //    njet[iCent+1] = 0;
  //    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  //      if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
  //        njet[iCent+1] += h_jet_counts[0][iPtJ][iCent][0]->GetBinContent (1);
  //    } // end loop over iPtJ
  //  } // end loop over iCent

  //  std::cout << "Number of pp jets: " << njet[0] << std::endl;
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
  //    std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << njet[iCent+1] << std::endl;

  //  std::cout << "Formatted for latex:" << std::endl;
  //  std::cout << (int) njet[0];
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
  //    std::cout << " & " << (int) njet[iCent+1];
  //  std::cout << std::endl << std::endl;


  //  float integral = 0;
  //  for (short iX = h_jet_pt_ref[0][0]->FindBin (trigpt); iX <= h_jet_pt_ref[0][0]->GetNbinsX (); iX++)
  //    integral += h_jet_pt_ref[0][0]->GetBinContent (iX) * h_jet_pt_ref[0][0]->GetBinWidth (iX);
  //  std::cout << "Average number of jets per trigger jet in pp: " << integral << std::endl;
  //  integral = 0;
  //  for (short iX = h_jet_pt[0][nZdcCentBins][0]->FindBin (trigpt); iX <= h_jet_pt[0][nZdcCentBins][0]->GetNbinsX (); iX++)
  //    integral += h_jet_pt[0][nZdcCentBins][0]->GetBinContent (iX) * h_jet_pt[0][nZdcCentBins][0]->GetBinWidth (iX);
  //  std::cout << "Average number of jets per trigger jet in p+Pb: " << integral << std::endl;
  //  std::cout << std::endl << std::endl;


  //  std::cout << "---------------" << std::endl << "JETS IN MC > " << trigpt << " GeV" << std::endl << "---------------" << std::endl;

  //  njet[0] = 0;
  //  for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  //    if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
  //      njet[0] += h_jet_counts_ref[1][iPtJ][0]->GetBinContent (1);
  //  } // end loop over iPtJ
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
  //    njet[iCent+1] = 0;
  //    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
  //      if (0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) >= trigpt && 0.5*(pTJBins[iPtJ]+pTJBins[iPtJ+1]) < maxpt)
  //        njet[iCent+1] += h_jet_counts[1][iPtJ][iCent][0]->GetBinContent (1);
  //    } // end loop over iPtJ
  //  } // end loop over iCent

  //  std::cout << "Number of pp jets: " << njet[0] << std::endl;
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
  //    std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << njet[iCent+1] << std::endl;

  //  std::cout << "Formatted for latex:" << std::endl;
  //  std::cout << (int) njet[0];
  //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++)
  //    std::cout << " & " << (int) njet[iCent+1];
  //  std::cout << std::endl;


  //  integral = 0;
  //  for (short iX = h_jet_pt_ref[1][0]->FindBin (trigpt); iX <= h_jet_pt_ref[1][0]->GetNbinsX (); iX++)
  //    integral += h_jet_pt_ref[1][0]->GetBinContent (iX) * h_jet_pt_ref[1][0]->GetBinWidth (iX);
  //  std::cout << "Average number of jets per trigger jet in pp: " << integral << std::endl;
  //  integral = 0;
  //  for (short iX = h_jet_pt[1][nZdcCentBins][0]->FindBin (trigpt); iX <= h_jet_pt[1][nZdcCentBins][0]->GetNbinsX (); iX++)
  //    integral += h_jet_pt[1][nZdcCentBins][0]->GetBinContent (iX) * h_jet_pt[1][nZdcCentBins][0]->GetBinWidth (iX);
  //  std::cout << "Average number of jets per trigger jet in p+Pb: " << integral << std::endl;
  //  std::cout << std::endl << std::endl;

  //  delete[] njet;
  //}



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (short iDType = 0; iDType < 2; iDType++) {

    const char* canvasName = Form ("c_jet_pt");

    const short iSamp = (iDType == 1 ? nSamps-1 : 0);

    TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
    c->cd ();

    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.12);
    c->SetLeftMargin (0.12);
    c->SetRightMargin (0.03);

    TH1D* h = nullptr;
    TGAE* g = nullptr;

    double ymin=2e-13;
    double ymax=1e3;

    c->SetLogx();
    c->SetLogy ();

    const double maxx = 400;
    h = new TH1D ("htemp", ";#it{p}_{T}^{jet} [GeV];(1 / N_{jet}) (dN_{jet} / d#it{p}_{T}^{jet}) [GeV^{-1}]", 1, pTJBins[0], maxx);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.036);
    h->GetXaxis ()->SetLabelSize (0.036);
    h->GetXaxis ()->SetTitleOffset (1.5);
    h->GetYaxis ()->SetTitleSize (0.036);
    h->GetYaxis ()->SetLabelSize (0.036);
    h->GetYaxis ()->SetTitleOffset (1.5);

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    {
      h = (TH1D*) h_jet_pt_ref[2*iDType][iSamp]->Clone ("htemp");
      h->Scale (std::pow (10, 3));
      myDrawHist (h, kBlack, 1, 2);
      SaferDelete (&h);

      h = h_jet_pt_ref[2*iDType][iSamp];
   
      g = make_graph (h);
      RecenterGraph (g);
      ScaleGraph (g, nullptr, std::pow (10, 3));
      myDraw (g, colorfulColors[0], kFullCircle, 1.4, 1, 3, "P", false);
      SaferDelete (&g);
    }

    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      h = (TH1D*) h_jet_pt_ref[2*iDType][iSamp]->Clone ("htemp");
      h->Scale (std::pow (10, 2-iCent));
      myDrawHist (h, kBlack, 1, 2);
      SaferDelete (&h);

      h = h_jet_pt[2*iDType][iCent][iSamp];

      g = make_graph (h);
      ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
      RecenterGraph (g);
      myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.4, 1, 3, "P", false);
      SaferDelete (&g);

    } // end loop over iCent

    myText (0.22, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.034);
    myText (0.61, 0.89, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.034);
    myText (0.61, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.034);

    mySimpleMarkerAndBoxAndLineText (0.27, 0.255, 1.4, 1001, colorfulSystColors[0], 1.0, colorfulColors[0], kFullCircle, 1.6, "#it{pp} (#times10^{3})", 0.028);

    for (short iCent = 0; iCent < nZdcCentBins; iCent++) {
      mySimpleMarkerAndBoxAndLineText (0.27 + (iCent >= 2 ? 0.34 : 0), 0.255-((iCent+1)%3)*0.035, 1.4, 1001, colorfulSystColors[iCent+1], 1.0, colorfulColors[iCent+1], kFullCircle, 1.6, Form ("%s %i-%i%% (#times10^{%i})", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent], 2-iCent), 0.028);
    } // end loop over iCent
    mySimpleMarkerAndBoxAndLineText (0.61, 0.15, 1.4, 1001, colorfulSystColors[nZdcCentBins+1], 1.0, colorfulColors[nZdcCentBins+1], kFullCircle, 1.6, Form ("All cent. (#times10^{%i})", 2-nZdcCentBins), 0.028);
    mySimpleMarkerAndBoxAndLineText (0.27, 0.15, 1.4, 1001, kWhite, 0.0, kBlack, kDot, 0.0, "#it{pp} (#it{scaled})", 0.028);

    c->RedrawAxis();

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtSpectrum%s.pdf", workPath.Data (), iDType == 1 ? "_mc" : ""));

  }



  for (short iDType = 0; iDType < 2; iDType++) {
    const char* canvasName = Form ("c_jet_trk_pt_jet_ratio_%s", iDType == 1 ? "mc" : "data");
    TCanvas* c = new TCanvas (canvasName, "", 1200, 800);
    c->Divide (3, 2);

    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      gPad->SetLogx ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];#it{p}+Pb / #it{pp}", 1, pTJBins[0], pTJBins[nPtJBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.5, 2.0);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (h_jet_pt_ratio[iDType][iCent][0], colorfulColors[nZdcCentBins-iCent], kFullCircle, 1.0, 1, 2, false);

      if (iCent < nZdcCentBins)
        myText (0.2, 0.865, colorfulColors[nZdcCentBins-iCent], Form ("#bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.2, 0.865, colorfulColors[0], "#bf{All centralities}", 0.05);
      if (iCent == nZdcCentBins) {
        myText (0.2, 0.80, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.05);
        myText (0.2, 0.74, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.05);
      }

    } // end loop over iCent

    c->cd ();
    myText (0.065, 0.971, kBlack, "#bf{#it{ATLAS}} Internal", 0.027);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtSpectrum_RatioSummary%s.pdf", workPath.Data (), iDType == 1 ? "_mc" : ""));
  } // end loop over iDType




  {
    const char* canvasName = "c_jet_trk_pt_jet_datamc_ratio";
    TCanvas* c = new TCanvas (canvasName, "", 1300, 700);
    c->Divide (4, 2);

    const short iSamp = nSamps-1;

    {
      c->cd (7);

      gPad->SetLogx ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];Data / MC", 1, pTJBins[0], pTJBins[nPtJBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.00, 1.8);
      //h->GetYaxis ()->SetRangeUser (0.10, 0.4);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (h_jet_pt_datamc_ratio_ref[0][iSamp], colorfulColors[0], kDot, 0.0, 1, 2, false, "PL");
      myDraw (h_jet_pt_datamc_ratio_ref[1][iSamp], colorfulColors[0], kDot, 0.0, 2, 2, false, "PL");
      //myDraw (h_jet_pt_datamc_ratio_ref[0][iSamp], colorfulColors[0], kOpenCircle, 1.0, 1, 2, false);
      //myDraw (h_jet_pt_datamc_ratio_ref[1][iSamp], colorfulColors[0], kFullCircle, 1.0, 1, 2, false);
      myDraw (f_jet_pt_datamc_ratio_ref[1][iSamp], colorfulColors[0], 1, 2);

      myText (0.24, 0.865, kBlack, "#bf{#it{pp}}", 0.05);
    }

    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
      c->cd (nZdcCentBins+1-iCent);

      gPad->SetLogx ();

      TH1D* h = new TH1D ("h", ";#it{p}_{T}^{jet} [GeV];Data / MC", 1, pTJBins[0], pTJBins[nPtJBins]);
      h->GetXaxis ()->SetMoreLogLabels ();
      h->GetYaxis ()->SetRangeUser (0.00, 1.8);
      //h->GetYaxis ()->SetRangeUser (0.10, 0.4);
      h->SetBinContent (1, 1);
      h->SetLineStyle (2);
      h->SetLineWidth (2);
      h->SetLineColor (kBlack);
      h->DrawCopy ("hist ][");
      SaferDelete (&h);

      myDraw (h_jet_pt_datamc_ratio[0][iCent][iSamp], colorfulColors[iCent+1], kDot, 0.0, 1, 2, false, "PL");
      myDraw (h_jet_pt_datamc_ratio[1][iCent][iSamp], colorfulColors[iCent+1], kDot, 0.0, 1, 2, false, "PL");
      //myDraw (h_jet_pt_datamc_ratio[0][iCent][iSamp], colorfulColors[iCent+1], kOpenCircle, 1.0, 1, 2, false);
      //myDraw (h_jet_pt_datamc_ratio[1][iCent][iSamp], colorfulColors[iCent+1], kFullCircle, 1.0, 1, 2, false);
      myDraw (f_jet_pt_datamc_ratio[1][iCent][iSamp], colorfulColors[iCent+1], 1, 2);

      if (iCent < nZdcCentBins)
        myText (0.24, 0.865, kBlack, Form ("#bf{#it{p}+Pb, ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.05);
      else
        myText (0.24, 0.865, kBlack, "#bf{#it{p}+Pb, 0-100%}", 0.05);

    } // end loop over iCent

    c->cd (8);
    myText (0.1, 0.84, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.07);
    myText (0.1, 0.75, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.07);
    myText (0.1, 0.66, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.07);
    myText (0.1, 0.57, kBlack, "Pythia8 JZ0-3", 0.07);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtSpectrum_DataMC_RatioSummary.pdf", workPath.Data ()));
  }




  /*{
    const char* canvasName = "c_jet_eta_phi_ref";

    TCanvas* c = new TCanvas (canvasName, "", 880, 800);

    c->SetRightMargin (0.15);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi_ref[0][0];

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    //myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetEtaPhiSpectrum_pp.pdf", workPath.Data ()));
  }



  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    const char* canvasName = Form ("c_jet_eta_phi_iCent%i", iCent);

    TCanvas* c = new TCanvas (canvasName, "", 880, 800);
    c->SetRightMargin (0.15);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi[0][iCent][0];

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    if (iCent < nZdcCentBins)
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    else
      myText (0.22, 0.85, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, All cent.", 0.032);
    //myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetEtaPhiSpectrum_pPb_iCent%i_%s.pdf", workPath.Data (), iCent));
  }*/




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // PLOT RESULTS FROM NITERS STUDY -- SUMS OF UNCERTAINTIES VS NITERS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
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
        if (xopt == -1 && (ydum > y || i == g_jet_pt_ref_unfIterUnc[iPtJInt]->GetN () - 1))
          xopt = x;
        else if (ydum < y)
          xopt = -1;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2}}", 1, 0, nIters1DMax);
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
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, -1);
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
        if (xopt == -1 && (ydum > y || i == g_jet_pt_unfIterUnc[iPtJInt][iCent]->GetN () - 1))
          xopt = x;
        else if (ydum < y)
          xopt = -1;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2}}", 1, 0, nIters1DMax);
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
        if (xopt == -1 && (ydum > y || i == g_jet_pt_ref_unfIterHybUnc[iPtJInt]->GetN () - 1))
          xopt = x;
        else if (ydum < y)
          xopt = -1;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n}", 1, 0, nIters1DMax);
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
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, -1);
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
        if (xopt == -1 && (ydum > y || i == g_jet_pt_unfIterHybUnc[iPtJInt][iCent]->GetN () - 1))
          xopt = x;
        else if (ydum < y)
          xopt = -1;
      }

      if (doLogY) ymin *= 0.2;
      else        ymin = 0;

      TH1D* h = new TH1D ("h", ";Iterations;#sqrt{#Sigma #sigma_{n}^{2} / n}", 1, 0, nIters1DMax);
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
      const int xopt_2d = GetJetSpectraNIters (false, iPtJInt, -1);
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
        if (xopt == -1 && (ydum > y || i == g_jet_pt_unfIterRelUnc[iPtJInt][iCent]->GetN () - 1))
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




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // PLOT RESULTS FROM NITERS STUDY -- AUXILIARY PLOTS (CONVERGENCE, REFOLDING, ...)
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
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
        ScaleGraph (g, (iIter == 0 ? h_jet_pt_ref[0][0] : h_jet_pt_ref_unf_nIters[iIter-1]));
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
        ScaleGraph (g, (iIter == 0 ? h_jet_pt[0][iCent][0] : h_jet_pt_unf_nIters[iCent][iIter-1]));
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
      gPad->SetLogx();

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

        TH1D* ht = (iIter == 0 ? h_jet_pt_ref[0][0] : h_jet_pt_ref_unf_nIters[iIter-1]);
        double den = ht->Integral (ht->FindBin (minJetPt+0.01), ht->FindBin (maxJetPt-0.01));

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
      gPad->SetLogx();

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

        TH1D* ht = (iIter == 0 ? h_jet_pt[0][iCent][0] : h_jet_pt_unf_nIters[iCent][iIter-1]);
        double den = ht->Integral (ht->FindBin (minJetPt+0.01), ht->FindBin (maxJetPt-0.01));

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
      gPad->SetLogx();

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

        //TH1D* ht = (iIter == 0 ? h_jet_pt_ref[0][0] : h_jet_pt_ref_unf_nIters[iIter-1]);
        //double den = ht->Integral (ht->FindBin (minJetPt+0.01), ht->FindBin (maxJetPt-0.01));

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
      gPad->SetLogx();

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

        //TH1D* ht = (iIter == 0 ? h_jet_pt[0][iCent][0] : h_jet_pt_unf_nIters[iCent][iIter-1]);
        //double den = ht->Integral (ht->FindBin (minJetPt+0.01), ht->FindBin (maxJetPt-0.01));

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
      gPad->SetLogx();

      g = new TGAE (nIters1DMax - nIters1DMin + 1);

      int nIter1p = -1;

      double den = h_jet_pt_ref[0][0]->Integral (h_jet_pt_ref[0][0]->FindBin (minJetPt+0.01), h_jet_pt_ref[0][0]->FindBin (maxJetPt-0.01));
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
      gPad->SetLogx();

      g = new TGAE (nIters1DMax - nIters1DMin + 1);

      int nIter1p = -1;

      double den = h_jet_pt[0][iCent][0]->Integral (h_jet_pt[0][iCent][0]->FindBin (minJetPt+0.01), h_jet_pt[0][iCent][0]->FindBin (maxJetPt-0.01));
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




  {
    const char* canvasName = "c_jet_pt_eta_jer_frac_ref";

    TCanvas* c = new TCanvas (canvasName, "", 1000, 800);
    c->SetRightMargin (0.19);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_pt_eta_jer_frac_ref[1][nSamps-1];

    h2->GetXaxis ()->SetTitle ("Reco. #it{p}_{T}^{jet} [GeV]");
    h2->GetYaxis ()->SetTitle ("Reco. #it{eta}^{jet}");
    h2->GetZaxis ()->SetTitle ("Fraction passing JER cut");

    h2->GetZaxis ()->SetTitleOffset (1.4 * h2->GetZaxis ()->GetTitleOffset ());

    h2->GetZaxis ()->SetRangeUser (0.92, 1);

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    //myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtEtaJERFrac_pp.pdf", workPath.Data ()));
  }



  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    const char* canvasName = Form ("c_jet_pt_eta_jer_frac_iCent%i", iCent);

    TCanvas* c = new TCanvas (canvasName, "", 1000, 800);
    c->SetRightMargin (0.19);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_pt_eta_jer_frac[1][iCent][nSamps-1];

    h2->GetXaxis ()->SetTitle ("Reco. #it{p}_{T}^{jet} [GeV]");
    h2->GetYaxis ()->SetTitle ("Reco. #it{eta}^{jet}");
    h2->GetZaxis ()->SetTitle ("Fraction passing JER cut");

    h2->GetZaxis ()->SetTitleOffset (1.4 * h2->GetZaxis ()->GetTitleOffset ());

    h2->GetZaxis ()->SetRangeUser (0.92, 1);

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    if (iCent < nZdcCentBins)
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    else
      myText (0.22, 0.85, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, All cent.", 0.032);
    //myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtEtaJERFrac_pPb_iCent%i.pdf", workPath.Data (), iCent));
  }


  return;
}


#endif
