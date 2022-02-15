#ifndef __JetHadronCorrelator_Make2DJetTrkCorr_C__
#define __JetHadronCorrelator_Make2DJetTrkCorr_C__

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
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;


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



void Make2DJetTrkCorr (const char* outFileTag) {//, const int nItersMax = 20) {

  const TString var = variations[0];

  const short nSamps = 6;

  TFile* inFile = nullptr;

  TH1D**    h_evt_counts_ref              = Get1DArray <TH1D*> (nSamps);
  TH1D***   h_jet_counts_ref              = Get2DArray <TH1D*> (nSamps, nPtJBins);

  TH1D**    h_jet_pt_ref                  = Get1DArray <TH1D*> (nSamps);
  TH2D**    h2_jet_pt_cov_ref             = Get1DArray <TH2D*> (nSamps);

  TH1D****  h_jet_trk_pt_ref              = Get3DArray <TH1D*> (nSamps, nPtJBins, nDir);
  TH2D****  h2_jet_trk_pt_cov_ref         = Get3DArray <TH2D*> (nSamps, nPtJBins, nDir);

  // now the 2D pTJet-pTch histograms
  TH2D***     h2_jet_trk_pt_ref           = Get2DArray <TH2D*> (nSamps, nDir);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/Make2DJetTrkCorr_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");

  for (short iSamp = 0; iSamp < nSamps; iSamp++) {

    const TString samp = GetSamp (1, iSamp);

    TH2D* h2_cov = nullptr;

    TString inFileName = Form ("%s/Histograms/All/JetsHists/Nominal/mc17_5TeV_hists_%s.root", rootPath.Data (), samp.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName.Data (), "read");
    outFile->cd ();

    h_evt_counts_ref[iSamp] = (TH1D*) inFile->Get ("h_evt_counts_mc17")->Clone  (Form ("h_evt_counts_ref_%s",   samp.Data ()));

    h_jet_pt_ref[iSamp]     = (TH1D*) inFile->Get ("h_jet_pt_mc17")->Clone      (Form ("h_jet_pt_ref_%s",       samp.Data ()));
    h2_cov                  = (TH2D*) inFile->Get ("h2_jet_pt_cov_mc17")->Clone (Form ("h2_jet_pt_cov_ref_%s",  samp.Data ()));

    CalcUncertainties (h_jet_pt_ref[iSamp], h2_cov, h_evt_counts_ref[iSamp]);

    h2_jet_pt_cov_ref[iSamp] = h2_cov;

    h_jet_pt_ref[iSamp]->Scale (h_evt_counts_ref[iSamp]->GetBinContent (2)); // convert distribution to total number of jets by un-scaling 1/N_evt factor
    h2_jet_pt_cov_ref[iSamp]->Scale (std::pow (h_evt_counts_ref[iSamp]->GetBinContent (2), 2));

    SaferDelete (&h_evt_counts_ref[iSamp]);


    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

      h_jet_counts_ref[iSamp][iPtJ] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mc17", pTJ.Data ()))->Clone (Form ("h_jet_counts_ref_%s_%s", pTJ.Data (), samp.Data ()));

      for (short iDir = 0; iDir < nDir; iDir++) {

        const TString dir = directions[iDir];

        h_jet_trk_pt_ref[iSamp][iPtJ][iDir] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_mc17",       dir.Data (), pTJ.Data ()))->Clone (Form ("h_jet_trk_pt_%s_ref_mc_%s_%s",      dir.Data (), pTJ.Data (), samp.Data ()));
        h2_cov                              = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_mc17",  dir.Data (), pTJ.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_ref_mc_%s_%s", dir.Data (), pTJ.Data (), samp.Data ()));

        CalcUncertainties (h_jet_trk_pt_ref[iSamp][iPtJ][iDir], h2_cov, h_jet_counts_ref[iSamp][iPtJ]);

        h2_jet_trk_pt_cov_ref[iSamp][iPtJ][iDir] = h2_cov;

        h_jet_trk_pt_ref[iSamp][iPtJ][iDir]->Scale (1, "width");
        //h_jet_trk_pt_ref[iSamp][iPtJ][iDir]->Scale (h_jet_pt_ref[iSamp]->GetBinContent (iPtJ+1), "width");
        //h2_jet_trk_pt_cov_ref[iSamp][iPtJ][iDir]->Scale (std::pow (h_jet_pt_ref[iSamp]->GetBinContent (iPtJ+1), 2));

      } // end loop over iDir

      SaferDelete (&h_jet_counts_ref[iSamp][iPtJ]);

    } // end loop over iPtJ

    inFile->Close ();
    SaferDelete (&inFile);
  }






  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // INTEGRATE OVER JET PT SELECTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iSamp = 0; iSamp < nSamps; iSamp++) {

    const TString samp = GetSamp (1, iSamp);

    for (short iDir = 0; iDir < nDir; iDir++) {

      const TString dir = directions[iDir];

      h2_jet_trk_pt_ref[iSamp][iDir] = new TH2D (Form ("h2_jet_trk_pt_%s_ref_mc_%s", dir.Data (), samp.Data ()), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{jet} [GeV];(1/N_{jet}) (dN_{ch}/d#it{p}_{T}^{ch}) [GeV^{-1}]", nPtChBins, pTChBins, nPtJBins, pTJBins);
      TH2D* h2 = h2_jet_trk_pt_ref[iSamp][iDir];
      h2->Sumw2 ();

      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        TH1D* h = h_jet_trk_pt_ref[iSamp][iPtJ][iDir];

        for (short iPtCh = 0; iPtCh < nPtChBins; iPtCh++) {
          h2->SetBinContent (iPtCh+1, iPtJ+1, h->GetBinContent (iPtCh+1));
          h2->SetBinError   (iPtCh+1, iPtJ+1, h->GetBinError   (iPtCh+1));
        } // end loop over iPtCh

      } // end loop over iPtJ

    } // end loop over iDir

  } // end loop over iSamp




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // FINALLY WRITE OUT EVERYTHING TO A SINGLE ROOT FILE WITH ALL THE RESULTS.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    outFile->cd ();

    for (short iSamp = 0; iSamp < nSamps; iSamp++) {

      h_jet_pt_ref[iSamp]->Write ();
      h2_jet_pt_cov_ref[iSamp]->Write ();

      for (short iDir = 0; iDir < nDir; iDir++) {

        h2_jet_trk_pt_ref[iSamp][iDir]->Write ();

        //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        //  h_jet_trk_pt_ref[iSamp][iPtJ][iDir]->Write ();
        //  h2_jet_trk_pt_cov_ref[iSamp][iPtJ][iDir]->Write ();

        //} // end loop over iPtJ

      } // end loop over iDir

    } // end loop over iSamp

    outFile->Close ();
  }
}


int main  (int argn, char** argv) {
  assert (argn >= 2); 
  Make2DJetTrkCorr (argv[1]);
}


#endif
