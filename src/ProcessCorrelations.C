#ifndef __JetHadronCorrelator_ProcessCorrelations_C__
#define __JetHadronCorrelator_ProcessCorrelations_C__

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


void ProcessCorrelations (const char* tag, const char* outFileTag) {//, const int nItersMax = 20) {

  TFile* inFile = nullptr;

  //const int nItersMin = 1;
  //const double* nItersVals = linspace (nItersMin, nItersMax, nItersMax-nItersMin);

  //const bool useJetWgts = true;

  TH1D***   h_evt_counts_ref              = Get2DArray <TH1D*> (2, nVar);
  TH1D****  h_jet_counts_ref              = Get3DArray <TH1D*> (2, nPtJBins, nVar);
  //TH1D***   h_evt_counts_ref_bkg          = Get2DArray <TH1D*> (2, nVar);
  TH1D****  h_jet_counts_ref_bkg          = Get3DArray <TH1D*> (2, nPtJBins, nVar);
  TH1D****  h_evt_counts                  = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH1D***** h_jet_counts                  = Get4DArray <TH1D*> (2, nPtJBins, nZdcCentBins+1, nVar);
  //TH1D****  h_evt_counts_bkg              = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH1D***** h_jet_counts_bkg              = Get4DArray <TH1D*> (2, nPtJBins, nZdcCentBins+1, nVar);

  TH1D***   h_jet_pt_ref                  = Get2DArray <TH1D*> (2, nVar);
  TH2D**    h2_jet_pt_cov_ref             = Get1DArray <TH2D*> (2);
  TH1D****  h_jet_pt                      = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH2D***   h2_jet_pt_cov                 = Get2DArray <TH2D*> (2, nZdcCentBins+1);

  TH1D***** h_jet_trk_pt_ref              = Get4DArray <TH1D*> (2, nPtJBins, nDir, nVar);
  TH2D****  h2_jet_trk_pt_cov_ref         = Get3DArray <TH2D*> (2, nPtJBins, nDir);
  TH1D***** h_jet_trk_dphi_ref            = Get4DArray <TH1D*> (2, nPtJBins, nPtChSelections, nVar);
  TH2D****  h2_jet_trk_dphi_cov_ref       = Get3DArray <TH2D*> (2, nPtJBins, nPtChSelections);

  TH1D***** h_jet_trk_pt_ref_bkg          = Get4DArray <TH1D*> (2, nPtJBins, nDir, nVar);
  TH2D****  h2_jet_trk_pt_cov_ref_bkg     = Get3DArray <TH2D*> (2, nPtJBins, nDir);
  TH1D***** h_jet_trk_dphi_ref_bkg        = Get4DArray <TH1D*> (2, nPtJBins, nPtChSelections, nVar);
  TH2D****  h2_jet_trk_dphi_cov_ref_bkg   = Get3DArray <TH2D*> (2, nPtJBins, nPtChSelections);

  TH1D******  h_jet_trk_pt                = Get5DArray <TH1D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH2D*****   h2_jet_trk_pt_cov           = Get4DArray <TH2D*> (2, nPtJBins, nDir, nZdcCentBins+1);
  TH1D******  h_jet_trk_dphi              = Get5DArray <TH1D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1, nVar);
  TH2D*****   h2_jet_trk_dphi_cov         = Get4DArray <TH2D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1);

  TH1D******  h_jet_trk_pt_bkg            = Get5DArray <TH1D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH2D*****   h2_jet_trk_pt_cov_bkg       = Get4DArray <TH2D*> (2, nPtJBins, nDir, nZdcCentBins+1);
  TH1D******  h_jet_trk_dphi_bkg          = Get5DArray <TH1D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1, nVar);
  TH2D*****   h2_jet_trk_dphi_cov_bkg     = Get4DArray <TH2D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1);

  TH1D*****   h_jet_trk_pt_ref_sig        = Get4DArray <TH1D*> (2, nPtJBins, nDir, nVar);
  TH2D****    h2_jet_trk_pt_cov_ref_sig   = Get3DArray <TH2D*> (2, nPtJBins, nDir);
  TH1D*****   h_jet_trk_dphi_ref_sig      = Get4DArray <TH1D*> (2, nPtJBins, nPtChSelections, nVar);
  TH2D****    h2_jet_trk_dphi_cov_ref_sig = Get3DArray <TH2D*> (2, nPtJBins, nPtChSelections);

  TH1D******  h_jet_trk_pt_sig            = Get5DArray <TH1D*> (2, nPtJBins, nDir, nZdcCentBins+1, nVar);
  TH2D*****   h2_jet_trk_pt_cov_sig       = Get4DArray <TH2D*> (2, nPtJBins, nDir, nZdcCentBins+1);
  TH1D******  h_jet_trk_dphi_sig          = Get5DArray <TH1D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1, nVar);
  TH2D*****   h2_jet_trk_dphi_cov_sig     = Get4DArray <TH2D*> (2, nPtJBins, nPtChSelections, nZdcCentBins+1);

  // now the pTJet-integrated histograms (e.g. > 30 GeV and > 60 GeV)
  TH1D*****   h_jetInt_trk_pt_ref                 = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt                     = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_pt_ref_bkg             = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt_bkg                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_pt_ref_sig             = Get4DArray <TH1D*> (2, 2, nDir, nVar);
  TH1D******  h_jetInt_trk_pt_sig                 = Get5DArray <TH1D*> (2, 2, nDir, nZdcCentBins+1, nVar);
  TH2D****    h2_jetInt_trk_pt_cov_ref_sig        = Get3DArray <TH2D*> (2, 2, nDir);
  TH2D*****   h2_jetInt_trk_pt_cov_sig            = Get4DArray <TH2D*> (2, 2, nDir, nZdcCentBins+1);

  TH1D*****   h_jetInt_trk_dphi_ref               = Get4DArray <TH1D*> (2, 2, nPtChSelections, nVar);
  TH1D******  h_jetInt_trk_dphi                   = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_dphi_ref_bkg           = Get4DArray <TH1D*> (2, 2, nPtChSelections, nVar);
  TH1D******  h_jetInt_trk_dphi_bkg               = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);
  TH1D*****   h_jetInt_trk_dphi_ref_sig           = Get4DArray <TH1D*> (2, 2, nPtChSelections, nVar);
  TH1D******  h_jetInt_trk_dphi_sig               = Get5DArray <TH1D*> (2, 2, nPtChSelections, nZdcCentBins+1, nVar);
  TH2D****    h2_jetInt_trk_dphi_cov_ref_sig      = Get3DArray <TH2D*> (2, 2, nPtChSelections);
  TH2D*****   h2_jetInt_trk_dphi_cov_sig          = Get4DArray <TH2D*> (2, 2, nPtChSelections, nZdcCentBins+1);


  //TGAE****  g_jetInt_trk_pt_ref_syst          = Get3DArray <TGAE*> (2, nDir, nVar);
  //TGAE****  g_jetInt_trk_pt_ref_bkg_syst      = Get3DArray <TGAE*> (2, nDir, nVar);
  TGAE****  g_jetInt_trk_pt_ref_sig_syst      = Get3DArray <TGAE*> (2, nDir, nVar);

  //TGAE***** g_jetInt_trk_pt_syst              = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);
  //TGAE***** g_jetInt_trk_pt_bkg_syst          = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_pt_sig_syst          = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, nVar);

  // for propagating uncertainties thru unfolding...
  TGAE****  g_jet_trk_pt_ref_sig_syst         = Get3DArray <TGAE*> (nPtJBins, nDir, nVar);
  TGAE***** g_jet_trk_pt_sig_syst             = Get4DArray <TGAE*> (nPtJBins, nDir, nZdcCentBins+1, nVar);

  //TGAE****  g_jetInt_trk_dphi_ref_syst        = Get3DArray <TGAE*> (2, nPtChSelections, nVar);
  //TGAE****  g_jetInt_trk_dphi_ref_bkg_syst    = Get3DArray <TGAE*> (2, nPtChSelections, nVar);
  TGAE****  g_jetInt_trk_dphi_ref_sig_syst    = Get3DArray <TGAE*> (2, nPtChSelections, nVar);

  //TGAE***** g_jetInt_trk_dphi_syst            = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);
  //TGAE***** g_jetInt_trk_dphi_bkg_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);
  TGAE***** g_jetInt_trk_dphi_sig_syst        = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, nVar);


  //TGAE****  g_jetInt_trk_pt_ref_systTot       = Get3DArray <TGAE*> (2, nDir, 3);
  //TGAE****  g_jetInt_trk_pt_ref_bkg_systTot   = Get3DArray <TGAE*> (2, nDir, 3);
  TGAE****  g_jetInt_trk_pt_ref_sig_systTot   = Get3DArray <TGAE*> (2, nDir, 3);

  //TGAE***** g_jetInt_trk_pt_systTot           = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
  //TGAE***** g_jetInt_trk_pt_bkg_systTot       = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_pt_sig_systTot       = Get4DArray <TGAE*> (2, nDir, nZdcCentBins+1, 3);
 
  //TGAE****  g_jetInt_trk_dphi_ref_systTot     = Get3DArray <TGAE*> (2, nPtChSelections, 3);
  //TGAE****  g_jetInt_trk_dphi_ref_bkg_systTot = Get3DArray <TGAE*> (2, nPtChSelections, 3);
  TGAE****  g_jetInt_trk_dphi_ref_sig_systTot = Get3DArray <TGAE*> (2, nPtChSelections, 3);

  //TGAE***** g_jetInt_trk_dphi_systTot         = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, 3);
  //TGAE***** g_jetInt_trk_dphi_bkg_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, 3);
  TGAE***** g_jetInt_trk_dphi_sig_systTot     = Get4DArray <TGAE*> (2, nPtChSelections, nZdcCentBins+1, 3);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/ProcessCorrelations_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");

  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];
      const TString varBkg = GetVarBkg (var);
      const short iVarBkg = (varBkg == var ? iVar : GetVarN (varBkg));

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      const bool hasRefBkg = (iDType != 1 && variationsWithNoppBkgd.count (var) == 0);
      const bool hasBkg = (variationsWithNopPbBkgd.count (var) == 0);

      TH2D* h2_cov = nullptr;

      {
        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s17_5TeV_hists.root", rootPath.Data (), tag, var.Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts_ref[iDType][iVar]  = (TH1D*) inFile->Get (Form ("h_evt_counts_%s17",   dType.Data ()))->Clone (Form ("h_evt_counts_ref_%s_%s",   dType.Data (), var.Data ()));

        h_jet_pt_ref[iDType][iVar]      = (TH1D*) inFile->Get (Form ("h_jet_pt_%s17",       dType.Data ()))->Clone (Form ("h_jet_pt_ref_%s_%s",       dType.Data (), var.Data ()));
        h2_cov                          = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s17",  dType.Data ()))->Clone (Form ("h2_jet_pt_cov_ref_%s_%s",  dType.Data (), var.Data ()));

        CalcUncertainties (h_jet_pt_ref[iDType][iVar], h2_cov, h_evt_counts_ref[iDType][iVar]);

        if (iVar == 0)  h2_jet_pt_cov_ref[iDType] = h2_cov;
        else            SaferDelete (&h2_cov);

        h_jet_pt_ref[iDType][iVar]->Scale (h_evt_counts_ref[iDType][iVar]->GetBinContent (2)); // convert distribution to total number of jets by un-scaling 1/N_evt factor
        if (iVar == 0) h2_jet_pt_cov_ref[iDType]->Scale (std::pow (h_evt_counts_ref[iDType][iVar]->GetBinContent (2), 2));
        //UnscaleWidth (h_jet_pt_ref[iDType][iVar]);

        SaferDelete (&h_evt_counts_ref[iDType][iVar]);
        //SaferDelete (&h2_jet_pt_cov_ref[iDType]);


        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          h_jet_counts_ref[iDType][iPtJ][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s17", pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_ref_%s_%s_%s", dType.Data (), pTJ.Data (), var.Data ()));

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]  = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_%s17",       dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_ref_%s_%s_%s",       dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_cov                                      = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_%s17",  dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_ref_%s_%s_%s",  dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar], h2_cov, h_jet_counts_ref[iDType][iPtJ][iVar]);

            h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]->Scale (h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1), "width");
            if (iVar == 0)  {
              h2_jet_trk_pt_cov_ref[iDType][iPtJ][iDir] = h2_cov;
              h2_cov->Scale (std::pow (h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1), 2));
            }
            else
              SaferDelete (&h2_cov);

            //SaferDelete (&h2_jet_trk_pt_cov_ref[iDType][iPtJ][iDir]);

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            const TString ptch = pTChSelections[iPtCh].Data ();

            h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_%s17",       ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_ref_%s_%s_%s",      ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_cov                                        = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_%s_%s17",  ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_cov_%s_ref_%s_%s_%s", ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar], h2_cov, h_jet_counts_ref[iDType][iPtJ][iVar]);

            h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar]->Scale (h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1));
            if (iVar == 0) {
              h2_jet_trk_dphi_cov_ref[iDType][iPtJ][iPtCh] = h2_cov;
              h2_cov->Scale (std::pow (h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1), 2));
            }
            else
              SaferDelete (&h2_cov);

            h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar]->Scale (1., "width");

            //SaferDelete (&h2_jet_trk_dphi_cov_ref[iDType][iPtJ][iPtCh]);

          } // end loop over iPtCh

          SaferDelete (&h_jet_counts_ref[iDType][iPtJ][iVar]);

        } // end loop over iPtJ

        inFile->Close ();
        SaferDelete (&inFile);
      }



      if (hasRefBkg && iVar == iVarBkg) {
        TString inFileName = Form ("%s/Histograms/%s/MixedHists/%s/%s17_5TeV_hists.root", rootPath.Data (), tag, var.Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          h_jet_counts_ref_bkg[iDType][iPtJ][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_%s17", pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_ref_bkg_%s_%s_%s", dType.Data (), pTJ.Data (), var.Data ()));

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]  = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_mixed_%s17",       dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_ref_bkg_%s_%s_%s",       dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_cov                                          = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_mixed_%s17",  dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_ref_bkg_%s_%s_%s",  dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar], h2_cov, h_jet_counts_ref_bkg[iDType][iPtJ][iVar]);

            h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]->Scale (h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1));
            if (iVar == 0) {
              h2_jet_trk_pt_cov_ref_bkg[iDType][iPtJ][iDir] = h2_cov;
              h2_cov->Scale (std::pow (h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1), 2));
            }
            else 
              SaferDelete (&h2_cov);

            h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]->Scale (1., "width");

            //SaferDelete (&h2_jet_trk_pt_cov_ref_bkg[iDType][iPtJ][iDir]);

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            const TString ptch = pTChSelections[iPtCh].Data ();

            h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_mixed_%s17",       ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_ref_bkg_%s_%s_%s",      ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_cov                                            = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_%s_mixed_%s17",  ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_cov_%s_ref_bkg_%s_%s_%s", ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar], h2_cov, h_jet_counts_ref_bkg[iDType][iPtJ][iVar]);

            h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar]->Scale (h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1));
            if (iVar == 0) {
              h2_jet_trk_dphi_cov_ref_bkg[iDType][iPtJ][iPtCh] = h2_cov;
              h2_cov->Scale (std::pow (h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1), 2));
            }
            else
              SaferDelete (&h2_cov);

            h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar]->Scale (1., "width");

            //SaferDelete (&h2_jet_trk_dphi_cov_ref_bkg[iDType][iPtJ][iPtCh]);

          } // end loop over iPtCh

          SaferDelete (&h_jet_counts_ref_bkg[iDType][iPtJ][iVar]);

        } // end loop over iPtJ

        inFile->Close ();
        SaferDelete (&inFile);

      } else {
        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
          h_jet_counts_ref_bkg[iDType][iPtJ][iVar] = h_jet_counts_ref_bkg[iDType][iPtJ][iVarBkg];

          for (short iDir = 0; iDir < nDir; iDir++) {
            h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar] = h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVarBkg];
          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
            h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar] = h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVarBkg];
          } // end loop over iPtCh
        } // end loop over iPtJ
      }



      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s16_5TeV_%s_hists.root", rootPath.Data (), tag, var.Data (), dType.Data (), cent.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s16",   dType.Data ()))->Clone (Form ("h_evt_counts_pPb_%s_%s_%s",  cent.Data (), dType.Data (), var.Data ()));

        h_jet_pt[iDType][iCent][iVar]     = (TH1D*) inFile->Get (Form ("h_jet_pt_%s16",       dType.Data ()))->Clone (Form ("h_jet_pt_pPb_%s_%s_%s",      cent.Data (), dType.Data (), var.Data ()));
        h2_cov                            = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s16",  dType.Data ()))->Clone (Form ("h2_jet_pt_cov_pPb_%s_%s_%s", cent.Data (), dType.Data (), var.Data ()));

        CalcUncertainties (h_jet_pt[iDType][iCent][iVar], h2_cov, h_evt_counts[iDType][iCent][iVar]);

        if (iVar == 0)  h2_jet_pt_cov[iDType][iCent] = h2_cov;
        else            SaferDelete (&h2_cov);

        h_jet_pt[iDType][iCent][iVar]->Scale (h_evt_counts[iDType][iCent][iVar]->GetBinContent (2)); // convert distribution to total number of jets by un-scaling 1/N_evt factor
        if (iVar == 0) h2_jet_pt_cov[iDType][iCent]->Scale (std::pow (h_evt_counts[iDType][iCent][iVar]->GetBinContent (2), 2));
        //UnscaleWidth (h_jet_pt[iDType][iCent][iVar]);

        SaferDelete (&h_evt_counts[iDType][iCent][iVar]);
        //SaferDelete (&h2_jet_pt_cov[iDType][iCent]);

        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

          h_jet_counts[iDType][iPtJ][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s16",  pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_pPb_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

          for (short iDir = 0; iDir < nDir; iDir++) {

            const TString dir = directions[iDir];

            h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_%s16",       dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_pPb_%s_%s_%s_%s",       dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_cov                                        = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_%s16",  dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_pPb_%s_%s_%s_%s",  dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar], h2_cov, h_jet_counts[iDType][iPtJ][iCent][iVar]);

            h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]->Scale (h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1), "width");
            if (iVar == 0) {
              h2_jet_trk_pt_cov[iDType][iPtJ][iDir][iCent] = h2_cov;
              h2_cov->Scale (std::pow (h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1), 2));
            }
            else
              SaferDelete (&h2_cov);

            //SaferDelete (&h2_jet_trk_pt_cov[iDType][iPtJ][iDir][iCent]);

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            const TString ptch = pTChSelections[iPtCh].Data ();

            h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]  = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_%s16",       ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_pPb_%s_%s_%s_%s",      ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            h2_cov                                            = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_%s_%s16",  ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_cov_%s_pPb_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            CalcUncertainties (h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar], h2_cov, h_jet_counts[iDType][iPtJ][iCent][iVar]);

            h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]->Scale (h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1));
            if (iVar == 0) {
              h2_jet_trk_dphi_cov[iDType][iPtJ][iPtCh][iCent] = h2_cov;
              h2_cov->Scale (std::pow (h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1), 2));
            }
            else
              SaferDelete (&h2_cov);

            h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]->Scale (1., "width");

            //SaferDelete (&h2_jet_trk_dphi_cov[iDType][iPtJ][iPtCh][iCent]);

          } // end loop over iPtCh

          SaferDelete (&h_jet_counts[iDType][iPtJ][iCent][iVar]);

        } // end loop over iPtJ

        inFile->Close ();
        SaferDelete (&inFile);
      } // end loop over iCent



      if (hasBkg && iVar == iVarBkg) {
        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

          TString inFileName = Form ("%s/Histograms/%s/MixedHists/%s/%s16_5TeV_%s_hists.root", rootPath.Data (), tag, var.Data (), dType.Data (), cent.Data ());
          std::cout << "Reading " << inFileName.Data () << std::endl;
          inFile = new TFile (inFileName.Data (), "read");
          outFile->cd ();

          for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

            const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

            h_jet_counts_bkg[iDType][iPtJ][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_mixed_%s16",  pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_counts_pPb_bkg_%s_%s_%s_%s", cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

            for (short iDir = 0; iDir < nDir; iDir++) {

              const TString dir = directions[iDir];

              h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_trk_pt_%s_%s_mixed_%s16",       dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_pt_%s_pPb_bkg_%s_%s_%s_%s",       dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
              h2_cov                                            = (TH2D*) inFile->Get (Form ("h2_jet_trk_pt_cov_%s_%s_mixed_%s16",  dir.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_pt_cov_%s_pPb_bkg_%s_%s_%s_%s",  dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

              CalcUncertainties (h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar], h2_cov, h_jet_counts_bkg[iDType][iPtJ][iCent][iVar]);

              h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar]->Scale (h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1), "width");
              if (iVar == 0) {
                h2_jet_trk_pt_cov_bkg[iDType][iPtJ][iDir][iCent] = h2_cov;
                h2_cov->Scale (std::pow (h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1), 2));
              }
              else
                SaferDelete (&h2_cov);

              //SaferDelete (&h2_jet_trk_pt_cov_bkg[iDType][iPtJ][iDir][iCent]);

            } // end loop over iDir

            for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

              const TString ptch = pTChSelections[iPtCh].Data ();

              h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar]  = (TH1D*) inFile->Get (Form ("h_jet_trk_dphi_%s_%s_mixed_%s16",       ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h_jet_trk_dphi_%s_pPb_bkg_%s_%s_%s_%s",      ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
              h2_cov                                                = (TH2D*) inFile->Get (Form ("h2_jet_trk_dphi_cov_%s_%s_mixed_%s16",  ptch.Data (), pTJ.Data (), dType.Data ()))->Clone (Form ("h2_jet_trk_dphi_cov_%s_pPb_bkg_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));

              CalcUncertainties (h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar], h2_cov, h_jet_counts_bkg[iDType][iPtJ][iCent][iVar]);

              h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar]->Scale (h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1), "width");
              if (iVar == 0) {
                h2_jet_trk_dphi_cov_bkg[iDType][iPtJ][iPtCh][iCent] = h2_cov;
                h2_cov->Scale (std::pow (h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1), 2));
              }
              else
                SaferDelete (&h2_cov);

              //SaferDelete (&h2_jet_trk_dphi_cov_bkg[iDType][iPtJ][iPtCh][iCent]);

            } // end loop over iPtCh

            //SaferDelete (&h_jet_counts_bkg[iDType][iPtJ][iCent]);

          } // end loop over iPtJ

          inFile->Close ();
          SaferDelete (&inFile);

        } // end loop over iCent
      } else {
        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
          for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
            h_jet_counts_bkg[iDType][iPtJ][iCent][iVar] = h_jet_counts_bkg[iDType][iPtJ][iCent][iVarBkg];

            for (short iDir = 0; iDir < nDir; iDir++) {
              h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar] = h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVarBkg];
            } // end loop over iDir

            for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
              h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar] = h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVarBkg];
            } // end loop over iPtCh
          } // end loop over iPtJ
        } // end loop over iCent
      }

    } // end loop over iVar

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CALCULATE SIGNAL YIELDS BY SUBTRACTING THE BACKGROUND COMPONENT
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      const bool hasRefBkg = true;//(iDType != 1 && variationsWithNoppBkgd.count (var) == 0);
      const bool hasBkg = true;//(variationsWithNopPbBkgd.count (var) == 0);

      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar] = (TH1D*) h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]->Clone (Form ("h_jet_trk_pt_%s_ref_sig_%s_%s_%s", dir.Data (), dType.Data (), pTJ.Data (), var.Data ()));
          if (hasRefBkg)
            h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar]->Add (h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar], -1);

          if (iVar == 0) {
            h2_jet_trk_pt_cov_ref_sig[iDType][iPtJ][iDir] = (TH2D*) h2_jet_trk_pt_cov_ref[iDType][iPtJ][iDir]->Clone (Form ("h2_jet_trk_pt_cov_%s_ref_sig_%s_%s", dir.Data (), dType.Data (), pTJ.Data ()));
            if (hasRefBkg)
              h2_jet_trk_pt_cov_ref_sig[iDType][iPtJ][iDir]->Add (h2_jet_trk_pt_cov_ref_bkg[iDType][iPtJ][iDir]);
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar] = (TH1D*) h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]->Clone (Form ("h_jet_trk_pt_%s_pPb_sig_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            if (hasBkg)
              h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar]->Add (h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar], -1);

            if (iVar == 0) {
              h2_jet_trk_pt_cov_sig[iDType][iPtJ][iDir][iCent] = (TH2D*) h2_jet_trk_pt_cov[iDType][iPtJ][iDir][iCent]->Clone (Form ("h2_jet_trk_pt_cov_%s_pPb_sig_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJ.Data ()));
              if (hasBkg)
                h2_jet_trk_pt_cov_sig[iDType][iPtJ][iDir][iCent]->Add (h2_jet_trk_pt_cov_bkg[iDType][iPtJ][iDir][iCent]);
            }

          } // end loop over iCent

        } // end loop over iDir


        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          h_jet_trk_dphi_ref_sig[iDType][iPtJ][iPtCh][iVar] = (TH1D*) h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar]->Clone (Form ("h_jet_trk_dphi_%s_ref_sig_%s_%s_%s", ptch.Data (), dType.Data (), pTJ.Data (), var.Data ()));
          if (hasRefBkg)
            h_jet_trk_dphi_ref_sig[iDType][iPtJ][iPtCh][iVar]->Add (h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar], -1);

          if (iVar == 0) {
            h2_jet_trk_dphi_cov_ref_sig[iDType][iPtJ][iPtCh] = (TH2D*) h2_jet_trk_dphi_cov_ref[iDType][iPtJ][iPtCh]->Clone (Form ("h2_jet_trk_dphi_cov_%s_ref_sig_%s_%s", ptch.Data (), dType.Data (), pTJ.Data ()));
            if (hasRefBkg)
              h2_jet_trk_dphi_cov_ref_sig[iDType][iPtJ][iPtCh]->Add (h2_jet_trk_dphi_cov_ref_bkg[iDType][iPtJ][iPtCh]);
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jet_trk_dphi_sig[iDType][iPtJ][iPtCh][iCent][iVar] = (TH1D*) h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]->Clone (Form ("h_jet_trk_dphi_%s_pPb_sig_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJ.Data (), var.Data ()));
            if (hasBkg)
              h_jet_trk_dphi_sig[iDType][iPtJ][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar], -1);

            if (iVar == 0) {
              h2_jet_trk_dphi_cov_sig[iDType][iPtJ][iPtCh][iCent] = (TH2D*) h2_jet_trk_dphi_cov[iDType][iPtJ][iPtCh][iCent]->Clone (Form ("h2_jet_trk_dphi_cov_%s_pPb_sig_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJ.Data ()));
              if (hasBkg)
                h2_jet_trk_dphi_cov_sig[iDType][iPtJ][iPtCh][iCent]->Add (h2_jet_trk_dphi_cov_bkg[iDType][iPtJ][iPtCh][iCent]);
            }

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iPtJ

    } // end loop over iVar

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // INTEGRATE OVER JET PT SELECTIONS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");
      const float minJetPt = (iPtJInt == 0 ? 30. : 60.);
      const float maxJetPt = 300;

      for (short iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;

        const bool hasRefBkg = (iDType != 1 && variationsWithNoppBkgd.count (var) == 0);
        const bool hasBkg = (variationsWithNopPbBkgd.count (var) == 0);

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          {
            h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]     = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_%s_%s_%s",     dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            if (hasRefBkg)
              h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_bkg_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_ref_sig_%s_%s_%s", dir.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);

            if (iVar == 0)
              h2_jetInt_trk_pt_cov_ref_sig[iDType][iPtJInt][iDir] = new TH2D (Form ("h2_jetInt_trk_pt_cov_%s_ref_sig_%s_%s", dir.Data (), dType.Data (), pTJInt.Data ()), "", nPtChBins, pTChBins, nPtChBins, pTChBins);

            h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]->Sumw2 ();
            if (hasRefBkg)
              h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]->Sumw2 ();
            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Sumw2 ();

            double totalJets = 0;//, totalJetsUF = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              const double nJets = h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1);
              h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]->Add (h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]);
              if (hasRefBkg)
                h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]->Add (h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]);
              h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Add (h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar]);

              if (iVar == 0)
                h2_jetInt_trk_pt_cov_ref_sig[iDType][iPtJInt][iDir]->Add (h2_jet_trk_pt_cov_ref_sig[iDType][iPtJ][iDir]);

              totalJets += nJets;

            } // end loop over iPtJ

            h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]->Scale (1./totalJets);
            if (hasRefBkg)
              h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]->Scale (1./totalJets);
            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Scale (1./totalJets);

            if (iVar == 0)
              h2_jetInt_trk_pt_cov_ref_sig[iDType][iPtJInt][iDir]->Scale (1./(totalJets*totalJets));
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]     = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_%s_%s_%s_%s",     dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            if (hasBkg)
              h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_bkg_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);
            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_pt_%s_pPb_sig_%s_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nPtChBins, pTChBins);

            if (iVar == 0)
              h2_jetInt_trk_pt_cov_sig[iDType][iPtJInt][iDir][iCent] = new TH2D (Form ("h2_jetInt_trk_pt_cov_%s_pPb_sig_%s_%s_%s", dir.Data (), cent.Data (), dType.Data (), pTJInt.Data ()), "", nPtChBins, pTChBins, nPtChBins, pTChBins);

            h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]->Sumw2 ();
            if (hasBkg)
              h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]->Sumw2 ();
            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Sumw2 ();

            double totalJets = 0;//, totalJetsUF = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              const double nJets = h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1);
              h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]->Add (h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]);
              if (hasBkg)
                h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]->Add (h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar]);
              h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Add (h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar]);

              if (iVar == 0)
                h2_jetInt_trk_pt_cov_sig[iDType][iPtJInt][iDir][iCent]->Add (h2_jet_trk_pt_cov_sig[iDType][iPtJ][iDir][iCent]);

              totalJets += nJets;

            } // end loop over iPtJ

            h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]->Scale (1./totalJets);
            if (hasBkg)
              h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]->Scale (1./totalJets);
            h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Scale (1./totalJets);

            if (iVar == 0)
              h2_jetInt_trk_pt_cov_sig[iDType][iPtJInt][iDir][iCent]->Scale (1./(totalJets*totalJets));

          } // end loop over iCent

        } // end loop over iDir


        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
  
          const TString ptch = pTChSelections[iPtCh].Data ();

          {
            h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]     = new TH1D (Form ("h_jetInt_trk_dphi_%s_ref_%s_%s_%s",     ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);
            if (hasRefBkg)
              h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar] = new TH1D (Form ("h_jetInt_trk_dphi_%s_ref_bkg_%s_%s_%s", ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);
            h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar] = new TH1D (Form ("h_jetInt_trk_dphi_%s_ref_sig_%s_%s_%s", ptch.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);

            if (iVar == 0)
              h2_jetInt_trk_dphi_cov_ref_sig[iDType][iPtJInt][iPtCh] = new TH2D (Form ("h2_jetInt_trk_dphi_cov_%s_ref_sig_%s_%s", ptch.Data (), dType.Data (), pTJInt.Data ()), "", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);

            h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]->Sumw2 ();
            if (hasRefBkg)
              h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]->Sumw2 ();
            h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]->Sumw2 ();

            double totalJets = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              const double nJets = h_jet_pt_ref[iDType][iVar]->GetBinContent (iPtJ+1);
              h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]->Add (h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar]);
              if (hasRefBkg)
                h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]->Add (h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar]);
              h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]->Add (h_jet_trk_dphi_ref_sig[iDType][iPtJ][iPtCh][iVar]);

              if (iVar == 0)
                h2_jetInt_trk_dphi_cov_ref_sig[iDType][iPtJInt][iPtCh]->Add (h2_jet_trk_dphi_cov_ref_sig[iDType][iPtJ][iPtCh]);

              totalJets += nJets;

            } // end loop over iPtJ

            h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]->Scale (1./totalJets);
            if (hasRefBkg)
              h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]->Scale (1./totalJets);
            h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]->Scale (1./totalJets);

            if (iVar == 0)
              h2_jetInt_trk_dphi_cov_ref_sig[iDType][iPtJInt][iPtCh]->Scale (1./(totalJets*totalJets));
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]     = new TH1D (Form ("h_jetInt_trk_dphi_%s_pPb_%s_%s_%s_%s",     ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);
            if (hasBkg)
              h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_dphi_%s_pPb_bkg_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);
            h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar] = new TH1D (Form ("h_jetInt_trk_dphi_%s_pPb_sig_%s_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data (), var.Data ()), "", nDPhiBins, dPhiBins);

            if (iVar == 0)
              h2_jetInt_trk_dphi_cov_sig[iDType][iPtJInt][iPtCh][iCent] = new TH2D (Form ("h2_jetInt_trk_dphi_cov_%s_pPb_sig_%s_%s_%s", ptch.Data (), cent.Data (), dType.Data (), pTJInt.Data ()), "", nDPhiBins, dPhiBins, nDPhiBins, dPhiBins);

            h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]->Sumw2 ();
            if (hasBkg)
              h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]->Sumw2 ();
            h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Sumw2 ();

            double totalJets = 0;
            for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

              if (0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1]) < minJetPt || maxJetPt < 0.5 * (pTJBins[iPtJ] + pTJBins[iPtJ+1])) continue;

              const double nJets = h_jet_pt[iDType][iCent][iVar]->GetBinContent (iPtJ+1);
              h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]);
              if (hasBkg)
                h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar]);
              h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Add (h_jet_trk_dphi_sig[iDType][iPtJ][iPtCh][iCent][iVar]);

              if (iVar == 0)
                h2_jetInt_trk_dphi_cov_sig[iDType][iPtJInt][iPtCh][iCent]->Add (h2_jet_trk_dphi_cov_sig[iDType][iPtJ][iPtCh][iCent]);

              totalJets += nJets;

            } // end loop over iPtJ

            h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]->Scale (1./totalJets);
            if (hasBkg)
              h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]->Scale (1./totalJets);
            h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Scale (1./totalJets);

            if (iVar == 0)
              h2_jetInt_trk_dphi_cov_sig[iDType][iPtJInt][iPtCh][iCent]->Scale (1./(totalJets*totalJets));

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iVar

    } // end loop over iPtJInt

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE TOTAL SYSTEMATIC UNCERTAINTIES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

    for (short iDir = 0; iDir < nDir; iDir++) {

      //g_jet_trk_pt_ref_syst[iPtJ][iDir][0]     = make_graph (h_jet_trk_pt_ref[0][iPtJ][iDir][0]);
      //g_jet_trk_pt_ref_bkg_syst[iPtJ][iDir][0] = make_graph (h_jet_trk_pt_ref_bkg[0][iPtJ][iDir][0]);
      g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][0] = make_graph (h_jet_trk_pt_ref_sig[0][iPtJ][iDir][0]);

      //ResetTGAEErrors (g_jet_trk_pt_ref_syst[iPtJ][iDir][0]);
      //ResetTGAEErrors (g_jet_trk_pt_ref_bkg_syst[iPtJ][iDir][0]);
      ResetTGAEErrors (g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][0]);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        //g_jet_trk_pt_syst[iPtJ][iDir][iCent][0]        = make_graph (h_jet_trk_pt[0][iPtJ][iDir][iCent][0]);
        //g_jet_trk_pt_bkg_syst[iPtJ][iDir][iCent][0]    = make_graph (h_jet_trk_pt_bkg[0][iPtJ][iDir][iCent][0]);
        g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][0]    = make_graph (h_jet_trk_pt_sig[0][iPtJ][iDir][iCent][0]);

        //ResetTGAEErrors (g_jet_trk_pt_syst[iPtJ][iDir][iCent][0]);
        //ResetTGAEErrors (g_jet_trk_pt_bkg_syst[iPtJ][iDir][iCent][0]);
        ResetTGAEErrors (g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][0]);

      } // end loop over iCent

    } // end loop over iDir

  } // end loop over iPtJInt

  for (short iPtJInt : {0, 1}) {

    for (short iDir = 0; iDir < nDir; iDir++) {

      //g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0]     = make_graph (h_jetInt_trk_pt_ref[0][iPtJInt][iDir][0]);
      //g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0] = make_graph (h_jetInt_trk_pt_ref_bkg[0][iPtJInt][iDir][0]);
      g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0] = make_graph (h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir][0]);

      //ResetTGAEErrors (g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0]);
      //ResetTGAEErrors (g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0]);
      ResetTGAEErrors (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0]);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        //g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0]        = make_graph (h_jetInt_trk_pt[0][iPtJInt][iDir][iCent][0]);
        //g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0]    = make_graph (h_jetInt_trk_pt_bkg[0][iPtJInt][iDir][iCent][0]);
        g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0]    = make_graph (h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent][0]);

        //ResetTGAEErrors (g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0]);
        //ResetTGAEErrors (g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0]);

      } // end loop over iCent

    } // end loop over iDir

    for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

      //g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0]     = make_graph (h_jetInt_trk_dphi_ref[0][iPtJInt][iPtCh][0]);
      //g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0] = make_graph (h_jetInt_trk_dphi_ref_bkg[0][iPtJInt][iPtCh][0]);
      g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0] = make_graph (h_jetInt_trk_dphi_ref_sig[0][iPtJInt][iPtCh][0]);

      //ResetTGAEErrors (g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0]);
      //ResetTGAEErrors (g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0]);
      ResetTGAEErrors (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0]);

      for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        //g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0]      = make_graph (h_jetInt_trk_dphi[0][iPtJInt][iPtCh][iCent][0]);
        //g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0]  = make_graph (h_jetInt_trk_dphi_bkg[0][iPtJInt][iPtCh][iCent][0]);
        g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0]  = make_graph (h_jetInt_trk_dphi_sig[0][iPtJInt][iPtCh][iCent][0]);

        //ResetTGAEErrors (g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0]);
        //ResetTGAEErrors (g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0]);
        ResetTGAEErrors (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0]);

      } // end loop over iCent

    } // end loop over iPtCh

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE SEPARATELY.
  // THEN CALCULATES SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE BY TAKING DIFFERENCES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iDType = 0; iDType < 2; iDType++) {

    for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

      for (short iVar = 1; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
          continue;

        for (short iDir = 0; iDir < nDir; iDir++) {

          g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar] = new TGAE ();

          CalcSystematics (g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar], h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][0],  h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar]);

          //if (variationsToSmooth.count (var) > 0) {
          //  SmoothSystematics (g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar], "[0]+[1]*pow(log(x),-1)+[2]*pow(log(x),-2)+[3]*pow(log(x),-3)");
          //}

          // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
          if (variationsThatDontCancelInSig.count (var) != 0) {
            TGAE* g_tot = new TGAE ();
            TGAE* g_bkg = new TGAE ();

            CalcSystematics (g_tot, h_jet_trk_pt_ref[iDType][iPtJ][iDir][0],      h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]);
            CalcSystematics (g_bkg, h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][0],  h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]);

            ResetTGAEErrors (g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar]);
            AddErrorsInQuadrature (g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar], g_tot, false, false);
            AddErrorsInQuadrature (g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar], g_bkg, false, false);

            delete g_tot;
            delete g_bkg;
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar]  = new TGAE ();

            CalcSystematics (g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar], h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][0], h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar]);

            //if (variationsToSmooth.count (var) > 0) {
            //  SmoothSystematics (g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar], "[0]+[1]*pow(log(x),-1)+[2]*pow(log(x),-2)+[3]*pow(log(x),-3)");
            //}

            // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
            if (variationsThatDontCancelInSig.count (var) != 0) {
              TGAE* g_tot = new TGAE ();
              TGAE* g_bkg = new TGAE ();

              CalcSystematics (g_tot, h_jet_trk_pt[iDType][iPtJ][iDir][iCent][0],     h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]);
              CalcSystematics (g_bkg, h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][0], h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar]);

              ResetTGAEErrors (g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar]);
              AddErrorsInQuadrature (g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar], g_tot, false, false);
              AddErrorsInQuadrature (g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar], g_bkg, false, false);

              delete g_tot;
              delete g_bkg;
            }

          } // end loop over iCent

        } // end loop over iDir

      } // end loop over iVar

    } // end loop over iPtJ

    for (short iPtJInt : {0, 1}) {

      for (short iVar = 1; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
          continue;

        for (short iDir = 0; iDir < nDir; iDir++) {

          g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar] = new TGAE ();

          CalcSystematics (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][0],  h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]);

          //if (variationsToSmooth.count (var) > 0) {
          //  SmoothSystematics (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], "[0]+[1]*pow(log(x),-1)+[2]*pow(log(x),-2)+[3]*pow(log(x),-3)");
          //}

          // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
          if (variationsThatDontCancelInSig.count (var) != 0) {
            TGAE* g_tot = new TGAE ();
            TGAE* g_bkg = new TGAE ();

            CalcSystematics (g_tot, h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][0],      h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]);
            CalcSystematics (g_bkg, h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][0],  h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]);

            ResetTGAEErrors (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar]);
            AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], g_tot, false, false);
            AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], g_bkg, false, false);

            delete g_tot;
            delete g_bkg;
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar] = new TGAE ();

            CalcSystematics (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]);

            //if (variationsToSmooth.count (var) > 0) {
            //  SmoothSystematics (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], "[0]+[1]*pow(log(x),-1)+[2]*pow(log(x),-2)+[3]*pow(log(x),-3)");
            //}

            // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
            if (variationsThatDontCancelInSig.count (var) != 0) {
              TGAE* g_tot = new TGAE ();
              TGAE* g_bkg = new TGAE ();

              CalcSystematics (g_tot, h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][0],     h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]);
              CalcSystematics (g_bkg, h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][0], h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]);

              ResetTGAEErrors (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar]);
              AddErrorsInQuadrature (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], g_tot, false, false);
              AddErrorsInQuadrature (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], g_bkg, false, false);

              delete g_tot;
              delete g_bkg;
            }

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar]  = new TGAE ();

          CalcSystematics (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar],  h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][0], h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]);

          // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
          if (variationsThatDontCancelInSig.count (var) != 0) {
            TGAE* g_tot = new TGAE ();
            TGAE* g_bkg = new TGAE ();

            CalcSystematics (g_tot, h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][0],     h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]);
            CalcSystematics (g_bkg, h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][0], h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]);

            ResetTGAEErrors (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar]);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], g_tot, false, false);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], g_bkg, false, false);

            delete g_tot;
            delete g_bkg;
          }

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar] = new TGAE ();

            CalcSystematics (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][0],  h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]);

            // allow some uncertainties to not cancel in the bkgd. subtracted (signal) yield by overwriting the current uncertainties
            if (variationsThatDontCancelInSig.count (var) != 0) {
              TGAE* g_tot = new TGAE ();
              TGAE* g_bkg = new TGAE ();

              CalcSystematics (g_tot, h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][0],      h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]);
              CalcSystematics (g_bkg, h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][0],  h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]);

              ResetTGAEErrors (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar]);
              AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], g_tot, false, false);
              AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], g_bkg, false, false);

              delete g_tot;
              delete g_bkg;
            }

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iVar

    } // end loop over iPtJInt

  } // end loop over iDType




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // SYSTEMATIC UNCERTAINTIES DERIVED IN MC MUST HAVE CENTRAL VALUES SET BY CENTRAL VALUES IN DATA
  // THE FINAL UNCERTAINTY IS ASSIGNED TO MATCH THE FRACTIONAL UNCERTAINTY IN MC
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iVar = 1; iVar < nVar; iVar++) {

    const TString var = variations[iVar];

    if (dataVariations.count (var) > 0 || mcVariations.count (var) == 0)
      continue; // skip variations already evaluated in data or that are not evaluated in MC

    //for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

    //  for (short iDir = 0; iDir < nDir; iDir++) {

    //    SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar],  h_jet_trk_pt_ref_sig[0][iPtJ][iDir][0]);

    //    for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

    //      SetCentralValuesKeepRelativeErrors (g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar], h_jet_trk_pt_sig[0][iPtJ][iDir][iCent][0]);

    //    } // end loop over iCent

    //  } // end loop over iDir

    //} // end loop over iPtJ

    for (short iPtJInt : {0, 1}) {

      //for (short iDir = 0; iDir < nDir; iDir++) {

      //  //SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar],      h_jetInt_trk_pt_ref[0][iPtJInt][iDir][0]);
      //  //SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar],  h_jetInt_trk_pt_ref_bkg[0][iPtJInt][iDir][0]);
      //  SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar],  h_jetInt_trk_pt_ref_sig[0][iPtJInt][iDir][0]);

      //  for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      //    //SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar],     h_jetInt_trk_pt[0][iPtJInt][iDir][iCent][0]);
      //    //SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_bkg[0][iPtJInt][iDir][iCent][0]);
      //    SetCentralValuesKeepRelativeErrors (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], h_jetInt_trk_pt_sig[0][iPtJInt][iDir][iCent][0]);

      //  } // end loop over iCent

      //} // end loop over iDir

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        //SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar],     h_jetInt_trk_dphi_ref[0][iPtJInt][iPtCh][0]);
        //SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar], h_jetInt_trk_dphi_ref_bkg[0][iPtJInt][iPtCh][0]);
        SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], h_jetInt_trk_dphi_ref_sig[0][iPtJInt][iPtCh][0]);

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          //SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar],      h_jetInt_trk_dphi[0][iPtJInt][iPtCh][iCent][0]);
          //SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar],  h_jetInt_trk_dphi_bkg[0][iPtJInt][iPtCh][iCent][0]);
          SetCentralValuesKeepRelativeErrors (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar],  h_jetInt_trk_dphi_sig[0][iPtJInt][iPtCh][iCent][0]);

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iPtJInt
    
  } // end loop over iVar




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // THESE GRAPHS STORE SUMMARY SYSTEMATIC UNCERTAINTIES FOR EACH CATEGORY: TRACKING, JETS, & MIXING.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

      for (short iDir = 0; iDir < nDir; iDir++) {

        //g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar]     = (TGAE*) g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0]->Clone ();
        //g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0]->Clone ();
        g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar] = (TGAE*) g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0]->Clone ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          //g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar]     = (TGAE*) g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0]->Clone ();
          //g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar] = (TGAE*) g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0]->Clone ();
          g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar] = (TGAE*) g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0]->Clone ();

        } // end loop over iCent

      } // end loop over iDir

      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

        //g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar]     = (TGAE*) g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0]->Clone ();
        //g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar] = (TGAE*) g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0]->Clone ();
        g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar] = (TGAE*) g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0]->Clone ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          //g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar]     = (TGAE*) g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          //g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar] = (TGAE*) g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0]->Clone ();
          g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar] = (TGAE*) g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0]->Clone ();

        } // end loop over iCent

      } // end loop over iPtCh

    } // end loop over iTotVar

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM EACH CATEGORY INTO ONE OF THREE GRAPHS
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    for (std::vector <TString> vgroup : variationGroups) {

      std::vector <int> iVars = {};
      for (TString s : vgroup) {
        const short iVar = GetVarN (s);
        if (0 <= iVar && iVar < nVar)
          iVars.push_back (iVar);
      }

      if (iVars.size () == 0)
        continue;

      else if (iVars.size () == 1) {

        const short iVar = iVars[0];
        const short iTotVar = GetTotVarN (GetTotVar (variations[iVar]));
        std::cout << "For var " << variations[iVar] << " assigning to " << totalVariations[iTotVar] << std::endl;

        for (short iDir = 0; iDir < nDir; iDir++) {

          //AddErrorsInQuadrature (g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar],      g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar], false, true);
          //AddErrorsInQuadrature (g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar], false, true);
          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar], false, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            //AddErrorsInQuadrature (g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar],     g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar], false, true);
            //AddErrorsInQuadrature (g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar], false, true);

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          //AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar],     g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar], false, true);
          //AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar], g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar], false, true);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar], g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar], false, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            //AddErrorsInQuadrature (g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar],      g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar], false, true);
            //AddErrorsInQuadrature (g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar], false, true);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar], false, true);

          } // end loop over iCent

        } // end loop over iPtCh

      }

      else {

        const short iTotVar = GetTotVarN (GetTotVar (variations[iVars[0]]));
        for (int iVar : iVars)
          std::cout << "For var " << variations[iVar] << " assigning to " << GetTotVar (variations[iVar]) << std::endl;

        for (short iDir = 0; iDir < nDir; iDir++) {

          //AddErrorsInQuadrature (g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar],      g_jetInt_trk_pt_ref_syst[iPtJInt][iDir],     &iVars, true);
          //AddErrorsInQuadrature (g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir], &iVars, true);
          AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar],  g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir], &iVars, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            //AddErrorsInQuadrature (g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar],     g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent],      &iVars, true);
            //AddErrorsInQuadrature (g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent],  &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar], g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent],  &iVars, true);

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          //AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar],     g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh],      &iVars, true);
          //AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar], g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh],  &iVars, true);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar], g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh],  &iVars, true);

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            //AddErrorsInQuadrature (g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar],      g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent],     &iVars, true);
            //AddErrorsInQuadrature (g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent], &iVars, true);
            AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar],  g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent], &iVars, true);

          } // end loop over iCent

        } // end loop over iPtCh

      }

    } // end loop over iVar

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM ALL SOURCES IN QUADRATURE, STORING RESULTS IN A SINGLE GRAPH
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (short iPtJInt : {0, 1}) {

    for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

      for (short iDir = 0; iDir < nDir; iDir++) {
    
        //AddErrorsInQuadrature (g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][0],      g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar]);
        //AddErrorsInQuadrature (g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][0],  g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar]);
        AddErrorsInQuadrature (g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][0],  g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar]);
    
        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
          //AddErrorsInQuadrature (g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][0],     g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar]);
          //AddErrorsInQuadrature (g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][0], g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][0], g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar]);
    
        } // end loop over iCent
    
      } // end loop over iDir
    
      for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {
    
        //AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][0],     g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar]);
        //AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][0], g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar]);
        AddErrorsInQuadrature (g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][0], g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar]);
    
        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    
          //AddErrorsInQuadrature (g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][0],      g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar]);
          //AddErrorsInQuadrature (g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][0],  g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar]);
          AddErrorsInQuadrature (g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][0],  g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar]);
    
        } // end loop over iCent
    
      } // end loop over iPtCh

    } // end loop over iTotVar

  } // end loop over iPtJInt




  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // FINALLY WRITE OUT EVERYTHING TO A SINGLE ROOT FILE WITH ALL THE RESULTS.
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  {
    outFile->cd ();

    for (short iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if (dataVariations.count (var) == 0 && mcVariations.count (var) == 0)
        continue;

      for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        const TString pTJ = Form ("%g-%gGeVJets", pTJBins[iPtJ], pTJBins[iPtJ+1]);

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          g_jet_trk_pt_ref_sig_syst[iPtJ][iDir][iVar]->Write (Form ("g_jet_trk_pt_%s_ref_sig_syst_%s_%s",  dir.Data (), pTJ.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            g_jet_trk_pt_sig_syst[iPtJ][iDir][iCent][iVar]->Write  (Form ("g_jet_trk_pt_%s_sig_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJ.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iDir

      } // end loop over iPtJ

      for (short iPtJInt : {0, 1}) {

        const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          //g_jetInt_trk_pt_ref_syst[iPtJInt][iDir][iVar]->Write     (Form ("g_jetInt_trk_pt_%s_ref_syst_%s_%s",      dir.Data (), pTJInt.Data (), var.Data ()));
          //g_jetInt_trk_pt_ref_bkg_syst[iPtJInt][iDir][iVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_bkg_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_pt_ref_sig_syst[iPtJInt][iDir][iVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_sig_syst_%s_%s",  dir.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_pt_syst[iPtJInt][iDir][iCent][iVar]->Write      (Form ("g_jetInt_trk_pt_%s_syst_%s_%s_%s",     dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            //g_jetInt_trk_pt_bkg_syst[iPtJInt][iDir][iCent][iVar]->Write  (Form ("g_jetInt_trk_pt_%s_bkg_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_pt_sig_syst[iPtJInt][iDir][iCent][iVar]->Write  (Form ("g_jetInt_trk_pt_%s_sig_syst_%s_%s_%s", dir.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          //g_jetInt_trk_dphi_ref_syst[iPtJInt][iPtCh][iVar]->Write      (Form ("g_jetInt_trk_dphi_%s_ref_syst_%s_%s",      ptch.Data (), pTJInt.Data (), var.Data ()));
          //g_jetInt_trk_dphi_ref_bkg_syst[iPtJInt][iPtCh][iVar]->Write  (Form ("g_jetInt_trk_dphi_%s_ref_bkg_syst_%s_%s",  ptch.Data (), pTJInt.Data (), var.Data ()));
          g_jetInt_trk_dphi_ref_sig_syst[iPtJInt][iPtCh][iVar]->Write  (Form ("g_jetInt_trk_dphi_%s_ref_sig_syst_%s_%s",  ptch.Data (), pTJInt.Data (), var.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_dphi_syst[iPtJInt][iPtCh][iCent][iVar]->Write     (Form ("g_jetInt_trk_dphi_%s_syst_%s_%s_%s",     ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            //g_jetInt_trk_dphi_bkg_syst[iPtJInt][iPtCh][iCent][iVar]->Write (Form ("g_jetInt_trk_dphi_%s_bkg_syst_%s_%s_%s", ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));
            g_jetInt_trk_dphi_sig_syst[iPtJInt][iPtCh][iCent][iVar]->Write (Form ("g_jetInt_trk_dphi_%s_sig_syst_%s_%s_%s", ptch.Data (), cent.Data (), pTJInt.Data (), var.Data ()));

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iPtJInt

    } // end loop over iVar


    for (short iPtJInt : {0, 1}) {

      const TString pTJInt = (iPtJInt == 0 ? "30GeV" : "60GeV");

      for (short iTotVar = 0; iTotVar < 3; iTotVar++) {

        const TString totVar = totalVariations[iTotVar];

        for (short iDir = 0; iDir < nDir; iDir++) {

          const TString dir = directions[iDir];

          //g_jetInt_trk_pt_ref_systTot[iPtJInt][iDir][iTotVar]->Write     (Form ("g_jetInt_trk_pt_%s_ref_%s_systTot_%s",      dir.Data (), totVar.Data (), pTJInt.Data ()));
          //g_jetInt_trk_pt_ref_bkg_systTot[iPtJInt][iDir][iTotVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_bkg_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_pt_ref_sig_systTot[iPtJInt][iDir][iTotVar]->Write (Form ("g_jetInt_trk_pt_%s_ref_sig_%s_systTot_%s",  dir.Data (), totVar.Data (), pTJInt.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_pt_systTot[iPtJInt][iDir][iCent][iTotVar]->Write      (Form ("g_jetInt_trk_pt_%s_%s_systTot_%s_%s",     dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            //g_jetInt_trk_pt_bkg_systTot[iPtJInt][iDir][iCent][iTotVar]->Write  (Form ("g_jetInt_trk_pt_%s_bkg_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_pt_sig_systTot[iPtJInt][iDir][iCent][iTotVar]->Write  (Form ("g_jetInt_trk_pt_%s_sig_%s_systTot_%s_%s", dir.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iDir

        for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

          const TString ptch = pTChSelections[iPtCh].Data ();

          //g_jetInt_trk_dphi_ref_systTot[iPtJInt][iPtCh][iTotVar]->Write      (Form ("g_jetInt_trk_dphi_%s_ref_%s_systTot_%s",      ptch.Data (), totVar.Data (), pTJInt.Data ()));
          //g_jetInt_trk_dphi_ref_bkg_systTot[iPtJInt][iPtCh][iTotVar]->Write  (Form ("g_jetInt_trk_dphi_%s_ref_bkg_%s_systTot_%s",  ptch.Data (), totVar.Data (), pTJInt.Data ()));
          g_jetInt_trk_dphi_ref_sig_systTot[iPtJInt][iPtCh][iTotVar]->Write  (Form ("g_jetInt_trk_dphi_%s_ref_sig_%s_systTot_%s",  ptch.Data (), totVar.Data (), pTJInt.Data ()));

          for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

            const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

            //g_jetInt_trk_dphi_systTot[iPtJInt][iPtCh][iCent][iTotVar]->Write     (Form ("g_jetInt_trk_dphi_%s_%s_systTot_%s_%s",      ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            //g_jetInt_trk_dphi_bkg_systTot[iPtJInt][iPtCh][iCent][iTotVar]->Write (Form ("g_jetInt_trk_dphi_%s_bkg_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));
            g_jetInt_trk_dphi_sig_systTot[iPtJInt][iPtCh][iCent][iTotVar]->Write (Form ("g_jetInt_trk_dphi_%s_sig_%s_systTot_%s_%s",  ptch.Data (), totVar.Data (), cent.Data (), pTJInt.Data ()));

          } // end loop over iCent

        } // end loop over iPtCh

      } // end loop over iTotVar

    } // end loop over iPtJInt



    for (short iDType = 0; iDType < 2; iDType++) {

      for (short iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;

        h_jet_pt_ref[iDType][iVar]->Write ();
        if (iVar == 0) h2_jet_pt_cov_ref[iDType]->Write ();

        for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          h_jet_pt[iDType][iCent][iVar]->Write ();
          if (iVar == 0) h2_jet_pt_cov[iDType][iCent]->Write ();

        } // end loop over iCent

      } // end loop over iVar


      for (short iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        const bool hasRefBkg = (iDType != 1 && variationsWithNoppBkgd.count (var) == 0);
        const bool hasBkg = (variationsWithNopPbBkgd.count (var) == 0);

        for (short iPtJInt : {0, 1}) {

          if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
            continue;

          for (short iDir = 0; iDir < nDir; iDir++) {

            h_jetInt_trk_pt_ref[iDType][iPtJInt][iDir][iVar]->Write ();
            if (hasRefBkg)
              h_jetInt_trk_pt_ref_bkg[iDType][iPtJInt][iDir][iVar]->Write ();
            h_jetInt_trk_pt_ref_sig[iDType][iPtJInt][iDir][iVar]->Write ();

            if (iVar == 0)
              h2_jetInt_trk_pt_cov_ref_sig[iDType][iPtJInt][iDir]->Write ();

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              h_jetInt_trk_pt[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              if (hasBkg)
                h_jetInt_trk_pt_bkg[iDType][iPtJInt][iDir][iCent][iVar]->Write ();
              h_jetInt_trk_pt_sig[iDType][iPtJInt][iDir][iCent][iVar]->Write ();

              if (iVar == 0)
                h2_jetInt_trk_pt_cov_sig[iDType][iPtJInt][iDir][iCent]->Write ();

            } // end loop over iCent

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            h_jetInt_trk_dphi_ref[iDType][iPtJInt][iPtCh][iVar]->Write ();
            if (hasRefBkg)
              h_jetInt_trk_dphi_ref_bkg[iDType][iPtJInt][iPtCh][iVar]->Write ();
            h_jetInt_trk_dphi_ref_sig[iDType][iPtJInt][iPtCh][iVar]->Write ();

            if (iVar == 0)
              h2_jetInt_trk_dphi_cov_ref_sig[iDType][iPtJInt][iPtCh]->Write ();

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              h_jetInt_trk_dphi[iDType][iPtJInt][iPtCh][iCent][iVar]->Write ();
              if (hasBkg)
                h_jetInt_trk_dphi_bkg[iDType][iPtJInt][iPtCh][iCent][iVar]->Write ();
              h_jetInt_trk_dphi_sig[iDType][iPtJInt][iPtCh][iCent][iVar]->Write ();

              if (iVar == 0)
                h2_jetInt_trk_dphi_cov_sig[iDType][iPtJInt][iPtCh][iCent]->Write ();

            } // end loop over iCent

          } // end loop over iPtCh

        } // end loop over iPtJInt


        for (short iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
            continue;

          for (short iDir = 0; iDir < nDir; iDir++) {

            h_jet_trk_pt_ref[iDType][iPtJ][iDir][iVar]->Write ();
            if (hasRefBkg)
              h_jet_trk_pt_ref_bkg[iDType][iPtJ][iDir][iVar]->Write ();
            h_jet_trk_pt_ref_sig[iDType][iPtJ][iDir][iVar]->Write ();

            if (iVar == 0)
              h2_jet_trk_pt_cov_ref_sig[iDType][iPtJ][iDir]->Write ();

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              h_jet_trk_pt[iDType][iPtJ][iDir][iCent][iVar]->Write ();
              if (hasBkg)
                h_jet_trk_pt_bkg[iDType][iPtJ][iDir][iCent][iVar]->Write ();
              h_jet_trk_pt_sig[iDType][iPtJ][iDir][iCent][iVar]->Write ();

              if (iVar == 0)
                h2_jet_trk_pt_cov_sig[iDType][iPtJ][iDir][iCent]->Write ();

            } // end loop over iCent

          } // end loop over iDir

          for (short iPtCh = 0; iPtCh < nPtChSelections; iPtCh++) {

            h_jet_trk_dphi_ref[iDType][iPtJ][iPtCh][iVar]->Write ();
            if (hasRefBkg)
              h_jet_trk_dphi_ref_bkg[iDType][iPtJ][iPtCh][iVar]->Write ();
            h_jet_trk_dphi_ref_sig[iDType][iPtJ][iPtCh][iVar]->Write ();

            for (short iCent = 0; iCent < nZdcCentBins+1; iCent++) {

              h_jet_trk_dphi[iDType][iPtJ][iPtCh][iCent][iVar]->Write ();
              if (hasBkg)
                h_jet_trk_dphi_bkg[iDType][iPtJ][iPtCh][iCent][iVar]->Write ();
              h_jet_trk_dphi_sig[iDType][iPtJ][iPtCh][iCent][iVar]->Write ();

            } // end loop over iCent

          } // end loop over iPtCh

        } // end loop over iPtJ

      } // end loop over iVar

    } // end loop over iDType

    outFile->Close ();
  }
}


int main  (int argn, char** argv) {
  assert (argn >= 3); 
  ProcessCorrelations (argv[1], argv[2]);//, atoi (argv[3]));
}


#endif
