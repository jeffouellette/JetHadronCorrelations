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

void ProcessJets (const char* tag, const char* outFileTag) {

  TFile* inFile = nullptr;

  TH1D***  h_evt_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_jet_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  TH1D**** h_evt_counts         = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH1D**** h_jet_counts         = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);

  TH1D*** h_jet_pt_ref          = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_jet_pt_cov_ref     = Get2DArray <TH2D*> (2, nVar);

  TH1D**** h_jet_pt             = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);
  TH2D**** h2_jet_pt_cov        = Get3DArray <TH2D*> (2, nZdcCentBins, nVar);

  TH1D**** h_jet_pt_ratio       = Get3DArray <TH1D*> (2, nZdcCentBins, nVar);

  TH2D***  h2_jet_eta_phi_ref   = Get2DArray <TH2D*> (2, nVar);
  TH2D**** h2_jet_eta_phi       = Get3DArray <TH2D*> (2, nZdcCentBins, nVar);

  TGAE**  g_jet_pt_ref_syst     = Get1DArray <TGAE*> (nVar);
  TGAE*** g_jet_pt_syst         = Get2DArray <TGAE*> (nZdcCentBins, nVar);

  TGAE*** g_jet_pt_ratio_syst   = Get2DArray <TGAE*> (nZdcCentBins, nVar);

  const int nAltPtJBins = (strcmp (tag, "30GeVJets") == 0 ? 13 : 16);
  double* altPtJBins = new double[nAltPtJBins+1];
  if (strcmp (tag, "30GeVJets") == 0) {
    altPtJBins[0] = pTJBins[0];
    altPtJBins[1] = pTJBins[2];
    altPtJBins[2] = pTJBins[4];
    altPtJBins[3] = pTJBins[6];
    altPtJBins[4] = pTJBins[8];
    altPtJBins[5] = pTJBins[10];
    altPtJBins[6] = pTJBins[12];
    altPtJBins[7] = pTJBins[14];
    altPtJBins[8] = pTJBins[18];
    altPtJBins[9] = pTJBins[24];
    altPtJBins[10] = pTJBins[30];
    altPtJBins[11] = pTJBins[40];
    altPtJBins[12] = pTJBins[50];
    altPtJBins[13] = pTJBins[60];
  }
  else {
    altPtJBins[0] = pTJBins[0];
    altPtJBins[1] = pTJBins[3];
    altPtJBins[2] = pTJBins[6];
    altPtJBins[3] = pTJBins[9];
    altPtJBins[4] = pTJBins[12];
    altPtJBins[5] = pTJBins[15];
    altPtJBins[6] = pTJBins[18];
    altPtJBins[7] = pTJBins[21];
    altPtJBins[8] = pTJBins[24];
    altPtJBins[9] = pTJBins[27];
    altPtJBins[10] = pTJBins[30];
    altPtJBins[11] = pTJBins[33];
    altPtJBins[12] = pTJBins[36];
    altPtJBins[13] = pTJBins[42];
    altPtJBins[14] = pTJBins[48];
    altPtJBins[15] = pTJBins[54];
    altPtJBins[16] = pTJBins[60];
  }
  


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/ProcessJets_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  for (int iDType = 0; iDType < 2; iDType++) {

    const TString dType = (iDType == 0 ? "data" : "mc");

    for (int iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
        continue;

      {
        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s17_5TeV_hists.root", rootPath.Data (), tag, var.Data (), dType.Data ());
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts_ref[iDType][iVar]    = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s17", tag, dType.Data ()))->Clone (Form ("h_evt_counts_ref_%s_%s", dType.Data (), var.Data ()));
        h_jet_counts_ref[iDType][iVar]    = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s17", tag, dType.Data ()))->Clone (Form ("h_jet_counts_ref_%s_%s", dType.Data (), var.Data ()));

        h_jet_pt_ref[iDType][iVar]        = (TH1D*) inFile->Get (Form ("h_jet_pt_%s_%s17", tag, dType.Data ()))->Clone (Form ("h_jet_pt_ref_%s_%s", dType.Data (), var.Data ()));
        h2_jet_pt_cov_ref[iDType][iVar]   = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s_%s17", tag, dType.Data ()))->Clone (Form ("h2_jet_pt_cov_ref_%s_%s", dType.Data (), var.Data ()));

        h2_jet_eta_phi_ref[iDType][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_%s17", tag, dType.Data ()))->Clone (Form ("h2_jet_eta_phi_ref_%s_%s", dType.Data (), var.Data ()));

        inFile->Close ();

        const double nEvts = h_evt_counts_ref[iDType][iVar]->GetBinContent (1); // total number of accepted evts
        const double nJets = h_jet_counts_ref[iDType][iVar]->GetBinContent (1); // total number of accepted jets

        //h_jet_pt_ref[iDType][iVar]->Rebin (4);
        //h2_jet_pt_cov_ref[iDType][iVar]->RebinX (4);
        //h2_jet_pt_cov_ref[iDType][iVar]->RebinY (4);

        CalcUncertainties (h_jet_pt_ref[iDType][iVar], h2_jet_pt_cov_ref[iDType][iVar], h_jet_counts_ref[iDType][iVar]);
        RebinSomeBins (&(h_jet_pt_ref[iDType][iVar]), nAltPtJBins, (double*) altPtJBins, true);

      }



      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/%s16_5TeV_iCent%i_hists.root", rootPath.Data (), tag, var.Data (), dType.Data (), iCent);
        std::cout << "Reading " << inFileName.Data () << std::endl;
        inFile = new TFile (inFileName.Data (), "read");
        outFile->cd ();

        h_evt_counts[iDType][iCent][iVar]   = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_%s16", tag, dType.Data ()))->Clone (Form ("h_evt_counts_pPb_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));
        h_jet_counts[iDType][iCent][iVar]   = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_%s16", tag, dType.Data ()))->Clone (Form ("h_jet_counts_pPb_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));

        h_jet_pt[iDType][iCent][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_pt_%s_%s16", tag, dType.Data ()))->Clone (Form ("h_jet_pt_pPb_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));
        h2_jet_pt_cov[iDType][iCent][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s_%s16", tag, dType.Data ()))->Clone (Form ("h2_jet_pt_cov_pPb_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));

        h2_jet_eta_phi[iDType][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_%s16", tag, dType.Data ()))->Clone (Form ("h2_jet_eta_phi_pPb_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));

        inFile->Close ();
    
        const double nEvts = h_evt_counts[iDType][iCent][iVar]->GetBinContent (1); // total number of accepted evts
        const double nJets = h_jet_counts[iDType][iCent][iVar]->GetBinContent (2); // total number of accepted jets

        //h_jet_pt[iDType][iCent][iVar]->Rebin (4);
        //h2_jet_pt_cov[iDType][iCent][iVar]->RebinX (4);
        //h2_jet_pt_cov[iDType][iCent][iVar]->RebinY (4);
    
        CalcUncertainties (h_jet_pt[iDType][iCent][iVar], h2_jet_pt_cov[iDType][iCent][iVar], h_jet_counts[iDType][iCent][iVar]);
        RebinSomeBins (&(h_jet_pt[iDType][iCent][iVar]), nAltPtJBins, (double*) altPtJBins, true);

      } // end loop over iCent



      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h_jet_pt_ratio[iDType][iCent][iVar] = (TH1D*) h_jet_pt[iDType][iCent][iVar]->Clone (Form ("h_jet_pt_ratio_iCent%i_%s_%s", iCent, dType.Data (), var.Data ()));
        h_jet_pt_ratio[iDType][iCent][iVar]->Divide (h_jet_pt_ref[iDType][iVar]);

      } // end loop over iCent

    } // end loop over iVar

  } // end loop over iDType



  {
    g_jet_pt_ref_syst[0] = make_graph (h_jet_pt_ref[0][0]);

    ResetTGAEErrors (g_jet_pt_ref_syst[0]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      g_jet_pt_syst[iCent][0]       = make_graph (h_jet_pt[0][iCent][0]);
      g_jet_pt_ratio_syst[iCent][0] = make_graph (h_jet_pt_ratio[0][iCent][0]);

      ResetTGAEErrors (g_jet_pt_syst[iCent][0]);
      ResetTGAEErrors (g_jet_pt_ratio_syst[iCent][0]);

    } // end loop over iCent
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // CREATE GRAPHS THAT WILL STORE SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE SEPARATELY.
  // THEN CALCULATES SYSTEMATIC UNCERTAINTIES FROM EACH SOURCE BY TAKING DIFFERENCES
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iDType = 0; iDType < 2; iDType++) {

    for (int iVar = 1; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) == 0))
        continue;

      g_jet_pt_ref_syst[iVar] = new TGAE ();

      CalcSystematics (g_jet_pt_ref_syst[iVar], h_jet_pt_ref[iDType][0], h_jet_pt_ref[iDType][iVar]);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        g_jet_pt_syst[iCent][iVar] = new TGAE ();
        g_jet_pt_ratio_syst[iCent][iVar] = new TGAE ();

        CalcSystematics (g_jet_pt_syst[iCent][iVar],        h_jet_pt[iDType][iCent][0],       h_jet_pt[iDType][iCent][iVar]);
        CalcSystematics (g_jet_pt_ratio_syst[iCent][iVar],  h_jet_pt_ratio[iDType][iCent][0], h_jet_pt_ratio[iDType][iCent][iVar]);

        // allow some uncertainties to not cancel in the p+Pb / pp ratio by overwriting the current uncertainties
        if (variationsThatDontCancelInRatio.count (var) != 0) {
          ResetTGAEErrors (g_jet_pt_ratio_syst[iCent][iVar]);
          AddRelErrorsInQuadrature (g_jet_pt_ratio_syst[iCent][iVar], g_jet_pt_ref_syst[iVar], false, true);
          AddRelErrorsInQuadrature (g_jet_pt_ratio_syst[iCent][iVar], g_jet_pt_syst[iCent][iVar], false, true);
        }

      } // end loop over iCent

    } // end loop over iVar

  } // end loop over iDType


  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // SYSTEMATIC UNCERTAINTIES DERIVED IN MC MUST HAVE CENTRAL VALUES SET BY CENTRAL VALUES IN DATA
  // THE FINAL UNCERTAINTY IS ASSIGNED TO MATCH THE FRACTIONAL UNCERTAINTY IN MC
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iVar = 1; iVar < nVar; iVar++) {

    const TString var = variations[iVar];

    if (dataVariations.count (var) > 0 || mcVariations.count (var) == 0)
      continue; // skip variations already evaluated in data or that are not evaluated in MC

    for (int iDir = 0; iDir < nDir; iDir++) {

      SetCentralValuesKeepRelativeErrors (g_jet_pt_ref_syst[iVar], h_jet_pt_ref[0][0]);

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        SetCentralValuesKeepRelativeErrors (g_jet_pt_syst[iCent][iVar],       h_jet_pt[0][iCent][0]);
        SetCentralValuesKeepRelativeErrors (g_jet_pt_ratio_syst[iCent][iVar], h_jet_pt_ratio[0][iCent][0]);

      } // end loop over iCent

    } // end loop over iDir

  } // end loop over iVar



  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  // ADD SYSTEMATIC UNCERTAINTIES FROM ALL SOURCES IN QUADRATURE, STORING RESULTS IN A SINGLE GRAPH
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  for (int iVar = 1; iVar < nVar; iVar++) {

    const TString var = variations[iVar];

    if (dataVariations.count (var) == 0 && mcVariations.count (var) == 0)
      continue;

    for (int iDir = 0; iDir < nDir; iDir++) {
  
      AddErrorsInQuadrature (g_jet_pt_ref_syst[0], g_jet_pt_ref_syst[iVar], false, true);
  
      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  
        AddErrorsInQuadrature (g_jet_pt_syst[iCent][0],       g_jet_pt_syst[iCent][iVar], false, true);
        AddErrorsInQuadrature (g_jet_pt_ratio_syst[iCent][0], g_jet_pt_ratio_syst[iCent][iVar], false, true);
  
      } // end loop over iCent
  
    } // end loop over iDir

  } // end loop over iVar
  


  {
    outFile->cd ();


    for (int iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      if (dataVariations.count (var) == 0 && mcVariations.count (var) == 0)
        continue;

      g_jet_pt_ref_syst[iVar]->Write (Form ("g_jet_pt_ref_syst_%s", var.Data ()));

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        g_jet_pt_syst[iCent][iVar]->Write (Form ("g_jet_pt_syst_pPb_iCent%i_%s", iCent, var.Data ()));
        g_jet_pt_ratio_syst[iCent][iVar]->Write (Form ("g_jet_pt_ratio_syst_iCent%i_%s", iCent, var.Data ()));

      } // end loop over iCent

    } // end loop over iVar


    for (int iDType = 0; iDType < 2; iDType++) {

      for (int iVar = 0; iVar < nVar; iVar++) {

        const TString var = variations[iVar];

        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;

        h_evt_counts_ref[iDType][iVar]->Write ();
        h_jet_counts_ref[iDType][iVar]->Write ();

        h_jet_pt_ref[iDType][iVar]->Write ();

        for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

          h_evt_counts[iDType][iCent][iVar]->Write ();
          h_jet_counts[iDType][iCent][iVar]->Write ();

          h_jet_pt[iDType][iCent][iVar]->Write ();
          h_jet_pt_ratio[iDType][iCent][iVar]->Write ();

        } // end loop over iCent

      } // end loop over iVar


      h2_jet_eta_phi_ref[iDType][0]->Write ();

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

        h2_jet_eta_phi[iDType][iCent][0]->Write ();

      } // end loop over iCent

    } // end loop over iDType

  
    outFile->Close ();
  }
}


#endif
