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

  TH1D**  h_evt_counts_ref = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_counts_ref = Get1DArray <TH1D*> (nVar);
  TH1D**  h_evt_counts_ref_bkg = Get1DArray <TH1D*> (nVar);
  TH1D**  h_jet_counts_ref_bkg = Get1DArray <TH1D*> (nVar);
  TH1D*** h_evt_counts = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_counts = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_evt_counts_bkg = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH1D*** h_jet_counts_bkg = Get2DArray <TH1D*> (nZdcCentBins, nVar);

  TH1D**  h_jet_pt_ref = Get1DArray <TH1D*> (nVar);
  TH2D**  h2_jet_pt_cov_ref = Get1DArray <TH2D*> (nVar);

  TH1D*** h_jet_pt = Get2DArray <TH1D*> (nZdcCentBins, nVar);
  TH2D*** h2_jet_pt_cov = Get2DArray <TH2D*> (nZdcCentBins, nVar);

  TH1D*** h_jet_pt_ratio = Get2DArray <TH1D*> (nZdcCentBins, nVar);

  TH2D**  h2_jet_eta_phi_ref = Get1DArray <TH2D*> (nVar);
  TH2D*** h2_jet_eta_phi = Get2DArray <TH2D*> (nZdcCentBins, nVar);

  TGAE**  g_jet_pt_ref_syst = Get1DArray <TGAE*> (nVar);
  TGAE*** g_jet_pt_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);

  TGAE*** g_jet_pt_ratio_syst = Get2DArray <TGAE*> (nZdcCentBins, nVar);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/PlotJets_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  for (int iVar = 0; iVar < nVar; iVar++) {

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), tag, variations[iVar].Data (), iCent);
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_counts[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", tag))->Clone (Form ("h_jet_counts_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_pt[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%s_data16", tag))->Clone (Form ("h_jet_pt_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h2_jet_pt_cov[iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s_data16", tag))->Clone (Form ("h2_jet_pt_cov_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
      h2_jet_eta_phi[iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_data16", tag))->Clone (Form ("h2_jet_eta_phi_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));

      inFile->Close ();
  
      const double nEvts = h_evt_counts[iCent][iVar]->GetBinContent (1); // total number of accepted evts
      const double nJets = h_jet_counts[iCent][iVar]->GetBinContent (2); // total number of accepted jets

      h_jet_pt[iCent][iVar]->Rebin (2);
      h2_jet_pt_cov[iCent][iVar]->RebinX (2);
      h2_jet_pt_cov[iCent][iVar]->RebinY (2);
  
      CalcUncertainties (h_jet_pt[iCent][iVar], h2_jet_pt_cov[iCent][iVar], nEvts);

      h_jet_pt[iCent][iVar]->Scale (nEvts / nJets);
      h2_jet_pt_cov[iCent][iVar]->Scale (pow (nEvts / nJets, 2));
    } // end loop over iCent



    {
      TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data17_5TeV_hists.root", rootPath.Data (), tag, variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts_ref[iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_ref_%s", variations[iVar].Data ()));
      h_jet_counts_ref[iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", tag))->Clone (Form ("h_jet_counts_ref_%s", variations[iVar].Data ()));
      h_jet_pt_ref[iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%s_data17", tag))->Clone (Form ("h_jet_pt_ref_%s", variations[iVar].Data ()));
      h2_jet_pt_cov_ref[iVar] = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s_data17", tag))->Clone (Form ("h2_jet_pt_cov_ref_%s", variations[iVar].Data ()));
      h2_jet_eta_phi_ref[iVar] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_data17", tag))->Clone (Form ("h2_jet_eta_phi_ref_%s", variations[iVar].Data ()));

      inFile->Close ();

      const double nEvts = h_evt_counts_ref[iVar]->GetBinContent (1); // total number of accepted evts
      const double nJets = h_jet_counts_ref[iVar]->GetBinContent (1); // total number of accepted jets

      h_jet_pt_ref[iVar]->Rebin (2);
      h2_jet_pt_cov_ref[iVar]->RebinX (2);
      h2_jet_pt_cov_ref[iVar]->RebinY (2);

      CalcUncertainties (h_jet_pt_ref[iVar], h2_jet_pt_cov_ref[iVar], nEvts);

      h_jet_pt_ref[iVar]->Scale (nEvts / nJets);
      h2_jet_pt_cov_ref[iVar]->Scale (pow (nEvts / nJets, 2));
    }



    //{
    //  inFile = new TFile (Form ("%s/Histograms/%s/MixedHists/Nominal/data16_5TeV_hists.root", rootPath.Data (), tag), "read");

    //  h_evt_counts_bkg[1] = (TH1D*) inFile->Get ("h_evt_counts_minbias_data16");
    //  h_jet_counts_bkg[1] = (TH1D*) inFile->Get ("h_jet_counts_minbias_data16");

    //  const double nEvts = h_evt_counts_bkg[1]->GetBinContent (1); // total number of accepted evts
    //  const double nJets = h_jet_counts_bkg[1]->GetBinContent (1); // total number of accepted jets
    //}


    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      h_jet_pt_ratio[iCent][iVar] = (TH1D*) h_jet_pt[iCent][iVar]->Clone (Form ("h_jet_pt_ratio_iCent%i_%s", iCent, variations[iVar].Data ()));
      h_jet_pt_ratio[iCent][iVar]->Divide (h_jet_pt_ref[iVar]);
    } // end loop over iCent
  } // end loop over iVar



  {
    g_jet_pt_ref_syst[0] = make_graph (h_jet_pt_ref[0]);

    ResetTGAEErrors (g_jet_pt_ref_syst[0]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      g_jet_pt_syst[iCent][0] = make_graph (h_jet_pt[iCent][0]);
      g_jet_pt_ratio_syst[iCent][0] = make_graph (h_jet_pt_ratio[iCent][0]);

      ResetTGAEErrors (g_jet_pt_syst[iCent][0]);
      ResetTGAEErrors (g_jet_pt_ratio_syst[iCent][0]);
    } // end loop over iCent
  }



  for (int iVar = 1; iVar < nVar; iVar++) {
    g_jet_pt_ref_syst[iVar] = new TGAE ();

    CalcSystematics (g_jet_pt_ref_syst[iVar], h_jet_pt_ref[0], h_jet_pt_ref[iVar]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      g_jet_pt_syst[iCent][iVar] = new TGAE ();

      g_jet_pt_ratio_syst[iCent][iVar] = new TGAE ();

      CalcSystematics (g_jet_pt_syst[iCent][iVar], h_jet_pt[iCent][0], h_jet_pt[iCent][iVar]);

      CalcSystematics (g_jet_pt_ratio_syst[iCent][iVar], h_jet_pt_ratio[iCent][0], h_jet_pt_ratio[iCent][iVar]);
    } // end loop over iCent
  } // end loop over iVar



/*
  // takes the maximum variation of the jet pT ES up/down variations
  {
    const int syst = 1;

    AddMaxSystematic (g_jet_pt_ref_syst[0], g_jet_pt_ref_syst[syst], g_jet_pt_ref_syst[syst+1]);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      AddMaxSystematic (g_jet_pt_syst[iCent][0], g_jet_pt_syst[iCent][syst], g_jet_pt_syst[iCent][syst+1]);

      AddMaxSystematic (g_jet_pt_ratio_ref_syst, g_jet_pt_ratio_syst[iCent][syst], g_jet_pt_ratio_syst[iCent][syst+1]);
    } // end loop over iCent
  }
*/



  {
    outFile->cd ();


    h_evt_counts_ref[0]->Write ();
    h_jet_counts_ref[0]->Write ();

    h_jet_pt_ref[0]->Write ();

    g_jet_pt_ref_syst[0]->Write ("g_jet_pt_ref_syst");

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      h_evt_counts[iCent][0]->Write ();
      h_jet_counts[iCent][0]->Write ();

      h_jet_pt[iCent][0]->Write ();
      h_jet_pt_ratio[iCent][0]->Write ();

      g_jet_pt_syst[iCent][0]->Write (Form ("g_jet_pt_syst_iCent%i", iCent));
      g_jet_pt_ratio_syst[iCent][0]->Write (Form ("g_jet_pt_ratio_syst_iCent%i", iCent));
    } // end loop over iCent


    for (int iVar = 0; iVar < nVar; iVar++) {
      h_jet_pt_ref[iVar]->Write ();

      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
        h_jet_pt[iCent][iVar]->Write ();
        h_jet_pt_ratio[iCent][iVar]->Write ();
      } // end loop over iCent
    } // end loop over iVar


    h2_jet_eta_phi_ref[0]->Write ();

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      h2_jet_eta_phi[iCent][0]->Write ();
    } // end loop over iCent

  
    outFile->Close ();
  }
}


#endif
