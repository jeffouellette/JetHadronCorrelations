#ifndef __JetHadronCorrelatorPlotJets_C__
#define __JetHadronCorrelatorPlotJets_C__

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

void ProcessJets (const char* tag, const char* outFileTag, const int iCent) {

  TFile* inFile = nullptr;

  TH1D*** h_evt_counts = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_jet_counts = Get2DArray <TH1D*> (2, nVar);
  //TH1D*** h_ljet_counts = Get2DArray <TH1D*> (2, nVar);
  //TH1D*** h_sljet_counts = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_evt_counts_bkg = Get2DArray <TH1D*> (2, nVar);
  TH1D*** h_jet_counts_bkg = Get2DArray <TH1D*> (2, nVar);
  //TH1D*** h_ljet_counts_bkg = Get2DArray <TH1D*> (2, nVar);
  //TH1D*** h_sljet_counts_bkg = Get2DArray <TH1D*> (2, nVar);

  TH1D*** h_jet_pt = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_jet_pt_cov = Get2DArray <TH2D*> (2, nVar);

  TH1D** h_jet_pt_ratio = Get1DArray <TH1D*> (nVar);

  TH2D*** h2_jet_eta_phi = Get2DArray <TH2D*> (2, nVar);


  TGAE*** g_jet_pt_syst = Get2DArray <TGAE*> (2, nVar);

  TGAE** g_jet_pt_ratio_syst = Get1DArray <TGAE*> (nVar);


  TString outFileName = outFileTag;
  outFileName.ReplaceAll (".root", "");
  outFileName = Form ("%s/Results/PlotJets_%s.root", rootPath.Data (), outFileName.Data ());
  std::cout << "Writing " << outFileName.Data () << std::endl;
  TFile* outFile = new TFile (outFileName.Data (), "recreate");


  for (int iVar = 0; iVar < nVar; iVar++) {
    {
      TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data16_5TeV_iCent%i_hists.root", rootPath.Data (), tag, variations[iVar].Data (), iCent);
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts[1][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_%s_pPb_%s", tag, variations[iVar].Data ()));
      h_jet_counts[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", tag))->Clone (Form ("h_jet_counts_%s_pPb_%s", tag, variations[iVar].Data ()));
      //h_ljet_counts[1][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_data16", tag))->Clone (Form ("h_evt_counts_%s_pPb_%s", tag, variations[iVar].Data ()));
      //h_sljet_counts[1][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_data16", tag))->Clone (Form ("h_sljet_counts_%s_pPb_%s", tag, variations[iVar].Data ()));
      h_jet_pt[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%s_data16", tag))->Clone (Form ("h_jet_pt_%s_pPb_%s", tag, variations[iVar].Data ()));
      h2_jet_pt_cov[1][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s_data16", tag))->Clone (Form ("h2_jet_pt_cov_%s_pPb_%s", tag, variations[iVar].Data ()));
      h2_jet_eta_phi[1][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_data16", tag))->Clone (Form ("h2_jet_eta_phi_%s_pPb_%s", tag, variations[iVar].Data ()));

      inFile->Close ();
  
      const double nEvts = h_evt_counts[1][iVar]->GetBinContent (1); // total number of accepted evts
      const double nJets = h_jet_counts[1][iVar]->GetBinContent (2); // total number of accepted jets

      h_jet_pt[1][iVar]->Rebin (2);
      h2_jet_pt_cov[1][iVar]->RebinX (2);
      h2_jet_pt_cov[1][iVar]->RebinY (2);
  
      CalcUncertainties (h_jet_pt[1][iVar], h2_jet_pt_cov[1][iVar], nEvts);

      h_jet_pt[1][iVar]->Scale (nEvts / nJets);
      h2_jet_pt_cov[1][iVar]->Scale (pow (nEvts / nJets, 2));
    }



    {
      TString inFileName = Form ("%s/Histograms/%s/JetsHists/%s/data17_5TeV_hists.root", rootPath.Data (), tag, variations[iVar].Data ());
      std::cout << "Reading " << inFileName.Data () << std::endl;
      inFile = new TFile (inFileName.Data (), "read");
      outFile->cd ();

      h_evt_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_%s_pp_%s", tag, variations[iVar].Data ()));
      h_jet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", tag))->Clone (Form ("h_jet_counts_%s_pp_%s", tag, variations[iVar].Data ()));
      //h_ljet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_data17", tag))->Clone (Form ("h_evt_counts_%s_pp_%s", tag, variations[iVar].Data ()));
      //h_sljet_counts[0][iVar] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_data17", tag))->Clone (Form ("h_sljet_counts_%s_pp_%s", tag, variations[iVar].Data ()));
      h_jet_pt[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%s_data17", tag))->Clone (Form ("h_jet_pt_%s_pp_%s", tag, variations[iVar].Data ()));
      h2_jet_pt_cov[0][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s_data17", tag))->Clone (Form ("h2_jet_pt_cov_%s_pp_%s", tag, variations[iVar].Data ()));
      h2_jet_eta_phi[0][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_data17", tag))->Clone (Form ("h2_jet_eta_phi_%s_pp_%s", tag, variations[iVar].Data ()));

      inFile->Close ();

      const double nEvts = h_evt_counts[0][iVar]->GetBinContent (1); // total number of accepted evts
      const double nJets = h_jet_counts[0][iVar]->GetBinContent (1); // total number of accepted jets

      h_jet_pt[0][iVar]->Rebin (2);
      h2_jet_pt_cov[0][iVar]->RebinX (2);
      h2_jet_pt_cov[0][iVar]->RebinY (2);

      CalcUncertainties (h_jet_pt[0][iVar], h2_jet_pt_cov[0][iVar], nEvts);

      h_jet_pt[0][iVar]->Scale (nEvts / nJets);
      h2_jet_pt_cov[0][iVar]->Scale (pow (nEvts / nJets, 2));

      //h_jet_pt[0][iVar]->Scale (1e-1);
      //h2_jet_pt_cov[0][iVar]->Scale (1e-2);
    }



    //{
    //  inFile = new TFile (Form ("%s/Histograms/%s/MixedHists/Nominal/data16_5TeV_hists.root", rootPath.Data (), tag), "read");

    //  h_evt_counts_bkg[1] = (TH1D*) inFile->Get ("h_evt_counts_minbias_data16");
    //  h_jet_counts_bkg[1] = (TH1D*) inFile->Get ("h_jet_counts_minbias_data16");
    //  h_ljet_counts_bkg[1] = (TH1D*) inFile->Get ("h_ljet_counts_minbias_data16");
    //  h_sljet_counts_bkg[1] = (TH1D*) inFile->Get ("h_sljet_counts_minbias_data16");

    //  const double nEvts = h_evt_counts_bkg[1]->GetBinContent (1); // total number of accepted evts
    //  const double nJets = h_jet_counts_bkg[1]->GetBinContent (1); // total number of accepted jets
    //  const double nLJets = h_ljet_counts_bkg[1]->GetBinContent (1); // total number of accepted leading jets
    //  const double nSLJets = h_sljet_counts_bkg[1]->GetBinContent (1); // total number of accepted subleading jets
    //}


    {
      h_jet_pt_ratio[iVar] = (TH1D*) h_jet_pt[1][iVar]->Clone (Form ("h_jet_pt_ratio_%s", variations[iVar].Data ()));
      h_jet_pt_ratio[iVar]->Divide (h_jet_pt[0][iVar]);
    }
  }



  {
    g_jet_pt_syst[0][0] = make_graph (h_jet_pt[0][0]);
    g_jet_pt_syst[1][0] = make_graph (h_jet_pt[1][0]);

    g_jet_pt_ratio_syst[0] = make_graph (h_jet_pt_ratio[0]);

    ResetTGAEErrors (g_jet_pt_syst[0][0]);
    ResetTGAEErrors (g_jet_pt_syst[1][0]);

    ResetTGAEErrors (g_jet_pt_ratio_syst[0]);
  }



  for (int iVar = 1; iVar < nVar; iVar++) {
    g_jet_pt_syst[0][iVar] = new TGAE ();
    g_jet_pt_syst[1][iVar] = new TGAE ();

    g_jet_pt_ratio_syst[iVar] = new TGAE ();

    CalcSystematics (g_jet_pt_syst[0][iVar], h_jet_pt[0][0], h_jet_pt[0][iVar]);
    CalcSystematics (g_jet_pt_syst[1][iVar], h_jet_pt[1][0], h_jet_pt[1][iVar]);

    CalcSystematics (g_jet_pt_ratio_syst[iVar], h_jet_pt_ratio[0], h_jet_pt_ratio[iVar]);
  }



  // takes the maximum variation of the jet pT ES up/down variations
  {
    const int syst = 1;

    AddMaxSystematic (g_jet_pt_syst[0][0], g_jet_pt_syst[0][syst], g_jet_pt_syst[0][syst+1]);
    AddMaxSystematic (g_jet_pt_syst[1][0], g_jet_pt_syst[1][syst], g_jet_pt_syst[1][syst+1]);

    AddMaxSystematic (g_jet_pt_ratio_syst[0], g_jet_pt_ratio_syst[syst], g_jet_pt_ratio_syst[syst+1]);
  }



  {
    outFile->cd ();

    h_evt_counts[0][0]->Write ("h_evt_counts_ref");
    h_jet_counts[0][0]->Write ("h_jet_counts_ref");
    //h_ljet_counts[0][0]->Write ("h_ljet_counts_ref");
    //h_sljet_counts[0][0]->Write ("h_sljet_counts_ref");
    h_evt_counts[1][0]->Write ("h_evt_counts");
    h_jet_counts[1][0]->Write ("h_jet_counts");
    //h_ljet_counts[1][0]->Write ("h_ljet_counts");
    //h_sljet_counts[1][0]->Write ("h_sljet_counts");

    h_jet_pt[0][0]->Write ("h_jet_pt_ref");
    h_jet_pt[1][0]->Write ("h_jet_pt");

    h_jet_pt_ratio[0]->Write ("h_jet_pt_ratio");

    g_jet_pt_syst[0][0]->Write ("g_jet_pt_ref_syst");
    g_jet_pt_syst[1][0]->Write ("g_jet_pt_syst");

    g_jet_pt_ratio_syst[0]->Write ("g_jet_pt_ratio_syst");


    for (int iVar = 1; iVar < nVar; iVar++) {
      h_jet_pt[0][iVar]->Write (Form ("h_jet_pt_ref_%s", variations[iVar].Data ()));
      h_jet_pt[1][iVar]->Write (Form ("h_jet_pt_%s", variations[iVar].Data ()));
  
      h_jet_pt_ratio[iVar]->Write (Form ("h_jet_pt_ratio_%s", variations[iVar].Data ()));
    }

    h2_jet_eta_phi[0][0]->Write ("h2_jet_eta_phi_ref");
    h2_jet_eta_phi[1][0]->Write ("h2_jet_eta_phi");
  
    outFile->Close ();
  }
}


#endif
