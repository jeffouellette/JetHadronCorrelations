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

using namespace JetHadronCorrelations;


TLine* l = new TLine ();
TLatex* tl = new TLatex ();


void PlotJets (const char* tag, const char* inFileTag) {

  TFile* inFile = nullptr;

  TH1D*  h_evt_counts_ref = nullptr;
  TH1D*  h_jet_counts_ref = nullptr;
  TH1D** h_evt_counts = Get1DArray <TH1D*> (nZdcCentBins);
  TH1D** h_jet_counts = Get1DArray <TH1D*> (nZdcCentBins);

  TH1D*   h_jet_pt_ref = nullptr;
  TH2D*   h2_jet_pt_cov_ref = nullptr;

  TH1D**  h_jet_pt = Get1DArray <TH1D*> (nZdcCentBins);
  TH2D**  h2_jet_pt_cov = Get1DArray <TH2D*> (nZdcCentBins);

  TH1D**  h_jet_pt_ratio = Get1DArray <TH1D*> (nZdcCentBins);

  TH2D*   h2_jet_eta_phi_ref = nullptr;
  TH2D**  h2_jet_eta_phi = Get1DArray <TH2D*> (nZdcCentBins);

  TGAE*   g_jet_pt_ref_syst = nullptr;
  TGAE**  g_jet_pt_syst = Get1DArray <TGAE*> (nZdcCentBins);

  TGAE**  g_jet_pt_ratio_syst = Get1DArray <TGAE*> (nZdcCentBins);

  TH1D**  h_jet_pt_ref_syst = Get1DArray <TH1D*> (nVar);
  TH1D*** h_jet_pt_syst = Get2DArray <TH1D*> (nZdcCentBins, nVar);

  TH1D*** h_jet_pt_ratio_syst = Get2DArray <TH1D*> (nZdcCentBins, nVar);


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotJets_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    h_evt_counts_ref = (TH1D*) inFile->Get ("h_evt_counts_ref_Nominal");
    h_jet_counts_ref = (TH1D*) inFile->Get ("h_jet_counts_ref_Nominal");

    h_jet_pt_ref = (TH1D*) inFile->Get ("h_jet_pt_ref_Nominal");

    g_jet_pt_ref_syst = (TGAE*) inFile->Get ("g_jet_pt_ref_syst");

    h2_jet_eta_phi_ref = (TH2D*) inFile->Get ("h2_jet_eta_phi_ref_Nominal");

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      h_evt_counts[iCent] = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_iCent%i_Nominal", iCent));
      h_jet_counts[iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_iCent%i_Nominal", iCent));

      h_jet_pt[iCent] = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_iCent%i_Nominal", iCent));
      h_jet_pt_ratio[iCent] = (TH1D*) inFile->Get (Form ("h_jet_pt_ratio_iCent%i_Nominal", iCent));

      g_jet_pt_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_pt_syst_pPb_iCent%i", iCent));
      g_jet_pt_ratio_syst[iCent] = (TGAE*) inFile->Get (Form ("g_jet_pt_ratio_syst_iCent%i", iCent));

      h2_jet_eta_phi[iCent] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_pPb_iCent%i_Nominal", iCent));
    } // end loop over iCent


    for (int iVar = 1; iVar < nVar; iVar++) {
      h_jet_pt_ref_syst[iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_%s", variations[iVar].Data ()));
      for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
        h_jet_pt_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_iCent%i_%s", iCent, variations[iVar].Data ()));
        h_jet_pt_ratio_syst[iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_ratio_iCent%i_%s", iCent, variations[iVar].Data ()));
      } // end loop over iCent
    } // end loop over iVar
  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_pt_iCent%i", iCent);

    TCanvas* c = new TCanvas (canvasName, "", 800, 920);
    c->cd ();

    const double fPad = 600./920.;
    TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0, 1-fPad, 1, 1.0);
    TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0, 0.0, 1, 1-fPad);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);

    uPad->SetLeftMargin (0.12);
    dPad->SetLeftMargin (0.12);
    uPad->SetRightMargin (0.03);
    dPad->SetRightMargin (0.03);

    uPad->Draw ();
    dPad->Draw ();

    TH1D* h = nullptr;
    TGAE* g = nullptr;

    double ymin=1e-7;
    double ymax=1e0;

    uPad->cd ();
    uPad->SetLogx();
    uPad->SetLogy ();

    h = (TH1D*) h_jet_pt_ref->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{jet} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fPad);
    h->GetXaxis ()->SetLabelSize (0.028/fPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fPad);
    h->GetYaxis ()->SetTitle ("(1 / N_{jet}) (dN_{jet} / d#it{p}_{T}^{jet})");
    h->GetYaxis ()->SetTitleSize (0.028/fPad);
    h->GetYaxis ()->SetLabelSize (0.028/fPad);
    h->GetYaxis ()->SetTitleOffset (1.8*fPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    //myDrawSyst (g_jet_pt_ref_syst, myBlue);
    h = h_jet_pt_ref;
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullCircle, 0.8);
    SaferDelete (&g);

    //myDrawSyst (g_jet_pt_syst[iCent], myRed);
    h = h_jet_pt[iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myRed, kFullCircle, 0.8);
    SaferDelete (&g);

    myText (0.50, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.028/fPad);
    myMarkerText (0.50, 0.72, myRed, kFullCircle, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.8, 0.028/fPad, true);
    myMarkerText (0.50, 0.66, myBlue, kFullCircle, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.8, 0.028/fPad, true);
    //myText (0.50, 0.60, kBlack, strcmp (tag, "30GeVJets") == 0 ? "Minimum bias trigger" : "J50 Trigger", 0.028/fPad);


    dPad->cd ();
    dPad->SetLogx();

    ymin=0.6;
    ymax=1.4;

    h = (TH1D*) h_jet_pt_ratio[iCent]->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{jet} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/(1-fPad));
    h->GetXaxis ()->SetLabelSize (0.028/(1-fPad));
    h->GetXaxis ()->SetTitleOffset (3.2*(1-fPad));
    h->GetYaxis ()->SetTitle ("#it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/(1-fPad));
    h->GetYaxis ()->SetLabelSize (0.028/(1-fPad));
    h->GetYaxis ()->SetTitleOffset (1.8*(1-fPad));
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    //myDrawSyst (g_jet_pt_ratio_syst[iCent], myBlue);
    h = h_jet_pt_ratio[iCent];
    g = make_graph (h);
    ResetXErrors (g);
    myDraw (g, myBlue, kFullCircle, 0.8);
    SaferDelete (&g);

    c->SaveAs (Form ("%s/Plots/JetPtSpectrum_%s.pdf", workPath.Data (), tag));
  }



  {
    const char* canvasName = "c_jet_eta_phi_ref";

    TCanvas* c = new TCanvas (canvasName, "", 880, 800);

    c->SetRightMargin (0.15);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi_ref;

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetEtaPhiSpectrum_pp_%s.pdf", workPath.Data (), tag));
  }



  for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
    const char* canvasName = Form ("c_jet_eta_phi_iCent%i", iCent);

    TCanvas* c = new TCanvas (canvasName, "", 880, 800);
    c->SetRightMargin (0.15);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi[1];

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetEtaPhiSpectrum_pPb_iCent%i_%s.pdf", workPath.Data (), iCent, tag));
  }


}


#endif
