#ifndef __JetHadronCorrelatorPlotter_C__
#define __JetHadronCorrelatorPlotter_C__

#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"

#include <Utilities.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include <iostream>
#include <math.h>

using namespace HadronYieldsAnalysis;

TColor* tcolor = new TColor ();
const Color_t myBlue = (Color_t) tcolor->GetColor ( 87, 132, 198);
const Color_t myPurple = (Color_t) tcolor->GetColor (130,  10, 130);
const Color_t myRed = (Color_t) tcolor->GetColor (255,  12,  73);

TLine* l = new TLine ();
TLatex* tl = new TLatex ();


void PlotPtCh () {

  SetupDirectories ("Data");

  TFile* inFile = new TFile ("./rootFiles/JetsData/Nominal/data16_5TeV_hists.root", "read");

  TH1D** h_jet_counts = new TH1D*[2];
  TH1D** h_ljet_counts = new TH1D*[2];
  TH1D** h_sljet_counts = new TH1D*[2];
  TH1D** h_jet_counts_bkg = new TH1D*[2];
  TH1D** h_ljet_counts_bkg = new TH1D*[2];
  TH1D** h_sljet_counts_bkg = new TH1D*[2];

  TH1D** h_jet_trk_pt_ns = new TH1D*[2];
  TH2D** h2_jet_trk_pt_ns_cov = new TH2D*[2];
  TH1D** h_ljet_trk_pt_ns = new TH1D*[2];
  TH2D** h2_ljet_trk_pt_ns_cov = new TH2D*[2];
  TH1D** h_sljet_trk_pt_ns = new TH1D*[2];
  TH2D** h2_sljet_trk_pt_ns_cov = new TH2D*[2];
  TH1D** h_jet_trk_pt_as = new TH1D*[2];
  TH2D** h2_jet_trk_pt_as_cov = new TH2D*[2];
  TH1D** h_ljet_trk_pt_as = new TH1D*[2];
  TH2D** h2_ljet_trk_pt_as_cov = new TH2D*[2];
  TH1D** h_sljet_trk_pt_as = new TH1D*[2];
  TH2D** h2_sljet_trk_pt_as_cov = new TH2D*[2];
  TH1D** h_jet_trk_pt_ns_bkg = new TH1D*[2];
  TH2D** h2_jet_trk_pt_ns_cov_bkg = new TH2D*[2];
  TH1D** h_ljet_trk_pt_ns_bkg = new TH1D*[2];
  TH2D** h2_ljet_trk_pt_ns_cov_bkg = new TH2D*[2];
  TH1D** h_sljet_trk_pt_ns_bkg = new TH1D*[2];
  TH2D** h2_sljet_trk_pt_ns_cov_bkg = new TH2D*[2];
  TH1D** h_jet_trk_pt_as_bkg = new TH1D*[2];
  TH2D** h2_jet_trk_pt_as_cov_bkg = new TH2D*[2];
  TH1D** h_ljet_trk_pt_as_bkg = new TH1D*[2];
  TH2D** h2_ljet_trk_pt_as_cov_bkg = new TH2D*[2];
  TH1D** h_sljet_trk_pt_as_bkg = new TH1D*[2];
  TH2D** h2_sljet_trk_pt_as_cov_bkg = new TH2D*[2];

  TH1D** h_jet_trk_pt_ns_sig = new TH1D*[2];
  TH1D** h_ljet_trk_pt_ns_sig = new TH1D*[2];
  TH1D** h_sljet_trk_pt_ns_sig = new TH1D*[2];
  TH1D** h_jet_trk_pt_as_sig = new TH1D*[2];
  TH1D** h_ljet_trk_pt_as_sig = new TH1D*[2];
  TH1D** h_sljet_trk_pt_as_sig = new TH1D*[2];

  TH1D* h_jet_trk_pt_ns_iaa = nullptr;
  TH1D* h_ljet_trk_pt_ns_iaa = nullptr;
  TH1D* h_sljet_trk_pt_ns_iaa = nullptr;
  TH1D* h_jet_trk_pt_as_iaa = nullptr;
  TH1D* h_ljet_trk_pt_as_iaa = nullptr;
  TH1D* h_sljet_trk_pt_as_iaa = nullptr;


  {
    h_jet_counts[1] = (TH1D*) inFile->Get ("h_jet_counts_jets_data16");
    h_ljet_counts[1] = (TH1D*) inFile->Get ("h_ljet_counts_jets_data16");
    h_sljet_counts[1] = (TH1D*) inFile->Get ("h_sljet_counts_jets_data16");
    h_jet_trk_pt_ns[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_jets_data16");
    h2_jet_trk_pt_ns_cov[1] = (TH2D*) inFile->Get ("h2_jet_trk_pt_ns_cov_jets_data16");
    h_ljet_trk_pt_ns[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns_jets_data16");
    h2_ljet_trk_pt_ns_cov[1] = (TH2D*) inFile->Get ("h2_ljet_trk_pt_ns_cov_jets_data16");
    h_sljet_trk_pt_ns[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns_jets_data16");
    h2_sljet_trk_pt_ns_cov[1] = (TH2D*) inFile->Get ("h2_sljet_trk_pt_ns_cov_jets_data16");
    h_jet_trk_pt_as[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_jets_data16");
    h2_jet_trk_pt_as_cov[1] = (TH2D*) inFile->Get ("h2_jet_trk_pt_as_cov_jets_data16");
    h_ljet_trk_pt_as[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_as_jets_data16");
    h2_ljet_trk_pt_as_cov[1] = (TH2D*) inFile->Get ("h2_ljet_trk_pt_as_cov_jets_data16");
    h_sljet_trk_pt_as[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_as_jets_data16");
    h2_sljet_trk_pt_as_cov[1] = (TH2D*) inFile->Get ("h2_sljet_trk_pt_as_cov_jets_data16");

    const double nEvts = h_jet_counts[1]->GetBinContent (1); // total number of accepted jets
    const double nLEvts = h_ljet_counts[1]->GetBinContent (1); // total number of accepted leading jets
    const double nSLEvts = h_sljet_counts[1]->GetBinContent (1); // total number of accepted subleading jets

    CalcUncertainties (h_jet_trk_pt_ns[1], h2_jet_trk_pt_ns_cov[1], nEvts);
    CalcUncertainties (h_ljet_trk_pt_ns[1], h2_ljet_trk_pt_ns_cov[1], nLEvts);
    CalcUncertainties (h_sljet_trk_pt_ns[1], h2_sljet_trk_pt_ns_cov[1], nSLEvts);
    CalcUncertainties (h_jet_trk_pt_as[1], h2_jet_trk_pt_as_cov[1], nEvts);
    CalcUncertainties (h_ljet_trk_pt_as[1], h2_ljet_trk_pt_as_cov[1], nLEvts);
    CalcUncertainties (h_sljet_trk_pt_as[1], h2_sljet_trk_pt_as_cov[1], nSLEvts);
  }



  {
    inFile = new TFile ("./rootFiles/JetsData/Nominal/data17_5TeV_hists.root", "read");

    h_jet_counts[0] = (TH1D*) inFile->Get ("h_jet_counts_jets_data17");
    h_ljet_counts[0] = (TH1D*) inFile->Get ("h_ljet_counts_jets_data17");
    h_sljet_counts[0] = (TH1D*) inFile->Get ("h_sljet_counts_jets_data17");
    h_jet_trk_pt_ns[0] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_jets_data17");
    h2_jet_trk_pt_ns_cov[0] = (TH2D*) inFile->Get ("h2_jet_trk_pt_ns_cov_jets_data17");
    h_ljet_trk_pt_ns[0] = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns_jets_data17");
    h2_ljet_trk_pt_ns_cov[0] = (TH2D*) inFile->Get ("h2_ljet_trk_pt_ns_cov_jets_data17");
    h_sljet_trk_pt_ns[0] = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns_jets_data17");
    h2_sljet_trk_pt_ns_cov[0] = (TH2D*) inFile->Get ("h2_sljet_trk_pt_ns_cov_jets_data17");
    h_jet_trk_pt_as[0] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_jets_data17");
    h2_jet_trk_pt_as_cov[0] = (TH2D*) inFile->Get ("h2_jet_trk_pt_as_cov_jets_data17");
    h_ljet_trk_pt_as[0] = (TH1D*) inFile->Get ("h_ljet_trk_pt_as_jets_data17");
    h2_ljet_trk_pt_as_cov[0] = (TH2D*) inFile->Get ("h2_ljet_trk_pt_as_cov_jets_data17");
    h_sljet_trk_pt_as[0] = (TH1D*) inFile->Get ("h_sljet_trk_pt_as_jets_data17");
    h2_sljet_trk_pt_as_cov[0] = (TH2D*) inFile->Get ("h2_sljet_trk_pt_as_cov_jets_data17");

    const double nEvts = h_jet_counts[0]->GetBinContent (1); // total number of accepted jets
    const double nLEvts = h_ljet_counts[0]->GetBinContent (1); // total number of accepted leading jets
    const double nSLEvts = h_sljet_counts[0]->GetBinContent (1); // total number of accepted subleading jets

    CalcUncertainties (h_jet_trk_pt_ns[0], h2_jet_trk_pt_ns_cov[0], nEvts);
    CalcUncertainties (h_ljet_trk_pt_ns[0], h2_ljet_trk_pt_ns_cov[0], nLEvts);
    CalcUncertainties (h_sljet_trk_pt_ns[0], h2_sljet_trk_pt_ns_cov[0], nSLEvts);
    CalcUncertainties (h_jet_trk_pt_as[0], h2_jet_trk_pt_as_cov[0], nEvts);
    CalcUncertainties (h_ljet_trk_pt_as[0], h2_ljet_trk_pt_as_cov[0], nLEvts);
    CalcUncertainties (h_sljet_trk_pt_as[0], h2_sljet_trk_pt_as_cov[0], nSLEvts);
  }



  {
    inFile = new TFile ("./rootFiles/MinBiasData/Nominal/data16_5TeV_hists.root", "read");

    h_jet_counts_bkg[1] = (TH1D*) inFile->Get ("h_jet_counts_minbias_data16");
    h_ljet_counts_bkg[1] = (TH1D*) inFile->Get ("h_ljet_counts_minbias_data16");
    h_sljet_counts_bkg[1] = (TH1D*) inFile->Get ("h_sljet_counts_minbias_data16");
    h_jet_trk_pt_ns_bkg[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_ns_minbias_data16");
    h2_jet_trk_pt_ns_cov_bkg[1] = (TH2D*) inFile->Get ("h2_jet_trk_pt_ns_cov_minbias_data16");
    h_ljet_trk_pt_ns_bkg[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_ns_minbias_data16");
    h2_ljet_trk_pt_ns_cov_bkg[1] = (TH2D*) inFile->Get ("h2_ljet_trk_pt_ns_cov_minbias_data16");
    h_sljet_trk_pt_ns_bkg[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_ns_minbias_data16");
    h2_sljet_trk_pt_ns_cov_bkg[1] = (TH2D*) inFile->Get ("h2_sljet_trk_pt_ns_cov_minbias_data16");
    h_jet_trk_pt_as_bkg[1] = (TH1D*) inFile->Get ("h_jet_trk_pt_as_minbias_data16");
    h2_jet_trk_pt_as_cov_bkg[1] = (TH2D*) inFile->Get ("h2_jet_trk_pt_as_cov_minbias_data16");
    h_ljet_trk_pt_as_bkg[1] = (TH1D*) inFile->Get ("h_ljet_trk_pt_as_minbias_data16");
    h2_ljet_trk_pt_as_cov_bkg[1] = (TH2D*) inFile->Get ("h2_ljet_trk_pt_as_cov_minbias_data16");
    h_sljet_trk_pt_as_bkg[1] = (TH1D*) inFile->Get ("h_sljet_trk_pt_as_minbias_data16");
    h2_sljet_trk_pt_as_cov_bkg[1] = (TH2D*) inFile->Get ("h2_sljet_trk_pt_as_cov_minbias_data16");

    const double nEvts = h_jet_counts_bkg[1]->GetBinContent (1); // total number of accepted jets
    const double nLEvts = h_ljet_counts_bkg[1]->GetBinContent (1); // total number of accepted leading jets
    const double nSLEvts = h_sljet_counts_bkg[1]->GetBinContent (1); // total number of accepted subleading jets

    CalcUncertainties (h_jet_trk_pt_ns_bkg[1], h2_jet_trk_pt_ns_cov_bkg[1], nEvts);
    CalcUncertainties (h_ljet_trk_pt_ns_bkg[1], h2_ljet_trk_pt_ns_cov_bkg[1], nLEvts);
    CalcUncertainties (h_sljet_trk_pt_ns_bkg[1], h2_sljet_trk_pt_ns_cov_bkg[1], nSLEvts);
    CalcUncertainties (h_jet_trk_pt_as_bkg[1], h2_jet_trk_pt_as_cov_bkg[1], nEvts);
    CalcUncertainties (h_ljet_trk_pt_as_bkg[1], h2_ljet_trk_pt_as_cov_bkg[1], nLEvts);
    CalcUncertainties (h_sljet_trk_pt_as_bkg[1], h2_sljet_trk_pt_as_cov_bkg[1], nSLEvts);
  }



  {
    h_jet_trk_pt_ns_sig[0] = (TH1D*) h_jet_trk_pt_ns[0]->Clone ("h_jet_trk_pt_ns_data17_sig");
    h_jet_trk_pt_ns_sig[1] = (TH1D*) h_jet_trk_pt_ns[1]->Clone ("h_jet_trk_pt_ns_data16_sig");
    h_jet_trk_pt_ns_sig[1]->Add (h_jet_trk_pt_ns_bkg[1], -1);

    h_ljet_trk_pt_ns_sig[0] = (TH1D*) h_ljet_trk_pt_ns[0]->Clone ("h_ljet_trk_pt_ns_data17_sig");
    h_ljet_trk_pt_ns_sig[1] = (TH1D*) h_ljet_trk_pt_ns[1]->Clone ("h_ljet_trk_pt_ns_data16_sig");
    h_ljet_trk_pt_ns_sig[1]->Add (h_ljet_trk_pt_ns_bkg[1], -1);

    h_sljet_trk_pt_ns_sig[0] = (TH1D*) h_sljet_trk_pt_ns[0]->Clone ("h_sljet_trk_pt_ns_data17_sig");
    h_sljet_trk_pt_ns_sig[1] = (TH1D*) h_sljet_trk_pt_ns[1]->Clone ("h_sljet_trk_pt_ns_data16_sig");
    h_sljet_trk_pt_ns_sig[1]->Add (h_sljet_trk_pt_ns_bkg[1], -1);

    h_jet_trk_pt_as_sig[0] = (TH1D*) h_jet_trk_pt_as[0]->Clone ("h_jet_trk_pt_as_data17_sig");
    h_jet_trk_pt_as_sig[1] = (TH1D*) h_jet_trk_pt_as[1]->Clone ("h_jet_trk_pt_as_data16_sig");
    h_jet_trk_pt_as_sig[1]->Add (h_jet_trk_pt_as_bkg[1], -1);

    h_ljet_trk_pt_as_sig[0] = (TH1D*) h_ljet_trk_pt_as[0]->Clone ("h_ljet_trk_pt_as_data17_sig");
    h_ljet_trk_pt_as_sig[1] = (TH1D*) h_ljet_trk_pt_as[1]->Clone ("h_ljet_trk_pt_as_data16_sig");
    h_ljet_trk_pt_as_sig[1]->Add (h_ljet_trk_pt_as_bkg[1], -1);

    h_sljet_trk_pt_as_sig[0] = (TH1D*) h_sljet_trk_pt_as[0]->Clone ("h_sljet_trk_pt_as_data17_sig");
    h_sljet_trk_pt_as_sig[1] = (TH1D*) h_sljet_trk_pt_as[1]->Clone ("h_sljet_trk_pt_as_data16_sig");
    h_sljet_trk_pt_as_sig[1]->Add (h_sljet_trk_pt_as_bkg[1], -1);


    h_jet_trk_pt_ns_iaa = (TH1D*) h_jet_trk_pt_ns_sig[1]->Clone ("h_jet_trk_pt_ns_data16_iaa");
    h_jet_trk_pt_ns_iaa->Divide (h_jet_trk_pt_ns_sig[0]);

    h_ljet_trk_pt_ns_iaa = (TH1D*) h_ljet_trk_pt_ns_sig[1]->Clone ("h_ljet_trk_pt_ns_data16_iaa");
    h_ljet_trk_pt_ns_iaa->Divide (h_ljet_trk_pt_ns_sig[0]);

    h_sljet_trk_pt_ns_iaa = (TH1D*) h_sljet_trk_pt_ns_sig[1]->Clone ("h_sljet_trk_pt_ns_data16_iaa");
    h_sljet_trk_pt_ns_iaa->Divide (h_sljet_trk_pt_ns_sig[0]);

    h_jet_trk_pt_as_iaa = (TH1D*) h_jet_trk_pt_as_sig[1]->Clone ("h_jet_trk_pt_as_data16_iaa");
    h_jet_trk_pt_as_iaa->Divide (h_jet_trk_pt_as_sig[0]);

    h_ljet_trk_pt_as_iaa = (TH1D*) h_ljet_trk_pt_as_sig[1]->Clone ("h_ljet_trk_pt_as_data16_iaa");
    h_ljet_trk_pt_as_iaa->Divide (h_ljet_trk_pt_as_sig[0]);

    h_sljet_trk_pt_as_iaa = (TH1D*) h_sljet_trk_pt_as_sig[1]->Clone ("h_sljet_trk_pt_as_data16_iaa");
    h_sljet_trk_pt_as_iaa->Divide (h_sljet_trk_pt_as_sig[0]);
  }


  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  {
    const char* canvasName = "c_jet_trk_pt_as";
    TCanvas* c = new TCanvas (canvasName, "", 1200, 1120);
    c->cd ();
    const double fuPad = 480./1120.;
    const double fdPad = 320./1120.;
    const double fcPad = 1.0 - fuPad - fdPad;
    const double fxPad = 0.42;
    TPad* ulPad = new TPad (Form ("%s_ulPad", canvasName), "", 0.0, 1.0-fuPad, fxPad+0.1, 1.0);
    TPad* clPad = new TPad (Form ("%s_clPad", canvasName), "", 0.0, fdPad, fxPad+0.1, 1.0-fuPad);
    TPad* dlPad = new TPad (Form ("%s_dlPad", canvasName), "", 0.0, 0.0, fxPad+0.1, fdPad);
    TPad* urPad = new TPad (Form ("%s_urPad", canvasName), "", fxPad+0.1, 1.0-fuPad, 1.0, 1.0);
    TPad* crPad = new TPad (Form ("%s_crPad", canvasName), "", fxPad+0.1, fdPad, 1.0, 1.0-fuPad);
    TPad* drPad = new TPad (Form ("%s_drPad", canvasName), "", fxPad+0.1, 0.0, 1.0, fdPad);

    ulPad->SetTopMargin (0.14);
    urPad->SetTopMargin (0.14);
    ulPad->SetBottomMargin (0);
    urPad->SetBottomMargin (0);
    clPad->SetTopMargin (0);
    crPad->SetTopMargin (0);
    clPad->SetBottomMargin (0);
    crPad->SetBottomMargin (0);
    dlPad->SetTopMargin (0);
    drPad->SetTopMargin (0);
    dlPad->SetBottomMargin (0.25);
    drPad->SetBottomMargin (0.25);

    ulPad->SetLeftMargin (0.1 / 0.52);
    clPad->SetLeftMargin (0.1 / 0.52);
    dlPad->SetLeftMargin (0.1 / 0.52);
    ulPad->SetRightMargin (0);
    clPad->SetRightMargin (0);
    dlPad->SetRightMargin (0);
    urPad->SetLeftMargin (0);
    crPad->SetLeftMargin (0);
    drPad->SetLeftMargin (0);
    urPad->SetRightMargin (0.03 / 0.48);
    crPad->SetRightMargin (0.03 / 0.48);
    drPad->SetRightMargin (0.03 / 0.48);
    ulPad->Draw ();
    clPad->Draw ();
    dlPad->Draw ();
    urPad->Draw ();
    crPad->Draw ();
    drPad->Draw ();

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    ulPad->cd (); 
    ulPad->SetLogx ();
    ulPad->SetLogy ();

    float ymin = 6e-3;
    float ymax = 1.5e2;

    h = (TH1D*) h_jet_trk_pt_ns[1]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ns[0];
    g = make_graph (h);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_ns[1];
    g = make_graph (h);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_ns_bkg[1];
    g = make_graph (h);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.28, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myText (0.28, 0.72, kBlack, "2016 #it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{0-20%}", 0.020/fuPad);
    myText (0.28, 0.66, kBlack, "2017 #it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);

    tl->SetTextAlign (22);
    tl->SetTextFont (42);
    tl->SetTextSize (0.022/fuPad);
    tl->DrawLatexNDC (0.1/0.52 + 0.5*(1.-0.1/0.52), 0.94, "#Delta#phi < #pi/8 (near-side)");


    urPad->cd (); 
    urPad->SetLogx ();
    urPad->SetLogy ();

    ymin = 6e-3;
    ymax = 1.5e2;

    h = (TH1D*) h_jet_trk_pt_as[1]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_as[0];
    g = make_graph (h);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_as[1];
    g = make_graph (h);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_as_bkg[1];
    g = make_graph (h);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myMarkerText (0.14, 0.78, myRed, kFullCircle, "#it{p}+Pb jet-tagged events (#it{all jets})", 0.8, 0.020/fuPad, true);
    myMarkerText (0.14, 0.72, myBlue, kFullCircle, "#it{pp} jet-tagged events (#it{all jets})", 0.8, 0.020/fuPad, true);
    myMarkerText (0.14, 0.66, kBlack, kOpenCircle, "#it{p}+Pb mixed events", 0.8, 0.020/fuPad);

    tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


    clPad->cd (); 
    clPad->SetLogx ();
    clPad->SetLogy ();

    ymin = 6e-3;
    ymax = 1.5e2;

    h = (TH1D*) h_jet_trk_pt_ns_sig[1]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ns_sig[0];
    g = make_graph (h);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_ns_sig[1];
    g = make_graph (h);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    g->SetMarkerColor (myRed);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    //myMarkerText (0.29, 0.84, myRed, kFullCircle, "#it{p}+Pb", 0.8, 0.020/fdPad, true);
    //myMarkerText (0.29, 0.75, myBlue, kFullCircle, "#it{pp} (no bkg.)", 0.8, 0.020/fdPad, true);
    myText (0.28, 0.86, kBlack, "Jet-hadron correlations", 0.020/fcPad);
    myText (0.28, 0.78, kBlack, "#it{p}_{T}^{jet} > 60 GeV", 0.020/fcPad);
    myText (0.28, 0.70, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy ();

    ymin = 6e-3;
    ymax = 1.5e2;

    h = (TH1D*) h_jet_trk_pt_as_sig[1]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_as_sig[0];
    g = make_graph (h);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_jet_trk_pt_as_sig[1];
    g = make_graph (h);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    g->SetMarkerColor (myRed);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    dlPad->cd (); 
    dlPad->SetLogx ();

    ymin = 0.53;
    ymax = 1.17;

    h = (TH1D*) h_jet_trk_pt_ns_iaa->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_ns_iaa;
    g = make_graph (h);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = 1.17;

    h = (TH1D*) h_jet_trk_pt_as_iaa->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_trk_pt_as_iaa;
    g = make_graph (h);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    c->SaveAs ("Plots/JetTagged_HadronYields_Central_comparison_PtCh.pdf"); 
  }



  {
    const char* canvasName = "c_ljet_trk_pt";
    TCanvas* c = new TCanvas (canvasName, "", 1200, 1120);
    c->cd ();
    const double fuPad = 480./1120.;
    const double fdPad = 320./1120.;
    const double fcPad = 1.0 - fuPad - fdPad;
    const double fxPad = 0.42;
    TPad* ulPad = new TPad (Form ("%s_ulPad", canvasName), "", 0.0, 1.0-fuPad, fxPad+0.1, 1.0);
    TPad* clPad = new TPad (Form ("%s_clPad", canvasName), "", 0.0, fdPad, fxPad+0.1, 1.0-fuPad);
    TPad* dlPad = new TPad (Form ("%s_dlPad", canvasName), "", 0.0, 0.0, fxPad+0.1, fdPad);
    TPad* urPad = new TPad (Form ("%s_urPad", canvasName), "", fxPad+0.1, 1.0-fuPad, 1.0, 1.0);
    TPad* crPad = new TPad (Form ("%s_crPad", canvasName), "", fxPad+0.1, fdPad, 1.0, 1.0-fuPad);
    TPad* drPad = new TPad (Form ("%s_drPad", canvasName), "", fxPad+0.1, 0.0, 1.0, fdPad);

    ulPad->SetTopMargin (0.14);
    urPad->SetTopMargin (0.14);
    ulPad->SetBottomMargin (0);
    urPad->SetBottomMargin (0);
    clPad->SetTopMargin (0);
    crPad->SetTopMargin (0);
    clPad->SetBottomMargin (0);
    crPad->SetBottomMargin (0);
    dlPad->SetTopMargin (0);
    drPad->SetTopMargin (0);
    dlPad->SetBottomMargin (0.25);
    drPad->SetBottomMargin (0.25);

    ulPad->SetLeftMargin (0.1 / 0.52);
    clPad->SetLeftMargin (0.1 / 0.52);
    dlPad->SetLeftMargin (0.1 / 0.52);
    ulPad->SetRightMargin (0);
    clPad->SetRightMargin (0);
    dlPad->SetRightMargin (0);
    urPad->SetLeftMargin (0);
    crPad->SetLeftMargin (0);
    drPad->SetLeftMargin (0);
    urPad->SetRightMargin (0.03 / 0.48);
    crPad->SetRightMargin (0.03 / 0.48);
    drPad->SetRightMargin (0.03 / 0.48);
    ulPad->Draw ();
    clPad->Draw ();
    dlPad->Draw ();
    urPad->Draw ();
    crPad->Draw ();
    drPad->Draw ();

    TH1D* h = nullptr; 
    TGAE* g = nullptr;

    ulPad->cd (); 
    ulPad->SetLogx ();
    ulPad->SetLogy ();

    float ymin = 6e-3;
    float ymax = 1.5e2;

    h = (TH1D*) h_ljet_trk_pt_ns[1]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_ljet_trk_pt_ns[0];
    g = make_graph (h);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_ljet_trk_pt_ns[1];
    g = make_graph (h);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_ljet_trk_pt_ns_bkg[1];
    g = make_graph (h);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.28, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.022/fuPad);
    myText (0.28, 0.72, kBlack, "2016 #it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{0-20%}", 0.020/fuPad);
    myText (0.28, 0.66, kBlack, "2017 #it{pp}, #sqrt{s} = 5.02 TeV", 0.020/fuPad);

    tl->SetTextAlign (22);
    tl->SetTextFont (42);
    tl->SetTextSize (0.022/fuPad);
    tl->DrawLatexNDC (0.1/0.52 + 0.5*(1.-0.1/0.52), 0.94, "#Delta#phi < #pi/8 (near-side)");


    urPad->cd (); 
    urPad->SetLogx ();
    urPad->SetLogy ();

    ymin = 6e-3;
    ymax = 1.5e2;

    h = (TH1D*) h_ljet_trk_pt_as[1]->Clone ("h");
    h->Reset ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{jet}) (dN_{ch} / d#it{p}_{T}^{ch}) [GeV^{-1}]");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_ljet_trk_pt_as[0];
    g = make_graph (h);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_ljet_trk_pt_as[1];
    g = make_graph (h);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_ljet_trk_pt_as_bkg[1];
    g = make_graph (h);
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myMarkerText (0.14, 0.78, myRed, kFullCircle, "#it{p}+Pb jet-tagged events (#it{leading jet})", 0.8, 0.020/fuPad, true);
    myMarkerText (0.14, 0.72, myBlue, kFullCircle, "#it{pp} jet-tagged events (#it{leading jet})", 0.8, 0.020/fuPad, true);
    myMarkerText (0.14, 0.66, kBlack, kOpenCircle, "#it{p}+Pb mixed events", 0.8, 0.020/fuPad);

    tl->DrawLatexNDC (0.5*(1-0.03/0.48), 0.94, "#Delta#phi > 7#pi/8 (away-side)");


    clPad->cd (); 
    clPad->SetLogx ();
    clPad->SetLogy (); 

    ymin = 6e-3;
    ymax = 1.5e2;

    h = (TH1D*) h_ljet_trk_pt_ns_sig[1]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_ljet_trk_pt_ns_sig[0];
    g = make_graph (h);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_ljet_trk_pt_ns_sig[1];
    g = make_graph (h);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    g->SetMarkerColor (myRed);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    //myMarkerText (0.29, 0.84, myRed, kFullCircle, "#it{p}+Pb", 0.8, 0.020/fdPad, true);
    //myMarkerText (0.29, 0.75, myBlue, kFullCircle, "#it{pp} (no bkg.)", 0.8, 0.020/fdPad, true);
    myText (0.28, 0.86, kBlack, "Leading jet-hadron correlations", 0.020/fcPad);
    myText (0.28, 0.78, kBlack, "#it{p}_{T}^{jet} > 60 GeV", 0.020/fcPad);
    myText (0.28, 0.70, kBlack, "|#eta_{ch} - #it{y}_{CoM}| < 2.035", 0.020/fcPad);


    crPad->cd (); 
    crPad->SetLogx ();
    crPad->SetLogy (); 

    ymin = 6e-3;
    ymax = 1.5e2;

    h = (TH1D*) h_ljet_trk_pt_as_sig[1]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("(Sig.+Bkg.) - Bkg.");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_ljet_trk_pt_as_sig[0];
    g = make_graph (h);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    h = h_ljet_trk_pt_as_sig[1];
    g = make_graph (h);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    g->SetMarkerColor (myRed);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    dlPad->cd (); 
    dlPad->SetLogx ();

    ymin = 0.53;
    ymax = 1.17;

    h = (TH1D*) h_ljet_trk_pt_ns_iaa->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_ljet_trk_pt_ns_iaa;
    g = make_graph (h);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);


    drPad->cd (); 
    drPad->SetLogx ();

    ymin = 0.53;
    ymax = 1.17;

    h = (TH1D*) h_ljet_trk_pt_as_iaa->Clone ("h");
    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fdPad);
    h->GetXaxis ()->SetLabelSize (0.028/fdPad);
    h->GetXaxis ()->SetTitleOffset (3.8*fdPad);
    //h->GetXaxis ()->SetLabelOffset (-0.05*fdPad);
    h->GetYaxis ()->SetTitle ("#it{I}_{#it{p}A} = #it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/fdPad);
    h->GetYaxis ()->SetLabelSize (0.028/fdPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fdPad);
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_ljet_trk_pt_as_iaa;
    g = make_graph (h);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    c->SaveAs ("Plots/LeadingJetTagged_HadronYields_Central_comparison_PtCh.pdf"); 
  }

}


#endif
