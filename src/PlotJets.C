#ifndef __JetHadronCorrelatorPlotJets_C__
#define __JetHadronCorrelatorPlotJets_C__

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

using namespace JetHadronCorrelations;

TColor* tcolor = new TColor ();
const Color_t myBlue = (Color_t) tcolor->GetColor ( 87, 132, 198);
const Color_t myPurple = (Color_t) tcolor->GetColor (130,  10, 130);
const Color_t myRed = (Color_t) tcolor->GetColor (255,  12,  73);

TLine* l = new TLine ();


void PlotJets (const char* tag) {

  SetupDirectories ("Data");

  TFile* inFile = new TFile (Form ("./rootFiles/%s/JetsHists/Nominal/data16_5TeV_hists.root", tag), "read");

  TH1D** h_evt_counts = new TH1D*[2];
  TH1D** h_jet_counts = new TH1D*[2];
  TH1D** h_ljet_counts = new TH1D*[2];
  TH1D** h_sljet_counts = new TH1D*[2];
  TH1D** h_evt_counts_bkg = new TH1D*[2];
  TH1D** h_jet_counts_bkg = new TH1D*[2];
  TH1D** h_ljet_counts_bkg = new TH1D*[2];
  TH1D** h_sljet_counts_bkg = new TH1D*[2];

  TH1D** h_jet_pt = new TH1D*[2];
  TH2D** h2_jet_pt_cov = new TH2D*[2];

  TH1D* h_jet_pt_ratio = nullptr;

  TH2D** h2_jet_eta_phi = new TH2D*[2];

  {
    inFile = new TFile (Form ("./rootFiles/%s/JetsHists/Nominal/data16_5TeV_hists.root", tag), "read");

    h_evt_counts[1] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data16", tag));
    h_jet_counts[1] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data16", tag));
    h_ljet_counts[1] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_data16", tag));
    h_sljet_counts[1] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_data16", tag));
    h_jet_pt[1] = (TH1D*) inFile->Get (Form ("h_jet_pt_%s_data16", tag));
    h2_jet_pt_cov[1] = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s_data16", tag));
    h2_jet_eta_phi[1] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_data16", tag));

    //const double nEvts = h_evt_counts[1]->GetBinContent (1); // total number of accepted evts
    //const double nJets = h_jet_counts[1]->GetBinContent (1); // total number of accepted jets
    //const double nLJets = h_ljet_counts[1]->GetBinContent (1); // total number of accepted leading jets
    //const double nSLJets = h_sljet_counts[1]->GetBinContent (1); // total number of accepted subleading jets

    CalcUncertainties (h_jet_pt[1], h2_jet_pt_cov[1], nEvts);
  }



  {
    inFile = new TFile (Form ("./rootFiles/%s/JetsHists/Nominal/data17_5TeV_hists.root", tag), "read");

    h_evt_counts[0] = (TH1D*) inFile->Get (Form ("h_evt_counts_%s_data17", tag));
    h_jet_counts[0] = (TH1D*) inFile->Get (Form ("h_jet_counts_%s_data17", tag));
    h_ljet_counts[0] = (TH1D*) inFile->Get (Form ("h_ljet_counts_%s_data17", tag));
    h_sljet_counts[0] = (TH1D*) inFile->Get (Form ("h_sljet_counts_%s_data17", tag));
    h_jet_pt[0] = (TH1D*) inFile->Get (Form ("h_jet_pt_%s_data17", tag));
    h2_jet_pt_cov[0] = (TH2D*) inFile->Get (Form ("h2_jet_pt_cov_%s_data17", tag));
    h2_jet_eta_phi[0] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_%s_data17", tag));

    //const double nEvts = h_evt_counts[0]->GetBinContent (1); // total number of accepted evts
    //const double nJets = h_jet_counts[0]->GetBinContent (1); // total number of accepted jets
    //const double nLJets = h_ljet_counts[0]->GetBinContent (1); // total number of accepted leading jets
    //const double nSLJets = h_sljet_counts[0]->GetBinContent (1); // total number of accepted subleading jets

    CalcUncertainties (h_jet_pt[0], h2_jet_pt_cov[0], nEvts);
  }



  //{
  //  inFile = new TFile (Form ("./rootFiles/%s/MixedHists/Nominal/data16_5TeV_hists.root", tag), "read");

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
    h_jet_pt_ratio = (TH1D*) h_jet_pt[1]->Clone ();
    h_jet_pt_ratio->Divide (h_jet_pt[0]);
  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  {
    TCanvas* c = new TCanvas ("c1", "", 800, 800);
    c->SetLogx ();
    c->SetLogy ();

    const double fuPad = 1.0;

    TH1D* h = nullptr;
    TGAE* g = nullptr;

    double ymin=5e-7;
    double ymax=8e-2;

    h = (TH1D*) h_jet_pt[0]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{jet} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fuPad);
    h->GetXaxis ()->SetLabelSize (0.028/fuPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fuPad);
    h->GetYaxis ()->SetTitle ("(1/N_{evt}) (dN_{jet} / d#it{p}_{T}^{jet})");
    h->GetYaxis ()->SetTitleSize (0.028/fuPad);
    h->GetYaxis ()->SetLabelSize (0.028/fuPad);
    h->GetYaxis ()->SetTitleOffset (2.1*fuPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = h_jet_pt[0];
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

    h = h_jet_pt[1];
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

    myText (0.58, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.040);
    myMarkerText (0.50, 0.82, myRed, kFullCircle, "2016 #it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{0-20%}", 0.8, 0.032, true);
    myMarkerText (0.50, 0.76, myBlue, kFullCircle, "2017 #it{pp}, #sqrt{s} = 5.02 TeV", 0.8, 0.032, true);

    c->SaveAs ("Plots/JetPtSpectrum.pdf");
  }



  {
    TCanvas* c = new TCanvas ("c2", "", 800, 800);
    c->SetRightMargin (0.18);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi[1];

    h2->Draw ("colz");

    c->SaveAs ("Plots/JetEtaPhiSpectrum.pdf");
  }


}


#endif
