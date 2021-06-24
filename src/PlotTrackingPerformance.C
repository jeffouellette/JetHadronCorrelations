#ifndef __PlotTrackingPerformance_C__
#define __PlotTrackingPerformance_C__

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>

#include <Utilities.h>

#include "CentralityDefs.h"
#include "Params.h"
#include "LocalUtilities.h"

using namespace JetHadronCorrelations;

const Color_t binColors[10] = {kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};


void PlotTrackingPerformance () {

  TFile* inFile = new TFile (Form ("%s/TrackingPerformance/Nominal/outFile.root", rootPath.Data ()), "read");

//  TH1D* h_truth_matching_prob[2] = {};
// = new TH1D (Form ("h_truth_matching_prob_%s", sys.Data ()), ";Truth matching prob.;N_{ch}^{rec}", 200, 0, 1);

  const int nFinerEtaTrkBins = 40;
  const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

  const double etaTrkBins[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
  const int nEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

//const double pTchBins[32] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 60, 80, 100};
//const int nPtChBins = 106;
  const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100};
  const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

  TH1D** h_truth_matching_prob = Get1DArray <TH1D*> (3);

  TH2D*** h2_truth_matched_reco_tracks = Get2DArray <TH2D*> (2, 3);
  TH2D*** h2_truth_tracks = Get2DArray <TH2D*> (2, 3);
  TH2D*** h2_efficiency = Get2DArray <TH2D*> (2, 3);

  TH2D*** h2_primary_reco_tracks = Get2DArray <TH2D*> (2, 3);
  TH2D*** h2_reco_tracks = Get2DArray <TH2D*> (2, 3);
  TH2D*** h2_purity = Get2DArray <TH2D*> (2, 3);


  TH1D**** h_truth_matched_reco_tracks = Get3DArray <TH1D*> (2, 3, nEtaTrkBins);
  TH1D**** h_truth_tracks = Get3DArray <TH1D*> (2, 3, nEtaTrkBins);
  TH1D**** h_efficiency = Get3DArray <TH1D*> (2, 3, nEtaTrkBins);


  TH1D**** h_primary_reco_tracks = Get3DArray <TH1D*> (2, 3, nEtaTrkBins);
  TH1D**** h_reco_tracks = Get3DArray <TH1D*> (2, 3, nEtaTrkBins);
  TH1D**** h_purity = Get3DArray <TH1D*> (2, 3, nEtaTrkBins);

  for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

    h_truth_matching_prob[iWP] = new TH1D (Form ("h_truth_matching_prob_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()), ";Truth matching prob.;N_{ch}^{rec}", 200, 0, 1);
    h_truth_matching_prob[iWP]->Sumw2 ();

    h2_truth_matched_reco_tracks[iWP] = new TH2D (Form ("h2_truth_matched_reco_tracks_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
    h2_truth_matched_reco_tracks[iWP]->Sumw2 ();

    h2_truth_tracks[iWP] = new TH2D (Form ("h2_truth_tracks_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
    h2_truth_tracks[iWP]->Sumw2 ();


    h2_primary_reco_tracks[iWP] = new TH2D (Form ("h2_primary_reco_tracks_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
    h2_primary_reco_tracks[iWP]->Sumw2 ();

    h2_reco_tracks[iWP] = new TH2D (Form ("h2_reco_tracks_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()), ";#eta;#it{p}_{T} [GeV]", nFinerEtaTrkBins, finerEtaTrkBins, nPtchBins, pTchBins);
    h2_reco_tracks[iWP]->Sumw2 ();

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      h_truth_matched_reco_tracks[iWP][iEta] = new TH1D (Form ("h_truth_matched_reco_tracks_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
      h_truth_matched_reco_tracks[iWP][iEta]->Sumw2 ();

      h_truth_tracks[iWP][iEta] = new TH1D (Form ("h_truth_tracks_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
      h_truth_tracks[iWP][iEta]->Sumw2 ();


      h_primary_reco_tracks[iWP][iEta] = new TH1D (Form ("h_primary_reco_tracks_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
      h_primary_reco_tracks[iWP][iEta]->Sumw2 ();

      h_reco_tracks[iWP][iEta] = new TH1D (Form ("h_reco_tracks_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta), ";#it{p}_{T} [GeV]", nPtchBins, pTchBins);
      h_reco_tracks[iWP][iEta]->Sumw2 ();
    }

  }


  //for (int iSys : {0, 1}) {
  for (int iSys : {0}) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iWP = 0; iWP < trackWPs.size (); iWP++) {

      h2_truth_matched_reco_tracks[iSys][iWP] = (TH2D*) inFile->Get (Form ("h2_truth_matched_reco_tracks_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()));
      h2_truth_tracks[iSys][iWP] = (TH2D*) inFile->Get (Form ("h2_truth_tracks_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()));

      h2_efficiency[iSys][iWP] = (TH2D*) h2_truth_matched_reco_tracks[iSys][iWP]->Clone (Form ("h2_efficiency_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()));
      BinomialDivide (h2_efficiency[iSys][iWP], h2_truth_matched_reco_tracks[iSys][iWP], h2_truth_tracks[iSys][iWP]);

      h2_primary_reco_tracks[iSys][iWP] = (TH2D*) inFile->Get (Form ("h2_primary_reco_tracks_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()));
      h2_reco_tracks[iSys][iWP] = (TH2D*) inFile->Get (Form ("h2_reco_tracks_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()));

      h2_purity[iSys][iWP] = (TH2D*) h2_primary_reco_tracks[iSys][iWP]->Clone (Form ("h2_purity_%s_%s", sys.Data (), trackWPNames[iWP].c_str ()));
      BinomialDivide (h2_purity[iSys][iWP], h2_primary_reco_tracks[iSys][iWP], h2_reco_tracks[iSys][iWP]);

      for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
        h_truth_matched_reco_tracks[iSys][iWP][iEta] = (TH1D*) inFile->Get (Form ("h_truth_matched_reco_tracks_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta));
        h_truth_tracks[iSys][iWP][iEta] = (TH1D*) inFile->Get (Form ("h_truth_tracks_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta));

        h_efficiency[iSys][iWP][iEta] = (TH1D*) h_truth_matched_reco_tracks[iSys][iWP][iEta]->Clone (Form ("h_efficiency_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta));
        BinomialDivide (h_efficiency[iSys][iWP][iEta], h_truth_matched_reco_tracks[iSys][iWP][iEta], h_truth_tracks[iSys][iWP][iEta]);

        h_primary_reco_tracks[iSys][iWP][iEta] = (TH1D*) inFile->Get (Form ("h_primary_reco_tracks_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta));
        h_reco_tracks[iSys][iWP][iEta] = (TH1D*) inFile->Get (Form ("h_reco_tracks_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta));

        h_purity[iSys][iWP][iEta] = (TH1D*) h_primary_reco_tracks[iSys][iWP][iEta]->Clone (Form ("h_purity_%s_%s_iEta%i", sys.Data (), trackWPNames[iWP].c_str (), iEta));
        BinomialDivide (h_purity[iSys][iWP][iEta], h_primary_reco_tracks[iSys][iWP][iEta], h_reco_tracks[iSys][iWP][iEta]);
      }
    }

  }



  const short iWP = 0;

  {
    TCanvas* c = new TCanvas ("c1", "", 800, 800);

    gPad->SetLogx ();

    gPad->SetBottomMargin (0.11);
    gPad->SetLeftMargin (0.11);
    gPad->SetRightMargin (0.04);
    gPad->SetTopMargin (0.04);

    TH1D* h = (TH1D*) h_efficiency[iWP][0][0]->Clone ("temp");
    h->Reset ();
    for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 1);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetYaxis ()->SetTitle ("Track Reco. Efficiency");

    double ymin = 0.5;
    double ymax = 1.06;

    h->SetLineWidth (0);

    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (24);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (24);

    h->DrawCopy ("hist");

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      h = h_efficiency[iWP][0][iEta];

      h->SetMarkerSize (0.8);
      h->SetMarkerColor (binColors[iEta%10]);
      h->SetLineColor (binColors[iEta%10]);

      h->DrawCopy ("e1 same");
    }

    myText (0.62, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.62, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.036);
    myText (0.62, 0.80, kBlack, "HITight tracks", 0.036);
    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      myLineText2 (0.68, 0.76-iEta*0.04, binColors[iEta%10], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.036);
    }

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencySummary.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c2", "", 800, 800);

    gPad->SetLogx ();

    gPad->SetBottomMargin (0.11);
    gPad->SetLeftMargin (0.11);
    gPad->SetRightMargin (0.04);
    gPad->SetTopMargin (0.04);

    TH1D* h = (TH1D*) h_purity[iWP][0][0]->Clone ("temp");
    h->Reset ();
    for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 1);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h->GetYaxis ()->SetTitle ("Track Purity");

    double ymin = 0.92;
    double ymax = 1.01;

    h->SetLineWidth (0);

    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (24);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (24);

    h->DrawCopy ("hist");

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      h = h_purity[iWP][0][iEta];

      h->SetMarkerSize (0.8);
      h->SetMarkerColor (binColors[iEta%10]);
      h->SetLineColor (binColors[iEta%10]);

      h->DrawCopy ("e1 same");
    }

    myText (0.62, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.62, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.036);
    myText (0.62, 0.80, kBlack, "HITight tracks", 0.036);
    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      myLineText2 (0.68, 0.76-iEta*0.04, binColors[iEta%10], kFullCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.036);
    }


    c->SaveAs (Form ("%s/Plots/TrackingPerformance/PuritySummary.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c3", "", 800, 800);

    gPad->SetLogy ();

    gPad->SetBottomMargin (0.11);
    gPad->SetLeftMargin (0.11);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.04);

    TH2D* h2 = (TH2D*) h2_efficiency[iWP][0]->Clone ("temp");
    h2->GetXaxis ()->SetTitle ("#eta");
    h2->GetYaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h2->GetZaxis ()->SetTitle ("Track Reco. Efficiency");

    double zmin = 0.5;
    double zmax = 1;

    h2->SetLineWidth (0);

    h2->GetZaxis ()->SetRangeUser (zmin, zmax);
    h2->GetYaxis ()->SetMoreLogLabels ();

    h2->GetXaxis ()->SetTitleFont (43);
    h2->GetXaxis ()->SetTitleSize (26);
    h2->GetYaxis ()->SetTitleFont (43);
    h2->GetYaxis ()->SetTitleSize (26);
    h2->GetZaxis ()->SetTitleFont (43);
    h2->GetZaxis ()->SetTitleSize (26);
    h2->GetZaxis ()->SetTitleOffset (1.4*h2->GetZaxis ()->GetTitleOffset ());
    h2->GetXaxis ()->SetLabelFont (43);
    h2->GetXaxis ()->SetLabelSize (24);
    h2->GetYaxis ()->SetLabelFont (43);
    h2->GetYaxis ()->SetLabelSize (24);
    h2->GetZaxis ()->SetLabelFont (43);
    h2->GetZaxis ()->SetLabelSize (24);

    h2->DrawCopy ("colz");
    SaferDelete (&h2);

    myText (0.62, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.62, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.036);
    myText (0.62, 0.80, kBlack, "HITight tracks", 0.036);

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencyMap_pp.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c4", "", 800, 800);

    gPad->SetLogy ();

    gPad->SetBottomMargin (0.11);
    gPad->SetLeftMargin (0.11);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.04);

    TH2D* h2 = (TH2D*) h2_purity[iWP][0]->Clone ("temp");
    h2->GetXaxis ()->SetTitle ("#eta");
    h2->GetYaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
    h2->GetZaxis ()->SetTitle ("Primary track fraction");

    double zmin = 0.8;
    double zmax = 1;

    h2->SetLineWidth (0);

    h2->GetZaxis ()->SetRangeUser (zmin, zmax);
    h2->GetYaxis ()->SetMoreLogLabels ();

    h2->GetXaxis ()->SetTitleFont (43);
    h2->GetXaxis ()->SetTitleSize (26);
    h2->GetYaxis ()->SetTitleFont (43);
    h2->GetYaxis ()->SetTitleSize (26);
    h2->GetZaxis ()->SetTitleFont (43);
    h2->GetZaxis ()->SetTitleSize (26);
    h2->GetZaxis ()->SetTitleOffset (1.4*h2->GetZaxis ()->GetTitleOffset ());
    h2->GetXaxis ()->SetLabelFont (43);
    h2->GetXaxis ()->SetLabelSize (24);
    h2->GetYaxis ()->SetLabelFont (43);
    h2->GetYaxis ()->SetLabelSize (24);
    h2->GetZaxis ()->SetLabelFont (43);
    h2->GetZaxis ()->SetLabelSize (24);

    h2->DrawCopy ("colz");
    SaferDelete (&h2);

    myText (0.62, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.62, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.036);
    myText (0.62, 0.80, kBlack, "HITight tracks", 0.036);

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/PurityMap_pp.pdf", workPath.Data ()));
  }

}

#endif
