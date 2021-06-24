#ifndef __PlotTrackMomentumResolution_C__
#define __PlotTrackMomentumResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

#include <Utilities.h>

#include "Params.h"
#include "CentralityDefs.h"
#include "LocalUtilities.h"

using namespace JetHadronCorrelations;

typedef TGraphAsymmErrors TGAE;

const int nFinerEtaTrkBins = 40;
const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

const double etaTrkBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const int nEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100};
const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

//const vector <int> systems = {0, 1};
const vector <int> systems = {0};

TH1D**** h_tms = nullptr;
TH2D** h2_avg_tms = nullptr;
TH2D** h2_avg_tmr = nullptr;
TH1D*** h_avg_tms = nullptr;
TH1D*** h_avg_tmr = nullptr;

const Color_t binColors[10] = {kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};
const Color_t colors[11] =    {kBlack, kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};


void PlotTrackMomentumResolution () {

  TFile* inFile = new TFile (Form ("%s/TrackMomentumResolution/Nominal/summary.root", rootPath.Data ()), "read");

  h2_avg_tms = new TH2D*[2];
  h2_avg_tmr = new TH2D*[2];
  h_avg_tms = new TH1D**[2];
  h_avg_tmr = new TH1D**[2];

  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    h2_avg_tms[iSys] = (TH2D*) inFile->Get (Form ("h2_avg_tms_%s", sys.Data ()));
    h2_avg_tmr[iSys] = (TH2D*) inFile->Get (Form ("h2_avg_tmr_%s", sys.Data ()));

    h_avg_tms[iSys] = new TH1D*[nEtaTrkBins+1];
    h_avg_tmr[iSys] = new TH1D*[nEtaTrkBins+1];
    for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {
      h_avg_tms[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_avg_tms_%s_iEta%i", sys.Data (), iEta));
      h_avg_tmr[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_avg_tmr_%s_iEta%i", sys.Data (), iEta));
    } // end loop over iEta
  } // end loop over iSys



  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();


  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_tms_%s", sys.Data ()), "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    {
      gPad->SetLogx ();

      TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("TMS [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 99.5;
      const double ymax = 100.5;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      tl->DrawLatex (0.5,  yoff, "0.5");
      tl->DrawLatex (0.7,  yoff, "0.7");
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (pTchBins[0], 100, pTchBins[nPtchBins], 100);
    }
     

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      TGAE* g = make_graph (h_avg_tms[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_avg_tms[iSys][nEtaTrkBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (30);
      tl->DrawLatexNDC (0.26, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.800, "Pythia8 dijets");
    }
    else if (iSys == 1) {
      tl->SetTextSize (30);
      tl->DrawLatexNDC (0.26, 0.845, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myMarkerTextNoLine (0.3, 0.39-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.3, 0.44, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/TrackMomentumResolution/TMS_%s.pdf", workPath.Data (), sys.Data ()));
  } // end loop over iSys



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_tmr_%s", sys.Data ()), "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    {
      gPad->SetLogx ();

      TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("TMR [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0;
      const double ymax = 20;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      tl->DrawLatex (0.5,  yoff, "0.5");
      tl->DrawLatex (0.7,  yoff, "0.7");
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
    }
     

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++) {
      TGAE* g = make_graph (h_avg_tmr[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_avg_tmr[iSys][nEtaTrkBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (30);
      tl->DrawLatexNDC (0.26, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.800, "Pythia8 dijets");
    }
    else if (iSys == 1) {
      tl->SetTextSize (30);
      tl->DrawLatexNDC (0.26, 0.845, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV");
    }

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myMarkerTextNoLine (0.5, 0.69-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.5, 0.74, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/TrackMomentumResolution/TMR_%s.pdf", workPath.Data (), sys.Data ()));
  }


}

#endif
