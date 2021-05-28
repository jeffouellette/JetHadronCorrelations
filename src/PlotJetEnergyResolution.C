#ifndef __PlotJetEnergyResolution_C__
#define __PlotJetEnergyResolution_C__

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

const int numFinerEtaBins = 56;
const double* finerEtaBins = linspace (-2.8, 2.8, numFinerEtaBins);

const double etaBins[] = {0, 0.3, 0.8, 1.2, 2.1, 2.8};
const int numEtaBins = sizeof (etaBins) / sizeof (etaBins[0]) - 1;

const double enJBins[] = {10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 360, 400, 450, 500};
const int numEnJBins = sizeof (enJBins) / sizeof (enJBins[0]) - 1;

//const vector <int> systems = {0, 1};
const vector <int> systems = {0};

TH2D** h2_r2_avg_jpts = nullptr;
TH2D** h2_r2_avg_jptr = nullptr;
TH1D*** h_r2_avg_jpts = nullptr;
TH1D*** h_r2_avg_jptr = nullptr;

TH2D** h2_r4_avg_jpts = nullptr;
TH2D** h2_r4_avg_jptr = nullptr;
TH1D*** h_r4_avg_jpts = nullptr;
TH1D*** h_r4_avg_jptr = nullptr;

TH2D** h2_r2_avg_jes = nullptr;
TH2D** h2_r2_avg_jer = nullptr;
TH1D*** h_r2_avg_jes = nullptr;
TH1D*** h_r2_avg_jer = nullptr;

TH2D** h2_r4_avg_jes = nullptr;
TH2D** h2_r4_avg_jer = nullptr;
TH1D*** h_r4_avg_jes = nullptr;
TH1D*** h_r4_avg_jer = nullptr;

TH2D** h2_r2_avg_jetacorr = nullptr;
TH2D** h2_r2_avg_jetares = nullptr;
TH1D*** h_r2_avg_jetacorr = nullptr;
TH1D*** h_r2_avg_jetares = nullptr;

TH2D** h2_r4_avg_jetacorr = nullptr;
TH2D** h2_r4_avg_jetares = nullptr;
TH1D*** h_r4_avg_jetacorr = nullptr;
TH1D*** h_r4_avg_jetares = nullptr;

const Color_t binColors[10] = {kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};
const Color_t colors[11] =    {kBlack, kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};


void PlotJetEnergyResolution () {

  TFile* inFile = new TFile (Form ("%s/JetEnergyResolution/Nominal/summary.root", rootPath.Data ()), "read");

  h2_r2_avg_jpts = new TH2D*[2];
  h2_r2_avg_jptr = new TH2D*[2];
  h_r2_avg_jpts = new TH1D**[2];
  h_r2_avg_jptr = new TH1D**[2];

  h2_r4_avg_jpts = new TH2D*[2];
  h2_r4_avg_jptr = new TH2D*[2];
  h_r4_avg_jpts = new TH1D**[2];
  h_r4_avg_jptr = new TH1D**[2];

  h2_r2_avg_jes = new TH2D*[2];
  h2_r2_avg_jer = new TH2D*[2];
  h_r2_avg_jes = new TH1D**[2];
  h_r2_avg_jer = new TH1D**[2];

  h2_r4_avg_jes = new TH2D*[2];
  h2_r4_avg_jer = new TH2D*[2];
  h_r4_avg_jes = new TH1D**[2];
  h_r4_avg_jer = new TH1D**[2];

  h2_r2_avg_jetacorr = new TH2D*[2];
  h2_r2_avg_jetares = new TH2D*[2];
  h_r2_avg_jetacorr = new TH1D**[2];
  h_r2_avg_jetares = new TH1D**[2];

  h2_r4_avg_jetacorr = new TH2D*[2];
  h2_r4_avg_jetares = new TH2D*[2];
  h_r4_avg_jetacorr = new TH1D**[2];
  h_r4_avg_jetares = new TH1D**[2];

  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    h2_r2_avg_jpts[iSys] = (TH2D*) inFile->Get (Form ("h2_r2_avg_jpts_%s", sys.Data ()));
    h2_r2_avg_jptr[iSys] = (TH2D*) inFile->Get (Form ("h2_r2_avg_jptr_%s", sys.Data ()));

    h2_r4_avg_jpts[iSys] = (TH2D*) inFile->Get (Form ("h2_r4_avg_jpts_%s", sys.Data ()));
    h2_r4_avg_jptr[iSys] = (TH2D*) inFile->Get (Form ("h2_r4_avg_jptr_%s", sys.Data ()));

    h2_r2_avg_jes[iSys] = (TH2D*) inFile->Get (Form ("h2_r2_avg_jes_%s", sys.Data ()));
    h2_r2_avg_jer[iSys] = (TH2D*) inFile->Get (Form ("h2_r2_avg_jer_%s", sys.Data ()));

    h2_r4_avg_jes[iSys] = (TH2D*) inFile->Get (Form ("h2_r4_avg_jes_%s", sys.Data ()));
    h2_r4_avg_jer[iSys] = (TH2D*) inFile->Get (Form ("h2_r4_avg_jer_%s", sys.Data ()));

    h2_r2_avg_jetacorr[iSys] = (TH2D*) inFile->Get (Form ("h2_r2_avg_jetacorr_%s", sys.Data ()));
    h2_r2_avg_jetares[iSys] = (TH2D*) inFile->Get (Form ("h2_r2_avg_jetares_%s", sys.Data ()));

    h2_r4_avg_jetacorr[iSys] = (TH2D*) inFile->Get (Form ("h2_r4_avg_jetacorr_%s", sys.Data ()));
    h2_r4_avg_jetares[iSys] = (TH2D*) inFile->Get (Form ("h2_r4_avg_jetares_%s", sys.Data ()));


    h_r2_avg_jpts[iSys] = new TH1D*[numEtaBins+1];
    h_r2_avg_jptr[iSys] = new TH1D*[numEtaBins+1];

    h_r4_avg_jpts[iSys] = new TH1D*[numEtaBins+1];
    h_r4_avg_jptr[iSys] = new TH1D*[numEtaBins+1];

    h_r2_avg_jes[iSys] = new TH1D*[numEtaBins+1];
    h_r2_avg_jer[iSys] = new TH1D*[numEtaBins+1];

    h_r4_avg_jes[iSys] = new TH1D*[numEtaBins+1];
    h_r4_avg_jer[iSys] = new TH1D*[numEtaBins+1];

    h_r2_avg_jetacorr[iSys] = new TH1D*[numEtaBins+1];
    h_r2_avg_jetares[iSys] = new TH1D*[numEtaBins+1];

    h_r4_avg_jetacorr[iSys] = new TH1D*[numEtaBins+1];
    h_r4_avg_jetares[iSys] = new TH1D*[numEtaBins+1];

    for (int iEta = 0; iEta <= numEtaBins; iEta++) {

      h_r2_avg_jpts[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r2_avg_jpts_%s_iEta%i", sys.Data (), iEta));
      h_r2_avg_jptr[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r2_avg_jptr_%s_iEta%i", sys.Data (), iEta));

      h_r4_avg_jpts[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r4_avg_jpts_%s_iEta%i", sys.Data (), iEta));
      h_r4_avg_jptr[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r4_avg_jptr_%s_iEta%i", sys.Data (), iEta));

      h_r2_avg_jes[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r2_avg_jes_%s_iEta%i", sys.Data (), iEta));
      h_r2_avg_jer[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r2_avg_jer_%s_iEta%i", sys.Data (), iEta));

      h_r4_avg_jes[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r4_avg_jes_%s_iEta%i", sys.Data (), iEta));
      h_r4_avg_jer[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r4_avg_jer_%s_iEta%i", sys.Data (), iEta));

      h_r2_avg_jetacorr[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r2_avg_jetacorr_%s_iEta%i", sys.Data (), iEta));
      h_r2_avg_jetares[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r2_avg_jetares_%s_iEta%i", sys.Data (), iEta));

      h_r4_avg_jetacorr[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r4_avg_jetacorr_%s_iEta%i", sys.Data (), iEta));
      h_r4_avg_jetares[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_r4_avg_jetares_%s_iEta%i", sys.Data (), iEta));

    } // end loop over iEta
  } // end loop over iSys



  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();


  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_r2_jpts_%s", sys.Data ()), "", 800, 800);
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

      TH1D* htemp = new TH1D ("htemp", "", 1, 30, enJBins[numEnJBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{E}_{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 96;
      const double ymax = 104;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
      tl->DrawLatex (200, yoff, "200");
      tl->DrawLatex (300, yoff, "300");
      tl->DrawLatex (500, yoff, "500");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (30, 100, enJBins[numEnJBins], 100);
    }
     

    for (int iEta = 0; iEta < numEtaBins; iEta++) {
      TGAE* g = make_graph (h_r2_avg_jpts[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_r2_avg_jpts[iSys][numEtaBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, "Anti-k_{T} HI Jets, R=0.2");
    }
    else if (iSys == 1) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.845, "Pythia8 #it{pp} + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < numEtaBins; iEta++)
      myMarkerTextNoLine (0.3, 0.39-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.3, 0.44, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JPTS_R2_%s.pdf", workPath.Data (), sys.Data ()));
  } // end loop over iSys



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_r4_jpts_%s", sys.Data ()), "", 800, 800);
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

      TH1D* htemp = new TH1D ("htemp", "", 1, 30, enJBins[numEnJBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{E}_{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 96;
      const double ymax = 104;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
      //tl->DrawLatex (10, yoff, "10");
      //tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
      tl->DrawLatex (200, yoff, "200");
      tl->DrawLatex (300, yoff, "300");
      tl->DrawLatex (500, yoff, "500");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (30, 100, enJBins[numEnJBins], 100);
    }
     

    for (int iEta = 0; iEta < numEtaBins; iEta++) {
      TGAE* g = make_graph (h_r4_avg_jpts[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_r4_avg_jpts[iSys][numEtaBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, "Anti-k_{T} HI Jets, R=0.4");
    }
    else if (iSys == 1) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.845, "Pythia8 #it{pp} + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < numEtaBins; iEta++)
      myMarkerTextNoLine (0.3, 0.39-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.3, 0.44, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JPTS_R4_%s.pdf", workPath.Data (), sys.Data ()));
  } // end loop over iSys



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_r2_jes_%s", sys.Data ()), "", 800, 800);
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

      TH1D* htemp = new TH1D ("htemp", "", 1, 30, enJBins[numEnJBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{E}_{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 96;
      const double ymax = 104;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
      //tl->DrawLatex (10, yoff, "10");
      //tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
      tl->DrawLatex (200, yoff, "200");
      tl->DrawLatex (300, yoff, "300");
      tl->DrawLatex (500, yoff, "500");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (30, 100, enJBins[numEnJBins], 100);
    }
     

    for (int iEta = 0; iEta < numEtaBins; iEta++) {
      TGAE* g = make_graph (h_r2_avg_jes[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_r2_avg_jes[iSys][numEtaBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, "Anti-k_{T} HI Jets, R=0.2");
    }
    else if (iSys == 1) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.845, "Pythia8 #it{pp} + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < numEtaBins; iEta++)
      myMarkerTextNoLine (0.3, 0.39-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.3, 0.44, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JES_R2_%s.pdf", workPath.Data (), sys.Data ()));
  } // end loop over iSys



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_r4_jes_%s", sys.Data ()), "", 800, 800);
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

      TH1D* htemp = new TH1D ("htemp", "", 1, 30, enJBins[numEnJBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{E}_{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 96;
      const double ymax = 104;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
      //tl->DrawLatex (10, yoff, "10");
      //tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
      tl->DrawLatex (200, yoff, "200");
      tl->DrawLatex (300, yoff, "300");
      tl->DrawLatex (500, yoff, "500");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (30, 100, enJBins[numEnJBins], 100);
    }
     

    for (int iEta = 0; iEta < numEtaBins; iEta++) {
      TGAE* g = make_graph (h_r4_avg_jes[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_r4_avg_jes[iSys][numEtaBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, "Anti-k_{T} HI Jets, R=0.4");
    }
    else if (iSys == 1) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.845, "Pythia8 #it{pp} + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < numEtaBins; iEta++)
      myMarkerTextNoLine (0.3, 0.39-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.3, 0.44, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JES_R4_%s.pdf", workPath.Data (), sys.Data ()));
  } // end loop over iSys



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_r2_jptr_%s", sys.Data ()), "", 800, 800);
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

      TH1D* htemp = new TH1D ("htemp", "", 1, 30, enJBins[numEnJBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{E}_{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0;
      const double ymax = 40;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
      //tl->DrawLatex (10, yoff, "10");
      //tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
      tl->DrawLatex (200, yoff, "200");
      tl->DrawLatex (300, yoff, "300");
      tl->DrawLatex (500, yoff, "500");
    }
     

    for (int iEta = 0; iEta < numEtaBins; iEta++) {
      TGAE* g = make_graph (h_r2_avg_jptr[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_r2_avg_jptr[iSys][numEtaBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, "Anti-k_{T} HI Jets, R=0.2");
    }
    else if (iSys == 1) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.845, "Pythia8 #it{pp} + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < numEtaBins; iEta++)
      myMarkerTextNoLine (0.5, 0.69-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.5, 0.74, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JPTR_R2_%s.pdf", workPath.Data (), sys.Data ()));
  }



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_r4_jptr_%s", sys.Data ()), "", 800, 800);
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

      TH1D* htemp = new TH1D ("htemp", "", 1, 30, enJBins[numEnJBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{E}_{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0;
      const double ymax = 40;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
      //tl->DrawLatex (10, yoff, "10");
      //tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
      tl->DrawLatex (200, yoff, "200");
      tl->DrawLatex (300, yoff, "300");
      tl->DrawLatex (500, yoff, "500");
    }
     

    for (int iEta = 0; iEta < numEtaBins; iEta++) {
      TGAE* g = make_graph (h_r4_avg_jptr[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_r4_avg_jptr[iSys][numEtaBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, "Anti-k_{T} HI Jets, R=0.4");
    }
    else if (iSys == 1) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.845, "Pythia8 #it{pp} + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < numEtaBins; iEta++)
      myMarkerTextNoLine (0.5, 0.69-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.5, 0.74, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JPTR_R4_%s.pdf", workPath.Data (), sys.Data ()));
  }



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_r2_jer_%s", sys.Data ()), "", 800, 800);
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

      TH1D* htemp = new TH1D ("htemp", "", 1, 30, enJBins[numEnJBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{E}_{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0;
      const double ymax = 40;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
      //tl->DrawLatex (10, yoff, "10");
      //tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
      tl->DrawLatex (200, yoff, "200");
      tl->DrawLatex (300, yoff, "300");
      tl->DrawLatex (500, yoff, "500");
    }
     

    for (int iEta = 0; iEta < numEtaBins; iEta++) {
      TGAE* g = make_graph (h_r2_avg_jer[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_r2_avg_jer[iSys][numEtaBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, "Anti-k_{T} HI Jets, R=0.2");
    }
    else if (iSys == 1) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.845, "Pythia8 #it{pp} + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < numEtaBins; iEta++)
      myMarkerTextNoLine (0.5, 0.69-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.5, 0.74, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JER_R2_%s.pdf", workPath.Data (), sys.Data ()));
  }



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_r4_jer_%s", sys.Data ()), "", 800, 800);
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

      TH1D* htemp = new TH1D ("htemp", "", 1, 30, enJBins[numEnJBins]);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{E}_{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0;
      const double ymax = 40;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.8,  yoff, "0.8");
      //tl->DrawLatex (10, yoff, "10");
      //tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
      //tl->DrawLatex (80, yoff, "80");
      tl->DrawLatex (100, yoff, "100");
      tl->DrawLatex (200, yoff, "200");
      tl->DrawLatex (300, yoff, "300");
      tl->DrawLatex (500, yoff, "500");
    }
     

    for (int iEta = 0; iEta < numEtaBins; iEta++) {
      TGAE* g = make_graph (h_r4_avg_jer[iSys][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_r4_avg_jer[iSys][numEtaBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (26);
    tl->DrawLatexNDC (0.26, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    if (iSys == 0) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
      tl->DrawLatexNDC (0.26, 0.810, "Anti-k_{T} HI Jets, R=0.4");
    }
    else if (iSys == 1) {
      tl->SetTextSize (26);
      tl->DrawLatexNDC (0.26, 0.845, "Pythia8 #it{pp} + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
    }
    for (int iEta = 0; iEta < numEtaBins; iEta++)
      myMarkerTextNoLine (0.5, 0.69-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.036);
    myMarkerTextNoLine (0.5, 0.74, colors[0], kFullCircle, "All #it{#eta}", 1.2, 0.036);

    c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JER_R4_%s.pdf", workPath.Data (), sys.Data ()));
  }


}

#endif
