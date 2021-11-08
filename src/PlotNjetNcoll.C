#ifndef __JetHadronCorrelatorPlotNjetNcoll_C__
#define __JetHadronCorrelatorPlotNjetNcoll_C__

#include <iostream>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TH1D.h>

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


void PlotNjetNcoll () {

  TFile* inFile = new TFile (Form ("%s/NjetNcoll.root", rootPath.Data ()), "read");

  TGAE** g_njet_ncoll = new TGAE*[2];
  TGAE** g_njet_ncoll_syst = new TGAE*[2];

  TGAE** g_njetoverncoll_ncoll = new TGAE*[2];
  TGAE** g_njetoverncoll_ncoll_syst = new TGAE*[2];

  for (short iPtJInt : {0, 1}) {
    g_njet_ncoll[iPtJInt] = (TGAE*) inFile->Get (Form ("g_njet_ncoll_%s", iPtJInt == 0 ? "30GeVJets" : "60GeVJets"));
    g_njet_ncoll_syst[iPtJInt] = (TGAE*) inFile->Get (Form ("g_njet_ncoll_syst_%s", iPtJInt == 0 ? "30GeVJets" : "60GeVJets"));
    g_njetoverncoll_ncoll[iPtJInt] = (TGAE*) inFile->Get (Form ("g_njetoverncoll_ncoll_%s", iPtJInt == 0 ? "30GeVJets" : "60GeVJets"));
    g_njetoverncoll_ncoll_syst[iPtJInt] = (TGAE*) inFile->Get (Form ("g_njetoverncoll_ncoll_syst_%s", iPtJInt == 0 ? "30GeVJets" : "60GeVJets"));
  }


  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();
  const Style_t shapes[] = {kFullCircle, kFullSquare, kOpenTriangleUp, kOpenTriangleDown};
  const Color_t colors[] = {(Color_t) TColor::GetColor (255, 102, 0), (Color_t) TColor::GetColor (0, 0, 153)};
  const Color_t systColors[] = {(Color_t) TColor::GetColor (255, 204, 0), (Color_t) TColor::GetColor (102, 102, 204)};

  {
    const char* canvasName = "c";

    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    //c->SetLogy ();

    TH1D* htemp = new TH1D ("htemp", ";#LTN_{coll}#GT;#it{R}_{#it{p}Pb}^{jet}", 1, 0, 18);
    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    xax->SetLabelFont (43);
    xax->SetLabelSize (32);

    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0;
    const double ymax = 2.6;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (0);
  
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (0, 1, 18, 1);

    tl->SetTextAlign (12);
    tl->SetTextFont (43);
    tl->SetTextSize (24);
    tl->DrawLatex (14.0, 1.10, "#LTN_{coll}#GT scaling");

    for (short iPtJInt : {1, 0}) {
      TGAE* g = g_njetoverncoll_ncoll_syst[iPtJInt];

      g->SetFillColorAlpha (systColors[iPtJInt], 1);
      g->Draw ("2");

      g = g_njetoverncoll_ncoll[iPtJInt];
      g->SetMarkerStyle (shapes[iPtJInt]);
      g->SetMarkerColor (colors[iPtJInt]);
      g->SetMarkerSize (1.6);
      g->SetLineColor (colors[iPtJInt]);
      g->SetLineWidth (3);
      g->Draw ("P");
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (12);
    tl->SetTextFont (43);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.55, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.55, 0.850, "#it{pp}, #sqrt{s} = 5.02 TeV");
    tl->DrawLatexNDC (0.55, 0.810, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV");
    //tl->DrawLatexNDC (0.55, 0.760, "#left|#eta_{CM}#right| < 2.0");
    tl->DrawLatexNDC (0.55, 0.760, "20% ZDC percentiles");
    //tl->DrawLatexNDC (0.55, 0.710, "FCal percentiles");

    mySimpleMarkerAndBoxAndLineText (0.62, 0.71, 1.4, 1001, systColors[0], 1.0, colors[0], kFullCircle, 1.6, "#it{p}_{T}^{jet} > 30 GeV", 0.028);
    mySimpleMarkerAndBoxAndLineText (0.62, 0.66, 1.4, 1001, systColors[1], 1.0, colors[1], kFullCircle, 1.6, "#it{p}_{T}^{jet} > 60 GeV", 0.028);

    c->SaveAs ("Plots/NjetOverNcoll_vs_Ncoll.pdf");
  }



  {
    const char* canvasName = "c2";

    TCanvas* c = new TCanvas (canvasName, "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    //c->SetLogy ();

    TH1D* htemp = new TH1D ("htemp", ";#LTN_{coll}#GT;#LTN_{coll}#GT #times #it{R}_{#it{p}Pb}^{jet}", 1, 0, 18);
    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    xax->SetLabelFont (43);
    xax->SetLabelSize (32);

    //yax->SetTitle ("(f_{evt} #times N_{jet} / L_{int}) / (f_{evt} #times N_{jet} / L_{int})_{#it{pp}}");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0;
    const double ymax = 19;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (0);
  
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (0, 0, 18, 18);

    tl->SetTextAlign (12);
    tl->SetTextFont (43);
    tl->SetTextSize (24);
    tl->SetTextAngle (43);
    tl->DrawLatex (14.5, 15.3, "#LTN_{coll}#GT scaling");
    tl->SetTextAngle (0);

    for (short iPtJInt : {1, 0}) {
      TGAE* g = g_njet_ncoll_syst[iPtJInt];
      g->SetFillColorAlpha (systColors[iPtJInt], 1);
      g->Draw ("2");

      g = g_njet_ncoll[iPtJInt];
      g->SetMarkerStyle (shapes[iPtJInt]);
      g->SetMarkerColor (colors[iPtJInt]);
      g->SetMarkerSize (1.6);
      g->SetLineColor (colors[iPtJInt]);
      g->SetLineWidth (3);
      g->Draw ("P");
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (12);
    tl->SetTextFont (43);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.25, 0.890, "#bf{#it{ATLAS}} Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.25, 0.850, "#it{pp}, #sqrt{s} = 5.02 TeV");
    tl->DrawLatexNDC (0.25, 0.810, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV");
    //tl->DrawLatexNDC (0.25, 0.760, "#left|#eta_{CM}#right| < 2.0");
    tl->DrawLatexNDC (0.25, 0.760, "20% ZDC percentiles");
    //tl->DrawLatexNDC (0.25, 0.710, "FCal percentiles");

    mySimpleMarkerAndBoxAndLineText (0.32, 0.71, 1.4, 1001, systColors[0], 1.0, colors[0], kFullCircle, 1.6, "#it{p}_{T}^{jet} > 30 GeV", 0.028);
    mySimpleMarkerAndBoxAndLineText (0.32, 0.66, 1.4, 1001, systColors[1], 1.0, colors[1], kFullCircle, 1.6, "#it{p}_{T}^{jet} > 60 GeV", 0.028);

    c->SaveAs ("Plots/Njet_vs_Ncoll.pdf");
  }

}


#endif
