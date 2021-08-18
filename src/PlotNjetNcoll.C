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

  TFile* f1 = new TFile ("rootFiles/Results/ProcessJets_60GeVJets.root", "read");
  TFile* f2 = new TFile ("rootFiles/Results/ProcessJets_30GeVJets.root", "read");
  TFile** files = new TFile*[2];
  files[0] = f1;
  files[1] = f2;
  const int nFiles = 2;

  TH1D**  h_jet_counts_ref  = Get1DArray <TH1D*> (nFiles);
  TH1D*** h_jet_counts      = Get2DArray <TH1D*> (nFiles, nZdcCentBins);

  TGAE** g_njet_ncoll = new TGAE*[nFiles];
  TGAE** g_njet_ncoll_syst = new TGAE*[nFiles];

  TGAE** g_njetoverncoll_ncoll = new TGAE*[nFiles];
  TGAE** g_njetoverncoll_ncoll_syst = new TGAE*[nFiles];

  for (int iFile = 0; iFile < nFiles; iFile++) {

    TFile* inFile = files[iFile];

    h_jet_counts_ref[iFile] = (TH1D*) inFile->Get (Form ("h_jet_counts_ref_data_Nominal"));
  
    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
  
      h_jet_counts[iFile][iCent] = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_iCent%i_data_Nominal", iCent));
  
    } // end loop over iCent

    g_njet_ncoll[iFile] = new TGAE ();
    g_njet_ncoll_syst[iFile] = new TGAE ();

    g_njetoverncoll_ncoll[iFile] = new TGAE ();
    g_njetoverncoll_ncoll_syst[iFile] = new TGAE ();

    double njetlumipp = 1;

    {
      const double njet = h_jet_counts_ref[iFile]->GetBinContent (1);
      const double ncoll = 1;
      const double sigma_ncoll = 0;
      const double lumi = (iFile == 0 ? 3.57487e6 : 2.65513e3) * 0.244; // multiply by Poisson probability of 1 with lambda=2.2 (average pileup rate for run 340718).

      njetlumipp = njet / lumi;

      g_njet_ncoll[iFile]->SetPoint (0, ncoll, 1);
      g_njet_ncoll[iFile]->SetPointEXhigh (0, 0);
      g_njet_ncoll[iFile]->SetPointEXlow (0, 0);
      g_njet_ncoll[iFile]->SetPointEYhigh (0, 0);
      g_njet_ncoll[iFile]->SetPointEYlow (0, 0);

      g_njet_ncoll_syst[iFile]->SetPoint (0, ncoll, 1);
      g_njet_ncoll_syst[iFile]->SetPointEXhigh (0, 0);
      g_njet_ncoll_syst[iFile]->SetPointEXlow (0, 0);
      g_njet_ncoll_syst[iFile]->SetPointEYhigh (0, 0);
      g_njet_ncoll_syst[iFile]->SetPointEYlow (0, 0);

      g_njetoverncoll_ncoll[iFile]->SetPoint (0, ncoll, 1);//njet / (lumi * ncoll));
      g_njetoverncoll_ncoll[iFile]->SetPointEXhigh (0, 0);
      g_njetoverncoll_ncoll[iFile]->SetPointEXlow (0, 0);
      g_njetoverncoll_ncoll[iFile]->SetPointEYhigh (0, 0);
      g_njetoverncoll_ncoll[iFile]->SetPointEYlow (0, 0);

      g_njetoverncoll_ncoll_syst[iFile]->SetPoint (0, ncoll, 1);
      g_njetoverncoll_ncoll_syst[iFile]->SetPointEXhigh (0, 0);
      g_njetoverncoll_ncoll_syst[iFile]->SetPointEXlow (0, 0);
      g_njetoverncoll_ncoll_syst[iFile]->SetPointEYhigh (0, 0);
      g_njetoverncoll_ncoll_syst[iFile]->SetPointEYlow (0, 0);

    }

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {

      const double njet = h_jet_counts[iFile][iCent]->GetBinContent (1);
      const double ncoll = zdcNcollValues[iCent];
      const double sigma_ncoll = zdcNcollErrors[iCent];
      const double lumi = (iFile == 0 ? 3.55545e2 : 25.1334) / 0.2;

      g_njet_ncoll[iFile]->SetPoint (iCent+1, ncoll, njet / (njetlumipp * lumi));
      g_njet_ncoll[iFile]->SetPointEXhigh (iCent+1, 0);
      g_njet_ncoll[iFile]->SetPointEXlow (iCent+1, 0);
      g_njet_ncoll[iFile]->SetPointEYhigh (iCent+1, std::sqrt (njet) / (njetlumipp * lumi));
      g_njet_ncoll[iFile]->SetPointEYlow (iCent+1, std::sqrt (njet) / (njetlumipp * lumi));

      g_njet_ncoll_syst[iFile]->SetPoint (iCent+1, ncoll, njet / (njetlumipp * lumi));
      g_njet_ncoll_syst[iFile]->SetPointEXhigh (iCent+1, sigma_ncoll);
      g_njet_ncoll_syst[iFile]->SetPointEXlow (iCent+1, sigma_ncoll);
      g_njet_ncoll_syst[iFile]->SetPointEYhigh (iCent+1, 0.4);
      g_njet_ncoll_syst[iFile]->SetPointEYlow (iCent+1, 0.4);

      g_njetoverncoll_ncoll[iFile]->SetPoint (iCent+1, ncoll, njet / (njetlumipp * lumi * ncoll));
      g_njetoverncoll_ncoll[iFile]->SetPointEXhigh (iCent+1, 0);
      g_njetoverncoll_ncoll[iFile]->SetPointEXlow (iCent+1, 0);
      g_njetoverncoll_ncoll[iFile]->SetPointEYhigh (iCent+1, std::sqrt (njet / std::pow (njetlumipp * lumi * ncoll, 2)));
      g_njetoverncoll_ncoll[iFile]->SetPointEYlow (iCent+1, std::sqrt (njet / std::pow (njetlumipp * lumi * ncoll, 2)));

      g_njetoverncoll_ncoll_syst[iFile]->SetPoint (iCent, ncoll, njet / (njetlumipp * lumi * ncoll));
      g_njetoverncoll_ncoll_syst[iFile]->SetPointEXhigh (iCent, sigma_ncoll);
      g_njetoverncoll_ncoll_syst[iFile]->SetPointEXlow (iCent, sigma_ncoll);
      g_njetoverncoll_ncoll_syst[iFile]->SetPointEYhigh (iCent, std::fabs (njet * sigma_ncoll / (njetlumipp * lumi * ncoll * ncoll)));
      g_njetoverncoll_ncoll_syst[iFile]->SetPointEYlow (iCent, std::fabs (njet * sigma_ncoll / (njetlumipp * lumi * ncoll * ncoll)));

    } // end loop over iCent

  }


  TLatex* tl = new TLatex ();

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

    TH1D* htemp = new TH1D ("htemp", ";#LTN_{coll}#GT;(1 / #LTN_{coll}#GT) #times (#sigma_{jet} / #sigma_{jet}^{#it{pp}})", 1, 0, 16);
    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);

    //yax->SetTitle ("f_{evt} #times N_{jet} / (#LTN_{coll}#GT #times L_{int}) [#mub]");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0;
    const double ymax = 3;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (0);
  
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    const Style_t shapes[] = {kFullCircle, kFullSquare, kOpenTriangleUp, kOpenTriangleDown};
    const Color_t colors[] = {(Color_t) TColor::GetColor (0, 0, 153), (Color_t) TColor::GetColor (255, 102, 0)};
    const Color_t systColors[] = {(Color_t) TColor::GetColor (102, 102, 204), (Color_t) TColor::GetColor (255, 204, 0)};
    for (int iFile = 0; iFile < nFiles; iFile++) {
      TGAE* g = g_njetoverncoll_ncoll_syst[iFile];

      g->SetFillColorAlpha (systColors[iFile], 1);
      g->Draw ("2");

      g = g_njetoverncoll_ncoll[iFile];
      g->SetMarkerStyle (shapes[iFile]);
      g->SetMarkerColor (colors[iFile]);
      g->SetMarkerSize (1.6);
      g->SetLineColor (colors[iFile]);
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

    mySimpleMarkerAndBoxAndLineText (0.62, 0.76, 1.4, 1001, systColors[0], 1.0, colors[0], kFullCircle, 1.6, "#it{p}_{T}^{jet} > 60 GeV", 0.028);
    mySimpleMarkerAndBoxAndLineText (0.62, 0.71, 1.4, 1001, systColors[1], 1.0, colors[1], kFullCircle, 1.6, "#it{p}_{T}^{jet} > 30 GeV", 0.028);

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

    TH1D* htemp = new TH1D ("htemp", ";#LTN_{coll}#GT;#sigma_{jet} / #sigma_{jet}^{#it{pp}}", 1, 0, 16);
    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);

    //yax->SetTitle ("(f_{evt} #times N_{jet} / L_{int}) / (f_{evt} #times N_{jet} / L_{int})_{#it{pp}}");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0;
    const double ymax = 24;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (0);
  
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    const Style_t shapes[] = {kFullCircle, kFullSquare, kOpenTriangleUp, kOpenTriangleDown};
    const Color_t colors[] = {(Color_t) TColor::GetColor (0, 0, 153), (Color_t) TColor::GetColor (255, 102, 0)};
    const Color_t systColors[] = {(Color_t) TColor::GetColor (102, 102, 204), (Color_t) TColor::GetColor (255, 204, 0)};
    for (int iFile = 0; iFile < nFiles; iFile++) {
      TGAE* g = g_njet_ncoll_syst[iFile];
      g->SetFillColorAlpha (systColors[iFile], 1);
      g->Draw ("2");

      g = g_njet_ncoll[iFile];
      g->SetMarkerStyle (shapes[iFile]);
      g->SetMarkerColor (colors[iFile]);
      g->SetMarkerSize (1.6);
      g->SetLineColor (colors[iFile]);
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

    mySimpleMarkerAndBoxAndLineText (0.62, 0.76, 1.4, 1001, systColors[0], 1.0, colors[0], kFullCircle, 1.6, "#it{p}_{T}^{jet} > 60 GeV", 0.028);
    mySimpleMarkerAndBoxAndLineText (0.62, 0.71, 1.4, 1001, systColors[1], 1.0, colors[1], kFullCircle, 1.6, "#it{p}_{T}^{jet} > 30 GeV", 0.028);

    c->SaveAs ("Plots/Njet_vs_Ncoll.pdf");
  }

}


#endif
