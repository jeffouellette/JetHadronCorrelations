#ifndef __PlotTrackingPerformanceWgtsComp_C__
#define __PlotTrackingPerformanceWgtsComp_C__

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

#include <ArrayTemplates.h>
#include <Utilities.h>
#include <MyStyle.h>
#include <MyColors.h>

#include "CentralityDefs.h"
#include "TreeVariables.h"
#include "Params.h"
#include "LocalUtilities.h"
#include "PiecewisePolynomialConstantFunc.h"

#include "AtlasStyle.h"

using namespace JetHadronCorrelations;

// Bin edge definitions: eta, pTch, pp/pPb.
const int nFinerEtaTrkBins = 40;
const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

//const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 70, 80, 100}; // old binning
const double pTchBins[] = {0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130}; // new binning, 113 elements or 112 bins
const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

                  //pi+, k+,  p+,   e-, mu-, sigma+, sigma-, xi-,  omega-, everyone
const short PIDs[] = {211, 321, 2212, 11, 13,  3222,   3112,   3312, 3334,   0};
const std::string partNames[] = {"#pi^{+}", "K^{+}", "p^{+}", "e^{-}", "#mu^{-}", "#Sigma^{+}", "#Sigma^{-}", "#Xi^{-}", "#Omega^{-}", "All h^{#pm}"};
const short nPIDs = sizeof (PIDs) / sizeof (PIDs[0]);

const vector <short> systems = {0, 1}; // 0 = pp, 1 = pPb
const short nSystems = systems.size ();


void PlotTrackingPerformanceWgtsComp () {

  TH2D***** h2_truth_matched_primary_tracks = Get4DArray <TH2D*> (2, nSystems, nPIDs, nMultBins);
  TH2D***** h2_truth_tracks                 = Get4DArray <TH2D*> (2, nSystems, nPIDs, nMultBins);
  TH2D***** h2_truth_tracks_wgt2            = Get4DArray <TH2D*> (2, nSystems, nPIDs, nMultBins);

  TH2D***** h2_efficiency                   = Get4DArray <TH2D*> (2, nSystems, nPIDs, nMultBins);

  TH1D****** h_truth_matched_primary_tracks  = Get5DArray <TH1D*> (2, nSystems, nPIDs, nMultBins, nEtaTrkBins);
  TH1D****** h_truth_tracks                  = Get5DArray <TH1D*> (2, nSystems, nPIDs, nMultBins, nEtaTrkBins);
  TH1D****** h_truth_tracks_wgt2             = Get5DArray <TH1D*> (2, nSystems, nPIDs, nMultBins, nEtaTrkBins);

  TH1D****** h_efficiency                    = Get5DArray <TH1D*> (2, nSystems, nPIDs, nMultBins, nEtaTrkBins);


  for (short iWgts : {0, 1}) {


    TString inFileName = Form ("%s/TrackingPerformance/%s/outFile.root", rootPath.Data (), iWgts == 0 ? "Nominal" : "MCFCalWeighted");
    std::cout << "Reading " << inFileName << std::endl;
    TFile* inFile = new TFile (inFileName, "read");

    for (short iSys : systems) {

      const TString sys = (iSys == 0 ? "pp" : "pPb");

      for (short iPID = 0; iPID < nPIDs; iPID++) {

        for (short iMult = 0; iMult < nMultBins; iMult++) {

          h2_truth_tracks[iWgts][iSys][iPID][iMult] = (TH2D*) inFile->Get (Form ("h2_truth_tracks_%s_PID%i_iMult%i", sys.Data (), PIDs[iPID], iMult));
          h2_truth_tracks_wgt2[iWgts][iSys][iPID][iMult] = (TH2D*) inFile->Get (Form ("h2_truth_tracks_wgt2_%s_PID%i_iMult%i", sys.Data (), PIDs[iPID], iMult));

          if (iMult == 0 || iMult == 1 || iMult == 2) {
            h2_truth_tracks[iWgts][iSys][iPID][iMult]->RebinY (4);
            h2_truth_tracks_wgt2[iWgts][iSys][iPID][iMult]->RebinY (4);
          }

          for (short iEta = 0; iEta < nEtaTrkBins; iEta++) {

            h_truth_tracks[iWgts][iSys][iPID][iMult][iEta] = (TH1D*) inFile->Get (Form ("h_truth_tracks_%s_PID%i_iMult%i_iEta%i", sys.Data (), PIDs[iPID], iMult, iEta));
            h_truth_tracks_wgt2[iWgts][iSys][iPID][iMult][iEta] = (TH1D*) inFile->Get (Form ("h_truth_tracks_wgt2_%s_PID%i_iMult%i_iEta%i", sys.Data (), PIDs[iPID], iMult, iEta));

            if (iMult == 0 || iMult == 1 || iMult == 2) {
              h_truth_tracks[iWgts][iSys][iPID][iMult][iEta]->Rebin (4);
              h_truth_tracks_wgt2[iWgts][iSys][iPID][iMult][iEta]->Rebin (4);
            }

          } // end loop over iEta

        } // end loop over iMult

        for (short iMult = 0; iMult < nMultBins; iMult++) {

          h2_truth_matched_primary_tracks[iWgts][iSys][iPID][iMult] = (TH2D*) inFile->Get (Form ("h2_truth_matched_primary_tracks_%s_PID%i_trk_TightPrimary_iMult%i", sys.Data (), PIDs[iPID], iMult));

          if (iMult == 0 || iMult == 1 || iMult == 2)
            h2_truth_matched_primary_tracks[iWgts][iSys][iPID][iMult]->RebinY (4);

          for (short iEta = 0; iEta < nEtaTrkBins; iEta++) {

            h_truth_matched_primary_tracks[iWgts][iSys][iPID][iMult][iEta] = (TH1D*) inFile->Get (Form ("h_truth_matched_primary_tracks_%s_PID%i_trk_TightPrimary_iMult%i_iEta%i", sys.Data (), PIDs[iPID], iMult, iEta));

            if (iMult == 0 || iMult == 1 || iMult == 2)
              h_truth_matched_primary_tracks[iWgts][iSys][iPID][iMult][iEta]->Rebin (4);

          } // end loop over iEta

        } // end loop over iMult

      } // end loop over iPID

    } // end loop over iSys

  } // end loop over iWgts


  for (short iWgts : {0, 1}) {

    for (short iSys : systems) {

      const TString sys = (iSys == 0 ? "pp" : "pPb");

      for (short iPID = 0; iPID < nPIDs; iPID++) {

        for (short iMult = 0; iMult < nMultBins; iMult++) {

          h2_efficiency[iWgts][iSys][iPID][iMult] = (TH2D*) h2_truth_matched_primary_tracks[iWgts][iSys][iPID][iMult]->Clone (Form ("h2_efficiency_%s_PID%i_trk_TightPrimary_iMult%i_%s", sys.Data (), PIDs[iPID], iMult, iWgts == 0 ? "noWgts" : "withWgts"));
          BinomialDivide (h2_efficiency[iWgts][iSys][iPID][iMult], h2_truth_matched_primary_tracks[iWgts][iSys][iPID][iMult], h2_truth_tracks[iWgts][iSys][iPID][iMult], h2_truth_tracks_wgt2[iWgts][iSys][iPID][iMult]);

          for (short iEta = 0; iEta < nEtaTrkBins; iEta++) {

            h_efficiency[iWgts][iSys][iPID][iMult][iEta] = (TH1D*) h_truth_matched_primary_tracks[iWgts][iSys][iPID][iMult][iEta]->Clone (Form ("h_efficiency_%s_PID%i_trk_TightPrimary_iMult%i_iEta%i_%s", sys.Data (), PIDs[iPID], iMult, iEta, iWgts == 0 ? "noWgts" : "withWgts"));
            BinomialDivide (h_efficiency[iWgts][iSys][iPID][iMult][iEta], h_truth_matched_primary_tracks[iWgts][iSys][iPID][iMult][iEta], h_truth_tracks[iWgts][iSys][iPID][iMult][iEta], h_truth_tracks_wgt2[iWgts][iSys][iPID][iMult][iEta]);

          } // end loop over iEta

        } // end loop over iMult

      } // end loop over iPID

    } // end loop over iSys

  } // end loop over iWgts



  TLatex* tl = new TLatex ();
  //TLine* l = new TLine ();


  //for (short iSys : systems) {
  {

    //const char* cname = Form ("c_eff_sum_vs_mult_%s", iSys == 0 ? "pp" : "pPb");
    const char* cname = Form ("c_eff_sum_vs_mult");

    const short iSys = 1;

    TCanvas* c = new TCanvas (cname, "", 800*(nMultBins-1), 800);
    c->Divide (3, 1);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    for (short iMult = 0; iMult < (nMultBins-1); iMult++) {

      c->cd (iMult+1);

      gPad->SetLeftMargin (lMargin);
      gPad->SetRightMargin (rMargin);
      gPad->SetBottomMargin (bMargin);
      gPad->SetTopMargin (tMargin);

      gPad->SetLogx ();

      TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
      htemp->SetBinContent (1, 1);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("Track Reco. Efficiency");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0.5;
      const double ymax = 1.06;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (1);
      htemp->SetLineStyle (2);
  
      htemp->DrawCopy ("hist");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (32);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
      //tl->DrawLatex (0.5,  yoff, "0.5");
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

      for (short iWgts : {0, 1})
        //for (short iEta = 0; iEta < nEtaTrkBins; iEta++)
        for (short iEta : {0, 3})
          myDraw (h_efficiency[iWgts][iSys][nPIDs-1][iMult][iEta], colors[iEta], iWgts == 0 ? kOpenCircle : kOpenSquare, 1.4);

      tl->SetTextColor (kBlack);
      tl->SetTextAlign (11);
      tl->SetTextSize (28);
      tl->DrawLatexNDC (0.22, 0.22, Form ("N_{ch}^{rec} = %i-%i", (int) std::ceil (multBins[iMult]), (int) std::floor (multBins[iMult+1])));

    } // end loop over iMult

    c->cd (1);
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.88, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.84, "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.80, "TightPrimary tracks");

    c->cd (3);
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.55, 0.376, "No weights");
    tl->DrawLatexNDC (0.635, 0.376, "With weights");
    //for (short iEta = 0; iEta < nEtaTrkBins; iEta++) {
    for (short iEta : {0, 3}) {
      myLineText2 (0.70, 0.326-(iEta/3)*0.050, colors[iEta], kOpenSquare, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.4, 0.04, true);
      myLineText2 (0.60, 0.326-(iEta/3)*0.050, colors[iEta], kOpenCircle, "", 1.4, 0.040, true);
    }

    //c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencySummary_vsMultiplicity_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
    c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencySummary_withFCalWeights.pdf", workPath.Data ()));
  } // end loop over iSys



  {
    const short iMult = nMultBins-1;
    const short iSys = 1;

    const char* cname = Form ("c_eff_sum_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 800, 800);

    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogx ();

    TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
    htemp->SetBinContent (1, 1);
  
    TAxis* xax = htemp->GetXaxis ();
    TAxis* yax = htemp->GetYaxis ();

    xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
    xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
    xax->SetLabelSize (0);

    yax->SetTitle ("Track Reco. Efficiency");
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
    const double ymin = 0.5;
    const double ymax = 1.06;
    yax->SetRangeUser (ymin, ymax);

    htemp->SetLineWidth (1);
    htemp->SetLineStyle (2);
  
    htemp->DrawCopy ("hist");
    SaferDelete (&htemp);

    tl->SetTextFont (43);
    tl->SetTextSize (32);
    tl->SetTextAlign (21);
    
    const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
    tl->DrawLatex (0.4,  yoff, "0.4");
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

    for (short iWgts : {0, 1})
      for (short iEta = 0; iEta < nEtaTrkBins; iEta++)
        myDraw (h_efficiency[iWgts][iSys][nPIDs-1][iMult][iEta], colors[iEta], iWgts == 0 ? kOpenCircle : kOpenSquare, 0.8, 1, 2, "PL");

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.88, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.84, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");
    tl->DrawLatexNDC (0.22, 0.80, "TightPrimary tracks");

    myText (0.50, 0.376, kBlack, "No weights", 0.026); 
    myText (0.64, 0.376, kBlack, "With weights", 0.026); 
    for (short iEta = 0; iEta < nEtaTrkBins; iEta++) {
      myLineText2 (0.60, 0.34-iEta*0.036, colors[iEta], kOpenCircle, "",  1.0, 0.026, true);
      myLineText2 (0.74, 0.34-iEta*0.036, colors[iEta], kOpenSquare, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.0, 0.026, true);
    }

    c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencySummary_%s_withFCalWeights.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  }




  //for (short iSys : systems) {

  //  const short iEta = 0;
  //  const short iMult = nMultBins-1;

  //  const char* cname = Form ("c_eff_com_%s", iSys == 0 ? "pp" : "pPb");

  //  TCanvas* c = new TCanvas (cname, "", 800, 800);

  //  const double lMargin = 0.15;
  //  const double rMargin = 0.04;
  //  const double bMargin = 0.15;
  //  const double tMargin = 0.04;

  //  c->SetLeftMargin (lMargin);
  //  c->SetRightMargin (rMargin);
  //  c->SetBottomMargin (bMargin);
  //  c->SetTopMargin (tMargin);

  //  c->SetLogx ();

  //  TH1D* htemp = new TH1D ("htemp", "", 1, pTchBins[0], pTchBins[nPtchBins]);
  //  htemp->SetBinContent (1, 1);
  //  
  //  TAxis* xax = htemp->GetXaxis ();
  //  TAxis* yax = htemp->GetYaxis ();
  //
  //  xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
  //  xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
  //  xax->SetLabelSize (0);
  //
  //  yax->SetTitle ("Track Reco. Efficiency");
  //  yax->SetLabelFont (43);
  //  yax->SetLabelSize (32);
  //  yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
  //  const double ymin = 0.0;
  //  const double ymax = 1.25;
  //  yax->SetRangeUser (ymin, ymax);
  //
  //  htemp->SetLineWidth (1);
  //  htemp->SetLineStyle (2);
  //
  //  htemp->DrawCopy ("hist");
  //  SaferDelete (&htemp);
  //
  //  tl->SetTextFont (43);
  //  tl->SetTextSize (32);
  //  tl->SetTextAlign (21);
  //  
  //  const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
  //  tl->DrawLatex (0.4,  yoff, "0.4");
  //  tl->DrawLatex (0.7,  yoff, "0.7");
  //  tl->DrawLatex (1,  yoff, "1");
  //  tl->DrawLatex (2,  yoff, "2");
  //  tl->DrawLatex (3,  yoff, "3");
  //  tl->DrawLatex (4,  yoff, "4");
  //  tl->DrawLatex (5,  yoff, "5");
  //  tl->DrawLatex (6,  yoff, "6");
  //  tl->DrawLatex (7,  yoff, "7");
  //  tl->DrawLatex (10, yoff, "10");
  //  tl->DrawLatex (20, yoff, "20");
  //  tl->DrawLatex (30, yoff, "30");
  //  tl->DrawLatex (40, yoff, "40");
  //  tl->DrawLatex (60, yoff, "60");
  //  //tl->DrawLatex (80, yoff, "80");
  //  tl->DrawLatex (100, yoff, "100");

  //  for (short iPID = 0; iPID < nPIDs; iPID++) {
  //    TGAE* g = make_graph (h_efficiency[iWgts][iSys][iPID][iMult][iEta]);
  //    TrimGraph (g, 0.4, 130);
  //    g->Print ();
  //    //for (int i = 0; i < g->GetN(); i++) {
  //    //  double x,y;
  //    //  g->GetPoint(i,x,y);
  //    //  std::cout << "i,x,y = " << i << ", " << x << ", " << y << std::endl;
  //    //}
  //    myDraw (g, systColors[iPID], kOpenCircle, 0.8);
  //  }
  //  TGAE* g = make_graph (h_efficiency[iWgts][iSys][nPIDs-1][iMult][iEta]);
  //  TrimGraph (g, 0.4, 130);
  //    g->Print ();
  //  myDraw (g, kBlack, kFullCircle, 0.8);

  //  tl->SetTextColor (kBlack);
  //  tl->SetTextAlign (11);

  //  tl->SetTextSize (28);
  //  tl->DrawLatexNDC (0.22, 0.89, "#bf{#it{ATLAS}} Simulation Internal");
  //  tl->SetTextSize (24);
  //  tl->DrawLatexNDC (0.22, 0.85, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");
  //  tl->DrawLatexNDC (0.22, 0.81, "TightPrimary tracks");

  //  for (short iPID = 0; iPID < nPIDs-1; iPID++)
  //    myLineText2 (0.25+(iPID>=nPIDs/2 ? 0.13 : 0), 0.42-(iPID%(nPIDs/2))*0.036, systColors[iPID], kOpenCircle, partNames[iPID].c_str (), 1.2, 0.026, true);
  //  myLineText2 (0.38, 0.42-((nPIDs-1)%(nPIDs/2))*0.036, kBlack, kFullCircle, partNames[nPIDs-1].c_str (), 1.2, 0.026, true);

  //  c->SaveAs (Form ("%s/Plots/TrackingPerformance/EfficiencyComposition_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  //} // end loop over iSys


}


int main  (int argn, char** argv) {
  SetAtlasStyle ();
  PlotTrackingPerformanceWgtsComp ();
  return 0;
}

#endif
