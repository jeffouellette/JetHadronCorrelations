#ifndef __PlotTrackMomentumResolution_C__
#define __PlotTrackMomentumResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

#include <ArrayTemplates.h>
#include <Utilities.h>
#include <MyStyle.h>
#include <MyColors.h>

#include "Params.h"
#include "CentralityDefs.h"
#include "LocalUtilities.h"
#include "TrackMomentumFit.h"

using namespace JetHadronCorrelations;

typedef TGraphAsymmErrors TGAE;

// Bin edge definitions: eta, pTch, pp/pPb.
const int nFinerEtaTrkBins = 40;
const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

//const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100};
//const double pTchBins[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100};
const double pTchBins[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 70, 80, 100};
const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

const vector <int> systems = {0, 1}; // 0 = pp, 1 = pPb
const int nSystems = systems.size ();

// Whether to fit with the double Gaussian (doesn't really work very well)
bool doDoubleGaussian = false;

// The TrackMomentumFit functor returns a double-sided Crystal ball function with either a single or double Gaussian core. It has 7 or 9 free parameters.
TrackMomentumFit tmf;
// The LogTrackMomentumFit functor returns the log of a double-sided Crystal ball function with either a single or double Gaussian core. It has 7 or 9 free parameters.
LogTrackMomentumFit ltmf;



void PlotTrackMomentumResolution () {

  TFile* inFile = new TFile (Form ("%s/TrackMomentumResolution/Nominal/summary_test.root", rootPath.Data ()), "read");

  TH2D** h2_avg_tms = Get1DArray <TH2D*> (2);
  TH2D** h2_avg_tmr = Get1DArray <TH2D*> (2);

  TH1D*** h_avg_tms = Get2DArray <TH1D*> (2, nEtaTrkBins+1);
  TH1D*** h_avg_tmr = Get2DArray <TH1D*> (2, nEtaTrkBins+1);

  TH1D**** h_tmr_integratedEta = Get3DArray <TH1D*> (2, nPtchBins, nEtaTrkBins+1);
  TH2D*** h2_ptreco_pttruth_integratedEta = Get2DArray <TH2D*> (nSystems, nEtaTrkBins+1);

  TF1**** f_tmr_fits = Get3DArray <TF1*> (2, nPtchBins, nEtaTrkBins+1);


  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    h2_avg_tms[iSys] = (TH2D*) inFile->Get (Form ("h2_avg_tms_%s", sys.Data ()));
    h2_avg_tmr[iSys] = (TH2D*) inFile->Get (Form ("h2_avg_tmr_%s", sys.Data ()));

    for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {

      h_avg_tms[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_avg_tms_%s_iEta%i", sys.Data (), iEta));
      h_avg_tmr[iSys][iEta] = (TH1D*) inFile->Get (Form ("h_avg_tmr_%s_iEta%i", sys.Data (), iEta));

      h2_ptreco_pttruth_integratedEta[iSys][iEta] = (TH2D*) inFile->Get (Form ("h2_ptreco_pttruth_integratedEta_%s_iEta%i", sys.Data (), iEta));

    } // end loop over iEta

    for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {
 
      for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {

        h_tmr_integratedEta[iSys][iPtch][iEta] = (TH1D*) inFile->Get (Form ("h_tmr_integratedEta_%s_iPtch%i_%s", sys.Data (), iPtch, (iEta == nEtaTrkBins ? "allEta" : Form ("iEta%i", iEta))));

        f_tmr_fits[iSys][iPtch][iEta] = (TF1*) inFile->Get (Form ("f_tmr_%s_iPtch%i_%s_integrated", sys.Data (), iPtch, (iEta == nEtaTrkBins ? "allEta" : Form ("iEta%i", iEta))));

      } // end loop over iEta

    } // end loop over iPtch

  } // end loop over iSys



  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();


  ltmf.DoDoubleGaussian (doDoubleGaussian);
  tmf.DoDoubleGaussian (doDoubleGaussian);


  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (Form ("c_tmr_fits_%s", sys.Data ()), "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    {
      TH1D* htemp = new TH1D ("htemp", "", 1, -0.4, 0.8);
      htemp->SetBinContent (1, 1);

      gPad->SetLogy ();
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("r = #it{p}_{T}^{truth} / #it{p}_{T}^{reco} - 1");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelFont (43);
      xax->SetLabelSize (32);

      yax->SetTitle ("dN / dr");
      yax->SetLabelFont (43);
      yax->SetLabelSize (32);
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 20;
      const double ymax = 2e9;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);
    }
     

    for (int iEta : {0, 3}) {

      for (int iPtch : {10, 30, 50}) {

        h_tmr_integratedEta[iSys][iPtch][iEta]->Rebin (4);
        h_tmr_integratedEta[iSys][iPtch][iEta]->Scale (1., "width");

        myDraw (h_tmr_integratedEta[iSys][iPtch][iEta], pastels[(iPtch-10)/20], iEta == 0 ? kFullCircle : kOpenSquare, 1.0);

        TF1* f_tmr = new TF1 (Form ("f_tmr_%s_iPtch%i_iEta%i", sys.Data (), iPtch, iEta), &tmf, -1, 4, tmf.ndf ());
        tmf.CopyParams (f_tmr_fits[iSys][iPtch][iEta], f_tmr);
        f_tmr->SetLineColor (pastels[(iPtch-10)/20]);
        f_tmr->SetLineStyle (iEta == 0 ? 1 : 2);
        f_tmr->Draw ("L same");

      } // end loop over iPtch

    } // end loop over iEta

    tl->SetTextColor (kBlack);
    tl->SetTextFont (43);
    tl->SetTextAlign (11);
    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.85, iSys == 0 ? "Pythia8 #it{pp}" : "Pythia8 + #it{p}+Pb overlay");
    tl->DrawLatexNDC (0.22, 0.81, iSys == 0 ? "#sqrt{s} = 5.02 TeV" : "#sqrt{s_{NN}} = 5.02 TeV");

    for (int iPtch : {10, 30, 50})
      for (int iEta : {0, 3})
        myLineText2 (0.60, 0.84-0.04*((iPtch-10)/20)-0.12*(iEta/3), pastels[(iPtch-10)/20], iEta == 0 ? kFullCircle : kOpenSquare, Form ("(#it{p}_{T}, |#it{#eta}|) #in (%g, %g GeV) #times (%g, %g)", pTchBins[iPtch], pTchBins[iPtch+1], etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.0, 0.024, true);

    c->SaveAs (Form ("%s/Plots/TrackMomentumResolution/TMR_Fits_%s.pdf", workPath.Data (), sys.Data ()));
  }



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
      const double ymin = -0.5;
      const double ymax = 0.5;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
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

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (pTchBins[0], 100, pTchBins[nPtchBins], 100);
    }
     

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myDraw (h_avg_tms[iSys][iEta], colors[iEta], kOpenCircle, 1.0);
    myDraw (h_avg_tms[iSys][nEtaTrkBins], kBlack, kFullCircle, 1.0);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.845, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myLineText2 (0.3, 0.40-0.04*iEta, colors[iEta], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.0, 0.026, true);
    myLineText2 (0.3, 0.44, kBlack, kFullCircle, "All #it{#eta}", 1.0, 0.026, true);

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
      const double ymax = 17;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
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
    }
     

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myDraw (h_avg_tmr[iSys][iEta], colors[iEta], kOpenCircle, 1.0);
    myDraw (h_avg_tmr[iSys][nEtaTrkBins], kBlack, kFullCircle, 1.0);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.845, iSys == 0 ? "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV" : "Pythia8 + #it{p}+Pb overlay, #sqrt{s_{NN}} = 5.02 TeV");

    for (int iEta = 0; iEta < nEtaTrkBins; iEta++)
      myLineText2 (0.5, 0.70-0.04*iEta, colors[iEta], kOpenCircle, Form ("%g < |#it{#eta}| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.0, 0.026, true);
    myLineText2 (0.5, 0.74, kBlack, kFullCircle, "All #it{#eta}", 1.0, 0.026, true);

    c->SaveAs (Form ("%s/Plots/TrackMomentumResolution/TMR_%s.pdf", workPath.Data (), sys.Data ()));
  }



  for (int iSys : systems) {
    const char* cname = Form ("c_ptresponse_all_%s", iSys == 0 ? "pp" : "pPb");

    TCanvas* c = new TCanvas (cname, "", 880, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.11;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetBottomMargin (bMargin);
    c->SetTopMargin (tMargin);

    c->SetLogx ();
    c->SetLogy ();
    c->SetLogz ();

    const int iEta = nEtaTrkBins;

    TH2D* h2 = (TH2D*) h2_ptreco_pttruth_integratedEta[iSys][iEta]->Clone ("temp");

    TAxis* xax = h2->GetXaxis ();
    TAxis* yax = h2->GetYaxis ();
    TAxis* zax = h2->GetZaxis ();

    // normalize distributions in response (along y-axis)
    for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
      double integral = 0;
      for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
        integral += h2->GetBinContent (iX, iY) * yax->GetBinWidth (iY);
      if (integral > 0) {
        for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
          h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / integral);
          h2->SetBinError (iX, iY, h2->GetBinError (iX, iY) / integral);
        } // end loop over iY
      }
    } // end loop over iX

    xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
    yax->SetTitle ("#it{p}_{T}^{reco} [GeV]");
    zax->SetTitle ("");

    const double zmin = 1e-5;
    const double zmax = 1;

    h2->SetLineWidth (0);

    xax->SetMoreLogLabels ();
    yax->SetMoreLogLabels ();
    zax->SetRangeUser (zmin, zmax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (32);
    yax->SetTitleFont (43);
    yax->SetTitleSize (32);
    zax->SetTitleFont (43);
    zax->SetTitleSize (32);
    zax->SetTitleOffset (1.1*zax->GetTitleOffset ());
    xax->SetLabelFont (43);
    xax->SetLabelSize (32);
    yax->SetLabelFont (43);
    yax->SetLabelSize (32);
    zax->SetLabelFont (43);
    zax->SetLabelSize (24);

    h2->DrawCopy ("colz");
    SaferDelete (&h2);

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (12);

    tl->SetTextSize (28);
    tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (24);
    tl->DrawLatexNDC (0.22, 0.850, iSys == 0 ? "Pythia8 #it{pp}" : "Pythia8 + #it{p}+Pb overlay");
    tl->DrawLatexNDC (0.22, 0.810, iSys == 0 ? "#sqrt{s} = 5.02 TeV" : "#sqrt{s_{NN}} = 5.02 TeV");

    c->SaveAs (Form ("%s/Plots/TrackMomentumResolution/FullMomentumResponse_%s.pdf", workPath.Data (), iSys == 0 ? "pp" : "pPb"));
  }


}

#endif
