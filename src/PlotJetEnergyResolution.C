#ifndef __PlotJetEnergyResolution_C__
#define __PlotJetEnergyResolution_C__

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

using namespace JetHadronCorrelations;

typedef TGraphAsymmErrors TGAE;


void PlotJetEnergyResolution () {

  const int nFinerEtaBins = 56;
  const double* finerEtaBins = linspace (-2.8, 2.8, nFinerEtaBins);
  
  const double etaBins[] = {-2.8, -2.1, -1.2, -0.8, -0.3, 0, 0.3, 0.8, 1.2, 2.1, 2.8};
  const int nEtaBins = sizeof (etaBins) / sizeof (etaBins[0]) - 1;
  
  const double pTJBins[] = {10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 360, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300};
  const int nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;
  
  const vector <int> systems = {0, 1}; // 0 = pp, 1 = pPb
  const int nSystems = systems.size ();
  
  const std::vector <JetRadius> radii = {JetRadius::R0p2, JetRadius::R0p4};
  const int nRadii = radii.size ();


  TFile* inFile = new TFile (Form ("%s/JetEnergyResolution/Nominal/summary.root", rootPath.Data ()), "read");

  TH2D*** h2_avg_jpts = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_avg_jptr = Get2DArray <TH2D*> (nSystems, nRadii);
  TH1D**** h_avg_jpts = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);
  TH1D**** h_avg_jptr = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);

  TH2D*** h2_avg_jes  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_avg_jer  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH1D**** h_avg_jes  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);
  TH1D**** h_avg_jer  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);

  TH2D*** h2_avg_jetacorr = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_avg_jetares  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH1D**** h_avg_jetacorr = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);
  TH1D**** h_avg_jetares  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);

  TH2D**** h2_jpts_integratedEta      = Get3DArray <TH2D*> (nSystems, nRadii, nEtaBins+1);
  TH2D**** h2_jes_integratedEta       = Get3DArray <TH2D*> (nSystems, nRadii, nEtaBins+1);
  TH2D**** h2_jetacorr_integratedEta  = Get3DArray <TH2D*> (nSystems, nRadii, nEtaBins+1);

  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iR = 0; iR < nRadii; iR++) {

      const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

      h2_avg_jpts[iSys][iR] = (TH2D*) inFile->Get (Form ("h2_r%i_avg_jpts_%s", r, sys.Data ()));
      h2_avg_jptr[iSys][iR] = (TH2D*) inFile->Get (Form ("h2_r2%i_avg_jptr_%s", r, sys.Data ()));

      h2_avg_jes[iSys][iR] = (TH2D*) inFile->Get (Form ("h2_r%i_avg_jes_%s", r, sys.Data ()));
      h2_avg_jer[iSys][iR] = (TH2D*) inFile->Get (Form ("h2_r%i_avg_jer_%s", r, sys.Data ()));

      h2_avg_jetacorr[iSys][iR] = (TH2D*) inFile->Get (Form ("h2_r%i_avg_jetacorr_%s", r, sys.Data ()));
      h2_avg_jetares[iSys][iR] = (TH2D*) inFile->Get (Form ("h2_r%i_avg_jetares_%s", r, sys.Data ()));

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        h_avg_jpts[iSys][iR][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_avg_jpts_%s_iEta%i", r, sys.Data (), iEta));
        h_avg_jptr[iSys][iR][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_avg_jptr_%s_iEta%i", r, sys.Data (), iEta));

        h_avg_jes[iSys][iR][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_avg_jes_%s_iEta%i", r, sys.Data (), iEta));
        h_avg_jer[iSys][iR][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_avg_jer_%s_iEta%i", r, sys.Data (), iEta));

        h_avg_jetacorr[iSys][iR][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_avg_jetacorr_%s_iEta%i", r, sys.Data (), iEta));
        h_avg_jetares[iSys][iR][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_avg_jetares_%s_iEta%i", r, sys.Data (), iEta));


        h2_jpts_integratedEta[iSys][iR][iEta] = (TH2D*) inFile->Get (Form ("h2_r%i_jpts_integratedEta_%s_iEta%i", r, sys.Data (), iEta));

        h2_jes_integratedEta[iSys][iR][iEta] = (TH2D*) inFile->Get (Form ("h2_r%i_jes_integratedEta_%s_iEta%i", r, sys.Data (), iEta));

        h2_jetacorr_integratedEta[iSys][iR][iEta] = (TH2D*) inFile->Get (Form ("h2_r%i_jetacorr_integratedEta_%s_iEta%i", r, sys.Data (), iEta));

      } // end loop over iEta

    } // end loop over iR

  } // end loop over iSys



  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();


  for (int iR = 0; iR < nRadii; iR++) {

    const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

    for (int iSys : systems) {

      const TString sys = (iSys == 0 ? "pp" : "pPb");

      // plot jet pT scale
      {
        TCanvas* c = new TCanvas (Form ("c_r%i_jpts_%s", r, sys.Data ()), "", 800, 800);
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

          TH1D* htemp = new TH1D ("htemp", "", 1, 20, pTJBins[nPtJBins]);
      
          TAxis* xax = htemp->GetXaxis ();
          TAxis* yax = htemp->GetYaxis ();

          xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
          xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
          xax->SetLabelSize (0);

          yax->SetTitle ("#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]");
          yax->SetLabelFont (43);
          yax->SetLabelSize (32);
          yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
          const double ymin = 92;
          const double ymax = 108;
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
          tl->DrawLatex (20, yoff, "20");
          tl->DrawLatex (30, yoff, "30");
          tl->DrawLatex (40, yoff, "40");
          tl->DrawLatex (50, yoff, "50");
          tl->DrawLatex (70, yoff, "70");
          tl->DrawLatex (100, yoff, "100");
          tl->DrawLatex (200, yoff, "200");
          tl->DrawLatex (300, yoff, "300");
          tl->DrawLatex (500, yoff, "500");
          tl->DrawLatex (800, yoff, "800");

          l->SetLineStyle (2);
          l->SetLineWidth (2);
          //l->SetLineColor (kPink-8);
          l->DrawLine (20, 100, pTJBins[nPtJBins], 100);
        }

        for (int iEta = 0; iEta < nEtaBins; iEta++)
          myDraw (h_avg_jpts[iSys][iR][iEta], systColors[iEta], kOpenCircle, 1.0);
        myDraw (h_avg_jpts[iSys][iR][nEtaBins], kBlack, kFullCircle, 1.0);

        tl->SetTextColor (kBlack);
        tl->SetTextAlign (12);

        tl->SetTextSize (26);
        tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
        if (iSys == 0)
          tl->DrawLatexNDC (0.22, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
        else if (iSys == 1)
          tl->DrawLatexNDC (0.22, 0.845, "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
        tl->DrawLatexNDC (0.22, 0.810, Form ("Anti-k_{T} HI Jets, R=0.%i", r));
        for (int iEta = 0; iEta < nEtaBins; iEta++)
          myLineText2 (0.3+(2*iEta>=nEtaBins ? 0.24 : 0), 0.328-0.032*(iEta%(nEtaBins/2)), systColors[iEta], kOpenCircle, Form ("%g < #it{#eta} < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.026, true);
        myLineText2 (0.3, 0.36, kBlack, kFullCircle, "|#it{#eta}| < 2.8", 1.2, 0.026, true);

        c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JPTS_R%i_%s.pdf", workPath.Data (), r, sys.Data ()));
      } // end loop over iSys


      // plot jet energy scale
      {

        TCanvas* c = new TCanvas (Form ("c_r%i_jes_%s", r, sys.Data ()), "", 800, 800);
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

          TH1D* htemp = new TH1D ("htemp", "", 1, 20, pTJBins[nPtJBins]);
      
          TAxis* xax = htemp->GetXaxis ();
          TAxis* yax = htemp->GetYaxis ();

          xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
          xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
          xax->SetLabelSize (0);

          yax->SetTitle ("#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]");
          yax->SetLabelFont (43);
          yax->SetLabelSize (32);
          yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
          const double ymin = 92;
          const double ymax = 108;
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
          tl->DrawLatex (20, yoff, "20");
          tl->DrawLatex (30, yoff, "30");
          tl->DrawLatex (40, yoff, "40");
          tl->DrawLatex (50, yoff, "50");
          tl->DrawLatex (70, yoff, "70");
          tl->DrawLatex (100, yoff, "100");
          tl->DrawLatex (200, yoff, "200");
          tl->DrawLatex (300, yoff, "300");
          tl->DrawLatex (500, yoff, "500");
          tl->DrawLatex (800, yoff, "800");

          l->SetLineStyle (2);
          l->SetLineWidth (2);
          //l->SetLineColor (kPink-8);
          l->DrawLine (20, 100, pTJBins[nPtJBins], 100);
        }

        for (int iEta = 0; iEta < nEtaBins; iEta++)
          myDraw (h_avg_jes[iSys][iR][iEta], systColors[iEta], kOpenCircle, 1.0);
        myDraw (h_avg_jes[iSys][iR][nEtaBins], kBlack, kFullCircle, 1.0);

        tl->SetTextColor (kBlack);
        tl->SetTextAlign (12);

        tl->SetTextSize (26);
        tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
        if (iSys == 0)
          tl->DrawLatexNDC (0.22, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
        else if (iSys == 1)
          tl->DrawLatexNDC (0.22, 0.845, "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
        tl->DrawLatexNDC (0.22, 0.810, Form ("Anti-k_{T} HI Jets, R=0.%i", r));
        for (int iEta = 0; iEta < nEtaBins; iEta++)
          myLineText2 (0.3+(2*iEta>=nEtaBins ? 0.24 : 0), 0.328-0.032*(iEta%(nEtaBins/2)), systColors[iEta], kOpenCircle, Form ("%g < #it{#eta} < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.026, true);
        myLineText2 (0.3, 0.36, kBlack, kFullCircle, "|#it{#eta}| < 2.8", 1.2, 0.026, true);

        c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JES_R%i_%s.pdf", workPath.Data (), r, sys.Data ()));
      }


      // plot jet pT resolution
      {
        TCanvas* c = new TCanvas (Form ("c_r%i_jptr_%s", r, sys.Data ()), "", 800, 800);
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

          TH1D* htemp = new TH1D ("htemp", "", 1, 20, pTJBins[nPtJBins]);
    
          TAxis* xax = htemp->GetXaxis ();
          TAxis* yax = htemp->GetYaxis ();

          xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
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
          tl->DrawLatex (20, yoff, "20");
          tl->DrawLatex (30, yoff, "30");
          tl->DrawLatex (40, yoff, "40");
          tl->DrawLatex (50, yoff, "50");
          tl->DrawLatex (70, yoff, "70");
          tl->DrawLatex (100, yoff, "100");
          tl->DrawLatex (200, yoff, "200");
          tl->DrawLatex (300, yoff, "300");
          tl->DrawLatex (500, yoff, "500");
          tl->DrawLatex (800, yoff, "800");
        }

        for (int iEta = 0; iEta < nEtaBins; iEta++)
          myDraw (h_avg_jptr[iSys][iR][iEta], systColors[iEta], kOpenCircle, 1.0);
        myDraw (h_avg_jptr[iSys][iR][nEtaBins], kBlack, kFullCircle, 1.0);

        tl->SetTextColor (kBlack);
        tl->SetTextAlign (12);

        tl->SetTextSize (26);
        tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
        if (iSys == 0)
          tl->DrawLatexNDC (0.22, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
        else if (iSys == 1)
          tl->DrawLatexNDC (0.22, 0.845, "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
        tl->DrawLatexNDC (0.22, 0.810, Form ("Anti-k_{T} HI Jets, R=0.%i", r));
        for (int iEta = 0; iEta < nEtaBins; iEta++)
          myLineText2 (0.5+(2*iEta>=nEtaBins ? 0.24 : 0), 0.718-0.032*(iEta%(nEtaBins/2)), systColors[iEta], kOpenCircle, Form ("%g < #it{#eta} < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.026, true);
        myLineText2 (0.5, 0.75, kBlack, kFullCircle, "|#it{#eta}| < 2.8", 1.2, 0.026, true);

        c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JPTR_R%i_%s.pdf", workPath.Data (), r, sys.Data ()));
      }



      // plot jet energy resolution
      {
        TCanvas* c = new TCanvas (Form ("c_r%i_jer_%s", r, sys.Data ()), "", 800, 800);
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

          TH1D* htemp = new TH1D ("htemp", "", 1, 20, pTJBins[nPtJBins]);
    
          TAxis* xax = htemp->GetXaxis ();
          TAxis* yax = htemp->GetYaxis ();

          xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
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
          tl->DrawLatex (20, yoff, "20");
          tl->DrawLatex (30, yoff, "30");
          tl->DrawLatex (40, yoff, "40");
          tl->DrawLatex (50, yoff, "50");
          tl->DrawLatex (70, yoff, "70");
          tl->DrawLatex (100, yoff, "100");
          tl->DrawLatex (200, yoff, "200");
          tl->DrawLatex (300, yoff, "300");
          tl->DrawLatex (500, yoff, "500");
          tl->DrawLatex (800, yoff, "800");
        }

        for (int iEta = 0; iEta < nEtaBins; iEta++)
          myDraw (h_avg_jer[iSys][iR][iEta], systColors[iEta], kOpenCircle, 1.0);
        myDraw (h_avg_jer[iSys][iR][nEtaBins], kBlack, kFullCircle, 1.0);

        tl->SetTextColor (kBlack);
        tl->SetTextAlign (12);

        tl->SetTextSize (26);
        tl->DrawLatexNDC (0.22, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
        if (iSys == 0)
          tl->DrawLatexNDC (0.22, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
        else if (iSys == 1)
          tl->DrawLatexNDC (0.22, 0.845, "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
        tl->DrawLatexNDC (0.22, 0.810, Form ("Anti-k_{T} HI Jets, R=0.%i", r));
        for (int iEta = 0; iEta < nEtaBins; iEta++)
          myLineText2 (0.5+(2*iEta>=nEtaBins ? 0.24 : 0), 0.718-0.032*(iEta%(nEtaBins/2)), systColors[iEta], kOpenCircle, Form ("%g < #it{#eta} < %g", etaBins[iEta], etaBins[iEta+1]), 1.2, 0.026, true);
        myLineText2 (0.5, 0.75, kBlack, kFullCircle, "|#it{#eta}| < 2.8", 1.2, 0.026, true);

        c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JER_R%i_%s.pdf", workPath.Data (), r, sys.Data ()));
      }


      for (int iEta = 0; iEta <= nEtaBins; iEta++) {
        const char* cname = Form ("c_jetResp_r%i_%s_iEta%i", r, sys.Data (), iEta);

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
        c->SetLogz ();

        TH2D* h2 = (TH2D*) h2_jpts_integratedEta[iSys][iR][iEta]->Clone ("temp");

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

        const double ymin = 0;
        const double ymax = 2.4;

        xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
        yax->SetTitle ("#it{p}_{T}^{reco} / #it{p}_{T}^{truth}");
        //zax->SetTitle ("(1 / N(#it{p}_{T}^{truth})) (dN / dr_{jet})");
        zax->SetTitle ("");

        //const double zmin = 0.8;
        //const double zmax = 1;

        h2->SetLineWidth (0);

        //zax->SetRangeUser (zmin, zmax);
        //xax->SetMoreLogLabels ();

        xax->SetTitleFont (43);
        xax->SetTitleSize (32);
        yax->SetTitleFont (43);
        yax->SetTitleSize (32);
        zax->SetTitleFont (43);
        zax->SetTitleSize (32);
        zax->SetTitleOffset (1.1*zax->GetTitleOffset ());
        xax->SetLabelFont (43);
        //xax->SetLabelSize (32);
        xax->SetLabelSize (0);
        yax->SetLabelFont (43);
        yax->SetLabelSize (32);
        zax->SetLabelFont (43);
        zax->SetLabelSize (24);

        h2->DrawCopy ("colz");
        SaferDelete (&h2);

        tl->SetTextFont (43);
        tl->SetTextSize (32);
        tl->SetTextAlign (21);
        
        const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
        tl->DrawLatex (20, yoff, "20");
        tl->DrawLatex (30, yoff, "30");
        tl->DrawLatex (40, yoff, "40");
        tl->DrawLatex (60, yoff, "60");
        tl->DrawLatex (100, yoff, "100");
        tl->DrawLatex (200, yoff, "200");
        tl->DrawLatex (400, yoff, "400");
        tl->DrawLatex (800, yoff, "800");

        tl->SetTextColor (kBlack);
        tl->SetTextAlign (12);

        tl->SetTextSize (26);
        tl->DrawLatexNDC (0.36, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
        if (iSys == 0)
          tl->DrawLatexNDC (0.36, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
        else if (iSys == 1)
          tl->DrawLatexNDC (0.36, 0.845, "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
        tl->DrawLatexNDC (0.36, 0.810, Form ("Anti-k_{T} HI Jets, R=0.%i", r));
        tl->DrawLatexNDC (0.36, 0.775, iEta == nEtaBins ? "|#it{#eta}| < 2.8" : Form ("%g < #it{#eta} < %g", etaBins[iEta], etaBins[iEta+1]));

        c->SaveAs (Form ("%s/Plots/JetEnergyResolution/JetResponse_R%i_%s_iEta%i.pdf", workPath.Data (), r, iSys == 0 ? "pp" : "pPb", iEta));
      }

    } // end loop over sys

  } // end loop over iR


}

#endif
