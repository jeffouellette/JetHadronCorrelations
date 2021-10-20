#ifndef __JetHadronCorrelatorPlotJets_C__
#define __JetHadronCorrelatorPlotJets_C__

#include <iostream>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLorentzVector.h>

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


const bool makeMCClosurePlots = true;


void PlotJets (const char* tag, const char* inFileTag) {

  TFile* inFile = nullptr;

  TH1D***  h_evt_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  TH1D***  h_jet_counts_ref     = Get2DArray <TH1D*> (2, nVar);
  TH1D**** h_evt_counts         = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH1D**** h_jet_counts         = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);

  TH1D*** h_jet_pt_ref          = Get2DArray <TH1D*> (2, nVar);
  TH2D*** h2_jet_pt_cov_ref     = Get2DArray <TH2D*> (2, nVar);

  TH1D**** h_jet_pt             = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);
  TH2D**** h2_jet_pt_cov        = Get3DArray <TH2D*> (2, nZdcCentBins+1, nVar);

  TH1D**** h_jet_pt_ratio       = Get3DArray <TH1D*> (2, nZdcCentBins+1, nVar);

  TH2D***  h2_jet_eta_phi_ref   = Get2DArray <TH2D*> (2, nVar);
  TH2D**** h2_jet_eta_phi       = Get3DArray <TH2D*> (2, nZdcCentBins+1, nVar);

  TGAE**  g_jet_pt_ref_syst     = Get1DArray <TGAE*> (nVar);
  TGAE*** g_jet_pt_syst         = Get2DArray <TGAE*> (nZdcCentBins+1, nVar);

  TGAE*** g_jet_pt_ratio_syst   = Get2DArray <TGAE*> (nZdcCentBins+1, nVar);


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/ProcessJets_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    for (int iDType = 0; iDType < 2; iDType++) {

      const TString dType = (iDType == 0 ? "data" : "mc");

      for (int iVar = 0; iVar < nVar; iVar++) {
  
        const TString var = variations[iVar];
  
        if ((iDType == 0 && dataVariations.count (var) == 0) || (iDType == 1 && mcVariations.count (var) + otherMCVariations.count (var) == 0))
          continue;
  
        h_evt_counts_ref[iDType][iVar]    = (TH1D*) inFile->Get (Form ("h_evt_counts_ref_%s_%s",    dType.Data (), var.Data ()));
        h_jet_counts_ref[iDType][iVar]    = (TH1D*) inFile->Get (Form ("h_jet_counts_ref_%s_%s",    dType.Data (), var.Data ()));
  
        h_jet_pt_ref[iDType][iVar]        = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_%s_%s",        dType.Data (), var.Data ()));
  
        h2_jet_eta_phi_ref[iDType][iVar]  = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_ref_%s_%s",  dType.Data (), var.Data ()));
  
        for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

          const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));
  
          h_evt_counts[iDType][iCent][iVar]   = (TH1D*) inFile->Get (Form ("h_evt_counts_pPb_%s_%s_%s",   cent.Data (), dType.Data (), var.Data ()));
          h_jet_counts[iDType][iCent][iVar]   = (TH1D*) inFile->Get (Form ("h_jet_counts_pPb_%s_%s_%s",   cent.Data (), dType.Data (), var.Data ()));
  
          h_jet_pt[iDType][iCent][iVar]       = (TH1D*) inFile->Get (Form ("h_jet_pt_pPb_%s_%s_%s",       cent.Data (), dType.Data (), var.Data ()));
          h_jet_pt_ratio[iDType][iCent][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_ratio_%s_%s_%s",     cent.Data (), dType.Data (), var.Data ()));
  
          h2_jet_eta_phi[iDType][iCent][iVar] = (TH2D*) inFile->Get (Form ("h2_jet_eta_phi_pPb_%s_%s_%s", cent.Data (), dType.Data (), var.Data ()));
  
        } // end loop over iCent

      } // end loop over iVar

    } // end loop over iDType


    for (int iVar = 0; iVar < nVar; iVar++) {

      const TString var = variations[iVar];

      g_jet_pt_ref_syst[iVar] = (TGAE*) inFile->Get (Form ("g_jet_pt_ref_syst_%s", var.Data ()));

      for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

        const TString cent = (iCent == nZdcCentBins ? "allCent" : Form ("iCent%i", iCent));

        g_jet_pt_syst[iCent][iVar]        = (TGAE*) inFile->Get (Form ("g_jet_pt_syst_pPb_%s_%s",   cent.Data (), var.Data ()));
        g_jet_pt_ratio_syst[iCent][iVar]  = (TGAE*) inFile->Get (Form ("g_jet_pt_ratio_syst_%s_%s", cent.Data (), var.Data ()));

      } // end loop over iCent

    } // end loop over iVar

  }



  {
    std::cout << "---------------" << std::endl << "JETS IN DATA" << std::endl << "---------------" << std::endl;
    std::cout << "Number of pp jets: " << h_jet_counts_ref[0][0]->GetBinContent (1) << std::endl;
    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++)
      std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << h_jet_counts[0][iCent][0]->GetBinContent (1) << std::endl;

    std::cout << "Formatted for latex:" << std::endl;
    std::cout << (int) h_jet_counts_ref[0][0]->GetBinContent (1);
    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++)
      std::cout << " & " << (int) h_jet_counts[0][iCent][0]->GetBinContent (1);
    std::cout << std::endl << std::endl;


    const float trigpt = (strcmp (tag, "30GeVJets") == 0 ? 30 : 60);
    float integral = 0;
    for (int iX = h_jet_pt_ref[0][0]->FindBin (trigpt); iX <= h_jet_pt_ref[0][0]->GetNbinsX (); iX++)
      integral += h_jet_pt_ref[0][0]->GetBinContent (iX) * h_jet_pt_ref[0][0]->GetBinWidth (iX);
    std::cout << "Average number of jets per trigger jet in pp: " << integral << std::endl;
    integral = 0;
    for (int iX = h_jet_pt[0][nZdcCentBins][0]->FindBin (trigpt); iX <= h_jet_pt[0][nZdcCentBins][0]->GetNbinsX (); iX++)
      integral += h_jet_pt[0][nZdcCentBins][0]->GetBinContent (iX) * h_jet_pt[0][nZdcCentBins][0]->GetBinWidth (iX);
    std::cout << "Average number of jets per trigger jet in p+Pb: " << integral << std::endl;
    std::cout << std::endl << std::endl;


    std::cout << "---------------" << std::endl << "JETS IN MC" << std::endl << "---------------" << std::endl;
    std::cout << "Number of pp jets: " << h_jet_counts_ref[1][0]->GetBinContent (1) << std::endl;
    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++)
      std::cout << "Number of p+Pb " << zdcCentPercs[iCent+1] << "\%-" << zdcCentPercs[iCent] << "\% jets: " << h_jet_counts[1][iCent][0]->GetBinContent (1) << std::endl;

    std::cout << "Formatted for latex:" << std::endl;
    std::cout << (int) h_jet_counts_ref[1][0]->GetBinContent (1);
    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++)
      std::cout << " & " << (int) h_jet_counts[1][iCent][0]->GetBinContent (1);
    std::cout << std::endl;


    integral = 0;
    for (int iX = h_jet_pt_ref[1][0]->FindBin (trigpt); iX <= h_jet_pt_ref[1][0]->GetNbinsX (); iX++)
      integral += h_jet_pt_ref[1][0]->GetBinContent (iX) * h_jet_pt_ref[1][0]->GetBinWidth (iX);
    std::cout << "Average number of jets per trigger jet in pp: " << integral << std::endl;
    integral = 0;
    for (int iX = h_jet_pt[1][nZdcCentBins][0]->FindBin (trigpt); iX <= h_jet_pt[1][nZdcCentBins][0]->GetNbinsX (); iX++)
      integral += h_jet_pt[1][nZdcCentBins][0]->GetBinContent (iX) * h_jet_pt[1][nZdcCentBins][0]->GetBinWidth (iX);
    std::cout << "Average number of jets per trigger jet in p+Pb: " << integral << std::endl;
    std::cout << std::endl << std::endl;
  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  for (int iDType = 0; iDType < 2; iDType++) {
  //for (int iDType = 0; iDType < 1; iDType++) {

    const char* canvasName = Form ("c_jet_pt");

    TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
    c->cd ();

    const double fPad = 600./1000.;
    TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0, 1-fPad, 1, 1.0);
    TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0, 0.0, 1, 1-fPad);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.20);

    uPad->SetLeftMargin (0.12);
    dPad->SetLeftMargin (0.12);
    uPad->SetRightMargin (0.03);
    dPad->SetRightMargin (0.03);

    uPad->Draw ();
    dPad->Draw ();

    TH1D* h = nullptr;
    TGAE* g = nullptr;

    double ymin=2e-9;
    double ymax=1e3;

    uPad->cd ();
    uPad->SetLogx();
    uPad->SetLogy ();

    const double maxx = (strcmp (tag, "30GeVJets") == 0 ? 90 : 200);
    h = new TH1D ("htemp", ";#it{p}_{T}^{jet} [GeV];(1 / N_{trig. jet}) (dN_{jet} / d#it{p}_{T}^{jet}) [GeV^{-1}]", 1, pTJBins[0], maxx);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.028/fPad);
    h->GetXaxis ()->SetLabelSize (0.028/fPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fPad);
    h->GetYaxis ()->SetTitleSize (0.028/fPad);
    h->GetYaxis ()->SetLabelSize (0.028/fPad);
    h->GetYaxis ()->SetTitleOffset (2.0*fPad);

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = (TH1D*) h_jet_pt_ref[iDType][0]->Clone ("htemp");
    h->Scale (std::pow (10, 3));
    myDrawHist (h, kBlack, 1, 2);
    SaferDelete (&h);

    h = h_jet_pt_ref[iDType][0];
    //g = (TGAE*) g_jet_pt_ref_syst[0]->Clone ();
    //SetCentralValuesKeepRelativeErrors (g, h);
    //ScaleGraph (g, nullptr, std::pow (10, 3));
    //myDrawSystFill (g, colorfulSystColors[0], 1, 1001);
    //SaferDelete (&g);
   
    g = make_graph (h);
    //ResetXErrors (g);
    ScaleGraph (g, nullptr, std::pow (10, 3));
    myDraw (g, colorfulColors[0], kFullCircle, 1.4, 1, 3, "P", false);
    SaferDelete (&g);

    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      h = (TH1D*) h_jet_pt_ref[iDType][0]->Clone ("htemp");
      h->Scale (std::pow (10, 2-iCent));
      myDrawHist (h, kBlack, 1, 2);
      SaferDelete (&h);

      h = h_jet_pt[iDType][iCent][0];
      //g = (TGAE*) g_jet_pt_syst[iCent][0]->Clone ();
      //SetCentralValuesKeepRelativeErrors (g, h);
      //ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
      //myDrawSystFill (g, colorfulSystColors[iCent+1], 1, 1001);
      //SaferDelete (&g);

      g = make_graph (h);
      //ResetXErrors (g);
      ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
      myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.4, 1, 3, "P", false);
      SaferDelete (&g);

    } // end loop over iCent

    myText (0.18, 0.89, kBlack, Form ("#bf{#it{ATLAS}} %sInternal", iDType == 0 ? "" : "Simulation "), 0.038);
    myText (0.64, 0.89, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.034);
    myText (0.64, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.034);
    myText (0.64, 0.79, kBlack, Form ("#bf{#it{p}_{T}^{leading} > %s GeV}", strcmp (tag, "30GeVJets") == 0 ? "30" : "60"), 0.034);
    //if (iDType == 0)
    //  myText (0.64, 0.74, kBlack, strcmp (tag, "30GeVJets") == 0 ? "HLT_mb_sptrk_L1MBTS_1" : "HLT_j50_ion_L1J10", 0.034);

    mySimpleMarkerAndBoxAndLineText (0.27, 0.16, 1.4, 1001, colorfulSystColors[0], 1.0, colorfulColors[0], kFullCircle, 1.6, "#it{pp} (#times10^{3})", 0.022/fPad);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      mySimpleMarkerAndBoxAndLineText (0.27 + (iCent >= 2 ? 0.4 : 0), 0.16-((iCent+1)%3)*0.04, 1.4, 1001, colorfulSystColors[iCent+1], 1.0, colorfulColors[iCent+1], kFullCircle, 1.6, Form ("%s %i-%i%% (#times10^{%i})", iDType == 0 ? "ZDC" : "FCal", zdcCentPercs[iCent+1], zdcCentPercs[iCent], 2-iCent), 0.022/fPad);
    } // end loop over iCent
    mySimpleMarkerAndBoxAndLineText (0.67, 0.04, 1.4, 1001, colorfulSystColors[nZdcCentBins+1], 1.0, colorfulColors[nZdcCentBins+1], kFullCircle, 1.6, Form ("All cent. (#times10^{%i})", 2-nZdcCentBins), 0.022/fPad);
    mySimpleMarkerAndBoxAndLineText (0.27, 0.04, 1.4, 1001, kWhite, 0.0, kBlack, kDot, 0.0, "#it{pp} (#it{scaled})", 0.022/fPad);

    uPad->RedrawAxis();


    dPad->cd ();
    dPad->SetLogx();

    ymin=0.00;
    ymax=2.00;

    h = new TH1D ("htemp", ";#it{p}_{T}^{jet} [GeV];#it{p}+Pb / #it{pp} + Const.", 1, pTJBins[0], maxx);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.028/(1-fPad));
    h->GetXaxis ()->SetLabelSize (0.028/(1-fPad));
    h->GetXaxis ()->SetTitleOffset (3.2*(1-fPad));
    h->GetYaxis ()->SetTitleSize (0.028/(1-fPad));
    h->GetYaxis ()->SetLabelSize (0.028/(1-fPad));
    h->GetYaxis ()->SetTitleOffset (2.0*(1-fPad));
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    double x, y;
    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const double offset = 0.25*(2-iCent);

      h = h_jet_pt_ratio[iDType][iCent][0];
      //g = (TGAE*) g_jet_pt_ratio_syst[iCent][0]->Clone ();
      //SetCentralValuesKeepRelativeErrors (g, h);
      //for (int i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  g->SetPoint (i, x, y+offset);
      //}
      //myDrawSystFill (g, colorfulSystColors[iCent+1], 1, 1001);
      //SaferDelete (&g);

      g = make_graph (h);
      //ResetXErrors (g);
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        g->SetPoint (i, x, y+offset);
      }
      myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.4, 1, 3, "P", false);
      SaferDelete (&g);

      l->SetLineWidth (3);
      l->SetLineColor (colorfulColors[iCent+1]);
      l->DrawLine (pTJBins[0], 1+offset, maxx, 1+offset);

      tl->SetTextAlign (offset > 0 ? 31 : 33);
      tl->SetTextFont (43);
      tl->SetTextSize (20);
      tl->SetTextColor (colorfulColors[iCent+1]);
      tl->DrawLatex (pTJBins[0] * std::exp (0.97 * std::log (maxx/pTJBins[0])), 1+offset+(offset > 0 ? 0.015:-0.015), Form ("#bf{#it{%s%g}}", offset == 0 ? "" : (offset > 0 ? "+ " : "#minus "), std::fabs (offset)));

    } // end loop over iCent

    dPad->RedrawAxis();

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtSpectrum_%s%s.pdf", workPath.Data (), tag, iDType == 1 ? "_mc" : ""));

  }



  if (makeMCClosurePlots) {
    const int iMCTruthLevel = GetVarN ("MCTruthLevel");

    const char* canvasName = Form ("c_jet_pt_mcclosure");

    TCanvas* c = new TCanvas (canvasName, "", 800, 1000);
    c->cd ();

    const double fPad = 600./1000.;
    TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0, 1-fPad, 1, 1.0);
    TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0, 0.0, 1, 1-fPad);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.20);

    uPad->SetLeftMargin (0.12);
    dPad->SetLeftMargin (0.12);
    uPad->SetRightMargin (0.03);
    dPad->SetRightMargin (0.03);

    uPad->Draw ();
    dPad->Draw ();

    TH1D* h = nullptr;
    TGAE* g = nullptr;

    double ymin=2e-9;
    double ymax=1e3;

    uPad->cd ();
    uPad->SetLogx();
    uPad->SetLogy ();

    const double maxx = (strcmp (tag, "30GeVJets") == 0 ? 90 : 200);
    h = new TH1D ("htemp", ";#it{p}_{T}^{jet} [GeV];(1 / N_{trig. jet}) (dN_{jet} / d#it{p}_{T}^{jet}) [GeV^{-1}]", 1, pTJBins[0], maxx);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.028/fPad);
    h->GetXaxis ()->SetLabelSize (0.028/fPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fPad);
    h->GetYaxis ()->SetTitleSize (0.028/fPad);
    h->GetYaxis ()->SetLabelSize (0.028/fPad);
    h->GetYaxis ()->SetTitleOffset (1.8*fPad);

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    h = (TH1D*) h_jet_pt_ref[1][iMCTruthLevel]->Clone ("htemp");
    h->Scale (std::pow (10, 3));
    myDrawHist (h, kBlack);
    SaferDelete (&h);

    h = h_jet_pt_ref[1][0];
    //g = (TGAE*) g_jet_pt_ref_syst[0]->Clone ();
    //SetCentralValuesKeepRelativeErrors (g, h);
    //ScaleGraph (g, nullptr, std::pow (10, 3));
    //myDrawSystFill (g, colorfulSystColors[0], 1, 1001);
    //SaferDelete (&g);
   
    g = make_graph (h);
    //ResetXErrors (g);
    ScaleGraph (g, nullptr, std::pow (10, 3));
    myDraw (g, colorfulColors[0], kFullCircle, 1.4, 1, 3, "P", false);
    SaferDelete (&g);

    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      h = (TH1D*) h_jet_pt[1][iCent][iMCTruthLevel]->Clone ("htemp");
      h->Scale (std::pow (10, 2-iCent));
      myDrawHist (h, kBlack);
      SaferDelete (&h);

      h = h_jet_pt[1][iCent][0];
      //g = (TGAE*) g_jet_pt_syst[iCent][0]->Clone ();
      //SetCentralValuesKeepRelativeErrors (g, h);
      //ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
      //myDrawSystFill (g, colorfulSystColors[iCent+1], 1, 1001);
      //SaferDelete (&g);

      g = make_graph (h);
      //ResetXErrors (g);
      ScaleGraph (g, nullptr, std::pow (10, 2-iCent));
      myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.4, 1, 3, "P", false);
      SaferDelete (&g);

    } // end loop over iCent

    myText (0.18, 0.89, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.038);
    myText (0.64, 0.89, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.034);
    myText (0.64, 0.84, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.034);
    myText (0.64, 0.79, kBlack, Form ("#bf{#it{p}_{T}^{leading} > %s GeV}", strcmp (tag, "30GeVJets") == 0 ? "30" : "60"), 0.034);

    mySimpleMarkerAndBoxAndLineText (0.27, 0.16, 1.4, 1001, colorfulSystColors[0], 1.0, colorfulColors[0], kFullCircle, 1.6, "#it{pp} (#times10^{3})", 0.022/fPad);

    for (int iCent = 0; iCent < nZdcCentBins; iCent++) {
      mySimpleMarkerAndBoxAndLineText (0.27 + (iCent >= 2 ? 0.4 : 0), 0.16-((iCent+1)%3)*0.04, 1.4, 1001, colorfulSystColors[iCent+1], 1.0, colorfulColors[iCent+1], kFullCircle, 1.6, Form ("FCal %i-%i%% (#times10^{%i})", zdcCentPercs[iCent+1], zdcCentPercs[iCent], 2-iCent), 0.022/fPad);
    } // end loop over iCent
    mySimpleMarkerAndBoxAndLineText (0.67, 0.04, 1.4, 1001, colorfulSystColors[nZdcCentBins+1], 1.0, colorfulColors[nZdcCentBins+1], kFullCircle, 1.6, Form ("All cent. (#times10^{%i})", 2-nZdcCentBins), 0.022/fPad);
    mySimpleMarkerAndBoxAndLineText (0.27, 0.04, 1.4, 1001, kWhite, 0.0, kBlack, kDot, 0.0, "MC Truth (scaled)", 0.022/fPad);

    uPad->RedrawAxis();


    dPad->cd ();
    dPad->SetLogx();

    ymin=0.00;
    ymax=2.00;

    h = new TH1D ("htemp", ";#it{p}_{T}^{jet} [GeV];Reco. / Truth + Const.", 1, pTJBins[0], maxx);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitleSize (0.028/(1-fPad));
    h->GetXaxis ()->SetLabelSize (0.028/(1-fPad));
    h->GetXaxis ()->SetTitleOffset (3.2*(1-fPad));
    h->GetYaxis ()->SetTitleSize (0.028/(1-fPad));
    h->GetYaxis ()->SetLabelSize (0.028/(1-fPad));
    h->GetYaxis ()->SetTitleOffset (1.8*(1-fPad));
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (0);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    double x, y;

    const double offset = 0.25*3;

    h = h_jet_pt_ref[1][0];
    //g = (TGAE*) g_jet_pt_ref_syst[0]->Clone ();
    //SetCentralValuesKeepRelativeErrors (g, h);
    //ScaleGraph (g, h_jet_pt_ref[1][iMCTruthLevel]);
    //for (int i = 0; i < g->GetN (); i++) {
    //  g->GetPoint (i, x, y);
    //  g->SetPoint (i, x, y+offset);
    //}
    //myDrawSystFill (g, colorfulSystColors[0], 1, 1001);
    //SaferDelete (&g);

    g = make_graph (h);
    ScaleGraph (g, h_jet_pt_ref[1][iMCTruthLevel]);
    //ResetXErrors (g);
    for (int i = 0; i< g->GetN (); i++) {
      g->GetPoint (i, x, y);
      g->SetPoint (i, x, y+offset);
    }
    myDraw (g, colorfulColors[0], kFullCircle, 1.4, 1, 3, "P", false);
    SaferDelete (&g);

    l->SetLineWidth (3);
    l->SetLineColor (colorfulColors[0]);
    l->DrawLine (pTJBins[0], 1+offset, maxx, 1+offset);

    tl->SetTextAlign (offset > 0 ? 31 : 33);
    tl->SetTextFont (43);
    tl->SetTextSize (20);
    tl->SetTextColor (colorfulColors[0]);
    tl->DrawLatex (pTJBins[0] * std::exp (0.97 * std::log (maxx/pTJBins[0])), 1+offset+(offset > 0 ? 0.015:-0.015), Form ("#bf{#it{%s%g}}", offset == 0 ? "" : (offset > 0 ? "+ " : "#minus "), std::fabs (offset)));

    for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {

      const double offset = 0.25*(2-iCent);

      h = h_jet_pt[1][iCent][0];
      //g = (TGAE*) g_jet_pt_syst[iCent][0]->Clone ();
      //SetCentralValuesKeepRelativeErrors (g, h);
      //ScaleGraph (g, h_jet_pt[1][iCent][iMCTruthLevel]);
      //for (int i = 0; i < g->GetN (); i++) {
      //  g->GetPoint (i, x, y);
      //  g->SetPoint (i, x, y+offset);
      //}
      //myDrawSystFill (g, colorfulSystColors[iCent+1], 1, 1001);
      //SaferDelete (&g);

      g = make_graph (h);
      ScaleGraph (g, h_jet_pt[1][iCent][iMCTruthLevel]);
      //ResetXErrors (g);
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        g->SetPoint (i, x, y+offset);
      }
      myDraw (g, colorfulColors[iCent+1], kFullCircle, 1.4, 1, 3, "P", false);
      SaferDelete (&g);

      l->SetLineWidth (3);
      l->SetLineColor (colorfulColors[iCent+1]);
      l->DrawLine (pTJBins[0], 1+offset, maxx, 1+offset);

      tl->SetTextAlign (offset > 0 ? 31 : 33);
      tl->SetTextFont (43);
      tl->SetTextSize (20);
      tl->SetTextColor (colorfulColors[iCent+1]);
      tl->DrawLatex (pTJBins[0] * std::exp (0.97 * std::log (maxx/pTJBins[0])), 1+offset+(offset > 0 ? 0.015:-0.015), Form ("#bf{#it{%s%g}}", offset == 0 ? "" : (offset > 0 ? "+ " : "#minus "), std::fabs (offset)));

    } // end loop over iCent

    dPad->RedrawAxis();

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetPtSpectrum_MCClosure_%s.pdf", workPath.Data (), tag));

  }



  {
    const char* canvasName = "c_jet_eta_phi_ref";

    TCanvas* c = new TCanvas (canvasName, "", 880, 800);

    c->SetRightMargin (0.15);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi_ref[0][0];

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.22, 0.85, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetEtaPhiSpectrum_pp_%s.pdf", workPath.Data (), tag));
  }



  for (int iCent = 0; iCent < nZdcCentBins+1; iCent++) {
    const char* canvasName = Form ("c_jet_eta_phi_iCent%i", iCent);

    TCanvas* c = new TCanvas (canvasName, "", 880, 800);
    c->SetRightMargin (0.15);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi[0][iCent][0];

    h2->Draw ("colz");

    myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    if (iCent < nZdcCentBins)
      myText (0.22, 0.85, kBlack, Form ("#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{ZDC %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);
    else
      myText (0.22, 0.85, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, All cent.", 0.032);
    myText (0.22, 0.81, kBlack, TString (tag).Contains ("30GeVJets") ? "#it{p}_{T}^{jet} > 30 GeV" : "#it{p}_{T}^{jet} > 60 GeV", 0.032);

    c->SaveAs (Form ("%s/Plots/JetDistributions/JetEtaPhiSpectrum_pPb_iCent%i_%s.pdf", workPath.Data (), iCent, tag));
  }


}


#endif
