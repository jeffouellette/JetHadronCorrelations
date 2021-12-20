#ifndef __PlotCentralityAnalysis2D_C__
#define __PlotCentralityAnalysis2D_C__

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>

#include <Utilities.h>
#include <MyColors.h>

#include <AggressiveAvocado.h>

#include "CentralityDefs.h"
#include "Params.h"
#include "LocalUtilities.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace JetHadronCorrelations;


void PlotCentralityAnalysis2D () { 

  SetAggressiveAvocadoStyle ();

  TFile* inFile = nullptr;


  int iZdc20Percent = 0;
  while (iZdc20Percent < nZdcCentBins && zdcCentPercs[iZdc20Percent] > 20) iZdc20Percent++;
  int iFcal20Percent = 0;
  while (iFcal20Percent < nFcalCentBins && fcalCentPercs[iFcal20Percent] > 100) iFcal20Percent++;


  TH2D* h2_jet_Pb_fcal_et_zdc_calibE = nullptr;
  TH2D* h2_mb_Pb_fcal_et_zdc_calibE = nullptr;

  TH2D* h2_ljet_pt_vs_Pb_fcal_et = nullptr; // correlation between leading jet pT and Pb side FCal Et without data overlay


  inFile = new TFile (Form ("%s/CentralityAnalysis/Nominal/allpPbRuns.root", rootPath.Data ()), "read");

  h2_jet_Pb_fcal_et_zdc_calibE = (TH2D*) inFile->Get ("h2_jet_Pb_fcal_et_zdc_calibE");//h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");
  h2_mb_Pb_fcal_et_zdc_calibE = (TH2D*) inFile->Get ("h2_mb_Pb_fcal_et_zdc_calibE");//h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");


  inFile = new TFile (Form ("%s/CentralityAnalysis/Nominal/mc16_so_all.root", rootPath.Data ()), "read");
  h2_ljet_pt_vs_Pb_fcal_et = (TH2D*) inFile->Get ("h2_ljet_pt_vs_Pb_fcal_et");


  {
    TCanvas* c = new TCanvas ("c_J50_zdc_fcalet_correlation", "", 1000, 812);

    gPad->SetLogz ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.16);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_zdc_calibE->Clone ("htemp");
    h->RebinX (4);
    h->RebinY (4);
    TProfile* h_px = h->ProfileX ("h_px");
    TProfile* h_py = h->ProfileX ("h_py");
    TGAE* g_px = TProfX2TGAE (h_px);
    TGAE* g_py = TProfY2TGAE (h_py);

    double x,y;
    for (int i = 0; i < g_py->GetN (); i++) {
      g_py->GetPoint (i, x, y);
      if (y > 127) g_py->RemovePoint (i--);
    }
    for (int i = 0; i < g_px->GetN (); i++) {
      g_px->GetPoint (i, x, y);
      if (x < -10 || x > 161) g_px->RemovePoint (i--);
    }

    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    h->GetYaxis ()->SetTitle ("Calibrated #Sigma#it{E}_{ZDC}^{Pb} [TeV]");
    h->GetZaxis ()->SetTitle ("Events");

    h->GetXaxis ()->SetTitleOffset (1.1 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.1 * h->GetYaxis ()->GetTitleOffset ());

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (30);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (30);
    h->GetZaxis ()->SetTitleFont (43);
    h->GetZaxis ()->SetTitleSize (30);

    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (26);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (26);
    h->GetZaxis ()->SetLabelFont (43);
    h->GetZaxis ()->SetLabelSize (26);

    h->DrawCopy ("colz");
    SaferDelete (&h);

    g_px->SetLineColor (kWhite);
    g_py->SetLineColor (kWhite);
    g_px->SetLineWidth (1);
    g_py->SetLineWidth (1);

    g_px->SetMarkerColor (kWhite);
    g_py->SetMarkerColor (kWhite);
    g_px->SetMarkerStyle (kDot);
    g_py->SetMarkerStyle (kDot);

    g_px->Draw ("P");
    g_py->Draw ("P");

    myText (0.14, 0.96, kWhite, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.45, 0.96, kWhite, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_zdc_fcalet_correlation.pdf", workPath.Data ()));

    SaferDelete (&h);
    SaferDelete (&h_px);
    SaferDelete (&h_py);
    SaferDelete (&g_px);
    SaferDelete (&g_py);
  }



  {
    TCanvas* c = new TCanvas ("MB_zdc_fcalet_correlation", "", 1000, 812);

    //gPad->SetLogz ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.16);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_zdc_calibE->Clone ("htemp");
    TProfile* h_px = h->ProfileX ("h_px");
    TProfile* h_py = h->ProfileY ("h_py");
    TGAE* g_px = TProfX2TGAE (h_px);
    TGAE* g_py = TProfY2TGAE (h_py);

    double x,y;
    for (int i = 0; i < g_py->GetN (); i++) {
      g_py->GetPoint (i, x, y);
      if (y > 127) g_py->RemovePoint (i--);
    }
    for (int i = 0; i < g_px->GetN (); i++) {
      g_px->GetPoint (i, x, y);
      if (x < -10 || x > 161) g_px->RemovePoint (i--);
    }

    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    h->GetYaxis ()->SetTitle ("Calibrated #Sigma#it{E}_{ZDC}^{Pb} [TeV]");
    h->GetZaxis ()->SetTitle ("Events");

    h->GetXaxis ()->SetTitleOffset (1.1 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.1 * h->GetYaxis ()->GetTitleOffset ());

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (30);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (30);
    h->GetZaxis ()->SetTitleFont (43);
    h->GetZaxis ()->SetTitleSize (30);

    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (26);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (26);
    h->GetZaxis ()->SetLabelFont (43);
    h->GetZaxis ()->SetLabelSize (26);

    h->DrawCopy ("colz");
    SaferDelete (&h);

    g_px->SetLineColor (kWhite);
    g_py->SetLineColor (kWhite);
    g_px->SetLineWidth (1);
    g_py->SetLineWidth (1);

    g_px->SetMarkerColor (kWhite);
    g_py->SetMarkerColor (kWhite);
    g_px->SetMarkerStyle (kDot);
    g_py->SetMarkerStyle (kDot);

    g_px->Draw ("P");
    g_py->Draw ("P");

    myText (0.14, 0.96, kWhite, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.47, 0.96, kWhite, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, MinBias Trigger", 0.032);
   
    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_zdc_fcalet_correlation.pdf", workPath.Data ()));

    SaferDelete (&h);
    SaferDelete (&h_px);
    SaferDelete (&h_py);
    SaferDelete (&g_px);
    SaferDelete (&g_py);
  }



  {
    TCanvas* c = new TCanvas ("c_ljet_pt_vs_Pb_fcal_et", "", 1000, 812);

    gPad->SetLogx ();
    gPad->SetLogz ();

    const double bMargin = 0.13;
    const double lMargin = 0.13;
    const double rMargin = 0.16;
    const double tMargin = 0.06;

    gPad->SetBottomMargin (bMargin);
    gPad->SetLeftMargin (lMargin);
    gPad->SetRightMargin (rMargin);
    gPad->SetTopMargin (tMargin);

    TH2D* h = (TH2D*) h2_ljet_pt_vs_Pb_fcal_et->Clone ("htemp");
    h->RebinX (2);

    // normalize distribution in response (along y-axis)
    for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
      double integral = 0;
      for (int iY = 1; iY <= h->GetNbinsY (); iY++)
        integral += h->GetBinContent (iX, iY);
      if (integral > 0) {
        for (int iY = 1; iY <= h->GetNbinsY (); iY++) {
          h->SetBinContent (iX, iY, h->GetBinContent (iX, iY) / integral);
          h->SetBinError (iX, iY, h->GetBinError (iX, iY) / integral);
        } // end loop over iY
      }
    } // end loop over iX

    const float zmin = 1e-4;
    const float zmax = h->GetMaximum ();

    TH1D* h_px = h->ProjectionX ("h_px");
    for (int iX = 1; iX <= h_px->GetNbinsX (); iX++) {
      TH1D* h_temp = h->ProjectionY ("h_temp", iX, iX);
      h_px->SetBinContent (iX, h_temp->GetBinCenter (h_temp->GetMaximumBin ()));
      h_px->SetBinError (iX, 0.5 * h_temp->GetBinWidth (h_temp->GetMaximumBin ()));
      SaferDelete (&h_temp);
    }
    //TProfile* h_px = h->ProfileX ("h_px");
    //TGAE* g_px = TProfX2TGAE (h->ProfileX ("h_px"));

    TF1* f_px = new TF1 ("f_ljet_pt_vs_mean_Pb_fcal_et", "[0]*log(x)+[1]", h->GetXaxis ()->GetBinLowEdge (1), h->GetXaxis ()->GetBinLowEdge (h->GetNbinsX ()+1));
    h_px->Fit (f_px, "RN0Q");

    //TGAE* g_px = TProfX2TGAE (h_px);
    TGAE* g_px = make_graph (h_px);


    h->GetXaxis ()->SetTitle ("Leading jet #it{p}_{T}^{truth} [GeV]");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, A} [GeV]");
    h->GetZaxis ()->SetTitle ("A.U.");

    h->GetXaxis ()->SetTitleOffset (1.1 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.1 * h->GetYaxis ()->GetTitleOffset ());

    h->GetZaxis ()->SetRangeUser (zmin, zmax);

    h->GetXaxis ()->SetMoreLogLabels ();

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (30);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (30);
    h->GetZaxis ()->SetTitleFont (43);
    h->GetZaxis ()->SetTitleSize (30);

    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (26);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (26);
    h->GetZaxis ()->SetLabelFont (43);
    h->GetZaxis ()->SetLabelSize (26);

    h->DrawCopy ("colz");

    myDraw (g_px, kWhite, kFullCircle, 0.8);

    myDraw (f_px, kOrange+2, 2, 2);

    myText (0.14, 0.96, kWhite, "#bf{#it{ATLAS}} Simulation Internal", 0.032);
    myText (0.45, 0.96, kWhite, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV, #it{y}_{CM} = -0.465", 0.032);

    myMarkerText (0.22, 0.86, kWhite, kFullCircle, "#color[0]{Mode of #Sigma#it{E}_{T}^{FCal,A}}", 1.2, 0.028, true);
    myText (0.18, 0.82, kWhite, "Linear fit", 0.028);
    myText (0.18, 0.78, kWhite, Form ("Slope: %.3f #pm %.3f", f_px->GetParameter (0), f_px->GetParError (0)), 0.028);
    myText (0.18, 0.74, kWhite, Form ("y-int: %.3f #pm %.3f", f_px->GetParameter (1), f_px->GetParError (1)), 0.028);


    TPad* subPad = new TPad ("c_allpPbRuns_zdc_subPad", "", 0.4, 0.6, 1.-rMargin, 1.-tMargin);

    subPad->SetBottomMargin (0.11);
    subPad->SetLeftMargin (0.09);
    subPad->SetRightMargin (0.00);
    subPad->SetTopMargin (0.00);

    subPad->SetFillColor (c->GetFillColor ());

    subPad->Draw ();

    subPad->cd ();
    subPad->SetLogx ();
    subPad->SetLogz ();

    TH2D* h_sub = (TH2D*) h->Clone ("h_sub");

    h_sub->GetYaxis ()->SetRangeUser (2, 12);
    h_sub->GetZaxis ()->SetRangeUser (zmin, zmax);

    h_sub->GetXaxis ()->SetTitleOffset (999);
    h_sub->GetYaxis ()->SetTitleOffset (999);

    h_sub->GetXaxis ()->SetLabelFont (43);
    h_sub->GetXaxis ()->SetLabelSize (16);
    h_sub->GetYaxis ()->SetLabelFont (43);
    h_sub->GetYaxis ()->SetLabelSize (16);

    h_sub->DrawCopy ("col");

    myDraw (g_px, kWhite, kFullCircle, 0.8);

    //((TF1*) f_px->Clone ())->Draw ("same");
    myDraw (f_px, kOrange+2, 2, 2);

    c->cd ();
    TLine* subPadBorder = new TLine ();
    subPadBorder->SetLineColor (kWhite);
    subPadBorder->DrawLineNDC (0.4, 0.6, 1.-rMargin, 0.6);
    subPadBorder->DrawLineNDC (1.-rMargin, 0.6, 1.-rMargin, 1.-tMargin);
    subPadBorder->DrawLineNDC (1.-rMargin, 1.-tMargin, 0.4, 1.-tMargin);
    subPadBorder->DrawLineNDC (0.4, 1.-tMargin, 0.4, 0.6);

    SaferDelete (&h);
    SaferDelete (&h_sub);
    SaferDelete (&h_px);
    SaferDelete (&g_px);
    SaferDelete (&f_px);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/Ljet_pt_vs_Pb_fcal_et.pdf", workPath.Data ()));
  }
}

#endif
