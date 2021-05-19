#ifndef __PlotJetSubtractedEnergy_C__
#define __PlotJetSubtractedEnergy_C__

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>

#include <Utilities.h>

#include "Params.h"
#include "LocalUtilities.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace JetHadronCorrelations;

TColor* tcolor = new TColor ();
const Color_t myBlue = (Color_t) tcolor->GetColor (45, 64, 245);
const Color_t myPurple = (Color_t) tcolor->GetColor (130,  10, 130);
const Color_t myRed = (Color_t) tcolor->GetColor (255,  12,  73);
const Color_t myGreen = (Color_t) tcolor->GetColor ( 54, 167,  80);
const Color_t myOrange = (Color_t) tcolor->GetColor (255,  68,   0);

const Color_t colors[] = {myBlue, myGreen, myOrange, myRed, myPurple};


void PlotJetSubtractedEnergy () { 

  TFile* inFile = new TFile (Form ("%s/JetSubtractedEnergy/Nominal/allRuns.root", rootPath.Data ()), "read");


  TH2D* h2_zdc_calibE_jet_subEt = (TH2D*) inFile->Get ("h2_zdc_calibE_jet_subEt");
  TH2D* h2_zdc_calibE_jet_subE  = (TH2D*) inFile->Get ("h2_zdc_calibE_jet_subE");
  TH2D* h2_Pb_fcal_et_jet_subEt = (TH2D*) inFile->Get ("h2_Pb_fcal_et_jet_subEt");
  TH2D* h2_Pb_fcal_et_jet_subE  = (TH2D*) inFile->Get ("h2_Pb_fcal_et_jet_subE");

  const int numCentBins = (DoFcalCentVar () ? numFcalCentBins : (DoFineFcalCentVar () ? numFineFcalCentBins : numZdcCentBins));
  TH1D** h_jet_subEt = new TH1D* [numCentBins];
  TH1D** h_jet_subEt_ratio = new TH1D* [numCentBins];
  TH1D** h_jet_subE = new TH1D* [numCentBins];
  TH1D** h_jet_subE_ratio = new TH1D* [numCentBins];

  TH1D*  h_jet_subEt_ref = (TH1D*) inFile->Get ("h_jet_subEt");
  h_jet_subEt_ref->Rebin (2);
  h_jet_subEt_ref->Scale (1./ h_jet_subEt_ref->Integral (), "width");
  TH1D*  h_jet_subE_ref = (TH1D*) inFile->Get ("h_jet_subE");
  h_jet_subE_ref->Rebin (2);
  h_jet_subE_ref->Scale (1./ h_jet_subE_ref->Integral (), "width");

  for (int iCent = 0; iCent < numCentBins; iCent++) {

    h_jet_subEt[iCent] = (TH1D*) inFile->Get (Form ("h_jet_subEt_iCent%i", iCent));
    h_jet_subEt[iCent]->Rebin (2);
    h_jet_subEt[iCent]->Scale (1./ h_jet_subEt[iCent]->Integral (), "width");
    h_jet_subEt_ratio[iCent] = (TH1D*) h_jet_subEt[iCent]->Clone (Form ("h_jet_subEt_ratio_iCent%i", iCent));
    h_jet_subEt_ratio[iCent]->Divide (h_jet_subEt_ref);

    h_jet_subE[iCent] = (TH1D*) inFile->Get (Form ("h_jet_subE_iCent%i", iCent));
    h_jet_subE[iCent]->Rebin (2);
    h_jet_subE[iCent]->Scale (1./ h_jet_subE[iCent]->Integral (), "width");
    h_jet_subE_ratio[iCent] = (TH1D*) h_jet_subE[iCent]->Clone (Form ("h_jet_subE_ratio_iCent%i", iCent));
    h_jet_subE_ratio[iCent]->Divide (h_jet_subE_ref);

  }




  {
    TCanvas* c = new TCanvas ("c_subEt_zdc", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_subEt_zdc_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_subEt_zdc_dPad", "", 0.0, 0.0, 1.0, 0.3);

    uPad->SetBottomMargin (0);
    uPad->SetLeftMargin (lMargin);
    uPad->SetRightMargin (rMargin);
    uPad->SetTopMargin (tMargin);

    dPad->SetBottomMargin (bMargin);
    dPad->SetLeftMargin (lMargin);
    dPad->SetRightMargin (rMargin);
    dPad->SetTopMargin (0);

    uPad->Draw ();
    dPad->Draw ();


    uPad->cd ();
    uPad->SetLogy ();

    TH1D* h = (TH1D*) h_jet_subEt_ref->Clone ("htemp");
    h->GetYaxis ()->SetTitle ("Normalized entries");

    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    double ymin = 5e-5;
    double ymax = 8e0;

    h->SetLineColor (kBlack);

    h->GetYaxis ()->SetRangeUser (ymin, ymax);

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (0);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (0);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (24);

    h->DrawCopy ("hist");
    SaferDelete (&h);

    for (int iCent = 0; iCent < numCentBins; iCent++) {

      h = (TH1D*) h_jet_subEt[iCent]->Clone ("htemp");

      h->SetLineColor (colors[iCent]);

      h->DrawCopy ("hist same");
      SaferDelete (&h);
    }

    myText (0.20, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.20, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.20, 0.820, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.20, 0.780, kBlack, "#it{p}_{T}^{jet} > 60 GeV", 0.032);
    myText (0.20, 0.740, kBlack, "#Sigma#it{E}_{ZDC}^{Pb} Percentiles", 0.032);

    myText (0.6, 0.860, kBlack, "#bf{#it{pp}}", 0.032);
    for (int iCent = 0; iCent < numCentBins; iCent++)
      myText (0.6, 0.82-iCent*0.04, colors[iCent], Form ("#bf{#it{p}+Pb, %i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 5e-2;
    ymax = 2e1;

    h = (TH1D*) h_jet_subEt_ratio[0]->Clone ("htemp");
    h->Reset ();
    for (int iX = 1; iX <= h->GetNbinsX (); iX++) h->SetBinContent (iX, 1);
    h->GetXaxis ()->SetTitle ("Underlying event \"size\" (#it{E}_{T}^{subtracted}) [GeV]");
    h->GetYaxis ()->SetTitle ("Ratio to #it{pp}");

    h->GetYaxis ()->SetRangeUser (ymin, ymax);

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (24);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (24);

    h->GetXaxis ()->SetTitleOffset (2.4 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());
    h->GetYaxis ()->CenterTitle ();

    h->SetLineStyle (2);
    h->SetLineWidth (1);
    h->DrawCopy ("hist");
    SaferDelete (&h);

    for (int iCent = 0; iCent < numCentBins; iCent++) {
      h = (TH1D*) h_jet_subEt_ratio[iCent]->Clone ("htemp");
      h->SetLineColor (colors[iCent]);

      h->DrawCopy ("hist same");
      SaferDelete (&h);
    }

    c->SaveAs (Form ("%s/Plots/JetSubtractedEnergy/subEt_zdc_comparisons.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("zdc_subEt_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.13);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_zdc_calibE_jet_subEt->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{ZDC}^{Pb} [TeV]");
    h->GetYaxis ()->SetTitle ("#it{E}_{T}^{subtracted} [GeV]");
    h->GetZaxis ()->SetTitle ("");

    h->GetXaxis ()->SetTitleOffset (1.2 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetZaxis ()->SetTitleFont (43);
    h->GetZaxis ()->SetTitleSize (26);

    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (22);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (22);
    h->GetZaxis ()->SetLabelFont (43);
    h->GetZaxis ()->SetLabelSize (22);

    h->GetXaxis ()->SetAxisColor (kWhite);
    h->GetXaxis ()->SetLabelColor (kWhite);
    h->GetXaxis ()->SetTitleColor (kWhite);
    h->GetYaxis ()->SetAxisColor (kWhite);
    h->GetYaxis ()->SetLabelColor (kWhite);
    h->GetYaxis ()->SetTitleColor (kWhite);
    h->GetZaxis ()->SetAxisColor (kWhite);
    h->GetZaxis ()->SetLabelColor (kWhite);
    h->GetZaxis ()->SetTitleColor (kWhite);

    h->DrawCopy ("colz");
    SaferDelete (&h);

    g_px->SetLineColor (kWhite);
    g_py->SetLineColor (kWhite);
    g_px->SetLineWidth (2);
    g_py->SetLineWidth (2);

    g_px->SetMarkerColor (kWhite);
    g_py->SetMarkerColor (kWhite);
    g_px->SetMarkerStyle (kDot);
    g_py->SetMarkerStyle (kDot);

    g_px->Draw ("P");
    g_py->Draw ("P");

    myText (0.14, 0.96, kWhite, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.44, 0.96, kWhite, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #it{p}_{T}^{jet} > 60 GeV", 0.032);
   
    c->SaveAs (Form ("%s/Plots/JetSubtractedEnergy/subEt_zdc_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("zdc_subE_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.13);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_zdc_calibE_jet_subE->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{ZDC}^{Pb} [TeV]");
    h->GetYaxis ()->SetTitle ("#it{E}^{subtracted} [GeV]");
    h->GetZaxis ()->SetTitle ("");

    h->GetXaxis ()->SetTitleOffset (1.2 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetZaxis ()->SetTitleFont (43);
    h->GetZaxis ()->SetTitleSize (26);

    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (22);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (22);
    h->GetZaxis ()->SetLabelFont (43);
    h->GetZaxis ()->SetLabelSize (22);

    h->GetXaxis ()->SetAxisColor (kWhite);
    h->GetXaxis ()->SetLabelColor (kWhite);
    h->GetXaxis ()->SetTitleColor (kWhite);
    h->GetYaxis ()->SetAxisColor (kWhite);
    h->GetYaxis ()->SetLabelColor (kWhite);
    h->GetYaxis ()->SetTitleColor (kWhite);
    h->GetZaxis ()->SetAxisColor (kWhite);
    h->GetZaxis ()->SetLabelColor (kWhite);
    h->GetZaxis ()->SetTitleColor (kWhite);

    h->DrawCopy ("colz");
    SaferDelete (&h);

    g_px->SetLineColor (kWhite);
    g_py->SetLineColor (kWhite);
    g_px->SetLineWidth (2);
    g_py->SetLineWidth (2);

    g_px->SetMarkerColor (kWhite);
    g_py->SetMarkerColor (kWhite);
    g_px->SetMarkerStyle (kDot);
    g_py->SetMarkerStyle (kDot);

    g_px->Draw ("P");
    g_py->Draw ("P");

    myText (0.14, 0.96, kWhite, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.44, 0.96, kWhite, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #it{p}_{T}^{jet} > 60 GeV", 0.032);
   
    c->SaveAs (Form ("%s/Plots/JetSubtractedEnergy/subE_zdc_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("fcalet_subEt_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.13);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_Pb_fcal_et_jet_subEt->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    h->GetYaxis ()->SetTitle ("#it{E}_{T}^{subtracted} [GeV]");
    h->GetZaxis ()->SetTitle ("");

    h->GetXaxis ()->SetTitleOffset (1.2 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetZaxis ()->SetTitleFont (43);
    h->GetZaxis ()->SetTitleSize (26);

    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (22);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (22);
    h->GetZaxis ()->SetLabelFont (43);
    h->GetZaxis ()->SetLabelSize (22);

    h->GetXaxis ()->SetAxisColor (kWhite);
    h->GetXaxis ()->SetLabelColor (kWhite);
    h->GetXaxis ()->SetTitleColor (kWhite);
    h->GetYaxis ()->SetAxisColor (kWhite);
    h->GetYaxis ()->SetLabelColor (kWhite);
    h->GetYaxis ()->SetTitleColor (kWhite);
    h->GetZaxis ()->SetAxisColor (kWhite);
    h->GetZaxis ()->SetLabelColor (kWhite);
    h->GetZaxis ()->SetTitleColor (kWhite);

    h->DrawCopy ("colz");
    SaferDelete (&h);

    g_px->SetLineColor (kWhite);
    g_py->SetLineColor (kWhite);
    g_px->SetLineWidth (2);
    g_py->SetLineWidth (2);

    g_px->SetMarkerColor (kWhite);
    g_py->SetMarkerColor (kWhite);
    g_px->SetMarkerStyle (kDot);
    g_py->SetMarkerStyle (kDot);

    g_px->Draw ("P");
    g_py->Draw ("P");

    myText (0.14, 0.96, kWhite, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.44, 0.96, kWhite, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #it{p}_{T}^{jet} > 60 GeV", 0.032);
   
    c->SaveAs (Form ("%s/Plots/JetSubtractedEnergy/subEt_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("fcalet_subE_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.13);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_Pb_fcal_et_jet_subE->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    h->GetYaxis ()->SetTitle ("#it{E}^{subtracted} [GeV]");
    h->GetZaxis ()->SetTitle ("");

    h->GetXaxis ()->SetTitleOffset (1.2 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetZaxis ()->SetTitleFont (43);
    h->GetZaxis ()->SetTitleSize (26);

    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (22);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (22);
    h->GetZaxis ()->SetLabelFont (43);
    h->GetZaxis ()->SetLabelSize (22);

    h->GetXaxis ()->SetAxisColor (kWhite);
    h->GetXaxis ()->SetLabelColor (kWhite);
    h->GetXaxis ()->SetTitleColor (kWhite);
    h->GetYaxis ()->SetAxisColor (kWhite);
    h->GetYaxis ()->SetLabelColor (kWhite);
    h->GetYaxis ()->SetTitleColor (kWhite);
    h->GetZaxis ()->SetAxisColor (kWhite);
    h->GetZaxis ()->SetLabelColor (kWhite);
    h->GetZaxis ()->SetTitleColor (kWhite);

    h->DrawCopy ("colz");
    SaferDelete (&h);

    g_px->SetLineColor (kWhite);
    g_py->SetLineColor (kWhite);
    g_px->SetLineWidth (2);
    g_py->SetLineWidth (2);

    g_px->SetMarkerColor (kWhite);
    g_py->SetMarkerColor (kWhite);
    g_px->SetMarkerStyle (kDot);
    g_py->SetMarkerStyle (kDot);

    g_px->Draw ("P");
    g_py->Draw ("P");

    myText (0.14, 0.96, kWhite, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.44, 0.96, kWhite, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #it{p}_{T}^{jet} > 60 GeV", 0.032);
   
    c->SaveAs (Form ("%s/Plots/JetSubtractedEnergy/subE_fcalet_correlation.pdf", workPath.Data ()));
  }




}

#endif
