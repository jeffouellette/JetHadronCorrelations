#ifndef __PlotZDCAnalysis_C__
#define __PlotZDCAnalysis_C__

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>

#include <Utilities.h>

#include <iostream>
#include <fstream>

using namespace std;

TGraphErrors* TProfY2TGE (TProfile* py) {
  TGraphErrors* g = new TGraphErrors ();
  for (int iX = 1; iX <= py->GetNbinsX (); iX++) {
    g->SetPoint (g->GetN (), py->GetBinContent (iX), py->GetBinCenter (iX));
    g->SetPointError (g->GetN ()-1, 0, py->GetBinWidth (iX) / 2.);
  }
  return g;
}

TGraphErrors* TProfX2TGE (TProfile* px) {
  TGraphErrors* g = new TGraphErrors ();
  for (int iX = 1; iX <= px->GetNbinsX (); iX++) {
    g->SetPoint (g->GetN (), px->GetBinCenter (iX), px->GetBinContent (iX));
    g->SetPointError (g->GetN ()-1, px->GetBinWidth (iX) / 2., 0);
  }
  return g;
}

void PlotZDCAnalysis () {

  TFile* inFile = new TFile ("rootFiles/ZDCAnalysis/Nominal/outFile.root", "read");

  int runs[6] = {312796, 312837, 312937, 312945, 312968, 314199};

  TH1D* h_jet_Pb_fcal_et[6] = {};
  TH1D* h_mb_Pb_fcal_et[6] = {};
  TH1D* h_jet_Pb_fcal_et_sum = nullptr;
  TH1D* h_mb_Pb_fcal_et_sum = nullptr;
  TH1D* h_ratio_Pb_fcal_et_sum = nullptr;

  TH1D* h_jet_Pb_zdc_calibE[6] = {};
  TH1D* h_mb_Pb_zdc_calibE[6] = {};
  TH1D* h_jet_Pb_zdc_calibE_sum = nullptr;
  TH1D* h_mb_Pb_zdc_calibE_sum = nullptr;
  TH1D* h_ratio_Pb_zdc_calibE_sum = nullptr; // ratio of jet-tagged / MB events

  TH2D* h2_jet_Pb_fcal_et_zdc_calibE[6] = {};
  TH2D* h2_mb_Pb_fcal_et_zdc_calibE[6] = {};
  TH2D* h2_jet_Pb_fcal_et_zdc_calibE_sum = nullptr;
  TH2D* h2_mb_Pb_fcal_et_zdc_calibE_sum = nullptr;

  for (int iRun = 0; iRun < 6; iRun++) {
    h_jet_Pb_fcal_et[iRun] = (TH1D*) inFile->Get (Form ("h_jet_Pb_fcal_et_run%i", runs[iRun]));
    h_mb_Pb_fcal_et[iRun] = (TH1D*) inFile->Get (Form ("h_mb_Pb_fcal_et_run%i", runs[iRun]));
    h_jet_Pb_zdc_calibE[iRun] = (TH1D*) inFile->Get (Form ("h_jet_Pb_zdc_calibE_run%i", runs[iRun]));
    h_mb_Pb_zdc_calibE[iRun] = (TH1D*) inFile->Get (Form ("h_mb_Pb_zdc_calibE_run%i", runs[iRun]));
    h2_jet_Pb_fcal_et_zdc_calibE[iRun] = (TH2D*) inFile->Get (Form ("h2_jet_Pb_fcal_et_zdc_calibE_run%i", runs[iRun]));
    h2_mb_Pb_fcal_et_zdc_calibE[iRun] = (TH2D*) inFile->Get (Form ("h2_mb_Pb_fcal_et_zdc_calibE_run%i", runs[iRun]));
  }

  h_jet_Pb_fcal_et_sum = (TH1D*) h_jet_Pb_fcal_et[0]->Clone ("h_jet_Pb_fcal_et_sum");
  h_mb_Pb_fcal_et_sum = (TH1D*) h_mb_Pb_fcal_et[0]->Clone ("h_mb_Pb_fcal_et_sum");
  h_jet_Pb_fcal_et_sum->Reset ();
  h_mb_Pb_fcal_et_sum->Reset ();
  h_jet_Pb_zdc_calibE_sum = (TH1D*) h_jet_Pb_zdc_calibE[0]->Clone ("h_jet_Pb_zdc_calibE_sum");
  h_mb_Pb_zdc_calibE_sum = (TH1D*) h_mb_Pb_zdc_calibE[0]->Clone ("h_mb_Pb_zdc_calibE_sum");
  h_jet_Pb_zdc_calibE_sum->Reset ();
  h_mb_Pb_zdc_calibE_sum->Reset ();
  h2_jet_Pb_fcal_et_zdc_calibE_sum = (TH2D*) h2_jet_Pb_fcal_et_zdc_calibE[0]->Clone ("h2_jet_Pb_fcal_et_zdc_calibE_sum");
  h2_mb_Pb_fcal_et_zdc_calibE_sum = (TH2D*) h2_mb_Pb_fcal_et_zdc_calibE[0]->Clone ("h2_mb_Pb_fcal_et_zdc_calibE_sum");
  h2_jet_Pb_fcal_et_zdc_calibE_sum->Reset ();
  h2_mb_Pb_fcal_et_zdc_calibE_sum->Reset ();
  for (int iRun = 0; iRun < 6; iRun++) {
    h_jet_Pb_fcal_et_sum->Add (h_jet_Pb_fcal_et[iRun]);
    h_mb_Pb_fcal_et_sum->Add (h_mb_Pb_fcal_et[iRun]);
    h_jet_Pb_zdc_calibE_sum->Add (h_jet_Pb_zdc_calibE[iRun]);
    h_mb_Pb_zdc_calibE_sum->Add (h_mb_Pb_zdc_calibE[iRun]);
    h2_jet_Pb_fcal_et_zdc_calibE_sum->Add (h2_jet_Pb_fcal_et_zdc_calibE[iRun]);
    h2_mb_Pb_fcal_et_zdc_calibE_sum->Add (h2_mb_Pb_fcal_et_zdc_calibE[iRun]);
  }

  h_jet_Pb_fcal_et_sum->Scale (1. / h_jet_Pb_fcal_et_sum->Integral (h_jet_Pb_fcal_et_sum->FindBin (49.03), h_jet_Pb_fcal_et_sum->GetNbinsX ()));
  h_mb_Pb_fcal_et_sum->Scale (1. / h_mb_Pb_fcal_et_sum->Integral (h_mb_Pb_fcal_et_sum->FindBin (49.03), h_mb_Pb_fcal_et_sum->GetNbinsX ()));
  h_jet_Pb_zdc_calibE_sum->Scale (1. / h_jet_Pb_zdc_calibE_sum->Integral (h_jet_Pb_zdc_calibE_sum->FindBin (65.9228), h_jet_Pb_zdc_calibE_sum->GetNbinsX ()));
  h_mb_Pb_zdc_calibE_sum->Scale (1. / h_mb_Pb_zdc_calibE_sum->Integral (h_mb_Pb_zdc_calibE_sum->FindBin (65.9228), h_mb_Pb_zdc_calibE_sum->GetNbinsX ()));

  h_ratio_Pb_fcal_et_sum = (TH1D*) h_jet_Pb_fcal_et_sum->Clone ("h_ratio_Pb_fcal_et_sum");
  h_ratio_Pb_fcal_et_sum->Divide (h_mb_Pb_fcal_et_sum); 
  h_ratio_Pb_zdc_calibE_sum = (TH1D*) h_jet_Pb_zdc_calibE_sum->Clone ("h_ratio_Pb_zdc_calibE_sum");
  h_ratio_Pb_zdc_calibE_sum->Divide (h_mb_Pb_zdc_calibE_sum); 

  //TH2D* h2_jet_Pb_fcal_et_zdc_calibE_run312796 = (TH2D*) inFile->Get ("h2_jet_Pb_fcal_et_zdc_calibE_run312796");
  //TH2D* h2_mb_Pb_fcal_et_zdc_calibE_run312796 = (TH2D*) inFile->Get ("h2_mb_Pb_fcal_et_zdc_calibE_run312796");
  //TH2D* h2_jet_Pb_fcal_et_zdc_calibE_run312968 = (TH2D*) inFile->Get ("h2_jet_Pb_fcal_et_zdc_calibE_run312796");
  //TH2D* h2_mb_Pb_fcal_et_zdc_calibE_run312968 = (TH2D*) inFile->Get ("h2_mb_Pb_fcal_et_zdc_calibE_run312796");


  for (int iRun = 0; iRun < 6; iRun++) {
    TCanvas* c = new TCanvas ("c1", "", 800, 800);

    gPad->SetLogy ();

    gPad->SetBottomMargin (0.11);
    gPad->SetLeftMargin (0.11);
    gPad->SetRightMargin (0.04);
    gPad->SetTopMargin (0.04);

    TH1D* h = h_jet_Pb_zdc_calibE[iRun];
    h->GetXaxis ()->SetTitle ("Calibrated #Sigma#it{E}_{ZDC}^{Pb} [TeV]");
    h->GetYaxis ()->SetTitle ("Normalized counts");

    double ymin = 5e-7;
    double ymax = 8e-2;

    h->SetLineColor (kBlue+3);

    h->Scale (1./h->Integral ());

    h->GetYaxis ()->SetRangeUser (ymin, ymax);


    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (24);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (24);

    h->DrawCopy ("hist");

    h = h_mb_Pb_zdc_calibE[iRun];

    h->Scale (1./h->Integral ());

    h->SetLineColor (kRed+1);

    double xq[17];
    double yq[17];
    TString percs[17];
    xq[0] = 0.00; percs[0] = "100%";
    xq[1] = 0.10; percs[1] = "90%";
    xq[2] = 0.20; percs[2] = "80%";
    xq[3] = 0.30; percs[3] = "70%";
    xq[4] = 0.40; percs[4] = "60%";
    xq[5] = 0.50; percs[5] = "50%";
    xq[6] = 0.60; percs[6] = "40%";
    xq[7] = 0.70; percs[7] = "30%";
    xq[8] = 0.80; percs[8] = "20%";
    xq[9] = 0.90; percs[9] = "10%";
    xq[10] = 0.95; percs[10] = "5%";
    xq[11] = 0.98; percs[11] = "2%";
    xq[12] = 0.99; percs[12] = "1%";
    xq[13] = 0.995; percs[13] = "0.5%";
    xq[14] = 0.998; percs[14] = "0.2%";
    xq[15] = 0.999; percs[15] = "0.1%";
    xq[16] = 1.00; percs[16] = "0%";
    h->GetQuantiles (17, yq, xq);
    for (int i = 0; i < 17; i++)
      cout << xq[i] << ", " << yq[i] << endl;

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 16; i++) {
      divs->DrawLine (yq[i], ymin, yq[i], h->GetBinContent (h->FindBin (yq[i])));
      tl->DrawLatex (yq[i]+0.15, exp (0.1*log(ymax/ymin)) * ymin, percs[i].Data ());
    }

    h->DrawCopy ("hist same");

    myText (0.20, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.20, 0.840, kBlack, Form ("Run %i (per. A)", runs[iRun]), 0.032);
    myText (0.65, 0.550, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.510, kRed+1, "MinBias Trigger", 0.032);

    TPad* subPad = new TPad ("c3_subPad", "", 0.6, 0.6, 1.-c->GetRightMargin (), 1.-c->GetTopMargin ());
    subPad->SetTopMargin (0);
    subPad->SetRightMargin (0);
    subPad->SetLeftMargin (0.10);
    subPad->SetBottomMargin (0.08);
//    subPad->SetBorderSize (1);
    subPad->SetLogy ();
    subPad->Draw ();
    subPad->cd ();

    h = h_jet_Pb_zdc_calibE[iRun];

    h->GetXaxis ()->SetTitleOffset (999);
    h->GetYaxis ()->SetTitleOffset (999);

    h->GetXaxis ()->SetRangeUser (0, 12);
    h->GetYaxis ()->SetRangeUser (4e-5, 2e-1);

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (18);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (18);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (14);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (14);

    h->DrawCopy ("hist");

    h = h_mb_Pb_zdc_calibE[iRun];
    h->DrawCopy ("hist same");

    c->cd ();
    TLine* subPadBorder = new TLine ();
    subPadBorder->DrawLineNDC (0.6, 0.6, 1.-c->GetRightMargin (), 0.6);
    subPadBorder->DrawLineNDC (1.-c->GetRightMargin (), 0.6, 1.-c->GetRightMargin (), 1.-c->GetTopMargin ());
    subPadBorder->DrawLineNDC (1.-c->GetRightMargin (), 1.-c->GetTopMargin (), 0.6, 1.-c->GetTopMargin ());
    subPadBorder->DrawLineNDC (0.6, 1.-c->GetTopMargin (), 0.6, 0.6);
   
    c->SaveAs (Form ("%s/Plots/ZDCAnalysis/run%i_zdc.pdf", workPath.Data (), runs[iRun]));
  }



  {
    TCanvas* c = new TCanvas ("c", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_dPad", "", 0.0, 0.0, 1.0, 0.3);
    TPad* subPad = new TPad ("c_subPad", "", 0.6, 0.6*0.7+0.3, 1.-rMargin, 1.-tMargin*0.7);

    uPad->SetBottomMargin (0);
    uPad->SetLeftMargin (lMargin);
    uPad->SetRightMargin (rMargin);
    uPad->SetTopMargin (tMargin);

    dPad->SetBottomMargin (bMargin);
    dPad->SetLeftMargin (lMargin);
    dPad->SetRightMargin (rMargin);
    dPad->SetTopMargin (0);

    subPad->SetBottomMargin (0.08);
    subPad->SetLeftMargin (0.10);
    subPad->SetRightMargin (0);
    subPad->SetTopMargin (0);

    uPad->Draw ();
    dPad->Draw ();
    subPad->Draw ();


    uPad->cd ();
    uPad->SetLogy ();

    TH1D* h = h_jet_Pb_zdc_calibE_sum;
    h->GetYaxis ()->SetTitle ("A.U. (normalized in 0-20%)");

    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    double ymin = 5e-7;
    double ymax = 8e-2;

    h->SetLineColor (kBlue+3);

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

    h = h_mb_Pb_zdc_calibE_sum;

    h->SetLineColor (kRed+1);

    double xq[17];
    double yq[17];
    TString percs[17];
    xq[0] = 0.00; percs[0] = "100%";
    xq[1] = 0.10; percs[1] = "90%";
    xq[2] = 0.20; percs[2] = "80%";
    xq[3] = 0.30; percs[3] = "70%";
    xq[4] = 0.40; percs[4] = "60%";
    xq[5] = 0.50; percs[5] = "50%";
    xq[6] = 0.60; percs[6] = "40%";
    xq[7] = 0.70; percs[7] = "30%";
    xq[8] = 0.80; percs[8] = "20%";
    xq[9] = 0.90; percs[9] = "10%";
    xq[10] = 0.95; percs[10] = "5%";
    xq[11] = 0.98; percs[11] = "2%";
    xq[12] = 0.99; percs[12] = "1%";
    xq[13] = 0.995; percs[13] = "0.5%";
    xq[14] = 0.998; percs[14] = "0.2%";
    xq[15] = 0.999; percs[15] = "0.1%";
    xq[16] = 1.00; percs[16] = "0%";
    h->GetQuantiles (17, yq, xq);
    for (int i = 0; i < 17; i++)
      cout << xq[i] << ", " << yq[i] << endl;

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 16; i++) {
      divs->DrawLine (yq[i], ymin, yq[i], h->GetBinContent (h->FindBin (yq[i])));
      tl->DrawLatex (yq[i]+0.20, exp (0.1*log(ymax/ymin)) * ymin, percs[i].Data ());
    }

    h->DrawCopy ("hist same");

    myText (0.20, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.20, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.20, 0.820, kBlack, "All runs", 0.032);
    myText (0.65, 0.550, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.510, kRed+1, "MinBias Trigger", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = h_ratio_Pb_zdc_calibE_sum;
    h->GetXaxis ()->SetTitle ("Calibrated #Sigma#it{E}_{ZDC}^{Pb} [TeV]");
    h->GetYaxis ()->SetTitle ("J50 / MB");

    h->GetXaxis ()->SetTitleOffset (2.4 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());
    h->GetYaxis ()->CenterTitle ();

    h->SetLineColor (kBlue+3);

    h->GetYaxis ()->SetRangeUser (ymin, ymax);

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (24);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (24);

    h->DrawCopy ("hist");

    divs->DrawLine (0, 1, 175, 1);


    subPad->cd ();
    subPad->SetLogy ();

    ymin = 4e-5;
    ymax = 2e-1;

    h = h_jet_Pb_zdc_calibE_sum;

    h->GetXaxis ()->SetTitleOffset (999);
    h->GetYaxis ()->SetTitleOffset (999);

    h->GetXaxis ()->SetRangeUser (0, 12);
    h->GetYaxis ()->SetRangeUser (ymin, ymax);

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (18);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (18);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (14);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (14);

    h->DrawCopy ("hist");

    h = h_mb_Pb_zdc_calibE_sum;
    h->DrawCopy ("hist same");

    c->cd ();
    TLine* subPadBorder = new TLine ();
    subPadBorder->DrawLineNDC (0.6, 0.6*0.7+0.3, 1.-rMargin, 0.6*0.7+0.3);
    subPadBorder->DrawLineNDC (1.-rMargin, 0.6*0.7+0.3, 1.-rMargin, 1.-tMargin*0.7);
    subPadBorder->DrawLineNDC (1.-rMargin, 1.-tMargin*0.7, 0.6, 1.-tMargin*0.7);
    subPadBorder->DrawLineNDC (0.6, 1.-tMargin*0.7, 0.6, 0.6*0.7+0.3);


    c->SaveAs (Form ("%s/Plots/ZDCAnalysis/allruns_zdc.pdf", workPath.Data ()));
  }




  {
    TCanvas* c = new TCanvas ("c2", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c2_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c2_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    TH1D* h = h_jet_Pb_fcal_et_sum;
    h->GetYaxis ()->SetTitle ("A.U. (normalized in 0-20%)");

    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    double ymin = 5e-7;
    double ymax = 1e0;

    h->SetLineColor (kBlue+3);

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

    h = h_mb_Pb_fcal_et_sum;

    h->SetLineColor (kRed+1);

    double xq[12];
    double yq[12];
    TString percs[12];
    xq[0] = 0.00;   yq[0] = 0;        percs[0] = "100%";
    xq[1] = 0.10;   yq[1] = 3.54;     percs[1] = "90%";
    xq[2] = 0.20;   yq[2] = 7.60;     percs[2] = "80%";
    xq[3] = 0.30;   yq[3] = 12.27;    percs[3] = "70%";
    xq[4] = 0.40;   yq[4] = 17.51;    percs[4] = "60%";
    xq[5] = 0.50;   yq[5] = 23.45;    percs[5] = "50%";
    xq[6] = 0.60;   yq[6] = 30.29;    percs[6] = "40%";
    xq[7] = 0.70;   yq[7] = 38.49;    percs[7] = "30%";
    xq[8] = 0.80;   yq[8] = 49.03;    percs[8] = "20%";
    xq[9] = 0.90;   yq[9] = 65.08;    percs[9] = "10%";
    xq[10] = 0.95;  yq[10] = 79.45;   percs[10] = "5%";
    xq[11] = 0.99;  yq[11] = 109.01;  percs[11] = "1%";
    for (int i = 0; i < 12; i++)
      cout << xq[i] << ", " << yq[i] << endl;

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 12; i++) {
      divs->DrawLine (yq[i], ymin, yq[i], h->GetBinContent (h->FindBin (yq[i])));
      tl->DrawLatex (yq[i]+0.20, exp (0.1*log(ymax/ymin)) * ymin, percs[i].Data ());
    }

    h->DrawCopy ("hist same");

    myText (0.65, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.65, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.820, kBlack, "All runs", 0.032);
    myText (0.65, 0.740, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.700, kRed+1, "MinBias Trigger", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = h_ratio_Pb_fcal_et_sum;
    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{Pb-FCal} [GeV]");
    h->GetYaxis ()->SetTitle ("J50 / MB");

    h->GetXaxis ()->SetTitleOffset (2.4 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());
    h->GetYaxis ()->CenterTitle ();

    h->SetLineColor (kBlue+3);

    h->GetYaxis ()->SetRangeUser (ymin, ymax);

    h->GetXaxis ()->SetTitleFont (43);
    h->GetXaxis ()->SetTitleSize (26);
    h->GetYaxis ()->SetTitleFont (43);
    h->GetYaxis ()->SetTitleSize (26);
    h->GetXaxis ()->SetLabelFont (43);
    h->GetXaxis ()->SetLabelSize (24);
    h->GetYaxis ()->SetLabelFont (43);
    h->GetYaxis ()->SetLabelSize (24);

    h->DrawCopy ("hist");

    divs->DrawLine (-30, 1, 220, 1);

    c->SaveAs (Form ("%s/Plots/ZDCAnalysis/allruns_fcal_et.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c3", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = h2_jet_Pb_fcal_et_zdc_calibE_sum;
    h->RebinY (4);
    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    h->GetYaxis ()->SetTitle ("Calibrated #Sigma#it{E}_{ZDC}^{Pb} [TeV]");
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

    h->Draw ("colz");


    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

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

    myText (0.15, 0.96, kWhite, "All runs, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/ZDCAnalysis/J50_zdc_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c4", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = h2_mb_Pb_fcal_et_zdc_calibE_sum;
    h->RebinY (4);

    h->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    h->GetYaxis ()->SetTitle ("Calibrated #Sigma#it{E}_{ZDC}^{Pb} [TeV]");
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

    h->Draw ("colz");


    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

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

    myText (0.15, 0.96, kWhite, "All runs, MinBias Trigger", 0.032);
   
    c->SaveAs (Form ("%s/Plots/ZDCAnalysis/MB_zdc_fcalet_correlation.pdf", workPath.Data ()));
  }
}

#endif
