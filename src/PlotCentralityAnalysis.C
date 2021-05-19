#ifndef __PlotCentralityAnalysis_C__
#define __PlotCentralityAnalysis_C__

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


void PlotCentralityAnalysis () { 

  TFile* inFile = new TFile (Form ("%s/CentralityAnalysis/Nominal/allRuns.root", rootPath.Data ()), "read");

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

  TH2D* h2_jet_Pb_fcal_et_zdc_calibE = nullptr;
  TH2D* h2_mb_Pb_fcal_et_zdc_calibE = nullptr;

  TH2D* h2_jet_Pb_fcal_et_Pb_q2 = nullptr;
  TH2D* h2_mb_Pb_fcal_et_Pb_q2 = nullptr;

  TH1D* h_jet_Pb_fcal_et_zdcSelected = nullptr;
  TH1D* h_mb_Pb_fcal_et_zdcSelected = nullptr;
  TH1D* h_ratio_Pb_fcal_et_zdcSelected = nullptr;

  TH3D* h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE");
  TH3D* h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE");

  for (int iRun = 0; iRun < 6; iRun++) {
    h_jet_Pb_fcal_et[iRun] = (TH1D*) inFile->Get (Form ("h_jet_Pb_fcal_et_run%i", runs[iRun]));
    h_mb_Pb_fcal_et[iRun] = (TH1D*) inFile->Get (Form ("h_mb_Pb_fcal_et_run%i", runs[iRun]));
    h_jet_Pb_zdc_calibE[iRun] = (TH1D*) inFile->Get (Form ("h_jet_Pb_zdc_calibE_run%i", runs[iRun]));
    h_mb_Pb_zdc_calibE[iRun] = (TH1D*) inFile->Get (Form ("h_mb_Pb_zdc_calibE_run%i", runs[iRun]));
  }

  h_jet_Pb_fcal_et_sum = (TH1D*) h_jet_Pb_fcal_et[0]->Clone ("h_jet_Pb_fcal_et_sum");
  h_mb_Pb_fcal_et_sum = (TH1D*) h_mb_Pb_fcal_et[0]->Clone ("h_mb_Pb_fcal_et_sum");
  h_jet_Pb_fcal_et_sum->Reset ();
  h_mb_Pb_fcal_et_sum->Reset ();
  h_jet_Pb_zdc_calibE_sum = (TH1D*) h_jet_Pb_zdc_calibE[0]->Clone ("h_jet_Pb_zdc_calibE_sum");
  h_mb_Pb_zdc_calibE_sum = (TH1D*) h_mb_Pb_zdc_calibE[0]->Clone ("h_mb_Pb_zdc_calibE_sum");
  h_jet_Pb_zdc_calibE_sum->Reset ();
  h_mb_Pb_zdc_calibE_sum->Reset ();
  for (int iRun = 0; iRun < 6; iRun++) {
    if (h_jet_Pb_fcal_et[iRun]) h_jet_Pb_fcal_et_sum->Add (h_jet_Pb_fcal_et[iRun]);
    if (h_mb_Pb_fcal_et[iRun]) h_mb_Pb_fcal_et_sum->Add (h_mb_Pb_fcal_et[iRun]);
    if (h_jet_Pb_zdc_calibE[iRun]) h_jet_Pb_zdc_calibE_sum->Add (h_jet_Pb_zdc_calibE[iRun]);
    if (h_mb_Pb_zdc_calibE[iRun]) h_mb_Pb_zdc_calibE_sum->Add (h_mb_Pb_zdc_calibE[iRun]);
  }

  h2_jet_Pb_fcal_et_zdc_calibE = (TH2D*) inFile->Get ("h2_jet_Pb_fcal_et_zdc_calibE");//h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");
  h2_mb_Pb_fcal_et_zdc_calibE = (TH2D*) inFile->Get ("h2_mb_Pb_fcal_et_zdc_calibE");//h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");

  h_jet_Pb_fcal_et_sum->Scale (1. / h_jet_Pb_fcal_et_sum->Integral (h_jet_Pb_fcal_et_sum->FindBin (fcalCentBins[1]), h_jet_Pb_fcal_et_sum->GetNbinsX ()));
  h_mb_Pb_fcal_et_sum->Scale (1. / h_mb_Pb_fcal_et_sum->Integral (h_mb_Pb_fcal_et_sum->FindBin (fcalCentBins[1]), h_mb_Pb_fcal_et_sum->GetNbinsX ()));
  h_jet_Pb_zdc_calibE_sum->Scale (1. / h_jet_Pb_zdc_calibE_sum->Integral (h_jet_Pb_zdc_calibE_sum->FindBin (zdcCentBins[1]), h_jet_Pb_zdc_calibE_sum->GetNbinsX ()));
  h_mb_Pb_zdc_calibE_sum->Scale (1. / h_mb_Pb_zdc_calibE_sum->Integral (h_mb_Pb_zdc_calibE_sum->FindBin (zdcCentBins[1]), h_mb_Pb_zdc_calibE_sum->GetNbinsX ()));

  h_ratio_Pb_fcal_et_sum = (TH1D*) h_jet_Pb_fcal_et_sum->Clone ("h_ratio_Pb_fcal_et_sum");
  h_ratio_Pb_fcal_et_sum->Divide (h_mb_Pb_fcal_et_sum); 
  h_ratio_Pb_zdc_calibE_sum = (TH1D*) h_jet_Pb_zdc_calibE_sum->Clone ("h_ratio_Pb_zdc_calibE_sum");
  h_ratio_Pb_zdc_calibE_sum->Divide (h_mb_Pb_zdc_calibE_sum); 

  h_jet_Pb_fcal_et_zdcSelected = h2_jet_Pb_fcal_et_zdc_calibE->ProjectionX ("h_jet_Pb_fcal_et_zdcSelected", h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[1]), h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->GetNbins ());
  h_mb_Pb_fcal_et_zdcSelected = h2_mb_Pb_fcal_et_zdc_calibE->ProjectionX ("h_mb_Pb_fcal_et_zdcSelected", h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[1]), h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->GetNbins ());

  h_jet_Pb_fcal_et_zdcSelected->Scale (1. / h_jet_Pb_fcal_et_zdcSelected->Integral (h_jet_Pb_fcal_et_zdcSelected->FindBin (fcalCentBins[1]), h_jet_Pb_fcal_et_zdcSelected->GetNbinsX ()));
  h_mb_Pb_fcal_et_zdcSelected->Scale (1. / h_mb_Pb_fcal_et_zdcSelected->Integral (h_mb_Pb_fcal_et_zdcSelected->FindBin (fcalCentBins[1]), h_mb_Pb_fcal_et_zdcSelected->GetNbinsX ()));

  h_ratio_Pb_fcal_et_zdcSelected = (TH1D*) h_jet_Pb_fcal_et_zdcSelected->Clone ("h_ratio_Pb_fcal_et_zdcSelected");
  h_ratio_Pb_fcal_et_zdcSelected->Divide (h_mb_Pb_fcal_et_zdcSelected);


  double xq[101];
  double yq[101];
  TString percs[101];
  for (int iPerc = 0; iPerc < 101; iPerc++) {
    xq[iPerc] = 0.00;
    xq[iPerc] = ((double)iPerc) / 100.;
    percs[iPerc] = Form ("%i%%", 100-iPerc);
  }


  //for (int iRun = 0; iRun < 6; iRun++) {
  //  TCanvas* c = new TCanvas ("c1", "", 800, 800);

  //  gPad->SetLogy ();

  //  gPad->SetBottomMargin (0.11);
  //  gPad->SetLeftMargin (0.11);
  //  gPad->SetRightMargin (0.04);
  //  gPad->SetTopMargin (0.04);

  //  TH1D* h = h_jet_Pb_zdc_calibE[iRun];
  //  h->GetXaxis ()->SetTitle ("Calibrated #Sigma#it{E}_{ZDC}^{Pb} [TeV]");
  //  h->GetYaxis ()->SetTitle ("Normalized counts");

  //  double ymin = 5e-7;
  //  double ymax = 8e-2;

  //  h->SetLineColor (kBlue+3);

  //  h->Scale (1./h->Integral ());

  //  h->GetYaxis ()->SetRangeUser (ymin, ymax);


  //  h->GetXaxis ()->SetTitleFont (43);
  //  h->GetXaxis ()->SetTitleSize (26);
  //  h->GetYaxis ()->SetTitleFont (43);
  //  h->GetYaxis ()->SetTitleSize (26);
  //  h->GetXaxis ()->SetLabelFont (43);
  //  h->GetXaxis ()->SetLabelSize (24);
  //  h->GetYaxis ()->SetLabelFont (43);
  //  h->GetYaxis ()->SetLabelSize (24);

  //  h->DrawCopy ("hist");

  //  h = h_mb_Pb_zdc_calibE[iRun];

  //  h->Scale (1./h->Integral ());

  //  h->SetLineColor (kRed+1);

  //  double xq[17];
  //  double yq[17];
  //  TString percs[17];
  //  xq[0] = 0.00; percs[0] = "100%";
  //  xq[1] = 0.10; percs[1] = "90%";
  //  xq[2] = 0.20; percs[2] = "80%";
  //  xq[3] = 0.30; percs[3] = "70%";
  //  xq[4] = 0.40; percs[4] = "60%";
  //  xq[5] = 0.50; percs[5] = "50%";
  //  xq[6] = 0.60; percs[6] = "40%";
  //  xq[7] = 0.70; percs[7] = "30%";
  //  xq[8] = 0.80; percs[8] = "20%";
  //  xq[9] = 0.90; percs[9] = "10%";
  //  xq[10] = 0.95; percs[10] = "5%";
  //  xq[11] = 0.98; percs[11] = "2%";
  //  xq[12] = 0.99; percs[12] = "1%";
  //  xq[13] = 0.995; percs[13] = "0.5%";
  //  xq[14] = 0.998; percs[14] = "0.2%";
  //  xq[15] = 0.999; percs[15] = "0.1%";
  //  xq[16] = 1.00; percs[16] = "0%";
  //  h->GetQuantiles (17, yq, xq);
  //  for (int i = 0; i < 17; i++)
  //    cout << xq[i] << ", " << yq[i] << endl;

  //  TLine* divs = new TLine ();
  //  TLatex* tl = new TLatex ();
  //  tl->SetTextAngle (-90);
  //  tl->SetTextAlign (11);
  //  tl->SetTextFont (43);
  //  tl->SetTextSize (14);
  //  divs->SetLineStyle (2);
  //  for (int i = 1; i < 16; i++) {
  //    divs->DrawLine (yq[i], ymin, yq[i], h->GetBinContent (h->FindBin (yq[i])));
  //    tl->DrawLatex (yq[i]+0.15, exp (0.1*log(ymax/ymin)) * ymin, percs[i].Data ());
  //  }

  //  h->DrawCopy ("hist same");

  //  myText (0.20, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
  //  myText (0.20, 0.840, kBlack, Form ("Run %i (per. A)", runs[iRun]), 0.032);
  //  myText (0.65, 0.550, kBlue+3, "J50 Trigger", 0.032);
  //  myText (0.65, 0.510, kRed+1, "MinBias Trigger", 0.032);

  //  TPad* subPad = new TPad ("c3_subPad", "", 0.6, 0.6, 1.-c->GetRightMargin (), 1.-c->GetTopMargin ());
  //  subPad->SetTopMargin (0);
  //  subPad->SetRightMargin (0);
  //  subPad->SetLeftMargin (0.10);
  //  subPad->SetBottomMargin (0.08);
////    subPad->SetBorderSize (1);
  //  subPad->SetLogy ();
  //  subPad->Draw ();
  //  subPad->cd ();

  //  h = h_jet_Pb_zdc_calibE[iRun];

  //  h->GetXaxis ()->SetTitleOffset (999);
  //  h->GetYaxis ()->SetTitleOffset (999);

  //  h->GetXaxis ()->SetRangeUser (0, 12);
  //  h->GetYaxis ()->SetRangeUser (4e-5, 2e-1);

  //  h->GetXaxis ()->SetTitleFont (43);
  //  h->GetXaxis ()->SetTitleSize (18);
  //  h->GetYaxis ()->SetTitleFont (43);
  //  h->GetYaxis ()->SetTitleSize (18);
  //  h->GetXaxis ()->SetLabelFont (43);
  //  h->GetXaxis ()->SetLabelSize (14);
  //  h->GetYaxis ()->SetLabelFont (43);
  //  h->GetYaxis ()->SetLabelSize (14);

  //  h->DrawCopy ("hist");

  //  h = h_mb_Pb_zdc_calibE[iRun];
  //  h->DrawCopy ("hist same");

  //  c->cd ();
  //  TLine* subPadBorder = new TLine ();
  //  subPadBorder->DrawLineNDC (0.6, 0.6, 1.-c->GetRightMargin (), 0.6);
  //  subPadBorder->DrawLineNDC (1.-c->GetRightMargin (), 0.6, 1.-c->GetRightMargin (), 1.-c->GetTopMargin ());
  //  subPadBorder->DrawLineNDC (1.-c->GetRightMargin (), 1.-c->GetTopMargin (), 0.6, 1.-c->GetTopMargin ());
  //  subPadBorder->DrawLineNDC (0.6, 1.-c->GetTopMargin (), 0.6, 0.6);
  // 
  //  c->SaveAs (Form ("%s/Plots/CentralityAnalysis/run%i_zdc.pdf", workPath.Data (), runs[iRun]));
  //}



  {
    TCanvas* c = new TCanvas ("c_allruns_zdc", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_allruns_zdc_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_allruns_zdc_dPad", "", 0.0, 0.0, 1.0, 0.3);
    TPad* subPad = new TPad ("c_allruns_zdc_subPad", "", 0.6, 0.6*0.7+0.3, 1.-rMargin, 1.-tMargin*0.7);

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

    TH1D* h = (TH1D*) h_jet_Pb_zdc_calibE_sum->Clone ("htemp");
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
    SaferDelete (&h);

    h = (TH1D*) h_mb_Pb_zdc_calibE_sum->Clone ("htemp");

    h->SetLineColor (kRed+1);

    double plot_xq[101];
    double plot_yq[17];
    TString plot_percs[17];
    plot_xq[0]  = 0.000; plot_percs[0]  = "100%";
    plot_xq[1]  = 0.100; plot_percs[1]  = "90%";
    plot_xq[2]  = 0.200; plot_percs[2]  = "80%";
    plot_xq[3]  = 0.300; plot_percs[3]  = "70%";
    plot_xq[4]  = 0.400; plot_percs[4]  = "60%";
    plot_xq[5]  = 0.500; plot_percs[5]  = "50%";
    plot_xq[6]  = 0.600; plot_percs[6]  = "40%";
    plot_xq[7]  = 0.700; plot_percs[7]  = "30%";
    plot_xq[8]  = 0.800; plot_percs[8]  = "20%";
    plot_xq[9]  = 0.900; plot_percs[9]  = "10%";
    plot_xq[10] = 0.950; plot_percs[10] = "5%";
    plot_xq[11] = 0.980; plot_percs[11] = "2%";
    plot_xq[12] = 0.990; plot_percs[12] = "1%";
    plot_xq[13] = 0.995; plot_percs[13] = "0.5%";
    plot_xq[14] = 0.998; plot_percs[14] = "0.2%";
    plot_xq[15] = 0.999; plot_percs[15] = "0.1%";
    plot_xq[16] = 1.000; plot_percs[16] = "0%";
    h->GetQuantiles (17, plot_yq, plot_xq);
    h->GetQuantiles (101, yq, xq);

    ofstream cutsfile;
    cutsfile.open (Form ("%s/aux/ZdcCentCuts.dat", workPath.Data ()));
    cutsfile << "CALIBRATED ENERGIES (over all runs)\n";
    for (int i = 0; i < 101; i++) {
      std::cout << xq[i] << ", " << yq[i] << std::endl;
      cutsfile << std::setw (6) << percs[i] << "\t" << yq[i] << "\n";
    }
    cutsfile.close ();

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 16; i++) {
      divs->DrawLine (plot_yq[i], ymin, plot_yq[i], h->GetBinContent (h->FindBin (plot_yq[i])));
      tl->DrawLatex (plot_yq[i]+0.20, exp (0.1*log(ymax/ymin)) * ymin, plot_percs[i].Data ());
    }

    h->DrawCopy ("hist same");
    SaferDelete (&h);

    myText (0.20, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.20, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.20, 0.820, kBlack, "All runs", 0.032);
    myText (0.65, 0.550, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.510, kRed+1, "MinBias Trigger", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = (TH1D*) h_ratio_Pb_zdc_calibE_sum->Clone ("htemp");
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
    SaferDelete (&h);

    divs->DrawLine (0, 1, 175, 1);


    subPad->cd ();
    subPad->SetLogy ();

    ymin = 4e-5;
    ymax = 2e-1;

    h = (TH1D*) h_jet_Pb_zdc_calibE_sum->Clone ("htemp");

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
    SaferDelete (&h);

    h = (TH1D*) h_mb_Pb_zdc_calibE_sum->Clone ("htemp");

    h->DrawCopy ("hist same");
    SaferDelete (&h);

    c->cd ();
    TLine* subPadBorder = new TLine ();
    subPadBorder->DrawLineNDC (0.6, 0.6*0.7+0.3, 1.-rMargin, 0.6*0.7+0.3);
    subPadBorder->DrawLineNDC (1.-rMargin, 0.6*0.7+0.3, 1.-rMargin, 1.-tMargin*0.7);
    subPadBorder->DrawLineNDC (1.-rMargin, 1.-tMargin*0.7, 0.6, 1.-tMargin*0.7);
    subPadBorder->DrawLineNDC (0.6, 1.-tMargin*0.7, 0.6, 0.6*0.7+0.3);


    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allruns_zdc.pdf", workPath.Data ()));
  }




  {
    TCanvas* c = new TCanvas ("c_allruns_fcal_et", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_allruns_fcal_et_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_allruns_fcal_et_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    TH1D* h = (TH1D*) h_jet_Pb_fcal_et_sum->Clone ("htemp");
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
    SaferDelete (&h);

    h = (TH1D*) h_mb_Pb_fcal_et_sum->Clone ("htemp");

    h->SetLineColor (kRed+1);

    double plot_xq[17];
    double plot_yq[17];
    TString plot_percs[17];
    plot_xq[0]  = 0.000; plot_percs[0]  = "100%";
    plot_xq[1]  = 0.100; plot_percs[1]  = "90%";
    plot_xq[2]  = 0.200; plot_percs[2]  = "80%";
    plot_xq[3]  = 0.300; plot_percs[3]  = "70%";
    plot_xq[4]  = 0.400; plot_percs[4]  = "60%";
    plot_xq[5]  = 0.500; plot_percs[5]  = "50%";
    plot_xq[6]  = 0.600; plot_percs[6]  = "40%";
    plot_xq[7]  = 0.700; plot_percs[7]  = "30%";
    plot_xq[8]  = 0.800; plot_percs[8]  = "20%";
    plot_xq[9]  = 0.900; plot_percs[9]  = "10%";
    plot_xq[10] = 0.950; plot_percs[10] = "5%";
    plot_xq[11] = 0.980; plot_percs[11] = "2%";
    plot_xq[12] = 0.990; plot_percs[12] = "1%";
    plot_xq[13] = 0.995; plot_percs[13] = "0.5%";
    plot_xq[14] = 0.998; plot_percs[14] = "0.2%";
    plot_xq[15] = 0.999; plot_percs[15] = "0.1%";
    plot_xq[16] = 1.000; plot_percs[16] = "0%";
    h->GetQuantiles (17, plot_yq, plot_xq);
    h->GetQuantiles (101, yq, xq);

    ofstream cutsfile;
    cutsfile.open (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()));
    for (int i = 0; i < 101; i++) {
      std::cout << xq[i] << ", " << yq[i] << std::endl;
      cutsfile << std::setw (6) << percs[i] << "\t" << yq[i] << "\n";
    }
    cutsfile.close ();

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 16; i++) {
      divs->DrawLine (plot_yq[i], ymin, plot_yq[i], h->GetBinContent (h->FindBin (plot_yq[i])));
      tl->DrawLatex (plot_yq[i]+0.20, exp (0.1*log(ymax/ymin)) * ymin, plot_percs[i].Data ());
    }

    h->DrawCopy ("hist same");
    SaferDelete (&h);

    myText (0.65, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.65, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.820, kBlack, "All runs", 0.032);
    myText (0.65, 0.740, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.700, kRed+1, "MinBias Trigger", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = (TH1D*) h_ratio_Pb_fcal_et_sum->Clone ("htemp");
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
    SaferDelete (&h);

    divs->DrawLine (-30, 1, 220, 1);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allruns_fcal_et.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_allruns_fcal_et_zdcSelected", "", 800, 800);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_allruns_fcal_et_zdcSelected_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_allruns_fcal_et_zdcSelected_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    TH1D* h = (TH1D*) h_jet_Pb_fcal_et_sum->Clone ("htemp");
    h->GetYaxis ()->SetTitle ("A.U. (normalized in 0-20%)");

    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    double ymin = 5e-9;
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
    SaferDelete (&h);

    h = (TH1D*) h_mb_Pb_fcal_et_sum->Clone ("htemp");

    h->SetLineColor (kRed+1);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    h = (TH1D*) h_jet_Pb_fcal_et_zdcSelected->Clone ("htemp");

    h->SetLineColor (kBlue+3);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    h = (TH1D*) h_mb_Pb_fcal_et_zdcSelected->Clone ("htemp");

    h->SetLineColor (kRed+1);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    myText (0.65, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.65, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.820, kBlack, "All runs", 0.032);
    myText (0.65, 0.740, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.700, kRed+1, "MinBias Trigger", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = (TH1D*) h_ratio_Pb_fcal_et_sum->Clone ("htemp");
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
    SaferDelete (&h);

    h = (TH1D*) h_ratio_Pb_fcal_et_zdcSelected->Clone ("htemp");

    h->SetLineColor (kBlue+3);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);

    TLine* tline = new TLine ();
    tline->SetLineStyle (2);
    tline->DrawLine (-30, 1, 220, 1);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allruns_fcal_et_zdcSelected.pdf", workPath.Data ()));
  }




}

#endif
