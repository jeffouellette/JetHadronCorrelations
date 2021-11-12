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
#include <MyColors.h>

#include "CentralityDefs.h"
#include "Params.h"
#include "LocalUtilities.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace JetHadronCorrelations;


void PlotCentralityAnalysis () { 

  TFile* inFile = nullptr;


  int iZdc20Percent = 0;
  while (iZdc20Percent < nZdcCentBins && zdcCentPercs[iZdc20Percent] > 20) iZdc20Percent++;
  int iFcal20Percent = 0;
  while (iFcal20Percent < nFcalCentBins && fcalCentPercs[iFcal20Percent] > 100) iFcal20Percent++;


  TH1D* h_jet_Pb_fcal_et[npPbRuns] = {};
  TH1D* h_mb_Pb_fcal_et[npPbRuns] = {};
  TH1D* h_jet_Pb_fcal_et_sum = nullptr;
  TH1D* h_mb_Pb_fcal_et_sum = nullptr;
  TH1D* h_ratio_Pb_fcal_et_sum = nullptr;

  TH1D* h_jet_p_fcal_et[nppRuns] = {};
  TH1D* h_mb_p_fcal_et[nppRuns] = {};
  TH1D* h_jet_p_fcal_et_sum = nullptr;
  TH1D* h_mb_p_fcal_et_sum = nullptr;
  TH1D* h_ratio_p_fcal_et_sum = nullptr;

  TH1D* h_jet_Pb_zdc_calibE[npPbRuns] = {};
  TH1D* h_mb_Pb_zdc_calibE[npPbRuns] = {};
  TH1D* h_jet_Pb_zdc_calibE_sum = nullptr;
  TH1D* h_mb_Pb_zdc_calibE_sum = nullptr;
  TH1D* h_ratio_Pb_zdc_calibE_sum = nullptr; // ratio of jet-tagged / MB events

  TH2D* h2_jet_Pb_fcal_et_zdc_calibE = nullptr;
  TH2D* h2_mb_Pb_fcal_et_zdc_calibE = nullptr;

  TH1D* h_jet_Pb_fcal_et_zdc_0t20 = nullptr;
  TH1D* h_mb_Pb_fcal_et_zdc_0t20 = nullptr;
  TH1D* h_ratio_Pb_fcal_et_zdc_0t20 = nullptr;
  TH1D* h_jet_Pb_fcal_et_zdc_20t40 = nullptr;
  TH1D* h_mb_Pb_fcal_et_zdc_20t40 = nullptr;
  TH1D* h_ratio_Pb_fcal_et_zdc_20t40 = nullptr;
  TH1D* h_jet_Pb_fcal_et_zdc_40t60 = nullptr;
  TH1D* h_mb_Pb_fcal_et_zdc_40t60 = nullptr;
  TH1D* h_ratio_Pb_fcal_et_zdc_40t60 = nullptr;
  TH1D* h_jet_Pb_fcal_et_zdc_60t80 = nullptr;
  TH1D* h_mb_Pb_fcal_et_zdc_60t80 = nullptr;
  TH1D* h_ratio_Pb_fcal_et_zdc_60t80 = nullptr;
  TH1D* h_jet_Pb_fcal_et_zdc_80t100 = nullptr;
  TH1D* h_mb_Pb_fcal_et_zdc_80t100 = nullptr;
  TH1D* h_ratio_Pb_fcal_et_zdc_80t100 = nullptr;

  TH3D* h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE = nullptr;
  TH3D* h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE = nullptr;

  //TGAE* g_ALICE_zna_calibE = new TGAE ();
  //TGAE* g_ALICE_zna_calibE_subpanel = new TGAE (); 


  TH1D* h_mb_Pb_fcal_et_sum_pPb_mc = nullptr;
  TH1D* h_mb_p_fcal_et_sum_pPb_mc = nullptr;
  TH1D* h_mb_Pb_fcal_et_sum_pPb_mc_corr = nullptr;
  TH1D* h_mb_p_fcal_et_sum_pPb_mc_corr = nullptr;
  TH1D* h_mb_Pb_fcal_et_sum_pp_mc = nullptr;
  TH1D* h_mb_p_fcal_et_sum_pp_mc = nullptr;


  inFile = new TFile (Form ("%s/CentralityAnalysis/Nominal/allppRuns.root", rootPath.Data ()), "read");

  for (int iRun = 0; iRun < nppRuns; iRun++) {
    h_jet_p_fcal_et[iRun] = (TH1D*) inFile->Get (Form ("h_jet_p_fcal_et_run%i", ppRuns[iRun]));
    h_mb_p_fcal_et[iRun] = (TH1D*) inFile->Get (Form ("h_mb_p_fcal_et_run%i", ppRuns[iRun]));
  }
  h_jet_p_fcal_et_sum = (TH1D*) h_jet_p_fcal_et[0]->Clone ("h_jet_p_fcal_et_sum");
  h_mb_p_fcal_et_sum = (TH1D*) h_mb_p_fcal_et[0]->Clone ("h_mb_p_fcal_et_sum");
  h_jet_p_fcal_et_sum->Reset ();
  h_mb_p_fcal_et_sum->Reset ();
  for (int iRun = 0; iRun < nppRuns; iRun++) {
    if (h_jet_p_fcal_et[iRun]) h_jet_p_fcal_et_sum->Add (h_jet_p_fcal_et[iRun]);
    if (h_mb_p_fcal_et[iRun]) h_mb_p_fcal_et_sum->Add (h_mb_p_fcal_et[iRun]);
  }
  //h_jet_p_fcal_et_sum->Scale (1., "width");
  //h_mb_p_fcal_et_sum->Scale (1., "width");
  //h_jet_p_fcal_et_sum->Scale (1./ h_jet_p_fcal_et_sum->Integral (), "width");
  //h_mb_p_fcal_et_sum->Scale (1./ h_mb_p_fcal_et_sum->Integral (), "width");

  h_ratio_p_fcal_et_sum = (TH1D*) h_jet_p_fcal_et_sum->Clone ("h_ratio_p_fcal_et_sum");
  h_ratio_p_fcal_et_sum->Divide (h_mb_p_fcal_et_sum); 


  inFile = new TFile (Form ("%s/CentralityAnalysis/Nominal/allpPbRuns.root", rootPath.Data ()), "read");


  h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE");
  h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE");

  for (int iRun = 0; iRun < npPbRuns; iRun++) {
    h_jet_Pb_fcal_et[iRun] = (TH1D*) inFile->Get (Form ("h_jet_Pb_fcal_et_run%i", pPbRuns[iRun]));
    h_mb_Pb_fcal_et[iRun] = (TH1D*) inFile->Get (Form ("h_mb_Pb_fcal_et_run%i", pPbRuns[iRun]));
    h_jet_Pb_zdc_calibE[iRun] = (TH1D*) inFile->Get (Form ("h_jet_Pb_zdc_calibE_run%i", pPbRuns[iRun]));
    h_mb_Pb_zdc_calibE[iRun] = (TH1D*) inFile->Get (Form ("h_mb_Pb_zdc_calibE_run%i", pPbRuns[iRun]));
  }

  h_jet_Pb_fcal_et_sum = (TH1D*) h_jet_Pb_fcal_et[0]->Clone ("h_jet_Pb_fcal_et_sum");
  h_mb_Pb_fcal_et_sum = (TH1D*) h_mb_Pb_fcal_et[0]->Clone ("h_mb_Pb_fcal_et_sum");
  h_jet_Pb_fcal_et_sum->Reset ();
  h_mb_Pb_fcal_et_sum->Reset ();
  h_jet_Pb_zdc_calibE_sum = (TH1D*) h_jet_Pb_zdc_calibE[0]->Clone ("h_jet_Pb_zdc_calibE_sum");
  h_mb_Pb_zdc_calibE_sum = (TH1D*) h_mb_Pb_zdc_calibE[0]->Clone ("h_mb_Pb_zdc_calibE_sum");
  h_jet_Pb_zdc_calibE_sum->Reset ();
  h_mb_Pb_zdc_calibE_sum->Reset ();
  for (int iRun = 0; iRun < npPbRuns; iRun++) {
    if (h_jet_Pb_fcal_et[iRun]) h_jet_Pb_fcal_et_sum->Add (h_jet_Pb_fcal_et[iRun]);
    if (h_mb_Pb_fcal_et[iRun]) h_mb_Pb_fcal_et_sum->Add (h_mb_Pb_fcal_et[iRun]);
    if (h_jet_Pb_zdc_calibE[iRun]) h_jet_Pb_zdc_calibE_sum->Add (h_jet_Pb_zdc_calibE[iRun]);
    if (h_mb_Pb_zdc_calibE[iRun]) h_mb_Pb_zdc_calibE_sum->Add (h_mb_Pb_zdc_calibE[iRun]);
  }

  h2_jet_Pb_fcal_et_zdc_calibE = (TH2D*) inFile->Get ("h2_jet_Pb_fcal_et_zdc_calibE");//h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");
  h2_mb_Pb_fcal_et_zdc_calibE = (TH2D*) inFile->Get ("h2_mb_Pb_fcal_et_zdc_calibE");//h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");

  //h_jet_Pb_fcal_et_sum->Scale (1., "width");
  //h_mb_Pb_fcal_et_sum->Scale (1., "width");
  //h_jet_Pb_zdc_calibE_sum->Scale (1., "width");
  //h_mb_Pb_zdc_calibE_sum->Scale (1., "width");
  //h_jet_Pb_fcal_et_sum->Scale (1. / h_jet_Pb_fcal_et_sum->Integral (h_jet_Pb_fcal_et_sum->FindBin (fcalCentBins[iFcal20Percent]), h_jet_Pb_fcal_et_sum->GetNbinsX ()), "width");
  //h_mb_Pb_fcal_et_sum->Scale (1. / h_mb_Pb_fcal_et_sum->Integral (h_mb_Pb_fcal_et_sum->FindBin (fcalCentBins[iFcal20Percent]), h_mb_Pb_fcal_et_sum->GetNbinsX ()), "width");
  //h_jet_Pb_zdc_calibE_sum->Scale (1. / h_jet_Pb_zdc_calibE_sum->Integral (h_jet_Pb_zdc_calibE_sum->FindBin (zdcCentBins[iZdc20Percent]), h_jet_Pb_zdc_calibE_sum->GetNbinsX ()), "width");
  //h_mb_Pb_zdc_calibE_sum->Scale (1. / h_mb_Pb_zdc_calibE_sum->Integral (h_mb_Pb_zdc_calibE_sum->FindBin (zdcCentBins[iZdc20Percent]), h_mb_Pb_zdc_calibE_sum->GetNbinsX ()), "width");

  h_ratio_Pb_fcal_et_sum = (TH1D*) h_jet_Pb_fcal_et_sum->Clone ("h_ratio_Pb_fcal_et_sum");
  h_ratio_Pb_fcal_et_sum->Divide (h_mb_Pb_fcal_et_sum); 
  h_ratio_Pb_zdc_calibE_sum = (TH1D*) h_jet_Pb_zdc_calibE_sum->Clone ("h_ratio_Pb_zdc_calibE_sum");
  h_ratio_Pb_zdc_calibE_sum->Divide (h_mb_Pb_zdc_calibE_sum); 
  
  h_jet_Pb_fcal_et_zdc_0t20 = h2_jet_Pb_fcal_et_zdc_calibE->ProjectionX ("h_jet_Pb_fcal_et_zdc_0t20", h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[4]), h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[5]));
  h_mb_Pb_fcal_et_zdc_0t20 = h2_mb_Pb_fcal_et_zdc_calibE->ProjectionX ("h_mb_Pb_fcal_et_zdc_0t20", h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[4]), h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[5]));
  h_jet_Pb_fcal_et_zdc_20t40 = h2_jet_Pb_fcal_et_zdc_calibE->ProjectionX ("h_jet_Pb_fcal_et_zdc_20t40", h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[3]), h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[4]));
  h_mb_Pb_fcal_et_zdc_20t40 = h2_mb_Pb_fcal_et_zdc_calibE->ProjectionX ("h_mb_Pb_fcal_et_zdc_20t40", h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[3]), h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[4]));
  h_jet_Pb_fcal_et_zdc_40t60 = h2_jet_Pb_fcal_et_zdc_calibE->ProjectionX ("h_jet_Pb_fcal_et_zdc_40t60", h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[2]), h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[3]));
  h_mb_Pb_fcal_et_zdc_40t60 = h2_mb_Pb_fcal_et_zdc_calibE->ProjectionX ("h_mb_Pb_fcal_et_zdc_40t60", h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[2]), h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[3]));
  h_jet_Pb_fcal_et_zdc_60t80 = h2_jet_Pb_fcal_et_zdc_calibE->ProjectionX ("h_jet_Pb_fcal_et_zdc_60t80", h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[1]), h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[2]));
  h_mb_Pb_fcal_et_zdc_60t80 = h2_mb_Pb_fcal_et_zdc_calibE->ProjectionX ("h_mb_Pb_fcal_et_zdc_60t80", h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[1]), h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[2]));
  h_jet_Pb_fcal_et_zdc_80t100 = h2_jet_Pb_fcal_et_zdc_calibE->ProjectionX ("h_jet_Pb_fcal_et_zdc_80t100", h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[0]), h2_jet_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[1]));
  h_mb_Pb_fcal_et_zdc_80t100 = h2_mb_Pb_fcal_et_zdc_calibE->ProjectionX ("h_mb_Pb_fcal_et_zdc_80t100", h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[0]), h2_mb_Pb_fcal_et_zdc_calibE->GetYaxis ()->FindBin (zdcCentBins[1]));

  h_jet_Pb_fcal_et_zdc_0t20->Scale (0.2 / h_jet_Pb_fcal_et_zdc_0t20->Integral (), "width");
  h_mb_Pb_fcal_et_zdc_0t20->Scale (0.2 / h_mb_Pb_fcal_et_zdc_0t20->Integral (), "width");
  h_jet_Pb_fcal_et_zdc_20t40->Scale (0.2 / h_jet_Pb_fcal_et_zdc_20t40->Integral (), "width");
  h_mb_Pb_fcal_et_zdc_20t40->Scale (0.2 / h_mb_Pb_fcal_et_zdc_20t40->Integral (), "width");
  h_jet_Pb_fcal_et_zdc_40t60->Scale (0.2 / h_jet_Pb_fcal_et_zdc_40t60->Integral (), "width");
  h_mb_Pb_fcal_et_zdc_40t60->Scale (0.2 / h_mb_Pb_fcal_et_zdc_40t60->Integral (), "width");
  h_jet_Pb_fcal_et_zdc_60t80->Scale (0.2 / h_jet_Pb_fcal_et_zdc_60t80->Integral (), "width");
  h_mb_Pb_fcal_et_zdc_60t80->Scale (0.2 / h_mb_Pb_fcal_et_zdc_60t80->Integral (), "width");
  h_jet_Pb_fcal_et_zdc_80t100->Scale (0.2 / h_jet_Pb_fcal_et_zdc_80t100->Integral (), "width");
  h_mb_Pb_fcal_et_zdc_80t100->Scale (0.2 / h_mb_Pb_fcal_et_zdc_80t100->Integral (), "width");

  h_ratio_Pb_fcal_et_zdc_0t20 = (TH1D*) h_jet_Pb_fcal_et_zdc_0t20->Clone ("h_ratio_Pb_fcal_et_zdc_0t20");
  h_ratio_Pb_fcal_et_zdc_0t20->Divide (h_mb_Pb_fcal_et_zdc_0t20);
  h_ratio_Pb_fcal_et_zdc_20t40 = (TH1D*) h_jet_Pb_fcal_et_zdc_20t40->Clone ("h_ratio_Pb_fcal_et_zdc_20t40");
  h_ratio_Pb_fcal_et_zdc_20t40->Divide (h_mb_Pb_fcal_et_zdc_20t40);
  h_ratio_Pb_fcal_et_zdc_40t60 = (TH1D*) h_jet_Pb_fcal_et_zdc_40t60->Clone ("h_ratio_Pb_fcal_et_zdc_40t60");
  h_ratio_Pb_fcal_et_zdc_40t60->Divide (h_mb_Pb_fcal_et_zdc_40t60);
  h_ratio_Pb_fcal_et_zdc_60t80 = (TH1D*) h_jet_Pb_fcal_et_zdc_60t80->Clone ("h_ratio_Pb_fcal_et_zdc_60t80");
  h_ratio_Pb_fcal_et_zdc_60t80->Divide (h_mb_Pb_fcal_et_zdc_60t80);
  h_ratio_Pb_fcal_et_zdc_80t100 = (TH1D*) h_jet_Pb_fcal_et_zdc_80t100->Clone ("h_ratio_Pb_fcal_et_zdc_80t100");
  h_ratio_Pb_fcal_et_zdc_80t100->Divide (h_mb_Pb_fcal_et_zdc_80t100);


  inFile = new TFile (Form ("%s/CentralityAnalysis/Nominal/mc16_all.root", rootPath.Data ()), "read");
  h_mb_Pb_fcal_et_sum_pPb_mc = (TH1D*) inFile->Get ("h_mb_Pb_fcal_et_run0");
  h_mb_p_fcal_et_sum_pPb_mc = (TH1D*) inFile->Get ("h_mb_p_fcal_et_run0");
  h_mb_Pb_fcal_et_sum_pPb_mc_corr = (TH1D*) inFile->Get ("h_mb_Pb_fcal_et_corr_run0");
  h_mb_p_fcal_et_sum_pPb_mc_corr = (TH1D*) inFile->Get ("h_mb_p_fcal_et_corr_run0");


  double xq[101];
  double yq[101];
  TString percs[101];
  for (int iPerc = 0; iPerc < 101; iPerc++) {
    xq[iPerc] = 0.00;
    xq[iPerc] = ((double)iPerc) / 100.;
    percs[iPerc] = Form ("%i%%", 100-iPerc);
  }


  //for (int iRun = 0; iRun < npPbRuns; iRun++) {
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
  //    std::cout << xq[i] << ", " << yq[i] << std::endl;

  //  TLine* divs = new TLine ();
  //  TLatex* tl = new TLatex ();
  //  tl->SetTextAngle (-90);
  //  tl->SetTextAlign (11);
  //  tl->SetTextFont (43);
  //  tl->SetTextSize (14);
  //  divs->SetLineStyle (2);
  //  for (int i = 1; i < 16; i++) {
  //    divs->DrawLine (yq[i], ymin, yq[i], h->GetBinContent (h->FindBin (yq[i])));
  //    tl->DrawLatex (yq[i]+0.15, std::exp (0.1*std::log (ymax/ymin)) * ymin, percs[i].Data ());
  //  }

  //  h->DrawCopy ("hist same");

  //  myText (0.20, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
  //  myText (0.20, 0.840, kBlack, Form ("Run %i (per. A)", pPbRuns[iRun]), 0.032);
  //  myText (0.65, 0.550, kBlue+3, "HLT_j50_L1J15", 0.032);
  //  myText (0.65, 0.510, kRed+1, "HLT_mb_sptrk_L1MBTS_1", 0.032);

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
  //  c->SaveAs (Form ("%s/Plots/CentralityAnalysis/run%i_zdc.pdf", workPath.Data (), pPbRuns[iRun]));
  //}



  {
    TCanvas* c = new TCanvas ("c_allpPbRuns_zdc", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_allpPbRuns_zdc_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_allpPbRuns_zdc_dPad", "", 0.0, 0.0, 1.0, 0.3);
    TPad* subPad = new TPad ("c_allpPbRuns_zdc_subPad", "", 0.6, 0.6*0.7+0.3, 1.-rMargin, 1.-tMargin*0.7);

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

    double ymin = 4e-6;
    double ymax = 2e1;


    TH1D* h = (TH1D*) h_mb_Pb_zdc_calibE_sum->Clone ("htemp");

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

    int ibin = 0;
    while (ibin < 17 && plot_percs[ibin++] != "20%");

    const double mbNorm = h->Integral (h->FindBin (plot_yq[ibin]), h->GetNbinsX ());
    h->Scale (1. / mbNorm, "width");

    ofstream cutsfile;
    cutsfile.open (Form ("%s/aux/ZdcCentCuts.dat", workPath.Data ()));
    for (int i = 0; i < 101; i++) {
      //std::cout << xq[i] << ", " << yq[i] << std::endl;
      cutsfile << std::setw (6) << percs[i] << "\t" << yq[i] << "\n";
    }
    cutsfile.close ();


    TH1D* hjet = (TH1D*) h_jet_Pb_zdc_calibE_sum->Clone ("htemp");

    const double j50Norm = hjet->Integral (hjet->FindBin (plot_yq[ibin]), hjet->GetNbinsX ());
    hjet->Scale (1. / j50Norm, "width");

    TAxis* xax = hjet->GetXaxis ();
    TAxis* yax = hjet->GetYaxis ();

    yax->SetTitle ("A.U. (normalized in 0-20%)");

    yax->SetTitleOffset (1.2 * yax->GetTitleOffset ());

    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (0);
    yax->SetTitleFont (43);
    yax->SetTitleSize (26);
    xax->SetLabelFont (43);
    xax->SetLabelSize (0);
    yax->SetLabelFont (43);
    yax->SetLabelSize (24);

    hjet->SetLineColor (kBlue+3);

    hjet->DrawCopy ("hist");
    SaferDelete (&hjet);
  

    TH1D* h_zna = nullptr;
    TH1D* h_zna_glau = nullptr;
    {
      std::ifstream f_alice (Form ("%s/aux/ALICE_ZNA.txt", workPath.Data ()));
      float xbins[1028] = {};
      float yvals[1027] = {};
      float yerrs[1027] = {};
      float x, y, ye;
      for (int i = 0; i < 1027; i++) {
        f_alice >> x >> y >> ye;
        xbins[i] = x;
        yvals[i] = y;
        yerrs[i] = ye;
      }
      f_alice >> x;
      xbins[1027] = x;

      h_zna = new TH1D ("h_alice_zna", ";E_{ZN} (TeV);Events (arb. units)", 1027, xbins);
      for (int i = 0; i < 1027; i++) {
        h_zna->SetBinContent (i+1, yvals[i]);
        h_zna->SetBinError (i+1, yerrs[i]);
      }

      double alice_xq[3];
      double alice_yq[3];
      alice_xq[0] = 0.000;
      alice_xq[1] = 0.800;
      alice_xq[2] = 1.000;
      h_zna->GetQuantiles (3, alice_yq, alice_xq);

      const double aliceNorm = 1e6;//h_zna->Integral (h_zna->FindBin (alice_yq[1]), h_zna->GetNbinsX ());
      std::cout << "aliceNorm = " << aliceNorm << std::endl;
      h_zna->Scale (1. / aliceNorm, "width");


      std::ifstream f_alice_glau (Form ("%s/aux/ALICE_ZNA_GLAU.txt", workPath.Data ()));
      for (int i = 0; i < 1027; i++) {
        f_alice_glau >> x >> y >> ye;
        xbins[i] = x;
        yvals[i] = y;
        yerrs[i] = ye;
      }
      f_alice_glau >> x;
      xbins[1027] = x;

      h_zna_glau = new TH1D ("h_alice_zna_glau", ";E_{ZN} (TeV);Events (arb. units)", 1027, xbins);
      for (int i = 0; i < 1027; i++) {
        h_zna_glau->SetBinContent (i+1, yvals[i]);
        h_zna_glau->SetBinError (i+1, yerrs[i]);
      }

      h_zna_glau->Scale (1. / aliceNorm, "width");
    }


    //const double alice_sf = 2.56/1.57;
    //{
    //  std::ifstream f_alice (Form ("%s/aux/ALICE_data.txt", workPath.Data ()));
    //  std::string dummyLine;
    //  std::getline (f_alice, dummyLine);
    //  double x, y, integral = 0;
    //  while (f_alice) {
    //    f_alice >> x >> y;
    //    if (x/(alice_sf*1e3) > zdcCentBins[4]) integral += y;
    //    g_ALICE_zna_calibE->SetPoint (g_ALICE_zna_calibE->GetN (), x/(alice_sf*1e3), y);
    //  }
    //  assert (integral > 0);
    //  ScaleGraph (g_ALICE_zna_calibE, nullptr, 1./integral);

    //  std::ifstream f_alice_subpanel (Form ("%s/aux/ALICE_data_subpanel.txt", workPath.Data ()));
    //  std::getline (f_alice_subpanel, dummyLine);
    //  while (f_alice_subpanel) {
    //    f_alice_subpanel >> x >> y;
    //    g_ALICE_zna_calibE_subpanel->SetPoint (g_ALICE_zna_calibE_subpanel->GetN (), x/(alice_sf*1e3), y);
    //  }
    //  assert (integral > 0);
    //  ScaleGraph (g_ALICE_zna_calibE_subpanel, nullptr, 1./integral);

    //}

    //TGAE* g = g_ALICE_zna_calibE;
    TGAE* g = make_graph (h_zna);
    g->SetMarkerColor (myGreen);
    g->SetLineColor (myGreen);
    g->SetMarkerStyle (kOpenCircle);
    g->Draw ("P");

    TH1D* h_g = h_zna_glau;
    h_g->SetLineColor (myLiteGreen);
    h_g->SetLineWidth (2);
    h_g->DrawCopy ("same hist");

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 16; i++) {
      divs->DrawLine (plot_yq[i], ymin, plot_yq[i], h->GetBinContent (h->FindBin (plot_yq[i])));
      tl->DrawLatex (plot_yq[i]+0.20, std::exp (0.1*std::log (ymax/ymin)) * ymin, plot_percs[i].Data ());
    }

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    myText (0.18, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.18, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.18, 0.820, kBlue+3, "HLT_j50_ion_L1J10", 0.032);
    myText (0.18, 0.780, kRed+1, "HLT_mb_sptrk_L1MBTS_1", 0.032);
    //myText (0.18, 0.740, kGreen+1, "ALICE 8.16 TeV, #Sigma#it{E}^{Pb}_{ZNA} #times 1.57/2.56", 0.032);
    myText (0.18, 0.740, myGreen, "ALICE, #Sigma#it{E}^{Pb}_{ZNA} (#times 10^{-6})", 0.032);
    myText (0.18, 0.700, myLiteGreen, "ALICE, SNM #otimes Glauber fit (#times 10^{-6})", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = (TH1D*) h_ratio_Pb_zdc_calibE_sum->Clone ("htemp");

    h->Scale (mbNorm / j50Norm);

    xax = h->GetXaxis ();
    yax = h->GetYaxis ();

    xax->SetTitle ("Calibrated #Sigma#it{E}_{ZDC}^{Pb} [TeV]");
    yax->SetTitle (Form ("#color[%i]{j50} / #color[%i]{mb_sptrk}", kBlue+3, kRed+1));

    xax->SetTitleOffset (2.4 * xax->GetTitleOffset ());
    yax->SetTitleOffset (1.2 * yax->GetTitleOffset ());
    yax->CenterTitle ();

    h->SetLineColor (kBlue+3);

    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (26);
    yax->SetTitleFont (43);
    yax->SetTitleSize (26);
    xax->SetLabelFont (43);
    xax->SetLabelSize (24);
    yax->SetLabelFont (43);
    yax->SetLabelSize (24);

    h->DrawCopy ("hist");
    SaferDelete (&h);

    divs->DrawLine (0, 1, 175, 1);


    subPad->cd ();
    subPad->SetLogy ();

    ymin = 3.2e-4;
    ymax = 1.6e0;


    h = (TH1D*) h_mb_Pb_zdc_calibE_sum->Clone ("htemp");
    h->Reset ();

    h->SetLineWidth (0);

    xax = h->GetXaxis ();
    yax = h->GetYaxis ();

    xax->SetTitleOffset (999);
    yax->SetTitleOffset (999);

    xax->SetRangeUser (0, 12);
    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (18);
    yax->SetTitleFont (43);
    yax->SetTitleSize (18);
    xax->SetLabelFont (43);
    xax->SetLabelSize (14);
    yax->SetLabelFont (43);
    yax->SetLabelSize (14);

    h->DrawCopy ("hist");
    SaferDelete (&h);


    h = (TH1D*) h_jet_Pb_zdc_calibE_sum->Clone ("htemp");

    h->Scale (1. / h->Integral (h->FindBin (plot_yq[ibin]), h->GetNbinsX ()), "width");

    h->SetLineColor (kBlue+3);

    h->DrawCopy ("hist same");
    SaferDelete (&h);

    //g = g_ALICE_zna_calibE_subpanel;
    g = make_graph (h_zna);
    g->SetMarkerColor (myGreen);
    g->SetLineColor (myGreen);
    g->SetMarkerStyle (kOpenCircle);
    g->Draw ("P");

    h_g->DrawCopy ("hist same");


    h = (TH1D*) h_mb_Pb_zdc_calibE_sum->Clone ("htemp");

    h->Scale (1. / h->Integral (h->FindBin (plot_yq[ibin]), h->GetNbinsX ()), "width");

    h->SetLineColor (kRed+1);

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
    TCanvas* c = new TCanvas ("c_allpPbRuns_fcal_et", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_allpPbRuns_fcal_et_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_allpPbRuns_fcal_et_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    double ymin = 5e-7;
    double ymax = 1e0;


    TH1D* h = (TH1D*) h_mb_Pb_fcal_et_sum->Clone ("htemp");

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

    int ibin = 0;
    while (ibin < 17 && plot_percs[ibin++] != "20%");

    ofstream cutsfile;
    cutsfile.open (Form ("%s/aux/FCalCentCuts.dat", workPath.Data ()));
    for (int i = 0; i < 101; i++) {
      //std::cout << xq[i] << ", " << yq[i] << std::endl;
      cutsfile << std::setw (6) << percs[i] << "\t" << yq[i] << "\n";
    }
    cutsfile.close ();

    const double mbNorm = h->Integral (h->FindBin (plot_yq[ibin]), h->GetNbinsX ());
    h->Scale (1. / mbNorm, "width");


    TH1D* hjet = (TH1D*) h_jet_Pb_fcal_et_sum->Clone ("htemp");

    const double j50Norm = hjet->Integral (hjet->FindBin (plot_yq[ibin]), hjet->GetNbinsX ());
    hjet->Scale (1. / j50Norm, "width");

    TAxis* xax = hjet->GetXaxis ();
    TAxis* yax = hjet->GetYaxis ();

    yax->SetTitle ("A.U. (normalized in 0-20%)");

    yax->SetTitleOffset (1.2 * yax->GetTitleOffset ());

    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (0);
    yax->SetTitleFont (43);
    yax->SetTitleSize (26);
    xax->SetLabelFont (43);
    xax->SetLabelSize (0);
    yax->SetLabelFont (43);
    yax->SetLabelSize (24);

    hjet->SetLineColor (kBlue+3);

    hjet->DrawCopy ("hist");
    SaferDelete (&hjet);


    h->DrawCopy ("hist same");

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 16; i++) {
      divs->DrawLine (plot_yq[i], ymin, plot_yq[i], h->GetBinContent (h->FindBin (plot_yq[i])));
      tl->DrawLatex (plot_yq[i]+0.20, std::exp (0.1*std::log (ymax/ymin)) * ymin, plot_percs[i].Data ());
    }

    SaferDelete (&h);


    myText (0.65, 0.890, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.65, 0.850, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.810, kBlack, "All runs", 0.032);
    myText (0.65, 0.770, kBlue+3, "HLT_j50_ion_L1J10", 0.032);
    myText (0.65, 0.730, kRed+1, "HLT_mb_sptrk_L1MBTS_1", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;


    h = (TH1D*) h_ratio_Pb_fcal_et_sum->Clone ("htemp");  

    h->Scale (mbNorm / j50Norm);

    xax = h->GetXaxis ();
    yax = h->GetYaxis ();

    xax->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    yax->SetTitle (Form ("#color[%i]{j50} / #color[%i]{mb_sptrk}", kBlue+3, kRed+1));

    xax->SetTitleOffset (2.4 * xax->GetTitleOffset ());
    yax->SetTitleOffset (1.2 * yax->GetTitleOffset ());
    yax->CenterTitle ();

    h->SetLineColor (kBlue+3);

    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (26);
    yax->SetTitleFont (43);
    yax->SetTitleSize (26);
    xax->SetLabelFont (43);
    xax->SetLabelSize (24);
    yax->SetLabelFont (43);
    yax->SetLabelSize (24);

    h->DrawCopy ("hist");
    SaferDelete (&h);

    divs->DrawLine (-30, 1, 220, 1);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allpPbRuns_fcal_et.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_pPb_fcal_et_dataMC", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_pPb_fcal_et_dataMC_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_pPb_fcal_et_dataMC_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    double ymin = 5e-7;
    double ymax = 1e0;


    TH1D* hmb = (TH1D*) h_mb_Pb_fcal_et_sum->Clone ("htemp");

    hmb->SetLineColor (kRed+1);

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
    hmb->GetQuantiles (17, plot_yq, plot_xq);

    int ibin = 0;
    while (ibin < 17 && plot_percs[ibin++] != "20%");

    const double mbNorm = hmb->Integral (hmb->FindBin (plot_yq[ibin]), hmb->GetNbinsX ());
    hmb->Scale (1. / mbNorm, "width");


    TH1D* hjet = (TH1D*) h_jet_Pb_fcal_et_sum->Clone ("htemp");

    const double j50Norm = hjet->Integral (hjet->FindBin (plot_yq[ibin]), hjet->GetNbinsX ());
    hjet->Scale (1. / j50Norm, "width");

    TAxis* xax = hjet->GetXaxis ();
    TAxis* yax = hjet->GetYaxis ();

    yax->SetTitle ("A.U. (normalized in 0-20%)");

    yax->SetTitleOffset (1.2 * yax->GetTitleOffset ());

    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (0);
    yax->SetTitleFont (43);
    yax->SetTitleSize (26);
    xax->SetLabelFont (43);
    xax->SetLabelSize (0);
    yax->SetLabelFont (43);
    yax->SetLabelSize (24);

    hjet->SetLineColor (kBlue+3);

    hjet->DrawCopy ("hist");


    TH1D* hmcc = (TH1D*) h_mb_Pb_fcal_et_sum_pPb_mc_corr->Clone ("hmc");
    const double mccNorm = hmcc->Integral (hmcc->FindBin (plot_yq[ibin]), hmcc->GetNbinsX ());
    hmcc->Scale (1. / mccNorm, "width");

    hmcc->SetLineColor (myLiteBlue);

    hmcc->DrawCopy ("hist same");


    TH1D* hmc = (TH1D*) h_mb_Pb_fcal_et_sum_pPb_mc->Clone ("hmc");
    const double mcNorm = hmc->Integral (hmc->FindBin (plot_yq[ibin]), hmc->GetNbinsX ());
    hmc->Scale (1. / mcNorm, "width");

    hmc->SetLineColor (myLiteRed);

    hmc->DrawCopy ("hist same");


    hmb->DrawCopy ("hist same");

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 16; i++) {
      divs->DrawLine (plot_yq[i], ymin, plot_yq[i], hmb->GetBinContent (hmb->FindBin (plot_yq[i])));
      tl->DrawLatex (plot_yq[i]+0.20, std::exp (0.1*std::log (ymax/ymin)) * ymin, plot_percs[i].Data ());
    }


    myText (0.65, 0.890, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.65, 0.850, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.810, kBlack, "All runs", 0.032);
    myText (0.65, 0.770, kBlue+3, "HLT_j50_ion_L1J10", 0.032);
    myText (0.65, 0.730, kRed+1, "HLT_mb_sptrk_L1MBTS_1", 0.032);
    myText (0.65, 0.690, myLiteRed, "Pythia + Overlay", 0.032);
    myText (0.65, 0.650, myLiteBlue, "Overlay only", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;


    TH1D* hjet_rat = (TH1D*) hjet->Clone ("hjet_ratio");  
    hjet_rat->Divide (hmc);
    

    xax = hjet_rat->GetXaxis ();
    yax = hjet_rat->GetYaxis ();

    xax->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    yax->SetTitle (Form ("Ratio to #color[%i]{MC+Overlay}", myLiteRed));

    xax->SetTitleOffset (2.4 * xax->GetTitleOffset ());
    yax->SetTitleOffset (1.2 * yax->GetTitleOffset ());
    yax->CenterTitle ();

    hjet_rat->SetLineColor (kBlue+3);

    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (26);
    yax->SetTitleFont (43);
    yax->SetTitleSize (26);
    xax->SetLabelFont (43);
    xax->SetLabelSize (24);
    yax->SetLabelFont (43);
    yax->SetLabelSize (24);

    hjet_rat->DrawCopy ("hist");

    TH1D* hmb_rat = (TH1D*) hmb->Clone ("hmb_ratio");
    hmb_rat->Divide (hmc);

    hmb_rat->DrawCopy ("hist same");

    TH1D* hmcc_rat = (TH1D*) hmcc->Clone ("hmcc_rat");
    hmcc_rat->Divide (hmc);

    hmcc_rat->DrawCopy ("hist same");
    

    //TFile* mcResampleFile = new TFile (Form ("%s/aux/MCResampling.root", workPath.Data ()), "recreate");
    //double a_mb = 0, a_jet = 0;
    //for (int iX = 1; iX <= hmb_rat->GetNbinsX (); iX++) {
    //  if (hmb_rat->GetBinCenter (iX) > 0 && hmb_rat->GetBinCenter (iX) < 150)
    //    a_mb = std::fmax (a_mb, hmb_rat->GetBinContent (iX));
    //}
    //for (int iX = 1; iX <= hjet_rat->GetNbinsX (); iX++) {
    //  if (hjet_rat->GetBinCenter (iX) > 0 && hjet_rat->GetBinCenter (iX) < 150)
    //    a_jet = std::fmax (a_jet, hjet_rat->GetBinContent (iX));
    //}
    //hmb_rat->Scale (1./a_mb);
    //hjet_rat->Scale (1./a_jet);
    //hmb_rat->Write ();
    //hjet_rat->Write ();
    //mcResampleFile->Write ();
    //mcResampleFile->Close ();

    divs->DrawLine (-30, 1, 220, 1);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/pPb_fcal_et_dataMC.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_allppRuns_fcal_et", "", 800, 1000);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_allppRuns_fcal_et_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_allppRuns_fcal_et_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    double ymin = 5e-8;
    double ymax = 1e0;


    TH1D* h = (TH1D*) h_mb_p_fcal_et_sum->Clone ("htemp");

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

    int ibin = 0;
    while (ibin < 17 && plot_percs[ibin++] != "20%");

    ofstream cutsfile;
    cutsfile.open (Form ("%s/aux/ppMixCuts.dat", workPath.Data ()));
    for (int i = 0; i < 101; i++) {
      //std::cout << xq[i] << ", " << yq[i] << std::endl;
      cutsfile << std::setw (6) << percs[i] << "\t" << yq[i] << "\n";
    }
    cutsfile.close ();

    const double mbNorm = h->Integral (h->FindBin (plot_yq[ibin]), h->GetNbinsX ());
    h->Scale (1. / mbNorm, "width");


    TH1D* hjet = (TH1D*) h_jet_p_fcal_et_sum->Clone ("htemp");

    const double j50Norm = hjet->Integral (hjet->FindBin (plot_yq[ibin]), hjet->GetNbinsX ());
    hjet->Scale (1. / j50Norm, "width");

    TAxis* xax = hjet->GetXaxis ();
    TAxis* yax = hjet->GetYaxis ();

    xax->SetTitle ("#Sigma#it{E}_{T}^{FCal} [GeV]");
    yax->SetTitle ("A.U. (normalized in 0-20%)");

    yax->SetTitleOffset (1.2 * yax->GetTitleOffset ());

    xax->SetRangeUser (-30, 140);
    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (26);
    yax->SetTitleFont (43);
    yax->SetTitleSize (26);
    xax->SetLabelFont (43);
    xax->SetLabelSize (24);
    yax->SetLabelFont (43);
    yax->SetLabelSize (24);

    hjet->SetLineColor (kBlue+3);

    hjet->DrawCopy ("hist");
    SaferDelete (&hjet);


    h->DrawCopy ("hist same");

    TLine* divs = new TLine ();
    TLatex* tl = new TLatex ();
    tl->SetTextAngle (-90);
    tl->SetTextAlign (11);
    tl->SetTextFont (43);
    tl->SetTextSize (14);
    divs->SetLineStyle (2);
    for (int i = 1; i < 16; i++) {
      divs->DrawLine (plot_yq[i], ymin, plot_yq[i], h->GetBinContent (h->FindBin (plot_yq[i])));
      tl->DrawLatex (plot_yq[i]+0.20, std::exp (0.1*std::log (ymax/ymin)) * ymin, plot_percs[i].Data ());
    }

    SaferDelete (&h);

 
    myText (0.65, 0.890, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.65, 0.850, kBlack, "#it{pp}, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.810, kBlack, "All runs", 0.032);
    myText (0.65, 0.770, kBlue+3, "HLT_j50_L1J15", 0.032);
    myText (0.65, 0.730, kRed+1, "HLT_mb_sptrk", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = (TH1D*) h_ratio_p_fcal_et_sum->Clone ("htemp");  

    h->Scale (mbNorm / j50Norm);

    xax = h->GetXaxis ();
    yax = h->GetYaxis ();

    xax->SetTitle ("#Sigma#it{E}_{T}^{FCal, A+C} [GeV]");
    yax->SetTitle (Form ("#color[%i]{j50} / #color[%i]{mb_sptrk}", kBlue+3, kRed+1));

    xax->SetTitleOffset (2.4 * xax->GetTitleOffset ());
    yax->SetTitleOffset (1.2 * yax->GetTitleOffset ());
    yax->CenterTitle ();

    h->SetLineColor (kBlue+3);

    xax->SetRangeUser (-30, 140);
    yax->SetRangeUser (ymin, ymax);

    xax->SetTitleFont (43);
    xax->SetTitleSize (26);
    yax->SetTitleFont (43);
    yax->SetTitleSize (26);
    xax->SetLabelFont (43);
    xax->SetLabelSize (24);
    yax->SetLabelFont (43);
    yax->SetLabelSize (24);

    h->DrawCopy ("hist");
    SaferDelete (&h);

    divs->DrawLine (-30, 1, 140, 1);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allppRuns_fcal_et.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_allpPbRuns_fcal_et_zdc_0t20", "", 800, 800);

    const double bMargin = 0.11;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    c->SetBottomMargin (bMargin);
    c->SetLeftMargin (lMargin);
    c->SetRightMargin (rMargin);
    c->SetTopMargin (tMargin);

    c->Draw ();


    c->cd ();
    //c->SetLogy ();

    TH1D* h = (TH1D*) h_mb_Pb_fcal_et_sum->Clone ("htemp");
    h->Scale (1./h->Integral (), "width");
    h->Rebin (2);
    h->Scale (0.5);

    TGAE* g = make_graph (h);

    //TH1D* h = (TH1D*) h_jet_Pb_fcal_et_sum->Clone ("htemp");
    g->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
    g->GetYaxis ()->SetTitle ("A.U.");

    g->GetXaxis ()->SetTitleOffset (0.8 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.2 * g->GetYaxis ()->GetTitleOffset ());

    //double ymin = 5e-9;
    //double ymax = 1e0;
    double ymin = 0;
    double ymax = 0.038;

    //g->SetLineColor (manyColors[0]);
    g->SetLineColor (kBlack);
    g->SetLineWidth (3);

    g->GetYaxis ()->SetRangeUser (ymin, ymax);

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (26);
    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (26);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (24);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (24);

    ((TGAE*) g->Clone ())->Draw ("AC");
    //h->DrawCopy ("hist");
    SaferDelete (&h);
    SaferDelete (&g);

    //h = (TH1D*) h_mb_Pb_fcal_et_zdc_0t20->Clone ("htemp");
    //h->Rebin (2);
    //h->Scale (0.5);
    //h->SetLineColor (manyColors[0]);
    //h->DrawCopy ("hist same");
    //SaferDelete (&h);
    //h = (TH1D*) h_mb_Pb_fcal_et_zdc_20t40->Clone ("htemp");
    //h->Rebin (2);
    //h->Scale (0.5);
    //h->SetLineColor (manyColors[2]);
    //h->DrawCopy ("hist same");
    //SaferDelete (&h);
    //h = (TH1D*) h_mb_Pb_fcal_et_zdc_40t60->Clone ("htemp");
    //h->Rebin (2);
    //h->Scale (0.5);
    //h->SetLineColor (manyColors[4]);
    //h->DrawCopy ("hist same");
    //SaferDelete (&h);
    //h = (TH1D*) h_mb_Pb_fcal_et_zdc_60t80->Clone ("htemp");
    //h->Rebin (2);
    //h->Scale (0.5);
    //h->SetLineColor (manyColors[6]);
    //h->DrawCopy ("hist same");
    //SaferDelete (&h);
    //h = (TH1D*) h_mb_Pb_fcal_et_zdc_80t100->Clone ("htemp");
    //h->Rebin (2);
    //h->Scale (0.5);
    //h->SetLineColor (manyColors[8]);
    //h->DrawCopy ("hist same");
    //SaferDelete (&h);

    h = (TH1D*) h_mb_Pb_fcal_et_zdc_0t20->Clone ("htemp");
    h->Rebin (2);
    h->Scale (0.5);
    g = make_graph (h);
    ResetTGAEErrors (g);
    ResetXErrors (g);
    myDraw (g, manyColors[0], kDot, 0, 1, 3, "C");
    SaferDelete (&g);

    h = (TH1D*) h_mb_Pb_fcal_et_zdc_20t40->Clone ("htemp");
    h->Rebin (2);
    h->Scale (0.5);
    g = make_graph (h);
    ResetTGAEErrors (g);
    ResetXErrors (g);
    myDraw (g, manyColors[2], kDot, 0, 1, 3, "C");
    SaferDelete (&g);

    h = (TH1D*) h_mb_Pb_fcal_et_zdc_40t60->Clone ("htemp");
    h->Rebin (2);
    h->Scale (0.5);
    g = make_graph (h);
    ResetTGAEErrors (g);
    ResetXErrors (g);
    myDraw (g, manyColors[4], kDot, 0, 1, 3, "C");
    SaferDelete (&g);

    h = (TH1D*) h_mb_Pb_fcal_et_zdc_60t80->Clone ("htemp");
    h->Rebin (2);
    h->Scale (0.5);
    g = make_graph (h);
    ResetTGAEErrors (g);
    ResetXErrors (g);
    myDraw (g, manyColors[6], kDot, 0, 1, 3, "C");
    SaferDelete (&g);

    h = (TH1D*) h_mb_Pb_fcal_et_zdc_80t100->Clone ("htemp");
    h->Rebin (2);
    h->Scale (0.5);
    g = make_graph (h);
    ResetTGAEErrors (g);
    ResetXErrors (g);
    myDraw (g, manyColors[8], kDot, 0, 1, 3, "C");
    SaferDelete (&g);


    myText (0.59, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.59, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.59, 0.820, kBlack, "HLT_mb_sptrk_L1MBTS_1", 0.032);
    myText (0.56, 0.740, kBlack,        "#bf{P(#Sigma#it{E}_{T}^{FCal,Pb})}", 0.026);
    myText (0.56, 0.700, manyColors[0], "#bf{0.2 #times P(#Sigma#it{E}_{T}^{FCal,Pb} | Zdc 0-20%)}", 0.026);
    myText (0.56, 0.660, manyColors[2], "#bf{0.2 #times P(#Sigma#it{E}_{T}^{FCal,Pb} | Zdc 20-40%)}", 0.026);
    myText (0.56, 0.620, manyColors[4], "#bf{0.2 #times P(#Sigma#it{E}_{T}^{FCal,Pb} | Zdc 40-60%)}", 0.026);
    myText (0.56, 0.580, manyColors[6], "#bf{0.2 #times P(#Sigma#it{E}_{T}^{FCal,Pb} | Zdc 60-80%)}", 0.026);
    myText (0.56, 0.540, manyColors[8], "#bf{0.2 #times P(#Sigma#it{E}_{T}^{FCal,Pb} | Zdc 80-100%)}", 0.026);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allruns_fcal_et_zdcBinned.pdf", workPath.Data ()));
  }




}

#endif
