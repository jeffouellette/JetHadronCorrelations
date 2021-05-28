#ifndef __PlotFCalJetBiasStudy_C__
#define __PlotFCalJetBiasStudy_C__

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>

#include <Utilities.h>

#include "Params.h"
#include "CentralityDefs.h"
#include "LocalUtilities.h"

using namespace std;
using namespace JetHadronCorrelations;


void PlotFCalJetBiasStudy () { 

  TFile* inFile = new TFile (Form ("%s/CentralityAnalysis/Nominal/allRuns.root", rootPath.Data ()), "read");


  TH3D* h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE");
  TH3D* h3_mb_Pb_fcal_et_Pb_q3_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_mb_Pb_fcal_et_Pb_q3_Pb_zdc_calibE");
  TH3D* h3_mb_Pb_fcal_et_Pb_q4_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_mb_Pb_fcal_et_Pb_q4_Pb_zdc_calibE");

  TH3D* h3_mb_Pb_fcal_et_p_q2_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_mb_Pb_fcal_et_p_q2_Pb_zdc_calibE");
  TH3D* h3_mb_Pb_fcal_et_p_q3_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_mb_Pb_fcal_et_p_q3_Pb_zdc_calibE");
  TH3D* h3_mb_Pb_fcal_et_p_q4_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_mb_Pb_fcal_et_p_q4_Pb_zdc_calibE");

  TH3D* h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE");
  TH3D* h3_jet_Pb_fcal_et_Pb_q3_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_jet_Pb_fcal_et_Pb_q3_Pb_zdc_calibE");
  TH3D* h3_jet_Pb_fcal_et_Pb_q4_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_jet_Pb_fcal_et_Pb_q4_Pb_zdc_calibE");

  TH3D* h3_jet_Pb_fcal_et_p_q2_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_jet_Pb_fcal_et_p_q2_Pb_zdc_calibE");
  TH3D* h3_jet_Pb_fcal_et_p_q3_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_jet_Pb_fcal_et_p_q3_Pb_zdc_calibE");
  TH3D* h3_jet_Pb_fcal_et_p_q4_Pb_zdc_calibE = (TH3D*) inFile->Get ("h3_jet_Pb_fcal_et_p_q4_Pb_zdc_calibE");


  TH2D* h2_jet_Pb_fcal_et_zdc_calibE = (TH2D*) h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");
  TH2D* h2_mb_Pb_fcal_et_zdc_calibE = (TH2D*) h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");

  TH2D* h2_jet_Pb_fcal_et_Pb_q2 = (TH2D*) h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_mb_Pb_fcal_et_Pb_q2 = (TH2D*) h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_jet_Pb_fcal_et_Pb_q3 = (TH2D*) h3_jet_Pb_fcal_et_Pb_q3_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_mb_Pb_fcal_et_Pb_q3 = (TH2D*) h3_mb_Pb_fcal_et_Pb_q3_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_jet_Pb_fcal_et_Pb_q4 = (TH2D*) h3_jet_Pb_fcal_et_Pb_q4_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_mb_Pb_fcal_et_Pb_q4 = (TH2D*) h3_mb_Pb_fcal_et_Pb_q4_Pb_zdc_calibE->Project3D ("xy");

  TH2D* h2_jet_Pb_fcal_et_p_q2 = (TH2D*) h3_jet_Pb_fcal_et_p_q2_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_mb_Pb_fcal_et_p_q2 = (TH2D*) h3_mb_Pb_fcal_et_p_q2_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_jet_Pb_fcal_et_p_q3 = (TH2D*) h3_jet_Pb_fcal_et_p_q3_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_mb_Pb_fcal_et_p_q3 = (TH2D*) h3_mb_Pb_fcal_et_p_q3_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_jet_Pb_fcal_et_p_q4 = (TH2D*) h3_jet_Pb_fcal_et_p_q4_Pb_zdc_calibE->Project3D ("xy");
  TH2D* h2_mb_Pb_fcal_et_p_q4 = (TH2D*) h3_mb_Pb_fcal_et_p_q4_Pb_zdc_calibE->Project3D ("xy");

  TH2D* h2_jet_Pb_zdc_calibE_p_q2 = (TH2D*) h3_jet_Pb_fcal_et_p_q2_Pb_zdc_calibE->Project3D ("zy");
  TH2D* h2_mb_Pb_zdc_calibE_p_q2 = (TH2D*) h3_mb_Pb_fcal_et_p_q2_Pb_zdc_calibE->Project3D ("zy");
  TH2D* h2_jet_Pb_zdc_calibE_p_q3 = (TH2D*) h3_jet_Pb_fcal_et_p_q3_Pb_zdc_calibE->Project3D ("zy");
  TH2D* h2_mb_Pb_zdc_calibE_p_q3 = (TH2D*) h3_mb_Pb_fcal_et_p_q3_Pb_zdc_calibE->Project3D ("zy");
  TH2D* h2_jet_Pb_zdc_calibE_p_q4 = (TH2D*) h3_jet_Pb_fcal_et_p_q4_Pb_zdc_calibE->Project3D ("zy");
  TH2D* h2_mb_Pb_zdc_calibE_p_q4 = (TH2D*) h3_mb_Pb_fcal_et_p_q4_Pb_zdc_calibE->Project3D ("zy");

  TH1D* h_jet_p_q2_PbFCalCentral = (TH1D*) h2_jet_Pb_fcal_et_p_q2->ProjectionX ("h_jet_p_q2_PbFCalCentral", h2_jet_Pb_fcal_et_p_q2->GetYaxis ()->FindBin (fcalCentBins[1]), h2_jet_Pb_fcal_et_p_q2->GetYaxis ()->GetNbins ());
  h_jet_p_q2_PbFCalCentral->Scale (1/h_jet_p_q2_PbFCalCentral->Integral (), "width");
  TH1D* h_mb_p_q2_PbFCalCentral = (TH1D*) h2_mb_Pb_fcal_et_p_q2->ProjectionX ("h_mb_p_q2_PbFCalCentral", h2_mb_Pb_fcal_et_p_q2->GetYaxis ()->FindBin (fcalCentBins[1]), h2_mb_Pb_fcal_et_p_q2->GetYaxis ()->GetNbins ());
  h_mb_p_q2_PbFCalCentral->Scale (1/h_mb_p_q2_PbFCalCentral->Integral (), "width");
  TH1D* h_ratio_p_q2_PbFCalCentral = (TH1D*) h_jet_p_q2_PbFCalCentral->Clone ("h_ratio_p_q2_PbFCalCentral");
  h_ratio_p_q2_PbFCalCentral->Divide (h_mb_p_q2_PbFCalCentral);
  TH1D* h_jet_p_q2_PbZdcCentral = (TH1D*) h2_jet_Pb_zdc_calibE_p_q2->ProjectionX ("h_jet_p_q2_PbZdcCentral", h2_jet_Pb_zdc_calibE_p_q2->GetYaxis ()->FindBin (fcalCentBins[1]), h2_jet_Pb_zdc_calibE_p_q2->GetYaxis ()->GetNbins ());
  h_jet_p_q2_PbZdcCentral->Scale (1/h_jet_p_q2_PbZdcCentral->Integral (), "width");
  TH1D* h_mb_p_q2_PbZdcCentral = (TH1D*) h2_mb_Pb_zdc_calibE_p_q2->ProjectionX ("h_mb_p_q2_PbZdcCentral", h2_mb_Pb_zdc_calibE_p_q2->GetYaxis ()->FindBin (fcalCentBins[1]), h2_mb_Pb_zdc_calibE_p_q2->GetYaxis ()->GetNbins ());
  h_mb_p_q2_PbZdcCentral->Scale (1/h_mb_p_q2_PbZdcCentral->Integral (), "width");
  TH1D* h_ratio_p_q2_PbZdcCentral = (TH1D*) h_jet_p_q2_PbZdcCentral->Clone ("h_ratio_p_q2_PbZdcCentral");
  h_ratio_p_q2_PbZdcCentral->Divide (h_mb_p_q2_PbZdcCentral);

  TH1D* h_jet_p_q3_PbFCalCentral = (TH1D*) h2_jet_Pb_fcal_et_p_q3->ProjectionX ("h_jet_p_q3_PbFCalCentral", h2_jet_Pb_fcal_et_p_q3->GetYaxis ()->FindBin (fcalCentBins[1]), h2_jet_Pb_fcal_et_p_q3->GetYaxis ()->GetNbins ());
  h_jet_p_q3_PbFCalCentral->Scale (1/h_jet_p_q3_PbFCalCentral->Integral (), "width");
  TH1D* h_mb_p_q3_PbFCalCentral = (TH1D*) h2_mb_Pb_fcal_et_p_q3->ProjectionX ("h_mb_p_q3_PbFCalCentral", h2_mb_Pb_fcal_et_p_q3->GetYaxis ()->FindBin (fcalCentBins[1]), h2_mb_Pb_fcal_et_p_q3->GetYaxis ()->GetNbins ());
  h_mb_p_q3_PbFCalCentral->Scale (1/h_mb_p_q3_PbFCalCentral->Integral (), "width");
  TH1D* h_ratio_p_q3_PbFCalCentral = (TH1D*) h_jet_p_q3_PbFCalCentral->Clone ("h_ratio_p_q3_PbFCalCentral");
  h_ratio_p_q3_PbFCalCentral->Divide (h_mb_p_q3_PbFCalCentral);
  TH1D* h_jet_p_q3_PbZdcCentral = (TH1D*) h2_jet_Pb_zdc_calibE_p_q3->ProjectionX ("h_jet_p_q3_PbZdcCentral", h2_jet_Pb_zdc_calibE_p_q3->GetYaxis ()->FindBin (fcalCentBins[1]), h2_jet_Pb_zdc_calibE_p_q3->GetYaxis ()->GetNbins ());
  h_jet_p_q3_PbZdcCentral->Scale (1/h_jet_p_q3_PbZdcCentral->Integral (), "width");
  TH1D* h_mb_p_q3_PbZdcCentral = (TH1D*) h2_mb_Pb_zdc_calibE_p_q3->ProjectionX ("h_mb_p_q3_PbZdcCentral", h2_mb_Pb_zdc_calibE_p_q3->GetYaxis ()->FindBin (fcalCentBins[1]), h2_mb_Pb_zdc_calibE_p_q3->GetYaxis ()->GetNbins ());
  h_mb_p_q3_PbZdcCentral->Scale (1/h_mb_p_q3_PbZdcCentral->Integral (), "width");
  TH1D* h_ratio_p_q3_PbZdcCentral = (TH1D*) h_jet_p_q3_PbZdcCentral->Clone ("h_ratio_p_q3_PbZdcCentral");
  h_ratio_p_q3_PbZdcCentral->Divide (h_mb_p_q3_PbZdcCentral);

  TH1D* h_jet_p_q4_PbFCalCentral = (TH1D*) h2_jet_Pb_fcal_et_p_q4->ProjectionX ("h_jet_p_q4_PbFCalCentral", h2_jet_Pb_fcal_et_p_q4->GetYaxis ()->FindBin (fcalCentBins[1]), h2_jet_Pb_fcal_et_p_q4->GetYaxis ()->GetNbins ());
  h_jet_p_q4_PbFCalCentral->Scale (1/h_jet_p_q4_PbFCalCentral->Integral (), "width");
  TH1D* h_mb_p_q4_PbFCalCentral = (TH1D*) h2_mb_Pb_fcal_et_p_q4->ProjectionX ("h_mb_p_q4_PbFCalCentral", h2_mb_Pb_fcal_et_p_q4->GetYaxis ()->FindBin (fcalCentBins[1]), h2_mb_Pb_fcal_et_p_q4->GetYaxis ()->GetNbins ());
  h_mb_p_q4_PbFCalCentral->Scale (1/h_mb_p_q4_PbFCalCentral->Integral (), "width");
  TH1D* h_ratio_p_q4_PbFCalCentral = (TH1D*) h_jet_p_q4_PbFCalCentral->Clone ("h_ratio_p_q4_PbFCalCentral");
  h_ratio_p_q4_PbFCalCentral->Divide (h_mb_p_q4_PbFCalCentral);
  TH1D* h_jet_p_q4_PbZdcCentral = (TH1D*) h2_jet_Pb_zdc_calibE_p_q4->ProjectionX ("h_jet_p_q4_PbZdcCentral", h2_jet_Pb_zdc_calibE_p_q4->GetYaxis ()->FindBin (fcalCentBins[1]), h2_jet_Pb_zdc_calibE_p_q4->GetYaxis ()->GetNbins ());
  h_jet_p_q4_PbZdcCentral->Scale (1/h_jet_p_q4_PbZdcCentral->Integral (), "width");
  TH1D* h_mb_p_q4_PbZdcCentral = (TH1D*) h2_mb_Pb_zdc_calibE_p_q4->ProjectionX ("h_mb_p_q4_PbZdcCentral", h2_mb_Pb_zdc_calibE_p_q4->GetYaxis ()->FindBin (fcalCentBins[1]), h2_mb_Pb_zdc_calibE_p_q4->GetYaxis ()->GetNbins ());
  h_mb_p_q4_PbZdcCentral->Scale (1/h_mb_p_q4_PbZdcCentral->Integral (), "width");
  TH1D* h_ratio_p_q4_PbZdcCentral = (TH1D*) h_jet_p_q4_PbZdcCentral->Clone ("h_ratio_p_q4_PbZdcCentral");
  h_ratio_p_q4_PbZdcCentral->Divide (h_mb_p_q4_PbZdcCentral);



  {
    TCanvas* c = new TCanvas ("c_J50_zdc_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_zdc_calibE->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

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

    myText (0.15, 0.96, kWhite, "All runs, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_zdc_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("MB_zdc_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_zdc_calibE->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

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

    myText (0.15, 0.96, kWhite, "All runs, MinBias Trigger", 0.032);
   
    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_zdc_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_J50_Pb_q2_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_Pb_q2->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{2}^{Pb}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_Pb_q2_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("MB_Pb_q2_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_Pb_q2->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{2}^{Pb}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, MinBias Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_Pb_q2_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_J50_Pb_q3_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_Pb_q3->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{3}^{Pb}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_Pb_q3_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("MB_Pb_q3_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_Pb_q3->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{3}^{Pb}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, MinBias Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_Pb_q3_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_J50_Pb_q4_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_Pb_q4->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{4}^{Pb}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_Pb_q4_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("MB_Pb_q4_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_Pb_q4->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{4}^{Pb}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, MinBias Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_Pb_q4_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_J50_p_q2_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_p_q2->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{2}^{p}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_p_q2_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("MB_p_q2_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_p_q2->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{2}^{p}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, MinBias Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_p_q2_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_J50_p_q3_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_p_q3->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{3}^{p}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_p_q3_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("MB_p_q3_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_p_q3->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{3}^{p}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, MinBias Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_p_q3_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_J50_p_q4_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_p_q4->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{4}^{p}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_p_q4_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("MB_p_q4_fcalet_correlation", "", 800, 800);

    gPad->SetLogz ();
    gPad->SetFillColor (kBlue+3);
    gPad->Update ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_p_q4->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

    h->GetXaxis ()->SetTitle ("#it{q}_{4}^{p}");
    h->GetYaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal, Pb} [GeV]");
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

    myText (0.15, 0.96, kWhite, "All runs, MinBias Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_p_q4_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_p_q2", "", 800, 800);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_p_q2_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_p_q2_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    TH1D* h = (TH1D*) h_jet_p_q2_PbFCalCentral->Clone ("htemp");

    h->GetYaxis ()->SetTitle ("A.U. (each normalized)");

    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    double ymin = 5e-4;
    double ymax = 1e1;

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

    h = (TH1D*) h_mb_p_q2_PbFCalCentral->Clone ("htemp");

    h->SetLineColor (kRed+1);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    h = (TH1D*) h_jet_p_q2_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kBlue+3);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    h = (TH1D*) h_mb_p_q2_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kRed+1);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    myText (0.65, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.65, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.780, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.740, kRed+1, "MinBias Trigger", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = (TH1D*) h_ratio_p_q2_PbFCalCentral->Clone ("htemp");
    h->GetXaxis ()->SetTitle ("#it{q}_{2}^{p}");
    h->GetYaxis ()->SetTitle ("J50 / MB");

    h->GetXaxis ()->SetTitleOffset (1.8 * h->GetXaxis ()->GetTitleOffset ());
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

    h = (TH1D*) h_ratio_p_q2_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kBlue+3);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);

    TLine* tline = new TLine ();
    tline->SetLineStyle (2);
    tline->DrawLine (-30, 1, 220, 1);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allruns_p_q2_Central.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_p_q3", "", 800, 800);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_p_q3_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_p_q3_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    TH1D* h = (TH1D*) h_jet_p_q3_PbFCalCentral->Clone ("htemp");

    h->GetYaxis ()->SetTitle ("A.U. (each normalized)");

    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    double ymin = 5e-4;
    double ymax = 1e1;

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

    h = (TH1D*) h_mb_p_q3_PbFCalCentral->Clone ("htemp");

    h->SetLineColor (kRed+1);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    h = (TH1D*) h_jet_p_q3_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kBlue+3);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    h = (TH1D*) h_mb_p_q3_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kRed+1);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    myText (0.65, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.65, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.780, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.740, kRed+1, "MinBias Trigger", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = (TH1D*) h_ratio_p_q3_PbFCalCentral->Clone ("htemp");
    h->GetXaxis ()->SetTitle ("#it{q}_{3}^{p}");
    h->GetYaxis ()->SetTitle ("J50 / MB");

    h->GetXaxis ()->SetTitleOffset (1.8 * h->GetXaxis ()->GetTitleOffset ());
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

    h = (TH1D*) h_ratio_p_q3_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kBlue+3);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);

    TLine* tline = new TLine ();
    tline->SetLineStyle (2);
    tline->DrawLine (-30, 1, 220, 1);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allruns_p_q3_Central.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c_p_q4", "", 800, 800);

    const double bMargin = 0.20;
    const double lMargin = 0.11;
    const double rMargin = 0.04;
    const double tMargin = 0.04;

    TPad* uPad = new TPad ("c_p_q4_uPad", "", 0.0, 0.3, 1.0, 1.0);
    TPad* dPad = new TPad ("c_p_q4_dPad", "", 0.0, 0.0, 1.0, 0.3);

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

    TH1D* h = (TH1D*) h_jet_p_q4_PbFCalCentral->Clone ("htemp");

    h->GetYaxis ()->SetTitle ("A.U. (each normalized)");

    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    double ymin = 5e-4;
    double ymax = 1e1;

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

    h = (TH1D*) h_mb_p_q4_PbFCalCentral->Clone ("htemp");

    h->SetLineColor (kRed+1);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    h = (TH1D*) h_jet_p_q4_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kBlue+3);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    h = (TH1D*) h_mb_p_q4_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kRed+1);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);


    myText (0.65, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.036);
    myText (0.65, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.65, 0.780, kBlue+3, "J50 Trigger", 0.032);
    myText (0.65, 0.740, kRed+1, "MinBias Trigger", 0.032);


    dPad->cd ();
    dPad->SetLogy ();

    ymin = 7e-2;
    ymax = 2e1;

    h = (TH1D*) h_ratio_p_q4_PbFCalCentral->Clone ("htemp");
    h->GetXaxis ()->SetTitle ("#it{q}_{4}^{p}");
    h->GetYaxis ()->SetTitle ("J50 / MB");

    h->GetXaxis ()->SetTitleOffset (1.8 * h->GetXaxis ()->GetTitleOffset ());
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

    h = (TH1D*) h_ratio_p_q4_PbZdcCentral->Clone ("htemp");

    h->SetLineColor (kBlue+3);
    h->SetLineStyle (2);

    h->DrawCopy ("hist same");
    SaferDelete (&h);

    TLine* tline = new TLine ();
    tline->SetLineStyle (2);
    tline->DrawLine (-30, 1, 220, 1);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/allruns_p_q4_Central.pdf", workPath.Data ()));
  }


}

#endif
