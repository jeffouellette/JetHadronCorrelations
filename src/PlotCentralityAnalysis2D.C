#ifndef __PlotCentralityAnalysis2D_C__
#define __PlotCentralityAnalysis2D_C__

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
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
  while (iZdc20Percent < numZdcCentBins && zdcCentPercs[iZdc20Percent] > 20) iZdc20Percent++;
  int iFcal20Percent = 0;
  while (iFcal20Percent < numFcalCentBins && fcalCentPercs[iFcal20Percent] > 100) iFcal20Percent++;


  TH2D* h2_jet_Pb_fcal_et_zdc_calibE = nullptr;
  TH2D* h2_mb_Pb_fcal_et_zdc_calibE = nullptr;

  inFile = new TFile (Form ("%s/CentralityAnalysis/Nominal/allpPbRuns.root", rootPath.Data ()), "read");

  h2_jet_Pb_fcal_et_zdc_calibE = (TH2D*) inFile->Get ("h2_jet_Pb_fcal_et_zdc_calibE");//h3_jet_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");
  h2_mb_Pb_fcal_et_zdc_calibE = (TH2D*) inFile->Get ("h2_mb_Pb_fcal_et_zdc_calibE");//h3_mb_Pb_fcal_et_Pb_q2_Pb_zdc_calibE->Project3D ("zx");


  {
    TCanvas* c = new TCanvas ("c_J50_zdc_fcalet_correlation", "", 1000, 812);

    gPad->SetLogz ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_jet_Pb_fcal_et_zdc_calibE->Clone ("htemp");
    h->RebinX (4);
    h->RebinY (4);
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

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

    h->GetXaxis ()->SetTitleOffset (1.2 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    //h->GetXaxis ()->SetTitleFont (43);
    //h->GetXaxis ()->SetTitleSize (26);
    //h->GetYaxis ()->SetTitleFont (43);
    //h->GetYaxis ()->SetTitleSize (26);
    //h->GetZaxis ()->SetTitleFont (43);
    //h->GetZaxis ()->SetTitleSize (26);

    //h->GetXaxis ()->SetLabelFont (43);
    //h->GetXaxis ()->SetLabelSize (22);
    //h->GetYaxis ()->SetLabelFont (43);
    //h->GetYaxis ()->SetLabelSize (22);
    //h->GetZaxis ()->SetLabelFont (43);
    //h->GetZaxis ()->SetLabelSize (22);

    //h->GetXaxis ()->SetAxisColor (kWhite);
    //h->GetXaxis ()->SetLabelColor (kWhite);
    //h->GetXaxis ()->SetTitleColor (kWhite);
    //h->GetYaxis ()->SetAxisColor (kWhite);
    //h->GetYaxis ()->SetLabelColor (kWhite);
    //h->GetYaxis ()->SetTitleColor (kWhite);
    //h->GetZaxis ()->SetAxisColor (kWhite);
    //h->GetZaxis ()->SetLabelColor (kWhite);
    //h->GetZaxis ()->SetTitleColor (kWhite);

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

    myText (0.15, 0.96, kWhite, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.50, 0.96, kWhite, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, J50 Trigger", 0.032);

    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/J50_zdc_fcalet_correlation.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("MB_zdc_fcalet_correlation", "", 1000, 812);

    gPad->SetLogz ();

    gPad->SetBottomMargin (0.13);
    gPad->SetLeftMargin (0.13);
    gPad->SetRightMargin (0.15);
    gPad->SetTopMargin (0.06);

    TH2D* h = (TH2D*) h2_mb_Pb_fcal_et_zdc_calibE->Clone ("htemp");
    TGraphErrors* g_px = TProfX2TGE (h->ProfileX ("h_px"));
    TGraphErrors* g_py = TProfY2TGE (h->ProfileY ("h_py"));

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

    h->GetXaxis ()->SetTitleOffset (1.2 * h->GetXaxis ()->GetTitleOffset ());
    h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

    //h->GetXaxis ()->SetTitleFont (43);
    //h->GetXaxis ()->SetTitleSize (26);
    //h->GetYaxis ()->SetTitleFont (43);
    //h->GetYaxis ()->SetTitleSize (26);
    //h->GetZaxis ()->SetTitleFont (43);
    //h->GetZaxis ()->SetTitleSize (26);

    //h->GetXaxis ()->SetLabelFont (43);
    //h->GetXaxis ()->SetLabelSize (22);
    //h->GetYaxis ()->SetLabelFont (43);
    //h->GetYaxis ()->SetLabelSize (22);
    //h->GetZaxis ()->SetLabelFont (43);
    //h->GetZaxis ()->SetLabelSize (22);

    //h->GetXaxis ()->SetAxisColor (kWhite);
    //h->GetXaxis ()->SetLabelColor (kWhite);
    //h->GetXaxis ()->SetTitleColor (kWhite);
    //h->GetYaxis ()->SetAxisColor (kWhite);
    //h->GetYaxis ()->SetLabelColor (kWhite);
    //h->GetYaxis ()->SetTitleColor (kWhite);
    //h->GetZaxis ()->SetAxisColor (kWhite);
    //h->GetZaxis ()->SetLabelColor (kWhite);
    //h->GetZaxis ()->SetTitleColor (kWhite);

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

    myText (0.15, 0.96, kWhite, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.50, 0.96, kWhite, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, MinBias Trigger", 0.032);
   
    c->SaveAs (Form ("%s/Plots/CentralityAnalysis/MB_zdc_fcalet_correlation.pdf", workPath.Data ()));
  }

}

#endif
