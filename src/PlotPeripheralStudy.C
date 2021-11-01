#ifndef __JetHadronCorrelator_PlotPeripheralStudy_C__
#define __JetHadronCorrelator_PlotPeripheralStudy_C__

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"

#include <ArrayTemplates.h>
#include <Utilities.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;


void PlotPeripheralStudy () {

  TString inFileName = Form ("%s/PeripheralStudy/MinBias.root", rootPath.Data ());
  TFile* inFile = new TFile (inFileName.Data (), "read");

  TH1D* h_sumGaps_Pb_0nXn = (TH1D*) inFile->Get ("h_sumGaps_Pb_0nXn");
  TH1D* h_sumGaps_Pb_0n0n = (TH1D*) inFile->Get ("h_sumGaps_Pb_0n0n");
  TH1D* h_cluster_sumGaps_Pb_0nXn = (TH1D*) inFile->Get ("h_cluster_sumGaps_Pb_0nXn");
  TH1D* h_cluster_sumGaps_Pb_0n0n = (TH1D*) inFile->Get ("h_cluster_sumGaps_Pb_0n0n");
  TH1D* h_sumGaps_Pb_0nXn_30GeVJet = (TH1D*) inFile->Get ("h_sumGaps_Pb_0nXn_30GeVJet");
  TH1D* h_sumGaps_Pb_0n0n_30GeVJet = (TH1D*) inFile->Get ("h_sumGaps_Pb_0n0n_30GeVJet");
  TH1D* h_cluster_sumGaps_Pb_0nXn_30GeVJet = (TH1D*) inFile->Get ("h_cluster_sumGaps_Pb_0nXn_30GeVJet");
  TH1D* h_cluster_sumGaps_Pb_0n0n_30GeVJet = (TH1D*) inFile->Get ("h_cluster_sumGaps_Pb_0n0n_30GeVJet");

  const float totEvts = h_sumGaps_Pb_0nXn->Integral () + h_sumGaps_Pb_0n0n->Integral ();
  const float totJetEvts = h_sumGaps_Pb_0nXn_30GeVJet->Integral () + h_sumGaps_Pb_0n0n_30GeVJet->Integral ();

  {
    TCanvas* c = new TCanvas ("c", "", 800, 800);
    gPad->SetLogy ();

    TH1D* h = h_sumGaps_Pb_0nXn;
    const float evtFrac_0nXn = h->Integral () / totEvts;
    const float upcFrac_0nXn = h->Integral (h->FindBin (2.5), h->GetNbinsX ()) / totEvts;//h->Integral ();
    h->Scale (1./totEvts, "width");
    h->GetYaxis ()->SetRangeUser (1e-3, 1e2);
    h->GetXaxis ()->SetTitle ("#Sigma_{Pb (#gamma)} #Delta#eta_{trk.+cl.}");
    h->GetYaxis ()->SetTitle ("(1/N_{tot}) (dN_{evt} / d#Sigma#Delta#eta)");
    h->SetLineColor (kRed);
    h->SetMarkerColor (kRed);
    h->DrawCopy ("e1");

    h = h_sumGaps_Pb_0nXn_30GeVJet;
    const float evtFrac_0nXn_30GeVJet = h->Integral () / totEvts;
    const float upcFrac_0nXn_30GeVJet = h->Integral (h->FindBin (2.5), h->GetNbinsX ()) / totJetEvts;//h->Integral ();
    h->Rebin (2);
    h->Scale (1./totJetEvts, "width");
    h->SetLineColor (kRed);
    h->SetMarkerColor (kRed);
    h->SetMarkerStyle (kOpenCircle);
    h->DrawCopy ("e1 same");

    h = h_sumGaps_Pb_0n0n;
    const float evtFrac_0n0n = h->Integral () / totEvts;
    const float upcFrac_0n0n = h->Integral (h->FindBin (2.5), h->GetNbinsX ()) / totEvts;//h->Integral ();
    h->Scale (1./totEvts, "width");
    h->SetLineColor (kBlue);
    h->SetMarkerColor (kBlue);
    h->DrawCopy ("e1 same");

    h = h_sumGaps_Pb_0n0n_30GeVJet;
    const float evtFrac_0n0n_30GeVJet = h->Integral () / totEvts;
    const float upcFrac_0n0n_30GeVJet = h->Integral (h->FindBin (2.5), h->GetNbinsX ()) / totJetEvts;//h->Integral ();
    h->Rebin (2);
    h->Scale (1./totJetEvts, "width");
    h->SetLineColor (kBlue);
    h->SetMarkerColor (kBlue);
    h->SetMarkerStyle (kOpenCircle);
    h->DrawCopy ("e1 same");

    //h_sumGaps_Pb_0nXn->DrawCopy ("e1 same");
    //h_sumGaps_Pb_0n0n->DrawCopy ("e1 same");
    //h_sumGaps_Pb_0nXn_30GeVJet->DrawCopy ("e1 same");
    //h_sumGaps_Pb_0n0n_30GeVJet->DrawCopy ("e1 same");

    TLine* l = new TLine ();
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (2.5, 1e-3, 2.5, 0.5);

    myText (0.26, 0.885, kBlack, "#bf{#it{ATLAS}} Internal", 0.030);
    myText (0.26, 0.850, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.030);
    myText (0.26, 0.815, kBlack, "HLT_mb_sptrk_L1MBTS_1", 0.030);
    myText (0.26, 0.780, kBlack, "#Sigma#it{E}_{ZDC}^{Pb} < 18.0971 TeV (#geq 80%)", 0.030);
    myText (0.33, 0.740, kBlack, "Selection", 0.030);
    myText (0.65, 0.740, kBlack, "Frac. w/ #Sigma#Delta#eta > 2.5", 0.030);
    myLineText2 (0.33, 0.710, kRed, kFullCircle, "\"0nXn\"", 1.2, 0.030);
    myLineText2 (0.33, 0.675, kBlue, kFullCircle, "\"0n0n\"", 1.2, 0.030);
    myLineText2 (0.33, 0.640, kRed, kOpenCircle, "\"0nXn\" + #it{p}_{T}^{lead} > 30 GeV", 1.2, 0.030);
    myLineText2 (0.33, 0.605, kBlue, kOpenCircle, "\"0n0n\" + #it{p}_{T}^{lead} > 30 GeV", 1.2, 0.030);
    myText (0.65, 0.705, kBlack, Form ("%.2f%%", upcFrac_0nXn*100), 0.030);
    myText (0.65, 0.670, kBlack, Form ("%.2f%%", upcFrac_0n0n*100), 0.030);
    myText (0.65, 0.635, kBlack, Form ("%.2f%%", upcFrac_0nXn_30GeVJet*100), 0.030);
    myText (0.65, 0.600, kBlack, Form ("%.2f%%", upcFrac_0n0n_30GeVJet*100), 0.030);

    c->SaveAs ("Plots/Peripheral_SumGaps.pdf");
  }

  {
    TCanvas* c = new TCanvas ("c2", "", 800, 800);
    gPad->SetLogy ();

    TH1D* h = h_cluster_sumGaps_Pb_0nXn;
    const float evtFrac_0nXn = h->Integral () / totEvts;
    const float upcFrac_0nXn = h->Integral (h->FindBin (2.5), h->GetNbinsX ()) / totEvts;//h->Integral ();
    h->Scale (1./totEvts, "width");
    h->GetYaxis ()->SetRangeUser (1e-3, 1e2);
    h->GetXaxis ()->SetTitle ("#Sigma_{Pb (#gamma)} #Delta#eta_{cl}");
    h->GetYaxis ()->SetTitle ("(1/N_{tot}) (dN_{evt} / d#Sigma#Delta#eta)");
    h->SetLineColor (kRed);
    h->SetMarkerColor (kRed);
    h->DrawCopy ("e1 ][");

    h = h_cluster_sumGaps_Pb_0nXn_30GeVJet;
    const float evtFrac_0nXn_30GeVJet = h->Integral () / totEvts;
    const float upcFrac_0nXn_30GeVJet = h->Integral (h->FindBin (2.5), h->GetNbinsX ()) / totJetEvts;//h->Integral ();
    h->Rebin (2);
    h->Scale (1./totJetEvts, "width");
    h->SetLineColor (kRed);
    h->SetMarkerColor (kRed);
    h->SetMarkerStyle (kOpenCircle);
    h->DrawCopy ("e1 same");

    h = h_cluster_sumGaps_Pb_0n0n;
    const float evtFrac_0n0n = h->Integral () / totEvts;
    const float upcFrac_0n0n = h->Integral (h->FindBin (2.5), h->GetNbinsX ()) / totEvts;//h->Integral ();
    h->Scale (1./totEvts, "width");
    h->SetLineColor (kBlue);
    h->SetMarkerColor (kBlue);
    h->DrawCopy ("e1 same");

    h = h_cluster_sumGaps_Pb_0n0n_30GeVJet;
    const float evtFrac_0n0n_30GeVJet = h->Integral () / totEvts;
    const float upcFrac_0n0n_30GeVJet = h->Integral (h->FindBin (2.5), h->GetNbinsX ()) / totJetEvts;//h->Integral ();
    h->Rebin (2);
    h->Scale (1./totJetEvts, "width");
    h->SetLineColor (kBlue);
    h->SetMarkerColor (kBlue);
    h->SetMarkerStyle (kOpenCircle);
    h->DrawCopy ("e1 same");

    //h_cluster_sumGaps_Pb_0nXn->DrawCopy ("e1 same");
    //h_cluster_sumGaps_Pb_0n0n->DrawCopy ("e1 same");
    //h_cluster_sumGaps_Pb_0nXn_30GeVJet->DrawCopy ("e1 same");
    //h_cluster_sumGaps_Pb_0n0n_30GeVJet->DrawCopy ("e1 same");

    TLine* l = new TLine ();
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->SetLineColor (kBlack);
    l->DrawLine (2.5, 1e-3, 2.5, 0.5);

    myText (0.26, 0.885, kBlack, "#bf{#it{ATLAS}} Internal", 0.030);
    myText (0.26, 0.850, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.030);
    myText (0.26, 0.815, kBlack, "HLT_mb_sptrk_L1MBTS_1", 0.030);
    myText (0.26, 0.780, kBlack, "#Sigma#it{E}_{ZDC}^{Pb} < 18.0971 TeV (#geq 80%)", 0.030);
    myText (0.33, 0.740, kBlack, "Selection", 0.030);
    myText (0.65, 0.740, kBlack, "Frac. w/ #Sigma#Delta#eta > 2.5", 0.030);
    myLineText2 (0.33, 0.710, kRed, kFullCircle, "\"0nXn\"", 1.2, 0.030);
    myLineText2 (0.33, 0.675, kBlue, kFullCircle, "\"0n0n\"", 1.2, 0.030);
    myLineText2 (0.33, 0.640, kRed, kOpenCircle, "\"0nXn\" + #it{p}_{T}^{lead} > 30 GeV", 1.2, 0.030);
    myLineText2 (0.33, 0.605, kBlue, kOpenCircle, "\"0n0n\" + #it{p}_{T}^{lead} > 30 GeV", 1.2, 0.030);
    myText (0.65, 0.705, kBlack, Form ("%.2f%%", upcFrac_0nXn*100), 0.030);
    myText (0.65, 0.670, kBlack, Form ("%.2f%%", upcFrac_0n0n*100), 0.030);
    myText (0.65, 0.635, kBlack, Form ("%.2f%%", upcFrac_0nXn_30GeVJet*100), 0.030);
    myText (0.65, 0.600, kBlack, Form ("%.2f%%", upcFrac_0n0n_30GeVJet*100), 0.030);

    c->SaveAs ("Plots/Peripheral_Cluster_SumGaps.pdf");
  }

  inFile->Close ();

}

#endif
