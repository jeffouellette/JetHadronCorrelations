#ifndef __JetHadronCorrelatorPlotJets_C__
#define __JetHadronCorrelatorPlotJets_C__

#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"

#include <ArrayTemplates.h>
#include <Utilities.h>
#include <MyStyle.h>
#include <MyColors.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include <iostream>
#include <math.h>

using namespace JetHadronCorrelations;


TLine* l = new TLine ();
TLatex* tl = new TLatex ();


void PlotJets (const char* inFileTag) {

  TFile* inFile = nullptr;

  TH1D** h_evt_counts = new TH1D*[2];
  TH1D** h_jet_counts = new TH1D*[2];
  //TH1D** h_ljet_counts = new TH1D*[2];
  //TH1D** h_sljet_counts = new TH1D*[2]; 
  TH1D** h_jet_pt = new TH1D*[2];
  TH2D** h2_jet_pt_cov = new TH2D*[2];

  TH1D* h_jet_pt_ratio = nullptr;

  TH2D** h2_jet_eta_phi = new TH2D*[2];

  TGAE** g_jet_pt_syst = new TGAE*[2];

  TGAE* g_jet_pt_ratio_syst = nullptr;


  TH1D*** h_jet_pt_syst = Get2DArray <TH1D*> (2, nVar);

  TH1D** h_jet_pt_ratio_syst = Get1DArray <TH1D*> (nVar);


  {
    TString inFileName = inFileTag;
    inFileName.ReplaceAll (".root", "");
    inFileName = Form ("%s/Results/PlotJets_%s.root", rootPath.Data (), inFileName.Data ());
    std::cout << "Reading " << inFileName.Data () << std::endl;
    inFile = new TFile (inFileName, "read");

    h_evt_counts[0] = (TH1D*) inFile->Get ("h_evt_counts_ref");
    h_jet_counts[0] = (TH1D*) inFile->Get ("h_jet_counts_ref");
    //h_ljet_counts[0] = (TH1D*) inFile->Get ("h_ljet_counts_ref");
    //h_sljet_counts[0] = (TH1D*) inFile->Get ("h_sljet_counts_ref");
    h_evt_counts[1] = (TH1D*) inFile->Get ("h_evt_counts");
    h_jet_counts[1] = (TH1D*) inFile->Get ("h_jet_counts");
    //h_ljet_counts[1] = (TH1D*) inFile->Get ("h_ljet_counts");
    //h_sljet_counts[1] = (TH1D*) inFile->Get ("h_sljet_counts");

    h_jet_pt[0] = (TH1D*) inFile->Get ("h_jet_pt_ref");
    h_jet_pt[1] = (TH1D*) inFile->Get ("h_jet_pt");

    h_jet_pt_ratio = (TH1D*) inFile->Get ("h_jet_pt_ratio");

    g_jet_pt_syst[0] = (TGAE*) inFile->Get ("g_jet_pt_ref_syst");
    g_jet_pt_syst[1] = (TGAE*) inFile->Get ("g_jet_pt_syst");

    g_jet_pt_ratio_syst = (TGAE*) inFile->Get ("g_jet_pt_ratio_syst");


    for (int iVar = 1; iVar < nVar; iVar++) {
      h_jet_pt_syst[0][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_ref_%s", variations[iVar].Data ()));
      h_jet_pt_syst[1][iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_%s", variations[iVar].Data ()));
  
      h_jet_pt_ratio_syst[iVar] = (TH1D*) inFile->Get (Form ("h_jet_pt_ratio_%s", variations[iVar].Data ()));
    }

    h2_jet_eta_phi[0] = (TH2D*) inFile->Get ("h2_jet_eta_phi_ref");
    h2_jet_eta_phi[1] = (TH2D*) inFile->Get ("h2_jet_eta_phi");
  }



  l->SetLineStyle (7);
  l->SetLineWidth (1);
  l->SetLineColor (kBlack);



  {
    const char* canvasName = "c_jet_pt";
    TCanvas* c = new TCanvas (canvasName, "", 800, 920);
    c->cd ();
    const double fPad = 600./920.;
    TPad* uPad = new TPad (Form ("%s_uPad", canvasName), "", 0, 1-fPad, 1, 1.0);
    TPad* dPad = new TPad (Form ("%s_dPad", canvasName), "", 0, 0.0, 1, 1-fPad);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);

    uPad->SetLeftMargin (0.12);
    dPad->SetLeftMargin (0.12);
    uPad->SetRightMargin (0.03);
    dPad->SetRightMargin (0.03);
    uPad->Draw ();
    dPad->Draw ();

    TH1D* h = nullptr;
    TGAE* g = nullptr;

    double ymin=1e-7;
    double ymax=1e0;

    uPad->cd ();
    uPad->SetLogy ();

    h = (TH1D*) h_jet_pt[0]->Clone ("h");
    h->Reset ();
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{jet} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/fPad);
    h->GetXaxis ()->SetLabelSize (0.028/fPad);
    h->GetXaxis ()->SetTitleOffset (2.1*fPad);
    h->GetYaxis ()->SetTitle ("(1 / N_{jet}) (dN_{jet} / d#it{p}_{T}^{jet})");
    h->GetYaxis ()->SetTitleSize (0.028/fPad);
    h->GetYaxis ()->SetLabelSize (0.028/fPad);
    h->GetYaxis ()->SetTitleOffset (1.8*fPad);

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_jet_pt_syst[0];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_pt[0];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myBlue);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    g = g_jet_pt_syst[1];
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myRed);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_pt[1];
    g = make_graph (h);
    ResetXErrors (g);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    g->SetMarkerColor (myRed);
    g->SetLineColor (myRed);
    g->SetLineWidth (2);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    myText (0.50, 0.78, kBlack, "#bf{#it{ATLAS}} Internal", 0.028/fPad);
    myMarkerText (0.50, 0.72, myRed, kFullCircle, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV, #bf{0-20%}", 0.8, 0.028/fPad, true);
    myMarkerText (0.50, 0.66, myBlue, kFullCircle, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.8, 0.028/fPad, true);
    //myText (0.50, 0.60, kBlack, strcmp (tag, "30GeVJets") == 0 ? "Minimum bias trigger" : "J50 Trigger", 0.028/fPad);


    dPad->cd ();

    h = (TH1D*) h_jet_pt_ratio->Clone ("h");

    //float avgy = 0;
    //for (int i = 1; i <= h->GetNbinsX (); i++) avgy += h->GetBinContent (i) / h->GetNbinsX ();
    //ymin = 0.3*avgy;
    //ymax = 3*avgy;
    ymin=0.7;
    ymax=1.3;

    h->Reset ();
    for (int i = 1; i <= h->GetNbinsX (); i++) h->SetBinContent (i, 1);
    h->GetXaxis ()->SetMoreLogLabels ();
    h->GetYaxis ()->SetRangeUser (ymin, ymax);
    h->GetXaxis ()->SetTitle ("#it{p}_{T}^{jet} [GeV]");
    h->GetXaxis ()->SetTitleSize (0.028/(1-fPad));
    h->GetXaxis ()->SetLabelSize (0.028/(1-fPad));
    h->GetXaxis ()->SetTitleOffset (3.2*(1-fPad));
    h->GetYaxis ()->SetTitle ("#it{p}+Pb / #it{pp}");
    h->GetYaxis ()->SetTitleSize (0.028/(1-fPad));
    h->GetYaxis ()->SetLabelSize (0.028/(1-fPad));
    h->GetYaxis ()->SetTitleOffset (1.8*(1-fPad));
    h->GetYaxis ()->CenterTitle ();

    h->SetLineWidth (1);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    g = g_jet_pt_ratio_syst;
    g->SetMarkerStyle (0);
    g->SetMarkerSize (0.);
    g->SetLineColor (myBlue);
    g->SetLineWidth (1);
    ((TGAE*) g->Clone ())->Draw ("5");
    SaferDelete (&g);
    h = h_jet_pt_ratio;
    g = make_graph (h);
    ResetXErrors (g);
    g->SetLineColor (myBlue);
    g->SetLineWidth (2);
    g->SetMarkerColor (myBlue);
    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (0.8);
    ((TGAE*) g->Clone ())->Draw ("p");
    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerColor (kBlack);
    g->SetLineWidth (0);
    ((TGAE*) g->Clone ())->Draw ("p");
    SaferDelete (&g);

    c->SaveAs (Form ("%s/Plots/JetPtSpectrum.pdf", workPath.Data ()));
  }



  {
    TCanvas* c = new TCanvas ("c2", "", 800, 800);
    c->SetRightMargin (0.18);
    c->SetLeftMargin (0.15);
    c->SetTopMargin (0.04);
    c->SetBottomMargin (0.15);

    TH2D* h2 = h2_jet_eta_phi[1];

    h2->Draw ("colz");

    c->SaveAs (Form ("%s/Plots/JetEtaPhiSpectrum.pdf", workPath.Data ()));
  }


}


#endif
