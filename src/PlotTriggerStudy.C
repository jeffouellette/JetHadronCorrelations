#ifndef __JetHadronCorrelator_PlotTriggerStudy_C__
#define __JetHadronCorrelator_PlotTriggerStudy_C__

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
#include <TLine.h>

#include <iostream>
#include <string.h>
#include <math.h>

using namespace JetHadronCorrelations;


void PlotTriggerStudy () {


  TString inFileName = Form ("%s/TriggerStudy.root", rootPath.Data ());
  std::cout << "Reading " << inFileName.Data () << std::endl;
  TFile* inFile = new TFile (inFileName.Data (), "read");

  TEfficiency* e_pPb_j50_eff = (TEfficiency*) inFile->Get ("e_pPb_j50_eff");


  {
    TCanvas* c = new TCanvas ("c", "", 800, 600);

    TH1D* h = new TH1D ("h", ";Leading jet #it{p}_{T}^{reco} [GeV];#varepsilon_{trigger}", 1, 20, 100);
    h->SetBinContent (1, 1);
    h->GetYaxis ()->SetRangeUser (0, 1.2);
    h->SetLineWidth (2);
    h->SetLineStyle (2);
    h->DrawCopy ("hist ][");
    SaferDelete (&h);

    TEfficiency* e = e_pPb_j50_eff;

    e->SetMarkerStyle (kOpenCircle);
    e->SetMarkerColor (kGreen+3);

    e->Draw ("e1 same");

    int bin = 1;
    while (bin < e->GetPassedHistogram ()->GetNbinsX () && e->GetEfficiency (bin) < 0.99) bin++;
    const float x = e->GetPassedHistogram ()->GetBinLowEdge (bin);

    TLine* l = new TLine ();
    l->SetLineWidth (2);
    l->SetLineStyle (3);
    l->SetLineColor (kGray+2);
    l->DrawLine (x, 0, x, 1);

    myText (0.56, 0.500, kBlack, "#bf{#it{ATLAS}} Internal", 0.040);
    myText (0.56, 0.450, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.040);
    myText (0.56, 0.400, kBlack, "anti-#it{k}_{T} R=0.4 HI jets", 0.040);
    myText (0.56, 0.350, kBlack, "Ref.: HLT_mb_sptrk_L1MBTS_1", 0.040);
    myLineText2 (0.60, 0.300, kGreen+3, kOpenCircle, "HLT_j50_ion_L1J10", 1.2, 0.040);
    
    c->SaveAs ("Plots/TriggerEfficiency.pdf");
  }

}

#endif
