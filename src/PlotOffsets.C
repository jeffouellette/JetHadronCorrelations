#ifndef __JetHadronCorrelatorPlotOffsets_C__
#define __JetHadronCorrelatorPlotOffsets_C__

#include <iostream>
#include <math.h>

#include <TColor.h>
#include <TLine.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

#include <ArrayTemplates.h>
#include <Utilities.h>
#include <MyStyle.h>
#include <MyColors.h>

#include "CentralityDefs.h"
#include "Params.h"
#include "OutTree.h"
#include "TreeVariables.h"
#include "LocalUtilities.h"
#include "Variations.h"

using namespace JetHadronCorrelations;

typedef TGraphErrors TGE;


TLine* l = new TLine ();
TLatex* tl = new TLatex ();


void PlotOffsets (const char* tag) {

  ifstream offsetsFile;
  offsetsFile.open (Form ("%s/aux/%s_IAAOffsets.dat", workPath.Data (), tag));

  string inTag, pTChSel, centStr, offpercerr;
  double offset = 0, offerr = 0;

  TGE** g_offsets = new TGE*[numZdcCentBins];
  for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
    g_offsets[iCent] = new TGE ();
  }

  while (!offsetsFile.eof ()) {
    offsetsFile >> inTag >> pTChSel >> centStr >> offset >> offerr >> offpercerr;

    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {

      if (TString (centStr.c_str ()) == TString (Form ("%i-%i%%", zdcCentPercs[iCent+1], zdcCentPercs[iCent]))) {
        g_offsets[iCent]->SetPoint (g_offsets[iCent]->GetN (), 0.5 * (pTChStrCuts[pTChSel].first + pTChStrCuts[pTChSel].second), -offset);
        g_offsets[iCent]->SetPointError (g_offsets[iCent]->GetN () -1, 0.5 * fabs (pTChStrCuts[pTChSel].first - pTChStrCuts[pTChSel].second), offerr);
        continue;
      }
    }
  }


  {
    TCanvas* c = new TCanvas ("c", "", 800, 800);
    gPad->SetLogx ();

    for (int iCent = 0; iCent < numZdcCentBins; iCent++) {
      g_offsets[iCent]->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      g_offsets[iCent]->GetYaxis ()->SetTitle ("#delta#it{I}_{pPb}");

      g_offsets[iCent]->GetYaxis ()->SetRangeUser (-0.5, 0.5);

      g_offsets[iCent]->SetLineColor (colors[iCent]);
      g_offsets[iCent]->SetLineWidth (2);
      g_offsets[iCent]->SetMarkerColor (colors[iCent]);

      g_offsets[iCent]->Draw (iCent == 0 ? "AP" : "P");
    }

    myText (0.60, 0.900, kBlack, "#bf{#it{ATLAS}} Internal", 0.032);
    myText (0.60, 0.860, kBlack, "#it{p}+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.032);
    myText (0.60, 0.820, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.032);
    myText (0.60, 0.780, kBlack, "#it{p}_{T}^{jet} > 60 GeV", 0.032);
    myText (0.60, 0.740, kBlack, "#Sigma#it{E}_{ZDC}^{Pb} Percentiles", 0.032);

    for (int iCent = 0; iCent < numZdcCentBins; iCent++)
      myText (0.6, 0.42-iCent*0.04, colors[iCent], Form ("#bf{%i-%i%%}", zdcCentPercs[iCent+1], zdcCentPercs[iCent]), 0.032);

    myText (0.2, 0.88, kBlack, "#bf{#it{I}_{pPb} #it{increases}}", 0.032);
    myText (0.2, 0.22, kBlack, "#bf{#it{I}_{pPb} #it{decreases}}", 0.032);

    c->SaveAs ("Plots/offsets.pdf");
  }

}


#endif
