
#include "include/TrackMomentumFit.h" 
#include "src/AnalyzeTrackMomentumResolution.C"

TF1* f = nullptr, *fdraw = nullptr;
TH1D* h = nullptr;
TFile* inFile = nullptr;
double stddev, inte;

TrackMomentumFit tmf;

void testFit () {

  ltmf.DoDoubleGaussian (false);
  tmf.DoDoubleGaussian (false);

  //inFile = new TFile ("rootFiles/TrackMomentumResolution/Nominal/summary.root", "read");
  inFile = new TFile ("rootFiles/TrackMomentumResolution/Nominal/allSamples.root", "read");
  
  h = (TH1D*) inFile->Get ("h_tmr_pp_iPtch4_iEta20");
  //h = (TH1D*) inFile->Get ("h_tmr_integratedEta_pp_iPtch0_iEta1");

  f = new TF1 ("f", &ltmf, -1.0, 4.0, ltmf.ndf ()); 

  DoLogFit (h, f);
  ////h->Rebin (2);


  //TGraphErrors* g = new TGraphErrors ();
  //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
  //  if (h->GetBinContent (iX) > 0) {
  //    g->SetPoint (g->GetN (), h->GetBinCenter (iX), std::log (h->GetBinContent (iX)));
  //    g->SetPointError (g->GetN () - 1, 0.5 * h->GetBinWidth (iX), std::fabs (h->GetBinError (iX) / h->GetBinContent (iX)));
  //    //h->SetBinError (iX, h->GetBinError (iX) / h->GetBinContent (iX));
  //    //h->SetBinContent (iX, std::log (h->GetBinContent (iX)));
  //  }
  //}
  //
  //stddev = h->GetStdDev ();
  //inte = h->Integral (h->FindBin (-stddev), h->FindBin (stddev));
  //f = new TF1 ("f", &ltmf, -1.0, 4.0, 9); 
  fdraw = new TF1 ("fdraw", &tmf, -1.0, 4.0, tmf.ndf ());
  
  //g->Draw ("AP");
  h->Draw ("e1");
  gPad->SetLogy();

  tmf.CopyParams (f, fdraw);

  fdraw->SetLineColor (kBlue+1);
  fdraw->Draw ("L same");

}
