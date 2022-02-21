#ifndef __AnalyzeJetEnergyResolution_C__
#define __AnalyzeJetEnergyResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TLatex.h>

#include <Utilities.h>
#include <MyColors.h>
#include <ArrayTemplates.h>

#include <AtlasStyle.h>

#include "Params.h"
#include "LocalUtilities.h"

using namespace JetHadronCorrelations;

typedef TGraphAsymmErrors TGAE;

// Number of sigma to use for second fit
const float nSigFit1 = 4.0;
const float nSigFit2 = 2.0;
const float nSigFit3 = 1.4;


// Minimum number of histogram entries to attempt a fit (400 is an arbitrary choice)
const int minEntriesForFit = 200;


/**
 * Wrapper function for performing a fit to a Gaussian TF1.
 */
TF1* DoGaussianFit (TH1D* h, double min = 0.0, double max = 1.4) {
  //std::cout << "Fitting " << h->GetName () << " to Gaussian." << std::endl;

  double mean = h->GetMean ();
  double sigma = h->GetStdDev ();
  //double amp = h->Integral (h->FindBin (mean-sigma), h->FindBin (mean+sigma)) / 0.68;
  double amp = h->Integral ();

  // Check we have "enough" entries to perform a fit.
  if (h->Integral (h->FindBin (mean-nSigFit1*sigma), h->FindBin (mean+nSigFit1*sigma)) < minEntriesForFit) {
    std::cout << "Not enough counts for fit! Will skip " << h->GetName () << std::endl;
    return nullptr;
  }
  //else
  //  std::cout << "Fitting " << h->GetName () << " to gaussian" << std::endl;

  TF1* fit = new TF1 ("fit", "gaus(0)", min, max);//std::fmax (mean-nSigFit1*sigma, 0.0), std::fmin (mean+nSigFit1*sigma, 2.4));

  fit->SetParameter (0, amp);
  fit->SetParameter (1, mean);
  fit->SetParameter (2, sigma);

  h->Fit (fit, "RN0Q");

  amp = fit->GetParameter (0);
  mean = fit->GetParameter (1);
  sigma = fit->GetParameter (2);

  delete fit;

  // Check we have "enough" entries to perform a fit.
  if (h->Integral (h->FindBin (mean-nSigFit2*sigma), h->FindBin (mean+nSigFit2*sigma)) < minEntriesForFit) {
    std::cout << "Not enough counts for fit! Will skip " << h->GetName () << std::endl;
    return nullptr;
  }
  //else
  //  std::cout << "Fitting " << h->GetName () << " to gaussian" << std::endl;

  fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-nSigFit2*sigma, min), std::fmin (mean+nSigFit2*sigma, max));

  fit->SetParameter (0, amp);
  fit->SetParameter (1, mean);
  fit->SetParameter (2, sigma);

  h->Fit (fit, "RN0Q");

  amp = fit->GetParameter (0);
  mean = fit->GetParameter (1);
  sigma = fit->GetParameter (2);

  delete fit;

  // Check we have "enough" entries to perform a fit.
  if (h->Integral (h->FindBin (mean-nSigFit3*sigma), h->FindBin (mean+nSigFit3*sigma)) < minEntriesForFit) {
    std::cout << "Not enough counts for fit! Will skip " << h->GetName () << std::endl;
    return nullptr;
  }
  //else
  //  std::cout << "Fitting " << h->GetName () << " to gaussian" << std::endl;

  fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-nSigFit3*sigma, min), std::fmin (mean+nSigFit3*sigma, max));

  fit->SetParameter (0, amp);
  fit->SetParameter (1, mean);
  fit->SetParameter (2, sigma);

  h->Fit (fit, "RN0Q");

  return fit;
}



void AnalyzeJetEnergyResolution () {

  SetAtlasStyle ();

  const int nFinerEtaBins = 56;
  const double* finerEtaBins = linspace (-2.8, 2.8, nFinerEtaBins);
  
  const double etaBins[] = {-2.8, -2.1, -1.2, -0.8, -0.3, 0, 0.3, 0.8, 1.2, 2.1, 2.8};
  const int nEtaBins = sizeof (etaBins) / sizeof (etaBins[0]) - 1;
  
  //const double pTJBins[] = {10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 360, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300};
  const double pTJBins[] = {8, 10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 190, 220, 250, 280, 310, 340, 370, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300};
  const int nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;

  const int nRespBins = 240;
  const double* respBins = linspace (0, 2.4, nRespBins);

  const int nEtaRespBins = 80;
  const double* etaRespBins = linspace (-0.2, 0.2, nEtaRespBins);
  
  const std::vector <int> systems = {0, 1}; // 0 = pp, 1 = pPb
  //const std::vector <int> systems = {0}; // 0 = pp, 1 = pPb // for testing
  const int nSystems = systems.size ();
  
  //const std::vector <JetRadius> radii = {JetRadius::R0p2, JetRadius::R0p4};
  const std::vector <JetRadius> radii = {JetRadius::R0p4}; // for testing
  const int nRadii = radii.size ();


  TFile* inFile = new TFile (Form ("%s/JetEnergyResolution/Nominal/allSamples.root", rootPath.Data ()), "read");

  TH1D***** h_jpts     = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nFinerEtaBins);
  TH1D***** h_jes      = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nFinerEtaBins);
  TH1D***** h_jetacorr = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nFinerEtaBins);

  TH2D*** h2_jet_ptreco_pttruth = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_jet_ereco_etruth   = Get2DArray <TH2D*> (nSystems, nRadii);

  for (int iSys : systems) {

    const std::string sys = (iSys == 0 ? "pp" : "pPb");

    for (int iR = 0; iR < nRadii; iR++) {

      const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

      h2_jet_ptreco_pttruth[iSys][iR]  = (TH2D*) inFile->Get (Form ("h2_r%i_jet_ptreco_pttruth_%s", r, sys.c_str ()));
      h2_jet_ereco_etruth[iSys][iR]    = (TH2D*) inFile->Get (Form ("h2_r%i_jet_ereco_etruth_%s", r, sys.c_str ()));

      for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        for (int iEta = 0; iEta < nFinerEtaBins; iEta++) {

          h_jpts[iSys][iR][iPtJ][iEta]      = (TH1D*) inFile->Get (Form ("h_r%i_jpts_%s_iPtJ%i_iEta%i", r, sys.c_str (), iPtJ, iEta));
          h_jes[iSys][iR][iPtJ][iEta]       = (TH1D*) inFile->Get (Form ("h_r%i_jes_%s_iPtJ%i_iEta%i", r, sys.c_str (), iPtJ, iEta));
          h_jetacorr[iSys][iR][iPtJ][iEta]  = (TH1D*) inFile->Get (Form ("h_r%i_jetacorr_%s_iPtJ%i_iEta%i", r, sys.c_str (), iPtJ, iEta));

        } // end loop over iEta

      } // end loop over iPtJ

    } // end loop over iR

  } // end loop over iSys


  TFile* outFile = new TFile (Form ("%s/JetEnergyResolution/Nominal/summary.root", rootPath.Data ()), "recreate");

  TH1D***** h_jes_integratedEta       = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nEtaBins);
  TH1D***** h_jpts_integratedEta      = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nEtaBins);
  TH1D***** h_jetacorr_integratedEta  = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nEtaBins);

  TH2D**** h2_jes_integratedEta       = Get3DArray <TH2D*> (nSystems, nRadii, nEtaBins);
  TH2D**** h2_jpts_integratedEta      = Get3DArray <TH2D*> (nSystems, nRadii, nEtaBins);
  TH2D**** h2_jetacorr_integratedEta  = Get3DArray <TH2D*> (nSystems, nRadii, nEtaBins);

  TH2D*** h2_avg_jpts = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_avg_jptr = Get2DArray <TH2D*> (nSystems, nRadii);
  TH1D**** h_avg_jpts = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);
  TH1D**** h_avg_jptr = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);

  TF1**** f_avg_jptr = Get3DArray <TF1*> (nSystems, nRadii, nEtaBins+1);

  TH2D*** h2_avg_jes  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_avg_jer  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH1D**** h_avg_jes  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);
  TH1D**** h_avg_jer  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);

  TF1**** f_avg_jer = Get3DArray <TF1*> (nSystems, nRadii, nEtaBins+1);

  TH2D*** h2_avg_jetacorr = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_avg_jetares  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH1D**** h_avg_jetacorr = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);
  TH1D**** h_avg_jetares  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);

  for (int iSys : systems) {

    const std::string sys = (iSys == 0 ? "pp" : "pPb");

    for (int iR = 0; iR < nRadii; iR++) {

      const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

      h2_avg_jpts[iSys][iR]     = new TH2D (Form ("h2_r%i_avg_jpts_%s", r, sys.c_str ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);
      h2_avg_jptr[iSys][iR]     = new TH2D (Form ("h2_r%i_avg_jptr_%s", r, sys.c_str ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);

      h2_avg_jes[iSys][iR]      = new TH2D (Form ("h2_r%i_avg_jes_%s", r, sys.c_str ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);
      h2_avg_jer[iSys][iR]      = new TH2D (Form ("h2_r%i_avg_jer_%s", r, sys.c_str ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);

      h2_avg_jetacorr[iSys][iR] = new TH2D (Form ("h2_r%i_avg_jetacorr_%s", r, sys.c_str ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};#LT#eta_{reco} - #eta_{truth}#GT", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);
      h2_avg_jetares[iSys][iR]  = new TH2D (Form ("h2_r%i_avg_jetares_%s", r, sys.c_str ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};#sigma#left[#eta_{reco} - #eta_{truth}#GT#right]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);


      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        h_avg_jpts[iSys][iR][iEta]      = new TH1D (Form ("h_r%i_avg_jpts_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]", nPtJBins, pTJBins);
        h_avg_jptr[iSys][iR][iEta]      = new TH1D (Form ("h_r%i_avg_jptr_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]", nPtJBins, pTJBins);

        h_avg_jes[iSys][iR][iEta]       = new TH1D (Form ("h_r%i_avg_jes_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]", nPtJBins, pTJBins);
        h_avg_jer[iSys][iR][iEta]       = new TH1D (Form ("h_r%i_avg_jer_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]", nPtJBins, pTJBins);

        h_avg_jetacorr[iSys][iR][iEta]  = new TH1D (Form ("h_r%i_avg_jetacorr_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#LT#eta_{reco} - #eta_{truth}#GT", nPtJBins, pTJBins);
        h_avg_jetares[iSys][iR][iEta]   = new TH1D (Form ("h_r%i_avg_jetares_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#sigma#left[#eta_{reco} - #eta_{truth}#GT#right]", nPtJBins, pTJBins);


        h2_jpts_integratedEta[iSys][iR][iEta] = new TH2D (Form ("h2_r%i_jpts_integratedEta_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", nPtJBins, pTJBins, nRespBins, respBins);
        h2_jpts_integratedEta[iSys][iR][iEta]->Sumw2 ();

        h2_jes_integratedEta[iSys][iR][iEta] = new TH2D (Form ("h2_r%i_jes_integratedEta_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#it{E}_{reco} / #it{E}_{truth};Counts", nPtJBins, pTJBins, nRespBins, respBins);
        h2_jes_integratedEta[iSys][iR][iEta]->Sumw2 ();

        h2_jetacorr_integratedEta[iSys][iR][iEta] = new TH2D (Form ("h2_r%i_jetacorr_integratedEta_%s_iEta%i", r, sys.c_str (), iEta), ";#it{p}_{T}^{truth} [GeV];#eta_{reco} - #eta_{truth};Counts", nPtJBins, pTJBins, nEtaRespBins, etaRespBins);
        h2_jetacorr_integratedEta[iSys][iR][iEta]->Sumw2 ();


        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          h_jpts_integratedEta[iSys][iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jpts_integratedEta_%s_iPtJ%i_iEta%i", r, sys.c_str (), iPtJ, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", nRespBins, respBins);
          h_jpts_integratedEta[iSys][iR][iPtJ][iEta]->Sumw2 ();

          h_jes_integratedEta[iSys][iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jes_integratedEta_%s_iPtJ%i_iEta%i", r, sys.c_str (), iPtJ, iEta), ";#it{E}_{reco} / #it{E}_{truth};Counts", nRespBins, respBins);
          h_jes_integratedEta[iSys][iR][iPtJ][iEta]->Sumw2 ();

          h_jetacorr_integratedEta[iSys][iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jetacorr_integratedEta_%s_iPtJ%i_iEta%i", r, sys.c_str (), iPtJ, iEta), ";#eta_{reco} - #eta_{truth};Counts", nEtaRespBins, etaRespBins);
          h_jetacorr_integratedEta[iSys][iR][iPtJ][iEta]->Sumw2 ();

        } // end loop over iPtJ

      } // end loop over iEta

    } // end loop over iR

  } // end loop over iSys


  TCanvas* c = new TCanvas ("c", "", 1200, 1200);
  c->Divide (4, 4);

  TCanvas* c2d = new TCanvas ("c2d", "", 880, 800);
  const double lMargin = 0.15;
  const double rMargin = 0.11;
  const double bMargin = 0.15;
  const double tMargin = 0.04;

  c2d->SetLeftMargin (lMargin);
  c2d->SetRightMargin (rMargin);
  c2d->SetBottomMargin (bMargin);
  c2d->SetTopMargin (tMargin);

  c2d->SetLogx();
  c2d->SetLogz();

  TLatex* tl = new TLatex ();

  int iCanvas = 1;

  for (int iSys : systems) {

    const std::string sys = (iSys == 0 ? "pp" : "pPb");

    for (int iR = 0; iR < nRadii; iR++) {

      const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

      for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          TH1D* h = h_jpts[iSys][iR][iPtJ][iFinerEta];

          // first add to eta-integrated plot
          const float binCenter = 0.5 * (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
          int iEta = 0;
          while (iEta < nEtaBins) {
            if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
              break;
            iEta++;
          }

          h_jpts_integratedEta[iSys][iR][iPtJ][iEta]->Add (h);
          h_jpts_integratedEta[iSys][iR][iPtJ][nEtaBins]->Add (h);

          //TF1* fit = DoGaussianFit (h);

          //if (fit) {
          //  const double mean = fit->GetParameter (1);
          //  const double sigma = fit->GetParameter (2);
          //  const double mean_err = fit->GetParError (1);
          //  const double sigma_err = fit->GetParError (2);

          //  delete fit;

          //  const double jpts = mean;
          //  const double jpts_err = mean_err;
          //  
          //  const double jptr = sigma / mean;
          //  const double jptr_err = std::fabs (jptr) * std::sqrt (std::pow (mean_err/mean, 2) + std::pow (sigma_err/sigma, 2));

          //  h2_avg_jpts[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jpts * 100);
          //  h2_avg_jpts[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jpts_err * 100);
          //  h2_avg_jptr[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jptr * 100);
          //  h2_avg_jptr[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jptr_err * 100);
          //}
          //else {
          //  h2_avg_jpts[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
          //  h2_avg_jpts[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          //  h2_avg_jptr[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
          //  h2_avg_jptr[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          //}

        } // end loop over iPtJ


        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          TH1D* h = h_jes[iSys][iR][iPtJ][iFinerEta];

          // first add to eta-integrated plot
          const float binCenter = 0.5 * (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
          int iEta = 0;
          while (iEta < nEtaBins) {
            if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
              break;
            iEta++;
          }

          h_jes_integratedEta[iSys][iR][iPtJ][iEta]->Add (h);
          h_jes_integratedEta[iSys][iR][iPtJ][nEtaBins]->Add (h);

          //TF1* fit = DoGaussianFit (h);

          //if (fit) {
          //  const double mean = fit->GetParameter (1);
          //  const double sigma = fit->GetParameter (2);
          //  const double mean_err = fit->GetParError (1);
          //  const double sigma_err = fit->GetParError (2);

          //  delete fit;

          //  const double jes = mean;
          //  const double jes_err = mean_err;
          //  
          //  const double jer = sigma / mean;
          //  const double jer_err = std::fabs (jer) * std::sqrt (std::pow (mean_err/mean, 2) + std::pow (sigma_err/sigma, 2));

          //  h2_avg_jes[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jes * 100);
          //  h2_avg_jes[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jes_err * 100);
          //  h2_avg_jer[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jer * 100);
          //  h2_avg_jer[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jer_err * 100);
          //}
          //else {
          //  h2_avg_jes[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
          //  h2_avg_jes[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          //  h2_avg_jer[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
          //  h2_avg_jer[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          //}

        } // end loop over iPtJ


        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          TH1D* h = h_jetacorr[iSys][iR][iPtJ][iFinerEta];

          // first add to eta-integrated plot
          const float binCenter = 0.5 * (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
          int iEta = 0;
          while (iEta < nEtaBins) {
            if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
              break;
            iEta++;
          }

          h_jetacorr_integratedEta[iSys][iR][iPtJ][iEta]->Add (h);
          h_jetacorr_integratedEta[iSys][iR][iPtJ][nEtaBins]->Add (h);

          //TF1* fit = DoGaussianFit (h, -0.2, 0.2);

          //if (fit) {
          //  const double mean = fit->GetParameter (1);
          //  const double sigma = fit->GetParameter (2);
          //  const double mean_err = fit->GetParError (1);
          //  const double sigma_err = fit->GetParError (2);

          //  delete fit;

          //  const double jetacorr = mean;
          //  const double jetacorr_err = mean_err;
          //  
          //  const double jetares = sigma;
          //  const double jetares_err = sigma_err;

          //  h2_avg_jetacorr[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jetacorr * 100);
          //  h2_avg_jetacorr[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jetacorr_err * 100);
          //  h2_avg_jetares[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jetares * 100);
          //  h2_avg_jetares[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jetares_err * 100);
          //}
          //else {
          //  h2_avg_jetacorr[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
          //  h2_avg_jetacorr[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          //  h2_avg_jetares[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
          //  h2_avg_jetares[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          //}

        } // end loop over iPtJ

      } // end loop over iFinerEta



      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        TH2D* h2 = h2_jpts_integratedEta[iSys][iR][iEta];
        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
          TH1D* h = h_jpts_integratedEta[iSys][iR][iPtJ][iEta];

          for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {

            h2->SetBinContent (iPtJ+1, iY, h2->GetBinContent (iPtJ+1, iY) + h->GetBinContent (iY));
            h2->SetBinError (iPtJ+1, iY, std::sqrt (std::pow (h2->GetBinError (iPtJ+1, iY), 2) + std::pow (h->GetBinError (iY), 2)));

          } // end loop over iY
  
        } // end loop over iPtJ


        h2 = h2_jes_integratedEta[iSys][iR][iEta];
        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
          TH1D* h = h_jes_integratedEta[iSys][iR][iPtJ][iEta];

          for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {

            h2->SetBinContent (iPtJ+1, iY, h2->GetBinContent (iPtJ+1, iY) + h->GetBinContent (iY));
            h2->SetBinError (iPtJ+1, iY, std::sqrt (std::pow (h2->GetBinError (iPtJ+1, iY), 2) + std::pow (h->GetBinError (iY), 2)));

          } // end loop over iY
  
        } // end loop over iPtJ


        h2 = h2_jetacorr_integratedEta[iSys][iR][iEta];
        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {
          TH1D* h = h_jetacorr_integratedEta[iSys][iR][iPtJ][iEta];

          for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {

            h2->SetBinContent (iPtJ+1, iY, h2->GetBinContent (iPtJ+1, iY) + h->GetBinContent (iY));
            h2->SetBinError (iPtJ+1, iY, std::sqrt (std::pow (h2->GetBinError (iPtJ+1, iY), 2) + std::pow (h->GetBinError (iY), 2)));

          } // end loop over iY
  
        } // end loop over iPtJ

      } // end loop over iEta



      const std::string plotFileName = "Plots/JetEnergyResolution/JetResponsePlots_R" + std::to_string (r) + "_" + sys + ".pdf";
      c->Print ((plotFileName + "[").c_str (), "pdf");

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        iCanvas = 1;
        c->Clear ("D");

        {
          c2d->cd ();
          TH2D* h2 = (TH2D*) h2_jpts_integratedEta[iSys][iR][iEta]->Clone ("htemp");

          TAxis* xax = h2->GetXaxis ();
          TAxis* yax = h2->GetYaxis ();
          TAxis* zax = h2->GetZaxis ();

          // normalize distributions in response (along y-axis)
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
            double integral = 0;
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              integral += h2->GetBinContent (iX, iY) * yax->GetBinWidth (iY);
            if (integral > 0) {
              for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
                h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / integral);
                h2->SetBinError (iX, iY, h2->GetBinError (iX, iY) / integral);
              } // end loop over iY
            }
          } // end loop over iX
          
          xax->SetTitleFont (43);
          xax->SetTitleSize (32);
          yax->SetTitleFont (43);
          yax->SetTitleSize (32);
          zax->SetTitleFont (43);
          zax->SetTitleSize (32);
          zax->SetTitleOffset (1.1*zax->GetTitleOffset ());
          xax->SetLabelFont (43);
          //xax->SetLabelSize (32);
          xax->SetLabelSize (0);
          yax->SetLabelFont (43);
          yax->SetLabelSize (32);
          zax->SetLabelFont (43);
          zax->SetLabelSize (24);
          
          zax->SetTitle ("");

          h2->DrawCopy ("colz");
          SaferDelete (&h2);

          tl->SetTextFont (43);
          tl->SetTextSize (32);
          tl->SetTextAlign (21);
  
          const double ymin = 0.0; 
          const double ymax = 2.4; 
          const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
          tl->DrawLatex (20, yoff, "20");
          tl->DrawLatex (30, yoff, "30");
          tl->DrawLatex (40, yoff, "40");
          tl->DrawLatex (60, yoff, "60");
          tl->DrawLatex (100, yoff, "100");
          tl->DrawLatex (200, yoff, "200");
          tl->DrawLatex (400, yoff, "400");
          tl->DrawLatex (800, yoff, "800");
  
          tl->SetTextAlign (12);
  
          tl->SetTextSize (26);
          tl->DrawLatexNDC (0.36, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
          if (iSys == 0)
            tl->DrawLatexNDC (0.36, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
          else if (iSys == 1)
            tl->DrawLatexNDC (0.36, 0.845, "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
          tl->DrawLatexNDC (0.36, 0.810, Form ("Anti-#it{k}_{T} HI Jets, R=0.%i", r));
          tl->DrawLatexNDC (0.36, 0.775, iEta == nEtaBins ? "|#it{#eta}| < 2.8" : Form ("%g < #it{#eta} < %g", etaBins[iEta], etaBins[iEta+1]));

          c2d->Print (plotFileName.c_str (), "pdf");
        }

        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          TH1D* h = h_jpts_integratedEta[iSys][iR][iPtJ][iEta];

          TF1* fit = DoGaussianFit (h);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

            c->cd (iCanvas);
            gPad->SetLogy ();
            h->GetXaxis ()->SetRangeUser (0.2, 1.8);
            h->GetXaxis ()->SetTitle (Form ("#it{p}_{T}^{reco} / #it{p}_{T}^{truth}, %g < #it{p}_{T}^{truth} < %g GeV", pTJBins[iPtJ], pTJBins[iPtJ+1]));
            h->SetLineColor (kBlack);
            h->SetMarkerStyle (kDot);
            h->DrawCopy ("e1");
            fit->SetLineColor (myBlue);
            fit->SetLineWidth (1);
            ((TF1*) fit->Clone ())->Draw ("same");
            myText (0.24, 0.8, kBlack, Form ("#chi^{2}/ndf = %.2f/%i", fit->GetChisquare (), fit->GetNDF ()), 0.06);
            if (iCanvas == 16 || iPtJ == nPtJBins-1) {
              c->Print (plotFileName.c_str (), "pdf");
              iCanvas = 1;
              c->Clear ("D");
            }
            else
              iCanvas++;

            delete fit;

            const double jpts = mean;
            const double jpts_err = mean_err;

            const double jptr = sigma / mean;
            const double jptr_err = std::fabs (jptr) * std::sqrt (std::pow (mean_err/mean, 2) + std::pow (sigma_err/sigma, 2));

            h_avg_jpts[iSys][iR][iEta]->SetBinContent (iPtJ+1, jpts * 100);
            h_avg_jpts[iSys][iR][iEta]->SetBinError (iPtJ+1, jpts_err * 100);
            h_avg_jptr[iSys][iR][iEta]->SetBinContent (iPtJ+1, jptr * 100);
            h_avg_jptr[iSys][iR][iEta]->SetBinError (iPtJ+1, jptr_err * 100);
          }
          else {
            h_avg_jpts[iSys][iR][iEta]->SetBinContent (iPtJ+1, -1);
            h_avg_jpts[iSys][iR][iEta]->SetBinError (iPtJ+1, 0);
            h_avg_jptr[iSys][iR][iEta]->SetBinContent (iPtJ+1, -1);
            h_avg_jptr[iSys][iR][iEta]->SetBinError (iPtJ+1, 0);
          }

        } // end loop over iPtJ


        iCanvas = 1;
        c->Clear ("D");

        {
          c2d->cd ();
          TH2D* h2 = (TH2D*) h2_jes_integratedEta[iSys][iR][iEta]->Clone ("htemp");

          TAxis* xax = h2->GetXaxis ();
          TAxis* yax = h2->GetYaxis ();
          TAxis* zax = h2->GetZaxis ();
  
          // normalize distributions in response (along y-axis)
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
            double integral = 0;
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              integral += h2->GetBinContent (iX, iY) * yax->GetBinWidth (iY);
            if (integral > 0) {
              for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
                h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / integral);
                h2->SetBinError (iX, iY, h2->GetBinError (iX, iY) / integral);
              } // end loop over iY
            }
          } // end loop over iX
          
          xax->SetTitleFont (43);
          xax->SetTitleSize (32);
          yax->SetTitleFont (43);
          yax->SetTitleSize (32);
          zax->SetTitleFont (43);
          zax->SetTitleSize (32);
          zax->SetTitleOffset (1.1*zax->GetTitleOffset ());
          xax->SetLabelFont (43);
          //xax->SetLabelSize (32);
          xax->SetLabelSize (0);
          yax->SetLabelFont (43);
          yax->SetLabelSize (32);
          zax->SetLabelFont (43);
          zax->SetLabelSize (24);
          
          zax->SetTitle ("");

          h2->DrawCopy ("colz");
          SaferDelete (&h2);

          tl->SetTextFont (43);
          tl->SetTextSize (32);
          tl->SetTextAlign (21);
  
          const double ymin = 0.0; 
          const double ymax = 2.4; 
          const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
          tl->DrawLatex (20, yoff, "20");
          tl->DrawLatex (30, yoff, "30");
          tl->DrawLatex (40, yoff, "40");
          tl->DrawLatex (60, yoff, "60");
          tl->DrawLatex (100, yoff, "100");
          tl->DrawLatex (200, yoff, "200");
          tl->DrawLatex (400, yoff, "400");
          tl->DrawLatex (800, yoff, "800");
  
          tl->SetTextAlign (12);
  
          tl->SetTextSize (26);
          tl->DrawLatexNDC (0.36, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
          if (iSys == 0)
            tl->DrawLatexNDC (0.36, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
          else if (iSys == 1)
            tl->DrawLatexNDC (0.36, 0.845, "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
          tl->DrawLatexNDC (0.36, 0.810, Form ("Anti-#it{k}_{T} HI Jets, R=0.%i", r));
          tl->DrawLatexNDC (0.36, 0.775, iEta == nEtaBins ? "|#it{#eta}| < 2.8" : Form ("%g < #it{#eta} < %g", etaBins[iEta], etaBins[iEta+1]));
  
          c2d->Print (plotFileName.c_str (), "pdf");
        }

        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          TH1D* h = h_jes_integratedEta[iSys][iR][iPtJ][iEta];

          TF1* fit = DoGaussianFit (h);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

            c->cd (iCanvas);
            gPad->SetLogy ();
            h->GetXaxis ()->SetRangeUser (0.2, 1.8);
            h->GetXaxis ()->SetTitle (Form ("#it{E}_{reco} / #it{E}_{truth}, %g < #it{p}_{T}^{truth} < %g GeV", pTJBins[iPtJ], pTJBins[iPtJ+1]));
            h->SetLineColor (kBlack);
            h->SetMarkerStyle (kDot);
            h->DrawCopy ("e1");
            fit->SetLineColor (myBlue);
            fit->SetLineWidth (1);
            ((TF1*) fit->Clone ())->Draw ("same");
            myText (0.24, 0.8, kBlack, Form ("#chi^{2}/ndf = %.2f/%i", fit->GetChisquare (), fit->GetNDF ()), 0.06);
            if (iCanvas == 16 || iPtJ == nPtJBins-1) {
              c->Print (plotFileName.c_str (), "pdf");
              iCanvas = 1;
              c->Clear ("D");
            }
            else
              iCanvas++;

            delete fit;

            const double jes = mean;
            const double jes_err = mean_err;

            const double jer = sigma / mean;
            const double jer_err = std::fabs (jer) * std::sqrt (std::pow (mean_err/mean, 2) + std::pow (sigma_err/sigma, 2));

            h_avg_jes[iSys][iR][iEta]->SetBinContent (iPtJ+1, jes * 100);
            h_avg_jes[iSys][iR][iEta]->SetBinError (iPtJ+1, jes_err * 100);
            h_avg_jer[iSys][iR][iEta]->SetBinContent (iPtJ+1, jer * 100);
            h_avg_jer[iSys][iR][iEta]->SetBinError (iPtJ+1, jer_err * 100);
          }
          else {
            h_avg_jes[iSys][iR][iEta]->SetBinContent (iPtJ+1, -1);
            h_avg_jes[iSys][iR][iEta]->SetBinError (iPtJ+1, 0);
            h_avg_jer[iSys][iR][iEta]->SetBinContent (iPtJ+1, -1);
            h_avg_jer[iSys][iR][iEta]->SetBinError (iPtJ+1, 0);
          }

        } // end loop over iPtJ


        iCanvas = 1;
        c->Clear ("D");

        {
          c2d->cd ();
          TH2D* h2 = (TH2D*) h2_jetacorr_integratedEta[iSys][iR][iEta]->Clone ("htemp");

          TAxis* xax = h2->GetXaxis ();
          TAxis* yax = h2->GetYaxis ();
          TAxis* zax = h2->GetZaxis ();
  
          // normalize distributions in response (along y-axis)
          for (int iX = 1; iX <= h2->GetNbinsX (); iX++) {
            double integral = 0;
            for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
              integral += h2->GetBinContent (iX, iY) * yax->GetBinWidth (iY);
            if (integral > 0) {
              for (int iY = 1; iY <= h2->GetNbinsY (); iY++) {
                h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) / integral);
                h2->SetBinError (iX, iY, h2->GetBinError (iX, iY) / integral);
              } // end loop over iY
            }
          } // end loop over iX

          xax->SetTitleFont (43);
          xax->SetTitleSize (32);
          yax->SetTitleFont (43);
          yax->SetTitleSize (32);
          zax->SetTitleFont (43);
          zax->SetTitleSize (32);
          zax->SetTitleOffset (1.1*zax->GetTitleOffset ());
          xax->SetLabelFont (43);
          //xax->SetLabelSize (32);
          xax->SetLabelSize (0);
          yax->SetLabelFont (43);
          yax->SetLabelSize (32);
          zax->SetLabelFont (43);
          zax->SetLabelSize (24);
          
          zax->SetTitle ("");

          h2->DrawCopy ("colz");
          SaferDelete (&h2);

          tl->SetTextFont (43);
          tl->SetTextSize (32);
          tl->SetTextAlign (21);
 
          const double ymin = -0.2; 
          const double ymax = 0.2; 
          const double yoff = ymin - 0.04 * (ymax-ymin) / (1.-tMargin-bMargin);
          tl->DrawLatex (20, yoff, "20");
          tl->DrawLatex (30, yoff, "30");
          tl->DrawLatex (40, yoff, "40");
          tl->DrawLatex (60, yoff, "60");
          tl->DrawLatex (100, yoff, "100");
          tl->DrawLatex (200, yoff, "200");
          tl->DrawLatex (400, yoff, "400");
          tl->DrawLatex (800, yoff, "800");
  
          tl->SetTextAlign (12);
  
          tl->SetTextSize (26);
          tl->DrawLatexNDC (0.36, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
          if (iSys == 0)
            tl->DrawLatexNDC (0.36, 0.850, "Pythia8 #it{pp}, #sqrt{s} = 5.02 TeV");
          else if (iSys == 1)
            tl->DrawLatexNDC (0.36, 0.845, "Pythia8 + #it{p}+Pb Overlay, #sqrt{s_{NN}} = 5.02 TeV");
          tl->DrawLatexNDC (0.36, 0.810, Form ("Anti-#it{k}_{T} HI Jets, R=0.%i", r));
          tl->DrawLatexNDC (0.36, 0.775, iEta == nEtaBins ? "|#it{#eta}| < 2.8" : Form ("%g < #it{#eta} < %g", etaBins[iEta], etaBins[iEta+1]));
  
          c2d->Print (plotFileName.c_str (), "pdf");
        }

        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          TH1D* h = h_jetacorr_integratedEta[iSys][iR][iPtJ][iEta];

          TF1* fit = DoGaussianFit (h, -0.2, 0.2);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

            c->cd (iCanvas);
            h->GetXaxis ()->SetRangeUser (-0.05, 0.05);
            h->GetXaxis ()->SetTitle (Form ("#eta_{reco} - #eta_{truth}, %g < #it{p}_{T}^{truth} < %g GeV", pTJBins[iPtJ], pTJBins[iPtJ+1]));
            h->SetLineColor (kBlack);
            h->SetMarkerStyle (kDot);
            h->DrawCopy ("e1");
            fit->SetLineColor (myBlue);
            fit->SetLineWidth (1);
            ((TF1*) fit->Clone ())->Draw ("same");
            myText (0.24, 0.8, kBlack, Form ("#chi^{2}/ndf = %.2f/%i", fit->GetChisquare (), fit->GetNDF ()), 0.06);
            if (iCanvas == 16 || iPtJ == nPtJBins-1) {
              c->Print (plotFileName.c_str (), "pdf");
              iCanvas = 1;
              c->Clear ("D");
            }
            else
              iCanvas++;

            delete fit;
  
            const double jetacorr = mean;
            const double jetacorr_err = mean_err;
  
            const double jetares = sigma;
            const double jetares_err = sigma_err;
  
            h_avg_jetacorr[iSys][iR][iEta]->SetBinContent (iPtJ+1, jetacorr);
            h_avg_jetacorr[iSys][iR][iEta]->SetBinError (iPtJ+1, jetacorr_err);
            h_avg_jetares[iSys][iR][iEta]->SetBinContent (iPtJ+1, jetares);
            h_avg_jetares[iSys][iR][iEta]->SetBinError (iPtJ+1, jetares_err);
          }
          else {
            h_avg_jetacorr[iSys][iR][iEta]->SetBinContent (iPtJ+1, -1);
            h_avg_jetacorr[iSys][iR][iEta]->SetBinError (iPtJ+1, 0);
            h_avg_jetares[iSys][iR][iEta]->SetBinContent (iPtJ+1, -1);
            h_avg_jetares[iSys][iR][iEta]->SetBinError (iPtJ+1, 0);
          }

        } // end loop over iPtJ

      } // end loop over iEta

      c->Print ((plotFileName + "]").c_str (), "pdf");


      // fit jet energy resolution to functional form
      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        {
          TF1* f = new TF1 (Form ("f_r%i_avg_jer_%s_iEta%i", r, sys.c_str (), iEta), "sqrt([0]*[0] + [1]*[1]/x + pow([2]/x, 2))", 12, pTJBins[nPtJBins]);
          f_avg_jer[iSys][iR][iEta] = f;

          TH1D* h = h_avg_jer[iSys][iR][iEta];
          h->Fit (f, "RN0Q");
        }

        {
          TF1* f = new TF1 (Form ("f_r%i_avg_jptr_%s_iEta%i", r, sys.c_str (), iEta), "sqrt([0]*[0] + [1]*[1]/x + pow([2]/x, 2))", 12, pTJBins[nPtJBins]);
          f_avg_jptr[iSys][iR][iEta] = f;

          TH1D* h = h_avg_jptr[iSys][iR][iEta];
          h->Fit (f, "RN0Q");
        }

      } // end loop over iEta

    } // end loop over iR

  } // end loop over iSys



  for (int iSys : systems) {

    for (int iR = 0; iR < nRadii; iR++) {

      h2_avg_jpts[iSys][iR]->Write ();
      h2_avg_jptr[iSys][iR]->Write ();

      h2_avg_jes[iSys][iR]->Write ();
      h2_avg_jer[iSys][iR]->Write ();

      h2_avg_jetacorr[iSys][iR]->Write ();
      h2_avg_jetacorr[iSys][iR]->Write ();

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        h_avg_jpts[iSys][iR][iEta]->Write ();
        h_avg_jptr[iSys][iR][iEta]->Write ();

        f_avg_jptr[iSys][iR][iEta]->Write ();

        h_avg_jes[iSys][iR][iEta]->Write ();
        h_avg_jer[iSys][iR][iEta]->Write ();

        f_avg_jer[iSys][iR][iEta]->Write ();

        h_avg_jetacorr[iSys][iR][iEta]->Write ();
        h_avg_jetares[iSys][iR][iEta]->Write ();

        h2_jpts_integratedEta[iSys][iR][iEta]->Write ();

        h2_jes_integratedEta[iSys][iR][iEta]->Write ();

        h2_jetacorr_integratedEta[iSys][iR][iEta]->Write ();

      } // end loop over iEta

      h2_jet_ptreco_pttruth[iSys][iR]->Write ();

      h2_jet_ereco_etruth[iSys][iR]->Write ();

    } // end loop oveer iR

  } // end loop over iSys

  outFile->Close ();
  
}


int main (int argn, char** argv) {
  AnalyzeJetEnergyResolution ();
  return 0;
}

#endif
