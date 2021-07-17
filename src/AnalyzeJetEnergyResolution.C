#ifndef __AnalyzeJetEnergyResolution_C__
#define __AnalyzeJetEnergyResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>

#include <Utilities.h>
#include <ArrayTemplates.h>

#include "Params.h"
#include "LocalUtilities.h"

using namespace JetHadronCorrelations;

typedef TGraphAsymmErrors TGAE;

// Number of sigma to use for second fit
const float nSigFit2 = 1.4;

// Minimum number of histogram entries to attempt a fit (400 is an arbitrary choice)
const int minEntriesForFit = 400;


/**
 * Wrapper function for performing a fit to a Gaussian TF1.
 */
TF1* DoGaussianFit (TH1D* h, double min = 0.0, double max = 2.4) {
  std::cout << "Fitting " << h->GetName () << " to Gaussian." << std::endl;
  // Check we have "enough" entries to perform a fit.
  if (h->Integral () < minEntriesForFit) {
    std::cout << "Not enough counts for fit! Will skip this histogram." << std::endl;
    return nullptr;
  }
  else
    std::cout << "Fitting " << h->GetName () << " to gaussian" << std::endl;


  TF1* fit = new TF1 ("fit", "gaus(0)", min, max);
  fit->SetParameter (0, h->Integral ());
  fit->SetParameter (1, h->GetMean ());
  fit->SetParameter (2, h->GetStdDev ());

  h->Fit (fit, "RN0Q");

  const double mean = fit->GetParameter (1);
  const double sigma = fit->GetParameter (2);

  delete fit;
  fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-nSigFit2*sigma, 0.0), std::fmin (mean+nSigFit2*sigma, 2.4));

  fit->SetParameter (0, h->Integral ());
  fit->SetParameter (1, mean);
  fit->SetParameter (2, sigma);

  h->Fit (fit, "RN0Q");

  return fit;
}



void AnalyzeJetEnergyResolution () {

  const int nFinerEtaBins = 56;
  const double* finerEtaBins = linspace (-2.8, 2.8, nFinerEtaBins);
  
  const double etaBins[] = {-2.8, -2.1, -1.2, -0.8, -0.3, 0, 0.3, 0.8, 1.2, 2.1, 2.8};
  const int nEtaBins = sizeof (etaBins) / sizeof (etaBins[0]) - 1;
  
  const double pTJBins[] = {10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 360, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300};
  const int nPtJBins = sizeof (pTJBins) / sizeof (pTJBins[0]) - 1;

  const int nRespBins = 240;
  const double* respBins = linspace (0, 2.4, nRespBins);

  const int nEtaRespBins = 80;
  const double* etaRespBins = linspace (-0.2, 0.2, nEtaRespBins);
  
  const std::vector <int> systems = {0, 1}; // 0 = pp, 1 = pPb
  const int nSystems = systems.size ();
  
  const std::vector <JetRadius> radii = {JetRadius::R0p2, JetRadius::R0p4};
  const int nRadii = radii.size ();


  TFile* inFile = new TFile (Form ("%s/JetEnergyResolution/Nominal/allSamples.root", rootPath.Data ()), "read");

  TH1D***** h_jpts     = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nFinerEtaBins);
  TH1D***** h_jes      = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nFinerEtaBins);
  TH1D***** h_jetacorr = Get4DArray <TH1D*> (nSystems, nRadii, nPtJBins, nFinerEtaBins);

  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iR = 0; iR < nRadii; iR++) {

      const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

      for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        for (int iEta = 0; iEta < nFinerEtaBins; iEta++) {

          h_jpts[iSys][iR][iPtJ][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_jpts_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta));
          h_jes[iSys][iR][iPtJ][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_jes_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta));
          h_jetacorr[iSys][iR][iPtJ][iEta] = (TH1D*) inFile->Get (Form ("h_r%i_jetacorr_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta));

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

  TH2D*** h2_avg_jes  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_avg_jer  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH1D**** h_avg_jes  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);
  TH1D**** h_avg_jer  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);

  TH2D*** h2_avg_jetacorr = Get2DArray <TH2D*> (nSystems, nRadii);
  TH2D*** h2_avg_jetares  = Get2DArray <TH2D*> (nSystems, nRadii);
  TH1D**** h_avg_jetacorr = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);
  TH1D**** h_avg_jetares  = Get3DArray <TH1D*> (nSystems, nRadii, nEtaBins+1);

  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iR = 0; iR < nRadii; iR++) {

      const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

      h2_avg_jpts[iSys][iR] = new TH2D (Form ("h2_r%i_avg_jpts_%s", r, sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);
      h2_avg_jptr[iSys][iR] = new TH2D (Form ("h2_r%i_avg_jptr_%s", r, sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);

      h2_avg_jes[iSys][iR] = new TH2D (Form ("h2_r%i_avg_jes_%s", r, sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);
      h2_avg_jer[iSys][iR] = new TH2D (Form ("h2_r%i_avg_jer_%s", r, sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);

      h2_avg_jetacorr[iSys][iR] = new TH2D (Form ("h2_r%i_avg_jetacorr_%s", r, sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#eta_{reco} - #eta_{truth}#GT", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);
      h2_avg_jetares[iSys][iR] = new TH2D (Form ("h2_r%i_avg_jetares_%s", r, sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma#left[#eta_{reco} - #eta_{truth}#GT#right]", nPtJBins, pTJBins, nFinerEtaBins, finerEtaBins);


      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        h_avg_jpts[iSys][iR][iEta] = new TH1D (Form ("h_r%i_avg_jpts_%s_iEta%i", r, sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]", nPtJBins, pTJBins);
        h_avg_jptr[iSys][iR][iEta] = new TH1D (Form ("h_r%i_avg_jptr_%s_iEta%i", r, sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]", nPtJBins, pTJBins);

        h_avg_jes[iSys][iR][iEta] = new TH1D (Form ("h_r%i_avg_jes_%s_iEta%i", r, sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]", nPtJBins, pTJBins);
        h_avg_jer[iSys][iR][iEta] = new TH1D (Form ("h_r%i_avg_jer_%s_iEta%i", r, sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]", nPtJBins, pTJBins);

        h_avg_jetacorr[iSys][iR][iEta] = new TH1D (Form ("h_r%i_avg_jetacorr_%s_iEta%i", r, sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#eta_{reco} - #eta_{truth}#GT", nPtJBins, pTJBins);
        h_avg_jetares[iSys][iR][iEta] = new TH1D (Form ("h_r%i_avg_jetares_%s_iEta%i", r, sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma#left[#eta_{reco} - #eta_{truth}#GT#right]", nPtJBins, pTJBins);


        h2_jpts_integratedEta[iSys][iR][iEta] = new TH2D (Form ("h2_r%i_jpts_integratedEta_%s_iEta%i", r, sys.Data (), iEta), ";#it{p}_{T}^{truth} [GeV];#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", nPtJBins, pTJBins, nRespBins, respBins);
        h2_jpts_integratedEta[iSys][iR][iEta]->Sumw2 ();

        h2_jes_integratedEta[iSys][iR][iEta] = new TH2D (Form ("h2_r%i_jes_integratedEta_%s_iEta%i", r, sys.Data (), iEta), ";#it{p}_{T}^{truth} [GeV];#it{E}_{reco} / #it{E}_{truth};Counts", nPtJBins, pTJBins, nRespBins, respBins);
        h2_jes_integratedEta[iSys][iR][iEta]->Sumw2 ();

        h2_jetacorr_integratedEta[iSys][iR][iEta] = new TH2D (Form ("h2_r%i_jetacorr_integratedEta_%s_iEta%i", r, sys.Data (), iEta), ";#it{p}_{T}^{truth} [GeV];#eta_{reco} - #eta_{truth};Counts", nPtJBins, pTJBins, nEtaRespBins, etaRespBins);
        h2_jetacorr_integratedEta[iSys][iR][iEta]->Sumw2 ();


        for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

          h_jpts_integratedEta[iSys][iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jpts_integratedEta_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", nRespBins, respBins);
          h_jpts_integratedEta[iSys][iR][iPtJ][iEta]->Sumw2 ();

          h_jes_integratedEta[iSys][iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jes_integratedEta_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta), ";#it{E}_{reco} / #it{E}_{truth};Counts", nRespBins, respBins);
          h_jes_integratedEta[iSys][iR][iPtJ][iEta]->Sumw2 ();

          h_jetacorr_integratedEta[iSys][iR][iPtJ][iEta] = new TH1D (Form ("h_r%i_jetacorr_integratedEta_%s_iPtJ%i_iEta%i", r, sys.Data (), iPtJ, iEta), ";#eta_{reco} - #eta_{truth};Counts", nEtaRespBins, etaRespBins);
          h_jetacorr_integratedEta[iSys][iR][iPtJ][iEta]->Sumw2 ();

        } // end loop over iPtJ

      } // end loop over iEta

    } // end loop over iR

  } // end loop over iSys



  for (int iSys : systems) {

    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iR = 0; iR < nRadii; iR++) {

      const int r = (radii[iR] == JetRadius::R0p4 ? 4 : 2);

      for (int iPtJ = 0; iPtJ < nPtJBins; iPtJ++) {

        for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

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

          TF1* fit = DoGaussianFit (h);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

            delete fit;

            const double jpts = mean;
            const double jpts_err = mean_err;
            
            const double jptr = sigma / mean;
            const double jptr_err = fabs (jptr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

            h2_avg_jpts[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jpts * 100);
            h2_avg_jpts[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jpts_err * 100);
            h2_avg_jptr[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jptr * 100);
            h2_avg_jptr[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jptr_err * 100);
          }
          else {
            h2_avg_jpts[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
            h2_avg_jpts[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
            h2_avg_jptr[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
            h2_avg_jptr[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          }

        } // end loop over iFinerEta

        for (int iEta = 0; iEta <= nEtaBins; iEta++) {

          TH1D* h = h_jpts_integratedEta[iSys][iR][iPtJ][iEta];

          TF1* fit = DoGaussianFit (h);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

            delete fit;

            const double jpts = mean;
            const double jpts_err = mean_err;

            const double jptr = sigma / mean;
            const double jptr_err = fabs (jptr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

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

        } // end loop over iEta



        for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

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

          TF1* fit = DoGaussianFit (h);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

            delete fit;

            const double jes = mean;
            const double jes_err = mean_err;
            
            const double jer = sigma / mean;
            const double jer_err = fabs (jer) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

            h2_avg_jes[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jes * 100);
            h2_avg_jes[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jes_err * 100);
            h2_avg_jer[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jer * 100);
            h2_avg_jer[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jer_err * 100);
          }
          else {
            h2_avg_jes[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
            h2_avg_jes[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
            h2_avg_jer[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
            h2_avg_jer[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          }

        } // end loop over iFinerEta

        for (int iEta = 0; iEta <= nEtaBins; iEta++) {

          TH1D* h = h_jes_integratedEta[iSys][iR][iPtJ][iEta];

          TF1* fit = DoGaussianFit (h);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

            delete fit;

            const double jes = mean;
            const double jes_err = mean_err;

            const double jer = sigma / mean;
            const double jer_err = fabs (jer) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

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

        } // end loop over iEta



        for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

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

          TF1* fit = DoGaussianFit (h, -0.2, 0.2);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

            delete fit;

            const double jetacorr = mean;
            const double jetacorr_err = mean_err;
            
            const double jetares = sigma;
            const double jetares_err = sigma_err;

            h2_avg_jetacorr[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jetacorr * 100);
            h2_avg_jetacorr[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jetacorr_err * 100);
            h2_avg_jetares[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, jetares * 100);
            h2_avg_jetares[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, jetares_err * 100);
          }
          else {
            h2_avg_jetacorr[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
            h2_avg_jetacorr[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
            h2_avg_jetares[iSys][iR]->SetBinContent (iPtJ+1, iFinerEta+1, -1);
            h2_avg_jetares[iSys][iR]->SetBinError (iPtJ+1, iFinerEta+1, 0);
          }

        } // end loop over iFinerEta

        for (int iEta = 0; iEta <= nEtaBins; iEta++) {

          TH1D* h = h_jetacorr_integratedEta[iSys][iR][iPtJ][iEta];

          TF1* fit = DoGaussianFit (h, -0.2, 0.2);

          if (fit) {
            const double mean = fit->GetParameter (1);
            const double sigma = fit->GetParameter (2);
            const double mean_err = fit->GetParError (1);
            const double sigma_err = fit->GetParError (2);

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

        } // end loop over iEta

      } // end loop over iPtJ

    } // end loop over iR

  } // end loop over iSys



  for (int iSys : systems) {

    for (int iR = 0; iR < nRadii; iR++) {

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

    } // end loop oveer iR

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

        h_avg_jes[iSys][iR][iEta]->Write ();
        h_avg_jer[iSys][iR][iEta]->Write ();

        h_avg_jetacorr[iSys][iR][iEta]->Write ();
        h_avg_jetares[iSys][iR][iEta]->Write ();

        h2_jpts_integratedEta[iSys][iR][iEta]->Write ();

        h2_jes_integratedEta[iSys][iR][iEta]->Write ();

        h2_jetacorr_integratedEta[iSys][iR][iEta]->Write ();

      } // end loop over iEta

    } // end loop oveer iR

  } // end loop over iSys

  outFile->Close ();
  
}

#endif
