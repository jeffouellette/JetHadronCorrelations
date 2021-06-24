#ifndef __AnalyzeJetEnergyResolution_C__
#define __AnalyzeJetEnergyResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>

#include <Utilities.h>

#include "Params.h"
#include "LocalUtilities.h"

using namespace JetHadronCorrelations;

typedef TGraphAsymmErrors TGAE;

const int nFinerEtaBins = 56;
const double* finerEtaBins = linspace (-2.8, 2.8, nFinerEtaBins);

const double etaBins[] = {0, 0.3, 0.8, 1.2, 2.1, 2.8};
const int nEtaBins = sizeof (etaBins) / sizeof (etaBins[0]) - 1;

const double enJBins[] = {10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 360, 400, 450, 500};
const int nEnJBins = sizeof (enJBins) / sizeof (enJBins[0]) - 1;

//const vector <int> systems = {0, 1};
const vector <int> systems = {0};

TH1D**** h_r2_jpts = nullptr;
TH1D**** h_r2_jes = nullptr;
TH1D**** h_r2_jetacorr = nullptr;
TH1D**** h_r4_jpts = nullptr;
TH1D**** h_r4_jes = nullptr;
TH1D**** h_r4_jetacorr = nullptr;

TH2D** h2_r2_avg_jpts = nullptr;
TH2D** h2_r2_avg_jptr = nullptr;
TH1D*** h_r2_avg_jpts = nullptr;
TH1D*** h_r2_avg_jptr = nullptr;
TH2D** h2_r4_avg_jpts = nullptr;
TH2D** h2_r4_avg_jptr = nullptr;
TH1D*** h_r4_avg_jpts = nullptr;
TH1D*** h_r4_avg_jptr = nullptr;

TH2D** h2_r2_avg_jes = nullptr;
TH2D** h2_r2_avg_jer = nullptr;
TH1D*** h_r2_avg_jes = nullptr;
TH1D*** h_r2_avg_jer = nullptr;
TH2D** h2_r4_avg_jes = nullptr;
TH2D** h2_r4_avg_jer = nullptr;
TH1D*** h_r4_avg_jes = nullptr;
TH1D*** h_r4_avg_jer = nullptr;

TH2D** h2_r2_avg_jetacorr = nullptr;
TH1D*** h_r2_avg_jetacorr = nullptr;
TH2D** h2_r2_avg_jetares = nullptr;
TH1D*** h_r2_avg_jetares = nullptr;
TH2D** h2_r4_avg_jetacorr = nullptr;
TH1D*** h_r4_avg_jetacorr = nullptr;
TH2D** h2_r4_avg_jetares = nullptr;
TH1D*** h_r4_avg_jetares = nullptr;

void AnalyzeJetEnergyResolution () {

  TFile* inFile = new TFile (Form ("%s/JetEnergyResolution/Nominal/allSamples.root", rootPath.Data ()), "read");

  const int nFinerEtaBins = 90;
  const double* finerEtaBins = linspace (-4.5, 4.5, nFinerEtaBins);

  const double enJBins[] = {10, 12, 15, 18, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 360, 400, 450, 500};
  const int nEnJBins = sizeof (enJBins) / sizeof (enJBins[0]) - 1;

  h_r2_jpts = new TH1D***[2];
  h_r2_jes = new TH1D***[2];
  h_r2_jetacorr = new TH1D***[2];

  h_r4_jpts = new TH1D***[2];
  h_r4_jes = new TH1D***[2];
  h_r4_jetacorr = new TH1D***[2];

  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    h_r2_jpts[iSys] = new TH1D**[nEnJBins];
    h_r2_jes[iSys] = new TH1D**[nEnJBins];
    h_r2_jetacorr[iSys] = new TH1D**[nEnJBins];

    h_r4_jpts[iSys] = new TH1D**[nEnJBins];
    h_r4_jes[iSys] = new TH1D**[nEnJBins];
    h_r4_jetacorr[iSys] = new TH1D**[nEnJBins];

    for (int iEnJ = 0; iEnJ < nEnJBins; iEnJ++) {

      h_r2_jpts[iSys][iEnJ] = new TH1D*[nFinerEtaBins];
      h_r2_jes[iSys][iEnJ] = new TH1D*[nFinerEtaBins];
      h_r2_jetacorr[iSys][iEnJ] = new TH1D*[nFinerEtaBins];

      h_r4_jpts[iSys][iEnJ] = new TH1D*[nFinerEtaBins];
      h_r4_jes[iSys][iEnJ] = new TH1D*[nFinerEtaBins];
      h_r4_jetacorr[iSys][iEnJ] = new TH1D*[nFinerEtaBins];

      for (int iEta = 0; iEta < nFinerEtaBins; iEta++) {

        h_r2_jpts[iSys][iEnJ][iEta] = (TH1D*) inFile->Get (Form ("h_r2_jpts_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta));
        h_r2_jes[iSys][iEnJ][iEta] = (TH1D*) inFile->Get (Form ("h_r2_jes_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta));
        h_r2_jetacorr[iSys][iEnJ][iEta] = (TH1D*) inFile->Get (Form ("h_r2_jetacorr_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta));

        h_r4_jpts[iSys][iEnJ][iEta] = (TH1D*) inFile->Get (Form ("h_r4_jpts_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta));
        h_r4_jes[iSys][iEnJ][iEta] = (TH1D*) inFile->Get (Form ("h_r4_jes_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta));
        h_r4_jetacorr[iSys][iEnJ][iEta] = (TH1D*) inFile->Get (Form ("h_r4_jetacorr_%s_iEnJ%i_iEta%i", sys.Data (), iEnJ, iEta));

      } // end loop over iEta
    } // end loop over iEnJ
  } // end loop over iSys


  TFile* outFile = new TFile (Form ("%s/JetEnergyResolution/Nominal/summary.root", rootPath.Data ()), "recreate");

  h2_r2_avg_jpts = new TH2D*[2];
  h2_r2_avg_jptr = new TH2D*[2];
  h_r2_avg_jpts = new TH1D**[2];
  h_r2_avg_jptr = new TH1D**[2];

  h2_r4_avg_jpts = new TH2D*[2];
  h2_r4_avg_jptr = new TH2D*[2];
  h_r4_avg_jpts = new TH1D**[2];
  h_r4_avg_jptr = new TH1D**[2];

  h2_r2_avg_jes = new TH2D*[2];
  h2_r2_avg_jer = new TH2D*[2];
  h_r2_avg_jes = new TH1D**[2];
  h_r2_avg_jer = new TH1D**[2];

  h2_r4_avg_jes = new TH2D*[2];
  h2_r4_avg_jer = new TH2D*[2];
  h_r4_avg_jes = new TH1D**[2];
  h_r4_avg_jer = new TH1D**[2];

  h2_r2_avg_jetacorr = new TH2D*[2];
  h2_r2_avg_jetares = new TH2D*[2];
  h_r2_avg_jetacorr = new TH1D**[2];
  h_r2_avg_jetares = new TH1D**[2];

  h2_r4_avg_jetacorr = new TH2D*[2];
  h2_r4_avg_jetares = new TH2D*[2];
  h_r4_avg_jetacorr = new TH1D**[2];
  h_r4_avg_jetares = new TH1D**[2];

  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    h2_r2_avg_jpts[iSys] = new TH2D (Form ("h2_r2_avg_jpts_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);
    h2_r2_avg_jptr[iSys] = new TH2D (Form ("h2_r2_avg_jptr_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);

    h2_r4_avg_jpts[iSys] = new TH2D (Form ("h2_r4_avg_jpts_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);
    h2_r4_avg_jptr[iSys] = new TH2D (Form ("h2_r4_avg_jptr_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);

    h2_r2_avg_jes[iSys] = new TH2D (Form ("h2_r2_avg_jes_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);
    h2_r2_avg_jer[iSys] = new TH2D (Form ("h2_r2_avg_jer_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);

    h2_r4_avg_jes[iSys] = new TH2D (Form ("h2_r4_avg_jes_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);
    h2_r4_avg_jer[iSys] = new TH2D (Form ("h2_r4_avg_jer_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);

    h2_r2_avg_jetacorr[iSys] = new TH2D (Form ("h2_r2_avg_jetacorr_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#eta_{reco} - #eta_{truth}#GT", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);
    h2_r2_avg_jetares[iSys] = new TH2D (Form ("h2_r2_avg_jetares_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma#left[#eta_{reco} - #eta_{truth}#GT#right]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);

    h2_r4_avg_jetacorr[iSys] = new TH2D (Form ("h2_r4_avg_jetacorr_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#LT#eta_{reco} - #eta_{truth}#GT", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);
    h2_r4_avg_jetares[iSys] = new TH2D (Form ("h2_r4_avg_jetares_%s", sys.Data ()), ";#it{E}_{truth} [GeV];#eta_{truth};#sigma#left[#eta_{reco} - #eta_{truth}#GT#right]", nEnJBins, enJBins, nFinerEtaBins, finerEtaBins);


    h_r2_avg_jpts[iSys] = new TH1D*[nEtaBins+1];
    h_r2_avg_jptr[iSys] = new TH1D*[nEtaBins+1];

    h_r4_avg_jpts[iSys] = new TH1D*[nEtaBins+1];
    h_r4_avg_jptr[iSys] = new TH1D*[nEtaBins+1];

    h_r2_avg_jes[iSys] = new TH1D*[nEtaBins+1];
    h_r2_avg_jer[iSys] = new TH1D*[nEtaBins+1];

    h_r4_avg_jes[iSys] = new TH1D*[nEtaBins+1];
    h_r4_avg_jer[iSys] = new TH1D*[nEtaBins+1];

    h_r2_avg_jetacorr[iSys] = new TH1D*[nEtaBins+1];
    h_r2_avg_jetares[iSys] = new TH1D*[nEtaBins+1];

    h_r4_avg_jetacorr[iSys] = new TH1D*[nEtaBins+1];
    h_r4_avg_jetares[iSys] = new TH1D*[nEtaBins+1];

    for (int iEta = 0; iEta <= nEtaBins; iEta++) {

      h_r2_avg_jpts[iSys][iEta] = new TH1D (Form ("h_r2_avg_jpts_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]", nEnJBins, enJBins);
      h_r2_avg_jptr[iSys][iEta] = new TH1D (Form ("h_r2_avg_jptr_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]", nEnJBins, enJBins);

      h_r4_avg_jpts[iSys][iEta] = new TH1D (Form ("h_r4_avg_jpts_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT [%]", nEnJBins, enJBins);
      h_r4_avg_jptr[iSys][iEta] = new TH1D (Form ("h_r4_avg_jptr_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma / #mu #left[#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#right] [%]", nEnJBins, enJBins);

      h_r2_avg_jes[iSys][iEta] = new TH1D (Form ("h_r2_avg_jes_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]", nEnJBins, enJBins);
      h_r2_avg_jer[iSys][iEta] = new TH1D (Form ("h_r2_avg_jer_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]", nEnJBins, enJBins);

      h_r4_avg_jes[iSys][iEta] = new TH1D (Form ("h_r4_avg_jes_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#it{E}_{reco} / #it{E}_{truth}#GT [%]", nEnJBins, enJBins);
      h_r4_avg_jer[iSys][iEta] = new TH1D (Form ("h_r4_avg_jer_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma / #mu #left[#it{E}_{reco} / #it{E}_{truth}#right] [%]", nEnJBins, enJBins);

      h_r2_avg_jetacorr[iSys][iEta] = new TH1D (Form ("h_r2_avg_jetacorr_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#eta_{reco} - #eta_{truth}#GT", nEnJBins, enJBins);
      h_r2_avg_jetares[iSys][iEta] = new TH1D (Form ("h_r2_avg_jetares_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma#left[#eta_{reco} - #eta_{truth}#GT#right]", nEnJBins, enJBins);

      h_r4_avg_jetacorr[iSys][iEta] = new TH1D (Form ("h_r4_avg_jetacorr_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#LT#eta_{reco} - #eta_{truth}#GT", nEnJBins, enJBins);
      h_r4_avg_jetares[iSys][iEta] = new TH1D (Form ("h_r4_avg_jetares_%s_iEta%i", sys.Data (), iEta), ";#it{E}_{truth} [GeV];#sigma#left[#eta_{reco} - #eta_{truth}#GT#right]", nEnJBins, enJBins);

    } // end loop over iEta
  } // end loop over iSys



  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iEnJ = 0; iEnJ < nEnJBins; iEnJ++) {

      TH1D** h_r2_jpts_integratedEta = new TH1D*[nEtaBins+1];
      for (int iEta = 0; iEta <= nEtaBins; iEta++) {
        h_r2_jpts_integratedEta[iEta] = new TH1D (Form ("h_r2_jpts_integratedEta_%s_iEta%i", sys.Data (), iEta), "#it{p}_{T}^{reco} / #it{p}_{T}^{truth}", 140, 0.3, 1.7);
        h_r2_jpts_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < nEtaBins) {
          if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
            break;
          iEta++;
        }

        h_r2_jpts_integratedEta[iEta]->Add (h_r2_jpts[iSys][iEnJ][iFinerEta]);
        h_r2_jpts_integratedEta[nEtaBins]->Add (h_r2_jpts[iSys][iEnJ][iFinerEta]);

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.3, 1.7);
        fit->SetParameter (0, h_r2_jpts[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r2_jpts[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, 0.3), std::fmin (mean+2*sigma, 1.7));

        fit->SetParameter (0, h_r2_jpts[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r2_jpts[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jpts = mean;
        const double jpts_err = mean_err;
        
        const double jptr = sigma / mean;
        const double jptr_err = fabs (jptr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h2_r2_avg_jpts[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jpts * 100);
        h2_r2_avg_jpts[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jpts_err * 100);
        h2_r2_avg_jptr[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jptr * 100);
        h2_r2_avg_jptr[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jptr_err * 100);
      }

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.3, 1.7);
        fit->SetParameter (0, h_r2_jpts_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r2_jpts_integratedEta[iEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, 0.3), std::fmin (mean+2*sigma, 1.7));

        fit->SetParameter (0, h_r2_jpts_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r2_jpts_integratedEta[iEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jpts = mean;
        const double jpts_err = mean_err;

        const double jptr = sigma / mean;
        const double jptr_err = fabs (jptr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h_r2_avg_jpts[iSys][iEta]->SetBinContent (iEnJ+1, jpts * 100);
        h_r2_avg_jpts[iSys][iEta]->SetBinError (iEnJ+1, jpts_err * 100);
        h_r2_avg_jptr[iSys][iEta]->SetBinContent (iEnJ+1, jptr * 100);
        h_r2_avg_jptr[iSys][iEta]->SetBinError (iEnJ+1, jptr_err * 100);

        delete h_r2_jpts_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_r2_jpts_integratedEta;



      TH1D** h_r4_jpts_integratedEta = new TH1D*[nEtaBins+1];
      for (int iEta = 0; iEta <= nEtaBins; iEta++) {
        h_r4_jpts_integratedEta[iEta] = new TH1D (Form ("h_r4_jpts_integratedEta_%s_iEta%i", sys.Data (), iEta), "#it{p}_{T}^{reco} / #it{p}_{T}^{truth}", 140, 0.3, 1.7);
        h_r4_jpts_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < nEtaBins) {
          if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
            break;
          iEta++;
        }

        h_r4_jpts_integratedEta[iEta]->Add (h_r4_jpts[iSys][iEnJ][iFinerEta]);
        h_r4_jpts_integratedEta[nEtaBins]->Add (h_r4_jpts[iSys][iEnJ][iFinerEta]);

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.3, 1.7);
        fit->SetParameter (0, h_r4_jpts[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r4_jpts[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, 0.3), std::fmin (mean+2*sigma, 1.7));

        fit->SetParameter (0, h_r4_jpts[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r4_jpts[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jpts = mean;
        const double jpts_err = mean_err;
        
        const double jptr = sigma / mean;
        const double jptr_err = fabs (jptr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h2_r4_avg_jpts[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jpts * 100);
        h2_r4_avg_jpts[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jpts_err * 100);
        h2_r4_avg_jptr[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jptr * 100);
        h2_r4_avg_jptr[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jptr_err * 100);
      }

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.3, 1.7);
        fit->SetParameter (0, h_r4_jpts_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r4_jpts_integratedEta[iEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, 0.3), std::fmin (mean+2*sigma, 1.7));

        fit->SetParameter (0, h_r4_jpts_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r4_jpts_integratedEta[iEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jpts = mean;
        const double jpts_err = mean_err;

        const double jptr = sigma / mean;
        const double jptr_err = fabs (jptr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h_r4_avg_jpts[iSys][iEta]->SetBinContent (iEnJ+1, jpts * 100);
        h_r4_avg_jpts[iSys][iEta]->SetBinError (iEnJ+1, jpts_err * 100);
        h_r4_avg_jptr[iSys][iEta]->SetBinContent (iEnJ+1, jptr * 100);
        h_r4_avg_jptr[iSys][iEta]->SetBinError (iEnJ+1, jptr_err * 100);

        delete h_r4_jpts_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_r4_jpts_integratedEta;



      TH1D** h_r2_jes_integratedEta = new TH1D*[nEtaBins+1];
      for (int iEta = 0; iEta <= nEtaBins; iEta++) {
        h_r2_jes_integratedEta[iEta] = new TH1D (Form ("h_r2_jes_integratedEta_%s_iEta%i", sys.Data (), iEta), "#it{E}_{reco} / #it{E}_{truth}", 140, 0.3, 1.7);
        h_r2_jes_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < nEtaBins) {
          if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
            break;
          iEta++;
        }

        h_r2_jes_integratedEta[iEta]->Add (h_r2_jes[iSys][iEnJ][iFinerEta]);
        h_r2_jes_integratedEta[nEtaBins]->Add (h_r2_jes[iSys][iEnJ][iFinerEta]);

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.3, 1.7);
        fit->SetParameter (0, h_r2_jes[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r2_jes[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, 0.3), std::fmin (mean+2*sigma, 1.7));

        fit->SetParameter (0, h_r2_jes[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r2_jes[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jes = mean;
        const double jes_err = mean_err;
        
        const double jer = sigma / mean;
        const double jer_err = fabs (jer) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h2_r2_avg_jes[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jes * 100);
        h2_r2_avg_jes[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jes_err * 100);
        h2_r2_avg_jer[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jer * 100);
        h2_r2_avg_jer[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jer_err * 100);
      }

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.3, 1.7);
        fit->SetParameter (0, h_r2_jes_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r2_jes_integratedEta[iEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, 0.3), std::fmin (mean+2*sigma, 1.7));

        fit->SetParameter (0, h_r2_jes_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r2_jes_integratedEta[iEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jes = mean;
        const double jes_err = mean_err;

        const double jer = sigma / mean;
        const double jer_err = fabs (jer) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h_r2_avg_jes[iSys][iEta]->SetBinContent (iEnJ+1, jes * 100);
        h_r2_avg_jes[iSys][iEta]->SetBinError (iEnJ+1, jes_err * 100);
        h_r2_avg_jer[iSys][iEta]->SetBinContent (iEnJ+1, jer * 100);
        h_r2_avg_jer[iSys][iEta]->SetBinError (iEnJ+1, jer_err * 100);

        delete h_r2_jes_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_r2_jes_integratedEta;



      TH1D** h_r4_jes_integratedEta = new TH1D*[nEtaBins+1];
      for (int iEta = 0; iEta <= nEtaBins; iEta++) {
        h_r4_jes_integratedEta[iEta] = new TH1D (Form ("h_r4_jes_integratedEta_%s_iEta%i", sys.Data (), iEta), "#it{E}_{reco} / #it{E}_{truth}", 140, 0.3, 1.7);
        h_r4_jes_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < nEtaBins) {
          if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
            break;
          iEta++;
        }

        h_r4_jes_integratedEta[iEta]->Add (h_r4_jes[iSys][iEnJ][iFinerEta]);
        h_r4_jes_integratedEta[nEtaBins]->Add (h_r4_jes[iSys][iEnJ][iFinerEta]);

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.3, 1.7);
        fit->SetParameter (0, h_r4_jes[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r4_jes[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, 0.3), std::fmin (mean+2*sigma, 1.7));

        fit->SetParameter (0, h_r4_jes[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r4_jes[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jes = mean;
        const double jes_err = mean_err;
        
        const double jer = sigma / mean;
        const double jer_err = fabs (jer) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h2_r4_avg_jes[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jes * 100);
        h2_r4_avg_jes[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jes_err * 100);
        h2_r4_avg_jer[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jer * 100);
        h2_r4_avg_jer[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jer_err * 100);
      }

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.3, 1.7);
        fit->SetParameter (0, h_r4_jes_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r4_jes_integratedEta[iEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, 0.3), std::fmin (mean+2*sigma, 1.7));

        fit->SetParameter (0, h_r4_jes_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r4_jes_integratedEta[iEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jes = mean;
        const double jes_err = mean_err;

        const double jer = sigma / mean;
        const double jer_err = fabs (jer) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h_r4_avg_jes[iSys][iEta]->SetBinContent (iEnJ+1, jes * 100);
        h_r4_avg_jes[iSys][iEta]->SetBinError (iEnJ+1, jes_err * 100);
        h_r4_avg_jer[iSys][iEta]->SetBinContent (iEnJ+1, jer * 100);
        h_r4_avg_jer[iSys][iEta]->SetBinError (iEnJ+1, jer_err * 100);

        delete h_r4_jes_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_r4_jes_integratedEta;



      TH1D** h_r2_jetacorr_integratedEta = new TH1D*[nEtaBins+1];
      for (int iEta = 0; iEta <= nEtaBins; iEta++) {
        h_r2_jetacorr_integratedEta[iEta] = new TH1D (Form ("h_r2_jetacorr_integratedEta_%s_iEta%i", sys.Data (), iEta), "#eta_{reco} - #eta_{truth}", 80, -0.2, 0.2);
        h_r2_jetacorr_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < nEtaBins) {
          if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
            break;
          iEta++;
        }

        h_r2_jetacorr_integratedEta[iEta]->Add (h_r2_jetacorr[iSys][iEnJ][iFinerEta]);
        h_r2_jetacorr_integratedEta[nEtaBins]->Add (h_r2_jetacorr[iSys][iEnJ][iFinerEta]);

        TF1* fit = new TF1 ("fit", "gaus(0)", -0.2, 0.2);
        fit->SetParameter (0, h_r2_jetacorr[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r2_jetacorr[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, -0.2), std::fmin (mean+2*sigma, 0.2));

        fit->SetParameter (0, h_r2_jetacorr[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r2_jetacorr[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jetacorr = mean;
        const double jetacorr_err = mean_err;
        
        const double jetares = sigma;
        const double jetares_err = sigma_err;

        h2_r2_avg_jetacorr[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jetacorr * 100);
        h2_r2_avg_jetacorr[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jetacorr_err * 100);
        h2_r2_avg_jetares[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jetares * 100);
        h2_r2_avg_jetares[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jetares_err * 100);
      }

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        TF1* fit = new TF1 ("fit", "gaus(0)", -0.2, 0.2);
        fit->SetParameter (0, h_r2_jetacorr_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r2_jetacorr_integratedEta[iEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, -0.2), std::fmin (mean+2*sigma, 0.2));

        fit->SetParameter (0, h_r2_jetacorr_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r2_jetacorr_integratedEta[iEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jetacorr = mean;
        const double jetacorr_err = mean_err;

        const double jetares = sigma;
        const double jetares_err = sigma_err;

        h_r2_avg_jetacorr[iSys][iEta]->SetBinContent (iEnJ+1, jetacorr);
        h_r2_avg_jetacorr[iSys][iEta]->SetBinError (iEnJ+1, jetacorr_err);
        h_r2_avg_jetares[iSys][iEta]->SetBinContent (iEnJ+1, jetares);
        h_r2_avg_jetares[iSys][iEta]->SetBinError (iEnJ+1, jetares_err);

        delete h_r2_jetacorr_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_r2_jetacorr_integratedEta;



      TH1D** h_r4_jetacorr_integratedEta = new TH1D*[nEtaBins+1];
      for (int iEta = 0; iEta <= nEtaBins; iEta++) {
        h_r4_jetacorr_integratedEta[iEta] = new TH1D (Form ("h_r4_jetacorr_integratedEta_%s_iEta%i", sys.Data (), iEta), "#eta_{reco} - #eta_{truth}", 80, -0.2, 0.2);
        h_r4_jetacorr_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < nFinerEtaBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtaBins[iFinerEta] + finerEtaBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < nEtaBins) {
          if (etaBins[iEta] < binCenter && binCenter < etaBins[iEta+1])
            break;
          iEta++;
        }

        h_r4_jetacorr_integratedEta[iEta]->Add (h_r4_jetacorr[iSys][iEnJ][iFinerEta]);
        h_r4_jetacorr_integratedEta[nEtaBins]->Add (h_r4_jetacorr[iSys][iEnJ][iFinerEta]);

        TF1* fit = new TF1 ("fit", "gaus(0)", -0.2, 0.2);
        fit->SetParameter (0, h_r4_jetacorr[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r4_jetacorr[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, -0.2), std::fmin (mean+2*sigma, 0.2));

        fit->SetParameter (0, h_r4_jetacorr[iSys][iEnJ][iFinerEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r4_jetacorr[iSys][iEnJ][iFinerEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jetacorr = mean;
        const double jetacorr_err = mean_err;
        
        const double jetares = sigma;
        const double jetares_err = sigma_err;

        h2_r4_avg_jetacorr[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jetacorr * 100);
        h2_r4_avg_jetacorr[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jetacorr_err * 100);
        h2_r4_avg_jetares[iSys]->SetBinContent (iEnJ+1, iFinerEta+1, jetares * 100);
        h2_r4_avg_jetares[iSys]->SetBinError (iEnJ+1, iFinerEta+1, jetares_err * 100);
      }

      for (int iEta = 0; iEta <= nEtaBins; iEta++) {

        TF1* fit = new TF1 ("fit", "gaus(0)", -0.2, 0.2);
        fit->SetParameter (0, h_r4_jetacorr_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_r4_jetacorr_integratedEta[iEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", std::fmax (mean-2*sigma, -0.2), std::fmin (mean+2*sigma, 0.2));

        fit->SetParameter (0, h_r4_jetacorr_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_r4_jetacorr_integratedEta[iEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double jetacorr = mean;
        const double jetacorr_err = mean_err;

        const double jetares = sigma;
        const double jetares_err = sigma_err;

        h_r4_avg_jetacorr[iSys][iEta]->SetBinContent (iEnJ+1, jetacorr);
        h_r4_avg_jetacorr[iSys][iEta]->SetBinError (iEnJ+1, jetacorr_err);
        h_r4_avg_jetares[iSys][iEta]->SetBinContent (iEnJ+1, jetares);
        h_r4_avg_jetares[iSys][iEta]->SetBinError (iEnJ+1, jetares_err);

        delete h_r4_jetacorr_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_r4_jetacorr_integratedEta;
    } // end loop over iEnJ
  } // end loop over iSys


  for (int iSys : systems) {
    h2_r2_avg_jpts[iSys]->Write ();
    h2_r2_avg_jptr[iSys]->Write ();

    h2_r4_avg_jpts[iSys]->Write ();
    h2_r4_avg_jptr[iSys]->Write ();

    h2_r2_avg_jes[iSys]->Write ();
    h2_r2_avg_jer[iSys]->Write ();

    h2_r4_avg_jes[iSys]->Write ();
    h2_r4_avg_jer[iSys]->Write ();

    h2_r2_avg_jetacorr[iSys]->Write ();
    h2_r2_avg_jetacorr[iSys]->Write ();

    h2_r4_avg_jetacorr[iSys]->Write ();
    h2_r4_avg_jetacorr[iSys]->Write ();

    for (int iEta = 0; iEta <= nEtaBins; iEta++) {
      h_r2_avg_jpts[iSys][iEta]->Write ();
      h_r2_avg_jptr[iSys][iEta]->Write ();

      h_r4_avg_jpts[iSys][iEta]->Write ();
      h_r4_avg_jptr[iSys][iEta]->Write ();

      h_r2_avg_jes[iSys][iEta]->Write ();
      h_r2_avg_jer[iSys][iEta]->Write ();

      h_r4_avg_jes[iSys][iEta]->Write ();
      h_r4_avg_jer[iSys][iEta]->Write ();

      h_r2_avg_jetacorr[iSys][iEta]->Write ();
      h_r2_avg_jetares[iSys][iEta]->Write ();

      h_r4_avg_jetacorr[iSys][iEta]->Write ();
      h_r4_avg_jetares[iSys][iEta]->Write ();
    } // end loop over iEta
  } // end loop over iSys

  outFile->Close ();
  
}

#endif
