#ifndef __AnalyzeTrackMomentumResolution_C__
#define __AnalyzeTrackMomentumResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>

#include <Utilities.h>

#include "Params.h"
#include "LocalUtilities.h"
#include "TrackMomentumFit.h"

using namespace JetHadronCorrelations;

typedef TGraphAsymmErrors TGAE;

const int nFinerEtaTrkBins = 40;
const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

const double etaTrkBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const int nEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100};
const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

//const vector <int> systems = {0, 1};
const vector <int> systems = {0};

TH1D**** h_tmr = nullptr;
TH2D** h2_avg_tms = nullptr;
TH2D** h2_avg_tmr = nullptr;
TH1D*** h_avg_tms = nullptr;
TH1D*** h_avg_tmr = nullptr;


void AnalyzeTrackMomentumResolution () {

  TFile* inFile = new TFile (Form ("%s/TrackMomentumResolution/Nominal/allSamples.root", rootPath.Data ()), "read");

  h_tmr = new TH1D***[2];
  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");
    h_tmr[iSys] = new TH1D**[nPtchBins];
    for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {
      h_tmr[iSys][iPtch] = new TH1D*[nFinerEtaTrkBins];
      for (int iEta = 0; iEta < nFinerEtaTrkBins; iEta++) {
        h_tmr[iSys][iPtch][iEta] = (TH1D*) inFile->Get (Form ("h_tmr_%s_iPtch%i_iEta%i", sys.Data (), iPtch, iEta));
      } // end loop over iEta
    } // end loop over iPtch
  } // end loop over iSys



  TFile* outFile = new TFile (Form ("%s/TrackMomentumResolution/Nominal/summary.root", rootPath.Data ()), "recreate");

  h2_avg_tms = new TH2D*[2];
  h2_avg_tmr = new TH2D*[2];
  h_avg_tms = new TH1D**[2];
  h_avg_tmr = new TH1D**[2];

  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");
    h2_avg_tms[iSys] = new TH2D (Form ("h2_avg_tms_%s", sys.Data ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMS [%]", nPtchBins, pTchBins, nFinerEtaTrkBins, finerEtaTrkBins);
    h2_avg_tmr[iSys] = new TH2D (Form ("h2_avg_tmr_%s", sys.Data ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMR [%]", nPtchBins, pTchBins, nFinerEtaTrkBins, finerEtaTrkBins);

    h_avg_tms[iSys] = new TH1D*[nEtaTrkBins+1];
    h_avg_tmr[iSys] = new TH1D*[nEtaTrkBins+1];
    for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {
      h_avg_tms[iSys][iEta] = new TH1D (Form ("h_avg_tms_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T}^{truth} [GeV];TMS [%]", nPtchBins, pTchBins);
      h_avg_tmr[iSys][iEta] = new TH1D (Form ("h_avg_tmr_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T}^{truth} [GeV];TMR [%]", nPtchBins, pTchBins);
    } // end loop over iEta
  } // end loop over iSys


  TrackMomentumFit tmf;

  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {

      TH1D** h_tmr_integratedEta = new TH1D*[nEtaTrkBins+1];
      for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {
        h_tmr_integratedEta[iEta] = new TH1D (Form ("h_tmr_integratedEta_%s_iEta%i", sys.Data (), iEta), "Track resolution, #it{p}_{T}^{truth} / #it{p}_{T}^{reco} - 1", 2000, -1.0, 4.0);
        h_tmr_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < nFinerEtaTrkBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtaTrkBins[iFinerEta] + finerEtaTrkBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < nEtaTrkBins) {
          if (etaTrkBins[iEta] < binCenter && binCenter < etaTrkBins[iEta+1])
            break;
          iEta++;
        }

        h_tmr_integratedEta[iEta]->Add (h_tmr[iSys][iPtch][iFinerEta]);
        h_tmr_integratedEta[nEtaTrkBins]->Add (h_tmr[iSys][iPtch][iFinerEta]);

        TF1* fit = new TF1 (Form ("f_tmr_%s_iPtch%i_iEta%i", sys.Data (), iPtch, iEta), &tmf, -1.0, 4.0, 9);

        TH1D* h = h_tmr[iSys][iPtch][iFinerEta];

        const double initRMS = h->GetRMS ();

        fit->SetParameter (0, -0.8);
        fit->SetParameter (1, 0.8);
        fit->SetParameter (2, 0);
        fit->SetParameter (3, 0.5 * 0.68 * h->Integral (h->FindBin (-initRMS), h->FindBin (initRMS)) / std::sqrt (2*M_PI*initRMS*initRMS));
        fit->SetParameter (4, 0.5 * 0.68 * h->Integral (h->FindBin (-initRMS), h->FindBin (initRMS)) / std::sqrt (2*M_PI*initRMS*initRMS));
        std::cout << "rms, int = " << initRMS << ", " << fit->GetParameter (3) << std::endl;
        fit->SetParameter (5, initRMS);
        fit->SetParameter (6, initRMS);
        fit->SetParameter (7, -2);
        fit->SetParameter (8, -2);

        fit->SetParLimits (0, -1., 0.0001);
        fit->SetParLimits (1, 0.0001, 4.);
        fit->SetParLimits (2, -0.5, 0.5);
        fit->SetParLimits (3, 0., 1e16);
        fit->SetParLimits (4, 0., 1e16);
        fit->SetParLimits (5, 0.0001, 2.);
        fit->SetParLimits (6, 0.0001, 2.);
        fit->SetParLimits (7, 0.0001, 10.);
        fit->SetParLimits (8, 0.0001, 10.);

        h->Fit (fit, "RN0");

        const double mean = fit->GetParameter (2);
        const double mean_err = fit->GetParError (1);
        const double sigma = std::sqrt (std::pow (fit->GetParameter (5), 2) + std::pow (fit->GetParameter (6), 2));
        const double sigma_err = std::sqrt (std::pow (fit->GetParameter (5) * fit->GetParError (5) / sigma, 2) + std::pow (fit->GetParameter (6) * fit->GetParError (6) / sigma, 2));

        const double tms = mean;
        const double tms_err = mean_err;
        
        const double tmr = sigma / mean;
        const double tmr_err = fabs (tmr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h2_avg_tms[iSys]->SetBinContent (iPtch+1, iFinerEta+1, tms * 100);
        h2_avg_tms[iSys]->SetBinError (iPtch+1, iFinerEta+1, tms_err * 100);
        h2_avg_tmr[iSys]->SetBinContent (iPtch+1, iFinerEta+1, tmr * 100);
        h2_avg_tmr[iSys]->SetBinError (iPtch+1, iFinerEta+1, tmr_err * 100);
      }

      for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {

        TF1* fit = new TF1 (Form ("f_tmr_%s_iPtch%i_iEta%i", sys.Data (), iPtch, iEta), &tmf, -1.0, 4.0, 9);

        TH1D* h = h_tmr_integratedEta[iEta];

        const double initRMS = h->GetRMS ();

        fit->SetParameter (0, -0.8);
        fit->SetParameter (1, 0.8);
        fit->SetParameter (2, 0);
        fit->SetParameter (3, 0.5 * 0.68 * h->Integral (h->FindBin (-initRMS), h->FindBin (initRMS)) / std::sqrt (2*M_PI*initRMS*initRMS));
        fit->SetParameter (4, 0.5 * 0.68 * h->Integral (h->FindBin (-initRMS), h->FindBin (initRMS)) / std::sqrt (2*M_PI*initRMS*initRMS));
        std::cout << "rms, int = " << initRMS << ", " << fit->GetParameter (3) << std::endl;
        fit->SetParameter (5, initRMS);
        fit->SetParameter (6, initRMS);
        fit->SetParameter (7, -2);
        fit->SetParameter (8, -2);

        fit->SetParLimits (0, -1., 0.0001);
        fit->SetParLimits (1, 0.0001, 4.);
        fit->SetParLimits (2, -0.5, 0.5);
        fit->SetParLimits (3, 0., 1e16);
        fit->SetParLimits (4, 0., 1e16);
        fit->SetParLimits (5, 0.0001, 2.);
        fit->SetParLimits (6, 0.0001, 2.);
        fit->SetParLimits (7, 0.0001, 10.);
        fit->SetParLimits (8, 0.0001, 10.);

        h->Fit (fit, "RN0");

        const double mean = fit->GetParameter (2);
        const double mean_err = fit->GetParError (1);
        const double sigma = std::sqrt (std::pow (fit->GetParameter (5), 2) + std::pow (fit->GetParameter (6), 2));
        const double sigma_err = std::sqrt (std::pow (fit->GetParameter (5) * fit->GetParError (5) / sigma, 2) + std::pow (fit->GetParameter (6) * fit->GetParError (6) / sigma, 2));

        const double tms = mean;
        const double tms_err = mean_err;
        
        const double tmr = sigma / mean;
        const double tmr_err = fabs (tmr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h_avg_tms[iSys][iEta]->SetBinContent (iPtch+1, tms * 100);
        h_avg_tms[iSys][iEta]->SetBinError (iPtch+1, tms_err * 100);
        h_avg_tmr[iSys][iEta]->SetBinContent (iPtch+1, tmr * 100);
        h_avg_tmr[iSys][iEta]->SetBinError (iPtch+1, tmr_err * 100);

        delete h_tmr_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_tmr_integratedEta;
    } // end loop over iPtch
  } // end loop over iSys


  for (int iSys : systems) {
    h2_avg_tms[iSys]->Write ();
    h2_avg_tmr[iSys]->Write ();

    for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {
      h_avg_tms[iSys][iEta]->Write ();
      h_avg_tmr[iSys][iEta]->Write ();
    } // end loop over iEta
  } // end loop over iSys

  outFile->Close ();
  
}

#endif
