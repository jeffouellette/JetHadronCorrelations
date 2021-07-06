#ifndef __AnalyzeTrackMomentumResolution_C__
#define __AnalyzeTrackMomentumResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>

#include <Utilities.h>
#include <ArrayTemplates.h>

#include "Params.h"
#include "LocalUtilities.h"
#include "TrackMomentumFit.h"

using namespace JetHadronCorrelations;

typedef TGraphAsymmErrors TGAE;


// Bin edge definitions: eta, pTch, pp/pPb.
const int nFinerEtaTrkBins = 40;
const double* finerEtaTrkBins = linspace (-2.5, 2.5, nFinerEtaTrkBins);

const double etaTrkBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const int nEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;

//const double pTchBins[] = {0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100};
const double pTchBins[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100};
const int nPtchBins = sizeof (pTchBins) / sizeof (pTchBins[0]) - 1;

//const vector <int> systems = {0, 1};
const vector <int> systems = {0}; // 0 = pp, 1 = pPb

TH1D**** h_tmr = nullptr;
TH2D** h2_avg_tms = nullptr;
TH2D** h2_avg_tmr = nullptr;
TH1D*** h_avg_tms = nullptr;
TH1D*** h_avg_tmr = nullptr;

TF1**** f_tmr_fits = nullptr;


// The LogTrackMomentumFit functor returns the log of a double-sided Crystal ball function with a double Gaussian core. It has 9 free parameters.
LogTrackMomentumFit ltmf;

// Minimum number of histogram entries to attempt a fit
const int minEntriesForFit = 1000;


/**
 * Wrapper function for performing a fit to a TF1 by taking to log of the entries in the histogram. (Sometimes better for fit stability.)
 * Note that the passed TF1 must be the log of the function you want to fit the histogram to. E.g. if the desired function is a gaussian, then the passed TF1 should be a quadratic.
 * Returns 0 if the fit was successful, 1 if there were not enough histogram entries, or 2 if the fit returned an error.
 */
int DoLogFit (TH1D* h, TF1* fit) {

  std::cout << "Fitting " << h->GetName () << " to " << fit->GetName () << std::endl;

  // Check we have "enough" entries to perform a fit. (1000 is an arbitrary choice)
  if (h->Integral () < minEntriesForFit) {
    std::cout << "Not enough counts for fit! Will skip this histogram." << std::endl;
    return 1;
  }

  // 2000 bins --> 1000 bins (dr = 0.005)
  int fStatus = -1;
  do {
    // Try reducing from 2000 bins
    if (fStatus != -1) {
      std::cout << "  --> Previous fit failed (fStatus = " << fStatus << ") --> rebinning by 2 and retrying." << std::endl;
      h->Rebin (2);
    }

    // Fit the log of the histogram (fit converges better this way)
    // Also drops points with 0 entries by storing non-zero points in a TGraphErrors
    TGraphErrors* g = new TGraphErrors ();
    for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
      if (h->GetBinContent (iX) > 0) {
        g->SetPoint (g->GetN (), h->GetBinCenter (iX), std::log (h->GetBinContent (iX)));
        g->SetPointError (g->GetN () - 1, 0.5 * h->GetBinWidth (iX), std::fabs (h->GetBinError (iX) / h->GetBinContent (iX)));
        //h->SetBinError (iX, h->GetBinError (iX) / h->GetBinContent (iX));
        //h->SetBinContent (iX, std::log (h->GetBinContent (iX)));
      }
    }

    const double initStdDev = h->GetStdDev ();
    const double integral = h->Integral (h->FindBin (-initStdDev), h->FindBin (initStdDev));

    ltmf.InitParams (fit, integral, initStdDev);

    //h->Fit (fit, "RN0QLL");
    //h->Fit (fit, "RN0Q");
    fStatus = (int) g->Fit (fit, "RN0"); // cannot use LL since graph is log(N)

    SaferDelete (&g);
  }
  while (fStatus != 0 && h->GetNbinsX () > 100);
  
  if (fStatus != 0) {
    std::cout << "  --> Fit failed! Please investigate further!" << std::endl;
    return 2;
  }
  else {
    std::cout << "  --> Fit succeeded!" << std::endl;
    return 0;
  }
}



/**
 * Fills an array with the track momentum scale (TMS) and its uncertainty, and the track momentum resolution (TMR) and its uncertainty, as extracted from a fitted TF1 response function.
 */
void GetFittedScaleAndResolution (TF1* f, double* res) {
  res[0] = (double) f->GetParameter (2);
  res[1] = (double) f->GetParError (2);

  const ldouble a1 = f->GetParameter (3);
  const ldouble a1e = std::fabs (f->GetParError (3));
  const ldouble a2 = f->GetParameter (4);
  const ldouble a2e = std::fabs (f->GetParError (4));
  const ldouble s1 = std::fabs (f->GetParameter (5));
  const ldouble s1e = std::fabs (f->GetParError (5));
  const ldouble s2 = std::fabs (f->GetParameter (6));
  const ldouble s2e = std::fabs (f->GetParError (6));

  std::cout << "  --> extracted fit parameters " << a1 << " +/- " << a1e << ", " << a2 << " +/- " << a2e << ", " << s1 << " +/- " << s1e << ", " << s2 << " +/- " << s2e << std::endl;

  // calculate the TMR from <x^2> where x is distributed according to the Gaussian core of f
  const ldouble norm = a1 * s1 + a2 * s2; // normalization factor
  const ldouble st = (norm != 0 ? std::sqrtl (std::fabs (a1*s1*s1*s1/norm + a2*s2*s2*s2/norm)) : 0);
  const ldouble ste = std::sqrtl (std::pow (a1e * a2 * s1*s2 * (s1*s1 - s2*s2), 2) +
                                  std::pow (a2e * a1 * s1*s2 * (s2*s2 - s1*s1), 2) + 
                                  std::pow (s1e * a1 * (2*a1*s1*s1*s1 + 3*a2*s1*s1*s2 - a2*s2*s2*s2), 2) + 
                                  std::pow (s2e * a2 * (2*a2*s2*s2*s2 + 3*a1*s1*s2*s2 - a1*s1*s1*s1), 2)) / (2*st * norm*norm);


  std::cout << "  --> fit st = " << st << " +/- " << ste << " (norm = " << norm << ")" << std::endl;
  
  res[2] = (double) st;
  res[3] = (double) ste; 
  return;
}


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

  // main results. Central values and widths of resolution functions in each bin.
  h2_avg_tms = Get1DArray <TH2D*> (2);
  h2_avg_tmr = Get1DArray <TH2D*> (2);

  h_avg_tms = Get2DArray <TH1D*> (2, nEtaTrkBins+1);
  h_avg_tmr = Get2DArray <TH1D*> (2, nEtaTrkBins+1);

  // eta-integrated histograms
  TH1D**** h_tmr_integratedEta = Get3DArray <TH1D*> (2, nPtchBins, nEtaTrkBins);
 
  // fit results
  f_tmr_fits = Get3DArray <TF1*> (2, nEtaTrkBins+1, nPtchBins);

  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");
    h2_avg_tms[iSys] = new TH2D (Form ("h2_avg_tms_%s", sys.Data ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMS [%]", nPtchBins, pTchBins, nFinerEtaTrkBins, finerEtaTrkBins);
    h2_avg_tmr[iSys] = new TH2D (Form ("h2_avg_tmr_%s", sys.Data ()), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMR [%]", nPtchBins, pTchBins, nFinerEtaTrkBins, finerEtaTrkBins);

    for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {
      h_avg_tms[iSys][iEta] = new TH1D (Form ("h_avg_tms_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T}^{truth} [GeV];TMS [%]", nPtchBins, pTchBins);
      h_avg_tmr[iSys][iEta] = new TH1D (Form ("h_avg_tmr_%s_iEta%i", sys.Data (), iEta), ";#it{p}_{T}^{truth} [GeV];TMR [%]", nPtchBins, pTchBins);
      for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {
        h_tmr_integratedEta[iSys][iPtch][iEta] = new TH1D (Form ("h_tmr_integratedEta_%s_iPtch%i_%s", sys.Data (), iPtch,  (iEta == nEtaTrkBins ? "allEta" : Form ("iEta%i", iEta))), "Track resolution, #it{p}_{T}^{truth} / #it{p}_{T}^{reco} - 1", 2000, -1.0, 4.0);
        h_tmr_integratedEta[iSys][iPtch][iEta]->Sumw2 ();
      } // end loop over iPtch
    } // end loop over iEta
  } // end loop over iSys

  // stores fitted TMS and TMR quantities
  double* fittedParams = new double[4];


  // Main loop over all histograms to do fits to resolution function
  for (int iSys : systems) {
    const TString sys = (iSys == 0 ? "pp" : "pPb");

    //for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {
    for (int iPtch = 0; iPtch < 2; iPtch++) {

      for (int iFinerEta = 0; iFinerEta < nFinerEtaTrkBins; iFinerEta++) {

        // First add to the relevant eta-integrated histograms
        const float binCenter = 0.5 * fabs (finerEtaTrkBins[iFinerEta] + finerEtaTrkBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < nEtaTrkBins) {
          if (etaTrkBins[iEta] < binCenter && binCenter < etaTrkBins[iEta+1])
            break;
          iEta++;
        }
        h_tmr_integratedEta[iSys][iPtch][iEta]->Add (h_tmr[iSys][iPtch][iFinerEta]);
        h_tmr_integratedEta[iSys][iPtch][nEtaTrkBins]->Add (h_tmr[iSys][iPtch][iFinerEta]);

        // Designate the histogram to fit
        TH1D* h = (TH1D*) h_tmr[iSys][iPtch][iFinerEta];

        // Define the TF1 fit function
        TF1* fit = new TF1 (Form ("f_tmr_%s_iPtch%i_iEta%i", sys.Data (), iPtch, iFinerEta), &ltmf, -1.0, 4.0, 9);

        // Perform the fit (defined in other function)
        const int fStatus = DoLogFit (h, fit);

        // Extract fit results if successful
        if (fStatus == 0) {
          GetFittedScaleAndResolution (fit, fittedParams);
          
          const double tms = (std::isfinite (fittedParams[0]) ? fittedParams[0] : 0);
          const double tms_err = (std::isfinite (fittedParams[1]) ? fittedParams[1] : 0);
          
          const double tmr = (std::isfinite (fittedParams[2]) ? fittedParams[2] : 0);
          const double tmr_err = (std::isfinite (fittedParams[3]) ? fittedParams[3] : 0);

          h2_avg_tms[iSys]->SetBinContent (iPtch+1, iFinerEta+1, tms * 100);
          h2_avg_tms[iSys]->SetBinError (iPtch+1, iFinerEta+1, tms_err * 100);
          h2_avg_tmr[iSys]->SetBinContent (iPtch+1, iFinerEta+1, tmr * 100);
          h2_avg_tmr[iSys]->SetBinError (iPtch+1, iFinerEta+1, tmr_err * 100);
        }

        // These fits are not stored (only eta-integrated ones are), though they could be in principle!
        SaferDelete (&fit);
      }

      for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {

        // Make a clone of the target histogram for fitting
        TH1D* h = h_tmr_integratedEta[iSys][iPtch][iEta];

        // Define the TF1 fit function
        TF1* fit = new TF1 (Form ("f_tmr_%s_iPtch%i_%s_integrated", sys.Data (), iPtch, (iEta == nEtaTrkBins ? "allEta" : Form ("iEta%i", iEta))), &ltmf, -1.0, 4.0, 9);

        // Perform the fit (defined in other function)
        const int fStatus = DoLogFit (h, fit);

        // Extract fit results if successful
        if (fStatus == 0) {
          GetFittedScaleAndResolution (fit, fittedParams);
          
          const double tms = (std::isfinite (fittedParams[0]) ? fittedParams[0] : 0);
          const double tms_err = (std::isfinite (fittedParams[1]) ? fittedParams[1] : 0);
          
          const double tmr = (std::isfinite (fittedParams[2]) ? fittedParams[2] : 0);
          const double tmr_err = (std::isfinite (fittedParams[3]) ? fittedParams[3] : 0);

          h_avg_tms[iSys][iEta]->SetBinContent (iPtch+1, tms * 100);
          h_avg_tms[iSys][iEta]->SetBinError (iPtch+1, tms_err * 100);
          h_avg_tmr[iSys][iEta]->SetBinContent (iPtch+1, tmr * 100);
          h_avg_tmr[iSys][iEta]->SetBinError (iPtch+1, tmr_err * 100);
        }

        // Store fit TF1 object
        f_tmr_fits[iSys][iEta][iPtch] = fit;
      } // end loop over iEta
    } // end loop over iPtch
  } // end loop over iSys


  for (int iSys : systems) {
    h2_avg_tms[iSys]->Write ();
    h2_avg_tmr[iSys]->Write ();

    for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {
      for (int iEta = 0; iEta < nFinerEtaTrkBins; iEta++) {
        h_tmr[iSys][iPtch][iEta]->Write ();
      } // end loop over iEta
    } // end loop over iPtch

    for (int iEta = 0; iEta <= nEtaTrkBins; iEta++) {
      h_avg_tms[iSys][iEta]->Write ();
      h_avg_tmr[iSys][iEta]->Write ();

      //for (int iPtch = 0; iPtch < nPtchBins; iPtch++) {
      for (int iPtch = 0; iPtch < 2; iPtch++) {
        h_tmr_integratedEta[iSys][iPtch][iEta]->Write ();
        f_tmr_fits[iSys][iEta][iPtch]->Write ();
      } // end loop over iPtch
    } // end loop over iEta
  } // end loop over iSys

  outFile->Close ();

}

#endif
