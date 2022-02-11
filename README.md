# JetHadronCorrelations
Analysis of jet-hadron correlations in 5.02 TeV p+Pb and pp with ATLAS run 2 data.

Notation: files ending in **.cxx** are meant to be compiled (see Makefile). Files ending in **.C** are ROOT macros.


## Before running the correlations code ##
There are several primary analysis routines for analyzing TTrees for auxiliary studies:

1. **TrackingPerformance.cxx** analyzes the tracking efficiency and primary fraction in MC, and outputs histograms per MC config to rootFiles/TrackingPerformance. Files must be hadd'ed to outFile.root, which further analyzed by PlotTrackingPerformance.cxx, which both creates plots as well as summary histograms and functions for applying these corrections, which are stored in aux/TrackingPerformance.root.

2. **TrackMomentumResolution.cxx** analyzes the track momentum resolution by comparing reco. and truth level transverse momenta in MC. Output is stored in rootFiles/TrackMomentumResolution. AnalyzeTrackMomentumResolution.C takes in the raw response distributions and performs fits to quantify the resolution and scale. Finally, PlotTrackMomentumResolution.C outputs summary plots of this study.

3. **JetEnergyResolution.cxx** analyzes the jet energy scale and resolution. Output is stored in rootFiles/JetEnergyResolution. AnalyzeJetEnergyResolution.C takes in the raw response distributions and performs the recursive Gaussian fits to quantify the JES and JER. Finally, PlotJetEnergyResolution.C outputs summary plots of this study.

4. **CentralityAnalysis.cxx** analyzes the aggregate centrality variables used in this analysis, including FCal total Et and ZDC energies. Output is stored in rootFiles/CentralityAnalysis. PlotCentralityAnalysis.C creates 1 dimensional plots of centrality variables, makes data/MC ratios for FCal reweighting, and outputs centrality category cuts. PlotCentralityAnalysis2D.C makes nice 2 dimensional plots of these variables, e.g. the FCal-ZDC energy correlation plot.


The primary TTree-to-TH1 routine which analyzes jet production and correlates jets with charged particles is **RunCorrelator.cxx**. This is designed to run over HTCondor, submission files can be found in run/RunCorrelator\*. Systematic flags, trigger flags, and data set configuration are all set via flags that are converted to instances of enum classes which make the code highly configurable. Output from RunCorrelator.cxx is stored in rootFiles/Histograms/(MC or trigger name). Output must be hadd'ed appropriately, which can be handled via HaddHistograms.sh. Note that this script should be modified to only hadd the configurations you want to add -- otherwise it will add everything together which can take a while...


## Producing inputs to the unfolding ##
To produce the response matrices, you must already have the per-jet weights as described below. The main code that produces the response matrices can be found in **MakeResponseMatrix.cxx**. It outputs to rootFiles/MakeResponseMatrices. Files must be hadd'ed together as usual before further processing.

Code that studies the post-unfolding closure of the analysis, as well as plots the response matrices, can be found under **PlotResponseMatrix.C**. This macro takes a boolean argument which, when true, studies the closure only when primary truth-matched reco. particles are considered. (This eliminates any possible effect entering from the primary fraction correction and/or the background subtraction.)  When false, it performs the nominal analysis on MC, including the background subtraction. Output plots are stored under Plots/Unfolding. Auxiliary information storing non-closure in the MC (for systematic variation purposes) is stored under aux/MCClosureHists.root.


## Analyzing and unfolding the results ##
The main analysis then proceeds in several modular steps, commands for which are stored in **ProcessHists.sh**.

1. First, **ProcessJets.C** and then **PlotJets.C** should be run to produce jet spectra and -- more importantly -- per-jet weights for MC events. The resultant weights are stored in aux/JetPtWeights.root. After deriving the jet weights, the MC can be fully processed with RunCorrelator.cxx. (It should have been previously processed to derive the weights, but should be reprocessed to incorporate these weights.)

2. Next, the macro **ProcessCorrelations.cxx** reads in all of the raw histograms with correlations data and (1) combines all systematic variations and collision configurations, (2) calculates uncertainties on data using covariance matrices, and (3) outputs the above into one file under rootFiles/Results/ProcessCorrelations\*.

3. After that, the jet-by-jet covariance matrices in charged particle transverse momentum are combined into one giga-matrix that is input to the unfolding using ProcessCovarianceMatrices.py. The results are stored in rootFiles/Results/ProcessCovarianceMatrices\* as .dat files and can be read into 2D histograms via LocalUtilities.cxx / GetCovarianceMatrix.

4. The unfolding convergence studies are performed by **ProcessNIters.C**. This macro takes as input the correlation histograms from step (2) and the covariance matrices from step (3) and performs the unfolding at a wide range of iterations values. **Warning**: the 2D unfolding takes an appreciable amount of time, and iterations beyond 20 are not recommended. A smaller number of iterations is suggested for testing purposes. The results from this study are all stored in rootFiles/Results/ProcessNIters\*. They can be analyzed and plotted using **PlotNIters.C**.

5. Finally, the main unfolding is performed on all the data (including systematic variations) in **ProcessUnfolding.cxx**. There are a few functions at the top of the source code which fix the number of iterations to be used in various histograms. These are just switch statements and can be easily edited. The output from this code is stored in rootFiles/Results/ProcessUnfolding\*.

6. After all of the above has been run, the main results can be studied. The primary plotting macro is **PlotPtCh.C**. This outputs the histograms in rootFiles/Results/ProcessCorrelations\* and rootFiles/Results/ProcessUnfolding\* (as well as auxiliary theory and data comparisons in aux/) to Plots/PtCh and Plots/Systematics/PtCh.



