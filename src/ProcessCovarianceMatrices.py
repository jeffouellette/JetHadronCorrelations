import sys
import os

from ROOT import TCanvas, TFile, TPaveText
from ROOT import gROOT, gBenchmark

from root_numpy import matrix

inFileName = os.getenv ("JETHADRONCORR_DATA_PATH") + "/rootFiles/Results/ProcessCorrelations_" + sys.argv[1] + ".root" # 1st argument is inFileTag
inFile = TFile (inFileName)

outFileName = ROOT.TString (sys.argv[2]) # 2nd argument is outFileTag
outFileName.ReplaceAll (".root", "")
outFileName = os.getenv ("JETHADRONCORR_DATA_PATH") + "/rootFiles/Results/ProcessCovarianceMatrices_" + outFileName + ".root"
outFile = TFile (outFileName)


dTypes = ["data", "mc"]
directions = ["ns", "perp", "as"]
pTJBins = [15, 20, 30, 45, 60, 90, 120, 160, 200, 240, 300, 350, 400]
pTChBins = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 50, 60, 75, 90, 120]

print ("Starting loop...")
for iDType, dType in enumerate (dTypes):
  print ("On " + dType)

  for iDir, direction in enumerate (directions):
    print ("On " + direction)

    h2_cov = []

    h2_bigcov = TH2D ("h2_jet_trk_pt_cov_" + direction + "_ref_sig_" + dType, "", len(pTJBins)*len(pTChBins), 0, len(pTJBins)*len(pTChBins), len(pTJBins)*len(pTChBins), 0, len(pTJBins)*len(pTChBins))

    for iPtJ,pTJ in enumerate (pTJBins[:-1]):

      h2_cov[iPtJ] = inFile.Get ("h2_jet_trk_pt_cov_" + direction + "_ref_sig_" + dType + "_" + pTJ + "-" + pTJBins[iPtJ+1] + "GeVJets")

      m_cov = TMatrixD (h2_cov[iPtJ].GetNbinsX ()+2, h2_cov[iPtJ].GetNbinsY ()+2, h2_cov[iPtJ].GetArray ());
      m_cov = matrix (m_cov)
      print (m_cov)
      m_cov = m_cov[1:-1][1:-1]
      print (m_cov)

    #done
  #done
#done
      



