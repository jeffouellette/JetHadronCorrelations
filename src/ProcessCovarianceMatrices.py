import sys
import os

import ROOT

import uproot
import pandas as pd
from scipy.linalg import block_diag
import numpy as np

def EmbedMatrix (A, B):
  C = B
  nb = B.shape[0]
  na = A.shape[0]
  lower = (nb) // 2 - (na // 2)
  upper = ((nb+1) // 2) + (na // 2)
  C[lower:upper, lower:upper] = A
  return C

def AddOverflow (cov):
  return EmbedMatrix (cov, np.zeros ((cov.shape[0]+2,cov.shape[1]+2)))


inFileName = os.getenv ("JETHADRONCORR_DATA_PATH") + "/rootFiles/Results/ProcessCorrelations_" + sys.argv[1] + ".root" # 1st argument is inFileTag
#inFile = TFile (inFileName)
inFile = uproot.open (inFileName)

outFileName = os.getenv ("JETHADRONCORR_DATA_PATH") + "/rootFiles/Results/ProcessCovarianceMatrices_" + sys.argv[2] + "/" # 2nd argument is outFileTag
#outFileName = os.getenv ("JETHADRONCORR_DATA_PATH") + "/rootFiles/Results/ProcessCovarianceMatrices_" + sys.argv[2] + ".root" # 2nd argument is outFileTag
##outFile = TFile (outFileName)
#outFile = uproot.recreate (outFileName)

# basic arrays...
dTypes = ["data", "mc"]
directions = ["ns", "perp", "as"] # angular integration regions; near-side, perpendicular, away-side
pTJBins = [15, 20, 30, 45, 60, 90, 120, 160, 200, 240, 300, 350, 400] # jet pT bins
pTChBins = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 50, 60, 75, 90, 120] # charged particle pT bins
centBins = ["iCent0", "iCent1", "iCent2", "iCent3", "iCent4", "allCent"]

numCovBins = (len(pTJBins)-1)*(len(pTChBins)-1)


print ("Starting loop...")
for iDType, dType in enumerate (dTypes):
  print ("On " + dType)

  for iDir, direction in enumerate (directions):
    print ("On " + direction)

    m_cov = []
    #m_cov.append (np.zeros ((len(pTChBins)+1,len(pTChBins)+1)))

    numCovBins = 0
    for iPtJ,pTJ in enumerate (pTJBins[:-1]):
      #print ("On" + str(pTJ))

      hname = "h2_jet_trk_pt_cov_" + direction + "_ref_sig_" + dType + "_" + str(pTJ) + "-" + str(pTJBins[iPtJ+1]) + "GeVJets"
      #m_cov.append (AddOverflow (inFile[hname].values ()))
      m_cov.append (inFile[hname].values ())

      if iDType + iDir + iPtJ == 0:
        df = pd.DataFrame (m_cov[-1])
        print (df)

      numCovBins += len(pTChBins) - 1 

    #done
    #m_cov.append (np.zeros ((len(pTChBins)+1,len(pTChBins)+1)))

    m_bigcov = block_diag (*m_cov)

    if iDType + iDir == 0:
      df = pd.DataFrame (m_bigcov)
      print (df)

    hname = "h2_jetAll_trk_pt_cov_" + direction + "_ref_sig_" + dType

    with open (outFileName + "/" + hname + ".txt", "w") as f:
      np.savetxt (f, m_bigcov)


    for iCent,cent in enumerate (centBins):

      m_cov = []

      #m_cov.append (np.zeros ((len(pTChBins)+1,len(pTChBins)+1)))

      numCovBins = 0
      for iPtJ,pTJ in enumerate (pTJBins[:-1]):
        #print ("On" + str(pTJ))

        hname = "h2_jet_trk_pt_cov_" + direction + "_pPb_sig_" + cent + "_" + dType + "_" + str(pTJ) + "-" + str(pTJBins[iPtJ+1]) + "GeVJets"
        #m_cov.append (AddOverflow (inFile[hname].values ()))
        m_cov.append (inFile[hname].values ())

        numCovBins += len(pTChBins) - 1 

      #done
      #m_cov.append (np.zeros ((len(pTChBins)+1,len(pTChBins)+1)))

      m_bigcov = block_diag (*m_cov)

      hname = "h2_jetAll_trk_pt_cov_" + direction + "_pPb_sig_" + cent + "_" + dType

      with open (outFileName + "/" + hname + ".txt", "w") as f:
        np.savetxt (f, m_bigcov)

    #done
  #done
#done
      



