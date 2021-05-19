#! /bin/bash

#declare -a systs=("None" "JetES5PercUpVar" "JetES5PercDownVar" "JetES2PercUpVar" "JetES2PercDownVar" "FcalCentVar" "FineFcalCentVar")
#declare -a systs=("None" "JetES5PercUpVar" "JetES5PercDownVar" "JetES2PercUpVar" "JetES2PercDownVar")
declare -a systs=("FineFcalCentVar")
declare -a trigs=("MinBias" "Jet50GeV")
#declare -a trigs=("Jet50GeV")

for syst in ${systs[@]}; do

  for trig in ${trigs[@]}; do
    condor_submit ProcessData.job trig=${trig} syst=${syst}
  done

done

