#! /bin/bash

#declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES2PercUpVar" "JetES2PercDownVar" "FcalCentVar" "FineFcalCentVar")
#declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES2PercUpVar" "JetES2PercDownVar" "FcalCentVar")
declare -a vardirs=("Nominal")
#declare -a minjpts=("15" "30" "60")
declare -a minjpts=("60")


for minjpt in ${minjpts[@]}; do

  sigdir=""
  if [[ minjpt -eq "30" ]] || [[ minjpt -eq "15" ]]; then
    sigdir="MinBiasData"
  elif [[ minjpt -eq "60" ]]; then
    sigdir="50GeVJetsData"
  elif [[ minjpt -eq "120" ]]; then
    sigdir="100GeVJetsData"
  fi

  for vardir in ${vardirs[@]}; do
    for cent in $(seq 0 4); do
    #for cent in $(seq 0 0); do
      condor_submit vardir=${vardir} minjpt=${minjpt} sigdir=${sigdir} cent=${cent} RunCorrelator.job
    done
  done

  #for cent in $(seq 0 99); do
  #for cent in $(seq 0 0); do
  #  condor_submit vardir=FineFcalCentVar minjpt=${minjpt} sigdir=${sigdir} cent=${cent} RunCorrelator.job
  #done

done
  
