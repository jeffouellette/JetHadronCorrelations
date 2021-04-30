#! /bin/bash

declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES5PercSmearVar" "JetES2PercUpVar" "JetES2PercDownVar" "JetES2PercSmearVar")
declare -a minjpts=("30" "60")

for minjpt in ${minjpts[@]}; do

  sigdir=""
  if [[ minjpt -eq "30" ]]; then
    sigdir="MinBiasData"
  elif [[ minjpt -eq "60" ]]; then
    sigdir="50GeVJetsData"
  fi

  for vardir in ${vardirs[@]}; do

    condor_submit vardir=${vardir} minjpt=${minjpt} sigdir=${sigdir} RunCorrelator.job

  done

done
  
