#! /bin/bash

#declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES2PercUpVar" "JetES2PercDownVar")
declare -a vardirs=("MixCatVar1" "MixCatVar2" "MixCatVar3")
#declare -a minjpts=("15" "30" "60")
#declare -a minjpts=("30" "60")
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
      condor_submit vardir=${vardir} minjpt=${minjpt} sigdir=${sigdir} cent=${cent} RunCorrelator_pPb.job
    done
    condor_submit vardir=${vardir} minjpt=${minjpt} sigdir=${sigdir} RunCorrelator_pp.job
  done


  #for cent in $(seq 0 4); do
  #  condor_submit vardir=FcalCentVar minjpt=${minjpt} sigdir=${sigdir} cent=${cent} RunCorrelator_pPb.job
  #done


  #for cent in $(seq 0 99); do
  #  condor_submit vardir=FineFcalCentVar minjpt=${minjpt} sigdir=${sigdir} cent=${cent} RunCorrelator_pPb.job
  #done

done
  
