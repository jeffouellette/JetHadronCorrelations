#! /bin/bash

declare -a systs=("TrkEffVar" "PrimFitVar")
#declare -a systs=("Nominal" "HITightVar" "HILooseVar" "TrkEffVar" "FakeRateVar" "PrimFitVar" "PartSpcUnc" "MixCatVar1" "MixCatVar2" "MixCatVar3")
#declare -a systs=("MixCatVar1" "MixCatVar2" "MixCatVar3")
declare -a jessysts=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar14" "JESVar15" "JESVar17" "JESVar18" "JESVar19")
#declare -a minjpts=("15" "30" "60")
declare -a minjpts=("30" "60")
#declare -a minjpts=("30")


for minjpt in ${minjpts[@]}; do

  trig=""
  ppjetstream=""
  pptracksstream="MinBias"
  pPbjetstream=""
  pPbtracksstream="MinBias"
  if [[ minjpt -eq "30" ]] || [[ minjpt -eq "15" ]]; then
    trig="MinBias"
    ppjetstream="MinBias"
    pPbjetstream="MinBias"
  elif [[ minjpt -eq "60" ]]; then
    trig="Jet50GeV"
    ppjetstream="Main"
    pPbjetstream="Main"
  elif [[ minjpt -eq "120" ]]; then
    trig="Jet100GeV"
    ppjetstream="Main"
    pPbjetstream="Main"
  fi


  for syst in ${systs[@]}; do
    #condor_submit syst=${syst} minjpt=${minjpt} trig=${trig} jetstream=${pPbjetstream} tracksstream=${pPbtracksstream} RunCorrelator_pPb.job
    condor_submit syst=${syst} minjpt=${minjpt} trig=${trig} jetstream=${ppjetstream} tracksstream=${pptracksstream} RunCorrelator_pp.job
  done

  #for jessyst in ${jessysts[@]}; do
  #  condor_submit syst=${jessyst} minjpt=${minjpt} RunCorrelator_MC.job
  #done
  #condor_submit syst=Nominal minjpt=${minjpt} RunCorrelator_MC_Mixing.job


  #condor_submit syst=FcalCentVar minjpt=${minjpt} trig=${trig} jetstream=${pPbjetstream} tracksstream=${pPbtracksstream} RunCorrelator_pPb.job

  #condor_submit syst=FineFcalCentVar minjpt=${minjpt} trig=${trig} jetstream=${pPbjetstream} tracksstream=${pPbtracksstream} RunCorrelator_pPb.job

done
  
