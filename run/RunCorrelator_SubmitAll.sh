#! /bin/bash

#declare -a systs=("Nominal")
#declare -a systs=("Nominal" "HITightVar" "HILooseVar" "TrkEffVar" "FakeRateVar" "PrimFitVar" "JetPrimFracVar" "PartSpcVar" "MixCatVar1" "MixCatVar2" "MixCatVar3" "MixCatVar4" "MixCatVar5")
declare -a systs=("HITightVar" "HILooseVar" "TrkEffVar" "FakeRateVar" "PrimFitVar" "JetPrimFracVar" "PartSpcVar" "MixCatVar1" "MixCatVar2" "MixCatVar3" "MixCatVar4" "MixCatVar5")
#declare -a systs=("Nominal" "HITightVar" "HILooseVar" "TrkEffVar" "FakeRateVar" "PrimFitVar" "JetPrimFracVar" "PartSpcVar" "MixCatVar1" "MixCatVar2" "MixCatVar3" "MixCatVar4" "MixCatVar5" "FcalCentVar" "FineFcalCentVar")

#declare -a mcsysts=("JESVar18 JESVar19")
#declare -a mcsysts=("Nominal" "MCRecoJetsTruthParts" "MCTruthJetsTruthParts")
#declare -a mcsysts=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar14" "JESVar15" "JESVar17" "JESVar18" "JESVar19" "MCRecoJetsTruthParts" "MCTruthJetsTruthParts")
declare -a mcsysts=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar15" "JESVar16" "JESVar17" "JESVar18" "JESVar19")

#declare -a mcmixsysts=("Nominal")
declare -a mcmixsysts=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar15" "JESVar16" "JESVar17" "JESVar18" "JESVar19")

declare -a trigs=("J50" "MinBias")
#declare -a trigs=("J50")


for trig in ${trigs[@]}; do

  for syst in ${systs[@]}; do
    condor_submit syst=${syst} trig=${trig} RunCorrelator_pPb.job
  done

done

for syst in ${systs[@]}; do
  condor_submit syst=${syst} RunCorrelator_pp.job
done



for syst in ${mcsysts[@]}; do
  condor_submit syst=${syst} RunCorrelator_MC.job
done

for syst in ${mcmixsysts[@]}; do
  condor_submit syst=${syst} RunCorrelator_MC_Mixing.job
done

