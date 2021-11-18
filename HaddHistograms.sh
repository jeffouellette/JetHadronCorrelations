#! /bin/bash

declare -a datavardirs=("Nominal")
#declare -a datavardirs=("Nominal" "HITightVar" "HILooseVar" "TrkEffVar" "FakeRateVar" "PrimFitVar" "JetPrimFracVar" "PartSpcVar" "MixCatVar1" "MixCatVar2" "MixCatVar3" "MixCatVar4" "MixCatVar5")
#declare -a datavardirs=("Nominal" "HITightVar" "HILooseVar" "TrkEffVar" "FakeRateVar" "PrimFitVar" "JetPrimFracVar" "PartSpcVar" "MixCatVar1" "MixCatVar3" "MixCatVar4" "MixCatVar5")
#declare -a datavardirs=("JetPrimFracVar")


#declare -a mcvardirs=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar14" "JESVar15" "JESVar16" "JESVar17" "JESVar18" "JESVar19" "JESVar20") # don't run JESVar16 or JESVar20, they are not applicable
declare -a mcvardirs=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar14" "JESVar15" "JESVar17" "JESVar18" "JESVar19")
#declare -a mcvardirs=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar14" "JESVar15" "JESVar17" "JESVar18" "JESVar19" "MCRecoJetsTruthParts" "MCTruthJetsTruthParts")
#declare -a mcvardirs=("Nominal" "MCTruthJetsTruthParts")
#declare -a mcvardirs=("MCRecoJetsTruthMatchedParts")

declare -a mcmixvardirs=("Nominal")
#declare -a mcmixvardirs=("Nominal" "MixCatVar2" "MixCatVar6")

declare -a sigmixdirs=("JetsHists" "MixedHists")
#declare -a sigmixdirs=("MixedHists")




for sigmix in ${sigmixdirs[@]}; do

  for vardir in ${datavardirs[@]}; do
  
    targetPath=rootFiles/Histograms/All/${sigmix}/${vardir}
    j50Path=rootFiles/Histograms/J50/${sigmix}/${vardir}
    mbPath=rootFiles/Histograms/MinBias/${sigmix}/${vardir}
  
    #hadd -f ${targetPath}/data17_5TeV_hists.root \
    #        ${j50Path}/340644_hists.root \
    #        ${j50Path}/340683_hists.root \
    #        ${j50Path}/340697_hists.root \
    #        ${j50Path}/340718_hists.root \
    #        ${j50Path}/340814_hists.root \
    #        ${j50Path}/340849_hists.root \
    #        ${j50Path}/340850_hists.root \
    #        ${j50Path}/340910_hists.root \
    #        ${j50Path}/340925_hists.root \
    #        ${j50Path}/340973_hists.root \
    #        ${j50Path}/341027_hists.root \
    #        ${j50Path}/341123_hists.root \
    #        ${j50Path}/341184_hists.root \
    #        ${mbPath}/340644_hists.root \
    #        ${mbPath}/340683_hists.root \
    #        ${mbPath}/340697_hists.root \
    #        ${mbPath}/340718_hists.root \
    #        ${mbPath}/340814_hists.root \
    #        ${mbPath}/340849_hists.root \
    #        ${mbPath}/340850_hists.root \
    #        ${mbPath}/340910_hists.root \
    #        ${mbPath}/340925_hists.root \
    #        ${mbPath}/340973_hists.root \
    #        ${mbPath}/341027_hists.root \
    #        ${mbPath}/341123_hists.root \
    #        ${mbPath}/341184_hists.root &

    #wait
  
    for cent in $(seq 0 4); do   
  
      hadd -f ${targetPath}/data16_5TeV_iCent${cent}_hists.root \
              ${j50Path}/312796*iCent${cent}_hists.root \
              ${j50Path}/312837*iCent${cent}_hists.root \
              ${j50Path}/312937*iCent${cent}_hists.root \
              ${j50Path}/312945*iCent${cent}_hists.root \
              ${j50Path}/312968*iCent${cent}_hists.root \
              ${j50Path}/314199*iCent${cent}_hists.root \
              ${mbPath}/312796*iCent${cent}_hists.root \
              ${mbPath}/312837*iCent${cent}_hists.root \
              ${mbPath}/312937*iCent${cent}_hists.root \
              ${mbPath}/312945*iCent${cent}_hists.root \
              ${mbPath}/312968*iCent${cent}_hists.root \
              ${mbPath}/314199*iCent${cent}_hists.root
  
    done
  
    hadd -f ${targetPath}/data16_5TeV_allCent_hists.root \
            ${targetPath}/data16_5TeV_iCent0_hists.root \
            ${targetPath}/data16_5TeV_iCent1_hists.root \
            ${targetPath}/data16_5TeV_iCent2_hists.root \
            ${targetPath}/data16_5TeV_iCent3_hists.root \
            ${targetPath}/data16_5TeV_iCent4_hists.root &

  done

done
  
  
  
  
#for vardir in ${mcvardirs[@]}; do
#
#  targetPath=rootFiles/Histograms/All/JetsHists/${vardir}
#  mcPath=rootFiles/Histograms/MC/JetsHists/${vardir}
#
#  hadd -f ${targetPath}/mc17_5TeV_hists.root \
#          ${mcPath}/pp17_JZ0*hists.root \
#          ${mcPath}/pp17_JZ1_a*hists.root \
#          ${mcPath}/pp17_JZ1_b*hists.root \
#          ${mcPath}/pp17_JZ2_a*hists.root \
#          ${mcPath}/pp17_JZ2_b*hists.root \
#          ${mcPath}/pp17_JZ3_a*hists.root \
#          ${mcPath}/pp17_JZ3_b*hists.root &
#
#  for cent in $(seq 0 4); do   
#
#    hadd -f ${targetPath}/mc16_5TeV_iCent${cent}_hists.root \
#            ${mcPath}/pPb16_JZ1*iCent${cent}_hists.root \
#            ${mcPath}/pPb16_JZ2*iCent${cent}_hists.root \
#            ${mcPath}/pPb16_JZ3*iCent${cent}_hists.root &
#
#  done
#
#  wait
#
#  hadd -f ${targetPath}/mc16_5TeV_allCent_hists.root \
#          ${targetPath}/mc16_5TeV_iCent0_hists.root \
#          ${targetPath}/mc16_5TeV_iCent1_hists.root \
#          ${targetPath}/mc16_5TeV_iCent2_hists.root \
#          ${targetPath}/mc16_5TeV_iCent3_hists.root \
#          ${targetPath}/mc16_5TeV_iCent4_hists.root & 
#
#done




#for vardir in ${mcmixvardirs[@]}; do
#
#  targetPath=rootFiles/Histograms/All/MixedHists/${vardir}
#  mcPath=rootFiles/Histograms/MC/MixedHists/${vardir}
#  
#  hadd -f ${targetPath}/mc17_5TeV_hists.root \
#          ${mcPath}/pp17_JZ0*hists.root \
#          ${mcPath}/pp17_JZ1_a*hists.root \
#          ${mcPath}/pp17_JZ1_b*hists.root \
#          ${mcPath}/pp17_JZ2_a*hists.root \
#          ${mcPath}/pp17_JZ2_b*hists.root \
#          ${mcPath}/pp17_JZ3_a*hists.root \
#          ${mcPath}/pp17_JZ3_b*hists.root &
#  
#  for cent in $(seq 0 4); do   
#  
#    hadd -f ${targetPath}/mc16_5TeV_iCent${cent}_hists.root \
#            ${mcPath}/pPb16_JZ1*iCent${cent}_hists.root \
#            ${mcPath}/pPb16_JZ2*iCent${cent}_hists.root \
#            ${mcPath}/pPb16_JZ3*iCent${cent}_hists.root &
#  
#  done
#  
#  wait
#  
#  hadd -f ${targetPath}/mc16_5TeV_allCent_hists.root \
#          ${targetPath}/mc16_5TeV_iCent0_hists.root \
#          ${targetPath}/mc16_5TeV_iCent1_hists.root \
#          ${targetPath}/mc16_5TeV_iCent2_hists.root \
#          ${targetPath}/mc16_5TeV_iCent3_hists.root \
#          ${targetPath}/mc16_5TeV_iCent4_hists.root & 
#
#done




#histpath=rootFiles/Histograms/${resdir}/JetsHists/FineFcalCentVar

#for cent in $(seq 0 99); do

#  hadd -f ${histpath}/data16_5TeV_iCent${cent}_hists.root \
#          ${histpath}/312796_iCent${cent}_hists.root \
#          ${histpath}/312837_iCent${cent}_hists.root \
#          ${histpath}/312937_iCent${cent}_hists.root \
#          ${histpath}/312945_iCent${cent}_hists.root \
#          ${histpath}/312968_iCent${cent}_hists.root \
#          ${histpath}/314199_iCent${cent}_hists.root

#done

#histpath=rootFiles/Histograms/${resdir}/MixedHists/FineFcalCentVar

#for cent in $(seq 0 99); do

#  hadd -f ${histpath}/data16_5TeV_iCent${cent}_hists.root \
#          ${histpath}/312796_iCent${cent}_hists.root \
#          ${histpath}/312837_iCent${cent}_hists.root \
#          ${histpath}/312937_iCent${cent}_hists.root \
#          ${histpath}/312945_iCent${cent}_hists.root \
#          ${histpath}/312968_iCent${cent}_hists.root \
#          ${histpath}/314199_iCent${cent}_hists.root
#done



wait
