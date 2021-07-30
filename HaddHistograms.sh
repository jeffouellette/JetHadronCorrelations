#! /bin/bash

declare -a datavardirs=("Nominal" "TrkEffVar" "FakeRateVar" "PrimFitVar" "PartSpcVar" "MixCatVar1" "MixCatVar2" "MixCatVar3")
#declare -a datavardirs=("Nominal" "HITightVar" "HILooseVar" "TrkEffVar" "FakeRateVar" "PrimFitVar" "PartSpcVar" "MixCatVar1" "MixCatVar2" "MixCatVar3")
#declare -a mcvardirs=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar14" "JESVar15" "JESVar16" "JESVar17" "JESVar18" "JESVar19" "JESVar20") # don't run JESVar16 or JESVar20, they are not applicable
declare -a mcvardirs=("Nominal" "JESVar0" "JESVar1" "JESVar2" "JESVar3" "JESVar4" "JESVar5" "JESVar6" "JESVar7" "JESVar8" "JESVar9" "JESVar10" "JESVar11" "JESVar12" "JESVar13" "JESVar14" "JESVar15" "JESVar17" "JESVar18" "JESVar19")
#declare -a minjpts=("30" "60")
declare -a minjpts=("60")

for minjpt in ${minjpts[@]}; do

  resdir="${minjpt}GeVJets"

  jzslice=""
  if [[ minjpt -eq "30" ]]; then
    jzslice="JZ1"
  elif [[ minjpt -eq "60" ]]; then
    jzslice="JZ2"
  fi
    

  #for vardir in ${datavardirs[@]}; do

  #  histpath=rootFiles/Histograms/${resdir}/JetsHists/${vardir}

  #  hadd -f ${histpath}/data17_5TeV_hists.root \
  #          ${histpath}/340644_hists.root \
  #          ${histpath}/340683_hists.root \
  #          ${histpath}/340697_hists.root \
  #          ${histpath}/340718_hists.root \
  #          ${histpath}/340814_hists.root \
  #          ${histpath}/340849_hists.root \
  #          ${histpath}/340850_hists.root \
  #          ${histpath}/340910_hists.root \
  #          ${histpath}/340925_hists.root \
  #          ${histpath}/340973_hists.root \
  #          ${histpath}/341027_hists.root \
  #          ${histpath}/341123_hists.root \
  #          ${histpath}/341184_hists.root

  #  for cent in $(seq 0 4); do   
 
  #    hadd -f ${histpath}/data16_5TeV_iCent${cent}_hists.root \
  #            ${histpath}/312796_iCent${cent}_hists.root \
  #            ${histpath}/312837_iCent${cent}_hists.root \
  #            ${histpath}/312937_iCent${cent}_hists.root \
  #            ${histpath}/312945_iCent${cent}_hists.root \
  #            ${histpath}/312968_iCent${cent}_hists.root \
  #            ${histpath}/314199_iCent${cent}_hists.root

  #  done


  #  histpath=rootFiles/Histograms/${resdir}/MixedHists/${vardir}

  #  hadd -f ${histpath}/data17_5TeV_hists.root \
  #          ${histpath}/340644_hists.root \
  #          ${histpath}/340683_hists.root \
  #          ${histpath}/340697_hists.root \
  #          ${histpath}/340718_hists.root \
  #          ${histpath}/340814_hists.root \
  #          ${histpath}/340849_hists.root \
  #          ${histpath}/340850_hists.root \
  #          ${histpath}/340910_hists.root \
  #          ${histpath}/340925_hists.root \
  #          ${histpath}/340973_hists.root \
  #          ${histpath}/341027_hists.root \
  #          ${histpath}/341123_hists.root \
  #          ${histpath}/341184_hists.root

  #  for cent in $(seq 0 4); do   

  #    hadd -f ${histpath}/data16_5TeV_iCent${cent}_hists.root \
  #            ${histpath}/312796_iCent${cent}_hists.root \
  #            ${histpath}/312837_iCent${cent}_hists.root \
  #            ${histpath}/312937_iCent${cent}_hists.root \
  #            ${histpath}/312945_iCent${cent}_hists.root \
  #            ${histpath}/312968_iCent${cent}_hists.root \
  #            ${histpath}/314199_iCent${cent}_hists.root

  #  done

  #done


  for vardir in ${mcvardirs[@]}; do

    histpath=rootFiles/Histograms/${resdir}/JetsHists/${vardir}

    hadd -f ${histpath}/mc17_5TeV_hists.root \
            ${histpath}/pp17_JZ1_a_hists.root \
            ${histpath}/pp17_JZ1_b_hists.root \
            ${histpath}/pp17_JZ2_a_hists.root \
            ${histpath}/pp17_JZ2_b_hists.root \
            ${histpath}/pp17_JZ3_a_hists.root \
            ${histpath}/pp17_JZ3_b_hists.root
            #${histpath}/pp17_${jzslice}_a_hists.root \
            #${histpath}/pp17_${jzslice}_b_hists.root

    for cent in $(seq 0 4); do   
 
      hadd -f ${histpath}/mc16_5TeV_iCent${cent}_hists.root \
              ${histpath}/pPb16_JZ1_iCent${cent}_hists.root \
              ${histpath}/pPb16_JZ2_iCent${cent}_hists.root \
              ${histpath}/pPb16_JZ3_iCent${cent}_hists.root
              #${histpath}/pPb16_${jzslice}_iCent${cent}_hists.root

    done

    #histpath=rootFiles/Histograms/${resdir}/MixedHists/${vardir}

    #hadd -f ${histpath}/mc17_5TeV_hists.root \
    #        ${histpath}/pp17_${jzslice}_a_hists.root \
    #        ${histpath}/pp17_${jzslice}_b_hists.root

    #for cent in $(seq 0 4); do   
 
    #  hadd -f ${histpath}/mc16_5TeV_iCent${cent}_hists.root \
    #          ${histpath}/pPb16_${jzslice}_iCent${cent}_hists.root

    #done

  done



  #histpath=rootFiles/Histograms/${resdir}/JetsHists/FcalCentVar

  #for cent in $(seq 0 4); do

  #  hadd -f ${histpath}/data16_5TeV_iCent${cent}_hists.root \
  #          ${histpath}/312796_iCent${cent}_hists.root \
  #          ${histpath}/312837_iCent${cent}_hists.root \
  #          ${histpath}/312937_iCent${cent}_hists.root \
  #          ${histpath}/312945_iCent${cent}_hists.root \
  #          ${histpath}/312968_iCent${cent}_hists.root \
  #          ${histpath}/314199_iCent${cent}_hists.root

  #done

  #histpath=rootFiles/Histograms/${resdir}/MixedHists/FcalCentVar

  #for cent in $(seq 0 4); do

  #  hadd -f ${histpath}/data16_5TeV_iCent${cent}_hists.root \
  #          ${histpath}/312796_iCent${cent}_hists.root \
  #          ${histpath}/312837_iCent${cent}_hists.root \
  #          ${histpath}/312937_iCent${cent}_hists.root \
  #          ${histpath}/312945_iCent${cent}_hists.root \
  #          ${histpath}/312968_iCent${cent}_hists.root \
  #          ${histpath}/314199_iCent${cent}_hists.root
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

done
