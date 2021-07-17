#! /bin/bash

#declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES5PercSmearVar" "JetES2PercUpVar" "JetES2PercDownVar" "JetES2PercSmearVar")
#declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES2PercUpVar" "JetES2PercDownVar")
declare -a vardirs=("Nominal")
#declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES2PercUpVar" "JetES2PercDownVar" "MixCatVar1" "MixCatVar2" "MixCatVar3")
#declare -a resdirs=("60GeVJets" "30GeVJets" "15GeVJets")
declare -a resdirs=("60GeVJets" "30GeVJets")
#declare -a resdirs=("30GeVJets")

for resdir in ${resdirs[@]}; do



  for vardir in ${vardirs[@]}; do

    histpath=rootFiles/Histograms/${resdir}/JetsHists/${vardir}

    hadd -f ${histpath}/data17_5TeV_hists.root \
            ${histpath}/340644_hists.root \
            ${histpath}/340683_hists.root \
            ${histpath}/340697_hists.root \
            ${histpath}/340718_hists.root \
            ${histpath}/340814_hists.root \
            ${histpath}/340849_hists.root \
            ${histpath}/340850_hists.root \
            ${histpath}/340910_hists.root \
            ${histpath}/340925_hists.root \
            ${histpath}/340973_hists.root \
            ${histpath}/341027_hists.root \
            ${histpath}/341123_hists.root \
            ${histpath}/341184_hists.root

    #for cent in $(seq 0 4); do   
 
    #  hadd -f ${histpath}/data16_5TeV_iCent${cent}_hists.root \
    #          ${histpath}/312796_iCent${cent}_hists.root \
    #          ${histpath}/312837_iCent${cent}_hists.root \
    #          ${histpath}/312937_iCent${cent}_hists.root \
    #          ${histpath}/312945_iCent${cent}_hists.root \
    #          ${histpath}/312968_iCent${cent}_hists.root
    #          #${histpath}/314199_iCent${cent}_hists.root

    #done

    histpath=rootFiles/Histograms/${resdir}/MixedHists/${vardir}

    hadd -f ${histpath}/data17_5TeV_hists.root \
            ${histpath}/340644_hists.root \
            ${histpath}/340683_hists.root \
            ${histpath}/340697_hists.root \
            ${histpath}/340718_hists.root \
            ${histpath}/340814_hists.root \
            ${histpath}/340849_hists.root \
            ${histpath}/340850_hists.root \
            ${histpath}/340910_hists.root \
            ${histpath}/340925_hists.root \
            ${histpath}/340973_hists.root \
            ${histpath}/341027_hists.root \
            ${histpath}/341123_hists.root \
            ${histpath}/341184_hists.root

    #for cent in $(seq 0 4); do   

    #  hadd -f ${histpath}/data16_5TeV_iCent${cent}_hists.root \
    #          ${histpath}/312796_iCent${cent}_hists.root \
    #          ${histpath}/312837_iCent${cent}_hists.root \
    #          ${histpath}/312937_iCent${cent}_hists.root \
    #          ${histpath}/312945_iCent${cent}_hists.root \
    #          ${histpath}/312968_iCent${cent}_hists.root
    #          #${histpath}/314199_iCent${cent}_hists.root

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
