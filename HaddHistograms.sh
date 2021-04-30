#! /bin/bash

#declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES5PercSmearVar" "JetES2PercUpVar" "JetES2PercDownVar" "JetES2PercSmearVar" "FcalCentVar" "FineFcalCentVar")
declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES5PercSmearVar" "JetES2PercUpVar" "JetES2PercDownVar" "JetES2PercSmearVar")
declare -a resdirs=("60GeVJets" "30GeVJets")

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
    
    hadd -f ${histpath}/data16_5TeV_hists.root \
            ${histpath}/312796_*_hists.root \
            ${histpath}/312837_*_hists.root \
            ${histpath}/312937_*_hists.root \
            ${histpath}/312945_*_hists.root \
            ${histpath}/312968_*_hists.root \
            ${histpath}/314199_*_hists.root
 
    histpath=rootFiles/Histograms/${resdir}/MixedHists/${vardir}

    hadd -f ${histpath}/data16_5TeV_hists.root \
            ${histpath}/312796_*_hists.root \
            ${histpath}/312837_*_hists.root \
            ${histpath}/312937_*_hists.root \
            ${histpath}/312945_*_hists.root \
            ${histpath}/312968_*_hists.root \
            ${histpath}/314199_*_hists.root

  done

done
