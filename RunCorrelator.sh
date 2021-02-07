#! /bin/bash

root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "jets_data16", "./rootFiles/JetsData/Nominal/data16_5TeV_hists.root", "./rootFiles/JetsData/Nominal/data16_5TeV.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340644_hists.root", "./rootFiles/JetsData/Nominal/340644.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340683_hists.root", "./rootFiles/JetsData/Nominal/340683.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340697_hists.root", "./rootFiles/JetsData/Nominal/340697.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340718_hists.root", "./rootFiles/JetsData/Nominal/340718.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340814_hists.root", "./rootFiles/JetsData/Nominal/340814.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340849_hists.root", "./rootFiles/JetsData/Nominal/340849.root")' &
wait

root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340850_hists.root", "./rootFiles/JetsData/Nominal/340850.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340910_hists.root", "./rootFiles/JetsData/Nominal/340910.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340925_hists.root", "./rootFiles/JetsData/Nominal/340925.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/340973_hists.root", "./rootFiles/JetsData/Nominal/340973.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/341027_hists.root", "./rootFiles/JetsData/Nominal/341027.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/341123_hists.root", "./rootFiles/JetsData/Nominal/341123.root")' &
root -b -q 'src/RunCorrelator.C("pp17_5TeV", "jets_data17", "./rootFiles/JetsData/Nominal/341184_hists.root", "./rootFiles/JetsData/Nominal/341184.root")' &
wait


#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent0_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent0.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent0.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent1_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent1.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent1.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent2_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent2.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent2.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent3_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent3.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent3.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent4_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent4.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent4.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent5_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent5.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent5.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent6_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent6.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent6.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent7_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent7.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent7.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent8_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent8.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent8.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent9_hists.root", "./rootFiles/JetsData/Nominal/312649/312649_iCent9.root", "./rootFiles/MinBiasData/Nominal/312649/312649_iCent9.root")' &
#wait


#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent0_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent0.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent0.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent1_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent1.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent1.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent2_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent2.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent2.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent3_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent3.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent3.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent4_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent4.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent4.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent5_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent5.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent5.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent6_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent6.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent6.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent7_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent7.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent7.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent8_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent8.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent8.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent9_hists.root", "./rootFiles/JetsData/Nominal/312796/312796_iCent9.root", "./rootFiles/MinBiasData/Nominal/312796/312796_iCent9.root")' &
#wait


#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent0_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent0.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent0.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent1_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent1.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent1.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent2_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent2.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent2.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent3_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent3.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent3.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent4_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent4.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent4.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent5_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent5.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent5.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent6_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent6.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent6.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent7_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent7.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent7.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent8_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent8.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent8.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent9_hists.root", "./rootFiles/JetsData/Nominal/312837/312837_iCent9.root", "./rootFiles/MinBiasData/Nominal/312837/312837_iCent9.root")' &
#wait


#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent0_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent0.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent0.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent1_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent1.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent1.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent2_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent2.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent2.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent3_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent3.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent3.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent4_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent4.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent4.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent5_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent5.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent5.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent6_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent6.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent6.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent7_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent7.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent7.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent8_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent8.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent8.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent9_hists.root", "./rootFiles/JetsData/Nominal/312937/312937_iCent9.root", "./rootFiles/MinBiasData/Nominal/312937/312937_iCent9.root")' &
wait


#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent0_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent0.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent0.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent1_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent1.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent1.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent2_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent2.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent2.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent3_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent3.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent3.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent4_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent4.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent4.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent5_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent5.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent5.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent6_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent6.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent6.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent7_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent7.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent7.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent8_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent8.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent8.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent9_hists.root", "./rootFiles/JetsData/Nominal/312945/312945_iCent9.root", "./rootFiles/MinBiasData/Nominal/312945/312945_iCent9.root")' &
#wait


#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent0_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent0.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent0.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent1_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent1.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent1.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent2_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent2.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent2.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent3_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent3.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent3.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent4_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent4.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent4.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent5_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent5.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent5.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent6_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent6.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent6.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent7_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent7.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent7.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent8_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent8.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent8.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent9_hists.root", "./rootFiles/JetsData/Nominal/312968/312968_iCent9.root", "./rootFiles/MinBiasData/Nominal/312968/312968_iCent9.root")' &
#wait


#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent0_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent0.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent0.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent1_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent1.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent1.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent2_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent2.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent2.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent3_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent3.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent3.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent4_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent4.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent4.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent5_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent5.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent5.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent6_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent6.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent6.root")' &
#root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent7_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent7.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent7.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent8_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent8.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent8.root")' &
root -b -q 'src/RunCorrelator.C("pPb16_5TeV", "minbias_data16", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent9_hists.root", "./rootFiles/JetsData/Nominal/314199/314199_iCent9.root", "./rootFiles/MinBiasData/Nominal/314199/314199_iCent9.root")' &
wait


#root -b -q 'src/Plotter.C+("./rootFiles/MinBiasData/Nominal/data16_5TeV_hists.root")'
