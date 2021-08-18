
#root -l -b -q 'src/ProcessCorrelations.C("60GeVJets", "60GeVJets")'
#root -l -b -q 'src/ProcessCorrelations.C("30GeVJets", "30GeVJets")'

#root -l -b -q 'src/ProcessFineFCalPtCh.C("60GeVJets", "60GeVJets")'
#root -l -b -q 'src/ProcessFineFCalPtCh.C("30GeVJets", "30GeVJets")'
#root -l -b -q 'src/CombineFineFCalPtCh.C("60GeVJets")'
#root -l -b -q 'src/CombineFineFCalPtCh.C("30GeVJets")'

#root -l -b -q 'src/ProcessJets.C("60GeVJets", "60GeVJets")'
#root -l -b -q 'src/ProcessJets.C("30GeVJets", "30GeVJets")'
  
#root -l -b -q 'src/PlotFineFCalPtCh.C("60GeVJets", "60GeVJets", "60GeVJets")'
#root -l -b -q 'src/PlotFineFCalPtCh.C("30GeVJets", "30GeVJets", "30GeVJets")'

#root -l -b -q 'src/PlotDPhi.C("60GeVJets", "60GeVJets")'
#root -l -b -q 'src/PlotDPhi.C("30GeVJets", "30GeVJets")'

#root -l -b -q 'src/PlotJets.C("60GeVJets", "60GeVJets")'
#root -l -b -q 'src/PlotJets.C("30GeVJets", "30GeVJets")'

#root -l -b -q 'src/CalculateIAAShifts.C("60GeVJets", "60GeVJets")'
#root -l -b -q 'src/CalculateIAAShifts.C("30GeVJets", "30GeVJets")'

root -l -b -q 'src/PlotPtCh.C("60GeVJets", "60GeVJets")'
root -l -b -q 'src/PlotPtCh.C("30GeVJets", "30GeVJets")'
