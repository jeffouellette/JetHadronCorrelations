#root -l -b -q 'src/ProcessCorrelations.C("All", "AllJets_NIters20", 20)'
#root -l -b -q 'src/ProcessCorrelations.C("All", "AllJets_NIters4", 4)' # same as below

#./bin/ProcessCorrelations.exe All AllJets_NIters4 4
#./bin/ProcessCorrelations.exe All AllJets

#root -l -b -q 'src/ProcessUnfolding.C("AllJets", "AllJets")'
#./bin/ProcessUnfolding.exe AllJets AllJets

#root -l -b -q 'src/ProcessFineFCalPtCh.C("60GeVJets", "60GeVJets")'
#root -l -b -q 'src/ProcessFineFCalPtCh.C("30GeVJets", "30GeVJets")'
#root -l -b -q 'src/CombineFineFCalPtCh.C("60GeVJets")'
#root -l -b -q 'src/CombineFineFCalPtCh.C("30GeVJets")'

#root -l -b -q 'src/ProcessJets.C("All", "AllJets")'

#root -l -b -q 'src/PlotFineFCalPtCh.C("60GeVJets", "60GeVJets", "60GeVJets")'
#root -l -b -q 'src/PlotFineFCalPtCh.C("30GeVJets", "30GeVJets", "30GeVJets")'

#root -l -b -q 'src/PlotDPhi.C("60GeVJets", "60GeVJets")'
#root -l -b -q 'src/PlotDPhi.C("30GeVJets", "30GeVJets")'

#root -l -b -q 'src/PlotJets.C("All", "AllJets")'

#root -l -b -q 'src/PlotNIters.C("AllJets_NIters20")'

root -l -b -q 'src/PlotPtCh.C("AllJets")'
