#./bin/ProcessCorrelations.exe All AllJets
#python src/ProcessCovarianceMatrices.py AllJets AllJets
#./bin/ProcessUnfolding.exe AllJets AllJets

root -l -b -q 'src/ProcessNIters.C("AllJets", "AllJets_NIters20", 20)'

#root -l -b -q 'src/ProcessJets.C("All", "AllJets")'

#root -l -b -q 'src/PlotJets.C("All", "AllJets")'

root -l -b -q 'src/PlotNIters.C("AllJets_NIters20", "AllJets_NIters20")'

#root -l -b -q 'src/PlotPtCh.C("AllJets", "AllJets")'
#root -l -b -q 'src/PlotDPhi.C("AllJets", "AllJets")'


#./bin/ProcessCorrelations.exe All AllJets_Test
#python src/ProcessCovarianceMatrices.py AllJets_Test AllJets
#./bin/ProcessUnfolding.exe AllJets AllJets_Test
#root -l -b -q 'src/ProcessNIters.C("AllJets_Test", "AllJets_NIters20_Test", 20)'
#root -l -b -q 'src/PlotNIters.C("AllJets_Test", "AllJets_NIters20_Test")'
#root -l -b -q 'src/ProcessNIters.C("AllJets_test", "AllJets_NIters4_Test", 4)'
#root -l -b -q 'src/PlotNIters.C("AllJets_Test", "AllJets_NIters4_Test")'
