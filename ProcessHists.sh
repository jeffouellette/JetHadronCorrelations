#root -l -b -q 'src/ProcessJets.C("All", "AllJets")'
#root -l -b -q 'src/PlotJets.C("All", "AllJets")'

#root -l -b -q 'src/PlotResponseMatrix.C+(false)'

#./bin/ProcessCorrelations.exe All AllJets
#python src/ProcessCovarianceMatrices.py AllJets AllJets
./bin/ProcessUnfolding.exe AllJets AllJets

#./bin/ProcessNIters.exe AllJets AllJets_NIters1 1 25
#./bin/ProcessNIters.exe AllJets AllJets_NIters20 20 25
#root -l -b -q 'src/PlotNIters.C("AllJets", "AllJets_NIters20", 20, 25)'

root -l -b -q 'src/PlotPtCh.C("AllJets", "AllJets")'
#root -l -b -q 'src/PlotDPhi.C("AllJets", "AllJets")'


# FOR TESTING
#./bin/ProcessCorrelations.exe All AllJets_Test
#python src/ProcessCovarianceMatrices.py AllJets_Test AllJets

#root -l -b -q 'src/ProcessNIters.C("AllJets_Test", "AllJets_NIters20", 20, 300)'
#root -l -b -q 'src/PlotNIters.C("AllJets_Test", "AllJets_NIters20")'
#root -l -b -q 'src/PlotNIters.C("AllJets_Test", "AllJets_NIters20", 20, 300)'

#root -l -b -q 'src/ProcessNIters.C("AllJets_Test", "AllJets_NIters20_JZ0123", 20)'
#root -l -b -q 'src/PlotNIters.C("AllJets_Test", "AllJets_NIters20_JZ0123")'

#root -l -b -q 'src/ProcessNIters.C("AllJets_Test", "AllJets_NIters4_Test", 4)'
#root -l -b -q 'src/PlotNIters.C("AllJets_Test", "AllJets_NIters4_Test")'

#./bin/ProcessUnfolding.exe AllJets_Test AllJets_Test

#root -l -b -q 'src/PlotPtCh.C("AllJets_Test", "AllJets_Test")'
