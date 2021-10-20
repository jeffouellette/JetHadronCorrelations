CXX=g++
CXXFLAGS=-O3 -g -Wall -fPIC -std=c++1y `root-config --cflags` -Iinclude -I${ROOT_UTILS_PATH}/include -I${ATLAS_PATH}/include -I${ROOUNFOLD_INCLUDE_DIR}
LDFLAGS=`root-config --glibs --ldflags` -Llib -L${ROOT_UTILS_PATH}/lib -L${ATLAS_PATH}/lib -L${ROOUNFOLD_LIBRARY} -lUtilities -lAtlasUtils -lAtlasStyle -lRooUnfold

libraries = LocalUtilities
algorithms = JetHadronSkimmer CentralityAnalysis JetSubtractedEnergy TrackingPerformance TrackMomentumResolution JetEnergyResolution
binaries = Process AnalyzeTrackMomentumResolution AnalyzeJetEnergyResolution RunCorrelator MakeResponseMatrix

.PHONY : libs algs bins directories clean

all : directories libs algs bins

libs : $(libraries)
algs : $(algorithms)
bins : $(binaries)

directories :
	mkdir -p bin lib obj

LocalUtilities : src/LocalUtilities.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o lib/lib$@.so src/$@.cxx

# Main needs to be compiled into binary
Process : $(libraries) $(algorithms) src/Process.cxx
	$(CXX) $(CXXFLAGS) src/Process.cxx $(LDFLAGS) $(libraries:%=-l%) $(algorithms:%=-l%) -o bin/Process.exe

RunCorrelator : $(libraries) src/RunCorrelator.cxx
	$(CXX) $(CXXFLAGS) src/RunCorrelator.cxx $(LDFLAGS) $(libraries:%=-l%) $(algorithms:%=-l%) -o bin/RunCorrelator.exe

MakeResponseMatrix : $(libraries) src/MakeResponseMatrix.cxx
	$(CXX) $(CXXFLAGS) src/MakeResponseMatrix.cxx $(LDFLAGS) $(libraries:%=-l%) $(algorithms:%=-l%) -o bin/MakeResponseMatrix.exe

AnalyzeTrackMomentumResolution : $(libraries) src/AnalyzeTrackMomentumResolution.C
	$(CXX) $(CXXFLAGS) src/AnalyzeTrackMomentumResolution.C $(LDFLAGS) $(libraries:%=-l%) $(algorithms:%=-l%) -o bin/AnalyzeTrackMomentumResolution.exe

AnalyzeJetEnergyResolution : $(libraries) src/AnalyzeJetEnergyResolution.C
	$(CXX) $(CXXFLAGS) src/AnalyzeJetEnergyResolution.C $(LDFLAGS) $(libraries:%=-l%) $(algorithms:%=-l%) -o bin/AnalyzeJetEnergyResolution.exe

obj/%.o : src/%.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $<

lib/lib%.so : %

# Algorithms need to be compiled into object files and then migrated to a shared library
% : obj/%.o $(libraries)
	$(CXX) -shared $(LDFLAGS) $(libraries:%=-l%) -o lib/lib$@.so $<

clean :
	rm -rf ./lib/*.so*
	rm -rf ./bin/*
	rm -rf ./obj/*
