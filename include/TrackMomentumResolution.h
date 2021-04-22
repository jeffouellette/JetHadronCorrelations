#ifndef __TrackMomentumResolution_h__
#define __TrackMomentumResolution_h__

namespace JetHadronCorrelations {

bool TrackMomentumResolution (const char* directory,
                              const int dataSet,
                              const char* inFileName = "",
                              const char* eventWeightsFileName = "");

} // end namespace

#endif
