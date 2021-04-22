#ifndef __TrackingPerformance_h__
#define __TrackingPerformance_h__

namespace JetHadronCorrelations {

bool TrackingPerformance (const char* directory,
                          const int dataSet,
                          const char* inFileName = "",
                          const char* eventWeightsFileName = "");

} // end namespace

#endif
