#ifndef __JetEnergyResolution_h__
#define __JetEnergyResolution_h__

namespace JetHadronCorrelations {

bool JetEnergyResolution (const char* directory,
                          const int dataSet,
                          const char* inFileName = "",
                          const char* eventWeightsFileName = "");

} // end namespace

#endif
