#ifndef __JetPtWeights_h__
#define __JetPtWeights_h__

namespace JetHadronCorrelations {

bool JetPtWeights (const char* directory,
                   const int dataSet,
                   const char* inFileName = "",
                   const char* eventWeightsFileName = "");

} // end namespace

#endif
