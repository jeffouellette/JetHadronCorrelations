#ifndef __JetHadronSkimmer_h__
#define __JetHadronSkimmer_h__

#include "TreeVariables.h"

namespace HadronYieldsAnalysis {

bool JetHadronSkimmer (const char* directory,
                       const int dataSet,
                       const char* inFileName = "");

} // end namespace

#endif
