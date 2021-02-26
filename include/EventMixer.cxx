#ifndef __JetHadronSkimmer_h__
#define __JetHadronSkimmer_h__

#include "TreeVariables.h"

namespace JetHadronCorrelations {

bool MeetsJetCuts (int iJ);

bool JetHadronSkimmer (const char* directory,
                       const int dataSet,
                       const char* inFileName = "");

} // end namespace

#endif
