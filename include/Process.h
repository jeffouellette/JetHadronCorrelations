#ifndef __Process_h__
#define __Process_h__

#include "Params.h"

using namespace JetHadronCorrelations;

namespace JetHadronCorrelations {

bool doHITightVar = false;
bool doPionsOnlyVar = false;
bool doWithPileupVar = false;
bool doJetES5PercUpVar = false;
bool doJetES5PercDownVar = false;
bool doJetES2PercUpVar = false;
bool doJetES2PercDownVar = false;

float crossSectionPicoBarns = 0;
float mcFilterEfficiency = 0;
int mcNumberEvents = 0;

CollisionSystem collisionSystem = CollisionSystem::pp15; // default is pp15
DataType dataType = DataType::Collisions; // default is collisions
TriggerType triggerType = TriggerType::Jet; // default is jet triggers

}

#endif
