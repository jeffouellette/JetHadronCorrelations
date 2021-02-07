#ifndef __run_h__
#define __run_h__

#include "Params.h"

using namespace HadronYieldsAnalysis;

namespace HadronYieldsAnalysis {

bool doHITightVar = false;
bool doPionsOnlyVar = false;
bool doWithPileupVar = false;

float crossSectionPicoBarns = 0;
float mcFilterEfficiency = 0;
int mcNumberEvents = 0;

CollisionSystem collisionSystem = CollisionSystem::pp15; // default is pp15
DataType dataType = DataType::Collisions; // default is collisions
TriggerType triggerType = TriggerType::Jet; // default is jet triggers

}

#endif
