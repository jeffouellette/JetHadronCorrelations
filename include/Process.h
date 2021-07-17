#ifndef __Process_h__
#define __Process_h__

#include "Params.h"

using namespace JetHadronCorrelations;

namespace JetHadronCorrelations {

float crossSectionPicoBarns = 0;
float mcFilterEfficiency = 0;
int mcNumberEvents = 0;

CollisionSystem collisionSystem = CollisionSystem::pp15; // default is pp15
DataType dataType = DataType::Collisions; // default is collisions
TriggerType triggerType = TriggerType::None; // default is no triggers
SystFlag systFlag = SystFlag::Nominal; // default is no systematics

}

#endif
