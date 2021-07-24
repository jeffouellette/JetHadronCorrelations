#ifndef __Trigger_h__
#define __Trigger_h__

/**
 * Implements a trigger struct.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 * Last modified: 7/23/2021
 */

#include <string>


/**
 * A trigger struct stores information about a trigger, mainly decision & prescale information.
 * Can be adapted to include momentum cuts, pseudorapidity cuts, efficiencies, etc.
 * Base struct contains a name, decision, and prescale.
 */
struct Trigger {
    
  public:
    std::string name = "";

    bool trigDecision = false;
    float trigPrescale = -1;

    Trigger (const string _name)  { name = _name; }
    Trigger (const Trigger* t)    { name = t->name; }

    ~Trigger () {}

};

#endif
