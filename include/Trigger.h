#ifndef __Trigger_h__
#define __Trigger_h__

/**
 * Implements a trigger struct.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 * Last modified: 11/12/2020
 */

#include <string>

using namespace std;

/**
 * A trigger struct stores information about a trigger, mainly decision & prescale information.
 * Can be adapted to include momentum cuts, pseudorapidity cuts, efficiencies, etc.
 * Base struct contains a name, decision, and prescale.
 */
struct Trigger {
    
  public:
    string name = "";

    bool trigDecision = false;
    float trigPrescale = -1;

    Trigger (const string _name)  { name = _name; }
    Trigger (const Trigger* t)    { name = t->name; }

    ~Trigger () {}

};

#endif
