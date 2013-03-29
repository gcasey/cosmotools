#ifndef MERGERTREEEVENT_H_
#define MERGERTREEEVENT_H_

#include "CosmologyToolsMacros.h"

#include <sstream>
#include <string>

/**
 * @brief Corresponding strings to event names
 */
static const char *EventName[] = {
  "(**SPLIT**)",
  "(**MERGE**)",
  "(**CONTINUATION**)",
  "(**BIRTH**)",
  "(**DEATH**)",
  "(**REBIRTH**)"
};

namespace cosmotk {

class MergerTreeEvent
{
public:

  enum
    {
    UNDEFINED = -1,
    SPLIT,
    MERGE,
    CONTINUATION,
    BIRTH,
    DEATH,
    REBIRTH,

    NUMBER_OF_EVENTS
    } EventTypes;


  /**
   * @brief Sets the user-supplied event in the given bitmask
   * @param mask bitmask where to set the given event
   * @param event the event to be set
   * @post MergerTreeEvent::IsEvent(mask,event)==true
   */
  static void SetEvent(unsigned char &mask, const int event)
  {
    assert("pre: event is out-of-range!" &&
           (event >= 0) && (event < MergerTreeEvent::NUMBER_OF_EVENTS) );
    assert("pre: insufficient number of bits to encode events!" &&
            (MergerTreeEvent::NUMBER_OF_EVENTS <= 8));
    mask |= (1 << event);
  }

  /**
   * @brief Unsets the user-supplied event from the given bitmask
   * @param mask bitmask to unset the given event from
   * @param event the event to unset
   * @post MergerTreeEvent::IsEvent(mask,event)==false
   */
  static void UnsetEvent(unsigned char &mask, const int event)
  {
    assert("pre: event is out-of-range!" &&
           (event >= 0) && (event < 8) );
    assert("pre: insufficient number of bits to encode events!" &&
            (MergerTreeEvent::NUMBER_OF_EVENTS <= 8));
    mask &= ~(1 << event);
  }

  /**
   * @brief Checks if the event in query is set in the given bitmask
   * @param mask bitmask to check
   * @param event the event to check
   * @return status true iff the event is checked, else, false.
   */
  static bool IsEvent(unsigned char &mask, const int event)
  {
    assert("pre: event is out-of-range!" &&
           (event >= 0) && (event < MergerTreeEvent::NUMBER_OF_EVENTS) );
    assert("pre: insufficient number of bits to encode events!" &&
            (MergerTreeEvent::NUMBER_OF_EVENTS <= 8));
    bool status = ( mask & (1 << event) )? true : false;
    return( status );
  }

  /**
   * @brief Resets all bits in the given mask to OFF
   * @param mask the bitmask that is being reset
   * @post MergerEvent::IsEvent(mask,i)==false \f$ \forall i \in \[0,7\] \f$
   */
  static void Reset(unsigned char &mask)
  {
    for(int i=0; i < 8; ++i)
      {
      MergerTreeEvent::UnsetEvent(mask,i);
      } // END for all bits of the mask
  }

  /**
   * @brief Returns a string represenation of the events encoded in the given
   * user-supplied event mask.
   * @param eventmask the encoded event mask in question.
   * @return s a string corresponding to the set of events that are encoded.
   */
  static std::string GetEventString(unsigned char eventmask)
  {
    std::ostringstream oss;
    for(int i=0; i < MergerTreeEvent::NUMBER_OF_EVENTS; ++i)
      {
      if( MergerTreeEvent::IsEvent(eventmask,i))
        {
        oss << EventName[ i ];
        } // END if event is set
      } // END for all possible events
    return( oss.str() );
  }

protected:
    MergerTreeEvent();
    ~MergerTreeEvent();

private:
    DISABLE_COPY_AND_ASSIGNMENT(MergerTreeEvent);
};

}


#endif
