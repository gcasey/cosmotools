#ifndef MERGERTREEEVENT_H_
#define MERGERTREEEVENT_H_

namespace cosmotk {

class MergerTreeEvent
{
public:
  enum
    {
    SPLIT = 0,
    MERGE,
    CONTINUATION,
    BIRTH,
    DEATH,
    UNDEFINED
    };
};

}


#endif
