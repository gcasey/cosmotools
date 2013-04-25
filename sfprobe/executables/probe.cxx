#include "probe.h"

//all we have to do is define the adapter here
#ifdef USEDAX
#  define DAX_DEVICE_ADAPTER DAX_DEVICE_ADAPTER_TBB
#endif

int main(int argc, char **argv)
{
  return probe(argc,argv);
}

