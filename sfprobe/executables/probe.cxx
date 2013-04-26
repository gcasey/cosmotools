
//setup this as the serial only version
#ifdef USEDAX
# undef USEDAX
#endif

#include "probe.h"

int main(int argc, char **argv)
{
  return probe(argc,argv);
}

