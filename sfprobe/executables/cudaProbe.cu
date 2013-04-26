//all we have to do is define the adapter here and how we use boost threads
#ifdef USEDAX
#  define BOOST_SP_DISABLE_THREADS
#  define DAX_DEVICE_ADAPTER DAX_DEVICE_ADAPTER_CUDA
#endif

#include "probe.h"

int main(int argc, char **argv)
{
  return probe(argc,argv);
}

