
//all we have to do is define the adapter here
#ifdef USEDAX
#  define DAX_DEVICE_ADAPTER DAX_DEVICE_ADAPTER_TBB
#endif

//needed to force the number of tbb threads to a give number
#include <tbb/task_scheduler_init.h>

#include "probe.h"

int main(int argc, char **argv)
{
  int num_cores = tbb::task_scheduler_init::automatic;

  //remove the last argument which is the number of cores to init
  if(argc < 6 || argc > 7)
    {
    std::cerr << "USAGE: ./probe <basedir> <rL> <ndim>";
    std::cerr << " <ndim2> <fringe> <tbb_cores (optional)>" << std::endl;
    }
  else if(argc == 7)
    {
    //we get to setup tbb number of cores
    num_cores = atoi(argv[6]);
    }

  //setup tbb with the user given number of cores, otherwise
  //we default to the automatic selection
  tbb::task_scheduler_init schedulerInit(num_cores);

  argc = 6; //trim argc to what probe expects
  int result = probe(argc,argv);
  return result;
}

