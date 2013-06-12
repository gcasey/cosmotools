#include "InSituAlgorithmInstantiator.h"

#include "InSituAlgorithm.h"
#include "LANLHaloFinderInSituAlgorithm.h"

#if defined(USEDIY) && defined(USEQHULL)
  #include "TessVoidFinderInSituAlgorithm.h"
#endif

#include <iostream>
#include <cassert>

namespace cosmotk
{

InSituAlgorithmInstantiator::InSituAlgorithmInstantiator()
{
  // TODO Auto-generated constructor stub

}

//-----------------------------------------------------------------------------
InSituAlgorithmInstantiator::~InSituAlgorithmInstantiator()
{
  // TODO Auto-generated destructor stub
}

//-----------------------------------------------------------------------------
InSituAlgorithm* InSituAlgorithmInstantiator::CreateInstance(
    std::string instanceName)
{
  InSituAlgorithm *tool = NULL;
  if(instanceName == "LANLHALOFINDER")
    {
    tool = new LANLHaloFinderInSituAlgorithm();
    }
#if defined(USEDIY) && defined(USEQHULL)
  else if(instanceName == "TESS")
    {
    tool = new TessVoidFinderInSituAlgorithm();
    }
#endif
  else
    {
    std::cerr << "ERROR: undefined analysis tool: " << instanceName;
    std::cerr << std::endl;
    std::cerr << "FILE: " << __FILE__ << std::endl;
    std::cerr << "LINE: " << __LINE__ << std::endl;
    }
  return( tool );
}

} /* namespace cosmotk */
