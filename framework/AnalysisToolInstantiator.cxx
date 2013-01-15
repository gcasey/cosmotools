#include "AnalysisToolInstantiator.h"

#include "AnalysisTool.h"
#include "LANLHaloFinderAnalysisTool.h"
#include "TessVoidFinderAnalysisTool.h"

#include <iostream>
#include <cassert>

namespace cosmotk
{

AnalysisToolInstantiator::AnalysisToolInstantiator()
{
  // TODO Auto-generated constructor stub

}

//-----------------------------------------------------------------------------
AnalysisToolInstantiator::~AnalysisToolInstantiator()
{
  // TODO Auto-generated destructor stub
}

//-----------------------------------------------------------------------------
AnalysisTool* AnalysisToolInstantiator::CreateInstance(
    std::string instanceName)
{
  AnalysisTool *tool = NULL;
  if(instanceName == "LANLHALOFINDER")
    {
    tool = new LANLHaloFinderAnalysisTool();
    }
  else if(instanceName == "TESS")
    {
    tool = new TessVoidFinderAnalysisTool();
    }
  else if(instanceName == "HALOTRACKER")
    {
    tool = new HaloTrackerAnalysisTool();
    }
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
