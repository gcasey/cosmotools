/**
 * @brief A singleton class used to instantiate AnalysisTools
 */
#ifndef ANALYSISTOOLINSTANTIATOR_H_
#define ANALYSISTOOLINSTANTIATOR_H_

#include "CosmologyToolsMacros.h"

namespace cosmotk
{

// Forward declarations
class AnalysisTool;

class AnalysisToolInstantiator
{
public:
  AnalysisToolInstantiator();
  virtual ~AnalysisToolInstantiator();

  /**
   * @brief Return an analysis tool instance corresponding to the given name.
   * @param instanceName the name of the analysis tool to instantiate
   * @return tool pointer to a concrete AnalysisTool instance
   * @post if (tool == NULL) it indicates that the instanceName was not found.
   */
  static AnalysisTool* CreateInstance(std::string instanceName);

private:
  DISABLE_COPY_AND_ASSIGNMENT(AnalysisToolInstantiator);
};

} /* namespace cosmotk */
#endif /* ANALYSISTOOLINSTANTIATOR_H_ */
