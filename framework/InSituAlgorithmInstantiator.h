/**
 * @brief A singleton class used to instantiate AnalysisTools
 */
#ifndef InSituAlgorithmINSTANTIATOR_H_
#define InSituAlgorithmINSTANTIATOR_H_

#include "CosmoToolsMacros.h"

namespace cosmotk
{

// Forward declarations
class InSituAlgorithm;

class InSituAlgorithmInstantiator
{
public:
  InSituAlgorithmInstantiator();
  virtual ~InSituAlgorithmInstantiator();

  /**
   * @brief Return an analysis tool instance corresponding to the given name.
   * @param instanceName the name of the analysis tool to instantiate
   * @return tool pointer to a concrete InSituAlgorithm instance
   * @post if (tool == NULL) it indicates that the instanceName was not found.
   */
  static InSituAlgorithm* CreateInstance(std::string instanceName);

private:
  DISABLE_COPY_AND_ASSIGNMENT(InSituAlgorithmInstantiator);
};

} /* namespace cosmotk */
#endif /* ANALYSISTOOLINSTANTIATOR_H_ */
