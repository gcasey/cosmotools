/**
 * @brief A class that provides simple functionality for timing
 */

#ifndef TASKTIMER_H_
#define TASKTIMER_H_

#include "CosmoToolsMacros.h"

namespace cosmotk {

class TaskTimer
{
public:
  TaskTimer();
  virtual ~TaskTimer();

  /**
   * @brief Starts a timer.
   */
  void StartTimer();

  /**
   * @brief Stops the timer.
   */
  void StopTimer();

  /**
   * @return ellapsedTime the ellapsed time.
   */
  double GetEllapsedTime();

protected:
  double StartTime;
  double EndTime;

private:
  DISABLE_COPY_AND_ASSIGNMENT(TaskTimer);
};

} /* namespace cosmotk */
#endif /* TASKTIMER_H_ */
