/**
 * @brief Implementation of TaskTimer
 */
#include "TaskTimer.h"

#include <mpi.h>

namespace cosmotk {

TaskTimer::TaskTimer()
{
	this->StartTime = 0.0;
	this->EndTime	= 0.0;
}

//------------------------------------------------------------------------------
TaskTimer::~TaskTimer()
{

}

//------------------------------------------------------------------------------
void TaskTimer::StartTimer()
{
	this->StartTime = MPI_Wtime();
	this->EndTime = this->StartTime;
}

//------------------------------------------------------------------------------
void TaskTimer::StopTimer()
{
	this->EndTime = MPI_Wtime();
}

//------------------------------------------------------------------------------
double TaskTimer::GetEllapsedTime()
{
	return( this->EndTime-this->StartTime );
}

} /* namespace cosmotk */
