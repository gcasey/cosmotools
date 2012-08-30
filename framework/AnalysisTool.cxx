#include "AnalysisTool.h"

#include "SimulationParticles.h"

#include <iostream>
#include <cassert>

namespace cosmotk
{

AnalysisTool::AnalysisTool()
{
  this->Name       = "";

  this->GenerateOutput = false;
  this->OutputFile = "";

  this->FrequencyType = EXPLICIT;
  this->ImplicitFrequency = 1;

  this->VisibilityStatus = true;
}

//-----------------------------------------------------------------------------
AnalysisTool::~AnalysisTool()
{
  this->ExplicitTimeSteps.clear();
}

//-----------------------------------------------------------------------------
void AnalysisTool::SetExplicitTimeSteps(INTEGER *tsteps, int N )
{
  assert("pre: timesteps array is NULL" && (tsteps != NULL) );

  for(int i=0; i < N; ++i)
    {
    this->ExplicitTimeSteps.insert( tsteps[i] );
    } // END for all timesteps
}

//-----------------------------------------------------------------------------
bool AnalysisTool::ShouldExecute(INTEGER ts)
{
  bool status = false;
  switch( this->FrequencyType )
    {
    case EXPLICIT:
      if(this->ExplicitTimeSteps.find(ts) != this->ExplicitTimeSteps.end() )
        {
        status = true;
        }
      else
        {
        status = false;
        }
      break;
    case IMPLICIT:
      if( (ts % this->ImplicitFrequency) == 0 )
        {
        status = true;
        }
      else
        {
        status = false;
        }
      break;
    default:
      std::cerr << "ERROR: Undefined frequency type!" << std::endl;
      status = false;
    }

  return( status );
}

} /* namespace cosmotk */
