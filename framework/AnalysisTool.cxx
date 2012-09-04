#include "AnalysisTool.h"

#include "SimulationParticles.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

namespace cosmotk
{

AnalysisTool::AnalysisTool()
{
  this->Name              = "";
  this->GenerateOutput    = false;
  this->OutputFile        = "";
  this->FrequencyType     = EXPLICIT;
  this->ImplicitFrequency = 1;
  this->VisibilityStatus  = true;
  this->Communicator      = MPI_COMM_NULL;
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

//-----------------------------------------------------------------------------
double AnalysisTool::GetDoubleParameter(std::string key)
{
  assert("pre: string parameter is empty!" && (!key.empty()) );
  assert("pre: parameter does not exist!" &&
          this->Parameters.find(key) != this->Parameters.end() );
  std::string s = this->Parameters[key];
  return(std::atof(s.c_str()));
}

//-----------------------------------------------------------------------------
int AnalysisTool::GetIntParameter(std::string key)
{
  assert("pre: string parameter is empty!" && (!key.empty()) );
  assert("pre: parameter does not exist!" &&
          this->Parameters.find(key) != this->Parameters.end() );
  std::string s = this->Parameters[key];
  return(std::atoi(s.c_str()));
}

//-----------------------------------------------------------------------------
std::set<int> AnalysisTool::GetIntListParameter(std::string key)
{
  assert("pre: string parameter is empty!" && (!key.empty()) );
  assert("pre: parameter does not exist!" &&
            this->Parameters.find(key) != this->Parameters.end() );

  std::set<int> intlist;
  std::string s = this->Parameters[key];
  std::istringstream iss( s );
  std::vector<std::string> tokens;
  std::copy(std::istream_iterator<std::string>(iss),
           std::istream_iterator<std::string>(),
           std::back_inserter< std::vector<std::string> >(tokens));

  for( int i=0; i < static_cast<int>(tokens.size()); ++i )
    {
    intlist.insert( this->GetIntParameter(tokens[i]));
    } // END for

  return(intlist);
}

//-----------------------------------------------------------------------------
bool AnalysisTool::GetBooleanParameter(std::string key)
{
  assert("pre: string parameter is empty!" && (!key.empty()) );
  assert("pre: parameter does not exist!" &&
            this->Parameters.find(key) != this->Parameters.end() );
  std::string s = this->Parameters[key];
  std::transform(s.begin(),s.end(),s.begin(),::toupper);
  if( (s == "YES") || (s == "TRUE") || (s == "ON") || (s == "ENABLED") )
    {
    return true;
    }
  return false;
}

//-----------------------------------------------------------------------------
std::string AnalysisTool::GetStringParameter(std::string key)
{
  assert("pre: string parameter is empty!" && (!key.empty()) );
  assert("pre: parameter does not exist!" &&
            this->Parameters.find(key) != this->Parameters.end() );
  return( this->Parameters[key] );
}

//-----------------------------------------------------------------------------
void AnalysisTool::ParseBasicParameters()
{
  assert( "pre: AnalysisTool parameters are empty" &&
          (!this->Parameters.empty()) );

  // STEP 0: Get frequency parameters
  this->FrequencyType = this->GetIntParameter("FREQUENCY_TYPE");
  switch( this->FrequencyType )
    {
    case EXPLICIT:
      this->ExplicitTimeSteps =
          this->GetIntListParameter("EXPLICIT_TIMESTEPS");
      break;
    case IMPLICIT:
      this->ImplicitFrequency = this->GetIntParameter("IMPLICIT_TIMESTEPS");
      break;
    default:
      std::cerr << "ERROR: Undefined frequency type!\n";
      std::cerr << "FILE: " << __FILE__ << std::endl;
      std::cerr << "LINE: " << __LINE__ << std::endl;
      return;
    } // END frequency type

  // STEP 1: Get I/O parameters
  this->GenerateOutput = this->GetBooleanParameter("WRITE_OUTPUT");
  if( this->GenerateOutput )
    {
    this->OutputFile = this->GetStringParameter("BASE_OUTPUT_FILE_NAME");
    }

  // STEP 2: Get visibility parameters
  this->VisibilityStatus = this->GetBooleanParameter("VISIBLE");
}

} /* namespace cosmotk */
