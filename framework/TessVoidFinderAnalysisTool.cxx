#include "TessVoidFinderAnalysisTool.h"

#include "tess.h"

#include <cassert>

namespace cosmotk
{

TessVoidFinderAnalysisTool::TessVoidFinderAnalysisTool()
{
  this->Name         = "TESS";
  this->Communicator = MPI_COMM_NULL;
  this->MinVol = this->MaxVol = -1.0;
  this->GhostFactor = 5.0;
}

//------------------------------------------------------------------------------
TessVoidFinderAnalysisTool::~TessVoidFinderAnalysisTool()
{
  // TODO: implement this
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::ParseParameters()
{
  // STEP 0: parse basic parameters, as defined by super-class
  this->ParseBasicParameters();

  // STEP 1: parse tess parameters
  this->GhostFactor =
      static_cast<REAL>(this->GetDoubleParameter("GHOST_ZONE_SIZE"));
  this->MinVol =
      static_cast<REAL>(this->GetDoubleParameter("MIN_VOL_THRESHOLD"));
  this->MaxVol =
      static_cast<REAL>(this->GetDoubleParameter("MAX_VOL_THRESHOLD"));
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::Execute(SimulationParticles *particles)
{
  assert("pre: input particles are NULL!" && (particles != NULL));
  assert("pre: MPI communicator is NULL!" &&
         (this->Communicator != MPI_COMM_NULL) );

  // STEP 0: short-circuit here
  if( particles->NumParticles == 0 )
    {
    return;
    }

  // STEP 1: parse the analysis tool parameters
  this->ParseParameters();

  // STEP 2: initialize tess, i.e., call tess_init()
  // TODO: implement this

  // STEP 3: execute tess on the given particle dataset, i.e., call tess()
  // TODO: implement this

  // STEP 4: finalize tess
  tess_finalize();
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::WriteOutput()
{
  if( !this->GenerateOutput )
    {
    return;
    }

  // TODO: write output from tess. If tess is finalized, how to we do I/O here?
}

//------------------------------------------------------------------------------
std::string TessVoidFinderAnalysisTool::GetInformation()
{
  return(this->GetBasicInformation());
}

} /* namespace cosmotk */
