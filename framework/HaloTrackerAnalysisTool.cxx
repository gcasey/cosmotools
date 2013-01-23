/**
 * @brief Implementation of HaloTracker
 */
#include "HaloTrackerAnalysisTool.h"

namespace cosmotk
{

namespace HaloTrackerInternals
{
  enum CenterFinder
  {
    AVERAGE = 0,
    CENTER_OF_MASS = 1,
    MBP = 2,
    MCP = 3,

    NUMBER_OF_CENTER_FINDER_METHODS
  };
} // END HaloTrackerInternals

HaloTrackerAnalysisTool::HaloTrackerAnalysisTool()
{
  this->Name               = "HALOTRACKER";
  this->LinkingLength      = 0.2;
  this->PMIN               = 250;
  this->CenterFinderMethod = HaloTrackerInternals::AVERAGE;
  this->ComputSODHalos     = false;
  this->RHO_C              = 2.77537e+11;
  this->INITIAL_SOD_MASS   = 1e+14;
  this->MIN_RADIUS_FACTOR  = 0.5;
  this->MAX_RADIUS_FACTOR  = 2;
  this->NUMBER_OF_BINS     = 20;
  this->FOF_SIZE_THRESHOLD = 500;
  this->Communicator       = MPI_COMM_NULL;

  this->HaloTracker = new cosmologytools::ForwardHaloTracker();
}

//------------------------------------------------------------------------------
HaloTrackerAnalysisTool::~HaloTrackerAnalysisTool()
{
  if( this->HaloTracker != NULL )
    {
    delete this->HaloTracker;
    }
}

//------------------------------------------------------------------------------
void HaloTrackerAnalysisTool::ParseParameters()
{
  // STEP 0: Parse basic parameters, as defined by the super-class
  this->ParseBasicParameters();

  // STEP 1: Parse halo-finder parameters
  this->LinkingLength = static_cast<REAL>(this->GetDoubleParameter("BB"));
  this->PMIN = this->GetIntParameter("PMIN");

  this->CenterFinderMethod = this->GetIntParameter("CENTER_FINDER_METHOD");
  assert(
    "pre: Invalid center finder method" &&
    (this->CenterFinderMethod >= 0) &&
    (this->CenterFinderMethod < HaloTrackerInternals::NUMBER_OF_CENTER_FINDER_METHODS));

  this->ComputSODHalos = this->GetBooleanParameter("COMPUTE_SOD_HALOS");

  if( this->ComputSODHalos )
    {
    this->RHO_C = static_cast<REAL>(this->GetDoubleParameter("RHO_C"));
    this->INITIAL_SOD_MASS =
        static_cast<REAL>(this->GetDoubleParameter("INITIAL_SOD_MASS"));
    this->MIN_RADIUS_FACTOR =
        static_cast<REAL>(this->GetDoubleParameter("MIN_RADIUS_FACTOR"));
    this->MAX_RADIUS_FACTOR =
        static_cast<REAL>(this->GetDoubleParameter("MAX_RADIUS_FACTOR"));
    this->NUMBER_OF_BINS = this->GetIntParameter("NUMBER_OF_BINS");
    this->FOF_SIZE_THRESHOLD = this->GetIntParameter("FOF_SIZE_THRESHOLD");
    }

  // STEP 2: Get Merger-Tree parameters
  this->MERGER_TREE_THRESHOLD =
      static_cast<INTEGER>(this->GetIntParameter("MERGER_TREE_THRESHOLD"));

  this->MERGER_TREE_FILE =
      this->GetStringParameter("MERGER_TREE_FILE");

  this->MERGER_TREE_FILE_FORMAT =
      this->GetIntParameter("MERGER_TREE_FILE_FORMAT");

}

//------------------------------------------------------------------------------
void HaloTrackerAnalysisTool::Execute(SimulationParticles *particles)
{
  assert("pre: input simulation particles instance is NULL!" &&
         (particles != NULL) );
  assert("pre: HaloTracker object is NULL" && (this->HaloTracker != NULL) );

  // STEP 0: parse the analysis tools parameters
  this->ParseParameters();

  // STEP 1: Setup Tracker
  this->HaloTracker->SetCommunicator(this->Communicator);
  this->HaloTracker->SetBoxLength(this->BoxLength);
  this->HaloTracker->SetLinkingLength(this->LinkingLength);
  this->HaloTracker->SetNG(this->NG);
  this->HaloTracker->SetNDIM(this->NDIM);
  this->HaloTracker->SetPMIN(this->PMIN);
  this->HaloTracker->SetMergerTreeFileName(this->MERGER_TREE_FILE);
  this->HaloTracker->SetMergerTreeThreshold(this->MERGER_TREE_THRESHOLD);

  // STEP 2: Register particles
  this->HaloTracker->RegisterParticles(
      particles->TimeStep,particles->RedShift,
      particles->X,particles->Y,particles->Z,
      particles->VX,particles->VY,particles->VZ,
      particles->Mass,
      particles->Potential,
      particles->GlobalIds,
      particles->Mask,
      particles->State,
      particles->NumParticles);

  // STEP 3: Track halos
  this->HaloTracker->TrackHalos();
}

//------------------------------------------------------------------------------
void HaloTrackerAnalysisTool::WriteOutput()
{
  if(this->GenerateOutput)
    {
    // TODO: implement this...The tracker should do a typical halo-dump here
    // since the merger-tree is only written at the end of the simulation.
    }
}

//------------------------------------------------------------------------------
std::string HaloTrackerAnalysisTool::GetInformation()
{
  return(this->GetBasicInformation());
}

} /* namespace cosmotk */
