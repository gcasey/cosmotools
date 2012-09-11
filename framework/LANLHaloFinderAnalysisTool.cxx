#include "LANLHaloFinderAnalysisTool.h"
#include "Partition.h"
#include <cassert>

#include <mpi.h>

namespace cosmotk
{

LANLHaloFinderAnalysisTool::LANLHaloFinderAnalysisTool()
{
  this->Name               = "LANLHaloFinder";
  this->LinkingLength      = 0.2;
  this->PMIN               = 250;
  this->CenterFinderMethod = AVERAGE;
  this->ComputSODHalos     = false;
  this->RHO_C              = 2.77537e+11;
  this->INITIAL_SOD_MASS   = 1e+14;
  this->MIN_RADIUS_FACTOR  = 0.5;
  this->MAX_RADIUS_FACTOR  = 2;
  this->NUMBER_OF_BINS     = 20;
  this->FOF_SIZE_THRESHOLD = 500;
  this->Communicator       = MPI_COMM_NULL;

  this->HaloFinder = new cosmologytools::CosmoHaloFinderP();
}

//-----------------------------------------------------------------------------
LANLHaloFinderAnalysisTool::~LANLHaloFinderAnalysisTool()
{
  if(this->HaloFinder != NULL)
    {
    delete this->HaloFinder;
    this->HaloFinder = NULL;
    }
}

//-----------------------------------------------------------------------------
void LANLHaloFinderAnalysisTool::ParseParameters()
{
  // STEP 0: Parse basic parameters, as defined by super-class
  this->ParseBasicParameters();

  // STEP 1: Parse LANL halo-finder parameters
  this->LinkingLength = static_cast<REAL>(this->GetDoubleParameter("BB"));
  this->PMIN = this->GetIntParameter("PMIN");

  this->CenterFinderMethod = this->GetIntParameter("CENTER_FINDER_METHOD");
  assert("pre: Invalid center finder method" &&
          (this->CenterFinderMethod >= 0) &&
          (this->CenterFinderMethod < NUMBER_OF_CENTER_FINDER_METHODS));

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
}

//-----------------------------------------------------------------------------
void LANLHaloFinderAnalysisTool::Execute(SimulationParticles *particles)
{
  assert("pre: halofinder is NULL!" && (this->HaloFinder != NULL));
  assert("pre: input particles are NULL!" && (particles != NULL) );
  assert("pre: MPI communicator is NULL!" &&
         (this->Communicator != MPI_COMM_NULL) );

  // STEP 0: Initialize partition
  cosmologytools::Partition::initialize(this->Communicator);

  // STEP 1: short-circuit here in the unlike event that we don't have any
  // particles at this time-step
  if( particles->NumParticles == 0 )
    {
    return;
    }

  // STEP 2: parse the analysis tool parameters
  this->ParseParameters();

  // STEP 3: Setup halo-finder
  if( this->GenerateOutput )
    {
    this->HaloFinder->setParameters(
        "",this->BoxLength,this->NG,this->NDIM,
        this->PMIN, this->LinkingLength);
    }
  else
    {
    this->HaloFinder->setParameters(
        this->OutputFile,this->BoxLength,this->NG,this->NDIM,
        this->PMIN, this->LinkingLength);
    }

  // STEP 4: Register the particles with the halo-finder
  // NOTE: cast this to long here since the halo-finder stores the total
  // number of particles in an ivar that is a long.
  long numParticles = static_cast<long>(particles->NumParticles);
  this->HaloFinder->setParticles(
      particles->X,particles->Y,particles->Z,
      particles->VX,particles->VY,particles->VZ,
      particles->Potential,particles->GlobalIds,
      particles->Mask,particles->State,
      numParticles);

  // STEP 5: Run the halo-finder at each rank
  this->HaloFinder->executeHaloFinder();

  // STEP 6: Merge results across ranks
  this->HaloFinder->collectHalos(true /*clearSerial*/);
  this->HaloFinder->mergeHalos();
}

//-----------------------------------------------------------------------------
void LANLHaloFinderAnalysisTool::WriteOutput()
{
  assert("pre: halofinder is NULL!" && (this->HaloFinder != NULL));
  assert("pre: MPI communicator is NULL!" &&
         (this->Communicator != MPI_COMM_NULL) );

  if( !this->GenerateOutput )
    {
    return;
    }

  this->HaloFinder->writeTaggedParticles(true /*clearSerial*/);
}

} /* namespace cosmotk */
