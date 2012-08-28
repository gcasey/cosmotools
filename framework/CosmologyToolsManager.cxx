#include "CosmologyToolsManager.h"

#include "SimulationParticles.h"
#include "AnalysisTool.h"

#include <iostream>
#include <cassert>


namespace cosmotk {

//------------------------------------------------------------------------------
CosmologyToolsManager::CosmologyToolsManager()
{
  this->ConfigurationFile = "";
  this->Communicator      = MPI_COMM_NULL;
  this->Rank  = this->NumRanks = -1;
  this->NDIM  = 0;

  this->Particles = new SimulationParticles();
}

//------------------------------------------------------------------------------
CosmologyToolsManager::~CosmologyToolsManager()
{
  this->Communicator = MPI_COMM_NULL;

  if( this->Particles != NULL )
    {
    delete this->Particles;
    }

  this->ClearAnalysisTools();
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::Initialize(MPI_Comm comm)
{
  assert("pre: invalid communicator" && (comm != MPI_COMM_NULL) );
  this->Communicator = comm;
  MPI_Comm_size(this->Communicator,&this->NumRanks);
  MPI_Comm_rank(this->Communicator,&this->Rank);
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::SetDomainParameters(
        REAL boxlength, INTEGER ghostoverlap, INTEGER NDIM)
{
  this->BoxLength    = boxlength;
  this->GhostOverlap = ghostoverlap;
  this->NDIM         = NDIM;
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::SetAnalysisConfiguration(char *configfile)
{
  assert("pre: configfile is NULL" && (configfile != NULL));
  this->ConfigurationFile = std::string( configfile );
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::SetParticles(
        INTEGER tstep, REAL redShift,
        POSVEL_T *px, POSVEL_T *py, POSVEL_T *pz,
        POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz,
        REAL *mass, POTENTIAL_T *potential,
        ID_T *tags, MASK_T *mask, STATUS_T *status,
        ID_T N)
{
  this->Particles->SetParticles(
        tstep,redShift,
        px,py,pz,
        vx,vy,vz,
        mass,
        potential,
        tags, mask, status,
        N);
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::ParseConfigurationFile()
{
  // TODO: implement this
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::ClearAnalysisTools()
{
  for( int i=0; i < this->AnalysisTools.size(); ++i)
    {
    if(this->AnalysisTools[i] != NULL)
      {
      delete this->AnalysisTools[i];
      }
    this->AnalysisTools[i] = NULL;
    }
  this->AnalysisTools.clear();
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::CoProcess()
{
  assert("pre: communicator is NULL!" && (this->Communicator!=MPI_COMM_NULL));
  assert("pre: simulation particles data-structure is NULL" &&
         (this->Particles != NULL));


  this->ParseConfigurationFile();
  this->CreateAnalysisTools();
  this->Barrier();

  if(this->Rank == 0 )
    {
    std::cout << "====================\n";
    std::cout << "COSMO TOOLS\n";
    std::cout << "Number of Analysis tools: " << this->AnalysisTools.size();
    std::cout << std::endl << std::endl;
    std::cout << "TIME STEP: " << this->Particles->TimeStep << std::endl;
    std::cout << "RED SHIFT: " << this->Particles->RedShift << std::endl;
    std::cout.flush();
    }


  for( int i=0; i < this->AnalysisTools.size(); ++i)
    {
    assert("pre: analysis tool is NULL" && (this->AnalysisTools[i] != NULL));

    if( this->AnalysisTools[i]->ShouldExecute(this->Particles->TimeStep) )
      {
      if( this->Rank == 0 )
        {
        std::cout << "Executing: " << this->AnalysisTools[i]->GetName();
        std::cout << "..." << std::endl;
        std::cout.flush();
        }
      this->Barrier();

      this->AnalysisTools[i]->Execute(this->Particles);

      if( this->Rank == 0 )
        {
        std::cout << "[DONE]\n";
        std::cout.flush();
        }
      this->Barrier();
      } // END if ShouldExecute


    } // END for all analysis tools
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::Finalize()
{
  this->Communicator      = MPI_COMM_NULL;
  this->ConfigurationFile = "";
  this->ClearAnalysisTools();
}

} // END cosmotk namespace
