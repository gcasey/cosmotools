#include "CosmologyToolsManager.h"

#include "AnalysisTool.h"
#include "LANLHaloFinderAnalysisTool.h"
#include "CosmologyToolsConfiguration.h"
#include "SimulationParticles.h"

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

  this->Particles     = new SimulationParticles();
  this->Configuration = new CosmologyToolsConfiguration();
}

//------------------------------------------------------------------------------
CosmologyToolsManager::~CosmologyToolsManager()
{
  this->Communicator = MPI_COMM_NULL;

  if( this->Particles != NULL )
    {
    delete this->Particles;
    }

  if( this->Configuration != NULL )
    {
    delete this->Configuration;
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
  assert( "pre: ConfigurationFile should not be empty!" &&
           !this->ConfigurationFile.empty());

  if( this->Configuration != NULL )
    {
    delete this->Configuration;
    }

  this->Configuration = new CosmologyToolsConfiguration();
  this->Configuration->SetCommunicator( this->Communicator );
  this->Configuration->SetConfigFile(this->ConfigurationFile);
  this->Configuration->ParseFile();
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
void CosmologyToolsManager::CreateAnalysisTools()
{
  assert("pre: NULL configuration!" && (this->Configuration != NULL) );

  for(int i=0; i < this->Configuration->GetNumberOfAnalysisTools(); ++i)
    {
    std::string toolName = this->Configuration->GetToolName( i );
    if( this->Configuration->IsToolEnabled(toolName) )
      {
      AnalysisTool *myAnalysisTool = NULL;

      Dictionary toolParameters;
      this->Configuration->GetToolParameters(toolName,toolParameters);

      if( toolName == "LANLHaloFinder" )
        {
        myAnalysisTool = new LANLHaloFinderAnalysisTool();
        myAnalysisTool->SetDomainParameters(
            this->BoxLength,this->GhostOverlap,this->NDIM);
        myAnalysisTool->SetParameters(toolParameters);
        }
      else
        {
        std::cerr << "ERROR: undefined analysis tool: " << toolName;
        std::cerr << std::endl;
        std::cerr << "FILE: " << __FILE__ << std::endl;
        std::cerr << "LINE: " << __LINE__ << std::endl;
        continue;
        }

      assert("pre: analysis tool is NULL!" && (myAnalysisTool != NULL));
      this->AnalysisTools.push_back(myAnalysisTool);
      } // END if tool is enabled
    } // END for all analysis configuration tools
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::CoProcess()
{
  assert("pre: communicator is NULL!" && (this->Communicator!=MPI_COMM_NULL));
  assert("pre: simulation particles data-structure is NULL" &&
         (this->Particles != NULL));


  // STEP 0: Parse the configuration file
  this->ParseConfigurationFile();

  // STEP 1: Create a vector of analysis tools
  this->ClearAnalysisTools();
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


  // STEP 2: Loop through all analysis tools. For each tool, check if it
  // should be executed according to the user-supplied frequency criteria
  // given at the configuration file and execute it. After the analysis is
  // executed, check if the output should be generated and if output is
  // requested, write it to disk.
  for( int i=0; i < this->AnalysisTools.size(); ++i)
    {
    assert("pre: analysis tool is NULL" && (this->AnalysisTools[i] != NULL));

    // Check if the analsys tools should be executed at this time-step
    if( this->AnalysisTools[i]->ShouldExecute(this->Particles->TimeStep) )
      {
      if( this->Rank == 0 )
        {
        std::cout << "Executing: " << this->AnalysisTools[i]->GetName();
        std::cout << "...";
        std::cout.flush();
        }
      this->Barrier();

      // Run the analysis
      this->AnalysisTools[i]->Execute(this->Particles);

      if( this->Rank == 0 )
        {
        std::cout << "[DONE]\n";
        std::cout.flush();
        }
      this->Barrier();

      // Check if we should write output for this analysis algorithm
      if( this->AnalysisTools[i]->GetGenerateOutput() == true )
        {
        if( this->Rank == 0 )
          {
          std::cout << "Generating output for "
                    << this->AnalysisTools[i]->GetName();
          std::cout << "...";
          std::cout.flush();
          }
        this->Barrier();

        this->AnalysisTools[i]->WriteOutput();

        if( this->Rank == 0 )
          {
          std::cout << "[DONE]\n";
          std::cout.flush();
          }
        this->Barrier();
        }

      // Check if the output should be visualized at this time-step
      if( (this->Configuration->GetVisualization()==true) &&
          (this->AnalysisTools[i]->IsVisible()) )
        {
        // TODO: update the visualization!
        }
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
