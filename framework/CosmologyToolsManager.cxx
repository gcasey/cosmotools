#include "CosmologyToolsManager.h"

#include "AnalysisTool.h"
#include "LANLHaloFinderAnalysisTool.h"
#include "CosmologyToolsConfiguration.h"
#include "SimulationParticles.h"

#include <iostream>
#include <cassert>

#define PRINTLN( str ) {               \
  if(this->Rank==0) {                  \
    std::cout str << std::endl;         \
    std::cout.flush();                  \
  }                                     \
  this->Barrier();                     \
}

#define PRINT( str ) {    \
  if(this->Rank==0) {     \
    std::cout str;         \
    std::cout.flush();     \
  }                        \
  this->Barrier();        \
}

namespace cosmotk {

//------------------------------------------------------------------------------
CosmologyToolsManager::CosmologyToolsManager()
{
  this->ConfigurationFile = "";
  this->Communicator      = MPI_COMM_NULL;
  this->Rank  = this->NumRanks = -1;
  this->NDIM  = 0;
  this->XYZPeriodic   = false;
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
        REAL boxlength, INTEGER ghostoverlap, INTEGER NDIM,
        bool XYZPeriodic)
{
  this->BoxLength    = boxlength;
  this->GhostOverlap = ghostoverlap;
  this->NDIM         = NDIM;
  this->XYZPeriodic  = XYZPeriodic;
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
  std::map<std::string,AnalysisTool*>::iterator iter;
  iter = this->AnalysisTools.begin();
  for(; iter != this->AnalysisTools.end(); ++iter)
    {
    assert("pre: AnalysisTool is NULL!" && (iter->second!=NULL) );
    delete iter->second;
    } // End for all analysis tools
  this->AnalysisTools.clear();
}

//------------------------------------------------------------------------------
bool CosmologyToolsManager::ToolExists(const std::string &toolName)
{
  if(this->AnalysisTools.find(toolName) != this->AnalysisTools.end() )
    {
    return true;
    }
  return false;
}

//------------------------------------------------------------------------------
AnalysisTool* CosmologyToolsManager::GetToolByName(
    const std::string &name)
{
  AnalysisTool *tool = NULL;
  if( this->ToolExists(name) )
    {
    tool = this->AnalysisTools[name];
    }
  else
    {
    if( name == "LANLHALOFINDER")
      {
      tool = new LANLHaloFinderAnalysisTool();
      }
    else
      {
      std::cerr << "ERROR: undefined analysis tool: " << name;
      std::cerr << std::endl;
      std::cerr << "FILE: " << __FILE__ << std::endl;
      std::cerr << "LINE: " << __LINE__ << std::endl;
      }
    }
  tool->SetCommunicator( this->Communicator );
  return( tool );
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::CreateAnalysisTools()
{
  assert("pre: NULL configuration!" && (this->Configuration != NULL) );

  AnalysisTool *tool = NULL;
  for(int i=0; i < this->Configuration->GetNumberOfAnalysisTools(); ++i)
    {
    std::string toolName = this->Configuration->GetToolName( i );

    Dictionary toolParameters;
    this->Configuration->GetToolParameters(toolName,toolParameters);

    tool = this->GetToolByName( toolName );
    assert("tool is NULL!" && (tool != NULL) );
    if( this->Configuration->IsToolEnabled(toolName))
      {
      tool->SetEnabled(true);
      }
    else
      {
      tool->SetEnabled(false);
      }
    tool->SetDomainParameters(
        this->BoxLength,this->GhostOverlap,this->NDIM);
    tool->SetParameters(toolParameters);
    this->AnalysisTools[ toolName ] = tool;
    } // END for all tools in the configuration file
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
  this->CreateAnalysisTools();
  this->Barrier();

  PRINTLN(
      << "\n\n====================\n"
      << "COSMO TOOLS\n"
      << "Number of Analysis tools: "
      << this->AnalysisTools.size()
      << std::endl << std::endl
      << "TIME STEP: " << this->Particles->TimeStep << std::endl
      << "RED SHIFT: " << this->Particles->RedShift );


  // STEP 2: Loop through all analysis tools. For each tool, check if it
  // should be executed according to the user-supplied frequency criteria
  // given at the configuration file and execute it. After the analysis is
  // executed, check if the output should be generated and if output is
  // requested, write it to disk.
  std::map<std::string,AnalysisTool*>::iterator toolIter;
  toolIter=this->AnalysisTools.begin();
  for( ;toolIter!=this->AnalysisTools.end(); ++toolIter)
    {
    std::string toolName = toolIter->first;
    AnalysisTool *tool   = toolIter->second;
    assert("pre: NULL tool!" && (tool != NULL) );

    PRINTLN(<< tool->GetInformation() );

    if( tool->ShouldExecute(this->Particles->TimeStep) )
      {
      PRINT(<< "Executing: " << tool->GetName() << "...");
      tool->Execute(this->Particles);
      PRINTLN(<< "[DONE]");

      if(tool->GetGenerateOutput() == true)
        {
        PRINT(<<"Writing output...");
        tool->WriteOutput();
        PRINTLN(<< "[DONE]");
        } // END if should write output

      if(this->Configuration->GetVisualization() &&
          tool->IsVisible())
        {
        // TODO: implement this
        } // END if the tool is visible
      } // END if the tool should execute

    } // END for all tools

  // STEP 3: Synchronize all ranks
  PRINTLN(<< "Finished co-processing @t=" << this->Particles->TimeStep);
  PRINTLN(<< "====================");
  this->Barrier();
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::Finalize()
{
  this->Communicator      = MPI_COMM_NULL;
  this->ConfigurationFile = "";
  this->ClearAnalysisTools();
}

} // END cosmotk namespace
