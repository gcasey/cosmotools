#include "CosmologyToolsManager.h"

#include "AnalysisTool.h"
#include "AnalysisToolInstantiator.h"
#include "CosmologyToolsConfiguration.h"
#include "SimulationParticles.h"


#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>

#define PRINTLN( str ) {        \
  if(this->Rank==0) {           \
    std::cout str << std::endl;  \
    std::cout.flush();           \
  }                              \
  this->Barrier();              \
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
  this->Rank              = this->NumRanks = -1;
  this->NDIM              = 0;
  this->ConfigurationIsRead = false;
  this->XYZPeriodic       = false;
  this->Particles         = new SimulationParticles();
  this->Configuration     = new CosmologyToolsConfiguration();
}

//------------------------------------------------------------------------------
CosmologyToolsManager::~CosmologyToolsManager()
{
   if( this->Particles != NULL )
    {
    delete this->Particles;
    }

  if( this->Configuration != NULL )
    {
    delete this->Configuration;
    }

  this->ClearAnalysisTools();

#ifdef ENABLESTATS
  this->FinalizeTimers();
#endif
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
bool CosmologyToolsManager::IsExecutionTimeStep(
        INTEGER tstep, REAL redShift)
{
  bool status = false;

  // STEP 0: Parse the configuration file
  this->StartTimer(tstep,"ReadConfiguration");
  this->ParseConfigurationFile();

  // STEP 1: Create a vector of analysis tools
  this->CreateAnalysisTools();
  this->EndTimer(tstep,"ReadConfiguration");
  this->Barrier();

  // STEP 2: Loop through all analsysis tools and see if anything needs to
  // be executed.
  std::map<std::string,AnalysisTool*>::iterator toolIter;
  toolIter=this->AnalysisTools.begin();
  for( ;toolIter!=this->AnalysisTools.end(); ++toolIter)
    {
    std::string toolName = toolIter->first;
    AnalysisTool *tool   = toolIter->second;
    assert("pre: NULL tool!" && (tool != NULL) );

    PRINTLN(<< tool->GetInformation() );

    if( tool->ShouldExecute(tstep) )
      {
      this->ConfigurationIsRead = true;
      status = true;
      break;
      } // END if the tool should execute

    } // END for all tools

   this->Barrier();
   return( status );
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
    tool = AnalysisToolInstantiator::CreateInstance(
        this->Configuration->GetToolClassInstance(name));
    assert("pre: Cannot create tool instance!" && (tool != NULL) );
    tool->SetCommunicator( this->Communicator );
    }

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
void CosmologyToolsManager::StartTimer(
      INTEGER tstep, std::string eventName)
{

#ifdef ENABLESTATS
  // STEP 0: Get the start time
  REAL startTime = MPI_Wtime();

  // STEP 1: Update event counter
  if( this->EventCounter.find(eventName) == this->EventCounter.end() )
    {
    this->EventCounter[ eventName ] = 1;
    this->GlobalTimers[ eventName ] = 0.0;
    }
  else
    {
    this->EventCounter[ eventName ]++;
    }

  // STEP 2: Store start time for this event at the given timestep
  std::ostringstream oss;
  oss << eventName << "-" << tstep;
  assert("pre: Event should not already exist!" &&
    (this->Timers.find(oss.str()) == this->Timers.end()) );

  this->Timers[oss.str()]=startTime;
#endif

}

//------------------------------------------------------------------------------
void CosmologyToolsManager::EndTimer(INTEGER tstep, std::string eventName)
{

#ifdef ENABLESTATS
  // STEP 0: Get the end time
  REAL endTime = MPI_Wtime();

  // STEP 1: Get the previously stored start time
  std::ostringstream oss;
  oss << eventName << "-" << tstep;
  assert("pre: Event does not exist!" &&
    (this->Timers.find(oss.str()) != this->Timers.end())  );
  REAL startTime = this->Timers[oss.str()];

  // STEP 2: Compute & store elapsed time
  REAL elapsed = endTime-startTime;
  this->Timers[oss.str()] = elapsed;

  // STEP 3: Update global timer across all timers
  assert("pre: Event does not exist in global timers" &&
          (this->GlobalTimers.find(eventName) != this->GlobalTimers.end()));
  this->GlobalTimers[eventName] += elapsed;
#endif

}

//------------------------------------------------------------------------------
void CosmologyToolsManager::GatherVector(
      std::vector<REAL> &snd, std::vector<REAL> &rcv)
{
  // STEP 0: Get the number of items each process will send
  int N = static_cast<int>( snd.size() );

  // STEP 1: Ensure each process sends the same number of items
  int NTOTAL = 0;
  MPI_Reduce(&N,&NTOTAL,1,MPI_INT,MPI_SUM,0,this->Communicator);
  assert("pre: Each proc must send the same number of data" &&
          NTOTAL==this->NumRanks*N);

  // STEP 2: Allocate rcv buffer at the root process
  if( this->Rank==0 )
    {
    rcv.resize(NTOTAL,0.0);
    }

  // STEP 3: Do gather
  MPI_Gather(&snd[0],N,MPI_REAL_T,
      &rcv[0],N,MPI_REAL_T,0,this->Communicator);
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::FinalizeTimers()
{
  // STEP 0: Print out the frequency of each event
  std::map<std::string,INTEGER>::iterator eventCounterIter=
      this->EventCounter.begin();
  for(; eventCounterIter != this->EventCounter.end(); ++eventCounterIter )
    {
    std::string eventName = eventCounterIter->first;
    INTEGER frequency = eventCounterIter->second;
    PRINTLN( << "EVENT: " << eventName << " RUN " << frequency << "times" );
    } // END for all events

  // STEP 2: Pack the discrete timers for this process
  std::vector< REAL > processTimers;
  processTimers.resize( this->Timers.size() );
  std::map<std::string,REAL>::iterator timersIter = this->Timers.begin();
  for(int idx=0; timersIter != this->Timers.end(); ++timersIter, ++idx )
    {
    processTimers[idx] = timersIter->second;
    } // END for all timers in this process

  // STEP 3: Pack the global timers for this process
  std::vector< REAL > processGlobalTimers;
  processGlobalTimers.resize( this->GlobalTimers.size() );
  std::map<std::string,REAL>::iterator globalTimersIter=
      this->GlobalTimers.begin();
  int idx=0;
  for(;globalTimersIter != this->GlobalTimers.end(); ++globalTimersIter,++idx)
    {
    processGlobalTimers[idx] = globalTimersIter->second;
    } // END for all global timers in this process


  // STEP 4: Gather vectors
  std::vector<REAL> DiscreteTimes;
  std::vector<REAL> TotalTimes;
  this->GatherVector(processTimers,DiscreteTimes);
  this->GatherVector(processGlobalTimers,TotalTimes);

  // STEP 4: Process 0 dumps the time statistics
  if( this->Rank == 0 )
    {
    std::map<std::string,REAL>::iterator iter;
    int NEVENTS = 0;

    // Write Total times per process
    std::ofstream ttimesFile;
    ttimesFile.open("TotalTimes.dat");
    assert("pre: Cannot write TotalTimes file!" && (ttimesFile.is_open()));

    ttimesFile << "ProcessID;";
    iter = this->GlobalTimers.begin();
    for(;iter != this->GlobalTimers.end(); ++iter)
      {
      ttimesFile << iter->first << ";";
      } // END for
    ttimesFile << std::endl;

    NEVENTS = this->GlobalTimers.size(); // number of events per proc.
    for(int rnk=0; rnk < this->NumRanks; ++rnk )
      {
      ttimesFile << rnk << ";";
      for(int event=0; event < NEVENTS; ++event )
        {
        ttimesFile << TotalTimes[rnk*NEVENTS+event] << ";";
        } // END for all events
      ttimesFile << std::endl;
      } // END for all ranks

    ttimesFile.close();

    // Write breakdown per time-step file
    std::ofstream dtimesFile;
    dtimesFile.open("DiscreteTimes.dat");
    assert("pre: Cannot write DiscreteTimes file!" && (dtimesFile.is_open()));

    dtimesFile << "ProcessID;";
    iter = this->Timers.begin();
    for(;iter != this->Timers.end(); ++iter)
      {
      dtimesFile << iter->first << ";";
      } // END for
    dtimesFile << std::endl;

    NEVENTS = this->Timers.size(); // number of events per proc.
    for(int rnk=0; rnk < this->NumRanks; ++rnk )
      {
      dtimesFile << rnk << ";";
      for( int event=0; event < NEVENTS; ++event )
        {
        dtimesFile << DiscreteTimes[rnk*NEVENTS+event] << ";";
        } // END for all events
      dtimesFile << std::endl;
      } // END for all ranks

    dtimesFile.close();

    } // END if rank==0

  this->Barrier();
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::CoProcess()
{
  assert("pre: communicator is NULL!" && (this->Communicator!=MPI_COMM_NULL));
  assert("pre: simulation particles data-structure is NULL" &&
         (this->Particles != NULL));


  if( !this->ConfigurationIsRead )
    {
    this->StartTimer(this->Particles->TimeStep,"ReadConfiguration");

    // STEP 0: Parse the configuration file
    this->ParseConfigurationFile();

    // STEP 1: Create a vector of analysis tools
    this->CreateAnalysisTools();
    this->EndTimer(this->Particles->TimeStep,"ReadConfiguration");
    this->Barrier();
    }

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
  this->StartTimer(this->Particles->TimeStep,"TotalInSituExecution");
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
      this->StartTimer(this->Particles->TimeStep,tool->GetName());
      tool->Execute(this->Particles);
      this->EndTimer(this->Particles->TimeStep,tool->GetName());
      PRINTLN(<< "[DONE]");

      if(tool->GetGenerateOutput() == true)
        {
        PRINT(<<"Writing output...");
        this->StartTimer(this->Particles->TimeStep,tool->GetName()+"IO");
        tool->WriteOutput();
        this->EndTimer(this->Particles->TimeStep,tool->GetName()+"IO");
        PRINTLN(<< "[DONE]");
        } // END if should write output

      if(this->Configuration->GetVisualization() &&
          tool->IsVisible())
        {
        this->StartTimer(this->Particles->TimeStep,tool->GetName()+"VIZ");
        // TODO: implement this
        this->EndTimer(this->Particles->TimeStep,tool->GetName()+"VIZ");
        } // END if the tool is visible
      } // END if the tool should execute

    } // END for all tools
  this->EndTimer(this->Particles->TimeStep,"TotalInSituExecution");

  // STEP 3: Synchronize all ranks
  PRINTLN(<< "Finished co-processing @t=" << this->Particles->TimeStep);
  PRINTLN(<< "====================");

  this->ConfigurationIsRead = false;
  this->Barrier();
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::Finalize()
{
//  this->Communicator      = MPI_COMM_NULL;
//  this->ConfigurationFile = "";
  this->ClearAnalysisTools();
}

} // END cosmotk namespace
