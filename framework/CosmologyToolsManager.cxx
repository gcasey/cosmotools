#include "CosmologyToolsManager.h"

#include "AnalysisTool.h"
#include "AnalysisToolInstantiator.h"
#include "CosmologyToolsConfiguration.h"
#include "SimulationParticles.h"

#include "diy.h"

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

unsigned char neigh_dirs[] = {
  DIY_X0,                   DIY_X1,
  DIY_Y0,                   DIY_Y1,
  DIY_Z0,                   DIY_Z1,
  DIY_X0 | DIY_Y0,          DIY_X1 | DIY_Y1,
  DIY_X0 | DIY_Y1,          DIY_X1 | DIY_Y0,
  DIY_Y0 | DIY_Z0,          DIY_Y1 | DIY_Z1,
  DIY_Y0 | DIY_Z1,          DIY_Y1 | DIY_Z0,
  DIY_Z0 | DIY_X0,          DIY_Z1 | DIY_X1,
  DIY_Z0 | DIY_X1,          DIY_Z1 | DIY_X0,
  DIY_X0 | DIY_Y0 | DIY_Z0, DIY_X1 | DIY_Y1 | DIY_Z1,
  DIY_X0 | DIY_Y0 | DIY_Z1, DIY_X1 | DIY_Y1 | DIY_Z0,
  DIY_X0 | DIY_Y1 | DIY_Z0, DIY_X1 | DIY_Y0 | DIY_Z1,
  DIY_X0 | DIY_Y1 | DIY_Z1, DIY_X1 | DIY_Y0 | DIY_Z0,
};


//
// Neighbors are enumerated so that particles can be attached to the correct
// neighbor, but these pairs must be preserved for the ParticleExchange.
// Every processor should be able to send and receive on every iteration of
// the exchange, so if everyone sends RIGHT and receives LEFT it works
//
// Do not change this pairing order.
//
enum NEIGHBOR
{
  X0,                   // Left face
  X1,                   // Right face

  Y0,                   // Bottom face
  Y1,                   // Top face

  Z0,                   // Front face
  Z1,                   // Back face

  X0_Y0,                // Left   bottom edge
  X1_Y1,                // Right  top    edge

  X0_Y1,                // Left   top    edge
  X1_Y0,                // Right  bottom edge

  Y0_Z0,                // Bottom front  edge
  Y1_Z1,                // Top    back   edge

  Y0_Z1,                // Bottom back   edge
  Y1_Z0,                // Top    front  edge

  Z0_X0,                // Front  left   edge
  Z1_X1,                // Back   right  edge

  Z0_X1,                // Front  right  edge
  Z1_X0,                // Back   left   edge

  X0_Y0_Z0,             // Left  bottom front corner
  X1_Y1_Z1,             // Right top    back  corner

  X0_Y0_Z1,             // Left  bottom back  corner
  X1_Y1_Z0,             // Right top    front corner

  X0_Y1_Z0,             // Left  top    front corner
  X1_Y0_Z1,             // Right bottom back  corner

  X0_Y1_Z1,             // Left  top    back  corner
  X1_Y0_Z0              // Right bottom front corner
};

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
void CosmologyToolsManager::ComputeRankNeighbors(
        int pos[3], int neighbor[26] )
{
  assert("pre: communicator is not cartesian!" &&
            this->IsCartesianCommunicator());

  int xpos = pos[0];
  int ypos = pos[1];
  int zpos = pos[2];

  // Face neighbors
  neighbor[X0] = this->GetRankByPosition(xpos-1, ypos, zpos);
  neighbor[X1] = this->GetRankByPosition(xpos+1, ypos, zpos);
  neighbor[Y0] = this->GetRankByPosition(xpos, ypos-1, zpos);
  neighbor[Y1] = this->GetRankByPosition(xpos, ypos+1, zpos);
  neighbor[Z0] = this->GetRankByPosition(xpos, ypos, zpos-1);
  neighbor[Z1] = this->GetRankByPosition(xpos, ypos, zpos+1);

  // Edge neighbors
  neighbor[X0_Y0] = this->GetRankByPosition(xpos-1, ypos-1, zpos);
  neighbor[X0_Y1] = this->GetRankByPosition(xpos-1, ypos+1, zpos);
  neighbor[X1_Y0] = this->GetRankByPosition(xpos+1, ypos-1, zpos);
  neighbor[X1_Y1] = this->GetRankByPosition(xpos+1, ypos+1, zpos);

  neighbor[Y0_Z0] = this->GetRankByPosition(xpos, ypos-1, zpos-1);
  neighbor[Y0_Z1] = this->GetRankByPosition(xpos, ypos-1, zpos+1);
  neighbor[Y1_Z0] = this->GetRankByPosition(xpos, ypos+1, zpos-1);
  neighbor[Y1_Z1] = this->GetRankByPosition(xpos, ypos+1, zpos+1);

  neighbor[Z0_X0] = this->GetRankByPosition(xpos-1, ypos, zpos-1);
  neighbor[Z0_X1] = this->GetRankByPosition(xpos+1, ypos, zpos-1);
  neighbor[Z1_X0] = this->GetRankByPosition(xpos-1, ypos, zpos+1);
  neighbor[Z1_X1] = this->GetRankByPosition(xpos+1, ypos, zpos+1);

  // Corner neighbors
  neighbor[X0_Y0_Z0] = this->GetRankByPosition(xpos-1, ypos-1, zpos-1);
  neighbor[X1_Y0_Z0] = this->GetRankByPosition(xpos+1, ypos-1, zpos-1);
  neighbor[X0_Y1_Z0] = this->GetRankByPosition(xpos-1, ypos+1, zpos-1);
  neighbor[X1_Y1_Z0] = this->GetRankByPosition(xpos+1, ypos+1, zpos-1);
  neighbor[X0_Y0_Z1] = this->GetRankByPosition(xpos-1, ypos-1, zpos+1);
  neighbor[X1_Y0_Z1] = this->GetRankByPosition(xpos+1, ypos-1, zpos+1);
  neighbor[X0_Y1_Z1] = this->GetRankByPosition(xpos-1, ypos+1, zpos+1);
  neighbor[X1_Y1_Z1] = this->GetRankByPosition(xpos+1, ypos+1, zpos+1);
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::GetBlockBounds(
        int decompSize[3], int pos[3], float min[3], float size[3])
{
  // NOTE: This method essentially computes the following from HACC
  // Domain::rL_local_alive(size);
  // Domain::corner_phys_alive(min);

  float m_grid2phys_pos = this->BoxLength/static_cast<float>(this->NDIM);

  for(int i=0; i < 3; ++i)
    {
    size[i] = this->BoxLength/static_cast<float>(decompSize[i]);
    float m_corner_grid_alive = pos[i]*(this->NDIM/decompSize[i]);
    min[i] = m_corner_grid_alive*m_grid2phys_pos;
    } // END for
}

//------------------------------------------------------------------------------
int CosmologyToolsManager::GetRankByPosition(int i, int j, int k)
{
  assert("pre: communicator is not cartesian!" &&
               this->IsCartesianCommunicator());
  int ijk[3];
  ijk[0]=i; ijk[1]=j; ijk[2]=k;
  int rnk;
  MPI_Cart_rank(this->Communicator,ijk,&rnk);
  return( rnk );
}

//------------------------------------------------------------------------------
bool CosmologyToolsManager::IsCartesianCommunicator()
{
  assert("pre: MPI communicator is NULL!" &&
         (this->Communicator != MPI_COMM_NULL) );

  int topology = 0;
  MPI_Topo_test(this->Communicator,&topology);
  return( (topology==MPI_CART) );
}

//------------------------------------------------------------------------------
void CosmologyToolsManager::SetupDIYDecomposition()
{
  assert("pre: MPI communicator is NULL!" &&
         (this->Communicator != MPI_COMM_NULL) );

  // STEP 0: We are dealing with a 3-D domain and with structured topology, so,
  // each block has 26 neighbors.
  const int DIMENSION        = 3;
  const int NUM_OF_NEIGHBORS = 26;

  // STEP 1: Ensure we are dealing with a cartesian communicator
  // TODO: In the long run we should have a way of dealing with non-cartesian
  // communicators.
  assert("pre: communicator is not cartesian!" &&
          this->IsCartesianCommunicator());
  if( this->IsCartesianCommunicator() )
    {
    std::cerr << "ERROR: cosmotools is expecting a cartesian communicator!\n";
    MPI_Abort(this->Communicator,-1);
    }

  // STEP 2: Get cartesian topology
  int myPosition[DIMENSION];  // the position of this rank
  int decomp_size[DIMENSION]; // the dimensions of the cartesian topology
  int periodicity[DIMENSION]; // periodicity for each direction
  MPI_Cart_get(
      this->Communicator,DIMENSION,decomp_size,periodicity,myPosition);

  // STEP 3: Compute number of blocks per process and total number of blocks
  // Currently, numBlocksPerProcess = 1 and totalNumberOfBlocks is given by
  // numBlocksPerProcess * numRanks.
  int numBlocksPerProcess = 1;
  int totalBlocks = 0;
  int numRanks = 0;
  int rank = 0;
  MPI_Comm_rank(this->Communicator,&rank);
  MPI_Comm_size(this->Communicator,&numRanks);
  totalBlocks = numBlocksPerProcess * numRanks;

  // STEP 4: Get the Block bounds
  bb_t bb;
  float min[DIMENSION], size[DIMENSION];
  this->GetBlockBounds(decomp_size,myPosition,min,size);
  for( int i=0; i < DIMENSION; ++i )
    {
    bb.min[i] = min[i];
    bb.max[i] = min[i] + size[i];
    }

  // STEP 5: data overall extents
  // assume all blocks are same size (as mine)
  // assume 0,0,0 is the overall data minimum corner
  float data_mins[DIMENSION], data_maxs[DIMENSION];
  for (int i = 0; i < DIMENSION; i++) {
    data_mins[i] = 0.0;
    data_maxs[i] = data_mins[i] + size[i] * decomp_size[i];
  }

  // STEP 6: Get neighbors
  int neigh_gids[NUM_OF_NEIGHBORS];
  this->ComputeRankNeighbors(myPosition,neigh_gids);
  gb_t **neighbors = new gb_t*[1];
  neighbors[0] = new gb_t[NUM_OF_NEIGHBORS];
  int num_neighbors[1] = { NUM_OF_NEIGHBORS };
  for (int i = 0; i < NUM_OF_NEIGHBORS; i++)
    {
    neighbors[0][i].gid       = neigh_gids[i];
    neighbors[0][i].proc      = neigh_gids[i];
    neighbors[0][i].neigh_dir = neigh_dirs[i];
    } // END for all neighbors

  // STEP 7: gids are trivial, only one block, my MPI rank
  int gids[1];
  gids[0] = rank;

  // STEP 8: Get wrap
  int wrap = (this->XYZPeriodic)? 1:0;

  // STEP 9: Describe the already decomposed domain in DIY
  DIY_Decomposed(numBlocksPerProcess,gids,&bb,
                  NULL,NULL,NULL,NULL,
                  neighbors,num_neighbors,wrap);

}

//------------------------------------------------------------------------------
void CosmologyToolsManager::Initialize(MPI_Comm comm)
{
  assert("pre: invalid communicator" && (comm != MPI_COMM_NULL) );
  this->Communicator = comm;
  MPI_Comm_size(this->Communicator,&this->NumRanks);
  MPI_Comm_rank(this->Communicator,&this->Rank);

  DIY_Init(3,NULL,1,this->Communicator);
  this->SetupDIYDecomposition();
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
  double startTime = MPI_Wtime();

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
  double endTime = MPI_Wtime();

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
      std::vector<double> &snd, std::vector<double> &rcv)
{
  // STEP 0: Get the number of items each process will send
  int N = static_cast<int>( snd.size() );

  // STEP 1: Ensure each process sends the same number of items
  int NTOTAL = 0;
  MPI_Allreduce(&N,&NTOTAL,1,MPI_INT,MPI_SUM,this->Communicator);
  assert("pre: Each proc must send the same number of data" &&
          NTOTAL==this->NumRanks*N);

  // STEP 2: Allocate rcv buffer at the root process
  if( this->Rank==0 )
    {
    rcv.resize(NTOTAL,0.0);
    }

  // STEP 3: Do gather
  MPI_Gather(&snd[0],N,MPI_DOUBLE,
      &rcv[0],N,MPI_DOUBLE,0,this->Communicator);
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
    PRINTLN( << "EVENT: " << eventName << " RUN " << frequency << " TIMES" );
    } // END for all events

  // STEP 2: Pack the discrete timers for this process
  std::vector< double > processTimers;
  processTimers.resize( this->Timers.size() );
  std::map<std::string,double>::iterator timersIter = this->Timers.begin();
  for(int idx=0; timersIter != this->Timers.end(); ++timersIter, ++idx )
    {
    processTimers[idx] = timersIter->second;
    } // END for all timers in this process

  // STEP 3: Pack the global timers for this process
  std::vector< double > processGlobalTimers;
  processGlobalTimers.resize( this->GlobalTimers.size() );
  std::map<std::string,double>::iterator globalTimersIter=
      this->GlobalTimers.begin();
  int idx=0;
  for(;globalTimersIter != this->GlobalTimers.end(); ++globalTimersIter,++idx)
    {
    processGlobalTimers[idx] = globalTimersIter->second;
    } // END for all global timers in this process


  // STEP 4: Gather vectors
  std::vector<double> DiscreteTimes;
  std::vector<double> TotalTimes;
  this->GatherVector(processTimers,DiscreteTimes);
  this->GatherVector(processGlobalTimers,TotalTimes);

  // STEP 4: Process 0 dumps the time statistics
  if( this->Rank == 0 )
    {
    std::map<std::string,double>::iterator iter;
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
  DIY_Finalize();
}

} // END cosmotk namespace
