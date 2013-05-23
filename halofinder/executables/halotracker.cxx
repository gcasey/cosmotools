/**
 * @brief C++ program to process out-of-core the halos and generate
 * a merger-tree.
 */

// C++ includes
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <vector>

// MPI include
#include <mpi.h>

// CosmologyTools includes
#include "CosmologyToolsMacros.h"
#include "ForwardHaloTracker.h"
#include "GenericIO.h"
#include "GenericIOMPIReader.h"
#include "GenericIOPosixReader.h"
#include "GenericIOReader.h"
#include "Halo.h"
#include "MPIUtilities.h"
#include "TaskTimer.h"

//==============================================================================
// Global variables
//==============================================================================
int rank;
int size;
MPI_Comm comm = MPI_COMM_WORLD;

cosmologytools::ForwardHaloTracker *HaloTracker;

// Command line parameters
std::string DataDir="";
std::string DataPrefix = "";
std::string TimeStepsFile = "";
std::string InDatFile = "";
int MergerTreeThreshold = -1;
bool SkipFofProperties = false;
bool Verify = false;
bool Synthetic = false;
bool NewLayout = false;
bool UsePosix  = false;

std::vector< int > timesteps;
std::vector<cosmotk::Halo> Halos;
std::map<std::string,int> Halo2Idx;

std::map<int,int> NumHalosAtTimeStep;

// These parameters are parsed from indat file
struct sim_parameters_t {
  REAL z_in;      // initial red-shift
  REAL z_fin;     // final red-shift
  INTEGER nsteps; // total number of time-steps
  REAL alpha;     // scale factor
  REAL omega_cdm;
  REAL deut;
  REAL hubble;
  REAL rl;
  INTEGER np;

  // computed values
  REAL p_in;
  REAL p_fin;
  REAL epsilon;
} SimulationParameters;

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

//==============================================================================
// Function prototypes
//==============================================================================
void CartCommInit(MPI_Comm comm);
void MPICart2DIYDecomposition(MPI_Comm comm);
int GetTotalNumberOfHalos();
void GetBlockBounds(
    int decompSize[3],int pos[3],float min[3],float dx[3]);
void ComputeRankNeighbors(int pos[3], int neighbor[26]);
int GetRankByPosition(int i, int j, int k);
void ParseArguments(int argc, char **argv);
void ParseSimulationParameters();
void ReadInAnalysisTimeSteps();
void ReadHalosAtTimeStep(int tstep);
int GetHaloIndex(int tstep,int haloTag);
REAL ComputeRedShift(const int tstep);
void CreateSyntheticHalo(
      const int tstep, const int haloIdx, ID_T start, ID_T end);
void CreateHalosAtTimeStep(const int tstep);
void WriteStatistics();
void GetFileNamesAtTimeStep(
        int tstep, std::string &fofFile, std::string &haloParticlesFile);
cosmotk::GenericIOReader* GetReader();

//------------------------------------------------------------------------------

/**
 * @brief Program main
 * @param argc the argument counter
 * @param argv the argument vector
 * @return rc return code
 */
int main(int argc, char **argv)
{
  // STEP 0: Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  cosmotk::MPIUtilities::Printf(comm,"- Initialize MPI...[DONE]\n");


  // STEP 1: Parse arguments
  ParseArguments(argc,argv);

  // STEP 2: Get cartesian communicator, needed for GenericIO
  CartCommInit(comm);
  cosmotk::MPIUtilities::Printf(
      comm,"- Initialize cartesian communicator...[DONE]\n");


  // STEP 3: Initialize DIY
  DIY_Init(3, NULL, 1, comm);
  cosmotk::MPIUtilities::Printf(comm,"- Initialized DIY...[DONE]\n");
  MPICart2DIYDecomposition(comm);
  cosmotk::MPIUtilities::Printf(
      comm,"- Prescribed decomposition to DIY...[DONE]\n");

  // STEP 4: Parse Simulation parameters
  ParseSimulationParameters();
  cosmotk::MPIUtilities::Printf(
      comm,"- Parse simulation parameters...[DONE]\n");


  // STEP 5: Get time-steps
  ReadInAnalysisTimeSteps();
  cosmotk::MPIUtilities::Printf(
      comm,"- Read in analysis time-steps...[DONE]\n");

  // STEP 6: Initialize timers, and the timers vector.
  // For each timestep, we measure the time to read the halos and the time
  // to track them.
  cosmotk::TaskTimer IOTimer,MergerTreeTimer;
  std::vector<double> timers;
  std::vector<double> walltimers;
  timers.resize(timesteps.size()*2,0.0);
  walltimers.resize(timesteps.size()*2,0.0);

  // STEP 7: Setup vector that holds memory usage information.
  // memUsage[i] holds the number of bytes at timestep i.
  std::vector<int> memUsage;
  memUsage.resize(timesteps.size(),0);

  // STEP 8: Loop through all time-steps and track halos
  HaloTracker = new cosmologytools::ForwardHaloTracker();
  HaloTracker->SetCommunicator( comm );
  HaloTracker->SetMergerTreeThreshold( MergerTreeThreshold );
  for(int t=0; t < timesteps.size(); ++t)
    {
    REAL z = ComputeRedShift(timesteps[t]);
    cosmotk::MPIUtilities::Printf(
        comm,"t=%d, redshift=%f\n", timesteps[t], z);

    cosmotk::MPIUtilities::Printf(comm,"Read halos...");
    IOTimer.StartTimer();
    ReadHalosAtTimeStep( timesteps[t] );
    IOTimer.StopTimer();
    timers[ t*2 ] = IOTimer.GetEllapsedTime();

    MPI_Allreduce(&timers[t*2],&walltimers[t*2],1,MPI_DOUBLE,MPI_MAX,comm);
    cosmotk::MPIUtilities::Printf(comm,"[DONE]\n");

    int numHalos = GetTotalNumberOfHalos();
    cosmotk::MPIUtilities::Printf(comm,"Number of halos=%d\n", numHalos );
    NumHalosAtTimeStep[ timesteps[t] ] = numHalos;

    cosmotk::MPIUtilities::Printf(comm,"Track halos...");
    MergerTreeTimer.StartTimer();
    HaloTracker->TrackHalos(timesteps[t],z,Halos);
    MergerTreeTimer.StopTimer();
    timers[ t*2+1 ] = MergerTreeTimer.GetEllapsedTime();

    MPI_Allreduce(
        &timers[t*2+1],&walltimers[t*2+1],1,MPI_DOUBLE,MPI_MAX,comm);

    cosmotk::MPIUtilities::Printf(comm,"[DONE]\n");

    // Update memory usage statistics
    int haloDataSize = cosmotk::Halo::GetHaloMetadataBytesize();
    memUsage[ t ]    = HaloTracker->GetTotalNumberOfHalos()*haloDataSize;
    memUsage[ t ] +=
        HaloTracker->GetTotalNumberOfHaloParticles()*sizeof(ID_T);
    memUsage[ t ] +=
        HaloTracker->GetHaloEvolutionTree()->GetTotalNumberOfBytes();


    cosmotk::MPIUtilities::Printf(
        comm,"\t - Processed timestep %d/%d SIM TSTEP=%d\n",
              t+1,timesteps.size(),timesteps[t]);
    Halos.clear();
    } // END for all time-step

  // STEP 9: Write the tree
  std::ostringstream oss;
  oss << "mtree-" << rank << ".ascii";
  std::ofstream ofs;
  ofs.open(oss.str().c_str());
  ofs << HaloTracker->GetHaloEvolutionTree()->ToString();
  ofs.close();
  //HaloTracker->WriteMergerTree("MergerTree");

  // STEP 10: Write statistics

  // STEP 10.1: timer statistics, each rank, writes each own
  oss.str("");
  oss.clear();
  oss << "Timing-" << rank << ".dat";
  ofs.open(oss.str().c_str());
  ofs << "I/O;MergerTree\n";
  for(int t=0; t < timesteps.size(); ++t)
    {
  ofs << timers[t*2] << ";" << timers[t*2+1] << std::endl;
    }
  ofs.close();

  // STEP 10.2: memory usage statistics, written by rank 0 only
  if( rank == 0 )
    {
    ofs.open("merger-tree-memory-usage.dat");
    ofs << "NumberOfBytesPerTimeStep\n";
    for(int t=0; t < timesteps.size(); ++t)
      {
      ofs << memUsage[t] <<  std::endl;
      } // END for all timesteps
    ofs.close();

    ofs.open("WallTiming.dat");
    ofs << "I/O;MergerTree\n";
    for(int t=0; t < timesteps.size(); ++t)
      {
      ofs << walltimers[t*2] << ";" << walltimers[t*2+1] << std::endl;
      }
    ofs.close();
    } // END if rank

  // STEP 11: Write statistics
  WriteStatistics();
  MPI_Barrier(comm);

  // STEP 11: Finalize
  delete HaloTracker;
  DIY_Finalize();
  cosmotk::MPIUtilities::Printf(comm,"Finalize DIY...[DONE]\n");
  cosmotk::MPIUtilities::Printf(comm,"Finalize MPI...[DONE]\n");
  MPI_Finalize();
  return 0;
}

//------------------------------------------------------------------------------
void CreateHalosAtTimeStep(const int tstep)
{
  switch( tstep )
    {
    case 0:
      CreateSyntheticHalo(tstep,0,0,5);
      CreateSyntheticHalo(tstep,1,6,10);
      CreateSyntheticHalo(tstep,2,11,20);
      break;
    case 1:
      CreateSyntheticHalo(tstep,0,0,10);
      CreateSyntheticHalo(tstep,1,25,30);
      break;
    case 2:
      CreateSyntheticHalo(tstep,0,0,10);
      CreateSyntheticHalo(tstep,1,25,27);
      CreateSyntheticHalo(tstep,2,28,30);
      break;
    default:
      assert("pre: invalid tstep for synthetic data-set!" && false);
    }
}

//------------------------------------------------------------------------------
cosmotk::GenericIOReader* GetReader()
{
  cosmotk::GenericIOReader* reader = NULL;
  if( UsePosix )
    {
    reader = new cosmotk::GenericIOPosixReader();
    }
  else
    {
    reader = new cosmotk::GenericIOMPIReader();
    }
  return( reader );
}

//------------------------------------------------------------------------------
void CreateSyntheticHalo(
      const int tstep, const int haloIdx, ID_T start, ID_T end)
{
  cosmotk::Halo halo;
  halo.Tag      = haloIdx;
  halo.TimeStep = tstep;
  halo.Redshift = ComputeRedShift( tstep );
  for(ID_T idx=start; idx <= end; ++idx)
    {
    halo.ParticleIds.insert( idx );
    }
  Halos.push_back( halo );
}

//------------------------------------------------------------------------------
int GetTotalNumberOfHalos()
{
  int totHalos = 0;
  int localNumHalos = Halos.size();
  MPI_Reduce(&localNumHalos,&totHalos,1,MPI_INT,MPI_SUM,0,comm);
  return( totHalos );
}

//------------------------------------------------------------------------------
void WriteStatistics()
{
  if( !Verify )
    {
    return;
    }

  std::ofstream ofs;
  if( rank == 0 )
    {
    ofs.open("HalosVsTreeNodes.txt");
    ofs << "Halos;TreeNodes;Zombies\n";
    }
  MPI_Barrier(comm);

  for(int t=0; t < timesteps.size(); ++t)
    {
    int tstep = timesteps[ t ];
    int numHalos = NumHalosAtTimeStep[ tstep ];
    int numNodes = HaloTracker->GetNumberOfMergerTreeNodes(tstep);
    int numZombies = HaloTracker->GetNumberOfMergerTreeZombieNodes(tstep);

    int totHalos = 0;
    int totNodes = 0;
    int totZombies = 0;
    MPI_Reduce(&numHalos,&totHalos,1,MPI_INT,MPI_SUM,0,comm);
    MPI_Reduce(&numNodes,&totNodes,1,MPI_INT,MPI_SUM,0,comm);
    MPI_Reduce(&numZombies,&totZombies,1,MPI_INT,MPI_SUM,0,comm);

    if( rank == 0 )
      {
      ofs << totHalos << ";" << totNodes << ";" << totZombies << std::endl;
      }
    } // END for all timesteps

  if( rank == 0 )
    {
    ofs.close();
    }
  MPI_Barrier(comm);
}

//------------------------------------------------------------------------------
void ParseArguments(int argc, char **argv)
{
  for(int i=1; i < argc; ++i )
    {
    if(strcmp(argv[i],"--timesteps")==0)
      {
      TimeStepsFile = std::string(argv[++i]);
      }
    else if(strcmp(argv[i],"--datadir")==0)
      {
      DataDir = std::string(argv[++i]);
      NewLayout = true;
      }
    else if(strcmp(argv[i],"--prefix")==0)
      {
      DataPrefix = std::string(argv[++i]);
      }
    else if(strcmp(argv[i],"--indat")==0)
      {
      InDatFile = std::string(argv[++i]);
      }
    else if(strcmp(argv[i],"--threshold")==0)
      {
      MergerTreeThreshold = atoi(argv[++i]);
      }
    else if(strcmp(argv[i],"--skip-fof")==0)
      {
      SkipFofProperties = true;
      }
    else if(strcmp(argv[i],"--synthetic")==0)
      {
      Synthetic = true;
      }
    else if(strcmp(argv[i],"--verify")==0)
      {
      Verify = true;
      }
    else if(strcmp(argv[i],"--use-posix")==0)
      {
      UsePosix = true;
      }
    else
      {
      std::cerr << "ERROR: invalid argument " << argv[i] << std::endl;
      MPI_Abort(comm,-1);
      }
    } // END for all arguments

  if( !Synthetic )
    {
    assert("pre: specify [--prefix <prefix>] arg" && (DataPrefix != "") );
    assert("pre: specify [--timesteps <timesteps.dat> arg" &&
            (TimeStepsFile!="") );
    assert("pre: specify [--indat <indat.dat>] arg" && (InDatFile != ""));
    }

  assert("pre: specify [--threshold <t>] arg (NOTE: t > 1)" &&
          (MergerTreeThreshold > 1) );

  MPI_Barrier(comm);
}

//------------------------------------------------------------------------------
void ParseSimulationParameters()
{
  if( Synthetic )
    {
    return;
    }

  // STEP 0: Parse raw data at rank 0 and broadcast to all ranks
  int buffSize  = 0;
  char *buffer = NULL;
  std::string indatfile;
  switch( rank )
    {
    case 0:
        {
        std::ifstream ifs;
        ifs.open(InDatFile.c_str());
        if( !ifs.is_open() )
          {
          std::cerr << "ERROR: Cannot open file: " << InDatFile << std::endl;
          std::cerr << "FILE: " << __FILE__ << std::endl;
          std::cerr << "LINE: " << __LINE__ << std::endl;
          MPI_Abort(comm,-1);
          }
        std::stringstream sstream;
        sstream.str(std::string(""));
        sstream << ifs.rdbuf();
        ifs.close();
        indatfile = sstream.str();

        buffSize = static_cast<int>(indatfile.size()+1);
        MPI_Bcast(&buffSize,1,MPI_INT,0,comm);
        MPI_Bcast(const_cast<char*>(indatfile.c_str()),
                    buffSize,
                    MPI_CHAR,
                    0,
                    comm);
        }
      break;
    default:
      MPI_Bcast(&buffSize,1,MPI_INT,0,comm);

      buffer = new char[buffSize];
      MPI_Bcast(buffer,buffSize,MPI_CHAR,0,comm);
      indatfile = std::string(buffer);
      delete [] buffer;
    } // END switch

  // STEP 1: Setup <key,value> dictionary
  Dictionary IndatParams;
  std::string line = "";
  std::stringstream parser(indatfile);
  while(std::getline(parser,line))
    {
    // Tokenize line
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter< std::vector<std::string>  >(tokens));

    if(tokens.size()==2)
      {
      IndatParams[ tokens[0] ] = tokens[1];
      }
    else
      {
      continue;
      }
    } // END while

  // STEP 2: Get indat arguments
  assert("pre: Did not find Z_IN" &&
   (IndatParams.find("Z_IN") != IndatParams.end()));
  assert("pre: Did not find Z_FIN" &&
   (IndatParams.find("Z_FIN") != IndatParams.end()));
  assert("pre: Did not find N_STEPS" &&
   (IndatParams.find("N_STEPS") != IndatParams.end()));
  assert("pre: Did not find ALPHA" &&
   (IndatParams.find("ALPHA") != IndatParams.end()));
  assert("pre: Did not find OMEGA_CDM" &&
   (IndatParams.find("OMEGA_CDM") != IndatParams.end()));
  assert("pre: Did not find DEUT" &&
   (IndatParams.find("DEUT") != IndatParams.end() ) );
  assert("pre: Did not find HUBBLE" &&
   (IndatParams.find("HUBBLE") != IndatParams.end()));
  assert("pre: Did not find RL" &&
   (IndatParams.find("RL") != IndatParams.end()));
  assert("pre: Did not find NP" &&
   (IndatParams.find("NP") != IndatParams.end()));

  SimulationParameters.z_in = std::atof(IndatParams["Z_IN"].c_str());
  SimulationParameters.z_fin = std::atof(IndatParams["Z_FIN"].c_str());
  SimulationParameters.nsteps = std::atoi(IndatParams["N_STEPS"].c_str());
  SimulationParameters.alpha = std::atof(IndatParams["ALPHA"].c_str());
  SimulationParameters.deut = std::atof(IndatParams["DEUT"].c_str());
  SimulationParameters.hubble = std::atof(IndatParams["HUBBLE"].c_str());
  SimulationParameters.rl = std::atof(IndatParams["RL"].c_str());
  SimulationParameters.np = std::atoi(IndatParams["NP"].c_str());
  SimulationParameters.omega_cdm =
      std::atof(IndatParams["OMEGA_CDM"].c_str());

  // STEP 3: Compute derived quantities
  REAL a_in  = 1.0/(1.0+SimulationParameters.z_in);
  REAL a_fin = 1.0/(1.0+SimulationParameters.z_fin);
  SimulationParameters.p_in  = std::pow(a_in,SimulationParameters.alpha);
  SimulationParameters.p_fin = std::pow(a_fin,SimulationParameters.alpha);
  REAL delta_p = SimulationParameters.p_fin - SimulationParameters.p_in;
  SimulationParameters.epsilon = delta_p/SimulationParameters.nsteps;

  // STEP 4: synchronize with all processes
  MPI_Barrier(comm);
}

//------------------------------------------------------------------------------
REAL ComputeRedShift(const int tstep)
{
  REAL z = 0.0;

  if( Synthetic )
    {
    z = static_cast<REAL>( tstep );
    }
  else
    {
    REAL base =
        SimulationParameters.p_in+(tstep+1)*SimulationParameters.epsilon;
    REAL power = (-1.0)/SimulationParameters.alpha;
    z = static_cast<REAL>( pow(base,power) ) - 1.0;
    }
  return( z );
}

//------------------------------------------------------------------------------
int GetHaloIndex(int tstep,int haloTag)
{
  std::string hashCode = cosmotk::Halo::GetHashCodeForHalo(haloTag,tstep);
  if( Halo2Idx.find(hashCode) != Halo2Idx.end() )
    {
    return( Halo2Idx[hashCode] );
    }
  else
    {
    cosmotk::Halo h;
    h.TimeStep = tstep;
    h.Redshift = ComputeRedShift(tstep);
    h.Tag      = haloTag;
    Halos.push_back( h );
    Halo2Idx[ hashCode ] = Halos.size()-1;
    return( Halos.size()-1 );
    }
}

//------------------------------------------------------------------------------
void GetFileNamesAtTimeStep(
        int tstep, std::string &fofFile, std::string &haloParticlesFile)
{
  std::ostringstream oss;
  oss.clear(); oss.str("");

  if( NewLayout )
    {
    oss << DataDir << "/STEP" << tstep << "/"
        << DataPrefix << "." << tstep << ".fofproperties";
    fofFile = oss.str();

    oss.clear(); oss.str("");
    oss << DataDir << "/STEP" << tstep << "/"
        << DataPrefix << "." << tstep << ".haloparticletags";
    haloParticlesFile = oss.str();
    } // END if
  else
    {
    oss << DataPrefix << "." << tstep << ".fofproperties";
    fofFile = oss.str();

    oss.clear(); oss.str("");
    oss << DataPrefix << "." << tstep << ".haloparticletags";
    haloParticlesFile = oss.str();
    } // END else
}

//------------------------------------------------------------------------------
void ReadHalosAtTimeStep(int tstep)
{
  assert("pre: halos at this tstep must be empty!" && (Halos.size()==0) );

  if( Synthetic )
    {
    CreateHalosAtTimeStep( tstep );
    }
  else
    {
    // STEP 0: Construct file names
    std::string fofPropertiesFile;
    std::string haloParticlesFile;
    GetFileNamesAtTimeStep(tstep,fofPropertiesFile,haloParticlesFile);


    // STEP 1: Open and read FOF properties file
    if( !SkipFofProperties )
      {
      cosmotk::GenericIOReader* FofPropertiesReader = GetReader();
      FofPropertiesReader->SetFileName(fofPropertiesFile);
      FofPropertiesReader->SetCommunicator(comm);
      FofPropertiesReader->OpenAndReadHeader();
      int nfof = FofPropertiesReader->GetNumberOfElements();
      ID_T *haloTags      = new ID_T[nfof];
      POSVEL_T *center_x = new POSVEL_T[nfof];
      POSVEL_T *center_y = new POSVEL_T[nfof];
      POSVEL_T *center_z = new POSVEL_T[nfof];
      POSVEL_T *mcx      = new POSVEL_T[nfof];
      POSVEL_T *mcy      = new POSVEL_T[nfof];
      POSVEL_T *mcz      = new POSVEL_T[nfof];
      POSVEL_T *halo_vx  = new POSVEL_T[nfof];
      POSVEL_T *halo_vy  = new POSVEL_T[nfof];
      POSVEL_T *halo_vz  = new POSVEL_T[nfof];
      POSVEL_T *halomass = new POSVEL_T[nfof];

      FofPropertiesReader->AddVariable("fof_halo_tag", haloTags);
      FofPropertiesReader->AddVariable("fof_halo_center_x", center_x);
      FofPropertiesReader->AddVariable("fof_halo_center_y", center_y);
      FofPropertiesReader->AddVariable("fof_halo_center_z", center_z);
      FofPropertiesReader->AddVariable("fof_halo_mean_x",mcx);
      FofPropertiesReader->AddVariable("fof_halo_mean_y",mcy);
      FofPropertiesReader->AddVariable("fof_halo_mean_z",mcz);
      FofPropertiesReader->AddVariable("fof_halo_mean_vx", halo_vx);
      FofPropertiesReader->AddVariable("fof_halo_mean_vy", halo_vy);
      FofPropertiesReader->AddVariable("fof_halo_mean_vz", halo_vz);
      FofPropertiesReader->AddVariable("fof_halo_mass", halomass);

      FofPropertiesReader->ReadData();

      for( int i=0; i < nfof; ++i )
        {
        int tag = haloTags[i];
        assert("pre: tag should not be negative!" &&
                (tag >= 0) );

//        std::cout << tstep << "\t";
//        std::cout << tag << "\t";
//        std::cout << center_x[i] << ", "
//                  << center_y[i] << ", "
//                  << center_z[i] << "\t";
//        std::cout << halo_vx[i] << ", ";
//        std::cout << halo_vy[i] << ", ";
//        std::cout << halo_vz[i] << std::endl;
//        std::cout.flush();

        int idx = GetHaloIndex(tstep,tag);
        Halos[idx].Center[0] = center_x[i];
        Halos[idx].Center[1] = center_y[i];
        Halos[idx].Center[2] = center_z[i];
        Halos[idx].MeanCenter[0] = mcx[i];
        Halos[idx].MeanCenter[1] = mcy[i];
        Halos[idx].MeanCenter[2] = mcz[i];
        Halos[idx].AverageVelocity[0] = halo_vx[i];
        Halos[idx].AverageVelocity[1] = halo_vy[i];
        Halos[idx].AverageVelocity[2] = halo_vz[i];
        Halos[idx].HaloMass = halomass[i];
        }

      delete [] haloTags;
      delete [] center_x;
      delete [] center_y;
      delete [] center_z;
      delete [] mcx;
      delete [] mcy;
      delete [] mcz;
      delete [] halo_vx;
      delete [] halo_vy;
      delete [] halo_vz;
      delete [] halomass;
      FofPropertiesReader->Close();
      delete FofPropertiesReader;
      } // END if skip FOF properties

    // STEP 2: Open and read HaloParticle tags
    cosmotk::GenericIOReader* HaloParticlesReader = GetReader();
    HaloParticlesReader->SetFileName(haloParticlesFile);
    HaloParticlesReader->SetCommunicator(comm);
    HaloParticlesReader->OpenAndReadHeader();
    int npart = HaloParticlesReader->GetNumberOfElements();
    ID_T *particleIds  = new ID_T[npart];
    ID_T *halo_tags    = new ID_T[npart];

    HaloParticlesReader->AddVariable("id",particleIds);
    HaloParticlesReader->AddVariable("fof_halo_tag",halo_tags);

    HaloParticlesReader->ReadData();

    for(int i=0; i < npart; ++i)
      {
//      std::cout << "tstep=" << tstep << "\t";
//      std::cout << "halo_tags: " << halo_tags[i] << std::endl;
//      std::cout.flush();
      int idx = GetHaloIndex(tstep,halo_tags[i]);
      Halos[idx].ParticleIds.insert(particleIds[i]);
      }

    // STEP 3: Close files
    HaloParticlesReader->Close();
    delete HaloParticlesReader;

    // STEP 4: De-allocate memory
    delete [] halo_tags;
    delete [] particleIds;
    }
}

//------------------------------------------------------------------------------
void ReadInAnalysisTimeSteps()
{
  if( Synthetic )
    {
    timesteps.resize( 3 );
    timesteps[ 0 ] = 0;
    timesteps[ 1 ] = 1;
    timesteps[ 2 ] = 2;
    } // END if synthetic
  else
    {
    int numTimeSteps;
    std::ifstream ifs;
    switch(rank)
      {
      case 0:
        ifs.open(TimeStepsFile.c_str());
        if( !ifs.is_open() )
          {
          std::cerr << "Cannot open file " << TimeStepsFile << std::endl;
          MPI_Abort(comm,-1);
          }

        // Read in time-steps
        for(int tstep; ifs >> tstep; timesteps.push_back(tstep) );

        // Broadcast send total number of time-steps in the file
        numTimeSteps = timesteps.size();
        MPI_Bcast(&numTimeSteps,1,MPI_INTEGER,0,comm);

        // Broadcast send the time-steps
        MPI_Bcast(&timesteps[0],numTimeSteps,MPI_INTEGER,0,comm);
        break;
      default:
        // Broad cast receive the total numer of time-steps
        MPI_Bcast(&numTimeSteps,1,MPI_INTEGER,0,comm);

        // Allocate time-steps vector
        timesteps.resize(numTimeSteps);

        // Broad cast receive all the time-steps
        MPI_Bcast(&timesteps[0],numTimeSteps,MPI_INTEGER,0,comm);
      } // END switch
    } // END else

  // Barrier synchronization
  MPI_Barrier(comm);
}

//------------------------------------------------------------------------------
int GetRankByPosition(int i, int j, int k)
{
  int ijk[3];
  ijk[0]=i; ijk[1]=j; ijk[2]=k;
  int rnk;
  MPI_Cart_rank(comm,ijk,&rnk);
  return( rnk );
}

//------------------------------------------------------------------------------
void GetBlockBounds(
    int decompSize[3], int pos[3], float min[3], float dx[3])
{
  float m_grid2phys_pos =
      SimulationParameters.rl/
      static_cast<float>(SimulationParameters.np);

  for(int i=0; i < 3; ++i)
    {
    dx[i] = SimulationParameters.rl/static_cast<float>(decompSize[i]);
    float m_corner_grid_alive =
        pos[i]*(SimulationParameters.np/decompSize[i]);
    min[i] = m_corner_grid_alive*m_grid2phys_pos;
    } // END for
}

//------------------------------------------------------------------------------
void ComputeRankNeighbors(int pos[3], int neighbor[26])
{

 int xpos = pos[0];
 int ypos = pos[1];
 int zpos = pos[2];

 // Face neighbors
 neighbor[X0] = GetRankByPosition(xpos-1, ypos, zpos);
 neighbor[X1] = GetRankByPosition(xpos+1, ypos, zpos);
 neighbor[Y0] = GetRankByPosition(xpos, ypos-1, zpos);
 neighbor[Y1] = GetRankByPosition(xpos, ypos+1, zpos);
 neighbor[Z0] = GetRankByPosition(xpos, ypos, zpos-1);
 neighbor[Z1] = GetRankByPosition(xpos, ypos, zpos+1);

 // Edge neighbors
 neighbor[X0_Y0] = GetRankByPosition(xpos-1, ypos-1, zpos);
 neighbor[X0_Y1] = GetRankByPosition(xpos-1, ypos+1, zpos);
 neighbor[X1_Y0] = GetRankByPosition(xpos+1, ypos-1, zpos);
 neighbor[X1_Y1] = GetRankByPosition(xpos+1, ypos+1, zpos);

 neighbor[Y0_Z0] = GetRankByPosition(xpos, ypos-1, zpos-1);
 neighbor[Y0_Z1] = GetRankByPosition(xpos, ypos-1, zpos+1);
 neighbor[Y1_Z0] = GetRankByPosition(xpos, ypos+1, zpos-1);
 neighbor[Y1_Z1] = GetRankByPosition(xpos, ypos+1, zpos+1);

 neighbor[Z0_X0] = GetRankByPosition(xpos-1, ypos, zpos-1);
 neighbor[Z0_X1] = GetRankByPosition(xpos+1, ypos, zpos-1);
 neighbor[Z1_X0] = GetRankByPosition(xpos-1, ypos, zpos+1);
 neighbor[Z1_X1] = GetRankByPosition(xpos+1, ypos, zpos+1);

 // Corner neighbors
 neighbor[X0_Y0_Z0] = GetRankByPosition(xpos-1, ypos-1, zpos-1);
 neighbor[X1_Y0_Z0] = GetRankByPosition(xpos+1, ypos-1, zpos-1);
 neighbor[X0_Y1_Z0] = GetRankByPosition(xpos-1, ypos+1, zpos-1);
 neighbor[X1_Y1_Z0] = GetRankByPosition(xpos+1, ypos+1, zpos-1);
 neighbor[X0_Y0_Z1] = GetRankByPosition(xpos-1, ypos-1, zpos+1);
 neighbor[X1_Y0_Z1] = GetRankByPosition(xpos+1, ypos-1, zpos+1);
 neighbor[X0_Y1_Z1] = GetRankByPosition(xpos-1, ypos+1, zpos+1);
 neighbor[X1_Y1_Z1] = GetRankByPosition(xpos+1, ypos+1, zpos+1);
}

//------------------------------------------------------------------------------
void MPICart2DIYDecomposition(MPI_Comm mycomm)
{
  assert("pre: MPI communicator is NULL!" &&
          (mycomm != MPI_COMM_NULL) );

  // STEP 0: We are dealing with a 3-D domain and with structured topology, so,
  // each block has 26 neighbors.
  const int DIMENSION        = 3;
  const int NUM_OF_NEIGHBORS = 26;

  // STEP 1: Get cartesian topology
  int myPosition[3];  // the position of this rank
  int decomp_size[3]; // the dimensions of the cartesian topology
  int periodicity[3]; // periodicity for each direction
  MPI_Cart_get(
     mycomm,DIMENSION,decomp_size,periodicity,myPosition);

  // STEP 2: Compute number of blocks per process and total number of blocks
  // Currently, numBlocksPerProcess = 1 and totalNumberOfBlocks is given by
  // numBlocksPerProcess * numRanks.
  int numBlocksPerProcess = 1;
  int totalBlocks = 0;
  int numRanks = size;
  totalBlocks = numBlocksPerProcess * numRanks;

  // STEP 3: Get the Block bounds
  bb_t bb;
  float min[DIMENSION], dx[DIMENSION];
  GetBlockBounds(decomp_size,myPosition,min,dx);
  for( int i=0; i < DIMENSION; ++i )
   {
   bb.min[i] = min[i];
   bb.max[i] = min[i] + dx[i];
   }

  // STEP 4: data overall extents
  // assume all blocks are same size (as mine)
  // assume 0,0,0 is the overall data minimum corner
  float data_mins[DIMENSION], data_maxs[DIMENSION];
  for (int i = 0; i < DIMENSION; i++) {
    data_mins[i] = 0.0;
    data_maxs[i] = data_mins[i] + dx[i] * decomp_size[i];
  }

  // STEP 5: Get neighbors
  int neigh_gids[NUM_OF_NEIGHBORS];
  ComputeRankNeighbors(myPosition,neigh_gids);
  gb_t **neighbors = new gb_t*[1];
  neighbors[0] = new gb_t[NUM_OF_NEIGHBORS];
  int num_neighbors[1] = { NUM_OF_NEIGHBORS };
  for (int i = 0; i < NUM_OF_NEIGHBORS; i++)
   {
   neighbors[0][i].gid       = neigh_gids[i];
   neighbors[0][i].proc      = neigh_gids[i];
   neighbors[0][i].neigh_dir = neigh_dirs[i];
   } // END for all neighbors

  // STEP 6: gids are trivial, only one block, my MPI rank
  int gids[1];
  gids[0] = rank;

  // STEP 7: Get wrap, i.e., domain is all-periodic
  int wrap = 1;

  // STEP 8: Describe the already decomposed domain in DIY
  DIY_Decomposed(numBlocksPerProcess,gids,&bb,
                 NULL,NULL,NULL,NULL,
                 neighbors,num_neighbors,wrap);
  delete [] neighbors[0];
  delete [] neighbors;
}

//------------------------------------------------------------------------------
void CartCommInit(MPI_Comm mycomm)
{
  int periodic[] = { 1, 1, 1 };
  int reorder = 0;
  int dims[3] = { 0, 0, 0 };

  // Get cartesian dimensions
  MPI_Dims_create(size,3,dims);

  // Create cartesian communicator
  MPI_Comm cartComm;
  MPI_Cart_create(mycomm,3,dims,periodic,reorder,&cartComm);

  // Reset global variables for rank and comm
  MPI_Comm_rank(cartComm,&rank);
  comm=cartComm;
}
