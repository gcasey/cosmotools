#include "mc3.h"

// voronoi analysis
#include "tess.h"
#include "diy.h"
void Tessellate(Particles *particles, int rank, int groupsize, int step, 
		int numsteps);
using namespace std;

// timing
double tess_times[MAX_TIMES];

//MAIN
int main(int argc, char* argv[])
{
  int step0 = 0;

  SimpleTimings::TimerRef t_total = SimpleTimings::getTimer("total");
  SimpleTimings::TimerRef t_init = SimpleTimings::getTimer("init");
  SimpleTimings::TimerRef t_stepr = SimpleTimings::getTimer("stepr");
  SimpleTimings::TimerRef t_step = SimpleTimings::getTimer("step");
  SimpleTimings::TimerRef t_xtra = SimpleTimings::getTimer("xtra");
  SimpleTimings::TimerRef t_sort = SimpleTimings::getTimer("sort");
  SimpleTimings::TimerRef t_map1 = SimpleTimings::getTimer("map1");

  if(argc < 2) {
usage:
    fprintf(stderr,"USAGE: mc3 <indat> <inBase|tfName> <outBase> <INIT|RECORD|BLOCK|COSMO|RESTART> <ROUND_ROBIN|ALL_TO_ALL|ONE_TO_ONE|restart_step>\n");
    fprintf(stderr,"-a <aliveDumpName>   : alive particle dumps\n");
    fprintf(stderr,"-A <smallDumpName>   : alive particle dumps from rank 0 only\n");
    fprintf(stderr,"-r <restartDumpName> : restart particle dumps\n");
    fprintf(stderr,"-f <refreshStepName> : steps for particle refresh\n");
    fprintf(stderr,"-o <analysisdat>     : config file for analysis\n");
    fprintf(stderr,"-s <staticDumpName>  : static time analysis dumps\n");
    fprintf(stderr,"-l <LCUpdateName>    : lightcone time updates\n");
    fprintf(stderr,"-h                   : halo outputs\n");
    fprintf(stderr,"-z                   : skewerz\n");
    fprintf(stderr,"-g                   : final grid output\n");
    fprintf(stderr,"-m                   : initialize MPI_Alltoall\n");
    fprintf(stderr,"-p <pkDumpName>      : P(k) dumps (+ initial, final, restarts)\n");
    fprintf(stderr,"-w                   : use white noise initializer\n");
    fprintf(stderr,"-I                   : use MPI IO for restart files and pseudo-outputs\n");
    fprintf(stderr,"-t <NXxNYxNZ>        : use 3D topology NX by NY by NZ\n");
    fprintf(stderr,"-V <step>            : the step number at which to produce viz. output\n");
    fprintf(stderr,"-v <step>            : the step number at which to produce viz. slab output\n");
    fprintf(stderr,"-C                   : write initial conditions to file(s)\n");
    exit(-1);
  }

  //INITIALIZE MPI
  MPI_Init(&argc, &argv);

  //sort command line options
  MC3Options options(argc, argv);

  string indatName;
  string inBase;
  string outBase;
  string dataType;
  string distributeType;

  int argvi = optind;
  if (!options.newParams()) {
    if (argc < 6) goto usage;
    indatName = argv[argvi++];
    inBase = argv[argvi++];
    outBase = argv[argvi++];
    dataType = argv[argvi++];
    distributeType = argv[argvi++];
  } else {
    string pName = argv[argvi++];
    options.readNewParams(pName);
    inBase = options.inBase();
    outBase = options.outBase();
    dataType = options.dataType();
    distributeType = options.distributeType();
  }

  //starting from restart file
  int restartQ = 0;
  if( dataType == "RESTART") {
    restartQ = 1;
    step0 = atoi(distributeType.c_str());
  }


  //need to convert Mpc to Mpc/h for true cosmo format
  int cosmoFormatQ=0;
  if( dataType == "COSMO") {
    cosmoFormatQ=1;
    dataType = "RECORD";
  }

  //INITIALIZE TOPOLOGY
  int tmpDims[DIMENSION];
  if(options.topologyQ()) {
    char *tmpDimsStr = (char *)options.topologyString().c_str();
    char *tmpDimsTok;
    tmpDimsTok = strtok(tmpDimsStr,"x");
    tmpDims[0] = atoi(tmpDimsTok);
    tmpDimsTok = strtok(NULL,"x");
    tmpDims[1] = atoi(tmpDimsTok);
    tmpDimsTok = strtok(NULL,"x");
    tmpDims[2] = atoi(tmpDimsTok);
    int nnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
    MY_Dims_init_3D(nnodes, DIMENSION, tmpDims);
  }
  Partition::initialize();
  int numranks = Partition::getNumProc();
  int rank = Partition::getMyProc();


  //READ INDAT FILE
  Basedata indat( indatName.c_str(), options.params() );


  //INITIALIZE GEOMETRY
  Domain::initialize(indat);


  //
  MC3Extras *extras = new MC3Extras(options, indat);


  //start some timers
  SimpleTimings::startTimer(t_total);
  SimpleTimings::startTimer(t_init);


  //(OPTIONALLY) INITIALIZE MPI_ALLTOALL
  if(options.initialAlltoallQ()) {
    if(rank==0) {
      printf("\nStarting MPI_Alltoall initialization\n");
      fflush(stdout);
    }
    initial_all_to_all(options.initialAlltoallQ());
    MPI_Barrier(MPI_COMM_WORLD);
  }


  //INSTANTIATE PARTICLES
  Particles particles(indat, options);

    
  if(!restartQ) {
    //START FROM THE BEGINNING
    loadParticles(indat, particles, 
		  inBase, dataType, distributeType, 
		  cosmoFormatQ, options);
  } else {
    //RESTARTING
    if(options.mpiio()) {
      particles.readRestart( create_outName( inBase + "." + MPI_RESTART_SUFFIX, step0).c_str() );
    } else {
      particles.readRestart( create_outName( create_outName( inBase + "." + RESTART_SUFFIX, step0), rank).c_str() );
    }
    step0++;
  }
  MPI_Barrier(MPI_COMM_WORLD);


  //WRITE OUT INITIAL CONDITIONS
  if (dataType == "INIT" && options.writeIC()) {
    if(options.mpiio()) {
      particles.writeAliveHCosmo( (outBase + "." + MPI_ALIVE_SUFFIX + ".IC").c_str(), 1.0/(1.0+indat.zin()));
    } else {
      particles.writeAliveHCosmo( create_outName( outBase + "." + ALIVE_SUFFIX + ".IC", rank).c_str(), 1.0/(1.0+indat.zin()));
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }  


  //MOVE PARTICLES TO CELL BEFORE ALLOCATING FFT MEMORY
  particles.shoveParticles();
  MPI_Barrier(MPI_COMM_WORLD);


  //LOCAL GRID
  GRID_T *rho_arr = particles.field();
  GRID_T *grad_phi_arr = rho_arr;
  GRID_T *phi_arr = rho_arr;


  //LOCAL COPIES OF GEOMETRIC INFO
  int ngla[DIMENSION], nglt[DIMENSION];
  Domain::ng_local_alive(ngla);
  Domain::ng_local_total(nglt);
  int Ngla = Domain::Ng_local_alive();
  int Nglt = Domain::Ng_local_total();
  int ngo = Domain::ng_overload();
  int ng = Domain::ng();


  //INITIALIZE GRID EXCHANGE
  GridExchange gexchange(nglt, ngo, ngo+1);


  //ALLOC POISSON SOLVER AND BUFFERS
  SolverDiscrete *solver = new SolverDiscrete(MPI_COMM_WORLD, ng);
  COMPLEX_T *fft_rho_arr, *fft_grad_phi_arr, *fft_phi_arr;
  poisson_alloc(&fft_rho_arr, &fft_grad_phi_arr, &fft_phi_arr);
  MPI_Barrier(MPI_COMM_WORLD);


  //TIMESTEPPER VARIABLES
  TimeStepper ts(indat.alpha(), indat.ain(), indat.afin(),
		 indat.nsteps(), indat.omega_matter(), indat.w_de(),
                 indat.wa_de());
  int64_t Np_local_alive, Np_global_alive;
  double rho_local_alive, rho_global_alive;


  MPI_Barrier(MPI_COMM_WORLD);


  //P(k) INITIAL
  if(!restartQ) {
    if(rank==0) {
      printf("P(k) initial\n");
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    map2_poisson_forward(particles, solver, rho_arr, fft_rho_arr);
    MPI_Barrier(MPI_COMM_WORLD);

    writePk(solver, outBase + "." + PK_SUFFIX + ".ini");
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0) printf("\n");
    MPI_Barrier(MPI_COMM_WORLD);
  }


  //DONE WITH INITIALIZATION
  SimpleTimings::stopTimerStats(t_init);
  if(rank==0) printf("\n");


  vector<int> Nplav;
  Nplav.reserve(ts.nsteps()+1);
  Nplav.push_back(particles.Np_local_alive());


  //TIMESTEPPER
  SimpleTimings::startTimer(t_stepr);

  //if restart, get timestepper variables up to speed
  for(int step = 0; step < step0; step++) {
    ts.advanceFullStep();
    extras->setStep(step, ts.aa());
  }
  MPI_Barrier(MPI_COMM_WORLD);


  //actual timestepping
  for(int step = step0; step < ts.nsteps(); step++) {
    SimpleTimings::startTimer(t_step);

    if(rank==0) {
      printf("STEP %d, pp = %f, a = %f, z = %f\n",
	     step, ts.pp(), ts.aa(), ts.zz());
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    //HALF STREAM
    SimpleTimings::startTimer(t_map1);
    particles.map1(ts.pp(), ts.tau2(), ts.adot());
    SimpleTimings::stopTimerStats(t_map1);
    MPI_Barrier(MPI_COMM_WORLD);


    //UPDATE TIME
    ts.advanceHalfStep();
    MPI_Barrier(MPI_COMM_WORLD);


    //POISSON FORWARD
    map2_poisson_forward(particles, solver, rho_arr, fft_rho_arr);
    MPI_Barrier(MPI_COMM_WORLD);


    //CHECKSUM DENSITY GRID
    rho_local_alive = sum_rho_alive(rho_arr);
    MPI_Allreduce(&rho_local_alive, &rho_global_alive, 1, 
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


    //POISSON BACKWARD GRADIENT
    map2_poisson_backward_gradient(particles, solver,
				   grad_phi_arr, fft_grad_phi_arr,
				   gexchange, ts, 1.0);
    MPI_Barrier(MPI_COMM_WORLD);


    //UPDATE TIME
    ts.advanceHalfStep();
    MPI_Barrier(MPI_COMM_WORLD);


    //HALF STREAM
    SimpleTimings::startTimer(t_map1);
    particles.map1(ts.pp(), ts.tau2(), ts.adot());
    SimpleTimings::stopTimerStats(t_map1);
    MPI_Barrier(MPI_COMM_WORLD);


    //CHECKSUM PARTICLES
    Np_local_alive = particles.Np_local_alive();
    Nplav.push_back(Np_local_alive);
    MPI_Allreduce(&Np_local_alive, &Np_global_alive, 1, 
		  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if(rank==0) {
      printf( "total alive density   = %f\n",rho_global_alive);
      cout << "total alive particles = " << Np_global_alive << endl;
      cout.flush();
    }
    MPI_Barrier(MPI_COMM_WORLD);




    //EXTRA STUFF
    SimpleTimings::startTimer(t_xtra);
    extras->setStep(step, ts.aa());


    if(extras->extrasStep()) {
      if(rank==0) {
	printf("EXTRAS: ");
	if(extras->vizStep() == step)printf("(viz. output) ");
	if(extras->staticStep())printf("(static output) ");
	if(extras->lcStep())printf("(light cone update) ");
	if(extras->aliveStep())printf("(alive particle output) ");
	if(extras->restartStep())printf("(restart dump) ");
	if(extras->pkStep())printf("(power spectrum) ");
	if(extras->refreshStep())printf("(overload particle refresh)");
	printf("\n");
	fflush(stdout);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    //FORWARD FFT SOLVE
    if(extras->fftfStep())
      map2_poisson_forward(particles, solver, rho_arr, fft_rho_arr);
    MPI_Barrier(MPI_COMM_WORLD);


    //P(k) DUMP
    if(extras->pkStep()) {
      writePk(solver, create_outName(outBase + "." + PK_SUFFIX, step) );
      MPI_Barrier(MPI_COMM_WORLD);
    }


    //BACKWARD POTENTIAL CALCULATION
    if(extras->fftbpotStep())
      map2_poisson_backward_potential(particles, solver, 
				      phi_arr, fft_phi_arr,
				      gexchange);
    MPI_Barrier(MPI_COMM_WORLD);


    //DROP FFT MEMORY
    int fftMemDropped = 0;
    if(extras->particleStep()) {
      fftMemDropped = 1;
      delete solver;
      poisson_free(&fft_rho_arr, &fft_grad_phi_arr, &fft_phi_arr);
      MPI_Barrier(MPI_COMM_WORLD);
    }


    if(extras->particleStep())
      extras->particleExtras(particles, indat, outBase, &Nplav);


    //RE-ALLOC FFT MEMORY
    if(fftMemDropped) {
      particles.shoveParticles();
      solver = new SolverDiscrete(MPI_COMM_WORLD, ng);
      poisson_alloc(&fft_rho_arr, &fft_grad_phi_arr, &fft_phi_arr);
      MPI_Barrier(MPI_COMM_WORLD);
    }

    SimpleTimings::stopTimerStats(t_xtra);
    SimpleTimings::stopTimerStats(t_step);

    if(rank==0) {
      printf("\n");
      fflush(stdout);
    }

//     voronoi analysis (only at a single (final) timestep for now)
//     if (step % 10 == 1 || step == ts.nsteps() - 1)
    if (step == ts.nsteps() - 1)
      Tessellate(&particles, rank, numranks, step, ts.nsteps());

    // test checkpoint / restart
    // if (step == ts.nsteps() - 1) {
    // if (step % 2 == 0) {
    // if (1) {

    //   // get particles
    //   vector<POSVEL_T> xx;
    //   vector<POSVEL_T> yy;
    //   vector<POSVEL_T> zz;
    //   vector<POSVEL_T> vx;
    //   vector<POSVEL_T> vy;
    //   vector<POSVEL_T> vz;
    //   vector<POSVEL_T> phi;
    //   vector<ID_T> pid;
    //   vector<MASK_T> mask;
    //   float anow = 1; // not worried about correct velocities
    //   particles.copyAliveIntoVectors(&xx,&yy,&zz,&vx,&vy,&vz, &phi, &pid, 
    // 				      &mask, anow);

    //   // write restart
    //   RestartIO io_wr(IO_WRITE_RESTART, (char *)"out.rs", MPI_COMM_WORLD);
    //   int64_t num_particles = xx.size();
    //   io_wr.WriteRestart(num_particles, &xx[0], &yy[0], &zz[0], &vx[0], &vy[0],
    // 			 &vz[0], &phi[0], &pid[0], &mask[0]);

    //   // define data to be read back in, allocated by ReadRestart()
    //   int r_num_particles;
    //   float *r_xx, *r_yy, *r_zz;
    //   float *r_vx, *r_vy, *r_vz;
    //   float *r_phi;
    //   int64_t *r_pid;
    //   uint16_t *r_mask;

    //   // read restart
    //   RestartIO io_rr(IO_READ_RESTART, (char *)"out.rs", MPI_COMM_WORLD);
    //   r_num_particles = io_rr.ReadRestart(r_xx, r_yy, r_zz, r_vx, r_vy, r_vz, 
    // 					  r_phi, r_pid, r_mask);

    //   // compare written and read restart data
    //   assert(num_particles == r_num_particles);
    //   for (int i = 0; i < num_particles; i++) {
    // 	assert(pid[i] == r_pid[i]);
    // 	assert(xx[i] == r_xx[i]);
    // 	assert(yy[i] == r_yy[i]);
    // 	assert(zz[i] == r_zz[i]);
    // 	assert(vx[i] == r_vx[i]);
    // 	assert(vy[i] == r_vy[i]);
    // 	assert(vz[i] == r_vz[i]);
    // 	assert(phi[i] == r_phi[i]);
    // 	assert(mask[i] == r_mask[i]);
    //   }

    //   fprintf(stderr, "I/O restart at step %d successful\n", step);

    // }

  } // end timestepper 


  SimpleTimings::stopTimerStats(t_stepr);


  //OUTPUT REST OF LIGHT CONE SKEWERS ACCUMULATED
  if(extras->lightconeQ() && extras->skewerQ()) {
    string outName = create_outName(create_outName(outBase+"."+LC_SKEWER_SUFFIX, ts.nsteps()), Partition::getMyProc());
    (extras->lcskewers())->WriteLocalSkewers(outName.c_str());
    (extras->lcskewers())->ClearSkewers();
  }


  //P(k) FINAL
  if(rank==0) {
    printf("P(k) final\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  map2_poisson_forward(particles, solver, rho_arr, fft_rho_arr);
  MPI_Barrier(MPI_COMM_WORLD);
  writePk(solver, outBase + "." + PK_SUFFIX + ".fin");  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("\n");
  MPI_Barrier(MPI_COMM_WORLD);

  
  //GRID OUTPUT
  if(extras->gridQ()) {
    //CIC ALREADY DONE FOR P(k)
    output_array_alive(rho_arr,create_outName(create_outName(outBase+"."+GRID_SUFFIX,ts.nsteps()),rank).c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);


  delete solver;
  poisson_free(&fft_rho_arr, &fft_grad_phi_arr, &fft_phi_arr);
  MPI_Barrier(MPI_COMM_WORLD);


  SimpleTimings::stopTimer(t_total);
  SimpleTimings::accumStats();


  // Shut down MPI
  Partition::finalize();
  MPI_Finalize();

  return 0;
}
//-------------------------------------------------------------------
//
// prepares info for calling tess and then calls it
//
// particles: hacc particle structure
// rank, groupsize: MPI usual
// step: current timestep
// numsteps: total number of timesteps
//
void Tessellate(Particles *particles, int rank, int groupsize, int step,
		int numsteps) {

  static bool first = true;

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

  if (rank == 0)
    fprintf(stderr, "\n*** Calling voronoi tessellation on step %d ***\n\n",
	    step);

  // get particles
  vector<POSVEL_T> xx;
  vector<POSVEL_T> yy;
  vector<POSVEL_T> zz;
  vector<POSVEL_T> vx;
  vector<POSVEL_T> vy;
  vector<POSVEL_T> vz;
  vector<POSVEL_T> phi;
  vector<ID_T> id;
  vector<MASK_T> mask;
  float anow = 1; // not worried about correct velocities
  particles->copyAliveIntoVectors(&xx,&yy,&zz,&vx,&vy,&vz, &phi, &id, 
				 &mask, anow);

  /* following must be malloc, not new, because tess will realloc them
     during its particle exchange, and it is a C library */
  float **pts = (float **)malloc(sizeof(float*)); // one block per process
  pts[0] = (float *)malloc(3 * xx.size() * sizeof(float));

  int num_pts[1]; // one block per process
  num_pts[0] = xx.size();
  for (int i = 0; i < num_pts[0]; i++) {
    pts[0][3 * i]     = xx[i];
    pts[0][3 * i + 1] = yy[i];
    pts[0][3 * i + 2] = zz[i];
  }

  // first timestep
  if (first) {

    // block bounds
    bb_t bb; // bounds
    float min[DIMENSION], size[DIMENSION];
    Domain::rL_local_alive(size);
    Domain::corner_phys_alive(min);
    for(int i = 0; i < DIMENSION; i++) {
      bb.min[i] = min[i];
      bb.max[i] = min[i] + size[i];
    }

    // decomposition size (number of blocks in each dimension)
    int decomp_size[DIMENSION];
    Partition::getDecompSize(decomp_size);

    // data overall extents
    // assume all blocks are same size (as mine)
    // assume 0,0,0 is the overall data minimum corner
    float data_mins[DIMENSION], data_maxs[DIMENSION];
    for (int i = 0; i < DIMENSION; i++) {
      data_mins[i] = 0.0;
      data_maxs[i] = data_mins[i] + size[i] * decomp_size[i];
    }

    // get neighbors
    // one block per process with 26 neighbors
    int neigh_gids[NUM_OF_NEIGHBORS]; // gids = process ranks in hacc 
    Partition::getNeighbors(neigh_gids);
    gb_t **neighbors = new gb_t*[1];
    neighbors[0] = new gb_t[NUM_OF_NEIGHBORS];
    int num_neighbors[1] = { NUM_OF_NEIGHBORS };
    for (int i = 0; i < NUM_OF_NEIGHBORS; i++) {
      neighbors[0][i].gid       = neigh_gids[i];
      neighbors[0][i].proc      = neigh_gids[i];
      neighbors[0][i].neigh_dir = neigh_dirs[i];
    }

    // guess at cell size and ghost factor
    float cell_size = 1.0;
    float ghost_factor = 2.0;

    // gids are trivial, only one block, my MPI rank
    int gids[1];
    gids[0] = rank;

    tess_init(1, groupsize, gids, &bb, neighbors, num_neighbors, 
	      cell_size, ghost_factor, data_mins, data_maxs, 1, .0001, -1.0, 
	      MPI_COMM_WORLD, tess_times);

    first = false;
    delete[] neighbors[0];
    delete[] neighbors;

  }

  // every timestep

  char buf[256];
  sprintf(buf, "/sandbox/vor-%d.out", step);
  tess(pts, num_pts, buf);

  // last timestep
  if (step == numsteps - 1)
    tess_finalize();

  // cleanup: free instead of delete, was malloc'ed instead of new'ed
  free(pts[0]);
  free(pts);

}
//-------------------------------------------------------------------
