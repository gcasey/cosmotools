#include "CosmologyTools.h"
#include "CosmologyToolsManager.h"
#include <iostream>
#include <cassert>

static cosmotk::CosmologyToolsManager *CosmoToolsManager = NULL;

void cosmotools_initialize(MPI_Comm *comm)
{
  CosmoToolsManager = new cosmotk::CosmologyToolsManager();
  assert(CosmoToolsManager != NULL);
  CosmoToolsManager->Initialize( *comm );
}

//------------------------------------------------------------------------------
void cosmotools_fortran_intialize(MPI_Fint *fcomm)
{
  MPI_Comm comm = MPI_Comm_f2c(*fcomm);
  cosmotools_initialize(&comm);
}

//------------------------------------------------------------------------------
void cosmotools_set_analysis_config(char *configfile)
{
  assert(CosmoToolsManager != NULL);
  CosmoToolsManager->SetAnalysisConfiguration(configfile);
}

//------------------------------------------------------------------------------
void cosmotools_set_domain_parameters(
        REAL *boxlength, INTEGER *ghostoverlap, INTEGER *NDIM,
        INTEGER *XYZPeriodic)
{
  assert(CosmoToolsManager != NULL);
  if( *XYZPeriodic == 1 )
    {
    CosmoToolsManager->SetDomainParameters(
        *boxlength,*ghostoverlap,*NDIM,true);
    }
  else
    {
    CosmoToolsManager->SetDomainParameters(
        *boxlength,*ghostoverlap,*NDIM,false);
    }

}

//------------------------------------------------------------------------------
void cosmotools_set_particles(
            INTEGER *tstep, REAL *redshift,
            POSVEL_T *px, POSVEL_T *py, POSVEL_T *pz,
            POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz,
            REAL *mass, POTENTIAL_T *potential,
            ID_T *tags, MASK_T *mask, STATUS_T *status,
            ID_T *N)
{
  assert(CosmoToolsManager != NULL);
  CosmoToolsManager->SetParticles(
      *tstep,*redshift,
      px,py,pz,
      vx,vy,vz,
      mass,potential,
      tags, mask, status, *N);
}

//------------------------------------------------------------------------------
void cosmotools_coprocess()
{
  assert(CosmoToolsManager != NULL);
  CosmoToolsManager->CoProcess();
}

//------------------------------------------------------------------------------
void cosmotools_finalize()
{
  assert(CosmoToolsManager != NULL);
  CosmoToolsManager->Finalize();
  delete CosmoToolsManager;
}
