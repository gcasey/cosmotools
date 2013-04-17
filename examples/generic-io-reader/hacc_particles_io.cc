#include "hacc_particles_io.h"

// C/C++ includes
#include <iostream>
#include <cassert>

// CosmoTools includes
#include "CosmologyToolsDefinitions.h"
#include "GenericIODefinitions.hpp"
#include "GenericIOMPIReader.h"
#include "GenericIOReader.h"
#include "GenericIOUtilities.h"

static cosmotk::GenericIOMPIReader *Reader = NULL;

void hacc_io_init(MPI_Comm *comm)
{
  Reader = new cosmotk::GenericIOMPIReader();
  Reader->SetCommunicator( *comm );
}


//------------------------------------------------------------------------------
void hacc_io_finit(MPI_Comm *fcomm)
{
  MPI_Comm comm = MPI_Comm_f2c(*fcomm);
  hacc_io_init(&comm);
}


//------------------------------------------------------------------------------
void hacc_io_set_file(char *file)
{
  assert("pre: Reader is NULL" && (Reader != NULL) );
  Reader->SetFileName(std::string(file));
}


//------------------------------------------------------------------------------
void hacc_io_get_num_elements(int *N)
{
  assert("pre: Reader is NULL" && (Reader != NULL) );

  Reader->OpenAndReadHeader();
  *N = Reader->GetNumberOfElements();
}


//------------------------------------------------------------------------------
void hacc_io_read_particles(
      float *x, float *y, float *z,
      float *vx, float *vy, float *vz,
      float *phi,
      int64_t *idx)
{
  // Register which variables to read
  Reader->AddVariable("x",x);
  Reader->AddVariable("y",y);
  Reader->AddVariable("z",z);
  Reader->AddVariable("vx",vx);
  Reader->AddVariable("vy",vy);
  Reader->AddVariable("vz",vz);
  Reader->AddVariable("phi",phi);
  Reader->AddVariable("id",idx);

  Reader->ReadData();
}

//------------------------------------------------------------------------------
void hacc_io_finalize()
{
  if( Reader != NULL )
    {
    Reader->Close();
    delete Reader;
    }
}
