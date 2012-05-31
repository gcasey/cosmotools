
#include "ParaViewInSituInterface.h"
#include "SimulationParticles.h"

// VTK includes
#include "vtkLiveInsituLink.h"
#include "vtkProcessModule.h"
#include "vtkSMSession.h"


// C/C++ includes
#include <cassert>


namespace cosmologytools
{

ParaViewInSituInterface::ParaViewInSituInterface()
{
  this->Port          = 2222;
  this->IsInitialized = false;
  this->PVLink        = vtkLiveInsituLink::New();
}

//------------------------------------------------------------------------------
ParaViewInSituInterface::~ParaViewInSituInterface()
{
  this->Finalize();
  if( this->PVLink != NULL )
    {
    this->PVLink->Delete();
    }
}

//------------------------------------------------------------------------------
void ParaViewInSituInterface::Initialize(int port)
{
  this->Port = port;
  this->PVLink->SetInsituPort( this->Port );

  vtkProcessModule *pm = vtkProcessModule::GetProcessModule();
  assert("pre: process module is NULL" && (pm != NULL) );

  vtkSMSession *session = vtkSMSession::SafeDownCast(pm->GetSession());
  assert("pre: PV session is NULL!" && (session != NULL) );

  // TODO: need to handle case where there is no running pvserver to connect
  this->PVLink->SimulationInitialize( session->GetSessionProxyManager() );
  this->IsInitialized = true;
}

//------------------------------------------------------------------------------
void ParaViewInSituInterface::Update(SimulationParticles *particles)
{
  assert("pre: input particles data-structure is NULL" && (particles != NULL));

  if(!this->IsInitialized)
    {
    return;
    }

  // For every new time-step update the simulation state before proceeding
  double redShift = static_cast<double>(particles->RedShift);
  this->PVLink->SimulationUpdate(redShift);

  // TODO: push the particles to a vtkPolyData and create a pipeline consisting
  // of the vtkPolyData. How does PVLink, know about what pipeline to send?

  this->PVLink->SimulationPostProcess(redShift);
}

//------------------------------------------------------------------------------
void ParaViewInSituInterface::Finalize()
{
  this->IsInitialized = false;
  // TODO: finalization code here!
}

} /* namespace cosmologytools */
