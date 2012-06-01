#include "ParaViewInSituInterface.h"
#include "SimulationParticles.h"

// VTK includes
#include "vtkCPDataDescription.h"
#include "vtkCPPipeline.h"
#include "vtkCPProcessor.h"
#include "vtkLiveInsituLink.h"
#include "vtkProcessModule.h"
#include "vtkSMSession.h"

// C/C++ includes
#include <iostream>
#include <cassert>


namespace cosmologytools
{

ParaViewInSituInterface::ParaViewInSituInterface()
{
  this->Port          = 2222;
  this->IsInitialized = false;
  this->PVLink        = vtkLiveInsituLink::New();
  this->PVCoProcessor = vtkCPProcessor::New();
}

//------------------------------------------------------------------------------
ParaViewInSituInterface::~ParaViewInSituInterface()
{
  this->Finalize();

  if( this->PVLink != NULL )
    {
    this->PVLink->Delete();
    }

  if( this->PVCoProcessor != NULL )
    {
    this->PVCoProcessor->Delete();
    }
}

//------------------------------------------------------------------------------
void ParaViewInSituInterface::Initialize(int port)
{
  // STEP 0: Initialize ParaView co-processors
  this->PVCoProcessor->Initialize();

  // STEP 1: Attach empty pipeline to co-processor
// Need to attach a pipeline with a particles source
//  vtkCPPipeline *pipeline = vtkCPPipeline::New();
//  this->PVCoProcessor->AddPipeline( pipeline );
//  pipeline->Delete();

  // STEP 2: Acquire PV session
  vtkProcessModule *pm = vtkProcessModule::GetProcessModule();
  assert("pre: process module is NULL" && (pm != NULL) );

  vtkSMSession *session = vtkSMSession::SafeDownCast(pm->GetSession());
  assert("pre: PV session is NULL!" && (session != NULL) );

  // STEP 3: Setup in-situ PV link on the user-supplied port
  // TODO: Need to deal with Visualization not being available on the
  // other end
  this->Port = port;
  this->PVLink->SetInsituPort( this->Port );
  this->PVLink->SimulationInitialize( session->GetSessionProxyManager() );
  this->IsInitialized = true;
}

//------------------------------------------------------------------------------
void ParaViewInSituInterface::Update(SimulationParticles *particles)
{
  assert("pre: input particles data-structure is NULL" && (particles != NULL));

  // STEP 0: Check if the interface has been innitialized
  if(!this->IsInitialized)
    {
    std::cerr << "WARNING: ParaView in-situ interface ";
    std::cerr << " has not be initialized.\n";
    std::cerr << __FILE__ << ":" << __LINE__;
    return;
    }

  // STEP 1: Notify visualization of new time-step
  // For every new time-step update the simulation state before proceeding
  double redShift = static_cast<double>(particles->RedShift);
  this->PVLink->SimulationUpdate( redShift );

  // STEP 2: Push data on to a pipeline
  vtkCPDataDescription *dataDescriptor = vtkCPDataDescription::New();
  dataDescriptor->AddInput("Particles");
  dataDescriptor->ForceOutputOn();

  // How does PVLink know about what pipeline to send?
  this->PVCoProcessor->GetPipeline(0)->CoProcess( dataDescriptor );

  // STEP 3: Push pipeline to remote visualization process
  this->PVLink->SimulationPostProcess( redShift );
}

//------------------------------------------------------------------------------
void ParaViewInSituInterface::Finalize()
{
  assert("pre: ParaView co-processor is NULL" && (this->PVCoProcessor!=NULL));
  this->IsInitialized = false;
  this->PVCoProcessor->Finalize();
}

} /* namespace cosmologytools */
