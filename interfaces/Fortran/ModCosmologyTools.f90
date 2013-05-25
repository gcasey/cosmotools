!-------------------------------------------------------------------------------
! Fortran interface for the CosmologyTools co-processing/co-visualization
! environment.
!
!-------------------------------------------------------------------------------
module ModCosmologyTools
    use mpi
    implicit none

    save

    !...List of available memory layouts
    integer :: ROWMAJOR = 0
    integer :: COLMAJOR = 1
    integer :: NUMMEMORYLAYOUT = 2

    !...List of available halo-finders
    integer :: LANLHALOFINDER=0
    integer :: NUMHALOFINDERS=1

contains

!-------------------------------------------------------------------------------
! Description:
! Initializes the ComsmologyTools co-processing/co-visualization environment
! with the given communicator.
! IN comm -- the MPI communicator
!-------------------------------------------------------------------------------
  subroutine ModCosmologyInit(comm)
      implicit none
      integer :: comm

      call CosmologyFinit( comm )
  end subroutine ModCosmologyInit

!-------------------------------------------------------------------------------
! Description:
! Enable in-situ visualization -- co-processed results are sent over TCP/IP
! to the visualization cluster
!-------------------------------------------------------------------------------
  subroutine ModCosmologyEnableVis()
      implicit none
      call CosmologyEnableVis()
  end subroutine ModCosmologyEnableVis

!-------------------------------------------------------------------------------
! Description:
! Disable in-situ visualization
!-------------------------------------------------------------------------------
  subroutine ModCosmologyDisableVis()
      implicit none
      call CosmologyDisableVis()
  end subroutine ModCosmologyDisableVis

!-------------------------------------------------------------------------------
! Description:
! Sets the memory layout to  use
!
! IN layout -- either 0(ROWMAJOR) or 1(COLMAJOR)
!-------------------------------------------------------------------------------
  subroutine ModCosmologySetMemoryLayout(layout)
      implicit none
      integer :: layout
      call CosmologySetMemoryLayout( layout )
  end subroutine ModCosmologySetMemoryLayout

!-------------------------------------------------------------------------------
! Description:
! Sets the halo tracker frequency
!
! IN frequency -- an integer that specifies how often the tracker is invoked.
!-------------------------------------------------------------------------------
  subroutine ModCosmologySetTrackerFrequency(frequency)
      implicit none
      integer :: frequency
      call CosmologySetTrackerFrequency( frequency )
  end subroutine ModCosmologySetTrackerFrequency

!-------------------------------------------------------------------------------
! Description:
! Sets the halo finder to use
!
! IN haloFinder -- an integer that specifies which halofinder to use.
!-------------------------------------------------------------------------------
  subroutine ModCosmologySetHaloFinder(haloFinder)
    implicit none
    integer :: haloFinder
    call CosmologySetHaloFinder(haloFinder)
  end subroutine ModCosmologySetHaloFinder

!-------------------------------------------------------------------------------
! Description:
! Sets the analysis timesteps the cosmology tools library will execute
!
! IN tsteps -- user-supplied array of timte-stesp
! IN N      -- the number of timesteps
!-------------------------------------------------------------------------------
  subroutine ModCosmologySetAnalysisTimesteps(tsteps,N)
    implicit none
    integer :: tsteps(:), N
    call CosmologySetAnalysisTimeSteps(tsteps,N)
  end subroutine ModCosmologySetAnalysisTimesteps

!-------------------------------------------------------------------------------
! Description:
! Sets the particles at the given timestep/redshift
!
! IN tstep     -- the current discrete timestep
! IN redshift  -- the redshift at the given timestep
! IN x         -- x-component of the particle position vector
! IN y         -- y-component of the particle position vector
! IN z         -- z-component of the particle position vector
! IN vx        -- x-component of the particles velocity vector
! IN vy        -- y-component of the particles velocity vector
! IN vz        -- z-component of the particles velocity vector
! IN mass      -- particle masses
! IN potential -- particle potential
! IN ids       -- the global IDs of each particle
! IN mask      -- particle mask
! IN state     -- particle state
! IN N         -- the total number of particles
!-------------------------------------------------------------------------------
  subroutine ModCosmologySetParticles(tstep,redshift,x,y,z,vx,vy,vz,mass,potential,ids,mask,state,N)
    implicit none
    integer :: tstep, N, ids(:), mask(:), state(:)
    real :: redshift, x(:), y(:), z(:), vx(:), vy(:), vz(:), mass(:), potential(:)

    call CosmologySetParticles(tstep,redshift,x,y,z,vx,vy,vz,mass,potential,ids,mask,state,N)
  end subroutine ModCosmologySetParticles

!-------------------------------------------------------------------------------
! Description:
! Acquires the halotag of each particle. Particles that are not in halos have
! a tag value of -1. Note, this method is to be used after FindHalos or
! TrackHalos is called.
!
! OUT tags -- user-supplied buffer of halo tags, the buffer must be prealloced.
!-------------------------------------------------------------------------------
  subroutine ModCosmologyGetHaloIds(tags)
    implicit none
    integer :: tags(:)
    call CosmologyGetHaloIds( tags )
  end subroutine ModCosmologyGetHaloIds

!-------------------------------------------------------------------------------
! Description:
! Acquires the halotag of each particle. Particles that are not in halos have
! a tag value of -1. Note, this method is to be used after FindHalos or
! TrackHalos is called.
!
! OUT tags -- user-supplied buffer of halo tags, the buffer must be prealloced.
!-------------------------------------------------------------------------------
  subroutine ModCosmologyGetSubHaloIds(tags)
    implicit none
    integer :: tags(:)
    call CosmologyGetSubHaloIds( tags )
  end subroutine ModCosmologyGetSubHaloIds

!-------------------------------------------------------------------------------
! Description:
! Calls the forward halo-tracker to track the halos at the prescribed frequency
! NOTE: this method must be called in combination with ModCosmologySetParticles
!-------------------------------------------------------------------------------
  subroutine ModCosmologyTrackHalos()
    implicit none
    call CosmologyTrackHalos()
  end subroutine ModCosmologyTrackHalos

!-------------------------------------------------------------------------------
! Description:
! Uses the prescribed halo-finder to find the halos of the given particle data
! at the given timestep.
!
! NOTE: this method must be called in combination with ModCosmologySetParticles
!-------------------------------------------------------------------------------
  subroutine ModCosmologyFindHalos()
    implicit none
    call CosmologyFindHalos()
  end subroutine ModCosmologyFindHalos

!-------------------------------------------------------------------------------
! Description:
! Finalizes the CosmologyTools environment
!-------------------------------------------------------------------------------
  subroutine ModCosmologyFinalize()
      implicit none
      call CosmologyFinalize()
  end subroutine ModCosmologyFinalize

end module ModCosmologyTools
