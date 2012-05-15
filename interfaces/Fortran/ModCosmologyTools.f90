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
! Enable in-situ visualization -- co-processed results are sent over TCP/IP
! to the visualization cluster
!-------------------------------------------------------------------------------
  subroutine ModCosmologyEnableVis()
      implicit none
      call CosmologyEnableVis()
  end subroutine ModCosmologyEnableVis

!-------------------------------------------------------------------------------
! Disable in-situ visualization
!-------------------------------------------------------------------------------
  subroutine ModCosmologyDisableVis()
      implicit none
      call CosmologyDisableVis()
  end subroutine ModCosmologyDisableVis

!-------------------------------------------------------------------------------
! Sets the memory layout to  use
! IN layout -- either 0(ROWMAJOR) or 1(COLMAJOR)
!-------------------------------------------------------------------------------
  subroutine ModCosmologySetMemoryLayout(layout)
      implicit none
      integer :: layout
      call CosmologySetMemoryLayout( layout )
  end subroutine ModCosmologySetMemoryLayout

!-------------------------------------------------------------------------------
! Sets the halo tracker frequency
! IN frequency -- an integer that specifies how often the tracker is invoked.
!-------------------------------------------------------------------------------
  subroutine ModCosmologySetTrackerFrequency(frequency)
      implicit none
      integer :: frequency
      call CosmologySetTrackerFrequency( frequency )
  end subroutine ModCosmologySetTrackerFrequency

!-------------------------------------------------------------------------------
! Sets the particles at the given timestep/redshift
! IN tstep    -- the current discrete timestep
! IN redshift -- the redshift at the given timestep
! IN x        -- x-component of the particle position vector
! IN y        -- y-component of the particle position vector
! IN z        -- z-component of the particle position vector
! IN vx       -- x-component of the particles velocity vector
! IN vy       -- y-component of the particles velocity vector
! IN vz       -- z-component of the particles velocity vector
! IN ids      -- the global IDs of each particle
! IN N        -- the total number of particles
!-------------------------------------------------------------------------------
  subroutine ModCosmologySetParticles(tstep,redshift,x,y,z,vx,vy,vz,ids,N)
    implicit none
    integer :: tstep, N, ids(:)
    real :: redshift, x(:), y(:), z(:), vx(:), vy(:), vz(:)

    call CosmologySetParticles(tstep,redshift,x,y,z,vx,vy,vz,ids,N)
  end subroutine ModCosmologySetParticles

!-------------------------------------------------------------------------------
! Finalizes the CosmologyTools environment
!-------------------------------------------------------------------------------
  subroutine ModCosmologyFinalize()
      implicit none
      call CosmologyFinalize()
  end subroutine ModCosmologyFinalize

end module ModCosmologyTools
