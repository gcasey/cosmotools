!-------------------------------------------------------------------------------
! Fortran interface for the CosmologyTools co-processing/co-visualization
! environment.
!
!-------------------------------------------------------------------------------
module ModCosmologyTools
    use mpi
    implicit none

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
! Finalizes the CosmologyTools environment
!-------------------------------------------------------------------------------
  subroutine ModCosmologyFinalize()
      implicit none
      call CosmologyFinalize()
  end subroutine ModCosmologyFinalize

end module ModCosmologyTools
