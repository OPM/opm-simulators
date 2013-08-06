!     AGMG expects an interface to MUMPS named like this
      subroutine DAGMG_MUMPS (id)

!       In here are the official defines for the parameter passed to
!       MUMPS. These must match the corresponding ones in dagmg.f90
        include 'dmumps_struc.h'
        type (DMUMPS_STRUC) :: id

!       Let the system library do all the heavy lifting
        call DMUMPS (id)
      end subroutine DAGMG_MUMPS
