! ************************************************************************
!
!               miniGhost: stencil computations with boundary exchange.
!                 Copyright (2012) Sandia Corporation
!
! Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
! license for use of this work by or on behalf of the U.S. Government.
!
! This library is free software; you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation; either version 2.1 of the
! License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
! USA
! Questions? Contact Richard F. Barrett (rfbarre@sandia.gov) or
!                    Michael A. Heroux (maherou@sandia.gov)
!
! ************************************************************************

MODULE DRIVER_BSPMA_MOD

CONTAINS
   SUBROUTINE CREATE_BACKFRONT_DATATYPE() 
   USE MG_OPTIONS_MOD
      INTEGER :: IERR
      INTEGER :: SIZE
      call MPI_TYPE_VECTOR(NVARS, (NY+2)*(NX+2), (NZ+2)*(NY+2)*(NX+2), MPI_DOUBLE, BACKFRONT_DATATYPE, IERR)
      IF(IERR /= MPI_SUCCESS) THEN
         WRITE(*,*) "create BACKFRONT_DATATYPE fails"
      END IF

      call MPI_TYPE_COMMIT(BACKFRONT_DATATYPE, IERR)
      IF(IERR /= MPI_SUCCESS) THEN
         WRITE(*,*) "commit BACKFRONT_DATATYPE fails"
      END IF
      call MPI_Type_size(BACKFRONT_DATATYPE, SIZE, IERR)
      !WRITE(*,*) "Create BACKFRONT_DATATYPE sucessfully = ", BACKFRONT_DATATYPE, ", SIZE = ", SIZE
      
   END SUBROUTINE CREATE_BACKFRONT_DATATYPE


   SUBROUTINE CREATE_EASTWEST_DATATYPE()
   USE MG_OPTIONS_MOD
      INTEGER :: IERR
      call MPI_TYPE_VECTOR(NVARS*(NY+2)*(NZ+2), 1, NX+2, MPI_DOUBLE, EASTWEST_DATATYPE, IERR)
      IF(IERR /= MPI_SUCCESS) THEN
         WRITE(*,*) "create EASTWEST_DATATYPE fails"
      END IF

      call MPI_TYPE_COMMIT(EASTWEST_DATATYPE, IERR)
      IF(IERR /= MPI_SUCCESS) THEN
         WRITE(*,*) "commit EASTWEST_DATATYPE fails"
      END IF
      !WRITE(*,*) "Create EASTWEST_DATATYPE sucessfully = ", EASTWEST_DATATYPE
   END SUBROUTINE CREATE_EASTWEST_DATATYPE


   SUBROUTINE CREATE_NORTHSOUTH_DATATYPE()
   USE MG_OPTIONS_MOD
      INTEGER :: IERR
      call MPI_TYPE_VECTOR(NVARS*(NZ+2), NX+2, (NX+2)*(NY+2), MPI_DOUBLE, NORTHSOUTH_DATATYPE, IERR)
      IF(IERR /= MPI_SUCCESS) THEN
         WRITE(*,*) "create NORTHSOUTH_DATATYPE fails"
      END IF

      call MPI_TYPE_COMMIT(NORTHSOUTH_DATATYPE, IERR)
      IF(IERR /= MPI_SUCCESS) THEN
         WRITE(*,*) "commit NORTHSOUTH_DATATYPE fails"
      END IF
      !WRITE(*,*) "Create NORTHSOUTH_DATATYPE sucessfully = ", NORTHSOUTH_DATATYPE
   END SUBROUTINE CREATE_NORTHSOUTH_DATATYPE


   SUBROUTINE MARK_ON_NODE_NEIGHBOR()
   USE MG_CONSTANTS_MOD
   USE MG_OPTIONS_MOD
      INTEGER :: MY_NODE, NEIGHBOR_NODE
      INTEGER :: CORE_PER_NODE = 36 ! Bebop broadwell node

      MY_NODE = MYPE / CORE_PER_NODE
      BREAK_DOWN_TIME(:) = 0.0
      ON_NODE(:) = -1
      ! Back boundary
      IF ( NEIGHBORS(BACK) /= -1 ) THEN
         NEIGHBOR_NODE = NEIGHBORS(BACK) / CORE_PER_NODE
         IF ( NEIGHBOR_NODE == MY_NODE) THEN
            ON_NODE(BACK) = 1
            ON_NODE_XY_CNT = 1
         ELSE
            ON_NODE(BACK) = 0
            OFF_NODE_XY_CNT = 1
         END IF
      END IF
   
      ! Front boundary
      IF ( NEIGHBORS(FRONT) /= -1 ) THEN
         NEIGHBOR_NODE = NEIGHBORS(FRONT) / CORE_PER_NODE
         IF ( NEIGHBOR_NODE == MY_NODE) THEN
            ON_NODE(FRONT) = 1
            ON_NODE_XY_CNT = 1
         ELSE
            ON_NODE(FRONT) = 0
            OFF_NODE_XY_CNT = 1
         END IF
      END IF
   
      ! East boundary
      IF ( NEIGHBORS(EAST) /= -1 ) THEN
         NEIGHBOR_NODE = NEIGHBORS(EAST) / CORE_PER_NODE
         IF ( NEIGHBOR_NODE == MY_NODE) THEN
            ON_NODE(EAST) = 1
            ON_NODE_YZ_CNT = 1
         ELSE
            ON_NODE(EAST) = 0
            OFF_NODE_YZ_CNT = 1
         END IF
      END IF
   
      ! West boundary
      IF ( NEIGHBORS(WEST) /= -1 ) THEN
         NEIGHBOR_NODE = NEIGHBORS(WEST) / CORE_PER_NODE
         IF ( NEIGHBOR_NODE == MY_NODE) THEN
            ON_NODE(WEST) = 1
            ON_NODE_YZ_CNT = 1
         ELSE
            ON_NODE(WEST) = 0
            OFF_NODE_YZ_CNT = 1
         END IF
      END IF
   
      ! North boundary
      IF ( NEIGHBORS(NORTH) /= -1 ) THEN
         NEIGHBOR_NODE = NEIGHBORS(NORTH) / CORE_PER_NODE
         IF ( NEIGHBOR_NODE == MY_NODE) THEN
            ON_NODE(NORTH) = 1
            ON_NODE_XZ_CNT = 1
         ELSE
            ON_NODE(NORTH) = 0
            OFF_NODE_XZ_CNT = 1
         END IF
      END IF
   
      ! South boundary
      IF ( NEIGHBORS(SOUTH) /= -1 ) THEN
         NEIGHBOR_NODE = NEIGHBORS(SOUTH) / CORE_PER_NODE
         IF ( NEIGHBOR_NODE == MY_NODE) THEN
            ON_NODE(SOUTH) = 1
            ON_NODE_XZ_CNT = 1
         ELSE
            ON_NODE(SOUTH) = 0
            OFF_NODE_XZ_CNT = 1
         END IF
      END IF
      ! WRITE(*,*) "Rank ", MYPE, " - ", ON_NODE
   END SUBROUTINE MARK_ON_NODE_NEIGHBOR

   SUBROUTINE DRIVER_BSPMA ( GRID, IERR )

   ! Purpose
   ! =======
   ! Driver for the bulk synchronous parallel with message aggregation
   ! communication strategy. Called by the Fortran driver subroutine MINI_GHOST.

   USE MG_CONSTANTS_MOD
   USE MG_CHECKPOINT_MOD
   USE MG_BSPMA_MOD
   USE MG_BSPMA_DIAGS_MOD
   USE MG_STENCIL_MOD
   USE MG_PROFILING_MOD
   USE MG_SUM_GRID_MOD

   IMPLICIT NONE

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER, INTENT(OUT) :: &
         IERR                       ! Return status

      REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1, 1:NVARS, 2), INTENT(INOUT) :: &
         GRID

      ! ------------------
      ! Local Declarations
      ! ------------------

      LOGICAL :: RESTART_FIRST_PASS=.FALSE.

      INTEGER ::                  &
         IVAR,                    &  ! Counter (over variables)
         ISPIKE,                  &  ! Counter (over NSPIKES)
         NUM_ERRORS,              &  ! Final check of the answer.
         STARTING_TSTEP=1,        &
         STARTING_SPIKE=1,        &
         TSTEP                       ! Time step counter.

      REAL(KIND=MG_REAL) ::       &
         ERROR_ITER,              &  ! Difference between new and old GRIDi sum.
         GSUM                        ! Global SUM across GRIDs.

      REAL(KIND=MG_REAL8) ::      &
         TIME_START,              &  ! Timing variable
         TIME_START_2,            &  ! Timing variable
         TIME_START_ALL              ! Timing variable

      INTEGER :: PID, EXCHG, LOCAL_COMM
      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0

      ! -------------------
      ! Begin time stepping
      ! -------------------

      CALL MG_CP_RESTART_GRID (GRID(0,0,0,1,IN_GRID_ID), RESTART_FIRST_PASS, STARTING_SPIKE, STARTING_TSTEP, GSUM, IERR )

      ! create MPI datatype for packing and unpacking
      call CREATE_BACKFRONT_DATATYPE()
      call CREATE_EASTWEST_DATATYPE()
      call CREATE_NORTHSOUTH_DATATYPE()
      ! PID = GETPID()
      ! WRITE(*,*) "PID ", PID
      ! call sleep(8)
      call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, LOCAL_COMM, IERR)
      call MARK_ON_NODE_NEIGHBOR()

      ! WARM UP
#if defined _MG_MPI
      ! CALL MPI_BARRIER(LOCAL_COMM, IERR)
      IF ( STENCIL == STENCIL_2D5PT .OR. STENCIL == STENCIL_3D7PT ) THEN
         CALL MG_BSPMA ( GRID(0,0,0,1,IN_GRID_ID), IERR )
         CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: MG_BSPMA', TSTEP )
      ELSE IF ( STENCIL == STENCIL_2D9PT .OR. STENCIL == STENCIL_3D27PT ) THEN
         CALL MG_BSPMA_DIAGS ( GRID, IERR )
         CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: MG_BSPMA_DIAG', TSTEP )
      ELSE
         IERR = -1
         CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: Unknown stencil for boundary exchange', STENCIL )
      END IF
#endif   
      
      IF ( STENCIL /= STENCIL_NONE ) THEN
         TIME_START = MG_TIMER()
         DO IVAR = 1, NVARS
            
            call MG_MPI_STENCIL(GRID(0,0,0,1,IN_GRID_ID), GRID(0,0,0,1,OUT_GRID_ID), IVAR, IERR)
            ! CALL MG_STENCIL ( GRID(0,0,0,IVAR), IVAR, IERR )
            CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA:MG_STENCIL', IVAR )
         END DO
      END IF

      EXCHG = IN_GRID_ID 
      IN_GRID_ID = OUT_GRID_ID
      OUT_GRID_ID = EXCHG 

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
      TIME_START_ALL = MG_TIMER()

      DO ISPIKE = STARTING_SPIKE, NSPIKES

         IF (RESTART_FIRST_PASS .EQV. .FALSE.) THEN
            CALL MG_INSERT_SPIKE ( GRID(0,0,0,1,IN_GRID_ID), ISPIKE, IERR )
         END IF

         DO TSTEP = STARTING_TSTEP, NTSTEPS
 100        FORMAT("Current iteration ", I0)
            IF(MYPE == ROOTPE) THEN
               WRITE(*, 100) TSTEP
            END IF
#if defined _MG_MPI
            ! CALL MPI_BARRIER(LOCAL_COMM, IERR)
            TIME_START = MG_TIMER()
            IF ( STENCIL == STENCIL_2D5PT .OR. STENCIL == STENCIL_3D7PT ) THEN
               CALL MG_BSPMA ( GRID(0,0,0,1,IN_GRID_ID), IERR )
               CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: MG_BSPMA', TSTEP )
            ELSE IF ( STENCIL == STENCIL_2D9PT .OR. STENCIL == STENCIL_3D27PT ) THEN
               CALL MG_BSPMA_DIAGS ( GRID, IERR )
               CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: MG_BSPMA_DIAG', TSTEP )
            ELSE
               IERR = -1
               CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: Unknown stencil for boundary exchange', STENCIL )
            END IF
            ! CALL MPI_BARRIER(LOCAL_COMM, IERR)
            MG_PERF%TIME_BSPMA_PE = MG_PERF%TIME_BSPMA_PE + MG_TIMER() - TIME_START
#endif   
            
            IF ( STENCIL /= STENCIL_NONE ) THEN
               TIME_START = MG_TIMER()
               DO IVAR = 1, NVARS
                  
                  call MG_MPI_STENCIL(GRID(0,0,0,1,IN_GRID_ID), GRID(0,0,0,1,OUT_GRID_ID), IVAR, IERR)
                  ! CALL MG_STENCIL ( GRID(0,0,0,IVAR), IVAR, IERR )
                  CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA:MG_STENCIL', IVAR )
                  

                  ! Reduction across GRID option
                  ! IF ( GRIDS_TO_SUM(IVAR) ) THEN
                  !    TIME_START_2 = MG_TIMER ()
                  !    CALL MG_SUM_GRID ( GRID(0,0,0,IVAR,OUT_GRID_ID), IVAR, GSUM, IERR )
                  !    CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA:MG_SUM_GRID', IVAR )
                  !    MG_PERF%TIME_SUMGRID_PE = MG_PERF%TIME_SUMGRID_PE + MG_TIMER () - TIME_START_2

                  !    IF ( MYPE == ROOTPE ) THEN
                  !       ERROR_ITER = ABS ( SOURCE_TOTAL(IVAR) - GSUM ) / SOURCE_TOTAL(IVAR)
                  !       IF ( ERROR_ITER > ERROR_TOL ) THEN
                  !          WRITE(*,99) TSTEP, ISPIKE, IVAR, ERROR_ITER, ERROR_TOL
                  !          IERR = -1
                  !          CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA', IVAR )
                  !       END IF
                  !       IF ( MOD ( TSTEP, REPORT_DIFFUSION ) == 0 ) THEN
                  !          WRITE(*,99) TSTEP, ISPIKE, IVAR, ERROR_ITER, ERROR_TOL
                  !       END IF
                  !    END IF
                  ! END IF
               END DO
               MG_PERF%TIME_STENCIL_PE = MG_PERF%TIME_STENCIL_PE + MG_TIMER() - TIME_START
            END IF

            EXCHG = IN_GRID_ID 
            IN_GRID_ID = OUT_GRID_ID
            OUT_GRID_ID = EXCHG 


            ! CALL MG_CP_CHECKPOINT ( GRID, TSTEP, ISPIKE, GSUM, IERR )

         END DO ! Time step

         STARTING_TSTEP=1
         RESTART_FIRST_PASS=.FALSE.

      END DO    ! SPIKE insertion.

      MG_PERF%TIME_WALL_PE = MG_TIMER() - TIME_START_ALL

      ! Final correctness check
      ! NUM_ERRORS = 0
      ! DO IVAR = 1, NVARS
      !    CALL MG_SUM_GRID ( GRID(0,0,0,IVAR,IN_GRID_ID), IVAR, GSUM, IERR )
      !    CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA:MG_SUM_GRID(Final check)', IVAR )
      !    MG_PERF%TIME_SUMGRID_PE = MG_PERF%TIME_SUMGRID_PE + MG_TIMER () - TIME_START_2

      !    IF ( MYPE == ROOTPE ) THEN
      !       ERROR_ITER = ABS ( SOURCE_TOTAL(IVAR) - GSUM ) / SOURCE_TOTAL(IVAR)
      !       IF ( ERROR_ITER > ERROR_TOL ) THEN
      !          WRITE(*,99) TSTEP, ISPIKE, IVAR, ERROR_ITER, ERROR_TOL
      !          NUM_ERRORS = NUM_ERRORS + 1
      !       END IF
      !    END IF
      ! END DO

      ! IF ( MYPE == ROOTPE ) THEN
      !    IF ( NUM_ERRORS == 0 ) THEN
      !       WRITE(*,*)
      !       WRITE(*,*) 'Computation within error tolerance.'
      !       WRITE(*,*)
      !    END IF
      ! END IF

      call MPI_Type_free(BACKFRONT_DATATYPE, IERR)
      call MPI_Type_free(EASTWEST_DATATYPE, IERR)
      call MPI_Type_free(NORTHSOUTH_DATATYPE, IERR)

      RETURN

 99   FORMAT ( 'Time step ', I4, ' for spike ', I3, ' for variable ', I4, &
               ' the error is ', 1PE12.5 '; error tolerance is ', 1PE12.5, '.' )

   END SUBROUTINE DRIVER_BSPMA

END MODULE DRIVER_BSPMA_MOD
