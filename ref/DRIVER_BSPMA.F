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
   !USE MPI

   IMPLICIT NONE

      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER, INTENT(OUT) :: &
         IERR                       ! Return status

      REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1, 1:NVARS), INTENT(INOUT) :: &
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
         STENCIL_TIME,            &
         TIME_START,              &  ! Timing variable
         TIME_START_2,            &  ! Timing variable
         TIME_START_ALL              ! Timing variable
      
      INTEGER :: RANK, NODE, DUMMY, PSIZE, TOTAL_NODE, LOCAL_RANK
      CHARACTER(LEN=64) LOG_FILE
     !  REAL(KIND=MG_REAL8) :: & 
     !     total_time
     !     comm_time

      ! ---------------------
      ! Executable Statements
      ! ---------------------

      IERR = 0

      ! -------------------
      ! Begin time stepping
      ! -------------------   
      CALL MPI_COMM_RANK(MPI_COMM_MG, RANK, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_MG, PSIZE, IERR)
      CALL MG_CP_RESTART_GRID (GRID, RESTART_FIRST_PASS, STARTING_SPIKE, STARTING_TSTEP, GSUM, IERR )
      CALL MPI_BARRIER(MPI_COMM_MG, IERR)
      TIME_START_ALL = MG_TIMER()
      !total_time = 0.0
      !comm_time = 0.0
      DO ISPIKE = STARTING_SPIKE, NSPIKES

         IF (RESTART_FIRST_PASS .EQV. .FALSE.) THEN
            CALL MG_INSERT_SPIKE ( GRID, ISPIKE, IERR )
         END IF

         DO TSTEP = STARTING_TSTEP, NTSTEPS
            IF(MYPE == ROOTPE) THEN
               110 FORMAT("Current Step ", I4, ", Total Step ", I4)
               WRITE(*,110) TSTEP, NTSTEPS
            END IF 
#if defined _MG_MPI
            TIME_START = MG_TIMER()
            IF ( STENCIL == STENCIL_2D5PT .OR. STENCIL == STENCIL_3D7PT ) THEN
               CALL MG_BSPMA ( GRID, IERR )
               CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: MG_BSPMA', TSTEP )
            ELSE IF ( STENCIL == STENCIL_2D9PT .OR. STENCIL == STENCIL_3D27PT ) THEN
               CALL MG_BSPMA_DIAGS ( GRID, IERR )
               CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: MG_BSPMA_DIAG', TSTEP )
            ELSE
               IERR = -1
               CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA: Unknown stencil for boundary exchange', STENCIL )
            END IF
            MG_PERF%TIME_BSPMA_PE = MG_PERF%TIME_BSPMA_PE + MG_TIMER() - TIME_START
#endif

            IF ( STENCIL /= STENCIL_NONE ) THEN
               TIME_START = MG_TIMER()
               DO IVAR = 1, NVARS
                  CALL MG_STENCIL ( GRID(0,0,0,IVAR), IVAR, IERR )
                  CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA:MG_STENCIL', IVAR )
               END DO 
               STENCIL_TIME = MG_TIMER() - TIME_START
               MG_PERF%TIME_STENCIL_PE = MG_PERF%TIME_STENCIL_PE + STENCIL_TIME

            END IF
            !WRITE(*,"(A,I0,A,I0,A,F0.3)") "RANK ", RANK, ", ITER ", TSTEP, ", STENCIL: ", STENCIL_TIME
            !CALL MG_CP_CHECKPOINT ( GRID, TSTEP, ISPIKE, GSUM, IERR )

         END DO ! Time step

         STARTING_TSTEP=1
         RESTART_FIRST_PASS=.FALSE.

      END DO    ! SPIKE insertion.

      MG_PERF%TIME_WALL_PE = MG_TIMER() - TIME_START_ALL

      ! Final correctness check
      NUM_ERRORS = 0

      
      NODE = RANK / 36
      TOTAL_NODE = PSIZE / 36
      LOCAL_RANK = MOD(RANK, 36)
      WRITE(LOG_FILE, "(A, I0, A, I0)") "timelog/NODE", TOTAL_NODE, "-", NODE
      
      IF(RANK == 0) THEN
         OPEN(UNIT=10, FILE=LOG_FILE)
         WRITE(10, 100) RANK, MG_PERF%TIME_BSPMA_PE, MG_PERF%TIME_STENCIL_PE, MG_PERF%TIME_BSPMA_PE + MG_PERF%TIME_STENCIL_PE
         CLOSE(10)
         CALL MPI_SEND(DUMMY, 1, MPI_INT, RANK + 1, 0, MPI_COMM_MG, IERR)
      ELSE
         CALL MPI_RECV(DUMMY, 1, MPI_INT, RANK - 1, 0, MPI_COMM_MG, MPI_STATUS_IGNORE, IERR)
         IF(LOCAL_RANK == 0) THEN
            OPEN(UNIT=10, FILE=LOG_FILE)
         ELSE   
            OPEN(UNIT=10, FILE=LOG_FILE, ACCESS = 'APPEND')
         END IF
         WRITE(10, 100) RANK, MG_PERF%TIME_BSPMA_PE, MG_PERF%TIME_STENCIL_PE, MG_PERF%TIME_BSPMA_PE + MG_PERF%TIME_STENCIL_PE, MG_PERF%TIME_WALL_PE
         CLOSE(10)
         IF(RANK /= PSIZE - 1) THEN
            CALL MPI_SEND(DUMMY, 1, MPI_INT, RANK + 1, 0, MPI_COMM_MG, IERR)
         END IF
      END IF
 100  FORMAT("RANK - ", I0, ", Comm: ", F0.3, ", Stencil: ", F0.3, ", Total: ", F0.3, ", Wall: ", F0.3)


      ! Final correctness check
      NUM_ERRORS = 0
      DO IVAR = 1, NVARS
         CALL MG_SUM_GRID ( GRID(0,0,0,IVAR), IVAR, GSUM, IERR )
         CALL MG_ASSERT ( IERR, 'DRIVER_BSPMA:MG_SUM_GRID(Final check)', IVAR )
         MG_PERF%TIME_SUMGRID_PE = MG_PERF%TIME_SUMGRID_PE + MG_TIMER () - TIME_START_2

         IF ( MYPE == ROOTPE ) THEN
            ERROR_ITER = ABS ( SOURCE_TOTAL(IVAR) - GSUM ) / SOURCE_TOTAL(IVAR)
            IF ( ERROR_ITER > ERROR_TOL ) THEN
               WRITE(*,99) TSTEP, ISPIKE, IVAR, ERROR_ITER, ERROR_TOL
               NUM_ERRORS = NUM_ERRORS + 1
            END IF
         END IF
      END DO
      IF ( MYPE == ROOTPE ) THEN
         IF ( NUM_ERRORS == 0 ) THEN
            WRITE(*,*)
            WRITE(*,*) 'Computation within error tolerance.'
            WRITE(*,*)
         END IF
      END IF

      RETURN

 99   FORMAT ( 'Time step ', I4, ' for spike ', I3, ' for variable ', I4, &
               ' the error is ', 1PE12.5 '; error tolerance is ', 1PE12.5, '.' )

   END SUBROUTINE DRIVER_BSPMA

END MODULE DRIVER_BSPMA_MOD
