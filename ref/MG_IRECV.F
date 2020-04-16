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

MODULE MG_IRECV_MOD

   ! Purpose
   ! =======
   ! Receiving of neighbor boundary data using MPI non-blocking functionality.

   USE MG_CONSTANTS_MOD
   USE MG_UTILS_MOD
   USE MG_PROFILING_MOD

   IMPLICIT NONE

CONTAINS
   
   SUBROUTINE MG_MPI_IRECV ( GRID, IERR )
   
      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER, INTENT(OUT) :: &
         IERR           ! Return status

      REAL(KIND=MG_REAL), DIMENSION(0:NX+1, 0:NY+1, 0:NZ+1, 1:NVARS), INTENT(INOUT) :: &
         GRID
      ! ------------------
      ! Local Declarations
      ! ------------------

      REAL(KIND=MG_REAL8) ::      &
         TIME_START,              &
         TIME_START_DIR

      INTEGER :: LOCAL
      ! ---------------------
      ! Executable Statements
      ! ---------------------
   
      IERR = 0

      IF ( NUMPES == 1 ) &
         RETURN

#if defined _MG_MPI

      TIME_START = MG_TIMER()

      ! --------------
      ! Post receives:
      ! --------------

      NUM_RECVS = 0

      COUNT_RECV_NORTH = RECV_BUFFER_NORTH_SIZE
      COUNT_RECV_SOUTH = RECV_BUFFER_SOUTH_SIZE
      COUNT_RECV_EAST  = RECV_BUFFER_EAST_SIZE
      COUNT_RECV_WEST  = RECV_BUFFER_WEST_SIZE
      COUNT_RECV_BACK  = RECV_BUFFER_BACK_SIZE
      COUNT_RECV_FRONT = RECV_BUFFER_FRONT_SIZE

      ! Back boundary
   
      IF ( NEIGHBORS(BACK) /= -1 ) THEN
         ! call MPIU_GET_IS_LOCAL(MPI_COMM_MG, NEIGHBORS(BACK), LOCAL)
         ! IF (LOCAL == 0) THEN
         !  ! inter-node comm
         !  INTER_NODE_SIZE = INTER_NODE_SIZE + COUNT_RECV_BACK
         ! ELSE
         !  ! WRITE(*, "('I am in intra comm to ', I0, ', comm size ', I0, ', local = ', I0)") NEIGHBORS(BACK), COUNT_RECV_BACK, LOCAL
         !  INTRA_NODE_SIZE = INTRA_NODE_SIZE + COUNT_RECV_BACK
         ! END IF
         COMM_DATA_SIZE = COMM_DATA_SIZE + NVARS*(NY+2)*(NX+2)
         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( GRID (0, 0, 0, 1), 1, BACKFRONT_DATATYPE, &
                          NEIGHBORS(BACK), MSG_TAGS(FRONT),              &
                          MPI_COMM_MG, MSG_REQS(BACK), IERR )
         MG_PERF%TIME_RECV_Z_PE = MG_PERF%TIME_RECV_Z_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_BACK
         IF ( COUNT_RECV_BACK > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_BACK
         END IF
         IF ( COUNT_RECV_BACK < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_BACK
         END IF

      END IF
   
      ! Front boundary
   
      IF ( NEIGHBORS(FRONT) /= -1 ) THEN
         ! call MPIU_GET_IS_LOCAL(MPI_COMM_MG, NEIGHBORS(FRONT), LOCAL)
         ! IF (LOCAL == 0) THEN
         !  ! inter-node comm
         !  INTER_NODE_SIZE = INTER_NODE_SIZE + COUNT_RECV_FRONT
         ! ELSE
         !  INTRA_NODE_SIZE = INTRA_NODE_SIZE + COUNT_RECV_FRONT
         ! END IF
         COMM_DATA_SIZE = COMM_DATA_SIZE + NVARS*(NY+2)*(NX+2)
         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( GRID (0, 0, NZ+1, 1), 1, BACKFRONT_DATATYPE,   &
                          NEIGHBORS(FRONT), MSG_TAGS(BACK),                &
                          MPI_COMM_MG, MSG_REQS(FRONT), IERR )
         MG_PERF%TIME_RECV_Z_PE = MG_PERF%TIME_RECV_Z_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_FRONT
         IF ( COUNT_RECV_FRONT > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_FRONT
         END IF
         IF ( COUNT_RECV_FRONT < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_FRONT
         END IF

      END IF
   
      ! East boundary
   
      IF ( NEIGHBORS(EAST) /= -1 ) THEN
         ! call MPIU_GET_IS_LOCAL(MPI_COMM_MG, NEIGHBORS(EAST), LOCAL)
         ! IF (LOCAL == 0) THEN
         !  ! inter-node comm
         !  INTER_NODE_SIZE = INTER_NODE_SIZE + COUNT_RECV_EAST
         ! ELSE
         !  INTRA_NODE_SIZE = INTRA_NODE_SIZE + COUNT_RECV_EAST
         ! END IF
         ! TIME_START_DIR = MG_TIMER ()
         COMM_DATA_SIZE = COMM_DATA_SIZE + NVARS*(NZ+2)*(NY+2)
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( GRID (NX+1, 0, 0, 1), 1, EASTWEST_DATATYPE, &
                          NEIGHBORS(EAST), MSG_TAGS(WEST),              &
                          MPI_COMM_MG, MSG_REQS(EAST), IERR )
         MG_PERF%TIME_RECV_X_PE = MG_PERF%TIME_RECV_X_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_EAST
         IF ( COUNT_RECV_EAST > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_EAST
         END IF
         IF ( COUNT_RECV_EAST < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_EAST
         END IF

      END IF
   
      ! West boundary
   
      IF ( NEIGHBORS(WEST) /= -1 ) THEN
         ! call MPIU_GET_IS_LOCAL(MPI_COMM_MG, NEIGHBORS(WEST), LOCAL)
         ! IF (LOCAL == 0) THEN
         !  ! inter-node comm
         !  INTER_NODE_SIZE = INTER_NODE_SIZE + COUNT_RECV_WEST
         ! ELSE
         !  INTRA_NODE_SIZE = INTRA_NODE_SIZE + COUNT_RECV_WEST
         ! END IF
         COMM_DATA_SIZE = COMM_DATA_SIZE + NVARS*(NZ+2)*(NY+2)
         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( GRID (0, 0, 0, 1), 1, EASTWEST_DATATYPE, &
                          NEIGHBORS(WEST), MSG_TAGS(EAST),                &
                          MPI_COMM_MG, MSG_REQS(WEST), IERR )
         MG_PERF%TIME_RECV_X_PE = MG_PERF%TIME_RECV_X_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_WEST
         IF ( COUNT_RECV_WEST > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_WEST
         END IF
         IF ( COUNT_RECV_WEST < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_WEST
         END IF

      END IF
   
      ! North boundary
   
      IF ( NEIGHBORS(NORTH) /= -1 ) THEN
         ! call MPIU_GET_IS_LOCAL(MPI_COMM_MG, NEIGHBORS(NORTH), LOCAL)
         ! IF (LOCAL == 0) THEN
         !  ! inter-node comm
         !  INTER_NODE_SIZE = INTER_NODE_SIZE + COUNT_RECV_NORTH
         ! ELSE
         !  INTRA_NODE_SIZE = INTRA_NODE_SIZE + COUNT_RECV_NORTH
         ! END IF
         COMM_DATA_SIZE = COMM_DATA_SIZE + NVARS*(NZ+2)*(NX+2)
         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( GRID (0, NY+1, 0, 1), 1, NORTHSOUTH_DATATYPE,   &
                          NEIGHBORS(NORTH), MSG_TAGS(SOUTH),                &
                          MPI_COMM_MG, MSG_REQS(NORTH), IERR )
         MG_PERF%TIME_RECV_Y_PE = MG_PERF%TIME_RECV_Y_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_NORTH
         IF ( COUNT_RECV_NORTH > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_NORTH
         END IF
         IF ( COUNT_RECV_NORTH < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_NORTH
         END IF

      END IF
   
      ! South boundary
   
      IF ( NEIGHBORS(SOUTH) /= -1 ) THEN
         ! call MPIU_GET_IS_LOCAL(MPI_COMM_MG, NEIGHBORS(SOUTH), LOCAL)
         ! IF (LOCAL == 0) THEN
         !  ! inter-node comm
         !  INTER_NODE_SIZE = INTER_NODE_SIZE + COUNT_RECV_SOUTH
         ! ELSE
         !  INTRA_NODE_SIZE = INTRA_NODE_SIZE + COUNT_RECV_SOUTH
         ! END IF
         COMM_DATA_SIZE = COMM_DATA_SIZE + NVARS*(NZ+2)*(NX+2)
         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( GRID (0, 0, 0, 1), 1, NORTHSOUTH_DATATYPE, &
                          NEIGHBORS(SOUTH), MSG_TAGS(NORTH),                &
                          MPI_COMM_MG, MSG_REQS(SOUTH), IERR )
         MG_PERF%TIME_RECV_Y_PE = MG_PERF%TIME_RECV_Y_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_SOUTH
         IF ( COUNT_RECV_SOUTH > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_SOUTH
         END IF
         IF ( COUNT_RECV_SOUTH < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_SOUTH
         END IF

      END IF

      MG_PERF%TIME_RECV_PE = MG_PERF%TIME_RECV_PE + MG_TIMER() - TIME_START

      MG_PERF%NUM_RECVS = MG_PERF%NUM_RECVS + NUM_RECVS

#endif _MPI

      RETURN

   END SUBROUTINE MG_MPI_IRECV
   

   SUBROUTINE MG_IRECV ( IERR )
   
      ! ---------------------
      ! Argument Declarations
      ! ---------------------

      INTEGER, INTENT(OUT) :: &
         IERR           ! Return status

      ! ------------------
      ! Local Declarations
      ! ------------------

      REAL(KIND=MG_REAL8) ::      &
         TIME_START,              &
         TIME_START_DIR

      ! ---------------------
      ! Executable Statements
      ! ---------------------
   
      IERR = 0

      IF ( NUMPES == 1 ) &
         RETURN

#if defined _MG_MPI

      TIME_START = MG_TIMER()

      ! --------------
      ! Post receives:
      ! --------------

      NUM_RECVS = 0

      COUNT_RECV_NORTH = RECV_BUFFER_NORTH_SIZE
      COUNT_RECV_SOUTH = RECV_BUFFER_SOUTH_SIZE
      COUNT_RECV_EAST  = RECV_BUFFER_EAST_SIZE
      COUNT_RECV_WEST  = RECV_BUFFER_WEST_SIZE
      COUNT_RECV_BACK  = RECV_BUFFER_BACK_SIZE
      COUNT_RECV_FRONT = RECV_BUFFER_FRONT_SIZE

      ! Back boundary
   
      IF ( NEIGHBORS(BACK) /= -1 ) THEN

         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( RECV_BUFFER_BACK, COUNT_RECV_BACK, MG_MPI_REAL, &
                          NEIGHBORS(BACK), MSG_TAGS(FRONT),              &
                          MPI_COMM_MG, MSG_REQS(BACK), IERR )
         MG_PERF%TIME_RECV_Z_PE = MG_PERF%TIME_RECV_Z_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_BACK
         IF ( COUNT_RECV_BACK > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_BACK
         END IF
         IF ( COUNT_RECV_BACK < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_BACK
         END IF

      END IF
   
      ! Front boundary
   
      IF ( NEIGHBORS(FRONT) /= -1 ) THEN

         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( RECV_BUFFER_FRONT, COUNT_RECV_FRONT, MG_MPI_REAL,   &
                          NEIGHBORS(FRONT), MSG_TAGS(BACK),                &
                          MPI_COMM_MG, MSG_REQS(FRONT), IERR )
         MG_PERF%TIME_RECV_Z_PE = MG_PERF%TIME_RECV_Z_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_FRONT
         IF ( COUNT_RECV_FRONT > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_FRONT
         END IF
         IF ( COUNT_RECV_FRONT < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_FRONT
         END IF

      END IF
   
      ! East boundary
   
      IF ( NEIGHBORS(EAST) /= -1 ) THEN

         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( RECV_BUFFER_EAST, COUNT_RECV_EAST, MG_MPI_REAL, &
                          NEIGHBORS(EAST), MSG_TAGS(WEST),              &
                          MPI_COMM_MG, MSG_REQS(EAST), IERR )
         MG_PERF%TIME_RECV_X_PE = MG_PERF%TIME_RECV_X_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_EAST
         IF ( COUNT_RECV_EAST > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_EAST
         END IF
         IF ( COUNT_RECV_EAST < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_EAST
         END IF

      END IF
   
      ! West boundary
   
      IF ( NEIGHBORS(WEST) /= -1 ) THEN

         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( RECV_BUFFER_WEST, COUNT_RECV_WEST, MG_MPI_REAL, &
                          NEIGHBORS(WEST), MSG_TAGS(EAST),                &
                          MPI_COMM_MG, MSG_REQS(WEST), IERR )
         MG_PERF%TIME_RECV_X_PE = MG_PERF%TIME_RECV_X_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_WEST
         IF ( COUNT_RECV_WEST > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_WEST
         END IF
         IF ( COUNT_RECV_WEST < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_WEST
         END IF

      END IF
   
      ! North boundary
   
      IF ( NEIGHBORS(NORTH) /= -1 ) THEN

         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( RECV_BUFFER_NORTH, COUNT_RECV_NORTH, MG_MPI_REAL,   &
                          NEIGHBORS(NORTH), MSG_TAGS(SOUTH),                &
                          MPI_COMM_MG, MSG_REQS(NORTH), IERR )
         MG_PERF%TIME_RECV_Y_PE = MG_PERF%TIME_RECV_Y_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_NORTH
         IF ( COUNT_RECV_NORTH > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_NORTH
         END IF
         IF ( COUNT_RECV_NORTH < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_NORTH
         END IF

      END IF
   
      ! South boundary
   
      IF ( NEIGHBORS(SOUTH) /= -1 ) THEN

         TIME_START_DIR = MG_TIMER ()
         NUM_RECVS = NUM_RECVS + 1
         CALL MPI_IRECV ( RECV_BUFFER_SOUTH, COUNT_RECV_SOUTH, MG_MPI_REAL, &
                          NEIGHBORS(SOUTH), MSG_TAGS(NORTH),                &
                          MPI_COMM_MG, MSG_REQS(SOUTH), IERR )
         MG_PERF%TIME_RECV_Y_PE = MG_PERF%TIME_RECV_Y_PE + MG_TIMER () - &
                                  TIME_START_DIR

         MG_PERF%RECV_COUNT  = MG_PERF%RECV_COUNT + COUNT_RECV_SOUTH
         IF ( COUNT_RECV_SOUTH > MG_PERF%RECV_COUNT_MAX ) THEN
            MG_PERF%RECV_COUNT_MAX = COUNT_RECV_SOUTH
         END IF
         IF ( COUNT_RECV_SOUTH < MG_PERF%RECV_COUNT_MIN ) THEN
            MG_PERF%RECV_COUNT_MIN = COUNT_RECV_SOUTH
         END IF

      END IF

      MG_PERF%TIME_RECV_PE = MG_PERF%TIME_RECV_PE + MG_TIMER() - TIME_START

      MG_PERF%NUM_RECVS = MG_PERF%NUM_RECVS + NUM_RECVS

#endif _MPI

      RETURN

   END SUBROUTINE MG_IRECV
END MODULE MG_IRECV_MOD
