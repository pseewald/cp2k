!--------------------------------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations                              ! 
!   Copyright (C) 2000 - 2017  CP2K developers group                                               ! 
!--------------------------------------------------------------------------------------------------!

! ************************************************************************************************** 
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT                               
!> \author AB (01.2018)                                                                              
! ************************************************************************************************** 

MODULE xas_tdp_utils
   
   USE dbcsr_api,                       ONLY : dbcsr_type, dbcsr_p_type, dbcsr_distribution_type, &
                                               dbcsr_get_info, dbcsr_distribution_get, &
                                               dbcsr_distribution_new, dbcsr_create, dbcsr_copy, &
                                               dbcsr_add_on_diag, dbcsr_get_block_p, dbcsr_print, &
                                               dbcsr_put_block, dbcsr_release, dbcsr_finalize,&
                                               dbcsr_distribution_release
   USE kinds,                           ONLY : dp
   USE qs_environment_types,            ONLY : qs_environment_type, get_qs_env
   USE xas_tdp_types,                   ONLY : donor_state_type

#include "./base/base_uses.f90"
   
   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "xas_tdp_utils"

   PUBLIC :: build_xas_tdp_matrix

CONTAINS

! **************************************************************************************************
!> \brief Builds the linear response tddft matrix to be diagonlized for excitation energies. The
!>        donor MOs are restricted to those of the given donor_state    
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \note At this point in time (01.2018), the Kernel is not taken into account, i.e. the matrix has
!>       the form | A  0 | and other subroutines take care of each submatrix
!>                | 0  A |
! **************************************************************************************************
   SUBROUTINE build_xas_tdp_matrix(donor_state, qs_env)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      
      CHARACTER(len=*), PARAMETER :: routineN = 'build_xas_tdp_matrix', &
                                     routineP = moduleN//":"//routineN
      
      INTEGER                                         :: handle
      TYPE(dbcsr_type), POINTER                       :: matrix_a
      TYPE(dbcsr_distribution_type), POINTER          :: submat_dist
      INTEGER, DIMENSION(:), POINTER                  :: submat_blk_size

      CALL timeset(routineN, handle)

      NULLIFY(matrix_a, submat_dist, submat_blk_size)

      CPASSERT(ASSOCIATED(donor_state))
      CPASSERT(ASSOCIATED(qs_env))

!  Get the submatrices distribution and block sizes
      CALL compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env) 

!  For now just testing with matrix A
      ALLOCATE(matrix_a)
      CALL build_matrix_a(matrix_a, submat_dist, submat_blk_size, donor_state, qs_env)
      CALL dbcsr_release(matrix_a)
      DEALLOCATE(matrix_a)

!  Clean-up
      CALL dbcsr_distribution_release(submat_dist)
      DEALLOCATE(submat_dist, submat_blk_size)

      CALL timestop(handle)

   END SUBROUTINE build_xas_tdp_matrix

! **************************************************************************************************
!> \breif This subroutine computes the row/column block structure as well as the dbcsr ditrinution
!>        for the submatrices making up the supr xas tddft matrix. They all share the same
!>        properties, which are base on the replication of the KS matrix
!> \param submat_dist the dbcsr distribution of the submatrices
!> \param submat_blk_size the row/col block sizes for the submatrices 
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
SUBROUTINE compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env)

      TYPE(dbcsr_distribution_type), POINTER          :: submat_dist
      INTEGER, DIMENSION(:), POINTER                  :: submat_blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_submat_dist_and_blk_size', &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: ndo_mo, group, nao, srow_dist, scol_dist, &
                                                         nblk_row, i
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks 
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size, row_dist, col_dist, &
                                                         row_dist_sub, col_dist_sub
      INTEGER, DIMENSION(:,:), POINTER                :: pgrid

      NULLIFY(matrix_ks, dbcsr_dist, row_blk_size, row_dist, col_dist, pgrid, col_dist_sub)
      NULLIFY(row_dist_sub)

!  The submatrices are indexed by M_{pi,qj}, where p,q label basis functions and i,j donor MOs
!  We consider them as being made of blocks of size nao x nao and i,j label those blocks 
!  Each of this block has the distribution and the block sizes of the KS matrix, that we simply
!  replicate ndo_mo*ndo_mo (# donor MOs) times

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices,1)
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, dbcsr_dist=dbcsr_dist) 
      CALL dbcsr_get_info(matrix_ks(1)%matrix, row_blk_size=row_blk_size)
      CALL dbcsr_distribution_get(dbcsr_dist, row_dist=row_dist, col_dist=col_dist, group=group,& 
                                  pgrid=pgrid) 
      nao = SUM(row_blk_size)
      nblk_row = SIZE(row_blk_size)
      srow_dist = SIZE(row_dist)
      scol_dist = SIZE(col_dist)

!  Creation if submatrix block size and col/row distribution
      ALLOCATE(submat_blk_size(ndo_mo*nblk_row))
      ALLOCATE(row_dist_sub(ndo_mo*srow_dist))
      ALLOCATE(col_dist_sub(ndo_mo*scol_dist))

      DO i = 1,ndo_mo
         submat_blk_size((i-1)*nblk_row+1:i*nblk_row) = row_blk_size
         row_dist_sub((i-1)*srow_dist+1:i*srow_dist) = row_dist         
         col_dist_sub((i-1)*scol_dist+1:i*scol_dist) = col_dist            
      END DO

!  Create the submatrix dbcsr distribution
      ALLOCATE(submat_dist)
      CALL dbcsr_distribution_new(submat_dist, group=group, pgrid=pgrid, row_dist=row_dist_sub, &
                                  col_dist=col_dist_sub)            

!  Clean-up
      DEALLOCATE(col_dist_sub, row_dist_sub)

   END SUBROUTINE compute_submat_dist_and_blk_size

! **************************************************************************************************
!> \brief Builds the so called matrix A, which is needed for the full xas_tdp matrix
!>        A_{pi,qj} = F_pq*delta_ij - epsilon_ij*delta_pq
!> \param matrix_a  pointer to a DBCSR matrix containing A
!> \param dist_a the dbcsr distribution of matrix A
!> \param blk_size_a the row/col block sizes of matrix A
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_matrix_a(matrix_a, dist_a, blk_size_a, donor_state, qs_env)

      TYPE(dbcsr_type), POINTER                       :: matrix_a
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(dbcsr_distribution_type), POINTER          :: dist_a
      INTEGER, DIMENSION(:), POINTER                  :: blk_size_a
      
      CHARACTER(len=*), PARAMETER :: routineN = 'build_matrix_a', routineP = moduleN//":"//routineN
      
      INTEGER                                         :: handle, ndo_mo, nao, i, nblk_row,&
                                                         group, j, iblk, jblk
      INTEGER, DIMENSION(:,:), POINTER                :: pgrid
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size 
      REAL(dp), DIMENSION(:,:), POINTER               :: energy_evals     
      TYPE(dbcsr_type)                                :: work_matrix        
      REAL(dp), DIMENSION(:), POINTER                 :: work_block 
      LOGICAL                                         :: found_block 

      CALL timeset(routineN, handle)

      NULLIFY(matrix_ks, dbcsr_dist, row_blk_size, work_block, pgrid, energy_evals)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices,1)
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_get_info(matrix_ks(1)%matrix, row_blk_size=row_blk_size)
      CALL dbcsr_distribution_get(dbcsr_dist, group=group, pgrid=pgrid)
      nao = SUM(row_blk_size)
      nblk_row = SIZE(row_blk_size)
      energy_evals => donor_state%energy_evals

!  Creating the matrix A with the right size and distribution (it is symmetric)
      CALL dbcsr_create(matrix=matrix_a, name="matrix_a", matrix_type="S", dist=dist_a, &
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)

!  Fill the matrix A, looping on blocks the size of matrix_ks, only iterate on upper triangle since
!  the matrix is symmetric
      DO i = 1,ndo_mo
         DO j = i,ndo_mo

            IF (i==j) THEN
               
               CALL dbcsr_copy(work_matrix, matrix_ks(1)%matrix, name="work_matrix")

            ELSE

               CALL dbcsr_create(matrix=work_matrix, row_blk_size=row_blk_size, dist=dbcsr_dist,& 
                                 col_blk_size=row_blk_size, name="work_matrix", matrix_type="S") 

            END IF

!           Substract on the diagonal only if non-zero eval in order to keep nzes to a minimum
            IF (energy_evals(i,j) .NE. 0.0_dp) THEN
               CALL dbcsr_add_on_diag(work_matrix,-energy_evals(i,j))
            END IF
            CALL dbcsr_finalize(work_matrix)

!           Loop over the blocks of the symmetric submatrix and copy one by one into A
            DO iblk = 1,nblk_row
               DO jblk = iblk,nblk_row
                  CALL dbcsr_get_block_p(work_matrix, iblk, jblk, work_block, found_block)

!                 Blocks filled with zero are dicsred
                  IF (found_block) THEN
                     PRINT*, "i,j: ", i,j                    
                     PRINT*, "iblk,jblk", iblk,jblk 
                     PRINT*, "block_size: ", SIZE(work_block)
                     
                     CALL dbcsr_put_block(matrix_a,(i-1)*nblk_row+iblk, (j-1)*nblk_row+jblk, &
                                          work_block)
!                    Need to fill the whole upper aprt of A, by making sure that non-diagonal
!                    sub-blocks are coiped twice in the non-diagonal upper part of A
                     IF ((iblk .NE. jblk) .AND. (i .NE. j)) THEN
                        CALL dbcsr_put_block(matrix_a,(i-1)*nblk_row+jblk, (j-1)*nblk_row+iblk, &
                                             work_block)
                     END IF   
                  END IF

                  NULLIFY(work_block)
               END DO
            END DO

            CALL dbcsr_release(work_matrix)

         END DO
      END DO

!  Test print
      CALL dbcsr_finalize(matrix_a)
      CALL dbcsr_print(matrix_a)
      CALL dbcsr_print(matrix_ks(1)%matrix)

      CALL timestop(handle)

   END SUBROUTINE build_matrix_a

END MODULE xas_tdp_utils

