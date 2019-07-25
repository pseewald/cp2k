!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility methods to build 3-center integral tensors of various types.
! **************************************************************************************************

MODULE qs_3c_tensors_types
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE distribution_2d_types,           ONLY: distribution_2d_create_prv => distribution_2d_create,&
                                              distribution_2d_release,&
                                              distribution_2d_type
   USE message_passing,                 ONLY: mp_cart_sub,&
                                              mp_comm_free,&
                                              mp_environ
   USE particle_types,                  ONLY: particle_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_iterator_p_type,&
                                              neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_3c_tensors_types'

   PUBLIC :: distribution_3d_type, neighbor_list_3c_type, neighbor_list_3c_iterator_type, &
             distribution_2d_create, distribution_3d_create, distribution_3d_destroy

   INTEGER, PARAMETER, PUBLIC :: symmetric_none = 0, symmetric_ij = 1, symmetric_jk = 2, symmetrik_ik = 3, symmetric_ijk = 4

   TYPE distribution_3d_type
      TYPE(distribution_2d_type), POINTER :: dist_2d_1 => NULL(), dist_2d_2 => NULL()
      INTEGER :: comm_3d, comm_2d_1, comm_2d_2
      LOGICAL :: owns_comm
   END TYPE distribution_3d_type

   TYPE neighbor_list_3c_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: ij_list => NULL(), jk_list => NULL()
      INTEGER :: sym
      TYPE(distribution_3d_type) :: dist_3d
      LOGICAL :: owns_dist
   END TYPE

   TYPE neighbor_list_3c_iterator_type
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iter_ij => NULL()
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iter_jk => NULL()
      INTEGER                                                    :: iter_level
      TYPE(neighbor_list_3c_type)                                :: ijk_nl
   END TYPE

CONTAINS
! **************************************************************************************************
!> \brief Create a 3d distribution
!> \param dist_3d 3d distribution object
!> \param dist1 distribution vector along 1st process grid dimension
!> \param dist2 distribution vector along 2nd process grid dimension
!> \param dist3 distribution vector along 3rd process grid dimension
!> \param nkind ...
!> \param particle_set ...
!> \param mp_comm_3d MPI communicator with a 3d cartesian topology
!> \param own_comm Whether mp_comm_3d should be owned by dist_3d (default false)
! **************************************************************************************************
   SUBROUTINE distribution_3d_create(dist_3d, dist1, dist2, dist3, nkind, particle_set, mp_comm_3d, own_comm)
      TYPE(distribution_3d_type)                         :: dist_3d
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist1, dist2, dist3
      INTEGER, INTENT(IN)                                :: nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: mp_comm_3d
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_comm

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_3d_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: comm_2d_1, comm_2d_2, handle, nproc_1, &
                                                            nproc_2
      INTEGER, DIMENSION(2)                              :: mp_coor_1, mp_coor_2, mp_dims_1, &
                                                            mp_dims_2

      CALL timeset(routineN, handle)

      IF (PRESENT(own_comm)) THEN
         IF (own_comm) dist_3d%comm_3d = mp_comm_3d
         dist_3d%owns_comm = own_comm
      ELSE
         dist_3d%owns_comm = .FALSE.
      ENDIF

      CALL mp_cart_sub(mp_comm_3d, [.TRUE., .TRUE., .FALSE.], comm_2d_1)
      CALL mp_cart_sub(mp_comm_3d, [.FALSE., .TRUE., .TRUE.], comm_2d_2)

      CALL mp_environ(nproc_1, mp_dims_1, mp_coor_1, comm_2d_1)
      CALL mp_environ(nproc_2, mp_dims_2, mp_coor_2, comm_2d_2)

      CPASSERT(mp_coor_1(2) == mp_coor_2(1))

      CALL distribution_2d_create(dist_3d%dist_2d_1, dist1, dist2, nkind, particle_set, comm_2d_1)
      CALL distribution_2d_create(dist_3d%dist_2d_2, dist2, dist3, nkind, particle_set, comm_2d_2)

      dist_3d%comm_2d_1 = comm_2d_1
      dist_3d%comm_2d_2 = comm_2d_2

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy a 3d distribution
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE distribution_3d_destroy(dist)
      TYPE(distribution_3d_type)                         :: dist

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_3d_destroy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL distribution_2d_release(dist%dist_2d_1)
      CALL distribution_2d_release(dist%dist_2d_2)
      CALL mp_comm_free(dist%comm_2d_1)
      CALL mp_comm_free(dist%comm_2d_2)
      IF (dist%owns_comm) CALL mp_comm_free(dist%comm_3d)

      NULLIFY (dist%dist_2d_1, dist%dist_2d_2)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create a 2d distribution. This mainly wraps distribution_2d_create
!>        for consistency with distribution_3d_create.
!> \param dist_2d 2d distribution object
!> \param dist1 distribution vector along 1st process grid dimension
!> \param dist2 distribution vector along 2nd process grid dimension
!> \param nkind ...
!> \param particle_set ...
!> \param mp_comm_2d MPI communicator with a 3d cartesian topology
! **************************************************************************************************
   SUBROUTINE distribution_2d_create(dist_2d, dist1, dist2, nkind, particle_set, mp_comm_2d)
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist1, dist2
      INTEGER, INTENT(IN)                                :: nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: mp_comm_2d

      INTEGER                                            :: iatom, ikind, n, natom, nproc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nparticle_local_col, nparticle_local_row
      INTEGER, DIMENSION(2)                              :: mp_coor, mp_dims
      INTEGER, DIMENSION(:, :), POINTER                  :: dist1_prv, dist2_prv
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER        :: local_particle_col, local_particle_row
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY (blacs_env, local_particle_col, local_particle_row, para_env)

      CALL mp_environ(nproc, mp_dims, mp_coor, mp_comm_2d)
      CALL cp_para_env_create(para_env, mp_comm_2d, owns_group=.FALSE.)
      CALL cp_blacs_env_create(blacs_env, para_env, &
                               grid_2d=mp_dims)

      CPASSERT(blacs_env%mepos(1) == mp_coor(1))
      CPASSERT(blacs_env%mepos(2) == mp_coor(2))

      CALL cp_para_env_release(para_env)

      natom = SIZE(particle_set)
      ALLOCATE (dist1_prv(natom, 2), dist2_prv(natom, 2))
      dist1_prv(:, 1) = dist1
      dist2_prv(:, 1) = dist2

      ALLOCATE (local_particle_col(nkind), local_particle_row(nkind))
      ALLOCATE (nparticle_local_row(nkind), nparticle_local_col(nkind))
      nparticle_local_row = 0; nparticle_local_col = 0

      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (dist1_prv(iatom, 1) == mp_coor(1)) nparticle_local_row(ikind) = nparticle_local_row(ikind)+1
         IF (dist2_prv(iatom, 1) == mp_coor(2)) nparticle_local_col(ikind) = nparticle_local_col(ikind)+1
      END DO

      DO ikind = 1, nkind
         n = nparticle_local_row(ikind)
         ALLOCATE (local_particle_row(ikind)%array(n))

         n = nparticle_local_col(ikind)
         ALLOCATE (local_particle_col(ikind)%array(n))
      END DO

      nparticle_local_row = 0; nparticle_local_col = 0
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (dist1_prv(iatom, 1) == mp_coor(1)) THEN
            nparticle_local_row(ikind) = nparticle_local_row(ikind)+1
            local_particle_row(ikind)%array(nparticle_local_row(ikind)) = iatom
         END IF
         IF (dist2_prv(iatom, 1) == mp_coor(2)) THEN
            nparticle_local_col(ikind) = nparticle_local_col(ikind)+1
            local_particle_col(ikind)%array(nparticle_local_col(ikind)) = iatom
         END IF
      END DO

      CALL distribution_2d_create_prv(dist_2d, row_distribution_ptr=dist1_prv, &
                                      col_distribution_ptr=dist2_prv, local_rows_ptr=local_particle_row, &
                                      local_cols_ptr=local_particle_col, blacs_env=blacs_env)

      CALL cp_blacs_env_release(blacs_env)
   END SUBROUTINE

END MODULE
