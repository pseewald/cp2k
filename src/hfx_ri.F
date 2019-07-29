!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI-methods for HFX
! **************************************************************************************************

MODULE hfx_ri

   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_dot, dbcsr_frobenius_norm, dbcsr_get_info, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_scalar, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_type_real_8
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, &
        dbcsr_t_copy_tensor_to_matrix, dbcsr_t_create, dbcsr_t_destroy, &
        dbcsr_t_distribution_destroy, dbcsr_t_distribution_new, dbcsr_t_distribution_type, &
        dbcsr_t_get_info, dbcsr_t_mp_environ_pgrid, dbcsr_t_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE hfx_types,                       ONLY: hfx_cyclic_dist,&
                                              hfx_ri_type
   USE iterate_matrix,                  ONLY: invert_hotelling,&
                                              matrix_sqrt_newton_schulz
   USE kinds,                           ONLY: dp
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_3c_tensors,                   ONLY: build_2c_integrals,&
                                              build_2c_neighbor_lists,&
                                              build_3c_integrals,&
                                              build_3c_neighbor_lists,&
                                              neighbor_list_2c_destroy,&
                                              neighbor_list_3c_destroy
   USE qs_3c_tensors_types,             ONLY: neighbor_list_3c_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf(qs_env, ri_data)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      INTEGER                                            :: natom, nkind, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_AO, sizes_RI
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      REAL(KIND=dp)                                      :: frob_matrix, frob_matrix_base
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_t_type), DIMENSION(1)                   :: t_2c_int
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int_1, t_3c_int_2
      TYPE(dbcsr_type), DIMENSION(1) :: matrix_tmp1, t_2c_int_mat, t_2c_op_pot, t_2c_op_pot_sqrt, &
         t_2c_op_pot_sqrt2, t_2c_op_pot_sqrt_inv, t_2c_op_RI, t_2c_op_RI_inv
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c_pot, nl_2c_RI
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (col_bsize, row_bsize, logger, para_env, dist_2d, basis_set_AO, basis_set_RI, &
               nl_2c_pot, nl_2c_RI, particle_set, qs_kind_set)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL get_qs_env(qs_env, natom=natom, nkind=nkind, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      distribution_2d=dist_2d, para_env=para_env)

      ! determine whether RI metric is same operator as used in HFX
      !same_op = op_ri == op_pot

      !IF(same_op .AND. op_pot == do_potential_truncated) THEN
      !   same_op = ABS(ri_data%r_cutoff_ri - hfx_rid_data%r_cutoff_pot) < 1.0E-16_dp
      !ENDIF

      !IF(same_op .AND. op_pot == do_potential_short) THEN
      !   same_op = ABS(ri_data%omega_ri - ri_data%omega_pot) < 1.0E-16_dp
      !ENDIF

      !CPASSERT(same_op)
      ALLOCATE (sizes_RI(natom), sizes_AO(natom))
      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, "RI_HFX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_RI)
      CALL basis_set_list_setup(basis_set_AO, "ORB", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_AO)

      !dft_control=dft_control)

      CALL dbcsr_t_create(t_3c_int_1(1, 1), "(RI | AO AO)", ri_data%dist1, [1], [2, 3], &
                          dbcsr_type_real_8, sizes_RI, sizes_AO, sizes_AO)

      CALL build_3c_neighbor_lists(nl_3c, basis_set_RI, basis_set_AO, basis_set_AO, ri_data%dist1_3d, "HFX_3c_nl", &
                                   qs_env, op_ij=ri_data%op_ri, x_range=ri_data%r_cutoff_ri)
      CALL build_3c_integrals(t_3c_int_1, ri_data%filter_eps, qs_env, nl_3c, basis_set_RI, basis_set_AO, basis_set_AO, &
                              op_ij=ri_data%op_ri, t_c_filename=ri_data%t_c_filename, &
                              r_cutoff=ri_data%r_cutoff_ri, omega=ri_data%omega_ri)

      CALL neighbor_list_3c_destroy(nl_3c)

      CALL build_2c_neighbor_lists(nl_2c_pot, basis_set_RI, basis_set_RI, "HFX_2c_nl_pot", qs_env, sym_ij=.FALSE., &
                                   op_ij=ri_data%op_pot, x_range=ri_data%r_cutoff_pot, &
                                   dist_2d=dist_2d)
      CALL build_2c_neighbor_lists(nl_2c_RI, basis_set_RI, basis_set_RI, "HFX_2c_nl_RI", qs_env, sym_ij=.FALSE., &
                                   op_ij=ri_data%op_ri, x_range=ri_data%r_cutoff_ri, &
                                   dist_2d=dist_2d)

      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(SIZE(sizes_RI)))
      ALLOCATE (col_bsize(SIZE(sizes_RI)))
      row_bsize(:) = sizes_RI
      col_bsize(:) = sizes_RI
      CALL dbcsr_create(t_2c_op_pot(1), "(R|P)", dbcsr_dist, dbcsr_type_no_symmetry, row_bsize, col_bsize)
      CALL dbcsr_create(t_2c_op_RI(1), template=t_2c_op_pot(1))
      CALL dbcsr_distribution_release(dbcsr_dist)

      CALL build_2c_integrals(t_2c_op_pot, ri_data%filter_eps, qs_env, nl_2c_pot, basis_set_RI, basis_set_RI, &
                              op_ij=ri_data%op_pot, &
                              t_c_filename=ri_data%t_c_filename, &
                              r_cutoff=ri_data%r_cutoff_pot, omega=ri_data%omega_pot)

      CALL neighbor_list_2c_destroy(nl_2c_pot)

      CALL build_2c_integrals(t_2c_op_RI, ri_data%filter_eps, qs_env, nl_2c_RI, basis_set_RI, basis_set_RI, &
                              op_ij=ri_data%op_ri, t_c_filename=ri_data%t_c_filename, &
                              r_cutoff=ri_data%r_cutoff_ri, omega=ri_data%omega_ri)

      CALL neighbor_list_2c_destroy(nl_2c_RI)

      ! todo: if problems, try other methods for matrix functions
      CALL dbcsr_create(t_2c_op_RI_inv(1), template=t_2c_op_RI(1), matrix_type=dbcsr_type_no_symmetry)
      CALL invert_hotelling(t_2c_op_RI_inv(1), t_2c_op_RI(1), ri_data%filter_eps)
      IF (ri_data%check_2c_inv) THEN
         CALL dbcsr_create(matrix_tmp1(1), template=t_2c_op_RI(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "N", 1.0_dp, t_2c_op_RI_inv(1), t_2c_op_RI(1), &
                             0.0_dp, matrix_tmp1(1), filter_eps=ri_data%filter_eps)
         frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp1(1))
         CALL dbcsr_add_on_diag(matrix_tmp1(1), -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(matrix_tmp1(1))
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) "Error for (inv(K)*K-I)", frob_matrix/frob_matrix_base
         ENDIF
         CALL dbcsr_release(matrix_tmp1(1))
      ENDIF
      CALL dbcsr_release(t_2c_op_RI(1))

      CALL dbcsr_create(t_2c_op_pot_sqrt(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(t_2c_op_pot_sqrt_inv(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_no_symmetry)
      CALL matrix_sqrt_newton_schulz(t_2c_op_pot_sqrt(1), t_2c_op_pot_sqrt_inv(1), t_2c_op_pot(1), &
                                     ri_data%filter_eps, ri_data%t2c_sqrt_order, ri_data%eps_lanczos, &
                                     ri_data%max_iter_lanczos)
      CALL dbcsr_release(t_2c_op_pot_sqrt_inv(1))
      IF (ri_data%check_2c_inv) THEN
         CALL dbcsr_create(matrix_tmp1(1), template=t_2c_op_pot(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(t_2c_op_pot_sqrt2(1), template=t_2c_op_pot_sqrt(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(t_2c_op_pot_sqrt2(1), t_2c_op_pot_sqrt(1))
         CALL dbcsr_multiply("N", "N", 1.0_dp, t_2c_op_pot_sqrt(1), t_2c_op_pot_sqrt2(1), &
                             0.0_dp, matrix_tmp1(1), filter_eps=ri_data%filter_eps)
         CALL dbcsr_add(matrix_tmp1(1), t_2c_op_pot(1), 1.0_dp, -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(matrix_tmp1(1))
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) "Error for (K^(1/2)*K^(1/2) - K)", frob_matrix
         ENDIF
         CALL dbcsr_release(matrix_tmp1(1))
         CALL dbcsr_release(t_2c_op_pot_sqrt2(1))
      ENDIF

      CALL dbcsr_create(t_2c_int_mat(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_multiply("N", "N", 1.0_dp, t_2c_op_RI_inv(1), t_2c_op_pot_sqrt(1), &
                          0.0_dp, t_2c_int_mat(1), filter_eps=ri_data%filter_eps)
      CALL dbcsr_release(t_2c_op_pot(1))
      CALL dbcsr_release(t_2c_op_RI_inv(1))
      CALL dbcsr_release(t_2c_op_pot_sqrt(1))

      CALL dbcsr_t_create(t_2c_int_mat(1), t_2c_int(1), name="(RI|RI)")
      CALL dbcsr_t_copy_matrix_to_tensor(t_2c_int_mat(1), t_2c_int(1))
      CALL dbcsr_release(t_2c_int_mat(1))

      CALL dbcsr_t_create(t_3c_int_1(1, 1), t_3c_int_2(1, 1), name="(RI | AO AO)")
      CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), t_2c_int(1), t_3c_int_1(1, 1), &
                            dbcsr_scalar(0.0_dp), t_3c_int_2(1, 1), &
                            contract_1=[2], notcontract_1=[1], &
                            contract_2=[1], notcontract_2=[2, 3], &
                            map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, move_data=.TRUE., &
                            unit_nr=unit_nr)
      CALL dbcsr_t_destroy(t_2c_int(1))
      CALL dbcsr_t_destroy(t_3c_int_1(1, 1))

      CALL dbcsr_t_copy(t_3c_int_2(1, 1), ri_data%t_3c_int(1, 1), order=[2, 1, 3], move_data=.TRUE.)
      CALL dbcsr_t_destroy(t_3c_int_2(1, 1))

      DEALLOCATE (basis_set_RI, basis_set_AO)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param ks_matrix_aux_fit_hfx ...
!> \param mo_coeff ...
!> \param rho_ao ...
!> \param do_admm ...
!> \param geometry_did_change ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks(qs_env, ri_data, ks_matrix, ehfx, ks_matrix_aux_fit_hfx, mo_coeff, rho_ao, &
                               do_admm, geometry_did_change, nspins)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: ks_matrix_aux_fit_hfx
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: rho_ao
      LOGICAL, INTENT(IN)                                :: do_admm, geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins

      INTEGER                                            :: ispin, n_mo_blocks
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsize_2, bsize_3, dist_MO
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: mo_bsize
      REAL(KIND=dp)                                      :: etmp
      TYPE(dbcsr_t_distribution_type)                    :: dist1, dist2
      TYPE(dbcsr_t_type)                                 :: ks_t, mo_coeff_t
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int_ctr, t_3c_int_ctr_opt, &
                                                            t_3c_int_ctr_opt_copy

      NULLIFY (mo_bsize)

      CPASSERT(SIZE(ks_matrix, 2) == 1)

      ! todo: mo localization

      ! assuming TREAT_LSD_IN_CORE is .FALSE.
      IF (geometry_did_change) CALL hfx_ri_pre_scf(qs_env, ri_data)

      !PRINT *, "checksum t3c tensor", dbcsr_t_checksum(ri_data%t_3c_int(1, 1))

      CALL dbcsr_t_mp_environ_pgrid(ri_data%pgrid, pdims, pcoord)

      CALL dbcsr_get_info(mo_coeff(1), nblkcols_total=n_mo_blocks, col_blk_size=mo_bsize)
      ALLOCATE (dist_MO(n_mo_blocks))
      CALL hfx_cyclic_dist(n_mo_blocks, pdims(1), mo_bsize, dist_MO)
      CALL dbcsr_t_distribution_new(dist1, ri_data%pgrid, [1], [2, 3], &
                                    dist_MO, ri_data%dist2_RI, ri_data%dist2_ao_2, dimsplit=2)
      CALL dbcsr_t_get_info(ri_data%t_3c_int(1, 1), blk_size_2=bsize_2, blk_size_3=bsize_3)
      CALL dbcsr_t_create(t_3c_int_ctr(1, 1), "(MO | RI AO)", dist1, [1], [2, 3], dbcsr_type_real_8, &
                          mo_bsize, bsize_2, bsize_3)
      CALL dbcsr_t_distribution_destroy(dist1)
      CALL dbcsr_t_distribution_new(dist2, ri_data%pgrid, [1, 2], [3], &
                                    dist_MO, ri_data%dist2_RI, ri_data%dist2_ao_2, dimsplit=1)
      CALL dbcsr_t_create(t_3c_int_ctr_opt(1, 1), "(MO RI | AO)", dist2, [1, 2], [3], dbcsr_type_real_8, &
                          mo_bsize, bsize_2, bsize_3)

      CALL dbcsr_t_create(t_3c_int_ctr_opt(1, 1), t_3c_int_ctr_opt_copy(1, 1))
      CALL dbcsr_t_distribution_destroy(dist2)
      CALL dbcsr_t_create(ks_matrix(1, 1)%matrix, ks_t)

      CALL dbcsr_t_create(mo_coeff(1), mo_coeff_t)

      DO ispin = 1, nspins
         ! todo: dbcsr_t_copy interface
         CALL dbcsr_t_copy_matrix_to_tensor(mo_coeff(ispin), mo_coeff_t)
         !PRINT *, "checksum mo coeff matrix", dbcsr_checksum(mo_coeff(ispin))
         !PRINT *, "checksum mo coeff tensor", dbcsr_t_checksum(mo_coeff_t)
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_t, ri_data%t_3c_int(1, 1), dbcsr_scalar(0.0_dp), t_3c_int_ctr(1, 1), &
                               contract_1=[1], notcontract_1=[2], &
                               contract_2=[1], notcontract_2=[2, 3], &
                               map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps)

         !PRINT *, "checksum contraction 1 out", dbcsr_t_checksum(t_3c_int_ctr(1,1))
         CALL dbcsr_t_copy(t_3c_int_ctr(1, 1), t_3c_int_ctr_opt(1, 1), move_data=.TRUE.)
         !PRINT *, "checksum opt", dbcsr_t_checksum(t_3c_int_ctr_opt(1,1))
         CALL dbcsr_t_copy(t_3c_int_ctr_opt(1, 1), t_3c_int_ctr_opt_copy(1, 1))
         !PRINT *, "checksum opt copy", dbcsr_t_checksum(t_3c_int_ctr_opt_copy(1,1))
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), t_3c_int_ctr_opt(1, 1), t_3c_int_ctr_opt_copy(1, 1), &
                               dbcsr_scalar(0.0_dp), ks_t, &
                               contract_1=[1, 2], notcontract_1=[3], &
                               contract_2=[1, 2], notcontract_2=[3], &
                               map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps)
         !PRINT *, "checksum ks tensor", ispin, dbcsr_t_checksum(ks_t)
         CALL dbcsr_t_copy_tensor_to_matrix(ks_t, ks_matrix(ispin, 1)%matrix, summation=.TRUE.)

      ENDDO

      CALL dbcsr_t_destroy(t_3c_int_ctr(1, 1))
      CALL dbcsr_t_destroy(t_3c_int_ctr_opt(1, 1))
      CALL dbcsr_t_destroy(t_3c_int_ctr_opt_copy(1, 1))
      CALL dbcsr_t_destroy(ks_t)
      CALL dbcsr_t_destroy(mo_coeff_t)

      !! ** Calculate the exchange energy, copied from integrate_four_center
      ehfx = 0.0_dp
      DO ispin = 1, nspins
         CALL dbcsr_dot(ks_matrix(ispin, 1)%matrix, rho_ao(ispin, 1)%matrix, &
                        etmp)
         ehfx = ehfx+0.5_dp*etmp
         IF (do_admm) THEN
            !for ADMMS, we need the exchange matrix k(d) for both spins
            CALL dbcsr_copy(ks_matrix_aux_fit_hfx(ispin)%matrix, ks_matrix(ispin, 1)%matrix, &
                            name="HF exch. part of matrix_ks_aux_fit for ADMMS")
         END IF
      ENDDO

      PRINT *, "RI HFX Energy", ehfx

   END SUBROUTINE

END MODULE
