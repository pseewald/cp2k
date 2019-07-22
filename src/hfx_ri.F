!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI-methods for HFX
! **************************************************************************************************

MODULE hfx_ri

   USE auto_basis,                      ONLY: create_ri_aux_basis_set
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type,&
                                              init_orb_basis_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_create, dbcsr_frobenius_norm, dbcsr_get_info, dbcsr_multiply, dbcsr_p_type, &
        dbcsr_release, dbcsr_scalar, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_real_8
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, &
        dbcsr_t_copy_tensor_to_matrix, dbcsr_t_create, dbcsr_t_distribution_destroy, &
        dbcsr_t_distribution_new, dbcsr_t_distribution_type, dbcsr_t_get_info, &
        dbcsr_t_mp_environ_pgrid, dbcsr_t_pgrid_create, dbcsr_t_pgrid_type, dbcsr_t_type
   USE hfx_types,                       ONLY: hfx_ri_type
   USE iterate_matrix,                  ONLY: invert_hotelling,&
                                              matrix_sqrt_newton_schulz
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: erfc_cutoff
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_3c_tensors,                   ONLY: build_2c_integrals,&
                                              build_2c_neighbor_lists,&
                                              build_3c_integrals,&
                                              build_3c_neighbor_lists,&
                                              distribution_3d_create
   USE qs_3c_tensors_types,             ONLY: distribution_3d_type,&
                                              neighbor_list_3c_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param hfx_ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_init(hfx_ri_data, qs_env)
      TYPE(hfx_ri_type), INTENT(OUT)                     :: hfx_ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: ikind, mp_comm_3d, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_AO_1, dist_AO_2, dist_RI, sizes_AO, &
                                                            sizes_RI
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: ri_aux_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL get_qs_env(qs_env, &
                      nkind=nkind, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      para_env=para_env)

      ! initialise RI basis
      DO ikind = 1, nkind
         NULLIFY (ri_aux_basis_set)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=ri_aux_basis_set, &
                          basis_type="RI_HFX")
         IF (.NOT. (ASSOCIATED(ri_aux_basis_set))) THEN
            ! Generate a default basis
            CALL create_ri_aux_basis_set(ri_aux_basis_set, qs_kind, dft_control%auto_basis_ri_aux)
            CALL init_orb_basis_set(ri_aux_basis_set)
            CALL init_interaction_radii_orb_basis(ri_aux_basis_set, dft_control%qs_control%eps_pgf_orb)
            CALL add_basis_set_to_container(qs_kind%basis_sets, ri_aux_basis_set, "RI_HFX")
         END IF
      END DO

      ! todo: read from ri section
      CALL erfc_cutoff(hfx_ri_data%eps_range, hfx_ri_data%omega_ri, hfx_ri_data%r_cutoff_ri)
      CALL erfc_cutoff(hfx_ri_data%eps_range, hfx_ri_data%omega_pot, hfx_ri_data%r_cutoff_pot)

      pdims = 0
      CALL mp_cart_create(para_env%group, 3, pdims, pcoord, mp_comm_3d)

      ALLOCATE (sizes_RI(natom), sizes_AO(natom))
      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, "RI_HFX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_RI)
      CALL basis_set_list_setup(basis_set_AO, "ORB", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_AO)

      ALLOCATE (dist_AO_1(natom))
      ALLOCATE (dist_AO_2(natom))
      ALLOCATE (dist_RI(natom))
      CALL cyclic_dist(natom, pdims(1), sizes_AO, dist_AO_1)
      CALL cyclic_dist(natom, pdims(2), sizes_AO, dist_AO_2)
      CALL cyclic_dist(natom, pdims(3), sizes_RI, dist_RI)
      CALL distribution_3d_create(qs_env, dist_3d, dist_ao_1, dist_ao_2, dist_RI, mp_comm_3d, own_comm=.TRUE.)

      CALL dbcsr_t_pgrid_create(para_env%group, pdims, pgrid, map1_2d=[1, 2], map2_2d=[3])
      hfx_ri_data%pgrid = pgrid

      hfx_ri_data%dist1_ao_1(:) = dist_ao_1
      hfx_ri_data%dist1_ao_2(:) = dist_ao_2
      hfx_ri_data%dist1_ri(:) = dist_RI
      hfx_ri_data%dist1_3d = dist_3d

      CALL cyclic_dist(natom, pdims(1), sizes_AO, dist_AO_1)
      CALL cyclic_dist(natom, pdims(3), sizes_AO, dist_AO_2)
      CALL cyclic_dist(natom, pdims(2), sizes_RI, dist_RI)

      hfx_ri_data%dist2_ao_1(:) = dist_ao_1
      hfx_ri_data%dist2_ao_2(:) = dist_ao_2
      hfx_ri_data%dist2_ri(:) = dist_RI

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param hfx_ri_data ...
!> \param t_3c_int ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf(qs_env, hfx_ri_data, t_3c_int)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(IN)                      :: hfx_ri_data
      TYPE(dbcsr_t_type), DIMENSION(1, 1), INTENT(OUT)   :: t_3c_int

      INTEGER                                            :: natom, nkind, op_pot, op_ri, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_AO, sizes_RI
      REAL(KIND=dp)                                      :: frob_matrix, frob_matrix_base
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_distribution_type)                    :: dist_1, dist_2
      TYPE(dbcsr_t_type), DIMENSION(1)                   :: t_2c_int
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int_1, t_3c_int_2
      TYPE(dbcsr_type), DIMENSION(1) :: matrix_tmp1, t_2c_int_mat, t_2c_op_pot, t_2c_op_pot_sqrt, &
         t_2c_op_pot_sqrt2, t_2c_op_pot_sqrt_inv, t_2c_op_RI, t_2c_op_RI_inv
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_AO, basis_set_RI
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c_pot, nl_2c_RI
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      op_pot = hfx_ri_data%op_pot
      op_ri = hfx_ri_data%op_ri
      ! determine whether RI metric is same operator as used in HFX
      !same_op = op_ri == op_pot

      !IF(same_op .AND. op_pot == do_potential_truncated) THEN
      !   same_op = ABS(hfx_ri_data%r_cutoff_ri - hfx_rid_data%r_cutoff_pot) < 1.0E-16_dp
      !ENDIF

      !IF(same_op .AND. op_pot == do_potential_short) THEN
      !   same_op = ABS(hfx_ri_data%omega_ri - hfx_ri_data%omega_pot) < 1.0E-16_dp
      !ENDIF

      !CPASSERT(same_op)
      ALLOCATE (sizes_RI(natom), sizes_AO(natom))
      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, "RI_HFX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_RI)
      CALL basis_set_list_setup(basis_set_AO, "ORB", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_AO)

      CALL get_qs_env(qs_env, natom=natom, nkind=nkind, qs_kind_set=qs_kind_set, para_env=para_env)
      !dft_control=dft_control)

      CALL dbcsr_t_distribution_new(dist_1, hfx_ri_data%pgrid, [1, 2], [3], &
                                    hfx_ri_data%dist1_ao_1, hfx_ri_data%dist1_ao_2, hfx_ri_data%dist1_RI, dimsplit=1)
      CALL dbcsr_t_create(t_3c_int_1(1, 1), "(munu|P)", dist_1, [1, 2], [3], dbcsr_type_real_8, sizes_AO, sizes_AO, sizes_RI)

      CALL build_3c_neighbor_lists(nl_3c, basis_set_AO, basis_set_AO, basis_set_RI, hfx_ri_data%dist1_3d, "HFX_3c_nl", &
                                   qs_env, op_jk=op_ri, x_range=hfx_ri_data%r_cutoff_ri, own_dist=.TRUE.)
      CALL build_3c_integrals(t_3c_int_1, hfx_ri_data%filter_eps, qs_env, nl_3c, basis_set_AO, basis_set_AO, &
                              basis_set_RI, op_jk=op_ri, t_c_filename=hfx_ri_data%t_c_filename, &
                              r_cutoff=hfx_ri_data%r_cutoff_ri, omega=hfx_ri_data%omega_ri)

      CALL dbcsr_t_distribution_destroy(dist_1)

      CALL build_2c_neighbor_lists(nl_2c_pot, basis_set_RI, basis_set_RI, "HFX_2c_nl_pot", qs_env, sym_ij=.FALSE., &
                                   op_ij=op_pot, x_range=hfx_ri_data%r_cutoff_pot)
      CALL build_2c_neighbor_lists(nl_2c_RI, basis_set_RI, basis_set_RI, "HFX_2c_nl_RI", qs_env, sym_ij=.FALSE., &
                                   op_ij=op_ri, x_range=hfx_ri_data%r_cutoff_pot)

      CALL build_2c_integrals(t_2c_op_pot, hfx_ri_data%filter_eps, qs_env, nl_2c_pot, basis_set_RI, basis_set_RI, &
                              op_ij=hfx_ri_data%op_pot, t_c_filename=hfx_ri_data%t_c_filename, &
                              r_cutoff=hfx_ri_data%r_cutoff_pot, omega=hfx_ri_data%omega_pot)

      CALL build_2c_integrals(t_2c_op_RI, hfx_ri_data%filter_eps, qs_env, nl_2c_RI, basis_set_RI, basis_set_RI, &
                              op_ij=hfx_ri_data%op_ri, t_c_filename=hfx_ri_data%t_c_filename, &
                              r_cutoff=hfx_ri_data%r_cutoff_ri, omega=hfx_ri_data%omega_ri)

      ! todo: if problems, try other methods for matrix functions
      CALL dbcsr_create(t_2c_op_RI_inv(1), template=t_2c_op_RI(1), matrix_type=dbcsr_type_no_symmetry)
      CALL invert_hotelling(t_2c_op_RI_inv(1), t_2c_op_RI(1), hfx_ri_data%filter_eps)
      IF (hfx_ri_data%check_2c_inv) THEN
         CALL dbcsr_create(matrix_tmp1(1), template=t_2c_op_RI(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "N", 1.0_dp, t_2c_op_RI_inv(1), t_2c_op_RI(1), &
                             0.0_dp, matrix_tmp1(1), filter_eps=hfx_ri_data%filter_eps)
         frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp1(1))
         CALL dbcsr_add_on_diag(matrix_tmp1(1), -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(matrix_tmp1(1))
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) "Error for (inv(K)*K-I)", frob_matrix/frob_matrix_base
         ENDIF
         CALL dbcsr_release(matrix_tmp1(1))
      ENDIF

      CALL dbcsr_create(t_2c_op_pot_sqrt(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(t_2c_op_pot_sqrt_inv(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_no_symmetry)
      CALL matrix_sqrt_newton_schulz(t_2c_op_pot_sqrt(1), t_2c_op_pot_sqrt_inv(1), t_2c_op_pot(1), &
                                     hfx_ri_data%filter_eps, hfx_ri_data%s_sqrt_order, hfx_ri_data%eps_lanczos, &
                                     hfx_ri_data%max_iter_lanczos)
      CALL dbcsr_release(t_2c_op_pot_sqrt_inv(1))
      IF (hfx_ri_data%check_2c_inv) THEN
         CALL dbcsr_create(matrix_tmp1(1), template=t_2c_op_pot(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(t_2c_op_pot_sqrt2(1), template=t_2c_op_pot_sqrt(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(t_2c_op_pot_sqrt2(1), t_2c_op_pot_sqrt(1))
         CALL dbcsr_multiply("N", "N", 1.0_dp, t_2c_op_pot_sqrt(1), t_2c_op_pot_sqrt2(1), &
                             0.0_dp, matrix_tmp1(1), filter_eps=hfx_ri_data%filter_eps)
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
                          0.0_dp, t_2c_int_mat(1), filter_eps=hfx_ri_data%filter_eps)
      CALL dbcsr_release(t_2c_op_pot(1))
      CALL dbcsr_release(t_2c_op_RI_inv(1))
      CALL dbcsr_release(t_2c_op_pot_sqrt(1))

      CALL dbcsr_t_create(t_2c_int_mat(1), t_2c_int(1), name="(P|R)")
      CALL dbcsr_t_copy_matrix_to_tensor(t_2c_int_mat(1), t_2c_int(1))
      CALL dbcsr_release(t_2c_int_mat(1))

      CALL dbcsr_t_create(t_3c_int_2(1, 1), t_3c_int_1(1, 1), name="(munu|R)")
      CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), t_3c_int_1(1, 1), t_2c_int(1), dbcsr_scalar(0.0_dp), t_3c_int_2(1, 1), &
                            contract_1=[3], notcontract_1=[1, 2], &
                            contract_2=[1], notcontract_2=[2], &
                            map_1=[1, 2], map_2=[3], filter_eps=hfx_ri_data%filter_eps, move_data=.TRUE., &
                            unit_nr=unit_nr)
      CALL dbcsr_t_destroy(t_2c_int(1))
      CALL dbcsr_t_destroy(t_3c_int_1(1, 1))

      CALL dbcsr_t_distribution_new(dist_2, hfx_ri_data%pgrid, [1, 2], [3], &
                                    hfx_ri_data%dist2_ao_1, hfx_ri_data%dist2_ri, hfx_ri_data%dist2_ao_2, dimsplit=1)
      CALL dbcsr_t_create(t_3c_int(1, 1), "(muR|nu)", dist_2, [1, 2], [3], dbcsr_type_real_8, sizes_AO, sizes_RI, sizes_AO)
      CALL dbcsr_t_copy(t_3c_int_2(1, 1), t_3c_int(1, 1), order=[1, 3, 2], move_data=.TRUE.)
      CALL dbcsr_t_destroy(t_3c_int_2(1, 1))
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param hfx_ri_data ...
!> \param ks_matrix ...
!> \param mo_coeff ...
!> \param t_3c_int ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks(hfx_ri_data, ks_matrix, mo_coeff, t_3c_int)
      TYPE(hfx_ri_type), INTENT(IN)                      :: hfx_ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int

      INTEGER                                            :: ispin, n_mo_blocks, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsize_1, bsize_2, dist_MO
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: mo_bsize
      TYPE(dbcsr_t_distribution_type)                    :: dist1, dist2
      TYPE(dbcsr_t_type)                                 :: ks_t, mo_coeff_t
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int_ctr, t_3c_int_ctr_opt, &
                                                            t_3c_int_ctr_opt_copy

      nspins = SIZE(ks_matrix, 1)
      CPASSERT(SIZE(ks_matrix, 2) == 1)

      ! todo: mo localization

      CALL dbcsr_t_mp_environ_pgrid(hfx_ri_data%pgrid, pdims, pcoord)

      CALL dbcsr_get_info(mo_coeff(1), nblkcols_total=n_mo_blocks, col_blk_size=mo_bsize)
      ALLOCATE (dist_MO(n_mo_blocks))
      CALL cyclic_dist(n_mo_blocks, pdims(3), mo_bsize, dist_MO)
      CALL dbcsr_t_distribution_new(dist1, hfx_ri_data%pgrid, [1, 2], [3], &
                                    hfx_ri_data%dist2_ao_1, hfx_ri_data%dist2_RI, dist_MO, dimsplit=1)
      CALL dbcsr_t_get_info(t_3c_int(1, 1), blk_size_1=bsize_1, blk_size_2=bsize_2)
      CALL dbcsr_t_create(t_3c_int_ctr(1, 1), "(muR|i)", dist1, [1, 2], [3], dbcsr_type_real_8, &
                          bsize_1, bsize_2, mo_bsize)
      CALL dbcsr_t_distribution_destroy(dist1)
      CALL dbcsr_t_distribution_new(dist2, hfx_ri_data%pgrid, [1], [2, 3], &
                                    hfx_ri_data%dist2_ao_1, hfx_ri_data%dist2_RI, dist_MO, dimsplit=2)
      CALL dbcsr_t_create(t_3c_int_ctr_opt(1, 1), "(mu|Ri)", dist2, [1], [2, 3], dbcsr_type_real_8, &
                          bsize_1, bsize_2, mo_bsize)

      CALL dbcsr_t_create(t_3c_int_ctr_opt(1, 1), t_3c_int_ctr_opt_copy(1, 1))
      CALL dbcsr_t_distribution_destroy(dist2)
      CALL dbcsr_t_create(ks_matrix(1, 1)%matrix, ks_t)

      CALL dbcsr_t_create(mo_coeff(1), mo_coeff_t)

      DO ispin = 1, nspins
         ! todo: dbcsr_t_copy interface
         CALL dbcsr_t_copy_matrix_to_tensor(mo_coeff(ispin), mo_coeff_t)
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), t_3c_int(1, 1), mo_coeff_t, dbcsr_scalar(0.0_dp), t_3c_int_ctr(1, 1), &
                               contract_1=[3], notcontract_1=[1, 2], &
                               contract_2=[1], notcontract_2=[2], &
                               map_1=[1, 2], map_2=[3], filter_eps=hfx_ri_data%filter_eps)

         CALL dbcsr_t_copy(t_3c_int_ctr(1, 1), t_3c_int_ctr_opt(1, 1), move_data=.TRUE.)
         CALL dbcsr_t_copy(t_3c_int_ctr_opt(1, 1), t_3c_int_ctr_opt_copy(1, 1))
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), t_3c_int_ctr_opt(1, 1), t_3c_int_ctr_opt_copy(1, 1), &
                               dbcsr_scalar(0.0_dp), ks_t, &
                               contract_1=[2, 3], notcontract_1=[1], &
                               contract_2=[2, 3], notcontract_2=[1], &
                               map_1=[1], map_2=[2], filter_eps=hfx_ri_data%filter_eps)
         CALL dbcsr_t_copy_tensor_to_matrix(ks_t, ks_matrix(ispin, 1)%matrix, summation=.TRUE.)

      ENDDO

      CALL dbcsr_t_destroy(t_3c_int_ctr(1, 1))
      CALL dbcsr_t_destroy(t_3c_int_ctr_opt(1, 1))
      CALL dbcsr_t_destroy(t_3c_int_ctr_opt_copy(1, 1))
      CALL dbcsr_t_destroy(ks_t)
      CALL dbcsr_t_destroy(mo_coeff_t)

   END SUBROUTINE

! **************************************************************************************************
!> \brief cyclic distribution of weighted elements todo: move this to somewhere else
!> \param nel ...
!> \param nbin ...
!> \param weights ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE cyclic_dist(nel, nbin, weights, dist)
      INTEGER, INTENT(IN)                                :: nel, nbin
      INTEGER, DIMENSION(nel), INTENT(IN)                :: weights
      INTEGER, DIMENSION(nel), INTENT(OUT)               :: dist

      INTEGER                                            :: ibin, iel, niter
      INTEGER, DIMENSION(nbin)                           :: occup

      occup(:) = 0
      ibin = 0
      DO iel = 1, nel
         niter = 0
         ibin = MOD(ibin+1, nbin)
         DO WHILE (occup(ibin+1)+weights(iel) .GE. MAXVAL(occup))
            IF (MINLOC(occup, DIM=1) == ibin+1) EXIT
            ibin = MOD(ibin+1, nbin)
            niter = niter+1
         ENDDO
         dist(iel) = ibin
         occup(ibin+1) = occup(ibin+1)+weights(iel)
      ENDDO

   END SUBROUTINE cyclic_dist

   !SUBROUTINE integrate_RI(qs_env, hfx_ri_data, ks_matrix, ehfx, hfx_ri_section, para_env, geometry_did_change, ispin)
   !   TYPE(qs_environment_type), POINTER                 :: qs_env
   !   TYPE(dbcsr_p_type)                                 :: ks_matrix
   !   REAL(KIND=dp), INTENT(OUT)                         :: ehfx
   !   TYPE(section_vals_type), POINTER                   :: hfx_ri_section
   !   TYPE(cp_para_env_type), POINTER                    :: para_env
   !   LOGICAL                                            :: geometry_did_change
   !   INTEGER, INTENT(IN)                                :: ispin
   !END SUBROUTINE

END MODULE
