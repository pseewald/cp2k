!--------------------------------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations                              ! 
!   Copyright (C) 2000 - 2018  CP2K developers group                                               ! 
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief This module deals with the RI GAPW part of the XAS TDDFPT calculations. For all excited
!>        atoms it computes and stores the soft basis based on the RI_XAS basis
!>        It then computes the matrix element <P|f(r)|Q> for P,Q in the RI basis and f(r) an
!>        analytical function of r. The GAPW procedure is simplified because there are no cross 
!>        atomic terms, and one does not need to sum over all atoms.
!>        Keeps all parameters (grid, pw, etc) from the gapw input
! **************************************************************************************************
MODULE xas_tdp_gapw
   USE ai_contraction_sphi,             ONLY: ab_contract
   USE ao_util,                         ONLY: exp_radius, transform_s2c
   USE atomic_kind_types,               ONLY: atomic_kind_type, get_atomic_kind
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: gto_basis_set_type, gto_basis_set_p_type, &
                                              allocate_gto_basis_set, get_gto_basis_set, &
                                              set_gto_basis_set, write_orb_basis_set
   USE cell_types,                      ONLY: cell_type, pbc
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type, cp_1d_r_p_type
   USE cp_control_types,                ONLY: qs_control_type, dft_control_type
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type, dbcsr_distribution_type,&
                                              dbcsr_create, dbcsr_print, dbcsr_put_block, &
                                              dbcsr_finalize, dbcsr_release, dbcsr_multiply, &
                                              dbcsr_type, dbcsr_release_p, &
                                              dbcsr_get_block_p
   USE gaussian_gridlevels,             ONLY: gridlevel_info_type, gaussian_gridlevel
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids, get_number_of_lebedev_grid,&
                                              init_lebedev_grids, lebedev_grid
   USE mathconstants,                   ONLY: dfac, pi
   USE mathlib,                         ONLY: get_pseudo_inverse_svd
   USE memory_utilities,                ONLY: reallocate
   USE qs_modify_pab_block,             ONLY: FUNC_AB
   USE orbital_pointers,                ONLY: nsoset, indso, ncoset, nco, nso, indco
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE paw_proj_set_types,              ONLY: paw_proj_set_type
   USE pw_env_methods,                  ONLY: pw_env_rebuild, pw_env_create
   USE pw_env_types,                    ONLY: pw_env_type, pw_env_get, pw_env_release
   USE pw_methods,                      ONLY: pw_zero, pw_transfer, pw_axpy, pw_integrate_function
   USE pw_pool_types,                   ONLY: pw_pool_p_type, pw_pools_create_pws, &
                                              pw_pools_give_back_pws, pw_pool_create_pw, &
                                              pw_pool_give_back_pw, pw_pool_type
   USE pw_types,                        ONLY: pw_p_type, pw_create, REALSPACE, REALDATA3D, &
                                              COMPLEXDATA1D, RECIPROCALSPACE
   USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
   USE qs_environment_types,            ONLY: qs_environment_type, get_qs_env
   USE qs_grid_atom,                    ONLY: grid_atom_type, allocate_grid_atom, &
                                              deallocate_grid_atom, create_grid_atom
   USE qs_harmonics_atom,               ONLY: create_harmonics_atom, get_maxl_CG, &
                                              harmonics_atom_type, allocate_harmonics_atom, &
                                              deallocate_harmonics_atom
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: qs_kind_type, get_qs_kind, get_qs_kind_set
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type, &
                                              deallocate_neighbor_list_set
   USE qs_o3c_methods,                  ONLY: calculate_o3c_integrals, contract12_o3c
   USE qs_o3c_types,                    ONLY: init_o3c_container, o3c_container_type, &
                                              release_o3c_container, o3c_iterate, &
                                              o3c_iterator_create, o3c_iterator_release, &
                                              o3c_iterator_type, get_o3c_iterator_info
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_type, qs_rho_get
   USE qs_rho_atom_types,               ONLY: rho_atom_type, allocate_rho_atom_set, &
                                              deallocate_rho_atom_set
   USE qs_vxc_atom,                     ONLY: gaVxcgb_noGC
   USE realspace_grid_types,            ONLY: realspace_grid_p_type, realspace_grid_type, &
                                              rs_grid_retain, rs_grid_zero, rs_pw_transfer, &
                                              rs_grid_release, rs2pw
   USE spherical_harmonics,             ONLY: clebsch_gordon, clebsch_gordon_deallocate, &
                                              clebsch_gordon_init
   USE soft_basis_set,                  ONLY: create_soft_basis
   USE task_list_types,                 ONLY: task_list_type, allocate_task_list,&
                                              deallocate_task_list
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE xc_functionals_utilities,        ONLY: set_util, calc_rs_pw

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "xas_tdp_gapw"

!  Some helper types
   TYPE grid_atom_p_type
      TYPE(grid_atom_type), POINTER                   :: grid_atom
   END TYPE grid_atom_p_type

   TYPE harmonics_atom_p_type
      TYPE(harmonics_atom_type), POINTER              :: harmonics_atom
   END TYPE harmonics_atom_p_type

! **************************************************************************************************
!> \brief a environment type that contains all the info needed for the RI GAPW calculations 
!> \param rho_atom the type containing the grid integral info for all excited atoms 
!> \param grid_atom_set the set of atomic grid for each kind
!> \param task_list the task list associated with the RI_XAS_SOFT basis for the pw integration
!> \param blk_size the block sizes of the dbcsr_matrix that will hold the integrals
!> \param sphi_s the spherical contraction coeffs fro the soft basis, in the same fornat as the 
!>               hard basis. TODO: use paw_porj_set if it turns out we need more
!> \param sri the RI overlap neighbor list. Atoms are only neighbor with themselves
!> \param coeff_h the coefficients to express the density in the hard RI basis for each atom
!> \param coeff_s the coefficients to express the density in the soft RI basis for each atom
!> \param pw_env a pw_env, based on the qs one, to integrate the smooth density
!> \note the rho_atom_set has the size of the number of atoms, even thogh some are not excited
!>       same for grid_atom_set, it has the same size as the number of kind, but some are empty
!>       same for harmonics_atom_set
! **************************************************************************************************
   TYPE ri_gapw_env_type
      INTEGER, DIMENSION(:), POINTER                  :: excited_kinds
      TYPE(rho_atom_type), DIMENSION(:) , POINTER     :: rho_atom_set
      TYPE(grid_atom_p_type), DIMENSION(:), POINTER   :: grid_atom_set
      TYPE(harmonics_atom_p_type), DIMENSION(:), &
                                            POINTER   :: harmonics_atom_set
      TYPE(task_list_type), POINTER                   :: task_list
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER     :: sphi_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
                                             POINTER  :: sri
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER     :: coeff_h
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER     :: coeff_s
      TYPE(pw_env_type), POINTER                      :: pw_env
   END TYPE ri_gapw_env_type

   PUBLIC :: ri_gapw_env_type

   PUBLIC :: init_ri_gapw_env, release_ri_gapw_env, integrate_ri_gapw_smooth, &
             integrate_ri_gapw_atoms

CONTAINS

! **************************************************************************************************
!> \brief Initializes a ri_gapw_env type given the qs_env
!> \param ri_gapw_env the ri_gapw_env to initialize
!> \param excited_kinds the indices of the excited kinds
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE init_ri_gapw_env(ri_gapw_env, excited_kinds, qs_env)

      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      INTEGER, DIMENSION(:), POINTER                  :: excited_kinds

      CHARACTER(len=*), PARAMETER :: routineN = 'init_ri_gapw_env', routineP = moduleN//":"//routineN

      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER               :: soft_basis, ri_basis
      REAL(dp)                                        :: rc, radius, eps_fit
      LOGICAL                                         :: paw_atom
      INTEGER                                         :: ikind, nkind, natom, nex_kinds, &
                                                         handle, nsgf, lmax_sphere, iset, &
                                                         nset, n, ipgf, ishell, lshell, icomin, &
                                                         icomax, isgfmin, isgfmax, isgf
      INTEGER, DIMENSION(:), POINTER                  :: lmax, npgf, nshell
      INTEGER, DIMENSION(:,:), POINTER                :: l, first_sgf, last_sgf
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_a 
      TYPE(qs_ks_env_type), POINTER                   :: ks_env
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      REAL(dp), DIMENSION(:,:), POINTER               :: sphi_h, sphi_s
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: radius_backup
      TYPE(qs_control_type), POINTER                  :: qs_control
      TYPE(dft_control_type), POINTER                 :: dft_control

      NULLIFY(qs_kind_set, soft_basis, ri_basis, basis_set_a, ks_env, particle_set, qs_control)
      NULLIFY(sphi_s, sphi_h, lmax, npgf, nshell, l, first_sgf, last_sgf, dft_control)

      CALL timeset(routineN, handle)

!  Initializing the type
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, ks_env=ks_env, natom=natom, &
                      particle_set=particle_set, dft_control=dft_control)

      qs_control => dft_control%qs_control
      nkind = SIZE(qs_kind_set)   
      nex_kinds = SIZE(excited_kinds)
      eps_fit = qs_control%gapw_control%eps_fit

      ALLOCATE(ri_gapw_env)
      CALL allocate_rho_atom_set(ri_gapw_env%rho_atom_set, natom)
      CALL allocate_task_list(ri_gapw_env%task_list)
      ALLOCATE(ri_gapw_env%excited_kinds(nex_kinds))
      ri_gapw_env%excited_kinds = excited_kinds
      ALLOCATE(ri_gapw_env%grid_atom_set(nkind))
      ALLOCATE(ri_gapw_env%harmonics_atom_set(nkind))
      ALLOCATE(ri_gapw_env%blk_size(natom))
      ALLOCATE(ri_gapw_env%sphi_s(nkind))
      ALLOCATE(ri_gapw_env%coeff_h(natom))
      ALLOCATE(ri_gapw_env%coeff_s(natom))

!  Allocate and initialize the atomic grids and harmonics
!  lmax_sphere should be at least maxlgto and max 2*maxlgto
      lmax_sphere = 60
      CALL init_ri_gapw_grid_harmo(lmax_sphere, ri_gapw_env, qs_env)

!  Go over all qs_kind and make sure the RI_XAS_SOFT basis is computed and stored
      DO ikind = 1,nkind
         IF (.NOT. ANY(excited_kinds == ikind)) CYCLE
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=soft_basis, basis_type="RI_XAS_SOFT")
         IF (.NOT. ASSOCIATED(soft_basis)) THEN
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS", & 
                             hard_radius=rc) 
            CALL allocate_gto_basis_set(soft_basis)
            CALL create_soft_basis(orb_basis=ri_basis, soft_basis=soft_basis, &
                                  eps_fit=eps_fit, rc=rc, paw_atom=paw_atom, &
                                  paw_type_forced=.TRUE., gpw_type_forced=.FALSE.) 
            CALL init_interaction_radii_orb_basis(soft_basis, qs_control%eps_pgf_orb)
            CALL add_basis_set_to_container(qs_kind_set(ikind)%basis_sets, &
                                            soft_basis, "RI_XAS_SOFT")
         END IF
         NULLIFY(soft_basis)

!     Compute the sphi_s coeff for the soft basis, in the same format as the sphi_h
         NULLIFY(sphi_s)
         CALL get_gto_basis_set(ri_basis, sphi=sphi_h, nsgf=nsgf, lmax=lmax, npgf=npgf, &
                                nshell=nshell, l=l, first_sgf=first_sgf, last_sgf=last_sgf, &
                                nset=nset)
         ALLOCATE(ri_gapw_env%sphi_s(ikind)%array(SIZE(sphi_h,1),SIZE(sphi_h,2)))
         sphi_s => ri_gapw_env%sphi_s(ikind)%array
         sphi_s = 0.0_dp
         DO iset = 1,nset
            n = ncoset(lmax(iset))
            DO ipgf = 1, npgf(iset)
               DO ishell = 1, nshell(iset)
                  lshell = l(ishell, iset)
                  icomin = ncoset(lshell-1)+1+n*(ipgf-1)
                  icomax = ncoset(lshell)+n*(ipgf-1)
                  isgfmin = first_sgf(ishell, iset)
                  isgfmax = last_sgf(ishell, iset)
                  radius = exp_radius(lshell, ri_basis%zet(ipgf, iset), eps_fit, 1.0_dp)
                  DO isgf = isgfmin, isgfmax
                     IF (radius > rc) sphi_s(icomin:icomax, isgf) = sphi_h(icomin:icomax, isgf)
                  END DO ! isgf
               END DO ! ishell
            END DO ! ipgf
         END DO ! iset

      END DO !ikind

!  Create a neighbor list based on the RI_XAS basis with the minimal image convention since
!  we only consider atom idependentely, we want to avoid the overlap with its image. Also this list
!  is such that atoms are only neighbor with themselves. Note: also valid for the RI_XAS_SOFT
      ALLOCATE(basis_set_a(nkind))
      ALLOCATE(radius_backup(nkind))
      CALL basis_set_list_setup(basis_set_a, "RI_XAS", qs_kind_set)
      DO ikind = 1,nkind
         IF (ASSOCIATED(basis_set_a(ikind)%gto_basis_set)) THEN
            CALL get_gto_basis_set(basis_set_a(ikind)%gto_basis_set, kind_radius=radius_backup(ikind))
            CALL set_gto_basis_set(basis_set_a(ikind)%gto_basis_set, kind_radius=0.1_dp) 
         END IF
      END DO
      CALL setup_neighbor_list(ri_gapw_env%sri, basis_set_a, qs_env=qs_env, mic=.TRUE.) 
      DO ikind = 1,nkind
         IF (ASSOCIATED(basis_set_a(ikind)%gto_basis_set)) THEN
            CALL set_gto_basis_set(basis_set_a(ikind)%gto_basis_set, kind_radius=radius_backup(ikind))
         END IF
      END DO

!  Build our own pw_env, based on the qs one, so that it accomodates the RI_XAS_SOFT basis
      CALL pw_env_create(ri_gapw_env%pw_env)
      CALL pw_env_rebuild(ri_gapw_env%pw_env, qs_env)

!  Create the task_list based on the neighbor list above. Keep input param
      CALL generate_qs_task_list(ks_env, ri_gapw_env%task_list, reorder_rs_grid_ranks=.TRUE., & 
                                 soft_valid=.FALSE., basis_type="RI_XAS_SOFT", &
                                 skip_load_balance_distributed=qs_control%skip_load_balance_distributed, &
                                 sab_orb_external=ri_gapw_env%sri, pw_env_external=ri_gapw_env%pw_env)

!  Create the block size for the integral matrix based on the size of the RI_XAS basis.
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=ri_gapw_env%blk_size, basis=basis_set_a)
      DEALLOCATE(basis_set_a)

!  Get the coefficients to expand the density in the RI_XAS and the RI_XAS_SOFT basis sets
      CALL calculate_density_coeffs(ri_gapw_env, qs_env)  

      CALL timestop(handle)

   END SUBROUTINE init_ri_gapw_env

! **************************************************************************************************
!> \brief Releases the ri_gapw_env type
!> \param ri_gapw_env the type to release
! **************************************************************************************************
   SUBROUTINE release_ri_gapw_env(ri_gapw_env)

      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env              

      CHARACTER(len=*), PARAMETER :: routineN = "release_ri_gapw_env", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: i

      IF (ASSOCIATED(ri_gapw_env%rho_atom_set)) THEN
         CALL deallocate_rho_atom_set(ri_gapw_env%rho_atom_set)
      END IF

      IF (ASSOCIATED(ri_gapw_env%grid_atom_set)) THEN
         DO i = 1,SIZE(ri_gapw_env%grid_atom_set) 
            IF (ASSOCIATED(ri_gapw_env%grid_atom_set(i)%grid_atom)) THEN
               CALL deallocate_grid_atom(ri_gapw_env%grid_atom_set(i)%grid_atom)
            END IF
         END DO
         DEALLOCATE(ri_gapw_env%grid_atom_set)
      END IF

      IF (ASSOCIATED(ri_gapw_env%harmonics_atom_set)) THEN
         DO i = 1,SIZE(ri_gapw_env%harmonics_atom_set)
            IF (ASSOCIATED(ri_gapw_env%harmonics_atom_set(i)%harmonics_atom)) THEN
               CALL deallocate_harmonics_atom(ri_gapw_env%harmonics_atom_set(i)%harmonics_atom)
            END IF
         END DO
         DEALLOCATE(ri_gapw_env%harmonics_atom_set)
      END IF

      IF (ASSOCIATED(ri_gapw_env%task_list)) THEN
         CALL deallocate_task_list(ri_gapw_env%task_list)
      END IF

      IF (ASSOCIATED(ri_gapw_env%excited_kinds)) THEN
         DEALLOCATE(ri_gapw_env%excited_kinds)
      END IF

      IF (ASSOCIATED(ri_gapw_env%blk_size)) THEN
         DEALLOCATE(ri_gapw_env%blk_size)
      END IF

      IF (ASSOCIATED(ri_gapw_env%sphi_s)) THEN
         DO i =1,SIZE(ri_gapw_env%sphi_s)
            IF (ASSOCIATED(ri_gapw_env%sphi_s(i)%array)) THEN
               DEALLOCATE(ri_gapw_env%sphi_s(i)%array)
            END IF
         END DO
         DEALLOCATE(ri_gapw_env%sphi_s)
      END IF

      IF (ASSOCIATED(ri_gapw_env%coeff_h)) THEN
         DO i = 1,SIZE(ri_gapw_env%coeff_h)
            IF (ASSOCIATED(ri_gapw_env%coeff_h(i)%array)) THEN
               DEALLOCATE(ri_gapw_env%coeff_h(i)%array)
            END IF
         END DO
         DEALLOCATE(ri_gapw_env%coeff_h)
      END IF

      IF (ASSOCIATED(ri_gapw_env%coeff_s)) THEN 
         DO i = 1,SIZE(ri_gapw_env%coeff_s)
            IF (ASSOCIATED(ri_gapw_env%coeff_s(i)%array)) THEN
               DEALLOCATE(ri_gapw_env%coeff_s(i)%array)
            END IF
         END DO
         DEALLOCATE(ri_gapw_env%coeff_s)
      END IF

      IF (ASSOCIATED(ri_gapw_env%sri)) THEN
         DO i = 1,SIZE(ri_gapw_env%sri)
            CALL deallocate_neighbor_list_set(ri_gapw_env%sri(i)%neighbor_list_set)
         END DO
         DEALLOCATE(ri_gapw_env%sri)
      END IF

      IF (ASSOCIATED(ri_gapw_env%pw_env)) THEN
         CALL pw_env_release(ri_gapw_env%pw_env)
      END IF

      DEALLOCATE(ri_gapw_env)

   END SUBROUTINE release_ri_gapw_env

! **************************************************************************************************
!> \brief Initializes the atomic grids and harmonics for the RI gapw calculation
!> \param lmax_sphere maximum L for the Clebsch-Gordan coeffs
!> \param ri_gapw_env the RI gapw environement 
!> \param qs_env ...
!> \note Largely inspired by init_rho_atom subroutine (which is sadly restricted to orbital basis set)
! **************************************************************************************************
   SUBROUTINE init_ri_gapw_grid_harmo(lmax_sphere, ri_gapw_env, qs_env)

      INTEGER                                         :: lmax_sphere
      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "init_ri_gapw_grid_harmo", &
                                     routineP = moduleN//":"//routineN

      REAL(dp), DIMENSION(:, :, :), POINTER           :: my_CG
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)          :: rga
      INTEGER                                         :: llmax, maxlgto, max_s_harm, max_s_set, &
                                                         lcleb, lc1, iso1, l1, m1, lc2, iso2, &
                                                         l2, m2, l1l2, mp, mm, lp, il, iso, &
                                                         ikind, quadrature, i, ll, na, la, nr, &
                                                         maxs, maxl
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER               :: tmp_basis
      TYPE(grid_atom_type), POINTER                   :: grid_atom
      TYPE(harmonics_atom_type), POINTER              :: harmonics
      TYPE(qs_control_type), POINTER                  :: qs_control
      TYPE(dft_control_type), POINTER                 :: dft_control

      NULLIFY(my_CG, qs_kind_set, tmp_basis, grid_atom, harmonics, qs_control, dft_control)

!  Initialization of some integer for the CG coeff generation
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control)
      CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="RI_XAS")
      qs_control => dft_control%qs_control

      llmax = MIN(lmax_sphere, 2*maxlgto)
      max_s_harm = nsoset(llmax)
      max_s_set = nsoset(maxlgto)
      lcleb = MAX(llmax, 2*maxlgto, 1)

!  Allocate and compute the CG coeffs (copied from init_rho_atom)
      CALL clebsch_gordon_init(lcleb)
      CALL reallocate(my_CG, 1, max_s_set, 1, max_s_set, 1, max_s_harm)

      ALLOCATE (rga(lcleb, 2))
      DO lc1 = 0, maxlgto
         DO iso1 = nsoset(lc1-1)+1, nsoset(lc1)
            l1 = indso(1, iso1)
            m1 = indso(2, iso1)
            DO lc2 = 0, maxlgto
               DO iso2 = nsoset(lc2-1)+1, nsoset(lc2)
                  l2 = indso(1, iso2)
                  m2 = indso(2, iso2)
                  CALL clebsch_gordon(l1, m1, l2, m2, rga)
                  IF (l1+l2 > llmax) THEN
                     l1l2 = llmax
                  ELSE
                     l1l2 = l1+l2
                  END IF
                  mp = m1+m2
                  mm = m1-m2
                  IF (m1*m2 < 0 .OR. (m1*m2 == 0 .AND. (m1 < 0 .OR. m2 < 0))) THEN
                     mp = -ABS(mp)
                     mm = -ABS(mm)
                  ELSE
                     mp = ABS(mp)
                     mm = ABS(mm)
                  END IF
                  DO lp = MOD(l1+l2, 2), l1l2, 2
                     il = lp/2+1
                     IF (ABS(mp) <= lp) THEN
                     IF (mp >= 0) THEN
                        iso = nsoset(lp-1)+lp+1+mp
                     ELSE
                        iso = nsoset(lp-1)+lp+1-ABS(mp)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 1)
                     ENDIF
                     IF (mp /= mm .AND. ABS(mm) <= lp) THEN
                     IF (mm >= 0) THEN
                        iso = nsoset(lp-1)+lp+1+mm
                     ELSE
                        iso = nsoset(lp-1)+lp+1-ABS(mm)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 2)
                     ENDIF
                  END DO
               ENDDO ! iso2
            ENDDO ! lc2
         ENDDO ! iso1
      ENDDO ! lc1
      DEALLOCATE (rga)
      CALL clebsch_gordon_deallocate()

!  Create the Lebedev grids and compute the spherical harmonics
      CALL init_lebedev_grids()
      quadrature = qs_control%gapw_control%quadrature
      
      DO i = 1,SIZE(ri_gapw_env%excited_kinds)
         ikind = ri_gapw_env%excited_kinds(i)

!        Allocate the grid and the harmonics for this kind
         CALL allocate_grid_atom(ri_gapw_env%grid_atom_set(ikind)%grid_atom)
         CALL allocate_harmonics_atom(ri_gapw_env%harmonics_atom_set(ikind)%harmonics_atom)

         NULLIFY(grid_atom, harmonics)
         grid_atom => ri_gapw_env%grid_atom_set(ikind)%grid_atom         
         harmonics => ri_gapw_env%harmonics_atom_set(ikind)%harmonics_atom 

!        Initialize some integers
         CALL get_qs_kind(qs_kind_set(ikind), ngrid_rad=nr, ngrid_ang=na)
         ll = get_number_of_lebedev_grid(n=na)
         na = lebedev_grid(ll)%n
         la = lebedev_grid(ll)%l
         grid_atom%ng_sphere = na
         grid_atom%nr = nr

         IF (llmax > la) THEN
            CPWARN("The Lebedev grid for XAS TDP RI GAPW is built for a too low angular momentum l")
            CPWARN("To fix this, either increase LEBEDEV_GRID or LMAXSPHERE.")
         END IF

         CALL create_grid_atom(grid_atom, nr, na, llmax, ll, quadrature)

         CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="RI_XAS")
         CALL get_gto_basis_set(gto_basis_set=tmp_basis, maxl=maxl)
         grid_atom%ng_sphere = na
         grid_atom%nr = nr
         maxs = nsoset(maxl)
         CALL create_harmonics_atom(harmonics, &
                                    my_CG, na, llmax, maxs, max_s_harm, ll, grid_atom%wa, &
                                    grid_atom%azi, grid_atom%pol) 
         CALL get_maxl_CG(harmonics, tmp_basis, llmax, max_s_harm)

      END DO

      CALL deallocate_lebedev_grids()
      DEALLOCATE(my_CG)

   END SUBROUTINE init_ri_gapw_grid_harmo

! **************************************************************************************************
!> \brief Compute the coefficients to project the density on the RI_XAS and the RI_XAS_SOFT basis sets
!> \param ri_gapw_env ...
!> \param qs_env ...
!> \note The density is n = sum_ab P_ab*phi_a*phi_b, the RI basis covers the products of orbital sgfs
!>       => n = sum_ab sum_cd P_ab (phi_a phi_b xi_c) S_cd^-1 xi_d
!>            = sum_d coeff_d xi_d , where xi are the RI basis func.
!>       The 3-center overlap is needed, as well as the inverse of the RI basis overlap. Note that 
!>       the RI basis functions are only centered on atoms of excited kinds (for now)
! **************************************************************************************************
   SUBROUTINE calculate_density_coeffs(ri_gapw_env, qs_env)

      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "calculate_ri_gapw_densities", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nkind, ilist, katom, kkind, ikind, natom, &
                                                         iat
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_ri_h, basis_set_orb, &
                                                         basis_set_ri_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &                                                
                                             POINTER  :: ac_list, sab_orb
      TYPE(o3c_container_type), POINTER               :: o3c_h, o3c_s
      TYPE(qs_rho_type), POINTER                      :: rho
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ao
      TYPE(o3c_iterator_type)                         :: o3c_iterator
      REAL(dp), DIMENSION(:,:), POINTER               :: tvec
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(qs_ks_env_type), POINTER                   :: ks_env
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER     :: sinv_s, sinv_h
      REAL(dp), DIMENSION(:,:), POINTER               :: block
      LOGICAL                                         :: found
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set 
      INTEGER, DIMENSION(:), POINTER                  :: nsgf
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: work

      NULLIFY(qs_kind_set, basis_set_ri_h, basis_set_ri_s, basis_set_orb, ac_list, rho, rho_ao)
      NULLIFY(o3c_h, o3c_s, tvec, matrix_s, ks_env, particle_set, nsgf)
      CPASSERT(ASSOCIATED(ri_gapw_env))

      CALL get_qs_env(qs_env, nkind=nkind, sab_orb=sab_orb, qs_kind_set=qs_kind_set, rho=rho, &
                      ks_env=ks_env, natom=natom, particle_set=particle_set)

!  TODO: Take also care of contribution from neighboring atoms by doing some kind of LRI scheme
!        The contribution from neighboring atoms will be smooth and should contribute the same way
!        to the hard and the soft densities (so that the sum/diff follow the rules)

!  Create the needed neighbor list and basis set lists. Neighbor list is the same for soft and hard
      ALLOCATE(basis_set_ri_h(nkind))
      ALLOCATE(basis_set_ri_s(nkind))
      ALLOCATE(basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri_h, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_ri_s, "RI_XAS_SOFT", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)
      CALL setup_neighbor_list(ac_list, basis_set_ri_h, basis_set_orb, qs_env) 

!  Create the inverse RI overlap matrices. Block diagonal because of only own neighbor list
!  Inverse block by block using SVD because not always positive definite with soft basis. 
!  TODO: need to be carful when parallelizing
      ALLOCATE(sinv_h(nkind))
      CALL build_overlap_matrix(ks_env=ks_env, matrix_s=matrix_s, basis_type_a="RI_XAS", &
                                basis_type_b="RI_XAS", sab_nl=ri_gapw_env%sri)
      DO ikind = 1,nkind
         found = .FALSE.
         CALL dbcsr_get_block_p(matrix_s(1)%matrix, ikind, ikind, block, found)
         IF (found) THEN
            ALLOCATE(sinv_h(ikind)%array(SIZE(block,1),SIZE(block,2)))
            CALL get_pseudo_inverse_svd(block, sinv_h(ikind)%array, 1.0E-8_dp)
         END IF
      END DO
      CALL dbcsr_release_p(matrix_s(1)%matrix)
      DEALLOCATE(matrix_s)

      ALLOCATE(sinv_s(nkind))
      CALL build_overlap_matrix(ks_env=ks_env, matrix_s=matrix_s, basis_type_a="RI_XAS_SOFT", & 
                                basis_type_b="RI_XAS_SOFT", sab_nl=ri_gapw_env%sri)
      DO ikind = 1,nkind
         found = .FALSE.
         CALL dbcsr_get_block_p(matrix_s(1)%matrix, ikind, ikind, block, found)
         IF (found) THEN
            ALLOCATE(sinv_s(ikind)%array(SIZE(block,1),SIZE(block,2)))
            CALL get_pseudo_inverse_svd(block, sinv_s(ikind)%array, 1.0E-8_dp)
         END IF
      END DO
      CALL dbcsr_release_p(matrix_s(1)%matrix)
      DEALLOCATE(matrix_s)

!  Create the hard and soft o3c 
      ALLOCATE(o3c_h, o3c_s)
      CALL init_o3c_container(o3c_h, 1, basis_set_orb, basis_set_orb, basis_set_ri_h, &
                              sab_orb, ac_list)
      CALL init_o3c_container(o3c_s, 1, basis_set_orb, basis_set_orb, basis_set_ri_s, &
                              sab_orb, ac_list)

!  Compute the overlaps
      CALL calculate_o3c_integrals(o3c_h)
      CALL calculate_o3c_integrals(o3c_s)

!  Contract the 3-center tensor with the density matrix
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL contract12_o3c(o3c_h, rho_ao)
      CALL contract12_o3c(o3c_s, rho_ao)

!  Allocate space for the projected density coefficients. Same nsgf for hard and soft by construction
      ALLOCATE(nsgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=nsgf, basis=basis_set_ri_h)
      DO iat = 1,natom
         ikind = particle_set(iat)%atomic_kind%kind_number
         IF (.NOT. ANY(ri_gapw_env%excited_kinds == ikind)) CYCLE
         ALLOCATE(ri_gapw_env%coeff_h(iat)%array(nsgf(iat)))
         ALLOCATE(ri_gapw_env%coeff_s(iat)%array(nsgf(iat)))
         ri_gapw_env%coeff_h(iat)%array = 0.0_dp
         ri_gapw_env%coeff_s(iat)%array = 0.0_dp
      END DO

!  Iterate over the o3c to access the contracted tvec = sum_ab P_ab*(abc). Hard then soft
!  TODO: careful whne parallelizing
      CALL o3c_iterator_create(o3c_h, o3c_iterator) 
      DO WHILE (o3c_iterate(o3c_iterator) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, katom=katom, kkind=kkind, tvec=tvec)
         IF (.NOT. ANY(ri_gapw_env%excited_kinds == kkind)) CYCLE

!        The coefficient for atom katom of kkind are tvec*S^-1
         ALLOCATE(work(nsgf(katom),1))
         CALL dgemm('N', 'N', nsgf(katom), 1, nsgf(katom), 1.0_dp, sinv_h(kkind)%array, nsgf(katom), &
                    tvec, nsgf(katom), 0.0_dp, work, nsgf(katom))
         ri_gapw_env%coeff_h(katom)%array(:) = ri_gapw_env%coeff_h(katom)%array(:) + work(:,1)
         DEALLOCATE(work)

      END DO
      CALL o3c_iterator_release(o3c_iterator)

      CALL o3c_iterator_create(o3c_s, o3c_iterator)
      DO WHILE (o3c_iterate(o3c_iterator) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, katom=katom, kkind=kkind, tvec=tvec)
         IF (.NOT. ANY(ri_gapw_env%excited_kinds == kkind)) CYCLE

!        The coefficient for atom katom of kkind are tvec*S^-1 
         ALLOCATE(work(nsgf(katom),1))
         CALL dgemm('N', 'N', nsgf(katom), 1, nsgf(katom), 1.0_dp, sinv_s(kkind)%array, nsgf(katom), &
                    tvec, nsgf(katom), 0.0_dp, work, nsgf(katom))
         ri_gapw_env%coeff_s(katom)%array(:) = ri_gapw_env%coeff_s(katom)%array(:) + work(:,1)
         DEALLOCATE(work)

      END DO
      CALL o3c_iterator_release(o3c_iterator)

!  some clean-up
      DO ilist = 1,SIZE(ac_list)
         CALL deallocate_neighbor_list_set(ac_list(ilist)%neighbor_list_set)
      END DO
      DO ikind = 1,nkind
         IF (ASSOCIATED(sinv_h(ikind)%array)) DEALLOCATE(sinv_h(ikind)%array)
         IF (ASSOCIATED(sinv_s(ikind)%array)) DEALLOCATE(sinv_s(ikind)%array)
      END DO
      DEALLOCATE(sinv_h, sinv_s)
      CALL release_o3c_container(o3c_h)
      CALL release_o3c_container(o3c_s)
      DEALLOCATE(basis_set_ri_h, basis_set_ri_s, basis_set_orb, o3c_h, o3c_s, nsgf, ac_list)

   END SUBROUTINE calculate_density_coeffs

! **************************************************************************************************
!> \brief Evaluates the density on a given atomic grid
!> \param n_h the array containing the hard density at each grid point
!> \param n_s the array containing the soft density at each grid point
!> \param atom_index the index of the atom in question
!> \param atom_kind the kind of the atom in question
!> \param ri_gape_env ...
!> \param qs_env ...
!> \note The density is expressed as n = sum_a coeff_a*xi_a. Knowing the coordinate of each grid
!>       grid point, one can simply evaluate xi_(r) 
! **************************************************************************************************
   SUBROUTINE put_density_on_atomic_grid(n_h, n_s, atom_index, atom_kind, ri_gapw_env, qs_env)

      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: n_h, n_s
      INTEGER, INTENT(IN)                             :: atom_index, atom_kind
      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "put_density_on_atomic_grid", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: na, nr, iset, ipgf, l, ico, iso, &
                                                         lx, ly, lz, ia, nset, isgf, start_c, &
                                                         start_s, sgfi
      TYPE(grid_atom_type), POINTER                   :: grid_atom
      TYPE(harmonics_atom_type), POINTER              :: harmonics
      TYPE(gto_basis_set_type), POINTER               :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      INTEGER, DIMENSION(:), POINTER                  :: lmax, npgf, nsgf_set
      REAL(dp), DIMENSION(:,:), POINTER               :: zet, sphi_h, sphi_s
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: sg, cg, sgf_h, sgf_s
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgf

      NULLIFY(grid_atom, ri_basis, qs_kind_set, harmonics, lmax, npgf, zet, nsgf_set, sphi_h, first_sgf)

!  Strategy: it makes sense to evaluate the so on the grid: first on the radial part (all at once)
!            then on the angular part (multiplying by the spherical harmonics). From there, one
!            should transform into co, and contract to sgf. Once we have the value of each sgf on 
!            the whole grid, one can take their weight

!  Generalities
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(atom_kind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, nset=nset, nsgf_set=nsgf_set,&
                             sphi=sphi_h, first_sgf=first_sgf)
      sphi_s => ri_gapw_env%sphi_s(atom_kind)%array

!  Get the grid, the harmonics and the info we need from it
      grid_atom => ri_gapw_env%grid_atom_set(atom_kind)%grid_atom
      harmonics => ri_gapw_env%harmonics_atom_set(atom_kind)%harmonics_atom
      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(harmonics))
      nr = grid_atom%nr
      na = grid_atom%ng_sphere

!  Warn if not enough spherical harmonics for accurate representation of the density
      IF (nsoset(MAXVAL(lmax)) > harmonics%max_iso_not0) THEN
         CPWARN("Not enough Spherical Harmonics for proper representation of the density")
      END IF

!  Loop over the sets and pgfs of the basis sets (same for hard and soft by construction)
      DO iset = 1,nset

         !allocate space to store the spherical and cartesian orbitals on the grid
         ALLOCATE(sg(npgf(iset)*nsoset(lmax(iset)),na,nr))
         ALLOCATE(cg(npgf(iset)*ncoset(lmax(iset)),na,nr))
         sg = 0.0_dp
         cg = 0.0_dp

         DO ipgf = 1,npgf(iset)
            start_s = (ipgf-1)*nsoset(lmax(iset))
            start_c = (ipgf-1)*ncoset(lmax(iset)) 

            !loop over the spherical gaussian orbitals 
            DO iso = 1,nsoset(lmax(iset))
               l = indso(1, iso)
               DO ia = 1,na
                  !compute: g = r**l * Y_lm * exp(-zet* r**2)
                  sg(start_s+iso,ia,1:nr) = EXP(-zet(ipgf, iset)*grid_atom%rad2(1:nr)) & 
                                            *grid_atom%rad(1:nr)**l *harmonics%slm(ia,iso)
               END DO !ia
            END DO !iso

            !transform into cartesian gaussian orbitals
            DO l = 0,lmax(iset)
               DO ico = 1,nco(l)
                  lx = indco(1, ico+ncoset(l-1))
                  ly = indco(2, ico+ncoset(l-1))
                  lz = indco(3, ico+ncoset(l-1))
                  DO iso = 1,nso(l)

                     cg(start_c+ncoset(l-1)+ico,:,:) = cg(start_c+ncoset(l-1)+ico,:,:) &
                                                     + sg(start_s+nsoset(l-1)+iso,:,:) &
                                                     * orbtramat(l)%s2c(iso,ico)       &
                              * SQRT(4.0_dp*pi/dfac(2*l+1)*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))

                  END DO !iso
               END DO !ico
            END DO !l

         END DO !ipgf

         !contract the co to get the sgf in both basis
         ALLOCATE(sgf_h(nsgf_set(iset),na,nr))
         ALLOCATE(sgf_s(nsgf_set(iset),na,nr))
         sgf_h = 0.0_dp
         sgf_s = 0.0_dp
         sgfi = first_sgf(1,iset)

         !TODO: can we use dgemm here ?
         DO isgf = 1,nsgf_set(iset)
            DO ico = 1,npgf(iset)*ncoset(lmax(iset))

               sgf_h(isgf,:,:) = sgf_h(isgf,:,:) + sphi_h(ico,sgfi+isgf-1)*cg(ico,:,:)
               sgf_s(isgf,:,:) = sgf_s(isgf,:,:) + sphi_s(ico,sgfi+isgf-1)*cg(ico,:,:)
               
            END DO !ico

            !put the sgf on the grid with the approriate coefficients and sum
            n_h(:,:,1) = n_h(:,:,1) + sgf_h(isgf,:,:)*ri_gapw_env%coeff_h(atom_index)%array(sgfi+isgf-1)
            n_s(:,:,1) = n_s(:,:,1) + sgf_s(isgf,:,:)*ri_gapw_env%coeff_s(atom_index)%array(sgfi+isgf-1)

         END DO !isgf

         DEALLOCATE(sg, cg, sgf_h, sgf_s)
      END DO !iset

   END SUBROUTINE put_density_on_atomic_grid

! **************************************************************************************************
!> \brief Integrate a funcion of r on the atomic grids, for the hard and the soft RI_XAS basis sets
!> \param int_mat_h the dbcsr_matrix in the usual block format where the hard integrals are stored
!> \param int_mat_s the dbcsr_matrix in the usual block format where the soft integrals are stored
!> \param singlet whether we want singlet or triplet excitations
!> \param ri_gapw_env ...
!> \param qs_env ...
!> \note test phase: we allocate/deallcoate everything locally
! **************************************************************************************************
   SUBROUTINE integrate_ri_gapw_atoms(int_mat_h, int_mat_s, singlet, ri_gapw_env, qs_env)

      TYPE(dbcsr_p_type)                              :: int_mat_h
      TYPE(dbcsr_p_type)                              :: int_mat_s
      LOGICAL, INTENT(IN)                             :: singlet
      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env  
      TYPE(qs_environment_type), POINTER              :: qs_env   

      CHARACTER(len=*), PARAMETER :: routineN = "integrate_ri_gapw_atoms", &
                                     routineP = moduleN//":"//routineN

      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      REAL(dp), DIMENSION(:, :, :), ALLOCATABLE       :: rho_h, rho_s
      REAL(dp), DIMENSION(:, :, :), POINTER           :: fxc_h, fxc_s
      INTEGER                                         :: natom, iat, ikind, na, nr, nsotot, nset, &
                                                         maxso, iset, jset, nsgf, ia, ir, sgfi, &
                                                         sgfj, maxco, ncoi, ncoj, ipgf, jpgf, &
                                                         starti_s, startj_s, starti_c, startj_c, &
                                                         ncotot, li, lj, n
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set 
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(grid_atom_type), POINTER                   :: grid_atom 
      TYPE(harmonics_atom_type), POINTER              :: harmonics
      TYPE(gto_basis_set_type), POINTER               :: ri_basis 
      TYPE(rho_atom_type), POINTER                    :: rho_atom
      REAL(dp), DIMENSION(:,:), POINTER               :: int_s_h, int_s_s, sphi_h, sphi_s
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: int_h, int_s, int_c_h, int_c_s
      INTEGER, DIMENSION(:), POINTER                  :: nsgf_set, npgf, lmax
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgf
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: rs_h, rs_s
      TYPE(dft_control_type), POINTER                 :: dft_control

      NULLIFY(dbcsr_dist, particle_set, qs_kind_set, grid_atom, harmonics, rho_atom, dft_control)
      NULLIFY(ri_basis, nsgf_set, first_sgf, lmax, npgf, sphi_h, sphi_s, int_s_h, int_s_s)
      NULLIFY(fxc_h, fxc_s)

!  Create the integral matrices based on the usual dbcsr_dist and the custom ri_gapw_env block sizes
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_create(matrix=int_mat_h%matrix, name="RI GAPW HARD INT", matrix_type='S', &         
                        dist=dbcsr_dist, row_blk_size=ri_gapw_env%blk_size, &                        
                        col_blk_size=ri_gapw_env%blk_size)
      CALL dbcsr_create(matrix=int_mat_s%matrix, name="RI GAPW SOFT INT", matrix_type='S', &         
                        dist=dbcsr_dist, row_blk_size=ri_gapw_env%blk_size, &                        
                        col_blk_size=ri_gapw_env%blk_size)

      CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set, natom=natom, &
                      dft_control=dft_control)

!  TODO: think carfully when parallelizing
!  Loop over the atoms of interest. 
      DO iat = 1, natom
         ikind = particle_set(iat)%atomic_kind%kind_number
         IF (.NOT. ANY(ri_gapw_env%excited_kinds == ikind)) CYCLE

         grid_atom => ri_gapw_env%grid_atom_set(ikind)%grid_atom
         harmonics => ri_gapw_env%harmonics_atom_set(ikind)%harmonics_atom
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")
         rho_atom => ri_gapw_env%rho_atom_set(iat)

!     Allocate the density, kernel and integral arrays
         na = grid_atom%ng_sphere
         nr = grid_atom%nr
         CALL get_gto_basis_set(ri_basis, nset=nset, maxso=maxso, maxco=maxco)
         nsotot = nset*maxso 

         ALLOCATE(rho_h(na, nr, 1)) ! 1 spin
         ALLOCATE(rho_s(na, nr, 1)) ! 1 spin

         ALLOCATE(fxc_h(na, nr, 1)) ! 1 spin                                                         
         ALLOCATE(fxc_s(na, nr, 1)) ! 1 spin

         NULLIFY(rho_atom%ga_Vlocal_gb_h)
         NULLIFY(rho_atom%ga_Vlocal_gb_s)
         ALLOCATE(rho_atom%ga_Vlocal_gb_h(1))
         ALLOCATE(rho_atom%ga_Vlocal_gb_s(1))
         NULLIFY(rho_atom%ga_Vlocal_gb_h(1)%r_coef)
         NULLIFY(rho_atom%ga_Vlocal_gb_s(1)%r_coef)
         ALLOCATE(rho_atom%ga_Vlocal_gb_h(1)%r_coef(nsotot, nsotot))
         ALLOCATE(rho_atom%ga_Vlocal_gb_s(1)%r_coef(nsotot, nsotot))

!     Put the density on the grid
         CALL put_density_on_atomic_grid(rho_h, rho_s, iat, ikind, ri_gapw_env, qs_env)

!     Compute the fxc kernel on the atomic grids. Use calc_rs_pw because it takes 3d arrays as arg
         CALL set_util(dft_control%qs_control%eps_rho_rspace)
         n = PRODUCT(SHAPE(rho_h))
         ALLOCATE(rs_h(n))
         ALLOCATE(rs_s(n))

         CALL calc_rs_pw(rho_h, rs_h, n)
         CALL calc_rs_pw(rho_s, rs_s, n)
         fxc_h = 0.0_dp
         fxc_s = 0.0_dp

         IF (singlet) THEN
            CALL pade_kernel_singlet(n, rho_h, rs_h, fxc_h, dft_control%qs_control%eps_rho_rspace)
            CALL pade_kernel_singlet(n, rho_s, rs_s, fxc_s, dft_control%qs_control%eps_rho_rspace)
         ELSE
            CALL pade_kernel_triplet(n, rho_h, rs_h, fxc_h, dft_control%qs_control%eps_rho_rspace)
            CALL pade_kernel_triplet(n, rho_s, rs_s, fxc_s, dft_control%qs_control%eps_rho_rspace)
         END IF

!     Take the grid weights into account
         DO ia = 1,na
            DO ir = 1,nr
               fxc_h(ia,ir,1) = fxc_h(ia,ir,1)*grid_atom%wa(ia)*grid_atom%wr(ir)
               fxc_s(ia,ir,1) = fxc_s(ia,ir,1)*grid_atom%wa(ia)*grid_atom%wr(ir)
            END DO
         END DO

!     Integrate
         CALL gaVxcgb_noGC(fxc_h, fxc_s, qs_kind_set(ikind), rho_atom, 1, ri_basis, harmonics, &
                           grid_atom)

         DEALLOCATE(rs_h, rs_s, fxc_h, fxc_s)

!     Turn the integrals (in spherical orbital representation) into something we can work with
         CALL get_gto_basis_set(ri_basis, nsgf=nsgf, nsgf_set=nsgf_set, first_sgf=first_sgf, &
                                npgf=npgf, lmax=lmax, sphi=sphi_h)

!     Initializing the arrays for the new integral representations
         NULLIFY(int_s_h, int_s_s)
!        integrals in co
         ncotot = nset*maxco
         ALLOCATE(int_c_h(ncotot,ncotot))
         ALLOCATE(int_c_s(ncotot,ncotot))
         int_c_h = 0.0_dp
         int_c_s = 0.0_dp
!        integrals in so
         int_s_h => rho_atom%ga_Vlocal_gb_h(1)%r_coef
         int_s_s => rho_atom%ga_Vlocal_gb_s(1)%r_coef
!        integrals in sgf
         ALLOCATE(int_h(nsgf,nsgf))
         ALLOCATE(int_s(nsgf,nsgf))
         int_h = 0.0_dp
         int_s = 0.0_dp

!     Switch from so integrals to co integrals. Subscript c => cartesian, s => spherical
!     It's symmetric, only compute upper triangle
         DO iset = 1,nset
            starti_s = (iset-1)*maxso+1
            starti_c = (iset-1)*maxco+1
            li = lmax(iset)
            DO ipgf = 1,npgf(iset)

               DO jset = iset,nset
                  startj_s = (jset-1)*maxso+1
                  startj_c = (jset-1)*maxco+1
                  lj = lmax(jset)
                  DO jpgf = 1,npgf(jset)

                     ! hard integrals from spherical to cartesian orbitals
                     CALL transform_s2c( &
                          int_s_h(starti_s:starti_s+nsoset(li)-1,startj_s:startj_s+nsoset(lj)-1), &
                          int_c_h(starti_c:starti_c+ncoset(li)-1,startj_c:startj_c+ncoset(lj)-1), &
                          MAX(li,lj), li, lj)

                     ! soft integrals from spherical to cartesian orbitals
                     CALL transform_s2c( & 
                          int_s_s(starti_s:starti_s+nsoset(li)-1,startj_s:startj_s+nsoset(lj)-1), &
                          int_c_s(starti_c:starti_c+ncoset(li)-1,startj_c:startj_c+ncoset(lj)-1), &
                          MAX(li,lj), li, lj)

                     startj_s = startj_s + nsoset(lj)
                     startj_c = startj_c + ncoset(lj)
                  END DO !jpgf
               END DO !jset

               starti_s = starti_s + nsoset(li)
               starti_c = starti_c + ncoset(li)
            END DO! ipgf
         END DO !iset

!     Contract into sgf 
         sphi_s => ri_gapw_env%sphi_s(ikind)%array
         DO iset = 1,nset
            starti_c = (iset-1)*maxco+1
            ncoi = npgf(iset)*ncoset(lmax(iset))
            sgfi = first_sgf(1,iset)

            DO jset = iset,nset
               startj_c = (jset-1)*maxco+1
               ncoj = npgf(jset)*ncoset(lmax(jset))
               sgfj = first_sgf(1,jset)

               ! hard integrals
               CALL ab_contract(int_h(sgfi:sgfi+nsgf_set(iset)-1,sgfj:sgfj+nsgf_set(jset)-1), &
                                int_c_h(starti_c:starti_c+ncoi-1,startj_c:startj_c+ncoj-1), &
                                sphi_h(:,sgfi:), sphi_h(:,sgfj:), ncoi, ncoj, &
                                nsgf_set(iset), nsgf_set(jset))

               ! soft integrals
               CALL ab_contract(int_s(sgfi:sgfi+nsgf_set(iset)-1,sgfj:sgfj+nsgf_set(jset)-1), & 
                                int_c_s(starti_c:starti_c+ncoi-1,startj_c:startj_c+ncoj-1), &
                                sphi_s(:,sgfi:), sphi_s(:,sgfj:), ncoi, ncoj, & 
                                nsgf_set(iset), nsgf_set(jset))

               ! deal with symmetries
               int_h(sgfj:sgfj+nsgf_set(jset)-1,sgfi:sgfi+nsgf_set(iset)-1) = &
                     TRANSPOSE(int_h(sgfi:sgfi+nsgf_set(iset)-1,sgfj:sgfj+nsgf_set(jset)-1))

               int_s(sgfj:sgfj+nsgf_set(jset)-1,sgfi:sgfi+nsgf_set(iset)-1) = &
                     TRANSPOSE(int_s(sgfi:sgfi+nsgf_set(iset)-1,sgfj:sgfj+nsgf_set(jset)-1))

            END DO !jset
         END DO !iset

!     Put the integrals for this atom as blocks in the dbcsr matrices
!     TODO: be carful when parallelizing
         CALL dbcsr_put_block(int_mat_h%matrix, iat, iat, int_h)
         CALL dbcsr_put_block(int_mat_s%matrix, iat, iat, int_s)

!     Clean-up 
         DEALLOCATE(rho_h, rho_s, int_h, int_s, int_c_h, int_c_s)

      END DO !iat

!  Finalize the matrix
      CALL dbcsr_finalize(int_mat_h%matrix)
      CALL dbcsr_finalize(int_mat_s%matrix)


   END SUBROUTINE integrate_ri_gapw_atoms

! **************************************************************************************************
!> \brief Collocates the fitted smooth density on a pw grid
!> \param rho_smooth_r the pw_p_type where the real smooth density over all space is stored
!> \param rho_smooth_g the pw_p_type where the imaginary smooth density is stored
!> \param ri_gapw_env ...
!> \param qs_env ...
!> \note This subroutine is HEAVILY inspired by "calculate_lri_rho_elec" form D. Golze
!>       As of now, use the qs pw_env. The reciprocal space density is just needed as a tool
! **************************************************************************************************
   SUBROUTINE put_density_on_pw_grid(rho_smooth_r, rho_smooth_g, ri_gapw_env, qs_env)

      TYPE(pw_p_type)                                 :: rho_smooth_r, rho_smooth_g
      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "put_density_on_pw_grid", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: igrid_level, maxco, ikind, nkind, na1, &
                                                         lmax_global, my_pos, group_size, iat, &
                                                         nseta, atom_a, m1, ipgf, dir, sgfa, ncoa, &
                                                         iset, offset
      REAL(dp)                                        :: eps_rho_rspace
      LOGICAL                                         :: map_consistent
      TYPE(pw_env_type), POINTER                      :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER          :: mgrid_gspace, mgrid_rspace
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER     :: pw_pools
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(realspace_grid_p_type), DIMENSION(:), &
                                             POINTER  :: rs_rho
      TYPE(realspace_grid_type), POINTER              :: rs_grid
      TYPE(cell_type), POINTER                        :: cell
      TYPE(gridlevel_info_type), POINTER              :: gridlevel_info
      TYPE(cube_info_type), DIMENSION(:), POINTER     :: cube_info
      TYPE(dft_control_type), POINTER                 :: dft_control
      TYPE(gto_basis_set_type), POINTER               :: soft_basis
      INTEGER, DIMENSION(:), POINTER                  :: la_max, atom_list, la_min, npgfa, &
                                                         nsgfa
      REAL(dp), DIMENSION(:,:), POINTER               :: pab, zeta, sphi_a, work
      TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
      REAL(dp), DIMENSION(3)                          :: ra  
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      REAL(dp), DIMENSION(:), POINTER                 :: coeffs
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgfa
      LOGICAL, DIMENSION(:), ALLOCATABLE              :: map_it
      INTEGER, DIMENSION(3)                           :: lb, location, tp, ub 

      NULLIFY(pw_env, mgrid_gspace, mgrid_rspace, pw_pools, qs_kind_set, rs_rho, rs_grid, cell)
      NULLIFY(gridlevel_info, cube_info, dft_control, soft_basis, pab, atomic_kind_set, coeffs)
      NULLIFY(particle_set, la_min, zeta, npgfa, sphi_a, first_sgfa, nsgfa, work)

!  Initialisation
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, cell=cell, &
                      dft_control=dft_control, nkind=nkind, atomic_kind_set=atomic_kind_set)
      pw_env => ri_gapw_env%pw_env

      cube_info => pw_env%cube_info 
      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      map_consistent = dft_control%qs_control%map_consistent
      gridlevel_info => pw_env%gridlevel_info

      ! *** set up the pw multi-grids *** !
      CALL pw_env_get(pw_env=pw_env, rs_grids=rs_rho, pw_pools=pw_pools)

      CALL pw_pools_create_pws(pw_pools, mgrid_rspace, &
                               use_data=REALDATA3D, &
                               in_space=REALSPACE)

      CALL pw_pools_create_pws(pw_pools, mgrid_gspace, &
                               use_data=COMPLEXDATA1D, &
                               in_space=RECIPROCALSPACE)

      ! *** set up the rs multi-grids *** !
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid)
         CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
      END DO

      !extract general basis info
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=maxco, basis_type="RI_XAS_SOFT")
      lmax_global = 0
      DO ikind = 1,nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=soft_basis, basis_type="RI_XAS_SOFT")
         IF (.NOT. ASSOCIATED(soft_basis)) CYCLE
         CALL get_gto_basis_set(soft_basis, lmax=la_max)
         lmax_global = MAX(MAXVAL(la_max), lmax_global)
      END DO

      !Get ready before we loop
      ALLOCATE(pab(maxco,1))
      offset = 0
      my_pos = mgrid_rspace(1)%pw%pw_grid%para%my_pos
      group_size = mgrid_rspace(1)%pw%pw_grid%para%group_size

      !Loop on kinds and atoms of kind. Only work on atoms of excited kinds
      DO ikind = 1,nkind

         IF (.NOT. ANY(ri_gapw_env%excited_kinds == ikind)) CYCLE

         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=soft_basis, basis_type="RI_XAS_SOFT")

         !get basis info for this kind
         CALL get_gto_basis_set(soft_basis, lmax=la_max, lmin=la_min, zet=zeta, nset=nseta, &
                                npgf=npgfa, sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)

         DO iat = 1,SIZE(atom_list)
            atom_a = atom_list(iat)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            coeffs => ri_gapw_env%coeff_s(atom_a)%array

            !deal with the grid distribution
            m1 = MAXVAL(npgfa(1:nseta))
            ALLOCATE(map_it(m1))
            DO iset = 1, nseta
               ! collocate this set locally? 
               map_it = .FALSE.
               DO ipgf = 1, npgfa(iset)
                  igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))
                  rs_grid => rs_rho(igrid_level)%rs_grid
                  IF (.NOT. ALL(rs_grid%desc%perd == 1)) THEN
                     DO dir = 1, 3
                        ! bounds of local grid (i.e. removing the 'wings'), if periodic
                        tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_grid%desc%npts(dir))
                        tp(dir) = MODULO(tp(dir), rs_grid%desc%npts(dir))
                        IF (rs_grid%desc%perd(dir) .NE. 1) THEN
                           lb(dir) = rs_grid%lb_local(dir)+rs_grid%desc%border
                           ub(dir) = rs_grid%ub_local(dir)-rs_grid%desc%border
                        ELSE
                           lb(dir) = rs_grid%lb_local(dir)
                           ub(dir) = rs_grid%ub_local(dir)
                        ENDIF
                        ! distributed grid, only map if it is local to the grid
                        location(dir) = tp(dir)+rs_grid%desc%lb(dir)
                     ENDDO
                     IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                         lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                         lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                        map_it(ipgf) = .TRUE.
                     ENDIF
                  ELSE
                     ! not distributed, just a round-robin distribution over the full set of CPUs
                     IF (MODULO(offset, group_size) == my_pos) map_it(ipgf) = .TRUE.
                  ENDIF
               END DO 
               offset = offset+1

               !Actually evaluating the density on the grid
               IF (ANY(map_it(1:npgfa(iset)))) THEN
                  sgfa = first_sgfa(1, iset)
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  m1 = sgfa+nsgfa(iset)-1
                  ALLOCATE(work(nsgfa(iset), 1))
                  work(1:nsgfa(iset), 1) = coeffs(sgfa:m1)
                  pab = 0.0_dp

                  CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), 1.0_dp, sphi_a(1, sgfa), &
                             SIZE(sphi_a, 1), work(1, 1), SIZE(work, 1), &
                             0.0_dp, pab(1, 1), SIZE(pab, 1))

                  DO ipgf = 1, npgfa(iset)
                     na1 = (ipgf-1)*ncoset(la_max(iset))
                     igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))
                     rs_grid => rs_rho(igrid_level)%rs_grid
                     IF (map_it(ipgf)) THEN

                        !take product of basis function with 1
                        CALL collocate_pgf_product_rspace(la_max=la_max(iset), &
                                                          zeta=zeta(ipgf, iset), &
                                                          la_min=la_min(iset), &
                                                          lb_max=0, zetb=0.0_dp, lb_min=0, &
                                                          ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                                          rab2=0.0_dp, scale=1._dp, &
                                                          pab=pab, o1=na1, o2=0, &
                                                          rsgrid=rs_grid, cell=cell, &
                                                          cube_info=cube_info(igrid_level), &
                                                          eps_rho_rspace=eps_rho_rspace, &
                                                          ga_gb_function=FUNC_AB, &
                                                          map_consistent=map_consistent, &
                                                          lmax_global=lmax_global)

                     END IF ! map_it(ipgf)
                  END DO !ipgf
                  DEALLOCATE(work)
               END IF !ANY(map_it)
            END DO !iset
            DEALLOCATE(map_it)
         END DO !iat
      END DO !ikind
      DEALLOCATE(pab)

      !Finalizing the density rho_smooth_r
      CALL pw_zero(rho_smooth_r%pw)
      CALL pw_zero(rho_smooth_g%pw)

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL pw_zero(mgrid_rspace(igrid_level)%pw)
         CALL rs_pw_transfer(rs=rs_rho(igrid_level)%rs_grid, & 
                             pw=mgrid_rspace(igrid_level)%pw, dir=rs2pw)
         CALL rs_grid_release(rs_rho(igrid_level)%rs_grid)
      END DO

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL pw_zero(mgrid_gspace(igrid_level)%pw)
         CALL pw_transfer(mgrid_rspace(igrid_level)%pw, &
                          mgrid_gspace(igrid_level)%pw)
         CALL pw_axpy(mgrid_gspace(igrid_level)%pw, rho_smooth_g%pw)
      END DO
      CALL pw_transfer(rho_smooth_g%pw, rho_smooth_r%pw)

      ! *** give back the multi-grids *** !
      CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace)
      CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace)

   END SUBROUTINE put_density_on_pw_grid

! **************************************************************************************************
!> \brief Integrate a function of r on the PW grid wrt the RI_XAS_SOFT basis
!> \param int_mat the matrix containing the integrals
!> \param singlet whther we consider singlet or triplet excitations
!> \param ri_gapw_env ...
!> \param qs_env ...
!> \note Only the PADE functional is available 
! **************************************************************************************************
   SUBROUTINE integrate_ri_gapw_smooth(int_mat, singlet, ri_gapw_env, qs_env)

      TYPE(dbcsr_p_type)                              :: int_mat
      LOGICAL                                         :: singlet
      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "integrate_ri_gapw_smooth", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: n
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      TYPE(pw_p_type)                                 :: rho_smooth_r, rho_smooth_g, fxc_r
      TYPE(pw_env_type), POINTER                      :: pw_env
      TYPE(pw_pool_type), POINTER                     :: auxbas_pw_pool
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: rs
      TYPE(dft_control_type), POINTER                 :: dft_control

      NULLIFY(dbcsr_dist, pw_env, auxbas_pw_pool, dft_control)

!  Create a dbcsr_matrix with the same distribution as the normal ones, but with the size of the 
!  orbital basis for non interesting atoms and the size of the RI basis for the ones we like
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, dft_control=dft_control)
      CALL dbcsr_create(matrix=int_mat%matrix, name="RI GAPW SMOOTH INT", matrix_type='S', &
                        dist=dbcsr_dist, row_blk_size=ri_gapw_env%blk_size, &
                        col_blk_size=ri_gapw_env%blk_size)

!  Create the smooth density pw_p_types here
      pw_env => ri_gapw_env%pw_env
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_smooth_r%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, fxc_r%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_smooth_g%pw, use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

!  Put the density on the pw grid
      CALL put_density_on_pw_grid(rho_smooth_r, rho_smooth_g, ri_gapw_env, qs_env)

!  Evaluate the functional on the grid
      CALL set_util(dft_control%qs_control%eps_rho_rspace)
      n = PRODUCT(SHAPE(rho_smooth_r%pw%cr3d))
      ALLOCATE(rs(n))
      CALL calc_rs_pw(rho_smooth_r%pw%cr3d, rs, n)
      CALL pw_zero(fxc_r%pw)

      IF (singlet) THEN
         CALL pade_kernel_singlet(n, rho_smooth_r%pw%cr3d, rs, fxc_r%pw%cr3d, &
                                  dft_control%qs_control%eps_rho_rspace)

      ELSE 
         CALL pade_kernel_triplet(n, rho_smooth_r%pw%cr3d, rs, fxc_r%pw%cr3d, &                         
                                  dft_control%qs_control%eps_rho_rspace)
      END IF


!  Before integration, need to multiply integrand by the volume element...
      fxc_r%pw%cr3d = fxc_r%pw%cr3d * fxc_r%pw%pw_grid%dvol

!  Integrate 
      CALL integrate_v_rspace(fxc_r, int_mat, qs_env=qs_env, calculate_forces=.FALSE., & 
                              basis_type="RI_XAS_SOFT", task_list_external=ri_gapw_env%task_list, &
                              pw_env_external=ri_gapw_env%pw_env)

!  Clean-up
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_smooth_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_smooth_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, fxc_r%pw)

   END SUBROUTINE integrate_ri_gapw_smooth

! **************************************************************************************************
!> \brief The actual pade kernel evaluation in the singlet case. Inspired from xc_pade.F
!> \param n the number of grid points
!> \param rho the array of the values of the density on the grid
!> \param rs the array of the values of rs on the grid
!> \param fxc the array of the values of fxc on the grid
!> \param eps_rho the density cutoff
!> \note TODO: OMP parallelize as in xc_pade.F
!>             In the singlet case, take d2exc_d2na + d2exc_dnadnb, can be shown to be the formula
!>             below, only dependent on n (spin-unpolarized)
! **************************************************************************************************
   SUBROUTINE pade_kernel_singlet(n, rho, rs, fxc, eps_rho)

      REAL(dp), DIMENSION(*), INTENT(IN)              :: rho, rs
      REAL(dp), DIMENSION(*), INTENT(INOUT)           :: fxc
      REAL(dp), INTENT(IN)                            :: eps_rho
      INTEGER                                         :: n

      CHARACTER(len=*), PARAMETER :: routineN = "pade_kernel_singlet", routineP=moduleN//":"//routineN 

      REAL(dp), PARAMETER  :: f13 = 1.0_dp/3.0_dp, &
                              f49 = 4.0_dp/9.0_dp

      REAL(dp), PARAMETER  :: a0 = 0.4581652932831429E+0_dp, &                                      
                              a1 = 0.2217058676663745E+1_dp, &
                              a2 = 0.7405551735357053E+0_dp, & 
                              a3 = 0.1968227878617998E-1_dp, &
                              b1 = 1.0000000000000000E+0_dp, &
                              b2 = 0.4504130959426697E+1_dp, &
                              b3 = 0.1110667363742916E+1_dp, & 
                              b4 = 0.2359291751427506E-1_dp 

      REAL(dp)                                        :: drs, d2rs, p, dp, d2p, q, dq, d2q, df, d2f
      INTEGER                                         :: ip

!     the functional has the form:               e = -f(rs) = -p(rs)/q(rs)
!     the first derivative wrt to rho (n):       de_dn = -df_drs*drs_dn
!     the seconde derivative wrt rho :           d2e_d2n = -d2f_d2rs*drs_dn**2 - df_drs*d2rs_d2n

!     rs has the form:                           rs = (3/4*pi)**1/3 * n**(-1/3)
!     the first derivative wrt n:                drs = -1/3 * rs * n**(-1)
!     the second derivative wrt n:               d2rs = 4/9 * rs * n**(-2)

!     We need p,q and the derivatives:           p = a0 + a1*rs + a2*rs**2 + a3*rs**3  
!                                                dp = a1 + 2*a2*rs + 3*a3*rs**2                     
!                                                d2p = 2*a2 + 6*a3*rs

!                                                q = b1*rs + b2*rs**2 + b3*rs**3 + b4*rs**4
!                                                dq = b1 + 2*b2*rs + 3*b3*rs**2 + 4*b4*rs**3
!                                                d2q = 2*b2 + 6*b3*rs + 12*b4*rs**2

!     Finally the derivatives of f wrt rs:       df = (dp*q - p*dq)/q**2
!                                                d2f = -2*dp*dq/q**2 + d2p/q 
!                                                      + p*( 2*dq**2/q**3 - d2q/q**2)

!     TODO: replace integer factors by double
      DO ip = 1,n
         IF (rho(ip) > eps_rho) THEN

            p = a0+(a1+(a2+a3*rs(ip))*rs(ip))*rs(ip)
            q = (b1+(b2+(b3+b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)  

            dp = a1 + 2*a2*rs(ip) + 3*a3*rs(ip)**2
            dq = b1 + 2*b2*rs(ip) + 3*b3*rs(ip)**2 + 4*b4*rs(ip)**3

            d2p = 2*a2 + 6*a3*rs(ip)
            d2q = 2*b2 + 6*b3*rs(ip) + 12*b4*rs(ip)**2

            df = (dp*q - p*dq)/q**2
            d2f = -2*dp*dq/q**2 + d2p/q + p*( 2*dq**2/q**3 - d2q/q**2)

            drs = -f13*rs(ip)/rho(ip)
            d2rs = f49*rs(ip)/rho(ip)**2

            !the factor 2 comes from the sum of the 2 derivatives
            fxc(ip) = -2*(d2f*drs**2 - df*d2rs)

         END IF
      END DO

   END SUBROUTINE pade_kernel_singlet

! **************************************************************************************************
!> \brief The actual pade kernel evaluation in the triplet case. Inspired from xc_pade.F
!> \param n the number of grid points
!> \param rho the array of the values of the density on the grid                                           
!> \param rs the array of the values of rs on the grid                                                     
!> \param fxc the array of the values of fxc on the grid                                                   
!> \param eps_rho the density cutoff                                                                       
!> \note TODO: OMP parallelize as in xc_pade.F                                                       
!>             In the triplet case, take d2exc_d2na - d2exc_dnadnb, can be shown to be the formula   
!>             below, only dependent on n (spin-unpolarized) 
! **************************************************************************************************
   SUBROUTINE pade_kernel_triplet(n, rho, rs, fxc, eps_rho)
      
      REAL(dp), DIMENSION(*), INTENT(IN)              :: rho, rs  
      REAL(dp), DIMENSION(*), INTENT(INOUT)           :: fxc  
      REAL(dp), INTENT(IN)                            :: eps_rho
      INTEGER                                         :: n

      CHARACTER(len=*), PARAMETER :: routineN = "pade_kernel_triplet", routineP=moduleN//":"//routineN

      REAL(dp), PARAMETER  :: a0 = 0.4581652932831429E+0_dp, &
                              a1 = 0.2217058676663745E+1_dp, & 
                              a2 = 0.7405551735357053E+0_dp, &
                              a3 = 0.1968227878617998E-1_dp, &
                              b1 = 1.0000000000000000E+0_dp, &
                              b2 = 0.4504130959426697E+1_dp, &
                              b3 = 0.1110667363742916E+1_dp, &
                              b4 = 0.2359291751427506E-1_dp

      REAL(dp), PARAMETER  :: da0 = 0.119086804055547E+0_dp, &
                              da1 = 0.6157402568883345E+0_dp, &
                              da2 = 0.1574201515892867E+0_dp, & 
                              da3 = 0.3532336663397157E-2_dp, &
                              db1 = 0.0000000000000000E+0_dp, &
                              db2 = 0.2673612973836267E+0_dp, &
                              db3 = 0.2052004607777787E+0_dp, &
                              db4 = 0.4200005045691381E-2_dp 

      REAL(dp), PARAMETER  :: d2f_d2xi = 4.0_dp/9.0_dp * 1.0_dp/(2.0_dp**(1.0_dp/3.0_dp)-1)

      REAL(dp)                                        :: p, dp, q, dq, de_df
      INTEGER                                         :: ip

!     We can express the functional as a function of fx =>        e = e(f) = -p(f)/q(f)
!     In spin polarized cases, one can show that we only need:    de_df
!                                                         and:    d2f_d2xi

!     The first derivitive of the functional has the form:        de_df = - (dp*q - p*dq)/q^2

!     We need p,q and their derivative wrt to f (at xi=0):        p = a0 + a1*rs + a2*rs^2 +a3*rs^3
!                                                                 dp = da0 +da1*rs +da2*rs^2 +da3*rs^3 

!                                                                 q = b1*rs + b2*rs^2 + b3*rs^3 +b4*rs^4
!                                                                 dq = db1*rs +db2*rs^2 +db3*rs^3 +db4*rs^4

!     The double derivative of f wrt xi (at xi=0):                d2f_d2xi = 4/9 * 1/(2^1/3 - 1)

!     In the end, we compute:                                     2/n^2 * de_df * d2f_d2xi

!  TODO: replace the integer factors by doubles
      DO ip = 1,n
         IF (rho(ip) > eps_rho) THEN 

            p = a0+(a1+(a2+a3*rs(ip))*rs(ip))*rs(ip)
            q = (b1+(b2+(b3+b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)   

            dp = da0+(da1+(da2+da3*rs(ip))*rs(ip))*rs(ip)
            dq = (db1+(db2+(db3+db4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

            de_df = -(dp*q - p*dq)/q**2

            fxc(ip) = 2/(rho(ip))**2 * de_df *d2f_d2xi

         END IF
      END DO

   END SUBROUTINE pade_kernel_triplet

END MODULE xas_tdp_gapw
