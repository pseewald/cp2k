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
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: gto_basis_set_type, gto_basis_set_p_type, &
                                              allocate_gto_basis_set, get_gto_basis_set
   USE cp_control_types,                ONLY: qs_control_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type, dbcsr_distribution_type,&
                                              dbcsr_create, dbcsr_print
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids, get_number_of_lebedev_grid,&
                                              init_lebedev_grids, lebedev_grid
   USE memory_utilities,                ONLY: reallocate
   USE orbital_pointers,                ONLY: nsoset, indso
   USE particle_types,                  ONLY: particle_type
   USE pw_types,                        ONLY: pw_p_type, pw_create, REALSPACE, REALDATA3D
   USE qs_environment_types,            ONLY: qs_environment_type, get_qs_env
   USE qs_grid_atom,                    ONLY: grid_atom_type, allocate_grid_atom, &
                                              deallocate_grid_atom, create_grid_atom
   USE qs_harmonics_atom,               ONLY: create_harmonics_atom, get_maxl_CG, &
                                              harmonics_atom_type, allocate_harmonics_atom, &
                                              deallocate_harmonics_atom
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: qs_kind_type, get_qs_kind
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type, &
                                              deallocate_neighbor_list_set
   USE qs_rho_types,                    ONLY: qs_rho_type, qs_rho_get
   USE qs_rho_atom_types,               ONLY: rho_atom_type, allocate_rho_atom_set, &
                                              deallocate_rho_atom_set
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
!> \param eps_rho the density cutoff in rreal space
!> \param rho_atom the type containing the grid integral info for all excited atoms 
!> \param grid_atom_set the set of atomic grid for each kind
!> \param task_list the task list associated with the RI_XAS_SOFT basis for the pw integration
!> \note the rho_atom_set has the size of the number of atoms, even thogh some are not excited
!>       same for grid_atom_set, it has the same size as the number of kind, but some are empty
!>       same for harmonics_atom_set
! **************************************************************************************************
   TYPE ri_gapw_env_type
      REAL(dp)                                        :: eps_rho
      INTEGER, DIMENSION(:), POINTER                  :: excited_kinds
      TYPE(rho_atom_type), DIMENSION(:) , POINTER     :: rho_atom_set
      TYPE(grid_atom_p_type), DIMENSION(:), POINTER   :: grid_atom_set
      TYPE(harmonics_atom_p_type), DIMENSION(:), &
                                            POINTER   :: harmonics_atom_set
      TYPE(task_list_type), POINTER                   :: task_list
   END TYPE ri_gapw_env_type

   PUBLIC :: ri_gapw_env_type

   PUBLIC :: init_ri_gapw_env, release_ri_gapw_env, integrate_ri_gapw_smooth, pade_kernel_pw_create

CONTAINS

! **************************************************************************************************
!> \brief Initializes a ri_gapw_env type given the qs_env and the qs_control
!> \param ri_gapw_env the ri_gapw_env to initialize
!> \param excited_kinds the indices of the excited kinds
!> \param qs_env ...
!> \param qs_control ...
! **************************************************************************************************
   SUBROUTINE init_ri_gapw_env(ri_gapw_env, excited_kinds, qs_env, qs_control)

      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(qs_control_type), POINTER                  :: qs_control
      INTEGER, DIMENSION(:), POINTER                  :: excited_kinds

      CHARACTER(len=*), PARAMETER :: routineN = 'init_ri_gapw_env', routineP = moduleN//":"//routineN

      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER               :: soft_basis, hard_basis
      REAL(dp)                                        :: rc
      LOGICAL                                         :: paw_atom
      INTEGER                                         :: ikind, nkind, ilist, natom, nex_kinds, &
                                                         handle
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_a 
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &                                                
                                             POINTER  :: ab_list
      TYPE(qs_ks_env_type), POINTER                   :: ks_env

      NULLIFY(qs_kind_set, soft_basis, hard_basis, basis_set_a, ab_list, ks_env)

      CALL timeset(routineN, handle)

!  Initializing the type
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, ks_env=ks_env, natom=natom)
      nkind = SIZE(qs_kind_set)   
      nex_kinds = SIZE(excited_kinds)

      ALLOCATE(ri_gapw_env)
      ri_gapw_env%eps_rho = qs_control%eps_rho_rspace
      CALL allocate_rho_atom_set(ri_gapw_env%rho_atom_set, natom)
      CALL allocate_task_list(ri_gapw_env%task_list)
      ALLOCATE(ri_gapw_env%excited_kinds(nex_kinds))
      ri_gapw_env%excited_kinds = excited_kinds
      ALLOCATE(ri_gapw_env%grid_atom_set(nkind))
      ALLOCATE(ri_gapw_env%harmonics_atom_set(nkind))

!  Allocate and initialize the atomic grids and harmonics
      CALL init_ri_gapw_grid_harmo(6, ri_gapw_env, qs_env, qs_control)

!  Go over all qs_kind and make sure the RI_XAS_SOFT basis is computed and stored
      DO ikind = 1,nkind
         IF (.NOT. ANY(excited_kinds == ikind)) CYCLE
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=soft_basis, basis_type="RI_XAS_SOFT")
         IF (.NOT. ASSOCIATED(soft_basis)) THEN
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=hard_basis, basis_type="RI_XAS", & 
                             hard_radius=rc) 
            CALL allocate_gto_basis_set(soft_basis)
            CALL create_soft_basis(orb_basis=hard_basis, soft_basis=soft_basis, &
                                  eps_fit=qs_control%gapw_control%eps_fit, rc=rc, &
                                  paw_atom=paw_atom, paw_type_forced=.TRUE., &
                                  gpw_type_forced=.FALSE.) 
            CALL init_interaction_radii_orb_basis(soft_basis, qs_control%eps_pgf_orb)
            CALL add_basis_set_to_container(qs_kind_set(ikind)%basis_sets, &
                                            soft_basis, "RI_XAS_SOFT")
         END IF
                                    
         NULLIFY(soft_basis)
      END DO

!  Create a neighbor list based on the RI_XAS_SOFT basis with the minimal image convention since
!  we only consider atom idependentely, we want to avoid the overlap with its image
! TODO: when smooth integral available, try making the neighbour list such that each atom is only
!       neighbor with itself and check if the block diagonal of the integral matrix changes
      ALLOCATE(basis_set_a(nkind))
      CALL basis_set_list_setup(basis_set_a, "RI_XAS_SOFT", qs_kind_set)
      CALL setup_neighbor_list(ab_list, basis_set_a, qs_env=qs_env, mic=.TRUE.) 

!  Create the task_list based on the neighbor list above. Keep input param
      CALL generate_qs_task_list(ks_env, ri_gapw_env%task_list, reorder_rs_grid_ranks=.TRUE., & 
                                 soft_valid=.FALSE., basis_type="RI_XAS_SOFT", &
                                 skip_load_balance_distributed=qs_control%skip_load_balance_distributed, &
                                 sab_orb_external=ab_list)
!  Clean up
      DO ilist = 1,SIZE(ab_list)
         CALL deallocate_neighbor_list_set(ab_list(ilist)%neighbor_list_set)
      END DO
      DEALLOCATE(ab_list)
      DEALLOCATE(basis_set_a)

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

      DEALLOCATE(ri_gapw_env)

   END SUBROUTINE release_ri_gapw_env

! **************************************************************************************************
!> \brief Initializes the atomic grids and harmonics for the RI gapw calculation
!> \param lmax_sphere maximum L for the Clebsch-Gordan coeffs
!> \param ri_gapw_env the RI gapw environement 
!> \param qs_env ...
!> \param qs_control ...
!> \note Largely inspired by init_rho_atom subroutine (which is sadly restricted to orbital basis set)
! **************************************************************************************************
   SUBROUTINE init_ri_gapw_grid_harmo(lmax_sphere, ri_gapw_env, qs_env, qs_control)

      INTEGER                                         :: lmax_sphere
      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(qs_control_type), POINTER                  :: qs_control

      CHARACTER(len=*), PARAMETER :: routineN = "init_ri_gapw_grid_harmo", &
                                     routineP = routineN//":"//moduleN

      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rga
      INTEGER                                            :: llmax, maxlgto, max_s_harm, max_s_set, &
                                                            lcleb, lc1, iso1, l1, m1, lc2, iso2, &
                                                            l2, m2, l1l2, mp, mm, lp, il, iso, &
                                                            ikind, quadrature, i, ll, na, la, nr, &
                                                            maxs, maxl
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: tmp_basis
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      NULLIFY(my_CG, qs_kind_set, tmp_basis, grid_atom, harmonics)

!  Initialization of some integer for the CG coeff generation
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      maxlgto = 0
      DO i = 1,SIZE(ri_gapw_env%excited_kinds)
         ikind = ri_gapw_env%excited_kinds(i)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="RI_XAS")
         CALL get_gto_basis_set(tmp_basis, maxl=il)
         maxlgto = MAX(maxlgto, il)
      END DO

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

         IF (llmax > la) CPABORT("Increase LMAXSPHERE for XAS_TDP GAPW")

         CALL create_grid_atom(grid_atom, nr, na, llmax, ll, quadrature)

         CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="RI_XAS")
         CALL get_gto_basis_set(gto_basis_set=tmp_basis, maxl=maxl)
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
!> \brief Integrate a function of r on the PW grid wrt the RI_XAS_SOFT basis
!> \param f_rspace the function of r to integrate
!> \param int_mat the matrix containing the integrals
!> \param ri_gapw_env the ri_gapw_env of the atom
!> \param qs_env ...
!> \note  Takes the default pw_env, (TODO) might want to make our own
!>        The matrix int_mat gets it's shape in this routine
! **************************************************************************************************
   SUBROUTINE integrate_ri_gapw_smooth(f_rspace, int_mat, ri_gapw_env, qs_env)

      TYPE(pw_p_type)                                 :: f_rspace
      TYPE(dbcsr_p_type)                              :: int_mat
      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "integrate_ri_gapw_smooth", &
                                     routineP = moduleN//":"//routineN

      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      INTEGER                                         :: nsgf, natom, iat, ikind
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set

      NULLIFY(dbcsr_dist, blk_size, particle_set, qs_kind_set)

!  Create a dbcsr_matrix with the same distribution as the normal ones, but with 1x1 block size for
!  non interesting atoms and the size of the RI basis for the ones we like
      CPASSERT(ASSOCIATED(int_mat%matrix))
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, natom=natom, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)
      ALLOCATE(blk_size(natom))
      DO iat = 1,natom
         ikind = particle_set(iat)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kind_set(ikind), basis_type="RI_XAS_SOFT", nsgf=nsgf)
         IF (ANY(ri_gapw_env%excited_kinds == ikind)) THEN
            blk_size(iat) = nsgf
         ELSE
            blk_size(iat) = 1
         END IF
      END DO
      CALL dbcsr_create(matrix=int_mat%matrix, name="RI GAPW SOFT INT", matrix_type='S', &
                        dist=dbcsr_dist, row_blk_size=blk_size, col_blk_size=blk_size)
      DEALLOCATE(blk_size)

!  Integrate 
      CALL integrate_v_rspace(f_rspace, int_mat, qs_env=qs_env, calculate_forces=.FALSE., & 
                              basis_type="RI_XAS_SOFT", task_list_external=ri_gapw_env%task_list)

   END SUBROUTINE integrate_ri_gapw_smooth

! **************************************************************************************************
!> \brief Computes the PADE fxc kernel on the real space grid for integration
!> \fxc_rspace the kernel on the real space grid
!> \rho_r the density on the real space grid
!> \note We assume that fxc will always be a function of r
! **************************************************************************************************
   SUBROUTINE pade_kernel_pw_create(fxc_rspace, rho_r, eps_rho)

      TYPE(pw_p_type)                                 :: fxc_rspace
      TYPE(pw_p_type)                                 :: rho_r
      REAL(dp)                                        :: eps_rho

      CHARACTER(len=*), PARAMETER :: routineN = "pade_kernel_pw_create", &
                                     routineP=routineN//":"//moduleN

      INTEGER                                         :: n
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: rs

!  Create the pw_p_type of fxc based on the grid of rho_r
      CALL pw_create(fxc_rspace%pw, rho_r%pw%pw_grid, use_data=REALDATA3D, in_space=REALSPACE)
      fxc_rspace%pw%cr3d = 0.0_dp

!  Get rs on the grid
      CALL set_util(eps_rho)
      n = PRODUCT(SHAPE(rho_r%pw%cr3d))
      ALLOCATE(rs(n))
      CALL calc_rs_pw(rho_r%pw%cr3d, rs, n)

!  It seems that all the distribution is taken care of by itself, so just work on the local data
!  We assume unpolarized spins PADE LDA, we do not multiply by rho (as in xc_pade)
      CALL pade_kernel_singlet(n, rho_r%pw%cr3d, rs, fxc_rspace%pw%cr3d, eps_rho)

   END SUBROUTINE pade_kernel_pw_create

! **************************************************************************************************
!> \brief The actual pade kernel evaluation in the singlet case. Inspired from xc_pade.F
!> \n the number of grid points
!> \rho the array of the values of the density on the grid
!> \rs the array of the values of rs on the grid
!> \fxc the array of the values of fxc on the grid
!> \eps_rho the density cutoff
!> \note TODO: OMP parallelize as in xc_pade.F
!>             In the singlet case, take d2exc_d2na + d2exc_dnadnb, can be shown to be the formula
!>             below, only dependent on n (spin-unpolarized)
! **************************************************************************************************
   SUBROUTINE pade_kernel_singlet(n, rho, rs, fxc, eps_rho)

      INTEGER, INTENT(IN)                             :: n
      REAL(dp), DIMENSION(*), INTENT(IN)              :: rho, rs
      REAL(dp), DIMENSION(*), INTENT(INOUT)           :: fxc
      REAL(dp), INTENT(IN)                            :: eps_rho

      CHARACTER(len=*), PARAMETER :: routineN = "pade_kernel_singlet", routineP=routineN//":"//moduleN 

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
!> \n the number of grid points                                                                      
!> \rho the array of the values of the density on the grid                                           
!> \rs the array of the values of rs on the grid                                                     
!> \fxc the array of the values of fxc on the grid                                                   
!> \eps_rho the density cutoff                                                                       
!> \note TODO: OMP parallelize as in xc_pade.F                                                       
!>             In the triplet case, take d2exc_d2na - d2exc_dnadnb, can be shown to be the formula   
!>             below, only dependent on n (spin-unpolarized) 
! **************************************************************************************************
   SUBROUTINE pade_kernel_triplet(n, rho, rs, fxc, eps_rho)
      
      INTEGER, INTENT(IN)                             :: n
      REAL(dp), DIMENSION(*), INTENT(IN)              :: rho, rs  
      REAL(dp), DIMENSION(*), INTENT(INOUT)           :: fxc  
      REAL(dp), INTENT(IN)                            :: eps_rho

      CHARACTER(len=*), PARAMETER :: routineN = "pade_kernel_triplet", routineP=routineN//":"//moduleN

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
