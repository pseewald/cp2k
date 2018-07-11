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
!  TODO: need to go over that and clean up some
   USE atomic_kind_types,               ONLY: atomic_kind_type, get_atomic_kind
   USE basis_set_container_types,       ONLY: add_basis_set_to_container, &
                                              remove_basis_from_container
   USE basis_set_types,                 ONLY: gto_basis_set_type, deallocate_gto_basis_set,&
                                              gto_basis_set_p_type, allocate_gto_basis_set
   USE dbcsr_api,                       ONLY: dbcsr_p_type, dbcsr_release, dbcsr_distribution_type,&
                                              dbcsr_distribution_get, dbcsr_create, dbcsr_print
   USE kinds,                           ONLY: dp
   USE cp_control_types,                ONLY: qs_control_type
   USE particle_types,                  ONLY: particle_type
   USE pw_types,                        ONLY: pw_p_type, pw_create, pw_release, REALSPACE, &
                                              REALDATA3D
   USE qs_environment_types,            ONLY: qs_environment_type, get_qs_env
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: qs_kind_type, get_qs_kind
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type, &
                                              deallocate_neighbor_list_set, &
                                              neighbor_list_iterate, & 
                                              neighbor_list_iterator_create, &                       
                                              neighbor_list_iterator_p_type, & 
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type, &
                                              get_iterator_info, &
                                              deallocate_neighbor_list_set
   USE qs_rho_types,                    ONLY: qs_rho_type, qs_rho_get
   USE qs_rho_atom_types,               ONLY: rho_atom_type, allocate_rho_atom_set, &
                                              deallocate_rho_atom_set
   USE soft_basis_set,                  ONLY: create_soft_basis
   USE task_list_types,                 ONLY: task_list_type, allocate_task_list,&
                                              deallocate_task_list
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE xc_functionals_utilities,        ONLY: set_util, calc_rs_pw

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "xas_tdp_gapw"

! **************************************************************************************************
!> \brief a environment type that contains all the info needed for the RI GAPW calculations 
!> \param eps_rho the density cutoff in rreal space
!> \param rho_atom the type containing the grid integral info for all excited atoms 
!> \param task_list the task list associated with the RI_XAS_SOFT basis for the pw integration
!> \note the rho_atom_set has the size of the number of atoms, even thogh some are not excited
! **************************************************************************************************
   TYPE ri_gapw_env_type
      REAL(dp)                                        :: eps_rho
      INTEGER, DIMENSION(:), POINTER                  :: excited_kinds
      TYPE(rho_atom_type), DIMENSION(:) , POINTER     :: rho_atom_set
      TYPE(task_list_type), POINTER                   :: task_list
   END TYPE ri_gapw_env_type

   PUBLIC :: ri_gapw_env_type

   PUBLIC :: init_ri_gapw_env, release_ri_gapw_env, integrate_ri_gapw_smooth

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
      INTEGER                                         :: ikind, nkind, ilist, natom
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_a 
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &                                                
                                             POINTER  :: ab_list
      TYPE(qs_ks_env_type), POINTER                   :: ks_env

      NULLIFY(qs_kind_set, soft_basis, hard_basis, basis_set_a, ab_list, ks_env)

!  Initializing the type
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, ks_env=ks_env, natom=natom)
      ALLOCATE(ri_gapw_env)
      ri_gapw_env%eps_rho = qs_control%eps_rho_rspace
      CALL allocate_rho_atom_set(ri_gapw_env%rho_atom_set, natom)
      CALL allocate_task_list(ri_gapw_env%task_list)
      ALLOCATE(ri_gapw_env%excited_kinds(SIZE(excited_kinds)))
      ri_gapw_env%excited_kinds = excited_kinds

!  Go over all qs_kind and make sure the RI_XAS_SOFT basis is computed and stored
      nkind = SIZE(qs_kind_set)   
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

   END SUBROUTINE init_ri_gapw_env

! **************************************************************************************************
!> \brief Releases the ri_gapw_env type
!> \param ri_gapw_env the type to release
! **************************************************************************************************
   SUBROUTINE release_ri_gapw_env(ri_gapw_env)

      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env              

      CHARACTER(len=*), PARAMETER :: routineN = "release_ri_gapw_env", &
                                     routineP = moduleN//":"//routineN

      IF (ASSOCIATED(ri_gapw_env%rho_atom_set)) THEN
         CALL deallocate_rho_atom_set(ri_gapw_env%rho_atom_set)
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
!> \brief Integrate a function of r on the PW grid wrt the RI_XAS_SOFT basis
!> \param f_rspace the function of r to integrate
!> \param int_mat the matrix containing the integrals
!> \param ri_gapw_env the ri_gapw_env of the atom
!> \param qs_env ...
!> \note  Takes the default pw_env, (TODO) might want to make our own
! **************************************************************************************************
   SUBROUTINE integrate_ri_gapw_smooth(ri_gapw_env, qs_env)

      TYPE(ri_gapw_env_type), POINTER                 :: ri_gapw_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "integrate_ri_gapw_smooth", &
                                     routineP = moduleN//":"//routineN

      TYPE(pw_p_type)                                 :: test_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER          :: rho_r
      TYPE(qs_rho_type), POINTER                      :: rho_struct
      TYPE(dbcsr_p_type)                              :: test_mat
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      INTEGER                                         :: nsgf, natom, iat, ikind
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set

      NULLIFY(dbcsr_dist, blk_size, test_mat%matrix, particle_set, qs_kind_set)
      NULLIFY(rho_r, rho_struct, test_rspace%pw)

!  Create a dbcsr_matrix with the same distribution as the normal ones, but with 1x1 block size for
!  non interesting atoms and the size of the RI basis for the ones we like
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
      ALLOCATE(test_mat%matrix)
      CALL dbcsr_create(matrix=test_mat%matrix, name="TESTMAT", matrix_type='S', dist=dbcsr_dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      DEALLOCATE(blk_size)

!  Generate a pw representation for the pade xc kernel
      CALL get_qs_env(qs_env, rho=rho_struct)
      CALL qs_rho_get(rho_struct, rho_r=rho_r)
      CALL pade_kernel_pw_create(test_rspace, rho_r(1), ri_gapw_env%eps_rho)

!  Integrate 
      CALL integrate_v_rspace(test_rspace, test_mat, qs_env=qs_env, calculate_forces=.FALSE., & 
                              basis_type="RI_XAS_SOFT", task_list_external=ri_gapw_env%task_list)

!  Clean-up
      CALL pw_release(test_rspace%pw)
      CALL dbcsr_release(test_mat%matrix)
      DEALLOCATE(test_mat%matrix)

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

!  Get rs on the grid
      CALL set_util(eps_rho)
      n = PRODUCT(SHAPE(rho_r%pw%cr3d))
      ALLOCATE(rs(n))
      CALL calc_rs_pw(rho_r%pw%cr3d, rs, n)

!  It seems that all the distribution is taken care of by itself, so just work on the local data
!  We assume unpolarized spins PADE LDA, we do not multiply by rho (as in xc_pade)
      CALL pade_kernel(n, rho_r%pw%cr3d, rs, fxc_rspace%pw%cr3d, eps_rho)

   END SUBROUTINE pade_kernel_pw_create

! **************************************************************************************************
!> \brief The actual pade kernel evaluation. Inspired from xc_pade.F
!> \n the number of grid points
!> \rho the array of the values of the density on the grid
!> \rs the array of the values of rs on the grid
!> \fxc the array of the values of fxc on the grid
!> \eps_rho the density cutoff
!> \note TODO: OMP parallelize as in xc_pade.F
! **************************************************************************************************
   SUBROUTINE pade_kernel(n, rho, rs, fxc, eps_rho)

      INTEGER, INTENT(IN)                             :: n
      REAL(dp), DIMENSION(*), INTENT(IN)              :: rho, rs
      REAL(dp), DIMENSION(*), INTENT(INOUT)           :: fxc
      REAL(dp)                                        :: eps_rho

      CHARACTER(len=*), PARAMETER :: routineN = "pade_kernel", routineP=routineN//":"//moduleN 

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

            fxc(ip) = -d2f*drs**2 - df*d2rs

         END IF
      END DO

   END SUBROUTINE

END MODULE xas_tdp_gapw
