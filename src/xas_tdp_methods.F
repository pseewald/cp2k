!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
 
! **************************************************************************************************
!> \brief Methods for X-Ray absorption spectroscopy (XAS) using TDDFPT  
!> \author AB (11.2017)
! **************************************************************************************************

MODULE xas_tdp_methods

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: allocate_sto_basis_set, create_gto_from_sto_basis, &
                                              deallocate_sto_basis_set, get_gto_basis_set, &
                                              gto_basis_set_type, set_sto_basis_set, srules,&
                                              sto_basis_set_type, deallocate_gto_basis_set,&
                                              init_orb_basis_set
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_type, cp_fm_get_submatrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE input_constants,                 ONLY: xas_tdp_by_index,&
                                              xas_tdp_by_kind,&
                                              xas_not_excited,&
                                              do_loc_none,&
                                              xas_1s_type,&
                                              xas_2s_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_get_subs_vals
   USE kinds,                           ONLY: dp
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: ptable
   USE qs_environment_types,            ONLY: qs_environment_type,&
                                              get_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind,&   
                                              qs_kind_type
   USE qs_loc_methods,                  ONLY: qs_loc_driver
   USE qs_loc_types,                    ONLY: localized_wfn_control_type,&
                                              qs_loc_env_create,&
                                              qs_loc_env_new_type
   USE qs_loc_utils,                    ONLY: qs_loc_control_init,&
                                              qs_loc_env_init,&
                                              set_loc_centers,&
                                              set_loc_wfn_lists
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE xas_methods,                     ONLY: calc_stogto_overlap
   USE xas_tdp_types,                   ONLY: xas_tdp_control_type,&
                                              read_xas_tdp_control,&
                                              xas_tdp_control_create,&
                                              xas_tdp_env_type,&
                                              xas_tdp_env_create,&
                                              xas_tdp_env_release,&
                                              set_xas_tdp_env,&
                                              donor_state_type,&
                                              set_donor_state
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_methods'

   PUBLIC :: xas_tdp

CONTAINS
                
! **************************************************************************************************
!> \brief Driver for XAS TDDFT calculations.
!> \param qs_env the inherited qs_environment
!> \param dft_control the set of input DFT keywords and values
!> \author AB
!> \note Empty for now...
! **************************************************************************************************
   SUBROUTINE xas_tdp(qs_env, dft_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, iat, istate, tmp_index,&
                                                            current_state_index, output_unit, itest
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      CHARACTER(len=2)                                   :: symbol_of_kind
      INTEGER, DIMENSION(:), POINTER                     :: atoms_of_kind
      TYPE(donor_state_type), POINTER                    :: current_state
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: dft_section, xas_tdp_section

      NULLIFY(xas_tdp_env, xas_tdp_control, atomic_kind_set, atoms_of_kind, current_state)

      CALL timeset(routineN, handle)

!  Logger initialization and XAS TDP banner printing
      NULLIFY (logger, dft_section, xas_tdp_section)                
      logger => cp_get_default_logger()
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")              
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")
      output_unit = cp_print_key_unit_nr(logger, xas_tdp_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".Log")                     

      IF (output_unit > 0) THEN
         WRITE(UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!                              XAS TDP                                      !", &
            "!    Starting TDDFPT deriven X-rays absorption spectroscopy calculations    !", &
            "!===========================================================================!"
      END IF

!  Safety checks
      IF (dft_control%lsd) THEN
         CPABORT("XAS TDP only implemented for closed shell ground state calculations")
      END IF
      
!  Initialization
      IF (output_unit > 0) THEN
         WRITE(UNIT=output_unit, FMT="(/,T3,A)") &
            "# Create and initialize the XAS_TDP environment"
      END IF 
      CALL xas_tdp_init(xas_tdp_env, xas_tdp_control, qs_env, dft_control)
      
!  Localization of the core orbitals if requested (used for better identification of donor states)
      IF (xas_tdp_env%qs_loc_env%do_localize) THEN
         IF (output_unit > 0) THEN
            WRITE(UNIT=output_unit, FMT="(/,T3,A,/)") &
               "# Localizing core orbitals for better identification"
         END IF
!        closed shell => myspin=1
         CALL qs_loc_driver(qs_env, xas_tdp_env%qs_loc_env, xas_tdp_control%print_loc_subsection,&
                           myspin=1)
      END IF
      
!  Assign lowest energy orbitals to excited atoms
      IF (xas_tdp_env%qs_loc_env%do_localize) THEN 
         CALL assign_loc_mos_to_ex_atoms(xas_tdp_env, xas_tdp_control, qs_env) 
      ELSE IF (.NOT. xas_tdp_env%qs_loc_env%do_localize) THEN
!        If no localization done, it is assumed that core orbitals are easily found with STO
!        overlap only. Then keep all possible links between MOs and donor states at this point
         xas_tdp_env%mos_of_ex_atoms(:,:) = 1
      END IF
!     for test purposes
      PRINT*, "Test for mos_of_ex_atoms:"
      DO iat = 1,xas_tdp_control%n_search
         PRINT*, xas_tdp_env%mos_of_ex_atoms(iat,:)
      END DO
      IF (output_unit > 0) THEN
         WRITE(UNIT=output_unit, FMT="(/,T3,A)") &
            "# Assign the ", xas_tdp_control%n_search,"  lowest energy MOs to excited atoms"
      END IF

!  Loop over donor states for calculation

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)      
      current_state_index = 1

!     Loop over atomic kinds
      DO ikind = 1,SIZE(atomic_kind_set)
         
         CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), element_symbol=symbol_of_kind,&
                              atom_list=atoms_of_kind)

!        Loop over atoms of kind
         DO iat = 1,SIZE(atoms_of_kind)
            
            IF (ANY(xas_tdp_env%ex_atom_indices .EQ. atoms_of_kind(iat))) THEN
               tmp_index = MINLOC(ABS(xas_tdp_env%ex_atom_indices-atoms_of_kind(iat)),1)  

!              Loop over states of excited atom of kind
               DO istate = 1,SIZE(xas_tdp_env%state_types,1)

                  IF (xas_tdp_env%state_types(istate, tmp_index) .NE. xas_not_excited) THEN
                    
                     current_state => xas_tdp_env%donor_states(current_state_index)
                     CALL set_donor_state(current_state,&
                                          at_index=xas_tdp_env%ex_atom_indices(tmp_index),&
                                          at_symbol=symbol_of_kind, kind_index=ikind,&
                                          state_type=xas_tdp_env%state_types(istate,tmp_index))

!                    Assign best fitting MO to current donnor core state
                     
                     PRINT*, "Assigning MO for donor state with properties: "
                     PRINT*, "Chemical symbol: ", current_state%at_symbol
                     PRINT*, "Kind index: ", current_state%kind_index 
                     PRINT*, "Atom index: ", current_state%at_index
                     PRINT*, "State type: ", current_state%state_type
   
                     CALL assign_mo_to_donor_state(current_state, xas_tdp_env,&
                                                   xas_tdp_control, qs_env)

                     PRINT*, "mos_of_ex_atoms after assignment"
                     DO itest = 1,xas_tdp_control%n_search
                        PRINT*, xas_tdp_env%mos_of_ex_atoms(itest,:) 
                     END DO

                     PRINT*, "Index of newly associated MO: ", current_state%mo_index
                     PRINT*, "Energy eval of MO: ", current_state%energy_eval

                     current_state_index = current_state_index + 1
                     NULLIFY(current_state)
                  END IF
               END DO
            END IF
         END DO
      END DO

!  Clean-up
      CALL xas_tdp_env_release(xas_tdp_env)
      IF (output_unit > 0) THEN
         WRITE(UNIT=output_unit, FMT="(/,T3,A,/,T3,A,/,T3,A,/)") &
            "!===========================================================================!", &
            "!     End of TDDFPT deriven X-rays absorption spectroscopy calculations     !", &
            "!===========================================================================!"   
      END IF

      CALL timestop(handle)

!  Structure would look like:
!     - Initialization of the xas_tdp_control
!     - Initialization of the xas_tdp_env (including allocation of the donor_state_type array)
!     - Loop over kinds, then atoms of kind, then states of atom of kind
!        - At this point initialize specific donor_state_types
!        - Do calculation for specific donor_state_type
!     - Collect and print results
!     - Release all types


   END SUBROUTINE xas_tdp
   
! **************************************************************************************************
!> \brief Overall control and  environment types initialization 
!> \param xas_tdp_env the environment type to initialize
!> \param xas_tdp_control the control type to initialize
!> \param qs_env the inherited qs environement type
!> \param dft_control the type in which xas_tdp_control lives
! **************************************************************************************************
   SUBROUTINE xas_tdp_init(xas_tdp_env, xas_tdp_control, qs_env, dft_control)
      
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type),    POINTER                 :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_init', routineP = moduleN//':'//routineN
      
      TYPE(section_vals_type), POINTER                   :: dft_section, xas_tdp_section,&
                                                            loc_section
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: at_kind_set
      INTEGER                                            :: nex_atoms, n_kinds, i, nat_of_kind,&
                                                            j, nex_kinds, k, n_donor_states,&
                                                            nspins, homo
      INTEGER, DIMENSION(2)                              :: n_mo, n_moloc
      INTEGER, DIMENSION(:), POINTER                     :: ex_atom_indices, ind_of_kind
      INTEGER, DIMENSION(:,:), POINTER                   :: state_types
      CHARACTER(len=2)                                   :: symbol
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY(dft_section, xas_tdp_section, at_kind_set, ex_atom_indices, state_types, ind_of_kind)
      NULLIFY(qs_loc_env, loc_section, mos)

!  XAS TDP control type initialization
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")                     
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")

      CALL xas_tdp_control_create(dft_control%xas_tdp_control)
      CALL read_xas_tdp_control(dft_control%xas_tdp_control,xas_tdp_section)   

      xas_tdp_control => dft_control%xas_tdp_control

!  XAS TDP environment type initialization
      CALL xas_tdp_env_create(xas_tdp_env)

!  Retrieving the excited atoms indices and correspondig state types
      IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

!        simply copy from xas_tdp_control
         nex_atoms = SIZE(xas_tdp_control%list_ex_atoms)
         CALL set_xas_tdp_env(xas_tdp_env, nex_atoms=nex_atoms)  
         ALLOCATE(xas_tdp_env%ex_atom_indices(nex_atoms))
         ALLOCATE(xas_tdp_env%state_types(SIZE(xas_tdp_control%state_types,1),nex_atoms))
         xas_tdp_env%ex_atom_indices = xas_tdp_control%list_ex_atoms
         xas_tdp_env%state_types = xas_tdp_control%state_types

      ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN
         
!        need to find out which atom of which kind is excited
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=at_kind_set)  
         n_kinds = SIZE(at_kind_set) 
         nex_atoms = 0

         DO i = 1,n_kinds
            CALL get_atomic_kind(atomic_kind=at_kind_set(i), element_symbol=symbol,&
                                 natom=nat_of_kind)
            IF (ANY(xas_tdp_control%list_ex_kinds .EQ. symbol)) nex_atoms = nex_atoms + nat_of_kind
         END DO

         ALLOCATE(xas_tdp_env%ex_atom_indices(nex_atoms))
         ALLOCATE(xas_tdp_env%state_types(SIZE(xas_tdp_control%state_types,1),nex_atoms))
         nex_kinds = SIZE(xas_tdp_control%list_ex_kinds)
         nex_atoms = 0

         DO i = 1,n_kinds
            CALL get_atomic_kind(atomic_kind=at_kind_set(i), element_symbol=symbol,&                 
                                                natom=nat_of_kind, atom_list=ind_of_kind)
            DO j = 1,nex_kinds
               IF (xas_tdp_control%list_ex_kinds(j) .EQ. symbol) THEN
                  xas_tdp_env%ex_atom_indices(nex_atoms+1:nex_atoms + nat_of_kind) = ind_of_kind
                  DO k = 1,SIZE(xas_tdp_control%state_types,1)
                     xas_tdp_env%state_types(k,nex_atoms+1:nex_atoms + nat_of_kind) = &
                                                                  xas_tdp_control%state_types(k,j)
                  END DO
                  nex_atoms = nex_atoms + nat_of_kind
               END IF
            END DO
         END DO

         CALL set_xas_tdp_env(xas_tdp_env, nex_atoms=nex_atoms)         

      END IF

!  Allocating memory for the array of donor states
      n_donor_states = COUNT(xas_tdp_env%state_types .NE. xas_not_excited)     
      ALLOCATE(xas_tdp_env%donor_states(n_donor_states))

!  Initializing the qs_loc_env from the LOCALIZE subsection of XAS_TDP (largely inpired by MI's XAS)
!  Closed shell ground state calculation is assumed (spin-unpolarized)
      CALL qs_loc_env_create(xas_tdp_env%qs_loc_env)
      qs_loc_env => xas_tdp_env%qs_loc_env
      loc_section => xas_tdp_control%loc_subsection
!     getting the number of (doubly occupied)  MOs 
      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mos(1)%mo_set, nmo=n_mo(1), homo=homo) 
      n_mo(2) = n_mo(1)
      nspins = 1
!     by default, all homo are localized
      IF (xas_tdp_control%n_search .LT. 0) xas_tdp_control%n_search = homo
      CALL qs_loc_control_init(qs_loc_env, loc_section, do_homo=.TRUE., &  
                               do_xas=.TRUE., nloc_xas=xas_tdp_control%n_search)
!     final qs_loc_env initialization
      IF (.NOT. qs_loc_env%do_localize) THEN
         qs_loc_env%localized_wfn_control%localization_method = do_loc_none
      ELSE IF (qs_loc_env%do_localize) THEN
         n_moloc = qs_loc_env%localized_wfn_control%nloc_states
         CALL set_loc_wfn_lists(qs_loc_env%localized_wfn_control, n_moloc, n_mo, nspins)  
         CALL set_loc_centers(qs_loc_env%localized_wfn_control, n_moloc, nspins)  
!        closed shell => myspin=1
         CALL qs_loc_env_init(qs_loc_env, qs_loc_env%localized_wfn_control, &    
                              qs_env, myspin=1, do_localize=qs_loc_env%do_localize) 
      END IF
      
!  Allocating memory for the array of excited atoms MOs. Worst case senario, all searched MOs are
!  associated to the same atom
   ALLOCATE(xas_tdp_env%mos_of_ex_atoms(xas_tdp_control%n_search, nex_atoms))
      

   END SUBROUTINE xas_tdp_init

!> *************************************************************************************************
!> \brief Assosciate previously localized lowest energy  MOs to each excited atoms. The procedure
!>        looks for MOs "centered" on the excited atom using the same method as the loc process. It
!>        then fills the mos_of_ex_atoms arrays of the xas_tdp_env. Only the xas_tdp_control%n_search
!>        lowest energy MOs are considered. Largely inspired by MI's implementation of XAS
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> *************************************************************************************************
   SUBROUTINE assign_loc_mos_to_ex_atoms(xas_tdp_env, xas_tdp_control, qs_env)
      
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'assign_loc_mos_to_ex_atoms',&
                                     routineP = moduleN//':'//routineN
      
      INTEGER, DIMENSION(:,:), POINTER                   :: mos_of_ex_atoms
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      INTEGER                                            :: iat, imo, n_search, nex_atoms, &
                                                            iat_memo, n_atoms, at_index
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, DIMENSION(3)                              :: perd_init
      REAL(dp), DIMENSION(3)                             :: wfn_center, at_pos, r_ac
      REAL(dp)                                           :: dist, dist_min
      
      NULLIFY(localized_wfn_control, mos_of_ex_atoms, cell, particle_set)

!  Initialization. mos_of_ex_atoms filled with -1, meaning no assigned state 
      mos_of_ex_atoms => xas_tdp_env%mos_of_ex_atoms
      mos_of_ex_atoms(:,:) = -1
      nex_atoms = xas_tdp_env%nex_atoms
      n_search = xas_tdp_control%n_search
      localized_wfn_control => xas_tdp_env%qs_loc_env%localized_wfn_control
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)
      n_atoms = SIZE(particle_set)

!     Temporarly impose periodic BCs in case of Berry's phase operator use for localization
      perd_init(1:3) = cell%perd(1:3)
      cell%perd(1:3) = 1

!  Loop over n_search lowest energy MOs and all atoms 
      DO imo = 1,n_search
!        retrieve MO wave function center coordinates. Closed-shell => centers_set(spin=1)
         wfn_center(1:3) = localized_wfn_control%centers_set(1)%array(1:3, imo) 
         iat_memo = 0

!        a large enough value to avoid surprises
         dist_min = 10000.0_dp
         DO iat = 1,n_atoms
            at_pos(1:3) = particle_set(iat)%r(1:3)
!           for test purposes
            PRINT*, "at_index: ", iat
            PRINT*, "particle_set(iat)%atom_index", particle_set(iat)%atom_index
            r_ac = pbc(at_pos, wfn_center, cell)
            dist = SQRT(r_ac(1)*r_ac(1) + r_ac(2)*r_ac(2) + r_ac(3)*r_ac(3))
         
!           keep memory of which atom is the closest to the wave function center
            IF (dist .LT. dist_min) THEN
               iat_memo = iat
               dist_min = dist
            END IF
         END DO

!        Verify that the closest atom is actually excited and assign the MO if so
         IF (ANY(xas_tdp_env%ex_atom_indices .EQ. iat_memo)) THEN
            at_index = MINLOC(ABS(xas_tdp_env%ex_atom_indices-iat_memo),1)
            mos_of_ex_atoms(imo, at_index) = 1
         END IF
      END DO

!  Go back to initial BCs
      cell%perd = perd_init

   END SUBROUTINE assign_loc_mos_to_ex_atoms

!> *************************************************************************************************
!> \brief Assign a core MO to a given donor_state, taking the type (1S, 2S, etc) into account.
!>        The projection on a representative Slater-type orbital is used as a indicator.
!> \param donor_state the donor_state to which a MO must be assigned
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> *************************************************************************************************
   SUBROUTINE assign_mo_to_donor_state(donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env        
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN ='assign_mo_to_donor_state', &
                                                      routineP = moduleN//':'//routineN

      REAL(dp), DIMENSION(:), POINTER                 :: zeta
      INTEGER, DIMENSION(:), POINTER                  :: lq, nq, first_sgf
      TYPE(sto_basis_set_type), POINTER               :: sto_basis_set
      TYPE(gto_basis_set_type), POINTER               :: sto_to_gto_basis_set, kind_basis_set
      INTEGER, DIMENSION(4,7)                         :: ne
      INTEGER                                         :: zval, i, l, nj, nsgf_sto, nsgf_kind, j,&
                                                         n_search, imo, iat, nao, my_mo, at_index,&
                                                         isgf
      REAL(dp)                                        :: zeff, max_overlap, overlap
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_fm_type), POINTER                       :: mo_coeff
      INTEGER, DIMENSION(:,:), POINTER                :: mos_of_ex_atoms
      REAL(dp), DIMENSION(:,:), POINTER               :: tmp_coeff, overlap_matrix
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      REAL(dp), DIMENSION(:), POINTER                 :: mo_evals
      
      NULLIFY(sto_basis_set, sto_to_gto_basis_set, qs_kind_set, kind_basis_set, lq, nq, zeta)
      NULLIFY(overlap_matrix, mos, mo_coeff, mos_of_ex_atoms, tmp_coeff, first_sgf, particle_set)
      NULLIFY(mo_evals)
      
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, mos=mos, particle_set=particle_set)

!  Contruction of a STO that fits the type of orbital we look for
      ALLOCATE(zeta(1))
      ALLOCATE(lq(1))
      ALLOCATE(nq(1))
!     Retrieving quantum numbers
      IF (donor_state%state_type .EQ. xas_1s_type) THEN
         nq(1) = 1
         lq(1) = 0
      ELSE IF (donor_state%state_type .EQ. xas_2s_type) THEN
         nq(1) = 2
         lq(1) = 0
      ELSE
         CPABORT("Procedure for required type not implemented")
      END IF

!     Getting the atomic number
      CALL get_qs_kind(qs_kind_set(donor_state%kind_index), zeff=zeff)
      zval = INT(zeff)

!     Electronic configureation (copied from MI's XAS)
      ne = 0
      DO l = 1,4
         nj = 2*(l-1)+1
         DO i = l,7
            ne(l,i) = ptable(zval)%e_conv(l-1)-2*nj*(i-l)  
            ne(l,i) = MAX(ne(l, i), 0)      
            ne(l,i) = MIN(ne(l, i), 2*nj)
         END DO
      END DO

!     computing zeta with the Slater sum rules
      zeta(1) = srules(zval, ne, nq(1), lq(1))

!     Allocating memory and initiate STO
      CALL allocate_sto_basis_set(sto_basis_set)
      CALL set_sto_basis_set(sto_basis_set, nshell=1, nq=nq, lq=lq, zet=zeta) 

!  Expanding the STO into (normalized) GTOs for later calculations
      CALL create_gto_from_sto_basis(sto_basis_set=sto_basis_set,&
                                     gto_basis_set=sto_to_gto_basis_set,&
                                     ngauss=xas_tdp_control%n_gauss)
      sto_to_gto_basis_set%norm_type = 2
      CALL init_orb_basis_set(sto_to_gto_basis_set)

!  Retrieving the atomic kind related GTO in which MOs are expanded
      CALL get_qs_kind(qs_kind_set(donor_state%kind_index), basis_set=kind_basis_set)

!  Allocating and computing the overlap between the two basis (they share the same center)
      CALL get_gto_basis_set(gto_basis_set=kind_basis_set, nsgf=nsgf_kind)  
      CALL get_gto_basis_set(gto_basis_set=sto_to_gto_basis_set, nsgf=nsgf_sto)  
      ALLOCATE(overlap_matrix(nsgf_sto,nsgf_kind))

!     Making use of MI's subroutine
      CALL calc_stogto_overlap(sto_to_gto_basis_set ,kind_basis_set, overlap_matrix)

!  Looping over the potential donor states to find the best fit
      mos_of_ex_atoms => xas_tdp_env%mos_of_ex_atoms
      n_search = xas_tdp_control%n_search
      at_index = donor_state%at_index
      iat = MINLOC(ABS(xas_tdp_env%ex_atom_indices-at_index),1)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, eigenvalues=mo_evals)
      ALLOCATE(first_sgf(SIZE(particle_set)))
      CALL get_particle_set(particle_set=particle_set, qs_kind_set=qs_kind_set, first_sgf=first_sgf)
      ALLOCATE(tmp_coeff(nao,1))

      max_overlap = 0.0_dp
      my_mo = 0
      DO imo = 1,n_search
         IF(mos_of_ex_atoms(imo,iat) > 0) THEN
!           Getting the relevant coefficients for the candidate state
            CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=tmp_coeff, start_row=1,&
                                     start_col=imo, n_rows=nao, n_cols=1, transpose=.FALSE.)
                                  
!           Computing the overlap
            overlap = 0.0_dp
            DO i = 1,nsgf_kind
               isgf = first_sgf(at_index)+i-1
               DO j = 1,nsgf_sto
                  overlap = overlap + overlap_matrix(j,i)*tmp_coeff(isgf,1)
               END DO
            END DO

            IF (overlap .GT. max_overlap) THEN
               my_mo = imo
               max_overlap = overlap
            END IF
         END IF
      END DO

!  Dealing with the result
      IF ((my_mo .GT. 0) .AND. (my_mo .LE. n_search)) THEN
!        Updating donor state        
         CALL set_donor_state(donor_state=donor_state, mo_index=my_mo)
         IF (.NOT. xas_tdp_env%qs_loc_env%do_localize) THEN
!           Canonical MO, no change in eigenvalues
            CALL set_donor_state(donor_state=donor_state, energy_eval=mo_evals(my_mo))
         ELSE IF (xas_tdp_env%qs_loc_env%do_localize) THEN
!           MO have been rotated, need epsilon_ii
!            CPABORT("Eigenvalues rotation not implemented yet")
            PRINT*, "Eigenvalue rotation not done yet"
         END IF
      ELSE 
         CPABORT("A core donor state could not be assigned a MO. Increasing NSEARCH might help.")
      END IF

         
!  Clean-up
      CALL deallocate_sto_basis_set(sto_basis_set)
      CALL deallocate_gto_basis_set(sto_to_gto_basis_set)
      DEALLOCATE(nq, lq, zeta, overlap_matrix, tmp_coeff, first_sgf)

   END SUBROUTINE assign_mo_to_donor_state

END MODULE xas_tdp_methods
