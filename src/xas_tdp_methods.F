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
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: xas_tdp_by_index,&
                                              xas_tdp_by_kind,&
                                              xas_not_excited
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_get_subs_vals
   USE qs_environment_types,            ONLY: qs_environment_type,&
                                              get_qs_env
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
                                                            current_state_index
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      CHARACTER(len=2)                                   :: symbol_of_kind
      INTEGER, DIMENSION(:), POINTER                     :: atoms_of_kind
      TYPE(donor_state_type), POINTER                    :: current_state

      NULLIFY(xas_tdp_env, xas_tdp_control, atomic_kind_set, atoms_of_kind, current_state)

      CALL timeset(routineN, handle)
      
!  Initialization
      CALL xas_tdp_init(xas_tdp_env, xas_tdp_control, qs_env, dft_control)
      
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
                                          at_symbol=symbol_of_kind,&
                                          state_type=xas_tdp_env%state_types(istate,tmp_index))
!                    Do calculations here...
                     
                     PRINT*, "Test for donor states: "
                     PRINT*, current_state%at_symbol
                     PRINT*, current_state%at_index
                     PRINT*, current_state%state_type

                     current_state_index = current_state_index + 1
                     NULLIFY(current_state)
                  END IF
               END DO
            END IF
         END DO
      END DO

!  Clean-up
      CALL xas_tdp_env_release(xas_tdp_env)

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
      
      TYPE(section_vals_type), POINTER                   :: dft_section, xas_tdp_section
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: at_kind_set
      INTEGER                                            :: nex_atoms, n_kinds, i, nat_of_kind,&
                                                            j, nex_kinds, k, n_donor_states
      INTEGER, DIMENSION(:), POINTER                     :: ex_atom_indices, ind_of_kind
      INTEGER, DIMENSION(:,:), POINTER                   :: state_types
      CHARACTER(len=2)                                   :: symbol

      NULLIFY(dft_section, xas_tdp_section, at_kind_set, ex_atom_indices, state_types, ind_of_kind)

!  XAS TDP control type initialization
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")                     
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")

      CALL xas_tdp_control_create(dft_control%xas_tdp_control)
      CALL read_xas_tdp_control(dft_control%xas_tdp_control,xas_tdp_section)   

      xas_tdp_control => dft_control%xas_tdp_control

!  XAS TDP environment type initialization
      CALL xas_tdp_env_create(xas_tdp_env)

!     Retrieving the excited atoms indices and correspondig state types
      IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

         nex_atoms = SIZE(xas_tdp_control%list_ex_atoms)
         CALL set_xas_tdp_env(xas_tdp_env, nex_atoms=nex_atoms)  
         ALLOCATE(xas_tdp_env%ex_atom_indices(nex_atoms))
         ALLOCATE(xas_tdp_env%state_types(SIZE(xas_tdp_control%state_types,1),nex_atoms))
         xas_tdp_env%ex_atom_indices = xas_tdp_control%list_ex_atoms
         xas_tdp_env%state_types = xas_tdp_control%state_types

      ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN
         
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

!     Allocating memory for the array of donor states
      n_donor_states = COUNT(xas_tdp_env%state_types .NE. xas_not_excited)     
      ALLOCATE(xas_tdp_env%donor_states(n_donor_states))

   END SUBROUTINE xas_tdp_init


END MODULE xas_tdp_methods
