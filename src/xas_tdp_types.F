!--------------------------------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations                              ! 
!   Copyright (C) 2000 - 2017  CP2K developers group                                               ! 
!--------------------------------------------------------------------------------------------------!

!> *************************************************************************************************
!> \brief Define XAS TDP control type and associated create, release, etc subroutines, as well as
!>        XAS TDP environment type and associated set, get, etc subroutines
!> \author AB (11.2017)
!> *************************************************************************************************
MODULE xas_tdp_types
   
   USE input_constants,                 ONLY: xas_1s_type,&
                                              xas_dip_vel,&
                                              xas_tdp_by_index,&
                                              xas_tdp_by_kind
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length
   USE memory_utilities,                ONLY: reallocate
                                     
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief Type containing control information for TDP XAS calculations
!> \param define_excited whether excited atoms are chosen by kind or index
!> \param dipole form whether the dipole moment is computed in the length or velocity representation
!> \param list_ex_atoms list of excited atom indices, kept empty if define_excited=by_kind
!> \param list_ex_kinds list of excited atom kinds, kept empty if define_excited=by_index
!> \param state_types columns correspond to the states to excite for each atom kind/index
!>                    the number of rows is the number of times the keyword is repeated
! **************************************************************************************************
   TYPE xas_tdp_control_type
      INTEGER                                 :: define_excited   
      INTEGER                                 :: dipole_form
      INTEGER, DIMENSION(:), POINTER          :: list_ex_atoms
      CHARACTER(len=default_string_length),&
                     DIMENSION(:), POINTER    :: list_ex_kinds
      INTEGER, DIMENSION(:,:), POINTER        :: state_types 
      
   END TYPE xas_tdp_control_type
   
!> *************************************************************************************************
!> \brief Type containing informations such as inputs and results for TDP XAS calculations
!> \param nex_atoms number of excited atoms
!> \param ex_atom_indices array containing the indices of the excited atoms
!> \param state_types columns correspond to the different donor states of each excited atom
!> *************************************************************************************************
   TYPE xas_tdp_env_type
      INTEGER                                 :: nex_atoms
      INTEGER, DIMENSION(:), POINTER          :: ex_atom_indices
      INTEGER, DIMENSION(:,:), POINTER        :: state_types
      TYPE(donor_state_type), DIMENSION(:),&
                        POINTER               :: donor_states
   END TYPE xas_tdp_env_type

!> *************************************************************************************************
!> \brief Type containing informations about a single donor state
!> \param index the index of the atom to which the state belongs
!> \param symbol the chemical symbol of the atom to which the state belongs
!> \paran state_type whether this is a 1s, 2s, etc state
!> *************************************************************************************************
TYPE donor_state_type
      INTEGER                                 :: index
      CHARACTER(len=2)                        :: symbol
      INTEGER                                 :: state_type
   END TYPE donor_state_type

   CHARACTER(len=*), PARAMETER, PRIVATE  :: moduleN = 'xas_tdp_types'

! *** Public data types ***
   PUBLIC :: xas_tdp_control_type, xas_tdp_env_type

! *** Public subroutines ***
   PUBLIC :: xas_tdp_control_create, xas_tdp_control_release, read_xas_tdp_control, &
             xas_tdp_env_create, xas_tdp_env_release, get_xas_tdp_env, set_xas_tdp_env

CONTAINS

! **************************************************************************************************
!> \brief Creates and initializes the xas_tdp_control_type
!> \param xas_tdp_control the type to initialize
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_create(xas_tdp_control)
      
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_control_create', &
                                                            routineP = moduleN//':'//routineN
      
      CPASSERT(.NOT. ASSOCIATED(xas_tdp_control))
      ALLOCATE (xas_tdp_control)
      
      xas_tdp_control%dipole_form = xas_dip_vel
      xas_tdp_control%define_excited = xas_tdp_by_index
      NULLIFY (xas_tdp_control%state_types)
      NULLIFY (xas_tdp_control%list_ex_atoms)
      NULLIFY (xas_tdp_control%list_ex_kinds)

      PRINT*, "Create TDP XAS control"

   END SUBROUTINE xas_tdp_control_create

! **************************************************************************************************
!> \brief Releases the xas_tdp_control_type
!> \param xas_tdp_control the type to release
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_release(xas_tdp_control)
      
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_control_release', &
                                                            routineP = moduleN//':'//routineN
      
      IF (ASSOCIATED(xas_tdp_control)) THEN
         IF (ASSOCIATED(xas_tdp_control%list_ex_atoms)) THEN
            DEALLOCATE (xas_tdp_control%list_ex_atoms)
         END IF
         IF (ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
            DEALLOCATE (xas_tdp_control%list_ex_kinds)
         END IF
         IF (ASSOCIATED(xas_tdp_control%state_types)) THEN
            DEALLOCATE(xas_tdp_control%state_types)
         END IF
         DEALLOCATE (xas_tdp_control)
      END IF

   PRINT*, "Release TDP XAS control"
      
   END SUBROUTINE xas_tdp_control_release

! **************************************************************************************************
!> \brief Reads the inputs and stores in xas_tdp_control_type
!> \param xas_tdp_control the type where inputs are stored
!> \param xas_tdp_section the section from which input are read
! **************************************************************************************************
   SUBROUTINE read_xas_tdp_control(xas_tdp_control, xas_tdp_section)

      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section
      
      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_xas_tdp_control', &
                                                            routineP = moduleN//':'//routineN

      INTEGER                                            :: irep,nrep,nexc                                                 
      CHARACTER(len=default_string_length),&
                     DIMENSION(:), POINTER               :: k_list
      INTEGER, DIMENSION(:), POINTER                     :: a_list
      INTEGER, DIMENSION(:), POINTER                     :: t_list

      NULLIFY (k_list,a_list)
      
      CALL section_vals_val_get(xas_tdp_section, "DEFINE_EXCITED", &
                                i_val=xas_tdp_control%define_excited)

      CALL section_vals_val_get(xas_tdp_section, "DIPOLE_FORM", i_val=xas_tdp_control%dipole_form)

      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            ALLOCATE (xas_tdp_control%list_ex_kinds(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN                   

            CALL section_vals_val_get(xas_tdp_section, "KIND_LIST", c_vals=k_list)
         
            IF (ASSOCIATED(k_list)) THEN
               nexc = SIZE(k_list)
               ALLOCATE (xas_tdp_control%list_ex_kinds(nexc))
               xas_tdp_control%list_ex_kinds = k_list
            END IF

         END IF
      END IF

      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_atoms)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN
            
            ALLOCATE (xas_tdp_control%list_ex_atoms(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            CALL section_vals_val_get(xas_tdp_section, "ATOM_LIST", i_vals=a_list)

            IF (ASSOCIATED(a_list)) THEN
               nexc = SIZE(a_list)
               CALL reallocate(xas_tdp_control%list_ex_atoms, 1, nexc)
               xas_tdp_control%list_ex_atoms = a_list
            END IF

         END IF
      END IF

      CALL section_vals_val_get(xas_tdp_section, "STATE_TYPES", n_rep_val=nrep)

      IF(.NOT. ASSOCIATED(xas_tdp_control%state_types)) THEN
         ALLOCATE (xas_tdp_control%state_types(nrep,nexc))
         DO irep = 1,nrep
            CALL section_vals_val_get(xas_tdp_section, "STATE_TYPES", i_rep_val=irep, i_vals=t_list)
            IF (SIZE(t_list) .NE. nexc) THEN
               CPABORT("The STATE_TYPES keywords do not have the correct number of entries.")
            END IF
            xas_tdp_control%state_types(irep,:) = t_list
         END DO
      END IF

      IF (ASSOCIATED(xas_tdp_control) .AND. ASSOCIATED(xas_tdp_section)) THEN
         PRINT*, "Read input for TDP XAS control:"
         PRINT*, "DEFINE_EXCITED: ", xas_tdp_control%define_excited
         PRINT*, "DIPOLE_FORM: ", xas_tdp_control%dipole_form
         PRINT*, "STATE_TYPES", xas_tdp_control%state_types
         PRINT*, "Excited atoms by index:"
         DO irep = 1, SIZE(xas_tdp_control%list_ex_atoms)
            PRINT*, xas_tdp_control%list_ex_atoms(irep)
         END DO
         PRINT*, "Excited kinds: "
         DO irep = 1, SIZE(xas_tdp_control%list_ex_kinds)
            PRINT*, xas_tdp_control%list_ex_kinds(irep)
         END DO
      END IF

   END SUBROUTINE read_xas_tdp_control

!> *************************************************************************************************
!> \brief Creates a TDP XAS environment type
!> \param xas_tdp_env the type to create
!> *************************************************************************************************
   SUBROUTINE xas_tdp_env_create(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_env_create', &
                                                            routineP = moduleN//':'//routineN

      ALLOCATE (xas_tdp_env)

      xas_tdp_env%nex_atoms = 1

      NULLIFY (xas_tdp_env%ex_atom_indices)
      NULLIFY (xas_tdp_env%state_types)
      NULLIFY (xas_tdp_env%donor_states)

   END SUBROUTINE xas_tdp_env_create

!> *************************************************************************************************
!> \brief Releases the TDP XAS environment type
!> \param xas_tdp_env the type to release
!> *************************************************************************************************
   SUBROUTINE xas_tdp_env_release(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env                              
                                                                                                           
      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_env_release', &         
                                                            routineP = moduleN//':'//routineN

      IF (ASSOCIATED(xas_tdp_env)) THEN
         IF (ASSOCIATED(xas_tdp_env%ex_atom_indices)) THEN
            DEALLOCATE(xas_tdp_env%ex_atom_indices)
         END IF
         IF (ASSOCIATED(xas_tdp_env%state_types)) THEN
            DEALLOCATE(xas_tdp_env%state_types)
         END IF
         IF (ASSOCIATED(xas_tdp_env%donor_states)) THEN
            CALL deallocate_donor_state_set(xas_tdp_env%donor_states)
         END IF
         DEALLOCATE (xas_tdp_env)
      END IF
   END SUBROUTINE xas_tdp_env_release

!> *************************************************************************************************
!> \brief Returns selected values contained in a TDP XAS environment type
!***************************************************************************************************
   SUBROUTINE get_xas_tdp_env(xas_tdp_env)
      
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      
      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_xas_env_type',&
                                                            routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(xas_tdp_env))

   END SUBROUTINE get_xas_tdp_env

!> *************************************************************************************************
!> \brief Sets values of selected variables within the TDP XAS environment type
!> *************************************************************************************************
   SUBROUTINE set_xas_tdp_env(xas_tdp_env, nex_atoms)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'set_xas_tdp_env',&
                                                            routineP = moduleN//':'//routineN
                              
      INTEGER, INTENT(IN), OPTIONAL                      :: nex_atoms

      CPASSERT(ASSOCIATED(xas_tdp_env))

      IF (PRESENT(nex_atoms)) xas_tdp_env%nex_atoms = nex_atoms
      
   END SUBROUTINE set_xas_tdp_env
   

!> *************************************************************************************************
!> \brief Deallocate a set of donor states
!> \param donor_state_set the set of donor states to deallocate
!> *************************************************************************************************
SUBROUTINE deallocate_donor_state_set(donor_state_set)
      TYPE(donor_state_type), DIMENSION(:) , POINTER   :: donor_state_set

      CHARACTER(len=*), PARAMETER                       :: routineN = 'release_excited_atomi_set',&
                                                           routineP = moduleN//':'//routineN

      IF (ASSOCIATED(donor_state_set)) THEN
!        Deallocate all other pointers here
         DEALLOCATE(donor_state_set)
      END IF

   END SUBROUTINE deallocate_donor_state_set

END MODULE xas_tdp_types
