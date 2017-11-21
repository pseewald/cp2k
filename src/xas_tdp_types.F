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
! **************************************************************************************************
   TYPE xas_tdp_control_type
      INTEGER                                 :: define_excited   
      INTEGER                                 :: dipole_form
      INTEGER                                 :: state_type
      INTEGER, DIMENSION(:), POINTER          :: list_ex_atoms
      CHARACTER(len=default_string_length),&
                     DIMENSION(:), POINTER    :: list_ex_kinds
   END TYPE xas_tdp_control_type
   
!> *************************************************************************************************
!> \brief Type containing informations such as inputs and results for TDP XAS calculations
!> \param nex_atoms number of excited atoms
!> \param ex_atoms list of indices of excited atoms
!> *************************************************************************************************
   TYPE xas_tdp_env_type
      INTEGER                                 :: nex_atoms
      INTEGER, DIMENSION(:), POINTER          :: ex_atoms
   END TYPE xas_tdp_env_type

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
      xas_tdp_control%state_type = xas_1s_type
      xas_tdp_control%define_excited = xas_tdp_by_index
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

      INTEGER                                            :: i                                                   
      CHARACTER(len=default_string_length),&
                     DIMENSION(:), POINTER               :: k_list
      INTEGER, DIMENSION(:), POINTER                     :: a_list

      NULLIFY (k_list,a_list)
      
      CALL section_vals_val_get(xas_tdp_section, "DEFINE_EXCITED", &
                                i_val=xas_tdp_control%define_excited)

      CALL section_vals_val_get(xas_tdp_section, "DIPOLE_FORM", i_val=xas_tdp_control%dipole_form)

      CALL section_vals_val_get(xas_tdp_section, "STATE_TYPE", i_val=xas_tdp_control%state_type)
      
      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            ALLOCATE (xas_tdp_control%list_ex_kinds(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN                   

            CALL section_vals_val_get(xas_tdp_section, "KIND_LIST", c_vals=k_list)
         
            IF (ASSOCIATED(k_list)) THEN
               ALLOCATE (xas_tdp_control%list_ex_kinds(SIZE(k_list)))
               DO i = 1, SIZE(k_list)
                  xas_tdp_control%list_ex_kinds(i) = k_list(i)
               END DO
            END IF

         END IF
      END IF

      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_atoms)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN
            
            ALLOCATE (xas_tdp_control%list_ex_atoms(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            CALL section_vals_val_get(xas_tdp_section, "ATOM_LIST", i_vals=a_list)

            IF (ASSOCIATED(a_list)) THEN
               CALL reallocate(xas_tdp_control%list_ex_atoms, 1, SIZE(a_list))
               DO i = 1, SIZE(a_list)
                  xas_tdp_control%list_ex_atoms(i) = a_list(i)
               END DO
            END IF

         END IF
      END IF

      IF (ASSOCIATED(xas_tdp_control) .AND. ASSOCIATED(xas_tdp_section)) THEN
         PRINT*, "Read input for TDP XAS control:"
         PRINT*, "DEFINE_EXCITED: ", xas_tdp_control%define_excited
         PRINT*, "DIPOLE_FORM: ", xas_tdp_control%dipole_form
         PRINT*, "STATE_TYPE", xas_tdp_control%state_type
         PRINT*, "Excited atoms by index:"
         DO i = 1, SIZE(xas_tdp_control%list_ex_atoms)
            PRINT*, xas_tdp_control%list_ex_atoms(i)
         END DO
         PRINT*, "Excited kinds: "
         DO i = 1, SIZE(xas_tdp_control%list_ex_kinds)
            PRINT*, xas_tdp_control%list_ex_kinds(i)
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

      NULLIFY (xas_tdp_env%ex_atoms)

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
         IF (ASSOCIATED(xas_tdp_env%ex_atoms)) THEN
            DEALLOCATE (xas_tdp_env%ex_atoms)
         END IF
         DEALLOCATE (xas_tdp_env)
      END IF
   END SUBROUTINE xas_tdp_env_release

!> *************************************************************************************************
!> \brief Returns selected values contained in a TDP XAS environment type
!***************************************************************************************************
   SUBROUTINE get_xas_tdp_env(xas_tdp_env, nex_atoms, ex_atoms)
      
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      INTEGER, INTENT(OUT), OPTIONAL                     :: nex_atoms
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: ex_atoms
      
      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_xas_env_type',&
                                                            routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(xas_tdp_env))

      IF (PRESENT(nex_atoms)) nex_atoms = xas_tdp_env%nex_atoms
      IF (PRESENT(ex_atoms)) ex_atoms => xas_tdp_env%ex_atoms
         

   END SUBROUTINE get_xas_tdp_env

!> *************************************************************************************************
!> \brief Sets values of selected variables within the TDP XAS environment type
!> *************************************************************************************************
   SUBROUTINE set_xas_tdp_env(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'set_xas_tdp_env',&
                                                            routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(xas_tdp_env))
      
   END SUBROUTINE set_xas_tdp_env
END MODULE xas_tdp_types
