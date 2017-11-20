!--------------------------------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations                              ! 
!   Copyright (C) 2000 - 2017  CP2K developers group                                               ! 
!--------------------------------------------------------------------------------------------------!

!> *************************************************************************************************
!> \brief Define XAS TDP control type and associated create, release, etc subroutines
!> \author AB (11.2017)
!> *************************************************************************************************
MODULE xas_tdp_control
   
   USE input_constants,                 ONLY: xas_1s_type,&
                                              xas_dip_vel
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
                                     
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief Type containing control information for TDP XAS calculations
! **************************************************************************************************
   TYPE xas_tdp_control_type
      INTEGER                                 :: dipole_form
      INTEGER                                 :: state_type
      INTEGER, DIMENSION(:), POINTER          :: list_ex_atoms
      CHARACTER(len=2), DIMENSION(:), POINTER :: list_ex_kinds
   END TYPE xas_tdp_control_type

   CHARACTER(len=*), PARAMETER, PRIVATE  :: moduleN = 'xas_tdp_control'

! *** Public data types ***
   PUBLIC :: xas_tdp_control_type

! *** Public subroutines ***
   PUBLIC :: xas_tdp_control_create, xas_tdp_control_release, read_xas_tdp_control

CONTAINS

! **************************************************************************************************
!> \brief Creates and initializes the xas_tdp_control_type
!> \param xas_tdp_control, the type to initialize
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_create(xas_tdp_control)
      
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_control_create', &
                                                            routineP = moduleN//':'//routineN
      
      CPASSERT(.NOT. ASSOCIATED(xas_tdp_control))
      ALLOCATE (xas_tdp_control)
      
      xas_tdp_control%dipole_form = xas_dip_vel
      xas_tdp_control%state_type = xas_dip_vel
      NULLIFY (xas_tdp_control%list_ex_atoms)
      NULLIFY (xas_tdp_control%list_ex_kinds)

      PRINT*, "Create TDP XAS control"

   END SUBROUTINE xas_tdp_control_create

! **************************************************************************************************
!> \brief Releases the xas_tdp_control_type
!> \param xas_tdp_control, the type to release
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
!> \param xas_tdp_control, the type where inputs are stored
!> \param xas_tdp_section, the section from which input are read
! **************************************************************************************************
   SUBROUTINE read_xas_tdp_control(xas_tdp_control, xas_tdp_section)

      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section
      
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'read_xas_tdp_control', &
                                                            routineP = moduleN//':'//routineN

      IF (ASSOCIATED(xas_tdp_control) .AND. ASSOCIATED(xas_tdp_section)) THEN
         PRINT*, "Read input for TDP XAS control"
      END IF

   END SUBROUTINE read_xas_tdp_control
END MODULE xas_tdp_control
