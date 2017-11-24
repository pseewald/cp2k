!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
 
! **************************************************************************************************
!> \brief Methods for X-Ray absorption spectroscopy (XAS) using TDDFPT  
!> \author AB (11.2017)
! **************************************************************************************************

MODULE xas_tdp_methods

   USE cp_control_types,                ONLY: dft_control_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_get_subs_vals
   USE qs_environment_types,            ONLY: qs_environment_type
   USE xas_tdp_types,                   ONLY: xas_tdp_control_type,&
                                              read_xas_tdp_control,&
                                              xas_tdp_control_create
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_tdp', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(section_vals_type), POINTER                   :: dft_section, xas_tdp_section

      CALL timeset(routineN, handle)
      
      NULLIFY(dft_section, xas_tdp_section)

      PRINT*, "Saying 'hi' from the ", routineP, " routine!"
      PRINT*, "Checking that the xas_tdp_control is created, read and released"

      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")                                 
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")             

      CALL xas_tdp_control_create(dft_control%xas_tdp_control)
      CALL read_xas_tdp_control(dft_control%xas_tdp_control,xas_tdp_section)


!  Structure would look like:
!     - Initialization of the xas_tdp_control
!     - Initialization of the xas_tdp_env (including each excited_atom type)
!     - Loop over excited_atoms array of xas_tdp_env for calculations
!     - Collect and print results
!     - Release all types

      CALL timestop(handle)
   END SUBROUTINE xas_tdp
END MODULE xas_tdp_methods
