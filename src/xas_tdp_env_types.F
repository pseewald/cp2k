!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! ************************************************************************************************** 
!> \brief create/release, get/set and such subroutines for TDP XAS environment types
!> \author AB (11.2017)                                                                             
! **************************************************************************************************
MODULE xas_tdp_env_types
   
#include "./base/base_uses.f90"
  
   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_env_types'

! *** Public data types ***

   PUBLIC :: xas_tdp_environment_type

! **************************************************************************************************
!> \param ...
!> \author AB (11.2017)
! **************************************************************************************************
   TYPE xas_tdp_environment_type

   END TYPE xas_tdp_environment_type

CONTAINS

   SUBROUTINE hi()
      PRINT*, "hi"
   END SUBROUTINE hi
END MODULE xas_tdp_env_types
