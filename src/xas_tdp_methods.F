!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
 
! **************************************************************************************************
!> \brief Methods for X-Ray absorption spectroscopy (XAS) using TDDFPT  
!> \author AB (11.2017)
! **************************************************************************************************

MODULE xas_tdp_methods

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_methods'

   PUBLIC :: xas_tdp

CONTAINS
                
! **************************************************************************************************
!> \breif Driver for XAS TDDFT calculations.
!> \author AB
!> \note Empty for now...
! **************************************************************************************************
   SUBROUTINE xas_tdp()

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_tdp', routineP = moduleN//':'//routineN
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      PRINT*, "Saying 'hi' from the ", routineP, " routine!"

      CALL timestop(handle)
   END SUBROUTINE xas_tdp
END MODULE xas_tdp_methods
