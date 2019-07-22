!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief 2-center electron repulsion integral routines based on libint2
!>        Currently available operators: Coulomb, Truncated Coulomb, Short Range (erfc)
!> todo: maybe merge with libint_3center
! **************************************************************************************************

MODULE libint_2center
   USE gamma,                           ONLY: fgamma => fgamma_0
   USE hfx_libint_wrapper,              ONLY: cp_libint_get_2eris,&
                                              cp_libint_set_params_eri,&
                                              cp_libint_t,&
                                              prim_data_f_size
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              t_c_g0_n
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_2center'

   PUBLIC :: eri_2center

   TYPE :: params_2c
      INTEGER                               :: m_max
      REAL(dp)                              :: ZetaInv, EtaInv, ZetapEtaInv, Rho
      REAL(dp), DIMENSION(3)                :: W
      REAL(dp), DIMENSION(prim_data_f_size) :: Fm
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Computes the 3-center electron repulsion integrals (ab|c) for a given set of cartesian
!>        gaussian orbitals
!> \param int_ab the integrals as array of cartesian orbitals (allocated before hand)
!> \param la_min ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param ra ...
!> \param lb_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rb ...
!> \param lib the libint_t object for evaluation (assume that it is initialized outside)
!> \param op which operator to use for the integrals
!> \param omega in case of short range erfc operator
!> \param r_cutoff the cutoff radius in case of truncated Coulomb
!> \note Prior to calling this routine, the cp_libint_t type passed as argument must be initialized,
!>       the libint library must be static initialized, and in case of truncated Coulomb operator,
!>       the latter must be initialized too
! **************************************************************************************************
   SUBROUTINE eri_2center(int_ab, la_min, la_max, npgfa, zeta, ra, &
                          lb_min, lb_max, npgfb, zetb, rb, &
                          lib, op, omega, r_cutoff)

      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: int_ab
      INTEGER, INTENT(IN)                                :: la_min, la_max, npgfa
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra
      INTEGER, INTENT(IN)                                :: lb_min, lb_max, npgfb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zetb
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rb
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      INTEGER, INTENT(IN)                                :: op
      REAL(dp), INTENT(IN), OPTIONAL                     :: omega, r_cutoff

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_2center', routineP = moduleN//':'//routineN

      INTEGER                                            :: a_mysize(1), a_offset, a_start, &
                                                            b_offset, b_start, i, ipgf, j, jpgf, &
                                                            li, lj, ncoa, ncob, p1, p2
      REAL(dp)                                           :: my_omega, my_r_cutoff, zeti, zetj
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      TYPE(params_2c), POINTER                           :: params

      NULLIFY (params, p_work)
      ALLOCATE (params)

      !Note: we want to compute all possible integrals based on the 3-centers (ab|c) before
      !      having to switch to (ba|c) (or the other way around) due to angular momenta in libint
      !      For a triplet of centers (k|ji), we can only compute integrals for which lj >= li

      my_r_cutoff = 0.0_dp
      IF (PRESENT(r_cutoff)) my_r_cutoff = r_cutoff

      my_omega = 0.0_dp
      IF (PRESENT(omega)) my_omega = omega

      !Looping over the pgfs
      DO ipgf = 1, npgfa
         zeti = zeta(ipgf)
         a_start = (ipgf-1)*ncoset(la_max)

         DO jpgf = 1, npgfb
            zetj = zetb(jpgf)
            b_start = (jpgf-1)*ncoset(lb_max)

            !start with all the (c|ba) integrals (stendard order) and keep to lb >= la
            CALL set_params(lib, ra, rb, la_max, lb_max, zeti, zetj, op, &
                            params_out=params, r_cutoff=my_r_cutoff, omega=my_omega)

            DO li = la_min, la_max
               a_offset = a_start+ncoset(li-1)
               ncoa = nco(li)
               DO lj = MAX(li, lb_min), lb_max
                  b_offset = b_start+ncoset(lj-1)
                  ncob = nco(lj)

                  a_mysize(1) = ncoa*ncob
                  CALL cp_libint_get_2eris(li, lj, lib, p_work, a_mysize)

                  DO j = 1, ncob
                     p1 = (j-1)*ncoa
                     DO i = 1, ncoa
                        p2 = p1+i
                        int_ab(a_offset+i, b_offset+j) = p_work(p2)
                     END DO
                  END DO

               END DO !lj
            END DO !li

            !swap centers 3 and 4 to compute (c|ab) with lb < la
            CALL set_params(lib, rb, ra, params_in=params)

            DO lj = lb_min, lb_max
               b_offset = b_start+ncoset(lj-1)
               ncob = nco(lj)
               DO li = MAX(lj+1, la_min), la_max
                  a_offset = a_start+ncoset(li-1)
                  ncoa = nco(li)

                  a_mysize(1) = ncoa*ncob
                  CALL cp_libint_get_2eris(lj, li, lib, p_work, a_mysize)

                  DO i = 1, ncoa
                     p1 = (i-1)*ncob
                     DO j = 1, ncob
                        p2 = p1+j
                        int_ab(a_offset+i, b_offset+j) = p_work(p2)
                     END DO
                  END DO

               END DO !li
            END DO !lj

         END DO !jpgf
      END DO !ipgf

      DEALLOCATE (params)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Sets the internals of the cp_libint_t object for integrals of type (k|j)
!> \param lib ..
!> \param rj ...
!> \param rk ...
!> \param lj_max ...
!> \param lk_max ...
!> \param zetj ...
!> \param zetk ...
!> \param op ...
!> \param r_cutoff ...
!> \param omega ...
!> \param params_in external parameters to use for libint
!> \param params_out returns the libint parameters computed based on the other arguments
!> \note The use of params_in and params_out comes from the fact that one might have to swap
!>       centers 1 and 2 because of angular momenta and pretty much all the parameters of libint
!>       remain the same upon such a change => might avoid recomputing things over and over again
! **************************************************************************************************
   SUBROUTINE set_params(lib, rj, rk, lj_max, lk_max, zetj, zetk, &
                         op, r_cutoff, omega, params_in, params_out)

      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rj, rk
      INTEGER, INTENT(IN), OPTIONAL                      :: lj_max, lk_max
      REAL(dp), INTENT(IN), OPTIONAL                     :: zetj, zetk
      INTEGER, INTENT(IN), OPTIONAL                      :: op
      REAL(dp), INTENT(IN), OPTIONAL                     :: r_cutoff, omega
      TYPE(params_2c), OPTIONAL, POINTER                 :: params_in, params_out

      CHARACTER(len=*), PARAMETER :: routineN = 'set_params', routineP = moduleN//':'//routineN

      INTEGER                                            :: l
      LOGICAL                                            :: use_gamma
      REAL(dp)                                           :: omega2, omega_corr, omega_corr2, prefac, &
                                                            R, T, tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: Fm
      TYPE(params_2c), POINTER                           :: params

      !Assume that one of params_in or params_out is present, and that in the latter case, all
      !other optinal arguments are here

      !The internal structure of libint2 is based on 4-center integrals
      !For 2-center, two of those are dummy centers
      !The integral is assumed to be (k|j) where the centers are ordered as:
      !k -> 1, j -> 3 and (the centers #2 & #4 are dummy centers)

      !If external parameters are given, just use them
      IF (PRESENT(params_in)) THEN
         params => params_in

         !If no external parameters to use, compute them
      ELSE
         params => params_out

         !Note: some variable of 4-center integrals simplify due to dummy centers:
         !      P -> rk, gammap -> zetk
         !      Q -> rj, gammaq -> zetj
         params%m_max = lj_max+lk_max
         params%ZetaInv = 1._dp/zetk; params%EtaInv = 1._dp/zetj
         params%ZetapEtaInv = 1._dp/(zetk+zetj)

         params%W = (zetk*rk+zetj*rj)*params%ZetapEtaInv
         params%Rho = zetk*zetj/(zetk+zetj)

         SELECT CASE (op)
         CASE (do_potential_coulomb)
            T = params%Rho*SUM((rj-rk)**2)
            prefac = 2._dp*pi/params%Rho*SQRT((pi*params%ZetapEtaInv)**3)
            params%Fm = 0.0_dp
            CALL fgamma(params%m_max, T, params%Fm)
            params%Fm = prefac*params%Fm
         CASE (do_potential_truncated)
            R = r_cutoff*SQRT(params%Rho)
            T = params%Rho*SUM((rj-rk)**2)
            prefac = 2._dp*pi/params%Rho*SQRT((pi*params%ZetapEtaInv)**3)

            CPASSERT(get_lmax_init() .GE. params%m_max) !check if truncated coulomb init correctly
            CALL t_c_g0_n(params%Fm, use_gamma, R, T, params%m_max)
            IF (use_gamma) CALL fgamma(params%m_max, T, params%Fm)
            params%Fm = prefac*params%Fm
         CASE (do_potential_short)
            T = params%Rho*SUM((rj-rk)**2)
            prefac = 2._dp*pi/params%Rho*SQRT((pi*params%ZetapEtaInv)**3)

            params%Fm = 0.0_dp
            CALL fgamma(params%m_max, T, params%Fm)

            omega2 = omega**2
            omega_corr2 = omega2/(omega2+params%Rho)
            omega_corr = SQRT(omega_corr2)
            T = T*omega_corr2
            ALLOCATE (Fm(prim_data_f_size))

            CALL fgamma(params%m_max, T, Fm)
            tmp = -omega_corr
            DO l = 1, params%m_max+1
               params%Fm(l) = params%Fm(l)+Fm(l)*tmp
               tmp = tmp*omega_corr2
            END DO
            params%Fm = prefac*params%Fm
         CASE DEFAULT
            CPABORT("Requested operator NYI")
         END SELECT

      END IF

      CALL cp_libint_set_params_eri(lib, rk, rk, rj, rj, params%ZetaInv, params%EtaInv, &
                                    params%ZetapEtaInv, params%Rho, rk, rj, params%W, &
                                    params%m_max, params%Fm)

   END SUBROUTINE set_params

END MODULE

