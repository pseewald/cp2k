!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI-methods for HFX
! **************************************************************************************************

MODULE hfx_ri


   CONTAINS

   SUBROUTINE hfx_ri_init(ri_data, qs_env)
      TYPE(hfx_ri_type), INTENT(OUT) :: ri_data
      TYPE(gto_basis_set_type), POINTER                  :: ri_aux_basis_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER :: ikind, nkind
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL get_qs_env(qs_env, &
                      nkind=nkind, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control)

      ! initialise RI basis
      DO ikind = 1, nkind
         NULLIFY (ri_aux_basis_set)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=ri_aux_basis_set, &
                          basis_type="RI_HFX")
         IF (.NOT. (ASSOCIATED(ri_aux_basis_set))) THEN
            ! Generate a default basis
            CALL create_ri_aux_basis_set(ri_aux_basis_set, qs_kind, dft_control%auto_basis_ri_aux)
            CALL init_orb_basis_set(ri_aux_basis_set)
            CALL init_interaction_radii_orb_basis(ri_aux_basis_set, dft_control%qs_control%eps_pgf_orb)
            CALL add_basis_set_to_container(qs_kind%basis_sets, ri_aux_basis_set, "RI_HFX")
         END IF
      END DO
   END SUBROUTINE

   SUBROUTINE hfx_ri_pre_scf(qs_env, hfx_ri_data, op_pot, op_ri, t_2c_int, t_3c_int)
      TYPE(hfx_ri_type), INTENT(IN) :: hfx_ri_data
      INTEGER, INTENT(IN) :: op_ri, op_pot
      TYPE(dbcsr_t_type), INTENT(OUT) :: t_3c_int, t_2c_int
      LOGICAL :: same_op
      TYPE(qs_environment_type), POINTER                 :: qs_env

      same_op = op_ri == op_pot
      CPASSERT(same_op)

      CALL build_3c_integrals(t_2c_int, hfx_ri_data%filter_eps, qs_env)
   END SUBROUTINE

   SUBROUTINE hfx_ri_update_ks()
   END SUBROUTINE

   !SUBROUTINE integrate_RI(qs_env, hfx_ri_data, ks_matrix, ehfx, hfx_ri_section, para_env, geometry_did_change, ispin)
   !   TYPE(qs_environment_type), POINTER                 :: qs_env
   !   TYPE(dbcsr_p_type)                                 :: ks_matrix
   !   REAL(KIND=dp), INTENT(OUT)                         :: ehfx
   !   TYPE(section_vals_type), POINTER                   :: hfx_ri_section
   !   TYPE(cp_para_env_type), POINTER                    :: para_env
   !   LOGICAL                                            :: geometry_did_change
   !   INTEGER, INTENT(IN)                                :: ispin
   !END SUBROUTINE

END MODULE
