module damage
	! preforms computations related to the damage part of the material model
	! eg: computation of effective stress vector and
	! computation of damaged elastoplastic tangent

	! load additional modules
	use constants
	use material_info
	use damage_var
	use damage_deriv

    implicit none

contains

	subroutine comp_damage_components(sigma, eps_p, delta_eps_p, alpha, l_char, &
										el_inc_flag, &
										alpha_d, omega, &
										C_ep, C_el_inv, &
										eps, eps_pl_old)
		! computes the nominal stress vector and
		! damaged elastoplastic tangent
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: eps_p						! plastic strain vector
		real(kind=dbl), dimension(6), intent(in) :: delta_eps_p					! plastic strain increment vector
		real(kind=dbl), dimension(:), intent(in) :: alpha						! internal variables vector
		logical, intent(in) :: el_inc_flag										! elastic increment flag, states if eincrement is elastic (true) or plastic (false)
		real(kind=dbl), intent(in) :: l_char									! characteristic element length
		real(kind=dbl), dimension(6), intent(in) :: eps, eps_pl_old				! current strains, plastic strains of last step
		real(kind=dbl), dimension(6,6), intent(in) :: C_el_inv					! inverse elasticity matrix

		! input/output variables
		real(kind=dbl), dimension(6), intent(inout) :: sigma					! stress vector (effective on entry, nominal on exit)
		real(kind=dbl), dimension(6,6), intent(inout) :: C_ep					! elastoplastic tangent (undamaged on entry, damaged on exit)
		real(kind=dbl), dimension(n_int_vars_dam), intent(inout) :: alpha_d		! internal damage variables vector
		real(kind=dbl), dimension(n_omega), intent(inout) :: omega				! damage variables vector

		! internal variables
		! --------------------------------------------------------------------------

		! update of damage variable
		! .
		! .
		! .
		! omega = ...

		! computation of nominal stress with updated damage variable
		! .
		! .
		! .
		! sigma = ...

		! computation of damaged elastoplastic tangent
		C_ep = comp_Cep_damaged(delta_eps_p, alpha_d, omega, &
								alpha, C_ep, C_el_inv, &
								sigma, &
								l_char, el_inc_flag)

	end subroutine comp_damage_components



	function comp_Cep_damaged(delta_eps_p, alpha_d, omega, &
								alpha, C_ep, C_el_inv, &
								sigma, &
								l_char, el_inc_flag)
		! computes damaged elastic plastic tangent
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: delta_eps_p					! plastic strain increment vector
		real(kind=dbl), dimension(:), intent(in) :: alpha_d						! internal damage variables vector
		real(kind=dbl), dimension(:), intent(in) :: omega						! damage variables vector
		real(kind=dbl), dimension(:), intent(in) :: alpha						! internal variables vector
		real(kind=dbl), dimension(6,6), intent(in) :: C_ep, C_el_inv			! elastic plastic tangent, inverse elasticity matrix
		real(kind=dbl), dimension(6), intent(in) :: sigma						! stress vector (effective)
		real(kind=dbl), intent(in) :: l_char									! characteristic element length
		logical, intent(in) :: el_inc_flag										! elastic increment flag, states if eincrement is elastic (true) or plastic (false)

		! output variable
		real(kind=dbl), dimension(6,6) :: comp_Cep_damaged						! damaged elatic plastic tangent

		! internal variables
		! --------------------------------------------------------------------------

	end function comp_Cep_damaged



end module damage
