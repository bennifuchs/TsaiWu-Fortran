module damage_var
	! contains routines for computation of damage part of the modified leon model

	! load additional modules
	use constants
	use derived_types

    implicit none

contains


	function comp_delta_alpha_d(delta_epsilon_p, alpha)
		! computes increment of internal scalar strain like softening variable,

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: delta_epsilon_p		! incremental plastic strains vector
		real(kind=dbl), dimension(:), intent(in) :: alpha				! internal variables vector (for continuum hardening and softening)

		! return variable
		real(kind=dbl) :: comp_delta_alpha_d					! increment of internal scalar strain like softening variable

		! internal variable
		! --------------------------------------------------------------------------

		! compuation of delta alpha_d
		! .
		! .
		! .
		! comp_delta_alpha_d = ...


	end function comp_delta_alpha_d



end module damage_var
