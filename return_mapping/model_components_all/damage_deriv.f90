module damage_deriv
	! contains routines for computation of the derivatives of the damage variables
	! with respect to the strain vector

	! load additional modules
	use constants
	use derived_types


    implicit none

contains

	function comp_DomegaDalpha_d(alpha_d, l_char)
		! computes derivative of scalar damage variable omega
		! with respect to internal strain-like damage variable alpha_d
		! domega/dalpha_d
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in) :: alpha_d							! internal strain type damage variable
		real(kind=dbl), intent(in) :: l_char							! characteristic element length

		! return variable
		real(kind=dbl) :: comp_DomegaDalpha_d					! domega/dalpha_d

		! internal variable
		! --------------------------------------------------------------------------

		! computation of domega/dalpha_d
		! .
		! .
		! .
		! comp_DomegaDalpha_d = ...


	end function comp_DomegaDalpha_d


end module damage_deriv






















