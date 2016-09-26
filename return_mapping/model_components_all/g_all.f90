module g_all
	! contains functions that compute the plastic potential functions of the material model
	! and one main (model-indepenent) function that calls all plastic potential functions of the model
	! and returns the function values as a vector

	use constants
	use material_info
	use derived_types

    implicit none

contains

	function g(sigma, alpha)
		! computes all plastic potential values for the given model
		! and returns them as real 1D array
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector

		! return variable
		real(kind=dbl), dimension(n_yf) :: g			! array for yield function values

		! internal variables
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of all plastic potentials (modified leon, ...)
		! returning plastic potential values in g vector
		! --------------------------------------

		g(1) = g_1(sigma, alpha)

	end function g



	function g_1(sigma, alpha)
		! computes first plastic potential function
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl) :: g_1

		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function g_1



end module g_all
