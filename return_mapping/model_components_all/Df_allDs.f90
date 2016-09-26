module Df_allDs
	! contains functions that compute the df/dsigma vectors of the material model
	! and one main (model-indepenent) function that calls all df/dsigma functions of the model
	! and returns the function values as a matrix

	! load additional modules
	use constants
	use material_info
	use derived_types
	use stress_routines
    implicit none

contains

	function DfDs(sigma, alpha)
		! computes the derivatives (gradient vectors) of the yield functions
		! with respect to sigma
		! in case of multi surface plasticity the gradients of all yield functions
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha		! internal variables vector

		! return variable
		real(kind=dbl), dimension(6,n_yf) :: 						DfDs	! results array

		! internal variable
		! --------------------------------------------------------------------------

		! computation of all gradient vector (modified leon, ...)
		! returning hardening function values in hardening_function vector

		DfDs(:,1) = Df_1Ds(sigma, alpha)

	end function DfDs




	function Df_1Ds(sigma, alpha)
		! compute df/dsigma vector of first yield function
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl), dimension(6) :: 				Df_1Ds

		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha		! internal variables vector

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function Df_1Ds



end module Df_allDs
