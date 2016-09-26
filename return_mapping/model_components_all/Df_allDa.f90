module Df_allDa
	! contains functions that compute the df/dalpha vectors of the material model
	! and one main (model-indepenent) function that calls all df/dalpha functions of the model
	! and returns the function values as a matrix

	! load additional modules
	use constants
	use material_info
	use derived_types
	use stress_routines
    implicit none

contains

	function DfDa(sigma, alpha)
		! computes the derivatives (gradient vectors) of the yield functions
		! with respect to the alpha vector
		! in case of multi surface plasticity the gradients of all yield functions
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha		! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha),n_yf) :: 			DfDa	! results array

		! --------------------------------------------------------------------------


		! --------------------------------------
		! computation of all gradient vectors
		! returning hardening function values in hardening_function vector
		! --------------------------------------

		DfDa(:,1) = Df_1Da(sigma, alpha)

	end function DfDa



	function Df_1Da(sigma, alpha)
		! gradient vector of first yield function

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha,1)) :: Df_1Da

		! internal variables
		! --------------------------------------------------------------------------

		! ...


	end function Df_1Da



end module Df_allDa
