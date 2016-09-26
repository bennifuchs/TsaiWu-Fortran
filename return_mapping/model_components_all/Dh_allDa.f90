module Dh_allDa
	! contains functions that compute the dh/dalpha matrices of the material model
	! and one main (model-indepenent) function that calls all dh/dalpha functions of the model
	! and returns the function values as a matrix

	! load additional modules
	use constants
	use material_info
	use derived_types

    implicit none

contains


	function DhDa(sigma, alpha)
		! computes the dh/dalpha matrix
		! in case of multi surface plasticity the dh/dalpha matrices of all yield surfaces
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 					sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 					alpha		! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha),size(alpha), n_yf) :: 	DhDa	! results array of rank 3

		! internal variable
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of all dh/dalpha matrices (modified leon, ...)
		! returning dh/dsigma matrices return array
		! --------------------------------------

		DhDa(:,:,1) = Dh_1Da(sigma, alpha)

	end function DhDa




	function Dh_1Da(sigma, alpha)
		! compute dh/dalpha matrix for the first hardening function
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 				sigma				! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 				alpha				! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha),size(alpha)) :: 	Dh_1Da

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function Dh_1Da


end module Dh_allDa
