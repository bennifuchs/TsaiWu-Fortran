module Dh_allDs
	! contains functions that compute the dh/dsigma matrices of the material model
	! and one main (model-indepenent) function that calls all dh/dsigma functions of the model
	! and returns the function values as a matrix

	use constants
	use material_info
	use derived_types

    implicit none

contains


	function DhDs(sigma, alpha)
		! computes the dh/dsigma matrix
		! in case of multi surface plasticity the dh/dsigma matrices of all yield surfaces
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha),6,n_yf) :: DhDs				! results array of rank 3

		! internal variable
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of all dh/dsigma matrices (modified leon, ...)
		! returning dh/dsigma matrices return array
		! --------------------------------------

		DhDs(:,:,1) = Dh_1Ds(sigma, alpha)

	end function DhDs




	function Dh_1Ds(sigma, alpha)

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma						! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha						! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha,1),6) :: Dh_1Ds

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function Dh_1Ds



end module Dh_allDs
