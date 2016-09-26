module DDg_allDDs
	! contains functions that compute the ddg/ddsigma matrices (hessian matrices) of the material model
	! and one main (model-indepenent) function that calls all ddg/ddsigma functions of the model
	! and returns the function values as a matrix

	! load additional modules
	use constants
	use material_info
	use derived_types
	use stress_routines
    implicit none

contains

	function DDgDDs(sigma, alpha)
		! computes the 2nd order derivatives (hessian matrices) of the plastic potentials
		! with respect to sigma
		! in case of multi surface plasticity the hessian matrices of all plastic potentials
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector

		! return variable
		real(kind=dbl), dimension(6,6,n_yf) :: 			DDgDDs		! results array

		! internal variable
		! --------------------------------------------------------------------------

		! --------------------------------------
		! returning hessian matrices in DDgDDsigma array
		! --------------------------------------

		DDgDDs(:,:,1) = DDg_1DDs(sigma, alpha)

	end function DDgDDs



	function DDg_1DDs(sigma, alpha)
		! computes ddg/ddsigma matrix (hessian) for the first plastic potential
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl), dimension(6,6) :: 				DDg_1DDs

		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha		! internal variables vector

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function DDg_1DDs



end module DDg_allDDs
