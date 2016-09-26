module DDg_allDsDa
	! contains functions that compute the ddg/dsigmadalpha matrices of the material model
	! and one main (model-indepenent) function that calls all ddg/dsigmadalpha functions of the model
	! and returns the function values as a matrix

	! load additional modules
	use constants
	use material_info
	use derived_types
	use stress_routines
    implicit none

contains

	function DDgDsDa(sigma, alpha)
		! computes the ddg/dsigmadalpha matrices of the plastic potentials
		! with respect to sigma and alpha
		! in case of multi surface plasticity the matrices of all plastic potentials
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 		sigma				! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 		alpha				! internal variables vector

		! return variable
		real(kind=dbl), dimension(6,size(alpha),n_yf) :: 	DDgDsDa		! results array

		! --------------------------------------------------------------------------


		! --------------------------------------
		! returning ddg/dsigmadalpha matrices in DDgDsigmaDalpha array
		! --------------------------------------

		DDgDsDa(:,:,1) = DDg_1DsDa(sigma, alpha)

	end function DDgDsDa



	function DDg_1DsDa(sigma, alpha)
		! computes the ddg/dsigmadalpha matrix for the material potential
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha		! internal variables vector

		! return variable
		real(kind=dbl), dimension(6,size(alpha)) :: 	DDg_1DsDa

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function DDg_1DsDa




end module DDg_allDsDa
