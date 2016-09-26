module Dg_allDs
	! contains functions that compute the derivatives of the plastic potentials dg/dsigma
	! of the material model
	! and one main (model-indepenent) function that calls all dg/dsigma functions of the model
	! and returns the derivative vectors as an array

	use constants
	use material_info
	use derived_types
	use yf_all, only: comp_tsaiwu_tensors

    implicit none

contains


	function DgDs(sigma, alpha, yield_flags)
		! computes the derivatives (gradient vectors) of the plastic potentials
		! with respect to sigma
		! in case of multi surface plasticity the gradients of all plastic potentials
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags			! active yield functions vector

		! return variable
		real(kind=dbl), dimension(6, n_yf) :: DgDs						! results array

		! internal variable
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of all gradient vector (modified leon, ...)
		! --------------------------------------

		if (yield_flags(1) .eqv. .false.) then		! check if first yield function is active and if not skip computation of derivatives
			! first yield function is active
			! compute gradient vector of first plastic potential (mod leon)
			! and assign it to the first line of DgDsigma array

			DgDs(:,1) = Dg_1Ds(sigma, alpha)		! first plastic potential

		end if

!		if (yield_flags(2) .eqv. .false.) then		! check if second yield function is active and if not skip computation of derivatives
!			! second yield function is active
!			! compute gradient vector of second plastic potential (mod leon cut off)
!			! and assign it to the second line of DgDsigma array
!			DgDs(:,2) = Dg_2Ds(sigma, alpha)
!
!		end if

	end function DgDs



	function Dg_1Ds(sigma, alpha)
		! computes the gradient vector dg/dsigma for the first plastic potential g
		! derivatives of the tsaiwu yield function (associated flow rule) for the present case

		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl), dimension(6) :: Dg_1Ds

		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma						! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha						! internal variables vector

		! internal variables
		real(kind=dbl), dimension(6) :: struc_vec				! transversely isotropic tsai-wu structural vector
		real(kind=dbl), dimension(6,6) :: struc_matrix			! transversely isotropic tsai-wu structural matrix
		! --------------------------------------------------------------------------

		! compute structural vector and matrix
		call comp_tsaiwu_tensors(struc_vec, struc_matrix)

		! compute tsaiwu gradient vector
		Dg_1Ds = struc_vec + 2._dbl*matmul(struc_matrix, sigma)

	end function Dg_1Ds



end module Dg_allDs
