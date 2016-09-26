module hf_all
	! contains functions that compute the hardening functions of the material model
	! and one main (model-indepenent) function that calls all hardening functions of the model
	! and returns the function values as a vector

	! load additional modules
	use constants
	use material_info
	use derived_types

    implicit none

contains


	function hf(sigma, alpha, yield_flags)
		! computes the hardening functions
		! in case of multi surface plasticity the hardening functions of all yield surfaces
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma				! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha				! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags				! active yield functions vector

		! return variable
		real(kind=dbl), dimension(size(alpha), n_yf) :: hf				! results array

		! internal variable
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of all gradient vector (modified leon, ...)
		! --------------------------------------

		if (yield_flags(1) .eqv. .false.) then			! first yield function is active

			! compute hardening function vector of first yield surface (mod leon)
			! and assign ist to first column of return hardening_function array
			hf(:,1) = hf_1(sigma, alpha)

		end if

!		if (yield_flags(2) .eqv. .false.) then			! second yield function is active
!
!			! compute hardening function vector of second yield surface (mod leon cut off)
!			! and assign it to second column of return hardening function array
!			hf(:,2) = hf_2(sigma, alpha)
!
!		end if


	end function hf



	function hf_1(sigma, alpha)
		! computes the first hardening function
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma						! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha						! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha)) :: hf_1

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function hf_1



end module hf_all
