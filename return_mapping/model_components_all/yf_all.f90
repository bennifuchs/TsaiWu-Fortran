module yf_all
	! contains functions that compute the yield functions of the material model
	! and one main (model-indepenent) function that calls all yield functions of the model
	! and returns the function values as a vector

	use constants
	use material_info
	use derived_types

    implicit none

contains


	function yf(sigma, alpha)
		! computes all yield function values for the given model
		! and returns them as real 1D array
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha							! internal variables vector

		! return variable
		real(kind=dbl), dimension(n_yf) :: yf					! array for yield function values

		! internal variables
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of each yield function (modified leon, ...)
		! returning yield functions in yield_function vector
		! --------------------------------------

		yf(1) = yf_1(sigma, alpha)

	end function yf


	function yf_ret_map(sigma, alpha, yield_flags)
		! computes yield function values for the given model
		! based on the entries of the yield flags vector
		! and returns them as real 1D array
		! if the yield function isn't active the corresponding yield function value won't be computed
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha			! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags			! active yield functions vector

		! return variable
		real(kind=dbl), dimension(n_yf) :: yf_ret_map				! array for yield function values

		! internal variables
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of each yield function (modified leon, ...)
		! --------------------------------------

		!write(*,*) size(yield_function)

		if (yield_flags(1) .eqv. .false.) then			! first yield function is active

			! call first yield function (modified leon)
			! and assign its value to first line of yield function return vector
			!yf_1 = yf_modleon( sigma, alpha, mat_data)
			yf_ret_map(1) = yf_1(sigma, alpha)

		end if

!		if (yield_flags(2) .eqv. .false.) then			! second yield function is active
!
!			! call second yield function (modified leon cut off)
!			! and assign its value to second line of yield function return vector
!			!yf_2 = yf_cutoff_modleon(sigma, alpha, mat_data)
!			yf_ret_map(2) = yf_2(sigma, alpha)
!
!		end if

	end function yf_ret_map



	function yf_1(sigma, alpha)
		! computes the first yield function
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl) :: yf_1

		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha		! internal variable vector

		! internal variables
		real(kind=dbl), dimension(6) :: struc_vec				! transversely isotropic tsai-wu structural vector
		real(kind=dbl), dimension(6,6) :: struc_matrix			! transversely isotropic tsai-wu structural matrix
		! --------------------------------------------------------------------------

		! compute structural vector and matrix
		call comp_tsaiwu_tensors(struc_vec, struc_matrix)

		! compute tsai-wu yield function
		yf_1 = dot_product(struc_vec, sigma) + dot_product(matmul(struc_matrix, sigma), sigma) - 1._dbl

	end function yf_1



	subroutine comp_tsaiwu_tensors(struc_vec, struc_matrix)
		! computes the tsai-wu structural vector and matrix
		! for the transversely isotropic case
		! from the model parameters imported via the derived_types module
		implicit none

		! --------------------------------------------------------------------------
		! return variables
		real(kind=dbl), dimension(:), intent(out) :: struc_vec
		real(kind=dbl), dimension(:,:), intent(out) :: struc_matrix
		! --------------------------------------------------------------------------

		struc_vec = (/ mat_pars%F2, mat_pars%F2, mat_pars%F3, &
						0._dbl, 0._dbl, 0._dbl /)

		struc_matrix = reshape( (/mat_pars%F22, mat_pars%F12, mat_pars%F23, 0._dbl, 0._dbl, 0._dbl, &
								  mat_pars%F12, mat_pars%F22, mat_pars%F23, 0._dbl, 0._dbl, 0._dbl, &
								  mat_pars%F23, mat_pars%F23, mat_pars%F33, 0._dbl, 0._dbl, 0._dbl, &
								  0._dbl, 0._dbl, 0._dbl, mat_pars%F44, 0._dbl, 0._dbl, &
								  0._dbl, 0._dbl, 0._dbl, 0._dbl, mat_pars%F44, 0._dbl, &
								  0._dbl, 0._dbl, 0._dbl, 0._dbl, 0._dbl, 2._dbl*(mat_pars%F22-mat_pars%F12)/), &
							   (/6,6/) )

	end subroutine comp_tsaiwu_tensors



end module yf_all
