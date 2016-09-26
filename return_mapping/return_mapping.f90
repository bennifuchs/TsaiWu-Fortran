module return_mapping
    ! contains functions which perform the return mapping step
    ! of an elastic-plastic constitutive model

    ! load additional modules
    use constants
    use aux_routines, only: pack_x_vector, &
    						unpack_x_vector, &
    						write_real_list
    use residuals, only: comp_R
    use jacobian, only: comp_jacobian, &
    					comp_jacobian_num
    use interface_definitions

    implicit none

contains

    subroutine multsurf_return_mapping(sigma_trial, alpha_old, &
    										sigma_new, alpha_new, delta_lambda_act, &
    										C_ep, &
    										C_el, &
											n_act_yf, yield_flags, &
											k, status)
		! performs multisurface return mapping step
		! based on an implicit integration scheme (implicit euler)
		! the nonlinear system of equations is solved with the newton method
		use material_info
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) ::						sigma_trial		! trial stress vector
		real(kind=dbl), dimension(:), intent(in) ::						alpha_old		! internal variables vector of the prior step
		real(kind=dbl), dimension(:,:), intent(in) ::					C_el			! elasticity matrix
		integer, intent(in) ::											n_act_yf		! number of currently active yield functions
		logical, dimension(:), intent(in) ::							yield_flags		! active yield functions vector

		! return variables
		real(kind=dbl), dimension(:), intent(out) ::	sigma_new				! updated stress vector
		real(kind=dbl), dimension(:), intent(out) ::			alpha_new		! updated internal variables vector
		real(kind=dbl), dimension(:), intent(out) ::				delta_lambda_act	! vector with plastic multipliers of active yield functions
		real(kind=dbl), dimension(:,:), intent(out) :: C_ep						! elastic-plastic tangent matrix
		integer, intent(out) ::											k				! counter variable for the number of itertations performed
		integer, intent(out) ::											status			! status variable

		! internal variables
		real(kind=dbl), dimension(6+n_int_vars+n_act_yf) :: x			! vector with variables which need to be solved (sigma, alpha, delta_lambda)
		real(kind=dbl), dimension(size(x)) ::							delta_x			! vector with increments of variables which need to be solved
		real(kind=dbl), allocatable, dimension(:) :: 					delta_x_rel_change	! vector with the relative change of the delta_x vector compared to the x vector
		real(kind=dbl) ::												delta_x_rel_change_norm		! norm of the relative change vector
		real(kind=dbl), dimension(size(x)) ::							R				! residual vector
		real(kind=dbl), dimension(size(x)) ::							S				! scaling vector for obtaining dimensionless residual vector
		real(kind=dbl) :: R_norm														! L2 norm of scaled (dimensionless) residual vector
		real(kind=dbl), dimension(size(x),size(x)) ::					Jac				! jacobian matrix
		real(kind=dbl), dimension(size(x),size(x)) ::					Jac_inv			! inverted jacobian matrix

		real(kind=dbl) ::		atol													! convergence tolerance for norm of resiual vector
		real(kind=dbl) ::		rtol													! convergence tolerance for norm of relative change of increment of x vector

		real(kind=dbl) ::		dnrm2													! datatype of BLAS level1 function for L2 norm of a vector

		! internal variables needed for solving the system of linear equations by LAPACK
		integer, dimension(size(x)) ::				lineq_piv_ind						! vector with pivot indices which define the premutation matrix (called IPIV in LAPACK)
		real(kind=dbl), dimension(size(x)) ::		lineq_resid							! residual vector (called WORK in LAPACK)
		real, dimension(size(x)*(size(x)+1)) ::		lineq_s_prec_vals					! vector containing single precision values (called SWORK in LAPACK)
		integer ::									lineq_iter							! status variable about the iterative refinement and when indicated the number of iterations
		integer ::									lapack_info							! exit status variable of the equation solving routine

		integer, dimension(n_act_yf) ::		act_yf_ind		! vector with indices of active yield functions, i.e. indices of entries of the yield flag vector which are false
		! --------------------------------------------------------------------------

		! --------------------------------------------------------------------------
		! INITIALIZATION OF VARIABLES FOR NEWTON ITERTATION
		! --------------------------------------------------------------------------

		! allocate delta_x_rel_change vector with right size depending on the presence of hardening
		if (hardening_present) then
			allocate(delta_x_rel_change(3))
		else
			allocate(delta_x_rel_change(2))
		end if

		! compute indices of active yield functions and store in act_yf_ind vector
		act_yf_ind = comp_act_yf_ind(n_act_yf, yield_flags)

		! initialize active delta lambdas vector with zeros
		delta_lambda_act = 0._dbl

		! initialize x vector with starting values of sigma, alpha and delta_lambda_act
		call pack_x_vector(sigma_trial, alpha_old, delta_lambda_act, x)

		! compute residual vector with initialized x vector
		R = comp_R(x, sigma_trial, alpha_old, &
					C_el, &
					n_act_yf, yield_flags, act_yf_ind)

!		write(*,*) 'initial R vector: ', R

		! compute scaling vector for the residual vector (in valentini2011 defined as scaling matrix)
		S(1:6) = 1._dbl/maxval(abs(sigma_trial))

		if (hardening_present) then
			! check if alpha_old vector is positive
			if (maxval(abs(alpha_old)) > 0._dbl) then
				! alpha_old vector is positive, set corresponding values to 1/max(abs(alpha_old))
				S(7:6+n_int_vars) = 1._dbl/maxval(abs(alpha_old))
			else
				! alpha_old vector is zero, set corresonding values in S to one
				S(7:6+n_int_vars) = 1._dbl
			end if
		end if

		S(7+n_int_vars:6+n_int_vars+n_act_yf) = 1._dbl

		! --------------------------------------------------------------------------
		! START OF NEWTON ITERTATION
		! --------------------------------------------------------------------------

		! start newton iteration
		! initialize k with zero
		k = 0
		inner_newton_loop: do

			! check if maximum number of iterations is reached
			max_iter_check: if (k > k_max) then

				! max number of allowed iterations reached
				! state error message and exit subroutine
				!write(*,*) 'Maxium number of allowed iterations reached. Increment does not converge!'
				status = -10
				return

			end if max_iter_check

			! check if jacobian matrix should be computed analytically or numerically
			if (jac_num_flag .eqv. .False.) then

				! compute jacobian matrix analytically based on the values of the x vector
				Jac = comp_jacobian(x, C_el, &
									n_act_yf, yield_flags, act_yf_ind)

			else if (jac_num_flag .eqv. .True.) then

				! compute jacobian matrix numerically
				Jac = comp_jacobian_num(x, sigma_trial, alpha_old, &
										C_el, &
										n_act_yf, yield_flags, act_yf_ind)

			end if

!			write(*,*) "Jacobian matrix: "
!			call write_real_list(Jac, size(Jac,1), size(Jac,2))

			! solve system of linear equations [J]*{delta_x} = {-R} for {delta_x}
			! call to LAPACK subroutine dgesv
			R = -1._dbl*R
			call dgesv( size(x), 1, &
						Jac, size(x), lineq_piv_ind, &
						R, size(x), &
						lapack_info )
			! call to LAPACK subroutine dsgesv
			! uses iterative refinement
!			call dsgesv( size(x), 1, &
!							J, size(x), lineq_piv_ind, &
!							-1._dbl*R, size(x), &
!							delta_x, size(x), &
!							lineq_resid, lineq_s_prec_vals, &
!							lineq_iter, lapack_info )

			! check if equations were solved successfully
			eq_check: if (lapack_info == 0) then

				! computaion was ok
				! assign results (which are written to the R vector on exit of the LAPCK routine)
				! to the delta_x vector
				delta_x = R

			else eq_check

				! error occurred
				! state error message and leave subroutine
				!write(*,*) 'Error while solving the system of linear equations!'
				status = -100
				return

			end if eq_check

			! update x vector
			x = x + delta_x

			! compute updated residual vector with updated x vector
			R = comp_R(x, sigma_trial, alpha_old, &
						C_el, &
						n_act_yf, yield_flags, act_yf_ind)

!			write(*,*) '		R: ', R

			! update k
			k = k+1

			! compute norm of updated scaled (dimensionless) residual vector
			! (call to BLAS level 1 function)
			R_norm = dnrm2( size(x), R*S, 1 )

!			write(*,*) '		R_norm: ', R_norm

			! compute the relative change of the delta_x vector compared to the x vector
			! (call to BLAS level 1 function for computations of L2 norms)
			delta_x_rel_change(1) = dnrm2( 6, &
											delta_x(1:6), &
											1) / &
									dnrm2( 6, &
											x(1:6), &
											1)

			if (hardening_present) then
				delta_x_rel_change(2) = dnrm2( n_int_vars, &
												delta_x(7:6+n_int_vars), &
												1 ) / &
										dnrm2( n_int_vars, &
												x(7:6+n_int_vars), &
												1 )

				delta_x_rel_change(3) = dnrm2( n_act_yf, &
												delta_x(7+n_int_vars: &
														6+n_int_vars+n_act_yf), &
												1 ) / &
										dnrm2( n_act_yf, &
												x(7+n_int_vars: &
													6+n_int_vars+n_act_yf), &
												1 )
			else
				delta_x_rel_change(2) = dnrm2( n_act_yf, &
								delta_x(7+n_int_vars: &
										6+n_int_vars+n_act_yf), &
								1 ) / &
						dnrm2( n_act_yf, &
								x(7+n_int_vars: &
									6+n_int_vars+n_act_yf), &
								1 )
			end if

			! compute the norm of the relative change vector
			! (call to BLAS level 1 function for computations of L2 norms)
			delta_x_rel_change_norm = dnrm2( size(delta_x_rel_change), delta_x_rel_change, 1 )
!			write(*,*) '		relative delta-x norm: ', delta_x_rel_change_norm

			! check the number of iterations and assign corresponding values to the tolerance variables
			n_iter_check: if (k < k_alt) then

				atol = atol_def
				rtol = rtol_def

			else if (k >= k_alt) then n_iter_check

				atol = atol_alt
				rtol = rtol_alt

			end if n_iter_check

			! check if both convergence criteria are met
			converg_check: if (R_norm < atol .and. delta_x_rel_change_norm < rtol) then

				! check was ok
				! compute jacobian with converged values (needed for C_ep)
				! check if jacobian matrix should be computed analytically or numerically
				if (jac_num_flag .eqv. .False.) then

					! compute jacobian analytically
					Jac = comp_jacobian(x, C_el, &
										n_act_yf, yield_flags, act_yf_ind)

				else if (jac_num_flag .eqv. .True.) then

					! compute jacobian matrix numerically
					Jac = comp_jacobian_num(x, sigma_trial, alpha_old, &
											C_el, &
											n_act_yf, yield_flags, act_yf_ind)

				end if

				! leave loop
				exit inner_newton_loop

			end if converg_check

		end do inner_newton_loop

!		write(*,*) 'no inner it: ', k
!		write(*,*) 'yf combination: ', yield_flags

		! change signs of first 6 lines of jacobian matrix
		Jac(1:6,:) = Jac(1:6,:)*(-1._dbl)

		! invert jacobian matrix
		call mat_inverse(Jac, Jac_inv, status)

		! check if inversion was ok
		if (status /= 0) then

			! inversion failed
			!write(*,*) 'Error in matrix inversion!'
			return

		end if

		! compute C_ep
		C_ep = matmul(Jac_inv(1:6,1:6),C_el)

		! unpack x vector
		call unpack_x_vector(x, sigma_new, alpha_new, delta_lambda_act)

		deallocate(delta_x_rel_change)

		! set status variable to zero and return
		status = 0
		return

    end subroutine multsurf_return_mapping


	subroutine mat_inverse(A, A_inv, status)
		! computes the inverse of a square double precision matrix
		! by applying two LAPACK routines
		! firstly LU decomposition
		! secondly computation of inverse matrix from LU matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variable
		real(kind=dbl), dimension(:,:), intent(in) :: A				! square matrix which is going to be inverted

		! return variable
		real(kind=dbl), dimension(size(A,1),size(A,1)), intent(out) :: A_inv		! inverted matrix
		integer, intent(out) :: status								! status variable

		! internal variables
		integer,dimension(size(A,1)) :: ipiv						! pivot indices vector used by LAPACK
		integer :: info												! status variable used by LAPACK
		real(kind=dbl), dimension(size(A)) :: work					! work vector used by LAPACK (blocksize!!) work(1) returns optimal lwork value
		integer :: lwork											! size of work vector used by LAPACK (should be n*optimal blocksize)
		! --------------------------------------------------------------------------

		! test if A is square
		if (size(A,1) /= size(A,2)) then

			! A is not square
			!write(*,*) 'Matrix is not square, cannot invert!'
			status = 10
			return

		end if

		! initialize A_inv with A
		A_inv = A

		! compute LU factorization of A
		! call to LAPACK subroutine dgetrf
		call dgetrf( size(A,1), size(A,2), &
						A_inv, size(A,1), &
						ipiv, info )

		! info variable check
		if (info /= 0) then

			! error in LU factorization
			!write(*,*) 'Error in LU factorization'
			status = -200
			return

		end if

		! compute inverse of A using the LU factorization
		! call to LAPACK subroutine dgetri
		lwork = size(A)
		call dgetri( size(A,1), &
						A_inv, size(A,1), &
						ipiv, work, lwork, &
						info)

		! info variable check
		if (info /= 0) then

			! error in inversion of A
			!write(*,*) 'Error in matrix inversion'
			status = -300
			return

		end if

		! everything was ok
		! A_inv was computed correctly
		! set status variable to zero
		status = 0

	end subroutine mat_inverse


	function comp_act_yf_ind(n_act_yf, yield_flags)
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		integer, intent(in) :: n_act_yf						! number of active yield functions
		logical, dimension(:), intent(in) :: yield_flags	! active yield functions vector

		! output variable
		integer, dimension(n_act_yf) ::		comp_act_yf_ind		! vector with indices of active yield functions, i.e. indices of entries of the yield flag vector which are false

		! internal variables
		integer :: i 		! counter variable
		integer :: j		! counter variable
		! --------------------------------------------------------------------------

		! compute indices of active yield functions and store in act_yf_ind vector
		j=1
		do i=1,size(yield_flags)
			if (yield_flags(i) .eqv. .false.) then
				comp_act_yf_ind(j) = i
				j = j+1
			end if
		end do

	end function comp_act_yf_ind

end module return_mapping
