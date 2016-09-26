module residuals
    ! contains functions for the computation of the residual vector
    ! which is needed for the newton iteration in the stress update algorithm
    ! (return mapping method) and additionally can be used for the numerical computation
    ! of the jacobian matrix needed in the newton iteration

    ! load additional modules
    use constants
	use material_info
	use aux_routines, only: unpack_x_vector
    use model_comp_act, only: Dg_actDsigma, &
    									hf_act, &
    									yf_act

    implicit none

contains

	function comp_R(x, sigma_trial, alpha_old, &
							C_el, &
							n_act_yf, yield_flags, act_yf_ind)
		! assembles the residual vector R
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	x				! x vector, containing stress vector, alpha vector and active delta_lambda vector
		real(kind=dbl), dimension(:), intent(in) ::		sigma_trial		! trial stress vector
		real(kind=dbl), dimension(:), intent(in) ::		alpha_old		! internal variables vector at the end of the prior load increment
		real(kind=dbl), dimension(:,:), intent(in) ::	C_el			! elasticity matrix
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(x)) ::			comp_R		! residual vector R

		! internal variables
		real(kind=dbl), dimension(6) :: 				sigma				! stress vector
		real(kind=dbl), dimension(n_int_vars) :: 		alpha				! internal variables vector
		real(kind=dbl), dimension(n_act_yf) :: 			delta_lambda_act	! active delta_lambdas vector
		integer :: i_start
		integer :: i_end
		! --------------------------------------------------------------------------

		! unpack x vector
		call unpack_x_vector(x, sigma, alpha, delta_lambda_act)

		! assemble residual vector
		i_start = 1
		i_end = 6
		comp_R(i_start:i_end) = R_sigma(sigma, alpha, delta_lambda_act, &
												sigma_trial, &
												C_el, &
												n_act_yf, yield_flags, act_yf_ind)

		i_start = i_end+1
		i_end = i_start+n_int_vars-1
		if (hardening_present) then
			comp_R(i_start:i_end) = R_alpha(sigma, alpha, delta_lambda_act, &
											alpha_old, &
											n_act_yf, yield_flags, act_yf_ind)
		end if


		i_start = i_end+1
		i_end = i_start+n_act_yf-1
		comp_R(i_start:i_end) = R_f(sigma, alpha, &
									n_act_yf, yield_flags, act_yf_ind)

	end function comp_R




	function R_sigma(sigma, alpha, delta_lambda_act, &
									sigma_trial, &
									C_el, &
									n_act_yf, yield_flags, act_yf_ind)
		! computes R_sigma part of the residual vector
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		real(kind=dbl), dimension(:), intent(in) ::		delta_lambda_act	! active delta lambda vector
		real(kind=dbl), dimension(:), intent(in) ::		sigma_trial		! trial stress vector
		real(kind=dbl), dimension(:,:), intent(in) :: 	C_el			! elasticity matrix
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(6) :: 	R_sigma

		! internal variables
		real(kind=dbl), dimension(6,n_act_yf) :: dg_act		! matrix with dg_act/dsigma vectors
		real(kind=dbl), dimension(6) ::		A				! auxiliary variable (vector) for sum( delta_lambda_act*dg_act/dsigma)
		integer :: i													! index variable
		! --------------------------------------------------------------------------

		! compute dg_act/dsigma matrix
		dg_act = Dg_actDsigma(sigma, alpha, yield_flags, act_yf_ind)

		! initialize A vector with zeros
		A = 0._dbl

		! loop over all active delta lambdas
		do i=1,n_act_yf

			! compute sum( delta_lambda_act*dg_act/dsigma)
			A = A + delta_lambda_act(i) * dg_act(:,i)

		end do

		! compute R_sigma
		R_sigma = -1._dbl*sigma + sigma_trial - matmul(C_el,A)

	end function R_sigma




	function R_alpha(sigma, alpha, delta_lambda_act, &
									alpha_old, &
									n_act_yf, yield_flags, act_yf_ind)
		! compute R_alpha part of the residual vector
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		real(kind=dbl), dimension(:), intent(in) ::		delta_lambda_act	! active delta lambda vector
		real(kind=dbl), dimension(:), intent(in) ::		alpha_old		! internal variables vector at the end of the prior load increment
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(alpha)) ::		R_alpha

		! internal variables
		real(kind=dbl), dimension(size(alpha),n_act_yf) :: h_act		! matrix with active hardening function vectors
		real(kind=dbl), dimension(size(alpha)) ::		A				! auxiliary variable (vector) for sum( delta_lambda_act*h_act)
		integer :: i													! index variable
		! --------------------------------------------------------------------------

		! compute active hardening functions matrix
		h_act = hf_act(sigma, alpha, yield_flags, act_yf_ind)

		! initialize A vector with zeros
		A = 0._dbl

		! loop over all active delta lambdas
		do i=1,n_act_yf

			! compute sum( delta_lambda_act*h_act)
			A = A + delta_lambda_act(i)*h_act(:,i)

		end do

		! compute R_alpha
		R_alpha = -1._dbl*alpha + alpha_old + A


	end function R_alpha




	function R_f(sigma, alpha, &
					n_act_yf, yield_flags, act_yf_ind)
		! computes R_f part of the residual vector
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(n_act_yf) ::			R_f

		! internal variables
		! --------------------------------------------------------------------------

		! compute R_f
		R_f = yf_act(sigma, alpha, yield_flags, act_yf_ind)

	end function R_f

end module residuals
