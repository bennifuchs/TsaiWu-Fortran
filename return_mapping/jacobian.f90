module jacobian
    ! module that contains functions for the
    ! computation of the jacobian matrix based on the active yield functions and active plastic potentials
    ! contains one main function that assembles the jacobian matrix
    ! and several local functions that compute sub matrices of the jacobian matrix
    !
    ! contains one additional function that computes the jacobian matrix numerically
    ! based on small disturbances of each variable of the residual vector

	! load additional modules
    use constants
	use material_info
	use aux_routines, only: unpack_x_vector
    use model_comp_act, only: DDg_actDDsigma, &
    									DDg_actDsigmaDalpha, &
    									Dg_actDsigma, &
    									Dh_actDsigma, &
    									Dh_actDalpha, &
    									hf_act, &
    									Df_actDsigma, &
    									Df_actDalpha
    use residuals, only: comp_R

    implicit none

contains

	function comp_jacobian(x, C_el, &
							n_act_yf, yield_flags, act_yf_ind)
		! assembles jacobian matrix, based on active yield functions and plastic potential
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	x				! x vector, containing stress vector, alpha vector and active delta_lambda vector
		real(kind=dbl), dimension(:,:), intent(in) ::	C_el			! elasticity matrix
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(x),size(x)) :: comp_jacobian

		! internal variables
		real(kind=dbl), dimension(6) :: 				sigma				! stress vector
		real(kind=dbl), dimension(n_int_vars) :: 		alpha				! internal variables vector
		real(kind=dbl), dimension(n_act_yf) :: 			delta_lambda_act	! active delta_lambdas vector
		! --------------------------------------------------------------------------

		! unpack x vector
		call unpack_x_vector(x, sigma, alpha, delta_lambda_act)

		! initialization of jacobian with zeros
		!construct_jacobian = 0._dbl


		! assign dR_sigma/dsigma to jacobian
		comp_jacobian(1:6, 1:6) = DR_sigmaDsigma(sigma, alpha, delta_lambda_act, &
														C_el, &
														n_act_yf, yield_flags, act_yf_ind)
		! assign dR_sigma/dalpha to jacobian
		if (hardening_present) then
			comp_jacobian(1:6, 6+1:6+n_int_vars) = DR_sigmaDalpha(sigma, alpha, delta_lambda_act, &
																		C_el, &
																		n_act_yf, yield_flags, act_yf_ind)
		end if

		! assign dR_sigma/ddelta_lambda to jacobian
		comp_jacobian(1:6, &
							6+n_int_vars+1:6+n_int_vars+n_act_yf) = DR_sigmaDdelta_lambda(sigma, alpha, &
																							C_el, &
																							n_act_yf, yield_flags, act_yf_ind)

		! assign dR_alpha/dsigma to jacobian
		if (hardening_present) then
			comp_jacobian(6+1:6+n_int_vars, &
								1:6) = DR_alphaDsigma(sigma, alpha, delta_lambda_act, &
														n_act_yf, yield_flags, act_yf_ind)
			! assign dR_alpha/dalpha to jacobian
			comp_jacobian(6+1:6+n_int_vars, &
								6+1:6+n_int_vars) = DR_alphaDalpha(sigma, alpha, delta_lambda_act, &
																	n_act_yf, yield_flags, act_yf_ind)
			! assign dR_alpha/ddelta_lambda to jacobian
			comp_jacobian(6+1:6+n_int_vars, &
								6+n_int_vars+1:6+n_int_vars+n_act_yf) = DR_alphaDdelta_lambda(sigma, alpha, &
																								n_act_yf, yield_flags, act_yf_ind)
		end if

		! assign dR_f/dsigma to jacobian
		comp_jacobian(6+n_int_vars+1:6+n_int_vars+n_act_yf, &
							1:6) = DR_fDsigma(sigma, alpha, &
												n_act_yf, yield_flags, act_yf_ind)
		! assign dR_f/dalpha to jacobian
		if (hardening_present) then
			comp_jacobian(6+n_int_vars+1:6+n_int_vars+n_act_yf, &
								6+1:6+n_int_vars) = DR_fDalpha(sigma, alpha, &
																n_act_yf, yield_flags, act_yf_ind)
		end if

		! assign dR_f/ddelta_lambda to jacobian
		comp_jacobian(6+n_int_vars+1:6+n_int_vars+n_act_yf, &
							6+n_int_vars+1:6+n_int_vars+n_act_yf) = DR_fDdelta_lambda(n_act_yf)

	end function comp_jacobian




	function DR_sigmaDsigma(sigma, alpha, delta_lambda_act, &
								C_el, &
								n_act_yf, yield_flags, act_yf_ind)
		! computes dR_sigma/dsigma matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		real(kind=dbl), dimension(:), intent(in) :: 	delta_lambda_act			! active plastic multipliers vector
		real(kind=dbl), dimension(:,:), intent(in) :: 	C_el			! elasticity matrix
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(6,6) :: 				DR_sigmaDsigma	! dR_sigma/dsigma

		! internal variables
		real(kind=dbl), dimension(6,6) ::				ident_m			! identity matrix
		real(kind=dbl), dimension(6,6,n_act_yf) :: 		DDg_act			! array with all ddg_act/ddsigma matrices
		real(kind=dbl), dimension(6,6) ::				A				! auxiliary matrix for sum(delta_lambda * ddg_act/ddsigma)
		integer ::										i				! index variable
		! --------------------------------------------------------------------------

		! compute ddg_act/ddsigma
		DDg_act = DDg_actDDsigma(sigma, alpha, yield_flags, act_yf_ind)

		! initialize A with zeros
		A = 0._dbl

		! loop over all active lambdas and ddg_act/ddsigma matrices
		do i = 1,n_act_yf

			! compute sum(delta_lambda_act_i * ddg_act_i/ddsigma)
			A = A + delta_lambda_act(i)*DDg_act(:,:,i)

		end do

		! initialize identity matrix
		ident_m = 0._dbl
		forall(i=1:6)
			ident_m(i,i) = 1._dbl
		end forall

		! compute dR_sigma/dsigma
		DR_sigmaDsigma = -1._dbl*ident_m - matmul(C_el,A)

	end function DR_sigmaDsigma




	function DR_sigmaDalpha(sigma, alpha, delta_lambda_act, &
								C_el, &
								n_act_yf, yield_flags, act_yf_ind)
		! computes dR_sigma/dalpha matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		real(kind=dbl), dimension(:), intent(in) :: 	delta_lambda_act			! active plastic multipliers vector
		real(kind=dbl), dimension(:,:), intent(in) :: 	C_el			! elasticity matrix
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(6,size(alpha)) ::	DR_sigmaDalpha	! dR_sigma/dalpha

		! internal variables
		real(kind=dbl), dimension(6,size(alpha),n_act_yf) :: DDg_act	! ddg_act/dsigmadalpha matrix
		real(kind=dbl), dimension(6,size(alpha)) ::	A				! auxiliary variable for sum(delta_lambda_act_i * ddg_act_i/dsigmadalpha)
		integer ::										i				! index variable
		! --------------------------------------------------------------------------

		! computation ddg_act/dsigmadalpha
		DDg_act = DDg_actDsigmaDalpha(sigma, alpha, yield_flags, act_yf_ind)

		! initialization of A with zeros
		A = 0._dbl

		! loop over all active delta_lambdas
		do i=1,n_act_yf

			! compute sum(delta_lambda_act_i * ddg_act_i/dsigmadalpha)
			A = A + delta_lambda_act(i)*DDg_act(:,:,i)

		end do

		! compute dR_sigma/dalpha
		DR_sigmaDalpha = matmul(-1._dbl*C_el,A)

	end function DR_sigmaDalpha




	function DR_sigmaDdelta_lambda(sigma, alpha, &
									C_el, &
									n_act_yf, yield_flags, act_yf_ind)
		! computes dR_sigma/ddelta_lambda matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		real(kind=dbl), dimension(:,:), intent(in) :: 	C_el			! elasticity matrix
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(6,n_act_yf) ::		DR_sigmaDdelta_lambda	! dR_sigma/ddelta_lambda

		! internal variables
		real(kind=dbl), dimension(6,n_act_yf) ::		Dg_act			! dg_act/dsigma
		integer ::										i				! index variable
		! --------------------------------------------------------------------------

		! compute dg_act/dsigma and transpose it ( matrix with dg_act/dsigma as column vectors)
		Dg_act = Dg_actDsigma(sigma, alpha, yield_flags, act_yf_ind)

		DR_sigmaDdelta_lambda = matmul(-1._dbl*C_el,Dg_act)

	end function DR_sigmaDdelta_lambda




	function DR_alphaDsigma(sigma, alpha, delta_lambda_act, &
								n_act_yf, yield_flags, act_yf_ind)
		! computes dR_alpha/dsigma matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		real(kind=dbl), dimension(:), intent(in) :: 	delta_lambda_act			! active plastic multipliers vector
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(alpha),6) ::	DR_alphaDsigma	! dR_alpha/dsigma

		! internal variables
		real(kind=dbl), dimension(size(alpha),6,n_act_yf) :: Dh_act	! dh_act/dsigma matrix
		integer ::										i				! index variable
		! --------------------------------------------------------------------------

		! computation of dh_act/dsigma matrix
		Dh_act = Dh_actDsigma(sigma, alpha, yield_flags, act_yf_ind)

		! initialization of dR_alpha/dsigma array with zeros
		DR_alphaDsigma = 0._dbl

		! loop over all active delta lambdas
		do i=1,n_act_yf

			! compute sum( delta_lambda_act * dh_act/dsigma )
			DR_alphaDsigma = DR_alphaDsigma + delta_lambda_act(i)*Dh_act(:,:,i)

		end do

	end function DR_alphaDsigma




	function DR_alphaDalpha(sigma, alpha, delta_lambda_act, &
								n_act_yf, yield_flags, act_yf_ind)
		! computes dR_alpha/dalpha matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		real(kind=dbl), dimension(:), intent(in) :: 	delta_lambda_act			! active plastic multipliers vector
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(alpha),size(alpha)) :: DR_alphaDalpha	! dR_alpha/dalpha

		! internal variables
		real(kind=dbl), dimension(size(alpha),size(alpha),n_act_yf) :: dh_act	! dh_act/dalpha matrix
		real(kind=dbl), dimension(size(alpha),size(alpha)) :: A					! auxiliary variable for sum( delta_lambda_act * dh_act/dalpha )
		real(kind=dbl), dimension(size(alpha),size(alpha)) :: ident_m			! identity matrix
		integer ::										i							! index variable
		! --------------------------------------------------------------------------

		! computation of dh_act/dalpha
		dh_act = Dh_actDalpha(sigma, alpha, yield_flags, act_yf_ind)

		! initialization of A with zeros
		A = 0._dbl

		! loop over all active delta lambdas
		do i = 1,n_act_yf

			! compute sum( delta_lambda_act * dh_act/dalpha )
			A = A + delta_lambda_act(i) * dh_act(:,:,i)

		end do

		! initialization of indentity matrix
		ident_m = 0._dbl
		forall(i=1:size(ident_m,1))
			ident_m(i,i) = 1._dbl
		end forall

		! dR_alpha/dalpha
		DR_alphaDalpha = -1._dbl*ident_m + A

	end function DR_alphaDalpha





	function DR_alphaDdelta_lambda(sigma, alpha, &
								n_act_yf, yield_flags, act_yf_ind)
		! computes dR_alpha/ddelta_lambda matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(alpha),n_act_yf) :: DR_alphaDdelta_lambda	! dR_alpha/ddelta_lambda matrix

		! internal variables
		real(kind=dbl), dimension(size(alpha),n_act_yf) :: h_act		! active hardening functions matrix
		! --------------------------------------------------------------------------

		! compute active hardening function matrix
		h_act = hf_act(sigma, alpha, yield_flags, act_yf_ind)

		! dR_alpha/ddelta_lambda
		DR_alphaDdelta_lambda = h_act

	end function DR_alphaDdelta_lambda




	function DR_fDsigma(sigma, alpha, &
								n_act_yf, yield_flags, act_yf_ind)
		! computes dR_f/dsigma matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(n_act_yf,6) ::		DR_fDsigma		! dR_f/dsigma

		! internal variables
		real(kind=dbl), dimension(6,n_act_yf) ::		df_act			! df_act/dsigma matrix
		! --------------------------------------------------------------------------

		! compute df_act/dsigma matrix and transpose it
		df_act = transpose( Df_actDsigma(sigma, alpha, yield_flags, act_yf_ind) )

		! dR_f/dsigma
		DR_fDsigma = df_act

	end function DR_fDsigma




	function DR_fDalpha(sigma, alpha, &
								n_act_yf, yield_flags, act_yf_ind)
		! computes dR_f/dalpha matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		integer, intent(in) :: 							n_act_yf		! number of active yield functions
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(n_act_yf,size(alpha)) :: DR_fDalpha	! DR_f/dalpha

		! internal variable
		real(kind=dbl), dimension(size(alpha),n_act_yf) :: df_act		! df_act/dalpha matrix
		! --------------------------------------------------------------------------

		! compute df_act/dalpha matrix and transpose it
		df_act = transpose( Df_actDalpha(sigma, alpha, yield_flags, act_yf_ind) )

		! dR_f/dalpha
		DR_fDalpha = df_act

	end function DR_fDalpha



	function DR_fDdelta_lambda(n_act_yf)
		! computes dR_f/ddelta_lambda matrix
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		integer, intent(in) :: 							n_act_yf		! number of active yield functions

		! return variable
		real(kind=dbl), dimension(n_act_yf,n_act_yf) :: DR_fDdelta_lambda	! DR_f/ddelta_lambda

		! internal variables
		! --------------------------------------------------------------------------

		! DR_f/ddelta_lmbda
		DR_fDdelta_lambda = 0._dbl

	end function DR_fDdelta_lambda




	function comp_jacobian_num(x, sigma_trial, alpha_old, &
											C_el, &
											n_act_yf, yield_flags, act_yf_ind)
		! assembles numerical jacobian matrix, based on active yield functions and plastic potential
		! by numerical differentiation of the residual vector
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
		real(kind=dbl), dimension(size(x),size(x)) :: comp_jacobian_num

		! internal variables
		real(kind=dbl), dimension(size(x)) :: h						! step width vector
		real(kind=dbl), dimension(size(x)), volatile :: xph						! positively disturbed x vector
		real(kind=dbl), dimension(size(x)), volatile :: xmh						! negatively disturbed x vector
		real(kind=dbl) :: dx_i											! numerically corrected step width
		integer :: i													! index variable
		! --------------------------------------------------------------------------

		! loop over all lines of the x vector
		do i=1,size(x)

			h = 0._dbl				! initialize h vector with zeros

			! compute step width for i-th x entry
			! based on the machine epsilon and the absolute value of the x entry (see wikipedia and fellin_ostermann_2002)
			h(i) = epsilon(x(i))**(1._dbl/3._dbl)*max(abs(x(i)),1._dbl)

			xph = x+h					! compute positively disturbed x vector
			xmh = x-h					! compute negatively disturbed x vector
			dx_i = xph(i) - xmh(i)		! compute numerically corrected step width

			! compute i-th column of jacobian matrix
			! by applying numerical differentiation to th residual vector
			comp_jacobian_num(:,i) = 1._dbl/dx_i * &
											( comp_R(xph, sigma_trial, alpha_old, &
															C_el, &
															n_act_yf, yield_flags, act_yf_ind) - &
											  comp_R(xmh, sigma_trial, alpha_old, &
															C_el, &
															n_act_yf, yield_flags, act_yf_ind) )
		end do

	end function comp_jacobian_num


end module jacobian
