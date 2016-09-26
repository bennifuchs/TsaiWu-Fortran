module multsurf_check
	! contains functions which call the multisurface returnmapping algorithm
	! and perform plausibility checks
	! with the results of the multisurface returnmapping algorithm (sigma, delta_lambda_active)

	! load additional modules
	use constants
	use material_info
	use return_mapping, only: multsurf_return_mapping, &
								comp_act_yf_ind
	use yf_all, only: yf
	use model_comp_act, only: Dg_actDsigma

    implicit none

contains

	subroutine multsurf_str_upd(sigma_trial, alpha_old, &
    										sigma_new, alpha_new, delta_lambda, &
    										delta_eps_pl, &
    										C_ep, &
    										C_el, &
											n_act_yf, yield_flags, &
											n_yf_c_tested, n_newton_it, &
											status)
		! performs multisurface return mapping step by
		! calling the multi_surface_return_mapping routine in the return_mapping module
		! and performs plausibility checks with the results of the multisurface return mapping routine
		!
		!If neccessary it starts a new return mapping step with a new set of yield functions
		! chosen based on the results of the plausibility checks
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: sigma_trial		! trial stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha_old		! internal variables vector of the prior step
		real(kind=dbl), dimension(:,:), intent(in) :: C_el			! elasticity matrix

		! return variables
		real(kind=dbl), dimension(:), intent(out) :: sigma_new				! updated stress vector
		real(kind=dbl), dimension(:), intent(out) :: alpha_new		! updated internal variables vector
		real(kind=dbl), dimension(:,:), intent(out) :: C_ep					! elastic-plastic tangent matrix
		integer, intent(inout) :: n_act_yf									! number of currently active yield functions
		logical, dimension(:), intent(inout) :: yield_flags				! active yield functions vector
		real(kind=dbl), dimension(:), intent(out) :: delta_lambda		! vector with plastic multipliers of all yield functions
		real(kind=dbl), dimension(:), intent(out) :: delta_eps_pl			! plastic strain increment vector
		integer, intent(out) :: n_yf_c_tested								! number of tested yield flags combinations
		integer, intent(out) :: n_newton_it									! number of newton iterations performed
		integer, intent(out) :: status										! status variable

		! internal variables
		real(kind=dbl), allocatable, dimension(:) :: delta_lambda_act		! vector with plastic multipliers of active yield functions
																			! varies in size, depending on the number of active yield functions
																			! therefore defined as allocatable
		logical, allocatable, dimension(:,:) :: yield_flags_tested			! array for the storage of already used yield flags vectors which are newly created
																			! after each failed plausibility check
																			! the newly generated yield flags vector is compared to the already used vectors
																			! and if it turns out that the vector was already used before,
																			! then search for a new yield surface combination is aborted and the stress update of increment failed
		integer :: n_max_yield_flags										! number of maximum unique yield_flag combinations
		integer :: k														! counter for the number of newton iterations performed in the return mapping step
		integer :: m														! counter for the number of already tested yield flags combinations
		logical :: plaus_status												! status variable which defines if plausibility check
																			! for results of returnmapping step was succesfull (true) or not (false)
		logical, allocatable, dimension(:) :: delta_lambda_act_bool			! vector with boolean values which define
																			! if the active plastic multipliers are >= 0 (True) or < 0 (False)
		logical, dimension(n_yf) :: yf_all_bool								! vector with boolean values which define
																			! if the yield function values are <= 0 (True) or > 0 (False)
		integer :: i														! index variable
		integer :: j														! second index variable
		integer :: n_act_yf_tmp												! variable for temporary storage of the updated number of active yield functions
		! --------------------------------------------------------------------------

		! compute number of maximum unique yield_flag combinations
		! done via the computation of number of neccesary digits for displaying a number
		! in a base-2 numeral system (i.e. binary system) (e.g. if there are four yield functions, yield_flags
		! is a n=4 components vector and the maximum number of unique active yield functions combinations
		! is given by 2**n-1 = 2**4-1 = 15
		! this is the number of columns needed for the yield_flags_tested array
		n_max_yield_flags = 2**n_yf - 1

		! allocate yield_flags_tested array
		! (lines: number of yield functions, columns: number of maximum unique yield_flag combinations)
		allocate(yield_flags_tested(n_yf, n_max_yield_flags))

		! initialize first column of yield_flags_tested array
		! with initial input yield_flags vector
		yield_flags_tested = .True.
		yield_flags_tested(:,1) = yield_flags

		allocate(delta_lambda_act(n_act_yf))				! allocate delta_lambda_act vector
		allocate(delta_lambda_act_bool(n_act_yf))			! and delta_lambda_act_bool vector, (size: n_act_yf)

		m = 0			! set counter for already tested yield flags combinations to zero

!		write(*,*) ''
!		write(*,*) 'start of multi surface loop: '
		! start multisurface stress update loop
		multsurf_loop: do

!			write(*,*) 'multi surface loop no: ', m
!			write(*,*) 'tested yield_flags: '
!			call write_bool_list(yield_flags_tested, n_yf, n_max_yield_flags)
!			write(*,*) yield_flags_tested

			! perform return mapping step with the initial set of active yield functions
			! defined by the initial input yield_flags vector
			call multsurf_return_mapping(sigma_trial, alpha_old, &
    										sigma_new, alpha_new, delta_lambda_act, &
    										C_ep, &
    										C_el, &
											n_act_yf, yield_flags, &
											k, status)

!			if (status /= 0) then
!				write(*,*) 'status after return mapping: ', status
!			end if

			! check for errors in return mapping step
			ret_mapping_stat_check: if (status /= 0) then

				! return mapping step didn' converge for the initial yield flags set
				! cycle through all yield functions and set the current yield function active (false)
				! while leaving all others inactive (true)
				! and try the return mapping step with this yield_flags combination

				! set number of active yield funcions to one
				n_act_yf = 1

				! adjust size of delta_lambda_act vector
				! and delta_lambda_act_bool vector
				deallocate(delta_lambda_act)
				allocate(delta_lambda_act(n_act_yf))

!				write(*,*) 'Start of first single yield function test loop:'

				! loop through each yield function an try return mapping
				call single_yf_test_loop(sigma_trial, alpha_old, &
											sigma_new, alpha_new, delta_lambda_act, &
											C_el, &
											C_ep, &
											yield_flags, &
											k, status)

!				write(*,*) 'single yield function test loop performed'

				if (status == 0) then
					exit multsurf_loop		! single yield function loop was successful, exit multisurface loop
				end if

				! single yield function loop was not successful
				! start yield function pair loop
				n_act_yf = 2		! set number of active yield funcions to two

				! adjust size of delta_lambda_act vector
				deallocate(delta_lambda_act)
				allocate(delta_lambda_act(n_act_yf))

!				write(*,*) 'Start of first yield function pair test loop: '

				! loop through all independent yield function pairs an try return mapping
				call yf_pair_test_loop(sigma_trial, alpha_old, &
											sigma_new, alpha_new, delta_lambda_act, &
											C_el, &
											C_ep, &
											yield_flags, &
											k, status)

!				write(*,*) 'yf_pair_test_loop performed'

				if (status == 0) then
					exit multsurf_loop		! yield function pair loop was successful, exit multisurface loop
				else
					deallocate(yield_flags_tested, &
								delta_lambda_act, delta_lambda_act_bool)	! some error occured, deallocate arrays

					if (status == -30) then			! check if failure in yf pair loop was - no addmissible active yield function pair found
						status = -40				! set error number for error in first yf pair loop
					end if

					return

				end if

			end if ret_mapping_stat_check

			! call routine for plausibility check of delta_lambda_act vector and active yield function vector
			call check_plaus(sigma_new, alpha_new, delta_lambda_act, &
										plaus_status, delta_lambda_act_bool, yf_all_bool, k)

			! update yield flag tested counter
			m = m+1

			! check results of plausibility check
			if (plaus_status .eqv. .true.) then

				! everything ok, exit multsurf loop
				exit multsurf_loop

			else

				! plausibility check failed
	            ! compute new yield_flags vector based on following rules:
	            !     if a delta_lambda entry is < 0-numerical tolerance (i.e boolean value is False)
	            !          => set corresponding yield function inactive by setting the yield_flag entry to True
	            !     if a yield function entry is > 0+numerical tolerance (i.e boolean value is False)
	            !          => set corresponding yield function active by setting the yield_flag entry to False
				yield_flags = comp_new_y_flags(yield_flags, delta_lambda_act_bool, yf_all_bool)

			end if

			! check if the new yield flag vector was already used before
			! loop through all already used yield flags vectors in yield_flags_tested array
			yield_flags_check_loop: do i=1,m

				! check if the i_th yield flag vector equals the new yield flag vector
				yf_comb_already_used: if (all(yield_flags_tested(:,i) .eqv. yield_flags)) then

					! new yield flag vector was already used before, without success
					! cycle through all yield functions and set the current yield function active (false)
					! while leaving all others inactive (true)
					! and try the return mapping step with this yield_flags combination

					! set number of active yield funcions to one
					n_act_yf = 1

					! adjust size of delta_lambda_act vector
					! and delta_lambda_act_bool vector
					deallocate(delta_lambda_act)
					allocate(delta_lambda_act(n_act_yf))

					! loop through each yield function an try return mapping
					call single_yf_test_loop(sigma_trial, alpha_old, &
												sigma_new, alpha_new, delta_lambda_act, &
												C_el, &
												C_ep, &
												yield_flags, &
												k, status)

					! check status of single yield function loop
					if (status /= 0) then
						deallocate(yield_flags_tested, &
									delta_lambda_act, delta_lambda_act_bool)	! single yield function loop failed, deallocate arrays and return
						return
					else
						exit multsurf_loop			! single yield function loop was successful, exit multisurface loop
					end if

				end if yf_comb_already_used

			end do yield_flags_check_loop

			! new yield flag vector was not used before
			! add new yield flags vector to the yield_flags_tested array
			yield_flags_tested(:,m+1) = yield_flags

			! compute the updated number of active yield functions
			! and store it in a temporary variable
			n_act_yf_tmp = n_yf - count(yield_flags)

			! check if size of the delta_lambda_act and delta_lambda_act_bool
			! vectors need to be adjusted
			if (n_act_yf_tmp /= n_act_yf) then

				! new and old sizes do not match
				! adjust size of delta_lambda_act vector
				! and delta_lambda_act_bool vector
				deallocate(delta_lambda_act)
				allocate(delta_lambda_act(n_act_yf_tmp))
				deallocate(delta_lambda_act_bool)
				allocate(delta_lambda_act_bool(n_act_yf_tmp))

			end if

			! update old size value with the new one
			n_act_yf = n_act_yf_tmp

			! return to the top of the loop and start a new return mapping step
			! with the newly generated yield_flags vector

		end do multsurf_loop

		! compute plastic epsilon increment
		delta_eps_pl = comp_delta_eps_pl(sigma_new, alpha_new, delta_lambda_act, &
											n_act_yf, yield_flags)

		! copy entries of delta_lambda_act vector to the right positions
		! of the delta_lambda vector with the help of the yield_flags vector
		delta_lambda = arrange_delta_lambda(delta_lambda_act, yield_flags)

		! save number of tested yield flags combinations and number of newton iterations
		n_yf_c_tested = m
		n_newton_it = k

		! set status variable to zero, deallocate arrays and return
		status = 0

		! deallocate arrays
		deallocate(yield_flags_tested, &
					delta_lambda_act, delta_lambda_act_bool)

		return

	end subroutine multsurf_str_upd




	subroutine check_plaus(sigma_tmp, alpha_tmp, delta_lambda_act_tmp, &
									plaus_status, delta_lambda_act_tmp_bool, yf_all_tmp_bool, &
									k)
		! checks if results from returnmapping step sigma_tmp, alpha_tmp and delta_lambda_act_tmp
		! fullfill following reqirements:
		! 1) all active delta_lambdas >= 0
		! 2) all yield function values <= 0
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: sigma_tmp						! stress vector from return mapping step
		real(kind=dbl), dimension(:), intent(in) :: alpha_tmp						! internal variables vector from return mapping step
		real(kind=dbl), dimension(:), intent(in) :: delta_lambda_act_tmp			! active plastic multipliers vector from return mapping step
		integer, intent(in) :: k													! number of newton iterations performed in the return mapping step

		! return variables
		logical, intent(out) :: plaus_status										! status variable of plausibility check, check ok = True, check failed = False
		logical, dimension(size(delta_lambda_act_tmp)), intent(out) :: delta_lambda_act_tmp_bool				! vector with boolean values which define
																					! if the active plastic multipliers are >= 0 (True) or < 0 (False)
		logical, dimension(n_yf), intent(out) :: yf_all_tmp_bool					! vector with boolean values which define
																					! if the yield function values are <= 0 (True) or > 0 (False)

		! internal variables
		real(kind=dbl), dimension(n_yf) :: yf_all_tmp
		real(kind=dbl) :: num_tol													! numerical tolerance for >= 0 and <= 0 tests
		! --------------------------------------------------------------------------

		! check the number k of newton iterations performed in the return mapping step
		! and based on this value adjust the numerical tolerance num_tol
		if (k < k_alt) then
			num_tol = yf_tol_def
		else
			num_tol = yf_tol_alt
		end if

		! compute all yield functions with sigma, alpha values from return mapping step
		yf_all_tmp = yf(sigma_tmp, alpha_tmp)

!		write(*,*) '		plausibility check: '
!		write(*,*) '		----------------------'
!		write(*,*) '		num_tol: ', num_tol
!		write(*,*) '		sigma: ', sigma_tmp
!		write(*,*) '		alpha: ', alpha_tmp
!		write(*,*) '		yf_all_tmp: ', yf_all_tmp
!		write(*,*) '		delta_lambda_act_tmp: ', delta_lambda_act_tmp
!		write(*,*) '		----------------------'

		! test if yield function values are <= 0 (True) or > 0 (False)
		! and assign boolean result values to yf_all_tmp_bool vector
		where (yf_all_tmp <= (0._dbl+num_tol))
			yf_all_tmp_bool = .true.
		elsewhere
			yf_all_tmp_bool = .false.
		end where

		! test if active delta_lambda values are >= 0 (True) or < 0 (False)
		! and assign boolean result values to delta_lambda_act_tmp_bool vector
!		where (delta_lambda_act_tmp >= (0._dbl-num_tol))
		where (.not.(delta_lambda_act_tmp < 0._dbl) )
			delta_lambda_act_tmp_bool = .true.
		elsewhere
			delta_lambda_act_tmp_bool = .false.
		end where

		! check if all boolean delta_lambda_act values
		! AND all boolean yield function values are true
		if (all(mask=yf_all_tmp_bool) .and. all(mask=delta_lambda_act_tmp_bool)) then

			! all values are true => plausibility check was successfull
			! set status to True
			plaus_status = .true.

		else

			! some values of delta_lambda_act_tmp_bool and/or
			! some values of yf_all_tmp_bool were false
			! => plausibility check failed, set status to False
			plaus_status = .false.

		end if

	end subroutine check_plaus




	function comp_new_y_flags(yield_flags_tmp, delta_lambda_act_tmp_bool, yf_all_tmp_bool)
        ! compute new yield_flags vector based on following rules:
        !     if a delta_lambda entry is < 0-numerical tolerance (i.e boolean value is False)
        !          => set corresponding yield function inactive by setting the yield_flag entry to True
        !     if a yield function entry is > 0+numerical tolerance (i.e boolean value is False)
        !          => set corresponding yield function active by setting the yield_flag entry to False
        implicit none

		! --------------------------------------------------------------------------
		! passed variables
		logical, dimension(:), intent(in) :: yield_flags_tmp					! yield flags vector which indicates which yield function is active (False) and which is inactive (True)
		logical, dimension(:), intent(in) :: delta_lambda_act_tmp_bool			! vector with boolean values which define
																				! if the active plastic multipliers are >= 0 (True) or < 0 (False)
		logical, dimension(:), intent(in) :: yf_all_tmp_bool					! vector with boolean values which define
																				! if the yield function values are <= 0 (True) or > 0 (False)

		! return variable
		logical, dimension(size(yield_flags_tmp)) :: comp_new_y_flags

		! internal variables
		integer :: i				! index for entries of yield_flags_tmp vector and boolean yield function vector
		integer :: j				! index for entries of delta_lambda_act_tmp_bool
		! --------------------------------------------------------------------------

		! initialize delta_lambda_act_tmp_bool index with one
		j = 1
		! loop through yield flags vector entries
		change_flags: do i=1,size(yield_flags_tmp)

			! check if i_th yield flag value is true, i.e. the yield function is initially inactive
			if (yield_flags_tmp(i) .eqv. .true.) then

				! yield function is inactive
				! set new yield flag for yield function, depending on the value of yf_all_tmp_bool
				! if yf_all_tmp_bool is false (yf > 0) => set yield function active, i.e yield_flag to false
				! if yf_all_tmp_bool is true (yf <= 0) => leave yield function inactive, i.e yield_flag to true
				comp_new_y_flags(i) = yf_all_tmp_bool(i)

			elseif (yield_flags_tmp(i) .eqv. .false.) then
				! yield function is active
				! set new yield flag for yield function, depending on the value of delta_lambda_act_tmp_bool
				! if delta_lambda_act_tmp_bool is true (delta_lambda >= 0) => leave yield function active, i.e yield_flag to false
				! if delta_lambda_act_tmp_bool is false (delta_lambda < 0) => set yield funtion inactive, i.e yield_flag to true
				comp_new_y_flags(i) = .not.delta_lambda_act_tmp_bool(j)
				j = j+1

			end if

		end do change_flags

	end function comp_new_y_flags




	function comp_delta_eps_pl(sigma, alpha, delta_lambda_act, &
								n_act_yf, yield_flags)
		! computes the plastic strain increment vector
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha		! internal variables vector
		real(kind=dbl), dimension(:), intent(in) :: delta_lambda_act	! active delta lambdas vector
		integer, intent(in) :: n_act_yf							! number of active yield functions
		logical, dimension(:), intent(in) :: yield_flags		! yield flags vector

		! return variable
		real(kind=dbl), dimension(size(sigma)) :: comp_delta_eps_pl		! plastic strain increment vector

		! internal variables
		real(kind=dbl), dimension(size(sigma), n_act_yf) :: Dg_act		! dg_act/dsigma matrix
		integer :: i											! counter

		integer, dimension(n_act_yf) ::		act_yf_ind		! vector with indices of active yield functions, i.e. indices of entries of the yield flag vector which are false
		! --------------------------------------------------------------------------

		act_yf_ind = comp_act_yf_ind(n_act_yf, yield_flags)

		! initialize delta_eps_pl with zero
		comp_delta_eps_pl = 0._dbl

		! compute the dg_act/dsigma matrix
		! and transpose it (dg_act_i/dsigma in columns)
		Dg_act = Dg_actDsigma(sigma, alpha, yield_flags, act_yf_ind)

		! loop over all active yield functions
		do i=1,n_act_yf

			! sum of (delta_lambda_act_i * dg_act_i/dsigma)
			comp_delta_eps_pl = comp_delta_eps_pl + delta_lambda_act(i)*Dg_act(:,i)

		end do

	end function comp_delta_eps_pl




	function arrange_delta_lambda(delta_lambda_act, yield_flags)
		! copy entries of delta_lambda_act vector to the right positions
		! of the delta_lambda vector with the help of the yield_flags vector
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: delta_lambda_act			! vector with plastic multpliers of active yield functions
		logical, dimension(:), intent(in) :: yield_flags						! yield flags vector which indicates which yield function is active and which is not

		! return variable
		real(kind=dbl), dimension(size(yield_flags)) :: arrange_delta_lambda	! vector with plastic multpliers of all yield functions
																				! for inactive yield functions the values are set to zero

		! internal variables
		integer :: i			! index variable for all delta lambdas
		integer :: j			! index variable for only active delta_lambdas
		! --------------------------------------------------------------------------

		! initialize delta_lambda vector with zeros
		arrange_delta_lambda = 0._dbl

		! initialize index variable for active delta_lambdas with one
		j = 1
		! loop through yield flags vector
		do i=1,size(yield_flags)

			! check i-th yield flag
			if (yield_flags(i) .eqv. .false.) then

				! yield function is active
				arrange_delta_lambda(i) = delta_lambda_act(j)
				! update index j
				j = j+1

			end if

		end do

	end function arrange_delta_lambda




	subroutine single_yf_test_loop(sigma_trial, alpha_old, &
									sigma_new, alpha_new, delta_lambda_act, &
									C_el, &
									C_ep, &
									yield_flags, &
									k, status)
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) ::						sigma_trial		! trial stress vector
		real(kind=dbl), dimension(:), intent(in) ::						alpha_old		! internal variables vector of the prior step
		real(kind=dbl), dimension(:,:), intent(in) ::					C_el			! elasticity matrix

		! return variables
		real(kind=dbl), dimension(6), intent(out) ::			sigma_new		! updated stress vector
		real(kind=dbl), dimension(n_int_vars), intent(out) ::	alpha_new		! updated internal variables vector
		real(kind=dbl), dimension(1), intent(out) :: 			delta_lambda_act! incremental plastic multiplier vector
		real(kind=dbl), dimension(6,6), intent(out) :: 			C_ep			! elastic-plastic tangent matrix
		logical, dimension(n_yf), intent(out) :: 				yield_flags		! active yield functions vector
		integer, intent(out) ::									status			! status variable
		integer, intent(out) ::									k				! counter for the number of newton iterations performed in the return mapping step


		! internal variables
		logical, dimension(1) :: 			delta_lambda_act_bool	! vector with boolean values which define
																	! if the active plastic multipliers are >= 0 (True) or < 0 (False)
		logical ::							plaus_status		! status variable which defines if plausibility check
		integer ::							n_act_yf			! number of currently active yield functions
		integer :: 							j					! number of tested yield functions
		logical, dimension(n_yf) ::			yf_all_bool			! vector with boolean values which define
		! --------------------------------------------------------------------------

		! set number of active yield funcions to one
		n_act_yf = 1

		single_yf_try_loop: do j=1,n_yf

			! set j-th yield function active an all others inactive
			yield_flags = .true.
			yield_flags(j) = .false.

!			write(*,*) '	single yf try loop no: ', j
!			write(*,*) '    single yf yield_flags: ', yield_flags

			! perform return mapping step with this yield flags vector
			call multsurf_return_mapping(sigma_trial, alpha_old, &
	    										sigma_new, alpha_new, delta_lambda_act, &
	    										C_ep, &
	    										C_el, &
												n_act_yf, yield_flags, &
												k, status)

!			write(*,*) '    status after return mapping: ', status

			! check for errors in return mapping step
			if (status /= 0) then
				cycle single_yf_try_loop	! errors occured, try with another active yield function
			end if

!			write(*,*) '	plausibility check entered'

			! perform plausibility check for the yield function values and delta_lambda_act values
			! obtained from the return mapping step
			call check_plaus(sigma_new, alpha_new, delta_lambda_act, &
										plaus_status, delta_lambda_act_bool, yf_all_bool, k)

!			write(*,*) '	plausibility check left'

			! check result of plausibility check
			if (plaus_status .eqv. .true.) then
				return		! everything ok, exit subroutine
			end if

		end do single_yf_try_loop

		! check if UMAT failed or single_yf_try_loop
		if (status /= 0) then
			return										! Umat failed
		else if (plaus_status .eqv. .false.) then
			status = -20								! no addmissible active yield function found
			return
		end if

	end subroutine single_yf_test_loop



	subroutine yf_pair_test_loop(sigma_trial, alpha_old, &
									sigma_new, alpha_new, delta_lambda_act, &
									C_el, &
									C_ep, &
									yield_flags, &
									k, status)
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) ::						sigma_trial		! trial stress vector
		real(kind=dbl), dimension(:), intent(in) ::						alpha_old		! internal variables vector of the prior step
		real(kind=dbl), dimension(:,:), intent(in) ::					C_el			! elasticity matrix

		! return variables
		real(kind=dbl), dimension(6), intent(out) ::			sigma_new		! updated stress vector
		real(kind=dbl), dimension(n_int_vars), intent(out) ::	alpha_new		! updated internal variables vector
		real(kind=dbl), dimension(2), intent(out) :: 		delta_lambda_act	! incremental plastic multiplier vector, length is two for two assumed active yfs
		real(kind=dbl), dimension(6,6), intent(out) :: 			C_ep			! elastic-plastic tangent matrix
		logical, dimension(n_yf), intent(out) :: 				yield_flags		! active yield functions vector
		integer, intent(out) ::									status			! status variable
		integer, intent(out) ::									k				! counter for the number of newton iterations performed in the return mapping step


		! internal variables
		logical, dimension(2) :: 		delta_lambda_act_bool	! vector with boolean values which define
																	! if the active plastic multipliers are >= 0 (True) or < 0 (False)
		logical ::							plaus_status		! status variable which defines if plausibility check
		integer ::							n_act_yf			! number of currently active yield functions
		integer :: 							j					! counter for first active yield function of test pair
		integer ::							i					! counter for second active yield function of test pair
		logical, dimension(n_yf) ::			yf_all_bool			! vector with boolean values which define
		integer :: l
		! --------------------------------------------------------------------------

		! set number of active yield funcions to one
		n_act_yf = 2

		l=1
		outer_yf_pair_try_loop: do j=1,n_yf-1		! loop through all yield functions, except the last one

			! set j-th yield function active and all others inactive
			yield_flags = .true.
			yield_flags(j) = .false.

			inner_yf_pair_try_loop: do i=j+1,n_yf		! loop through remaining yfs, starting from the (j+1)-th one

				yield_flags(i) = .false.		! set i-th yield function as second one of the pair active

!				write(*,*) '	yf pair try loop no: ', l
!				write(*,*) '    yf pair yield_flags: ', yield_flags

				! perform return mapping step with this yield flags vector
				call multsurf_return_mapping(sigma_trial, alpha_old, &
		    										sigma_new, alpha_new, delta_lambda_act, &
		    										C_ep, &
		    										C_el, &
													n_act_yf, yield_flags, &
													k, status)

!				write(*,*) '    status after return mapping: ', status

				! check for errors in return mapping step
				if (status /= 0) then
					cycle inner_yf_pair_try_loop	! errors occured, try with another active yield function
				end if

				! perform plausibility check for the yield function values and delta_lambda_act values
				! obtained from the return mapping step
				call check_plaus(sigma_new, alpha_new, delta_lambda_act, &
											plaus_status, delta_lambda_act_bool, yf_all_bool, k)

				! check result of plausibility check
				if (plaus_status .eqv. .true.) then
					return		! everything ok, exit subroutine
				end if
				l = l+1
			end do inner_yf_pair_try_loop

		end do outer_yf_pair_try_loop

		! check if UMAT failed or single_yf_try_loop
		if (status /= 0) then
			return										! Umat failed
		else if (plaus_status .eqv. .false.) then
			status = -30								! no addmissible active yield function pair found
			return
		end if

	end subroutine yf_pair_test_loop

end module multsurf_check
