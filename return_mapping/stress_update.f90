module stress_update
	! contains the umat subroutine for the multi surface elastic-plastic damage model for concrete/rock mass
	implicit none

contains

	subroutine umat_tsaiwu(stress, statev, ddsdde, sse, spd, scd, &
	                 rpl, ddsddt, drplde, drpldt, &
	                 stran , dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
	                 ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt, &
	                 celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
		! computes updated stresses based on a multi surface elastic-plastic damage model for concrete/rock mass
		! the return mapping algorithm for the update of the stresses and internal variables is based on
		! an implicit euler integration scheme, i.e. a system of nonlinear equations must be solved with
		! the newton's method

		! load additional modules
		use constants
		use material_info
		use derived_types
		use aux_routines, only: conv_to_3D, &
								unpack_state_vars_new, &
								pack_state_vars_new, &
								conv_to_plane_cond, &
								unpack_material_data, &
								material_setup
		use stress_routines, only: calc_C_elastic, &
									calc_inverse_C_elastic, &
									calc_transv_istr_C_elastic, &
									calc_inv_transv_istr_C_elastic
		use yf_all, only: yf
		use multsurf_check, only: multsurf_str_upd
		use damage, only: comp_damage_components, &
								comp_Cep_damaged

		implicit none

		! ##########################################################################
		! DECLARATION OF VARIABLES
		! ##########################################################################

		! --------------------------------------------------------------
		! passed variables
		! --------------------------------------------------------------
		integer, intent(in) :: ndi			! Number of direct stress components at this point
		integer, intent(in) :: nshr			! Number of engineering shear stress components at this point
		integer, intent(in) :: ntens		! Size of the stress or strain component array (NDI + NSHR)
		integer, intent(in) :: nstatv		! Number of solution-dependent state variables that are associated with this material type
		integer, intent(in) :: nprops		! User-defined number of material constants associated with this user material

		integer, intent(in) :: noel			! element number
		integer, intent(in) :: npt			! integration point number
		integer, intent(in) :: layer		! Layer number (for composite shells and layered solids)
		integer, intent(in) :: kspt			! Section point number within the current layer
		integer, intent(in) :: kstep		! Step number
		integer, intent(in) :: kinc			! increment number

		real(kind=dbl), intent(in), dimension(ntens) :: stran		! An array containing the total strains at the beginning of the increment
!f2py	depend(ntens) :: stran
		real(kind=dbl), intent(in), dimension(ntens) :: dstran		! Array of strain increments
!f2py	depend(ntens) :: dstran
		real(kind=dbl), intent(in), dimension(2) :: time			! time(1): Value of step time at the beginning of the current increment or frequency
																	! time(2): Value of total time at the beginning of the current increment
		real(kind=dbl), intent(in) :: dtime							! time increment
		real(kind=dbl), intent(in) :: temp							! Temperature at the start of the increment
		real(kind=dbl), intent(in) :: dtemp							! Increment of temperature
		real(kind=dbl), intent(in) :: predef						! Array of interpolated values of predefined field variables at this point at the start of the increment, based on the values read in at the nodes
		real(kind=dbl), intent(in) :: dpred							! Array of increments of predefined field variables
		character(len=80), intent(in) :: cmname						! User-defined material name, left justified
		real(kind=dbl), intent(in), dimension(nprops) :: props		! User-specified array of material constants associated with this user material
!f2py	depend(nprops) :: props
		real(kind=dbl), intent(in), dimension(3) :: coords			! An array containing the coordinates of this point
		real(kind=dbl), intent(in), dimension(3,3) :: drot			! Rotation increment matrix. This matrix represents the increment of rigid body rotation of the basis system in which the components of stress (STRESS) and strain (STRAN) are stored
		real(kind=dbl), intent(in) :: celent						! Characteristic element length, which is a typical length of a line across an element for a first-order element; it is half of the same typical length for a second-order element
		real(kind=dbl), intent(in), dimension(3,3) :: dfgrd0		! Array containing the deformation gradient at the beginning of the increment
		real(kind=dbl), intent(in), dimension(3,3) :: dfgrd1		! Array containing the deformation gradient at the end of the increment

		! --------------------------------------------------------------
		! output variables in all situations
		! --------------------------------------------------------------
		real(kind=dbl), intent(out), dimension(ntens,ntens) :: ddsdde		! Jacobian matrix of the constitutive model ddelta_sigma/ddelta_lambda
!f2py	depend(ntens) :: ddsdde
		real(kind=dbl), intent(inout), dimension(ntens) :: stress			! This array is passed in as the stress tensor at the beginning of the increment and must be updated in this routine to be the stress tensor at the end of the increment
!f2py	depend(ntens) :: stress
!f2py	intent(in,out) :: stress
		real(kind=dbl), intent(inout), dimension(nstatv) :: statev			! An array containing the solution-dependent state variables. These are passed in as the values at the beginning of the increment
!f2py	depend(nstatv) :: statev
!f2py	intent(in,out) :: statev
		real(kind=dbl), intent(inout) :: sse								! Specific elastic strain energy, passed in as the values at the start of the increment and should be updated to the corresponding specific energy values at the end of the increment
!f2py	intent(in,out) :: sse
		real(kind=dbl), intent(inout) :: spd								! plastic dissipation, passed in as the values at the start of the increment and should be updated to the corresponding specific energy values at the end of the increment
!f2py	intent(in,out) :: spd
		real(kind=dbl), intent(inout) :: scd								! “creep” dissipation, passed in as the values at the start of the increment and should be updated to the corresponding specific energy values at the end of the increment
!f2py	intent(in,out) :: scd

		! output variables Only in a fully coupled thermal-stress or a coupled thermal-electrical-structural analysis
		real(kind=dbl), intent(out) :: rpl							! Volumetric heat generation per unit time at the end of the increment caused by mechanical working of the material
		real(kind=dbl), intent(out), dimension(ntens) :: ddsddt		! Variation of the stress increments with respect to the temperature
!f2py	depend(ntens) :: ddsddt
		real(kind=dbl), intent(out), dimension(ntens) :: drplde		! Variation of RPL with respect to the strain increments
!f2py	depend(ntens) :: drplde
		real(kind=dbl), intent(out) :: drpldt						! Variation of RPL with respect to the temperature

		! output variables Only in a geostatic stress procedure or a coupled pore fluid diffusion/stress analysis for pore pressure cohesive elements
		!real(kind=dbl), intent(out) :: rpl							! RPL is used to indicate whether or not a cohesive element is open to the tangential flow of pore fluid

		! variables that can be updated
		real(kind=dbl), intent(out) :: pnewdt						! Ratio of suggested new time increment to the time increment being used. This variable allows you to provide input to the automatic time incrementation algorithms in Abaqus/Standard.


		! --------------------------------------------------------------
		! internal variables
		! --------------------------------------------------------------
		real(kind=dbl), dimension(6) :: sigma						! working variable of stress vector
		real(kind=dbl), dimension(6) :: sigma_trial					! working variable of trial stress vector
		real(kind=dbl), dimension(6) :: sigma_trial_old				! trial stress vector of the last step
		real(kind=dbl), dimension(6) :: eps							! working variable of strain vector
		real(kind=dbl), dimension(6) :: delta_eps					! working variable of strain increment vector
		real(kind=dbl), dimension(6) :: eps_new						! working variable of updated strain vector

		integer :: n_act_yf											! number of active yield functions

		real(kind=dbl), dimension(nprops-4) :: mat_data				! internal material parameters vector (just pure parameters, without n_yf, ... variables)

		real(kind=dbl), dimension(6) :: eps_pl, eps_pl_old			! plastic strains vector
		real(kind=dbl), dimension(6) :: eps_el						! elastic strains vector
		real(kind=dbl), dimension(6) :: delta_eps_pl				! plastic strain increments vector, from the last step on entry of UMAT (via state vars)
																	! and to be updated in UMAT for the current increment

		real(kind=dbl), allocatable, dimension(:) :: delta_lambda	! plastic multiplier vector
		real(kind=dbl), allocatable, dimension(:) :: alpha_old		! internal variables vector at the end of last step
		real(kind=dbl), allocatable, dimension(:) :: alpha_new		! updated internal variables vector

		logical, allocatable, dimension(:) :: yield_flags			! yield flags vector, defines which yield function is active and which not
		real(kind=dbl), allocatable, dimension(:) :: f_trial		! yield functions values for the trial state
		logical :: elastic_inc_flag										! flag which indicates if the increment is elastic (true) or not (false), needed for more efficient damage computation

		real(kind=dbl), dimension(6,6) :: C_el, C_el_inv			! elasticity matrix (sig=C_el*eps) and inverse elasticity matrix (eps=C_el_inv*sig)
		real(kind=dbl) :: E1										! elasticity module, in case of transversal isotropy: elasticity module parallel to the structure (isotropic) plane
		real(kind=dbl) :: nu1										! in case of transversal isotropy: poisson ration perpendicular to the structure (isotropic) plane
		real(kind=dbl) :: E2										! in case of transversal isotropy: elasticity module perpendicular to the structure (isotropic) plane
		real(kind=dbl) :: G2										! in case of transversal isotropy: shear modulus for shear loading in/parallel to the structure (isotropic) plane
		real(kind=dbl) :: nu2										! in case of transversal isotropy: poisson ration in the structure (isotropic) plane

		real(kind=dbl), dimension(6,6) :: C_ep						! elastic-plastic tangent

		real(kind=dbl), allocatable, dimension(:) :: alpha_dam		! internal damage variable vector at the end of last step
		real(kind=dbl), allocatable, dimension(:) :: omega			! damage variable vector, at the end of the last step on entry of UMAT (via state vars)
																	! and to be updated in UMAT for the current increment

		integer :: n_yf_c_tested									! number of tested yield flags combinations
		integer :: n_newton_it										! number of performed newton iterations
		integer :: status_var										! status variable
		real(kind=dbl), dimension(2) :: comp_time					! array for time stamps of beginning and ending of execution time
		! ##########################################################################

!		write(*,*) 'Fortran Code entered: beginning of UMAT'

		! get beginning time stamp
		call CPU_TIME(comp_time(1))


		! ##########################################################################
		! DATA UNPACKING AND PREPARING
		! ##########################################################################

		! set proposed time increment factor to one
		pnewdt = 1._dbl

		! check for number of stress and strain components
		! and if neccessary convert them to array suitable
		! for 3D return mapping
		call conv_to_3D(ndi, nshr, &
	                         stress, stran, dstran, &
	                         sigma, eps, delta_eps, &
	                         status_var)

	    ! check conversion status
	    if (status_var /= 0) then
	    	! an error occured in conversion
	    	! set error flag at the end of the state variables vector
	    	! to the value of status_var and return
	    	statev(nstatv) = real(status_var)
	    	return
	    end if

		! perform setup of the material model
		! and store relevant values in the variables defined in the material_info module
		call material_setup(props(1:4), status_var)
		if (status_var /= 0) then
	    	! material setup failed due to invalid model components definition
	    	statev(nstatv) = real(status_var)
	    	return
	    end if

		! unpack material data frpm props array
		! and store the parameters in the mat_pars variable
		! of the derived type material_data stored in the module derived_types
		call unpack_material_data(props(5:))

		! check plausibility of internal material_data vector, needs to be adjusted for anisotropic model/material params
!		call check_parameters(mat_data, celent, status_var)
!		if (status_var /= 0) then
!	    	! material parameters wrong
!	    	! set error flag at the end of the state variables vector
!	    	! to the value of status_var and return
!	    	statev(nstatv) = real(status_var)
!			return
!		end if

		! unpack first element of statevars vector
		n_act_yf = int(statev(1))

		! allocate vectors
		allocate(delta_lambda(n_yf), stat=status_var)
		allocate(yield_flags(n_yf), stat=status_var)

		! allocate internal hardening variables vector alpha
		! depending on the presence of a hardening law
		if (hardening_present) then
			allocate(alpha_old(n_int_vars), stat=status_var)	! hardening law present with one or more internal hardening variables
		else
			allocate(alpha_old(1), stat=status_var)				! no hardening law, hence allocate one entry for dummy hardening variable, which remains zero
		end if

		! allocate alpha_d and omega vectors only in case of active damage flag
		damage_check2: if (damage_present) then
			allocate(alpha_dam(n_int_vars_dam), stat=status_var)	! damage law present with one or more alpha_d and one or more omega
			allocate(omega(n_omega), stat=status_var)				! alloctate alpha_d and omega vectors
		else
			allocate(alpha_dam(1), stat=status_var)			! no damage law, hence allocate one entry for dummy alpha_d and dummy omega, which remain zero
			allocate(omega(1), stat=status_var)
		end if damage_check2

		! allocation check
		if (status_var /= 0) then
			status_var = 30
			! an error occured in allocation of the vectors
	    	! set error flag at the end of the state variables vector
	    	! to the value of status_var and return
	    	statev(nstatv) = real(status_var)
			return
		end if

		! unpack state variables with converged values of last step
		! elastic-plastic tangent of last increment is assigned to C_ep
		call unpack_state_vars_new(statev, &
									eps_pl, eps_el, &
									delta_lambda, alpha_old, yield_flags, &
									alpha_dam, delta_eps_pl, omega, &
									C_ep)

		! extract elastic material parameters from the props vector
		! and compute elasticity matrix
		if (anisotropic_elastic_flag .eqv. .False.) then
			! isotropic elasticity
			E1 = mat_pars%E1
			nu1 = mat_pars%nu1
			C_el = calc_C_elastic(E1, nu1)	! compute isotropic elasticity matrix with E and nu
			C_el_inv = calc_inverse_C_elastic(E1,nu1) ! compute inverse isotropic elasticity matrix with e and nu
		else
			! transverse isotropy
			E1 = mat_pars%E1
			nu1 = mat_pars%nu1
			E2 = mat_pars%E2
			G2 = mat_pars%G2
			nu2 = mat_pars%nu2
			C_el = calc_transv_istr_C_elastic(E1, E2, G2, nu1, nu2)				! compute transversal isotropic elasticity matrix with E1, nu1, E2, G2 and nu2
			C_el_inv = calc_inv_transv_istr_C_elastic(E1, E2, G2, nu1, nu2)		! compute inverse transversal isotropic elasticity matrix with E1, nu1, E2, G2 and nu2
		end if


		! ##########################################################################
		! check if umat is called for computation of C_ep only
		! ##########################################################################

		! check if umat is called for computation of C_ep only
	    ! This step is important because at the start of the first increment
	    ! of a new step Abaqus calls the UMAT (update_sigma) with a zero strain increment as input value
	    ! to obtain a first approximation for the elastoplastic tangent C_ep.
		! check if strain increment vector equals zero
		Cep_only_check: if (all(dstran == 0._dbl)) then

			! zero strain increment -> only elastic plastic tangent matrix from last increment needs to be returned
			! convert stress vector and elastoplastic tangent back to input dimensions
			! and assign them to the return variables stress and ddsdde
			if (kinc == 1 .and. kstep == 1) then
				call conv_to_plane_cond(ndi, nshr, &
								        sigma, C_el, &
								        stress, ddsdde, &
								        status_var)
			else
				call conv_to_plane_cond(ndi, nshr, &
								        sigma, C_ep, &
								        stress, ddsdde, &
								        status_var)
			end if

			! deallocate arrays
			deallocate(delta_lambda, alpha_old, yield_flags, alpha_dam, omega)

	    	! set error flag at the end of the state variables vector
	    	! to the value of status_var and return
	    	statev(nstatv) = real(status_var)

			! get ending time stamp
			call CPU_TIME(comp_time(2))
			! compute execution time and write to std-out / .msg file
			!write(*,*) 'execution time: ', comp_time(2)-comp_time(1)

			! exit umat
			return

		end if Cep_only_check
		! ##########################################################################



		! ##########################################################################
		! computation of trial state
		! ##########################################################################
		! computation of new strain vector
		eps_new = eps + delta_eps

		! computation of sigma trial from old stresses and plastic strain increment of last load increment
		! the nominal stresses are converted to effective stresses
		damage_check5: if (damage_present .and. damage_flag) then
			sigma_trial = sigma/(1._dbl-omega(1)) + matmul(C_el, delta_eps)
		else
			sigma_trial = sigma + matmul(C_el, delta_eps)
		end if damage_check5

!		write(*,*) '			sigma_trial:', sigma_trial

		! allocate f_trial vector
		allocate(f_trial(n_yf))

		! evaluation of all yield functions with sigma trial
		f_trial = yf(sigma_trial, alpha_old)

		write(*,*) ''
		write(*,*) 'f_trial: ', f_trial

		! compute yield flags vector and number of active yield functions
		! based on the values of f_trial
!		where ( f_trial <= (0._dbl+1.0D-8) )
!			yield_flags = .true.		! inactive yield function
!		elsewhere
!			yield_flags = .false.		! active yield function
!		end where

		yield_flags = .true.
		if (maxval(f_trial) > 0._dbl) then
			yield_flags(maxloc(f_trial)) = .false.
		end if

!		write(*,*) 'initial yield flags: ', yield_flags

		! compute number of active yield functions
		n_act_yf = n_yf - count(yield_flags)
		! ##########################################################################



		! ##########################################################################
		! check of trial state
		! ##########################################################################

		elastic_step_check: if (all(yield_flags)) then

			! all yield function values are <= 0, i.e all yield flags are true
			! => elastic step => updated stress vector is sigma trial
			! and elasto-plastic tangent C_ep equals elasticity matrix C_el
			! set elastic increment flag to true (needed for more efficient damage computation)
			elastic_inc_flag = .true.

			! set delta_lambda vector to zero
			delta_lambda = 0._dbl

			! set plastic strain increment vector to zero
			delta_eps_pl = 0._dbl

			! set elastic epsilon to new epsilon
			eps_el = eps_new

			! set values for return variables
			sse = 0._dbl		! Specific elastic strain energy
			spd = 0._dbl		! plastic dissipation
			scd = 0._dbl		! “creep” dissipation

			! set number of tested yield flags combinations
			! and number of performed newton iterations to zero
			n_yf_c_tested = 0
			n_newton_it = 0

			! check if damage is activated
			damage_check6: if (damage_present .and. damage_flag) then

				! compute nominal stresses and damaged elasticity matrix
				call comp_damage_components(sigma_trial, eps_pl, delta_eps_pl, alpha_old, celent, &
											elastic_inc_flag, &
											alpha_dam, omega, &
											C_el, C_el_inv, &
											eps_new, eps_pl)

			end if damage_check6

			! damage is deactivated, hence omit damage computations and
			! pack state variables for case of deactivated damage
			call pack_state_vars_new(statev, nstatv, &
										n_act_yf, &
										eps_pl, eps_el, &
										delta_lambda, alpha_old, yield_flags, &
										alpha_dam, delta_eps_pl, omega, &
										C_el, &
										n_yf_c_tested, n_newton_it)


			! convert stress vector and C_el back to input dimensions
			call conv_to_plane_cond(ndi, nshr, &
							        sigma_trial, C_el, &
							        stress, ddsdde, &
							        status_var)

			! deallocate arrays
			deallocate(delta_lambda, alpha_old, yield_flags, f_trial, alpha_dam, omega)


	    	! set error flag at the end of the state variables vector
	    	! to the value of status_var and return
	    	statev(nstatv) = real(status_var)

			! get ending time stamp
			call CPU_TIME(comp_time(2))
			! compute execution time and write to std-out / .msg file
			!write(*,*) 'execution time: ', comp_time(2)-comp_time(1)

			! exit umat
			return

		else elastic_step_check

			! set elastic increment flag to false (needed for more efficient damage computation)
			elastic_inc_flag = .false.

		end if elastic_step_check
		! ##########################################################################



		! ##########################################################################
		! start of multisurface stress update
		! ##########################################################################

		! some yield function values are > 0
		! allocate alpha_new
		if (hardening_present) then
			allocate(alpha_new(n_int_vars))
		else
			allocate(alpha_new(1))
		end if

		! perform multi surface stress update
		call multsurf_str_upd(sigma_trial, alpha_old, &
    										sigma, alpha_new, delta_lambda, &
    										delta_eps_pl, &
    										C_ep, &
    										C_el, &
											n_act_yf, yield_flags, &
											n_yf_c_tested, n_newton_it, &
											status_var)

		! check if stress update finished without errors
		return_mapping_check: if (status_var /= 0) then

	    	! an error occured in multi surface stress update

			! check if damage is activated
			damage_check8: if (damage_present .and. damage_flag) then

				! compute nominal stresses and damaged elasticity matrix
				call comp_damage_components(sigma_trial, eps_pl, delta_eps_pl, alpha_new, celent, &
											elastic_inc_flag, &
											alpha_dam, omega, &
											C_el, C_el_inv, &
											eps_new, eps_pl)

			end if damage_check8

			! damage is deactivated, hence omit damage computations and
			! pack state variables for case of deactivated damage
			call pack_state_vars_new(statev, nstatv, &
										n_act_yf, &
										eps_pl, eps_el, &
										delta_lambda, alpha_new, yield_flags, &
										alpha_dam, delta_eps_pl, omega, &
										C_el, &
										n_yf_c_tested, n_newton_it)


			! deallocate arrays
			deallocate(delta_lambda, alpha_old, yield_flags, &
						f_trial, alpha_new, alpha_dam, omega)

	    	! set error flag at the end of the state variables vector
	    	! to the value of status_var and return
	    	statev(nstatv) = real(status_var)

	    	! set proposed time increment
	    	pnewdt = 0.250_dbl

			! get ending time stamp
			call CPU_TIME(comp_time(2))
			! compute execution time and write to std-out / .msg file
			!write(*,*) 'execution time: ', comp_time(2)-comp_time(1)

	    	return

	    end if return_mapping_check
		! ##########################################################################



		! ##########################################################################
		! postprocessing and perparation of data for output variables
		! ##########################################################################

		! compute new plastic and elastic strains
		eps_pl_old = eps_pl				! store old plastic strains before plastic strain update, needed for numerical d_alpha_d/depsilon
		eps_pl = eps_pl + delta_eps_pl
		eps_el = eps_new - eps_pl

!		write(*,*) 'eps_pl: ', eps_pl
!		write(*,*) 'delta_eps_pl: ', delta_eps_pl

		! check if damage is activated
		damage_check10: if (damage_present .and. damage_flag) then
			! compute damage components
			call comp_damage_components(sigma, eps_pl, delta_eps_pl, alpha_new, celent, &
										elastic_inc_flag, &
										alpha_dam, omega, &
										C_ep,  C_el_inv, &
										eps_new, eps_pl_old)
		end if damage_check10

		! damage is deactivated, hence omit damage computations and
		! pack state variables for case of deactivated damage
		call pack_state_vars_new(statev, nstatv, &
									n_act_yf, &
									eps_pl, eps_el, &
									delta_lambda, alpha_new, yield_flags, &
									alpha_dam, delta_eps_pl, omega, &
									C_el, &
									n_yf_c_tested, n_newton_it)


		! set values for return variables
		sse = 0._dbl		! Specific elastic strain energy
		spd = 0._dbl		! plastic dissipation
		scd = 0._dbl		! “creep” dissipation

		! convert stress and strain vector back to input dimensions
		call conv_to_plane_cond(ndi, nshr, &
						        sigma, C_ep, &
						        stress, ddsdde, &
						        status_var)

		! deallocate arrays
		deallocate(delta_lambda, alpha_old, yield_flags, &
					f_trial, alpha_new, alpha_dam, omega)

    	! set error flag at the end of the state variables vector
    	! to the value of status_var and return
    	statev(nstatv) = real(status_var)
		! ##########################################################################

		! get ending time stamp
		call CPU_TIME(comp_time(2))
		! compute execution time and write to std-out / .msg file
		!write(*,*) 'execution time: ', comp_time(2)-comp_time(1)

	end subroutine umat_tsaiwu

end module stress_update
