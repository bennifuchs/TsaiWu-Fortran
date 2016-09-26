module material_info
    implicit none

    integer, save :: n_yf				! number of yield functions
	integer, save :: n_int_vars			! number of internal hardening variables
	integer, save :: n_int_vars_dam		! number of internal softening variables
	integer, save :: n_omega			! number of damage variables

	logical, save :: hardening_present			! True if hardening is present, False if no hardening law is detected

	logical, save :: damage_present		! True if damage law is detected, false if no damage law is detected

end module material_info
module constants
	! module containing constants
    implicit none

	! kind number of double precision float
    integer, parameter :: dbl = kind(1.0D0)

	! Pi
	real(kind=dbl), parameter :: Pi = &
			 3.14159265358979323846264338327950288_dbl

	! iteration number for switching tolerances of inner newton iteration
	integer, parameter :: k_alt = 25

	! maximum allowed iterations of inner newton iteration
	integer, parameter :: k_max = 50

	! tolerances for the inner newton iteration
	! for k < k_alt
	real(kind=dbl), parameter :: atol_def = 1.0D-9
	real(kind=dbl), parameter :: rtol_def = 1.0D-6

	! tolerances for the inner newton iteration
	! for k > k_alt
	real(kind=dbl), parameter :: atol_alt = 1.0D-8
	real(kind=dbl), parameter :: rtol_alt = 1.0D-4

	! tolerances for yield function values in plausibility check
	real(kind=dbl), parameter :: yf_tol_def = 1.0D-7		! for k in inner newton iteration < k_alt
	real(kind=dbl), parameter :: yf_tol_alt = 1.0D-6		! for k in inner newton iteration > k_alt

	! maximum value of damage variable
	real(kind=dbl), parameter :: omega_max = 0.99_dbl

	! flag for switch between numerical/analytical computation of jacobian
	logical, parameter :: jac_num_flag = .True.			! computation via numerical differentiation of R-vector: .True.
														! computation via analytical derivatives: .False.

!	! flag to check if model formulation contains a damage law
!	logical, save :: damage_present				! True if damage law is detected, false if no damage law is detected

	! flag for activation/deactivation of damage computation
	logical, parameter :: damage_flag = .False.			! damage ist taken into account: .True.
														! damage is neglected: .False.

	! flag for switch between numerical/analytical computation of damaged elastoplastic tangent
	logical, parameter :: Cep_dam_num_flag = .False.	! Cep_dam is computed analytically: .False.
														! C_ep_dam is computed numerically: .True.

	! flag for switch between isotropic/transversal isotropic elasticity
	logical, parameter :: anisotropic_elastic_flag = .True.				! isotropic elasticity: .False.
																		! transversal isotropic elasticity: .True.

end module constants
module derived_types
    use constants
    implicit none

	! save defined data values between different loads of the module
	save

	! derived data type for storage of material/model parameters
	type :: material_data
		real(kind=dbl) :: E1	! elasticity modulus for loading in foliation plane
		real(kind=dbl) :: nu1	! poisson ratio for lateral in plane strains due to in plane loading
		real(kind=dbl) :: E2	! elasticity modulus perpendicular to foliation plane
		real(kind=dbl) :: G2	! shear modulus in planes perependicular to the foliation plane
		real(kind=dbl) :: nu2	! poisson ratio for lateral in plane strains due to loading perpendicular to the foliation plane
		real(kind=dbl) :: F2	! independent components of transversely isotropic tsai-wu structural vector [F2, F2, F3, 0, 0, 0]
		real(kind=dbl) :: F3
		real(kind=dbl) :: F22	! independent components of transversely isotropic tsau-wu structural matrix
		real(kind=dbl) :: F33	! | F22 F12 F23  0   0       0     |
		real(kind=dbl) :: F44	! |     F22 F23  0   0       0     |
		real(kind=dbl) :: F12	! |         F33  0   0       0     |
		real(kind=dbl) :: F23	! |             F44  0       0     |
	end type material_data		! |  symm.          F44      0     |
								! |                     2(F22-F12) |


	type(material_data) :: mat_pars

end module derived_types
module interface_definitions

! 	use constants

    interface

		! BLAS function for computation of L2 norm
    	function dnrm2( n, x, incx )
       		use constants
    		implicit none

    		! input variables
    		integer,intent(in) ::						n
    		real(kind=dbl), dimension(:), intent(in) :: x
    		integer, intent(in) ::						incx

    		! output variables
    		real(kind=dbl) :: 							dnrm2

    	end function dnrm2


		! LAPACK subroutine for solving system of linear equations
		subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
			use constants
			implicit none

			! input variables
			integer, intent(in) :: n
			integer, intent(in) :: nrhs
			integer, intent(in) :: lda
			integer, intent(in) :: ldb

			! output variables
			real(kind=dbl), dimension(lda,n), intent(inout) :: a
			integer, dimension(n), intent(out) :: ipiv
			real(kind=dbl), dimension(lda,nrhs), intent(inout) :: b
			integer, intent(out) :: info

		end subroutine dgesv


		! LAPACK subroutine for LU factorization
		subroutine dgetrf(m, n, a, lda, ipiv, info)
			use constants
			implicit none

			! input variables
			integer, intent(in) :: m
			integer, intent(in) :: n
			integer, intent(in) :: lda

			! output variables
			real(kind=dbl), dimension(lda,n), intent(inout) :: a
			integer, dimension(n), intent(out) :: ipiv
			integer, intent(out) :: info

		end subroutine dgetrf


		! LAPACK subroutine for computation of inverse matrix based on LU factorization
		subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
			use constants
			implicit none

			! input variables
			integer, intent(in) :: n
			integer, intent(in) :: lda
			integer, intent(in) :: lwork
			integer, dimension(n), intent(in) :: ipiv

			! output variables
			real(kind=dbl), dimension(lda,n), intent(inout) :: a
			real(kind=dbl), dimension(lwork), intent(out) :: work
			integer, intent(out) :: info

		end subroutine dgetri


    end interface

end module interface_definitions
module aux_routines

	use constants
	use material_info

	implicit none

	! Module containing auxiliary routines used in various parts of the Umat

contains



	subroutine conv_to_3D(n_dir, n_shr, &
	                         sigma_in, epsilon_in, delta_epsilon_in, &
	                         sigma_out, epsilon_out, delta_epsilon_out, &
	                         err)
		implicit none
	    ! function that converts stress and strain vectors from plane stress or plane strain analyses
	    ! to 3D stress and strain vectors
	    ! by assigning the values of the plane stress / plane strain vectors to the right positions of the
	    ! 3D stress and strain vectors

		! --------------------------------------------------------------------------
		! passed variables
	    integer, intent(in) :: n_dir            						! number of direct stress/strain components
	    integer, intent(in) :: n_shr            						! number of shear stress/strain components
	    real(kind=dbl), intent(in), dimension(:) :: sigma_in    			! input stress vector
	    real(kind=dbl), intent(in), dimension(:) :: epsilon_in    		! input strain vector
	    real(kind=dbl), intent(in), dimension(:) :: delta_epsilon_in  	! input strain increment vector

		! return variables
	    real(kind=dbl), intent(out), dimension(6) :: sigma_out    				! output stress vector
	    real(kind=dbl), intent(out), dimension(6) :: epsilon_out    				! output strain vector
	    real(kind=dbl), intent(out), dimension(6) :: delta_epsilon_out    		! output strain increment vector

		integer, intent(out) :: err				! error variable

	    ! internal variables
		! --------------------------------------------------------------------------


	    dim_check: if ( (n_dir==2) .and. (n_shr==1) ) then
	        ! plane stress element (sig11, sig22, sig12)
			! cannot be converted

!	        ! move stress and strain components in the right position of the 3D vectors
!	        sigma_out = (/ sigma_in(1), sigma_in(2), 0._dbl, sigma_in(3), 0._dbl, 0._dbl /)
!	        epsilon_out = (/ epsilon_in(1), epsilon_in(2), 0._dbl, epsilon_in(3), 0._dbl, 0._dbl /)
!	        delta_epsilon_out = (/ delta_epsilon_in(1), delta_epsilon_in(2), 0._dbl, delta_epsilon_in(3), 0._dbl, 0._dbl /)

	        ! set error variable to 20
	        err = 20

	    else if ( (n_dir==3) .and. (n_shr==1) ) then
	        ! plane strain element (sig11, sig22, sig33, sig12)

	        ! move stress and strain components in the right position of the 3D vectors
	        sigma_out = (/ sigma_in(1), sigma_in(2), sigma_in(3), sigma_in(4), 0._dbl, 0._dbl /)
	        epsilon_out = (/ epsilon_in(1), epsilon_in(2), epsilon_in(3), epsilon_in(4), 0._dbl, 0._dbl /)
	        delta_epsilon_out = (/ delta_epsilon_in(1), delta_epsilon_in(2), delta_epsilon_in(3), &
	        							delta_epsilon_in(4), 0._dbl, 0._dbl /)

			! set error variable zero
			err = 0

	    else if ( (n_dir==3) .and. (n_shr==3) ) then
	        ! 3D element (sig11, sig22, sig33, sig12, sig13, sig23)

	        ! do nothing
	        sigma_out = sigma_in
	        epsilon_out = epsilon_in
	        delta_epsilon_out = delta_epsilon_in

	        ! set error variable zero
	        err = 0

		else
			! set output arrays to zero
			sigma_out = 0._dbl
			epsilon_out = 0._dbl
			delta_epsilon_out = 0._dbl

			! set error variable to 20
			err = 20

		end if dim_check


	end subroutine conv_to_3D




	subroutine conv_to_plane_cond(n_dir, n_shr, &
						                     sigma_in, C_ep_in, &
						                     sigma_out, C_ep_out, &
						                     err)
	    ! function that converts stress vector and elastic-plastic matrix from 3D stress and strain analyses
	    ! to plane stress or plane strain vectors
	    ! by assigning the values of the 3D stress and strain vectors to the right positions of the
	    ! plane stress / plane strain vectors
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
	    integer, intent(in) :: n_dir            						! number of direct stress/strain components
	    integer, intent(in) :: n_shr            						! number of shear stress/strain components
	    real(kind=dbl), intent(in), dimension(:) :: sigma_in    		! input stress vector (3D)
	    real(kind=dbl), intent(in), dimension(:,:) :: C_ep_in    		! input elastic plastic matrix (3D)

		! return variables
	    real(kind=dbl), intent(out), dimension(n_dir+n_shr) :: sigma_out    			! output stress vector (plane conditions)
	    real(kind=dbl), intent(out), dimension(n_dir+n_shr,n_dir+n_shr) :: C_ep_out    	! output strain vector (plane conditions)

		integer, intent(out) :: err				! error variable

	    ! internal variables
		! --------------------------------------------------------------------------

	    dim_check: if ( (n_dir==2) .and. (n_shr==1) ) then
	        ! plane stress element (sig11, sig22, sig12)
			! cannot be converted

	        ! set error variable to 21
	        err = 21

	    else if ( (n_dir==3) .and. (n_shr==1) ) then
	        ! plane strain element (sig11, sig22, sig33, sig12)

	        ! move stress and strain components in the right position of the 3D vectors
	        sigma_out = (/ sigma_in(1), sigma_in(2), sigma_in(3), sigma_in(4) /)
	        C_ep_out = C_ep_in(1:4,1:4)

			! set error variable zero
			err = 0

	    else if ( (n_dir==3) .and. (n_shr==3) ) then
	        ! 3D element (sig11, sig22, sig33, sig12, sig13, sig23)

	        ! do nothing
	        sigma_out = sigma_in
			C_ep_out = C_ep_in

	        ! set error variable zero
	        err = 0

		else
			! set output arrays to zero
			sigma_out = 0._dbl
			C_ep_out = 0._dbl

			! set error variable to 21
			err = 21

		end if dim_check

	end subroutine conv_to_plane_cond



	subroutine unpack_state_vars_new(state_vars, &
									eps_pl, eps_el, &
									delta_lambda, alpha, yield_flags, &
									alpha_dam, delta_eps_pl, omega, &
									C_ep_old)
		! unpacks the state variables vector which is passed from ABAQUS to the umat
		! and which contains the state variables at the end of the last step
		! this version is intended for activated damage
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in), dimension(:) :: state_vars		! state variables array

		! output variables
		real(kind=dbl), intent(out), dimension(6) :: eps_pl			! plastic strain vector at the end of the last step (eps_pl n)
		real(kind=dbl), intent(out), dimension(6) :: eps_el			! elastic strain vector at the end of the last step (eps_el n)

		real(kind=dbl), intent(out), dimension(n_yf) :: delta_lambda		! delta Lambda vector at the end of the last step (delta_lambda n)
		real(kind=dbl), intent(out), dimension(n_int_vars) :: alpha			! internal variables vector at the end of the last step (alpha n)
		logical, intent(out), dimension(n_yf) :: yield_flags				! active yield functions vector of the last step (yield_flags n)
!		real(kind=dbl), intent(out), dimension(n_int_vars) :: alpha_old		! internal variables vector at the beginning of the last step (alpha n-1)
!		real(kind=dbl), intent(out), dimension(6) :: sigma_trial			! trial stress vector of the last step (sigma_trial n-1)

		real(kind=dbl), intent(out), dimension(n_int_vars_dam) :: alpha_dam	! internal damage variables vector at the end of the last step (alpha n)
		real(kind=dbl), intent(out), dimension(6) :: delta_eps_pl			! plastic strain increment vector computed in the last step (delta_eps_pl n = eps_pl n - eps_pl n-1)
		real(kind=dbl), intent(out), dimension(n_omega) :: omega			! damage variables vector at the end of the last step (omega n)

		real(kind=dbl), intent(out), dimension(6,6) :: C_ep_old				! elastic-plastic tangent matrix of the last step (C_ep n)

		! internal variables
		integer :: i_start			! start index
		integer :: i_end			! end index
		! --------------------------------------------------------------------------

		! unpack remaining elements of statevars vector:

		! plastic strains at end of last step
		i_start = 1+1
		i_end = i_start+6-1
		eps_pl = state_vars(i_start:i_end)

		! elastic strains at end of last step
		i_start = i_end+1
		i_end = i_start+6-1
		eps_el = state_vars(i_start:i_end)

		! plastic multipliers at end of last step
		i_start = i_end+1
		i_end = i_start+n_yf-1
		delta_lambda = state_vars(i_start:i_end)

		! internal hardening variables alpha at end of last step
		i_start = i_end+1
		i_end = i_start+n_int_vars-1
		if (hardening_present) then
			alpha = state_vars(i_start:i_end)
		else
			alpha = 0._dbl
		end if

		! yield flags at end of last step
		i_start = i_end+1
		i_end = i_start+n_yf-1
		yield_flags = conv_to_logical( int(state_vars(i_start:i_end)) )

		! internal softening variables alpha_d at end of last step
		i_start = i_end+1
		i_end = i_start+n_int_vars_dam-1
		if (damage_present) then
			alpha_dam = state_vars(i_start:i_end)
		else
			alpha_dam = 0._dbl
		end if

		! plastic strain incerement at end of last step
		i_start = i_end+1
		i_end = i_start+6-1
		delta_eps_pl = state_vars(i_start:i_end)

		! damage variables omega at end of last step
		i_start = i_end+1
		i_end = i_start+n_omega-1
		if (damage_present) then
			omega = state_vars(i_start:i_end)
		else
			omega = 0._dbl
		end if

		! elastic-plastic tangent matrix at end of last step
		i_start = i_end+1
		i_end = i_start+36-1
		C_ep_old = reshape(state_vars(i_start:i_end), (/6,6/))

	end subroutine unpack_state_vars_new



	subroutine pack_state_vars_new(state_vars, n_state_vars, &
									n_act_yf, &
									eps_pl, eps_el, &
									delta_lambda, alpha, yield_flags, &
									alpha_dam, delta_eps_pl, omega, &
									C_ep_old, &
									n_yf_c_tested, n_newton_it)
		! packs the state variables vector which is returned from the umat to ABAQUS
		! and which contains the state variables at the end of the current step
		! this version is intended for activated damage
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		integer, intent(in) :: n_state_vars				! number of state variables (length of state variables vector)
		integer, intent(in) :: n_act_yf					! number of active yield functions
		real(kind=dbl), intent(in), dimension(:) :: eps_pl		! plastic strain vector
		real(kind=dbl), intent(in), dimension(:) :: eps_el		! elastic strain vector
		real(kind=dbl), intent(in), dimension(:) :: delta_lambda	! delta Lambda vector
		real(kind=dbl), intent(in), dimension(:) :: alpha			! internal variables vector
		logical, intent(in), dimension(:) :: yield_flags			! active yield functions vector
!		real(kind=dbl), intent(in), dimension(:) :: alpha_old		! internal variables vector at the beginning of the step
!		real(kind=dbl), intent(in), dimension(:) :: sigma_trial		! trial stress vector

		real(kind=dbl), intent(in), dimension(:) :: alpha_dam		! internal damage variable vector
		real(kind=dbl), intent(in), dimension(:) :: delta_eps_pl	! plastic strain increment vector
		real(kind=dbl), intent(in), dimension(:) :: omega			! damage variable vector

		real(kind=dbl), intent(in), dimension(:,:) :: C_ep_old		! elastic-plastic tangent matrix

		integer, intent(in) :: n_yf_c_tested						! number of tested yield flags combinations
		integer, intent(in) :: n_newton_it							! number of performed newton iterations

		! return variables
		real(kind=dbl), intent(out), dimension(n_state_vars) :: state_vars			! state variables array

		! internal variables
		integer :: i_start			! start index
		integer :: i_end			! end index
		! --------------------------------------------------------------------------

		! pack state variables array

		! 1. number of active yield functions
		i_start = 1
		i_end = i_start+1-1
		state_vars(i_start:i_end) = real(n_act_yf)

		! 2. plastic strain vector
		i_start = i_end+1
		i_end = i_start+6-1
		state_vars(i_start:i_end) = eps_pl

		! 3. elastic strain vector
		i_start = i_end+1
		i_end = i_start+6-1
		state_vars(i_start:i_end) = eps_el

		! 4. incremental plastic multiplier
		i_start = i_end+1
		i_end = i_start+n_yf-1
		state_vars(i_start:i_end) = delta_lambda

		! 5. internal hardening variable (if model contains hardening law)
		i_start = i_end+1
		i_end = i_start+n_int_vars-1
		if (hardening_present) then
			state_vars(i_start:i_end) = alpha
		end if

		! 6. yield flags
		i_start = i_end+1
		i_end = i_start+n_yf-1
		state_vars(i_start:i_end) = real( conv_to_int(yield_flags) )

		! 7. internal softeing variable (if model contains damage law)
		i_start = i_end+1
		i_end = i_start+n_int_vars_dam-1
		if (damage_present) then
			state_vars(i_start:i_end) = alpha_dam
		end if

		! 7. plastic strain increment vector
		i_start = i_end+1
		i_end = i_start+6-1
		state_vars(i_start:i_end) = delta_eps_pl

		! 8. damage variable (if model contains damage law)
		i_start = i_end+1
		i_end = i_start+n_omega-1
		if (damage_present) then
			state_vars(i_start:i_end) = omega
		end if

		! 9. elasto-plastic tangent matrix
		i_start = i_end+1
		i_end = i_start+36-1
		state_vars(i_start:i_end) = reshape(C_ep_old, (/ 36 /) )

		! 10. number of tested yield function combinations
		i_start = i_end+1
		i_end = i_start+1-1
		state_vars(i_start:i_end) = real(n_yf_c_tested)

		! 11. number of newton iterations
		i_start = i_end+1
		i_end = i_start+1-1
		state_vars(i_start:i_end) = real(n_newton_it)

	end subroutine pack_state_vars_new



	subroutine pack_x_vector(sigma, alpha, delta_lambda_act, x)
		! packs the sigma, alpha and delta_lambda_act vectors in one large vector
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma
		real(kind=dbl), dimension(:), intent(in) :: alpha
		real(kind=dbl), dimension(:), intent(in) :: delta_lambda_act

		! return variable
		real(kind=dbl), dimension(6+n_int_vars+size(delta_lambda_act)), intent(out) :: x
!f2py	depend(n_int_vars,delta_lambda_act) :: x

		! internal variable
		integer :: i_start
		integer :: i_end
		! --------------------------------------------------------------------------

		! pack sigma, alpha and delta_lambda_act vector
		i_start = 1
		i_end = 6
		x(i_start:i_end) = sigma

		i_start = i_end+1
		i_end = i_start+n_int_vars-1
		if (hardening_present) then
			x(i_start:i_end) = alpha
		end if

		i_start = i_end+1
		i_end = i_start+size(delta_lambda_act)-1
		x(i_start:i_end) = delta_lambda_act

	end subroutine pack_x_vector



	subroutine unpack_x_vector(x, sigma, alpha, delta_lambda_act)
		! unpacks the sigma, alpha and delta_lambda_act vectors from the large x vector
		! not suitable for call from python because alpha is returned as assumed shape array
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: x

		! return variable
		real(kind=dbl), dimension(6), intent(out) :: sigma
		real(kind=dbl), dimension(:), intent(out) :: alpha
		real(kind=dbl), dimension(size(x)-6-n_int_vars), intent(out) :: delta_lambda_act

		! internal variables
		integer :: n_act_yf
		integer :: i_start
		integer :: i_end
		! --------------------------------------------------------------------------
		! compute number of active yield functions
		n_act_yf = size(x)-6-n_int_vars

		! unpack sigma, alpha and delta_lambda_act vector
		i_start = 1
		i_end = 6
		sigma = x(i_start:i_end)

		i_start = i_end+1
		i_end = i_start+n_int_vars-1
		if (hardening_present) then
			alpha = x(i_start:i_end)
		else
			alpha = 0._dbl
		end if

		i_start = i_end+1
		i_end = i_start+n_act_yf-1
		delta_lambda_act = x(i_start:i_end)

	end subroutine unpack_x_vector



	subroutine unpack_material_data(mat_data_vec)
		! unpacks material data vector and return entries in
		! the variable mat_pars of the derived data type 'material_data'
		! this variable is defined in the module derived_types
		use derived_types
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: mat_data_vec

		! internal variables
		! --------------------------------------------------------------------------

		! unpack remaining parameters
		mat_pars%E1 = mat_data_vec(1)
		mat_pars%nu1 = mat_data_vec(2)
		mat_pars%E2 = mat_data_vec(3)
		mat_pars%G2 = mat_data_vec(4)
		mat_pars%nu2 = mat_data_vec(5)
		mat_pars%F2 = mat_data_vec(6)
		mat_pars%F3 = mat_data_vec(7)
		mat_pars%F22 = mat_data_vec(8)
		mat_pars%F33 = mat_data_vec(9)
		mat_pars%F44 = mat_data_vec(10)
		mat_pars%F12 = mat_data_vec(11)
		mat_pars%F23 = mat_data_vec(12)

	end subroutine unpack_material_data



!	subroutine check_parameters(mat_data, l_char, status_var)
!		! checks plausibility of model/material parameters
!
!		use derived_types
!		implicit none
!
!		! --------------------------------------------------------------------------
!		! passed variables
!		real(kind=dbl), dimension(:), intent(in) :: mat_data		! material data vector
!		real(kind=dbl), intent(in) :: l_char						! characteristic element length
!
!		! return variable
!		integer, intent(out) :: status_var							! exit status variable
!
!		! internal variables
!		real(kind=dbl) :: E				! elasticity modulus
!		real(kind=dbl) :: nu			! poisson ratio
!		real(kind=dbl) :: m0			! friction parameter
!		real(kind=dbl) :: fcu			! uniaxial compressive ultimate strength
!		real(kind=dbl) :: ecc			! eccentricity parameter
!		real(kind=dbl) :: fcy			! uniaxial compressive yield strenght
!		real(kind=dbl) :: GfI			! mode I fracture energy
!		real(kind=dbl) :: Ad			! model parameter for shape of damage evolution function (Ad > 0, i.e. positive)
!		real(kind=dbl) :: Bd			! model parameter for shape of damage evolution function (Bd > 1.)
!
!		real(kind=dbl) :: ftu			! uniaxial tensile strength, needs to be computed from fcu, e and m0
!		logical :: test_res				! logical value for result of plausibility check
!		! --------------------------------------------------------------------------
!
!		E = mat_pars%E1
!		nu = mat_pars%nu1
!		m0 = mat_pars%m0
!		fcu = mat_pars%fcu
!		ecc = mat_pars%e
!		fcy = mat_pars%fcy
!		GfI = mat_pars%GfI
!		Ad = mat_pars%Ad
!		Bd = mat_pars%Bd
!
!		! computation of ftu
!		ftu = fcu/(6._dbl*ecc)*(-m0*(ecc+1._dbl)+sqrt(m0**2*(ecc+1._dbl)**2+36._dbl*ecc**2))
!
!		! perform checks and create assign result to boolean variable
!		test_res = ( (E > 0._dbl) .and. &
!						(nu > 0._dbl .and. nu < 0.5_dbl) .and. &
!						(fcu > fcy) .and. &
!						(m0 > 0._dbl .and. ecc > 0.5_dbl .and. ecc <= 1._dbl) .and. &
!						(GfI > 0._dbl) .and. &
!						(l_char < (2._dbl*E*GfI)/ftu**2) .and. &
!						(Ad > 0._dbl .and. Bd > 1._dbl) )
!
!		! set status variable according to the value of the boolean variable and return
!		if (test_res) then
!			status_var = 0
!			return
!		else
!			status_var = 60
!			return
!		end if
!
!	end subroutine check_parameters



	subroutine conv_param(fc, ft, mu, tau)
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in) :: fc			! uniaxial compressive strength
		real(kind=dbl), intent(in) :: ft			! uniaxial tensile strength

		! return variables
		real(kind=dbl), intent(out) :: mu			! friction angle
		real(kind=dbl), intent(out) :: tau		! cohesion
		! --------------------------------------------------------------------------

		! computation of mu and tau
		mu = -(ft-fc)/(fc+ft) * sqrt(2.)
		tau = 2./sqrt(3.) * fc*ft/(fc+ft)

	end subroutine conv_param



	function conv_to_logical(int_array)
		! function that converts an 1D-integer array to a 1D-logical array based on following rules:
		! int_array(i) = 0	=>	.TRUE.
		! int_array(i) = 1	=>	.FALSE.
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		integer, dimension(:), intent(in) :: int_array		! input array

		! return variables
		logical, dimension(size(int_array)) :: conv_to_logical			! output array

		! internal variables
		integer :: i										! index variable
		! --------------------------------------------------------------------------

		! loop over all entries of int_array
		do i = 1,size(int_array)

			! check the value of the integer array entry
			if (int_array(i) == 0) then
				conv_to_logical(i) = .true.
			else if (int_array(i) == 1) then
				conv_to_logical(i) = .false.
			end if

		end do

	end function conv_to_logical



	function conv_to_int(logical_array)
		! function that converts an 1D-logical array to a 1D-integer array based on following rules:
		! logical_array(i) == .TRUE.	=>	0
		! int_array(i) == .FALSE.	=>	1
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		logical, dimension(:), intent(in) :: logical_array		! input array

		! return variables
		integer, dimension(size(logical_array)) :: conv_to_int			! output array

		! internal variables
		integer :: i										! index variable
		! --------------------------------------------------------------------------

		! loop over all entries of logical_array
		do i = 1,size(logical_array)

			! check the value of the logical array entry
			if (logical_array(i) .eqv. .true.) then
				conv_to_int(i) = 0
			else if (logical_array(i) .eqv. .false.) then
				conv_to_int(i) = 1
			end if

		end do

	end function conv_to_int


	subroutine material_setup(mat_info, err)
		! performs identification of the material models components
		! and sets respective variables in the material_info module
		implicit none

		! input variable
		real(kind=dbl), dimension(4), intent(in) :: mat_info	! array containing the basic information about the material model components

		! output variable
		integer, intent(out) :: err

		! unpack the elements of mat_info vector and stare them in the variables defined in the material_info module
		n_yf = int(mat_info(1))
		n_int_vars = int(mat_info(2))
		n_int_vars_dam = int(mat_info(3))
		n_omega = int(mat_info(4))

		! check the values of n_yf, n_int_vars, n_int_vars_dam and n_omega to be within the valid range
		if ( (n_yf < 1) .or. (n_int_vars < 0) .or. &
			 (n_int_vars_dam < 0) .or. (n_omega < 0) ) then
			! input error for n_yf or n_int_vars
	    	! set error variable to 40 and return
			err = 40
			return
		end if

		! check number if internal hardening variables and set the hardening_present flag to the appropriate value
		if (n_int_vars == 0) then
			hardening_present = .false.
		else
			hardening_present = .true.
		end if

		! check number of damage variables  and set the damage_present flag to the appropriate value
		if ( (n_int_vars_dam == 0) .or. (n_omega == 0) ) then
			damage_present = .false.
		else
			damage_present = .true.
		end if

!		write(*,*) 'material setup entered: '
!		write(*,*) 'n_yf', n_yf
!		write(*,*) 'n_int_vars', n_int_vars
!		write(*,*) 'n_int_vars_dam', n_int_vars_dam
!		write(*,*) 'n_omega', n_omega
!		write(*,*) 'hardening present: ', hardening_present
!		write(*,*) 'damage present: ', damage_present



	end subroutine material_setup


end module aux_routines
module math_routines
	! module for various mathematical functions and subroutines

	! load additional modules
	use constants

    implicit none

contains

	function calc_rot_matrix_zyz(alpha, beta, gam)
		! computes rotation matrix (with direction cosines) for coordinate transformation of vector
		! following the zyz convention:
		! 1. rotation around z axis,
		! 2. rotation around rotated y axis,
		! 3. rotation around rotated z axis,
		! for further informations see euler angles on wikipedia
		! for definition of anisotropy (eg transverse isotropy) assume:
		! y axis to show into northern direction,
		! x axis to show into eastern direction,
		! z axis to show into skywards direction,
		! alpha to be positive in eastern direction,
		! beta to be positive in downwards direction,
		! gamma to be zero
		! for further information see Wittke: 'Rock Mechanics Based on an Anisotropic Jointed Rock Model', p. 44

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in) :: alpha		! rotation around z axis
		real(kind=dbl), intent(in) :: beta		! rotation around rotated y axis
		real(kind=dbl), intent(in) :: gam		! rotation around rotated z axis, in case if transverse isotropy zero

		! return variable
		real(kind=dbl), dimension(3,3) :: calc_rot_matrix_zyz		! rotation matrix with direction cosines between old an rotated coordinate basis

		! internal variables
		real(kind=dbl), dimension(3,3) :: Tz1	! rotation matrix for first rotation around z axis
		real(kind=dbl), dimension(3,3) :: Ty2	! rotation matrix for second rotation around rotated y axis
		real(kind=dbl), dimension(3,3) :: Tz3	! rotation matrix for third rotation around rotated z axis
		real(kind=dbl), dimension(3,3) :: Tzyz	! rotation matrix all three above combined into one: Tz1*Ty2*Tz3
		! --------------------------------------------------------------------------

!		Tz1 = reshape( (/ cos(-alpha), -sin(-alpha), 0._dbl, sin(-alpha), cos(-alpha), 0._dbl, 0._dbl, 0._dbl, 1._dbl /), (/3,3/) )
!		Ty2 = reshape( (/ cos(beta), 0._dbl, sin(beta), 0._dbl, 1._dbl, 0._dbl, -sin(beta), 0._dbl, cos(beta) /), (/3,3/) )
!		Tz3 = reshape( (/ cos(-gam), -sin(-gam), 0._dbl, sin(-gam), cos(-gam), 0._dbl, 0._dbl, 0._dbl, 1._dbl /), (/3,3/) )
		Tz1 = reshape( (/ cos(-gam), -sin(-gam), 0._dbl, sin(-gam), cos(-gam), 0._dbl, 0._dbl, 0._dbl, 1._dbl /), (/3,3/) )
		Ty2 = reshape( (/ cos(beta), 0._dbl, sin(beta), 0._dbl, 1._dbl, 0._dbl, -sin(beta), 0._dbl, cos(beta) /), (/3,3/) )
		Tz3 = reshape( (/ cos(-alpha), -sin(-alpha), 0._dbl, sin(-alpha), cos(-alpha), 0._dbl, 0._dbl, 0._dbl, 1._dbl /), (/3,3/) )
		calc_rot_matrix_zyz = matmul(Tz1, matmul(Ty2, Tz3))

!		Tzyz = reshape( (/ (-sin(-gam)*sin(-alpha)+cos(-gam)*cos(beta)*cos(-alpha)), &
!							(-sin(-gam)*cos(-alpha)-cos(-gam)*cos(beta)*sin(-alpha)), &
!							(cos(-gam)*sin(beta)), &
!							(cos(-gam)*sin(-alpha)+sin(-gam)*cos(beta)*cos(-alpha)), &
!							(cos(-gam)*cos(alpha)-sin(-gam)*cos(beta)*sin(-alpha)), &
!							(sin(-gam)*sin(beta)), &
!							(-sin(beta)*cos(-alpha)), &
!							(sin(beta)*sin(-alpha)), &
!							(cos(beta)) /), (/3,3/) )
!
!		calc_rot_matrix_zyz = Tzyz

!		write(*,*) calc_rot_matrix_zyz
!		write(*,*) Tzyz
!		write(*,*) (calc_rot_matrix_zyz - Tzyz)

	end function calc_rot_matrix_zyz



	elemental function mcauly(x)
		! computes <x> (mcauly brackets
		! if x < 0   =>  <x> = 0
		! if x >= 0  =>  <x> = x
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl) :: mcauly

		! passed variable
		real(kind=dbl), intent(in) :: x
		! --------------------------------------------------------------------------

		if (x < 0._dbl) then

			mcauly = 0._dbl

		else

			mcauly = x

		end if

	end function mcauly



	elemental function heaviside(x)
		! computes heaviside function H(x)
		! if x < 0  => H(x) = 0
		! if x >= 0  => H(x) = 1
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl) :: heaviside

		! passed variable
		real(kind=dbl), intent(in) :: x
		! --------------------------------------------------------------------------

		if (x < 0._dbl) then

			heaviside = 0._dbl

		else

			heaviside = 1._dbl

		end if

	end function heaviside


end module math_routines
module stress_routines
	! module containing various routines for computation of stresses and stress invariants

	! load additional modules
	use constants
	use math_routines, only: calc_rot_matrix_zyz

	implicit none


	! definition and initialization of constant arrays

	! deviatoric operator
	real(kind=dbl), dimension(6,6), parameter :: I_dev = &
			 reshape( (/ 2._dbl/3._dbl, -1._dbl/3._dbl, -1._dbl/3._dbl, 0._dbl, 0._dbl, 0._dbl, &
						-1._dbl/3._dbl, 2._dbl/3._dbl, -1._dbl/3._dbl, 0._dbl, 0._dbl, 0._dbl, &
						-1._dbl/3._dbl, -1._dbl/3._dbl, 2._dbl/3._dbl, 0._dbl, 0._dbl, 0._dbl, &
						0._dbl, 0._dbl, 0._dbl, 1._dbl, 0._dbl, 0._dbl, &
						0._dbl, 0._dbl, 0._dbl, 0._dbl, 1._dbl, 0._dbl, &
						0._dbl, 0._dbl, 0._dbl, 0._dbl, 0._dbl, 1._dbl /), (/ 6, 6 /) )

	! dsigma_m/dsigma
	real(kind=dbl), dimension(6), parameter :: Dsigma_mDsigma = (/ 1._dbl/3._dbl, 1._dbl/3._dbl, 1._dbl/3._dbl, &
																	 0._dbl, 0._dbl, 0._dbl /)

	! numerical zero
	real(kind=dbl), parameter :: num_zero = 1.0d-10

	real(kind=dbl), dimension(6), parameter :: eps_eng_scale = (/ 1._dbl, 1._dbl, 1._dbl, 0.5_dbl, 0.5_dbl, 0.5_dbl /)	!  scaling vector for strain vector to convert engineering shear strain components

contains


	function calc_I1(sigma)
		implicit none
		! computes first invariant of stress tensor

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_I1

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector
		! --------------------------------------------------------------------------

		! computation of first invariant
		calc_I1 = sum(sigma(1:3))

	end function calc_I1




	function calc_I2(sigma)
		implicit none
		! computes second invariant of stress tensor

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_I2

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: s11
		real(kind=dbl) :: s22
		real(kind=dbl) :: s33
		real(kind=dbl) :: s12
		real(kind=dbl) :: s13
		real(kind=dbl) :: s23
		! --------------------------------------------------------------------------

		! unpacking of stress vector
		s11 = sigma(1)
		s22 = sigma(2)
		s33 = sigma(3)
		s12 = sigma(4)
		s13 = sigma(5)
		s23 = sigma(6)

		! compuation of second invariant
		calc_I2 = s11*s22 + s22*s33 + s11*s33 - s12**2 - s13**2 - s23**2

	end function calc_I2




	function calc_I3(sigma)
		implicit none
		! computes third invariant of stress tensor

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_I3

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: s11
		real(kind=dbl) :: s22
		real(kind=dbl) :: s33
		real(kind=dbl) :: s12
		real(kind=dbl) :: s13
		real(kind=dbl) :: s23
		! --------------------------------------------------------------------------

		! unpacking of stress vector
		s11 = sigma(1)
		s22 = sigma(2)
		s33 = sigma(3)
		s12 = sigma(4)
		s13 = sigma(5)
		s23 = sigma(6)

		! compuation of second invariant
		calc_I3 = s11*s22*s33 + 2._dbl*s12*s23*s13 - s12**2*s33 - s23**2*s11 - s13**2*s22

	end function calc_I3




	function calc_s_mean(sigma)
		implicit none
		! computes mean stress

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_s_mean

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector
		! --------------------------------------------------------------------------

		! computation of mean stress
		calc_s_mean = 1._dbl/3._dbl * sum( sigma(1:3) )

	end function calc_s_mean




	function calc_s_dev(sigma)
		implicit none
		! computes deviatoric stress vector

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6) :: calc_s_dev

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector

		! internal variable
		real(kind=dbl) :: s_mean								! mean stress
		! --------------------------------------------------------------------------

		! compuation of s_mean
		s_mean = calc_s_mean(sigma)

		! computation of s_dev
		calc_s_dev = sigma - (/ s_mean, s_mean, s_mean, 0._dbl, 0._dbl, 0._dbl /)

	end function calc_s_dev




	function calc_J1(s_dev)
		implicit none
		! computes first invariant of deviatoric stress vector

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_J1

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: s_dev		! deviatoric stress vector
		! --------------------------------------------------------------------------

		! computation of J1
		calc_J1 = sum( s_dev(1:3) )

	end function calc_J1




	function calc_J2(sigma)
		implicit none
		! computes second invariant of deviatoric stress vector

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_J2

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector
		! --------------------------------------------------------------------------

		! computation J2
		calc_J2 = 1._dbl/3._dbl*calc_I1(sigma)**2 - calc_I2(sigma)

		! test if J2 < 0 (J2 must be positive)
		if (calc_J2 < 0._dbl) then
			calc_J2 = 0._dbl
		end if

	end function calc_J2




	function calc_J3(sigma)
		implicit none
		! computes third invariant of deviatoric stress vector

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_J3

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector
		! --------------------------------------------------------------------------

		! computation J3
		calc_J3 = 2._dbl/27._dbl*calc_I1(sigma)**3 - 1._dbl/3._dbl*calc_I1(sigma)*calc_I2(sigma) + calc_I3(sigma)

	end function calc_J3



	function calc_DI1Dsigma(sigma)
		implicit none
	    ! computes the gradient of the first invariant of the stress vector with respect to sigma
	    ! dI1/dsigma

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6) :: calc_DI1Dsigma

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector
		! --------------------------------------------------------------------------

		! dI1/dsigma
		calc_DI1Dsigma = (/ 1._dbl, 1._dbl, 1._dbl, 0._dbl, 0._dbl, 0._dbl /)

	end function calc_DI1Dsigma



	function calc_DI2Dsigma(sigma)
		implicit none
	    ! computes the gradient of the second invariant of the stress vector with respect to sigma
	    ! dI2/dsigma

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6) :: calc_DI2Dsigma

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector
		! --------------------------------------------------------------------------

		! dI2/dsigma
		calc_DI2Dsigma = (/ sigma(2) + sigma(3), &
						sigma(1) + sigma(3), &
						sigma(2) + sigma(1), &
						-2._dbl*sigma(4), &
						-2._dbl*sigma(5), &
						-2._dbl*sigma(6) /)

	end function calc_DI2Dsigma



	function calc_DI3Dsigma(sigma)
		implicit none
	    ! computes the gradient of the third invariant of the stress vector with respect to sigma
	    ! dI3/dsigma

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6) :: calc_DI3Dsigma

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: s11
		real(kind=dbl) :: s22
		real(kind=dbl) :: s33
		real(kind=dbl) :: s12
		real(kind=dbl) :: s13
		real(kind=dbl) :: s23
		! --------------------------------------------------------------------------

		! unpacking of stress vector
		s11 = sigma(1)
		s22 = sigma(2)
		s33 = sigma(3)
		s12 = sigma(4)
		s13 = sigma(5)
		s23 = sigma(6)

		! dI3/dsigma
		calc_DI3Dsigma = (/ s22*s33 - s23**2, &
						s11*s33 - s13**2, &
						s11*s22 - s12**2, &
						2._dbl*(s23*s13 - s12*s33), &
						2._dbl*(s12*s23 - s13*s22), &
						2._dbl*(s12*s13 - s23*s11) /)

	end function calc_DI3Dsigma



	function calc_DJ2Dsigma(sigma)
		implicit none
	    ! computes the gradient of the second invariant of the deviatoric stress vector with respect to sigma
	    ! dJ2/dsigma

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6) :: calc_DJ2Dsigma

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: I1
		! --------------------------------------------------------------------------

		!return s_dev
		!calc_DJ2Dsigma = calc_s_dev(sigma)

		calc_DJ2Dsigma = 2._dbl/3._dbl*calc_I1(sigma)*calc_DI1Dsigma(sigma) - calc_DI2Dsigma(sigma)

	end function calc_DJ2Dsigma



	function calc_DJ3Dsigma(sigma)
		implicit none
	    ! computes the gradient if the third invariant of the deviatoric stress vector with respect to sigma
	    ! dJ3/dsigma

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6) :: calc_DJ3Dsigma

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: dJ3dI1		! dJ3/dI1
		real(kind=dbl) :: dJ3dI2		! dJ3/dI2
		real(kind=dbl) :: dJ3dI3		! dJ3/dI3
		! --------------------------------------------------------------------------

		! dJ3/dI1
		dJ3dI1 = 6._dbl/27._dbl*calc_I1(sigma)**2 - 1._dbl/3._dbl*calc_I2(sigma)

		! dJ3/dI2
		dJ3dI2 = -1._dbl/3._dbl*calc_I1(sigma)

		! dJ3/dI2
		dJ3dI3 = 1._dbl

	    ! dJ3/dsigma
	    calc_DJ3Dsigma = dJ3dI1*calc_DI1Dsigma(sigma) + dJ3dI2*calc_DI2Dsigma(sigma) + dJ3dI3*calc_DI3Dsigma(sigma)

	end function calc_DJ3Dsigma



	function calc_rho(sigma)
		implicit none
		! computes the Haigh Westergaard coordinate rho

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_rho

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: J2									! second invariant of deviatoric stress vector
		! --------------------------------------------------------------------------

		! computation of J2
		J2 = calc_J2(sigma)

		! computation of rho
		calc_rho = sqrt(2._dbl*J2)

	end function calc_rho



	function calc_theta(sigma)
		implicit none
		! computes the Haigh Westergaard coordinate theta

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl) :: calc_theta

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: J2									! second invariant of deviatoric stress vector
		real(kind=dbl) :: J3									! third invariant of deviatoric stress vector
		! --------------------------------------------------------------------------

		! computation of J2 and J3
		J2 = calc_J2(sigma)
		J3 = calc_J3(sigma)

		! computation of theta
		! test if both J2 and J3 are zero, or J2 is negative, i.e a stress state on the hydrostatic axis
		if ( (abs(J3) < epsilon(1._dbl) .and. sqrt(J2)**3 < epsilon(1._dbl) ) .or. &

				(J2 <= 0._dbl) ) then

			! set theta to zero
			calc_theta = 0._dbl

		else if ( 3._dbl*sqrt(3._dbl)/2._dbl * J3/sqrt(J2)**3 > 1._dbl ) then

			! set theta to zero
			calc_theta = 0._dbl

		else if ( 3._dbl*sqrt(3._dbl)/2._dbl * J3/sqrt(J2)**3 < -1._dbl ) then

			calc_theta = Pi/3._dbl

		else

			calc_theta = 1._dbl/3._dbl * acos( 3._dbl*sqrt(3._dbl)/2._dbl * J3/sqrt(J2)**3 )

		end if

	end function calc_theta



	function calc_DrhoDsigma(sigma)
		implicit none
	    ! computes the gradient of the Haigh Weestergaard coordinate rho with respect to sigma ( drho/dsigma )

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6) :: calc_DrhoDsigma

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: rho
		! --------------------------------------------------------------------------

		! compute rho
		rho = calc_rho(sigma)

		! droh/dsigma
		if (rho < epsilon(rho)) then	! numerically equivalent to zero
			calc_DrhoDsigma = huge(rho)*1.0D-10
		else
		    calc_DrhoDsigma = 1._dbl/rho * calc_DJ2Dsigma(sigma)
		end if

	end function calc_DrhoDsigma



	function calc_DthetaDsigma(sigma)
		implicit none
	    ! computes the gradient of the Haigh Weestergaard coordinate theta with respect to sigma ( dtheta/dsigma )

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6) :: calc_DthetaDsigma

		! passed variable
		real(kind=dbl), dimension(6), intent(in) :: sigma		! stress vector

		! internal variables
		real(kind=dbl) :: dtheta_dJ2				! dtheta/dJ2
		real(kind=dbl) :: dtheta_dJ3				! dtheta/dJ3
		! --------------------------------------------------------------------------

		! check if cos(3*theta) is outside the interval ]-1,1[
		if ( .not. (cos(3._dbl*calc_theta(sigma)) < 1._dbl) .or. &
			 .not. (cos(3._dbl*calc_theta(sigma)) > -1._dbl) ) then

			! cos(3*theta) is on the borders of or outside the interval,
			! set dtheta/dJ2 and dtheta/dJ3 to zero
			dtheta_dJ2 = 0.0_dbl
			dtheta_dJ3 = 0.0_dbl

		else

			! dtheta/dJ2
			dtheta_dJ2 = sqrt(27._dbl)/4._dbl * calc_J3(sigma)/ &
							( sqrt(calc_J2(sigma))**5 * sqrt(1._dbl - cos(3._dbl*calc_theta(sigma))**2) )

			! dtheta/dJ3
			dtheta_dJ3 = -sqrt(3._dbl)/2._dbl * 1._dbl/ &
							( sqrt(calc_J2(sigma))**3 * sqrt(1._dbl - cos(3._dbl*calc_theta(sigma))**2) )

		end if

	    ! DthetaDsigma
	    calc_DthetaDsigma = dtheta_dJ2*calc_DJ2Dsigma(sigma) + dtheta_dJ3*calc_DJ3Dsigma(sigma)

	end function calc_DthetaDsigma



	function calc_C_elastic(E,nu)
		implicit none
		! computes elasticity matrix according to the generalized hooke's law
		! sigma = C * epsilon
	    ! see "Concepts and Application of Finite Element Analysis" p. 79, Cook R. et al, 2002, John Wiley & Sons
	    ! engineering shear strains are assumed -> tau_xy = 2*epsilon_xy

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6,6) :: calc_C_elastic

		! passed variables
		real(kind=dbl), intent(in) :: E		! elasticity module
		real(kind=dbl), intent(in) :: nu		! poisson ratio

		! internal variables
		real(kind=dbl) :: a
		real(kind=dbl) :: b
		real(kind=dbl) :: c
		! --------------------------------------------------------------------------

		! initialization of internal variables
		a = E*(1._dbl-nu)/( (1._dbl+nu)*(1._dbl-2._dbl*nu) )
		b = nu/(1._dbl-nu)
		c = (1._dbl-2._dbl*nu)/( 2._dbl*(1._dbl-nu) )

		! initialization of C_el with zeros
		calc_C_elastic = 0._dbl

		! initialization of non zero components of C_el
		calc_C_elastic(1,1) = a
		calc_C_elastic(2,2) = a
		calc_C_elastic(3,3) = a
		calc_C_elastic(4,4) = a*c
		calc_C_elastic(5,5) = a*c
		calc_C_elastic(6,6) = a*c
		calc_C_elastic(1,2) = a*b
		calc_C_elastic(1,3) = a*b
		calc_C_elastic(2,1) = a*b
		calc_C_elastic(2,3) = a*b
		calc_C_elastic(3,1) = a*b
		calc_C_elastic(3,2) = a*b

	end function calc_C_elastic




	function calc_inverse_C_elastic(E,nu)
		implicit none
	    ! computes inverse elasticity matrix C^(-1) according to the generalized hooke's law
	    ! epsilon = C^(-1) * sigma
	    ! see "Festigkeitslehre" p. 86; Mang H. and G. Hofstetter; 2004, SpringerWienNewYork
	    ! engineering shear strains are assumed -> tau_xy = 2*epsilon_xy

		! --------------------------------------------------------------------------
		! type declaration of return variable
		real(kind=dbl), dimension(6,6) :: calc_inverse_C_elastic

		! passed variables
		real(kind=dbl), intent(in) :: E		! elasticity module
		real(kind=dbl), intent(in) :: nu		! poisson ratio

		! internal variables
		real(kind=dbl) :: c
		real(kind=dbl) :: G
		! --------------------------------------------------------------------------

		! initialization of internal variables
		c = -nu/E
		G = E/( 2._dbl*(1._dbl+nu) )

		! initialization of C_el with zeros
		calc_inverse_C_elastic = 0._dbl

		! initialization of non zero components of C_el
		calc_inverse_C_elastic(1,1) = 1._dbl/E
		calc_inverse_C_elastic(2,2) = 1._dbl/E
		calc_inverse_C_elastic(3,3) = 1._dbl/E
		calc_inverse_C_elastic(4,4) = 1._dbl/G
		calc_inverse_C_elastic(5,5) = 1._dbl/G
		calc_inverse_C_elastic(6,6) = 1._dbl/G
		calc_inverse_C_elastic(1,2) = c
		calc_inverse_C_elastic(1,3) = c
		calc_inverse_C_elastic(2,1) = c
		calc_inverse_C_elastic(2,3) = c
		calc_inverse_C_elastic(3,1) = c
		calc_inverse_C_elastic(3,2) = c

	end function calc_inverse_C_elastic


! CORRECTED TRANSVERSE ELASTICITY MATRIX, ERROR IN WITTKE'S BOOK !!!!!!
	function calc_transv_istr_C_elastic(E1, E2, G2, nu1, nu2)
		! computes transverse isotropic elasticity matrix,
		! defined according to 'Rock Mechanics Based on an Anisotropic Jointed Rock Model' by Walter Wittke, p. 42 ff
		! C_el(1,2) = C_el(2,1) is wrong in the book, here the corrected expressions are implemented
		! C_el(3,3) is wrong in the book, here the corrected expressions are implemented
		! sigma = C_el * epsilon
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in) :: E1		! elasticity constant parallel to the structure (isotropic) plane
		real(kind=dbl), intent(in) :: E2		! elasticity constant perpendicular to the structure (isotropic) plane
		real(kind=dbl), intent(in) :: G2		! shear modulus for shear loading in/parallel to the structure (isotropic) plane
		real(kind=dbl), intent(in) :: nu1		! poisson ration perpendicular to the structure (isotropic) plane
		real(kind=dbl), intent(in) :: nu2		! poisson ration in the structure (isotropic) plane

		! return variable
		real(kind=dbl), dimension(6,6) :: calc_transv_istr_C_elastic		! return variable
!f2py	intent(out) :: calc_transv_istr_C_elastic

		! internal variables
		real(kind=dbl) :: n		! auxiliary variable
		real(kind=dbl) :: m		! auxiliary variable
		! --------------------------------------------------------------------------

		! initialize auxiliary variables
		n = E1/E2
		m = 1._dbl-nu1-2._dbl*n*nu2**2

		! initialize elasticity matrix with zeros
		calc_transv_istr_C_elastic = 0._dbl

		! compute nonzero components of elasticity matrix
		calc_transv_istr_C_elastic(1,1) = E1*(1._dbl-n*nu2**2)/((1._dbl+nu1)*m)
		calc_transv_istr_C_elastic(1,2) = E1*(nu1+n*nu2**2)/((1._dbl+nu1)*m)		! E1*(1._dbl+n*nu2**2)/((1._dbl+nu1)*m) in Wittke's book which is wrong
		calc_transv_istr_C_elastic(2,1) = E1*(nu1+n*nu2**2)/((1._dbl+nu1)*m)		! E1*(1._dbl+n*nu2**2)/((1._dbl+nu1)*m) in Wittke's book which is wrong
		calc_transv_istr_C_elastic(2,2) = E1*(1._dbl-n*nu2**2)/((1._dbl+nu1)*m)
		calc_transv_istr_C_elastic(1:2,3) = E1*nu2/m
		calc_transv_istr_C_elastic(3,1:2) = E1*nu2/m
		calc_transv_istr_C_elastic(3,3) = E2*(1._dbl-nu1)/m							! E1*(1._dbl-nu1)/m in Wittke's book which is wrong
		calc_transv_istr_C_elastic(4,4) = E1/(2._dbl*(1._dbl+nu1))
		calc_transv_istr_C_elastic(5,5) = G2
		calc_transv_istr_C_elastic(6,6) = G2

	end function calc_transv_istr_C_elastic



!	function calc_transv_istr_C_elastic(E1, E2, G2, nu1, nu2)
!		! computes transverse isotropic elasticity matrix,
!		! defined according to Maple Computation
!		! sigma = C_el * epsilon
!		implicit none
!
!		! --------------------------------------------------------------------------
!		! passed variables
!		real(kind=dbl), intent(in) :: E1		! elasticity constant parallel to the structure (isotropic) plane
!		real(kind=dbl), intent(in) :: E2		! elasticity constant perpendicular to the structure (isotropic) plane
!		real(kind=dbl), intent(in) :: G2		! shear modulus for shear loading in/parallel to the structure (isotropic) plane
!		real(kind=dbl), intent(in) :: nu1		! poisson ration perpendicular to the structure (isotropic) plane
!		real(kind=dbl), intent(in) :: nu2		! poisson ration in the structure (isotropic) plane
!
!		! return variable
!		real(kind=dbl), dimension(6,6) :: calc_transv_istr_C_elastic		! return variable
!
!		! internal variables
!		! --------------------------------------------------------------------------
!
!		! initialize elasticity matrix with zeros
!		calc_transv_istr_C_elastic = 0._dbl
!
!		! compute nonzero components of elasticity matrix
!		calc_transv_istr_C_elastic(1,1) = E1*(-E2+nu2**2*E1)/((E2*nu1-E2+2._dbl*nu2**2*E1)*(1._dbl+nu1))
!		calc_transv_istr_C_elastic(1,2) = -E1*(E2*nu1+nu2**2*E1)/((E2*nu1-E2+2._dbl*nu2**2*E1)*(1._dbl+nu1))
!		calc_transv_istr_C_elastic(1,3) = -nu2*E1*E2/(E2*nu1-E2+2*nu2**2*E1)
!		calc_transv_istr_C_elastic(2,1) = -E1*(E2*nu1+nu2**2*E1)/((E2*nu1-E2+2._dbl*nu2**2*E1)*(1._dbl+nu1))
!		calc_transv_istr_C_elastic(2,2) = E1*(-E2+nu2**2*E1)/((E2*nu1-E2+2._dbl*nu2**2*E1)*(1._dbl+nu1))
!		calc_transv_istr_C_elastic(2,3) = -nu2*E1*E2/(E2*nu1-E2+2._dbl*nu2**2*E1)
!		calc_transv_istr_C_elastic(3,1) = -nu2*E1*E2/(E2*nu1-E2+2._dbl*nu2**2*E1)
!		calc_transv_istr_C_elastic(3,2) = -nu2*E1*E2/(E2*nu1-E2+2._dbl*nu2**2*E1)
!		calc_transv_istr_C_elastic(3,3) = (nu1-1._dbl)*E2**2/(E2*nu1-E2+2._dbl*nu2**2*E1)
!		calc_transv_istr_C_elastic(4,4) = 1._dbl/2._dbl*E1/(1._dbl+nu1)
!		calc_transv_istr_C_elastic(5,5) = G2
!		calc_transv_istr_C_elastic(6,6) = G2
!
!	end function calc_transv_istr_C_elastic



	function calc_inv_transv_istr_C_elastic(E1, E2, G2, nu1, nu2)
		! computes transverse isotropic elasticity matrix,
		! defined according to 'Rock Mechanics Based on an Anisotropic Jointed Rock Model' by Walter Wittke, p. 42 ff
		! epsilon = C_el^(⁻1) * sigma
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in) :: E1		! elasticity constant parallel to the structure (isotropic) plane
		real(kind=dbl), intent(in) :: E2		! elasticity constant perpendicular to the structure (isotropic) plane
		real(kind=dbl), intent(in) :: G2		! shear modulus for shear loading in/parallel to the structure (isotropic) plane
		real(kind=dbl), intent(in) :: nu1		! poisson ration perpendicular to the structure (isotropic) plane
		real(kind=dbl), intent(in) :: nu2		! poisson ration in the structure (isotropic) plane

		! return variable
		real(kind=dbl), dimension(6,6) :: calc_inv_transv_istr_C_elastic		! return variable

		! internal variables
		! --------------------------------------------------------------------------

		! initialize elasticity matrix with zeros
		calc_inv_transv_istr_C_elastic = 0._dbl

		! compute nonzero components of elasticity matrix
		calc_inv_transv_istr_C_elastic(1,1) = 1._dbl/E1
		calc_inv_transv_istr_C_elastic(1,2) = -nu1/E1
		calc_inv_transv_istr_C_elastic(2,1) = -nu1/E1
		calc_inv_transv_istr_C_elastic(2,2) = 1._dbl/E1
		calc_inv_transv_istr_C_elastic(1:2,3) = -nu2/E2
		calc_inv_transv_istr_C_elastic(3,1:2) = -nu2/E2
		calc_inv_transv_istr_C_elastic(3,3) = 1._dbl/E2
		calc_inv_transv_istr_C_elastic(4,4) = 2._dbl*(1._dbl+nu1)/E1
		calc_inv_transv_istr_C_elastic(5,5) = 1._dbl/G2
		calc_inv_transv_istr_C_elastic(6,6) = 1._dbl/G2

	end function calc_inv_transv_istr_C_elastic



	function calc_sig_princ_hw(s_mean, rho, theta)
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl), dimension(3) :: calc_sig_princ_hw

		! passed variables
		real(kind=dbl), intent(in) :: s_mean
		real(kind=dbl), intent(in) :: rho
		real(kind=dbl), intent(in) :: theta

		! internal variables
		real(kind=dbl), dimension(3) :: s_mean_vec
		real(kind=dbl), dimension(3) :: theta_vec
		! --------------------------------------------------------------------------

		s_mean_vec = (/ s_mean, s_mean, s_mean /)

		theta_vec = (/ cos(theta), -sin(Pi/6._dbl-theta), -sin(Pi/6._dbl+theta) /)

		calc_sig_princ_hw = s_mean_vec + sqrt(2._dbl/3._dbl)*rho*theta_vec

	end function calc_sig_princ_hw



	function calc_sig_princ_hw_sig(sigma)
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl), dimension(3) :: calc_sig_princ_hw_sig

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma

		! internal variables
		real(kind=dbl) :: s_mean
		real(kind=dbl) :: rho
		real(kind=dbl) :: theta
		real(kind=dbl), dimension(3) :: s_mean_vec
		real(kind=dbl), dimension(3) :: theta_vec
		! --------------------------------------------------------------------------

		! test if shear stress components are zero
		if ( (abs(sigma(4)) < epsilon(1._dbl)) .and. &
				(abs(sigma(5)) < epsilon(1._dbl)) .and. &
				(abs(sigma(6)) < epsilon(1._dbl)) ) then

			! shear stresses are zero and therefore
			! the components of sigma already resemble
			! the principal stresses
			calc_sig_princ_hw_sig = sort_sig_princ(sigma(1:3))
			return

		end if

		s_mean = calc_s_mean(sigma)
		rho = calc_rho(sigma)
		theta = calc_theta(sigma)

		s_mean_vec = (/ s_mean, s_mean, s_mean /)

		theta_vec = (/ cos(theta), -sin(Pi/6._dbl-theta), -sin(Pi/6._dbl+theta) /)

		calc_sig_princ_hw_sig = s_mean_vec + sqrt(2._dbl/3._dbl)*rho*theta_vec

	end function calc_sig_princ_hw_sig



	function calc_Dsigma_princDsigma(sigma)
		! computes derivatives of principal stress vector
		! with respect to the stress vector via HW-coordinates
		! dsigma_princ/dsigma
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl), dimension(3,6) :: calc_Dsigma_princDsigma	! dsigma_princ/dsigma

		! passed variables
		real(kind=dbl), intent(in), dimension(6) :: sigma			! stress vector

		! internal variables
		real(kind=dbl) :: s_mean									! HW-coordinate
		real(kind=dbl) :: rho										! HW-coordinate
		real(kind=dbl) :: theta										! HW-coordinate
		real(kind=dbl), dimension(3,1) :: Dsig_princDs_mean			! dsigma_princ/ds_mean initialized as rank 2 variable (matrix with 1 column),
																	! so that it can be transposed and used in matrix multiplications
																	! (not possible with rank 1 arrays, i.e. vectors)

		real(kind=dbl), dimension(3,1) :: Dsig_princDrho			! dsigma_princ/drho initialized as rank 2 variable (matrix with 1 column),
																	! so that it can be transposed and used in matrix multiplications
																	! (not possible with rank 1 arrays, i.e. vectors)

		real(kind=dbl), dimension(3,1) :: Dsig_princDtheta			! dsigma_princ/dtheta initialized as rank 2 variable (matrix with 1 column),
																	! so that it can be transposed and used in matrix multiplications
																	! (not possible with rank 1 arrays, i.e. vectors)

		real(kind=dbl), dimension(6,1) :: Ds_meanDsig				! ds_mean/dsigma initialized as rank 2 variable (matrix with 1 column),
																	! so that it can be transposed and used in matrix multiplications
																	! (not possible with rank 1 arrays, i.e. vectors)

		real(kind=dbl), dimension(6,1) :: DrhoDsig					! drho/dsigma initialized as rank 2 variable (matrix with 1 column),
																	! so that it can be transposed and used in matrix multiplications
																	! (not possible with rank 1 arrays, i.e. vectors)

		real(kind=dbl), dimension(6,1) :: dthetaDsig				! dtheta/dsigma initialized as rank 2 variable (matrix with 1 column),
																	! so that it can be transposed and used in matrix multiplications
																	! (not possible with rank 1 arrays, i.e. vectors)
		! --------------------------------------------------------------------------

		! computation of HW-coordinates
		s_mean = calc_s_mean(sigma)
		rho = calc_rho(sigma)
		theta = calc_theta(sigma)

		! computation of dsigma_princ/ds_mean
		Dsig_princDs_mean = reshape( (/ 1._dbl, 1._dbl, 1._dbl /), (/3,1/) )

		! computation of dsigma_princ/drho
		Dsig_princDrho = sqrt(2._dbl/3._dbl) * &
							reshape( (/ cos(theta), -sin(Pi/6._dbl-theta), -sin(Pi/6._dbl+theta) /), &
									(/3,1/) )

		! computation of dsigma_princ/dtheta
		Dsig_princDtheta = sqrt(2._dbl/3._dbl) * rho * &
							reshape( (/ -sin(theta), cos(Pi/6._dbl-theta), -cos(Pi/6._dbl+theta) /), &
									(/3,1/) )

		! computation of ds_mean/dsigma and assignment to
		! rank 2 (matrix) variable Ds_meanDsig
		Ds_meanDsig = reshape( Dsigma_mDsigma, (/6,1/) )

		! computation of drho/dsigma and assignment to
		! rank 2 (matrix) variable DrhoDsig
		DrhoDsig = reshape( calc_DrhoDsigma(sigma), (/6,1/) )

		! computation of dtheta/dsigma and assignment to
		! rank 2 (matrix) variable DthetaDsig
		DthetaDsig = reshape( calc_DthetaDsigma(sigma), (/6,1/) )

		! computation of dsigma_princ/dsigma with chain rule
		! dsigma_princ/dsigma = (dsigma_princ/ds_mean * ds_mean/dsigma) +
		! 						+ (dsigma_princ/drho * drho/dsigma) +
		! 						+ (dsigma_princ/dtheta * dtheta/dsigma)
		calc_Dsigma_princDsigma = matmul( Dsig_princDs_mean, transpose(Ds_meanDsig) ) + &
									matmul( Dsig_princDrho, transpose(DrhoDsig) ) + &
									matmul( Dsig_princDtheta, transpose(DthetaDsig) )

	end function calc_Dsigma_princDsigma



	function sort_sig_princ(sig_princ)
		! sorts the entries of the principal stress vector,
		! depending on their size in descending order
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl), dimension(3) :: sort_sig_princ

		! input variable
		real(kind=dbl), intent(in), dimension(3) :: sig_princ

		! internal variables
		integer, dimension(1) :: ind_max	! index of maximum sigma component
		integer, dimension(1) :: ind_min	! index of minimum sigma component
		integer, dimension(1) :: ind_mid	! index of intermediate sigma component
		! --------------------------------------------------------------------------

		! find index of maximum value
		ind_max = maxloc(sig_princ)
		! find index of minimum value
		ind_min = minloc(sig_princ)

		if (ind_max(1) == ind_min(1)) then
			! all components have the same size
			! return principal stress vector in unchanged order
			sort_sig_princ = sig_princ
			return
		else
			! find index of intermediate value
			ind_mid = 6 - ind_max - ind_min

			! assign maximum value to first entry of array
			sort_sig_princ(1) = sig_princ(ind_max(1))
			! assign intermediate value to second entry of array
			sort_sig_princ(2) = sig_princ(ind_mid(1))
			! assign minimum value to third(last) entry of array
			sort_sig_princ(3) = sig_princ(ind_min(1))
		end if

	end function sort_sig_princ



	function calc_rot_matrix_epsilon(alpha, beta, gam)
		! computes rotation matrix for coordinate transformation of the strain vector,
		! with shear strains in engineering notation
		! according to 'Concepts and Applications of Finite Element Analysis' by Robert D. Cook et al, p. 274, ff
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in) :: alpha		! rotation around z axis
		real(kind=dbl), intent(in) :: beta		! rotation around rotated y axis
		real(kind=dbl), intent(in) :: gam		! rotation around rotated z axis, in case if transverse isotropy zero

		! output variable
		real(kind=dbl), dimension(6,6) :: calc_rot_matrix_epsilon	! output rotation matrix

		! internal variables
		real(kind=dbl), dimension(3,3) :: T						! rotation matrix with direction cosines of angles between original and rotated coordinate axes
		real(kind=dbl) :: lx, mx, nx							! components of rotated x-base vector (first line of T)
		real(kind=dbl) :: ly, my, ny							! components of rotated y-base vector (second line of T)
		real(kind=dbl) :: lz, mz, nz							! components of rotated z-base vector (third line of T)
		real(kind=dbl), dimension(3,3) :: T11, T12, T21, T22	! submatrices of rotation matrix
		! --------------------------------------------------------------------------

		! computation of rotation matrix out of the three euler rotation angle defined by zyz convention
		T = calc_rot_matrix_zyz(alpha, beta, gam)

		! assign values of T to components
		lx = T(1,1)
		mx = T(1,2)
		nx = T(1,3)
		ly = T(2,1)
		my = T(2,2)
		ny = T(2,3)
		lz = T(3,1)
		mz = T(3,2)
		nz = T(3,3)

		! fill submatrices
		T11 = T**2
		T12(:,1) = T(:,1)*T(:,2)
		T12(:,2) = T(:,1)*T(:,3)
		T12(:,3) = T(:,2)*T(:,3)
		T21(1,:) = 2._dbl*T(1,:)*T(2,:)
		T21(2,:) = 2._dbl*T(1,:)*T(3,:)
		T21(3,:) = 2._dbl*T(2,:)*T(3,:)
		T22 = reshape((/ (lx*my+mx*ly), (lx*mz+mx*lz), (ly*mz+my*lz), &
							(lx*ny+nx*ly), (lx*nz+nx*lz), (ly*nz+ny*lz), &
							(mx*ny+nx*my), (mx*nz+nx*mz), (my*nz+ny*mz) /), (/3,3/))

		! fill rotation matrix with submatrices
		calc_rot_matrix_epsilon(1:3,1:3) = T11
		calc_rot_matrix_epsilon(1:3,4:6) = T12
		calc_rot_matrix_epsilon(4:6,1:3) = T21
		calc_rot_matrix_epsilon(4:6,4:6) = T22

	end function calc_rot_matrix_epsilon



	function calc_rot_matrix_sigma(alpha, beta, gam)
		! computes rotation matrix for coordinate transformation of the stress vector,
		! according to 'Concepts and Applications of Finite Element Analysis' by Robert D. Cook et al, p. 274, ff
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in) :: alpha		! rotation around z axis
		real(kind=dbl), intent(in) :: beta		! rotation around rotated y axis
		real(kind=dbl), intent(in) :: gam		! rotation around rotated z axis, in case if transverse isotropy zero

		!output variable
		real(kind=dbl), dimension(6,6) :: calc_rot_matrix_sigma		! output stress rotation matrix

		! internal variables
		real(kind=dbl), dimension(6,6) :: T_epsilon		! rotation matrix for strain vector
		! --------------------------------------------------------------------------

		! compute rotation matrix for strain vector
		T_epsilon = calc_rot_matrix_epsilon(alpha, beta, gam)

		! change submatrices of T_epsilon and assign them to output stress rotation matrix
		calc_rot_matrix_sigma(1:3,1:3) = T_epsilon(1:3,1:3)
		calc_rot_matrix_sigma(1:3,4:6) = 2._dbl*T_epsilon(1:3,4:6)
		calc_rot_matrix_sigma(4:6,1:3) = 0.5_dbl*T_epsilon(4:6,1:3)
		calc_rot_matrix_sigma(4:6,4:6) = T_epsilon(4:6,4:6)

	end function calc_rot_matrix_sigma



	function rotate_sigma(sigma, alpha, beta, gam)
		! computes and returns rotated stress vector
		implicit none

		! --------------------------------------------------------------------------
		! input variables
		real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector
		real(kind=dbl), intent(in) :: alpha						! first euler rotation angle, around z-axis, eastwards positive (zyz convention)
		real(kind=dbl), intent(in) :: beta						! second euler rotation angle, around rotated y-axis, downwards positive (zyz convention)
		real(kind=dbl), intent(in) :: gam						! third euler rotation angle, around rotated z-axis, eastwards positive (zyz convention)

		! return variable
		real(kind=dbl), dimension(6) :: rotate_sigma			! rotated stress vector

		! internal variables
		real(kind=dbl), dimension(6,6) :: T_sig					! rotation matrix for rotation of stress vector
		! --------------------------------------------------------------------------

		! computation of stress rotation matrix with the three euler rotation angle defined by zyz convention
		T_sig = calc_rot_matrix_sigma(alpha, beta, gam)

		! compute rotated stress vector
		rotate_sigma = matmul(T_sig, sigma)

	end function rotate_sigma



	function rotate_epsilon(eps, alpha, beta, gam)
		! computes and returns rotated strain vector
		implicit none

		! --------------------------------------------------------------------------
		! input variables
		real(kind=dbl), intent(in), dimension(6) :: eps			! strain vector
		real(kind=dbl), intent(in) :: alpha						! first euler rotation angle, around z-axis, eastwards positive (zyz convention)
		real(kind=dbl), intent(in) :: beta						! second euler rotation angle, around rotated y-axis, downwards positive (zyz convention)
		real(kind=dbl), intent(in) :: gam						! third euler rotation angle, around rotated z-axis, eastwards positive (zyz convention)

		! return variable
		real(kind=dbl), dimension(6) :: rotate_epsilon			! rotated strain vector

		! internal variables
		real(kind=dbl), dimension(6,6) :: T_eps					! rotation matrix for rotation of strain vector
		! --------------------------------------------------------------------------

		! computation of strain rotation matrix with the three euler rotation angle defined by zyz convention
		T_eps = calc_rot_matrix_epsilon(alpha, beta, gam)

		! compute rotated stress vector
		rotate_epsilon = matmul(T_eps, eps)

	end function rotate_epsilon



	function rotate_material_matrix(C, alpha, beta, gam)
		! computes and returns rotated material tensor/matrix (C_elastic or C_elastoplastic),
		! according to: 'Rock Mechanics Based on an Anisotropic Jointed Rock Model' by Walter Wittke, p. 45 ff
		!				'Concepts and Applications of Finite Element Analysis' by Robert D. Cook et al, p. 275, ff
		implicit none

		! --------------------------------------------------------------------------
		! input variables
		real(kind=dbl), intent(in), dimension(6,6) :: C			! material matrix
		real(kind=dbl), intent(in) :: alpha						! first euler rotation angle, around z-axis, eastwards positive (zyz convention)
		real(kind=dbl), intent(in) :: beta						! second euler rotation angle, around rotated y-axis, downwards positive (zyz convention)
		real(kind=dbl), intent(in) :: gam						! third euler rotation angle, around rotated z-axis, eastwards positive (zyz convention)

		! return variable
		real(kind=dbl), dimension(6,6) :: rotate_material_matrix			! rotated material matrix
!f2py	intent(out) :: rotate_material_matrix

		! internal variables
		real(kind=dbl), dimension(6,6) :: T_eps					! strain rotation matrix for rotation of material matrix
		real(kind=dbl), dimension(6,6) :: T_sigma				! stress rotation matrix for rotation of material matrix
		! --------------------------------------------------------------------------

		! computation of strain rotation matrix
		T_eps = calc_rot_matrix_epsilon(-alpha, -beta, -gam)

		! computation of stress rotation matrix
		T_sigma = calc_rot_matrix_sigma(alpha, beta, gam)

		! compute rotated material matrix
		rotate_material_matrix = matmul(T_sigma, matmul(C,T_eps))

	end function rotate_material_matrix



	function rotate_inv_material_matrix(C_inv, alpha, beta, gam)
		! computes and returns rotated material tensor/matrix (inverse C_elastic or inverse C_elastoplastic),
		! according to: 'Rock Mechanics Based on an Anisotropic Jointed Rock Model' by Walter Wittke, p. 45 ff
		!				'Concepts and Applications of Finite Element Analysis' by Robert D. Cook et al, p. 275, ff
		implicit none

		! --------------------------------------------------------------------------
		! input variables
		real(kind=dbl), intent(in), dimension(6,6) :: C_inv			! material matrix
		real(kind=dbl), intent(in) :: alpha						! first euler rotation angle, around z-axis, eastwards positive (zyz convention)
		real(kind=dbl), intent(in) :: beta						! second euler rotation angle, around rotated y-axis, downwards positive (zyz convention)
		real(kind=dbl), intent(in) :: gam						! third euler rotation angle, around rotated z-axis, eastwards positive (zyz convention)

		! return variable
		real(kind=dbl), dimension(6,6) :: rotate_inv_material_matrix			! rotated inverse material matrix

		! internal variables
		real(kind=dbl), dimension(6,6) :: T_sig					! rotation matrix for rotation of inverse material matrix
		real(kind=dbl), dimension(6,6) :: T_eps					! rotation matrix for rotation of inverse material matrix
		! --------------------------------------------------------------------------

		! computation of rotation matrix
		!T_sig = calc_rot_matrix_sigma(alpha, beta, gam)
		T_eps = calc_rot_matrix_epsilon(alpha,beta, gam)

		! compute rotated inverse material matrix
		!rotate_inv_material_matrix = matmul(transpose(T_sig), matmul(C_inv,T_sig))
		rotate_inv_material_matrix = matmul(T_eps, matmul(C_inv,transpose(T_eps)))

	end function rotate_inv_material_matrix



	function calc_gen_l_vec(sigma)
		! computes the normed generalized loading vector l
		! according to "On failure criteria for anisotropic cohesive-frictional materials" by S. Pietruszczak and Z. Mroz, 2000
		implicit none

		! --------------------------------------------------------------------------
		! input variables
		real(kind=dbl), intent(in), dimension(6) :: sigma			! stress vector

		! return variable
		real(kind=dbl), dimension(3) :: calc_gen_l_vec		! generalized loading vector
!f2py	intent(out) :: calc_gen_loading_vec

		! internal variables
		real(kind=dbl), dimension(3) :: L		! loading vector
		! --------------------------------------------------------------------------

		! compute L
!		L(1) = sqrt(sigma(1)**2+sigma(4)**2+sigma(5)**2)		! s11, s12, s13
!		L(2) = sqrt(sigma(2)**2+sigma(4)**2+sigma(6)**2)		! s22, s12, s23
!		L(3) = sqrt(sigma(3)**2+sigma(5)**2+sigma(6)**2)		! s33, s13, s23
		if (.not. all(sigma .eq. 0._dbl)) then
			L = calc_gen_l_vec_unnormed(sigma)			! not all stress components are zero, proceed to compute L

			! compute generalized loading vector
			calc_gen_l_vec = 1._dbl/norm2(L) * L
		else

			calc_gen_l_vec = 0._dbl		! all stress components are zero, set l to zero and return

		end if

	end function calc_gen_l_vec



	function calc_gen_l_vec_unnormed(sigma)
	        ! computes unnormed generalized loading vector L
	        implicit none

	        ! --------------------------------------------------------------------------
	        ! passed variables
			real(kind=dbl), intent(in), dimension(6) :: sigma		! stress vector

	        ! return variable
	        real(kind=dbl), dimension(3) :: calc_gen_l_vec_unnormed	! L
! f2py		intent(out) :: calc_gen_loading_vec_unnormed

	        ! internal variables
	        ! --------------------------------------------------------------------------

		! compute L
		calc_gen_l_vec_unnormed(1) = sqrt(sigma(1)**2+sigma(4)**2+sigma(5)**2)		! s11, s12, s13
		calc_gen_l_vec_unnormed(2) = sqrt(sigma(2)**2+sigma(4)**2+sigma(6)**2)		! s22, s12, s23
		calc_gen_l_vec_unnormed(3) = sqrt(sigma(3)**2+sigma(5)**2+sigma(6)**2)		! s33, s13, s23

	    end function calc_gen_l_vec_unnormed


end module stress_routines
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
module Df_allDs
	! contains functions that compute the df/dsigma vectors of the material model
	! and one main (model-indepenent) function that calls all df/dsigma functions of the model
	! and returns the function values as a matrix

	! load additional modules
	use constants
	use material_info
	use derived_types
	use stress_routines
    implicit none

contains

	function DfDs(sigma, alpha)
		! computes the derivatives (gradient vectors) of the yield functions
		! with respect to sigma
		! in case of multi surface plasticity the gradients of all yield functions
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha		! internal variables vector

		! return variable
		real(kind=dbl), dimension(6,n_yf) :: 						DfDs	! results array

		! internal variable
		! --------------------------------------------------------------------------

		! computation of all gradient vector (modified leon, ...)
		! returning hardening function values in hardening_function vector

		DfDs(:,1) = Df_1Ds(sigma, alpha)

	end function DfDs




	function Df_1Ds(sigma, alpha)
		! compute df/dsigma vector of first yield function
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl), dimension(6) :: 				Df_1Ds

		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha		! internal variables vector

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function Df_1Ds



end module Df_allDs
module Df_allDa
	! contains functions that compute the df/dalpha vectors of the material model
	! and one main (model-indepenent) function that calls all df/dalpha functions of the model
	! and returns the function values as a matrix

	! load additional modules
	use constants
	use material_info
	use derived_types
	use stress_routines
    implicit none

contains

	function DfDa(sigma, alpha)
		! computes the derivatives (gradient vectors) of the yield functions
		! with respect to the alpha vector
		! in case of multi surface plasticity the gradients of all yield functions
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha		! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha),n_yf) :: 			DfDa	! results array

		! --------------------------------------------------------------------------


		! --------------------------------------
		! computation of all gradient vectors
		! returning hardening function values in hardening_function vector
		! --------------------------------------

		DfDa(:,1) = Df_1Da(sigma, alpha)

	end function DfDa



	function Df_1Da(sigma, alpha)
		! gradient vector of first yield function

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha,1)) :: Df_1Da

		! internal variables
		! --------------------------------------------------------------------------

		! ...


	end function Df_1Da



end module Df_allDa
module g_all
	! contains functions that compute the plastic potential functions of the material model
	! and one main (model-indepenent) function that calls all plastic potential functions of the model
	! and returns the function values as a vector

	use constants
	use material_info
	use derived_types

    implicit none

contains

	function g(sigma, alpha)
		! computes all plastic potential values for the given model
		! and returns them as real 1D array
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector

		! return variable
		real(kind=dbl), dimension(n_yf) :: g			! array for yield function values

		! internal variables
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of all plastic potentials (modified leon, ...)
		! returning plastic potential values in g vector
		! --------------------------------------

		g(1) = g_1(sigma, alpha)

	end function g



	function g_1(sigma, alpha)
		! computes first plastic potential function
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl) :: g_1

		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function g_1



end module g_all
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
module Dh_allDs
	! contains functions that compute the dh/dsigma matrices of the material model
	! and one main (model-indepenent) function that calls all dh/dsigma functions of the model
	! and returns the function values as a matrix

	use constants
	use material_info
	use derived_types

    implicit none

contains


	function DhDs(sigma, alpha)
		! computes the dh/dsigma matrix
		! in case of multi surface plasticity the dh/dsigma matrices of all yield surfaces
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha),6,n_yf) :: DhDs				! results array of rank 3

		! internal variable
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of all dh/dsigma matrices (modified leon, ...)
		! returning dh/dsigma matrices return array
		! --------------------------------------

		DhDs(:,:,1) = Dh_1Ds(sigma, alpha)

	end function DhDs




	function Dh_1Ds(sigma, alpha)

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma						! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha						! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha,1),6) :: Dh_1Ds

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function Dh_1Ds



end module Dh_allDs
module Dh_allDa
	! contains functions that compute the dh/dalpha matrices of the material model
	! and one main (model-indepenent) function that calls all dh/dalpha functions of the model
	! and returns the function values as a matrix

	! load additional modules
	use constants
	use material_info
	use derived_types

    implicit none

contains


	function DhDa(sigma, alpha)
		! computes the dh/dalpha matrix
		! in case of multi surface plasticity the dh/dalpha matrices of all yield surfaces
		! of the model must be computed in this subroutine
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 					sigma		! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 					alpha		! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha),size(alpha), n_yf) :: 	DhDa	! results array of rank 3

		! internal variable
		! --------------------------------------------------------------------------

		! --------------------------------------
		! computation of all dh/dalpha matrices (modified leon, ...)
		! returning dh/dsigma matrices return array
		! --------------------------------------

		DhDa(:,:,1) = Dh_1Da(sigma, alpha)

	end function DhDa




	function Dh_1Da(sigma, alpha)
		! compute dh/dalpha matrix for the first hardening function
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: 				sigma				! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 				alpha				! internal variables vector

		! return variable
		real(kind=dbl), dimension(size(alpha),size(alpha)) :: 	Dh_1Da

		! internal variables
		! --------------------------------------------------------------------------

		! ...

	end function Dh_1Da


end module Dh_allDa
module model_comp_act
	! modules that contains functions which filter out the inactive model components
	! ( yield function values, derivative vectors, derivative matrices) and return only the active components
	! as 2D arrays

	! load additional modules
	use constants
	use material_info

	use yf_all, only: yf_ret_map
	use Df_allDs, only: DfDs
	use Df_allDa, only: DfDa
	use Dg_allDs, only: DgDs
	use DDg_allDDs, only: DDgDDs
	use DDg_allDsDa, only: DDgDsDa
	use hf_all, only: hf
	use Dh_allDs, only: DhDs
	use Dh_allDa, only: DhDa

    implicit none

contains

	function yf_act(sigma, alpha, yield_flags, act_yf_ind)
		! returns array/vector with active yield function values

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(:), intent(in) :: 	sigma			! stress vector
		real(kind=dbl), dimension(:), intent(in) :: 	alpha			! internal variables vector
		logical, dimension(:), intent(in) :: 			yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(act_yf_ind)) :: 			yf_act		! active yield functions vector

		! internal variables
		real(kind=dbl), dimension(size(yield_flags)) ::	 yf_all		! all yield functions vector
		! --------------------------------------------------------------------------

		! compute all yield functions
		yf_all = yf_ret_map(sigma, alpha, yield_flags)

		! filter active yield functions via vector with indices of active yield functions
		yf_act = yf_all(act_yf_ind)

	end function yf_act



	function Df_actDsigma(sigma, alpha, yield_flags, act_yf_ind)
		! computes array with df_act/dsigma vectors of active yield functions

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(6,size(act_yf_ind)) :: Df_actDsigma				! array with df/dsigma vectors of active yield functions

		! internal variables
		real(kind=dbl), dimension(6,n_yf) :: Df_all			! array with df/dsigma vectors of all yield functions
		! --------------------------------------------------------------------------

		! compute all df/dsigma vectors
		Df_all = DfDs(sigma, alpha)

		! filter active df/dsigma vectors via vector with indices of active yield functions
		Df_actDsigma = Df_all(:,act_yf_ind)

	end function Df_actDsigma




	function Df_actDalpha(sigma, alpha, yield_flags, act_yf_ind)
		! computes vector with df_act/dalpha vectors of active yield functions

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags					! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(alpha), size(act_yf_ind)) :: Df_actDalpha	! array with df/dalpha vectors of active yield functions

		real(kind=dbl), dimension(size(alpha), n_yf) :: Df_all	! array with df/dalpha vectors of all yield functions
		! --------------------------------------------------------------------------

		! compute all df/dalpha vectors
		Df_all = DfDa(sigma, alpha)

		! filter active df/dalpha vectors via vector with indices of active yield functions
		Df_actDalpha = Df_all(:,act_yf_ind)

	end function Df_actDalpha




	function Dg_actDsigma(sigma, alpha, yield_flags, act_yf_ind)
		! computes vector with dg_act/dalpha vectors of active plastic potentials

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(6, size(act_yf_ind)) :: Dg_actDsigma	! array with dg/dsigma vectors of active yield functions

		! internal variables
		real(kind=dbl), dimension(6, n_yf) :: Dg_all	! array with dg/dsigma vectors of all plastic potentials
		! --------------------------------------------------------------------------

		! compute all df/dalpha vectors
		Dg_all = DgDs(sigma, alpha, yield_flags)

		! filter active df/dalpha vectors via vector with indices of active yield functions
		Dg_actDsigma = Dg_all(:,act_yf_ind)


	end function Dg_actDsigma




	function DDg_actDDsigma(sigma, alpha, yield_flags, act_yf_ind)
		! computes array with ddg_act/ddsigma matrices of active plastic potentials

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(6,6,size(act_yf_ind)) :: DDg_actDDsigma			! array, with ddg/ddsigma (hessians) arrays of active yield functions

		! internal variables
		real(kind=dbl), dimension(6,6,n_yf) :: DDg_all	! array with ddg/ddsigma (hessians) matrices of all plastic potentials
		! --------------------------------------------------------------------------

		! compute all ddg/ddsigma matrices
		DDg_all = DDgDDs(sigma, alpha)

		! filter active ddg/ddsigma matrices via vector with indices of active yield functions
		DDg_actDDsigma = DDg_all(:,:,act_yf_ind)

	end function DDg_actDDsigma




	function DDg_actDsigmaDalpha(sigma, alpha, yield_flags, act_yf_ind)
		! computes array with ddg_act/dsigmadalpha matrices of active plastic potentials

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(6,size(alpha), size(act_yf_ind)) :: DDg_actDsigmaDalpha		! array with ddg/dsigmadalpha matrices of active plastic potentials

		! internal variables
		real(kind=dbl), dimension(6,size(alpha),n_yf) :: DDg_all	! array with ddg/dsigmadalpha matrices of all plastic potentials
		! --------------------------------------------------------------------------

		! compute all ddg/dsigmadalpha matrices
		DDg_all = DDgDsDa(sigma, alpha)

		! filter active ddg/dsigmadalpha matrices via vector with indices of active yield functions
		DDg_actDsigmaDalpha = DDg_all(:,:,act_yf_ind)

	end function DDg_actDsigmaDalpha




	function hf_act(sigma, alpha, yield_flags, act_yf_ind)
		! computes array with hardening function vectors of active yield functions

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(alpha),size(act_yf_ind)) :: hf_act	! array with hardening vectors of active yield functions

		! internal variables
		real(kind=dbl), dimension(size(alpha),n_yf) :: hardening_f_all	! array with hardening vectors of all yield functions
		! --------------------------------------------------------------------------

		! compute all hardening function vectors
		hardening_f_all = hf(sigma, alpha, yield_flags)

		! filter active hardening function vectors via vector with indices of active yield functions
		hf_act = hardening_f_all(:,act_yf_ind)

	end function hf_act



	function Dh_actDsigma(sigma, alpha, yield_flags, act_yf_ind)
		! computes array with dh/dsigma matrices of active yield functions

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(alpha),6, size(act_yf_ind)) :: Dh_actDsigma		! array with dh/dsigma matrices of active yield functions

		! internal variables
		real(kind=dbl), dimension(size(alpha),6, n_yf) :: dh_all	! array with dh/dsigma matrices of all yield functions
		! --------------------------------------------------------------------------

		! compute all dh/dsigma matrices
		dh_all = DhDs(sigma, alpha)

		! filter active dh/dsigma matrices via vector with indices of active yield functions
		Dh_actDsigma = dh_all(:,:,act_yf_ind)

	end function Dh_actDsigma



	function Dh_actDalpha(sigma, alpha, yield_flags, act_yf_ind)
		! computes array with dh/dalpha matrices of active yield functions

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: sigma					! stress vector
		real(kind=dbl), dimension(:), intent(in) :: alpha					! internal variables vector
		logical, dimension(:), intent(in) :: yield_flags		! active yield functions vector
		integer, dimension(:), intent(in) ::			act_yf_ind		! vector with indices of active yield functions

		! return variable
		real(kind=dbl), dimension(size(alpha),size(alpha), size(act_yf_ind)) :: Dh_actDalpha		! array with dh/dalpha matrices of active yield functions

		! internal variables
		real(kind=dbl), dimension(size(alpha),size(alpha), n_yf) :: dh_all	! array with dh/dalpha matrices of all yield functions
		! --------------------------------------------------------------------------

		! compute all dh/dalpha matrices
		dh_all = DhDa(sigma, alpha)

		! filter active dh/dalpha matrices via vector with indices of active yield functions
		Dh_actDalpha = dh_all(:,:,act_yf_ind)

	end function Dh_actDalpha


end module model_comp_act
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
module return_mapping
    ! contains functions which perform the return mapping step
    ! of an elastic-plastic constitutive model

    ! load additional modules
    use constants
    use aux_routines, only: pack_x_vector, &
    						unpack_x_vector
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
module damage_var
	! contains routines for computation of damage part of the modified leon model

	! load additional modules
	use constants
	use derived_types

    implicit none

contains


	function comp_delta_alpha_d(delta_epsilon_p, alpha)
		! computes increment of internal scalar strain like softening variable,

		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: delta_epsilon_p		! incremental plastic strains vector
		real(kind=dbl), dimension(:), intent(in) :: alpha				! internal variables vector (for continuum hardening and softening)

		! return variable
		real(kind=dbl) :: comp_delta_alpha_d					! increment of internal scalar strain like softening variable

		! internal variable
		! --------------------------------------------------------------------------

		! compuation of delta alpha_d
		! .
		! .
		! .
		! comp_delta_alpha_d = ...


	end function comp_delta_alpha_d



end module damage_var
module damage_deriv
	! contains routines for computation of the derivatives of the damage variables
	! with respect to the strain vector

	! load additional modules
	use constants
	use derived_types


    implicit none

contains

	function comp_DomegaDalpha_d(alpha_d, l_char)
		! computes derivative of scalar damage variable omega
		! with respect to internal strain-like damage variable alpha_d
		! domega/dalpha_d
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), intent(in) :: alpha_d							! internal strain type damage variable
		real(kind=dbl), intent(in) :: l_char							! characteristic element length

		! return variable
		real(kind=dbl) :: comp_DomegaDalpha_d					! domega/dalpha_d

		! internal variable
		! --------------------------------------------------------------------------

		! computation of domega/dalpha_d
		! .
		! .
		! .
		! comp_DomegaDalpha_d = ...


	end function comp_DomegaDalpha_d


end module damage_deriv






















module damage
	! preforms computations related to the damage part of the material model
	! eg: computation of effective stress vector and
	! computation of damaged elastoplastic tangent

	! load additional modules
	use constants
	use material_info
	use damage_var
	use damage_deriv

    implicit none

contains

	subroutine comp_damage_components(sigma, eps_p, delta_eps_p, alpha, l_char, &
										el_inc_flag, &
										alpha_d, omega, &
										C_ep, C_el_inv, &
										eps, eps_pl_old)
		! computes the nominal stress vector and
		! damaged elastoplastic tangent
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: eps_p						! plastic strain vector
		real(kind=dbl), dimension(6), intent(in) :: delta_eps_p					! plastic strain increment vector
		real(kind=dbl), dimension(:), intent(in) :: alpha						! internal variables vector
		logical, intent(in) :: el_inc_flag										! elastic increment flag, states if eincrement is elastic (true) or plastic (false)
		real(kind=dbl), intent(in) :: l_char									! characteristic element length
		real(kind=dbl), dimension(6), intent(in) :: eps, eps_pl_old				! current strains, plastic strains of last step
		real(kind=dbl), dimension(6,6), intent(in) :: C_el_inv					! inverse elasticity matrix

		! input/output variables
		real(kind=dbl), dimension(6), intent(inout) :: sigma					! stress vector (effective on entry, nominal on exit)
		real(kind=dbl), dimension(6,6), intent(inout) :: C_ep					! elastoplastic tangent (undamaged on entry, damaged on exit)
		real(kind=dbl), dimension(n_int_vars_dam), intent(inout) :: alpha_d		! internal damage variables vector
		real(kind=dbl), dimension(n_omega), intent(inout) :: omega				! damage variables vector

		! internal variables
		! --------------------------------------------------------------------------

		! update of damage variable
		! .
		! .
		! .
		! omega = ...

		! computation of nominal stress with updated damage variable
		! .
		! .
		! .
		! sigma = ...

		! computation of damaged elastoplastic tangent
		C_ep = comp_Cep_damaged(delta_eps_p, alpha_d, omega, &
								alpha, C_ep, C_el_inv, &
								sigma, &
								l_char, el_inc_flag)

	end subroutine comp_damage_components



	function comp_Cep_damaged(delta_eps_p, alpha_d, omega, &
								alpha, C_ep, C_el_inv, &
								sigma, &
								l_char, el_inc_flag)
		! computes damaged elastic plastic tangent
		implicit none

		! --------------------------------------------------------------------------
		! passed variables
		real(kind=dbl), dimension(6), intent(in) :: delta_eps_p					! plastic strain increment vector
		real(kind=dbl), dimension(:), intent(in) :: alpha_d						! internal damage variables vector
		real(kind=dbl), dimension(:), intent(in) :: omega						! damage variables vector
		real(kind=dbl), dimension(:), intent(in) :: alpha						! internal variables vector
		real(kind=dbl), dimension(6,6), intent(in) :: C_ep, C_el_inv			! elastic plastic tangent, inverse elasticity matrix
		real(kind=dbl), dimension(6), intent(in) :: sigma						! stress vector (effective)
		real(kind=dbl), intent(in) :: l_char									! characteristic element length
		logical, intent(in) :: el_inc_flag										! elastic increment flag, states if eincrement is elastic (true) or plastic (false)

		! output variable
		real(kind=dbl), dimension(6,6) :: comp_Cep_damaged						! damaged elatic plastic tangent

		! internal variables
		! --------------------------------------------------------------------------

	end function comp_Cep_damaged



end module damage
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

!		write(*,*) ''
!		write(*,*) 'f_trial: ', f_trial

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
subroutine umat(stress, statev, ddsdde, sse, spd, scd, &
	                 rpl, ddsddt, drplde, drpldt, &
	                 stran , dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
	                 ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt, &
	                 celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

	use constants
	use stress_update, only: umat_tsaiwu

	implicit none

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
	real(kind=dbl), intent(in), dimension(ntens) :: dstran		! Array of strain increments
	real(kind=dbl), intent(in), dimension(2) :: time			! time(1): Value of step time at the beginning of the current increment or frequency
																! time(2): Value of total time at the beginning of the current increment
	real(kind=dbl), intent(in) :: dtime							! time increment
	real(kind=dbl), intent(in) :: temp							! Temperature at the start of the increment
	real(kind=dbl), intent(in) :: dtemp							! Increment of temperature
	real(kind=dbl), intent(in) :: predef						! Array of interpolated values of predefined field variables at this point at the start of the increment, based on the values read in at the nodes
	real(kind=dbl), intent(in) :: dpred							! Array of increments of predefined field variables
	character(len=80), intent(in) :: cmname						! User-defined material name, left justified
	real(kind=dbl), intent(in), dimension(nprops) :: props		! User-specified array of material constants associated with this user material
	real(kind=dbl), intent(in), dimension(3) :: coords			! An array containing the coordinates of this point
	real(kind=dbl), intent(in), dimension(3,3) :: drot			! Rotation increment matrix. This matrix represents the increment of rigid body rotation of the basis system in which the components of stress (STRESS) and strain (STRAN) are stored
	real(kind=dbl), intent(in) :: celent						! Characteristic element length, which is a typical length of a line across an element for a first-order element; it is half of the same typical length for a second-order element
	real(kind=dbl), intent(in), dimension(3,3) :: dfgrd0		! Array containing the deformation gradient at the beginning of the increment
	real(kind=dbl), intent(in), dimension(3,3) :: dfgrd1		! Array containing the deformation gradient at the end of the increment

	! --------------------------------------------------------------
	! output variables in all situations
	! --------------------------------------------------------------
	real(kind=dbl), intent(out), dimension(ntens,ntens) :: ddsdde		! Jacobian matrix of the constitutive model ddelta_sigma/ddelta_lambda
	real(kind=dbl), intent(inout), dimension(ntens) :: stress			! This array is passed in as the stress tensor at the beginning of the increment and must be updated in this routine to be the stress tensor at the end of the increment
	real(kind=dbl), intent(inout), dimension(nstatv) :: statev			! An array containing the solution-dependent state variables. These are passed in as the values at the beginning of the increment
	real(kind=dbl), intent(inout) :: sse								! Specific elastic strain energy, passed in as the values at the start of the increment and should be updated to the corresponding specific energy values at the end of the increment
	real(kind=dbl), intent(inout) :: spd								! plastic dissipation, passed in as the values at the start of the increment and should be updated to the corresponding specific energy values at the end of the increment
	real(kind=dbl), intent(inout) :: scd								! “creep” dissipation, passed in as the values at the start of the increment and should be updated to the corresponding specific energy values at the end of the increment

	! output variables Only in a fully coupled thermal-stress or a coupled thermal-electrical-structural analysis
	real(kind=dbl), intent(out) :: rpl							! Volumetric heat generation per unit time at the end of the increment caused by mechanical working of the material
	real(kind=dbl), intent(out), dimension(ntens) :: ddsddt		! Variation of the stress increments with respect to the temperature
	real(kind=dbl), intent(out), dimension(ntens) :: drplde		! Variation of RPL with respect to the strain increments
	real(kind=dbl), intent(out) :: drpldt						! Variation of RPL with respect to the temperature

	! output variables Only in a geostatic stress procedure or a coupled pore fluid diffusion/stress analysis for pore pressure cohesive elements
	!real(kind=dbl), intent(out) :: rpl							! RPL is used to indicate whether or not a cohesive element is open to the tangential flow of pore fluid

	! variables that can be updated
	real(kind=dbl), intent(out) :: pnewdt						! Ratio of suggested new time increment to the time increment being used. This variable allows you to provide input to the automatic time incrementation algorithms in Abaqus/Standard.
	! --------------------------------------------------------------

	call umat_tsaiwu(stress, statev, ddsdde, sse, spd, scd, &
	                 rpl, ddsddt, drplde, drpldt, &
	                 stran , dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
	                 ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt, &
	                 celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

end subroutine umat
