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



	SUBROUTINE write_real_list (matrix, m_zeilen, n_spalten)
		IMPLICIT NONE
		! subroutine that writes a real 2D array to the standard io

		! Definition der calling-Parameter
		INTEGER, INTENT(IN) :: m_zeilen								! Anzahl der Zeilen (eingehender Wert)
		INTEGER, INTENT(IN) :: n_spalten							! Anzahl der Spalten (eingehender Wert)
		REAL(kind=dbl), DIMENSION(m_zeilen,n_spalten), INTENT(IN) :: matrix	! Matrix mit m Zeilen und n Spalten (eingehender Wert)

		! interne Parameter
		INTEGER :: i		! Zeilenzähler
		INTEGER :: j		! Spaltenzähler
		CHARACTER(len=50) :: form		! Formatbeschreiber

		! Schreiben des Inhaltes des Formatbeschreibers
		WRITE (form,10) n_spalten
		10 FORMAT ( '(' I3, '(1X, E14.7))')			! I3 wird im Programmablauf durch die Zahl ersetzt, die  '(I3(1X, F10.3)'
																		! der Spaltenanzahl entspricht

		! Schreiben der Matrix im richtigen Format

		WRITE (*,form) ((matrix(i,j), j=1,n_spalten), i=1,m_zeilen)

	END SUBROUTINE write_real_list


end module aux_routines
