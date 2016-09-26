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
