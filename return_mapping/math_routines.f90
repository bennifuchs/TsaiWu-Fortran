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
