module math_routines_lib
	! module for various mathematical functions and subroutines

	! load additional modules
	!use constants_lib
	use iso_c_binding, only: dbl => c_double

    implicit none

contains

	subroutine calc_rot_matrix_zyz(T_rot, alpha, beta, gam) bind(c)
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
		real(kind=dbl), dimension(3,3), intent(out) :: T_rot		! rotation matrix with direction cosines between old an rotated coordinate basis

		! internal variables
		real(kind=dbl), dimension(3,3) :: Tz1	! rotation matrix for first rotation around z axis
		real(kind=dbl), dimension(3,3) :: Ty2	! rotation matrix for second rotation around rotated y axis
		real(kind=dbl), dimension(3,3) :: Tz3	! rotation matrix for third rotation around rotated z axis
		real(kind=dbl), dimension(3,3) :: Tzyz	! rotation matrix all three above combined into one: Tz1*Ty2*Tz3
		! --------------------------------------------------------------------------

		Tz1 = reshape( (/ cos(-gam), -sin(-gam), 0._dbl, sin(-gam), cos(-gam), 0._dbl, 0._dbl, 0._dbl, 1._dbl /), (/3,3/) )
		Ty2 = reshape( (/ cos(beta), 0._dbl, sin(beta), 0._dbl, 1._dbl, 0._dbl, -sin(beta), 0._dbl, cos(beta) /), (/3,3/) )
		Tz3 = reshape( (/ cos(-alpha), -sin(-alpha), 0._dbl, sin(-alpha), cos(-alpha), 0._dbl, 0._dbl, 0._dbl, 1._dbl /), (/3,3/) )
		T_rot = transpose(matmul(Tz1, matmul(Ty2, Tz3) ))


	end subroutine calc_rot_matrix_zyz



	function mcauly(x) bind(c)
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


	elemental function mcauly_elemental(x)
		! computes <x> (mcauly brackets
		! if x < 0   =>  <x> = 0
		! if x >= 0  =>  <x> = x
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl) :: mcauly_elemental

		! passed variable
		real(kind=dbl), intent(in) :: x
		! --------------------------------------------------------------------------

		if (x < 0._dbl) then

			mcauly_elemental = 0._dbl

		else

			mcauly_elemental = x

		end if

	end function mcauly_elemental


	function heaviside(x) bind(c)
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


	elemental function heaviside_elemental(x)
		! computes heaviside function H(x)
		! if x < 0  => H(x) = 0
		! if x >= 0  => H(x) = 1
		implicit none

		! --------------------------------------------------------------------------
		! return variable
		real(kind=dbl) :: heaviside_elemental

		! passed variable
		real(kind=dbl), intent(in) :: x
		! --------------------------------------------------------------------------

		if (x < 0._dbl) then

			heaviside_elemental = 0._dbl

		else

			heaviside_elemental = 1._dbl

		end if

	end function heaviside_elemental

end module math_routines_lib
