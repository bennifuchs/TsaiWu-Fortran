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
