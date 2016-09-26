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
