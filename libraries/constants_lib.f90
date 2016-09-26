module constants_lib
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

	! flag for activation/deactivation of damage computation
	logical, parameter :: damage_flag = .True.			! damage ist taken into account: .True.
														! damage is neglected: .False.

	! flag for switch between numerical/analytical computation of damaged elastoplastic tangent
	logical, parameter :: Cep_dam_num_flag = .False.	! Cep_dam is computed analytically: .False.
														! C_ep_dam is computed numerically: .True.

	! flag for switch between new/old plastic potential
	integer, parameter :: plast_pot = 3					! 1 ... plastic potential according to Unteregger (old version 2012)
														! 2 ... plastic potential according to Grassl, Jirasek, Valentini
														! 3 ... plastic potential according to Unteregger (new version 2014)

	! flag for switch between isotropic/transversal isotropic elasticity
	logical, parameter :: anisotropic_elastic_flag = .True.				! isotropic elasticity: .False.
																		! transversal isotropic elasticity: .True.
end module constants_lib

