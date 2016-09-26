module material_info
    implicit none

    integer, save :: n_yf				! number of yield functions
	integer, save :: n_int_vars			! number of internal hardening variables
	integer, save :: n_int_vars_dam		! number of internal softening variables
	integer, save :: n_omega			! number of damage variables

	logical, save :: hardening_present			! True if hardening is present, False if no hardening law is detected

	logical, save :: damage_present		! True if damage law is detected, false if no damage law is detected

end module material_info
