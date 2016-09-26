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
