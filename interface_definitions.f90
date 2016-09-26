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
