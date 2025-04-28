! We use the Newton Raphson method to solve the equation coth(A) - 1/A = exp(-s_ab) as it was
! derived in [Nanbu, Yonemura, 1998].  

module root_finder
    implicit none
    contains
        function f(a1,a2) result(res)
            real(8), intent(in) :: a1,a2
            real(8)             :: res
            res = 1./tanh(a1) - 1./a1 - exp(-a2)
        end function f

        function df(a1) result(outp)
            real(8), intent(in) :: a1
            real(8)             :: outp
            outp = -1./sinh(a1)**2 + 1./a1**2
        end function df

        real function Newton_Raphson(t)
            real(8), intent(in) :: t
            real(8) :: A0 = 0.01, tol = 1e-6
            integer :: max_iteration = 100
            integer :: i
            real(8) :: fval, dfval
            do i=1,max_iteration
                fval = f(A0,t)
                dfval = df(A0)

                if (abs(fval) < tol) then
                    Newton_Raphson = A0
                    return
                end if
                A0 = A0 - fval/dfval
            end do
            Newton_Raphson = A0
        end function Newton_Raphson
end module root_finder



!program Newton_test
!    use root_finder
!    implicit none
!    real :: t = 0.001
!    real :: A 
!    A = Newton_Raphson(t)
!    print *, 'A = ', A
!end program