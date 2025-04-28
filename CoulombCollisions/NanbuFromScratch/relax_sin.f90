program relax_sin
    implicit none
    integer :: Number_particles = 200000, i_particle, i_it
    real :: polar_angle, azimuthal_angle
    real :: chi, psi, ran1, ran2, pi = atan(-1.0)
    real :: sum = 0.0
    real, dimension(3,3) :: I
    I = 0.0d0  ! Initialize the matrix with zeros
    do i_it = 1, 3
        I(i_it, i_it) = 1.0d0  ! Set diagonal elements to 1
    end do

    do i_particle = 1, Number_particles
        call random_number(ran1)
        call random_number(ran2)

        azimuthal_angle = 2*pi*ran1
        polar_angle = atan((pi/180.)/(2*sqrt(ran2)))

    enddo


    contains
        function A(theta, phi)
            implicit none
            real(8), dimension(3,3) :: A
            real(8), dimension(3,3) :: arrayTMP
            real(8) :: theta, phi

            arrayTMP(1,1) = cos(theta)*cos(phi)**2+sin(phi)**2
            arrayTMP(1,2) = (cos(theta)-1.)*sin(phi)*cos(phi)
            arrayTMP(1,3) = -sin(theta)*cos(phi)

            arrayTMP(2,1) = (cos(theta)-1.)*sin(phi)*cos(phi)
            arrayTMP(2,2) = cos(theta)*sin(phi)**2+cos(phi)**2
            arrayTMP(2,3) = -sin(theta)*cos(phi)

            arrayTMP(3,1) = sin(theta)*cos(phi)
            arrayTMP(3,2) = sin(theta)*cos(phi)
            arrayTMP(3,3) = cos(theta)
            A=arrayTMP

        end function A

        function C(theta, phi)
            implicit none
            real(8), dimension(3,3) :: C
            real(8), dimension(3,3) :: arrayTMP
            real(8) :: theta, phi
            arrayTMP = 0.
            
            arrayTMP(1,3) = -cos(phi)
            arrayTMP(2,3) = -sin(phi)
            arrayTMP(3,1) = cos(phi)
            arrayTMP(3,2) = sin(phi)

            C=arrayTMP

        end function C

        function D(theta, phi)
            implicit none
            real(8), dimension(3,3) :: D
            real(8), dimension(3,3) :: arrayTMP
            real(8) :: theta, phi
            arrayTMP = 0.
            
            arrayTMP(1,1) = cos(phi)**2
            arrayTMP(1,2) = cos(phi)*sin(phi)

            arrayTMP(2,1) = cos(phi)*sin(phi)
            arrayTMP(2,2) = sin(phi)**2

            arrayTMP(3,3) = 1

            D=arrayTMP

        end function D

end program relax_sin