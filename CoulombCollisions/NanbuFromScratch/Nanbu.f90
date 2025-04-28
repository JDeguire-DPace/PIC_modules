program Nanbu
    implicit none
    integer :: i, j, k, n
    integer, parameter :: number_particles = 10**5
    integer, parameter :: number_microCollision = 3000
    real(8) :: random_polar, random_azimuthal
    real(8) :: theta_average = 1.017e-3
    real(8) :: s,A, mean_v_x, mean_v_y, mean_v_z
    integer, parameter :: n_samples = 100000            ! Number of samples to generate
    real(8), parameter :: pi = 3.141592653589793d0     ! Pi constant
    
    real(8) :: x, y, fx, fmax, u1, u2                  ! Temporary variables for sampling
    real(8), dimension(n_samples) :: samples           ! Array to store accepted samples
    integer :: count, index           
    real(8), dimension(number_particles) :: v_x,v_y,v_z,chi,phi
    open(unit=20, file="vel.txt", status="replace")
    open(unit=10, file="samples.txt", status="old")
    do i = 1, n_samples
        read(10,*) samples(i)
    end do
    close(10)
    v_x = 0.0
    v_y = 0.0
    v_z = 1.0
    chi = 0.0
    phi = 0.0
    do i=0,10000
        mean_v_z = 0.0
        mean_v_x = 0.0
        mean_v_y = 0.0

        do j=1,number_particles
            call random_number(random_azimuthal)
            call random_number(random_polar)

            index = int(random_polar * dble(n_samples)) + 1
            if (index > n_samples) index = n_samples

            random_polar = samples(index)
            random_azimuthal = 2.0 * pi * random_azimuthal
            chi(j) = chi(j)+random_polar
            phi(j) = phi(j)+random_azimuthal
            v_x(j) = cos(phi(j)) * sin(chi(j))
            v_y(j) = sin(phi(j)) * sin(chi(j))
            v_z(j) = cos(chi(j))
            mean_v_x = mean_v_x + v_x(j)**2/number_particles
            mean_v_y = mean_v_y + (v_y(j)**2+v_x(j)**2)/number_particles
            mean_v_z = mean_v_z + v_z(j)/number_particles
            
        end do
        write(20,'(I5,3F15.6)') i, mean_v_x, mean_v_y, abs(mean_v_z)
    end do
    close(20)
end program Nanbu