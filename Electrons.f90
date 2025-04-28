module solution_transcendental_nanbu
    implicit none
    contains
        real(kind=8) function f1(A, t)
            implicit none
            real(kind=8), intent(in) :: A, t
            ! Function to solve: coth(A) - 1/A - e^t
            f1 = 1.0 / dtanh(A) - 1.d0 / A - dexp(-t)
        end function f1

        real(kind=8) function f2(A)
            implicit none
            real(kind=8), intent(in) :: A
            ! Derivative of f(A): -csch^2(A) + 1/A^2
            f2 = -1.0 / dsinh(A)**2 + 1.d0 / A**2
        end function f2

        real(kind=8) function NewtRaph(t)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8) :: A0 = 10.d0, tol = 1.e-8
            integer :: max_iter = 50
            real(kind=8) :: A, f_val, df_val
            integer :: i

            A = A0
            do i = 1, max_iter
                f_val = f1(A,t)
                df_val = f2(A)

                if (abs(f_val) < tol) then
                    NewtRaph = A  ! Converged
                    return
                end if

                A = A - f_val / df_val  ! Newton step

                if (A <= 0.d0) then  ! Ensure A remains positive
                    A = A0 / 2.d0  ! Reset to a safer value
                end if
            end do

            NewtRaph = A  ! Return the last value if not converged
        end function NewtRaph


        real(8) function analytical_approx(t)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8) :: A_tau

            ! Analytical approximation for small tau
            if (t <= 0.01d0) then
                A_tau = 1.d0 / max(t, 1.d-10)
            else if (t <= 0.1d0) then
                A_tau = (84.708d0 * t + 166.4288d0) / (166.4681d0 * t - 0.0003d0)
            else if (t <= 0.6d0) then
                A_tau = 1.2528d0 * t**(-0.9208d0) + 0.1115d0 * t
            else if (t < 2.d0) then
                A_tau = (-50.8179d0 * t + 181.7747d0) / (88.0921d0 * t + 20.5056d0)
            else
                A_tau = 3.d0 * dexp(-min(t, 2000.d0))
            end if

            analytical_approx = A_tau
        end function analytical_approx

end module solution_transcendental_nanbu





module Part
    ! Module defining the particle type
    implicit none
    real(kind=8), parameter :: pi = 4.d0*datan(1.d0)
    type :: particle
        real(kind=8) :: mass, charge, energy, density
        real(kind=8), dimension(3) :: velocity
    end type
end module Part







program CoulombCollisions
    use Part
    use solution_transcendental_nanbu

    implicit none

    ! Physical constants and parameters
    real(kind=8), parameter :: density_electrons = 1.d18, k_B = 1.380649d-23
    real(kind=8), parameter :: e = dexp(1.d0), epsilon0 = 8.854187817d-12
    real(kind=8), parameter :: mass_electrons = 9.10938356d-31, elementary_charge = 1.6021766209d-19
    integer, parameter :: Number_electrons = 10000000 ! Reduce the number of particles to accelerate the simulation
    

    ! Derived parameters
    real(kind=8), parameter :: charge_electrons = -1.d0*elementary_charge
    real(kind=8), parameter :: EnergyTotaleElectrons = 1.1d0 * elementary_charge
    real(kind=8), parameter :: Initial_Tp = 3.d0*EnergyTotaleElectrons / 3.3d0  ! Initial temperature in the z-direction
    real(kind=8), parameter :: Initial_Tz = 1.3d0*Initial_Tp 
    real(kind=8), parameter :: Initial_energy_electrons = (1.d0/3.d0) * Initial_Tz + (2.d0/3.d0) * Initial_Tp
    real(kind=8), parameter :: density = density_electrons
    integer, parameter :: Number_collision_el_el = (Number_electrons) / 2.d0
    real(kind=8) :: E_parallel, E_perp
    real(kind=8), parameter :: LogCoulomb = dlog(4.d0*pi*(epsilon0*EnergyTotaleElectrons)**1.5d0/(dsqrt(density_electrons)*abs(charge_electrons)**3))
    real(kind=8), parameter :: tau0 =  3.d0*pi*dsqrt(2.d0*pi)*epsilon0**2*mass_electrons**2*(1.1*1.602d-19)**1.5/(charge_electrons**4*density_electrons*LogCoulomb)
    real(kind=8), parameter :: Delta_t = 0.02d0*tau0
    
    real(kind=8), parameter :: t_final = 12.d0*tau0
    
    ! Variables
    type(particle), dimension(Number_electrons) :: array_electrons
    integer, dimension(Number_electrons) :: index_electrons
    integer :: i, j, i_particle, k_particle, i_tstep, number_time_step = t_final/Delta_t
    real(kind=8) :: time, Temp_elec, vz, vp, Tp, Tz
    integer :: i_seed = 123456789
    real(kind=8) :: random_Gauss_elec(6)
    real(kind=8), dimension(3) :: tot_velec, tot_velec2 
    print *, Initial_Tp, Initial_Tz, Initial_energy_electrons, LogCoulomb
    ! Open output file
    open(10, file="Pre55_B.txt", status="replace") ! Adapt title to compare different dt's
    ! Initialize index arrays
    do i = 1, Number_electrons
        index_electrons(i) = i
    end do

    ! Initialize particles
    array_electrons%mass = mass_electrons
    array_electrons%charge = charge_electrons
    array_electrons%energy = Initial_energy_electrons
    array_electrons%density = density_electrons
    Temp_elec = Initial_energy_electrons

    ! Assign initial velocities to particles following Maxwell-Boltzmann distribution (Nanbu & Yonemura, 1998)
    do i = 1, Number_electrons
        do j = 1, 6
            call rran2(random_Gauss_elec(j), i_seed)
        end do
        array_electrons(i)%velocity = maxwell_boltzmann(dsqrt(3.d0 * Initial_Tz / mass_electrons), dsqrt(3.d0 * Initial_Tp / mass_electrons), random_Gauss_elec)
    enddo

    ! Time-stepping loop
    do i_tstep = 0, number_time_step
        time = i_tstep * Delta_t
        call random_permutation(Number_electrons, index_electrons, i_seed)
        if (i_tstep == 0) then
            array_electrons%energy = Initial_energy_electrons
        else
            array_electrons%energy = Temp_elec
        end if
        ! Randomize particle indices for random partner selection
        

        ! Electron-electron collisions
        do i_particle = 1, Number_collision_el_el
            if (i_particle > Number_electrons) exit ! Safety check on i_particle
            k_particle = i_particle + Number_collision_el_el
            if (k_particle > Number_electrons) k_particle = k_particle - Number_electrons ! Safety check on k_particle
            call collision(array_electrons(index_electrons(i_particle)), array_electrons(index_electrons(k_particle))) ! Collision subroutine
        end do

        ! Compute total velocity vectors and kinetic energies
        vz = 0.d0
        vp = 0.d0
        tot_velec = 0.d0
        tot_velec2 = 0.d0

        do i = 1, Number_electrons
            tot_velec = tot_velec + array_electrons(i)%velocity / Number_electrons
            tot_velec2 = tot_velec2 + array_electrons(i)%velocity**2 / Number_electrons
        end do

        ! Compute parallel and perpendicular kinetic energies
        E_parallel = 0.5d0 * mass_electrons * tot_velec2(1)  ! Parallel kinetic energy (x-direction)
        E_perp = 0.5d0 * mass_electrons * (tot_velec2(2) + tot_velec2(3))  ! Perpendicular kinetic energy (y and z directions)

        ! Compute temperatures
        Tp = (mass_electrons / (3.d0 * elementary_charge)) * (0.5d0*(tot_velec2(1)+tot_velec2(2)) - 0.5d0*(tot_velec(1)**2+tot_velec(2)**2))
        Tz = (mass_electrons / (3.d0 * elementary_charge)) * (tot_velec2(3) - tot_velec(3)**2)
        Temp_elec = (1.d0 / 3.d0) * Tz + (2.d0 / 3.d0) * Tp

        ! Write time, temperatures, and kinetic energies to the output file
        write(10, '(6EN14.4)') time/tau0, Temp_elec, Tp/Temp_elec, Tz/Temp_elec, E_parallel, E_perp
    end do

    close(10)



contains

    subroutine rran2(r, seed)
        ! Generates a random number between 0 and 1 using a seed
        implicit none
        real(kind=8), intent(inout) :: r  ! Output random number
        integer, intent(inout) :: seed   ! Random seed
        integer :: k
        integer, parameter :: ia = 16807, im = 2147483647, iq = 127773, ir = 2836
        real(kind=8), parameter :: am = 1.d0 / real(im, 8)

        ! Compute the next random number using the LCG formula
        k = seed / iq
        seed = ia * (seed - k * iq) - ir * k
        if (seed < 0) seed = seed + im

        ! Normalize the random number to the range [0, 1]
        r = am * real(seed, 8)
    end subroutine rran2





    ! Function adapted from LEPIC -> load_gauss() 
    function maxwell_boltzmann(v_parallel, v_perp, rnd) result(velocity)
        ! Generates a 3D velocity vector based on parallel and perpendicular energy distributions
        implicit none
        real(kind=8), intent(in) :: v_parallel  ! Parallel velocity magnitude
        real(kind=8), intent(in) :: v_perp      ! Perpendicular velocity magnitude
        real(kind=8), intent(in) :: rnd(6)      ! Array of random numbers
        real(kind=8), dimension(3) :: velocity  ! Resulting 3D velocity vector

        ! Random angle for the perpendicular velocity distribution

        ! Assign velocity components
        velocity(1) = v_perp * dsqrt(-2.d0 * dlog(rnd(1))) * dcos(2.d0 * pi * rnd(2))  ! Parallel component in the x-direction
        velocity(2) = v_perp * dsqrt(-2.d0 * dlog(rnd(3))) * dcos(2.d0 * pi * rnd(4))  ! Perpendicular component in the y-direction
        velocity(3) = v_parallel * dsqrt(-2.d0 * dlog(rnd(5))) * dcos(2.d0 * pi * rnd(6))  ! Perpendicular component in the z-direction
    end function maxwell_boltzmann



    real(kind=8) function reduced_mass(m1, m2)
        ! Computes the reduced mass of two particles
        implicit none
        real(kind=8), intent(in) :: m1, m2
        reduced_mass = m1 * m2 / (m1 + m2)
    end function reduced_mass



    subroutine collision(particle1, particle2)
        ! Handles collisions between two particles
        implicit none
        type(particle), intent(inout) :: particle1, particle2
        real(kind=8), dimension(3) :: g, h
        real(kind=8) :: g_perp, angle_azimuthal, rand1, A_tau, argument_acos, angle_diffusion
        real(kind=8) :: LAMBDA, debye_length, average_g2, A, tau

        ! Nanbu's method
        g = particle1%velocity - particle2%velocity
        debye_length = dsqrt(epsilon0*Temp_elec*elementary_charge / (density_electrons*charge_electrons**2)) ! Debye length
        average_g2 = 3.d0*elementary_charge * (particle1%energy/particle1%mass + particle2%energy/particle2%mass)
        !LAMBDA = 4.d0*pi*epsilon0*reduced_mass(particle1%mass, particle2%mass)*average_g2*debye_length/abs(particle1%charge*particle2%charge)
        LAMBDA = 4.d0*pi*(epsilon0*particle1%energy)**1.5d0/(dsqrt(density)*abs(charge_electrons)**3)

        A = density / (4.d0 * pi) * (particle1%charge * particle2%charge / (epsilon0 * reduced_mass(particle1%mass, particle2%mass)))**2 * dlog(LAMBDA)
        tau = A * Delta_t / max(1.d-28, dsqrt(sum(g**2))**3)
        g_perp = dsqrt(max(0.d0, g(2)**2 + g(3)**2))
        call rran2(angle_azimuthal, i_seed)
        call rran2(rand1, i_seed)
        
        
        ! Either use the numerical Newton Raphson solution (NewtRaph(tau)) or the approximation analytical solution (analytical_approx(tau)) 
        A_tau = analytical_approx(tau)
        angle_azimuthal = 2.d0 * pi * angle_azimuthal
        h(1) = g_perp * dcos(angle_azimuthal)
        h(2) = -(g(1) * g(2) * dcos(angle_azimuthal) + dsqrt(max(0.d0, sum(g**2))) * g(3) * dsin(angle_azimuthal)) / max(1.d-28, g_perp)
        h(3) = -(g(1) * g(3) * dcos(angle_azimuthal) - dsqrt(max(0.d0, sum(g**2))) * g(2) * dsin(angle_azimuthal)) / max(1.d-28, g_perp)
        argument_acos = (A_tau+dlog(rand1+(1.d0-rand1)*dexp(-2.d0*A_tau))) / A_tau
        !if (argument_acos<-1.d0 .or. argument_acos>1.d0) print *, "argument_acos = ", argument_acos

        argument_acos = max(-1.d0, min(1.d0, argument_acos))
        angle_diffusion = dacos(argument_acos)
        !angle_diffusion = max(0.d0, min(pi, argument_acos))
        particle1%velocity = particle1%velocity - reduced_mass(particle1%mass, particle2%mass)/particle1%mass * (g * (1.d0 - dcos(angle_diffusion)) + h * dsin(angle_diffusion))
        particle2%velocity = particle2%velocity + reduced_mass(particle1%mass, particle2%mass)/particle2%mass * (g * (1.d0 - dcos(angle_diffusion)) + h * dsin(angle_diffusion))
    end subroutine collision



    ! Ensure the randomness of the collision partners
    subroutine random_permutation(n, indices, seed)
        ! Generates a random permutation of indices
        implicit none
        integer, intent(in) :: n
        integer, intent(inout) :: indices(n)  ! Array to be permuted
        integer, intent(inout) :: seed        ! Random seed
        integer :: i, j, temp
        real(kind=8) :: r
        do i = n, 2, -1
            call rran2(r, seed)  ! Generate a random number between 0 and 1
            j = 1 + int(r * i)   ! Generate a random index in the range [1, i]
    
            ! Swap indices(i) and indices(j)
            temp = indices(i)
            indices(i) = indices(j)
            indices(j) = temp
        end do
    end subroutine random_permutation



end program CoulombCollisions