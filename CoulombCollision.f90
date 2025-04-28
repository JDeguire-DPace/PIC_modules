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
                A_tau = 3.d0 * dexp(-min(t, 40.d0))
            end if

            analytical_approx = A_tau
        end function analytical_approx

end module solution_transcendental_nanbu

module rand_generator
    ! Module for random number generation using a linear congruential generator
    implicit none
    integer, parameter :: ia = 16807, im = 2147483647, iq = 127773, ir = 2836
    real(kind=8), parameter :: am = 1.d0 / real(im,8)
contains
    subroutine rran2(r, seed)
        ! Generates a random number between 0 and 1 using a seed
        implicit none
        real(kind=8), intent(out) :: r
        integer, intent(inout) :: seed
        integer :: k

        k = seed / iq
        seed = ia * (seed - k * iq) - ir * k
        if (seed < 0) seed = seed + im
        r = am * real(seed,8)
    end subroutine rran2
end module rand_generator





module Part
    ! Module defining the particle type
    implicit none
    real(kind=8), parameter :: pi = 4.d0*datan(1.d0)
    type :: particle
        real(kind=8) :: mass, charge, flow_velocity, energy, density
        real(kind=8), dimension(3) :: velocity
    end type
end module Part







program CoulombCollisions
    use Part
    use rand_generator
    use solution_transcendental_nanbu

    implicit none

    ! Physical constants and parameters
    real(kind=8), parameter :: e = dexp(1.d0)
    integer, parameter :: Number_electrons = 10**5, Number_ions = 10**5, Z_ion = 1 ! Reduce the number of particles to accelerate the simulation
    real(kind=8), parameter :: Delta_t = 1.0d-11, density_ions = 1.d21 ! Change Delta_t to check convergence
    real(kind=8), parameter :: mass_electrons = 9.10938356d-31, elementary_charge = 1.6021766209d-19
    real(kind=8), parameter :: density_electrons = real(Z_ion,8) * density_ions, k_B = 1.380649d-23, t_final = 2.d-5
    real(kind=8), parameter :: epsilon0 = 8.854187817d-12

    ! Derived parameters
    real(kind=8) :: mass_ion = 5.d0 * mass_electrons
    real(kind=8) :: charge_electrons = -1.d0*elementary_charge
    real(kind=8) :: charge_ions = real(Z_ion,8) * elementary_charge
    real(kind=8), parameter :: Initial_energy_electrons = 1000.d0 * elementary_charge
    real(kind=8), parameter :: Initial_energy_ions = 100.d0 * elementary_charge
    real(kind=8) :: Initial_flow_velocity_elec = dsqrt(Initial_energy_electrons / mass_electrons)
    real(kind=8) :: Initial_flow_velocity_ion = 0.d0
    real(kind=8), parameter :: density = density_electrons + density_ions
    integer, parameter :: Number_collision_el_ion = Number_electrons * density_ions / density
    integer, parameter :: Number_collision_el_el = (Number_electrons - Number_collision_el_ion) / 2.d0
    integer, parameter :: Number_collision_ion_ion = (Number_ions - Number_collision_el_ion) / 2.d0

    ! Variables
    type(particle), dimension(Number_electrons) :: array_electrons
    type(particle), dimension(Number_ions) :: array_ions
    integer, dimension(Number_electrons) :: index_electrons
    integer, dimension(Number_ions) :: index_ions
    integer :: i, j, i_particle, k_particle, i_tstep, number_time_step = t_final/Delta_t
    real(kind=8) :: time
    real(kind=8), dimension(3) :: tot_velec, tot_vion, tot_velec2, tot_vion2
    real(kind=8) :: averageV_ele, averageV2_ele, averageV2_ion, averageV_ion, Temp_elec, Temp_ion, Flow_ion, Flow_elec
    integer :: i_seed = 123456789
    real(kind=8) :: random_Gauss_ion(3), random_Gauss_elec(3)

    ! Open output file
    open(10, file="Coulomb_11.txt", status="replace") ! Adapt title to compare different dt's

    ! Initialize index arrays
    do i = 1, Number_electrons
        index_electrons(i) = i
    end do
    do i = 1, Number_ions
        index_ions(i) = i
    end do

    ! Initialize particles
    array_electrons%mass = mass_electrons
    array_electrons%charge = charge_electrons
    array_ions%mass = mass_ion
    array_ions%charge = charge_ions
    array_electrons%flow_velocity = Initial_flow_velocity_elec
    array_ions%flow_velocity = Initial_flow_velocity_ion
    array_electrons%energy = Initial_energy_electrons
    array_ions%energy = Initial_energy_ions

    Temp_elec = Initial_energy_electrons
    Temp_ion = Initial_energy_ions
    Flow_ion = Initial_flow_velocity_ion
    Flow_elec = Initial_flow_velocity_elec

    array_electrons%density = density_electrons
    array_ions%density = density_ions



    ! Assign initial velocities to particles following Maxwell-Boltzmann distribution (Nanbu & Yonemura, 1998)
    do i = 1, Number_electrons
        do j = 1, 3
            call rran2(random_Gauss_elec(j), i_seed)
        end do
        array_electrons(i)%velocity = maxwell_boltzmann(Initial_flow_velocity_elec, dsqrt(3.d0 * Initial_energy_electrons / mass_electrons), random_Gauss_elec)
    enddo
    do i = 1, Number_ions    
        do j = 1, 3
            call rran2(random_Gauss_ion(j), i_seed)
        end do
        array_ions(i)%velocity = maxwell_boltzmann(Initial_flow_velocity_ion, dsqrt(3.d0 * Initial_energy_ions / mass_ion), random_Gauss_ion)
    end do




    ! Time-stepping loop
    do i_tstep = 0, number_time_step
        time = i_tstep * Delta_t
        if (i_tstep == 0) then
            array_electrons%flow_velocity = Initial_flow_velocity_elec
            array_ions%flow_velocity = Initial_flow_velocity_ion
            array_electrons%energy = Initial_energy_electrons
            array_ions%energy = Initial_energy_ions
        else
            array_electrons%flow_velocity = Flow_elec
            array_ions%flow_velocity = Flow_ion
            array_electrons%energy = Temp_elec
            array_ions%energy = Temp_ion
        end if
        ! Reset averages
        averageV_ele = 0.d0
        averageV2_ele = 0.d0
        averageV2_ion = 0.d0
        averageV_ion = 0.d0

        ! Randomize particle indices for random partner selection
        call random_permutation(Number_electrons, index_electrons)
        call random_permutation(Number_ions, index_ions)

        !----------------------------------------------------------
        !                   Collision loops
        !
        ! 1- Ion-electron collisions
        ! 2- Ion-ion collisions
        ! 3- Electron-electron collisions
        !----------------------------------------------------------

        ! Ion-electron collisions
        do i_particle = 1, Number_collision_el_ion
            
            if (i_particle > Number_electrons .or. i_particle > Number_ions) exit ! Safety check

            ! Adding the velocity of the electron and ion to the averages
            averageV_ele = averageV_ele + (1.d0 / Number_electrons) * dsqrt(max(0.d0, sum(array_electrons(i_particle)%velocity**2)))
            averageV_ion = averageV_ion + (1.d0 / Number_ions) * dsqrt(max(0.0, sum(array_ions(i_particle)%velocity**2)))
            averageV2_ele = averageV2_ele + (1.d0 / Number_electrons) * sum(array_electrons(i_particle)%velocity**2)
            averageV2_ion = averageV2_ion + (1.d0 / Number_ions) * sum(array_ions(i_particle)%velocity**2)
            
            call collision(array_electrons(index_electrons(i_particle)), array_ions(index_ions(i_particle))) ! Collision subroutine
        end do

        ! Ion-ion collisions
        do i_particle = Number_collision_el_ion + 1, Number_collision_el_ion + Number_collision_ion_ion
            
            if (i_particle > Number_ions) exit ! Safety check on i_particle
            k_particle = i_particle + Number_collision_ion_ion
            if (k_particle > Number_ions) k_particle = k_particle - Number_ions ! Safety check on k_particle

            ! Adding the velocity of the 2 ions to the averages. We divide by 2 because we have 2 ions 
            averageV_ion = averageV_ion + (1.d0 / 2.d0 / Number_ions) * (dsqrt(sum(array_ions(i_particle)%velocity**2)) + dsqrt(sum(array_ions(k_particle)%velocity**2)))
            averageV2_ion = averageV2_ion + (1.d0 / 2.d0 / Number_ions) * (sum(array_ions(i_particle)%velocity**2) + sum(array_ions(k_particle)%velocity**2))
            
            call collision(array_ions(index_ions(i_particle)), array_ions(index_ions(k_particle))) ! Collision subroutine
        end do

        ! Electron-electron collisions
        do i_particle = Number_collision_el_ion + 1, Number_collision_el_ion + Number_collision_el_el

            if (i_particle > Number_electrons) exit ! Safety check on i_particle
            k_particle = i_particle + Number_collision_el_el
            if (k_particle > Number_electrons) k_particle = k_particle - Number_electrons ! Safety check on k_particle

            ! Adding the velocity of the 2 electrons to the averages. We divide by 2 because we have 2 electrons
            averageV_ele = averageV_ele + (1.d0 / 2.d0 / Number_electrons) * (dsqrt(sum(array_electrons(i_particle)%velocity**2)) + dsqrt(sum(array_electrons(k_particle)%velocity**2)))
            averageV2_ele = averageV2_ele + (1.d0 / 2.d0 / Number_electrons) * (sum(array_electrons(i_particle)%velocity**2) + sum(array_electrons(k_particle)%velocity**2))
            
            call collision(array_electrons(index_electrons(i_particle)), array_electrons(index_electrons(k_particle))) ! Collision subroutine
        end do

        ! Compute total velocity vectors
        tot_velec = 0.d0
        tot_vion = 0.d0
        tot_velec2 = 0.d0
        tot_vion2 = 0.d0


        do i = 1, Number_electrons
            tot_velec = tot_velec + array_electrons(i)%velocity
            tot_velec2 = tot_velec2 + array_electrons(i)%velocity**2
        end do

        tot_velec = tot_velec / Number_electrons
        tot_velec2 = tot_velec2 / Number_electrons

    

        do i = 1, Number_ions
            tot_vion = tot_vion + array_ions(i)%velocity
            tot_vion2 = tot_vion2 + array_ions(i)%velocity**2
        end do
        tot_vion = tot_vion / Number_ions
        tot_vion2 = tot_vion2 / Number_ions

        ! Compute temperatures and write to file
        Temp_elec = mass_electrons * (sum(tot_velec2) - sum(tot_velec**2)) / (3.d0 * elementary_charge)
        Temp_ion = mass_ion * (sum(tot_vion2) - sum(tot_vion**2)) / (3.d0 * elementary_charge)
        Flow_elec = dsqrt(max(0.0, sum(tot_velec**2))) / Initial_flow_velocity_elec
        Flow_ion = dsqrt(max(0.0, sum(tot_vion**2))) / Initial_flow_velocity_elec

        write(10, '(5EN14.4)') time, Temp_elec, Temp_ion, Flow_elec, Flow_ion ! Write to file
    end do

    close(10)



contains



    ! Function adapted from LEPIC -> load_gauss() 
    function maxwell_boltzmann(flow_speed, v_th, rnd) result(velocity)
        ! Generates a 3D velocity vector following a Maxwell-Boltzmann distribution
        implicit none
        real(kind=8), intent(in) :: flow_speed, v_th
        real(kind=8), intent(in) :: rnd(3)
        real(kind=8), dimension(3) :: velocity
        real(kind=8) :: theta, phi, vp

        vp = v_th * dsqrt(-dlog(max(1.0d-12, rnd(1))))
        theta = 2.0d0 * pi * rnd(2)
        phi = acos(2.0d0 * rnd(3) - 1.0d0)

        velocity(1) = flow_speed + vp * sin(phi) * cos(theta)
        velocity(2) = vp * sin(phi) * sin(theta)
        velocity(3) = vp * cos(phi)
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
        real(kind=8) :: CoulombLog, debye_length, A, tau, b_min

        ! Nanbu's method
        g = particle1%velocity - particle2%velocity
        debye_length = dsqrt(epsilon0*(particle1%energy+particle2%energy) / (density_ions*charge_ions**2 + density_electrons*charge_electrons**2)) ! Debye length
        b_min = abs(particle1%charge * particle2%charge) / (4.d0 * pi * epsilon0 * (particle1%energy+particle2%energy)) ! Minimum impact parameter
        CoulombLog = dlog(debye_length/b_min)
        A = (density / (4.d0 * pi)) * (particle1%charge * particle2%charge / (epsilon0 * reduced_mass(particle1%mass, particle2%mass)))**2 * CoulombLog
        tau = A * Delta_t / max(1.d-18, dsqrt(sum(g**2))**3)
        g_perp = dsqrt(max(0.d0, g(2)**2 + g(3)**2))
        call rran2(angle_azimuthal, i_seed)
        call rran2(rand1, i_seed)
        
        
        ! Either use the numerical Newton Raphson solution (NewtRaph(tau)) or the approximation analytical solution (analytical_approx(tau)) 
        A_tau = NewtRaph(tau)

        angle_azimuthal = 2.d0 * pi * angle_azimuthal
        h(1) = g_perp * dcos(angle_azimuthal)
        h(2) = -(g(1) * g(2) * dcos(angle_azimuthal) + dsqrt(max(0.d0, sum(g**2))) * g(3) * dsin(angle_azimuthal)) / max(1.d-18, g_perp)
        h(3) = -(g(1) * g(3) * dcos(angle_azimuthal) - dsqrt(max(0.d0, sum(g**2))) * g(2) * dsin(angle_azimuthal)) / max(1.d-18, g_perp)

        argument_acos = (A_tau+dlog(rand1+(1.d0-rand1)*dexp(-2.d0*A_tau))) / A_tau
        angle_diffusion = dacos(argument_acos)
        
        particle1%velocity = particle1%velocity - particle2%mass / (particle1%mass + particle2%mass) * (g * (1.d0 - dcos(angle_diffusion)) + h * dsin(angle_diffusion))
        particle2%velocity = particle2%velocity + particle1%mass / (particle1%mass + particle2%mass) * (g * (1.d0 - dcos(angle_diffusion)) + h * dsin(angle_diffusion))
    end subroutine collision



    ! Ensure the randomness of the collision partners
    subroutine random_permutation(n, indices)
        ! Generates a random permutation of indices
        implicit none
        integer, intent(in) :: n
        integer, intent(inout) :: indices(n)
        integer :: i, j, temp
        real(kind=8) :: r

        do i = n, 2, -1
            call rran2(r, i_seed)
            j = 1 + int(r * i)
            temp = indices(i)
            indices(i) = indices(j)
            indices(j) = temp
        end do
    end subroutine random_permutation




end program CoulombCollisions