module rand_generator
    ! Module for random number generation
    implicit none
    integer, parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
    real(kind=8), parameter :: am = 1.0_16 / im
  
  contains
    ! Subroutine rran2 generates a random number between 0 and 1
    ! It modifies the seed to ensure different random numbers are generated each time
    subroutine rran2(r,seed)
        implicit none
        real(kind=8), intent(out) :: r
        integer, intent(inout) :: seed
        integer :: k
    
        k = seed / iq
        seed = ia * (seed - k * iq) - ir * k
        if (seed < 0) seed = seed + im
        r = am * seed
    end subroutine rran2
  
  end module rand_generator




module Part
    ! Creation of the type particle
    ! A particle has a mass, a charge and a velocity
    implicit none
    real(kind=8), parameter :: pi = dacos(-1.d0)
    type :: particle
        real(kind=8) :: mass, charge
        real(kind=8), dimension(3) :: velocity
    end type
end module Part


program CoulombCollisions
    use Part
    use rand_generator

    implicit none
    
    ! Parameters
    real(kind=8), parameter :: e = dexp(1.d0)
    integer, parameter :: Number_alpha = 10**5, Number_beta = 10**5, Z_ion = 1
    real(kind=8), parameter :: Delta_t = 1.0d-7, density_ions = 1.d21, CoulombLog = 15.9d0
    real(kind=8), parameter :: mass_electrons = 9.109d-31, charge_electrons = 1.602d-19
    real(kind=8), parameter :: density_electrons = Z_ion * density_ions, k_B = 1.380649d-23, t_final = 2.d-5
    real(kind=8), parameter :: epsilon0 = 8.854d-12

    ! Derived parameters
    real(kind=8) :: mass_ion = 5. * mass_electrons
    real(kind=8) :: charge_ions = Z_ion * charge_electrons
    real(kind=8), parameter :: Initial_energy_electrons = 1000 * charge_electrons
    real(kind=8), parameter :: Initial_energy_ions = 100 * charge_electrons
    real(kind=8) :: Initial_flow_velocity_elec = dsqrt(2 * Initial_energy_electrons / mass_electrons)
    real(kind=8) :: Initial_flow_velocity_ion = 0.d0
    integer, parameter :: Number_collision_el_ion = Number_alpha * density_ions / (density_electrons + density_ions)
    integer, parameter :: Number_collision_el_el = (Number_alpha - Number_collision_el_ion) / 2.d0
    integer, parameter :: Number_collision_ion_ion = (Number_beta - Number_collision_el_ion) / 2.d0

    ! Variables
    type(particle), dimension(Number_alpha) :: array_electrons
    type(particle), dimension(Number_beta) :: array_ions
    integer, dimension(Number_alpha) :: index_electrons
    integer, dimension(Number_beta) :: index_ions
    integer :: i, j, i_tstep, number_time_step = int(t_final / Delta_t),k
    real(kind=8) :: density = density_electrons + density_ions, time
    real(kind=8), dimension(3) :: tot_velec, tot_vion, tot_velec2, tot_vion2
    real(kind=8) :: averageV_ele, averageV2_ele, averageV2_ion, averageV_ion, Temp_elec, Temp_ion, sum_vel, sum_vion
    integer :: i_seed = 123456789
    real(kind=8) :: random_Gauss(3)
    

    ! Open output file
    open(10, file="CeciEstUnTest7_1.txt", status="replace")

    ! Initialize indices
    do i = 1, Number_alpha
        index_electrons(i) = i
    end do
    do i = 1, Number_beta
        index_ions(i) = i
    end do

    ! Initialize particles
    array_electrons%mass = mass_electrons
    array_electrons%charge = charge_electrons
    array_ions%mass = mass_ion
    array_ions%charge = charge_ions

    do i = 1, Number_alpha
        do j = 1, 3
            call rran2(random_Gauss(j), i_seed)
        end do
        array_electrons(i)%velocity = maxwell_boltzmann(Initial_flow_velocity_elec, dsqrt(3.d0*Initial_energy_electrons/mass_electrons), random_Gauss)
        do j = 1, 3
            call rran2(random_Gauss(j), i_seed)
        end do
        array_ions(i)%velocity = maxwell_boltzmann(Initial_flow_velocity_ion, dsqrt(3.d0*Initial_energy_ions/mass_ion), random_Gauss)
    end do

    ! Time-stepping loop
    do i_tstep = 0, number_time_step
        time = i_tstep * Delta_t

        ! Reset averages: averageV=<v>, averageV2=<v^2>
        averageV_ele = 0.d0
        averageV2_ele = 0.d0
        averageV2_ion = 0.d0
        averageV_ion = 0.d0

        ! Randomize particle indices
        call random_permutation(Number_alpha, index_electrons)
        call random_permutation(Number_beta, index_ions)


        ! Ion-electron collisions
        do i = 1, Number_collision_el_ion
            if (i>Number_alpha .or. i>Number_beta) exit
            ! Compute averages
            averageV_ele = averageV_ele + (1.d0 / Number_alpha) * dsqrt(max(0.d0, sum(array_electrons(i)%velocity**2)))
            averageV_ion = averageV_ion + (1.d0 / Number_beta) * dsqrt(max(0.0, sum(array_ions(i)%velocity**2)))
            averageV2_ele = averageV2_ele + (1.d0 / Number_alpha) * sum(array_electrons(i)%velocity**2)
            averageV2_ion = averageV2_ion + (1.d0 / Number_beta) * sum(array_ions(i)%velocity**2)
            call collision(array_electrons(index_electrons(i)), array_ions(index_ions(i)))
        enddo


        ! Ion-ion collisions
        do i = Number_collision_el_ion + 1, Number_collision_el_ion + Number_collision_ion_ion
            if (i>Number_beta) exit
            k = i + Number_collision_ion_ion
            if (k>Number_beta) k = k-Number_beta
            averageV_ion = averageV_ion + (1.d0 / 2.d0 / Number_beta) * (dsqrt(sum(array_ions(i)%velocity**2)) +  dsqrt(sum(array_ions(k)%velocity**2)))
            averageV2_ion = averageV2_ion + (1.d0/ 2.d0 / Number_beta) * (sum(array_ions(i)%velocity**2) + sum(array_ions(k)%velocity**2))
            call collision(array_ions(index_ions(i)), array_ions(index_ions(k)))
        end do

        
        
        ! Electron-electron collisions
        do i = Number_collision_el_ion + 1, Number_collision_el_ion + Number_collision_el_el  
            if (i>Number_alpha) exit
            k = i + Number_collision_el_el
            if (k>Number_alpha) k = k-Number_alpha
            averageV_ele = averageV_ele + (1.d0 / 2.d0 / Number_alpha) * (dsqrt(sum(array_electrons(i)%velocity**2)) +  dsqrt(sum(array_electrons(k)%velocity**2)))
            averageV2_ele = averageV2_ele + (1.d0/ 2.d0 / Number_alpha) * (sum(array_electrons(i)%velocity**2) + sum(array_electrons(k)%velocity**2))
            call collision(array_electrons(index_electrons(i)), array_electrons(index_electrons(k)))   
        end do

        tot_velec = 0.d0
        tot_vion = 0.d0
        tot_velec2 = 0.d0
        tot_vion2 = 0.d0


        ! Compute total velocity vectors (average velocity vector, parallelized)
        do i = 1, Number_alpha
            tot_velec = tot_velec + array_electrons(i)%velocity
            tot_velec2 = tot_velec2 + array_electrons(i)%velocity**2
        end do
        tot_velec = tot_velec / Number_alpha
        tot_velec2 = tot_velec2 / Number_alpha

        
        do i = 1, Number_beta
            tot_vion = tot_vion + array_ions(i)%velocity
            tot_vion2 = tot_vion2 + array_ions(i)%velocity**2
        end do
        tot_vion = tot_vion / Number_beta
        tot_vion2 = tot_vion2 / Number_beta

        Temp_elec = mass_electrons * (sum(tot_velec2)-sum(tot_velec**2)) / (3.d0 * 1.6d-19)
        Temp_ion = mass_ion * (sum(tot_vion2)-sum(tot_vion**2)) / (3.d0 * charge_electrons)
        sum_vel = dsqrt(max(0.0, sum(tot_velec**2)))
        sum_vion = dsqrt(max(0.0, sum(tot_vion**2)))
        write(10, '(5EN14.4)') time, Temp_elec, Temp_ion, sum_vel / Initial_flow_velocity_elec, sum_vion / Initial_flow_velocity_elec
    end do

    close(10)

    contains

    function maxwell_boltzmann(flow_speed, v_th, rnd) result(velocity)

        implicit none
        real(kind=8), intent(in) :: flow_speed, v_th
        real(kind=8), intent(in) :: rnd(3)  ! Array of random numbers
        real(kind=8), dimension(3) :: velocity  ! Resulting 3D velocity vector
        real(kind=8) :: theta, phi, vp
        
        ! Gaussian loading
        vp = v_th * dsqrt(-dlog(max(1.0d-12, rnd(1))))  ! Magnitude of velocity
        
        ! Random angles for spherical coordinates
        theta = 2.0d0 * pi * rnd(2)             ! Azimuthal angle
        phi = acos(2.0d0 * rnd(3) - 1.0d0)      ! Polar angle
        
        ! Calculate 3D velocity components
        velocity(1) = flow_speed + vp * sin(phi) * cos(theta)  ! Add flow velocity to x-component
        velocity(2) = vp * sin(phi) * sin(theta)
        velocity(3) = vp * cos(phi)
        
        return
    end function maxwell_boltzmann

    real(kind=8) function reduced_mass(m1, m2)
        real(kind=8), intent(in) :: m1, m2
        reduced_mass = m1 * m2 / (m1 + m2)
    end function reduced_mass

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
        real(kind=8) :: A0 = 10.d0, tol = 1.e-6
        integer :: max_iter = 10
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


    subroutine collision(particle1, particle2)
        implicit none
        type(particle), intent(inout) :: particle1, particle2
        real(kind=8), dimension(3) :: g, h
        real(kind=8) :: g_perp, angle_azimuthal, rand1, A_tau, argument_acos, angle_diffusion
        real(kind=8) :: A, tau

        g = particle1%velocity-particle2%velocity  
            
        A = density/(4.*pi) * (particle1%charge*particle2%charge/(epsilon0*reduced_mass(particle1%mass,particle2%mass)))**2 * CoulombLog    
        tau = A * Delta_t / max(1d-18, dsqrt(sum(g**2))**3)
        g_perp = dsqrt(max(0.0, g(2)**2 + g(3)**2))
        call rran2(angle_azimuthal,i_seed)
        call rran2(rand1,i_seed)

        ! Use either the NewtonRaphson method or the approximated analytical solution

        ! A_tau = NewtRaph(tau)
        if (tau <= 0.01d0) then
            A_tau = 1.d0 / max(tau, 0.0000000001d0)
        else if (tau <= 0.1d0) then
            A_tau = (84.708d0 * tau + 166.4288d0) / (166.4681d0 * tau - 0.0003d0)
        else if (tau <= 0.6d0) then
            A_tau = 1.2528d0 * tau**(-0.9208d0) + 0.1115d0 * tau
        else if (tau < 2.d0) then
            A_tau = (-50.8179d0 * tau + 181.7747d0) / (88.0921d0 * tau + 20.5056d0)
        else
            A_tau = 3.d0 * dexp(-min(tau, 40.d0))
        end if


        angle_azimuthal = 2.d0 * pi * angle_azimuthal
        h(1) = g_perp * dcos(angle_azimuthal)
        h(2) = -(g(1) * g(2) * dcos(angle_azimuthal) + dsqrt(max(0.0, sum(g**2))) * g(3) * dsin(angle_azimuthal)) / max(1d-18, g_perp)
        h(3) = -(g(1) * g(3) * dcos(angle_azimuthal) - dsqrt(max(0.0, sum(g**2))) * g(2) * dsin(angle_azimuthal)) / max(1d-18, g_perp)

        argument_acos = dlog(max(1d-18, dexp(-A_tau) + 2.0 * rand1 * dsinh(A_tau))) / A_tau
        argument_acos = max(-1.d0, min(1.d0, argument_acos))
        angle_diffusion = dacos(argument_acos)
        particle1%velocity = particle1%velocity - particle2%mass/(particle1%mass+particle2%mass) * (g * (1. - dcos(angle_diffusion)) + h * dsin(angle_diffusion))
        particle2%velocity = particle2%velocity + particle1%mass/(particle1%mass+particle2%mass) * (g * (1. - dcos(angle_diffusion)) + h * dsin(angle_diffusion))

    end subroutine collision

    subroutine random_permutation(n, indices)
        implicit none
        integer, intent(in) :: n
        integer, intent(inout) :: indices(n)
        integer :: i, j, temp
        real(kind=8) :: r
    
        do i = n, 2, -1
            ! call random_number(r)
            call rran2(r,i_seed)
            j = 1 + int(r * i)  ! Generate random index in range [1, i]
    
            ! Swap indices(i) and indices(j)
            temp = indices(i)
            indices(i) = indices(j)
            indices(j) = temp
        end do
    end subroutine random_permutation

end program CoulombCollisions