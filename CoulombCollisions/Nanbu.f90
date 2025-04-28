module Part
    real, parameter :: pi = acos(-1.0)
    type :: particle
        real(8) :: mass, charge
        real(8), dimension(3) :: velocity
    end type

    contains

    function maxwell_boltzmann(flow_speed, e0, mass) result(array_v)
        implicit none
        real(8), intent(in) :: flow_speed, e0, mass
        real(8), dimension(3) :: array_v
        real(8) :: v_th, x1, x2, gauss1, gauss2, gauss3

        ! Correct thermal velocity (v_rms)
        v_th = sqrt(2.0 * e0 / mass)

        ! Generate Gaussian-distributed random numbers using Box-Muller transform
        call random_number(x1)
        call random_number(x2)
        gauss1 = sqrt(-2.0 * log(x1)) * cos(2.0 * pi * x2)
        gauss2 = sqrt(-2.0 * log(x1)) * sin(2.0 * pi * x2)
        call random_number(x1)
        call random_number(x2)
        gauss3 = sqrt(-2.0 * log(x1)) * cos(2.0 * pi * x2)

        ! Assign Maxwell-Boltzmann velocity components
        array_v(1) = flow_speed + v_th * gauss1
        array_v(2) = v_th * gauss2
        array_v(3) = v_th * gauss3
    end function maxwell_boltzmann
end module Part

subroutine random_permutation(n, indices)
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: indices(n)
    integer :: i, j, temp
    real :: r

    do i = n, 2, -1
        call random_number(r)
        j = 1 + int(r * i)  ! Generate random index in range [1, i]

        ! Swap indices(i) and indices(j)
        temp = indices(i)
        indices(i) = indices(j)
        indices(j) = temp
    end do
end subroutine random_permutation

program CoulombCollisions
    use Part
    implicit none

    ! Parameters
    real(8), parameter :: e = exp(1.0)
    integer, parameter :: Number_alpha = 10**5, Number_beta = 10**5, Z_ion = 1
    real(8), parameter :: Delta_t = 1.e-7, density_ions = 1.e21, CoulombLog = 15.9
    real(8), parameter :: mass_electrons = 9.109e-31, charge_electrons = 1.602e-19
    real(8), parameter :: density_electrons = Z_ion * density_ions, k_B = 1.380649e-23, t_final = 8.e-5
    real(8), parameter :: Initial_energy_electrons = 1000 * charge_electrons
    real(8), parameter :: Initial_energy_ions = 100 * charge_electrons
    real(8), parameter :: epsilon0 = 8.854e-12

    ! Derived parameters
    real(8) :: mass_ion = 5 * mass_electrons
    real(8) :: charge_ions = Z_ion * charge_electrons
    real(8) :: Initial_flow_velocity_elec = sqrt(2 * Initial_energy_electrons / mass_electrons)
    real(8) :: Initial_flow_velocity_ion = 0.0
    integer, parameter :: Number_collision_el_ion = Number_alpha * density_ions / (density_electrons + density_ions)
    integer, parameter :: Number_collision_el_el = (Number_alpha - Number_collision_el_ion) / 2
    integer, parameter :: Number_collision_ion_ion = (Number_beta - Number_collision_el_ion) / 2

    ! Variables
    type(particle), dimension(Number_alpha) :: array_electrons
    type(particle), dimension(Number_beta) :: array_ions
    integer, dimension(Number_alpha) :: index_electrons
    integer, dimension(Number_beta) :: index_ions
    integer :: i, i_tstep, k, number_time_step = int(t_final / Delta_t)
    real(8) :: density = density_electrons + density_ions, time
    real(8), dimension(3) :: g, h, tot_velec, tot_vion
    integer, dimension(6, Number_alpha) :: resj
    real(8) :: rand1, angle_diffusion, angle_azimuthal, g_perp, A, tau, A_tau, yo
    real(8) :: averageV_ele, averageV2_ele, averageV2_ion, averageV_ion, Temp_elec, Temp_ion, sum_vel, sum_vion

    ! Open output file
    open(10, file="OnTestLa.txt", status="replace")

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
        array_electrons(i)%velocity = maxwell_boltzmann(Initial_flow_velocity_elec, Initial_energy_electrons, mass_electrons)
        array_ions(i)%velocity = maxwell_boltzmann(Initial_flow_velocity_ion, Initial_energy_ions, mass_ion)
    end do

    write(10, '(5EN14.4)') 0.0, 1000.0, 100.0, 1.0, 0.0

    ! Time-stepping loop
    do i_tstep = 1, number_time_step
        time = i_tstep * Delta_t

        ! Reset averages
        averageV_ele = 0.0
        averageV2_ele = 0.0
        averageV2_ion = 0.0
        averageV_ion = 0.0

        tot_velec = 0.0
        tot_vion = 0.0

        ! Compute total velocities
        do i = 1, Number_alpha
            tot_velec = tot_velec + (1.0 / Number_alpha) * array_electrons(i)%velocity
            tot_vion = tot_vion + (1.0 / Number_beta) * array_ions(i)%velocity
        end do

        ! Randomize particle indices
        call random_permutation(Number_alpha, index_electrons)
        call random_permutation(Number_beta, index_ions)

        resj(i_tstep + 1, :) = index_electrons

        ! Collision loop
        do i = 1, int(3 * Number_alpha / 4)
            ! Compute averages
            averageV_ele = averageV_ele + (1.0 / Number_alpha) * sqrt(sum(array_electrons(i)%velocity**2))
            averageV_ion = averageV_ion + (1.0 / Number_beta) * sqrt(sum(array_ions(i)%velocity**2))
            averageV2_ele = averageV2_ele + (1.0 / Number_alpha) * sum(array_electrons(i)%velocity**2)
            averageV2_ion = averageV2_ion + (1.0 / Number_beta) * sum(array_ions(i)%velocity**2)

            ! Electron-ion collisions
            if (i <= Number_collision_el_ion) then
                g = array_electrons(index_electrons(i))%velocity - array_ions(index_ions(i))%velocity
                A = density / (4.0 * pi) * (charge_electrons * charge_ions / (epsilon0 * reduced_mass(mass_electrons, mass_ion)))**2 * CoulombLog
                tau = A * Delta_t / max(1e-12, sqrt(sum(g**2))**3)
                g_perp = sqrt(g(2)**2+g(3)**2)
                call random_number(angle_azimuthal)
                call random_number(rand1)
                A_tau = f2(tau)
                if (tau<=0.02) A_tau = f1(tau)
                if (tau>=2.) A_tau = f3(tau)
                
                ! A_tau = Newton_Raphson(tau)
                angle_azimuthal = 2*pi*angle_azimuthal

                h(1) = g_perp*cos(angle_azimuthal)
                h(2) = -(g(1)*g(2)*cos(angle_azimuthal)+sqrt(sum(g**2))*g(3)*sin(angle_azimuthal))/g_perp
                h(3) = -(g(1)*g(3)*cos(angle_azimuthal)-sqrt(sum(g**2))*g(2)*sin(angle_azimuthal))/g_perp
                yo = log(exp(-A_tau)+2.*rand1*sinh(A_tau))/A_tau
                yo = max(-1.0, min(1.0, yo))

                angle_diffusion = acos(yo)
                
                array_electrons(index_electrons(i))%velocity = array_electrons(index_electrons(i))%velocity - mass_ion/(mass_ion+mass_electrons) * (g*(1.-cos(angle_diffusion))+h*sin(angle_diffusion))
                array_ions(index_ions(i))%velocity = array_ions(index_ions(i))%velocity + mass_electrons/(mass_ion+mass_electrons) * (g*(1.-cos(angle_diffusion))+h*sin(angle_diffusion))
            else
                

                ! Ions
                k = index_ions(i+ Number_collision_ion_ion) 
                g = array_ions(index_ions(i))%velocity-array_ions(k)%velocity
                
                A = density/(4.*pi) * (charge_ions*charge_ions/(epsilon0*reduced_mass(mass_ion,mass_ion)))**2 * CoulombLog
                
                tau = A * Delta_t / (max(1e-12, sqrt(g(1)**2 + g(2)**2 + g(3)**2)**3))
                g_perp = sqrt(g(2)**2+g(3)**2)
                call random_number(angle_azimuthal)
                call random_number(rand1)
                
                
                
                A_tau = f2(tau)
                if (tau<=0.02) A_tau = f1(tau)
                if (tau>=2.) A_tau = f3(tau)

                ! A_tau = Newton_Raphson(tau)
                angle_azimuthal = 2*pi*angle_azimuthal
                h(1) = g_perp*cos(angle_azimuthal)
                h(2) = -(g(1)*g(2)*cos(angle_azimuthal)+sqrt(sum(g**2))*g(3)*sin(angle_azimuthal))/g_perp
                h(3) = -(g(1)*g(3)*cos(angle_azimuthal)-sqrt(sum(g**2))*g(2)*sin(angle_azimuthal))/g_perp

                yo = log(exp(-A_tau)+2.*rand1*sinh(A_tau))/A_tau
                yo = max(-1.0, min(1.0, yo))

                angle_diffusion = acos(yo)
                
                array_ions(index_ions(i))%velocity = array_ions(index_ions(i))%velocity - 0.5 * (g * (1. - cos(angle_diffusion)) + h * sin(angle_diffusion))
                array_ions(k)%velocity = array_ions(k)%velocity + 0.5 * (g * (1. - cos(angle_diffusion)) + h * sin(angle_diffusion))

                ! Electrons
                
                k = index_electrons(i + Number_collision_el_el)
                g = array_electrons(index_electrons(i))%velocity-array_electrons(k)%velocity
                
                A = density/(4.*pi) * (charge_electrons*charge_electrons/(epsilon0*reduced_mass(mass_electrons,mass_electrons)))**2 * CoulombLog
                
                tau = A * Delta_t / (max(1e-12, sqrt(g(1)**2 + g(2)**2 + g(3)**2)**3))
                g_perp = sqrt(g(2)**2+g(3)**2)
                call random_number(angle_azimuthal)
                call random_number(rand1)
                A_tau = f2(tau)
                if (tau<=0.02) A_tau = f1(tau)
                if (tau>=2.) A_tau = f3(tau)
                ! A_tau = Newton_Raphson(tau)
                angle_azimuthal = 2*pi*angle_azimuthal
                h(1) = g_perp*cos(angle_azimuthal)
                h(2) = -(g(1)*g(2)*cos(angle_azimuthal)+sqrt(sum(g**2))*g(3)*sin(angle_azimuthal))/g_perp
                h(3) = -(g(1)*g(3)*cos(angle_azimuthal)-sqrt(sum(g**2))*g(2)*sin(angle_azimuthal))/g_perp

                yo = log(exp(-A_tau)+2.*rand1*sinh(A_tau))/A_tau
                yo = max(-1.0, min(1.0, yo))
                
                angle_diffusion = acos(yo)
                
                array_electrons(index_electrons(i))%velocity = array_electrons(index_electrons(i))%velocity - 0.5 * (g * (1. - cos(angle_diffusion)) + h * sin(angle_diffusion))
                array_electrons(k)%velocity = array_electrons(k)%velocity + 0.5 * (g * (1. - cos(angle_diffusion)) + h * sin(angle_diffusion))

            end if
        
        end do

        ! Write results
        Temp_elec = mass_electrons * (averageV2_ele - averageV_ele**2) / (3.0 * 1.6e-19)
        Temp_ion = mass_ion * (averageV2_ion - averageV_ion**2) / (3.0 * charge_electrons)
        sum_vel = sqrt(sum(tot_velec**2))
        sum_vion = sqrt(sum(tot_vion**2))
        write(10, '(5EN14.4)') time, Temp_elec, Temp_ion, sum_vel / Initial_flow_velocity_elec, sum_vion / Initial_flow_velocity_elec
    end do

    close(10)

    contains

    real(8) function reduced_mass(m1, m2)
        real(8), intent(in) :: m1, m2
        reduced_mass = m1 * m2 / (m1 + m2)
    end function reduced_mass

    function f1(tau) result(output)
        real(8), intent(in) :: tau
        real(8) :: output
        output = 1.0 / tau
    end function f1

    function f2(tau) result(output)
        real(8), intent(in) :: tau
        real(8) :: output
        output = (1.0 / tau + 3.0 * exp(-tau)) / 2.0
    end function f2

    function f3(tau) result(output)
        real(8), intent(in) :: tau
        real(8) :: output
        output = 3.0 * exp(-min(tau, 40.0))
    end function f3

end program CoulombCollisions