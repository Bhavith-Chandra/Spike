from .constants import g, pi, epsilon_0, k, R, Avogadro_number
import numpy as np
class AdvancedPhysicsMathChem:

    @staticmethod
    def weight(mass):
        return mass * g

    @staticmethod
    def density(mass, vol):
        return mass / vol if vol != 0 else ValueError("Division by zero encountered")

    @staticmethod
    def volume(length, width, height):
        return length * width * height

    @staticmethod
    def area(length, width):
        return length * width

    @staticmethod
    def surface_area_sphere(radius):
        return 4 * pi * radius**2

    @staticmethod
    def volume_cylinder(radius, height):
        return pi * radius**2 * height

    @staticmethod
    def kinetic_energy(mass, velocity):
        return 0.5 * mass * velocity**2

    @staticmethod
    def velocity_from_height_gravity(height):
        return (2 * g * height)**0.5

    @staticmethod
    def moles_to_mass(moles, molar_mass):
        return moles * molar_mass

    @staticmethod
    def mass_to_moles(mass, molar_mass):
        return mass / molar_mass if molar_mass != 0 else ValueError("Division by zero encountered")

    @staticmethod
    def gravitational_force(m1, m2, distance):
        return (g * m1 * m2) / distance**2

    @staticmethod
    def electric_force(charge1, charge2, distance):
        return k * (charge1 * charge2) / distance**2

    @staticmethod
    def gravitational_potential_energy(mass, height):
        return mass * g * height

    @staticmethod
    def electric_potential_energy(charge, voltage):
        return charge * voltage

    @staticmethod
    def ideal_gas_law(pressure, volume, temperature):
        return (pressure * volume) / (R * temperature)

    @staticmethod
    def capacitance_parallel_plate(capacitance_0, distance):
        return epsilon_0 * (capacitance_0 / distance)

    @staticmethod
    def capacitance_cylinder(radius, length):
        return 2 * pi * epsilon_0 * length / (np.log(2 * radius / R))

    @staticmethod
    def ideal_gas_moles(pressure, volume, temperature):
        return (pressure * volume) / (R * temperature)

    @staticmethod
    def number_of_particles_in_gas(pressure, volume, temperature):
        return (pressure * volume) / (R * temperature) * Avogadro_number

    @staticmethod
    def torque_magnetic_dipole(magnetic_moment, magnetic_field, angle):
        return magnetic_moment * magnetic_field * np.sin(angle)

    @staticmethod
    def angular_velocity_rotational_kinetic_energy(inertia, angular_velocity):
        return 0.5 * inertia * angular_velocity**2

    @staticmethod
    def gravitational_potential_energy_orbit(mass, celestial_mass, radius):
        return (-g * mass * celestial_mass) / radius

    @staticmethod
    def work_done_gas_isothermal_expansion(n, R, temperature_initial, volume_initial, volume_final):
        return n * R * temperature_initial * np.log(volume_final / volume_initial)

    @staticmethod
    def electrostatic_force_between_charged_spheres(charge1, charge2, distance):
        return k * (charge1 * charge2) / distance**2

    @staticmethod
    def dipole_potential_electric_field(p, distance, angle):
        return (k * p * np.cos(angle)) / distance**2

    @staticmethod
    def half_life_radioactive_decay(initial_amount, decay_constant):
        return np.log(2) / decay_constant

    @staticmethod
    def heat_transfer_conduction(thermal_conductivity, area, temperature_difference, thickness):
        return (thermal_conductivity * area * temperature_difference) / thickness

    @staticmethod
    def reaction_rate_arrhenius(A, activation_energy, temperature):
        return A * np.exp(-activation_energy / (R * temperature))