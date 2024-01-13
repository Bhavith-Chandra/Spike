from .constants import pi, epsilon0, vonklitzing_constant, plancks_constant, \
    compton_wavelength
from .constants import finestructure_constant, Coulombs_constant, elementary_charge, \
    speed_of_light
import math
import numpy as np
from .Mass import Mass
c = 3 * (10**8)

class Electromagnetism:

    def __init__(self, current=1, voltage=1, resistance=1, area=1, length=1, distance=1,
                 charge=1, force=1, magnetic_flux=1, initial_wavelength=1,
                 scattering_angle=1, magnetic_flux_density=1, magnetic_field_strength=1,
                 velocity=1, charge_density=1, electron_density=1, induced_emf=1,
                 time=1, number_of_turns=1, inductance=1, current_rate_of_change=1,
                 charge1=1, charge2=1, principal_quantum_number=1, momentum=1,
                 fermi_level=1, mass=1, elementary_charge1=1, change_in_magnetic_flux=1,
                 magnetic_field=1, gravitational_field_strength=1, mu0=1, n1=1, n2=1,
                 angle_of_incidence=1, hv=1, threshold_frequency=1, work_function=1,
                 stopping_potential=1, initial_quantity=1, decay_constant=1,
                 mass_defect=1, final_quantity=1) -> None:
        self.current = current
        self.voltage = voltage
        self.resistance = resistance
        self.area = area
        self.length = length
        self.distance = distance
        self.charge = charge
        self.force = force
        self.magnetic_flux = magnetic_flux
        self.initial_wavelength = initial_wavelength
        self.scattering_angle = scattering_angle
        self.magnetic_flux_density = magnetic_flux_density
        self.magnetic_field_strength = magnetic_field_strength
        self.velocity = velocity
        self.charge_density = charge_density
        self.electron_density = electron_density
        self.induced_emf = induced_emf
        self.time = time
        self.number_of_turns = number_of_turns
        self.inductance = inductance
        self.current_rate_of_change = current_rate_of_change
        self.charge1 = charge1
        self.charge2 = charge2
        self.principal_quantum_number = principal_quantum_number
        self.momentum = momentum
        self.fermi_level = fermi_level
        self.mass = mass
        self.elementary_charge1 = elementary_charge1
        self.change_in_magnetic_flux = change_in_magnetic_flux
        self.magnetic_field = magnetic_field
        self.gravitational_field_strength = gravitational_field_strength
        self.mu0 = mu0
        self.n1 = n1
        self.n2 = n2
        self.angle_of_incidence = angle_of_incidence
        self.hv = hv
        self.threshold_frequency = threshold_frequency
        self.work_function = work_function
        self.stopping_potential = stopping_potential
        self.initial_quantity = initial_quantity
        self.decay_constant = decay_constant
        self.mass_defect = mass_defect
        self.final_quantity = final_quantity

    def resistivity(self):
        return (self.resistance * self.area) / self.length

    def electric_field(self):
        return self.voltage / self.distance

    def electric_potential_energy(self):
        return self.charge * self.voltage

    def magnetic_field_strength(self):
        return self.force / (self.length * self.current)

    def magnetic_flux_density(self):
        return self.magnetic_flux / self.area

    def magnetic_flux(self):
        return self.magnetic_flux_density * self.area

    def magnetic_force(self):
        return self.magnetic_field_strength * self.length * self.current

    def lorentz_force(self):
        return self.charge * self.magnetic_field_strength * self.velocity

    def hall_voltage(self):
        return (self.magnetic_field_strength * self.current) / (self.charge_density *
                                                                self.electron_density)

    def faradays_law(self):
        return self.induced_emf * self.time / self.number_of_turns

    def self_inductance(self):
        return self.inductance * self.current_rate_of_change

    def mutual_inductance(self):
        return self.induced_emf / self.current_rate_of_change

    def coulombs_law(self):
        return (1 / (4 * pi * epsilon0)) * ((self.charge1 * self.charge2) /
                                            self.distance ** 2)

    def capacitance(self):
        return self.charge / self.voltage

    def electric_field_strength(self):
        return (1 / (4 * pi * epsilon0)) * (self.charge / self.distance ** 2)

    def electric_flux_density(self):
        return self.charge / self.area

    def magnetic_flux_quantum(self):
        return self.magnetic_flux / (vonklitzing_constant / (2 * pi))

    def de_broglie_wavelength(self):
        return plancks_constant / self.momentum

    def compton_wavelength_change(self):
        return self.initial_wavelength - \
               (compton_wavelength * (1 - math.cos(self.scattering_angle)))

    def bohr_orbit_radius(self):
        return (finestructure_constant ** 2) * (Coulombs_constant *
                (plancks_constant ** 2)) / ((pi * elementary_charge ** 2) *
                (self.principal_quantum_number ** 2) * (speed_of_light ** 2))

    def fermi_energy(self):
        return self.fermi_level * elementary_charge

    def de_broglie_wavelength_matter_wave(self):
        return plancks_constant / (self.mass * self.velocity )

    def number_of_ions(self):
        return self.charge / self.elementary_charge1

    def acceleration_due_to_gravity(self):
        return self.gravitational_field_strength * self.mass

    def magnetic_field_inside_a_solenoid(self):
        return (self.mu0 * self.number_of_turns * self.current) / self.length

    def magnetic_field_around_a_wire(self, radius):
        return (self.mu0 * self.current) / (2 * pi * radius)

    def induced_emf(self):
        return -self.change_in_magnetic_flux / self.time

    def maxwell_equation(self):
        return np.cross(self.magnetic_field, self.magnetic_flux_density) - \
               Coulombs_constant * self.magnetic_field

    def snells_law(self):
        return (self.n1 * math.sin(self.angle_of_incidence)) / self.n2

    def critical_angle(self):
        return math.asin(self.n2 / self.n1)

    def photoelectric_effect_work_function(self):
        return self.hv - self.threshold_frequency

    def photoelectric_effect_max_velocity(self):
        return math.sqrt((2 * elementary_charge *
                        (self.stopping_potential - self.work_function)) / Mass.electron)

    def decay_law(self):
        return self.initial_quantity * math.exp(-self.decay_constant * self.time)

    def half_life(self):
        return math.log(2) / self.decay_constant

    def nuclear_binding_energy(self):
        return self.mass_defect * speed_of_light ** 2

    def radioactive_decay(self):
        return -self.initial_quantity * math.log(self.final_quantity /
                                                 self.initial_quantity) / self.time

    def einstein_mass_energy_equivalence(self):
        return self.mass * speed_of_light ** 2

    def coulomb_law(self):
        return (1 / (4 * pi * epsilon0)) * ((self.charge1 * self.charge2) /
                                            self.distance ** 2)
    def mutual_inductance_magnetic_materials(self, permeability, area):
        """
        Calculate the mutual inductance for magnetic materials.
        Formula: Mutual Inductance = (permeability * area) / length
        """
        return (permeability * area) / self.length


    def skin_depth(self, conductivity, frequency, magnetic_permeability):
        """
        Calculate the skin depth in a conducting material.
        Formula: Skin Depth = sqrt(2 / (pi * frequency * conductivity * magnetic_permeability))
        """
        return np.sqrt(2 / (pi * frequency * conductivity * magnetic_permeability))

    def eddy_current_loss(self, resistivity, magnetic_field_frequency, magnetic_field_amplitude, thickness):
        """
        Calculate eddy current losses in a conducting material.
        Formula: Eddy Current Loss = (pi**2 * resistivity * magnetic_field_frequency**2 * magnetic_field_amplitude**2 * thickness**2) / (6 * conductivity)
        """
        return (pi**2 * resistivity * magnetic_field_frequency**2 * magnetic_field_amplitude**2 * thickness**2) / (6 * self.conductivity)

    def magnetic_circuit_reluctance(self, length, permeability, area):
        """
        Calculate the magnetic circuit reluctance.
        Formula: Reluctance = length / (permeability * area)
        """
        return length / (permeability * area)

    def magnetic_circuit_flux(self, reluctance, magnetomotive_force):
        """
        Calculate the magnetic circuit flux.
        Formula: Flux = magnetomotive_force / reluctance
        """
        return magnetomotive_force / reluctance

    def magnetic_circuit_series_connection(self, *reluctances):
        """
        Calculate the total reluctance for components in series in a magnetic circuit.
        Formula: Total Reluctance = sum(reluctances)
        """
        return sum(reluctances)

    def magnetic_circuit_parallel_connection(self, *reluctances):
        """
        Calculate the total reluctance for components in parallel in a magnetic circuit.
        Formula: Total Reluctance = 1 / sum(1/reluctances)
        """
        return 1 / np.sum(1 / np.array(reluctances))

    def wave_equation_speed_of_light(self):
        """
        Calculate the wave equation for the speed of light.
        Formula: (1 / c**2) * (d^2(E)/dt^2 - nabla^2(E)) = 0
        """
        return (1 / c**2) * (self.electric_field_rate_of_change_of_rate_of_change - 
                            np.gradient(np.gradient(self.electric_field)))

    def dirac_delta_function(self, position, a):
        """
        Calculate the Dirac delta function.
        Formula: delta(x - position) = 1 / (2 * a) if -a < x - position < a else 0
        """
        return 1 / (2 * a) if -a < self.x - position < a else 0
    def magnetic_flux_linkage(self, magnetic_field, area):
        """
        Calculate the magnetic flux linkage.
        Formula: Magnetic Flux Linkage = magnetic_field * area * cos(theta)
        """
        return magnetic_field * area * np.cos(self.scattering_angle)

    def electromagnetic_wave_speed(self, permittivity, permeability):
        """
        Calculate the speed of electromagnetic waves in a medium.
        Formula: Wave Speed = 1 / sqrt(permittivity * permeability)
        """
        return 1 / np.sqrt(permittivity * permeability)

    def poynting_vector(self, electric_field, magnetic_field):
        """
        Calculate the Poynting vector.
        Formula: S = electric_field x magnetic_field
        """
        return np.cross(electric_field, magnetic_field)

    def magnetic_dipole_moment(self, current, area, loops):
        """
        Calculate the magnetic dipole moment of a current loop.
        Formula: Magnetic Dipole Moment = current * area * loops
        """
        return current * area * loops

    def magnetic_susceptibility(self, magnetization, magnetic_field):
        """
        Calculate the magnetic susceptibility of a material.
        Formula: Magnetic Susceptibility = magnetization / magnetic_field
        """
        return magnetization / magnetic_field

    def magnetic_circuit_mmf(self, magnetic_field, length):
        """
        Calculate the magnetomotive force (MMF) in a magnetic circuit.
        Formula: MMF = magnetic_field * length
        """
        return magnetic_field * length

    def magnetic_circuit_flux_density(self, magnetic_flux, area):
        """
        Calculate the magnetic flux density in a magnetic circuit.
        Formula: Magnetic Flux Density = magnetic_flux / area
        """
        return magnetic_flux / area

    def hertzian_dipole_radiated_power(self, current, frequency, length):
        """
        Calculate the radiated power of a Hertzian dipole antenna.
        Formula: Radiated Power = (current * frequency)^2 * (length / (6 * pi * speed_of_light))
        """
        return (current * frequency)**2 * (length / (6 * pi * speed_of_light))
    def inductive_reactance(self, inductance, frequency):
        """
        Calculate the inductive reactance in an AC circuit.
        Formula: Inductive Reactance = 2 * pi * frequency * inductance
        """
        return 2 * pi * frequency * inductance

    def skin_effect_depth(self, conductivity, frequency, permeability):
        """
        Calculate the skin effect depth in a conducting material.
        Formula: Skin Effect Depth = sqrt(2 / (pi * frequency * conductivity * permeability))
        """
        return np.sqrt(2 / (pi * frequency * conductivity * permeability))

    def gyromagnetic_ratio(self, charge, mass):
        """
        Calculate the gyromagnetic ratio in magnetic resonance.
        Formula: Gyromagnetic Ratio = charge / (2 * mass)
        """
        return charge / (2 * mass)

    def magnetic_susceptibility_ferromagnetic(self, magnetization, magnetic_field_strength):
        """
        Calculate the magnetic susceptibility for ferromagnetic materials.
        Formula: Magnetic Susceptibility = (magnetization - magnetic_field_strength) / magnetic_field_strength
        """
        return (magnetization - magnetic_field_strength) / magnetic_field_strength

    def self_inductance_toroid(self, permeability, N, radius, height):
        """
        Calculate the self-inductance of a toroidal coil.
        Formula: Self-Inductance = (permeability * N**2 * pi * radius**2) / height
        """
        return (permeability * N**2 * pi * radius**2) / height

    def magnetic_dipole_moment_orbital_motion(self, charge, angular_momentum):
        """
        Calculate the magnetic dipole moment due to orbital motion of charged particle.
        Formula: Magnetic Dipole Moment = (charge * angular_momentum) / (2 * mass)
        """
        return (charge * angular_momentum) / (2 * Mass.electron)

    def magnetic_hall_effect(self, electric_field, magnetic_field, charge_carrier_density, thickness):
        """
        Calculate the Hall voltage in the Hall effect.
        Formula: Hall Voltage = (electric_field * thickness) / (magnetic_field * charge_carrier_density)
        """
        return (electric_field * thickness) / (magnetic_field * charge_carrier_density)

    def magnetic_moment_spin(self, gyromagnetic_ratio, spin_angular_momentum):
        """
        Calculate the magnetic moment due to spin angular momentum.
        Formula: Magnetic Moment = gyromagnetic_ratio * spin_angular_momentum
        """
        return gyromagnetic_ratio * spin_angular_momentum

    def magnetic_flux_transformer(self, magnetic_flux, turns_ratio):
        """
        Calculate the magnetic flux in the secondary coil of a transformer.
        Formula: Magnetic Flux (Secondary) = turns_ratio * Magnetic Flux (Primary)
        """
        return turns_ratio * magnetic_flux

    def magnetic_circuit_permeance(self, reluctance):
        """
        Calculate the permeance in a magnetic circuit.
        Formula: Permeance = 1 / Reluctance
        """
        return 1 / reluctance

    def magnetic_circuit_relay(self, flux_density, core_area):
        """
        Calculate the relay magnetomotive force (MMF) in a magnetic circuit.
        Formula: MMF = Flux Density * Core Area
        """
        return flux_density * core_area
    def magnetic_moment_orbital_motion(self, current, area, loops):
        """
        Calculate the magnetic moment due to orbital motion of current loops.
        Formula: Magnetic Moment = current * area * loops
        """
        return current * area * loops

    def magnetic_torque(self, magnetic_moment, magnetic_field, angle):
        """
        Calculate the torque experienced by a magnetic dipole in a magnetic field.
        Formula: Torque = magnetic_moment * magnetic_field * sin(angle)
        """
        return magnetic_moment * magnetic_field * np.sin(angle)

    def magnetic_flux_cylinder_coil(self, magnetic_field, radius, length, turns):
        """
        Calculate the magnetic flux through a cylindrical coil.
        Formula: Magnetic Flux = magnetic_field * pi * radius**2 * length * turns
        """
        return magnetic_field * pi * radius**2 * length * turns

    def magnetic_inductive_coupling(self, magnetic_flux_primary, turns_primary, turns_secondary):
        """
        Calculate the induced voltage in a secondary coil due to magnetic inductive coupling.
        Formula: Induced Voltage (Secondary) = -d/dt (magnetic_flux_primary * turns_primary) * turns_secondary
        """
        return -np.gradient(magnetic_flux_primary * turns_primary) * turns_secondary

    def magnetic_circuit_permanent_magnet(self, remanence, magnetic_length, area):
        """
        Calculate the magnetic flux in a permanent magnet.
        Formula: Magnetic Flux = remanence * area
        """
        return remanence * area

    def magnetic_solenoid_self_inductance(self, permeability, n, length, area):
        """
        Calculate the self-inductance of a solenoid.
        Formula: Self-Inductance = (permeability * n**2 * area) / length
        """
        return (permeability * n**2 * area) / length

    def magnetic_hysteresis_loss(self, coercivity, frequency, magnetic_field_amplitude, volume):
        """
        Calculate the hysteresis loss in a magnetic material.
        Formula: Hysteresis Loss = (2 * pi * coercivity * frequency * magnetic_field_amplitude)**2 * volume
        """
        return (2 * pi * coercivity * frequency * magnetic_field_amplitude)**2 * volume

    def magnetic_circuit_parallel_inductance(self, *inductances):
        """
        Calculate the total inductance for components in parallel in a magnetic circuit.
        Formula: Total Inductance = 1 / sum(1/inductances)
        """
        return 1 / np.sum(1 / np.array(inductances))
    
    def magnetic_energy_stored(self, inductance, current):
        """
        Calculate the magnetic energy stored in an inductor.
        Formula: Magnetic Energy = 0.5 * inductance * current**2
        """
        return 0.5 * inductance * current**2

    def magnetic_magnetic_dipole_energy(self, magnetic_moment, magnetic_field):
        """
        Calculate the potential energy of a magnetic dipole in a magnetic field.
        Formula: Magnetic Dipole Energy = -magnetic_moment * magnetic_field * cos(theta)
        """
        return -magnetic_moment * magnetic_field

    def magnetic_circuit_series_inductance(self, *inductances):
        """
        Calculate the total inductance for components in series in a magnetic circuit.
        Formula: Total Inductance = sum(inductances)
        """
        return np.sum(inductances)

    def magnetic_circuit_parallel_reluctance(self, *reluctances):
        """
        Calculate the total reluctance for components in parallel in a magnetic circuit.
        Formula: Total Reluctance = 1 / sum(1/reluctances)
        """
        return 1 / np.sum(1 / np.array(reluctances))

    def magnetic_circuit_transient_response(self, inductance, resistance, time_constant):
        """
        Calculate the transient response of a magnetic circuit.
        Formula: Transient Response = (1 - exp(-t / time_constant)) * (inductance / resistance)
        """
        return (1 - np.exp(-self.time / time_constant)) * (inductance / resistance)

    def magnetic_induction(self, magnetic_flux, area, orientation):
        """
        Calculate the magnetic induction (magnetic flux density) in a material.
        Formula: Magnetic Induction = magnetic_flux / (area * cos(theta))
        """
        return magnetic_flux / (area * np.cos(orientation))

    def magnetic_electric_field_rate_of_change(self, magnetic_flux, time):
        """
        Calculate the electric field rate of change induced by a changing magnetic flux.
        Formula: Electric Field Rate of Change = -d/dt (magnetic_flux) / time
        """
        return -np.gradient(magnetic_flux) / time

    def magnetic_relay_force(self, mmf, air_gap_length, permeability):
        """
        Calculate the force exerted by a magnetic relay.
        Formula: Force = (mmf**2 * permeability * area) / (2 * air_gap_length)
        """
        return (mmf**2 * permeability * self.area) / (2 * air_gap_length)

    def magnetic_circuit_mutual_inductance(self, magnetic_flux_secondary, current_primary):
        """
        Calculate the mutual inductance in a magnetic circuit.
        Formula: Mutual Inductance = magnetic_flux_secondary / current_primary
        """
        return magnetic_flux_secondary / current_primary

    def magnetic_circuit_magnetic_potential(self, magnetic_field, permeability, length):
        """
        Calculate the magnetic potential in a magnetic circuit.
        Formula: Magnetic Potential = magnetic_field * length / permeability
        """
        return magnetic_field * length / permeability

    def magnetic_reluctance(self, length, permeability, area):
        """
        Calculate the reluctance in a magnetic circuit.
        Formula: Reluctance = length / (permeability * area)
        """
        return length / (permeability * area)

    def magnetic_energy_stored_capacitor(self, capacitance, voltage):
        """
        Calculate the magnetic energy stored in a capacitor.
        Formula: Magnetic Energy = 0.5 * capacitance * voltage**2
        """
        return 0.5 * capacitance * voltage**2

    def magnetic_flux_mutual_inductance(self, mutual_inductance, current_secondary):
        """
        Calculate the magnetic flux induced in the secondary coil.
        Formula: Magnetic Flux (Secondary) = mutual_inductance * current_secondary
        """
        return mutual_inductance * current_secondary

    def magnetic_energy_density_inductor(self, inductance, current):
        """
        Calculate the energy density in an inductor.
        Formula: Magnetic Energy Density = 0.5 * inductance * current**2
        """
        return 0.5 * inductance * current**2

    def magnetic_circuit_flux_linkage(self, magnetic_flux, turns):
        """
        Calculate the flux linkage in a magnetic circuit.
        Formula: Flux Linkage = magnetic_flux * turns
        """
        return magnetic_flux * turns

    def magnetic_magnetic_monopole(self, magnetic_flux, surface_area):
        """
        Calculate the magnetic charge (hypothetical magnetic monopole).
        Formula: Magnetic Charge = magnetic_flux / (4 * pi * surface_area)
        """
        return magnetic_flux / (4 * pi * surface_area)
     