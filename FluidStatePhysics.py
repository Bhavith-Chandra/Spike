from .constants import pi
import math


class Mechanics:

    def init(self, initial_velocity=1, acceleration=1, time=1, final_velocity=1,
                 mass=1, velocity=1, force=1, displacement=1, angle=0, work=1,
                 radius=1, period=1, lever_arm=1, angular_displacement=1,
                 angular_velocity=1):
        self.initial_velocity = initial_velocity
        self.acceleration = acceleration
        self.time = time
        self.final_velocity = final_velocity
        self.mass = mass
        self.velocity = velocity
        self.force = force
        self.displacement = displacement
        self.angle = angle
        self.work = work
        self.radius = radius
        self.period = period
        self.lever_arm = lever_arm
        self.angular_displacement = angular_displacement
        self.angular_velocity = angular_velocity

    def velocity(self):
        return self.initial_velocity + self.acceleration * self.time

    def displacement(self):
        return self.initial_velocity * self.time + 0.5 * self.acceleration * \
               self.time ** 2

    def acceleration(self):
        return (self.final_velocity - self.initial_velocity) / self.time

    def uniform_accelerated_motion(self):
        displacement = self.initial_velocity * self.time + 0.5 * self.acceleration * \
                       self.time ** 2
        final_velocity = self.initial_velocity + self.acceleration * self.time
        return displacement, final_velocity

    def force(self):
        return self.mass * self.acceleration

    def work(self):
        return self.force * self.displacement * math.cos(math.radians(self.angle))

    def kinetic_energy(self):
        return 0.5 * self.mass * self.velocity ** 2

    def potential_energy(self, height, gravitational_field_strength):
        return self.mass * height * gravitational_field_strength

    def power(self):
        return self.work / self.time

    def momentum(self):
        return self.mass * self.velocity

    def impulse(self):
        return self.force * self.time

    def circular_velocity(self):
        return 2 * pi * self.radius / self.period

    def centripetal_acceleration(self):
        return self.velocity ** 2 / self.radius

    def torque(self):
        return self.force * self.lever_arm

    def angular_velocity(self):
        return self.angular_displacement / self.time

    def angular_acceleration(self):
        return self.angular_velocity / self.time
    
    def gravitational_force(self):
        gravitational_constant = 6.67430e-11
        return (gravitational_constant * self.mass1 * self.mass2) / self.distance**2

    def elastic_potential_energy(self):
        return 0.5 * self.spring_constant * self.displacement**2

    def angular_momentum(self):
        return self.radius * self.momentum()

    def angular_impulse(self):
        return self.torque() * self.time

    def rotational_kinetic_energy(self):
        return 0.5 * self.radius*2 * self.angular_velocity*2

    def gravitational_potential_energy(self):
        gravitational_constant = 6.67430e-11
        return (-gravitational_constant * self.mass * self.mass_of_earth) / self.radius

    def simple_pendulum_period(self):
        return 2 * pi * math.sqrt(self.length / self.gravitational_field_strength)
    
    def electric_force(self):
        electric_constant = 8.9875e9
        return (electric_constant * self.charge1 * self.charge2) / self.distance**2

    def electric_potential_energy(self):
        return self.charge * self.electric_potential

    def electric_power(self):
        return self.electric_current * self.electric_potential

    def magnetic_force(self):
        magnetic_constant = 4e-7 * pi
        return (magnetic_constant * self.charge * self.velocity * self.magnetic_field) / math.sin(math.radians(self.angle))

    def magnetic_torque(self):
        return self.magnetic_moment * self.magnetic_field * math.sin(math.radians(self.angle))

    def magnetic_flux(self):
        return self.magnetic_field * self.area * math.cos(math.radians(self.angle))

    def resistive_power(self):
        return self.current**2 * self.resistance

    def capacitive_reactance(self):
        return 1 / (2 * pi * self.frequency * self.capacitance)

    def inductive_reactance(self):
        return 2 * pi * self.frequency * self.inductance

    def impedance(self):
        return math.sqrt(self.resistance**2 + (self.inductive_reactance - self.capacitive_reactance)**2)

    def power_factor(self):
        return math.cos(math.atan((self.inductive_reactance - self.capacitive_reactance) / self.resistance))

    def doppler_effect_frequency(self):
        return self.source_frequency * (self.speed_of_sound + self.observer_speed) / (self.speed_of_sound - self.source_speed)

    def escape_velocity(self):
        return math.sqrt(2 * self.gravitational_constant * self.mass / self.radius)

    def fluid_pressure(self):
        return self.density * self.gravitational_field_strength * self.depth

    def fluid_flow_rate(self):
        return self.cross_sectional_area * self.velocity

    def Bernoulli_equation(self):
        return self.pressure + 0.5 * self.density * self.velocity**2 + self.density * self.gravitational_field_strength * self.height

    def Doppler_effect_velocity(self):
        return self.source_speed * (self.speed_of_sound + self.observer_speed) / (self.speed_of_sound - self.source_speed)

    def electric_circuit_power(self):
        return self.voltage * self.electric_current * math.cos(math.radians(self.power_factor_angle))

    def torque_due_to_gravity(self):
        return self.mass * self.gravitational_field_strength * self.lever_arm

    def magnetic_induction(self):
        return self.magnetic_flux / self.area

    def specific_heat_capacity(self):
        return self.heat_added / (self.mass * (self.final_temperature - self.initial_temperature))

    def projectile_motion_range(self):
        return (self.initial_velocity**2 * math.sin(2 * math.radians(self.angle))) / self.gravitational_field_strength

    def Lorentz_force(self):
        magnetic_force = self.charge * (self.velocity * self.magnetic_field)
        electric_force = self.charge * self.electric_field
        return math.sqrt(magnetic_force*2 + electric_force*2)

    def quantum_mechanics_probability_density(self):
        wavefunction_squared = abs(self.wavefunction)**2
        return wavefunction_squared

    def chaos_theory_lyapunov_exponent(self):
        return 1 / self.time * abs(self.chaos_function_derivative)

    def neural_network_activation(self):
        return 1 / (1 + math.exp(-self.weighted_sum))

    def game_theory_nash_equilibrium(self):
        return self.strategy_optimal_response()

    def information_theory_entropy(self):
        return -sum(p * math.log2(p) for p in self.probabilities)

    def elastic_collision_velocity(self):
        velocity_after_collision_1 = (self.velocity1 * (self.mass1 - self.mass2) + 2 * self.mass2 * self.velocity2) / (self.mass1 + self.mass2)
        velocity_after_collision_2 = (self.velocity2 * (self.mass2 - self.mass1) + 2 * self.mass1 * self.velocity1) / (self.mass1 + self.mass2)
        return velocity_after_collision_1, velocity_after_collision_2

    def fluid_static_pressure(self):
        return self.density * self.gravitational_field_strength * self.depth

    def fluid_velocity_from_torricelli(self):
        return math.sqrt(2 * self.gravitational_field_strength * self.height_difference)

    def dynamic_equilibrium_tension(self):
        return self.mass * (self.gravity_acceleration + self.circular_velocity**2 / self.radius)

    def escape_speed(self):
        return math.sqrt(2 * self.gravitational_constant * self.mass / self.distance_from_center)

    def Hookes_law_stress(self):
        return self.elastic_modulus * self.strain

    def fluid_viscosity_shear_rate(self):
        return self.shear_stress / self.viscosity

    def escape_velocity_orbit(self):
        return math.sqrt(self.gravitational_constant * self.mass / self.distance_from_center)

    def wave_interference_beats(self):
        beat_frequency = abs(self.frequency_source1 - self.frequency_source2)
        return beat_frequency

    def fluid_surface_tension(self):
        return self.surface_tension / self.radius

    def rotational_kinetic_energy_hoop(self):
        return 0.5 * self.mass * self.radius*2 * self.angular_velocity*2

    def moment_of_inertia_disk(self):
        return 0.5 * self.mass * self.radius**2

    def magnetic_force_on_current(self):
        return self.current * self.magnetic_field * self.length

    def magnetic_dipole_moment(self):
        return self.current * self.area_normal_to_current

    def projectile_motion_height(self):
        return (self.initial_velocity**2 * math.sin(math.radians(2 * self.angle))) / (2 * self.gravitational_field_strength)

    def electric_potential_energy_capacitor(self):
        return 0.5 * self.capacitance * self.electric_potential**2

    def rotational_dynamics_torque_moment_of_inertia(self):
        angular_acceleration = self.torque / self.moment_of_inertia_rotational
        return angular_acceleration

    def fluid_flow_rate_poisson_equation(self):
        poiseuille_constant = 8 * self.viscosity * self.length / (math.pi * self.radius**4)
        return poiseuille_constant * (self.pressure_difference / self.viscosity)

    def harmonic_oscillator_amplitude(self):
        angular_frequency = 2 * pi * self.frequency_harmonic_oscillator
        return self.initial_amplitude * math.exp(-self.damping_ratio * angular_frequency * self.time)

    def thermodynamic_work(self):
        return self.pressure * self.volume_change

    def fluid_dynamics_boundary_layer_thickness(self):
        reynolds_number = self.density * self.velocity * self.characteristic_length / self.viscosity
        return 0.37 * self.characteristic_length / math.sqrt(reynolds_number)

    def fluid_pressure_hydrostatic(self):
        return self.density * self.gravitational_field_strength * self.depth

    def magnetic_flux_density(self):
        magnetic_flux = self.magnetic_flux_density * self.area
        return magnetic_flux

    def harmonic_oscillator_velocity(self):
        angular_frequency = 2 * pi * self.frequency_harmonic_oscillator
        return self.initial_amplitude * angular_frequency * math.exp(-self.damping_ratio * angular_frequency * self.time)

    def electric_current_density(self):
        return self.electric_current / self.cross_sectional_area

    def potential_energy_spring(self):
        return 0.5 * self.spring_constant * self.displacement**2

    def fluid_dynamics_bernoulli_principle(self):
        return self.pressure + 0.5 * self.density * self.velocity**2 + self.density * self.gravitational_field_strength * self.height

    def electric_circuit_resistance_series(self):
        total_resistance = sum(self.resistances)
        return total_resistance

    def elastic_potential_energy_compression(self):
        return 0.5 * self.spring_constant * self.compression**2

    def fluid_viscosity_critical_speed(self):
        return math.sqrt(self.radius * self.gravitational_field_strength / self.density)

    def work_energy_theorem(self):
        kinetic_energy_final = 0.5 * self.mass * self.final_velocity**2
        kinetic_energy_initial = 0.5 * self.mass * self.initial_velocity**2
        return kinetic_energy_final - kinetic_energy_initial

    def fluid_dynamics_velocity_head(self):
        return self.velocity**2 / (2 * self.gravitational_field_strength)