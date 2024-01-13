from .constants import Coulombs_constant, permeability_of_free_space, pi, elementary_charge
import math

class Electrostatics:

    def __init__(self, distance=1, charge=1, charge_1=1, charge_2=1, voltage=1, time=1,
                 current=1, force=1, velocity=1, length=1, charge_enclosed=1,
                 epsilon_0=1, emf=1, electric_field=1, magnetic_field=1,
                 resistance=1, charge_carrier_density=1, thickness=1,
                 cross_sectional_area=1, magnetic_permeability=1) -> None:
        self.distance = distance
        self.charge = charge
        self.c_1 = charge_1
        self.c_2 = charge_2
        self.voltage = voltage
        self.time = time
        self.current = current
        self.force = force
        self.velocity = velocity
        self.length = length
        self.charge_enclosed = charge_enclosed
        self.epsilon_0 = epsilon_0
        self.emf = emf
        self.electric_field1 = electric_field
        self.magnetic_field1 = magnetic_field
        self.r = resistance
        self.charge_carrier_density = charge_carrier_density
        self.thickness = thickness
        self.cross_sectional_area = cross_sectional_area
        self.mu_0 = magnetic_permeability

    def electric_force(self):
        return Coulombs_constant * self.c_1 * self.c_2 / self.distance ** 2

    def electric_field(self):
        return Coulombs_constant * self.charge / self.distance ** 2

    def electric_potential(self):
        return Coulombs_constant * self.charge / self.distance

    def capacitance(self):
        return self.charge / self.voltage

    def electric_current(self):
        return self.charge / self.time

    def resistance(self):
        return self.voltage / self.current

    def ohms_law(self):
        return self.voltage / self.current

    def coulombs_law(self):
        return self.electric_field1 * self.charge

    def gauss_law(self):
        return self.charge_enclosed / self.epsilon_0

    def faradays_law(self):
        return self.emf * self.time

    def magnetic_field(self):
        return self.force / (self.charge * self.velocity)

    def lorentz_force(self):
        return self.charge * self.velocity * self.magnetic_field1

    def electric_potential_energy(self):
        return Coulombs_constant * self.c_1 * self.c_2 / self.distance

    def work_done(self):
        return self.charge * self.voltage

    def potential_difference(self):
        return self.voltage / self.distance

    def electric_field_due_to_point_charge(self):
        return Coulombs_constant * self.c_1 / self.distance ** 2

    def electric_flux(self):
        return self.electric_field1 * self.cross_sectional_area

    def dielectric_strength(self):
        return self.voltage / self.thickness

    def electric_displacement_field(self):
        return self.electric_field1 + self.charge_carrier_density / self.epsilon_0

    def polarization_density(self):
        return self.charge_carrier_density

    def capacitance_parallel_plate(self):
        return (self.epsilon_0 * self.cross_sectional_area) / self.distance

    def capacitance_cylindrical_capacitor(self):
        return (2 * pi * self.epsilon_0 * self.length) / (math.log(self.distance / self.thickness))

    def capacitance_spherical_capacitor(self):
        return (4 * pi * self.epsilon_0 * self.distance) / (1 / self.length + 1 / (self.distance - self.thickness))

    def energy_density_electric_field(self):
        return 0.5 * self.epsilon_0 * self.electric_field1 ** 2
    
    def electric_dipole_moment(self):
        return self.charge * self.distance

    def potential_energy_of_dipole(self):
        return -self.electric_dipole_moment() * self.electric_field1

    def torque_on_dipole(self):
        return self.electric_dipole_moment() * self.electric_field1 * math.sin(math.radians(90))

    def potential_due_to_point_charge(self):
        return Coulombs_constant * self.c_1 / self.distance

    def electric_flux_through_surface(self):
        return self.electric_field1 * self.cross_sectional_area * math.cos(math.radians(180))

    def image_charge(self):
        return -self.charge

    def electric_field_inside_conductor(self):
        return 0

    def potential_of_conductor(self):
        return self.electric_field1 * self.distance

    def induced_charge_on_conductor(self):
        return -self.charge

    def electric_field_inside_dielectric(self):
        return self.electric_field1 / self.epsilon_0

    def potential_due_to_dipole(self):
        return Coulombs_constant * self.electric_dipole_moment() / self.distance ** 2

    def energy_stored_in_capacitor(self):
        return 0.5 * self.capacitance() * self.voltage ** 2
    
    def electric_field_due_to_line_charge(self):
        return (2 * Coulombs_constant * self.charge) / (self.distance * math.sqrt(math.pi))

    def potential_due_to_line_charge(self):
        return (2 * Coulombs_constant * self.charge) / math.sqrt(math.pi * self.distance)

    def electric_potential_energy_of_system(self):
        return 0.5 * Coulombs_constant * self.c_1 * self.c_2 / self.distance

    def electric_field_due_to_plane_charge(self):
        return (self.charge_carrier_density * self.epsilon_0) / 2

    def potential_due_to_plane_charge(self):
        return (self.charge_carrier_density * self.epsilon_0 * self.distance) / 2

    def energy_density_magnetic_field(self):
        return 0.5 * self.mu_0 * self.magnetic_field1 ** 2

    def inductance(self):
        return (self.mu_0 * self.cross_sectional_area) / self.length

    def self_inductance_of_coil(self):
        return (self.mu_0 * self.cross_sectional_area * self.length) / math.log(self.length / self.thickness)

    def mutual_inductance(self, other_coil):
        return (self.mu_0 * self.cross_sectional_area * other_coil.length) / math.log(other_coil.length / other_coil.thickness)
    
    def electric_field_due_to_ring_of_charge(self):
        return (Coulombs_constant * self.charge) / (2 * math.pi * self.distance)

    def potential_due_to_ring_of_charge(self):
        return (Coulombs_constant * self.charge) / (self.distance)

    def electric_field_due_to_spherical_shell(self):
        return (Coulombs_constant * self.charge) / (4 * math.pi * self.epsilon_0 * self.distance ** 2)

    def potential_due_to_spherical_shell(self):
        return (Coulombs_constant * self.charge) / (self.epsilon_0 * self.distance)

    def electric_field_due_to_non_uniformly_charged_rod(self):
        return (Coulombs_constant * self.charge_carrier_density * self.length) / (2 * self.epsilon_0)

    def potential_due_to_non_uniformly_charged_rod(self):
        return (Coulombs_constant * self.charge_carrier_density * self.length ** 2) / (4 * self.epsilon_0)

    def electric_field_due_to_disk_of_charge(self):
        return (Coulombs_constant * self.charge) / (2 * self.epsilon_0 * self.distance)

    def potential_due_to_disk_of_charge(self):
        return (Coulombs_constant * self.charge) / (self.epsilon_0 * self.distance)

    def electric_field_due_to_infinite_line_of_charge(self):
        return (Coulombs_constant * self.charge_carrier_density) / (2 * math.pi * self.distance)

    def potential_due_to_infinite_line_of_charge(self):
        return (Coulombs_constant * self.charge_carrier_density) / (2 * math.pi * self.epsilon_0 * self.distance)

    def energy_stored_in_inductor(self):
        return 0.5 * self.inductance() * self.current ** 2
    
    def electric_field_due_to_point_dipole(self):
        return (2 * Coulombs_constant * self.electric_dipole_moment()) / self.distance**3

    def potential_due_to_point_dipole(self):
        return (Coulombs_constant * self.electric_dipole_moment()) / self.distance**2

    def electric_field_inside_dielectric_sphere(self):
        return (self.charge / (4 * math.pi * self.epsilon_0 * self.distance**3)) * (3 - self.distance)

    def potential_due_to_dielectric_sphere(self):
        return (self.charge / (4 * math.pi * self.epsilon_0 * self.distance))

    def electric_field_due_to_charged_conducting_cylinder(self):
        return (self.charge / (2 * math.pi * self.epsilon_0 * self.distance))

    def potential_due_to_charged_conducting_cylinder(self):
        return (self.charge / (2 * math.pi * self.epsilon_0) * math.log(self.distance))

    def electric_field_due_to_charged_conducting_sphere(self):
        return (self.charge / (4 * math.pi * self.epsilon_0 * self.distance**2))

    def potential_due_to_charged_conducting_sphere(self):
        return (self.charge / (4 * math.pi * self.epsilon_0 * self.distance))

    def energy_stored_in_capacitor_with_dielectric(self, dielectric_constant):
        return 0.5 * dielectric_constant * self.epsilon_0 * self.capacitance() * self.voltage**2

    def energy_density_in_electric_field_with_dielectric(self, dielectric_constant):
        return 0.5 * dielectric_constant * self.epsilon_0 * self.electric_field1**2
    
    def capacitance_of_parallel_plate_capacitor_with_dielectric(self, dielectric_constant):
        # Calculates the capacitance of a parallel plate capacitor with a dielectric material between the plates.
        return (dielectric_constant * self.epsilon_0 * self.cross_sectional_area) / self.distance

    def energy_stored_in_inductor_with_magnetic_core(self, magnetic_core_permeability):
        # Considers the presence of a magnetic core within an inductor.
        return 0.5 * (self.mu_0 * magnetic_core_permeability) * self.inductance() * self.current**2
    
    def energy_stored_in_parallel_plate_capacitor(self):
        # Calculates the energy stored in a parallel plate capacitor.
        return 0.5 * self.capacitance() * self.voltage**2

    def magnetic_flux_through_loop(self):
        # Calculates the magnetic flux through a loop in a magnetic field.
        return self.magnetic_field1 * self.cross_sectional_area

    def energy_stored_in_magnetic_field_of_inductor(self):
        # Calculates the energy stored in the magnetic field of an inductor.
        return 0.5 * self.inductance() * self.magnetic_field1**2

    def energy_stored_in_magnetic_field_of_solenoide(self):
        # Calculates the energy stored in the magnetic field of an ideal solenoid.
        return 0.5 * self.inductance() * self.current**2

    