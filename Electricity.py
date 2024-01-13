from .constants import Coulombs_constant, plancks_constant
import math
import numpy as np

class Electricity:
    """
        This class is the Electricity Suit. it contains the functions
        functions for Electricity calculations.
    """
    def __init__(self, voltage=1, current=1, resistance=1,
                 q1=1, q2=1, r=1) -> None:
        self.V = voltage
        self.C = current
        self.R = resistance
        self.k = Coulombs_constant
        self.q1 = q1
        self.q2 = q2
        self.r = r

    def force_electrostatics(self) -> float:
        return (self.k * self.q1 * self.q2) / (self.r ** 2)

    def resistance(self) -> float:
        if self.C != 0 and self.C is not None:
            return self.V / self.C
        else:
            # Handle division by zero
            raise ValueError("Current cannot be zero for resistance calculation.")

    def current(self) -> float:
        if self.R != 0 and self.R is not None:
            return self.V / self.R
        else:
            # Handle division by zero
            raise ValueError("Resistance cannot be zero for current calculation.")

    def voltage(self) -> float:
        return self.C * self.R

    def power(self) -> float:
        return self.V * self.C

    def ohms_law(self):
        if self.R != 0 and self.R is not None:
            return self.V / self.R
        else:
            # Handle division by zero
            raise ValueError("Resistance cannot be zero for ohms law calculation.")

    def capacitance(self) -> float:
        """
        Calculate the capacitance of a capacitor.
        Formula: Capacitance = q / V
        """
        if self.V != 0 and self.V is not None:
            return self.q1 / self.V
        else:
            # Handle division by zero
            raise ValueError("Voltage cannot be zero for capacitance calculation.")

    def energy_stored_in_capacitor(self) -> float:
        """
        Calculate the energy stored in a capacitor.
        Formula: Energy Stored = 0.5 * q * V
        """
        return 0.5 * self.q1 * self.V

    def inductance(self) -> float:
        """
        Calculate the inductance of an inductor.
        Formula: Inductance = (V * dt) / di
        """
        if self.C != 0 and self.C is not None:
            return (self.V * self.r) / self.C
        else:
            # Handle division by zero
            raise ValueError("Current cannot be zero for inductance calculation.")

    def magnetic_energy_stored_in_inductor(self) -> float:
        """
        Calculate the magnetic energy stored in an inductor.
        Formula: Magnetic Energy Stored = 0.5 * L * I^2
        """
        return 0.5 * self.inductance() * self.C**2

    def skin_effect_depth(self, frequency, conductivity, permeability) -> float:
        """
        Calculate the skin effect depth in a conductor.
        Formula: Skin Effect Depth = sqrt(2 / (pi * frequency * conductivity * permeability))
        """
        return math.sqrt(2 / (math.pi * frequency * conductivity * permeability))

    def power_factor(self) -> float:
        """
        Calculate the power factor in an AC circuit.
        Formula: Power Factor = cos(theta), where theta is the phase angle between voltage and current
        """
        if self.V != 0 and self.V is not None:
            return math.cos(math.atan(self.C / self.V))
        else:
            # Handle division by zero
            raise ValueError("Voltage cannot be zero for power factor calculation.")

    def complex_power(self) -> complex:
        """
        Calculate the complex power in an AC circuit.
        Formula: Complex Power = V * conj(I)
        """
        return self.V * complex(self.C, -self.R)
    
    def energy_transmission_loss(self, distance, resistivity, cross_sectional_area) -> float:
        """
        Calculate energy transmission loss in a conductor.
        Formula: Transmission Loss = (resistivity * distance * self.C**2) / cross_sectional_area
        """
        return (resistivity * distance * self.C**2) / cross_sectional_area

    def impedance_matching(self, source_impedance, load_impedance) -> float:
        """
        Calculate the reflection coefficient for impedance matching.
        Formula: Reflection Coefficient = (load_impedance - source_impedance) / (load_impedance + source_impedance)
        """
        return (load_impedance - source_impedance) / (load_impedance + source_impedance)

    def skin_depth(self, frequency, permeability, conductivity) -> float:
        """
        Calculate the skin depth in a conductor.
        Formula: Skin Depth = sqrt(2 / (pi * frequency * permeability * conductivity))
        """
        return np.sqrt(2 / (np.pi * frequency * permeability * conductivity))

    def power_factor_correction(self, original_power_factor, desired_power_factor) -> float:
        """
        Calculate the required reactive power for power factor correction.
        Formula: Required Reactive Power = P * tan(arccos(desired_power_factor)) - P * tan(arccos(original_power_factor))
        """
        original_power = self.power()
        return original_power * (np.tan(np.arccos(desired_power_factor)) - np.tan(np.arccos(original_power_factor)))

    def thomson_effect(self, temperature_difference, current_density) -> float:
        """
        Calculate the Thomson effect in a conductor.
        Formula: Thomson Effect = (self.C * temperature_difference) / (self.C * current_density)
        """
        if self.C != 0 and self.C is not None and current_density != 0 and current_density is not None:
            return (self.C * temperature_difference) / (self.C * current_density)
        else:
            # Handle division by zero
            raise ValueError("Current density and self capacitance cannot be zero for Thomson effect calculation.")
        
    def magnetic_cooling_power(self, magnetic_field_strength, temperature_difference, thermal_conductivity) -> float:
        """
        Calculate the magnetic cooling power in a material.
        Formula: Magnetic Cooling Power = (self.C * magnetic_field_strength * temperature_difference) / thermal_conductivity
        """
        if self.C != 0 and self.C is not None and thermal_conductivity != 0 and thermal_conductivity is not None:
            return (self.C * magnetic_field_strength * temperature_difference) / thermal_conductivity
        else:
            # Handle division by zero
            raise ValueError("Thermal conductivity and self capacitance cannot be zero for magnetic cooling power calculation.")

    def drift_velocity_derivation(self, electric_field, charge_density, mobility) -> float:
        """
        Derive the drift velocity equation in terms of electric field, charge density, and mobility.
        Formula: Drift Velocity = electric_field / (charge_density * mobility)
        """
        if charge_density != 0 and charge_density is not None and mobility != 0 and mobility is not None:
            return electric_field / (charge_density * mobility)
        else:
            # Handle division by zero
            raise ValueError("Charge density and mobility cannot be zero for drift velocity derivation.")

    def skin_effect_loss_derivative(self, frequency, conductivity, voltage) -> float:
        """
        Derive the skin effect loss equation in terms of frequency, conductivity, and voltage.
        Formula: Skin Effect Loss Derivative = (sqrt(pi * frequency * conductivity) * voltage) / (2 * sqrt(2))
        """
        return (np.sqrt(np.pi * frequency * conductivity) * voltage) / (2 * np.sqrt(2))

    def magnetic_flux_density_derivative(self, magnetic_flux, area) -> float:
        """
        Derive the magnetic flux density equation in terms of magnetic flux and area.
        Formula: Magnetic Flux Density Derivative = magnetic_flux / area
        """
        if area != 0 and area is not None:
            return magnetic_flux / area
        else:
            # Handle division by zero
            raise ValueError("Area cannot be zero for magnetic flux density derivation.")
    
    def thermoelectric_power(self, temperature_difference, Seebeck_coefficient, electrical_conductivity) -> float:
        """
        Calculate the thermoelectric power (Seebeck effect) in a material.
        Formula: Thermoelectric Power = Seebeck_coefficient * temperature_difference / electrical_conductivity
        """
        if electrical_conductivity != 0 and electrical_conductivity is not None:
            return Seebeck_coefficient * temperature_difference / electrical_conductivity
        else:
            # Handle division by zero
            raise ValueError("Electrical conductivity cannot be zero for thermoelectric power calculation.")

    def hall_effect_voltage(self, magnetic_field, current, charge_carrier_density) -> float:
        """
        Calculate the Hall effect voltage in a conductor.
        Formula: Hall Effect Voltage = magnetic_field * current / (charge_carrier_density * self.C)
        """
        if charge_carrier_density != 0 and charge_carrier_density is not None:
            return (magnetic_field * current) / (charge_carrier_density * self.C)
        else:
            # Handle division by zero
            raise ValueError("Charge carrier density cannot be zero for Hall effect voltage calculation.")
    
    def poynting_vector(self, electric_field, magnetic_field) -> np.ndarray:
        """
        Calculate the Poynting vector in an electromagnetic field.
        Formula: Poynting Vector = electric_field x magnetic_field
        """
        return np.cross(electric_field, magnetic_field)

    def transmission_line_impedance(self, capacitance_per_unit_length, inductance_per_unit_length) -> float:
        """
        Calculate the characteristic impedance of a transmission line.
        Formula: Characteristic Impedance = sqrt(inductance_per_unit_length / capacitance_per_unit_length)
        """
        if capacitance_per_unit_length != 0 and capacitance_per_unit_length is not None:
            return np.sqrt(inductance_per_unit_length / capacitance_per_unit_length)
        else:
            # Handle division by zero
            raise ValueError("Capacitance per unit length cannot be zero for characteristic impedance calculation.")

    def reflection_coefficient(self, load_impedance, transmission_line_impedance) -> float:
        """
        Calculate the reflection coefficient in a transmission line.
        Formula: Reflection Coefficient = (load_impedance - transmission_line_impedance) / (load_impedance + transmission_line_impedance)
        """
        return (load_impedance - transmission_line_impedance) / (load_impedance + transmission_line_impedance)

    def cyclotron_frequency(self, magnetic_field_strength, charge) -> float:
        """
        Calculate the cyclotron frequency for charged particles in a magnetic field.
        Formula: Cyclotron Frequency = (charge * magnetic_field_strength) / (2 * pi * self.mass)
        """
        if charge != 0 and charge is not None:
            return (charge * magnetic_field_strength) / (2 * np.pi * self.mass)
        else:
            # Handle division by zero
            raise ValueError("Charge cannot be zero for cyclotron frequency calculation.")

    def electron_mobility(self, electrical_conductivity, charge_density, electric_field) -> float:
        """
        Calculate the electron mobility in a material.
        Formula: Electron Mobility = electrical_conductivity / (charge_density * electric_field)
        """
        if charge_density != 0 and charge_density is not None and electric_field != 0 and electric_field is not None:
            return electrical_conductivity / (charge_density * electric_field)
        else:
            # Handle division by zero
            raise ValueError("Charge density and electric field cannot be zero for electron mobility calculation.")

    def electrostriction_constant(self, applied_voltage, electric_field, strain) -> float:
        """
        Calculate the electrostriction constant in a material.
        Formula: Electrostriction Constant = (applied_voltage / (electric_field * strain)) * (1 / self.C)
        """
        if electric_field != 0 and electric_field is not None and strain != 0 and strain is not None:
            return (applied_voltage / (electric_field * strain)) * (1 / self.C)
        else:
            # Handle division by zero
            raise ValueError("Electric field and strain cannot be zero for electrostriction constant calculation.")

    def electron_heat_capacity(self, temperature) -> float:
        """
        Calculate the electron heat capacity in a material.
        Formula: Electron Heat Capacity = self.C * temperature
        """
        return self.C * temperature

    def mutual_inductance_coupling_coefficient(self, mutual_inductance, self_inductance_1, self_inductance_2) -> float:
        """
        Calculate the mutual inductance coupling coefficient between two coils.
        Formula: Coupling Coefficient = mutual_inductance / sqrt(self_inductance_1 * self_inductance_2)
        """
        if self_inductance_1 != 0 and self_inductance_1 is not None and self_inductance_2 != 0 and self_inductance_2 is not None:
            return mutual_inductance / np.sqrt(self_inductance_1 * self_inductance_2)
        else:
            # Handle division by zero
            raise ValueError("Self-inductance values cannot be zero for coupling coefficient calculation.")

    def skin_effect_depth(self, frequency, conductivity, permeability) -> float:
        """
        Calculate the skin effect depth in a conductor with magnetic permeability.
        Formula: Skin Effect Depth = sqrt(2 / (pi * frequency * conductivity * permeability))
        """
        return np.sqrt(2 / (np.pi * frequency * conductivity * permeability))

    def quantum_tunneling_probability(self, barrier_width, particle_mass, energy) -> float:
        """
        Calculate the probability of quantum tunneling through a potential barrier.
        Formula: Tunneling Probability = exp(-2 * barrier_width * sqrt(2 * particle_mass * energy) / h_bar)
        """
        h_bar = plancks_constant / (2 * np.pi)
        return np.exp(-2 * barrier_width * np.sqrt(2 * particle_mass * energy) / h_bar)

    