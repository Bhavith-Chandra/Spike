import math


class Charge:
    electron = 1.602 * (math.pow(10, -19))
    elementary_charge = 1.602176634e-19  # Elementary charge in Coulombs

    # Common ions and their charges
    sodium_ion = 1 * elementary_charge  # Na+
    magnesium_ion = 2 * elementary_charge  # Mg2+
    aluminum_ion = 3 * elementary_charge  # Al3+
    chloride_ion = -1 * elementary_charge  # Cl-
    sulfate_ion = -2 * elementary_charge  # SO4^2-
    ammonium_ion = 1 * elementary_charge  # NH4+

    # Exotic particles and their charges
    muon = -1 * elementary_charge  # μ-
    positron = 1 * elementary_charge  # e+

    def __init__(self, value):
        """
        Initialize a Charge object with a custom charge value.

        Parameters:
            value (float): The charge value in Coulombs.
        """
        self.value = value

    def is_positive(self):
        """
        Check if the charge is positive.

        Returns:
            bool: True if the charge is positive, False otherwise.
        """
        return self.value > 0

    def is_negative(self):
        """
        Check if the charge is negative.

        Returns:
            bool: True if the charge is negative, False otherwise.
        """
        return self.value < 0

    def __repr__(self):
        """
        Return a string representation of the Charge object.

        Returns:
            str: String representation.
        """
        return f"Charge({self.value})"
    
        elementary_charge = 1.602176634e-19  # Elementary charge in Coulombs

    # Common ions and their charges
    sodium_ion = 1 * elementary_charge  # Na+
    magnesium_ion = 2 * elementary_charge  # Mg2+
    aluminum_ion = 3 * elementary_charge  # Al3+
    chloride_ion = -1 * elementary_charge  # Cl-
    sulfate_ion = -2 * elementary_charge  # SO4^2-
    ammonium_ion = 1 * elementary_charge  # NH4+

    # Exotic particles and their charges
    muon = -1 * elementary_charge  # μ-
    positron = 1 * elementary_charge  # e+

    def __init__(self, value):
        """
        Initialize a Charge object with a custom charge value.

        Parameters:
            value (float): The charge value in Coulombs.
        """
        self.value = value

    def is_positive(self):
        """
        Check if the charge is positive.

        Returns:
            bool: True if the charge is positive, False otherwise.
        """
        return self.value > 0

    def is_negative(self):
        """
        Check if the charge is negative.

        Returns:
            bool: True if the charge is negative, False otherwise.
        """
        return self.value < 0

    def calculate_force(self, other_charge, distance):
        """
        Calculate the electrostatic force between two charges.

        Parameters:
            other_charge (Charge): Another Charge object.
            distance (float): Distance between the charges.

        Returns:
            float: Electrostatic force between the charges.
        """
        k = 8.9875e9  # Coulomb's constant
        return k * abs(self.value * other_charge.value) / (distance ** 2)

    def __repr__(self):
        """
        Return a string representation of the Charge object.

        Returns:
            str: String representation.
        """
        return f"Charge({self.value})"
    
    def is_positive(self):
        """
        Check if the charge is positive.

        Returns:
            bool: True if the charge is positive, False otherwise.
        """
        return self.value > 0

    def is_negative(self):
        """
        Check if the charge is negative.

        Returns:
            bool: True if the charge is negative, False otherwise.
        """
        return self.value < 0

    def calculate_force(self, other_charge, distance):
        """
        Calculate the electrostatic force between two charges.

        Parameters:
            other_charge (Charge): Another Charge object.
            distance (float): Distance between the charges.

        Returns:
            float: Electrostatic force between the charges.
        """
        k = 8.9875e9  # Coulomb's constant
        return k * abs(self.value * other_charge.value) / (distance ** 2)

    def electric_field_strength(self, distance):
        """
        Calculate the electric field strength at a certain distance from the charge.

        Parameters:
            distance (float): Distance from the charge.

        Returns:
            float: Electric field strength.
        """
        k = 8.9875e9  # Coulomb's constant
        return k * abs(self.value) / (distance ** 2)

    def potential_energy(self, other_charge, distance):
        """
        Calculate the potential energy between two charges.

        Parameters:
            other_charge (Charge): Another Charge object.
            distance (float): Distance between the charges.

        Returns:
            float: Potential energy between the charges.
        """
        k = 8.9875e9  # Coulomb's constant
        return k * (self.value * other_charge.value) / distance

    def __repr__(self):
        """
        Return a string representation of the Charge object.

        Returns:
            str: String representation.
        """
        return {self.value}