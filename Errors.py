from .constants import *
import numpy as np

class Errors:

    def __init__(self):
        pass

    @staticmethod
    def relative_error(a, b):
        """
        Calculate the relative error between two values.
        Formula: abs(a - b) / max(abs(a), abs(b))
        """
        try:
            return abs(a - b) / max(abs(a), abs(b))
        except ZeroDivisionError:
            return 0  # Avoid division by zero

    @staticmethod
    def error_muldiv(a, b, c, d):
        """
        Calculate the error for a multiplication and division operation.
        Formula: ((a / c) + (b / d))
        """
        try:
            return (a / c) + (b / d)
        except ZeroDivisionError:
            return ValueError("Division by zero encountered")

    @staticmethod
    def error_addsub(a, b):
        """
        Calculate the error for an addition and subtraction operation.
        Formula: a + b
        """
        return a + b

    @staticmethod
    def percentage_error(observed, expected):
        """
        Calculate the percentage error between observed and expected values.
        Formula: ((observed - expected) / expected) * 100
        """
        try:
            return ((observed - expected) / expected) * 100
        except ZeroDivisionError:
            return ValueError("Division by zero encountered")

    @staticmethod
    def absolute_error(data):
        """
        Calculate the absolute error for a list of values.
        Formula: abs(sum(data) / len(data))
        """
        try:
            return abs(np.sum(data) / len(data))
        except ZeroDivisionError:
            return ValueError("Division by zero encountered")

    @staticmethod
    def mean_absolute_error(data):
        """
        Calculate the mean absolute error for a list of values.
        Formula: sum(abs(data)) / len(data)
        """
        try:
            return np.sum(np.abs(data)) / len(data)
        except ZeroDivisionError:
            return ValueError("Division by zero encountered")

    # Additional error-related functions can be added here for a more comprehensive error handling system.

    @staticmethod
    def weighted_mean_absolute_error(data, weights):
        """
        Calculate the weighted mean absolute error for a list of values and corresponding weights.
        Formula: sum(weights * abs(data)) / sum(weights)
        """
        try:
            return np.sum(weights * np.abs(data)) / np.sum(weights)
        except ZeroDivisionError:
            return ValueError("Division by zero encountered")

    @staticmethod
    def rms_error(data):
        """
        Calculate the root mean square error for a list of values.
        Formula: sqrt(sum(data**2) / len(data))
        """
        try:
            return np.sqrt(np.sum(np.square(data)) / len(data))
        except ZeroDivisionError:
            return ValueError("Division by zero encountered")