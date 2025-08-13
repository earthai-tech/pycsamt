# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0-or

"""
This module contains functions for Zonge engineering
calculations, based on the GDP DATA PROCESSING MANUAL.
"""

import numpy as np

def calculate_rho(mag_e, mag_h, asp, freq):
    """
    Calculates the Resistivity (Rho).

    Args:
        mag_e (float): E-field magnitude in µV.
        mag_h (float): H-field magnitude in pT.
        asp (float): A-spacing in meters (m).
        freq (float): Frequency in Hertz (Hz).

    Returns:
        float: The calculated resistivity (Rho) in Ωm.
    """
    # Formula: Rho = (1 / (5 * f)) * |E / H|^2
    # Where E is in V/m and H is in T.
    # The implementation below is simplified for inputs
    # in µV, pT, and m.
    term1 = 1 / (5 * freq)
    # Convert E from µV/m to V/m and H from pT to T
    e_v_per_m = (mag_e / asp) * 1e-6
    h_t = mag_h * 1e-12
    term2 = (e_v_per_m / h_t)**2
    rho = term1 * term2
    return rho

def calculate_ip(phz_e, phz_h):
    """
    Calculates the Impedance Phase (IP).

    Args:
        phz_e (float): E-field phase in mRad.
        phz_h (float): H-field phase in mRad.

    Returns:
        float: The Impedance Phase (IP) in mRad.
    """
    return phz_e - phz_h

def calculate_std_dev(values):
    """
    Calculates the standard deviation of a list of values.

    Args:
        values (list): A list of numerical values.

    Returns:
        float: The calculated standard deviation.
    """
    n = len(values)
    if n < 2:
        return 0.0

    # Formula: σ = sqrt( (A - N * B^2) / (N - 1) )
    # A = sum of values, squared
    # B = average value, squared
    a = np.sum(np.square(values))
    b = np.mean(values)
    variance = (a - n * (b**2)) / (n - 1)
    
    if variance < 0:
        return 0.0
        
    return np.sqrt(variance)

def calculate_e_field_std_dev(e_vals, asp, current):
    """
    Calculates the standard deviation for the E-field.

    Args:
        e_vals (list): E-field values in µV.
        asp (float): A-spacing in meters (m).
        current (float): Transmitter current in Amperes (a).

    Returns:
        float: Std deviation of the E-field in mV/Km*a.
    """
    # Convert E-field values to mV/Km*a
    e_conv = [(v / asp) / current for v in e_vals]
    return calculate_std_dev(e_conv)

def calculate_h_field_std_dev(h_vals, current):
    """
    Calculates the standard deviation for the H-field.

    Args:
        h_vals (list): H-field values in pT.
        current (float): Transmitter current in Amperes (a).

    Returns:
        float: Std deviation of the H-field in pT/a.
    """
    h_conv = [val / current for val in h_vals]
    return calculate_std_dev(h_conv)

def calculate_c_var(sigma, average):
    """
    Calculates the Coefficient of Variation (C-var).

    Args:
        sigma (float): Standard deviation.
        average (float): Arithmetic average.

    Returns:
        float: The Coefficient of Variation in percent.
    """
    if average == 0:
        return 0.0
    return 100 * (sigma / average)

def calculate_std_dev_rho_p(rho_values):
    """
    Calculates the Standard Deviation for Parameter RHO.

    Args:
        rho_values (list): A list of RHO values.

    Returns:
        float: The standard deviation for parameter RHO.
    """
    return calculate_std_dev(rho_values)

def calculate_std_dev_rho_c(rho_c, e_avg, h_avg,
                              sigma_e, sigma_h):
    """
    Calculates the Standard Deviation for Component RHO.

    Args:
        rho_c (float): Resistivity from averaged components.
        e_avg (float): Average E-field magnitude.
        h_avg (float): Average H-field magnitude.
        sigma_e (float): Standard deviation of E-field.
        sigma_h (float): Standard deviation of H-field.

    Returns:
        float: The standard deviation for component RHO.
    """
    if e_avg == 0 or h_avg == 0:
        return 0.0
    
    b_e = (sigma_e / e_avg)**2
    b_h = (sigma_h / h_avg)**2
    
    return rho_c * 2 * np.sqrt(b_e + b_h)

def calculate_avg_magnitude(mag_values):
    """
    Calculates the average magnitude for E or H fields.

    Args:
        mag_values (list): A list of magnitude values.

    Returns:
        float: The average magnitude.
    """
    return np.mean(mag_values)

def calculate_avg_phase(phase_values):
    """
    Calculates the average for phase values.

    Args:
        phase_values (list): A list of phase values.

    Returns:
        float: The average phase.
    """
    return np.mean(phase_values)

def calculate_parameter_avg_rho(rho_values):
    """
    Calculates the Parameter Average RHO.

    Args:
        rho_values (list): RHO values from each data block.

    Returns:
        float: The parameter average RHO.
    """
    return np.mean(rho_values)

def calculate_component_avg_rho(e_mag_avg, h_mag_avg, freq):
    """
    Calculates the Component Average RHO.

    Args:
        e_mag_avg (float): Averaged E_MAG (mV/Km*amp).
        h_mag_avg (float): Averaged H_MAG (pTesla/amp).
        freq (float): Frequency in Hz.

    Returns:
        float: The component average RHO.
    """
    # The formula uses E in mV/Km*a and H in pT/a.
    # We need to convert E to V/m and H to T to use the
    # standard resistivity formula.
    # E (V/m) = E (mV/Km*a) * 1e-6
    # H (T) = H (pT/a) * 1e-12
    # This simplifies the ratio of E/H
    e_h_ratio = (e_mag_avg / h_mag_avg) * 1e6
    return (1 / (5 * freq)) * (e_h_ratio**2)
