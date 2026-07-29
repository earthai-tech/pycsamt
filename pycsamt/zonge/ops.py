# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

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
    term2 = (e_v_per_m / h_t) ** 2
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


def calculate_std_dev_rho_c(rho_c, e_avg, h_avg, sigma_e, sigma_h):
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

    b_e = (sigma_e / e_avg) ** 2
    b_h = (sigma_h / h_avg) ** 2

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


def calculate_magnetic_induction(h_mag, rho):
    """
    Calculates the Magnetic Induction (M) from the magnetic field
    amplitude and resistivity.

    Args:
        h_mag (float): Magnetic field magnitude in nT.
        rho (float): Resistivity in Ωm.

    Returns:
        float: The calculated magnetic induction (M) in nT·m.

    .. math::
        M = \frac{H}{\rho}
        \text{Where H is the magnetic field magnitude in nT and }
        \rho \text{ is the resistivity in Ωm.}
    """
    return h_mag / rho


def calculate_apparent_resistivity(e_mag, h_mag, geometric_factor=1.0):
    r"""
    Calculates the apparent resistivity (Rho Apparent) from the E-field
    and H-field magnitudes using the geometric factor.

    Args:
        e_mag (float): Electric field magnitude in V/m.
        h_mag (float): Magnetic field magnitude in T.
        geometric_factor (float): Geometric factor, typically 1.0 for
            vertical dipole configuration (default is 1.0).

    Returns:
        float: The calculated apparent resistivity (Rho Apparent) in Ωm.

    .. math::
        \rho_a = \frac{5 \cdot E}{H} \cdot \text{Geometric Factor}
        \text{Where E is in V/m and H is in T.}
    """
    return (5 * e_mag / h_mag) * geometric_factor


def calculate_snr(signal_values, noise_values):
    r"""
    Calculates the Signal-to-Noise Ratio (SNR) for the given signal and
    noise values.

    Args:
        signal_values (list): A list of signal values.
        noise_values (list): A list of noise values.

    Returns:
        float: The signal-to-noise ratio (SNR).

    .. math::
        \text{SNR} = \frac{\text{Signal Mean}}{\text{Noise Std Dev}}
        \text{Where Signal Mean is the average signal value, and}
        \text{Noise Std Dev is the standard deviation of the noise values.}
    """
    signal_mean = np.mean(signal_values)
    noise_std_dev = calculate_std_dev(noise_values)
    if noise_std_dev == 0:
        return np.inf  # Avoid division by zero
    return signal_mean / noise_std_dev


def calculate_phase_error(phz_e, phz_h):
    r"""
    Calculates the phase error based on the difference between the E-field
    and H-field phases.

    Args:
        phz_e (float): E-field phase in degrees or radians.
        phz_h (float): H-field phase in degrees or radians.

    Returns:
        float: The calculated phase error in degrees or radians.

    .. math::
        \text{Phase Error} = \left| \phi_E - \phi_H \right|
        \text{Where }\phi_E\text{ is the E-field phase and }\\
            \phi_H\text{ is the H-field phase.}
    """
    return np.abs(phz_e - phz_h)


def propagate_resistivity_error(rho, e_avg, h_avg, sigma_e, sigma_h):
    r"""
    Propagates the error in resistivity based on the standard deviations
    of E-field and H-field.

    Args:
        rho (float): The calculated resistivity (Rho).
        e_avg (float): Average E-field magnitude.
        h_avg (float): Average H-field magnitude.
        sigma_e (float): Standard deviation of the E-field.
        sigma_h (float): Standard deviation of the H-field.

    Returns:
        float: The propagated error in resistivity.

    .. math::
        \text{Error in Resistivity} = \rho \cdot \sqrt{
            \left( \frac{\sigma_E}{E_{\text{avg}}} \right)^2 +
            \left( \frac{\sigma_H}{H_{\text{avg}}} \right)^2 }
        \text{Where }\sigma_E\text{ and }\sigma_H\\\
            text{ are the standard deviations of E-field and H-field.}
    """
    e_h_ratio = (e_avg / h_avg) * 1e6  # noqa
    error_in_rho = rho * np.sqrt((sigma_e / e_avg) ** 2 + (sigma_h / h_avg) ** 2)
    return error_in_rho


def calculate_avg_amplitude(field_values):
    r"""
    Calculates the average amplitude for E-field or H-field values.

    Args:
        field_values (list): A list of E-field or H-field values.

    Returns:
        float: The average amplitude.

    .. math::
        \text{Avg Amplitude} = \frac{1}{N} \sum_{i=1}^N |x_i|
        \text{Where } x_i\text{ represents the individual field values.}
    """
    return np.mean(np.abs(field_values))


def calculate_relative_error(rho, sigma_rho):
    """
    Calculates the relative error in resistivity.

    Args:
        rho (float): The calculated resistivity.
        sigma_rho (float): The standard deviation of resistivity.

    Returns:
        float: The relative error in resistivity as a percentage.

    .. math::
        \text{Relative Error} = \frac{\\sigma_{\rho}}{\rho} \times 100
        \text{Where }\\sigma_{\rho}\\\
            text{ is the standard deviation of resistivity.}
    """
    if rho == 0:
        return 0.0
    return (sigma_rho / rho) * 100


def calculate_magnitude_ratio(e_mag, h_mag):
    """
    Calculates the magnitude ratio between the E-field and H-field.

    Args:
        e_mag (float): Electric field magnitude in V/m.
        h_mag (float): Magnetic field magnitude in T.

    Returns:
        float: The magnitude ratio (E/H).

    .. math::
        \text{Magnitude Ratio} = \frac{E_{\text{mag}}}{H_{\text{mag}}}
        \text{Where } E_{\text{mag}} \text{ is the electric field magnitude and }
        H_{\text{mag}} \text{ is the magnetic field magnitude.}
    """
    return e_mag / h_mag


def calculate_resistivity_phase(rho, phase_e, phase_h):
    """
    Calculates the resistivity phase based on resistivity and phase
    differences between E-field and H-field.

    Args:
        rho (float): The resistivity value in Ωm.
        phase_e (float): E-field phase in radians.
        phase_h (float): H-field phase in radians.

    Returns:
        float: The resistivity phase in radians.

    .. math::
        \text{Resistivity Phase} = \text{atan2}(E_{\text{phase}} - H_{\text{phase}})
        \text{Where E}_{\text{phase}}\text{ and H}_\\
            {\text{phase}}\text{ are the phases of E-field and H-field.}
    """
    return np.arctan2(phase_e - phase_h, rho)


def calculate_frequency_dependent_resistivity(e_mag, h_mag, freq):
    """
    Calculates frequency-dependent resistivity (Rho) based on the E-field
    and H-field magnitudes.

    Args:
        e_mag (float): Electric field magnitude in V/m.
        h_mag (float): Magnetic field magnitude in T.
        freq (float): Frequency in Hz.

    Returns:
        float: The calculated frequency-dependent resistivity (Rho) in Ωm.

    .. math::
        \rho_{\text{freq}} = \frac{E_{\text{mag}}}\\
            {H_{\text{mag}}} \\cdot \frac{1}{f}
        \text{Where } f \text{ is the frequency in Hz.}
    """
    return (e_mag / h_mag) / freq


def calculate_rho_correction(rho, e_std, h_std, e_avg, h_avg):
    """
    Corrects resistivity values based on the standard deviations of E-field
    and H-field, and their average values.

    Args:
        rho (float): Resistivity in Ωm.
        e_std (float): Standard deviation of E-field.
        h_std (float): Standard deviation of H-field.
        e_avg (float): Average E-field magnitude.
        h_avg (float): Average H-field magnitude.

    Returns:
        float: The corrected resistivity value.

    .. math::
        \rho_{\text{corr}} = \rho \\cdot \\left\\
            ( 1 + \frac{\\sigma_E}{E_{\text{avg}}} +\\
             \frac{\\sigma_H}{H_{\text{avg}}} \right)
        \text{Where } \\sigma_E \text{ and }\\
            \\sigma_H \text{ are the standard deviations of E-field and H-field.}
    """
    return rho * (1 + (e_std / e_avg) + (h_std / h_avg))


def calculate_averaged_magnitude(values):
    """
    Calculates the average of magnitudes, useful for processing both
    E-field and H-field magnitudes.

    Args:
        values (list): A list of magnitude values.

    Returns:
        float: The average magnitude.

    .. math::
        \text{Avg Magnitude} = \frac{1}{N} \\sum_{i=1}^N |x_i|
        \text{Where } x_i\text{ represents individual field values.}
    """
    return np.mean(np.abs(values))


def calculate_conductivity(rho):
    """
    Calculates the conductivity (Sigma) from the resistivity.

    Args:
        rho (float): Resistivity in Ωm.

    Returns:
        float: The calculated conductivity (Sigma) in S/m.

    .. math::
        \\sigma = \frac{1}{\rho}
        \text{Where } \rho \text{ is the resistivity in Ωm.}
    """
    return 1 / rho


def calculate_error_propagation_amplitude(e_std, h_std, rho_std, e_avg, h_avg, rho):
    """
    Propagates error for the amplitude based on the standard deviations
    of the E-field, H-field, and resistivity.

    Args:
        e_std (float): Standard deviation of the E-field.
        h_std (float): Standard deviation of the H-field.
        rho_std (float): Standard deviation of resistivity.
        e_avg (float): Average value of the E-field.
        h_avg (float): Average value of the H-field.
        rho (float): Resistivity value.

    Returns:
        float: The propagated error for amplitude.

    .. math::
        \text{Error in Amplitude} = \\sqrt{ \\left( \frac{\\sigma_E}\\
                                                 {E_{\text{avg}}} \right)^2
        + \\left( \frac{\\sigma_H}{H_{\text{avg}}}\\
                \right)^2 + \\left( \frac{\\sigma_{\rho}}{\rho} \right)^2 }
        \text{Where }\\sigma_E\text{, }\\sigma_H\text{, and }\\
        \\sigma_{\rho}\text{ are the standard deviations\\
                           for E-field, H-field, and resistivity.}
    """
    return np.sqrt((e_std / e_avg) ** 2 + (h_std / h_avg) ** 2 + (rho_std / rho) ** 2)


def calculate_e_field_error(e_vals, asp, current):
    """
    Calculates the error for the E-field based on field values, A-spacing,
    and current.

    Args:
        e_vals (list): E-field values in µV.
        asp (float): A-spacing in meters (m).
        current (float): Transmitter current in Amperes (a).

    Returns:
        float: The calculated error for the E-field.

    .. math::
        \text{E-field Error} = \frac{1}{N} \\sum_{i=1}^N\\
            \\left| E_{\text{val}} - \frac{E_{\text{avg}}}{\text{current}} \right|
        \text{Where } E_{\text{val}} \text{ are the individual E-field values and }
        E_{\text{avg}} \text{ is the average E-field.}
    """
    e_conv = [(v / asp) / current for v in e_vals]
    e_avg = np.mean(e_conv)
    return np.mean(np.abs(np.array(e_vals) - e_avg))
