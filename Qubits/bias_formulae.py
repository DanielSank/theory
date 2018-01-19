import numpy as np

import labrad.units as U

"""
EC = e^2 / 2C
EJ = I_c Phi_0 / 2pi
"""

PI = np.pi
TWO_PI = 2 * np.pi
EIGHT_PI = 8 * np.pi
h = 6.62E-34 * U.J * U.s
HBAR = h / TWO_PI
ELECTRON_CHARGE = 1.6E-19 * U.C
R_K = h / (ELECTRON_CHARGE**2)
PHI0 = h / (2 * ELECTRON_CHARGE)
FIFTY_OHM = 50 * U.Ohm


def dBm_to_power(value_dBm):
    return 1*U.mW * 10**(value_dB/10.0)


def I_c_to_L(I_c):
    """Convert junction critical current to inductance."""
    return PHI0 / (TWO_PI * I_c)


def voltage_rms_at_open(power_available, Z_0):
    """RMS voltage at an open circuit, given an available power from source.

    You can also use this to convert power spectral density to voltage spectral
    density.

    Args:
        power_available (Value[W]): Available power from source.
        Z_0 (Value[Ohm]): Impedance of line.
    """
    return 2 * (power_available * Z_0)**0.5


def current_rms_at_short(power_available, Z_0):
    """RMS current at a short circuit, given an available power from source.

    You can also use this to convert power spectral density to voltage spectral
    density.

    Args:
        power_available (Value[W]): Available power from source.
        Z_0 (Value[Ohm]): Impedance of line.
    """
    return 2 * (power_available / Z_0)**0.5


def current_dc_at_short(V_source_available, attenuation, Z_0):
    return 2 * V_source_available / np.sqrt(attenuation) / Z_0


class Junction(object):
    """A Josephson junction
    
    Attributes:
            I_c (Value[A]): Critical current
            E (Value[J]): Junction energy scale
            f (Value[Hz]): Frequency scale
            L (Value[H]): Inductance
    """
    def __init__(self, I_c):
        self.I_c = I_c
        self.E = I_c * PHI0 / TWO_PI
        self.f = self.E / h
        self.L = PHI0 / (TWO_PI * I_c)


class Capacitor(object):
    """A capacitor

    Attributes:
        C (Value[F]): Capacitance
        E (Value[J]): Energy scale
        f (Value[Hz]): Frequency scale
    """
    def __init__(self, C):
        self.C = C
        self.E = ELECTRON_CHARGE**2 / (2 * C)
        self.f = self.E / h


class Transmon(object):
    """A transmon qubit
    
    Attributes:
            f_10_max (Value[Hz]): Unbiased frequency
            eta (Value[Hz]): Anharmonicity (negative!)
            capacitor (Capacitor): The qubit's capacitor
            junction (Junction): The qubit's junction
            impedance (Value[Ohm]): The qubit's characteristic impedance
            z (float): Qubit impedance divided by (R_K / 8 pi)
    """

    def __init__(self, f_10_max, f_10_min, eta):
        self.f_10_max = f_10_max
        self.eta = eta

        self.capacitor = Capacitor(-1 / (2 * R_K * eta))
        self.junction = Junction(
                TWO_PI * PHI0 * self.capacitor.C * (f_10_max - eta)**2)
        self.impedance = (self.junction.L / self.capacitor.C)**(0.5)
        self.z = self.impedance / (R_K / EIGHT_PI)

        self.f_plasma_0 = f_10_max + self.capacitor.f
        self.frustration_max = self.frustration_for_frequency(
                f_10_max, f_10_min, eta)

    @staticmethod
    def frustration_for_frequency(f_10_max, f_10, anharmonicity):
        """
        Returns (float): SQUID frustration Phi / Phi0
        """
        f_c = -anharmonicity
        return (1.0/PI) * np.arccos(((f_10 + f_c)/(f_10_max + f_c))**2)

    def Q_capacitive(self, C_c, R_ext):
        return (self.capacitor.C / C_c)**2 * (self.impedance / R_ext)

    def Q_inductive(self, M_c, R_ext):
        return (self.junction.L / M_c)**2 * (R_ext / self.impedance)

    def gamma_vs_Q(self, S_noise_power_available, Q):
        return S_noise_power_available / (Q * HBAR)

    def pi_length_vs_Q(self, signal_power_available, Q):
        return (PI/2) * (HBAR * Q / signal_power_available)**0.5

    def df_10_dPhi(self, f_10):
        """Returns df_10 / dPhi assuming symmetric junctions"""
        f = f_10 / self.f_10_max
        return -self.f_10_max / PHI0 * (PI/2) * (1.0/f) * (1 - f**4)**0.5

    def t_phi(
            self,
            f_10,
            noise_power_spectral_density_available,
            M_SQUID,
            Z_0):
        domega_10_dPhi = (2*PI) * self.df_10_dPhi(f_10)
        S_I = current_rms_at_short(
                noise_power_spectral_density_available,
                Z_0)**2
        S_flux = S_I * (M_SQUID**2)
        return 4 * (domega_10_dPhi**-2) / S_flux
