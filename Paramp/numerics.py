import matplotlib.pyplot as plt
import numpy as np


PI = np.pi


def L(omega, Q, omega_0):
    gamma = omega_0 / (2 * Q)
    return omega_0**2 - omega**2 + 1.0j*2*omega*gamma


def alpha(omega, omega_p, omega_0, Q, A):
    omega_i = omega_p - omega
    numerator = np.conjugate(L(omega_i, Q, omega_0))
    denominator = L(omega, Q, omega_0) * np.conjugate(L(omega_i, Q, omega_0)) - (A * omega_0**2 / 2)**2
    return numerator/denominator


def plot():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    omega_p = 2 * PI * 8.1
    omega_0 = 2 * PI * 4
    Q = 5.6

    omegas = np.linspace(2*PI*3.25, 2*PI*5, 100)


    for A in [0, 0.1, 0.2, 0.3, 0.4, 0.7]:
        ax.semilogy(
            omegas/(2*PI),
            np.abs(alpha(omegas, omega_p, omega_0, Q, A)),
            label='A = {}'.format(A),)
        ax.grid(True)
        ax.legend(loc='upper right')
        ax.set_xlabel('Frequency')
        ax.set_ylabel('|alpha| / J')
