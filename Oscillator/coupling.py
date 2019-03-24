import numpy as np
import matplotlib.pyplot as plt

def normal_frequencies_rwa(omega_a, omega_b, g):
    delta = omega_a - omega_b
    omega_plus = (omega_a + omega_b)/2 + np.sqrt(g**2 + (delta / 2)**2)
    omega_minus = (omega_a + omega_b)/2 - np.sqrt(g**2 + (delta / 2)**2)
    return omega_plus, omega_minus


def normal_frequencies(omega_a, omega_b, g_plus, g_minus):
    delta = omega_a - omega_b
    s = omega_a + omega_b
    omega_plus - (1 / np.sqrt(2)) * np.sqrt(
        2 * g_minus**2 - 2 * g_plus**2 + omega_a**2 + omega_b**2 + inner_term)

    omega_minus = (1 / np.sqrt(2)) * np.sqrt(
        2 * g_minus**2 - 2 * g_plus**2 + omega_a**2 + omega_b**2 - inner_term)

    return omega_plus, omega_minus


def omega_minus_rwa(omega_a, omega_b, g):
    delta = omega_a - omega_b
    return (omega_a + omega_b)/2 - np.sqrt(g**2 + (delta / 2)**2)


def plot_rwa(omega_a, omega_b, gs, plot_keywords):
    _, ax = plt.subplots()
    for g, linestyle in zip(gs, ['--', '-']):
        omega_plus, omega_minus = normal_frequencies_rwa(omega_a, omega_b, g)
        ax.plot(
            omega_a,
            omega_plus,
            color='r',
            linestyle=linestyle,
            label=r'$\omega$ +, g={}'.format(g),
            **plot_keywords,)
        ax.plot(
            omega_a,
            omega_minus,
            color='b',
            linestyle=linestyle,
            label=r'$\omega -$, g={}'.format(g),
            **plot_keywords,)
    return ax


def go_rwa(plot_keywords, legend=True):
    ax = plot(
        np.linspace(9.5, 10.5),
        10,
        [0.02, 0.1],
        plot_keywords,)
    plt.grid()
    plt.xlabel(r"$\omega_b'$")
    plt.ylabel(r'$\omega_\pm$')
    if legend:
        ax.legend(loc='upper left')
    return ax
