import numpy as np
import matplotlib.pyplot as plt
pi = np.pi


def square_window_k(k):
    return np.sin(pi * k) / (pi * k)


def leakage_data(xi, k_max, window=square_window_k):

    ks = np.linspace(-k_max, k_max, 2*k_max + 1)
    k_continuous = np.linspace(-k_max, k_max, 1000)

    fourier_coefficients = window(ks)
    fourier_transform = window(k_continuous)
    #commensurate frequency

    fourier_coefficients_shifted = window(ks - xi)
    fourier_transform_shifted = window(k_continuous - xi)
    # non-commensurate frequency

    return (ks,
            k_continuous,
            (fourier_coefficients, fourier_transform),
            (fourier_coefficients_shifted, fourier_transform_shifted))


def plot(xi, k_max):
    marker_size = 18

    ks, k_continuous, unshifted_data, shifted_data = leakage_data(xi, k_max)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(k_continuous, shifted_data[1], 'b-', linewidth=3,
            label=r'$\xi={}$'.format(xi))
    ax.plot(ks, shifted_data[0], 'b.', markersize=marker_size)
    ax.plot(k_continuous, unshifted_data[1], 'k--', linewidth=3, alpha=0.7,
            label=r'$\xi={}$'.format(0))
    ax.plot(ks, unshifted_data[0], 'k.', markersize=marker_size, alpha=0.7)
    ax.plot(0, 1, 'k.', markersize=marker_size, alpha=0.7)

    ax.grid()
    ax.legend()


def plot_squared(xi, k_max):
    marker_size = 16

    ks, k_continuous, unshifted_data, shifted_data = leakage_data(xi, k_max)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.semilogy(k_continuous, np.abs(shifted_data[1])**2, 'b-', linewidth=2)
    ax.semilogy(ks, np.abs(shifted_data[0])**2, 'b.', markersize=marker_size)
    ax.semilogy(k_continuous, np.abs(unshifted_data[1])**2, 'k--',
            linewidth=2, alpha=0.7)
    ax.semilogy(ks, np.abs(unshifted_data[0])**2, 'k.',
            markersize=marker_size, alpha=0.7)

    plt.grid()
