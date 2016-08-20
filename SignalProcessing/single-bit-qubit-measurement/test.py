from __future__ import division

import numpy as np
import scipy.special as ss
import scipy.stats as scipystats
import matplotlib.pyplot as plt
import pyle.analysis.signal_processing as pasp
#import pyle.analysis.fluxqubit.numerics as pafn


PI = np.pi


def white_noise_modulus_pdf(r, sigma, n):
    """PDF of modulus of a white noise DFT coefficient.

    Note that there is no assumption about the distribution of the noise in the
    time domain, so long as the samples are uncorrelated.

    Args:
        r (float): Modulus at which to compute the PDF.
        sigma (float): Standard deviation of the time domain signal from whence
            came the DFT coefficients.
        n (int): Number of points in the time domain signal.

    Returns:
        (float): Probability density at the provided modulus value.
    """
    q = sigma**2 / n
    return (2 / q) * r * np.exp(- r**2 / q)


def white_noise_modulus_cdf(r, sigma, n):
    """Integral of previous function."""
    q = sigma**2 / n
    return 1 - np.exp(-r**2 / q)


def white_noise_modulus_squared_pdf(r_squared, sigma, n):
    """Probability density of the mod square of a white noise DFT coefficient.

    Args:
        r_squared (float): Mod square at which to compute the PDF.
        sigma (float): Standard deviation of the time domain signal from whence
            came the DFT coefficients.
        n (int): Number of points in the time domain signal

    Returns:
        (float): Probability density at the provided mod square value.
    """
    q = sigma**2 / n
    return (1 / q) * np.exp(-r_squared / q)


def white_noise_product_modulus_pdf(r, sigma, n):
    """PDF of modulus of product of two white noise DFT coefficients.

    Suppose we have two white noise time sequences each with time domain
    standard deviation equal to sigma:
        a = [b_0, a_1, ..., a_{n-1}] and
        b = [a_0, b_1, ..., b_{n-1}] .

    Now we form the discrete Fourier transforms of each sequence
        A = [A_0, A_k, ..., A_n]
        B = [B_0, B_k, ..., B_n]
    where
        A_k = (1/n) \sum_{m=0}^{n-1} a_m \exp(-i 2 \pi m k / n) .

    Now we multiply the DFT's together to form
        C = [A_0 B_0^*, A_1 B_1^*, ..., A_{n-1} B_{n-1}^*]

    This function tells you the probability distribution P(r) where r = |C_k|.
    Note that since we're dealing with white noise the distribution is
    independent of k.

    Args:
        r (float): Modulus at which you want the probability density.
        sigma (float): Standard deviation of the time domain signal from whence
            came the Fourier amplitudes.
        n (int): Number of points in the time domain signal.

    Returns:
        (float): Probability density at the provided modulus value.
    """
    return (2 * n / sigma**2)**2 * r * ss.kn(0, 2 * r * n / sigma**2)


def white_noise_product_modulus_mean(sigma, n):
    return (np.pi/4) * sigma**2


def test_white_noise_modulus_pdf(n):
    """Test the pdf for the modulus of white noise DFT coefficient.

    We test with noise whose time domain distribution is 1/2 probability of
    +1 and -1. This is not Gaussian noise. The statistics of the modulus of the
    DFT coefficients does not depend on the time domain distribution of the
    noise.

    Args:
        n (int): Number of points in time domain.

    Returns:
    """
    x = (np.random.binomial(1, 0.5, n) * 2) - 1
    # Uncorrelated random samples, each +/- 1 with 1/2 probability of each.

    sigma = 1.0
    # The distribution of x has deviation of 1.

    ft_x = np.fft.fft(x) / n
    r_x = np.abs(ft_x)

    h, bin_edges = np.histogram(r_x, bins=500)
    bin_width = bin_edges[1] - bin_edges[0]
    bin_centers = np.convolve(
            bin_edges, np.array([0.5, 0.5]), mode='valid')

    expected = white_noise_modulus_pdf(
            bin_centers,
            sigma,
            n) * n * bin_width

    plt.semilogy(bin_centers, h, '.', markersize=8)
    plt.semilogy(bin_centers, expected, 'r-')


def test_white_noise_modulus_squared_pdf(n, sigma=1, plot_func=plt.semilogy):
    """Like the previous function, but for the mod square."""
    x = np.random.normal(0, sigma, n)
    ft_x = np.fft.fft(x) / n
    r_x_squared = np.abs(ft_x)**2

    h, bin_edges = np.histogram(r_x_squared, bins=500)
    bin_width = bin_edges[1] - bin_edges[0]
    bin_centers = np.convolve(
            bin_edges, np.array([0.5, 0.5]), mode='valid')

    expected = white_noise_modulus_squared_pdf(
            bin_centers,
            sigma,
            n) * n * bin_width

    plot_func(bin_centers, h, '.', markersize=8)
    plot_func(bin_centers, expected, 'r-')


def test_white_noise_product_modulus_pdf(n):
    """Test pdf of the modulus of the product of white noise DFT coefficients.

    Args:
        n (int): Number of time domain samples of each of the two sequences we
            use to generate the DFT coefficients which are multiplied together.
    """
    x = (np.random.binomial(1, 0.5, n) * 2) - 1
    y = (np.random.binomial(1, 0.5, n) * 2) - 1

    ft_x = np.fft.fft(x) / n
    ft_y = np.fft.fft(y) / n

    sigma = 1.0

    c = ft_x * ft_y.conjugate()

    h, bin_edges = np.histogram(np.abs(c), bins=500)
    bin_centers = np.convolve(bin_edges, np.array([0.5, 0.5]), mode='valid')
    bin_width = bin_edges[1] - bin_edges[0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(bin_centers, h, '.', markersize=12)
    r = np.linspace(0, max(bin_centers), 500)
    ax.semilogy(
            r,
            white_noise_product_modulus_pdf(r, sigma, n) * bin_width * n,
            color='r',
            linewidth=2)


def test_gaussian(n):
    """Show DFT of binary white noise has same stats as Gaussian white noise.

    Args:
        n (int): Number of time domain samples.

    We form a sequence of random values each +1 or -1 with equal probability.
    
    """
    data = np.random.binomial(1, 0.5, n)
    data = (data * 2) - 1
    ft = np.fft.fft(data) / n

    re = np.real(ft)
    h, bin_edges = np.histogram(re, bins=500)
    bin_centers = np.convolve(bin_edges, np.array([0.5, 0.5]), mode='valid')
    bin_width = bin_edges[1] - bin_edges[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(bin_centers, h/(n*1.0), '.', markersize=12)
    x = np.linspace(min(bin_centers), max(bin_centers), 500)
    sigma = 1.0 / np.sqrt(2 * n)
    ax.semilogy(
            x,
            (1.0/np.sqrt(2*np.pi*sigma**2))*np.exp(
                    -x**2 / (2 * sigma**2))*bin_width,
            color='r',
            linewidth=2)


def test_theory(n, dt, plotfunc=plt.loglog):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Sequence of uncorrelated +1 and -1.
    x = (2 * np.random.binomial(1, 0.5, n)) - 1

    cs, frequencies = pasp.self_cross_spectrum(x, dt)

    frequencies_averaged, cs_averaged, bins = pasp.log_average(
            frequencies[1:],
            cs[1:],
            bins_per_decade=10)

    sigma = 1

    # Simulation csd, not averaged
    plotfunc(frequencies,
             np.abs(cs),
             '.',
             markersize=8,
             label='no averaging - simulated')

    # Theoretical csd (mean), not averaged
    plotfunc(frequencies_averaged,
             np.zeros(len(frequencies_averaged)) + (2 * (2 * dt) * (np.pi/4)*sigma**2),
             color='g',
             linewidth=4,
             label='no averaging - theory')

    # Simulation csd, averaged
    plotfunc(frequencies_averaged,
             np.abs(cs_averaged),
             color='r',
             marker='o',
             linewidth=0,
             markersize=12,
             label='Log averaged - simulated')

    # Theoretical csd, averaged
    psd_factor = 4 * dt
    # Converts between a mod-squared DFT (with the 1/n normalization) and a
    # power spectral density where the sampling interval is dt.

    plotfunc(frequencies_averaged,
             psd_factor * white_noise_product_modulus_mean(sigma, n) / np.sqrt(bins),
             color='r',
             linewidth=4,
             label='log averaged - theory')

    ax.grid(which='both')
    lgnd = ax.legend(loc='lower left', numpoints=1, prop={'size':30})
    ax.set_xlabel('Frequency', fontsize=32)
    ax.set_ylabel('Cross spectrum', fontsize=32)
    ax.tick_params(axis='both', which='major', labelsize=32)

