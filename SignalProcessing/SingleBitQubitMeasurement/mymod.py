from __future__ import division

from typing import TypeVar
from typing import cast
import numpy as np
import scipy
import scipy.special as ss
import matplotlib.pyplot as plt
import pyle.analysis.signal_processing as pasp


FloatOrArray = TypeVar("FloatOrArray", float, np.ndarray)

PI = np.pi
TWOPI = 2 * PI


def white_noise_modulus_pdf(r: FloatOrArray, sigma: float, n: int) -> FloatOrArray:
    """PDF of modulus of the modulus of the normalized DFT of white noise.

    Suppose we have a sequence of random numbers, each of which is indepdently
    distributed from a distribution with standard deviation sigma. For example,
        x = np.random.normal(0, 1, n).
    Now we construct the normalized discrete Fourier transform
        ft = (1 / n) np.fft.fft(x).
    The numbers ft are complex random numbers. The probability density of the
    modulii of these numbers is given by this function.

    Note that there is no assumption about the distribution of the noise in the
    time domain (i.e. they don't have to be Gaussian), so long as the samples are
    uncorrelated.

    Args:
        r: Modulus at which to compute the PDF.
        sigma: Standard deviation of the time domain signal from whence came the
            normalized DFT coefficients.
        n: Number of points in the time domain signal.

    Returns:
        Probability density at the provided modulus value.
    """
    q = sigma**2 / n
    return (2 / q) * r * np.exp(- r**2 / q)


def white_noise_modulus_pdf_demo() -> None:
    """Show that white_noise_modulus_pdf is correct"""
    sigma = 1
    n = 1_000_000
    x = np.random.normal(0, sigma, n)
    modulus = np.abs(np.fft.fft(x) / n)
    pdf, bin_edges = np.histogram(modulus, bins=100, density=True)
    bin_centers = np.convolve(bin_edges, np.array([0.5, 0.5]), mode="valid")

    _, ax = cast(tuple[plt.Figure, plt.Axes], plt.subplots())
    ax.semilogy(
            bin_centers,
            pdf,
            color="tab:blue",
            marker=".",
            markersize=10,
            linewidth=0,
            label="Simulation",
    )
    (sigma_fit, n_fit), _ = scipy.optimize.curve_fit(
            f=white_noise_modulus_pdf,
            xdata=bin_centers,
            ydata=pdf,
            p0=(1, n),
    )
    ax.semilogy(
            bin_centers,
            white_noise_modulus_pdf(bin_centers, sigma=sigma_fit, n=n_fit),
            markersize=0,
            linewidth=2,
            color="tab:orange",
            label="Theory",
    )
    ax.set_xlabel("Modulus of normalized Fourier coefficient")
    ax.set_ylabel("PDF")
    ax.legend()
    ax.grid()


def white_noise_modulus_cdf(r: FloatOrArray, sigma: float, n: int) -> FloatOrArray:
    """Integral of previous function."""
    q = sigma**2 / n
    return 1 - np.exp(-r**2 / q)


def white_noise_modulus_squared_pdf(r_squared: FloatOrArray, sigma: float, n: int) -> FloatOrArray:
    """PDF of modulus square of the normalized DFT of white noise.

    Suppose we have a sequence of random numbers, each of which is indepdently
    distributed from a distribution with standard deviation sigma. For example,
        x = np.random.normal(0, 1, n).
    Now we construct the normalized discrete Fourier transform
        ft = (1 / n) np.fft.fft(x).
    The numbers ft are complex random numbers. The probability density of the
    modulii squared of these numbers is given by this function.

    Note that there is no assumption about the distribution of the noise in the
    time domain (i.e. they don't have to be Gaussian), so long as the samples are
    uncorrelated.

    Args:
        r_squared: Modulus squared at which to compute the PDF.
        sigma: Standard deviation of the time domain signal from whence came the
            normalized DFT coefficients.
        n: Number of points in the time domain signal.

    Returns:
        Probability density at the provided modulus squared value.
    """
    q = sigma**2 / n
    return (1 / q) * np.exp(-r_squared / q)


def white_noise_modulus_squared_pdf_demo() -> None:
    """Show that white_noise_modulus_squared_pdf is correct"""
    sigma = 1
    n = 1_000_000
    x = np.random.normal(0, sigma, n)
    modulus_squared = np.abs(np.fft.fft(x) / n)**2
    pdf, bin_edges = np.histogram(modulus_squared, bins=100, density=True)
    bin_centers = np.convolve(bin_edges, np.array([0.5, 0.5]), mode="valid")

    _, ax = cast(tuple[plt.Figure, plt.Axes], plt.subplots())
    ax.semilogy(
            bin_centers,
            pdf,
            marker=".",
            linewidth=0,
            markersize=10,
            color="tab:blue",
            label="Simulation",
    )
    (sigma_fit, n_fit), _ = scipy.optimize.curve_fit(
            f=white_noise_modulus_squared_pdf,
            xdata=bin_centers,
            ydata=pdf,
            p0=(1, n),
    )
    ax.semilogy(
            bin_centers,
            white_noise_modulus_squared_pdf(bin_centers, sigma=sigma_fit, n=n_fit),
            markersize=0,
            linewidth=2,
            color="tab:orange",
            label="Theory",
    )
    ax.set_xlabel("Modulus squared of normalized Fourier coefficient")
    ax.set_ylabel("PDF")
    ax.legend()
    ax.grid()


def white_noise_product_modulus_pdf(r: FloatOrArray, sigma: float, n: int) -> FloatOrArray:
    """PDF of modulus of product of two normalized white noise DFT coefficients.

    Suppose we have two white noise time sequences
        a = [b_0, a_1, ..., a_{n-1}] and
        b = [a_0, b_1, ..., b_{n-1}]
    each with time domain standard deviation equal to sigma.

    We form the normalized discrete Fourier transforms of each sequence
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
        r: Modulus at which you want the probability density.
        sigma: Standard deviation of the time domain signal from whence came the
            Fourier amplitudes.
        n: Number of points in the time domain signal.

    Returns:
        Probability density at the provided modulus value.
    """
    return (2 * n / sigma**2)**2 * r * ss.kn(0, 2 * r * n / sigma**2)


def white_noise_product_modulus_pdf_demo(n: int) -> None:
    """Show that white_noise_product_modulus_pdf is correct.

    Args:
        n: Number of time domain samples of each of the two sequences we use to
            generate the DFT coefficients which are multiplied together.
    """
    x = (np.random.binomial(1, 0.5, n) * 2) - 1
    y = (np.random.binomial(1, 0.5, n) * 2) - 1

    ft_x = np.fft.fft(x) / n
    ft_y = np.fft.fft(y) / n

    sigma = 1.0

    cross_spectrum = ft_x * ft_y.conjugate()

    h, bin_edges = np.histogram(np.abs(cross_spectrum), bins=500, density=True)
    bin_centers = np.convolve(bin_edges, np.array([0.5, 0.5]), mode="valid")

    _, ax = cast(tuple[plt.Figure, plt.Axes], plt.subplots())
    ax.semilogy(bin_centers, h, '.', markersize=12, label="Simulation")
    r = np.linspace(0, max(bin_centers), 500)
    ax.semilogy(
            r,
            white_noise_product_modulus_pdf(r, sigma, n),
            color="tab:red",
            linewidth=2,
            label="Theory",
    )
    ax.set_xlabel("Modulus of product of white noise DFT's")
    ax.set_ylabel("PDF")
    ax.legend()
    ax.grid()


def white_noise_product_modulus_mean(sigma: float) -> float:
    return (np.pi/4) * sigma**2


def test_gaussian(n: int):
    """Compare stats of normalized DFT of binary and Gaussian white noise.

    Args:
        n: Number of time domain samples.
    """
    data = (2 * np.random.binomial(1, 0.5, n)) - 1
    ft = np.fft.fft(data) / n

    re = np.real(ft)
    pdf, bin_edges = np.histogram(re, bins=500, density=True)
    bin_centers = np.convolve(bin_edges, np.array([0.5, 0.5]), mode='valid')
    bin_width = bin_edges[1] - bin_edges[0]

    _, ax = cast(tuple[plt.Figure, plt.Axes], plt.subplots())

    ax.semilogy(bin_centers, pdf, '.', markersize=12, color="tab:blue", label="Simulation")

    x = np.linspace(min(bin_centers), max(bin_centers), 500)
    sigma_re = 1.0 / np.sqrt(2 * n)  # Standard deviation of real part of ft
    ax.semilogy(
            x,
            (1.0 / np.sqrt(TWOPI * sigma_re**2)) * np.exp(-x**2 / (2 * sigma_re**2)),
            color="tab:red",
            linewidth=2,
            label="Theory",
    )
    ax.set_xlabel("Real part of normalized DFT")
    ax.set_ylabel("PDF")
    ax.legend()
    ax.grid()


def test_theory(n: int, dt: float, bins_per_decade: int, plotfunc=plt.loglog):
    _, ax = cast(tuple[plt.Figure, plt.Axes], plt.subplots())

    # Sequence of uncorrelated +1 and -1.
    x = (2 * np.random.binomial(1, 0.5, n)) - 1

    cs, frequencies = pasp.self_cross_spectrum(x, dt)
    sigma = 1

    frequencies_averaged, cs_averaged, bins = pasp.log_average(
            frequencies[1:],
            cs[1:],
            bins_per_decade=bins_per_decade)

    # Simulation csd, not averaged
    indices = np.unique(
            np.logspace(0, np.log10(len(frequencies)-1), 1000).astype(int))
    plotfunc(
            frequencies[indices],
            np.abs(cs[indices]),
            marker='o',
            markersize=12,
            alpha=0.4,
            label='no averaging - simulated',
    )

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
             alpha=1.0,
             label='Log averaged - simulated')

    # Theoretical csd, averaged
    psd_factor = 4 * dt
    # Converts between a mod-squared DFT (with the 1/n normalization) and a
    # power spectral density where the sampling interval is dt.

    plotfunc(
            frequencies_averaged,
            psd_factor * white_noise_product_modulus_mean(sigma) / np.sqrt(bins),
            color='r',
            linewidth=4,
            label='log averaged - theory',
    )

    ax.grid(which='both')
    ax.legend(loc='lower left', numpoints=1, prop={'size':30})
    ax.set_xlabel('Frequency', fontsize=32)
    ax.set_ylabel('Cross spectrum', fontsize=32)
    ax.tick_params(axis='both', which='major', labelsize=32)
