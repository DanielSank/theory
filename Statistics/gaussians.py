from __future__ import division

import numpy as np
import matplotlib.pyplot as plt


def compare(samples, pdf, pdf_args=(), bins=100, plotter=plt.plot):
    """Compare a pdf against a numerical simulation.

    Args:
        samples (ndarray): Set of numbers produced by a random process.
        pdf (function): The pdf we think describes the samples. This function
            takes at least one argument, which is the value of the random
            process, and returns the probability density for that value.
        pdf_args (tuple): Additional arguments passed to pdf.
        bins (int): Number of histogram bins to use when binning the samples.
        plotter (function): Plot function used to compare the samples against
            the pdf.
    """
    n = len(samples)
    counts, edges = np.histogram(samples, bins=bins)
    bin_width = edges[1] - edges[0]
    bin_centers = (edges[1:] + edges[0:-1]) / 2
    plotter(bin_centers, counts / bin_width / n, 'b.')
    plotter(bin_centers, pdf(bin_centers, *pdf_args), 'r-')


def gaussian_pdf(x, sigma):
    """PDF of a one dimensional Gaussian distribution."""
    return (1/np.sqrt(2*np.pi*sigma**2)) * np.exp(-x**2 / (2 * sigma**2))


def gaussian_squared_pdf(a, sigma):
    """PDF of the square of a one dimensional Gaussian distribution."""
    return (1/np.sqrt(2*np.pi*sigma**2)) * np.exp(-a/(2*sigma**2))/np.sqrt(a)


def mod_z_squared_pdf(a, sigma):
    """PDF of mod squared of a two dimensional Gaussian distributed variable.

    Suppose X and Y are both Gaussian random variables with stdandard deviation
    sigma. Now define Z = X + iY. This function gives the distribution of |Z|^2.
    """
    return (1/(2*sigma**2)) * np.exp(-a / (2 * sigma**2))


def mod_z_pdf(a, sigma):
    return (1/sigma**2) * a * np.exp(-a**2 / (2 * sigma**2))

