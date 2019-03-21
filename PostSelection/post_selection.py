import matplotlib.pyplot as plt
import numpy as np

from scipy.special import erf, erfc

RT2 = np.sqrt(2)
two_rt_two = 2 * np.sqrt(2)

def kept_data_fraction(snr, r):
    return 1 - 0.5 * (
        erf(np.sqrt(snr)/2 + (r / two_rt_two)) +
        erf(-np.sqrt(snr)/2 + (r / two_rt_two))
    )

def error_raw(snr):
    return 0.5 * erfc(np.sqrt(snr) / 2)

def error_post_selected(snr, r):
    return 0.5 * erfc(np.sqrt(snr)/2 + r / two_rt_two)

def fidelity_gain(snr, r):
    return (1 - error_post_selected(snr, r)) / (1 - error_raw(snr))


def go():
    _, ax = plt.subplots()
    rs = np.linspace(1, 6)
    for snr in [4, 8, 12, 16]:
        fidelity_gain = (1 - error_post_selected(snr, rs)) / (1 - error_raw(snr))
        stats_loss = np.sqrt(kept_data_fraction(snr, rs))
        plt.plot(
            rs,
            fidelity_gain * stats_loss,
            label='SNR={}'.format(snr),)
    ax.legend(loc='lower left')
    plt.grid()
    plt.xlabel('r')
    plt.ylabel('supremacy SNR gain')

def compare(snrs):
    _, axes = plt.subplots(2, 1)
    rs = np.linspace(0, 6)
    for snr in snrs:
        axes[0].plot(
            rs,
            fidelity_gain(snr, rs) - 1,
            label='SNR={}'.format(snr),)
        axes[0].set_xlabel('r')
        axes[0].set_ylabel('fidelity fractional increase')
        axes[0].grid()

        axes[1].plot(
            rs,
            np.sqrt(kept_data_fraction(snr, rs)),
            label='SNR={}'.format(snr),)
        axes[1].set_xlabel('r')
        axes[1].set_ylabel('fraction of data kept')
        axes[1].grid()
        axes[1].legend(loc='lower left')
