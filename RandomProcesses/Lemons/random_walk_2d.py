import matplotlib.pyplot as plt
import numpy as np


def simulate_walk(number_of_steps):
    dx = np.random.normal(size=number_of_steps)
    dy = np.random.normal(size=number_of_steps)

    x = np.cumsum(dx)
    y = np.cumsum(dy)

    return x, y


def plot_walk(x, y, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.plot(
        x,
        y,
        color='b',
        marker='.',
        markersize=8,
        linewidth=2,)
    ax.grid(True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    return ax
