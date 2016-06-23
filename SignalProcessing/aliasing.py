import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


pi = np.pi


LABEL_FONT_SIZE = 20
TICK_WIDTH = 3
TICK_LENGTH = 6
TICK_LABEL_FONT_SIZE = 40


def make_presentable(
        ax,
        xlabel=None,
        ylabel=None,
        label_font_size=LABEL_FONT_SIZE,
        tick_width=TICK_WIDTH,
        tick_length=TICK_LENGTH,
        tick_color='k',
        tick_label_font_size=TICK_LABEL_FONT_SIZE,
        x_ticks=None,
        y_ticks=None,
        grid_line_width=2,
        grid_color='k',
        y_limits=None,
        x_limits=None,
        title=''):
    """Format plots for inclusion in technical documentation.

    Args:
        ax (matplotlib axis object): Axis onto which the plot is drawn.
        xlabel (string): x axis label. If None (the default) no label is used.
        ylabel (string): y axis label. If None (the default) no label is used.
        label_font_size (int): Font size for axis labels.
        tick_width (int): Width of tick marks.
        tick_length (int): Length of tick marks.
        tick_color (string): Tick color.
        tick_label_font_size (int): Size of tick labels.
        x_ticks (iterable of float): Values at which to draw x ticks.
        y_ticks (iterable of float): Values at which to draw y ticks.
        grid_line_width (int): Width of grid lines (default 2).
        grid_color (string): Grid line color (default 'k'=black).
        x_limits (tuple of float): (left, right) edges of plot.
        y_limits (tuple of float): (bottom, top) edge of plot.

    Nothing is returned. The axis object is updated.
    """
    for axis, ticks, label in [(ax.xaxis, x_ticks, xlabel), (ax.yaxis, y_ticks, ylabel)]:
        if label is not None:
            axis.set_label_text(label)
        # ticks
        if ticks is not None:
            axis.set_ticks(ticks)
            axis.set_ticklabels([str(t) for t in ticks])
        axis.set_tick_params(length=tick_length, width=tick_width,
            color=tick_color, labelsize=tick_label_font_size)
        # Label font
        label = axis.get_label()
        label.set_fontsize(label_font_size)
    # grid
    ax.grid(linewidth=grid_line_width, color=grid_color)
    # x and y limits
    if x_limits is not None:
        left, right = x_limits
        ax.set_xlim(left=left, right=right)
    if y_limits is not None:
        bottom, top = y_limits
        ax.set_ylim(bottom=bottom, top=top)
    ax.set_title(title, fontsize=label_font_size)


def plot(f_a, f_b, n, dot_size=40, line_width=9, x_buffer=0.2):
    result = points(f_a, f_b, n, x_buffer=x_buffer)
    xs, ys_a, ys_b, xs_continuous, ys_a_continuous, ys_b_continuous = result

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(xs_continuous, ys_a_continuous.real, 'b-', linewidth=line_width)
    ax.plot(xs_continuous, ys_a_continuous.imag, 'b--', linewidth=line_width)
    ax.plot(xs_continuous, ys_b_continuous.real, 'r-', linewidth=line_width)
    ax.plot(xs_continuous, ys_b_continuous.imag, 'r--', linewidth=line_width)

    ax.plot(xs, ys_a.real, 'k.', markersize=dot_size)
    ax.plot(xs, ys_a.imag, 'ks', markersize=dot_size/2)
    ax.plot(xs, ys_b.real, 'k.', markersize=dot_size)
    ax.plot(xs, ys_b.imag, 'ks', markersize=dot_size/2)

    make_presentable(ax, x_limits=(-n-x_buffer, n+x_buffer),
                     y_limits=(-1.1, 1.1))

    return ax


def points(f_a, f_b, n, x_buffer):
    xs = np.linspace(-n, n, (2*n) + 1)
    ys_a = np.exp(2j * pi * f_a * xs)
    ys_b = np.exp(2j * pi * f_b * xs)

    xs_continuous = np.linspace(-n-x_buffer, n+x_buffer, 1000)
    ys_a_continuous = np.exp(2j * pi * f_a * xs_continuous)
    ys_b_continuous = np.exp(2j * pi * f_b * xs_continuous)
    return xs, ys_a, ys_b, xs_continuous, ys_a_continuous, ys_b_continuous


def plot_3d(fa, fb, n, dot_size=24, line_width=3):
    x, za, zb, x_cont, za_cont, zb_cont = points_3d(fa, fb, n)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(za.real, za.imag, x, 'k.', markersize=dot_size)
    ax.plot(zb.real, zb.imag, x, 'k.', markersize=dot_size)
    ax.plot(za_cont.real, za_cont.imag, x_cont, 'b-', linewidth=line_width,
            alpha=0.5)
    ax.plot(zb_cont.real, zb_cont.imag, x_cont, 'r-', linewidth=line_width,
            alpha=0.5)

    return ax


def points_3d(fa, fb, n):
    x = np.linspace(-n, n, (2*n) + 1)
    za = np.exp(2j * pi * fa * x)
    zb = np.exp(2j * pi * fa * x)

    x_cont = np.linspace(-n, n, 1000)
    za_cont = np.exp(2j * pi * fa * x_cont)
    zb_cont = np.exp(2j * pi * fb * x_cont)

    return x, za, zb, x_cont, za_cont, zb_cont
