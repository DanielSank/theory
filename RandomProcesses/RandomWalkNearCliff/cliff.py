import numpy as np
import matplotlib.pyplot as plt


LABEL_FONT_SIZE = 40
TICK_WIDTH = 2
TICK_LENGTH = 4
TICK_LABEL_FONT_SIZE = int(LABEL_FONT_SIZE * 0.75)
LEGEND_FONT_SIZE = int(LABEL_FONT_SIZE * 0.6)


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
        title='',
        legend=True,
        legend_loc='upper left',
        legend_font_size=LEGEND_FONT_SIZE):
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
    for axis, ticks, label in [(ax.xaxis, x_ticks, xlabel),
                               (ax.yaxis, y_ticks, ylabel)]:
        if label is not None:
            axis.set_label_text(label)
        # ticks
        if ticks is not None:
            axis.set_ticks(ticks)
            axis.set_ticklabels([str(t) for t in ticks])
        axis.set_tick_params(length=tick_length, width=tick_width,
            color=tick_color, labelsize=tick_label_font_size)
        # Label font
        axis.get_label().set_fontsize(label_font_size)
    # grid
    ax.grid(linewidth=grid_line_width, color=grid_color)
    # x and y limits
    if x_limits is not None:
        left, right = x_limits
        ax.set_xlim(left=left, right=right)
    if y_limits is not None:
        bottom, top = y_limits
        ax.set_ylim(bottom=bottom, top=top)
    if legend:
        ax.legend(prop={'size': legend_font_size}, loc=legend_loc)

    ax.set_title(title, fontsize=label_font_size)


def P_minus_one(p, z):
    """Free generating function for one leftward step"""
    q = 1 - p
    f = np.sqrt(1 - 4*p*q*z**2)
    return (1 - f)/(2*q*z * f)


def P_zero(p, z):
    """Free generating function for zero steps, i.e. loop-back"""
    q = 1 - p
    f = np.sqrt(1 - 4*p*q*z**2)
    return 1 / f


def F_minus_one(p, z):
    """Absorbing generating function for one leftward step"""
    return P_minus_one(p, z) / P_zero(p, z)


def plot_F_minus_one(ps, zs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for z in zs:
        plt.plot(ps, F_minus_one(ps, z), label='$z={}$'.format(z),
                 linewidth=4)
    make_presentable(ax,
                     xlabel='$p$',
                     ylabel='Absorption probability')
    return ax


def plot_P(ps, zs):
    colors = ['b', 'g', 'r', 'k', 'y']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for z, color in zip(zs, colors):
        plt.semilogy(ps, P_minus_one(ps, z),
                     color+'-',
                     label='$z={}$'.format(z),
                     linewidth=4)
        plt.semilogy(ps, P_zero(ps, z),
                     color+'--',
                     linewidth=4)
        make_presentable(ax,
                         xlabel='$p$',
                         ylabel='Generating function value')
