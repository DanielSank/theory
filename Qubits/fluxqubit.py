import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

import functools
import numpy as np
import scipy.interpolate as interpolate
import scipy.optimize as opt
import itertools

import labrad
import labrad.units as U
from labrad.units import Value, ValueArray as VA

us, GHz, Ohm, pH = (U.Unit(s) for s in ['us', 'GHz', 'Ohm', 'pH'])


# Physical constants

pi = np.pi
PHI_0 = Value(2E-15, 'V*s')
HBAR = Value(1.05E-34, 'J*s')
R_K = Value(25.8E3, 'Ohm')


# Data management

NUM_BETAS = 10
NUM_PHI_POINTS = 200

DATA_NAMES = ['Z_LC', 'beta', 'f_LC', 'phi_x', 'f_10', 'f_21',
    'phi_left', 'phi_right', 'phi_max', 'phi_10', '|0>', '|1>']

DATA_TYPES = ['f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
              'f8',
              '(%d,)f8'%(NUM_PHI_POINTS,), '(%d,)f8'%(NUM_PHI_POINTS,)]

LABELS = { # data key, label, unit string, unit label
    'beta': (r'$\beta$', '', ''),
    'f_LC': (r'$f_{LC}$', 'GHz', 'GHz'),
    'Z_LC': (r'$Z_{LC}$', 'Ohm', r'$\Omega$'),
    'phi_x': (r'$\phi_x$', '', '')
    }


# Plot styling

LABEL_FONT_SIZE = 36
TICK_LABEL_FONT_SIZE = 32
SUBPLOT_LABEL_FONT_SIZE = 32
TICK_WIDTH = 4
TICK_LENGTH = 4
CB_TICK_LABEL_FONT_SIZE = 18
MARKER_SIZE = 8
LEGEND_FONT_SIZE = 22
LINE_WIDTH = 4


# Default analysis parameters

FREQUENCIES = VA([20], 'GHz')
FREQUENCY = FREQUENCIES[0]

IMPEDANCES = VA([50, 150], 'Ohm')
IMPEDANCE = IMPEDANCES[0]

RESONATOR_FREQUENCY = 8*GHz
RESONATOR_IMPEDANCE = 50*Ohm


# plotting/analysis utilities

def make_presentable(ax, xlabel=None, ylabel=None,
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
                    x_limits=None):
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
    for axis, ticks, label in zip((ax.xaxis, ax.yaxis), (x_ticks, y_ticks),
                                  (xlabel, ylabel)):
        if label is not None:
            axis.set_label_text(label)
        # ticks
        if ticks is not None:
            axis.set_ticks(ticks)
            axis.set_ticklabels([str(t) for t in ticks])
        axis.set_tick_params(length=tick_length ,width=tick_width,
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


def collapse(z):
    """Get a sorted copy of a collection with all duplicates removed"""
    s = set(z)
    l = np.array([x for x in s])
    l.sort()
    return l


def plot_iter(ax, data, const_names, indep, deps,
    plot_type='linear', args=None):
    """Generate plots while holding several parameters constant.

    Args
        ax: The axis onto which the plot is drawn.
        data (np.ndarray): Array of numerical data. Columns are named according
            to DATA_NAMES.
        const_names (iterable of string): Names of variables held constant.
        indep ():
        deps (list of dict):
        plot_type (string): Axis scaling. Options are 'linear' (the default),
            'semilogx', 'semilogy', and 'loglog'.

    This function draws several curves on a single axis. The number of curves is
    equal to the number of combinations of constants.
    """
    constants = (collapse(data[name]) for name in const_names)
    labels = set()
    colors = ['b', 'r', 'g', 'k', 'c', 'y']
    # Loop over all combinations of constants. For each combination we create
    # unique legend label.
    for i, c in enumerate(itertools.product(*constants)):
        legend_label = ''
        dat = data[:]  # copy data
        for val, name in zip(c, const_names):
            dat = dat[dat[name] == val]
            # Add a section to the legend label indicating constant's value
            label_name, unit, unit_text = LABELS[name]
            if label_name == r'$f_{LC}$':
                label_name = r'$\omega_{LC} / 2\pi$'
            legend_label += r'%s$=$%.2f %s  '%(label_name, val, unit_text)
        for dep in deps:
            d = dep.copy()
            func = d.pop('func')
            d['label'] = d.get('label','') + '  ' + legend_label
            # Hack to avoid labelling each branch of a multi-branch plot
            if d['label'] in labels:
                d.pop('label')
            else:
                labels.add(d['label'])

            d.setdefault('linewidth', 4)
            plotfunc = {'linear': ax.plot,
                        'semilogy': ax.semilogy,
                        'semilogx': ax.semilogx,
                        'loglog': ax.loglog}[plot_type]
            if not 'color' in d.keys():
                d['color'] = colors[i]
            if args is None:
                plotfunc(indep(dat), func(dat), **d)
            else:
                plotfunc(indep(dat), func(dat, *args), **d)


# Simulation utilities

def diagonalize_H_over_params(Z_LCs, betas, f_LCs, phi_xs):
    """Diagonalize the flux qubit Hamiltonian over a set of parameters

    Args:
        Z_LCs (ValueArray[Ohm]): Harmonic limit impedances
        betas (np.ndarray): betas
        f_LCs (ValueArray[GHz]): Harmonic limit frequencies
        phi_xs (np.ndarray): Dimensionless tilt fluxes.

    Returns (phi, data). phi is the discretized  dimensionless flux
    phi = 2 pi Phi / Phi_0
    a.k.a the superconducting phase, which we use as the coordinate for the flux
    qubit Hamiltonian. data is a structured ndarray with the following columns:
        Z_LC: Impedance of device in harmonic limit.
        beta
        f_LC: Frequency of well in harmonic limit.
        phi_x: Dimensionless tilt flux.
        f_10: Frequency difference between |1> and |0>
        f_21: Frequency difference between |2> and |1>
        phi_left: phi value at left well minimum
        phi_right: phi value at right well minimum
        phi_max: phi value at peak of interwell barrier
        matrix_phi_10
        |0>: Ground state wave function
        |1>: First excited state wave function
    """
    num_rows = len(Z_LCs) * len(betas) * len(f_LCs) * len(phi_xs)
    data = np.zeros((num_rows,), dtype={'names':DATA_NAMES,
        'formats':DATA_TYPES})

    for i, (beta, Z_LC, f_LC, phi_x) in \
        enumerate(itertools.product(betas, Z_LCs, f_LCs, phi_xs)):
        omega_LC = 2*pi*f_LC
        E_L = Inductor.L_to_E_L(Z_LC / omega_LC)
        phi, potential, kinetic = hamiltonian_position_space(
            beta,
            Z_LC,
            phi_x=phi_x
        )
        vals, vecs = diagonalize_H(potential + kinetic)

        # Point |0> upward
        if vecs[:,0][NUM_PHI_POINTS/2] < 0:
            vecs[:,0] = -vecs[:,0]

        # Find well minima
        if beta > 1:
            phi_left, phi_right = minima(beta, phi_x)
        else:
            phi_left, phi_right = (np.nan, np.nan)

        # Find position of barrier peak
        if beta > 1:
            phi_max = maximum(beta, phi_x)
        else:
            phi_max = np.nan

        f_10, f_21 = np.diff(vals)[0:2] * (E_L / (2*pi*HBAR))
        phi_10 = np.abs(np.sum(vecs[:,0] * phi * vecs[:,1]))
        # TODO: Find a better way to get the units right
        data[i] = (
            Z_LC['Ohm'],
            beta,
            f_LC['GHz'],
            phi_x,
            f_10['GHz'],
            f_21['GHz'],
            phi_left,
            phi_right,
            phi_max,
            phi_10,
            vecs[:,0],
            vecs[:,1]
        )
    return phi, data


def hamiltonian_position_space(beta, Z_LC, phi_x=0):
    """The phi (phase) basis Hamiltonian for the flux qubit in units of E_L.

    Args:
        beta (float): Ratio of linear to junction inductance (see below)
        Z_LC (Value[Ohm]): Harmonic limit impedance.
        phi_x (float): Dimensionless tilt flux (default 0).

    Returns (phi, potential, kinetic). phi is an array of the discrete phase
        points over which the Hamiltonian is discretized. potential is the
        (diagonal) matrix representing the potential energy evaluated at the
        points given by phi. kinetic is the matrix representing the potential
        energy evaluated over the points given by phi.

    The computation is done by discretizing the Hamiltonian. From basic circuit
    theory you can find the Hamiltonian in the phase (phi) basis:
    
    H/E_L = (1/2)(phi - phi_x)^2 + beta(cos(phi) - 1) + (1/2) 8 (E_C/E_J)q^2

    where [phi, q] = i
    E_L = (Phi_0/2pi)^2 / L
    E_C = e^2 / 2C
    beta = L / L_{J_0}
    L_{J_0} is the zero bias inductance of the junction.

    From the commutation relation [phi, q] = i we deduce that in the phi basis
    q => -i d/dphi. Then q^2 = -(d/dphi)^2. Discretizing with lattice spacing
    dphi gives

    (d/dphi)^2 => [f(phi + dphi) - 2f(phi) + f(phi - dphi)] / dphi^2

    As a matrix, this has -2/dphi^2 along the diagonal and 1/dphi^2 on the first
    off diagonals.
    """

    E_C_over_E_L = 2 * (2*pi * Z_LC / R_K)**2

    phi = np.linspace(-pi, pi, NUM_PHI_POINTS) #includes endpoints
    dphi = phi[1] - phi[0]

    q_squared = -(np.diag([-2]*NUM_PHI_POINTS) +
                  np.diag([1]*(NUM_PHI_POINTS-1), 1) +
                  np.diag([1]*(NUM_PHI_POINTS-1), -1)) / dphi**2

    charge_term = 4 * E_C_over_E_L * q_squared
    inductor_term = (0.5) * (np.diag((phi - phi_x)**2))
    junction_term = beta * (np.diag(np.cos(phi)) - np.eye(NUM_PHI_POINTS))

    potential = inductor_term + junction_term
    kinetic = charge_term

    return phi, potential, kinetic


def diagonalize_H(H):
    """Get ordered eigenvalues and eigenvectors of a Hamiltonian

    Args:
        H (np.ndarray): Matrix representation of the Hamiltonian.

    Returns (eigenvalues, eigenvectors):
        eigenvalues is a np.ndarray of the eigenvalues in ascending order.
        eigenvectors is a np.ndarray of the eigenvectors. It's shape is the same
        as H. Each column represents an eigenvector, and the columns are in the
        same order as their corresponding eigenvalues.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    #Sort eigenvalues and eigenvectors by increasing energy
    inds = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[inds]
    eigenvectors = eigenvectors[:,inds]
    return eigenvalues, eigenvectors


# Analytic formulae (see John's fluxmon write-up)

def potential_over_E_L(phi, phi_x, beta):
    """Flux qubit potential divided by E_L

    Args:
        phi (float): dimensionless flux (i.e. phase) points at which to
            evaluate the potential.
        phi_x (float): dimensionless tilt flux
        beta (float): ratio of linear and junction inductances

    Returns: Value of the flux qubit potential at the specified parameters, in
        units of E_L.
    """
    return ((phi - phi_x)**2 / 2) + beta * (np.cos(phi) - 1)


def potential_minima_approx(beta):
    """Approximate position of the well minima for zero tilt flux

    Args:
        beta (float): ratio of linear and junction inductances.

    Returns:
        Dimensionless flux (i.e. phase) position of well minima.
    """
    return np.sqrt(6*((beta-1)/beta) + (9.0/5)*((beta-1)/beta)**2)


def barrier_height_over_E_L_approx(beta):
    """Approximate height of interwell barrier in units of E_L

    Args:
        beta (float): Ratio of linear and junction inductances.

    Returns (float):
        Approximate height of interwell barrier, in units of E_L.
    """
    return (3.0/2) * (beta - 1)**2 / beta


def well_frequency_over_f_LC_approx(beta):
    """Approximate well frequency in units of f_LC

    Args:
        beta (float): Ratio of linear inductance to junction inductancel

    Returns (float):
        Approximate value of well frequency divided by f_LC.
    """
    return np.sqrt(2) * np.sqrt(beta - 1)


# Numerically computed properties

def barrier_top_objective_func(phi, phi_x, beta):
    """Function which attains a minimum at the interwell barrier peak.

    phi (float): Dimensionless flux (i.e. phase)
    phi_x (float): Dimensionless tilt flux.
    beta (float): Ratio of linear to junction inductance

    Returns minus the potential energy (in units of E_L).
    """
    return -potential_over_E_L(phi, phi_x, beta)


def minima(beta, phi_x):
    """Numerically find left and right minima of the potential.

    Args:
        beta (float): Ratio of linear to junction inductance.
        phi_x (float): Dimensionless tilt flux.

    Returns: (left minimum, right minimum), the locations of the well minima in
        dimensionless flux (i.e. phase) units.

    This function sometimes finds the zero derivative point at phi=0, which is a
    local max, not a min. Adding bounds to the opt.minimize call might fix this.
    """
    beta = float(beta)
    phi_x = float(phi_x)
    right = opt.minimize(potential_over_E_L,
        potential_minima_approx(beta),# bounds=((0.03, None)),
        args=(phi_x, beta))
    left = opt.minimize(potential_over_E_L,
        -potential_minima_approx(beta),# bounds=((None, -0.03)),
        args=(phi_x, beta))
    return left.x[0], right.x[0]


def maximum(beta, phi_x):
    """Find the position of the interwell barrier peak.

    Args:
        beta (float): ratio of linear to junction inductance.
        phi_x (float): dimensionless tilt flux.

    Returns (float):
        Location of barrier peak in dimensionless flux (i.e. phase) units.
    """
    beta = float(beta)
    phi_x = float(phi_x)
    return opt.minimize(barrier_top_objective_func, 0,
        args=(phi_x, beta)).x[0]


def tunnel_time_times_f_LC(f_LC, Z_LC, phi_well, phi_barrier, beta, phi_x):
    """The tunnelling time out of a well, normalized to f_LC

    Args:
        f_LC (Value[Hz]): Harmonic limit frequency of the qubit.
        Z_LC (Value[Ohm]): Harmonic limit impedance of the qubit.
        phi_well (float): Dimensionless flux (i.e. phase) position of the well.
        phi_barrier (float): Dimensionless flux (i.e. phase) position of the
            interwell barrier peak.
        beta (float): Linear to junction inductance ratio.
        phi_x (float): Dimensionless tilt flux.

    Returns (float):
        Tunneling time multiplied by f_LC.

    The formula used here to compute the tunneling time assumes escape from the
    metastable well into a continuum of states. This is certainly _not_ the case
    in a real flux qubit, so the results here could be completely bogus.
    """

    L = Z_LC / (2*pi * f_LC)
    E_L = Inductor.L_to_E_L(L)
    a = 1.3
    b = 3.77

    well_bottom = potential_over_E_L(phi_well, phi_x, beta)
    well_top = potential_over_E_L(phi_barrier, phi_x, beta)
    barrier_height = E_L * (well_top - well_bottom)

    omega_well = 2*pi * f_LC * well_frequency_over_f_LC(phi_well, beta)
    omega_barrier = 2*pi * f_LC * \
        np.sqrt(-1.0 + beta * np.cos(phi_barrier))

    A = b * barrier_height / (HBAR * omega_barrier)
    # Eq. (14) John's writeup
    Delta = a * HBAR * omega_well * np.sqrt(A) * np.exp(-A)
    return (2*pi*HBAR / Delta) * f_LC


def number_of_states(f_LC, Z_LC, phi_well, phi_max, phi_x, beta):
    """Number of bounds states in a well.

    Args:
        f_LC (Value[Hz]): Harmonic limit frequency.
        Z_LC (Value[Ohm]): Harmonic limit impedance.
        phi_well (float): Dimensionless flux (i.e. phase) position of well.
        phi_max (float): Dimensionless flux (i.e. phase) position of interwell
            barrier peak.
        phi_x (float): Dimensionless tilt flux.
        beta (float): Ratio of linear to junction inductance.

    Returns (float):
        Approximate number of bound states in the well.
    """
    L = Z_LC / (2*pi*f_LC)
    E_L = Inductor.L_to_E_L(L)
    potential_bottom = potential_over_E_L(phi_well, phi_x, beta)
    potential_top = potential_over_E_L(phi_max, phi_x, beta)
    well_height_over_E_L = potential_top - potential_bottom
    f_well = f_LC * well_frequency_over_f_LC(phi_well, beta)
    well_energy_over_E_L = HBAR * 2*pi * f_well / E_L
    return well_height_over_E_L / well_energy_over_E_L


def well_frequency_over_f_LC(phi, beta):
        """Compute the well frequency from known position (phi), beta and f_LC.

        Args:
            phi: The bottom point of the well.
            beta: The beta for the potential.

        Returns:
            The frequency of the well.

        For a harmonic system with Hamiltonian
        H = (1/2)a u^2 + (1/2)b v^2,    [u,v]= i gamma
        the resonance frequency is
        hbar omega = gamma (a b)**1/2.
        Writing the potential energy as U = (1/2)a u^2, we can re-express the
        resonance frequency equation as
        hbar omega = gamma (U'' b)**1/2.

        For the flux qubit, the Hamiltonian is
        H = (1/2) E_L (phi - phi_x)**2
            + E_L beta (cos(phi) - 1)
            + (1/2) 8 E_c q^2
        with [phi, q] = i.
        Therefore, we have b = 8 E_C, gamma = 1, and
        U'' = E_L (1 - beta cos(phi)).

        Putting it together gives
        hbar omega = (8 E_C E_L (1 - beta cos(phi)))**1/2
                   = hbar omega_LC (1 - beta cos(phi))**1/2.
        or
        f = f_LC (1 - beta cos(phi))**1/2.
        """

        return np.sqrt(1 - beta * np.cos(phi))


# Dependent curves suitable for passing to plot_iter

def dep_omega_10_over_omega_LC(d):
    """|0> -> |1> frequency normalized to harmonic limit frequency

    Args:
        d (np.ndarray): Structured data array with named columns.

    Returns:
        |0> -> |1> frequency divided by f_LC.
    """
    return d['f_10'] / d['f_LC']


def dep_tunnelling_time_times_f_LC(which, d):
    """Tunneling time normalized to f_LC

    Args:
        d (np.ndarray): Structured data array with named columns.

    Returns (float):
        Tunneling time multiplied by f_LC.
    """
    f_LC = d['f_LC'] * U.Unit(LABELS['f_LC'][1])
    Z_LC = d['Z_LC'] * U.Unit(LABELS['Z_LC'][1])
    phi_well = d[which]
    phi_barrier = d['phi_max']
    return tunnel_time_times_f_LC(f_LC, Z_LC, phi_well, phi_barrier,
        d['beta'], d['phi_x'])


def dep_number_of_states(which, d):
    """Number of states in a well

    Args:
        which (string): 'phi_left' or 'phi_right', indicating which well to
            analyze.
        d (np.ndarray): Structured data array with named columns.

    Returns (float):
        Number of states in the well.
    """
    f_LC = d['f_LC'] * U.Unit(LABELS['f_LC'][1])
    Z_LC = d['Z_LC'] * U.Unit(LABELS['Z_LC'][1])
    phi = d[which]
    return number_of_states(f_LC, Z_LC, phi, d['phi_max'], d['phi_x'],
        d['beta'])


def dep_well_frequency_over_f_LC(which, d):
    """Oscillation frequency of a well

    Args:
        which (string): 'phi_left' or 'phi_right', indicating which well to
            analyze.
        d (np.ndarray): Structured data array with named columns.

    Returns (float):
        Oscillation frequency of the well divideded by harmonic limit frequency.
    """
    phi = d[which]
    beta = d['beta']
    return well_frequency_over_f_LC(phi, beta)


def dep_well_frequency_over_f_LC_approx(d):
    """Approximate well frequency with zero tilt flux.

    Args:
        d (np.ndarray): Structured data array with named columns.

    Returns (float):
        Well frequency divided by harmonic limit frequency.
    """
    return well_frequency_over_f_LC_approx(d['beta'])


def dep_coupling_factor_no_s(eta_c, N, d):
    """Coupling factor for two qubits, ignoring the self-inductance adjustment

    Args:
        eta_c (float): Coupling efficiency, i.e. the ratio M/L where M is the
            mutual inductance and L is the self inductance.
        N (int): Number of qubits you wish to couple to a single qubit.
        d (np.ndarray): Structured data array, with named columns.

    Returns (float):
        The dimensionless part of the qubit-qubit coupling. Multiply by omega_LC
        to get a coupling frequency.

    The qubit-qubit coupling is written
    J / hbar = [s (eta_c / N) (R_K / (8 pi Z_LC)) <0|phi|1>^2] omega_LC
    The part in [..] is the "coupling factor", a dimensionless parameter. The s
    factor is defined as (1 - M^2 / (L1 L2))^-1 and is near unity, so we ignore
    it.

    See the annealer-handbook for a full development of the coupling strength
    formula.
    """
    Z_LC = d['Z_LC'] * U.Unit(LABELS['Z_LC'][1])
    f_LC = d['f_LC'] * U.Unit(LABELS['f_LC'][1])

    phi_10 = d['phi_10']
    return (eta_c / N) * (R_K / (8 * np.pi * Z_LC)) * phi_10**2


# Data production. These functions actually analyze the Hamiltonian and create
# plots of the results. Each one groups several related plots together.

def frequency(f_LCs=FREQUENCIES, Z_LCs=IMPEDANCES):
    """Generate plots of relevant frequency scales for the flux qubit

    Args:
        f_LCs (ValueArray[Hz]): Harmonic limit frequencies.
        Z_LCs (ValueArray[Ohm]): Harmonic limit impedances.
    """

    phi, data = diagonalize_H_over_params(Z_LCs,
        np.linspace(0, 3.5, 50),
        f_LCs,
        VA([0], '')
        )

    # omega_10/omega_LC and J/omega_LC versus beta
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'phi_x']
    deps = [
        {
            'func': dep_omega_10_over_omega_LC,
            'label': r'$\omega_{10} / \omega_{LC}$'
        },
        {
            'func': functools.partial(dep_coupling_factor_no_s, 0.1, 10),
            'linestyle': '--', 'label': r'$J/\omega_{LC}$'
        }
    ]
    def indep(d):
        return d['beta']
    plot_iter(ax, data, constants, indep, deps, plot_type='semilogy')
    #ax.set_title(r"Normalized frequency and coupling versus $\beta$")
    ax.legend(loc='lower right', prop={'size': LEGEND_FONT_SIZE})
    make_presentable(ax, xlabel=r'$\beta$')
    ax.set_ylim(bottom=0.002, top=1)
    ax.set_xlim(left=0, right=2.0)

    # relative anharmonicity versus frequency
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'phi_x']
    def dep(d):
        return (d['f_21'] - d['f_10']) / d['f_10']
    deps = [{'func': dep}]
    indep = dep_omega_10_over_omega_LC
    plot_iter(ax, data, constants, indep, deps, plot_type='semilogy')
    ax.axhline(0.05, linewidth=3, color='k')
    ax.text(0.2, 0.03, r'5% anharmonicity', size=LEGEND_FONT_SIZE)
    ax.set_title("Relative anharmonicity versus frequency")
    make_presentable(ax, xlabel=r'$\omega_{10}/\omega_{LC}$',
        ylabel=r'Relative anharmonicity: $(\omega_{21} - \omega_{10})/\omega_{10}$')
    ax.legend(loc='upper right', prop={'size': LEGEND_FONT_SIZE})
    ax.set_ylim(bottom=0.001, top=10)
    ax.set_xlim(left=0.1, right=0.9)

    # Matrix element
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'phi_x']
    def dep(d):
        return d['phi_10']
    deps = [{'func': dep}]
    def indep(d):
        return d['beta']
    plot_iter(ax, data, constants, indep, deps, plot_type='semilogy')
    ax.set_title(r'Matrix element versus $\beta$')
    ax.legend(loc='upper left', prop={'size': LEGEND_FONT_SIZE})
    make_presentable(ax, xlabel=r'$\beta$',
        ylabel=r'$\langle 0 | \phi | 1 \rangle$')
    ax.set_xlim(left=0, right=2.0)

    # Sensitivity to phi_x
    phi, data = diagonalize_H_over_params(Z_LCs,
        [1.2],
        f_LCs,
        VA(np.linspace(0, 0.2), '')
    )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'beta']
    deps = [{'func': dep_omega_10_over_omega_LC}]
    plot_iter(ax, data, constants, lambda d,c: d['phi_x'], deps)
    ax.set_title(r'$\omega_{10} / \omega_{LC}$ versus $\phi_x$')
    ax.legend(loc='lower right', prop={'size': LEGEND_FONT_SIZE})
    make_presentable(ax, xlabel=r'$\phi_x$',
        ylabel=r'$\omega_{10}/\omega_{LC}$')


def measurement(f_LCs=FREQUENCIES, Z_LCs=IMPEDANCES):
    """Plots relevant to the measurement process

    Args:
        f_LCs (ValueArray[Hz]): Harmonic limit frequencies of the qubit.
        Z_LCs (ValueArray[Ohm]): Harmonic limit impedances of the qubit.
    """

    # Well frequencies versus beta at phi_x = 0
    betas = np.linspace(1.1, 2.5, 20)
    phi_xs = [0]
    phi, data = diagonalize_H_over_params(Z_LCs, betas, f_LCs, phi_xs)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'phi_x']
    deps = [{'func': functools.partial(dep_well_frequency_over_f_LC, 'phi_left')},
            {'func': dep_well_frequency_over_f_LC_approx,
             'label': 'approx', 'linestyle': '--', 'linewidth': 2}
    ]
    def indep(d):
        return d['beta']
    plot_iter(ax, data, constants, indep, deps)
    ax.set_title(r"Relative well frequency")
    make_presentable(ax, xlabel=r'$\beta$',
        ylabel=r'$\omega_{w} / \omega_{LC}$')
    ax.legend(loc='lower right', prop={'size': LEGEND_FONT_SIZE})

    # Tunnel time versus well frequency at phi_x = 0
    phi_xs = [0]
    phi, data = diagonalize_H_over_params(Z_LCs, betas, f_LCs, phi_xs)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'phi_x']
    deps = [{'func': functools.partial(dep_tunnelling_time_times_f_LC, 'phi_left'),
        'alpha': 0.4}]
    def indep(d):
        return dep_well_frequency_over_f_LC('phi_left', d)
    plot_iter(ax, data, constants, indep, deps, plot_type='semilogy')
    ax.set_title("Tunnel time versus well frequency")
    make_presentable(ax, xlabel=r'$f_w / f_{LC}$',
        ylabel=r'$T \times f_{LC}$')
    ax.legend(loc='upper left')

    # Data for varied phi_x
    betas = [2.0, 2.25, 2.5]
    phi_xs = np.linspace(0, 0.2, 20)
    del data
    phi, data = diagonalize_H_over_params(Z_LCs, betas, f_LCs, phi_xs)

    # Number of states versus beta and phi_x
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'beta']
    deps = [
        {
            'func': functools.partial(dep_number_of_states, 'phi_left'),
            'label': '',
            'alpha': 1.0
        },
        {
            'func': functools.partial(dep_number_of_states, 'phi_right'),
            'label': '',
            'alpha': 1.0
        }
    ]
    def indep(d):
        return d['phi_x']
    plot_iter(ax, data, constants, indep, deps)
    ax.set_title(r'Number of states versus $\phi_x$')
    make_presentable(ax, xlabel=r'$\phi_x$',
        ylabel=r'Number of states')
    ax.legend(loc='upper right', prop={'size': LEGEND_FONT_SIZE})

    # Well frequencies and well frequency differences versus phi_x
    def fit_curves(which, d):
        # Exact well frequency at phi_x = 0
        initial = well_frequency_over_f_LC(d[which][0], d['beta'])
        # Approximate values of phi_m and f_well / f_LC from e.g. John's writeup
        phi_m = potential_minima_approx(d['beta'])
        f_w_over_f_LC = well_frequency_over_f_LC_approx(d['beta'])
        # Slope computed from analytic (approximate) calculation
        slope = {'phi_left': -1, 'phi_right': 1}[which] * (phi_m / 2) * f_w_over_f_LC**-3
        return initial + d['phi_x'] * slope

    def frequency_diff(d):  # Frequency difference between wells
        f_left = well_frequency_over_f_LC(d['phi_left'], d['beta'])
        f_right = well_frequency_over_f_LC(d['phi_right'], d['beta'])
        return np.abs(f_left - f_right)

    # Well frequencies
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'beta']
    deps = [
        {
            'func': functools.partial(dep_well_frequency_over_f_LC, 'phi_left'),
            'label': '', 'alpha': 0.5
        },
        {
            'func': functools.partial(dep_well_frequency_over_f_LC, 'phi_right'),
            'label': '', 'alpha': 0.5
        },
        {
            'func': functools.partial(fit_curves, 'phi_left'),
            'linestyle': '--', 'color': 'k', 'linewidth': 2
        },
        {
            'func': functools.partial(fit_curves, 'phi_right'),
            'linestyle': '--', 'color': 'k', 'linewidth': 2
        }
    ]
    def indep(d):
        return d['phi_x']
    plot_iter(ax, data, constants, indep, deps)
    make_presentable(ax, xlabel=r'$\phi_x$', ylabel=r'$f_w / f_{LC}$')
    ax.set_title(r'Left and right well frequencies versus $\phi_x$')
    ax.legend(loc='lower left', prop={'size': LEGEND_FONT_SIZE})

    # Well frequency differences versus phi_x
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'beta']
    def indep(d):
        return d['phi_x']
    deps = [{'func': frequency_diff, 'alpha': 0.5, 'linestyle': '--'}]
    plot_iter(ax, data, constants, indep, deps)
    make_presentable(ax, xlabel=r'$\phi_x$',
        ylabel=r'$\delta \omega_w / \omega_{LC}$')
    ax.set_title(r'Well frequency difference versus $\phi_x$')
    ax.legend(loc='upper left', prop={'size': LEGEND_FONT_SIZE})

    # Number of states versus well frequency difference
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'beta']
    deps = [
        {'func': functools.partial(dep_number_of_states, 'phi_left')},
        {'func': functools.partial(dep_number_of_states, 'phi_right')}
    ]
    plot_iter(ax, data, constants, frequency_diff, deps)
    make_presentable(ax, xlabel=r'$\delta \omega_w / \omega_{LC}$',
        ylabel='Number of states')
    ax.legend(loc='lower left')

    # Tunneling time versus well frequency difference
    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'beta']
    deps = [
        {'func': functools.partial(dep_tunnelling_time_times_f_LC, 'phi_left')},
    ]
    plot_iter(ax, data, constants, frequency_diff, deps, plot_type='semilogy')
    make_presentable(ax, xlabel=r'$\delta \omega_w / \omega_{LC}$',
        ylabel=r'$T \times f_{LC}$')
    ax.legend(loc='lower left', prop={'size': LEGEND_FONT_SIZE})


def f_10_and_f_well_versus_beta():
    """Plot f_10 and well frequency against beta for phi_x=0"""
    Z_LCs = VA([IMPEDANCE['Ohm']], 'Ohm')
    f_LCs = VA([FREQUENCY['GHz']], 'GHz')
    betas = np.linspace(0.1, 1.8, 30)
    phi_xs = [0]
    phi, data = diagonalize_H_over_params(Z_LCs, betas, f_LCs, phi_xs)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    names = ['f_LC', 'Z_LC', 'phi_x']
    def well_freq(d):
        f_LC = d['f_LC'] * U.Unit(LABELS['f_LC'][1])
        Z_LC = d['Z_LC'] * U.Unit(LABELS['Z_LC'][1])
        betas = d['beta']
        f_well = f_LC * well_frequency_over_f_LC(d['phi_left'], betas)
        return f_well[LABELS['f_LC'][1]] / d['f_LC']
    def f_10(d):
        return d['f_10'] / d['f_LC']
    def indep(d):
        return d['beta']
    deps = [
        {'func': f_10, 'label': r'$\omega_{10} / \omega_{LC}$'},
        {'func': well_freq, 'label': r'$\omega_w / \omega_{LC}$', 'color': 'g'}
    ]
    plot_iter(ax, data, names, indep, deps, plot_type='semilogy')
    make_presentable(ax, xlabel=r'$\beta$', ylabel='Frequency',
        y_limits=(0.1, 1.4))
    ax.legend(loc='lower left', prop={'size': LEGEND_FONT_SIZE})
    ax.set_ylim


def freezing():
    """Time for qubit to tunnel out of metastable well"""
    betas = np.linspace(1.1, 2, 20)
    phi_xs = [0]
    phi, data = diagonalize_H_over_params(IMPEDANCES, betas, FREQUENCIES, phi_xs)

    # Plot tunnelling time versus beta
    fig = plt.figure()
    ax = fig.add_subplot(111)
    names = ['f_LC', 'Z_LC', 'phi_x']

    deps = [{'func': functools.partial(dep_tunnelling_time_times_f_LC, 'phi_left')}]
    def indep(d):
        return d['beta']
    plot_iter(ax, data, names, indep, deps, plot_type='semilogy')
    ax.set_title(r"Tunneling time")
    make_presentable(ax, xlabel=r'$\beta$', ylabel=r'$T \times f_{LC}$')
    ax.legend(loc='lower right', prop={'size':LEGEND_FONT_SIZE})

    # Tunnelling time versus frequency/f_LC
    fig = plt.figure()
    ax = fig.add_subplot(111)
    names = ['f_LC', 'Z_LC', 'phi_x']
    def indep(d):
        return d['f_10'] / d['f_LC']

    def dep(d):
        f_LC = d['f_LC'] * U.Unit(LABELS['f_LC'][1])
        t = dep_tunnelling_time_times_f_LC(d)
        return t

    deps = [{'func': functools.partial(dep_tunnelling_time_times_f_LC, 'phi_left')}]
    plot_iter(ax, data, names, indep, deps, plot_type='loglog')
    ax.set_title(r"Tunneling time versus $f_{10}$")
    make_presentable(ax, xlabel=r'$f_{10}/f_{LC}$', ylabel=r'Tunnel time $\times f_{LC}$')
    ax.legend(loc='upper right')


def bias_damping():
    """Bias line limitation on qubit T1"""
    betas = [1.0, 1.5]
    phi_xs = [0]
    phi, data = diagonalize_H_over_params([150*Ohm], betas, FREQUENCIES, phi_xs)
    M = VA(np.linspace(0.01, 2), 'pH')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    constants = ['f_LC', 'Z_LC', 'phi_x', 'beta']

    def T1_us(d, M):
        f_LC = d['f_LC'][0] * U.Unit(LABELS['f_LC'][1])
        Z_LC = d['Z_LC'][0] * U.Unit(LABELS['Z_LC'][1])
        f_10 = d['f_10'][0] * U.Unit(LABELS['f_LC'][1])
        phi_10 = d['phi_10'][0]
        R = 50 * Ohm
        Q = 4 * np.pi * (R / R_K) * (Z_LC / (2 * np.pi * f_LC * M))**2 / (phi_10**2)
        return (Q / (2 * np.pi * f_10))['us']

    deps = [{'func': T1_us}]

    def indep(d):
        return M['pH']

    plot_iter(ax, data, constants, indep, deps, plot_type='loglog', args=(M,))
    ax.set_title(r"$T_1$ from bias line versus $M$")
    make_presentable(ax, xlabel=r'$M$[pH]', ylabel=r'$T_1$ [us]')
    ax.legend(loc='upper left', prop={'size': LEGEND_FONT_SIZE})


def g_versus_kappa(Delta, delta_f_w, f_r):
    """Compute qubit-resonator coupling g needed versus resonator Q

    Args:
        Delta (Value[Hz]): Qubit-resonator detuning.
        delta_f_w (Value[Hz]): Frequency difference between flux qubit wells
            during readout phase.
        f_r (Value[Hz]): Resonator frequency (i.e. average of the two well
            frequencies).
    """
    #Convert to angular frequency
    Delta = 2*pi * Delta
    delta_omega_w = 2*np.pi * delta_f_w
    omega_r = 2*np.pi * f_r

    kappas = 1.0 / ValueArray(np.logspace(0, 3), 'us')
    Qs = omega_r / kappas
    gs = ((Delta**2 * kappas)/(2 * delta_omega_w))**0.5  # Angular frequency

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogx(Qs, gs['MHz']/(2*np.pi), linewidth=4)
    ax.grid()
    make_presentable(ax, xlabel=r'$Q_r$', ylabel=r'$g$ [MHz]')


class Inductor(object):
    """
    Represents an inductor and provides properties useful for qubit analysis.
    """

    def __init__(self, L):
        assert Value(1.0, 'H').isCompatible(L)
        self._L = L

    @staticmethod
    def L_to_E_L(L):
        """Convert inductance to energy scale"""
        return (PHI_0 / (2*pi))**2 / L

    @property
    def L(self):
        return self._L

    @property
    def E_L(self):
        """Inductor energy scale"""
        return self.L_to_E_L(self._L)

    @staticmethod
    def potential(L, phi, phi_ext):
        """Potential energy of an inductor"""
        return 0.5 * Inductor.L_to_E_L(L) * (phi - phi_ext)**2

    def potential_func(self):
        """Get a function representing this inductor's potential energy"""
        return functools.partial(Inductor.potential, self.L)


class Junction(object):
    """
    Represents a junction and provides properties useful for qubit analysis.
    """

    def __init__(self, I_c):
        assert Value(1.0, 'A').isCompatible(I_c)
        self._I_c = I_c

    @staticmethod
    def I_c_to_E_J(I_c):
        """Convert from critical current to junction energy scale."""
        return Phi_0 * I_c / (2*pi)

    @property
    def I_c(self):
        return self._I_c

    @property
    def E_J(self):
        return self.I_c_to_E_J(self.I_c)

    @staticmethod
    def potential(I_c, phi):
        return -1 * Junction.I_c_to_E_J(I_c) * np.cos(phi)

    def potential_func(self):
        """Get a function representing this junction's potential energy"""
        return functools.partial(Junction.potential, self.I_c)
