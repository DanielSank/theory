from __future__ import division

# This file processes SPICE simulations associated with a purcell filter
# for fast qubit readout. Only data analysis is done here.
# LtSpice.litspice and LtSpice.rawfile invoke SPICE and extract data from the
# files written by SPICE

import numpy as np
import re
import matplotlib.pyplot as plt
import scipy.signal
import LtSpice.ltspice as spice

from LtSpice.rawfile import load_spice_raw, spicefile_get_columns

#f_res = 6.75e9
f_q = 6.0e9
C_q = 85e-15
#C_res = np.pi / (4*2*np.pi*f_res*50)
eta = 240e6
h = 6.626e-34


# Left over parsing functions that haven't been moved into LtSpice.rawfile

def pair2cplx(s):
    a,b = s.split(',')
    return float(a)+1j*float(b)

def parse_dBangle(s):
    match_float = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
    rexp = r'\((%s)dB,(%s).*\)' % (match_float, match_float)
    mo = re.match(rexp, s)

# Analysis
    
def fwhm(x, threshold=0.5):
    z = x[:,1] + 1j*x[:,2]
    f = x[:,0]
    half_max = np.max(np.abs(z))*threshold
    mask = np.abs(z)<half_max
    tmp = np.diff(mask)
    return np.diff(f[np.nonzero(tmp)])

def coupling_cap(g, fResonator, CResonator, fQubit):
    '''Coupling capacitance between qubit and resonator
    
    g: qubit-resonator coupling (Frequency dimensions)
    fResonator: Frequency of resonator
    fQubit: frequency of Qubit
    CResonator: Equivalent capacitance of resonator
    
    Returns coupling capacitance in SI units
    '''
    C_g = 2*g * np.sqrt(C_q*CResonator) / (np.sqrt(fQubit*fResonator))
    return C_g

def g_eff(f, Cc):
    '''
    Calculate g for a given coupling cap as a function of qubit frequency.
    '''
    g = np.sqrt(f_res*f)*Cc/np.sqrt(C_q*C_res)
    return g

def PurcellSim(x, g0):
    Cc = coupling_cap(g0)
    f = x[:,0]
    z = x[:,1] + 1j*x[:,2] + 1/(1j*2*np.pi * f * Cc)
    Qex = abs(z.imag / z.real)
    Q = Qex * C_q/Cc
    T1 = Q / (2*np.pi*f)
    return f,T1,Q,Qex

def PurcellFormula(f, g0, kappa):
    '''
    The simple and incorrect (g/delta)^2 kappa model for the purcell effect.  This is
    only valid in for a single mode, which is not a good approximation for us.
    '''
    Cc = coupling_cap(g0)
    g = g_eff(f, Cc)
    delta = f-f_res
    damping = (g/delta)**2 * kappa
    T1 = 1/damping
    return T1

def analyze_purcell(fname, g_nom, fResonator, CResonator, fQubit_nom):
    '''Estimate the purcell losses as a function of qubit frequency.
    
    fname - str: Name of .raw file containing result of SPICE simulation.
    g_nom - real: Value of resonator-qubit coupling. Because g is frequency
                  but coupling capacitance is fixed we take a "nominal" g here
                  along with a nominal qubit frequency to compute the fixed
                  coupling capacitance.
    fResonator - real: Frequency of readout resonator
    CResonator - real: Equivalent capacitance of readout resonator
    fQubit_nom - real: Nominal qubit frequency. This is the frequency of the
                       qubit corresponding to g_nom
    
    To generate an appropraite .raw file run a small-signal AC simulation with 
    the input resonator drive off, but with an amplitude '1' current source
    where the qubit coupling cap would go.
    '''
    data = load_spice_raw(fname)
    C_g = coupling_cap(g_nom, fResonator, CResonator, fQubit_nom)
    x = spicefile_get_columns(data, 'V(vres)')
    f = x[:,0].real
    Z_g = 1/(1j*2*np.pi*f*C_g)
    z = x[:,1] + Z_g
    Qex = abs(z.imag/z.real)
    Q = Qex * C_q/C_g #C_q is globally defined qubit capacitance
    T1 = Q / (2*np.pi*f)
    #T1_simple = PurcellFormula(f, g0, kappa)
    return f, T1
    
def analyze_purcell_v2(fname):
    data = load_spice_raw(fname)
    x = spicefile_get_columns(data, ['V(vqubit)','I(v2)'])
    Y = x[:,2] / x[:,1]
    parallelResistance = 1.0 / np.abs(Y.real)
    f = x[:,0]
    T1 = C_q * parallelResistance
    return f, T1
    
def analyze_readout(fname, g0=None, f_q=f_q, plot=False):
    '''
    Analyzes readout parameters from a spice simulation.

    Run a small signal AC analysis with the drive voltage source
    on, and sweep over the frequency range that includes the resonator.  If g0 is not given, 
    it will be calculated assuming 2*chi = 2*kappa.
    '''
    data = load_spice_raw(fname)
    x = spicefile_get_columns(data, ['V(vres)', 'V(out)'])
    f =x[:,0].real
    vres = x[:,1]
    vout = x[:,2]
    #Find resonator resonance frequency
    idx = np.argmax(np.abs(vres))
    fres = f[idx]
    #Find FWHM of resonator
    mask = np.abs(vres) > np.max(np.abs(vres))/np.sqrt(2) # Find the 1/sqrt(2) amplitude points
    fa = f[np.min(np.nonzero(mask))]
    fb = f[np.max(np.nonzero(mask))]
    kappa = (fb-fa)/2.0
    print "fres: %g kappa_a: %g kappa_b: %g kappa: %g" % (fres, fres-fa, fb-fres, kappa)
    delta = fres - f_q
    #If no g0 is provided, use \chi=\kappa to compute g_0
    #\chi = (g^2/\Delta)(1/(1+\Delta/\eta))
    if g0 is None:
        g0 = np.sqrt(delta*kappa*(1+delta/eta))
        print "optimal g0: %g" % g0
    ncrit = (delta/g0)**2
    chi = (g0**2 / delta)/(1+(delta/eta))
    nmax = ncrit/2
    print "g0: %f, delta: %s, ncrit: %s, chi: %g" % (g0, delta, ncrit, chi)
    f1 = fres - chi
    f2 = fres + chi
    print "f1: %s, f2: %s" % (f1, f2)
    v_out1 = vout[np.argmin(np.abs(f-f1))]
    v_out2 = vout[np.argmin(np.abs(f-f2))]
    v_res1 = vres[np.argmin(np.abs(f-f1))]
    v_res2 = vres[np.argmin(np.abs(f-f2))]
    print "v_out1: %s, v_out2: %s, v_res1: %s, v_res2: %s" % (v_out1, v_out2, v_res1, v_res2)
    vmax = np.sqrt(4*h*f_res**2*50*nmax)
    v_drive = vmax / max([abs(v_res1), abs(v_res2)])
    v_amp = max([abs(v_out1), abs(v_out2)])*v_drive
    print "Pin: %f dBm Pamp: %f dBm" % (10*np.log10(1000 * v_drive**2/50), 
                                        10*np.log10(1000 * v_amp**2/50))
    print "Pamp %f photon/ns" % ( 1e-9*(v_amp**2/50) / (h*fres))
    separation = abs(v_drive * (v_out1 - v_out2))
    print "separation: %f photon/ns" % (1e-9*separation**2/50/(h*fres),)
    t, y1, demod1 = TF2sine_response(f, vout, f1)
    _, y2, demod2 = TF2sine_response(f, vout, f2)
    tmask = np.logical_and(t>0, t<1000e-9)
    if plot:
        plt.figure()
        plt.plot(t[tmask], demod1[tmask].real, label="f1 real")
        plt.plot(t[tmask], demod1[tmask].imag, label="f1 imag")
        plt.plot(t[tmask], demod2[tmask].real, label="f2 real")
        plt.plot(t[tmask], demod2[tmask].imag, label="f2 imag")
        plt.legend()
        plt.grid()
        plt.figure()
        plt.plot(t[tmask], abs(demod1[tmask]-demod2[tmask]), label="separation")
        plt.grid()
    return {'g0':g0, 'ncrit':ncrit}
        
def TF2sine_response(f, tf, f0):
    '''
    Take the frequency space response function and calculate the time domain response to a
    sine wave that turns on at time zero.

    For convenience, we also return a demodulated signal to show the envelope ring-up time.
    '''
    df = f[1]-f[0]
    prepad = int(min(f) // df)
    npts = prepad + len(f)
    postpad = 2**int(np.ceil(np.log2(npts)))-npts + 1
    npts = prepad + postpad + len(f)
    sr = 2*(df * npts)
    tmp = np.zeros(npts, dtype=complex)
    tmp[prepad:prepad+len(f)] = tf
    tf_f = tmp
    
    tmp = np.arange(npts-1)
    x = np.sin(2.0*np.pi*f0*tmp/sr)
    x = np.hstack((0*x, x))
    x_f = np.fft.rfft(x)
    y_f = x_f * tf_f
    y_t = np.fft.irfft(y_f)
    max_t = 1/(2.0*df)
    t = np.linspace(-max_t, max_t, num=len(y_t))
    ref = np.exp(-1j*2*np.pi*f0*tmp/sr)
    ref = np.hstack((0*ref, ref))
    demod = ref*y_t
    filter_pts = int(5e-9*sr) # 5 nanosecond
    filter_win = scipy.signal.kaiser(filter_pts, 14)
    demod = np.convolve(demod, filter_win, mode='same')
    return t,y_t, demod

def T1_theory(dx, f_q, C_q, C_g, Qf, Qrr, order='all'):
    Zq = 1.0/(2*np.pi*f_q*C_q)
    f_r = f_q / (1 + dx)
    if order == 'all':
        Qq = (np.pi/4) * Qf**2 * Qrr * (C_q/C_g)**2 * \
             (Zq/50.0) * ((2*dx + dx**2)/(1+dx))**4
    else:
        raise RuntimeError("order must be 'all'")
    return   (1.0/(2*np.pi*f_q))*Qq
    
# User functions to analyze the filtered measurement system

### BEGIN bad documentation region - functions may have ERRORS!

def purcell_vs_Qrr(netlistFilename, outFilename, Qrrs, sweepValues,
                   f_res=6.75E9):
    colors=['b','g','y','r','k']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for Qrr, color in zip(Qrrs,colors):
        delta = np.abs(f_q-f_res)
        kappa = f_res/(2*Qrr)
        g0 = np.sqrt(delta*kappa*(1+delta/eta)) #For matched kappa=chi
        print g0
        #Analyze purcell effect
        userParameters = {'f0': f_res, 'qc':Qrr, 'v_readout':0, 'i_qubit':1}
        allParameters = spice.processWithParameters(netlistFilename,
            outFilename, userParameters, sweepValues)
        #f_res = allParameters['params']['f0']
        C_res = allParameters['params']['cres']
        Ck = allParameters['params']['c_coupling']
        f, T1 = analyze_purcell(outFilename+'.raw', g0, f_res, C_res, f_q)
        _delta = f-f_res
        #Make the theory line
        Qf = allParameters['params']['qbp']
        Qrr = allParameters['params']['qc']
        C_g = coupling_cap(g0, f_res, C_res, f_q)
        dx = (f-f_res)/f_res
        T1_thy = T1_theory(dx, f_q, C_q, C_g, Qf, Qrr)
        legendLabel = 'Qrr = %d, Cg=%.2g, Ck=%.2g, g0=%.2g'%(Qrr, C_g, Ck, g0)
        ax.semilogy((_delta)/1E9, T1*1E6, color, linewidth=5, alpha=0.5, label=legendLabel)
        ax.semilogy(_delta/1E9, T1_thy*1E6, color, linewidth=2, label='Theory')
    plt.legend(loc='lower left', numpoints=1)
    plt.grid()
    plt.title(r'Purcell $T_1$ limit vs. $\Delta$ for various $Q_{rr}, Q_F=%d$'%Qf)
    plt.xlabel(r'Detuning $\Delta$ [GHz]')
    plt.ylabel('Qubit T1 [us]')

def purcell_vs_frrDetuning(netlistFilename, outFilename, detunings, sweepValues,
                           f_res=6.75E9, Qrr=1000):
    raise RuntimeError("Several lines in this function look wrong")
    colors=['b','g','y','r','k']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for detuning, color in zip(detunings,colors):
        delta = np.abs(f_q - f_res + detuning)
        chi = f_res/(2*Qrr)
        g0 = np.sqrt(delta*chi*(1+delta/eta))
        print g0
        #Analyze purcell effect
        userParameters = {'qc':Qrr, 'f0':f_res+detuning, 'v_readout':0, 'i_qubit':1}
        allParameters = spice.processWithParameters(netlistFilename, outFilename, userParameters, sweepValues)
        #f_res = allParameters['params']['f0']
        C_res = allParameters['params']['cres']
        Ck = allParameters['params']['c_coupling']
        f, T1 = analyze_purcell(outFilename+'.raw', g0, f_res, C_res, f_q)
        _delta = f-(f_res+detuning)
        Qf = allParameters['params']['qbp']
        Qrr = allParameters['params']['qc']
        C_g = coupling_cap(g0, f_res, C_res, f_q)
        legendLabel = 'resonator error = %.3g MHz'%(detuning/1E6)
        ax.semilogy((_delta[0::50])/1E9, T1[0::50]*1E6, color+'.', markersize=8, label=legendLabel)
    plt.legend(loc='lower left', numpoints=1)
    plt.grid()
    plt.title(r'Purcell $T_1$ limit vs. readout resonator frequency error for $Q_F=%d$, $Q_{rr}=%d$'%(Qf,Qrr))
    plt.xlabel(r'Detuning $\Delta$ [GHz]')
    plt.ylabel('Qubit T1 [us]')
    
def purcell_vs_filterQ(netlistFilename, outFilename, Qbps, sweepValues,
                       f_res=6.75E9, Qrr=1500):
    colors=['b','g','y','r','k']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for Qbp,color in zip(Qbps,colors):
        delta = np.abs(f_q-f_res)
        kappa = f_res/(2*Qrr)
        g0 = np.sqrt(delta*kappa*(1+delta/eta))
        userParameters = {'f0': f_res, 'qbp':Qbp, 'v_readout':0, 'i_qubit':1, 'qc':Qrr}
        #Run with the current parameters
        allParameters = spice.processWithParameters(netlistFilename, outFilename, userParameters, sweepValues)
        f_res = allParameters['params']['f0']
        C_res = allParameters['params']['cres']
        Ck = allParameters['params']['c_coupling']
        Qf = allParameters['params']['qbp']
        Qrr = allParameters['params']['qc']
        f, T1 = analyze_purcell(outFilename+'.raw', g0, f_res, C_res, f_q)
        delta = f-f_res
        #Make the theory line
        C_g = coupling_cap(g0, f_res, C_res, f_q)
        dx = (f-f_res)/f_res
        T1_thy = T1_theory(dx, f_q, C_q, C_g, Qf, Qrr)
        #Plot
        legendLabel = 'Qf = %d, Cg=%.2g, Ck=%.2g, g0=%.2g'%(Qf, C_g, Ck, g0)
        ax.semilogy((delta[0::50])/1E9, T1[0::50]*1E6, color+'.', markersize=8, label='$Q_F$ = '+str(Qbp))
        ax.semilogy(delta/1E9, T1_thy*1E6, color, linewidth=1, label='Theory')
    plt.legend(loc='lower left', numpoints=1)
    plt.grid()
    plt.title(r'Purcell $T_1$ limit vs. $\Delta$ for various $Q_F, Q_{rr}=%d$'%Qrr)
    plt.xlabel('Detuning [GHz]')
    plt.ylabel('Qubit T1 [us]')
    
### END bad documentation region
    
def purcell_vs_Delta(netlistFilename, outFilename, QFilter,
                     QResonators, sweepValues, fResonator=6.75E9):
    """Compute purcell limited T1 for a range of detunings resonator Qs
    
    Written 14 November 2013 - Daniel Sank
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colors = ['b','g','y','r','k']
    for QResonator, color in zip(QResonators, colors):
        #Choose C_g
        #Coupling is frequency dependent, but coupling capacitance is fixed.
        #We use a nominal qubit frequency of f_q (global variable) and whatever
        #fResonator is passed into this function. This gives a nominal Delta.
        #With the nominal Delta and a kappa computed from the desired
        #QResonator we find a nominal value of g from which we can find C_g.
        Delta_nom = np.abs(f_q - fResonator)
        kappa = 2*np.pi*fResonator / QResonator
        #Assume matched kappa=chi
        chi_nom = 0.5 * (kappa/(2.0*np.pi))
        g_nom = np.sqrt(Delta_nom * chi_nom * (1 + Delta_nom/eta))
        userParameters = {'f0':fResonator, 'qc': QResonator, 'qbp':QFilter,
                          'v_readout': 0, 'i_qubit':1}
        allParameters = spice.processWithParameters(netlistFilename,
            outFilename, userParameters, sweepValues)
        CResonator = allParameters['params']['cres']
        CKappa = allParameters['params']['c_coupling']
        #Compute T1 vs frequency. We use globally nominal qubit frequency where
        #analyze_purcell wants a nominal qubit frequency. This is ok because we
        #used the same nominal frequency when we computed g_nom.
        f, T1 = analyze_purcell(outFilename+'.raw', g_nom, fResonator,
                                CResonator, f_q)
        Delta = f - fResonator
        #Make theory line
        assert allParameters['params']['qbp'] == QFilter
        assert allParameters['params']['qc'] == QResonator
        C_g = coupling_cap(g_nom, fResonator, CResonator, f_q)
        dx = (f - fResonator) / fResonator
        T1_thy_all = T1_theory(dx, f_q, C_q, C_g, QFilter, QResonator,
                               order='all')
        #Add line to plot
        #legendLabel = 'SPICE - Qr = %d, Cg=%.2g, Ck=%.2g, g_nom=%.2g'%\
        #    (QResonator, C_g, CKappa, g_nom)
        legendLabel = r'SPICE, $Q_r =$ %d'%QResonator
        ax.semilogy((Delta)/1E9, T1*1E6, color+'-', linewidth=5, alpha=0.5,
                    label=legendLabel)
        #ax.semilogy(Delta/1E9, T1_thy_1*1E6, color+'-', linewidth=2,
        #            label='Theory - first order')
        ax.semilogy(Delta/1E9, T1_thy_all*1E6, color+'--', linewidth=2,
                    label='Theory')
    plt.legend(loc='lower left', numpoints=1)
    plt.grid()
    plt.title(r'Purcell $T_1$ limit vs. $\Delta$ for various $Q_{r}, Q_F=%d$'\
        %QFilter)
    plt.xlabel(r'Detuning $\Delta$ [GHz]')
    plt.ylabel('Qubit T1 [us]')
    return ax

def purcell_vs_Delta_v2(netlistFilename, outFilename, QFilter,
                        QResonators, sweepValues, fResonator=6.75E9,
                        ax=None):
    """Compute purcell limited T1 for a range of detunings resonator Qs
    
    Here we use a much simplified analysis with the spice file
    bpfilter_pi4_qubitT1.asc
    
    Written 3 December 2013 - Daniel Sank
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    colors = ['b','g','y','r','k']
    for QResonator, color in zip(QResonators, colors):
        #Choose C_g
        #Coupling is frequency dependent, but coupling capacitance is fixed.
        #We use a nominal qubit frequency of f_q (global variable) and whatever
        #fResonator is passed into this function. This gives a nominal Delta.
        #With the nominal Delta and a kappa computed from the desired
        #QResonator we find a nominal value of g from which we can find C_g.
        Delta_nom = np.abs(f_q - fResonator)
        kappa = 2*np.pi*fResonator / QResonator
        #Assume matched kappa=chi
        chi_nom = 0.5 * (kappa/(2*np.pi))
        g_nom = np.sqrt(Delta_nom * chi_nom * (1 + Delta_nom/eta))
        CResonator = 1.0 / (8 * fResonator * 50)
        C_g = coupling_cap(g_nom, fResonator, CResonator, f_q)
        userParameters = {'f0':fResonator, 'qc': QResonator, 'qbp':QFilter,
                          'c_g':C_g,
                          'v_readout': 0, 'i_qubit':0, 'v_qubit':1}
        allParameters = spice.processWithParameters(netlistFilename,
            outFilename, userParameters, sweepValues)
        CKappa = allParameters['params']['c_coupling']
        #Compute T1 vs frequency. We use globally nominal qubit frequency where
        #analyze_purcell wants a nominal qubit frequency. This is ok because we
        #used the same nominal frequency when we computed g_nom.
        f, T1 = analyze_purcell_v2(outFilename+'.raw')
        Delta = f - fResonator
        #Make theory line
        assert allParameters['params']['cres'] == CResonator
        assert allParameters['params']['qbp'] == QFilter
        assert allParameters['params']['qc'] == QResonator
        dx = (f - fResonator) / fResonator
        T1_thy_all = T1_theory(dx, f_q, C_q, C_g, QFilter, QResonator,
                               order='all')
        #Add line to plot
        riseTime = (1.0/kappa)*1E9
        legendLabel = r'SPICE, $\kappa_r^{-1} =$ %d ns'%(riseTime,)
        ax.semilogy((Delta)/1E9, T1*1E6, color+'-', linewidth=5, alpha=0.5,
                    label=legendLabel)
        ax.semilogy(Delta/1E9, T1_thy_all*1E6, color+'--', linewidth=2,
                    label='Theory')
    plt.legend(loc='lower left', numpoints=1, prop={'size':22})
    plt.grid()
    plt.title(r'Purcell $T_1$ limit vs. $\Delta$ for various $Q_{r}, Q_F=%d$'\
        %QFilter)
    plt.xlabel(r'Detuning $\Delta$ [GHz]')
    plt.ylabel('Qubit T1 [us]')
    return ax
    
def resonatorQ_vs_detuningFromFilter(netlistFilename, outFilename,
    fFilter, QF, detunings, Qr_target, sweepValues):
    """Find Q_r vs (f_F - f_r)
    
    If the measurement resonator is not in the center of the filter passband,
    we expect its quality factor Q_r to be higher than if it were in the
    center. We quantify this here.
    
    We plot Q_r against frequency for resonators whoe resonance frequency f_r
    is detuned from the "ideal" point at the filter's resonance f_F. What we
    call "Q_r" in the plot is literally
    
    Energy in the measurement resonator / Energy lost per radian 
    
    which is the definition of Q_r. Note that this ratio can be computed at any
    frequency. What we actually call "Q" is this ratio _at the resonance
    frequency of the mode under consideration_, in this case f_r. So, to read
    the resulting plot, you want to look at Q_r at the frequency of the detuned
    measurement resonator.
    
    fFilter: design resonance frequency of filter
    QF: design Q of filter
    detunings - iterable of Value: range of frequency differences between
                filter and measurement resontaor
    Qr_target: design Q of measurement resonator
    """
    colors = ['b', 'g', 'r', 'y', 'k']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for detuning, c in zip(detunings, colors):
        f0 = fFilter + detuning
        userParameters = {'f0': f0['Hz'], 'qbp':QF, 'qc':Qr_target,
                          'fbp':fFilter['Hz'],
                          'v_readout': 0, 'i_qubit':1}
        allParameters = spice.processWithParameters(netlistFilename,
            outFilename, userParameters, sweepValues)
        data = load_spice_raw(outFilename+'.raw')
        data = spicefile_get_columns(data, ['V(vres)', 'V(out)'])
        f = data[:,0].real
        #Q is defined as energy stored / energy lost per radian
        #For a \lambda/4 resonator this comes out to
        # Q = (pi/4)(R/Z_0) (Vres/VR)^2
        #where R is the resistor damping the system, Z_0 is the line impedance
        #of the the resonator, Vres is the voltage on the voltage amplitude at
        #the antinode of the resonator, and VR is the voltage amplitude on the
        #resistor.
        Q = (np.pi/4)*(50.0/50.0)*np.abs((data[:,1]/data[:,2]))**2
        ax.semilogy(f/1E9 - fFilter['GHz'], Q, c+'-', linewidth=4,
                    label=str(detuning['GHz'])+' GHz', alpha=0.8)
    ax.grid(which='major')
    ax.grid(which='minor')
    ax.legend(loc='lower left')
    ax.set_xlabel('Detuning [GHz]')
    ax.set_ylabel(r'$Q_r$')
    return ax