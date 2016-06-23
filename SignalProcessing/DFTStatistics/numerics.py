import numpy as np
import matplotlib.pyplot as plt

import labrad.units as units

ns, MHz = (units.Unit(s) for s in ['ns', 'MHz'])

ADC_SAMPLE_TIME = 2*ns

def gaussianWithSigma(sigma):
    """Get a Gaussian function with a specified sigma"""
    def f(n, m):
        return np.exp(-(n-m)**2 / float((2 * sigma**2)))
    return f

def covariance(corrFunc, N, i, q):
    """
    Compute the covariance matrix for a specific correlation function.
    
    i and q choose real and imaginary quadratures.
    """
    n, m, k, l = np.indices((N,N,N,N))
    funcs = {'i': np.cos, 'q': np.sin}
    mat = corrFunc(n,m)*funcs[i](2*np.pi*n*k/N)*\
                        funcs[q](2*np.pi*m*l/N)
    return np.sum(np.sum(mat, axis=0), axis=0)/float(N**2)
    
def iqCovarance(N, sigma):
    frequencies = np.arange(0, N)
    corr = gaussianWithSigma(sigma)
    datRR = covariance(corr, N, 'r', 'r')
    datRI = covariance(corr, N, 'r', 'i')
    datII = covariance(corr, N, 'i', 'i')
    resultRR = np.diag(datRR, 0)
    resultRI = np.diag(datRI, 0)
    resultII = np.diag(datII, 0)
    return resultRR, resultRI, resultII

def sigmas(N, sigmas):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for sigma in sigmas:
        corr = gaussianWithSigma(sigma)
        result = covariance(corr, N, 'i', 'i')
        ax.semilogy(np.arange(N)[1:], np.diag(result, 0)[1:],
                    '-', linewidth=4, label=str(sigma))
    ax.legend(loc = 'lower right')
    ax.grid()
    
def plot(resultRR, resultRI, resultII):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(frequencies[1:], resultRR[1:], 'b.', markersize=8, label="RR")
    ax.semilogy(frequencies[1:], resultRI[1:], 'r.', markersize=8, label="RI")
    ax.semilogy(frequencies[1:], resultII[1:], 'g.', markersize=8, label="II")
    ax.legend(loc='upper right')
    ax.grid()

### Asynchronous demodulation case
"""
In the lab we compute the frequency response of the system with an asynchronous
method. We do not actually compute a DFT. We fix a time T (or number of
samples) AND a frequency f and then compute the sum
sum_{n=0}^{N-1} cos(2 pi f t_n) mySignal(t_n)
where t_n = T * n/N.
"""

def covarianceAsyncDemod(N, quad1, quad2, freqs, corrFunc):
    """
    N - int: Number of sample points.
    quad1 - 'i' or 'q': Chooses cosine or sine.
    quad2 - 'i' or 'q': Chooses cosine or sine.
    freqs - iterable: demodulation frequencies, in sample units.
    sigma - scalar: width of correlation function in sample units.
    corrFunc - function: Time correlation function, in sample units.
    """
    funcs = {'i': np.cos, 'q': np.sin}
    n, m = np.indices((N, N))
    result = []
    #Loop over frequencies of interest. Not sure how to vectorize this.
    for freq in freqs:
        mat = corrFunc(n, m)*funcs[quad1](2*np.pi*n*freq)*\
                             funcs[quad2](2*np.pi*m*freq)
        result.append(np.sum(np.sum(mat, axis=0), axis=0)/float(N**2))
    return np.array(result)

def covarianceVsSigma(N, quad1, quad2, freqs, rolloffFrequencies,
                      adcSampleTime):
    """
    N - int: Number of sample points.
    freqs - iterable of Value[Hz]: Set of demodulation frequencies with
            physical units.
    rolloffFrequencies - iterable of Value[Hz]: Gaussian filter 3dB
                         frequencies in physical units. By 3dB frequency I mean
                         the frequency at which the transmitted voltage
                         amplitude is down by 3dB.
    adcSampleTime - Value[s]: Time between ADC samples, in physical units.
    """
    #Convert frequencies to sample units
    freqs_samples = [f['GHz']*adcSampleTime['ns'] for f in freqs]
    #Convert rolloff frequencies to sigma in sample units
    sigmas_freq = [f / (2*np.sqrt(2 * np.log(2))) for f in rolloffFrequencies]
    sigmas_ns = [1.0/(2*np.pi*s['GHz']) for s in sigmas_freq]
    sigmas_samples = [s/adcSampleTime['ns'] for s in sigmas_ns]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for sigma, rof in zip(sigmas_samples, rolloffFrequencies):
        corrFunc = gaussianWithSigma(sigma)
        result = covarianceAsyncDemod(N, quad1, quad2, freqs_samples, corrFunc)
        ax.semilogy([f['MHz'] for f in freqs][1:], result[1:]/result[1],
                    '-', linewidth=4, label=str(rof))
    
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Noise variancs [au]')
    ax.grid()
    ax.legend(loc='lower left')
    return ax