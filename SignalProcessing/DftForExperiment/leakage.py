# We plot the Discrete Fourier Transform of a sinusoid that makes 12.2
# oscillations on the aquisition window. We also plot an analytic formula.

# Note that the largest bin is number 12, which is the nearest to the actual
# frequency of 12.2.

import numpy as np
from matplotlib import pyplot as plt

from numpy import pi,sin,linspace

MARKERSIZE=30
plt.rcParams['font.size']=30

N=20
xi = 12.2
I=0+1.0j

n = np.linspace(0,N-1,N)
ks = n

theory = ((sin(pi*(xi-ks))/sin(pi*(xi-ks)/N))**2)/N**2

signal = np.exp(I*2*pi*xi*n/N)
ft = np.fft.fft(signal)
noncommensurateS = ft*np.conjugate(ft)/N**2


plt.semilogy(ks,theory,'r',linewidth=5,label='Expected for xi=12.2')
plt.semilogy(ks,noncommensurateS,'.',markersize=MARKERSIZE, label='Numerical for xi=12.2')
plt.xlabel('Bin number')
plt.ylabel('Normalized Power')
plt.grid()
plt.show()
