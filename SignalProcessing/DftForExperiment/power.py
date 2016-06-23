import numpy as np
from numpy import pi,cos,sin,exp

from matplotlib import pyplot as plt

plt.rcParams['font.size']=35

nus = np.linspace(-3,3,1000)
curve = 0.5*(1+sin(4*pi*nus)/(4*pi*nus))

plt.plot(nus,curve,linewidth=5)
plt.grid()
plt.xlabel('Xi')
plt.ylabel('Power')
plt.show()
