from InelasticAnalysisLib import RateVsEnergy
from matplotlib import pyplot as plt
import numpy as np


ER = np.linspace(0.1,500,2000)

n0 = RateVsEnergy(ER, 1000, 132, 54, 1, 1, 1e-40, 0., 0.3)
n100 = RateVsEnergy(ER, 1000, 132, 54, 1, 1, 1e-40, 100., 0.3)
n200 = RateVsEnergy(ER, 1000, 132, 54, 1, 1, 1e-40, 200., 0.3)
n300 = RateVsEnergy(ER, 1000, 132, 54, 1, 1, 1e-40, 300., 0.3)

#mask = ER>100
#n100[mask] = n100[mask]/2
#n200[mask] = n200[mask]/2
#n300[mask] = n300[mask]/2
#n0[mask] = n0[mask]/2

scalingFactor = 60 * 60 * 24

plt.plot(ER,n0   * scalingFactor)
plt.plot(ER,n100 * scalingFactor)
plt.plot(ER,n200 * scalingFactor)
plt.plot(ER,n300 * scalingFactor)

plt.figure(1)
plt.axis([0,500,1e-12,1])
plt.semilogy()
plt.show()

plt.figure(2)
n0 = RateVsEnergy(ER, 10000, 132, 54, 1, 1, 1e-40, 0., 0.3)
n100 = RateVsEnergy(ER, 10000, 132, 54, 1, 1, 1e-40, 100., 0.3)
n200 = RateVsEnergy(ER, 10000, 132, 54, 1, 1, 1e-40, 200., 0.3)
n300 = RateVsEnergy(ER, 10000, 132, 54, 1, 1, 1e-40, 300., 0.3)

plt.plot(ER,n0   * scalingFactor)
plt.plot(ER,n100 * scalingFactor)
plt.plot(ER,n200 * scalingFactor)
plt.plot(ER,n300 * scalingFactor)

plt.axis([0,500,1e-12,1])
plt.semilogy()
plt.show()


