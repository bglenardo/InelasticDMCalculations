'''
This function interpolates the efficiency curves so that they
have the same dimensions as the RateVsEnergy arrays. 

The efficiency calculations were done using the high-energy
libNEST model developed for high-E NR LUX analyses.
'''

import numpy as np

def InterpolateEfficiency(eMin, eMax, nPoints, eff_txt_file):

  Energy_vec = np.linspace(eMin,eMax,nPoints)
  Eff_vec = np.zeros(nPoints)

  eff = np.genfromtxt(eff_txt_file,delimiter=',')

  for i in range(0,nPoints):
      if eff[-2,0] <= Energy_vec[i]:
          Eff_vec[i] = 0.
          continue
      mask = eff[:,0] < Energy_vec[i]
#      print(mask)
#      print(Energy_vec[i])
#      print(eff[:,0][mask][-1])
#      print(eff[:,0][np.logical_not(mask)][0])

      slope = (eff[:,1][np.logical_not(mask)][0] - eff[:,1][mask][-1]) / \
              (eff[:,0][np.logical_not(mask)][0] - eff[:,0][mask][-1])
      Eff_vec[i] = eff[:,1][mask][-1] + slope*(Energy_vec[i] - eff[:,0][mask][-1])

  return Eff_vec


