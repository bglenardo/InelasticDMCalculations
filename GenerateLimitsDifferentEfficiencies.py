import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from InelasticAnalysisLib import CutAndCount90PercentConfXSec_Efficiency

delta = np.linspace(0,550,150)
mdens = 0.3 # GeV/cm^3/c^2
fn = 1.
fp = 1.
Mdm = 10000.

LUX_exp = 3.35e4
LUX_Emin = 1.
LUX_Emax = 400.
nevents = 0

xsec_LUX_2_100 = CutAndCount90PercentConfXSec_Efficiency(LUX_Emin, LUX_Emax, \
                 LUX_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents,'Efficiencies/vals_2_100.txt')
xsec_LUX_2_200 = CutAndCount90PercentConfXSec_Efficiency(LUX_Emin, LUX_Emax, \
                 LUX_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents,'Efficiencies/vals_2_200.txt')
xsec_LUX_2_300 = CutAndCount90PercentConfXSec_Efficiency(LUX_Emin, LUX_Emax, \
                 LUX_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents,'Efficiencies/vals_2_300.txt')
xsec_LUX_2_350 = CutAndCount90PercentConfXSec_Efficiency(LUX_Emin, LUX_Emax, \
                 LUX_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents,'Efficiencies/vals_2_350.txt')
xsec_LUX_2_400 = CutAndCount90PercentConfXSec_Efficiency(LUX_Emin, LUX_Emax, \
                 LUX_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents,'Efficiencies/vals_2_400.txt')

plt.figure(1)
plt.semilogy()

plt.plot(delta,xsec_LUX_2_100,'-',linewidth=2,color='Black',label='S1c = 2-100')
plt.plot(delta,xsec_LUX_2_200,'-',linewidth=2,color='DarkBlue',label='S1c = 2-200')
plt.plot(delta,xsec_LUX_2_300,'-',linewidth=2,color='Blue',label='S1c = 2-300')
plt.plot(delta,xsec_LUX_2_350,'-',linewidth=2,color='CadetBlue',label='S1c = 2-350')
plt.plot(delta,xsec_LUX_2_400,'-',linewidth=2,color='MediumSeaGreen',label='S1c = 2-400')

plt.xlabel('Mass splitting (keV)');
plt.ylabel('Cross section (cm^2)');
plt.legend(loc='upper left')
plt.axis([0.,400.,1e-45,1e-36])
plt.show()
plt.savefig('InelasticDMProjection_m_10TeV_differentEfficiencies.pdf')


