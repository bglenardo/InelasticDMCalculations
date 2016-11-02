import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from InelasticAnalysisLib import CutAndCount90PercentConfXSec

delta = np.linspace(0,550,150)
mdens = 0.3 # GeV/cm^3/c^2
fn = 1.
fp = 1.
Mdm = 10000.

LUX_exp = 1.4e4
LUX_Emin = 1.
LUX_Emax = 30.
LUX_Emax_best = 500.
LUX_Emax_NEST = 256.
nevents = 0

PICO_exp = 1.3e3
PICO_Emin = 7
PICO_Emax = 1000
nevents_PICO = 0

CRESST_exp = 0.052e3
CRESST_Emin = 30
CRESST_Emax = 120
nevents_CRESST = 4

XENON100_exp = 100.9*48
XENON100_Emin = 8.4
XENON100_Emax = 44.6
nevents_XENON100 = 3


xsec_LUX_1_30 = CutAndCount90PercentConfXSec(LUX_Emin, LUX_Emax, \
                 LUX_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents)
xsec_LUX_best = CutAndCount90PercentConfXSec(LUX_Emin, LUX_Emax_best, \
                 LUX_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents)
xsec_LUX_NEST = CutAndCount90PercentConfXSec(LUX_Emin, LUX_Emax_NEST, \
                 LUX_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents)

xsec_PICO = CutAndCount90PercentConfXSec(PICO_Emin, PICO_Emax, \
                 PICO_exp, Mdm, 127, 53, fn, fp, delta, mdens, nevents_PICO)

xsec_CRESST = CutAndCount90PercentConfXSec(CRESST_Emin, CRESST_Emax, \
                 CRESST_exp, Mdm, 184, 74, fn, fp, delta, mdens, nevents_CRESST)

xsec_XENON100 = CutAndCount90PercentConfXSec( XENON100_Emin, XENON100_Emax,\
                 XENON100_exp, Mdm, 131, 54, fn, fp, delta, mdens, nevents_XENON100)
plt.figure(1)
plt.semilogy()

plt.plot(delta,xsec_LUX_1_30,'-b',linewidth=2,label="LUX, E = 1-30 keVr")
plt.plot(delta,xsec_LUX_NEST,'--b',linewidth=2, label="LUX, E=1-256 keVr")
plt.plot(delta,xsec_LUX_best,':b',linewidth=2,label="LUX, E=1-500 keVr")

plt.plot(delta,xsec_PICO,'-g',linewidth=2,label="PICO (current)")
plt.plot(delta,xsec_CRESST,'-r',linewidth=2,label="CRESST (current)")
plt.plot(delta,xsec_XENON100,'-m',linewidth=2, label="XENON100 (current)")

plt.xlabel('Mass splitting (keV)');
plt.ylabel('Cross section (cm^2)');
plt.legend(loc='lower right')
plt.axis([0.,550.,1e-45,1e-32])
plt.show()
plt.savefig('InelasticDMProjection_m_10TeV.pdf')


