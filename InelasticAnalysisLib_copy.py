import numpy as np
import scipy as sp


################################################################

def ER(Mdm, Mxe, v, theta_lab, delta, plus_or_minus):
   
   mu = Mdm*Mxe/(Mdm+Mxe)
   c2 = np.cos(theta_lab)**2
   first_term = mu*v**2*c2 - delta
   second_term = np.sqrt(mu*v**2*c2)*np.sqrt(mu*v**2*c2 - 2*delta)

   return mu/Mxe * ( first_term + plus_or_minus*second_term )


###########################################################################


def VelocityDist(v, ve, theta_earth, v0, v_esc):
   
   norm = np.pi**(3/2) * v0**3 * ( sp.special.erf(v_esc/v_0) - 2*v_esc/(np.pi**(1/2)*v0)*np.exp(-v_esc**2/v0**2) )

   return np.exp( - (v**2 + v_e**2 + 2*ve*v*np.cos(theta_earth))/v0**2 )/norm

###########################################################################

def HelmFF( ER, A ):

   c = 1.23*A**(1/3) - 0.6 # fm
   a = 0.52 # fm
   s = 1 # fm
   m_N = 0.932 # GeV/c^2
 
   rsq = c**2 + (7/3) * np.pi**2 * a**2 - 5 * s**2

   q = np.sqrt(2 * A * m_N * ER*1e-6 )/0.197

   j1 = np.sin(q*np.sqrt(rsq))/(q**2 * rsq) - np.cos(q*np.sqrt(rsq))/q*np.sqrt(rsq)

   return (3*j1)**2/(q**2 * rsq) * np.exp(-(q*s)**2)

################################################################

def IntegrateVelocityDist(ER, Mdm, A, delta, ve, v0, v_esc):

   m_N = 0.939 # neutron mass in GeV/c^2

   Mxe = m_N*A * 1e6 # keV/c^2
   Mdm = Mdm * 1e6 # keV/c^2

   mu = Mxe*Mdm/(Mxe + Mdm) # keV/c^2

   v_min = 1/np.sqrt(2*ER*Mxe) * (ER*Mxe/mu + delta) # units of c

   v_min = v_min * 3e5 # convert to units of km/s

 
   if (v_esc - ve) > v_min:

      # Compute the first integral
      v_int = np.linspace(v_min,v_esc-ve,1000000)
      dv = v_int[1] - v_int[0]
      array = v_int * velocityDistCosInt(v_int,1,ve,v0,v_esc)
      term1 = sum(array)*dv

      # Compute the second integral
      v_int = np.linspace(v_esc-ve,v_esc+ve,1000000)
      dv = v_int[1] - v_int[0]
      cstar = (v_esc**2 - v_int**2 - ve**2)/(2*ve*v_int)
      array = v_int * velocityDistCosInt(v_int,cstar,ve,v0,v_esc)     
      term2 = sum(array)*dv
      
      return (term1 + term2)*2*np.pi
   
   if ((v_esc - ve) < v_min) and (v_min < (v_esc + ve)):
     
      v_int = np.linspace(v_esc-ve,v_esc+ve,1000000)
      dv = v_int[1] - v_int[0]
      cstar = (v_esc**2 - v_int**2 - ve**2)/(2*ve*v_int)
      array = v_int * velocityDistCosInt(v_int,cstar,ve,v0,v_esc)

      return sum(array)*dv*2*np.pi

###########################################################################

def velocityDistCosInt(v, x, ve, v0, v_esc):
  
   norm = np.pi**(3/2) * v0**3 * ( sp.special.erf(v_esc/v0) - 2*v_esc/(np.pi**(1/2)*v0)*np.exp(-v_esc**2/v0**2) )

   prefactor = v0/(2*norm*ve*v)*(np.exp(-x) - np.exp(1))

   return prefactor * np.exp( -(v**2 + ve**2)/v0**2 )


###########################################################################

def Multiply(a,b):
   return a*b

###########################################################################

def RateVsEnergy(ER, Mdm, A, Z, fn, fp, xsec, mdens): 

   m_N = 0.939 # mass of neutron in GeV/c^2  
   Mxe = m_N * A

   mu = Mxe * Mdm / (Mxe + Mdm)

   ndens = mdens/Mdm

   xsec_prefactor = xsec * Mxe/(2 * mu**2) * (Z*fp + (A-Z)*fn)**2/fn**2 * HelmFF(ER,A)

     ## STOPPED HERE - NEED TO FIGURE OUT UNITS FOR EQ. 6
