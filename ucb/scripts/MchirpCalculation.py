

import numpy as np
import random
from scipy.optimize import fsolve
import argparse

TSUN=4.9169e-6 ## this is the solar mass in seconds as defined in Constants.h, but it doesn't agree with lalsuite

p = argparse.ArgumentParser()

p.add_argument("--samples", type=str, required=True,help="File containing the posterior samples")
p.add_argument("--margq", action="store_true",help="Marginalize over mass ratio. Use this option if fdot is computed to 1PN")

args = p.parse_args()

## Convert from m1, m2 to Mc, eta
def m1_m2_to_mchirp_eta(m1,m2):
    Mchirp = np.power(m1*m2,3./5)*np.power(m1+m2,-1./5)
    eta=(m2*m1)/(m1+m2)/(m1+m2)
    return Mchirp,eta

# Convert from Mc, eta to m1, m2
def mchirp_q_to_m1_m2(mchirp,q):
    factor = mchirp*np.power(1+q, 1.0/5.0)
    m1 = factor*np.power(q, -3.0/5.0)
    m2 = factor*np.power(q, 2.0/5.0)
    return m1,m2

# 0PN fdot term (function of the chirp mass)
def fdot0PN(f0,m1,m2):
    Mc,eta=m1_m2_to_mchirp_eta(m1,m2)
    
    x=np.pi*Mc*f0
    
    return 96./5/np.pi/Mc/Mc*np.power(x,11./3)

# 1PN fdot term (function of the chirp mass and the mass ratio)
def fdot1PN(f0,m1,m2):
    Mc,eta=m1_m2_to_mchirp_eta(m1,m2)
    
    x=np.pi*Mc*f0
    
    return -2./35/np.pi/Mc/Mc*np.power(eta,-2./5)*(743+924*eta)*np.power(x,13./3)

# 1.5PN fdot term (function of the chirp mass and the mass ratio, and the spins)
# We do not use it as it doesn't seem to matter
def fdot1p5PN(f0,m1,m2):
    Mc,eta=m1_m2_to_mchirp_eta(m1,m2)
    
    x=np.pi*Mc*f0
    
    ## let's assume maximal spins. we find below that even in this case this term doesn't matter
    chi1=-1
    chi2=-1
    eta=(m1*m2)/(m1+m2)/(m1+m2)
    chieff=(m1*chi1+m2*chi2)/(m1+m2)
    chia=(chi1-chi2)/2
    delta=(m1-m2)/(m1+m2)
    
    beta = 1./3*((113-76*eta)/4.*chieff+76./4*delta*eta*chia)
    
    return -96./5/np.pi/Mc/Mc*np.power(eta,-3./5)*(beta-4*np.pi)*np.power(x,14./3)

# 0PN fddot term (function of the chirp mass)    
def fddot0PN(f0,m1,m2):
    Mc,eta=m1_m2_to_mchirp_eta(m1,m2)
    
    x=np.pi*Mc*f0
    
    return 33792./25/np.pi/Mc/Mc/Mc*np.power(x,19./3)

# t-domain phase using all the above terms
def PsiFull(t,f0,m1,m2):
    Psi0=2*np.pi*f0*t
    Psidot=np.pi*(fdot0PN(f0,m1,m2)+fdot1PN(f0,m1,m2))*t*t
    Psiddot=0.5*np.pi*fddot0PN(f0,m1,m2)*t*t*t
    
    return Psi0+Psidot+Psiddot

# t-domain phase using only the 0PN terms
def Psi0PN(t,f0,m1,m2):
    Psi0=2*np.pi*f0*t
    Psidot=np.pi*(fdot0PN(f0,m1,m2))*t*t
    Psiddot=0.5*np.pi*fddot0PN(f0,m1,m2)*t*t*t
    
    return Psi0+Psidot+Psiddot

# Compute Mc from f0 and 0PN fdot
def McfromFdot(f0,fdot):
    den=96./5*np.power(np.pi,8./3)*np.power(f0,11./3)
    
    return np.power(fdot/den,3./5)

# Compute Mc from f0 and fddot
def McfromFddot(f0,fddot):
    den=33792./25*np.power(np.pi,16./3)*np.power(f0,19./3)
    
    return np.power(fddot/den,3./10)

def func_opt(Mc,f0,fdot,q):
    m1,m2=mchirp_q_to_m1_m2(Mc,q)
    
    return np.abs(fdot-(fdot0PN(f0,m1,m2)+fdot1PN(f0,m1,m2)))

# Compute Mc from f0 and 0PN+1PN fdot
def McfromFdotFull(f0,fdot,q):
    InitialGuess = McfromFdot(f0,fdot)
    Mcsol = fsolve(func_opt,InitialGuess,args=(f0,fdot,q))[0]## we need to do this numerically
    
    return Mcsol


ChainData=np.genfromtxt(args.samples)
margq=1

f0samp=ChainData[:,0]
fdotsamp=ChainData[:,1]
fddotsamp=ChainData[:,8]

McFromFddot=McfromFddot(f0samp,fddotsamp)/TSUN

if (args.margq):
    np.random.seed(1)
    qmarg=np.random.uniform(0.1,1,len(f0samp))

    McFromFdot=[McfromFdotFull(f0samp[i],fdotsamp[i],qmarg[i])/TSUN for i in range(len(qmarg))]
else:
    McFromFdot=McfromFdot(f0samp,fdotsamp)/TSUN


with open('Mchirpchain.dat', 'w') as f:
    for Mc in McFromFdot:
        f.write("%0.10f\n" % Mc)






