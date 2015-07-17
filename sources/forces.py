#Import of several modules

import math
import time
import numpy as np
from scipy import integrate



#================================= Long range forces ==================================

def vdw (z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc):
    if (z + zc) < a0:
        force = 0
    else:
        force = -H*Rt/(6*(z + zc)*(z + zc)) 
    return force

def dlvo (z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc):
    epsilon = 8.854187*10**-12
    if (z + zc) < a0:
        force = 4*math.pi*Rt/(epsilon*epsilonr)*sigmat*sigmas*landadeb*np.exp(-(a0)/landadeb)
    else:
        force = -H*Rt/(6*(z + zc)*(z + zc)) + 4*math.pi*Rt/(epsilon*epsilonr)*sigmat*sigmas*landadeb*np.exp(-(z + zc)/landadeb)

    return force


#================================= Short range forces ==================================


def hertz (z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc):
    if (z + zc) < a0:
        force = 4.0/3.0*E*math.sqrt(Rt)*math.pow(-(z + zc - a0),3.0/2.0)
    else:
        force = 0

    return force

def dmt (z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc):
    if (z + zc) < a0:
        force = 4.0/3.0*E*math.sqrt(Rt)*math.pow(-(z + zc - a0),3.0/2.0) -H*Rt/(6*a0*a0)
    else:
        force = 0

    return force   

def tatara (z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc):
    if ((z + zc) < a0) and nc > 0:
        force = alfa*math.pow(-(z + zc - a0),1.5)/math.pow(2,1.5) + 3*alfa**2*(z + zc - a0)**2/(8*nc) + 15*alfa**3/(math.pow(2,11/2.0)*nc**2)*math.pow(-(z + zc - a0),2.5) -H*Rt/(6*a0*a0)
    else:
        force = 0

    return force   

def jkr (z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc):
    if (z + zc) < a0:
        force = 0
    else:
        force = 0

    return force 


#================================= Viscosity forces ==================================

def visco (z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc):
    if (z + zc) < a0:
        force = -eta*math.sqrt(-Rt*(z + zc))*v
    else:
        force = 0

    return force

#================================= Driven forces ==================================

def driven (F0,f0,t):
    if f0.size == 1:
        force = F0*math.cos(f0*2*math.pi*t)
    else:
        force = F0[1]*math.cos(f0[1]*2*math.pi*t) + F0[2]*math.cos(f0[2]*2*math.pi*t)

    return force




    return force

f_vector = [dlvo,vdw,hertz,dmt,tatara,visco]

def sumatory(option_mask,z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc):
    force = sum([f(z,v,t,f0,k,Q,Rt,a0,E,H,mu,m,D,F0,zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc) for f,mask in zip(f_vector,option_mask) if mask]) + driven(F0,f0,t)

    return force
