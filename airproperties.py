#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:19:21 2020


The function rhoair computes the air density according to the gas perfect law.

The function etaair computes the dynamical viscosity of dry air according to the
reference [1].

[1] Kadoya, K. and Matsunaga, N. and Nagashima, A., "Viscosity and Thermal
Conductivity of Dry Air in the Gaseous Phase", J. Phys. Chem. Ref. Data,
vol. 14(4), pp. 947--970 (1985).

@author: franck
"""
import numpy as np
from scipy.constants import R

# ---------------
# Dry air density
# ---------------

def rhoair(P,T):
    """
    Determination of the density of air from the perfect gas law.

    Parameters
    ----------
    P : Float
        Pressure of air in Pa.
    T : Float
        Temperarure of air in K.

    Returns
    -------
    Float
        density of air in kg/m^3.

    """
    Mair=2.89647e-2
    rair=R/Mair
    return P/(rair*T)
# end rhoair(P,T)

# ---------------------------
# Dry air dynamical viscosity
# ---------------------------

def etaair(P,T):
    H=6.1609e-6
    Tstart=132.5
    rhostart=314.3
    A1=0.128517
    A05=2.60661
    A=np.array([-1.00000,-0.197846,-0.709661,0.662534,0.00770147])
    B=np.array([0.465601,1.26469,-0.511425,0.2746])
    Tr=T/Tstart
    rhor=rhoair(P,T)/rhostart
    eta0=A1*Tr+A05*np.sqrt(Tr)
    for i in range(5):
        eta0+=A[i]*Tr**(-i)
    #end for
    deta=0.
    for i in range(4):
        deta+=B[i]*rhor**(i+1)
    #end for
    return H*(eta0+deta)
# end etaair(P,T)

# ---------------------------
# Dry air kinematic viscosity
# ---------------------------

def nuair(P,T):
    """
    Determination of the kinematic viscosity of aire

    Parameters
    ----------
    P : Float
        Pressure of air in Pa.
    T : Float
        Temperarure of air in K.

    Returns
    -------
    Float
        Kinematic viscosity in m^2/s.

    """
    return etaair(P,T)/rhoair(P,T)
# end nuair(P,T)

# ----------------------------
# Dry air thermal conductivity
# ----------------------------

def lambdaair(P,T):
    H=25.9778e-3
    Tstart=132.5
    rhostart=314.3
    A1=0.239503
    A05=0.00649768
    A=np.array([1.,-1.92615,2.00383,-1.07553,0.229414])
    B=np.array([0.402287,0.356603,-0.163159,0.138059,-0.0201725])
    Tr=T/Tstart
    rhor=rhoair(P,T)/rhostart
    lambda0=A1*Tr+A05*np.sqrt(Tr)
    for i in range(5):
        lambda0+=A[i]*Tr**(-i)
    #end for
    dlambda=0.
    for i in range(5):
        dlambda+=B[i]*rhor**(i+1)
    #end for
    return H*(lambda0+dlambda)
# end lambdaair(P,T)

# Determination of the specific heat capacity of dry air from N2 and O2 data

def CpN2(T):
    return 3.5+3390.**2*np.exp(-3390./T)/(T*(1.-np.exp(-3390./T)))**2
#end CpN2(T)

def CpO2(T):
    A1=2270.**2*np.exp(-2270./T)/(T*(1.-np.exp(-2270./T)))**2
    A2=(2.*11390.**2*np.exp(-11390./T)+18990.**2*np.exp(-18990./T))/(T**2*(3.+2.*np.exp(-11390./T)+np.exp(-18990./T)))
    A3=(2.*11390.*np.exp(-11390./T)+18990.*np.exp(-18990./T))**2/(T*(3.+2.*np.exp(-11390./T)+np.exp(-18990./T)))**2
    return 3.5+A1+A2-A3
# end CpO2(T)

def Cpair(T):
    xO2=0.21
    xN2=1.-xO2
    MAir=xO2*31.9988+xN2*28.0134
    return R*(xO2*CpO2(T)+xN2*CpN2(T))*1.e3/MAir
# end Cpair(T)

# Thermal diffusivity of air

def kappaair(P,T):
    return lambdaair(P,T)/(rhoair(P,T)*Cpair(T))
# end kappaair(P,T)
