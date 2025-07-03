#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:57:40 2025

Computation of the fluid properties of argon using a database generated from
Nist web site.

The properties are determined for one atmosphere, 101325 Pa and the temperature
ranges in [300;2000] K.

@author: fpigeonneau
"""

import numpy as np
from scipy.constants import R

def etaAr(T):
    """
    Dynamic viscosity of the argon

    Parameters
    ----------
    T : Float
        Temperature (K).

    Returns
    -------
    Float
        Dynamic viscosity in Pa.s.

    """
    return -1.91994324e-05+2.37689280e-06*np.sqrt(T)
#end etaAr

def lambdaAr(T):
    """
    Thermal conductivity of the argon

    Parameters
    ----------
    T : Float
        Temperature (K).

    Returns
    -------
    Float
        Thermal conductivity (W/mK)

    """
    return -0.01468708+0.00184834*np.sqrt(T)
#end lambdaAr

def rhoAr(T):
    """
    Density of argon at P=101325. Pa.

    Parameters
    ----------
    T : Float
        Temperure (K).

    Returns
    -------
    Float
        Density (kg/m^3).

    """
    P=101325.0
    MAr=0.039948
    rAr=R/MAr
    return P/(rAr*T)-0.00192367+5.87475e-07*T+0.036457/np.sqrt(T)
#end rhoAr

def nuAr(T):
    """
    Kinematic viscosity of argon.

    Parameters
    ----------
    T : Float
        Temperature (K).

    Returns
    -------
    Float
        Kinematic viscosity of argon (m^2/s).

    """
    return etaAr(T)/rhoAr(T)
#end nuAr

def kappaAr(T):
    """
    Kinematic viscosity of argon.
    
    Parameters
    ----------
    T : Float
       Temperature (K).
    
    Returns
    -------
    Float
        Thermal diffusivity of argon (m^2/s).
          
    """
    
    Cp=520.4814619883042
    return lambdaAr(T)/(rhoAr(T)*Cp)
#end kappaAr
