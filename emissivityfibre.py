#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 18:42:03 2021

This function determines the spherical emmisivity of a SiO2 glass fiber.
    
The absorption of the fused silica is obtained from the data of the literature
in the bandwidth of the visible and the close infrared.

The overall emissivity is obtained by integration over the wavelength of the
absorption of the fused silica.

    Input parameters:
        T: Temperarture [K]
        R: fiber radius [m]
        wavelength: wavelength of the spectrum [m]
        kappa: absorption coefficient [1/m]
    
    Output variables:
        emiss: spherical emmisivity of the fiber
    
@author: franck
"""

# -------------------------------------------------------------
# Importation of the physical constants used in the computation
# -------------------------------------------------------------

from scipy.constants import c
from scipy.constants import k
from scipy.constants import h

# -------------------------------
# Importation of the numpy module
# -------------------------------

import numpy as np

# -----------------------------------------
# Beginning of the function emissivityfibre
# -----------------------------------------

def emissivityfibre(T,R,wavelength,kappa):
    # Determination of zeta
    zeta=h*c/(k*T*wavelength)
    
    # Determination of the spectral emmisivity of the fibre with a radius R
    epsiloncyl=0.91*(1.-np.exp(-2.*R*kappa))
    
    # Determination of the average of the emmisivity weighted by the Planck's law
    fintegrand=15.*zeta**3*epsiloncyl/(np.pi**4*(np.exp(zeta)-1.))
    emiss=-np.trapz(fintegrand,zeta)

    # Return of the solution
    return emiss
# -----------------------------------
# End of the function emissivityfibre
# -----------------------------------
