#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 15:00:26 2021

This function determines the radiative flux absorbed by the fiber. The general
form of this flux is given in Myers [1]. This flux is obtained by integration of
a product of three quantities: the radiative flux emitted by the wall of the liner,
a geometry factor integrated over the azimuthal angle and the emissivity of the
fibre at the temperature of the wall of the liner.

Input variables:
    
    N: Number of element over the length of the liner
    Twall: profile of the wall temperature normalised by Ts
    Phiwall: profile of the radiative flux normalised by sigma*Ts**4
    b: the interior radius of the liner normalised by the total length of the liner
    Ts: softneting temperature [K]
    beta: ratio of the radius of the preform to the length of the liner
    R0: radius of the preform [m]
    R: profile of the radius of the fibre normalised by R0
    wavelength: array of the wavelength of the absorption of the fused silica [m]
    kappa: array of the absorption of the fused silica [m^-1]
    
Output variable:
    
    qa: absorbed flux by the fibre normalised by sigma*Ts**4
    
@author: F. Pigeonneau, Mines Paris - CEMEF (CFL)
"""

# -------------------------------
# Importation of the module numpy
# -------------------------------
import numpy as np

# --------------------------------------------------------
# Importation of our own functions useful for the function
# --------------------------------------------------------

from finitediff1order import FiniteDiff1order
from emissivityfibre import emissivityfibre
from viewfactor import ViewFactor

# ---------------------------------------------
# Beginning of the function RadiativeAbsorption
# ---------------------------------------------

def RadiativeAbsorption(N,Twall,Phiwall,b,Ts,beta,R0,R,wavelength,kappa):
    
    # Size zwall
    z=np.linspace(0.,1.,N)
    dz=z[1]-z[0]
    
    # Determination of the derivative of R respect to z
    dRdz=FiniteDiff1order(dz,R)
    
    # Determination of qa
    qa=np.zeros(N)
    
    for i in range(N):
        # Determination of the absortivity of the fiber
        emiss=np.zeros(N)
        for j in range(N):
            emiss[j]=emissivityfibre(Ts*Twall[j],R0*R[i],wavelength,kappa)
        #end for
    
        # Determination of the view factor F
        F=np.zeros(N)
        for j in range(N):
            F[j]=ViewFactor(z[i],z[j],R[i],dRdz[i],b,beta)
        #end for
    
        # Determination of the absorption flux
        qa[i]=np.trapz(Phiwall*F*emiss,z)
    #end for
    
    # Return of the radiative absorption flux
    return qa
# ---------------------------------------
# End of the function RadiativeAbsorption
# ---------------------------------------
