#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 16:46:06 2021

This function determines the radius and the temperature of the fibre. The coupled
problem is composed of two equations corresponding to the radius solution of the
momentum equation in the steady-state regime and the energy balance. A Picard
scheme is used to solve the problem.

Input parameters:
    
    z: array corresponding of the axial coordinate norlized by the total length
    Rini: Profile of the initial (guess) solution of the fibre radius
    R0: Radius of the preform [m]
    L: total length of the liner [m]
    b: interior radius of the liner [m]
    Ts: softneting temperature of the fused silica [K]
    alpha: ratio of the fiber at the exit to the radius of the preform
    beta: ratio of the radius of the preform to the length of the liner
    Pe: Péclet number
    Bo: Boltzmann number
    numeps: numerical accuracy required in the Picard algorithm
    Amu: coefficient of the VFT law of the viscosity
    Bmu: coefficient of the VFT law of the viscosity
    Tmu: coefficient of the VFT law of the viscosity
    wavelength: array of the wavelength of the absorption of the fused silica [m]
    kappa: array of the absorption of the fused silica [m^-1]
    Phiwall: profile of the radiative flux normalised by sigma*Ts**4
    Twall: profile of the temperature of the wall of the liner normalised by Ts

Output variables:
    
    R: radius of the fibre from the inlet to the exit of the liner normalized by R0
    T: temperature of the fibre from the inlet to the exit of the liner normalized by Ts
    
@author: F. Pigeonneau, Mines-Paris - CEMEF (CFL)
"""

# -----------------------------
# Importation of useful modules
# -----------------------------

import numpy as np
from scipy.optimize import fsolve

# --------------------------------
# Importation of our own functions
# --------------------------------

from fiberradius import FiberRadius
from radiativeabsorption import RadiativeAbsorption
from thermalsolverfinitediff import ThermalSolverFiniteDiff
from biotfiber import BiotFiber

# ---------------------------------------------
# Beginning of the function lubricationsolution
# ---------------------------------------------

def lubricationsolution(N,kmax,numeps,a0,b,Ts,P,alpha,beta,Pe,Birad,aBiconv,Amu,Bmu,Tmu,\
                        wavelength,kappa,Phiwall,Twall):
    
    # Initialisation of the error
    L2normRT=1.e0
    
    # Initiliation of the temperature
    T=Twall
    Tnm1=T
    
    # Determination of the radius
    a=FiberRadius(N,alpha,Amu,Bmu,Tmu,Ts*T)
    anm1=a
    
    # Determination of the radiative absorption of the fiber
    qa=RadiativeAbsorption(N,Twall,Phiwall,b,Ts,beta,a0,a,wavelength,kappa)
    
    # Determination of the Biot number due to the convection
    Biconv=BiotFiber(aBiconv,Ts,P,Twall,a)
    
    # Determination of R and T by a Picard scheme
    k=1
    while (L2normRT>numeps and k<kmax):
        
        # Determination of the temperature
        
        T=fsolve(ThermalSolverFiniteDiff,Tnm1,args=(N,a0,b,Ts,alpha,beta,Pe,Birad,\
                 Biconv,wavelength,kappa,Twall,a,qa))
        
        # Determination of the radius
        a=FiberRadius(N,alpha,Amu,Bmu,Tmu,Ts*T)
        
        # Determination of the radiative absorption of the fiber
        qa=RadiativeAbsorption(N,Twall,Phiwall,b,Ts,beta,a0,a,wavelength,kappa)

        # Determination of the difference between the two iterations
        amean=np.mean(a)
        Tmean=np.mean(T)
        stda=np.std(a-anm1)
        stdT=np.std(T-Tnm1)
        L2normRT=np.sqrt((stda/amean)**2+(stdT/Tmean)**2)
        
        print('For k=',k,' L2normRT=',L2normRT)
        
        # Save for the next iteration
        anm1=a
        Tnm1=T
        k+=1
        
    # end while
    
    # ------------------------------
    # Return to the output variables
    # ------------------------------
    return a,T,qa,Biconv

# ---------------------------------------
# End of the function lubricationsolution
# ---------------------------------------
