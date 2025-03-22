#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 14:48:05 2021

This function determines the discretisation of the thermal equation which must
be balanced, i.e. equal to zero. The output is used in the nonlinear Newton 
solver to determine the temperature in the fibre.

Input variables:
    
    T: temperature of the fibre, solution of the Newton solver
    N: number of nodes of the segment [0;1]
    a0: radius of the preform [m]
    b: interior radius of the liner [-]
    Ts: softening temperature [K]
    alpha: ratio of the fiber at the exit to the radius of the preform
    beta: ratio of the radius of the preform to the length of the liner
    Pe: Péclet number
    Bo: Boltzmann number
    wavelength: array of the wavelength of the absorption of the fused silica [m]
    kappa: array of the absorption of the fused silica [m^-1]
    Twall: profile of the temperature of the wall of the liner normalised by Ts
    a: radius of the fibre normalized by R0
    qa: absorbed flux by the fibre normalised by sigma*Ts**4
    
Output variable:
    
    f: function corresponding to the discretisation form of the thermal equation
    which must be equal to zero.
    
References:
    
[1] Lu, Z., Blanc, W. and Pigeonneau, F. Prediction of the heating and coolind 
    rates of the fibre during the drawing, TechReport, 2021.

@author: F. Pigeonneau, Mines-Paristech, CEMEF (CFL)
"""

# ---------------------------------
# Importation of the useful modules
# ---------------------------------
import numpy as np

# --------------------------------
# Importation of our own functions
# --------------------------------

from emissivityfibre import emissivityfibre

# ---------------------------------------
# Beginning of the function ThermalSolver
# ---------------------------------------

def ThermalSolverFiniteDiff(T,N,a0,b,Ts,alpha,beta,Pe,Birad,Biconv,wavelength,kappa,\
                            Twall,a,qa):
    # Set of the step size in z
    dz=1./np.float64(N-1)
    
    # Determination of the square radius in each element
    Scell=np.zeros(N-1)
    for i in range(N-1):
        Scell[i]=0.5*(a[i]**2+a[i+1]**2)
    #end for
    
    # Determination of the radiative emission of the fiber
    qe=np.zeros(N)
    for i in range(N):
        qe[i]=T[i]**4*emissivityfibre(Ts*T[i],a0*a[i],wavelength,kappa)
    #end for
    
    f=np.zeros(N)
    f[0]=T[0]-Twall[0]
    c1=Pe*alpha**2/dz
    invdz2=1./dz**2
    
    # Determination of the discretisation of the heat transfer equation
    for i in range(1,N-1):
        f[i]=c1*(T[i]-T[i-1])-invdz2*(Scell[i-1]*T[i-1]-(Scell[i-1]+Scell[i])*T[i]+Scell[i]*T[i+1])+\
             2.*Birad*a[i]*(qe[i]-qa[i])+2.*Biconv[i]*a[i]*(T[i]-Twall[i])
    #end for
    #f[N-1]=3.*T[N-1]-4*T[N-2]+T[N-3]
    i=N-1
    f[i]=c1*(T[i]-T[i-1])-invdz2*(1.5*a[i]**2*(1.5*T[i]-2.*T[i-1]+0.5*T[i-2])-\
         a[i-1]**2*(T[i]-T[i-2])+0.5*a[i-2]**2*(T[i-1]-T[i-3]))+2.*Birad*a[i]*(qe[i]-qa[i])+\
         2.*Biconv[i]*a[i]*(T[i]-Twall[i])
    
    # Return of f
    return f
# -------------------------------------------
# End of the function ThermalSolverFiniteDiff
# -------------------------------------------






