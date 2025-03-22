#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 12:57:27 2021

This function determines the radiative flux and the temperature of the liner
according the method described by Myers [1].

The problem is normalized as described in the document [2]: temperature is reduced
by the softening temperature and the radiative powers are divided by sigma*Ts^4.

The integral equation is solved by a trapezoidal quadrature.

Input parameters:
    
    z: axial position normalized by the total length of the liner
    b: radius of the liner [m]
    Le: length of the inlet domain non heated [m]
    Lg: length of the graphite domain heated [m]
    Ls: length of the exit domain non heated [m]
    Te: Temperature at each extremity normalised by Ts
    Tmax: Temperature maximum of the wall normalised by Ts
    epswall: emissivity of the wall of the liner
    
Output variables:
    Phiwall: normalised radiative flux emitted from the wall of the liner
    Twall: wall temperature of the liner normalised by Ts

References:

[1] Myers, M. R. A model for unsteady analysis of preform drawing AIChE J.,
    1989, 35, 592-602.
[2] Lu, Z., Blanc, W. and Pigeonneau, F. Prediction of the heating and coolind 
    rates of the fibre during the drawing, TechReport, 2021.

@author: F. Pigeonneau, Mines Paris - CEMEF
"""

# -----------------------------
# Importation of useful modules
# -----------------------------

import numpy as np
import matplotlib.pyplot as plt
from radiativekernels import K

# -----------------------------------------
# Function determining the source of energy
# -----------------------------------------

def S(z,Le,Lg,n):
    if (z<Le):
        s=0.
    else:
        if (z<Le+Lg):
            x=-1.+2.*(z-Le)/Lg
            s=(0.5*(1.+np.cos(np.pi*x)))**n
        else:
            s=0.
        # end if
    # end if
    return s
#end def

# ----------------------------------------
# Beginning of the function radiativeliner
# ----------------------------------------

def radiativeliner(N,b,Le,Lg,Ls,Te,Tmax,epswall,n,Plot,SaveFig):

    # Normalisation of the longitudinal dimension by the diameter of the liner
    Le=Le/(2.*b)
    Lg=Lg/(2.*b)
    Ls=Ls/(2.*b)
    L=Le+Lg+Ls

    # Normalisation of the length of the liner by the diameter of the liner
    z=np.linspace(0.,L,N)
    dz=z[1]-z[0]

    # ----------------------------------
    # Determination of the linear system
    # ----------------------------------
    
    A=np.zeros((N,N))
    B=np.zeros(N)

    # Computation of the coefficients of A of the first row, i=0
    # ----------------------------------------------------------
    i=0

    # First column, j=0
    j=0
    A[i,j]=1.-0.5*dz

    # For  j=1 to N-2
    for j in range(1,N-1):
        A[i,j]=-dz*K(np.abs(z[i]-z[j]))
    # end for

    # Case for j=N-1
    j=N-1
    A[i,j]=-0.5*dz*K(np.abs(z[i]-z[j]))

    # Determination of the source term
    B[i]=S(z[i],Le,Lg,n)

    # Computation of A and B for i=1 à N-2
    # ------------------------------------

    for i in range(1,N-1):
        # Case for j=0
        j=0
        A[i,j]=-0.5*dz*K(np.abs(z[i]-z[j]))
    
        # Case for j=1 à N-2
        for j in range(1,N-1):
            A[i,j]=np.float64(i==j)-dz*K(np.abs(z[i]-z[j]))
        # end for
    
        # Case for j=N-1
        j=N-1
        A[i,j]=-0.5*dz*K(np.abs(z[i]-z[j]))
    
        # Determination of the source term
        B[i]=S(z[i],Le,Lg,n)
    # end for

    # Computation of the coefficients of A of the last row, i=N-1
    # -----------------------------------------------------------
    i=N-1

    # Case for j=0
    j=0
    A[i,j]=-0.5*dz*K(np.abs(z[i]-z[j]))

    # Case for j=1 to N-2
    for j in range(1,N-1):
        A[i,j]=-dz*K(np.abs(z[i]-z[j]))
    # end for

    # Case for j=N-1
    j=N-1
    A[i,j]=1.-0.5*dz

    # Determination of the source term
    B[i]=S(z[i],Le,Lg,n)

    # -----------------------------
    # Solution of the linear system
    # -----------------------------
    
    Gamma=np.linalg.solve(A,B)
    
    # Determination of the maximum of the heat flux
    # ---------------------------------------------

    S0=(Tmax**4-Te**4)*epswall/(1.-epswall+epswall*np.max(Gamma))
    print('S0=',S0)
    
    # Computation of the wall temperature
    # -----------------------------------
    Twall=np.zeros(N)
    for i in range(N):
        Twall[i]=(S0*((1.-epswall)*B[i]/epswall+Gamma[i])+Te**4)**0.25
    # end for
    
    if (Plot):
        plt.figure()
        plt.plot(z,Gamma,'k',linewidth=2)
        plt.xlabel(r'$\tilde{z}$',fontsize=14)
        plt.ylabel(r'$\bar{\Gamma}$',fontsize=14)
        plt.twinx()
        plt.plot(z,B,'b',linewidth=2)
        plt.ylabel(r'$\bar{S}$',fontsize=14,color='b')
        plt.yticks(color='b')
        if (SaveFig):
            plt.savefig('GammaSvstildez.png',dpi=300,bbox_inches ='tight')
        #endif
    #endif
    
    # Return of the solution
    return Te**4+S0*Gamma, Twall
# ----------------------------------
# End of the function radiativeliner
# ----------------------------------
