#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 15:42:17 2021

This function determines the radius of the fiber according the solution valid
in the steady-state regime with a viscosity given by a VFT law.

@author: F. Pigeonneau, Mines-Paris - CEMEF (CFL)
"""

import numpy as np
from etanorm import etanorm

def FiberRadius(N,alpha,Amu,Bmu,Tmu,T):
    
    # Détermination of z
    z=np.linspace(0.,1,N)
    
    # Determination of the inverse of the viscosity
    inveta=np.ones(N)
    for i in range(N):
        inveta[i]=1./etanorm(Amu,Bmu,Tmu,T[i])
    #end for
    
    # Determination of the constant of integration
    A=np.log(alpha)/np.trapz(inveta,z)
    
    a=np.zeros(N)
    for i in range(N):
        a[i]=np.exp(A*np.trapz(inveta[0:i+1],z[0:i+1]))
    #end for
    
    # Return of the fiber radius
    return a
# end FiberRadius
