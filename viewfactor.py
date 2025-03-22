#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 16:30:44 2021

This function compute the view factor according to the formulation given by
Myers [1]. To do the computation, a change of variable is achieved in such a way

x=cos(theta).

The function integrate of scipy module is used to compute the integrate.

Input parameters:
    R: radius of the fiber
    dRdz: the first derivative repect to z of R
    z: axial position along the fibre 
    beta: ratio R0/L
    b: ratio of the liner
    zwall: axial position along the wall of the liner

Output variable:
    F: geometry factor
    
@author: fpigeonneau
"""

# ------------------------------
# Importation of usefull modules
# ------------------------------

import numpy as np
from scipy import integrate


def fvf(x,x0,beta,b,R,z,zwall):
    return (x-x0)*(1.-beta*R*x/b)/((1.+beta**2*R**2/b**2+((z-zwall)/b)**2-2.*\
                   beta*R*x/b)**2*np.sqrt(np.abs(1.-x**2)))

    
#def fvf(x,theta0,beta,b,R,z,zwall):
#    return (np.cos(x)-np.cos(theta0))*(1.-beta*R*np.cos(x)/b)/(1.+beta**2*R**2/b**2+((z-zwall)/b)**2-\
#                                                               2.*beta*R*np.cos(x)/b)**2
#end fvf

# ------------------------------------
# Beginning of the function ViewFactor
# ------------------------------------
    
def ViewFactor(z,zwall,R,dRdz,b,beta):
            
    # Inferior bound of the integral
    x0=(beta/b)*(R-(z-zwall)*dRdz)

    # Determination of the integral
    I,err=integrate.quad(fvf,x0,1.,args=(x0,beta,b,R,z,zwall),limit=200)
    
    # Determination of the view factor
    F=2.*I/(np.pi*b*np.sqrt(1.+(beta*dRdz)**2))
    
    # Return of the solution
    return F
# ------------------------------
# End of the function ViewFactor
# ------------------------------
