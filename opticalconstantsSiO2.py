#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 08:59:14 2024

@author: fpigeonneau
"""

# -----------------------------
# Importation of useful modules
# -----------------------------

import numpy as np
from scipy import special

def SellmeierSiO2n(lmbd):
    """
    Determination of the refractive index of silica from the Sellmeier relation.
    This equation is valid for lambda<6.7 mum.

    Parameters
    ----------
    lmbd : Float
        Wavelength of light in m.

    Returns
    -------
    Float
        Refractive index of SiO2.

    """
    
    return np.sqrt(1.+0.6961663*lmbd**2/(lmbd**2-0.0684043**2)+\
                   0.4079426*lmbd**2/(lmbd**2-0.1162414**2)+\
                   0.8974794*lmbd**2/(lmbd**2-9.896161**2))
# end SellmeierSiO2n

def Gauss1(x,alpha,eta0,sig):
    eta=1.e-2/x # en m^-1
    return alpha*(np.exp(-(eta-eta0)**2/sig**2)-np.exp(-(eta+eta0)**2/sig**2))
#end Gauss1

def Gauss2(x,alpha,eta0,sig):
    eta=1.e-2/x # en m^-1
    x1=(eta+eta0)/sig
    x2=(eta-eta0)/sig
    return 2.*alpha*(special.dawsn(x1)-special.dawsn(x2))/np.sqrt(np.pi)
#end Gauss2

def permittivity(lmbd):
    """
    Determination of the permittivity of silica according the model provided by
    Kitamura et al. (2007).
    
    Parameters
    ----------
    lmbd : Float
        Wavelength of light.
    
    Returns
    -------
    epsprime : Float
        Real part of the permittivity.
    epssecond : Float
        Imaginary part of the permittivity.
    
    """
    alphaSiO2=np.array([3.7998,0.46089,1.2520,7.8147,1.0313,5.3757,6.3305,1.2948])
    etaSiO20=np.array([1089.7,1187.7,797.78,1058.2,446.13,443.00,465.80,1026.7])
    sigmaSiO2=np.array([31.454,100.46,91.601,63.153,275.111,45.220,22.680,232.14])/(2*np.sqrt(np.log(2.)))
    
    epsprime=2.1232
    epssecond=0
    N=np.size(alphaSiO2)
    for i in range(N):
        epsprime+=Gauss2(lmbd,alphaSiO2[i],etaSiO20[i],sigmaSiO2[i])
        epssecond+=Gauss1(lmbd,alphaSiO2[i],etaSiO20[i],sigmaSiO2[i])
    #end for
    
    return epsprime,epssecond
#end permittivity

def refractiveindexkappa(lmbd):
    """
    Determination the refractive index and absorption coefficient of the silica
    from the permittivity.

    Parameters
    ----------
    lmbd : Float
        Wavelength of light.

    Returns
    -------
    n : Float
        Refractive index of silica.
    kappa : Float
        Absorption coeeficient of silica.

    """
    epsprime, epssecond=permittivity(lmbd)
    n=np.sqrt(0.5*(epsprime+np.sqrt(epsprime**2+epssecond**2)))
    absorp=np.sqrt(0.5*(np.sqrt(epsprime**2+epssecond**2)-epsprime))
    kappa=4.*np.pi*absorp/lmbd
    return n,kappa
#end refractiveindexkappa


def nkappaSiO2(lmbd,lmbdSiO2,kappaSiO2):
    """
    This function determines the refractive index and the absorption coefficient.
    
    For the absorption coeefienct, when lambda<7.65 mum, the data of Khashan
    & Nassif (2001) and Drummond (1936) are used. For larger wavelength, the 
    fitting of Kitamura et al. (2007) is used.
    
    For the refractive index, the Sellmeir expression is used whan lambda <6.5 mum
    otherwise the fitting of Kitamura et al. (2007) is used.
    
    Parameters
    ----------
    lmbd : Array of float
        Wavelength of light for which n and kappa are determined.
    lmbdSiO2 : Array of float
        Wavelength of data of Drummond and Khashan & Nassif.
    kappaSiO2 : Array of float
        kappa of data of Drummond and Khashan & Nassif..

    Returns
    -------
    n : Array of float size N
        Refractive index of a silica.
    kappa :Array of float size N 
        Coefficient d'absorption of the silica.
    """
    # Size of lmbd
    # ------------
    
    N=np.size(lmbd)
    
    # Initialisation of n and kappa
    # -----------------------------
    n=np.zeros(N)
    kappa=np.zeros(N)
    
    # Determination of n and kappa
    # ----------------------------
    
    for i in range(N):
        # Determination of the kappa
        if (lmbd[i]<7.65e-6):
            kappa[i]=np.interp(lmbd[i],lmbdSiO2,kappaSiO2)
        else:
            n[i],kappa[i]=refractiveindexkappa(lmbd[i])
        #end if
        
        # Determination of the refractive index
        if (lmbd[i]<6.5e-6):
            n[i]=SellmeierSiO2n(lmbd[i]*1.e6)
        else:
            n[i],kappa2=refractiveindexkappa(lmbd[i])
        #end if
    #end for
    
    # Return of n and kappa
    # ---------------------
    
    return n,kappa
#end nkappaSiO2
