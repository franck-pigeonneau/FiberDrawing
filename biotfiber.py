#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 08:24:57 2025

In this function, the Biot number due to the convection is computed according
to the law proposed by Kase and Matsuo (1965).

@author: fpigeonneau
"""

# ---------------------------------
# Importation of the useful modules
# ---------------------------------
import numpy as np

# --------------------------------
# Importation of our own functions
# --------------------------------

from argonproperties import nuAr
from argonproperties import lambdaAr
from argonproperties import rhoAr

# ---------------------------------------
# Beginning of the function ThermalSolver
# ---------------------------------------

def BiotFiber(aBiconv,Ts,P,b,T,a):
    """
    This function determines the Biot number using the heat transfer coefficient
    of Kase and Matsuo (1965) taking account the counter flow of gas in the furnace.

    Parameters
    ----------
    aBiconv : Float
        Prefactor of the Biot number resulting of the normalization.
    alpha : Float
        Ratio of a0/aL.
    beta : Float
        Aspect ratio a0/L.
    UAr0 : Float
        Normalized velocity of gas flow rate.
    Ts : Float
        Softening temperature of the glass.
    P : Float
        Absolute pressure in the liner.
    b : Float
        Normalized inner radius of the liner.
    T : Float
        Normalized temperature in the liner.
    a : Float
        Normalized fiber radius.

    Returns
    -------
    Float
        Biot number due to convective exchange.

    """
        
    # Determination of the Biot number
    return (aBiconv/a**(4./3.))*(lambdaAr(Ts*T)/lambdaAr(Ts))*(nuAr(Ts)/nuAr(Ts*T))**(1./3.)
# end BiotFiber
