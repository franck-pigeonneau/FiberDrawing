#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:50:06 2022

Determination of the normalised viscosity.

@author: franck
"""

from glassproperties import muVFT

def etanorm(Amu,Bmu,Tmu,T):
    etas=10.**6.65
    eta=muVFT(T,Amu,Bmu,Tmu)/etas
    if (eta>223872.11385683378):
        eta=223872.11385683378
    #endif
    return eta
#end etanorm
