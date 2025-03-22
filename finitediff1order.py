#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 17:16:25 2021

Determination of the spatial derivative of the first order by a finite difference
scheme at the second order.

@author: fpigeonneau
"""

import numpy as np

def FiniteDiff1order(h,F):
    N=np.size(F)
    dFdx=np.zeros(N)
    dFdx[0]=0.5*(-3.*F[0]+4.*F[1]-F[2])/h
    for i in range(1,N-1):
        dFdx[i]=0.5*(F[i+1]-F[i-1])/h
    #end for
    dFdx[N-1]=0.5*(3.*F[N-1]-4.*F[N-2]+F[N-3])/h
    
    # Return of the solution
    return dFdx
# end FiniteDiff1order
