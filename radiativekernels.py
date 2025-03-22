#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 12:44:13 2021

@author: fpigeonneau
"""
import numpy as np

def F(x):
    return (0.5+x**2)/np.sqrt(1.+x**2)-x
# end F(x)
    
def K(x):
    return 1.-(x**3+1.5*x)/(1.+x**2)**1.5
# end K(x)
