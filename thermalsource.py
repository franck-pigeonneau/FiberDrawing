#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 14:32:02 2024

@author: fpigeonneau
"""

import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(-1,1,200)

n=10.

plt.figure()
plt.plot(x,(0.5*(1.+np.cos(np.pi*x)))**n)
plt.xlabel(r'$x$')
plt.ylabel(r'$S$')
