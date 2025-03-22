# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:52:02 2022

@author: Zhuorui Lu
"""

import numpy as np
from matplotlib import pyplot as plt

#df = pd.read_excel("C:\\Users\\Zhuorui Lu\\.spyder-py3\\Thermal profile\\result.xlsx", sheet_name='T1800 V21.4')
#x1 = df['z*L']
#x2 = df['a0*a']

# Measure of the temperature in the fiber tower of InPhyNi
zexpe= np.array([8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13.,13.5,14.,14.5,15.,15.5])
Texpe= np.array([1618.,1705.,1760.,1815.,1845.,1890.,1915.,1942.,1947.,1950.,1947.,1938.,1920.,1885.,1850.])

# Numerical solution
A=np.loadtxt('UL21Tmax2223.dat')
znum=A[:,0]
Twnum=A[:,7]

plt.figure()
plt.plot(zexpe,Texpe,'ko',label='Exp. data')
plt.plot(znum*1.e2,Twnum-273.15,'k',label='Num. sol.')
plt.xlabel(r'$z$ (cm)',fontsize=14)
plt.ylabel(r'$T$ (°C)',fontsize=14)
plt.legend(loc='upper right',fontsize=12,frameon=False)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('Tvszcompaexpenum.png',dpi=300,bbox_inches ='tight')
