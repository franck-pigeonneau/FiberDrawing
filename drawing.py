#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 19:37:21 2021

This program determines the profil of the fiber heating in the drawing furnage.
A lubrication model is used in steady-state regime. The heating of the fiber is
taken into account. 

The program is decomposed in various steps:
    
    1. Input of physical and numerical parameters
    2. Computation of the radiative flux and the liner temperature
    3. Determination of a, T, qa, Biconv
    4. Heating rate of the fiber
    5. Plotting
    6. Save of the results

The model is described in detail in ref. [1].

To easy change the glass properties, the glass liner chatacteristics, the 
working conditions are gathered in csv files. Here, these files are given for
our own applications.

The modules required to run this program are:
    
    - numpy;
    - matplotlib;
    - pandas;
    - os;
    - scipy.

References:

[1] F. Pigeonneau, Z. Lu & W. Blanc (2025). Thermal and mechanic behaviors of 
optic silica glass fiber during the drawing process, Int. J. Heat Mass Transfer,
under review.

List of variables
-----------------

AVFT: prefactor of the VFT law [Pa.s]
BVFT: coefficient of the activation energy of the viscosity [K]
TVFT: Temperature of the VFT law [K]
P: Atmospheric pressure [Pa]
etas: Dynamic viscosity at Ts [Pa.s]
Ts: Softening temperature [K]
rho: density [kg/m3]
Cp: specific heat capacity [J/kg/K]
lmbd: phonic thermal conductivity [W/m/K]
b: interior radius of the liner [m]
Linlet: Length of the entering part without heating [m]
Lheating: Length of the heating area in graphite [m]
Lexit: Length of the exit part without heating [m]
L: Total length of the liner [m]
epswall: wall emissivity of the liner
UL: drawing velocity [m/s]
a0: Radius of the fiber at the inlet [m]
aL: Radius of the fiber at the exit [m]
Te: Normalized temperature at the extremities, i.e. divided by Ts
Tmax: Maximum of normalized temperature of the liner
alpha: aspect ratio of the fibre (aL/a0)
a0sL: ratio of a0 to L
Pe: Péclet number
Birad: Radiative Biot number
aBiconv: Pre-factor of convective Biot number
N: Number of the nodes on the segment [0,1] corresponding to the normalized z
    coordinate
numeps: numerical accuracy asked in the Picard scheme
z: discrete z axis coordinate over the segment [0,1]
dz: length of the element between z[i] and z[i+1]
Phiwall: Radiative flux on the wall of the liner normalized by sigma*Ts**4
Twall: Normalized temperature of the wall of the liner
a: fiber radius
T: fiber temperature
qa: absorbed radiative fluw
Biconv: Convective Biot number
U: fiber velocity
dTdz: Temperature derivative respect to z
dTdt: Temperature rate
dadz: Fiber radius derivative respect to z
F: Drawing force
Fcap: Capillary force
t: Residence time of the glass in the liner

@author: F. Pigeonneau, Mines Paris - PSL - CEMEF (CFL)
"""

# ---------------------------------
# Importation of the useful modules
# ---------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# ----------------------------------
# Importation of physical constantes
# ----------------------------------

from scipy.constants import sigma
from scipy.constants import g
from scipy.constants import physical_constants

# --------------------------------
# Importation of our own functions
# --------------------------------

from glassproperties import Tsoft
from glassproperties import Tglass
from glassproperties import muVFT
from airproperties import nuair
from airproperties import lambdaair
from radiativeliner import radiativeliner
from lubricationsolution import lubricationsolution
from finitediff1order import FiniteDiff1order
from emissivityfibre import emissivityfibre

# ------------------------
# Constants of the problem
# ------------------------

P=physical_constants['standard atmosphere'][0]
etas=10.**6.65

# ---------------------------------------------
# 1. Input of physical and numerical parameters
# ---------------------------------------------

# Working conditions
Plot=True
SaveFig=True
dworking=pd.read_csv('working.csv',index_col=0)
UL=dworking['UL'].values[0]/60. #m/s
print('UL=',UL,'m/s')
Tmax=dworking['Tmax'].values[0]+273.15 # K
print('Tmax=',Tmax,'K')
a0=dworking['a0'].values[0] # m
print('a0=',a0,'m')
aL=dworking['aL'].values[0] # m
print('aL=',aL,'m')

# Glass liner
dliner=pd.read_csv('glassliner.csv',index_col=0)
b=dliner['b'].values[0]
print('b=',b,'m')
L=dliner['L'].values[0]
print('L=',L,'m')
Linlet=dliner['Linlet'].values[0]
print('Linlet=',Linlet,'m')
Lexit=dliner['Lexit'].values[0]
print('Lexit=',Lexit,'m')
Lheating=L-(Linlet+Lexit)
print('Lheating=',Lheating,'m')
epswall=dliner['epswall'].values[0]
print('epswall=',epswall)

# Glass properties
name='SiO2'
dglass=pd.read_csv('glass.csv',index_col=0)
iglass=np.argwhere(dglass.index==name)[0][0]
rho=dglass['rho'].values[iglass]
print('rho=',rho)
Cp=dglass['Cp'].values[iglass]
print('Cp=',Cp)
lmbd=dglass['lmbd'].values[iglass]
print('lmbd=',lmbd)
gamma=dglass['gamma'].values[iglass]
print('gamma=',gamma)
AVFT=dglass['AVFT'].values[iglass]
print('AVFT=',AVFT)
BVFT=dglass['BVFT'].values[iglass]
print('BVFT=',BVFT)
TVFT=dglass['TVFT'].values[iglass]
print('TVFT=',TVFT)
dkappa=pd.read_csv('kappaSiO2vslmbd.csv')
wavelength=dkappa['lambda'].values
kappa=dkappa['kappa'].values

# Softening & glass transition temperature
Ts=Tsoft(AVFT,BVFT,TVFT)
print('Ts=',Ts)
Tg=Tglass(AVFT,BVFT,TVFT)
print('Tg=',Tg)

# Air properties at T=Ts
nuairTs=nuair(P,Ts)
kairTs=lambdaair(P,Ts)

# Numerical parameters
dnum=pd.read_csv('numerics.csv',index_col=0)
N=dnum['N'].values[0]
numeps=dnum['numeps'].values[0]
kmax=dnum['kmax'].values[0]
beta=dnum['beta'].values[0]
dz=1./np.float64(N-1)
z=np.linspace(0,1.,N)

# Temperature normalised by Ts
Te=dnum['Te'].values[0]
Tmax/=Ts

# Dimensionless numbers
# ---------------------

# Radius aspect ratio
alpha=aL/a0
# Reynolds number
Re=rho*UL*L/etas
print('Re=',Re)
# Capillary number
Ca=etas*UL*a0/(gamma*L)
print('Ca=',Ca)
# Galilei number
Ga=rho*g*L**2/(etas*UL)
print('Ga=',Ga)
# Péclet number
Pe=L*UL*rho*Cp/lmbd
print('Pe=',Pe)
# Radiative Biot number
Birad=sigma*Ts**3*L**2/(lmbd*a0)
print('Birad=',Birad)
# Pre-factor of convective Biot number
aBiconv=0.21*(L/a0)**2*alpha**(2./3.)*(kairTs/lmbd)*(2.*a0*UL/nuairTs)**(1./3.)

# Dimensionless dimensions
# ------------------------

a0sL=a0/L
Linlet/=L
Lheating/=L
Lexit/=L
b/=L

# --------------------------------------------------------------
# 2. Computation of the radiative flux and the liner temperature
# --------------------------------------------------------------

Phiwall,Twall=radiativeliner(N,b,Linlet,Lheating,Lexit,Te,Tmax,epswall,beta,Plot,SaveFig)

# ------------------------------------
# 3. Determination of a, T, qa, Biconv
# ------------------------------------

a,T,qa,Biconv=lubricationsolution(N,kmax,numeps,a0,b,Ts,P,alpha,a0sL,Pe,Birad,
                                  aBiconv,AVFT,BVFT,TVFT,wavelength,kappa,
                                  Phiwall,Twall)

# ----------------------------
# 4. Heating rate of the fiber
# ----------------------------

# 4.1 Determination of the velocity
# ---------------------------------

U=alpha**2/a**2

# 4.2 Determination of the heating rate
# -------------------------------------

dTdz=FiniteDiff1order(dz,T)
dTdt=Ts*UL*dTdz*U/L

# 4.3 Determination of the tension
# --------------------------------

dUdz=FiniteDiff1order(dz,U)
F=3.*a0**2*(UL/L)*muVFT(Ts*T,AVFT,BVFT,TVFT)*dUdz*np.pi*a**2
print('F=',F[int(N/2)])

# Determination of the amplitude of the inertia and capillary forces
dadz=FiniteDiff1order(dz,a)
Fcap=-2.*gamma*a0*a*np.pi/np.sqrt(1.+(a0*dadz/L)**2)
print('max Fcap : ', np.max(np.abs(Fcap)))

# Time residence
t=np.zeros(N)
for i in range(1,N):
    t[i]=t[i-1]+0.5*(1./U[i-1]+1./U[i])*dz
#end for

# -----------
# 5. Plotting
# -----------

directory='UL'+str(int(60.*UL))+'Tmax'+str(int(Tmax*Ts))
if (not os.path.exists(directory)):
    os.mkdir(directory)
#end if

if (Plot):
    # Determination of time residence
    t=np.zeros(N)
    for i in range(1,N):
        t[i]=t[i-1]+0.5*(1./U[i-1]+1./U[i])*dz
    #end for
    
    fig1,ax1=plt.subplots()
    ax1.plot(z*L,Ts*T,'k',linewidth=2)
    ax1.plot(z*L,Ts*Twall,'b',linewidth=2)
    ax1.set_xlabel('$z$ (m)')
    ax1.set_ylabel('$T$ (K)')
    ax1.legend((r'Fibre',r'Wall'),loc=0)
    
    fig2,ax2=plt.subplots()
    ax2.plot(z*L,a0*a*1.e3,'k')
    ax2.set_xlabel('$z$ (m)')
    ax2.set_ylabel('$a$ (mm)')
    
    fig3,ax3=plt.subplots()
    ax3.plot(z*L,dTdt,'k')
    ax3.set_xlabel('$z$ (m)')
    ax3.set_ylabel(r'$\frac{dT}{dt}$ (K/s)')

    fig4,ax4=plt.subplots()
    ax4.plot(z[:-3]*L,F[:-3],'k',label='Num. sol.')
    ax4.set_xlabel('$z$ (m)')
    ax4.set_ylabel(r'$F$ (N)')
    ax4.set_ylim((0.8*np.min(F),2*np.min(F)))
    
    # Computation of the radiative and convective thermal source terms
    qe=np.zeros(N)
    for i in range(N):
        qe[i]=T[i]**4*emissivityfibre(Ts*T[i],a0*a[i],wavelength,kappa)
    #end for
    phirad=-2.*Birad*a*(qe-qa)
    # Determination of the Biot number due to the convection
    phiconv=-2.*Biconv*a*(T-Twall)
    
    fig5,ax5=plt.subplots()
    ax5.plot(z,phirad,'k',label=r'Radiative flux',linewidth=2)
    ax5.plot(z,phiconv,'r',label=r'Convective flux',linewidth=2)
    ax5.set_xlabel(r'$\bar{z}$')
    ax5.set_ylabel(r'$\varphi$')
    ax5.legend(loc=0)
    
    fig6,ax6=plt.subplots()
    ax6.semilogx(t[1:]*L/UL,Ts*T[1:],'k',linewidth=2)
    ax6.set_xlabel(r'$t$ (s)')
    ax6.set_ylabel(r'$T$ (K)')
    
    if (SaveFig):
        fig1.savefig(directory+'/Tvsz.png',dpi=300,bbox_inches ='tight')
        fig2.savefig(directory+'/avsz.png',dpi=300,bbox_inches ='tight')
        fig3.savefig(directory+'/dTdtvsz.png',dpi=300,bbox_inches ='tight')
        fig4.savefig(directory+'/Fvsz.png',dpi=300,bbox_inches ='tight')
        fig5.savefig(directory+'/phitermvsz.png',dpi=300,bbox_inches ='tight')
        fig6.savefig(directory+'/Tvst.png',dpi=300,bbox_inches ='tight')
    #end if
#end if

# ----------------------
# 6. Save of the results
# ----------------------

A=np.transpose(np.array([z*L,t*L/UL,a0*a,Ts*T,UL*U,dUdz*UL/L,dTdt,F,Twall*Ts,Phiwall*sigma*Ts**4]))
data=pd.DataFrame(A,columns=np.array(['z','t','a','T','U','dUdz','dTdt','F','Twall','Phiwall']))
data.to_csv(directory+'/results.csv')

# Cooling fiber after the exit of the furnace
cooling=False
if (cooling):
    hair=300.
    Tair=293.15
    lout=rho*Cp*aL*UL/hair
    zg=rho*Cp*aL*UL*np.log((Ts*T[N-1]-Tair)/(Tg-Tair))/(2.*hair)
    Nout=50
    zout=np.linspace(0.,zg,Nout)
    tout=zout/UL
    Tout=Tair+(Ts*T[N-1]-Tair)*np.exp(-2.*zout/lout)
    dToutdt=-2.*(Ts*T[N-1]-Tair)*np.exp(-2.*zout/lout)*UL/lout
    
    # Overall figure
    
    ztotal=L*z
    ztotal=np.append(ztotal,zout[1:]+L)
    Ttotal=Ts*T
    Ttotal=np.append(Ttotal,Tout[1:])
    ttotal=t*L/UL
    ttotal=np.append(ttotal,tout[1:]+t[N-1]*L/UL)
    
    dTtotaldt=dTdt
    dTtotaldt=np.append(dTtotaldt,dToutdt[1:])
    
    Tmic=1750.+273.15
    iloc=np.where(Ttotal>Tmic)[0]
    
    plt.figure()
    plt.plot(ttotal,Ttotal,'k',linewidth=2)
    if (np.size(iloc)>0):
        plt.plot(ttotal[iloc],Tmic*np.ones(np.size(iloc)),'r',linewidth=2)
        plt.annotate('$t$='+str(np.round(ttotal[iloc[0]],2)),(ttotal[iloc[0]]-300,Ttotal[iloc[0]]+10))
        plt.annotate('$t$='+str(np.round(ttotal[iloc[-1]],2)),(ttotal[iloc[-1]]-200,Ttotal[iloc[-1]]+10))
        plt.annotate(r'$T_{\mathrm{bi}}$='+str(Tmic),(np.mean(ttotal[iloc])-300,Tmic-40))
    #endif
    plt.xlabel(r'$t$ (s)',fontsize=14)
    plt.ylabel(r'$T$ (K)',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    file='Tvst'+str(int(60.*UL))+'Tmax'+str(int(Tmax*Ts))
    plt.savefig(file+'.png',dpi=300,bbox_inches ='tight')
    
    tTmaxfib=ttotal[Ttotal==np.max(Ttotal)]
    deltat=ttotal[np.size(ttotal)-1]-tTmaxfib[0]
    print('deltat=',deltat,' s')
#end if

plt.show()
