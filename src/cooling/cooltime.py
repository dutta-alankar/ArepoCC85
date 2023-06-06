# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 23:37:21 2021

@author: alankar
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from decimal import Decimal

def fexp(number):
    (sign, digits, exponent) = Decimal(number).as_tuple()
    return len(digits) + exponent - 1

def fman(number):
    return Decimal(number).scaleb(-fexp(number)).normalize()

Msun = 2e33
yr = 365*24*60**2
mp = 1.6726219e-24
kB = 1.380649e-16
pc = 3.086e18
kpc = 1e3*pc
UNIT_VELOCITY = 1e5

X, Y, Z = 0.7154, 0.2703, 0.0143
Z *= 10**-0.2
X = 1.-(Y+Z)
mu = 1/(2*X+0.75*Y+0.5625*Z)
mue = 2./(1+X)
mui = 1/(1/mu-1/mue)
muH = 1/X

coolT = np.loadtxt('./cooltable_townsend.dat')
coolT = coolT[:-2,:]
cool = interpolate.interp1d(coolT[:,0], coolT[:,1], fill_value='extrapolate', kind='linear')

gamma = 5/3.
P = 1e3*kB

T = np.logspace(3,6, 1000)
Tmix = np.sqrt(1e3*1e6)

plt.figure(figsize=(13,10))
plt.loglog(T, cool(T), 'tab:purple', linewidth=5)
plt.vlines(Tmix, np.min(cool(T))/10, np.max(cool(T))*10, colors='tab:gray', linestyles='--', linewidth=5, label=r'$T_{mix}$')
plt.grid()
plt.ylim(ymin=np.min(cool(T)), ymax=np.max(cool(T))*5)
plt.text(Tmix+3000, 10**-23.8, r'$T_{\rm mix}$', fontsize=22)
plt.ylabel(r'Cooling function $\Lambda(T)$ [$erg cm^{3} s^{-1}$]', size=25 )
plt.xlabel(r'Gas temperature [$K$]', size=25 )
plt.tick_params(axis='both', which='major', labelsize=22)
plt.tick_params(axis='both', which='minor', labelsize=22)
plt.show()

tcool_ib = gamma/(gamma-1)*(1/P)/(cool(T)/(kB*T)**2)/(1e6*yr)
tcool_ib *= (mue*mui/mu**2)

plt.figure(figsize=(13,10))
plt.loglog(T, tcool_ib, 'tab:purple', linewidth=5)
plt.vlines(Tmix, np.min(tcool_ib)/10, np.max(tcool_ib)*10, colors='tab:gray', linestyles='--', linewidth=5, label=r'$T_{mix}$')
plt.grid()
plt.ylim(ymin=np.min(tcool_ib), ymax=np.max(tcool_ib))
plt.text(Tmix+3000, 8e-2, r'$T_{\rm mix}$', fontsize=22)
plt.ylabel(r'Isobaric cooling time [$Myr$] ($p/k_B=10^3~Kcm^{-3}$)', size=25 )
plt.xlabel(r'Gas temperature [$K$]', size=25 )
plt.tick_params(axis='both', which='major', labelsize=22)
plt.tick_params(axis='both', which='minor', labelsize=22)
plt.show()
