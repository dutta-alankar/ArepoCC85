# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 14:56:15 2023

@author: Alankar
"""

import numpy as np
from scipy.interpolate import interp1d
import sys
#sys.argv[1] will have rTld

gamma = 5/3.
data = np.loadtxt('CC85_steady-prof_gamma_%.3f.txt'%gamma)

rhoTld = interp1d(data[:,0], data[:,1], fill_value='extrapolate')
prsTld = interp1d(data[:,0], data[:,2], fill_value='extrapolate')
velTld = interp1d(data[:,0], data[:,3], fill_value='extrapolate')

def momTerm(rTld): #Assumption: rTld>=1
    r = np.arange(0., rTld, 1e-3)
    global rhoTld, prsTld, velTld
    prsTld_at_rTld = prsTld(rTld)

    rhoTld = rhoTld(r)
    prsTld = prsTld(r)
    velTld = velTld(r)

    intg = np.trapz(r**2*np.gradient(prsTld, r), r)

    return prsTld_at_rTld - (1./rTld**2)*intg

with open(f'{sys.argv[2]}.txt', 'w') as file:
    file.write( '%f\n'%momTerm(float(sys.argv[1])) )
