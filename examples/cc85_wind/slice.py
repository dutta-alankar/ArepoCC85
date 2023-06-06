# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:00:45 2022

@author: alankar
"""

import arepo_run as arun
import gadget
import matplotlib.pylab as plt
import numpy as np
#import cmasher #This gives you additional, nice colour tables

xc, yc, zc = 25000, 25000, 25000
boxsize = [5000,50000]

num  = 43 # Snapshot number

o    = arun.Run(snappath='./output/',snapbase="snap_")

for num in range(42,48):
    s    = o.loadSnap(snapnum=num)

    s.plot_Aslice(value='rho', axes=[0,2], cmap='inferno', colorbar=True,
    box = boxsize, proj=False, proj_fact=0.01, center=[xc,yc,zc], res=1024,
    cblabel=r'Density [code units]',logplot=True, contour=False)

    plt.xlabel(r'$x \,\rm [pc]$', fontsize=18)
    plt.ylabel(r'$z \,\rm [pc]$', fontsize=18)

    plt.savefig('./output/plots/density.%04d.png'%num)

#plt.show()
