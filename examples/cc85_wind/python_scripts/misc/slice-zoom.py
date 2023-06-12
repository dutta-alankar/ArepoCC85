# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:00:45 2022

@author: alankar
"""

import arepo_run as arun
import gadget
import matplotlib.pylab as plt
import numpy as np
import os
import sys
#import cmasher #This gives you additional, nice colour tables

os.system('mkdir -p ./output/plots/')
start = int(sys.argv[1])
stop = int(sys.argv[2])+1

xc, yc, zc = 25000, 25000, 25000

boxsize = [100,500]

print('Without cells')

o    = arun.Run(snappath='./output/',snapbase="snap_")
for num in range(start, stop):
    s    = o.loadSnap(snapnum=num)

    plt.figure(figsize=(8,13), num=num)
    pc = s.plot_Aslice(value='rho', axes=[0,1], cmap='inferno', colorbar=True,
    box = boxsize, proj=False, proj_fact=0.01, center=[xc,yc+500,zc], res=1024,
    cblabel=r'Density [code units]', logplot=True, contour=False, newfig=False)
    pc.set_clim(vmin=3e-3, vmax=3)

    plt.xlabel(r'$x \,\rm [pc]$', fontsize=18)
    plt.ylabel(r'$y \,\rm [pc]$', fontsize=18)
    plt.hlines(25285, 24990, 25010, colors='tab:green', linestyles='solid', linewidth=3)

    plt.savefig('./output/plots/density-noCont.%04d.png'%num)
    plt.close()

    plt.figure(figsize=(8,13), num=num)
    pc = s.plot_Aslice(value='pres', axes=[0,1], cmap='inferno', colorbar=True,
    box = boxsize, proj=False, proj_fact=0.01, center=[xc,yc+500,zc], res=1024,
    cblabel=r'Pressure [code units]', logplot=True, contour=False, newfig=False)
    pc.set_clim(vmin=20, vmax=2e4)

    plt.xlabel(r'$x \,\rm [pc]$', fontsize=18)
    plt.ylabel(r'$y \,\rm [pc]$', fontsize=18)
    plt.hlines(25285, 24990, 25010, colors='tab:green', linestyles='solid', linewidth=3)

    plt.savefig('./output/plots/pressure-noCont.%04d.png'%num)
    plt.close()
