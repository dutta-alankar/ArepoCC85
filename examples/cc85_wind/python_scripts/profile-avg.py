# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 10:09:53 2023

@author: alankar
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys

# dir = './output-relaxed'
# dir = './output-relaxedSteady'
dir = './output'
base = 'snapshot'

# constants
mu = 0.61
mp = 1.676e-24
Msun = 2e33
yr = 365*24*60**2
pc = 3.0857e18
kB = 1.38e-16
gamma = 5/3.
xH    = 0.76

# wind params
AGNWindMdot = 0.5*(Msun/yr)
AGNWindEdot = 2.e41
AGNWindPdot = 0.
Rinj        = 200*pc
rho0 = AGNWindMdot**(3/2)*AGNWindEdot**(-1/2)*Rinj**(-2)
P0   = AGNWindMdot**(1/2)*AGNWindEdot**(1/2)*Rinj**(-2)
v0   = AGNWindMdot**(-1/2)*AGNWindEdot**(1/2)

# units and constants
UnitVelocity_in_cm_per_s = 1.e5
UnitLength_in_cm         = pc
UnitMass_in_g            = 4.924e31

UnitTime_in_s = UnitLength_in_cm/UnitVelocity_in_cm_per_s
UnitU_in_cm_per_s_sq   = UnitVelocity_in_cm_per_s**2
UnitDensity_in_g_per_cm3 = UnitMass_in_g/UnitLength_in_cm**3

os.system(f'mkdir -p {dir}/plots/profile')

# free parameters
boxsize   = 50000                              # pc
boxhres   = 10000                              # pc
lrpoints  = 10                                 # low resolution fill of the cell
res       = 128                                # roughly res cels per side at high resolution
nside     = 6                                  # HEALPIX base resolution nside
bndryR    = Rinj/UnitLength_in_cm              # radius of boundary
size_fac  = 1.                                 # thickness of boundary layers in units of half the grid size
SpherePos = [0.5, 0.5, 0.5]                    # position of the sphere centre
SpherePos = np.array(SpherePos)*boxsize
pert      = 0.1                                # magnitude of mesh perturbation in terms of dx
contrast  = 10000.0

cc85 = np.loadtxt('CC85_steady-prof_gamma_1.667.txt')

for num in range(0,int(sys.argv[1])+1):
    print(num, end='\r')
    data = {}
    rTill = 5000 # pc
    with h5py.File(f'{dir}/{base}_{num:03d}.hdf5', 'r') as hdf:
        # generate array
        count = np.array(hdf['/PartType0/AGNFlag']).shape[0]
        data['pos']       = np.zeros((count,3), dtype=np.float32)
        data['velocity']  = np.zeros((count,3), dtype=np.float32)
        data['density']   = np.zeros(count, dtype=np.float32)
        data['mass']      = np.zeros(count, dtype=np.float32)
        data['pressure']  = np.zeros(count, dtype=np.float32)

        # assign the solution
        data['pos'][:,:]      = np.array(hdf['/PartType0/Coordinates'])
        data['pressure'][:]   = np.array(hdf['/PartType0/Pressure'])
        data['density'][:]    = np.array(hdf['/PartType0/Density'])
        data['mass'][:]       = np.array(hdf['/PartType0/Masses'])
        data['velocity'][:,:] = np.array(hdf['/PartType0/Velocities'])

        rr      = np.sqrt((data['pos'][:,0] - SpherePos[0])**2 + \
                        (data['pos'][:,1] - SpherePos[1])**2 + \
                        (data['pos'][:,2] - SpherePos[2])**2)

        ind     = np.where( np.logical_not(np.logical_and(rr>=150, rr<=rTill)) )
        remove  = np.size(ind)             # number of cells that need removal
        # print('Remove %d cell/s'%remove)

        newposx = np.delete(data['pos'][:,0],ind)
        newposy = np.delete(data['pos'][:,1],ind)
        newposz = np.delete(data['pos'][:,2],ind)
        data['pos'] = np.vstack( (newposx, newposy, newposz) ).T

        newvelx = np.delete(data['velocity'][:,0],ind)
        newvely = np.delete(data['velocity'][:,1],ind)
        newvelz = np.delete(data['velocity'][:,2],ind)
        data['velocity'] = np.vstack( (newvelx, newvely, newvelz) ).T

        data['pressure'] = np.delete(data['pressure'][:],ind)
        data['density']  = np.delete(data['density'][:],ind)
        data['mass']     = np.delete(data['mass'][:],ind)
        data['radius']   = np.delete(rr[:],ind)

    # print('Cutout: ', np.min(data['radius']), np.max(data['radius']) )
    # print('Shapes: ', rr.shape, data['radius'].shape )
    # print('Cutout: ', data['radius'])

    data['density']  = data['density']*UnitDensity_in_g_per_cm3
    data['velocity'] = data['velocity']*UnitVelocity_in_cm_per_s
    data['radius']   = data['radius']*UnitLength_in_cm
    data['pressure'] = data['pressure']*UnitDensity_in_g_per_cm3*UnitVelocity_in_cm_per_s**2
    data['mass']     = data['mass']*UnitMass_in_g

    # Perform 2D histogram on the data
    hist, x_edges, y_edges = np.histogram2d(np.log10(data['radius']/Rinj), np.log10(data['density']/rho0),
                    bins = (100,100), density=True)

    # Find the percentiles of the histogram
    percentiles = [16, 50, 84]
    percentile_values = np.percentile(hist, percentiles)

    # Generate the coordinates for the percentile lines
    x_coords, y_coords = np.meshgrid( 0.5*(x_edges[1:]+x_edges[:-1]), 0.5*(y_edges[1:]+y_edges[:-1]) )

    plt.figure(figsize=(13,10))
    '''
    for i in range(len(percentiles)):
        y_line = y_coords[np.where(hist>percentile_values[i])[0][0],:]
        x_line = np.ones(len(y_line))*x_coords[np.where(hist>percentile_values[i])[0][0],0]
        plt.plot(x_line,y_line,color='tab:red', linestyle='dashed')
    '''
    plt.pcolormesh(x_edges, y_edges, hist.T, cmap='viridis', norm='log')
    plt.colorbar()
    plt.plot(np.log10(cc85[:,0]), np.log10(cc85[:,1]), linewidth=3, color='black')


    for i in range(len(percentiles)):
        plt.contour(x_coords, y_coords, hist.T, levels=[percentile_values[i]],
                    colors='tab:red', linestyles='dashed')

    plt.xlim(0., 1.35)
    plt.ylim(-3.9, -0.8)
    plt.xlabel(r'Distance [$pc$] (log)', size=28)
    plt.ylabel(r'Normalized density' ,size=28)
    ax = plt.gca()
    # Get artists and labels for legend and chose which ones to display
    # handles, labels = ax.get_legend_handles_labels()
    # from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
    # leg = plt.legend([(p1, p2, p3), (p4, p5, p6), (p1,p4), (p2,p5), (p3,p6)],
    #                 [labels[0], labels[3],
    #                 r'$R_{cl}/\Delta$ = %d'%8,
    #                 r'$R_{cl}/\Delta$ = %d'%16,
    #                 r'$R_{cl}/\Delta$ = %d'%32], numpoints=1,
    #                 handler_map={tuple: HandlerTuple(ndivide=None)},
    #                 loc='lower left', ncol=1, fancybox=True, fontsize=25)
    plt.tick_params(axis='both', which='major', length=12, width=3, labelsize=24)
    plt.tick_params(axis='both', which='minor', length=8, width=2, labelsize=22)
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{dir}/plots/profile/dens_snap_{num:03d}.png') #, transparent=True)
    # plt.show()
    plt.close()

    plt.figure(figsize=(13,10))

    # Perform 2D histogram on the data
    hist, x_edges, y_edges = np.histogram2d(np.log10(data['radius']/Rinj), np.log10(data['pressure']/P0),
                    bins = (100,100), density=True)

    # Find the percentiles of the histogram
    percentiles = [16, 50, 84]
    percentile_values = np.percentile(hist, percentiles)

    # Generate the coordinates for the percentile lines
    x_coords, y_coords = np.meshgrid( 0.5*(x_edges[1:]+x_edges[:-1]), 0.5*(y_edges[1:]+y_edges[:-1]) )

    plt.pcolormesh(x_edges, y_edges, hist.T, cmap='viridis', norm='log')
    plt.colorbar()

    for i in range(len(percentiles)):
        plt.contour(x_coords, y_coords, hist.T, levels=[percentile_values[i]],
                    colors='tab:red', linestyles='dashed')
    plt.plot(np.log10(cc85[:,0]), np.log10(cc85[:,2]), linewidth=3, color='black')

    plt.xlim(0., 1.35)
    plt.ylim(-7.9, -1.1)
    plt.xlabel(r'Distance [$pc$] (log)', size=28)
    plt.ylabel(r'Normalized pressure' ,size=28)
    ax = plt.gca()
    # Get artists and labels for legend and chose which ones to display
    # handles, labels = ax.get_legend_handles_labels()
    # from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
    # leg = plt.legend([(p1, p2, p3), (p4, p5, p6), (p1,p4), (p2,p5), (p3,p6)],
    #                 [labels[0], labels[3],
    #                 r'$R_{cl}/\Delta$ = %d'%8,
    #                 r'$R_{cl}/\Delta$ = %d'%16,
    #                 r'$R_{cl}/\Delta$ = %d'%32], numpoints=1,
    #                 handler_map={tuple: HandlerTuple(ndivide=None)},
    #                 loc='lower left', ncol=1, fancybox=True, fontsize=25)
    plt.tick_params(axis='both', which='major', length=12, width=3, labelsize=24)
    plt.tick_params(axis='both', which='minor', length=8, width=2, labelsize=22)
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{dir}/plots/profile/pres_snap_{num:03d}.png') #, transparent=True)
    # plt.show()
    plt.close()

    plt.figure(figsize=(13,10))

    # Perform 2D histogram on the data
    hist, x_edges, y_edges = np.histogram2d(np.log10(data['radius']/Rinj), np.sqrt(np.sum(data['velocity']**2, axis=1))/v0,
                    bins = (100,100), density=True)

    # Find the percentiles of the histogram
    percentiles = [16, 50, 84]
    percentile_values = np.percentile(hist, percentiles)

    # Generate the coordinates for the percentile lines
    x_coords, y_coords = np.meshgrid( 0.5*(x_edges[1:]+x_edges[:-1]), 0.5*(y_edges[1:]+y_edges[:-1]) )

    plt.pcolormesh(x_edges, y_edges, hist.T, cmap='viridis', norm='log')
    plt.colorbar()

    for i in range(len(percentiles)):
        plt.contour(x_coords, y_coords, hist.T, levels=[percentile_values[i]],
                    colors='tab:red', linestyles='dashed')
    plt.plot(np.log10(cc85[:,0]), cc85[:,3], linewidth=3, color='black')

    plt.xlim(0., 1.35)
    plt.ylim(0., 2.0)
    plt.xlabel(r'Distance [$pc$] (log)', size=28)
    plt.ylabel(r'Normalized velocity' ,size=28)
    ax = plt.gca()
    # Get artists and labels for legend and chose which ones to display
    # handles, labels = ax.get_legend_handles_labels()
    # from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
    # leg = plt.legend([(p1, p2, p3), (p4, p5, p6), (p1,p4), (p2,p5), (p3,p6)],
    #                 [labels[0], labels[3],
    #                 r'$R_{cl}/\Delta$ = %d'%8,
    #                 r'$R_{cl}/\Delta$ = %d'%16,
    #                 r'$R_{cl}/\Delta$ = %d'%32], numpoints=1,
    #                 handler_map={tuple: HandlerTuple(ndivide=None)},
    #                 loc='lower left', ncol=1, fancybox=True, fontsize=25)
    plt.tick_params(axis='both', which='major', length=12, width=3, labelsize=24)
    plt.tick_params(axis='both', which='minor', length=8, width=2, labelsize=22)
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{dir}/plots/profile/vel_snap_{num:03d}.png') #, transparent=True)
    # plt.show()
    plt.close()

    # Mass
    plt.figure(figsize=(13,10))

    # Perform 2D histogram on the data
    hist, x_edges, y_edges = np.histogram2d(np.log10(data['radius']/Rinj), np.log10(data['mass']/Msun),
                    bins = (100,100), density=True)

    # Find the percentiles of the histogram
    percentiles = [16, 50, 84]
    percentile_values = np.percentile(hist, percentiles)

    # Generate the coordinates for the percentile lines
    x_coords, y_coords = np.meshgrid( 0.5*(x_edges[1:]+x_edges[:-1]), 0.5*(y_edges[1:]+y_edges[:-1]) )

    plt.pcolormesh(x_edges, y_edges, hist.T, cmap='viridis', norm='log')
    plt.colorbar()

    for i in range(len(percentiles)):
        plt.contour(x_coords, y_coords, hist.T, levels=[percentile_values[i]],
                    colors='tab:red', linestyles='dashed')

    plt.xlim(0., 1.35)
    # plt.ylim(0., 2.0)
    plt.xlabel(r'Distance [$pc$] (log)', size=28)
    plt.ylabel(r'Cell Mass [$M_{\odot}$] (log)' ,size=28)
    ax = plt.gca()
    # Get artists and labels for legend and chose which ones to display
    # handles, labels = ax.get_legend_handles_labels()
    # from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
    # leg = plt.legend([(p1, p2, p3), (p4, p5, p6), (p1,p4), (p2,p5), (p3,p6)],
    #                 [labels[0], labels[3],
    #                 r'$R_{cl}/\Delta$ = %d'%8,
    #                 r'$R_{cl}/\Delta$ = %d'%16,
    #                 r'$R_{cl}/\Delta$ = %d'%32], numpoints=1,
    #                 handler_map={tuple: HandlerTuple(ndivide=None)},
    #                 loc='lower left', ncol=1, fancybox=True, fontsize=25)
    plt.tick_params(axis='both', which='major', length=12, width=3, labelsize=24)
    plt.tick_params(axis='both', which='minor', length=8, width=2, labelsize=22)
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{dir}/plots/profile/mass_snap_{num:03d}.png') #, transparent=True)
    # plt.show()
    plt.close()

    # Volume
    plt.figure(figsize=(13,10))

    # Perform 2D histogram on the data
    hist, x_edges, y_edges = np.histogram2d(np.log10(data['radius']/Rinj),
                                            np.log10(data['mass'])-np.log10(data['density'])-3*np.log10(pc),
                                            bins = (100,100), density=True)

    # Find the percentiles of the histogram
    percentiles = [16, 50, 84]
    percentile_values = np.percentile(hist, percentiles)

    # Generate the coordinates for the percentile lines
    x_coords, y_coords = np.meshgrid( 0.5*(x_edges[1:]+x_edges[:-1]), 0.5*(y_edges[1:]+y_edges[:-1]) )

    plt.pcolormesh(x_edges, y_edges, hist.T, cmap='viridis', norm='log')
    plt.colorbar()

    for i in range(len(percentiles)):
        plt.contour(x_coords, y_coords, hist.T, levels=[percentile_values[i]],
                    colors='tab:red', linestyles='dashed')

    plt.xlim(0., 1.35)
    # plt.ylim(0., 2.0)
    plt.xlabel(r'Distance [$pc$] (log)', size=28)
    plt.ylabel(r'Cell Volume [$pc^3$] (log)' ,size=28)
    ax = plt.gca()
    # Get artists and labels for legend and chose which ones to display
    # handles, labels = ax.get_legend_handles_labels()
    # from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
    # leg = plt.legend([(p1, p2, p3), (p4, p5, p6), (p1,p4), (p2,p5), (p3,p6)],
    #                 [labels[0], labels[3],
    #                 r'$R_{cl}/\Delta$ = %d'%8,
    #                 r'$R_{cl}/\Delta$ = %d'%16,
    #                 r'$R_{cl}/\Delta$ = %d'%32], numpoints=1,
    #                 handler_map={tuple: HandlerTuple(ndivide=None)},
    #                 loc='lower left', ncol=1, fancybox=True, fontsize=25)
    plt.tick_params(axis='both', which='major', length=12, width=3, labelsize=24)
    plt.tick_params(axis='both', which='minor', length=8, width=2, labelsize=22)
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{dir}/plots/profile/vol_snap_{num:03d}.png') #, transparent=True)
    # plt.show()
    plt.close()
