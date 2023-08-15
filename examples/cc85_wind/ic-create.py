# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 01:04:40 2021

@author: alankar

Masses are interpreted as densities.
Make sure corresponding Arepo flag is on in Config.sh
"""

import numpy as np
import healpy as hp
#import matplotlib.pylab as plt
import sys

from gadget import gadget_write_ics
from const import KB

sys.path.append('./arepo-snap-util')

#background parameters
nH_bg    = 300.
temp_bg  = 2.e4
mu_bg    = 0.61           # mean particle weight
gamma    = 5./3

#units and constants
unit_l   = 3.0857e21
unit_v   = 1.65e6
unit_m   = 1.9885e43
m_p      = 1.673e-24
xH       = 0.76

unit_t   = unit_l/unit_v
unit_u   = unit_v**2
unit_rho = unit_m/unit_l**3

#free parameters
boxsize  = 10
res      = 128                                # roughly res cels per side
nside    = 6                                  # HEALPIX base resolution nside
bndryR   = boxsize/res                        # radius of boundary
size_fac = 1.                                 # thickness of boundary layers in units of half the grid size
SpherePos= [boxsize/2, boxsize/2, boxsize/2]  # position of the sphere centre
pert     = 0.1                                # magnitude of mesh perturbation in terms of dx
contrast = 1e-3

#other quantities (do not change)
ndim     = 3
rho_bg   = m_p * nH_bg / xH
cs       = np.sqrt(gamma * KB * temp_bg / (mu_bg * m_p))
dgrid    = boxsize/res
ncell    = res**ndim
dx       = boxsize / ncell**(1./ndim)
npix     = hp.pixelfunc.nside2npix(nside)
volbndry = (size_fac*dgrid/2.)**3
vol      = dx**3

#generate arrays
data = {}
data['pos']       = np.zeros((ncell + 2*npix,3), dtype=np.float32)
data['vel']       = np.zeros((ncell + 2*npix,3), dtype=np.float32)
data['mass']      = np.zeros(ncell  + 2*npix, dtype=np.float32)
data['u']         = np.zeros(ncell  + 2*npix, dtype=np.float32)
data['type']      = np.zeros(ncell  + 2*npix, dtype=np.float32)
data['flga']      = np.zeros(ncell  + 2*npix, dtype=int)
data['count']     = ncell + 2*npix

rho               = rho_bg
u                 = 1./(gamma - 1)/gamma * cs**2

# Add boundary cells
ipix   = np.arange(npix)
xx,yy,zz = hp.pixelfunc.pix2vec(nside, ipix)

# Inner Layer
data['pos'][ncell:ncell+npix,0]           = (bndryR-size_fac*dgrid/2.) * xx + SpherePos[0]
data['pos'][ncell:ncell+npix,1]           = (bndryR-size_fac*dgrid/2.) * yy + SpherePos[1]
data['pos'][ncell:ncell+npix,2]           = (bndryR-size_fac*dgrid/2.) * zz + SpherePos[2]
data['flga'][ncell:ncell+npix]            = 1
data['mass'][ncell:ncell+npix]            = rho / unit_rho *contrast
data['u'][ncell:ncell+npix]               = u / unit_u

#Outer Layer
data['pos'][ncell+npix:ncell+2*npix,0]    = (bndryR+size_fac*dgrid/2.) * xx + SpherePos[0]
data['pos'][ncell+npix:ncell+2*npix,1]    = (bndryR+size_fac*dgrid/2.) * yy + SpherePos[1]
data['pos'][ncell+npix:ncell+2*npix,2]    = (bndryR+size_fac*dgrid/2.) * zz + SpherePos[2]
data['flga'][ncell+npix:ncell+2*npix]     = 2
data['mass'][ncell+npix:ncell+2*npix]     = rho / unit_rho *contrast
data['u'][ncell+npix:ncell+2*npix]        = u / unit_u

#Fill the entire box with normal AREPO cells
points   = np.linspace( 0, 1., endpoint=True, num=int(np.round(ncell**(1./ndim)))+1, dtype=np.float32 )
cell_pos = 0.5 * (points[1:]+points[:-1]) * boxsize
cells_x, cells_y, cells_z  = np.meshgrid(cell_pos,cell_pos,cell_pos)

data['pos'][:ncell,0] = cells_x.flatten()
data['pos'][:ncell,1] = cells_y.flatten()
data['pos'][:ncell,2] = cells_z.flatten()

# Perturb Cartesian grid
theta  = (np.random.rand(ncell).astype(np.float32)) * np.pi
phi    = (np.random.rand(ncell).astype(np.float32)) * np.pi * 2

pert_x = np.sin(theta) * np.cos(phi)
pert_y = np.sin(theta) * np.sin(phi)
pert_z = np.cos(theta)

data['pos'][:ncell,0]  += pert * dx * pert_x
data['pos'][:ncell,1]  += pert * dx * pert_y
data['pos'][:ncell,2]  += pert * dx * pert_z
data['mass'][:ncell]    = rho / unit_rho
#data['mass'][data['flga']>0]  = rho / unit_rho * 1.e-3
data['u'][:ncell]       = u / unit_u

# Remove normal AREPO cells from volume bounded by boundary sphere
rr      = np.sqrt((data['pos'][:,0] - SpherePos[0])**2 + (data['pos'][:,1] - SpherePos[1])**2 + (data['pos'][:,2] - SpherePos[2])**2)
flag    = data['flga']
ind     = np.where(np.logical_and( rr <= (bndryR+size_fac*dgrid), (flag == 0)) )
number  = np.size(ind)             # number of cells that need removal

newposx = np.delete(data['pos'][:,0],ind)
newposy = np.delete(data['pos'][:,1],ind)
newposz = np.delete(data['pos'][:,2],ind)

newvelx = np.delete(data['vel'][:,0],ind)
newvely = np.delete(data['vel'][:,1],ind)
newvelz = np.delete(data['vel'][:,2],ind)

newmass = np.delete(data['mass'],ind)
newtype = np.delete(data['type'],ind)
newu    = np.delete(data['u'],ind)
newflga = np.delete(data['flga'],ind)

rr = np.delete(rr,ind)

# Update Data and Save
ncell = ncell - number
data['pos']       = np.zeros((ncell + 2*npix,3))
data['vel']       = np.zeros((ncell + 2*npix,3))
data['mass']      = np.zeros(ncell  + 2*npix)
data['u']         = np.zeros(ncell  + 2*npix)
data['type']      = np.zeros(ncell  + 2*npix)
data['flga']      = np.zeros(ncell  + 2*npix, dtype=int)
data['count']     = data['count'] - number

data['pos'][:,0]  = newposx
data['pos'][:,1]  = newposy
data['pos'][:,2]  = newposz
data['vel'][:,0]  = newvelx
data['vel'][:,1]  = newvely
data['vel'][:,2]  = newvelz
data['mass']      = newmass
data['u']         = newu
data['type']      = newtype
data['flga']      = newflga

print('ReferenceGasPartMass: ', rho / unit_rho*(boxsize/res)**3)
gadget_write_ics("homogeneous_nH3e2_%d_nside%d_boxsize_%d.dat.ic"%(res,nside,boxsize), data, format='hdf5' )
