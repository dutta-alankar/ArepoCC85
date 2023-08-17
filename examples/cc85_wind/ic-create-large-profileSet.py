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
#from mpl_toolkits import mplot3d
import sys
from gadget import gadget_write_ics
from scipy.interpolate import interp1d

#sys.path.append('./arepo-snap-util')

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

data = np.loadtxt('CC85_steady-prof_gamma_%.3f.txt'%gamma)

rhoTld = interp1d(data[:,0], data[:,1], fill_value='extrapolate')
prsTld = interp1d(data[:,0], data[:,2], fill_value='extrapolate')
velTld = interp1d(data[:,0], data[:,3], fill_value='extrapolate')

#units and constants
UnitVelocity_in_cm_per_s = 1.e5
UnitLength_in_cm         = pc
UnitMass_in_g            = 4.924e31

UnitTime_in_s = UnitLength_in_cm/UnitVelocity_in_cm_per_s
UnitU_in_cm_per_s_sq   = UnitVelocity_in_cm_per_s**2
UnitDensity_in_g_per_cm3 = UnitMass_in_g/UnitLength_in_cm**3

#free parameters
boxsize   = 50000                              # pc
boxhres   = 10000                              # pc
lrpoints  = 10                                 # low resolution fill of the cell
res       = 512                                # roughly res cels per side at high resolution
nside     = 6                                  # HEALPIX base resolution nside
bndryR    = Rinj/pc                            # radius of boundary
size_fac  = 1.                                 # thickness of boundary layers in units of half the grid size
SpherePos = [0.5, 0.5, 0.5]                    # position of the sphere centre
SpherePos = np.array(SpherePos)*boxsize
pert      = 0.1                                # magnitude of mesh perturbation in terms of dx
contrast  = 10000.0

#other quantities (do not change)
ndim     = 3
ncell    = res**ndim
npix     = hp.pixelfunc.nside2npix(nside)

#generate arrays
data = {}
data['pos']       = np.zeros((ncell + 2*npix,3), dtype=np.float32)
data['vel']       = np.zeros((ncell + 2*npix,3), dtype=np.float32)
data['mass']      = np.zeros(ncell  + 2*npix, dtype=np.float32)
data['u']         = np.zeros(ncell  + 2*npix, dtype=np.float32)
data['type']      = np.zeros(ncell  + 2*npix, dtype=np.float32)
data['flga']      = np.zeros(ncell  + 2*npix, dtype=int)
data['count']     = ncell + 2*npix

#Fill the entire box with normal AREPO cells
npoints   = int(np.round(ncell**(1./ndim)))+1

edgesx = np.linspace( Rinj/pc + SpherePos[0], boxhres/2 + SpherePos[0], endpoint=True, num=npoints//2, dtype=np.float32 )
edgesx = np.hstack( (np.linspace(-boxhres/2 + SpherePos[0], -Rinj/pc + SpherePos[0], endpoint=True, num=npoints-len(edgesx), dtype=np.float32),edgesx) )
dgridx = edgesx[1]-edgesx[0]
cell_posx = 0.5 * (edgesx[1:]+edgesx[:-1])

edgesy = np.linspace( Rinj/pc + SpherePos[1], boxhres/2 + SpherePos[1], endpoint=True, num=npoints//2, dtype=np.float32 )
edgesy = np.hstack( (np.linspace(-boxhres/2 + SpherePos[1], -Rinj/pc + SpherePos[1], endpoint=True, num=npoints-len(edgesy), dtype=np.float32),edgesy) )
dgridy = edgesy[1]-edgesy[0]
cell_posy = 0.5 * (edgesy[1:]+edgesy[:-1])

edgesz = np.linspace( Rinj/pc + SpherePos[2], boxhres/2 + SpherePos[2], endpoint=True, num=npoints//2, dtype=np.float32 )
edgesz = np.hstack( (np.linspace(-boxhres/2 + SpherePos[2], -Rinj/pc + SpherePos[2], endpoint=True, num=npoints-len(edgesz), dtype=np.float32),edgesz) )
dgridz = edgesz[1]-edgesz[0]
cell_posz = 0.5 * (edgesz[1:]+edgesz[:-1])

volbndry = (size_fac/2.)**3 * dgridx*dgridy*dgridz
dx       = dgridx
vol      = dx**3
cells_x, cells_y, cells_z  = np.meshgrid(cell_posx,cell_posy,cell_posz, indexing='ij')
print("Length Points: ", len(edgesx))
print("Cell edges: ", edgesx)

#Fill the entire box with normal AREPO cells highres
data['pos'][:ncell,0] = cells_x.flatten()
data['pos'][:ncell,1] = cells_y.flatten()
data['pos'][:ncell,2] = cells_z.flatten()

# Add boundary cells
ipix   = np.arange(npix)
xx,yy,zz = hp.pixelfunc.pix2vec(nside, ipix)
closeness = 4

dgrid = min(dgridx, dgridy, dgridz)/closeness
# Inner Layer
print("Inner layer: %.2f"%(bndryR-size_fac*dgrid/2.))
data['pos'][ncell:ncell+npix,0]           = (bndryR-size_fac*dgrid/2.) * xx + SpherePos[0]
data['pos'][ncell:ncell+npix,1]           = (bndryR-size_fac*dgrid/2.) * yy + SpherePos[1]
data['pos'][ncell:ncell+npix,2]           = (bndryR-size_fac*dgrid/2.) * zz + SpherePos[2]
data['flga'][ncell:ncell+npix]            = 1
# data['mass'][ncell:ncell+npix]            = rho0 / UnitDensity_in_g_per_cm3
# data['u'][ncell:ncell+npix]               = (1/(gamma-1))* (P0/rho0) /  UnitU_in_cm_per_s_sq

#Outer Layer
print("Outer layer: %.2f"%(bndryR+size_fac*dgrid/2.))
data['pos'][ncell+npix:ncell+2*npix,0]    = (bndryR+size_fac*dgrid/8.) * xx + SpherePos[0]
data['pos'][ncell+npix:ncell+2*npix,1]    = (bndryR+size_fac*dgrid/8.) * yy + SpherePos[1]
data['pos'][ncell+npix:ncell+2*npix,2]    = (bndryR+size_fac*dgrid/8.) * zz + SpherePos[2]
data['flga'][ncell+npix:ncell+2*npix]     = 2
# data['mass'][ncell+npix:ncell+2*npix]     = rho0 / UnitDensity_in_g_per_cm3
# data['u'][ncell+npix:ncell+2*npix]        = (1/(gamma-1))* (P0/rho0) /  UnitU_in_cm_per_s_sq

'''
fig = plt.figure(figsize = (10, 10))
ax = plt.axes(projection ="3d")
ax.scatter3D((bndryR-size_fac*dgrid/2.) * xx + SpherePos[0], (bndryR-size_fac*dgrid/2.) * yy + SpherePos[1], (bndryR-size_fac*dgrid/2.) * zz + SpherePos[2])
plt.show()
'''

# Perturb Cartesian grid
theta  = (np.random.rand(ncell).astype(np.float32)) * np.pi
phi    = (np.random.rand(ncell).astype(np.float32)) * np.pi * 2

pert_x = np.sin(theta) * np.cos(phi)
pert_y = np.sin(theta) * np.sin(phi)
pert_z = np.cos(theta)

data['pos'][:ncell,0]  += pert * dx * pert_x
data['pos'][:ncell,1]  += pert * dx * pert_y
data['pos'][:ncell,2]  += pert * dx * pert_z
# data['mass'][:ncell]    = 1e-4*rho0 / UnitDensity_in_g_per_cm3
# data['mass'][data['flga']>0]  = rho / unit_rho * 1.e-3
# data['u'][:ncell]       = (1/(gamma-1))* (P0/rho0) /  UnitU_in_cm_per_s_sq

rr      = np.sqrt((data['pos'][:,0] - SpherePos[0])**2 + \
                  (data['pos'][:,1] - SpherePos[1])**2 + \
                  (data['pos'][:,2] - SpherePos[2])**2)

data['mass'][:] = rho0*rhoTld(rr/(Rinj/pc)) / UnitDensity_in_g_per_cm3
data['u'][:]    = (1/(gamma-1))* ((P0*prsTld(rr/(Rinj/pc)))/(rho0*rhoTld(rr/(Rinj/pc)))) /  UnitU_in_cm_per_s_sq
vel = (v0*velTld(rr/(Rinj/pc))) / UnitVelocity_in_cm_per_s

data['vel'][:,0]  = vel * (data['pos'][:,0] - SpherePos[0])/rr
data['vel'][:,1]  = vel * (data['pos'][:,1] - SpherePos[1])/rr
data['vel'][:,2]  = vel * (data['pos'][:,2] - SpherePos[2])/rr

# Remove normal AREPO cells from volume bounded by boundary sphere
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

# Update Data
data['count']     = data['count'] - number
data['pos']       = np.zeros((data['count'],3))
data['vel']       = np.zeros((data['count'],3))
data['mass']      = np.zeros(data['count'])
data['u']         = np.zeros(data['count'])
data['type']      = np.zeros(data['count'])
data['flga']      = np.zeros(data['count'], dtype=int)

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

# Fill Outer region with low resolution
# Assumes uniform box
ncell_lowres = int(lrpoints**ndim)
lowres_data = {}
lowres_data['pos']       = np.zeros((ncell_lowres,3), dtype=np.float32)
lowres_data['vel']       = np.zeros((ncell_lowres,3), dtype=np.float32)
lowres_data['mass']      = np.zeros(ncell_lowres, dtype=np.float32)
lowres_data['u']         = np.zeros(ncell_lowres, dtype=np.float32)
lowres_data['type']      = np.zeros(ncell_lowres, dtype=np.float32)
lowres_data['flga']      = np.zeros(ncell_lowres, dtype=int)
lowres_data['count']     = ncell_lowres

npoints   = int(np.round(ncell_lowres**(1./ndim)))+1
edges = np.linspace( -boxsize/2 + SpherePos[0], boxsize/2 + SpherePos[0], endpoint=True, num=npoints, dtype=np.float32 )
cell_pos = 0.5 * (edges[1:]+edges[:-1])
cells_x, cells_y, cells_z  = np.meshgrid(cell_pos,cell_pos,cell_pos, indexing='ij')
print("Length Points: ", len(edges))
print("Cell edges: ", edges)

# Fill the entire box with normal AREPO cells lowres
lowres_data['pos'][:ncell_lowres,0] = cells_x.flatten()
lowres_data['pos'][:ncell_lowres,1] = cells_y.flatten()
lowres_data['pos'][:ncell_lowres,2] = cells_z.flatten()

# Perturb Cartesian grid
theta  = (np.random.rand(ncell_lowres).astype(np.float32)) * np.pi
phi    = (np.random.rand(ncell_lowres).astype(np.float32)) * np.pi * 2

pert_x = np.sin(theta) * np.cos(phi)
pert_y = np.sin(theta) * np.sin(phi)
pert_z = np.cos(theta)

dx = edges[1] - edges[0]
pert     = 0.001                                # magnitude of mesh perturbation in terms of dx
lowres_data['pos'][:ncell_lowres,0]  += pert * dx * pert_x
lowres_data['pos'][:ncell_lowres,1]  += pert * dx * pert_y
lowres_data['pos'][:ncell_lowres,2]  += pert * dx * pert_z
# lowres_data['mass'][:ncell_lowres]    = 1e-4*rho0 / UnitDensity_in_g_per_cm3
# lowres_data['u'][:ncell_lowres]       = (1/(gamma-1))* (P0/rho0) /  UnitU_in_cm_per_s_sq

rr      = np.sqrt((lowres_data['pos'][:ncell_lowres,0] - SpherePos[0])**2 + \
                  (lowres_data['pos'][:ncell_lowres,1] - SpherePos[1])**2 + \
                  (lowres_data['pos'][:ncell_lowres,2] - SpherePos[2])**2)


lowres_data['mass'][:ncell_lowres] = rho0*rhoTld(rr/(Rinj/pc)) / UnitDensity_in_g_per_cm3
lowres_data['u'][:ncell_lowres]    = (1/(gamma-1))* ((P0*prsTld(rr/(Rinj/pc)))/(rho0*rhoTld(rr/(Rinj/pc)))) /  UnitU_in_cm_per_s_sq
vel  = (v0*velTld(rr/(Rinj/pc))) / UnitVelocity_in_cm_per_s

data['vel'][:ncell_lowres,0]  = vel * (lowres_data['pos'][:ncell_lowres,0] - SpherePos[0])/rr
data['vel'][:ncell_lowres,1]  = vel * (lowres_data['pos'][:ncell_lowres,1] - SpherePos[1])/rr
data['vel'][:ncell_lowres,2]  = vel * (lowres_data['pos'][:ncell_lowres,2] - SpherePos[2])/rr

# Remove lowres cells from highres region
flag    = lowres_data['flga']
ind     = np.where(
                np.logical_and(np.fabs(cells_x.flatten())<boxhres/2,
                               np.fabs(cells_y.flatten())<boxhres/2,
                               np.fabs(cells_z.flatten())<boxhres/2 ) )
number  = np.size(ind)             # number of cells that need removal

newposx = np.hstack((data['pos'][:,0],np.delete(lowres_data['pos'][:,0],ind)))
newposy = np.hstack((data['pos'][:,1],np.delete(lowres_data['pos'][:,1],ind)))
newposz = np.hstack((data['pos'][:,2],np.delete(lowres_data['pos'][:,2],ind)))

newvelx = np.hstack((data['vel'][:,0],np.delete(lowres_data['vel'][:,0],ind)))
newvely = np.hstack((data['vel'][:,1],np.delete(lowres_data['vel'][:,1],ind)))
newvelz = np.hstack((data['vel'][:,2],np.delete(lowres_data['vel'][:,2],ind)))

newmass = np.hstack((data['mass'],np.delete(lowres_data['mass'],ind)))
newtype = np.hstack((data['type'],np.delete(lowres_data['type'],ind)))
newu    = np.hstack((data['u'],np.delete(lowres_data['u'],ind)))
newflga = np.hstack((data['flga'],np.delete(lowres_data['flga'],ind)))

# Update Data and Save
ncell_lowres = ncell_lowres - number
data['count']     = data['count'] + ncell_lowres
data['pos']       = np.zeros((data['count'],3))
data['vel']       = np.zeros((data['count'],3))
data['mass']      = np.zeros(data['count'])
data['u']         = np.zeros(data['count'])
data['type']      = np.zeros(data['count'])
data['flga']      = np.zeros(data['count'], dtype=int)

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

def write_xmf(fileName):
    import h5py
    with open(fileName+'.xmf','w') as f:
        f.write('<?xml version=\"1.0\" ?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n')
        f.write('<Xdmf Version=\"2.0\">\n')
        f.write(' <Domain>')
        hdf = h5py.File(fileName+'.hdf5','r')
        numParts = 1 #for the time being, only gas particles are there
        for part in range(numParts):
            prec = 8
            f.write( '  <Grid Name=\"PartType%d\" GridType=\"Uniform\">\n'%part)
            f.write('   <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n'%len(np.array(hdf['PartType%d/Masses'%part])))
            f.write('   <Geometry GeometryType=\"XYZ\">\n')
            f.write('    <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n'%(len(np.array(hdf['PartType%d/Masses'%part])), prec))
            f.write('     %s:/PartType%d/Coordinates\n'%('./'+fileName+'.hdf5',part))
            f.write('    </DataItem>\n')
            f.write('   </Geometry>\n')

            f.write('   <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Node\">\n'%'Velocities')
            f.write('    <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n'%(len(np.array(hdf['PartType%d/Masses'%part])), prec))
            f.write('     %s:/PartType%d/%s\n'%('./'+fileName+'.hdf5',part,'Velocities'))
            f.write('    </DataItem>\n')
            f.write('   </Attribute>\n')

            f.write('   <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n'%'Masses')
            f.write('    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n'%(len(np.array(hdf['PartType%d/Masses'%part])), prec))
            f.write('     %s:/PartType%d/%s\n'%('./'+fileName+'.hdf5',part,'Masses'))
            f.write('    </DataItem>\n')
            f.write('   </Attribute>\n')

            f.write('   <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n'%'InternalEnergy')
            f.write('    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n'%(len(np.array(hdf['PartType%d/Masses'%part])), prec))
            f.write('     %s:/PartType%d/%s\n'%('./'+fileName+'.hdf5',part,'InternalEnergy'))
            f.write('    </DataItem>\n')
            f.write('   </Attribute>\n')

            f.write('   <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n'%'AGNFlag')
            f.write('    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n'%(len(np.array(hdf['PartType%d/Masses'%part])), prec))
            f.write('     %s:/PartType%d/%s\n'%('./'+fileName+'.hdf5',part,'AGNFlag'))
            f.write('    </DataItem>\n')
            f.write('   </Attribute>\n')

            f.write('   <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n'%'ParticleIDs')
            f.write('    <DataItem Dimensions=\"%d\" NumberType=\"Integer\" Precision=\"%d\" Format=\"HDF\">\n'%(len(np.array(hdf['PartType%d/Masses'%part])), prec))
            f.write('     %s:/PartType%d/%s\n'%('./'+fileName+'.hdf5',part,'ParticleIDs'))
            f.write('    </DataItem>\n')
            f.write('   </Attribute>\n')

            f.write('  </Grid>\n')
            f.write(' </Domain>\n')
            f.write('</Xdmf>')
            hdf.close()

fileName = "homogeneous_res_%d_nside%d_boxsize_%d.dat.ic"%(res,nside,boxsize)
print('ReferenceGasPartMass: ', rho0 / UnitDensity_in_g_per_cm3 *(boxhres/res)**3)
gadget_write_ics(fileName, data, format='hdf5' )
write_xmf(fileName)
