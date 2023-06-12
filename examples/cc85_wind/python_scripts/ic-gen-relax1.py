# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 01:04:40 2021

@author: alankar

Masses are interpreted as densities. 
Make sure corresponding Arepo flag is on in Config.sh

Read the relaxed grid and re-fill the field values from CC85 solution.
"""

import numpy as np
import h5py
import sys
from gadget import gadget_write_ics
from scipy.interpolate import interp1d

# sys.path.append('./arepo-snap-util')

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

# units and constants
UnitVelocity_in_cm_per_s = 1.e5
UnitLength_in_cm         = pc
UnitMass_in_g            = 4.924e31

UnitTime_in_s = UnitLength_in_cm/UnitVelocity_in_cm_per_s
UnitU_in_cm_per_s_sq   = UnitVelocity_in_cm_per_s**2
UnitDensity_in_g_per_cm3 = UnitMass_in_g/UnitLength_in_cm**3

# mesh relaxed data
relaxedField = h5py.File('./mrelax/snap_%03d.hdf5'%96, 'r')

# free parameters
boxsize   = 50000                              # pc
boxhres   = 10000                              # pc
lrpoints  = 10                                 # low resolution fill of the cell
res       = 128                                # roughly res cels per side at high resolution
nside     = 6                                  # HEALPIX base resolution nside
bndryR    = Rinj/pc                            # radius of boundary
size_fac  = 1.                                 # thickness of boundary layers in units of half the grid size
SpherePos = [0.5, 0.5, 0.5]                    # position of the sphere centre
SpherePos = np.array(SpherePos)*boxsize
pert      = 0.1                                # magnitude of mesh perturbation in terms of dx
contrast  = 10000.0

# generate array
count = np.array(relaxedField['/PartType0/AGNFlag']).shape[0]
data = {}
data['pos']       = np.zeros((count,3), dtype=np.float32)
data['vel']       = np.zeros((count,3), dtype=np.float32)
data['mass']      = np.zeros(count, dtype=np.float32)
data['u']         = np.zeros(count, dtype=np.float32)
data['type']      = np.zeros(count, dtype=np.float32)
data['flga']      = np.zeros(count, dtype=int)
data['count']     = count

# assign the solution
data['pos'][:,:]  = np.array(relaxedField['/PartType0/Coordinates'])
data['flga'][:]   = np.array(relaxedField['/PartType0/AGNFlag'])

rr      = np.sqrt((data['pos'][:,0] - SpherePos[0])**2 + \
                  (data['pos'][:,1] - SpherePos[1])**2 + \
                  (data['pos'][:,2] - SpherePos[2])**2)

data['mass'][:] = rho0*rhoTld(rr/(Rinj/pc)) / UnitDensity_in_g_per_cm3
data['u'][:]    = (1/(gamma-1))* ((P0*prsTld(rr/(Rinj/pc)))/(rho0*rhoTld(rr/(Rinj/pc)))) /  UnitU_in_cm_per_s_sq
vel = (v0*velTld(rr/(Rinj/pc))) / UnitVelocity_in_cm_per_s
sintheta = np.sqrt(data['pos'][:,0]**2 + data['pos'][:,1]**2)/np.sqrt(data['pos'][:,0]**2 + data['pos'][:,1]**2 + data['pos'][:,2]**2)
costheta = data['pos'][:,2]/np.sqrt(data['pos'][:,0]**2 + data['pos'][:,1]**2 + data['pos'][:,2]**2)
cosphi = data['pos'][:,0]/np.sqrt(data['pos'][:,0]**2 + data['pos'][:,1]**2)
sinphi = data['pos'][:,1]/np.sqrt(data['pos'][:,0]**2 + data['pos'][:,1]**2)

data['vel'][:,0]  = vel * sintheta * cosphi 
data['vel'][:,1]  = vel * sintheta * sinphi 
data['vel'][:,2]  = vel * costheta 

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

fileName = "homogeneous_res_%d_nside%d_boxsize_%d_relaxed.dat.ic"%(res,nside,boxsize)
gadget_write_ics(fileName, data, format='hdf5' )
write_xmf(fileName)
