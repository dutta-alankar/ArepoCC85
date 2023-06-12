# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 01:04:40 2021

@author: alankar

Masses are interpreted as densities. 
Make sure corresponding Arepo flag is on in Config.sh
Takes in output from a steady wind simulation and introduces the cloud
"""

import numpy as np
import h5py
import sys
import healpy as hp
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
relaxedField = h5py.File('./mrelax/steady/snap_%03d.hdf5'%98, 'r')

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

Rcl  = 10 # pc
Rini = 300 # pc
chi  = 100

Rcl  = Rcl*pc/UnitLength_in_cm
Rini = Rini*pc/UnitLength_in_cm
RclBydcell = 16

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
data['vel'][:,:]  = np.array(relaxedField['/PartType0/Velocities'])
data['mass'][:]   = np.array(relaxedField['/PartType0/Density'])
data['u'][:]      = np.array(relaxedField['/PartType0/InternalEnergy'])

# Remove normal AREPO cells from volume inside the cloud
extra = 3*RclBydcell
layers = RclBydcell+extra

SpherePos = [0.5, 0.5, 0.5]                    # position of the sphere centre
SpherePos = np.array(SpherePos)*boxsize

SpherePos[1] = SpherePos[1] + Rini

rrCl      = np.sqrt((data['pos'][:,0] - SpherePos[0])**2 + \
                    (data['pos'][:,1] - SpherePos[1])**2 + \
                    (data['pos'][:,2] - SpherePos[2])**2)
# print(np.min(rrCl), np.max(rrCl))
flag    = data['flga']
ind     = np.where( rrCl <= (Rcl+extra*Rcl/RclBydcell) ) 
number  = np.size(ind)             # number of cells that need removal
print('Removed %d cell/s.'%number)
count_wnd = count - number

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

rrCl = np.delete(rrCl,ind)

#other quantities (do not change)
nside    = 3                # HEALPIX base resolution nside
ndim     = 3
npix     = hp.pixelfunc.nside2npix(nside)

# Add cells of the cloud
ipix   = np.arange(npix)
xx,yy,zz = hp.pixelfunc.pix2vec(nside, ipix) 

# generate array for the cloud
count = layers*npix
data_cloud = {}
data_cloud['pos']       = np.zeros((count,3), dtype=np.float32)
data_cloud['vel']       = np.zeros((count,3), dtype=np.float32)
data_cloud['mass']      = np.zeros(count, dtype=np.float32)
data_cloud['u']         = np.zeros(count, dtype=np.float32)
data_cloud['type']      = np.zeros(count, dtype=np.float32)
data_cloud['flga']      = np.zeros(count, dtype=int)
data_cloud['count']     = count

for layer in range(layers):
    SpherePos = [0.5, 0.5, 0.5]                    # position of the sphere centre
    SpherePos = np.array(SpherePos)*boxsize
    
    data_cloud['pos'][layer*npix:(layer+1)*npix,0]           = (layer*(Rcl/RclBydcell)+0.5*(Rcl/RclBydcell)) * xx + SpherePos[0]
    data_cloud['pos'][layer*npix:(layer+1)*npix,1]           = (layer*(Rcl/RclBydcell)+0.5*(Rcl/RclBydcell)) * yy + SpherePos[1] + Rini
    data_cloud['pos'][layer*npix:(layer+1)*npix,2]           = (layer*(Rcl/RclBydcell)+0.5*(Rcl/RclBydcell)) * zz + SpherePos[2]
    
    # print(layer)
    # print('x_avg = %.2f'%np.average(data_cloud['pos'][layer*npix:(layer+1)*npix,0]))
    # print('y_avg = %.2f'%np.average(data_cloud['pos'][layer*npix:(layer+1)*npix,1]))
    # print('z_avg = %.2f'%np.average(data_cloud['pos'][layer*npix:(layer+1)*npix,2]))
    
    xw = data_cloud['pos'][:,0][layer*npix:(layer+1)*npix] - SpherePos[0]
    yw = data_cloud['pos'][:,1][layer*npix:(layer+1)*npix] - SpherePos[1]
    zw = data_cloud['pos'][:,2][layer*npix:(layer+1)*npix] - SpherePos[2]
        
    rr_wind = np.sqrt(xw**2 + yw**2 + zw**2)
    
    # print('rr_avg = %.2f'%np.average(rr_wind))
    
    data_cloud['flga'][layer*npix:(layer+1)*npix]            = 0
    data_cloud['mass'][layer*npix:(layer+1)*npix]            = rho0*rhoTld(rr_wind/(Rinj/UnitLength_in_cm)) / UnitDensity_in_g_per_cm3
    data_cloud['u'][layer*npix:(layer+1)*npix]               = (1/(gamma-1))* ((P0*prsTld(rr_wind/(Rinj/UnitLength_in_cm)))/(rho0*rhoTld(rr_wind/(Rinj/UnitLength_in_cm)))) /  UnitU_in_cm_per_s_sq
    
    if layer<RclBydcell: # fill with cloud
        data_cloud['mass'][layer*npix:(layer+1)*npix]            = chi*data_cloud['mass'][layer*npix:(layer+1)*npix]
        data_cloud['u'][layer*npix:(layer+1)*npix]               = data_cloud['u'][layer*npix:(layer+1)*npix]/chi
        
    else: # fill with wind
        vel  = (v0*velTld(rr_wind/(Rinj/UnitLength_in_cm))) / UnitVelocity_in_cm_per_s
        
        sintheta = np.sqrt(xw**2 + yw**2)/rr_wind
        costheta = zw/rr_wind
        cosphi = xw/np.sqrt(xw**2 + yw**2)
        sinphi = yw/np.sqrt(xw**2 + yw**2)

        data_cloud['vel'][layer*npix:(layer+1)*npix,0]  = vel * sintheta * cosphi 
        data_cloud['vel'][layer*npix:(layer+1)*npix,1]  = vel * sintheta * sinphi 
        data_cloud['vel'][layer*npix:(layer+1)*npix,2]  = vel * costheta 
        
total_count = layers*npix + count_wnd

totalposx = np.hstack((newposx,data_cloud['pos'][:,0]))
totalposy = np.hstack((newposy,data_cloud['pos'][:,1]))
totalposz = np.hstack((newposz,data_cloud['pos'][:,2]))

totalvelx = np.hstack((newvelx,data_cloud['vel'][:,0]))
totalvely = np.hstack((newvely,data_cloud['vel'][:,1]))
totalvelz = np.hstack((newvelz,data_cloud['vel'][:,2]))

totalmass = np.hstack((newmass,data_cloud['mass']))
totaltype = np.hstack((newtype,data_cloud['type']))
totalu    = np.hstack((newu,data_cloud['u']))
totalflga = np.hstack((newflga,data_cloud['flga']))

print('total: %d==%d'%(total_count, totalposx.shape[0]))

# Update Data and Save
data['count']     = total_count
data['pos']       = np.zeros((data['count'],3))
data['vel']       = np.zeros((data['count'],3))
data['mass']      = np.zeros(data['count'])
data['u']         = np.zeros(data['count'])
data['type']      = np.zeros(data['count'])
data['flga']      = np.zeros(data['count'], dtype=int)

data['pos'][:,0]  = totalposx
data['pos'][:,1]  = totalposy
data['pos'][:,2]  = totalposz
data['vel'][:,0]  = totalvelx
data['vel'][:,1]  = totalvely
data['vel'][:,2]  = totalvelz
data['mass']      = totalmass
data['u']         = totalu
data['type']      = totaltype
data['flga']      = totalflga

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

fileName = "homogeneous_res_%d_nside%d_boxsize_%d_cloud.dat.ic"%(res,nside,boxsize)
gadget_write_ics(fileName, data, format='hdf5' )
write_xmf(fileName)
