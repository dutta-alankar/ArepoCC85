""" @package ./examples/shocktube_1d/check
Code that checks results of 1d shocktube problem

created by Rainer Weinberger, last modified: 19.02.2019
"""

""" load libraries """
import sys    # needed for exit codes
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os      # file specific calls
import matplotlib.pyplot as plt    ## needs to be active for plotting!
from scipy.interpolate import interp1d
from tqdm.auto import tqdm



cool_loc = "./cooltable_townsend.dat"
LAMBDA = np.loadtxt(cool_loc)
LAMBDA = interp1d(LAMBDA[:,0], LAMBDA[:,1])

mp = 1.6726e-24
kB = 1.3807e-16
Msun = 2e33
yr = 365*24*60**2
kpc = 3.086e21

X, Y, Z = 0.7154, 0.2703, 0.0143
Z *= 10**-0.2
X = 1.-(Y+Z)
mu = 1/(2*X+0.75*Y+0.5625*Z)
mue = 2./(1+X)
mui = 1/(1/mu-1/mue)
muH = 1/X

n = 0.1 #CGS

ne = n * (mu/mue)
ni = n * (mu/mui)
Tini = 1e7
gamma = 5/3.

UnitVelocity_in_cm_per_s    =              1e5
UnitLength_in_cm            =              3.086e18
UnitMass_in_g               =              4.91e31

tcool = n*kB*Tini/(gamma-1)/(ne*ni*LAMBDA(Tini))/(UnitLength_in_cm/UnitVelocity_in_cm_per_s)
print("mu %f   mue %f   mui%f   lam %e     tcool %f"%(mu, mue, mui, LAMBDA(Tini), tcool))

Dtype = np.float64  # double precision: np.float64, for single use np.float32
## open initial conditiions to get parameters

i_file = -1
pressure = []
temperature = []
time = []
cell = 5
while True:
    i_file += 1

    ## read in data
    directory = "./output/"
    filename = "snap_%03d.hdf5" % (i_file)
    print("try to open "+directory+filename)
    ## open hdf5 file
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break
    print("analyzing "+ filename)

    ## get data from snapshot
    time.append( np.float( data["Header"].attrs["Time"] ))
    position = np.array(data["PartType0"]["Coordinates"], dtype = np.float64)[cell]
    density = np.array(data["PartType0"]["Density"], dtype = np.float64)[cell]
    vel = np.array(data["PartType0"]["Velocities"], dtype = np.float64)[cell]
    internalEnergy = np.array(data["PartType0"]["InternalEnergy"], dtype = np.float64)[cell]
    ## convert to more useful data structure
    pressure.append( density*internalEnergy*(gamma-1)*((UnitMass_in_g/UnitLength_in_cm**3)*UnitVelocity_in_cm_per_s**2))
    temperature.append( internalEnergy*(gamma-1)*(mu*mp)/kB *UnitVelocity_in_cm_per_s**2)
    data.close()

pressure = np.array(pressure)
temperature = np.array(temperature)
time = np.array(time) * (UnitLength_in_cm/UnitVelocity_in_cm_per_s)

dPdt= np.gradient(pressure,time) #CGS
dT = np.gradient(temperature, time)

fig=plt.figure(figsize=(8,5))
factor = 1.0
plt.plot(temperature,-factor*dPdt/(ne*ni*(gamma-1)), 'o', c='tab:blue', label=r'$\rm \frac{-dP/dt}{(\gamma -1) n_e n_i}$', ms=3)
temperature = np.logspace(np.log10(np.min(temperature)), np.log10(np.max(temperature)), 300)
plt.plot(temperature,LAMBDA(temperature), c='tab:red', label=r'$\rm \Lambda (T)$', lw=5, alpha=0.5)
plt.legend(loc='lower right', prop={'size': 20})

plt.xlabel(r'Temperature (K)',fontsize=12)
plt.ylabel(r'Cooling function',fontsize=12)
plt.grid()
plt.yscale('log')
plt.xscale('log')
plt.ylim(ymin=3e-25, ymax=1.4e-21)
plt.savefig('LAMBDA_CHECK.png')
plt.show()
