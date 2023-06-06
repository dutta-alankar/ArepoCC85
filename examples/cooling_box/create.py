""" @package examples/shocktube_1d/create.py
Code that creates 1d shocktube problem initial conditions
supposed to be as siple as possible and self-contained

created by Rainer Weinberger, 19.02.2019
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format

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

UnitVelocity_in_cm_per_s    =              1e5
UnitLength_in_cm            =              3.086e18
UnitMass_in_g               =              4.91e31

simulation_directory = str(sys.argv[1])
print("examples/scool-test/create.py: creating ICs in directory " + simulation_directory)

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(1.0)
NumberOfCells = IntType(10)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

## Cooling test
n = 0.1 #CGS
density_0 = FloatType(n*mu*mp)
velocity_0 = FloatType(0.0)
temperature_0 = FloatType(1e7)
pressure_0 = FloatType(density_0/(mu*mp)*kB*temperature_0)

gamma = FloatType(5./3)  ## note: this has to be consistent with the parameter settings for Arepo!
gamma_minus1 = gamma - FloatType(1.0)
uthermal_0 = pressure_0 / gamma_minus1 / density_0

""" set up grid: uniform 1d grid """
## spacing
dx = Boxsize / FloatType(NumberOfCells)
## position of first and last cell
pos_first, pos_last = FloatType(0.5) * dx, Boxsize - FloatType(0.5) * dx

## set up grid
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = np.linspace(pos_first, pos_last, NumberOfCells, dtype=FloatType)
Volume = np.full(NumberOfCells, dx, dtype=FloatType)

""" set up hydrodynamical quantitites """
## Gas state
Density = np.full(Pos.shape[0], density_0, dtype=FloatType) / (UnitMass_in_g/UnitLength_in_cm**3)
Velocity = np.zeros(Pos.shape, dtype=FloatType)
Velocity[:,0] = velocity_0 / UnitVelocity_in_cm_per_s
Uthermal = np.full(Pos.shape[0], uthermal_0, dtype=FloatType) / UnitVelocity_in_cm_per_s**2


## mass instead of density needed for input
Mass = Density * Volume

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")    # create particle group for gas cells

## write header entries
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype=IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType) )
header.attrs.create("MassTable", np.zeros(6, dtype=IntType) )
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", Boxsize)
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
header.attrs.create("UnitVelocity_in_cm_per_s", UnitVelocity_in_cm_per_s)
header.attrs.create("UnitLength_in_cm", UnitLength_in_cm)
header.attrs.create("UnitMass_in_g", UnitMass_in_g)
if Pos.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

## write cell data
part0.create_dataset("ParticleIDs", data=np.arange(1, NumberOfCells+1) )
part0.create_dataset("Coordinates", data=Pos)
part0.create_dataset("Masses", data=Mass)
part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)

## close file
IC.close()

""" normal exit """
sys.exit(0)
