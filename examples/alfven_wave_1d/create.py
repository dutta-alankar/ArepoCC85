""" @package examples/alfven_wave_1d/create.py
Code that creates 1d Alfven wave test problem;
supposed to be as simple and self-contained as possible

created by Alessandro Stenghel and Federico Marinacci,
last modified 13.7.2020 -- comments welcome
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/alfven_wave_1d/create.py: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/ics.hdf5'

## ensure calculations happen with predefined precision
FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32 # integer type

## computational domain
Boxsize = FloatType(1.0)
if len(sys.argv) > 3:
  NumberOfCells = IntType(sys.argv[3])
else:
  NumberOfCells = IntType(32)

## initial state
density_0 = FloatType(1.0)
velocity_0 = FloatType(0.0)
pressure_0 = FloatType(1.0)
gamma = FloatType(5.0) / FloatType(3.0)
gamma_minus_one = gamma - FloatType(1.0)
delta = FloatType(1e-6)    # relative velocity perturbation
uthermal_0 = pressure_0 / density_0 / gamma_minus_one
bfield_0= FloatType(1.0)
k_z = 2*np.pi
omega = bfield_0*k_z/np.sqrt(density_0)

""" set up grid: uniform 1d grid """
## spacing
dx = Boxsize / FloatType(NumberOfCells)

## position of first and last cell
pos_first, pos_last = FloatType(0.5) * dx, Boxsize - FloatType(0.5) * dx

## set up grid
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = np.linspace(pos_first, pos_last, NumberOfCells, dtype=FloatType)
Volume = np.full(NumberOfCells, dx, dtype=FloatType)

""" set up magnetohydrodynamical quantitites """
## set up unperturbed system; density, velocity and specific internal energy
Density = np.full(Pos.shape[0], density_0, dtype=FloatType)
Velocity = np.zeros(Pos.shape, dtype=FloatType)
Uthermal = np.full(Pos.shape[0], uthermal_0, dtype=FloatType)
Bfield = np.zeros(Pos.shape, dtype=FloatType)

## perturbations
Velocity[:,0] = velocity_0
Velocity[:,1] = delta*np.sin(k_z*Pos[:,0])
Velocity[:,2] = delta*np.cos(k_z*Pos[:,0])
Uthermal *= (Density / density_0)**gamma_minus_one
Bfield[:,0] = bfield_0
Bfield[:,1] = -k_z*bfield_0/omega*Velocity[:,1]
Bfield[:,2] = -k_z*bfield_0/omega*Velocity[:,2]

## mass instead of density needed for input
Mass = Density * Volume

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")  # create particle group for gas cells

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
part0.create_dataset("MagneticField", data=Bfield*np.sqrt(4*np.pi)) #conversion Lorentz System->Gaussian System

## close file
IC.close()

""" normal exit """
sys.exit(0)
