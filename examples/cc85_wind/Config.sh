#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#--------------------------------------- Basic operation mode of code
NTYPES=6                        # number of particle types
#PERIODIC
#REFLECTIVE_X=2 #=2            # if set to 2, the boundary is inflow/outflow
#REFLECTIVE_Y=2 #=2
#REFLECTIVE_Z=2 #=2

#COOLING
GAMMA=1.66667
#ISOTHERM_EQS
#USE_ENTROPY_FOR_COLD_FLOWS
#ENTROPY_MACH_THRESHOLD=1.1
#INPUT_IN_DOUBLEPRECISION

#----------------------------------------MPI/Threading Hybrid
#NUM_THREADS=4                           # use OpenMP, with the given number of threads per MPI task
#IMPOSE_PINNING
#IMPOSE_PINNING_OVERRIDE_MODE
#GENERIC_ASYNC                           # enables asynchronous communication scheme

#--------------------------------------- Mesh Type
#AMR
#VORONOI

#--------------------------------------- Riemann solver
#VARIABLE_GAMMA
#RIEMANN_HLL
#RIEMANN_HLLC
#RIEMANN_ROSUNOV
#RIEMANN_HLLD
#RIEMANN_GAMMA

#AMR_CONNECTIONS
#AMR_GRADIENTS
#AMR_REDUCE_DOMAIN_DECOMPOISTION

#--------------------------------------- Reconstruction
#TVD_SLOPE_LIMITER
#TVD_SLOPE_LIMITER_VANLEER
#TVD_SLOPE_LIMITER_SUPERBEE
#TVD_SLOPE_LIMITER_ALBADA
#TVD_SLOPE_LIMITER_MINBEE
#TVD_SLOPE_LIMITER_MINMOD
#TVD_SLOPE_LIMITER_MC
#GRADIENT_LIMITER_DUFFELL
#DISABLE_TIME_EXTRAPOLATION              # use only when you know exactly what you are doing. activating this option will make your results wrong but can tell you about the behaviour of your code
#DISABLE_SPATIAL_EXTRAPOLATION           # use only when you know exactly what you are doing. activating this option will make your results wrong but can tell you about the behaviour of your code
#NO_SCALAR_GRADIENTS                     # disables time and spatial extrapolation for passive scalar fields
#GRADIENTS_GREEN_GAUSS                   # original (now depreciated) gradient estimate, reduced hydro scheme to first order

#--------------------------------------- Mesh motion and regularization
#VORONOI_STATIC_MESH
#VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION  # for VORONOI_STATIC_MESH force domain decomposition if there exist non-gas particles
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE
#REGULARIZE_MESH_LLOYD
#OUTPUT_MESH_FACE_ANGLE
#STICKY_POINTS_ON_REFLECTIVE_SURFACE     # if reflective boundaries are used, allows points to move only tangentially at boundary

#--------------------------------------- Time integration options
#FORCE_EQUAL_TIMESTEPS    # this chooses a variable but global timestep
TREE_BASED_TIMESTEPS     # non-local timestep criterion (take 'signal speed' into account)
#DECOUPLE_TIMESTEPS       # allows different timebins for gravity and hydro. use only WITHOUT FORCE_EQUAL_TIMESTEPS
#MUSCL_HANCOCK           # original (now depreciated) time integration scheme, only first order
#RUNGE_KUTTA_FULL_UPDATE

#--------------------------------------- Image generation
#VORONOI_MESHOUTPUT                      # 2D and 3D mesh output
#VORONOI_FREQUENT_IMAGES                 # creates images with frequency 'TimeBetweenImages' given in parameterfile, independent of snapshots

#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS
#REFINEMENT_SPLIT_MOST_DISTANCE_NEIGHBOUR
#REFINEMENT_MERGE_PAIRS
REFINEMENT_VOLUME_LIMIT
#REFINEMENT_HIGH_RES_GAS

#--------------------------------------- Mesh-relaxing or mesh-adding (this will not carry out a simulation)
#MESHRELAX                     # this keeps the mass constant and only regularizes the mesh
#MESHRELAX_DENSITY_IN_INPUT
#ADDBACKGROUNDGRID=16

#--------------------------------------- Things that are always recommended
#AUTO_SWAP_ENDIAN_READIC                # Enables automatic ENDIAN swapping for reading ICs
#CHUNKING                 # will calculated the gravity force in interleaved blocks. This can reduce imbalances in case multiple iterations due to insufficient buffer size need to be done

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
OUTPUT_IN_DOUBLEPRECISION                 # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
#OUTPUT_COORDINATES_IN_DOUBLEPRECISION    # will always output coordinates in double precision
#NGB_TREE_DOUBLEPRECISION                 # if this is enabled, double precision is used for the neighbor node extension

#-------------------------------------------- Things for special behaviour
#VORONOI_DYNAMIC_UPDATE          # keeps track of mesh connectivity, which speeds up mesh construction
VORONOI_MESH_KEEP_DT_AND_DTC    # keeps DTC and DT in memory, i.e. for anisotropic transport solvers
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#MPI_HYPERCUBE_ALLGATHERV       # some MPI-libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPISENDRECV_CHECKSUM
#NOTREERND
#BITS_PER_DIMENSION=42      # Peano-Hilbert order
OVERRIDE_PEANOGRID_WARNING


#--------------------------------------- Output/Input options
#UPDATE_GRADIENTS_FOR_OUTPUT
#OUTPUT_PRESSURE_GRADIENT
#OUTPUT_DENSITY_GRADIENT
#OUTPUT_VELOCITY_GRADIENT
#OUTPUT_BFIELD_GRADIENT
#OUTPUT_VERTEX_VELOCITY
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE
#OUTPUT_VOLUME
#OUTPUT_CENTER_OF_MASS
#OUTPUT_SURFACE_AREA
OUTPUT_PRESSURE
#OUTPUTPOTENTIAL
#OUTPUTACCELERATION
#OUTPUTTIMESTEP
#OUTPUT_SOFTENINGS            # output particle softenings
#OUTPUTGRAVINTERACTIONS       # output gravitatational interactions (from the tree) of particles
HAVE_HDF5                     # needed when HDF5 I/O support is desired
HDF5_FILTERS                  # activate snapshot compression and checksum for HDF5 output
OUTPUT_XDMF                   #writes an .xmf file for each snapshot, which can be read by visit (with the hdf5 snapshot)
#OUTPUTCOOLRATE                # outputs cooling rate, and conduction rate if enabled
#OUTPUT_DIVVEL                 # output  velocity divergence
#OUTPUT_CURLVEL                 # output  velocity curl
#OUTPUT_COOLHEAT               # output actual energy loss/gain in cooling/heating routine
OUTPUT_VORTICITY
#OUTPUT_CELL_SPIN
#MEASURE_DISSIPATION_RATE      # measures and outputs dissipation rate. Note: requires USE_ENTROPY_FOR_COLD_FLOWS, even though it will then always use the thermal energy update
#OUTPUT_MACHNUM                # output maximum mach number of a cell
#OUTPUT_TASK
#OUTPUT_ENTROPY
#OUTPUT_CSND
READ_MASS_AS_DENSITY_IN_INPUT  # Reads the mass field in the IC as density

#--------------------------------------- Testing and Debugging options
DEBUG                         # enables core-dumps
#VERBOSE                       # reports readjustments of buffer sizes
HOST_MEMORY_REPORTING         # reports after start-up the available system memory by analyzing /proc/meminfo

#--------------------------------------- Tracer
PASSIVE_SCALARS=1

#-------------------------------------- Navier-Stokes Terms
#GLOBAL_VISCOSITY 	   #needs dynamic and bulk coefficients
#USE_KINEMATIC_VISCOSITY    #needs only one input parameter
#ALPHA_VISCOSITY=2            #for accretion disks
#LOCAL_VISCOSITY=1          #=1 Sutherland viscosity/ =2 Spitzer viscosity
#THERMAL_CONDUCTION
#TRACER_DIFFUSION           #requires TRACER_FIELD switched on

#-------------------------------------- Special Boundaries within domain
#SPECIAL_BOUNDARY          #Main Switch
#BOUNDARY_FLAG
AGNWIND_FLAG
REFINE_AGNWIND
#CLOUD_PRESENT
#AGNWIND_DUTYCYCLE
#BRITNI_ACCRETION
#STICKYFLAGS
#OUTPUT_STICKYFLAGS
