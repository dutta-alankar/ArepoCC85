#!/bin/sh
module purge
module load gcc/10 openmpi/4 hdf5-mpi/1.12.0 gsl/2.4 fftw-mpi/3.3.10

export LD_LIBRARY_PATH="/mpcdf/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/lib:/mpcdf/soft/SLE_15/packages/skylake/gsl/gcc_10-10.3.0/2.4/lib:/mpcdf/soft/SLE_15/packages/skylake/fftw/gcc_10-10.3.0-openmpi_4-4.0.7/3.3.10/lib:/mpcdf/soft/SLE_15/packages/skylake/hdf5/gcc_10-10.3.0-openmpi_4-4.0.7/1.12.0/lib:$LD_LIBRARY_PATH"

source /freya/ptmp/mpa/adutt/ArepoCC85/.venv/bin/activate
python slice.py $1 $2
python slice-zoom.py $1 $2
