# Using Anaconda

IMPACT-T is available through conda-forge and can be installed via:
```bash
conda create -n impact
source activate impact # or conda activate impact
# For non-MPI
conda install -c conda-forge impact-t

# For OpenMPI
conda install -c conda-forge impact-t=*=mpi_openmpi*

# For MPICH
conda install -c conda-forge impact-t=*=mpi_mpich*
```
# Compiling The Code

## Unix

### Single Processor Code:

```shell script
mkdir build
cd build
cmake ../src
make
make install
```

### Multi Processor Code:

```shell script
mkdir build
cd build
cmake ../src -DUSE_MPI=ON
make
make install
```

## Windows

### Single Processor Code:

```shell script
mkdir build
cd build
cmake -G "NMake Makefiles" ../src
cmake --build .
cmake --install
```

### Multi Processor Code:

**Not Tested**


## Testing

```shell script
cd examples
cd Sample1
# Execute non-MPI
ImpactTexe

# Execute MPI
mpirun -n <cores> ImpactTexe-mpi
```

# Using at NERSC

If you decide to use IMPACT-T at NERSC you can do so via Anaconda or building from source.
The instructions for Anaconda are the same as above.

## Compiling the code

### For Haswell
```bash
module load openmpi # if using OpenMPI otherwise skip for MPICH
mkdir build
cd build
cmake ../src -DUSE_MPI=ON
```

### For KNL
```bash
module swap craype-haswell craype-mic-knl
module load openmpi # if using OpenMPI otherwise skip for MPICH
mkdir build
cd build
cmake ../src -DUSE_MPI=ON
```

## Running at NERSC

There is one caveat with the conda-forge installed MPI version of the code.
Instead of running with `srun -n <cores> ` you must use `mpirun -n <cores>`.

