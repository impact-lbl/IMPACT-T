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
After these steps, the IMPACT-T executable `ImpactTexe` or `ImpactTexe-mpi`, respectively, will be in your [PATH](https://en.wikipedia.org/wiki/PATH_(variable)) environment variable and is thus ready to use like any regular command-line command.

# Compiling The Code

If you are new to CMake, [this short tutorial](https://hsf-training.github.io/hsf-training-cmake-webpage/) from the HEP Software foundation is the perfect place to get started with it.

If you just want to use CMake to build the project, jump into sections *1. Introduction*, *2. Building with CMake* and *9. Finding Packages*.

## Unix

### Single Processor Code:

```shell script
# inside the IMPACT-T src/ directory:
cmake -S . -B build
cmake --build build
# the executable in now in build/bin/

# this command needs sudo if you install into system paths:
cmake --build build --target install
```

### Multi Processor Code:

```shell script
# inside the IMPACT-T src/ directory:
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
cmake --build build --target install
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
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
# find the executable in build/bin/
```

### For KNL
```bash
module swap craype-haswell craype-mic-knl
module load openmpi # if using OpenMPI otherwise skip for MPICH
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
# find the executable in build/bin/
```

## Running at NERSC

There is one caveat with the conda-forge installed MPI version of the code.
Instead of running with `srun -n <cores> ` you must use `mpirun -n <cores>`.
