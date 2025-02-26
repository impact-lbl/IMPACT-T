
# IMPACT-T

IMPACT-T is a fully three-dimensional program to track relativistic charged particles taking into account space charge forces, short-range longitudinal and transverse wakefields, coherent synchrotron radiation (CSR) wakefield in accelerators. IMPACT-T code can run on both massive parallel supercomputers and single processor computers such as Windows PC, Mac, and Linux system. It is one of the few codes used in the photoinjector community that has a parallel implementation, making it very useful for high statistics simulations of beam halos and beam diagnostics. It has a comprehensive set of beamline elements, and furthermore allows arbitrary overlap of their fields, which gives the IMPACT-T a capability to model both the standing wave structure and traveling wave structure. It includes mean-field space-charge solvers based on an integrated Green function to efficiently and accurately treat beams with large aspect ratio, and a shifted Green function to efficiently treat image charge effects of a cathode. It is also unique in its inclusion of energy binning in the space-charge calculation to model beams with large energy spread. It also has a direct N-body solver to calculate stochastic space-charge forces. IMPACT-T has a flexible data structure that allows particles to be stored in containers with common characteristics; for photoinjector simulations the containers represent multiple slices, but in other applications they could correspond, e.g., to particles of different species. Together, all these features make IMPACT-T a powerful and versatile tool for modeling beams in photoinjectors and other systems.

Main contact: Ji Qiang (jqiang@lbl.gov), Lawrence Berkeley National Laboratory

Citation:
J. Qiang, S. Lidia, R. D. Ryne, C. Limborg-Deprey, "A Three-Dimensional Quasi-Static Model for High Brightnees Beam Dynamics Simulation",
Phys. Rev. Special Topics - Accel. Beams 9, 044204, (2006).

# Installation using Anaconda

Information about Anaconda, including install instructions, can be found on the [Conda Reference](https://docs.conda.io/projects/conda/en/latest/) website.

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

## Using conda to compile IMPACT-T

`conda-forge` has all of the necessary compilers and dependencies to build IMPACT-T from source.

Create a build environment like so:

```bash
conda create -n impactt-build -c conda-forge compilers cmake openmpi
conda activate impactt-build
```

Then to build the non-parallel version:

```bash
cmake -S src/ -B build-single
cmake --build build-single -j 4
ls build-single/ImpactTexe
```

And to build the MPI-parallelized version with OpenMPI:

```bash
cmake -S src/ -B build-mpi -DUSE_MPI=ON
cmake --build build-mpi -j 4
ls build-mpi/ImpactTexe-mpi
```

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
If you like to install IMPACT-T into another directory than the default, pass to the `cmake -S . -B build` line the additional argument `-DCMAKE_INSTALL_PREFIX=/your/custom/install/path`.

### Multi Processor Code:

```shell script
# inside the IMPACT-T src/ directory:
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
cmake --build build --target install
```

## Windows

For Windows it will be necessary to use `NMake` to read and execute the generated makefiles.

`NMake` is command-line tool included with Visual Studio that builds projects based on commands that are contained in a description file.

More information on `NMake` can be found on the [NMAKE Reference](https://docs.microsoft.com/en-us/cpp/build/reference/nmake-reference?view=msvc-160) website.

### Single Processor Code:

```shell script
cmake -S . -B build -G "NMake Makefiles"
cmake --build build
cmake --build build --target install
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
module load cmake
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
# find the executable in build/bin/
```

### For KNL
```bash
module swap craype-haswell craype-mic-knl
module load openmpi # if using OpenMPI otherwise skip for MPICH
module load cmake
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
# find the executable in build/bin/
```

## Running at NERSC

There is one caveat with the conda-forge installed MPI version of the code.
Instead of running with `srun -n <cores> ` you must use `mpirun -n <cores>`.
