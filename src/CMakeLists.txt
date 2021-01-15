# set minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Uncomment for VERBOSE Makefiles. Useful for debugging.
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# project name and language
project(ImpactT LANGUAGES Fortran)

# BUILD_SHARED_LIBS is a global flag offered by CMake
# to toggle the behavior of add_library
# With OFF it will create a static library
set(BUILD_SHARED_LIBS OFF)

option(USE_MPI "Activate MPI build" OFF)


list(APPEND _sources
    DataStruct/NumConst.f90
    DataStruct/PhysConst.f90
    DataStruct/Pgrid.f90
	DataStruct/Data.f90
	Func/Timer.f90
	Func/Transpose.f90
	Func/Fldmger.f90
	Func/Ptclmger.f90
	Func/FFT.f90
    Appl/BPM.f90
    Appl/CCL.f90
    Appl/CCDTL.f90
    Appl/DTL.f90
    Appl/SC.f90
    Appl/Multipole.f90
	Appl/DriftTube.f90
	Appl/Quadrupole.f90
	Appl/ConstFoc.f90
	Appl/SolRF.f90
	Appl/Sol.f90
	Appl/Dipole.f90
	Appl/EMfld.f90
	Appl/EMfldCart.f90
	Appl/EMfldCyl.f90
	Appl/EMfldAna.f90
	Appl/BeamLineElem.f90
	Appl/CompDom.f90
	Appl/BeamBunch.f90
	Appl/Depositor.f90
	Appl/Ranger.f90
	Appl/Field.f90
	Appl/Distribution.f90
	Contrl/Input.f90
	Contrl/Output.f90
	Contrl/AccSimulator.f90
	Contrl/main.f90
)

set(OUTPUT_NAME "ImpactTexe")

if(USE_MPI)
    find_package(MPI REQUIRED)
    add_definitions(${MPI_Fortran_COMPILE_FLAGS} -DUSE_MPI)
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    set(OUTPUT_NAME "ImpactTexe-mpi")
else()
    configure_file(mpif_single.h mpif.h COPYONLY)
    list(APPEND _sources mpistub.f90)
endif()

add_executable(${OUTPUT_NAME} ${_sources})
target_link_libraries(${OUTPUT_NAME} ${MPI_Fortran_LIBRARIES})

# Enables preprocessing
target_compile_options(${OUTPUT_NAME} PRIVATE -cpp)

install(TARGETS ${OUTPUT_NAME} RUNTIME DESTINATION bin)
