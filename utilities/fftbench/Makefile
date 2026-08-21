# Builds fftbench, which compares the two IMPACT-T FFT backends
# (bundled FFTPACK 5.1 vs FFTW3) for per-function agreement and speed.
#
# Usage:
#   make [FC=gfortran] [FFTW_HOME=/opt/homebrew]
#   ./fftbench             # default sizes 64 128 256 512
#   ./fftbench 64 90 128   # explicit transform sizes (even numbers)
#
# The main build compiles exactly one FFT implementation and neither is a
# standalone backend module: FFT.f90 keeps the FFTPACK primitives inside
# module FFTclass (mixed with the MPI/Transpose-dependent high-level
# routines), and FFT_fftw.f90 keeps the FFTW primitives in module
# FFTbackendclass next to its own copy of FFTclass. To exercise both
# backends from one binary, the two sets of primitives are lifted into
# their own dependency-free modules, renamed so they can link together.
# They are generated from the canonical sources below, so they never drift.

FC ?= gfortran
FFLAGS ?= -O3
FFTW_HOME ?= $(shell pkg-config --variable=prefix fftw3 2>/dev/null || echo /usr/local)
FFT_SRC = ../../src/Func

fftbench: fftpack_backend.f90 fftw_backend.f90 fftbench.f90
	$(FC) $(FFLAGS) -I$(FFTW_HOME)/include -o $@ $^ -L$(FFTW_HOME)/lib -lfftw3

# FFTPACK: four1/realft/sinft/cosft1 and the bundled FFTPACK 5.1 library run
# from `subroutine four1` to the end of module FFTclass and use nothing
# outside that block, so wrap that span in a standalone module.
fftpack_backend.f90: $(FFT_SRC)/FFT.f90
	awk 'BEGIN{e=0} \
	     /^subroutine four1\(/&&!e{e=1;print "      module fftpackbackendclass";print "      contains"} \
	     e&&/end module FFTclass/{print "      end module fftpackbackendclass";exit} \
	     e' $< > $@

# FFTW: take only module FFTbackendclass from FFT_fftw.f90 (skip its FFTclass,
# which pulls in Timer/Transpose/MPI) and rename it so both backends coexist.
fftw_backend.f90: $(FFT_SRC)/FFT_fftw.f90
	sed -n '/^      module FFTbackendclass/,/^      end module FFTbackendclass/p' $< \
	  | sed 's/FFTbackendclass/fftwbackendclass/g' > $@

clean:
	rm -f fftbench fftpack_backend.f90 fftw_backend.f90 *.mod *.o

.PHONY: clean
