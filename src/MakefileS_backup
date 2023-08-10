#*****************************************************
#  General Makefile
#
#*****************************************************

#**************************************************************************
# Macros defining the Fortran, C/C++ compiler and linker.

CC = ifort
LINK = ifort
FFLAGS = -O3

#**************************************************************************
# List of .o files that EXENAME depends on.  Edit as appropriate for MP.

OBJS = \
	mpistub.o DataStruct/NumConst.o DataStruct/PhysConst.o DataStruct/Pgrid.o \
	DataStruct/Data.o \
	Func/Timer.o Func/Transpose.o Func/Fldmger.o Func/Ptclmger.o  \
	Func/FFT.o \
        Appl/BPM.o Appl/CCL.o Appl/CCDTL.o Appl/DTL.o Appl/SC.o Appl/Multipole.o \
	Appl/DriftTube.o Appl/Quadrupole.o Appl/ConstFoc.o Appl/SolRF.o \
	Appl/Sol.o Appl/Dipole.o Appl/EMfld.o Appl/EMfldCart.o Appl/EMfldCyl.o Appl/EMfldAna.o \
	Appl/BeamLineElem.o \
	Appl/CompDom.o Appl/BeamBunch.o Appl/Depositor.o Appl/Ranger.o Appl/Field.o Appl/Distribution.o \
	Contrl/Input.o Contrl/Output.o Contrl/AccSimulator.o Contrl/main.o

OBJS2 = \
	mpistub.o NumConst.o PhysConst.o Pgrid.o Data.o \
        Timer.o Transpose.o Fldmger.o Ptclmger.o FFT.o \
	BPM.o CCL.o CCDTL.o DTL.o SC.o DriftTube.o Quadrupole.o Multipole.o ConstFoc.o \
	SolRF.o Sol.o Dipole.o EMfld.o EMfldCart.o EMfldCyl.o EMfldAna.o \
	BeamLineElem.o CompDom.o BeamBunch.o Depositor.o Ranger.o Field.o Distribution.o \
	Input.o Output.o AccSimulator.o main.o	
#**************************************************************************
# Change this line if you don't like 'a.out'.

EXENAME = ImpactTexe

#************************************************************************
# disable predefined suffixes and define your own set of allowed suffixes
 .SUFFIXES:
 .SUFFIXES: .o .f .F .c .f90 .F90

#*************************************************************************
# inference rules (how to compile object files that have no explicit rules)
#  $* = name part of target
#  $@ = full target name
#  $< = dependent name

.f90.o:
	$(CC) -c $(FFLAGS) $<

#**************************************************************************
# Rules for building EXENAME from OBJS and OBJS from your source.

$(EXENAME):  $(OBJS) 
	$(LINK) -o $(EXENAME) $(OBJS2) 

#************************************************************************
# if you wish to compile a certain object with different flags
# or in some special way, then specify the target & dependency explicitly
# the following line say Timer.o is depended on Timer.f90
#Timer.o: Timer.f90
#	$(CC) -c -O3 Timer.f90

	cp  AccSimulator.o main.o Input.o Output.o Contrl
	cp  BPM.o CCL.o CCDTL.o DTL.o SC.o DriftTube.o Quadrupole.o \
	    ConstFoc.o BeamLineElem.o BeamBunch.o Field.o CompDom.o \
	    Distribution.o SolRF.o Sol.o Dipole.o Multipole.o EMfld.o \
	    EMfldCart.o EMfldCyl.o EMfldAna.o Depositor.o Ranger.o Appl
	cp  Timer.o Transpose.o Fldmger.o Ptclmger.o FFT.o Func	
	cp  NumConst.o PhysConst.o Data.o Pgrid.o DataStruct 
#***********************************************************************
clean:
	-rm *.o $(EXENAME) *.mod */*.o
