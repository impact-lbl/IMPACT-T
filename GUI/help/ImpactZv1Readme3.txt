******************* Copyright Notice ********************************
IMPACT-Z  v 1.7 Copyright (c) 2012, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy).
All rights reserved.

If you have questions about your rights to use or distribute this
software, please contact Berkeley Lab's Technology Transfer Department
at  TTD@lbl.gov .

NOTICE.  This software was developed under partial funding from the U.S.
Department of Energy.  As such, the U.S. Government has been granted for
itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
worldwide license in the Software to reproduce, prepare derivative works,
and perform publicly and display publicly.  Beginning five (5) years after
the date permission to assert copyright is obtained from the U.S.
Department of Energy, and subject to any subsequent five (5) year
renewals, the U.S. Government is granted for itself and others acting
on its behalf a paid-up, nonexclusive, irrevocable, worldwide license
in the Software to reproduce, prepare derivative works, distribute
copies to the public, perform publicly and display publicly, and to
permit others to do so.

******************** End of Copyright Notice****************************
------------
AUTHORS
------------
Ji Qiang, Rober D. Ryne
Ms 71 J, ATAP
1 cyclotron Rd.
Lawrence Berkeley National Laboratory
Berkeley, CA 94720
Email: jqiang@lbl.gov
Tel: 510-495-2608
------------
CONTRIBUTORS
------------

------------
INTRODUCTION
------------
IMPACT-Z is a 3D parallel/serial Particle-In-Cell (PIC) code based on
multi-layer object-oriented design.
The present version of IMPACT-Z can treat intense beams propagating
through drifts, magnetic quadrupoles, magnetic soleniods, and rf cavities
using linear map integrator and drifts, magnetic dipoles, quadrupoles,
soleniods, and rf cavities using nonlinear Lorentz integrator.
It has a novel treatment of rf cavities, in which the gap transfer maps
are computed during the simulations by reading in Superfish rf fields.
The goal is to avoid time-consuming (and unnecessary) fine-scale
integration of millions of particles through the highly z-dependent
cavity fields. Instead, fine-scale integration is used to compute the
maps (which involve a small number of terms), and the maps are applied
to particles. If you are familiar with magnetic optics, then you will
recognize that this is analogous to the technique used to simulate
beam transport through magnets with fringe fields.

The version of IMPACT-Z currently
has a 3D space-charge model that assumes several types of
boundary conditions: 1) 3D open boundary conditions,
2) a solver that uses open boundary conditions transversely and
periodic boundary conditions longitudinally, and 3) solvers for
round and rectangular conducting pipes (where the longitudinal
boundary condition is either open or periodic).
It does not include the calculation of CSR wakefield through
the bending magnet since this version is primary used for
modeling the mult-charge state ion beams.
The error studies can include the field errors, misalignment
errors, and rotation errors.
---------------------------------------------------

To run the IMPACTz code on Windows PC, one needs to open a
DOS window using the "cmd" command. Then, one needs to type in
the ImpactTz executable name to run the code. All input files should be
under the same directory as the IMPACT executable code.

To run the IMPACTz code, the user needs to prepare some input files.
Note: when download those files from the web, make sure that there
is NO extra file extension that is automatically added to the file.

-----------
There are a few utility programs under:
http://portal.nersc.gov/project/m669/IMPACT/ilinac/utilities/
-----------
Tr2Impact.f90: takes the Twiss parameters of the initial
distribution and converts those into the initial distribution
parameters used in the IMPACT input file.
Engscan.py: does a single cavity phase scan.
RFcoef.f90: prepares Fourier expansion coefficients of RF or solenoid
field to be used in the simulation (with Lorentz nonlinear integrator.)
phaseOptZ.py: sets up the RF driven phase of the cavities with the
user specified design phases.
Theta: sets up the initial driven phase for rf cavities.
-----------
-----------
There are some sample examples under:
http://portal.nersc.gov/project/m669/IMPACT/ilinac/examples/

The IMPACTz has been used to model the Spallation Neutron
Source linac (room temperature and superconducting versions), the
Accelerator Production of Tritium linac, the quadrupole channel
for the LEDA beam halo experiment, the CERN
Superconducting Proton Linac, the JPAC linac, the RIA linac driver, etc.
-----------

IMPACT-Z INPUT FILES
------------------
The main IMPACT input file is called test.in.
Here is a sample with some of my comments in ()
WARNING: The / seen at the end of each line below is important:
the read statement in IMPACT that reads all the beamline element
info wants to read in a lot of numbers.
In a fortran free-format read statement, you can make it ignore
the rest by using / at the end of the line.
So, if you neglect this, the code will crash.
Also, any line starting with ! will be treated as a comment line.
!=========================================================================
                        INPUT FILE test.in
!=========================================================================
1 1                     (processor layout: NOTE WELL: the product of
                         these numbers must equal the number of
                         processors that you run on!)

6 10000000 1 0 1        (6D problem with 10M particles, "1" for linear map
                         integrator ("2" for nonlinear Lorentz integrator),
                         "0" for no error studies ("1" for error studies),
                         "1" for standard output ("2" for 90%, 95%, 99%
                         emittance, radius output).)
64 64 64 1 0.014 0.014 0.10  (64x64x64 space-charge grid, "1" for 3D open
                             all mesh numbers have to a power of 2, i.e. 2^n;
                             "2" for transverse open, longitudinal periodic,
                             (the 3rd mesh number has to 2^n+1)
                             "3" for transverse finite, longitudinal open
                             round pipe (in this case, the 2nd mesh number has
                             to be 2^n+1), "4" for transverse finite,
                             longitudinal periodic round pipe (in this case,
                             the 2nd and the 3rd mesh number has to be
                             2^n+1), "5" for
                             transverse finite, longitudinal open rectangular
                             pipe (in this case, the 1st and the 2nd mesh
                             number have to be 2^n+1),"6" for transverse finite,
                             longitudinal periodic rectangular pipe (in this
                             case, all mesh number have to be 2^n+1)),
                             x pipe width
                             "0.014" m, y pipe width "0.014" m, period length
                             "0.10" m.
3 0 0 4                 (Input distribution is "type 3",
                        Here is what the distributions mean:
                        1=rectangle in 6D phase space (for single charge)
                        2=6D Gaussian (single charge)
                        3=6D Waterbag  (single charge)
                        4=Semi-Gaussian [uniform position,Gaussian momentum,sc]
                        5=KV transverse, uniform longitudinal (sc).
                        16=6D Waterbag for multiple charge state (mc).
                        17=6D Gaussian for multiple charge state (mc).
                        19=Read coords/momenta from external file "partcl.data"
                        [this is
                        undergoing tests in Testing directory. More later.])
                        "0" means no restart ("1" means to restart from
                        some point after stop), "0" means no sub-cycle
                        ("1" means sub-cycle).
                        "4" for # of charge states = 4.
16000 16000 16000 16000 (# of particle list for each charge state)
0.0 0.0 0.0 0.0         (current for each charge state in unit of A)
3.97e-10 4.01e-10 4.06e-10 4.10e-10 (q_i/m_i for each charge state. Here, we normalized each charge
                                     by the mass of reference particle so that the reference particle
                                     has 1 AMU mass, but less charge. For example, for a proton, it is
                                     1./938.23e6.)
2.60529297811E-2, 3.6536731692E-4,  -0.95073401 1.0 1.0 0.0 0.0
1.69213266215E-2, 2.3573085744E-4,  -0.67283857 1.0 1.0 0.0 0.0
6.14232727818E-2, 1.7583221689E-4,  0.148895386 1.0 1.0 0.0 0.0
(Explanation of preceding 3 lines: They are
 sigmax    lambdax      mux     mismatchx  mismatchpx offsetX  offsetPx
 sigmay    lambday      muy     mismatchy  mismatchpy offsetY  offsetPy
 sigmaz    lambdaz      muz     mismatchz  mismatchE offsetPhase offsetEnergy
         The distribution in each dimension is a quadratic form.
         For example, in x, it is a function of
         x^2/sigmax^2 + 2*x*px*mux/(sigmax*lambdax)+px^2/lambdax^2
         Some day the code may be changed to work with Courant-Snyder
         parameters. All those quantities are dimensionless IMPACT internal units.)

0.070 120.e6 939.294e6 1.0 352.2e6 0.0 (current averaged over rf period=70 mA,
                                    initial kinetic energy=120MeV
                                    particle mass=939.294 MeV/c^2
                                    1.0=charge
                                    scaling frequency=352.2 MHz, 0.0 is
                                    initial phase of the reference particle)
!-----------------------------------------------------------------------------------
!lattice beam line element description starting below:
!-----------------------------------------------------------------------------------
0.770000  4  20    0    0.014  /    (drift, length=0.77m,
                   ^                 4 "steps" where each step
                 drift               through the beamline element
                                     consists of a half-step +
                                     a space-charge kick + a half-step.
                                     Each half-step involves
                                     computing a map for that
                                     half-element, computed by
                                     numerical integration with
                                     20 "map steps", pipe radus is 0.014 m.)
!
0.30      4  20    1 -5.67  0. 0.014  0. 0. 0. 0. 0. /(magnetic quad, length=0.3m
                   ^                        4 "steps," 20 "map steps"
                 quad                       gradient=-5.67 Tesla/m, input
                                            gradient file ID, radius=0.014m,
                                            x misalignment error=0.0m, y
                                            misalignment error=0.0m, rotation
                                            error x, y, z = 0.0, 0., 0. rad.)
!
0.30      4  20    2  9.8696 9.8696 9.8696 0.014 / (3D constant focusing,
                                            length=0.3m, 4 "steps",
                   ^                        " 20 "map steps",kx0^2=9.8696,
                 constant focusing          ky0^2=9.8696,kz0^2=9.8696,
                                            radius=0.014m. Note, it does not work
                                            for Lorentz integrator option.)
!
0.30      4  20    3  5.67  0. 0.014  0. 0. 0. 0. 0. 0. /(magnetic solenoid,length=0.3m
                   ^                        4 "steps," 20 "map steps"
                 solenoid                   Bz0=5.67 Tesla, input
                                            field file ID 0., radius=0.014m,
                                            x misalignment error=0.0m, y
                                            misalignment error=0.0m, rotation
                                            error x, y, z=0.0, 0., 0. rad.
                                            Note: the map integrator only
                                            works for v1.7 and later.
                                            For the Lorentz integrator, the
                                            length includes two linear fringe regions
                                            and a flat top region. Here, the
                                            length of the fringe region is
                                            defined by the 2*radius.
                                            The total length = effective length +
                                            2*radius.)
!
1.48524   10 20   4   0.1  0.0  150.  0.014  0.0 0.0 0.0 0. 0. /
                  ^                  (magnetic dipole, length=1.48524m,
                dipole               10 "steps", 20 "map steps"
                                     bending angle 0.1rad, k1=0.0, input
                                     switch 150., radius=0.014m,0.0=entr.
                                     pole face angle(rad),
                                     0.0=exit pole face angle, 0.0=curvature of
                                     entr. face, 0.0=curv. Of exit face,
                                     0.0=integrated fringe field.
!
1.48524   10 20   101   1.0 700.0e6 30. 1.0 0.014 0.01 1.0 0.01 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 /
                  ^              (RF cavity, length=1.48524m,
                 DTL             10 "steps", 20 "map steps"
                                 field scaling=1.0, RF frequency=700.0e6,
                                 driven phase=30.0 degree, input field
                                 ID=1.0, radius=0.014m, quad 1 length=0.01m,
                                 quad 1 gradient=1.0T/m, quad 2 length=0.01m,
                                 quad 2 gradient, x misalignment
                                 error=0.0m, y misalignment error=0.0m,
                                 rotation error x, y, z=0.0, 0., 0. rad for
                                 quad, x displacement, y displacement,
                                 rotation error x, y, z for RF field.
!
1.48524   10 20   102   1.0 700.0e6 30. 1.0 0.014 0. 0. 0. 0. 0. /
                  ^              (RF cavity, length=1.48524m,
                 CCDTL             10 "steps", 20 "map steps"
                                 field scaling=1.0, RF frequency=700.0e6,
                                 driven phase=30.0 degree, input field
                                 ID=1.0, radius=0.014m, x misalignment
                                 error=0.0m, y misalignment error=0.0m,
                                 rotation error x, y, z=0.0, 0., 0. rad)
!
1.48524   10 20   103   1.0 700.0e6 30. 1.0 0.014 0. 0. 0. 0. 0. /
                  ^              (RF cavity, length=1.48524m,
                 CCL             10 "steps", 20 "map steps"
                                 field scaling=1.0, RF frequency=700.0e6,
                                 driven phase=30.0 degree, input field
                                 ID=1.0, radius=0.014m, x misalignment
                                 error=0.0m, y misalignment error=0.0m,
                                 rotation error x, y, z=0.0, 0., 0. rad)
!
1.48524   10 20   104   1.0 700.0e6 30. 1.0 0.014 0. 0. 0. 0. 0. /
                  ^              (RF cavity, length=1.48524m,
                  SC             10 "steps", 20 "map steps"
                                 field scaling=1.0, RF frequency=700.0e6,
                                 driven phase=30.0 degree, input field
                                 ID=1.0, radius=0.014m, x misalignment
                                 error=0.0m, y misalignment error=0.0m,
                                 rotation error x, y, z=0.0, 0., 0. rad)
!
1.48524   10 20   105   1.0 700.0e6 30. 1.0 0.014 0. 0. 0. 0. 0. 1.0 /
                  ^              (Solenoid with RF cavity, length=1.48524m,
                  SolRF           10 "steps", 20 "map steps"
                                 field scaling=1.0, RF frequency=700.0e6,
                                 driven phase=30.0 degree, input field
                                 ID=1.0, radius=0.014m, x misalignment
                                 error=0.0m, y misalignment error=0.0m,
                                 rotation error x, y, z=0.0, 0., 0. rad,
                                 Bz0=1.0 Tesla)
!
1.48524   10 20   110   1.0 700.0e6 30. 1.0 0.014 0.014 0. 0. 0. 0. 0. 1.0 1.0 /
                  ^              (user defined RF cavity, length=1.48524m,
                  EMfld          10 "steps", 20 "map steps"
                                 field scaling=1.0, RF frequency=700.0e6,
                                 driven phase=30.0 degree, input field
                                 ID=1.0, radius=0.014m, x misalignment
                                 error=0.0m, y misalignment error=0.0m,
                                 rotation error x, y, z=0.0, 0., 0. rad),
                                 1.0 (using discrete data only, 2.0 using both
                                 discrete data and analytical function, other
                                 using analytical function only), 2.0 (field
                                 in Cartesian coordinate, 1.0 in Cylindrical
                                 coordinate)
!
0. 0 0 -1 /   (shift the centroid of beam to the axis)
!
0. 0 N -2 0.0 1 /   (write the full particle set on fort.N
        ^     ^ NOTE WELL: N must not equal 5 or 6 (Fortran code)
      write     NOTE WELL: This prints the data set with 1 as sample
                frequency. NOTE: those particle data are dimensionless in
                IMPACT internal unit. The normalization constants are
                described in the following Partcl.data.
!
0. 0 0 -3 0.014 0.02 0.02 0.02 0.02 0.02 0.02 /
                       (write the accumulated density along R, X, Y into file
                        Xprof.data, Yprof.data, RadDens.data, radius=0.014m,
                        xmax=0.02m,pxmax=0.02 mc,ymax=0.02m,pymax=0.02 mc,
                        zmax=0.02 rad,pzmax=0.02 mc^2, if no frame range
                        is specified, the program will use the maximum
                        amplitude at given location as the frame.)
!
0. 0 0 -4 0.014 0.02 0.02 0.02 0.02 0.02 0.02 /
                       (write the density along R, X, Y into file
                       Xprof2.data, Yprof2.data, RadDens2.data, radius=0.014m,
                       xmax=0.02m,pxmax=0.02 mc,ymax=0.02m,pymax=0.02 mc,
                       zmax=0.02 rad,pzmax=0.02 mc^2, if no frame range
                       is specified, the program will use the maximum
                       amplitude at given location as the frame.)
!
0. 0 0 -5 0.014 0.02 0.02 0.02 0.02 0.02 0.02 /
                       (write the 2D projections of 6D distribution in
                       XY.data, XPx.data, XZ.data, YPy.data, YZ.data,
                       ZPz.data. radius=0.014m,
                       xmax=0.02,pxmax=0.02 mc,ymax=0.02,pymax=0.02 mc,
                       zmax=0.02 rad,pzmax=0.02 mc^2, if no frame range
                       is specified, the program will use the maximum
                       amplitude at given location as the frame.)
!
0. 0 0 -6 0.014 0.02 0.02 0.02 0.02 2 0.02 /
                       (write the 3D density into the file fort.8.
                       radius=0.014m,
                       xmax=0.02m,pxmax=0.02 mc,ymax=0.02m,pymax=0.02 mc,
                       zmax=2 degree,pzmax=0.02 mc^2, if no frame range
                       is specified, the program will use the maximum
                       amplitude at given location as the frame.)
!
0. 0 0 -7 /
                       (write the 6D phase space information and local
                       computation domain into files ph0000sN and gm0000sN.
                       This function is used for restart purpose.
!
0. 0 0 -21 0.014 0.02 0.02 0.02 0.02 0.02 0.02 /
                      (shift the beam centroid in 6D phase space.
                      radius=0.014m (not used),
                      xshift=0.02m,pxshift=0.02rad,yshift=0.02m,pyshift=0.02rad,
                      zshift=0.02deg,pzshift=0.02MeV.)
!
0.  1 1 -99 /
         ^      (halt execution at this point in the input file.
        halt     this is useful if you have a big file want to run
                part-way through it without deleting a lot of lines.
!
=====================================================================
The other data files are called rfdataN, e.g. rfdata1, rfdata52,...
The first column of this file is "z," the second column is "Ez"
on-axis, and the third column is d/dz(Ez on-axis). This is what
is needed to compute the linear transfer map for the rf cavities.
These formats are used for the MAP integrator. For the NONLINEAR Lorentz
integrator, the rfdataN contains the Fourier expansion coefficients for the Ez
on the axis. These coefficients can be generated from the discrete on-axis
Ez data using the code Rfcoef. Rfcoef will take rfdata.in as input and output
Fourier expansion coefficients into file rfdataX. Here, X needs to be replaced by N for different type of rf cavities.

The initial distribution type 19 for read in distribution that is stored in the file "partcl.data".
The first line of the partcl.data is:
# of particles, 0.0 0.0
The rest of lines are particle information.
There are 9 columns for each particle:
X, Px, y, Pz, t, Pt, charge/mass, charge for each particle, id.
Here:
------------------
1st col: x normalized by c/omega
2nd col: Px normalized by mc, here this column has to be divided by gamma beta
         to convert to unit radian
3rd col: y normalized by c/omega
4th col: Py normalized by mc, here this column has to be divided by gamma beta
         to convert to unit radian
5th col: phase (radian)
6th col: kinetic energy deviation (E0 -E) normalized by mc^2
7th col: q/m, for an electron, q = -1, m = 0.511001e6
8th col: charge per macro-particle
9th col: particle id
--------------------

IMPACT-Z OUTPUT FILES
-------------------
File fort.18: reference particle information
----------------------------------------------------------------
1st col: distance (m)
2nd col: absolute phase (radian)
3rd col: gamma
4th col: kinetic energy (MeV)
5th col: beta
6th col: Rmax (m) R is measured from the axis of pipe
----------------------------------------------------------------
File fort.24 (for X) fort.25 (for Y), fort.26 (for Z): RMS size information
----------------------------------------------------------------
1st col: z distance (m)
2nd col: centroid location (m)
3rd col: RMS size (m)
4th col: Centroid momentum (rad for fort.24 and fort.25, MeV for fort.26)
5th col: RMS momentum (rad for fort.24 and fort.25, MeV for fort.26)
6th col: Twiss parameter, alpha
7th col: normalized RMS emittance (m-rad for transverse and degree-MeV for for
t.26)
----------------------------------------------------------------
File fort.27: maximum amplitude information
----------------------------------------------------------------
1st col: z distance (m)
2nd col: Max. X (m)
3rd col: Max. Px (rad)
4th col: Max. Y (m)
5th col: Max. Py (rad)
6th col: Max. Phase (degree)
7th col: Max. Energy deviation (MeV)
----------------------------------------------------------------
File fort.28: load balance and loss diagnostic
----------------------------------------------------------------
1st col: z distance (m)
2nd col: min # of particles on a PE
3rd col: max # of particles on a PE
4th col: total # of particles in the bunch
----------------------------------------------------------------
File fort.29: cubic root of 3rd moments of the beam distribution
----------------------------------------------------------------
1st col: z distance (m)
2nd col: X (m)
3rd col: Px (rad)
4th col: Y (m)
5th col: Py  (rad)
6th col: phase (degree)
7th col: Energy deviation (MeV)
----------------------------------------------------------------
File fort.30: square root, square root of 4th moments of the beam distribution
----------------------------------------------------------------
1st col: z distance (m)
2nd col: X (m)
3rd col: Px (rad)
4th col: Y (m)
5th col: Py (rad)
6th col: phase (degree)
7th col: Energy deviation (MeV)
----------------------------------------------------------------
File fort.32: number of particles for each charge state
----------------------------------------------------------------
1st col: z distance (m)
2nd col: number of particles for each charge state
----------------------------------------------------------------
Lastly, if you use the "-2" type code (described in the
explanation of test.in) to print the data at various locations,
the resulting data will be in columns of the form
x,px,y,py,t,pt, q/m, charge/per particle, id
------------------------------------
1st col: x normalized by c/omega
2nd col: Px normalized by mc, here this column has to be divided by gamma beta
         to convert to unit radian
3rd col: y normalized by c/omega
4th col: Py normalized by mc, here this column has to be divided by gamma beta
         to convert to unit radian
5th col: phase (radian)
6th col: kinetic energy deviation (E0 -E) normalized by mc^2
7th col: q/m, for an electron, q = -1, m = 0.511001e6
8th col: charge per macro-particle
9th col: particle id
!

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!Theta code
!
! This code is used to find the driven phase theta0 (in cos(wt+theta0))
! and field amplitude (E0) to achieve the designed energy gain
! and 0 current phase advance in one lattice period. This code also
! sets up the lattice structure for the macroparticle beam dynamics
! simulation.
! Here, the synchronous phase is defined at either the center of
! RF gap (isynflag!=0) or at the electrical center (isynflag=0).
! The code uses the longitudinal coordinate, z, as independent variable.
!
!--------------------------------------------------------------------------
1)Input file: sim.dat
-----------------------------------------------------------
402.5e6 88.e6 939.294e6 0.0 -1.0 2 0 2 7 88.78 0.001 0.03 75.0 0.001 0.02
                                        scaling frequency, input kinetic
                                        energi (eV), particle mass (eV),
                                        initial phase a(1) (rad), charge
                                        number, idata (1 for discrete RF
                                        data, 2 for analytical function,
                                        0(1) without(with) error study ),
                                        2(1) given energy gain (given syn.
                                        phase), 7 for number of elements per
                                        period, 88.78 (longitudinal 0 current
                                        phase advance), 0.001 (convergence
                                        tolerance), 0.03 (E field amplitude
                                        increase step), 75.0 (transverse 0
                                        current phase advance), 0.001
                                        (convergence tolerance), 0.02 (B field
                                         amplitude increase step).
0.1364384E+00   40  2   0   10.000000E+00 /  length (m), # of steps, # of
                        ^                   substeps, beam line element type,
                       drift                radius.
                                            focusing strength (T/m), RF field
                                            scaling constant, RF frequency
                                            (Hz), energy gain (MeV),
                                            synchronous phase (degree),
                                            isynflag (location of synchronous,
                                            0 at the gap center, 1 at the
                                            electrical center
0.8000000E-01   40  2   1   0.305570E+02 0.0 10.0 / length (m), # of steps,
                        ^                   # of substeps, beam line element
                       quad                 type,focusing strength (T/m),field
                                            file ID, radius=10.0m.
0.8000000E-01   40  2   2   0.3E+02 0.3E+2 0.3E+2 10.0 /length, # of steps,
                        ^                   # of substeps, beam line element
                   Uniform focusing         type,kx0^2,ky0^2,kz0^2,radius
0.8000000E-01   40  2   3   0.305570E+02 0.0 10.0 / length (m), # of steps,
                        ^                   # of substeps, beam line element
                   solenoid                 type,focusing field Bz (T),field
                                            file ID, radius=10.0m.
0.8000000E-01   40  2   4   0.305570E+02 0.0 0.0 10.0 /length, # of steps,
                        ^                    # of substeps, beam line element
                      dipole                 type, Bx(T), By(T), field
                                            file ID, radius=10.0m.
0.6060880E+00  400  2   103  0.1000E+07 805.e6 -30.0 2.0 1.0 10.0 1.0 /
                        ^                  length (m), # of steps, # of
                       CCL                 sub steps, beam line element type
                                           (CCDTL, CCL, SC),
                                           scaling constant, RF freq.,
                                           synchronous phase (degree), energy
                                           gain (MeV), field file ID, radius,
                                           isynflag (location of synchronous,
                                           >=1 at the gap center, 0 at the
                                           electrical center).
0.6060880E+00  400  2   101  0.1000E+07 805.e6 -30.0 2.0 1.0 10.0 1.0 0.01, 1.0 0.01 1.0 /
                        ^                  length (m), # of steps, # of
                        DTL                sub steps, beam line element type,
                                           scaling constant, RF freq.,
                                           synchronous phase (degree), energy
                                           gain (MeV), field file ID, radius,
                                           isynflag (location of synchronous,
                                           >=1 at the gap center, 0 at the
                                           electrical center), quad 1 length,
                                           quad 1 strength, quad 2 length,
                                           quad 2 strength.
0.6060880E+00  400  2   105  0.1000E+07 805.e6 -30.0 2.0 1.0 10.0 1.0 2.0 /
                        ^                  length (m), # of steps, # of
                      Sol+RF               sub steps, beam line element type,
                                           scaling constant, RF freq.,
                                           synchronous phase (degree), energy
                                           gain (MeV), field file ID, radius,
                                           isynflag (location of synchronous,
                                           >=1 at the gap center, 0 at the
                                           electrical center), Bz (T).
! -----------------------------------
2) Input file: errorlimit.dat
! -----------------------------------
1.0 1.0 1.0 1.0 1.0 1.0 1.0   Quad x, y, z displacement errors, pitch, yaw,
                              roll errors, and gradient error
1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0  RF cavity x, y, z displacement errors, pitch,                                 yaw, roll errors, field, and phase errors
! -----------------------------------
3) Input file: rfdatax: x is the field file ID.
Warning: Ez' has to be accurate. Discrete data for Ez and Ez' is recomended.
! -----------------------------------
z, ez, ez'   (for discrete Ez(0,z) data)
---
2.0          (for Fourier coeficients of Ez(0,z) data, first line is F0,
0.0          second line is F1 cos term,
0.0          third line is F1 sin term,
1.0          4th line is F2 cos term,
0.0          5th line is F2 sin term). The line number has to be odd.
! -----------------------------------
Output files:
!----------------------------------
fort.24: output lattice layout for IMPACT input file test.in
fort.18: output z,phase wt,gamma,energy,beta,wt+theta0.
fort.19: output z,g0,ez,ez'.
fort.68: output gap or cell #,z,phase at cent. of gap or elect. field.
fort.69: output gap or cell #,z,phase at cent. of gap or elect. field.
phase.adv: 0 current phase advance period.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!Rfcoef code
! This code is used to find Fourier expansion coefficients for the RF data, Ez on axis. The off-axis field will be reconstructed based these coefficients and axis symmetry assumption.
! Input: rfdata.in: Ez vs. z
! Output: rfdata.out: reconstruct the field from the Fourier coefficients.
!         rfdatax: is the one to be used for the IMPACT-Z
!         simulation.
