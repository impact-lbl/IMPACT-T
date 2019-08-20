IMPACT-T Documentation  
====================

IMPACT-T: A 3D Parallel Particle Tracking Code in Time Domain
-------------------------------------------------------------------------------------

IMPACT-T is a fully three-dimensional program to track relativistic charged particles taking into account space charge forces, short-range longitudinal and transverse wakefields, coherent synchrotron radiation (CSR) wakefield in accelerators. IMPACT-T code can run on both massive parallel supercomputers and single processor computers such as Windows PC, Mac, and Linux system. It is one of the few codes used in the photoinjector community that has a parallel implementation, making it very useful for high statistics simulations of beam halos and beam diagnostics. It has a comprehensive set of beamline elements, and furthermore allows arbitrary overlap of their fields, which gives the IMPACT-T a capability to model both the standing wave structure and traveling wave structure. It includes mean-field space-charge solvers based on an integrated Green function to efficiently and accurately treat beams with large aspect ratio, and a shifted Green function to efficiently treat image charge effects of a cathode. It is also unique in its inclusion of energy binning in the space-charge calculation to model beams with large energy spread. It also has a direct N-body solver to calculate stochastic space-charge forces. IMPACT-T has a flexible data structure that allows particles to be stored in containers with common characteristics; for photoinjector simulations the containers represent multiple slices, but in other applications they could correspond, e.g., to particles of different species. Together, all these features make IMPACT-T a powerful and versatile tool for modeling beams in photoinjectors and other systems.

Here is the link to the home page of IMPACT-T: <https://amac.lbl.gov/~jiqiang/IMPACT-T/index.html>

Here is the link to the GitHub of IMPACT-T: <https://github.com/impact-lbl/IMPACT-T>

This is the license statement: 
---------------------------------------
> *** License Agreement ***
>> IMPACT-T, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved."
>>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
>> 1.		Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
>> 2.		Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
>> 3.		Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

>> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

>> You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a  non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.

This is the README file: 
---------------------------------
> *** Copyright Notice ***
>> "IMPACT-T" Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

>> If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

>> NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
> ****************************
> V2.0
>> Note: 
>> 1.		The current version of the code is for serial single processor computer with Fortran90 compiler. To run the code on a parall computer with MPI, the user has to comment out the line "use mpistub" in Contrl/Input.f90, DataStruct/Data.f90, DataStruct/Pgrid.f90, DataStruct/PhysConst.f90, and Func/Timer.f90. The user also has to remove the mpif.h file under the Appl, Control, DataStruct, and Func directories. The user also has to modify the Makefile to remove the mpistub.o inside the file and to use the appropriate parallel F90 compiler such as mpif90.
>> 2.		The phaseOpt.py is used to find the driven phase of a RF cavity with initial design phase. This code needs to be modified for each input ImpactT.in file in order to use it correctly.
>> 3.		The subroutines in FFT.f90: realft, four1, and sinft, can be replaced with functions from the Numerical Recipe or some equavilent 1D FFT functions.

