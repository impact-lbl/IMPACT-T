!----------------------------------------------------------------
!*** Copyright Notice ***
!
!"IMPACT-T" Copyright (c) 2016, The Regents of the University of California,
!through Lawrence Berkeley National Laboratory (subject to receipt of any
!required approvals from the U.S. Dept. of Energy).  All rights reserved.
!
!If you have questions about your rights to use or distribute this software,
!please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
!NOTICE.  This Software was developed under funding from the U.S. Department of
!Energy and the U.S. Government consequently retains certain rights. As such,
!the U.S. Government has been granted for itself and others acting on its behalf
!a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce,
!distribute copies to the public, prepare derivative works, and perform publicly
!and display publicly, and to permit other to do so.
!
!****************************
!
! AccSimulatorclass: Linear accelerator simulator class in CONTROL layer.
! Version: 2.1
! Author: Ji Qiang
! Description: This class defines functions to set up the initial beam 
!              particle distribution, field information, computational
!              domain, beam line element lattice and run the dynamics
!              simulation through the system.
! Comments:
!----------------------------------------------------------------
      module AccSimulatorclass
        use Pgrid2dclass
        use CompDomclass
        use FieldQuantclass
        use BeamLineElemclass
        use Ptclmgerclass
        use BeamBunchclass
        use Timerclass
        use Inputclass
        use Outputclass
        use Dataclass
        use PhysConstclass
        use NumConstclass
        use Distributionclass
        use Rangerclass
        use Depositorclass
        implicit none
        !# of phase dim., num. total and local particles, int. dist. 
        !and restart switch, error study switch, substep for space-charge
        !switch,# of time step
        integer :: Dim, Flagdist,Rstartflg,Flagerr,&
                            Flagsubstep,ntstep 
        integer, dimension(Nbunchmax) :: Np, Nplocal
        !# of num. total x, total and local y mesh pts., type of BC, 
        !# of beam elems, type of integrator.
        !FlagImage: switch flag for image space-charge force calculation: "1" for yes, 
        !otherwise for no. 
        integer :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal,Flagbc,&
                            Nblem,Flagmap,Flagdiag,FlagImage

        !# of processors in column and row direction.
        integer :: npcol, nprow

        !initial # of bunches/bins
        integer :: Nbunch

        !beam current, kin. energy, part. mass, charge, ref. freq., period length, 
        !time step size 
        double precision :: Bcurr,Bkenergy,Bmass,Bcharge,Bfreq,&
                                     Perdlen,dt,xrad,yrad

        !conts. in init. dist.
        integer, parameter :: Ndistparam = 21
        double precision, dimension(Ndistparam) :: distparam

        !2d logical processor array.
        type (Pgrid2d) :: grid2d

        !beam particle object and array.
        type (BeamBunch), dimension(Nbunchmax) :: Ebunch

        !beam charge density and field potential arrays.
        type (FieldQuant) :: Potential

        !geometry object.
        type (CompDom) :: Ageom

        !overlaped external field data array
        type (fielddata), dimension(Maxoverlap) :: fldmp

        !maximum e- emission time
        double precision :: temission
        !number of steps for emission
        integer :: Nemission

        !distance after that to turn off image space-charge
        double precision :: zimage

        !restart time and step
        double precision :: tend,dtlessend
        integer :: iend,ibchend,nfileout,ioutend,itszend,isteerend,isloutend

        !beam line element array.
        type (BPM),target,dimension(Nbpmmax) :: beamln0
        type (DriftTube),target,dimension(Ndriftmax) :: beamln1
        type (Quadrupole),target,dimension(Nquadmax) :: beamln2
        type (DTL),target,dimension(Ndtlmax) :: beamln3
        type (CCDTL),target,dimension(Nccdtlmax) :: beamln4
        type (CCL),target,dimension(Ncclmax) :: beamln5
        type (SC),target,dimension(Nscmax) :: beamln6
        type (ConstFoc),target,dimension(Ncfmax) :: beamln7
        type (SolRF),target,dimension(Nslrfmax) :: beamln8
        type (Sol),target,dimension(Nslmax) :: beamln9
        type (Dipole),target,dimension(Ndipolemax) :: beamln10
        type (EMfld),target,dimension(Ncclmax) :: beamln11
        type (EMfldCart),target,dimension(Ncclmax) :: beamln12
        type (EMfldCyl),target,dimension(Ncclmax) :: beamln13
        type (EMfldAna),target,dimension(Ncclmax) :: beamln14
        type (Multipole),target,dimension(Nquadmax) :: beamln15
        type (BeamLineElem),dimension(Nblemtmax)::Blnelem
        !longitudinal position of each element (min and max).
        double precision, dimension(2,Nblemtmax)::zBlnelem
        !beam line element period.
        interface construct_AccSimulator
          module procedure init_AccSimulator
        end interface
      contains
        !set up objects and parameters.
        subroutine init_AccSimulator(time)
        implicit none
        include 'mpif.h'
        integer :: i,test1,test2,j
        integer :: myid,myidx,myidy,ierr,inb,jstp
        integer, allocatable, dimension(:) :: bnseg,bmpstp,bitype
        double precision, allocatable, dimension(:) :: blength,val1,val0,&
        val2, val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24
        double precision :: time
        double precision :: t0,zz
        double precision :: z,phsini
        double precision, dimension(2) :: tmpdr 
        double precision, dimension(5) :: tmpcf 
        double precision, dimension(9) :: tmpbpm 
        double precision, dimension(9) :: tmpquad
        double precision, dimension(10) :: tmpdipole 
        double precision, dimension(11) :: tmprf
        double precision, dimension(12) :: tmpslrf
        double precision, dimension(13) :: tmp13
        double precision, dimension(25) :: tmpdtl
        integer :: iqr,idr,ibpm,iccl,iccdtl,idtl,isc,icf,islrf,isl,idipole,&
          iemfld,iemfldcart,iemfldcyl,iemfldana,ib,imult
        integer, allocatable, dimension(:) :: seedarray
        real*8 rancheck
        integer :: seedsize

        !start up MPI.
        call init_Input(time)

        ! initialize Timer.
        call construct_Timer(0.0d0)

        call starttime_Timer(t0)

        !Flagmap = 0

!-------------------------------------------------------------------
! get all global input parameters.
        call in_Input(Dim,Np(1),Nx,Ny,Nz,Flagbc,Flagdist,Rstartflg,&
              Flagmap,distparam,Ndistparam,Bcurr,Bkenergy,Bmass,Bcharge,&
        Bfreq,xrad,yrad,Perdlen,Nblem,npcol,nprow,Flagerr,Flagdiag,&
        Flagsubstep,phsini,dt,ntstep,Nbunch,FlagImage,Nemission,&
        temission,zimage)
 
!        print*,"Np: ",Np,dt,ntstep,Nx,Ny,Nz
!        print*,"Bcurr: ",Bcurr,Bkenergy,Bmass,Bcharge,Bfreq
!-------------------------------------------------------------------
! construct 2D logical processor Cartesian coordinate
        call construct_Pgrid2d(grid2d,MPI_COMM_WORLD,nprow,npcol)
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        if(myid.eq.0) then
          !print*,"Start simulation:"
          print*,"!-----------------------------------------------------------"
          print*,"! IMPACT-T Parallel Beam Dynamics Tracking Code: 2.1 beta version"
          print*,"! Copyright of The Regents of the University of California"
          print*,"!-----------------------------------------------------------"
        endif

        !construct Constants class.
        call constructT_PhysConst(dt,Bfreq)

!-------------------------------------------------------------------
! construct computational domain CompDom class and get local geometry 
! information on each processor.
        !if(Rstartflg.eq.1) then
        !  call ingeom_Output(1500,z,inb,jstp,nprow,npcol,Ageom,Nx,Ny,Nz,&
        !                    myidx,myidy)
        !  if(myid.eq.0) print*,"rstart at: ",z,inb,jstp
        !  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !else
          !xrad = 0.1363243029*0.2
          call construct_CompDom(Ageom,distparam,Ndistparam,Flagdist,&
               Nx,Ny,Nz,grid2d,nprow,npcol,Flagbc,xrad,yrad,Perdlen)
        !endif

!-------------------------------------------------------------------
! initialize Data class.
        do i = 1, Maxoverlap
          call initt_Data(fldmp(i))
        enddo

!-------------------------------------------------------------------
! construct BeamBunch class.
        call construct_BeamBunch(Ebunch(1),Bcurr,Bkenergy,Bmass,Bcharge,&
                            Np(1),phsini)

        !initial value for the output file number
        nfileout = 40
!-------------------------------------------------------------------
! construct beam line elements.
        allocate(blength(Nblem),bnseg(Nblem),bmpstp(Nblem))
        allocate(bitype(Nblem))
        allocate(val0(Nblem))
        allocate(val1(Nblem),val2(Nblem),val3(Nblem),val4(Nblem))
        allocate(val5(Nblem),val6(Nblem),val7(Nblem),val8(Nblem))
        allocate(val9(Nblem),val10(Nblem),val11(Nblem),val12(Nblem))
        allocate(val13(Nblem),val14(Nblem),val15(Nblem),val16(Nblem))
        allocate(val17(Nblem),val18(Nblem),val19(Nblem),val20(Nblem))
        allocate(val21(Nblem),val22(Nblem),val23(Nblem),val24(Nblem))

        call in_Input(Nblem,blength,bnseg,bmpstp,bitype,val0,val1,&
        val2,val3,&
        val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24)

        iccl = 0
        iccdtl = 0
        idtl = 0
        isc = 0
        idr = 0
        iqr = 0
        ibpm = 0
        icf = 0
        islrf = 0
        isl = 0
        idipole = 0
        iemfld = 0
        iemfldcart = 0
        iemfldcyl = 0
        iemfldana = 0
        imult = 0
        zz = 0.0
        !If we allow the starting edge "zz" of one beam line element
        !to be inside the preceding beam, this allows the beam line
        !elment overlaping. 
        do i = 1, Nblem
          zBlnelem(1,i) = val0(i)
          if(bitype(i).lt.0) then
            ibpm = ibpm + 1
            call construct_BPM(beamln0(ibpm),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpbpm(1) = val0(i)
            tmpbpm(2) = val1(i)
            tmpbpm(3) = val2(i)
            tmpbpm(4) = val3(i)
            tmpbpm(5) = val4(i)
            tmpbpm(6) = val5(i)
            tmpbpm(7) = val6(i)
            tmpbpm(8) = val7(i)
            tmpbpm(9) = val8(i)
            call setparam_BPM(beamln0(ibpm),tmpbpm)
            Blnelem(i) = assign_BeamLineElem(beamln0(ibpm))
            !reset the output file number from BPM type "-3".
            if(bitype(i).eq.(-3)) then
              nfileout = bmpstp(i)
            endif
          else if(bitype(i).eq.0) then
            idr = idr + 1
            call construct_DriftTube(beamln1(idr),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpdr(1) = val0(i)
            tmpdr(2) = val1(i)
            call setparam_DriftTube(beamln1(idr),tmpdr)
            Blnelem(i) = assign_BeamLineElem(beamln1(idr))
          else if(bitype(i).eq.1) then
            iqr = iqr + 1
            call construct_Quadrupole(beamln2(iqr),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmprf(1) = val0(i)
            tmprf(2) = val1(i)
            tmprf(3) = val2(i)
            tmprf(4) = val3(i)
            tmprf(5) = val4(i)
            tmprf(6) = val5(i)
            tmprf(7) = val6(i)
            tmprf(8) = val7(i)
            tmprf(9) = val8(i)
            tmprf(10) = val9(i)
            tmprf(11) = val10(i)
            call setparam_Quadrupole(beamln2(iqr),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln2(iqr))
          else if(bitype(i).eq.2) then
            icf = icf + 1
            call construct_ConstFoc(beamln7(icf),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpcf(1) = val0(i)
            tmpcf(2) = val1(i)
            tmpcf(3) = val2(i)
            tmpcf(4) = val3(i)
            tmpcf(5) = val4(i)
            call setparam_ConstFoc(beamln7(icf),tmpcf)
            Blnelem(i) = assign_BeamLineElem(beamln7(icf))
          else if(bitype(i).eq.3) then
            isl = isl + 1
            call construct_Sol(beamln9(isl),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpquad(1) = val0(i)
            tmpquad(2) = val1(i)
            tmpquad(3) = val2(i)
            tmpquad(4) = val3(i)
            tmpquad(5) = val4(i)
            tmpquad(6) = val5(i)
            tmpquad(7) = val6(i)
            tmpquad(8) = val7(i)
            tmpquad(9) = val8(i)
            call setparam_Sol(beamln9(isl),tmpquad)
            Blnelem(i) = assign_BeamLineElem(beamln9(isl))
          else if(bitype(i).eq.4) then
            idipole = idipole + 1
            call construct_Dipole(beamln10(idipole),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpdipole(1) = val0(i)
            tmpdipole(2) = val1(i)
            tmpdipole(3) = val2(i)
            tmpdipole(4) = val3(i)
            tmpdipole(5) = val4(i)
            tmpdipole(6) = val5(i)
            tmpdipole(7) = val6(i)
            tmpdipole(8) = val7(i)
            tmpdipole(9) = val8(i)
            tmpdipole(10) = val9(i)
            call setparam_Dipole(beamln10(idipole),tmpdipole)
            Blnelem(i) = assign_BeamLineElem(beamln10(idipole))
          else if(bitype(i).eq.5) then
            imult = imult + 1
            call construct_Multipole(beamln15(imult),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpdipole(1) = val0(i)
            tmpdipole(2) = val1(i)
            tmpdipole(3) = val2(i)
            tmpdipole(4) = val3(i)
            tmpdipole(5) = val4(i)
            tmpdipole(6) = val5(i)
            tmpdipole(7) = val6(i)
            tmpdipole(8) = val7(i)
            tmpdipole(9) = val8(i)
            tmpdipole(10) = val9(i)
            call setparam_Multipole(beamln15(imult),tmpdipole)
            Blnelem(i) = assign_BeamLineElem(beamln15(imult))
          else if(bitype(i).eq.101) then
            idtl = idtl + 1
            call construct_DTL(beamln3(idtl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpdtl(1) = val0(i)
            tmpdtl(2) = val1(i) 
            tmpdtl(3) = val2(i) 
            tmpdtl(4) = val3(i) 
            tmpdtl(5) = val4(i) 
            tmpdtl(6) = val5(i) 
            tmpdtl(7) = val6(i) 
            tmpdtl(8) = val7(i) 
            tmpdtl(9) = val8(i) 
            tmpdtl(10) = val9(i) 
            tmpdtl(11) = val10(i) 
            tmpdtl(12) = val11(i) 
            tmpdtl(13) = val12(i) 
            tmpdtl(14) = val13(i) 
            tmpdtl(15) = val14(i) 
            tmpdtl(16) = val15(i) 
            tmpdtl(17) = val16(i) 
            tmpdtl(18) = val17(i) 
            tmpdtl(19) = val18(i) 
            tmpdtl(20) = val19(i) 
            tmpdtl(21) = val20(i) 
            tmpdtl(22) = val21(i) 
            tmpdtl(23) = val22(i) 
            tmpdtl(24) = val23(i) 
            tmpdtl(25) = val24(i) 
            call setparam_DTL(beamln3(idtl),tmpdtl)
            Blnelem(i) = assign_BeamLineElem(beamln3(idtl))
          else if(bitype(i).eq.102) then
            iccdtl = iccdtl + 1
            call construct_CCDTL(beamln4(iccdtl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = val0(i)
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i)
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_CCDTL(beamln4(iccdtl),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln4(iccdtl))
          else if(bitype(i).eq.103) then
            iccl = iccl + 1
            call construct_CCL(beamln5(iccl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = val0(i)
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i) 
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_CCL(beamln5(iccl),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln5(iccl))
          else if(bitype(i).eq.104) then
            isc = isc + 1
            call construct_SC(beamln6(isc),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = val0(i)
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i)
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_SC(beamln6(isc),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln6(isc))
          else if(bitype(i).eq.105) then
            islrf = islrf + 1
            call construct_SolRF(beamln8(islrf),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpslrf(1) = val0(i)
            tmpslrf(2) = val1(i) 
            tmpslrf(3) = val2(i) 
            tmpslrf(4) = val3(i) 
            tmpslrf(5) = val4(i) 
            tmpslrf(6) = val5(i) 
            tmpslrf(7) = val6(i) 
            tmpslrf(8) = val7(i) 
            tmpslrf(9) = val8(i) 
            tmpslrf(10) = val9(i) 
            tmpslrf(11) = val10(i) 
            tmpslrf(12) = val11(i) 
            call setparam_SolRF(beamln8(islrf),tmpslrf)
            Blnelem(i) = assign_BeamLineElem(beamln8(islrf))
          else if(bitype(i).eq.110) then
            iemfld = iemfld + 1
            call construct_EMfld(beamln11(iemfld),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmp13(1) = val0(i)
            tmp13(2) = val1(i) 
            tmp13(3) = val2(i) 
            tmp13(4) = val3(i) 
            tmp13(5) = val4(i) 
            tmp13(6) = val5(i) 
            tmp13(7) = val6(i) 
            tmp13(8) = val7(i) 
            tmp13(9) = val8(i) 
            tmp13(10) = val9(i) 
            tmp13(11) = val10(i) 
            tmp13(12) = val11(i) 
            tmp13(13) = val12(i) 
            call setparam_EMfld(beamln11(iemfld),tmp13)
            Blnelem(i) = assign_BeamLineElem(beamln11(iemfld))
          else if(bitype(i).eq.111) then
            iemfldcart = iemfldcart + 1
            call construct_EMfldCart(beamln12(iemfldcart),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmp13(1) = val0(i)
            tmp13(2) = val1(i)
            tmp13(3) = val2(i)
            tmp13(4) = val3(i)
            tmp13(5) = val4(i)
            tmp13(6) = val5(i)
            tmp13(7) = val6(i)
            tmp13(8) = val7(i)
            tmp13(9) = val8(i)
            tmp13(10) = val9(i)
            tmp13(11) = val10(i)
            call setparam_EMfldCart(beamln12(iemfldcart),tmp13)
            Blnelem(i) = assign_BeamLineElem(beamln12(iemfldcart))
          else if(bitype(i).eq.112) then
            iemfldcyl = iemfldcyl + 1
            call construct_EMfldCyl(beamln13(iemfldcyl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmp13(1) = val0(i)
            tmp13(2) = val1(i)
            tmp13(3) = val2(i)
            tmp13(4) = val3(i)
            tmp13(5) = val4(i)
            tmp13(6) = val5(i)
            tmp13(7) = val6(i)
            tmp13(8) = val7(i)
            tmp13(9) = val8(i)
            tmp13(10) = val9(i)
            tmp13(11) = val10(i)
            call setparam_EMfldCyl(beamln13(iemfldcyl),tmp13)
            Blnelem(i) = assign_BeamLineElem(beamln13(iemfldcyl))
          else if(bitype(i).eq.113) then
            iemfldana = iemfldana + 1
            call construct_EMfldAna(beamln14(iemfldana),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmp13(1) = val0(i)
            tmp13(2) = val1(i)
            tmp13(3) = val2(i)
            tmp13(4) = val3(i)
            tmp13(5) = val4(i)
            tmp13(6) = val5(i)
            tmp13(7) = val6(i)
            tmp13(8) = val7(i)
            tmp13(9) = val8(i)
            tmp13(10) = val9(i)
            tmp13(11) = val10(i)
            call setparam_EMfldAna(beamln14(iemfldana),tmp13)
            Blnelem(i) = assign_BeamLineElem(beamln14(iemfldana))
          else
          endif 
          zz = val0(i) + blength(i)
          zBlnelem(2,i) = zz
        enddo
!-------------------------------------------------------------------
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) print*,"pass setting up lattice..."

        deallocate(blength,bnseg,bmpstp,bitype)
        deallocate(val0)
        deallocate(val1,val2,val3,val4,val5,val6,val7,val8,val9)
        deallocate(val10,val11,val12,val13,val14,val15,val16)
        deallocate(val17,val18,val19,val20,val21,val22,val23)
        deallocate(val24)

!-------------------------------------------------------------------
! sample initial particle distribution.
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!set random number
        call random_seed(SIZE=seedsize)
        allocate(seedarray(seedsize))
        do i = 1, seedsize
            seedarray(i) = 10.0d0 + myid*Dim*20+i*1.0d0*myid
        enddo
        call random_seed(PUT=seedarray)
        do i = 1, 3000
           call random_number(rancheck)
        enddo
        deallocate(seedarray)
        print*,"check randomness: ",myid,rancheck

        !tend = 0.0d0
        tend = phsini
        dtlessend = 1.0d0
        isloutend = 0
        iend = 0
        ioutend = 0
        itszend = 0
        isteerend = 0
        if(Rstartflg.eq.1) then
          call inpoint_Output(nfileout+myid,Ebunch(1),tend,iend,ibchend,nprow,npcol,&
          Ageom,Nx,Ny,Nz,myidx,myidy,Np(1),ioutend,itszend,isteerend,isloutend,&
          dtlessend)
          if(myid.eq.0) print*,"restart at: ",tend,iend,ib
        else
          !call sample_Dist(Ebunch(1),distparam,Ndistparam,Flagdist,Ageom,grid2d,Flagbc)
          ib = 1
          call sample_Dist(Ebunch(1),distparam,Ndistparam,Flagdist,Ageom,grid2d,Flagbc,ib,Nbunch)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) print*,"pass generating initial distribution..."

!-------------------------------------------------------------------
! construct FieldQuant class objects.
        call construct_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d)

!        print*,"start to read in the inputs for the other bunches/bins:.."
        !generate the particle distribution for the other bunches/bins
        do ib = 2, Nbunch
          call in_Input(Dim,Np(ib),Nx,Ny,Nz,Flagbc,Flagdist,Rstartflg,&
              Flagmap,distparam,Ndistparam,Bcurr,Bkenergy,Bmass,Bcharge,&
          Bfreq,xrad,yrad,Perdlen,Nblem,npcol,nprow,Flagerr,Flagdiag,&
          Flagsubstep,phsini,dt,ntstep,Nbunch,FlagImage,ib)
          call construct_BeamBunch(Ebunch(ib),Bcurr,Bkenergy,Bmass,Bcharge,&
                            Np(ib),phsini)
          if(Rstartflg.eq.1) then
             call inpoint_Output(ib*nfileout+myid,Ebunch(ib),tend,iend,ibchend,&
                nprow,npcol,&
                Ageom,Nx,Ny,Nz,myidx,myidy,Np(ib),ioutend,itszend,isteerend,&
                isloutend,dtlessend)
          else
            !call sample_Dist(Ebunch(ib),distparam,Ndistparam,Flagdist,Ageom,grid2d,Flagbc)
            call sample_Dist(Ebunch(ib),distparam,Ndistparam,Flagdist,Ageom,grid2d,Flagbc,ib,Nbunch)
          endif
        enddo

        !get local particle number and mesh number on each processor.
        do ib = 1, Nbunch
          call getnpt_BeamBunch(Ebunch(ib),Nplocal(ib))
        enddo

        t_init = t_init + elapsedtime_Timer(t0)
        !this is just for finding the driven phase of each cavity
        !tend = phsini

        end subroutine init_AccSimulator

        !Run beam dynamics simulation through accelerator.
        subroutine run_AccSimulator()
        implicit none
        include 'mpif.h'
        integer :: i,j,bnseg,bmpstp,bitype
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,npx,npy,&
                   totnp,ierr,nylcr,nzlcr
        integer :: nbal,ibal,istep,ifile
        integer :: ibalend,istepend,nfile
        integer, dimension(3) :: lcgrid
        integer, allocatable, dimension(:) :: lctabnmx,lctabnmy
        integer, allocatable, dimension(:,:,:) :: temptab
        double precision :: z0,z,tau1,tau2,blength,t0,zend
        double precision, allocatable, dimension(:,:) :: lctabrgx, lctabrgy
        double precision, dimension(6) :: lcrange,range,ptrange,sgcenter,grange
        double precision, dimension(3) :: msize
        double precision :: hy,hz,ymin,zmin,piperad,zedge,hzwake
        double precision :: tmp1,tmp2,tmp3,tmp4,rfile
        double precision, allocatable, dimension(:,:,:) :: chgdens,tmppot
        double precision, allocatable, dimension(:,:,:) :: exg,eyg,ezg
        double precision, allocatable, dimension(:,:,:) :: bxg,byg,bzg
        double precision, allocatable, dimension(:,:,:) :: besscoef
        double precision, allocatable, dimension(:,:) :: bessnorm,gml
        integer, allocatable, dimension(:) :: modth,pydisp
        integer :: nmod,k
        !double precision :: sumtest, sumtest2, sumtest3
        double precision, dimension(8) :: drange
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: realSamplePeriod
        integer :: nsubstep,integerSamplePeriod
        double precision :: zcent,distance,blnLength,dzz
        integer, allocatable, dimension(:,:) :: idrfile
        integer :: ibend,ibstart,isw,ibinit,ibendold,iifile,ii,ibinitold,idbeamln
        double precision :: zmax,t,dtless,zshift,gammazavg,curr
        integer :: tmpflag,ib,ibb,ibunch,inib,nplctmp,nptmp,nptottmp
        double precision, allocatable, dimension(:) :: gammaz
        double precision, allocatable, dimension(:,:) :: brange
        double precision :: dGspread
        integer, dimension(Maxoverlap) :: tmpfile
        double precision :: tmpcur,totchrg,r0
        integer :: flagpt2pt,flagpos, flagcathode !//switch for point-to-point calculation and switch for space-charge with z>0
        double precision :: zz,zorgin,zorgin2,vref,gamin,gam
        double precision, dimension(6) :: ptref
        integer :: idbd,idbend,flagbctmp
        integer :: npttmplc,npttmp
        double precision :: ztmp1,ztmp2,deltaz
        integer :: ipt,iptnew
        double precision :: betazini
        !//switch for output the phase space plot: 0 - no output, 1 - output 
        !//6D particle phase space at given time, 2 - output for restart 
        integer :: flagphout,iout
        !//switch for changing time step size: 0 - no change, 1 - change
        !// time step size after tszstart
        integer :: flagtimesz
        double precision :: trstart
        integer, dimension(100) :: nsamp,nfileouttmp,nfileslout,nslout
        double precision, dimension(101) :: tszstart,dtnewsize
        double precision, dimension(101) :: tphout,tsteer,xoffset,yoffset,&
                                     pxoffset,pyoffset,zoffset,pzoffset,tslout
        integer :: flagazmuth
        !//the following parameters are defined for the wakefield
        !//calculation: only short range transverse dipole, longitudinal monopole wakes
        !//are included. This follows K. Bane's paper (SLAC-PUB-9663), March 2003, for
        !//forward traveling wave structures using three parameters, iris radius a, gap g, 
        !//period L. For backward traveling wave structures, the wakes are hardwired inside
        !//the code following the report: P. Craievich, T. Weiland, I. Zagorodnov, "The
        !//short-range wakefields in the BTW accelerating structure of the ELETTRA linac,"
        !//ST/M-04/02.
        integer :: flagwake,kz,jy,ix,jadd,kadd,kst,jst,nz2,isteer,iizz
        double precision, allocatable, dimension(:) :: ztable,zdisp,&
            denszlc,densz,exwake,eywake,ezwake,xwakelc,xwakez,ywakelc,ywakez,&
            denszp,denszpp,ans,csg
        double precision, allocatable, dimension(:,:) :: sendensz,&
                                                         recvdensz
        double precision :: xx,yy,t3dstart,rr,tmger,tmpwk,tstop
        double precision, dimension(100) :: aawk,ggwk,lengwk,twkstart,twkend
        double precision :: aawk1,ggwk1,lengwk1,tbtwstart
        double precision, allocatable, dimension(:,:)  :: tmppts
        double precision, allocatable, dimension(:,:,:)  :: tmpptsb
        integer :: iwk,nwk,flagbtw
        !//for the 1D csr wake field
        integer :: iizz1,flagcsr
        double precision :: eeff,ssll,r0bend,kkcsr,Bybend
        !//for N body
        integer  :: nplocal0,np0
        integer :: itsz
        real*8 :: bendlen,zz1,zz2,zwkmin,poscent,zfmin,zfmax,coeftol
        integer :: ldsg,nlsg,nrsg,npsg,msg,ncoefreal,iz
        !//for readin wake function
        integer :: ndatawk,ndatawkmax,itmp
        real*8, allocatable, dimension(:,:) :: wklong2d,wktran2d
        real*8, allocatable, dimension(:) :: wklong,wktran
        integer :: flagstep,flagspc,itspc
        double precision, dimension(101) :: tspcstart,vspc
        real*8, dimension(3) :: tmpfld
        real*8, dimension(4) :: pos
        integer :: nplocalmax,islout
!------------------------------------------------
!for longitudinal mesh adjustment
        real*8 :: zrandi,zadjust,zadjmax
!for first order emission model
        real*8 :: tmpzz,dtmp,recpgam
        integer :: totnpts
        real*8 :: qchg 
!for collimator
        integer :: icol
        real*8, dimension(101) :: tcol,xradmin,xradmax,yradmin,yradmax
!for instant applying linear transfer matrix to the beam
        real*8, dimension(101) :: tmap
        integer :: imap
        real*8 :: tmpx,tmppx,tmpy,tmppy,tmpz,tmppz
        real*8, dimension(6,6,101) :: rmt
!for instant applying the uncorrelated energy heater
        real*8, dimension(101) :: theat,siggambetz
        real*8 :: rdnorm,twopi,x1,x2
        integer :: iheat
!for instant applying the rotation along the z-axis
        real*8, dimension(101) :: tzrot,zrotang
        integer :: izrot
!for switch integrator: -18
        integer :: flaginteg = 0
        real*8, dimension(101) :: tinteg
!--------------
! dielectric wakefield module implemented by Daniel Mihalcea:
! PRST-AB, 15, 081304, (2012).
! for DWA:
        integer :: iwrite,dflagwake,dndatawk,dnwk,dfound
        double precision, allocatable, dimension(:,:,:) ::dexwake,deywake,&
                   dezwake,dbxwake,dbywake,dbzwake
        double precision, dimension(100) :: dtype,da,db,dLx,dtwkstart&
                                           ,dtwkend,deps,dnkx,dnky
        double precision, dimension(50,50,2) :: dky,amp
        double precision, dimension(50) :: dkx
        double precision :: dtype1,da1,db1,deps1,dLx1,dnkx1,dnky1,hx,xmi,ymi


        twopi = 4*asin(1.0d0)
        !zadjmax = 0.15 !30% increase of z domain
        zadjmax = 0.0d0 !0% increase of z domain

        ndatawkmax = 1000

        flagspc = 1
        flagbtw = 0 !initial no backward tws wakefield
        flagcsr = 1
        flagazmuth = 1
        if(Nemission.gt.0) then
          flagcathode = 1
        else !no cathode model
          flagcathode = 0
          FlagImage = 0
        endif
        !flagazmuth = 0
!-------------------------------------------------------------------
! prepare initial parameters, allocate temporary array.
        !ibalend = 0
        ibalend = 1
        istepend = 0
        zend = 0.0
!warning: in the point-to-point space charge force model, the
!choice of the cut-off radius, r0, for each particle has to be very careful.
!otherwise, it can end up with untrustable results.
        r0 = 0.0/Scxlt
        !flagpt2pt = 1 !"0" not point-2-point, "1" use point-2-point
        flagpt2pt = 0 !"0" not point-2-point, "1" use point-2-point

        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        allocate(lctabnmx(0:npx-1))
        allocate(lctabnmy(0:npy-1))
        allocate(lctabrgx(2,0:npx-1))
        allocate(lctabrgy(2,0:npy-1))
        allocate(temptab(2,0:npx-1,0:npy-1))

        nbal = 100
        ibal = ibalend
        istep = istepend
        z = zend

        !assign initial storage for charge density,Ex,Ey,Ez,Bx,By,Bz,image
        !charge potential. 
        allocate(chgdens(1,1,1))
        allocate(exg(1,1,1))
        allocate(eyg(1,1,1))
        allocate(ezg(1,1,1))
        allocate(bxg(1,1,1))
        allocate(byg(1,1,1))
        allocate(bzg(1,1,1))
        allocate(tmppot(1,1,1))

!!! DWA
        allocate(dexwake(1,1,1))
        allocate(deywake(1,1,1))
        allocate(dezwake(1,1,1))
        allocate(dbxwake(1,1,1))
        allocate(dbywake(1,1,1))
        allocate(dbzwake(1,1,1))
!!!!

        !assign storage used for wakefield calculation
        allocate(ztable(0:npx-1))
        allocate(zdisp(0:npx-1))

        allocate(denszlc(Nz)) 
        allocate(densz(Nz)) 
        allocate(denszp(Nz)) 
        allocate(denszpp(Nz)) 
        allocate(ans(Nz)) 
        allocate(csg(Nz)) 
        allocate(xwakelc(Nz))
        allocate(xwakez(Nz))
        allocate(ywakelc(Nz))
        allocate(ywakez(Nz))
        allocate(recvdensz(Nz,2))
        allocate(sendensz(Nz,2))
        allocate(exwake(Nz))
        allocate(eywake(Nz))
        allocate(ezwake(Nz))
        allocate(wklong2d(ndatawkmax,100)) 
        allocate(wktran2d(ndatawkmax,100)) 
        allocate(wktran(ndatawkmax)) 
        allocate(wklong(ndatawkmax)) 
        wklong2d = 0.0
        wklong = 0.0
        wktran2d = 0.0
        wktran = 0.0
        exwake = 0.0
        eywake = 0.0
        ezwake = 0.0
        flagwake = 0
        aawk1 = 0.05
        ggwk1 = 0.05
        lengwk1 = 0.1

!!! DWA
        dflagwake=0
        dnwk=0
!!! end DWA

        !no phase space print out
        flagphout = 0
        tphout = 1.0e12 
        tslout = 1.0e12
        !no change of the time step size
        flagtimesz = 0
        tszstart = 1.0e12
        tspcstart = 1.0e12
        dtnewsize = 1.0e12
        trstart = 1.0e12
        t3dstart = 1.0e12
        twkstart = 1.0e12
        twkend = 1.0e12
        tsteer = 1.0e12 
        tmger = 1.0e12
        tstop = 1.0e12
        tcol = 1.0e12
        tmap = 1.0e12
        theat = 1.0e12
        tzrot = 1.0e12

        izrot = 0
        iheat = 0
        imap = 0       
        islout = 0
        iout = 0
        itsz = 0
        itspc = 0
        nwk = 0
        isteer = 0
        icol = 0

        !idrfile is used to store the <element type>, <external data file name>,
        !and <id> for the internal data storage of each beamline element  
        allocate(idrfile(3,Nblem))
        idrfile = 1
        idrfile(2,:) = -10
        do i = 1, Nblem
          call getparam_BeamLineElem(Blnelem(i),blength,bnseg,bmpstp,&
                                     bitype)
          idrfile(1,i) = bitype
          !get external file id for each rf beam line element.
          if(bitype.gt.100) then
            call getparam_BeamLineElem(Blnelem(i),5,rfile)
            idrfile(2,i) = int(rfile + 0.1)
          !get external file for solenoid.
          else if(bitype.eq.3) then
            call getparam_BeamLineElem(Blnelem(i),3,rfile)
            idrfile(2,i) = int(rfile + 0.1)
          else if(bitype.eq.4) then
            call getparam_BeamLineElem(Blnelem(i),4,rfile)
            idrfile(2,i) = int(rfile + 0.1)
          endif
          if(bitype.eq.(-1)) then
            isteer = isteer + 1
            if(isteer.gt.100) then
              print*,"The maximum steering location is 100!!!"
              isteer = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),2,tsteer(isteer))
            call getparam_BeamLineElem(Blnelem(i),3,xoffset(isteer))
            call getparam_BeamLineElem(Blnelem(i),4,pxoffset(isteer))
            call getparam_BeamLineElem(Blnelem(i),5,yoffset(isteer))
            call getparam_BeamLineElem(Blnelem(i),6,pyoffset(isteer))
            call getparam_BeamLineElem(Blnelem(i),7,zoffset(isteer))
            call getparam_BeamLineElem(Blnelem(i),8,pzoffset(isteer))
            xoffset(isteer) = xoffset(isteer)/Scxlt
            yoffset(isteer) = yoffset(isteer)/Scxlt
            zoffset(isteer) = zoffset(isteer)/Scxlt
          endif
          if(bitype.eq.(-2)) then
            flagphout = 1
            iout = iout + 1
            if(iout.gt.100) then
              print*,"The maximum phase space output is 100!!!"
              iout = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),3,tphout(iout))
            nfileouttmp(iout) = bmpstp 
            nsamp(iout) = bnseg
          endif
          if(bitype.eq.(-3)) then
            flagphout = 2
            call getparam_BeamLineElem(Blnelem(i),3,trstart)
            nfileout = bmpstp 
!            nsamp = bnseg
            print*,"trstart: ",trstart
          endif
          if(bitype.eq.(-4)) then
            flagtimesz = 1
            itsz = itsz + 1
            if(itsz.gt.100) then
              print*,"The maximum time step size change is 100!!!"
              itsz = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),3,tszstart(itsz))
            call getparam_BeamLineElem(Blnelem(i),4,dtnewsize(itsz))
            !switch for different integrators
            call getparam_BeamLineElem(Blnelem(i),5,tinteg(itsz))
          endif
          if(bitype.eq.(-5)) then
            call getparam_BeamLineElem(Blnelem(i),3,t3dstart)
          endif
          if(bitype.eq.(-6)) then
            nwk = nwk + 1
            if(nwk.gt.100) then
              print*,"The maximum wakefield cavity is 100!!!"
              nwk = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),3,twkstart(nwk))
            call getparam_BeamLineElem(Blnelem(i),4,twkend(nwk))
            call getparam_BeamLineElem(Blnelem(i),5,aawk(nwk))
            call getparam_BeamLineElem(Blnelem(i),6,ggwk(nwk))
            call getparam_BeamLineElem(Blnelem(i),7,lengwk(nwk))
            !read in the wake function from external file: fort.bmpstp
            if(bnseg.gt.0) then
              !count the number of data points in wake function
              aawk(nwk) = 0.0
              do itmp = 1, ndatawkmax
                read(bmpstp,*,end=111)wklong2d(itmp,nwk),wktran2d(itmp,nwk)
                aawk(nwk) = aawk(nwk) + 1
              enddo
111           continue
              close(bmpstp)
            else
            !use the analytical formulae for wakefield.
              if(aawk(nwk).gt.100) then !BTWS
                flagbtw = 1
              else if(aawk(nwk).lt.0.0) then
                if(aawk(nwk).gt.-10) then !1.3 GHz
                  flagbtw = 2
                else !3.9 GHz
                  flagbtw = 3
                endif
              else
                flagbtw = 0
              endif

            endif
            if(aawk(nwk).gt.100) then
              tbtwstart = ggwk(nwk)
            else
              tbtwstart = 1.d20
            endif
          endif
          if(bitype.eq.(-7)) then
            call getparam_BeamLineElem(Blnelem(i),3,tmger)
          endif
          if(bitype.eq.(-8)) then
            itspc = itspc + 1
            if(itspc.gt.100) then
              print*,"The maximum space-charge change is 100!!!"
              itspc = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),2,vspc(itspc))
            call getparam_BeamLineElem(Blnelem(i),3,tspcstart(itspc))
          endif
          !output slice information (current,emittance,energy spread,etc) at given location.
          if(bitype.eq.(-9)) then
            islout = islout + 1
            if(islout.gt.100) then
              print*,"The maximum slice output is 100!!!"
              islout = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),3,tslout(islout))
            nfileslout(islout) = bmpstp
            nslout(islout) = bnseg
          endif

!collimator information
          if(bitype.eq.(-11)) then
            icol = icol + 1
            if(icol.gt.100) then
              print*,"The maximum # of collimator is 100!!!"
              icol = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),2,tcol(icol))
            call getparam_BeamLineElem(Blnelem(i),3,xradmin(icol))
            call getparam_BeamLineElem(Blnelem(i),4,xradmax(icol))
            call getparam_BeamLineElem(Blnelem(i),5,yradmin(icol))
            call getparam_BeamLineElem(Blnelem(i),6,yradmax(icol))
          endif

!transfer matrix information
          if(bitype.eq.(-12)) then
            imap = imap + 1
            if(imap.gt.100) then
              print*,"The maximum # of transfer matrix is 100!!!"
              imap = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),1,tmap(imap))
            open(11,file="linearmap.in",status="old")
            do ii = 1, (imap-1)*6
              read(11,*)
            enddo
            do ii = 1, 6
              read(11,*)rmt(ii,1,imap),rmt(ii,2,imap),rmt(ii,3,imap),&
                       rmt(ii,4,imap),rmt(ii,5,imap),rmt(ii,6,imap)
            enddo
            close(11)
          endif

!!! DWA
          if(bitype.eq.(-13)) then
            dnwk = dnwk + 1
            if(dnwk.gt.100) then
              print*,"The maximum dielectric wakefield cavity is 100!!!"
              dnwk = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),1,dtwkstart(dnwk))
            call getparam_BeamLineElem(Blnelem(i),2,dtwkend(dnwk))
            call getparam_BeamLineElem(Blnelem(i),3,dtype(dnwk))
            call getparam_BeamLineElem(Blnelem(i),4,da(dnwk))
            call getparam_BeamLineElem(Blnelem(i),5,db(dnwk))
            call getparam_BeamLineElem(Blnelem(i),6,deps(dnwk))
            call getparam_BeamLineElem(Blnelem(i),7,dLx(dnwk))
            call getparam_BeamLineElem(Blnelem(i),8,dnkx(dnwk))
            call getparam_BeamLineElem(Blnelem(i),9,dnky(dnwk))
          endif
!!! end DWA

          if(bitype.eq.(-15)) then !turn on pt2pt for SC calculation
            call getparam_BeamLineElem(Blnelem(i),3,r0)
            flagpt2pt = 1
          endif

          if(bitype.eq.(-16)) then !instant gamma beta_z heating
            iheat = iheat + 1
            call getparam_BeamLineElem(Blnelem(i),1,theat(iheat))
            call getparam_BeamLineElem(Blnelem(i),2,siggambetz(iheat))
          endif

          if(bitype.eq.(-17)) then !instant rotation along z-axis
            izrot = izrot + 1
            call getparam_BeamLineElem(Blnelem(i),1,tzrot(izrot))
            call getparam_BeamLineElem(Blnelem(i),2,zrotang(izrot))
          endif

          if(bitype.eq.(-99)) then
            call getparam_BeamLineElem(Blnelem(i),3,tstop)
            print*,"tstop: ",tstop
          endif
        enddo

!        Ebunch(1)%Pts1 = 1.0d0
!        call phaseinadios(Ebunch(1))

        dtless = dtlessend !dimensionless time step size.
        t = tend
        distance = 0.0d0
        tmpfile = 0
        !length of the total beamline
        blnLength = zBlnelem(2,Nblem)
        ibinitold = 1
        ibendold = 0
        iifile = 0
        islout = isloutend
        imap = 0
        iheat = 0
        izrot = 0

        allocate(gammaz(Nbunch))
        allocate(brange(12,Nbunch))
        !count the total current and # of particles and local # of particles for each 
        !bunch or bin
        curr = 0.0d0
        do ib = 1, Nbunch
          curr = curr + Ebunch(ib)%Current
          Nplocal(ib) = Ebunch(ib)%Nptlocal
          Np(ib) = Ebunch(ib)%Npt
        enddo
        nplocal0 = Nplocal(1)
        np0 = Np(1)
        if(Flagdiag.eq.1) then
          !output the moments from the average of all bunches at fixed t.
          call diagnostic1avg_Output(t,Ebunch,Nbunch)
        else if(Flagdiag.eq.2) then
          !output the moments from the average of all bunches at fixed z.
          call diagnostic1avgZ_Output(t,Ebunch,Nbunch)
        else if(Flagdiag.eq.3) then
          !output the moments from the average of all bunches at fixed z.
          call diagnostic1avgZtest_Output(t,Ebunch,Nbunch)
        else
        endif
!        ibunch = 0
!        iout = 0
!        itsz = 0
!        isteer = 0
        itspc = 0
        iout = ioutend
        itsz = itszend
        isteer = isteerend
!        ibunch = 1
        ibunch = Nbunch
        tmpcur = curr
        flagpos = 1
        flagstep = 1
        zorgin = 0.0d0
        idbd = 1
        icol = 0
        !particles behind the cathode will use this beta for emission.
        !betazini = sqrt(1.0-1.0/(1.0+distparam(21)**2))
        betazini = sqrt(1.0d0-1.0d0/(1.0d0+Bkenergy/Bmass)**2)
        !output initial phase space distribution 
        do ib = 1, Nbunch
          call phase_Output(40+ib-1,Ebunch(ib),1)
          qchg = Ebunch(ib)%Current/Scfreq
          call sliceprocdep_Output(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                   Nz,qchg,Ebunch(ib)%Mass,60+ib-1)
        enddo
        dzz = betazini*Clight*dtless*Dt
        zmin = 0.0d0
        call MPI_BARRIER(comm2d,ierr)
        !print*,"iout: ",iout,tphout(1)
!----------------------------------------------------------------------
! start looping through ntstep time step.
        !iend is the time step number from last simulation (used in
        !restart function).
        do i = iend+1, ntstep
          if(myid.eq.0) then
            print*,"i,t,<z>: ",i,t,distance
          endif

          !switch on the backward traveling wave structure
          if(t.ge.tbtwstart) then
            flagbtw = 1
          endif

          !steering the beam centroid to the given X, Px,Y, Py, Z, Pz values at given location
          !if(t.le.tsteer(isteer+1) .and. (t+dtless*Dt).ge.tsteer(isteer+1)) then
          if(distance.le.tsteer(isteer+1) .and. (distance+dzz).ge.tsteer(isteer+1)) then
            isteer = isteer + 1
            do ib = 1, Nbunch
              !//find the range and center information of each bunch/bin
              call singlerange(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                             ptrange,sgcenter)
              do ipt = 1, Nplocal(ib)
                 Ebunch(ib)%Pts1(1,ipt) = Ebunch(ib)%Pts1(1,ipt) - sgcenter(1) + xoffset(isteer)
                 Ebunch(ib)%Pts1(2,ipt) = Ebunch(ib)%Pts1(2,ipt) - sgcenter(2) + pxoffset(isteer)
                 Ebunch(ib)%Pts1(3,ipt) = Ebunch(ib)%Pts1(3,ipt) - sgcenter(3) + yoffset(isteer)
                 Ebunch(ib)%Pts1(4,ipt) = Ebunch(ib)%Pts1(4,ipt) - sgcenter(4) + pyoffset(isteer)
                 Ebunch(ib)%Pts1(5,ipt) = Ebunch(ib)%Pts1(5,ipt) - sgcenter(5) + zoffset(isteer)
                 Ebunch(ib)%Pts1(6,ipt) = Ebunch(ib)%Pts1(6,ipt) - sgcenter(6) + pzoffset(isteer)
              enddo
            enddo
          endif

          !apply instant linear transfer matrix to the beam
          if(distance.le.tmap(imap+1) .and. (distance+dzz).ge.tmap(imap+1)) then
            imap = imap + 1
            do ib = 1, Nbunch

            !//find the range and center information of each bunch/bin
              call singlerange(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                             ptrange,sgcenter)

              do ipt = 1, Nplocal(ib)
                 tmpx = Ebunch(ib)%Pts1(1,ipt)  - sgcenter(1)
                 tmppx = Ebunch(ib)%Pts1(2,ipt) - sgcenter(2)
                 tmpy = Ebunch(ib)%Pts1(3,ipt)  - sgcenter(3)
                 tmppy = Ebunch(ib)%Pts1(4,ipt) - sgcenter(4)
                 tmpz = Ebunch(ib)%Pts1(5,ipt)  - sgcenter(5)
                 tmppz = Ebunch(ib)%Pts1(6,ipt) - sgcenter(6)
                 Ebunch(ib)%Pts1(1,ipt)=tmpx*rmt(1,1,imap)+tmppx*rmt(1,2,imap)+&
                    tmpy*rmt(1,3,imap)+tmppy*rmt(1,4,imap)+&
                    tmpz*rmt(1,5,imap)+tmppz*rmt(1,6,imap)+ sgcenter(1)
                 Ebunch(ib)%Pts1(2,ipt)=tmpx*rmt(2,1,imap)+tmppx*rmt(2,2,imap)+&
                    tmpy*rmt(2,3,imap)+tmppy*rmt(2,4,imap)+&
                    tmpz*rmt(2,5,imap)+tmppz*rmt(2,6,imap)+ sgcenter(2)
                 Ebunch(ib)%Pts1(3,ipt)=tmpx*rmt(3,1,imap)+tmppx*rmt(3,2,imap)+&
                    tmpy*rmt(3,3,imap)+tmppy*rmt(3,4,imap)+&
                    tmpz*rmt(3,5,imap)+tmppz*rmt(3,6,imap)+ sgcenter(3)
                 Ebunch(ib)%Pts1(4,ipt)=tmpx*rmt(4,1,imap)+tmppx*rmt(4,2,imap)+&
                    tmpy*rmt(4,3,imap)+tmppy*rmt(4,4,imap)+&
                    tmpz*rmt(4,5,imap)+tmppz*rmt(4,6,imap)+ sgcenter(4)
                 Ebunch(ib)%Pts1(5,ipt)=tmpx*rmt(5,1,imap)+tmppx*rmt(5,2,imap)+&
                    tmpy*rmt(5,3,imap)+tmppy*rmt(5,4,imap)+&
                    tmpz*rmt(5,5,imap)+tmppz*rmt(5,6,imap)+ sgcenter(5)
                 Ebunch(ib)%Pts1(6,ipt)=tmpx*rmt(6,1,imap)+tmppx*rmt(6,2,imap)+&
                    tmpy*rmt(6,3,imap)+tmppy*rmt(6,4,imap)+&
                    tmpz*rmt(6,5,imap)+tmppz*rmt(6,6,imap)+ sgcenter(6)
              enddo
            enddo
          endif

          !instant heating
          if(distance.le.theat(iheat+1) .and. (distance+dzz).ge.theat(iheat+1)) then
            iheat = iheat + 1
            do ib = 1, ibunch
              do ipt = 1, Nplocal(ibunch)
                call random_number(x2)
                call random_number(x1)
                if(x1.eq.0.0d0) x1 = 1.0d-10
                rdnorm = sqrt(-2.0d0*log(x1))*cos(twopi*x2)
                Ebunch(ib)%Pts1(6,ipt) = Ebunch(ib)%Pts1(6,ipt) +  &
                                         siggambetz(iheat)*rdnorm
              enddo
            enddo
          endif

          !instant z-rotation
          if(distance.le.tzrot(izrot+1) .and. (distance+dzz).ge.tzrot(izrot+1)) then
            izrot = izrot + 1
            do ib = 1, ibunch
              do ipt = 1, Nplocal(ibunch)
                 tmpx = Ebunch(ib)%Pts1(1,ipt)
                 tmppx = Ebunch(ib)%Pts1(2,ipt)
                 tmpy = Ebunch(ib)%Pts1(3,ipt)
                 tmppy = Ebunch(ib)%Pts1(4,ipt)
                 Ebunch(ib)%Pts1(1,ipt) = tmpx*cos(zrotang(izrot))+tmpy*sin(zrotang(izrot))
                 Ebunch(ib)%Pts1(2,ipt) = tmppx*cos(zrotang(izrot))+tmppy*sin(zrotang(izrot))
                 Ebunch(ib)%Pts1(3,ipt) = -tmpx*sin(zrotang(izrot))+tmpy*cos(zrotang(izrot))
                 Ebunch(ib)%Pts1(4,ipt) = -tmppx*sin(zrotang(izrot))+tmppy*cos(zrotang(izrot))
              enddo
            enddo
          endif

          if(i.le.Nemission) then !first Nemission steps for emission
            dtless = temission/Nemission/Dt
          else
            if(flagstep .eq. 1) then
              if(iend.le.Nemission) then
                dtless = 1.0d0
              else
                dtless = dtlessend
              endif
            endif
          endif
          if(zmin.gt.0.0) then
            nbal = 50
            !nbal = 50000
            !nbal = 500
            !nbal = 5
          else
            nbal = 100000
            !nbal = 100
          endif

          !change time step size at tszstart
          !if(t.le.tszstart(itsz+1) .and. (t+dtless*Dt).ge.tszstart(itsz+1)) then
          if(distance.le.tszstart(itsz+1) .and. (distance+dzz).ge.tszstart(itsz+1)) then
            itsz = itsz + 1
            dtless = dtnewsize(itsz)/Dt
            flagstep = 0
            !change the type of the integrator
            flaginteg = tinteg(itsz) + 0.000001
          endif

!          print*,"flatinteg: ",flaginteg,tinteg(itsz),dtless,itsz0

          !change charge space at tspcstart
          !if(t.le.tspcstart(itspc+1) .and. (t+dtless*Dt).ge.tspcstart(itspc+1)) then
          if(distance.le.tspcstart(itspc+1) .and. (distance+dzz).ge.tspcstart(itspc+1)) then
            itspc = itspc + 1
            if(vspc(itspc).gt.0.0) then
              flagspc = 1
            else
              flagspc = -1
            endif
          endif

          !turn off 2D azmuthal symmetry deposition and use full 3d deposition.
          !if(t.ge.t3dstart) flagazmuth = 0
          if(distance.ge.t3dstart) flagazmuth = 0

          !turn on the wakefield between twkstart and twkend.
          do iwk = 1, nwk
            !if(t.ge.twkstart(iwk) .and. t.le.twkend(iwk)) then
            if(distance.ge.twkstart(iwk) .and. distance.le.twkend(iwk)) then
              flagwake = 1
              aawk1 = aawk(iwk)
              ggwk1 = ggwk(iwk)
              lengwk1 = lengwk(iwk)
              if(ggwk1.le.0.0) then 
                ndatawk = aawk(iwk) + 0.1
                wklong(:) = wklong2d(:,iwk)
                wktran(:) = wktran2d(:,iwk)
              endif
              exit
            else
              flagwake = 0
            endif
          enddo

!!! DWA
          !turn on the DWA between dtwkstart and dtwkend.
          dfound=0
          do iwk = 1, dnwk
            !if(t.ge.twkstart(iwk) .and. t.le.twkend(iwk)) then
           if(dfound.eq.0)then
            if(distance.ge.dtwkstart(iwk) .and. distance.le.dtwkend(iwk))then
              dflagwake = 1
              dtype1=dtype(iwk)
              da1=da(iwk)
              db1=db(iwk)
              deps1=deps(iwk)
              dLx1=dLx(iwk)
              dnkx1=dnkx(iwk)
              dnky1=dnky(iwk)
              if(dtype1.eq.0)then
                call dispEqSlab_FieldQuant(da1,db1,deps1,dLx1,dnkx1,dnky1,&
                     dkx,dky)
                call ampSlab_FieldQuant(da1,db1,deps1,dnkx1,dnky1,dkx,dky,&
                     amp)
              else
                call dispEqCyl_FieldQuant(da1,db1,deps1,dnkx1,dnky1,dky)
                call ampCyl_FieldQuant(da1,db1,deps1,dnkx1,dnky1,dky,amp)
              endif
              dfound=1
            else
              dflagwake = 0
            endif
           endif
          enddo
!!! end DWA

          !//update particle positions using velocity for half step.
          if(flagcathode.eq.1) then
            do ib = 1, Nbunch   
              call drifthalf_BeamBunch(Ebunch(ib),t,dtless,betazini)
            enddo   
          else
            do ib = 1, Nbunch   
              call drifthalforg_BeamBunch(Ebunch(ib),t,dtless)
            enddo   
          endif
          t = t + 0.5*dtless*Dt

          !ibunch is total number of bunches within the effective computational domain
          !only bunch/bin with zmax>0 is counted as effective bunch/bin 
          if(ibunch.lt.Nbunch) then
            !ibunch = 0
            do ib = 1, Nbunch
              !//find the range and center information of each bunch/bin
              call singlerange(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                             ptrange,sgcenter)
              if(ptrange(6).gt.0.0) then
                !ibunch = ibunch + 1
              endif 
            enddo
          endif

          !only the bunch with zmax>0 is counted as an effective bunch 
          do ib = 1, ibunch
            !//find the range and center of each bunch/bin
            call singlerange(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                             ptrange,sgcenter)
            !Ebunch(ib)%refptcl(5) = sgcenter(5) + zorgin
            Ebunch(ib)%refptcl(5) = sgcenter(5) 
            gammaz(ib) = sqrt(1.0+sgcenter(6)**2)!//gammaz from <gamma_i betaz_i>
            !for test
            !gammaz(ib) = Bkenergy/Bmass+1.0d0
            Ebunch(ib)%refptcl(6) = -gammaz(ib)
            do inib = 1, 6
              brange(inib,ib) = ptrange(inib)
              brange(inib+6,ib) = sgcenter(inib)
            enddo
            !the longitudinal range for space-charge calculation has to be > 0
!            if(flagpos.eq.1) then
!J.Q. 10/27/08
              if(ptrange(5).gt.0) then
                brange(5,ib) = ptrange(5)
              else
                brange(5,ib) = 0.0
              endif

!            else
!              brange(5,ib) = ptrange(5)
!            endif
          enddo
          !here, ptrange(5) is the zminlc of the last bunch
          !this is the criterion for all bunches out of the cathode
          !if(ptrange(5).gt.0.0  .or. (flagcathode.eq.0) ) flagpos = 0
          !get the global computational domain and center for all effective bunches/bins 
!          print*,"brange: ",brange,Np
          call globalrange(brange,grange,gammazavg,zcent,Np,ibunch) 
          if(grange(5).gt.0.0d0  .or. (flagcathode.eq.0) ) flagpos = 0
          if(grange(5).gt.0.0d0) flagcathode = 0
!          print*,"zcent: ",zcent

!          print*,"gamz: ",gammaz(1)
          !get the distance of the center of all effective bunches/bins
          distance = zcent*Scxlt
          dzz = sqrt(1.0d0-1.0d0/gammazavg**2)*Clight*dtless*Dt
          nptottmp = sum(Np)
          !exit if the beam is outside the beamline
          if(distance.gt.blnLength .or. distance.gt.tstop .or. nptottmp.lt.1) then
            exit
          endif

          !check the particles outside the computational domain
          do ib = 1, ibunch
            call lost_BeamBunch(Ebunch(ib),xrad,yrad,Perdlen,zcent,&
                                nplctmp,nptmp)
            Nplocal(ib) = nplctmp
            Np(ib) = nptmp
!!            print*,"npt: ",ib,myid,nplctmp,nptmp
          enddo 

          totnpts = sum(Np)
          if(totnpts.le.0) exit

          !find the beginning and end beam line elements that the effective
          !bunch particles occupy
          !zmin = grange(5)*Scxlt + zorgin
          !zmax = grange(6)*Scxlt + zorgin
          zmin = grange(5)*Scxlt 
          zmax = grange(6)*Scxlt
!          if(zmax.le.0.0) goto 1000 !if no particles emittted pass field calcuation

          if(myid.eq.0) print*,"zmin,zmax: ",zmin,zmax
          !using the following way, we can find the id of the min element
          !and the max. element between which the bunch stays. It works even
          !if there are more than one element overlaps each other. However,
          !the beam line element has to be aranged so that zmin(i+1)>=zmin(i). 
          idbeamln = 1
          do ii = 1, Nblem
            if( (zmin.ge.zBlnelem(1,ii)) .and. &
                (zmin.le.zBlnelem(2,ii)) ) then
               idbeamln = ii
               exit !exit the loop from the first element id containing zmin
            endif
          enddo
          ibinit = idbeamln
          !exit the loop from the last element id containing zmax
          do ii = ibinit, Nblem
            if( (zmax.ge.zBlnelem(1,ii)) .and. &
                (zmax.le.zBlnelem(2,ii)) ) then
               idbeamln = ii 
            endif
          enddo
          ibend = idbeamln
          !print*,"zblnelem: ",i,zmin,zmax,zBlnelem(1,1:3),zBlnelem(2,1:3)
          
          !In the following, we will find the id in the internal global data array
          !and assign it to each beam line element.
          !This information is used to get the external field at given position.
          !This should also work for the overlaped beam line element.
!          if(i.eq.1) then
!           ibstart = max(ibendold,ibinit) 
!         else
!           ibstart = max(ibendold,ibinit) + 1
!         endif
          if(ibendold.ge.ibinit) then
            ibstart = ibendold + 1
          else
            ibstart = ibinit
          endif

!          print*,"ibinit: ",ibinit,ibend,ibstart,ibend,ibendold
          idbend = 0
          do ii = ibstart,ibend
            !for element type > 100 or solenoid
            if((idrfile(1,ii).gt.100).or.(idrfile(1,ii).eq.3).or.&
               (idrfile(1,ii).eq.4)) then 
              isw = 1
              !check whether the new element is in the old elements range,
              !which already been read in. 
              do j = ibinitold,ibendold
                if(idrfile(2,ii).ne.idrfile(2,j)) then
                else
                  isw = 0
                  exit
                endif
              enddo
              !check whether the new element file has been read in by the
              !old element.
              do j = 1, mod(iifile-1,Maxoverlap)+1
                if(idrfile(2,ii).ne.tmpfile(j)) then
                else
                  isw = 0
                  idrfile(3,ii) = j
                  exit
                endif
              enddo
              !accumulate the external file 1d in the global data array
!              if(isw.eq.1) iifile = iifile + 1
!              !assign the id in global data array to each beam line element
!              idrfile(3,ii) = mod(iifile-1,Maxoverlap)+1
              if(isw.eq.1) then
                iifile = iifile + 1
                !assign the id in global data array to each beam line element
                idrfile(3,ii) = mod(iifile-1,Maxoverlap)+1
              endif
              tmpfile(idrfile(3,ii)) = idrfile(2,ii)
!              print*,"idrfile: ",ii,idrfile(2,ii),idrfile(3,ii),iifile,Maxiifile
              if(isw.eq.1) then !readin new data
              if(idrfile(1,ii).eq.3) then !numerical data for solenoid.
                call read2tsol_Data(fldmp(idrfile(3,ii)),idrfile(2,ii))
              else if(idrfile(1,ii).lt.110) then !Fcoef coefficient description field
                call read1t_Data(fldmp(idrfile(3,ii)),idrfile(2,ii))
              else if(idrfile(1,ii).eq.111) then !3D Cartesian coordinate of field
                call read3t_Data(fldmp(idrfile(3,ii)),idrfile(2,ii))
              else if(idrfile(1,ii).eq.112) then !azmuthal symmetric cylindrical coordinate
                call read2t_Data(fldmp(idrfile(3,ii)),idrfile(2,ii))
              else if(idrfile(1,ii).eq.113) then !analytical function description of field
                if(idrfile(2,ii).eq.10) then
                  call read1t_Data(fldmp(idrfile(3,ii)),idrfile(2,ii))
                else if(idrfile(2,ii).gt.0) then
                  call read4t_Data(fldmp(idrfile(3,ii)),idrfile(2,ii))
                endif
              else if(idrfile(1,ii).eq.4) then
                call read1t_Data(fldmp(idrfile(3,ii)),idrfile(2,ii))
              else
                print*,"wrong element type to read in...."
                stop
              endif
              endif
            endif
            if(idrfile(1,ii).eq.4) then
               !this is used to avoid inifinite loop at the exit of the bend
               !where some particles behind the reference particle are still
               !inside the "bend". here, the "bend" is defined as a short section
               !of drift + bend + a short section of drift. 
               if(idbd.ne.ii) then
                 idbend = 1
                 idbd = ii
               else
                 idbend = 0
               endif
            else
               idbend = 0
            endif
          enddo 
          ibinitold = ibinit
          ibendold = ibend

          if(idbend.ne.1) then !not bending magnet

          if(zmax.le.0.0) goto 1000 !if no particles emittted pass field calcuation

          !calculate space-charge force for curr > 0 
          !if((curr.gt.0.0).and.(flagspc.eq.1)) then !solve the Poission equation 
          if((curr.gt.0.0)) then !solve the Poission equation 
            if(flagpt2pt.ne.1) then !using mesh solver

            call random_number(zrandi)
            zadjust = zadjmax*zrandi
            grange(5) = grange(5)-(grange(6)-grange(5))*zadjust
            grange(6) = grange(6)+(grange(6)-grange(5))*zadjust

            ! get new boundary from the range of beam particles.
            if(Flagbc.eq.4) then
            else
              !the use of grange instead of ptrange is due to multiple bunch
              call update_CompDom(Ageom,grange,grid2d,Flagbc)
            endif
            call getlcmnum_CompDom(Ageom,lcgrid)
            Nxlocal = lcgrid(1)
            if(npy.gt.1) then
              Nylocal = lcgrid(2) + 2
            else
              Nylocal = lcgrid(2)
            endif
            if(npx.gt.1) then
              Nzlocal = lcgrid(3) + 2
            else
              Nzlocal = lcgrid(3)
            endif

            call getlcrange_CompDom(Ageom,lcrange)
            !move all effective particles to their local processor 
            if(totnp.gt.1) then
              do ib = 1, ibunch
                nplctmp = Nplocal(ib)
                call ptsmv2_Ptclmger(Ebunch(ib)%Pts1,nplctmp,grid2d,Pdim,&
                                Nplcmax,lcrange)
                Nplocal(ib) = nplctmp
                ! assign new 'Nplocal' local particles on each processor.
                call setnpt_BeamBunch(Ebunch(ib),Nplocal(ib))
              enddo
            endif

            !assign the storage for potential and charge density
            call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)
            deallocate(chgdens)
            allocate(chgdens(Nxlocal,Nylocal,Nzlocal))

!-------------------------------------------------------------------
! start load balance. (at the location of new space charge calculation)
! we have used the charge density from the first bunch as the density to
! determine the balance boundary. This is not exact right for multiple bunches.
            if((mod(ibal,nbal).eq.0).and.(totnp.gt.1) ) then
              call MPI_BARRIER(comm2d,ierr)
              if(myid.eq.0) then 
                print*," load balance! "
              endif
              ! deposit particles onto grid to obtain charge density of bunch 1.
              if(flagazmuth.eq.1) then
                call chgdenstest_Depositor(Ebunch(1),chgdens,Ageom,grid2d,&
                   gammaz(1),Flagbc,Perdlen,zcent,flagpos,Ny)
              else
                call chgdens_Depositor(Ebunch(1),chgdens,Ageom,grid2d,&
                   gammaz(1),Flagbc,Perdlen,zcent,flagpos)
              endif

              call getlctabnm_CompDom(Ageom,temptab)
              lctabnmx(0:npx-1) = temptab(1,0:npx-1,0)
              lctabnmy(0:npy-1) = temptab(2,0,0:npy-1)
              call getmsize_CompDom(Ageom,msize)
              hy = msize(2)
              hz = msize(3) 
              call getrange_CompDom(Ageom,range)
              ymin = range(3)
              zmin = range(5)
              !find the new local computational domain geometry so that the
              !# of particles are roughly balanced with each local domain
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call balance_CompDom(chgdens,lctabnmx,&
                lctabrgx,npx,npy,commrow,commcol, &
                lcgrid(1),lcgrid(2),lcgrid(3),Ny,Nz,hz,zmin)
                call setlctab_CompDom(Ageom,lctabnmx,lctabrgx,&
                                     npx,npy,myidx,myidy)
              else
                call balance_CompDom(chgdens,lctabnmx,&
                lctabnmy,lctabrgx,lctabrgy,npx,npy,commrow,commcol, &
                lcgrid(1),lcgrid(2),lcgrid(3),Ny,Nz,hy,hz,ymin,zmin)
                call setlctab_CompDom(Ageom,lctabnmx,lctabnmy,lctabrgx,&
                                     lctabrgy,npx,npy,myidx,myidy)
              endif
              call getlcmnum_CompDom(Ageom,lcgrid)
              Nxlocal = lcgrid(1) 
              if(npy.gt.1) then
                Nylocal = lcgrid(2) + 2
              else
                Nylocal = lcgrid(2) 
              endif
              if(npx.gt.1) then
                Nzlocal = lcgrid(3) + 2
              else
                Nzlocal = lcgrid(3) 
              endif
              call getlcrange_CompDom(Ageom,lcrange)
              !//move particles around using the updated local geometry
              do ib = 1, ibunch
                nplctmp = Nplocal(ib) 
                ! pass particles to local space domain on new processor.
!                print*,"before particle move",myid,lcrange
                if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                  print*,"not available in this version yet....."
                  stop
                else
                  call ptsmv2_Ptclmger(Ebunch(ib)%Pts1,nplctmp,grid2d,Pdim,&
                                 Nplcmax,lcrange)
                endif
!                print*,"after particle move",myid
                ! assign new 'Nplocal' local particles on each processor.
                call setnpt_BeamBunch(Ebunch(ib),nplctmp)
                Nplocal(ib) = nplctmp
              enddo
              !prepare for charge density and Poisson solver
              deallocate(chgdens)
              allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
              call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,npy)
            endif
            ibal = ibal + 1
!end load balance.
!-------------------------------------------------------------------

            !//initialize the space-charge fields.
            deallocate(exg)
            deallocate(eyg)
            deallocate(ezg)
            deallocate(bxg)
            deallocate(byg)
            deallocate(bzg)
            allocate(exg(Nxlocal,Nylocal,Nzlocal))
            allocate(eyg(Nxlocal,Nylocal,Nzlocal))
            allocate(ezg(Nxlocal,Nylocal,Nzlocal))
            allocate(bxg(Nxlocal,Nylocal,Nzlocal))
            allocate(byg(Nxlocal,Nylocal,Nzlocal))
            allocate(bzg(Nxlocal,Nylocal,Nzlocal))
            exg = 0.0
            eyg = 0.0
            ezg = 0.0
            bxg = 0.0
            byg = 0.0
            bzg = 0.0

!!! DWA
            deallocate(dexwake)
            deallocate(deywake)
            deallocate(dezwake)
            deallocate(dbxwake)
            deallocate(dbywake)
            deallocate(dbzwake)
            allocate(dexwake(Nxlocal,Nylocal,Nzlocal))
            allocate(deywake(Nxlocal,Nylocal,Nzlocal))
            allocate(dezwake(Nxlocal,Nylocal,Nzlocal))
            allocate(dbxwake(Nxlocal,Nylocal,Nzlocal))
            allocate(dbywake(Nxlocal,Nylocal,Nzlocal))
            allocate(dbzwake(Nxlocal,Nylocal,Nzlocal))
            dexwake = 0.0
            deywake = 0.0
            dezwake = 0.0
            dbxwake = 0.0
            dbywake = 0.0
            dbzwake = 0.0
!!! end DWA

            if(npx.gt.1) then
                nzlcr = Nzlocal-2
                kadd = 1
            else
                nzlcr = Nzlocal
                kadd = 0
            endif
            if(npy.gt.1) then
                nylcr = Nylocal-2
                jadd = 1
            else
                nylcr = Nylocal
                jadd = 0
            endif

            deallocate(tmppot)
            allocate(tmppot(Nxlocal,Nylocal,Nzlocal))
            call getrange_CompDom(Ageom,range)

            !for wakefield calculation
            !call getmsize_CompDom(Ageom,msize)
            !hz = msize(3)
           
            !//sum up the space-charge fields from all bunches/bins
            do ib = 1, ibunch
              ! deposit particles onto grid to obtain charge density of each bunch/bin.
              if(flagazmuth.eq.1) then
                call chgdenstest_Depositor(Ebunch(ib),chgdens,Ageom,grid2d,&
                            gammaz(ib),Flagbc,Perdlen,zcent,flagpos,Ny)
              else
                call chgdens_Depositor(Ebunch(ib),chgdens,Ageom,grid2d,&
                                     gammaz(ib),Flagbc,Perdlen,zcent,flagpos)
              endif

              !includes wakefield
              if(flagwake.eq.1) then
                hzwake = (grange(6)-grange(5))*1.0000001/(Nz-1) !avoid over index
                xwakelc = 0.0
                ywakelc = 0.0
                denszlc = 0.0
                xwakez = 0.0
                ywakez = 0.0
                densz = 0.0
                !call getrange_CompDom(Ageom,range)
                do ipt = 1, Nplocal(ib)
                  iizz = (Ebunch(ib)%Pts1(5,ipt)-grange(5))/hzwake + 1
                  eeff = ((grange(5)-Ebunch(ib)%Pts1(5,ipt))+iizz*hzwake)/hzwake
                  iizz1 = iizz + 1
                  denszlc(iizz) = denszlc(iizz) + eeff
                  denszlc(iizz1) = denszlc(iizz1) + 1.0d0 - eeff
                  xwakelc(iizz) = xwakelc(iizz) + Ebunch(ib)%Pts1(1,ipt)*eeff
                  xwakelc(iizz1) = xwakelc(iizz1) + Ebunch(ib)%Pts1(1,ipt)*(1.-eeff)
                  ywakelc(iizz) = ywakelc(iizz) + Ebunch(ib)%Pts1(3,ipt)*eeff
                  ywakelc(iizz1) = ywakelc(iizz1) + Ebunch(ib)%Pts1(3,ipt)*(1.-eeff)
                enddo
                call MPI_ALLREDUCE(denszlc,densz,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,commcol,ierr) 
                call MPI_ALLREDUCE(xwakelc,xwakez,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,commcol,ierr) 
                call MPI_ALLREDUCE(ywakelc,ywakez,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,commcol,ierr) 

                do kz = 1, Nz
                  sendensz(kz,1) = xwakez(kz)
                enddo
                do kz = 1, Nz
                  sendensz(kz,2) = ywakez(kz)
                enddo

                Nz2 = 2*Nz
                call MPI_ALLREDUCE(sendensz,recvdensz,Nz2,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,commrow,ierr) 
 
                do kz = 1, Nz
                  xwakez(kz) = recvdensz(kz,1)
                enddo
                do kz = 1, Nz
                  ywakez(kz) = recvdensz(kz,2)
                enddo

                do kz = 1, Nz
                  sendensz(kz,1) = densz(kz)
                enddo
                call MPI_ALLREDUCE(sendensz,recvdensz,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,commrow,ierr) 

                do kz = 1, Nz
                  densz(kz) = recvdensz(kz,1)
                enddo

                !get the line charge density along z
                do kz = 1, Nz
                  recvdensz(kz,1) = densz(kz)*Ebunch(ib)%Current/Scfreq/Np(ib)/(hzwake*Scxlt)*&
                                    Ebunch(ib)%Charge/abs(Ebunch(ib)%Charge)
                enddo
                do kz = 1, Nz
                  if(densz(kz).ne.0.0) then
                    xwakez(kz) = xwakez(kz)/densz(kz)*Scxlt
                    ywakez(kz) = ywakez(kz)/densz(kz)*Scxlt
                  else
                    xwakez(kz) = 0.0
                    ywakez(kz) = 0.0
                  endif
                  recvdensz(kz,2) = ywakez(kz)
                enddo

                if(ggwk1.gt.0.0) then
                  call wakefield_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                     hzwake,aawk1,ggwk1,lengwk1,flagbtw)
                else
                  call wakereadin_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                     hzwake,lengwk1,ndatawk,wklong,wktran)
                endif

              endif

!!! DWA
              !includes dielectric wakefields
              if(dflagwake.eq.1) then
!!                hzwake = (grange(6)-grange(5))*1.0000001/(Nz-1) !avoid over index
               call getmsize_CompDom(Ageom,msize)
               hx = msize(1)*Scxlt
               hy = msize(2)*Scxlt
               hz = msize(3)*Scxlt
               xmi=range(1)*Scxlt
               ymi=range(3)*Scxlt
               call dwakefield_FieldQuant(Nx,Ny,Nz,hx,hy,hz,chgdens,dexwake,&
               deywake,dezwake,dbxwake,dbywake,dbzwake,dtype1,deps1,dnkx1,&
               dnky1,dkx,dky,amp,gammaz(ib),xmi,ymi)
              endif
!!! end DWA

              !// solve 3D Poisson's equation for each bunch/bin
              !turn off image space-charge after zimage due to the screen of pipe
              if(distance.gt.zimage) FlagImage = 0
              !zshift = -(brange(5,ib)+brange(6,ib))*Scxlt
              zshift = -(range(5)+range(6))*gammaz(ib)*Scxlt
              if((Flagbc.eq.1) .and. (flagspc.eq.1)) then
                ! solve Poisson's equation using 3D isolated boundary condition.
                ! the image space-charge potential with respect to z = 0 is calculated
                ! if the FlagImage is on
                !no space-charge, only wake field
                call update3Otnew_FieldQuant(Potential,chgdens,Ageom,&
                grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,gammaz(ib),&
                tmppot,FlagImage,zshift)
              else
                print*,"no such boundary condition type!!!"
!                stop
                Potential%FieldQ = 0.0d0
              endif

              !find the E and B fields in the lab frame from the effective bunch/bin  
              tmpflag = 0

              call gradEB_FieldQuant(Nxlocal,Nylocal,Nzlocal,&
              Potential%FieldQ,Ageom,grid2d,Flagbc,gammaz(ib),tmpflag,&
              exg,eyg,ezg,bxg,byg,bzg)

              !find the E and B fields in the lab frame from the image potential of 
              !effective bunch/bin  
              if(FlagImage.eq.1) then
                call gradEB_FieldQuant(Nxlocal,Nylocal,Nzlocal,&
                tmppot,Ageom,grid2d,Flagbc,gammaz(ib),FlagImage,&
                exg,eyg,ezg,bxg,byg,bzg)
              endif
              !add the E field from the wake to the space-charge field.
              if(flagwake.eq.1) then
                !print*,"include wake: "
                call getlctabnm_CompDom(Ageom,temptab)
                ztable(0:npx-1) = temptab(1,0:npx-1,0)
                zdisp(0) = 0
                do kz = 1, npx-1
                  zdisp(kz) = zdisp(kz-1)+ztable(kz-1)
                enddo
                do kz = 1,Nzlocal
                  kst = zdisp(myidx)+kz-kadd !get the global index
                  if(kst.eq.0) then
                    kst = 1
                  else if(kst.eq.Nz+1) then
                    kst = Nz
                  else
                  endif
                  do jy = 1, Nylocal
                    do ix = 1, Nxlocal
                      exg(ix,jy,kz) = exg(ix,jy,kz) + exwake(kst)
                      eyg(ix,jy,kz) = eyg(ix,jy,kz) + eywake(kst)
                      ezg(ix,jy,kz) = ezg(ix,jy,kz) + ezwake(kst)
                    enddo
                  enddo  
                enddo

              endif
!!! DWA
              if(dflagwake.eq.1) then
                do kz=1,Nzlocal
                 do jy=1,Nylocal
                  do ix=1,Nxlocal
                   exg(ix,jy,kz) = exg(ix,jy,kz) + dexwake(ix,jy,kz)
                   eyg(ix,jy,kz) = eyg(ix,jy,kz) + deywake(ix,jy,kz)
                   ezg(ix,jy,kz) = ezg(ix,jy,kz) + dezwake(ix,jy,kz)
                   bxg(ix,jy,kz) = bxg(ix,jy,kz) + dbxwake(ix,jy,kz)
                   byg(ix,jy,kz) = byg(ix,jy,kz) + dbywake(ix,jy,kz)
                   bzg(ix,jy,kz) = bzg(ix,jy,kz) + dbzwake(ix,jy,kz)
                  enddo
                 enddo
                enddo
 
                if(zmin.gt.dtwkstart(1).and.zmax.lt.dtwkend(1).and.&
                 iwrite.eq.1)then
                  do kz=1,Nzlocal
                    write(73,*)(kz-1)*hz,dexwake(Nx,Ny,kz),deywake(Nx,Ny&
                    ,kz),dezwake(Nx,Ny,kz),dbxwake(Nx,Ny,kz),&
                    dbywake(Nx,Ny,kz),dbzwake(Nx,Ny,kz)
                    write(72,*)(kz-1)*hz,exg(Nx,Ny,kz),eyg(Nx,Ny&
                    ,kz),ezg(Nx,Ny,kz),bxg(Nx,Ny,kz),&
                    byg(Nx,Ny,kz),bzg(Nx,Ny,kz)
                    if(kz.eq.Nzlocal)then
                     iwrite=0
                    endif
                  enddo
                endif
 
              endif
!!! end DWA

            enddo

            !interpolate the space-charge fields using CIC + external fields to
            do ib = 1, ibunch
              if(flaginteg.eq.2) then
                print*,"not available in this version yet"
                stop
              else if(flaginteg.eq.3) then
                print*,"not available in this version yet"
                stop
              else
                !Boris's 2nd order integrator
                call kick2t_BeamBunch(Nplocal(ib),Nxlocal,Nylocal,Nzlocal,&
                Ebunch(ib)%Pts1,exg,eyg,ezg,bxg,byg,bzg,Ageom,npx,npy,myidx,&
                myidy,t,Ebunch(ib)%Charge,Ebunch(ib)%Mass,dtless,Blnelem,&
                zBlnelem,idrfile,Nblem,ibinit,ibend,fldmp,Flagerr)
              endif
            enddo

            else ! point-2-point space charge force

            do ib = 1, ibunch
              totchrg = Ebunch(ib)%Current/Bfreq
              call kickpt2ptImg_BeamBunch(nplocal0,Ebunch(ib)%Pts1,t,&
                 Ebunch(ib)%Charge,&
                 Ebunch(ib)%Mass,dtless,Blnelem,zBlnelem,idrfile,Nblem,&
                 ibinit,ibend,fldmp,totchrg,r0,np0,Np(ib))
            enddo

            endif
          else !no space-charge fields, only external fields
            do ib = 1, ibunch
              call scatter20t_BeamBunch(Nplocal(ib),Ebunch(ib)%Pts1,t,&
                 Ebunch(ib)%Charge,&
                 Ebunch(ib)%Mass,dtless,Blnelem,zBlnelem,idrfile,Nblem,&
                 ibinit,ibend,fldmp,Flagerr)
            enddo
          endif

1000      continue
          !//update particle positions using new velocity for half step.
          if(flagcathode.eq.1) then
            do ib = 1, Nbunch
              call drifthalf_BeamBunch(Ebunch(ib),t,dtless,betazini)
              call driftemission_BeamBunch(Ebunch(ib),t,dtless,betazini)
            enddo
            !first order emission model with random emission time within each step
            do ib = 1, Nbunch
              do ipt = 1, Nplocal(ib)
                tmpzz = Ebunch(ib)%Pts1(5,ipt)-dtless*betazini
                recpgam = 1.0/sqrt(1.0+Ebunch(ib)%Pts1(2,ipt)**2+&
                          Ebunch(ib)%Pts1(4,ipt)**2+&
                                 Ebunch(ib)%Pts1(6,ipt)**2)
 
                if(tmpzz.le.0.0d0 .and. Ebunch(ib)%Pts1(5,ipt).ge.0.0d0) then
                  dtmp = dtless*Ebunch(ib)%Pts1(5,ipt)/(dtless*betazini)
                  Ebunch(ib)%Pts1(1,ipt) = Ebunch(ib)%Pts1(1,ipt)+dtmp*Ebunch(ib)%Pts1(2,ipt)*recpgam
                  Ebunch(ib)%Pts1(3,ipt) = Ebunch(ib)%Pts1(3,ipt)+dtmp*Ebunch(ib)%Pts1(4,ipt)*recpgam
                  Ebunch(ib)%Pts1(5,ipt) = dtmp*Ebunch(ib)%Pts1(6,ipt)*recpgam
 
                endif
              enddo
            enddo
          else
            do ib = 1, Nbunch
              call drifthalforg_BeamBunch(Ebunch(ib),t,dtless)
            enddo
          endif
          t = t + 0.5*dtless*Dt
 
          !output for every 5 steps
          if(mod(i,5).eq.0) then

          if(Flagdiag.eq.1) then
            call diagnostic1avg_Output(t,Ebunch,Nbunch)
          else if(Flagdiag.eq.2) then
            call diagnostic1avgZ_Output(t,Ebunch,Nbunch)
          else if(Flagdiag.eq.3) then
            call diagnostic1avgZtest_Output(t,Ebunch,Nbunch)
          else if(Flagdiag.eq.105) then
            nplocalmax = 0
            do ib = 1, Nbunch
              if(nplocalmax.le.Nplocal(ib)) then
                nplocalmax = Nplocal(ib)
              endif
            enddo
            allocate(tmpptsb(6,nplocalmax,Nbunch))
            do ib = 1, Nbunch
              do ipt = 1, Nplocal(ib)
                tmpptsb(:,ipt,ib) = Ebunch(ib)%Pts1(:,ipt)
              enddo
            enddo
            do ib = 1, Nbunch
              do ipt = 1, Nplocal(ib)
                pos(1)= Ebunch(ib)%Pts1(1,ipt)*Scxlt
                pos(2)= Ebunch(ib)%Pts1(3,ipt)*Scxlt
                zz = Ebunch(ib)%Pts1(5,ipt)*Scxlt
                pos(3)= zz
                pos(4) = t
                do ii = 1, Nblem
                  if( (zz.ge.zBlnelem(1,ii)) .and. &
                      (zz.le.zBlnelem(2,ii)) ) then
                    call getparam_BeamLineElem(Blnelem(ii),blength,bnseg,&
                                  bmpstp,bitype)
                    call getparam_BeamLineElem(Blnelem(ii),12,rfile)
                         
                    if((bitype.eq.105).and.(abs(rfile).gt.0.0)) then
                      call getvecAt_SolRF(pos,tmpfld,Blnelem(ii)%pslrf,&
                      fldmp(idrfile(3,ii)))
                      Ebunch(ib)%Pts1(2,ipt) = Ebunch(ib)%Pts1(2,ipt)+&
                        tmpfld(1)*Clight/Ebunch(ib)%Mass*Ebunch(ib)%Charge 
                      Ebunch(ib)%Pts1(4,ipt) = Ebunch(ib)%Pts1(4,ipt)+&
                        tmpfld(2)*Clight/Ebunch(ib)%Mass*Ebunch(ib)%Charge 
                    endif
                  endif
                enddo
              enddo
            enddo
            call diagnostic1avg_Output(t,Ebunch,Nbunch)
            do ib = 1, Nbunch
              do ipt = 1, Nplocal(ib)
                Ebunch(ib)%Pts1(:,ipt)=tmpptsb(:,ipt,ib)
              enddo
            enddo
            deallocate(tmpptsb)
          endif

          endif

          !test the rebin function
          dGspread = 1.0

          else !into bending magnet

            zorgin = zBlnelem(1,idbd) 
            zorgin2 = zBlnelem(2,idbd) 
            if(fldmp(idrfile(3,idbd))%Fcoeft(1) > 0) then
              flagcsr = 1
            else
              flagcsr = 0
            endif
            gamin = fldmp(idrfile(3,idbd))%Fcoeft(2)
           
            zz1 = fldmp(idrfile(3,idbd))%Fcoeft(21) !effective starting location of bend
            zz2 = fldmp(idrfile(3,idbd))%Fcoeft(22) !effective end location of bend

            bendlen = zz2 - zz1

            if(Nbunch>1) then
              ibb = Nbunch/2
            else
              ibb = 1
            endif
!            gamin = -Ebunch(ibb)%refptcl(6)
            vref = sqrt(1.d0-1.d0/gamin**2)
            !zz = Ebunch(ibb)%refptcl(5)*Scxlt - 0.5*vref*dtless*Scxlt
            call getparam_BeamLineElem(Blnelem(idbd),3,Bybend)
            r0bend = Ebunch(1)%Mass/Clight*sqrt(gamin**2-1.0d0)/abs(Bybend)
            kkcsr = -2.0/(3*r0bend)**(1.0/3.0)
            !here, we need to make sure that the reference particle has the right
            !designed energy so that the reference trajectory through the bend
            !is correct.
            do ib = 1, Nbunch
              call convEntr_BeamBunch(Ebunch(ib),zorgin,gamin)
            enddo
            zz = Ebunch(ibb)%refptcl(5)*Scxlt - 0.5*vref*dtless*Scxlt

            !//back half step
            do ib = 1, Nbunch
              call driftbackhalf_BeamBunch(Ebunch(ib),t,dtless)
            enddo
            t = t - 0.5*dtless*Dt

            Flagbctmp = 1
            deallocate(exg)
            deallocate(eyg)
            deallocate(ezg)
            deallocate(bxg)
            deallocate(byg)
            deallocate(bzg)
            allocate(exg(Nxlocal,Nylocal,Nzlocal))
            allocate(eyg(Nxlocal,Nylocal,Nzlocal))
            allocate(ezg(Nxlocal,Nylocal,Nzlocal))
            allocate(bxg(Nxlocal,Nylocal,Nzlocal))
            allocate(byg(Nxlocal,Nylocal,Nzlocal))
            allocate(bzg(Nxlocal,Nylocal,Nzlocal))
            do ii = 1, ntstep

              do ib = 1, Nbunch
                call drifthalfBd_BeamBunch(Ebunch(ib),t,dtless)
              enddo
              t = t + 0.5*dtless*Dt

              do ib = 1, Nbunch
                ptref = Ebunch(ib)%refptcl
                gammazavg = sqrt(1.0+ptref(2)**2+ptref(6)**2)

                !//go to the local "ptref" coordinates
                call rottoT_BeamBunch(Ebunch(ib),ptref,ptrange,poscent)

                if(Bcurr.gt.1.0e-10) then

                  call update_CompDom(Ageom,ptrange,grid2d,Flagbctmp)
                  call getlcmnum_CompDom(Ageom,lcgrid)
                  Nxlocal = lcgrid(1)
                  if(npy.gt.1) then
                    Nylocal = lcgrid(2) + 2
                  else
                    Nylocal = lcgrid(2)
                  endif
                  if(npx.gt.1) then
                    Nzlocal = lcgrid(3) + 2
                  else
                    Nzlocal = lcgrid(3)
                  endif
                  call getlcrange_CompDom(Ageom,lcrange)
                  !move all effective particles to their local processor
                  if(totnp.gt.1) then
                    nplctmp = Nplocal(ib)
                    call ptsmv2_Ptclmger(Ebunch(ib)%Pts1,nplctmp,grid2d,Pdim,&
                                Nplcmax,lcrange)
                    Nplocal(ib) = nplctmp
                    ! assign new 'Nplocal' local particles on each processor.
                    call setnpt_BeamBunch(Ebunch(ib),Nplocal(ib))
                  endif

                  !assign the storage for potential and charge density
                  call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)
                  deallocate(chgdens)
                  allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
                  ! deposit particles onto grid to obtain charge density of each bunch/bin.
                  call chgdens_Depositor(Ebunch(ib),chgdens,Ageom,grid2d,&
                                     gammazavg,Flagbc,Perdlen,zcent,flagpos)

                  if(npx.gt.1) then
                    nzlcr = Nzlocal-2
                  else
                    nzlcr = Nzlocal
                  endif
                  if(npy.gt.1) then
                    nylcr = Nylocal-2
                  else
                    nylcr = Nylocal
                  endif
                  !// solve 3D Poisson's equation for each bunch/bin
                  deallocate(tmppot)
                  allocate(tmppot(Nxlocal,Nylocal,Nzlocal))
                  zshift = 0.0
                  tmpflag = 0
                  if(Flagbc.eq.1) then
                    ! solve Poisson's equation using 3D isolated boundary condition.
                    call update3Otnew_FieldQuant(Potential,chgdens,Ageom,&
                    grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,gammazavg,&
                    tmppot,tmpflag,zshift)
                  else
                    print*,"no such boundary condition type!!!"
                    stop
                  endif

                  exg = 0.0
                  eyg = 0.0
                  ezg = 0.0
                  bxg = 0.0
                  byg = 0.0
                  bzg = 0.0
                  !find the E and B fields in the lab frame from the effective bunch/bin
                  if(flagcsr.eq.1) then
                    call getrange_CompDom(Ageom,grange)
                    hzwake = (grange(6)-grange(5))*1.0000001/(Nz-1) !avoid over index
                    denszlc = 0.0
                    do ipt = 1, Nplocal(ib)
                      iizz = (Ebunch(ib)%Pts1(5,ipt)-grange(5))/hzwake + 1
                      eeff = ((grange(5)-Ebunch(ib)%Pts1(5,ipt))+iizz*hzwake)/hzwake
                      iizz1 = iizz + 1
                      denszlc(iizz) = denszlc(iizz) + eeff
                      denszlc(iizz1) = denszlc(iizz1) + 1.0d0 - eeff
                    enddo
                    densz = 0.0
                    call MPI_ALLREDUCE(denszlc,densz,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,commcol,ierr)
                    do kz = 1, Nz
                      sendensz(kz,1) = densz(kz)
                    enddo
                    call MPI_ALLREDUCE(sendensz,recvdensz,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,commrow,ierr)
                    do kz = 1, Nz
                      densz(kz) = recvdensz(kz,1)
                    enddo
                    !get the line charge density along z
                    hzwake = hzwake*Scxlt
                    do kz = 1, Nz
                      densz(kz) = recvdensz(kz,1)*Ebunch(ib)%Current/Scfreq/Np(ib)/hzwake*&
                                    Ebunch(ib)%Charge/abs(Ebunch(ib)%Charge)
                    enddo

!no filtering is needed with IGF 
!-------------

                    zwkmin = (grange(5)-poscent)*Scxlt + (zz - zz1)

                    !This includes both transient and steady state csr wake
                    !using integrated Green's function method
                    ezwake = 0.0
                    denszp = 0.0d0
                    denszpp = 0.0d0
                    call csrwakeTrIGF_FieldQuant(Nz,r0bend,zwkmin,hzwake,&
                                   bendlen,densz,denszp,denszpp,gamin,ezwake)

                    call getlctabnm_CompDom(Ageom,temptab)
                    ztable(0:npx-1) = temptab(1,0:npx-1,0)
                    zdisp(0) = 0
                    do kz = 1, npx-1
                       zdisp(kz) = zdisp(kz-1)+ztable(kz-1)
                    enddo
                    do kz = 1,Nzlocal
                       kst = zdisp(myidx)+kz-kadd !get the global index
                       if(kst.eq.0) then
                         kst = 1
                       else if(kst.eq.Nz+1) then
                         kst = Nz
                       else
                       endif
                       do jy = 1, Nylocal
                         do ix = 1, Nxlocal
                           ezg(ix,jy,kz) = ezg(ix,jy,kz) + ezwake(kst)
                         enddo 
                       enddo
                    enddo

                  endif

                  call kick2tBd_BeamBunch(Nplocal(ib),Nxlocal,Nylocal,Nzlocal,&
                  Ebunch(ib)%Pts1,exg,eyg,ezg,bxg,byg,bzg,Ageom,npx,npy,myidx,&
                  myidy,t,Ebunch(ib)%Charge,Ebunch(ib)%Mass,dtless,Blnelem,&
                  zBlnelem,idrfile,Nblem,idbd,fldmp,Ebunch(ib)%refptcl)
                else 
                  call kick2tBd0_BeamBunch(Nplocal(ib),Ebunch(ib)%Pts1,&
                  t,Ebunch(ib)%Charge,Ebunch(ib)%Mass,dtless,Blnelem,&
                  zBlnelem,idrfile,Nblem,idbd,fldmp,Ebunch(ib)%refptcl)
                endif

                call rotbackT_BeamBunch(Ebunch(ib),ptref)

              enddo

              do ib = 1, Nbunch
                call drifthalfBd_BeamBunch(Ebunch(ib),t,dtless)
              enddo
              t = t + 0.5*dtless*Dt

              gam = sqrt(1.0+Ebunch(ibb)%refptcl(2)**2+Ebunch(ibb)%refptcl(6)**2)
              vref = sqrt((Ebunch(ibb)%refptcl(2)/gam)**2+(Ebunch(ibb)%refptcl(6)/gam)**2)
              zz = zz + vref*dtless*Scxlt

              call diagnostic1avgB_Output(t,Ebunch,Nbunch)
              print*,"zz: ",zz,zorgin2

              if(zz.gt.(zorgin2-zorgin)) then
                exit
              endif
            enddo

            zorgin = zorgin2 
            do ib =  1, Nbunch
              !call convExit_BeamBunch(Ebunch(ib))
              call convExit_BeamBunch(Ebunch(ib),zorgin2)
            enddo
          endif
          if(distance.le.tphout(iout+1) .and. &
            (distance+dzz).ge.tphout(iout+1)) then

            !convert mechanic momentum into mechanic momentum in solenoid
            if(Flagdiag.eq.105) then 

            nplocalmax = 0
            do ib = 1, Nbunch
              if(nplocalmax.le.Nplocal(ib)) then
                nplocalmax = Nplocal(ib)
              endif
            enddo
            allocate(tmpptsb(6,nplocalmax,Nbunch))
            do ib = 1, Nbunch
              do ipt = 1, Nplocal(ib)
                tmpptsb(:,ipt,ib) = Ebunch(ib)%Pts1(:,ipt)
              enddo
            enddo

            do ib = 1, Nbunch
              do ipt = 1, Nplocal(ib)
                pos(1)= Ebunch(ib)%Pts1(1,ipt)*Scxlt
                pos(2)= Ebunch(ib)%Pts1(3,ipt)*Scxlt
                zz = Ebunch(ib)%Pts1(5,ipt)*Scxlt
                pos(3)= zz
                pos(4) = t
                do ii = 1, Nblem
                  if( (zz.ge.zBlnelem(1,ii)) .and. &
                      (zz.le.zBlnelem(2,ii)) ) then
                    call getparam_BeamLineElem(Blnelem(ii),blength,bnseg,&
                                  bmpstp,bitype)
                    call getparam_BeamLineElem(Blnelem(ii),12,rfile)
                    if((bitype.eq.105).and.(abs(rfile).gt.0.0)) then
                      call getvecAt_SolRF(pos,tmpfld,Blnelem(ii)%pslrf,&
                      fldmp(idrfile(3,ii)))
                      Ebunch(ib)%Pts1(2,ipt) = Ebunch(ib)%Pts1(2,ipt)+&
                        tmpfld(1)*Clight/Ebunch(ib)%Mass*Ebunch(ib)%Charge
                      Ebunch(ib)%Pts1(4,ipt) = Ebunch(ib)%Pts1(4,ipt)+&
                        tmpfld(2)*Clight/Ebunch(ib)%Mass*Ebunch(ib)%Charge
                    endif
                  endif
                enddo
              enddo
            enddo

            iout = iout + 1
            do ib = 1, Nbunch
              call phase_Output(nfileouttmp(iout)+ib-1,Ebunch(ib),nsamp(iout))
            enddo

            do ib = 1, Nbunch
              do ipt = 1, Nplocal(ib)
                Ebunch(ib)%Pts1(:,ipt)=tmpptsb(:,ipt,ib)
              enddo
            enddo
            deallocate(tmpptsb)

            else

            iout = iout + 1
            do ib = 1, Nbunch
              call phase_Output(nfileouttmp(iout)+ib-1,Ebunch(ib),nsamp(iout))
            enddo

            endif

          endif

          if(distance.le.tslout(islout+1) .and. &
            (distance+dzz).ge.tslout(islout+1)) then
 
            islout = islout + 1
            do ib = 1, Nbunch
              qchg = Ebunch(ib)%Current/Scfreq
              call sliceprocdep_Output(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                   nslout(islout),qchg,Ebunch(ib)%Mass,nfileslout(islout)+ib-1)
            enddo
          endif

!          if(t.le.trstart .and. (t+dtless*Dt).ge.trstart) then
          if(distance.le.trstart .and. (distance+dzz).ge.trstart) then
            do ib = 1, Nbunch
              call outpoint_Output(ib*nfileout+myid,Ebunch(ib),t,&
                        i,ib,npx,npy,Ageom,iout,itsz,isteer,islout,dtless)
            enddo
          endif

          !collimation
          if(distance.le.tcol(icol+1) .and. (distance+dzz).ge.tcol(icol+1)) then
            icol = icol + 1
            do ib = 1, ibunch
              call lostXY_BeamBunch(Ebunch(ib),xradmin(icol),xradmax(icol),&
                   yradmin(icol),yradmax(icol),nplctmp,nptmp)
              Nplocal(ib) = nplctmp
              Np(ib) = nptmp
            enddo
          endif
         
          !meager multiple bins into one bin
          !if(t.le.tmger .and. (t+dtless*Dt).ge.tmger) then
          if(distance.le.tmger .and. (distance+dzz).ge.tmger) then
             nplctmp = 0
             nptottmp = 0
             tmpcur = 0.0
             do ib = 1, Nbunch
               nplctmp = nplctmp+Nplocal(ib)
               nptottmp = nptottmp+Np(ib)
               tmpcur = tmpcur + Ebunch(ib)%Current
             enddo
             allocate(tmppts(6,Nplocal(1)))
             do ipt = 1, Nplocal(1)
               tmppts(:,ipt) = Ebunch(1)%Pts1(:,ipt)
             enddo
             deallocate(Ebunch(1)%Pts1)
             allocate(Ebunch(1)%Pts1(6,nplctmp))
             do ipt = 1, Nplocal(1)
               Ebunch(1)%Pts1(:,ipt) = tmppts(:,ipt) 
             enddo
             iptnew = Nplocal(1)
             do ib = 2, Nbunch
               do ipt = 1, Nplocal(ib)
                 iptnew = iptnew + 1
                 Ebunch(1)%Pts1(:,iptnew) = Ebunch(ib)%Pts1(:,ipt) 
               enddo
             enddo 
             Nbunch = 1
             ibunch = Nbunch
             Nplocal(1) = nplctmp
             Np(1) = nptottmp
             Ebunch(1)%Npt = Np(1)
             Ebunch(1)%Nptlocal = Nplocal(1)
             Ebunch(1)%Current = tmpcur
             deallocate(tmppts)
          endif
        enddo

! final output.
        call MPI_BARRIER(comm2d,ierr)
        !drift back half time step
        do ib = 1, Nbunch   
          call drifthalf_BeamBunch(Ebunch(ib),t,-dtless,betazini)
          do ipt = 1, Nplocal(ib)
            deltaz = blnLength/Scxlt - Ebunch(ib)%Pts1(5,ipt)
            Ebunch(ib)%Pts1(1,ipt) = Ebunch(ib)%Pts1(1,ipt)+&
                     Ebunch(ib)%Pts1(2,ipt)/Ebunch(ib)%Pts1(6,ipt)*deltaz
            Ebunch(ib)%Pts1(3,ipt) = Ebunch(ib)%Pts1(3,ipt)+&
                     Ebunch(ib)%Pts1(4,ipt)/Ebunch(ib)%Pts1(6,ipt)*deltaz
          enddo
        enddo   
        !output six 2-D phase projections.
        !call phase2dold_Output(30,Ebunch,Np)
        !output all particles in 6d phase space at given location blnLength.
        do ib = 1, Nbunch
          i = 50+ib-1 
          call phase_Output(i,Ebunch(ib),1)
          qchg = Ebunch(ib)%Current/Scfreq
          call sliceprocdep_Output(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                   Nz,qchg,Ebunch(ib)%Mass,70+ib-1)
        enddo
!        call dens2d_Output(1,50,Ebunch(1),Np(1),0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,&
!                           0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0)
!        call phaserd2_Output(20,Ebunch,Np,0.0,0.0,0.0,0.0,0.0)
        t_integ = t_integ + elapsedtime_Timer(t0)
        call showtime_Timer()

        deallocate(lctabnmx,lctabnmy)
        deallocate(lctabrgx,lctabrgy)
        deallocate(temptab)
        deallocate(chgdens)
        deallocate(tmppot)
        deallocate(exg)
        deallocate(eyg)
        deallocate(ezg)
        deallocate(bxg)
        deallocate(byg)
        deallocate(bzg)
        deallocate(idrfile)
        deallocate(gammaz)
        deallocate(brange)
        if(Flagbc.eq.3) then
          deallocate(besscoef)
          deallocate(bessnorm)
          deallocate(gml)
          deallocate(pydisp)
          deallocate(modth)
        endif
        deallocate(ztable)
        deallocate(zdisp)
        deallocate(denszlc)
        deallocate(densz)
        deallocate(denszp)
        deallocate(denszpp)
        deallocate(ans)
        deallocate(csg)
        deallocate(recvdensz)
        deallocate(exwake)
        deallocate(eywake)
        deallocate(ezwake)
        deallocate(sendensz)
        deallocate(xwakelc)
        deallocate(xwakez)
        deallocate(ywakelc)
        deallocate(ywakez)
        deallocate(wklong2d)
        deallocate(wktran2d)
        deallocate(wktran)
        deallocate(wklong)

!!! DWA
        deallocate(dexwake)
        deallocate(deywake)
        deallocate(dezwake)
        deallocate(dbxwake)
        deallocate(dbywake)
        deallocate(dbzwake)
!!! end DWA

        end subroutine run_AccSimulator

        subroutine rebin_Utility(this,Nbunch,ibunch,dGspread)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: dGspread
        type (BeamBunch), dimension(:), intent(inout) :: this
        integer, intent(in) :: Nbunch
        integer, intent(out) :: ibunch
        double precision, allocatable, dimension(:,:,:) :: tmppts
        double precision :: gambetminlc,gambetmaxlc,gambetmin,&
               gambetmax,gami,del1,deltaG,del1rel,hgam,gammamin,&
               gammamax
        integer, allocatable, dimension(:) :: nbinlc,nbingl
        integer :: ibin,i,j,k,ib,innp,ierr
     
        gambetminlc = 1.0e10
        gambetmaxlc = -1.0e10
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          do i = 1, innp
            if(gambetminlc.gt.this(ib)%Pts1(6,i)) then
              gambetminlc = this(ib)%Pts1(6,i)
            endif
            if(gambetmaxlc.lt.this(ib)%Pts1(6,i)) then
              gambetmaxlc = this(ib)%Pts1(6,i)
            endif
          enddo
        enddo 
        call MPI_ALLREDUCE(gambetminlc,gambetmin,1,MPI_DOUBLE_PRECISION,&
                        MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(gambetmaxlc,gambetmax,1,MPI_DOUBLE_PRECISION,&
                        MPI_MAX,MPI_COMM_WORLD,ierr)
      
        !spread associated with the gamma betaz
        gammamin = sqrt(1+gambetmin**2)
        gammamax = sqrt(1+gambetmax**2)
        deltaG = gammamax - gammamin
        print*,"gammamin: ",gammamin,gammamax,deltaG
        !determine the # of bins for the given relative energy spread
        do i = 1, Nbunchmax 
          ibunch = i
          del1 = deltaG/ibunch
          del1rel = del1/(gammamin+0.5*del1)
          if(del1rel.lt.dGspread) then
            exit
          endif 
        enddo

        hgam = deltaG/ibunch

        allocate(nbinlc(ibunch))
        nbinlc = 0
        !find the # of local particles for each bin
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          do i = 1, innp
            gami = sqrt(1.0+this(ib)%Pts1(6,i)**2)
            ibin = (gami-gammamin)/hgam + 1
            if(ibin.ge.ibunch) ibin = ibunch
            nbinlc(ibin) = nbinlc(ibin) + 1
          enddo
        enddo

        !//copy the particles to the new bin stored by a temporary array
        do i = 1, ibunch
          allocate(tmppts(6,nbinlc(i),i))
        enddo
        nbinlc = 0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          do i = 1, innp
            gami = sqrt(1.0+this(ib)%Pts1(6,i)**2)
            ibin = (gami-gammamin)/hgam + 1
            if(ibin.ge.ibunch) ibin = ibunch
            nbinlc(ibin) = nbinlc(ibin) + 1
            do j = 1, 6
              tmppts(j,nbinlc(ibin),ibin) = this(ib)%Pts1(j,i)
            enddo
          enddo
        enddo

        !release the memory held by Ebunch
        do i = 1, Nbunch
          deallocate(Ebunch(i)%Pts1)
          allocate(Ebunch(i)%Pts1(6,1))
        enddo

        !//copy the particles back to Ebunch with new bin
        do i = 1, ibunch
          allocate(Ebunch(i)%Pts1(6,nbinlc(i)))
          do j = 1, nbinlc(i)
            do k = 1, 6
              Ebunch(i)%Pts1(k,j) = tmppts(k,j,i)
            enddo 
          enddo
          Ebunch(i)%Nptlocal = nbinlc(i)
        enddo 

        deallocate(tmppts)

        allocate(nbingl(ibunch))
        call MPI_ALLREDUCE(nbinlc,nbingl,ibunch,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)

        do i = 1, ibunch
          Ebunch(i)%Npt = nbingl(i)
        enddo

        deallocate(nbinlc)
        deallocate(nbingl)

        if(ibunch.gt.Nbunch) then
          do i = Nbunch + 1, ibunch
            Ebunch(i)%Current = Ebunch(1)%Current 
            Ebunch(i)%Mass = Ebunch(1)%Mass
            Ebunch(i)%Charge = Ebunch(1)%Charge
            Ebunch(i)%refptcl = Ebunch(1)%refptcl
          enddo
        endif

        end subroutine rebin_Utility

        subroutine destruct_AccSimulator(time)
        implicit none
        include 'mpif.h'
        double precision :: time
        integer :: ib
 
!        print*,"before  destruct1: "
        do ib = 1, Maxoverlap
!          print*,"inside: ",ib
!          call destructt_Data(fldmp(ib))
        enddo
        print*,"before  destruct2: "
        do ib = 1, Nbunch
          call destruct_BeamBunch(Ebunch(ib))
        enddo
        print*,"before  destruct3: "
        call destruct_FieldQuant(Potential)
        print*,"before  destruct4: "
        call destruct_CompDom(Ageom)

        call end_Output(time)

        end subroutine destruct_AccSimulator

      end module AccSimulatorclass
