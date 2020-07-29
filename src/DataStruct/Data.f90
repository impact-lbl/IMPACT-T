!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! Dataclass: Field data class in DATA STRUCTURE layer.
! Version: 1.0
! Author: Ji Qiang 
! Description: This class stores the rf cavity data Ez, Ez', Ez'' on the 
!              axis; Fourier coefficients of Ez on the axis; Ez(r,z),
!              Er(r,z), Htheta(r,z) on the r-z grid plane; and Ex(x,y,z),
!              Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z), Bz(x,y,z) on
!              uniform x, y, z grid, and Br(r,z) and Bz(r,z) on the r-z grid.
! Comments: 
!----------------------------------------------------------------
      module Dataclass
        use NumConstclass
        use mpistub
!        save
!-----------------------------------------------------------------------
! using the x-y-z field data (Ex,Ey,Ez,Bx,By,Bz) directly.
        !number of grid points along x, y, and z direction.
        integer :: NxIntvRfg = 1
        integer :: NyIntvRfg = 1
        integer :: NzIntvRfg = 1
        !range in x, y, and zdirections.
        double precision :: XmaxRfg,XminRfg,YmaxRfg,YminRfg,ZmaxRfg,ZminRfg
        ! discrete Ex(x,y,z), Ey(x,y,z), Ez(x,y,z) and Bx(x,y,z), By(x,y,z), and 
        ! Bz(x,y,z) rf data. Here, the grid is uniform in x, y and z.
        double precision,allocatable,dimension(:,:,:) :: &
               Exgrid,Eygrid,Ezgrid,Bxgrid,Bygrid,Bzgrid
!-----------------------------------------------------------------------
! using the r-z field data (Er,Ez,Htheta) directly.
        !number of grid points along r direction.
        integer :: NrIntvRf = 1
        !number of grid points along z direction.
        integer :: NzIntvRf = 1
        !range in r and z directions.
        double precision :: RmaxRf,RminRf,ZmaxRf,ZminRf
        ! discrete Ez(r,z), Er(r,z) and Htheta(r,z) rf data. Here, the grid
        ! is uniform in r and z.
        double precision,allocatable,dimension(:,:) :: &
               ezdata,erdata,btdata
!-----------------------------------------------------------------------
! using only on-axis field data and its derivities.
        !initial number of grid points on the axis.
        integer, parameter :: Ndataini = 5000
        ! discrete Ez(0,z), Ez'(0,z), Ez''(0,z) rf data.
        double precision,dimension(Ndataini) :: zdat,edat,epdat,eppdat
!----------------------------------------------------------------------
! using the Fourier coefficients
        !number of Fourier expansion coefficients.
        integer, parameter :: NcoefF = 401
        double precision,dimension(NcoefF) :: Fcoef
        !Fcoef(1): constant
        !Fcoef(2k): cosine term
        !Fcoef(2k+1): sine term
!----------------------------------------------------------------------
        ! practical number of grid data on the axis or Fourier coefficients.
        integer :: Ndata
!------------------------------------------------------------------------
!
!       data storage for the t code
        type fielddata
          ! using the x-y-z field data (Ex,Ey,Ez,Bx,By,Bz) 
          !directly for beam line element i.
          !number of grid points along x, y, and z direction.
          integer :: NxIntvRfgt 
          integer :: NyIntvRfgt
          integer :: NzIntvRfgt
          !range in x, y, and zdirections.
          double precision  :: XmaxRfgt,XminRfgt,&
               YmaxRfgt,YminRfgt,ZmaxRfgt,ZminRfgt
          ! discrete Ex(x,y,z,i), Ey(x,y,z,i), Ez(x,y,z,i) and Bx(x,y,z,i), 
          ! By(x,y,z,i), and
          ! Bz(x,y,z,i) rf data. Here, the grid is uniform in x, y and z.
          double complex,pointer,dimension(:,:,:) :: &
               Exgridt,Eygridt,Ezgridt,Bxgridt,Bygridt,Bzgridt
          ! using the r-z field data (Er,Ez,Htheta) directly.
          !number of grid points along r direction.
          integer :: NrIntvRft
          !number of grid points along z direction.
          integer :: NzIntvRft 
          !range in r and z directions.
          double precision :: RmaxRft,RminRft,&
              ZmaxRft,ZminRft
          ! discrete Ezt(r,z,i), Ert(r,z,i) and 
          ! Btheta(r,z,i) rf data. Here, the grid
          ! is uniform in r and z.
          double precision,pointer,dimension(:,:) :: &
               ezdatat,erdatat,btdatat,brdatat,bzdatat
          ! using the analytical function for field data 
          ! (Ex,Ey,Ez,Bx,By,Bz) for beam line element i directly.
          !number of coefficients in Ex.
          integer :: Ncex 
          !number of coefficients in Ey.
          integer :: Ncey 
          !number of coefficients in Ez.
          integer :: Ncez 
          !number of coefficients in Bx.
          integer :: Ncbx
          !number of coefficients in By.
          integer :: Ncby
          !number of coefficients in Bz.
          integer :: Ncbz 
          !coefficients for analytical function description of the rf data.
          double precision,pointer,dimension(:) :: &
              coefex,coefey,coefez,coefbx,coefby,coefbz
          ! using the Fourier coefficients
          !integer, parameter :: NcoefFt = 401
          ! practical number of grid data on the axis or Fourier coefficients.
          double precision,dimension(401) :: Fcoeft
          !store the discrete data of ez,ez',ez'',ez''' on axis. 
          !here, maximum data point is 5000 for E field or B field.
          double precision,dimension(4,10002) :: Fcoeftdata
          integer :: Ndatat
        end type fielddata
      contains

!-------------------------------------------------------------------------
! the following are used to store the field data in position z domain, i.e.
! used by Impact-Z code.
        !Initialize the data storage arrays.
        subroutine init_Data()
        implicit none
        include 'mpif.h' 
        integer :: i,j

        NzIntvRf = 1
        NrIntvRf = 1
        allocate(ezdata(NzIntvRf+1,NrIntvRf+1))
        allocate(erdata(NzIntvRf+1,NrIntvRf+1))
        allocate(btdata(NzIntvRf+1,NrIntvRf+1))

        do j = 1, NrIntvRf+1
          do i = 1, NzIntvRf+1
            ezdata(i,j) = 0.0
            erdata(i,j) = 0.0
            btdata(i,j) = 0.0
          enddo
        enddo

        do i = 1, Ndataini
          zdat(i) = 0.0
          edat(i) = 0.0
          epdat(i) = 0.0
          eppdat(i) = 0.0
        enddo

        Ndata = 1
        RminRf = 0.0
        RmaxRf = 1.0
        ZminRf = 0.0
        ZmaxRf = 1.0

!initialization of Ex,Ey,Ez,Bx,By,Bz
        NxIntvRfg = 1
        NyIntvRfg = 1
        NzIntvRfg = 1
        allocate(Exgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Eygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Ezgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bxgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bzgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        Exgrid = 0.0
        Eygrid = 0.0
        Ezgrid = 0.0
        Bxgrid = 0.0
        Bygrid = 0.0
        Bzgrid = 0.0
        XminRfg = 0.0
        XmaxRfg = 1.0
        YminRfg = 0.0
        YmaxRfg = 1.0
        ZminRfg = 0.0
        ZmaxRfg = 1.0

        end subroutine init_Data

        subroutine destruct_Data()
        implicit none
        include 'mpif.h' 

        deallocate(ezdata)
        deallocate(erdata)
        deallocate(btdata)
        deallocate(Exgrid)
        deallocate(Eygrid)
        deallocate(Ezgrid)
        deallocate(Bxgrid)
        deallocate(Bygrid)
        deallocate(Bzgrid)
         
        end subroutine destruct_Data

        !read in the discrete field data Ez(0,z), Ez'(0,z), Ez''(0,z) 
        !distribution along axis zdat from files "rfdatax or rfdataxx 
        !or rfdataxxx".
        subroutine read1_Data(ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1
        character*7 name1
        character*8 name2
        character*9 name3

        name1 = 'rfdatax'
        name2 = 'rfdataxx'
        name3 = 'rfdataxxx'

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
!            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          n = 0
50        continue
            read(14,*,end=77)tmp1,tmp2,tmp3,tmp4
            n = n + 1
            zdat(n) = tmp1
            edat(n) = tmp2
            epdat(n) = tmp3
            eppdat(n) = tmp4
            !write(15,100)zdat(n),edat(n),epdat(n),eppdat(n)
            !divided by 100 is due to the unit in rf data is cm.
          goto 50
77        continue
          close(14)
!          close(15)
          Ndata = n
          zdat1 = zdat(1)
          do i = 1, Ndata
            zdat(i) = zdat(i) - zdat1
          enddo
        endif
100     format(6(1x,1pe15.8))

        call MPI_BCAST(Ndata,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(zdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(edat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(epdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(eppdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndata
        end subroutine read1_Data

        ! read in discrete Ez(r,z), Er(r,z) and Btheta(r,z) rf data from
        ! files "1Tx.T7 or 1Txx.T7 or 1Txxx.T7. Here, the grid
        ! is uniform in r and z.
        subroutine read2_Data(ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0
        character*6 name1
        character*7 name2
        character*8 name3

        name1 = '1Tx.T7'
        name2 = '1Txx.T7'
        name3 = '1Txxx.T7'
   
        mu0 = 4*2*asin(1.0)*1.0e-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(3:3) = char(ifile+48)
            open(14,file=name1,status='old')
            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(3:3) = char(ii+48)
            name2(4:4) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(3:3) = char(ii+48)
            name3(4:4) = char(kk+48)
            name3(5:5) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          ! the input range units are cm
          read(14,*,end=33)tmp1,tmp2,tmpint
          ZminRf = tmp1/100.0
          ZmaxRf = tmp2/100.0
          NzIntvRf = tmpint
          if(tmpint.ne.NzIntvRf) then
            print*,"input data wrong in Z: ",NzIntvRf,tmpint
            stop
          endif
          ! the input range units are cm
          read(14,*,end=33)tmp1
          read(14,*,end=33)tmp1,tmp2,tmpint
          RminRf = tmp1/100.0
          RmaxRf = tmp2/100.0
          NrIntvRf = tmpint
          if(tmpint.ne.NrIntvRf) then
            print*,"input data wrong in R: ",NrIntvRf,tmpint
            stop
          endif
          !print*,"Nz: ",NzIntvRf,ZminRf,ZmaxRf
          !print*,"Nr: ",NrIntvRf,RminRf,RmaxRf
        endif
33      continue
        call MPI_BCAST(NrIntvRf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NzIntvRf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        deallocate(ezdata)
        deallocate(erdata)
        deallocate(btdata)
        allocate(ezdata(NzIntvRf+1,NrIntvRf+1))
        allocate(erdata(NzIntvRf+1,NrIntvRf+1))
        allocate(btdata(NzIntvRf+1,NrIntvRf+1))

        if(myrank.eq.0) then
          n = 0
50        continue
            if(mod(n,2).eq.0) then
              read(14,*,end=77)tmp1,tmp2,tmp3
              nn = n/2+1
              j  = (nn-1)/(NzIntvRf+1) + 1
              i = mod((nn-1),NzIntvRf+1) + 1
              ezdata(i,j) = tmp1
              erdata(i,j) = tmp2
              n = n + 1
              write(15,100)float(i-1),ezdata(i,j)
            else
              read(14,*,end=77)tmp1
              nn = (n+1)/2
              j  = (nn-1)/(NzIntvRf+1) + 1
              i = mod((nn-1),NzIntvRf+1) + 1
              btdata(i,j) = tmp1
              n = n + 1
            endif
          goto 50
77        continue
          close(14)
          close(15)
          Ndatalc = n/2
          !print*,"Ndata in 0: ",Ndatalc
        endif
100     format(2(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZminRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZmaxRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(RminRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(RmaxRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ezdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(erdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(btdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndatalc
        end subroutine read2_Data

        !readin the Fourier coefficients of the RF field from files 
        !"rfdatax or rfdataxx or rfdataxxx".
        subroutine read3_Data(ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1
        character*7 name1
        character*8 name2
        character*9 name3

        name1 = 'rfdatax'
        name2 = 'rfdataxx'
        name3 = 'rfdataxxx'

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          n = 0
50        continue
            read(14,*,end=77)tmp1
            n = n + 1
            Fcoef(n) = tmp1
          goto 50
77        continue
          close(14)
          Ndata = n
        endif

        call MPI_BCAST(Ndata,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Fcoef,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        !enforce the total # of coefficients to be an odd number.
        if(mod(Ndata,2).eq.0) then
          Ndata = Ndata + 1
          Fcoef(Ndata) = 0.0
        endif

        !print*,"Ndata: ",Ndata
        end subroutine read3_Data

        ! read in discrete Ex(x,y,z), Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z),
        ! Bz(x,y,z) rf data from
        ! files "1Tx.T7 or 1Txx.T7 or 1Txxx.T7. Here, the grid
        ! is uniform in x, y and z.
        subroutine read4_Data(ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint,k
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0,tmp5,tmp6
        character*6 name1
        character*7 name2
        character*8 name3

        name1 = '1Tx.T7'
        name2 = '1Txx.T7'
        name3 = '1Txxx.T7'
   
        mu0 = 4*2*asin(1.0)*1.0e-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(3:3) = char(ifile+48)
            open(14,file=name1,status='old')
            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(3:3) = char(ii+48)
            name2(4:4) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(3:3) = char(ii+48)
            name3(4:4) = char(kk+48)
            name3(5:5) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          ! the input range units are m
          read(14,*,end=33)tmp1,tmp2,tmpint
          XminRfg = tmp1
          XmaxRfg = tmp2
          NxIntvRfg = tmpint
          read(14,*,end=33)tmp1,tmp2,tmpint
          YminRfg = tmp1
          YmaxRfg = tmp2
          NyIntvRfg = tmpint
          read(14,*,end=33)tmp1,tmp2,tmpint
          ZminRfg = tmp1
          ZmaxRfg = tmp2
          NzIntvRfg = tmpint
          !print*,"Nz: ",NzIntvRf,ZminRf,ZmaxRf
          !print*,"Nr: ",NrIntvRf,RminRf,RmaxRf
        endif
33      continue
        call MPI_BCAST(NxIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NyIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NzIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(XminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(XmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(YminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(YmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        deallocate(Exgrid)
        deallocate(Eygrid)
        deallocate(Ezgrid)
        deallocate(Bxgrid)
        deallocate(Bygrid)
        deallocate(Bzgrid)
        allocate(Exgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Eygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Ezgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bxgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bzgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        Exgrid = 0.0
        Eygrid = 0.0
        Ezgrid = 0.0
        Bxgrid = 0.0
        Bygrid = 0.0
        Bzgrid = 0.0

        if(myrank.eq.0) then
          n = 0
50        continue
              read(14,*,end=77)tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
              n = n+1
              k = (n-1)/((NxIntvRfg+1)*(NyIntvRfg+1))+1
              j = (n-1-(k-1)*(NxIntvRfg+1)*(NyIntvRfg+1))/(NxIntvRfg+1) + 1
              i = n - (k-1)*(NxIntvRfg+1)*(NyIntvRfg+1) - (j-1)*(NxIntvRfg+1)
              Exgrid(i,j,k) = tmp1
              Eygrid(i,j,k) = tmp2
              Ezgrid(i,j,k) = tmp3
              Bxgrid(i,j,k) = tmp4
              Bygrid(i,j,k) = tmp5
              Bzgrid(i,j,k) = tmp6
          goto 50
77        continue
          close(14)
          Ndatalc = n
          !print*,"Ndata in 0: ",Ndatalc
        endif
100     format(2(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Exgrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Eygrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Ezgrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bxgrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bygrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bzgrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndata
        end subroutine read4_Data

!-------------------------------------------------------------------------
! the following are used for storing the field data in time t domain 
! used by Impact-T code.
        subroutine initt_Data(this)
        implicit none
        include 'mpif.h' 
        type (fielddata), intent(inout) :: this
        integer :: i,j

!initialization of Er, Ez, Btheta
        this%NzIntvRft = 1
        this%NrIntvRft = 1
        allocate(this%ezdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%erdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%btdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%brdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%bzdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        this%ezdatat = 0.0
        this%erdatat = 0.0
        this%btdatat = 0.0
        this%brdatat = 0.0
        this%bzdatat = 0.0

        this%Ndatat = 1
        this%RminRft = 0.0
        this%RmaxRft = 1.0
        this%ZminRft = 0.0
        this%ZmaxRft = 1.0

!initialization of Ex,Ey,Ez,Bx,By,Bz
        this%NxIntvRfgt = 1
        this%NyIntvRfgt = 1
        this%NzIntvRfgt = 1
        allocate(this%Exgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Eygridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Ezgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bxgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bygridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bzgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        this%Exgridt = 0.0
        this%Eygridt = 0.0
        this%Ezgridt = 0.0
        this%Bxgridt = 0.0
        this%Bygridt = 0.0
        this%Bzgridt = 0.0

        this%XminRfgt = 0.0
        this%XmaxRfgt = 1.0
        this%YminRfgt = 0.0
        this%YmaxRfgt = 1.0
        this%ZminRfgt = 0.0
        this%ZmaxRfgt = 1.0
! initialization of the coefficients for analytical discription of Ex,Ey,Ez,Bx,By,Bz.
        this%Ncex = 1
        this%Ncey = 1
        this%Ncez = 1
        this%Ncbx = 1
        this%Ncby = 1
        this%Ncbz = 1
        allocate(this%coefex(this%Ncex))
        allocate(this%coefey(this%Ncey))
        allocate(this%coefez(this%Ncez))
        allocate(this%coefbx(this%Ncbx))
        allocate(this%coefby(this%Ncby))
        allocate(this%coefbz(this%Ncbz))
        this%coefex = 0.0
        this%coefey = 0.0
        this%coefez = 0.0
        this%coefbx = 0.0
        this%coefby = 0.0
        this%coefbz = 0.0

        end subroutine initt_Data

        subroutine destructt_Data(this)
        implicit none
        include 'mpif.h' 
        type (fielddata), intent(inout) :: this

        print*,"d1: "
        deallocate(this%ezdatat)
        deallocate(this%erdatat)
        deallocate(this%btdatat)
        deallocate(this%brdatat)
        deallocate(this%bzdatat)
        print*,"d2: "
        deallocate(this%Exgridt)
        print*,"d20: "
        deallocate(this%Eygridt)
        print*,"d21: "
        deallocate(this%Ezgridt)
        print*,"d22: "
        deallocate(this%Bxgridt)
        print*,"d23: "
        deallocate(this%Bygridt)
        print*,"d24: "
        deallocate(this%Bzgridt)
        print*,"d3: "
        deallocate(this%coefex)
        deallocate(this%coefey)
        deallocate(this%coefez)
        deallocate(this%coefbx)
        deallocate(this%coefby)
        deallocate(this%coefbz)
        print*,"d4: "
         
        end subroutine destructt_Data

        !readin the Fourier coefficients of the RF field from files 
        !"rfdatax or rfdataxx or rfdataxxx".
        subroutine read1t_Data(this,ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        type (fielddata), intent(inout) :: this
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1
        character*7 name1
        character*8 name2
        character*9 name3

        name1 = 'rfdatax'
        name2 = 'rfdataxx'
        name3 = 'rfdataxxx'

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          print*,"name: ",name1,name2,name3
          n = 0
50        continue
            read(14,*,end=77)tmp1
            n = n + 1
            this%Fcoeft(n) = tmp1
          goto 50
77        continue
          close(14)
          this%Ndatat = n
        endif

        call MPI_BCAST(this%Ndatat,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%Fcoeft(1),this%Ndatat,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        !enforce the total # of coefficients to be an odd number.
        if(mod(this%Ndatat,2).eq.0) then
          this%Ndatat = this%Ndatat + 1
          this%Fcoeft(this%Ndatat) = 0.0
        endif

        print*,"Ndata: ",this%Ndatat
        end subroutine read1t_Data

        !readin the Fourier coefficients of the RF field from files 
        !"rfdatax or rfdataxx or rfdataxxx".
        subroutine read1tdata_Data(this,ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        type (fielddata), intent(inout) :: this
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,i1,ntmp
        double precision :: tmp1,tmp2,tmp3
        character*7 name1
        character*8 name2
        character*9 name3

        name1 = 'rfdatax'
        name2 = 'rfdataxx'
        name3 = 'rfdataxxx'

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          print*,"name: ",name1,name2,name3

          !initialization
          this%Fcoeftdata = 0.0d0
          ntmp = 1

          !read in Ez RF data
          !read in # of data, starting location of the field, ending location 
          read(14,*)tmp1,tmp2,tmp3
          this%Fcoeftdata(1,1) = tmp1
          this%Fcoeftdata(2,1) = tmp2
          this%Fcoeftdata(3,1) = tmp3
          ntmp = tmp1 + 0.1
          do i = 2,ntmp+1
            read(14,*)this%Fcoeftdata(1,i),this%Fcoeftdata(2,i),&
                      this%Fcoeftdata(3,i),this%Fcoeftdata(4,i)
          enddo
          !read in static Bz data
          !read in # of data, starting location of the field, ending location 
          read(14,*)tmp1,tmp2,tmp3
          i1 = ntmp + 2
          this%Fcoeftdata(1,i1) = tmp1
          this%Fcoeftdata(2,i1) = tmp2
          this%Fcoeftdata(3,i1) = tmp3
          ntmp = tmp1 + 0.1
          do i = i1+1,i1+ntmp
            read(14,*)this%Fcoeftdata(1,i),this%Fcoeftdata(2,i),&
                      this%Fcoeftdata(3,i),this%Fcoeftdata(4,i)
          enddo

          i1 = i1 + ntmp
          close(14)

          this%Ndatat = i1
          ntmp = i1*4
        endif

        call MPI_BCAST(this%Ndatat,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        ntmp = this%Ndatat*4
        call MPI_BCAST(this%Fcoeftdata(1,1),ntmp,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        print*,"Ndata: ",this%Ndatat
        end subroutine read1tdata_Data

        !J.Q. changed the format of inputs (7/28/2020) to the new
        !Superfilsh RF cavity output.
        ! read in discrete Ez(r,z), Er(r,z) and Btheta(r,z) rf data from
        ! files "1Tx.T7 or 1Txx.T7 or 1Txxx.T7. Here, the grid
        ! is uniform in r and z.
        subroutine read2t_Data(this,ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        type (fielddata), intent(inout) :: this
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0
        character*6 name1
        character*7 name2
        character*8 name3
        double precision :: tmpmax

        name1 = '1Tx.T7'
        name2 = '1Txx.T7'
        name3 = '1Txxx.T7'
   
        mu0 = 4*2*asin(1.0d0)*1.0d-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(3:3) = char(ifile+48)
            open(14,file=name1,status='old')
            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(3:3) = char(ii+48)
            name2(4:4) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(3:3) = char(ii+48)
            name3(4:4) = char(kk+48)
            name3(5:5) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          ! the input range units are cm
          read(14,*,end=33)tmp1,tmp2,tmpint
          this%ZminRft = tmp1/100.0
          this%ZmaxRft = tmp2/100.0
          this%NzIntvRft = tmpint
          if(tmpint.ne.this%NzIntvRft) then
            print*,"input data wrong in Z: ",this%NzIntvRft,tmpint
            stop
          endif
          ! the input range units are cm
          read(14,*,end=33)tmp1
          read(14,*,end=33)tmp1,tmp2,tmpint
          this%RminRft = tmp1/100.0
          this%RmaxRft = tmp2/100.0
          this%NrIntvRft = tmpint
          if(tmpint.ne.this%NrIntvRft) then
            print*,"input data wrong in R: ",this%NrIntvRft,tmpint
            stop
          endif
          !print*,"Nz: ",NzIntvRft(ib),ZminRft(ib),ZmaxRft(ib)
          !print*,"Nr: ",NrIntvRft(ib),RminRft(ib),RmaxRft(ib)
        endif
33      continue
        call MPI_BCAST(this%NrIntvRft,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%NzIntvRft,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        deallocate(this%ezdatat)
        deallocate(this%erdatat)
        deallocate(this%btdatat)
        allocate(this%ezdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%erdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%btdatat(this%NzIntvRft+1,this%NrIntvRft+1))

        if(myrank.eq.0) then
          tmpmax = -1.0e10
          n = 1
50        continue
              read(14,*,end=77)tmp1,tmp2,tmp3,tmp4
              j  = (n-1)/(this%NzIntvRft+1) + 1
              i = mod((n-1),this%NzIntvRft+1) + 1
              this%ezdatat(i,j) = tmp1*1.0e6
              this%erdatat(i,j) = tmp2*1.0e6
              !convert from H (A/m) to Tesla
              this%btdatat(i,j) = tmp4*mu0
              n = n + 1
              write(15,100)float(i-1),this%ezdatat(i,j),this%erdatat(i,j)
          goto 50
77        continue
          close(14)
          close(15)
          Ndatalc = n - 1
          !print*,"Ndata in 0: ",Ndatalc
          do i = 1, this%NzIntvRft+1
            !find the max. of e field on axis
            tmp3 = sqrt(this%ezdatat(i,1)**2+this%erdatat(i,1)**2) 
            if(tmpmax.lt.abs(tmp3)) tmpmax = abs(tmp3)
          enddo
          print*,"maximum E field on axis: ",tmpmax
        endif
100     format(3(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ZminRft,1,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ZmaxRft,1,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%RminRft,1,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%RmaxRft,1,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ezdatat(1,1),Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%erdatat(1,1),Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%btdatat(1,1),Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndatalc
        end subroutine read2t_Data

        ! read in discrete Ez(r,z), Er(r,z) and Btheta(r,z) rf data from
        ! files "1Tx.T7 or 1Txx.T7 or 1Txxx.T7. Here, the grid
        ! is uniform in r and z.
        subroutine read2told_Data(this,ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        type (fielddata), intent(inout) :: this
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0
        character*6 name1
        character*7 name2
        character*8 name3
        double precision :: tmpmax

        name1 = '1Tx.T7'
        name2 = '1Txx.T7'
        name3 = '1Txxx.T7'
   
        mu0 = 4*2*asin(1.0)*1.0e-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(3:3) = char(ifile+48)
            open(14,file=name1,status='old')
            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(3:3) = char(ii+48)
            name2(4:4) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(3:3) = char(ii+48)
            name3(4:4) = char(kk+48)
            name3(5:5) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          ! the input range units are cm
          read(14,*,end=33)tmp1,tmp2,tmpint
          this%ZminRft = tmp1/100.0
          this%ZmaxRft = tmp2/100.0
          this%NzIntvRft = tmpint
          if(tmpint.ne.this%NzIntvRft) then
            print*,"input data wrong in Z: ",this%NzIntvRft,tmpint
            stop
          endif
          ! the input range units are cm
          read(14,*,end=33)tmp1
          read(14,*,end=33)tmp1,tmp2,tmpint
          this%RminRft = tmp1/100.0
          this%RmaxRft = tmp2/100.0
          this%NrIntvRft = tmpint
          if(tmpint.ne.this%NrIntvRft) then
            print*,"input data wrong in R: ",this%NrIntvRft,tmpint
            stop
          endif
          !print*,"Nz: ",NzIntvRft(ib),ZminRft(ib),ZmaxRft(ib)
          !print*,"Nr: ",NrIntvRft(ib),RminRft(ib),RmaxRft(ib)
        endif
33      continue
        call MPI_BCAST(this%NrIntvRft,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%NzIntvRft,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        deallocate(this%ezdatat)
        deallocate(this%erdatat)
        deallocate(this%btdatat)
        allocate(this%ezdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%erdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%btdatat(this%NzIntvRft+1,this%NrIntvRft+1))

        if(myrank.eq.0) then
          tmpmax = -1.0e10
          n = 0
50        continue
            if(mod(n,2).eq.0) then
              read(14,*,end=77)tmp1,tmp2,tmp3
              nn = n/2+1
              j  = (nn-1)/(this%NzIntvRft+1) + 1
              i = mod((nn-1),this%NzIntvRft+1) + 1
              this%ezdatat(i,j) = tmp1*1.0e6
              this%erdatat(i,j) = tmp2*1.0e6
              n = n + 1
              write(15,100)float(i-1),this%ezdatat(i,j),this%erdatat(i,j)
            else
              read(14,*,end=77)tmp1
              nn = (n+1)/2
              j  = (nn-1)/(this%NzIntvRft+1) + 1
              i = mod((nn-1),this%NzIntvRft+1) + 1
              !convert from H (A/m) to Tesla
              this%btdatat(i,j) = tmp1*mu0
              n = n + 1
            endif
          goto 50
77        continue
          close(14)
          close(15)
          Ndatalc = n/2
          !print*,"Ndata in 0: ",Ndatalc
          do i = 1, this%NzIntvRft+1
            !find the max. of e field on axis
            tmp3 = sqrt(this%ezdatat(i,1)**2+this%erdatat(i,1)**2) 
            if(tmpmax.lt.abs(tmp3)) tmpmax = abs(tmp3)
          enddo
          print*,"maximum E field on axis: ",tmpmax
        endif
100     format(3(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ZminRft,1,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ZmaxRft,1,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%RminRft,1,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%RmaxRft,1,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ezdatat(1,1),Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%erdatat(1,1),Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%btdatat(1,1),Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndatalc
        end subroutine read2told_Data

        ! read in analytical function coefficients of 
        ! Ex(x,y,z), Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z),
        ! Bz(x,y,z) rf data from
        !"rfdatax or rfdataxx or rfdataxxx".
        subroutine read4t_Data(this,ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        type (fielddata), intent(inout) :: this
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
        character*7 name1
        character*8 name2
        character*9 name3

        name1 = 'rfdatax'
        name2 = 'rfdataxx'
        name3 = 'rfdataxxx'

        print*,"inside read4t: ",ifile
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif
        endif


          if(myrank.eq.0) then
            read(14,*,end=33)this%Ncex,this%Ncey,this%Ncez,&
                             this%Ncbx,this%Ncby,this%Ncbz
          endif
33        continue
          call MPI_BCAST(this%Ncex,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%Ncey,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%Ncez,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%Ncbx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%Ncby,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(this%Ncbz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          deallocate(this%coefex)
          deallocate(this%coefey)
          deallocate(this%coefez)
          deallocate(this%coefbx)
          deallocate(this%coefby)
          deallocate(this%coefbz)
          allocate(this%coefex(this%Ncex))
          allocate(this%coefey(this%Ncey))
          allocate(this%coefez(this%Ncez))
          allocate(this%coefbx(this%Ncbx))
          allocate(this%coefby(this%Ncby))
          allocate(this%coefbz(this%Ncbz))

        n = 0
        if(myrank.eq.0) then
50        continue
            read(14,*,end=77)tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
            n = n + 1
            this%coefex(n) = tmp1
            this%coefey(n) = tmp2
            this%coefez(n) = tmp3
            this%coefbx(n) = tmp4
            this%coefby(n) = tmp5
            this%coefbz(n) = tmp6
          goto 50
77        continue
          close(14)
          this%Ndatat = n
        endif

        call MPI_BCAST(this%Ndatat,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%coefex(1),this%Ndatat,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%coefey(1),this%Ndatat,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%coefez(1),this%Ndatat,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%coefbx(1),this%Ndatat,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%coefby(1),this%Ndatat,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%coefbz(1),this%Ndatat,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        print*,"Ndata: ",this%Ndatat
        end subroutine read4t_Data

        ! read in discrete Bz(r,z), Br(r,z) data from
        ! files "1Tx.T7 or 1Txx.T7 or 1Txxx.T7 for solenoid. Here, the grid
        ! is uniform in r and z.
        subroutine read2tsol_Data(this,ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        type (fielddata), intent(inout) :: this
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0
        character*6 name1
        character*7 name2
        character*8 name3
        double precision :: tmpmax

        name1 = '1Tx.T7'
        name2 = '1Txx.T7'
        name3 = '1Txxx.T7'
   
        mu0 = 4*2*asin(1.0)*1.0e-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(3:3) = char(ifile+48)
            open(14,file=name1,status='old')
            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(3:3) = char(ii+48)
            name2(4:4) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(3:3) = char(ii+48)
            name3(4:4) = char(kk+48)
            name3(5:5) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          ! the input range units are cm
          read(14,*,end=33)tmp1,tmp2,tmpint
          this%RminRft = tmp1/100.0
          this%RmaxRft = tmp2/100.0
          this%NrIntvRft = tmpint
          if(tmpint.ne.this%NrIntvRft) then
            print*,"input data wrong in R: ",this%NrIntvRft,tmpint
            stop
          endif
          ! the input range units are cm
          read(14,*,end=33)tmp1,tmp2,tmpint
          this%ZminRft = tmp1/100.0
          this%ZmaxRft = tmp2/100.0
          this%NzIntvRft = tmpint
          if(tmpint.ne.this%NzIntvRft) then
            print*,"input data wrong in Z: ",this%NzIntvRft,tmpint
            stop
          endif
          !print*,"Nz: ",NzIntvRft(ib),ZminRft(ib),ZmaxRft(ib)
          !print*,"Nr: ",NrIntvRft(ib),RminRft(ib),RmaxRft(ib)
        endif
33      continue
        call MPI_BCAST(this%NrIntvRft,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%NzIntvRft,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        deallocate(this%bzdatat)
        deallocate(this%brdatat)
        allocate(this%bzdatat(this%NrIntvRft+1,this%NzIntvRft+1))
        allocate(this%brdatat(this%NrIntvRft+1,this%NzIntvRft+1))

        if(myrank.eq.0) then
          tmpmax = -1.0e10
          n = 0
50        continue
              read(14,*,end=77)tmp1,tmp2
              nn = n+1
              j  = (nn-1)/(this%NrIntvRft+1) + 1
              i = mod((nn-1),this%NrIntvRft+1) + 1
              this%brdatat(i,j) = tmp1
              this%bzdatat(i,j) = tmp2
!              write(15,100)float(i-1),bzdatat(i,j)
              n = n + 1
          goto 50
77        continue
          close(14)
          Ndatalc = n
          !print*,"Ndata in 0: ",Ndatalc
          do i = 1, this%NzIntvRft+1 
            tmp3 = sqrt(this%brdatat(1,i)**2+this%bzdatat(1,i)**2)
            if(tmpmax.lt.tmp3) then
              tmpmax = tmp3
            endif
            write(15,100)float(i-1),this%brdatat(1,i),this%brdatat(2,i),&
                         this%bzdatat(1,i),this%bzdatat(2,i)
          enddo
          print*,"maximum B field on axix: ",tmpmax
          close(15)
        endif
100     format(5(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ZminRft,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ZmaxRft,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%RminRft,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%RmaxRft,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%bzdatat(1,1),Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%brdatat(1,1),Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndatalc
        end subroutine read2tsol_Data

        ! read in discrete complex Ex(x,y,z), Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z),
        ! Bz(x,y,z) rf data from
        ! files "1Tx.T7 or 1Txx.T7 or 1Txxx.T7. Here, the grid
        ! is uniform in x, y and z.
        ! This field can be used to model the traveling wave structure
        subroutine read3t_Data(this,ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        type (fielddata), intent(inout) :: this
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint,k
        double complex :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0,tmp5,tmp6
        real*8 :: tmp11,tmp22
        character*6 name1
        character*7 name2
        character*8 name3

        name1 = '1Tx.T7'
        name2 = '1Txx.T7'
        name3 = '1Txxx.T7'
   
        mu0 = 4*2*asin(1.0)*1.0e-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(3:3) = char(ifile+48)
            open(14,file=name1,status='old')
            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(3:3) = char(ii+48)
            name2(4:4) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(3:3) = char(ii+48)
            name3(4:4) = char(kk+48)
            name3(5:5) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          ! the input range units are m
          read(14,*,end=33)tmp11,tmp22,tmpint
          this%XminRfgt = tmp11
          this%XmaxRfgt = tmp22
          this%NxIntvRfgt = tmpint
          read(14,*,end=33)tmp11,tmp22,tmpint
          this%YminRfgt = tmp11
          this%YmaxRfgt = tmp22
          this%NyIntvRfgt = tmpint
          read(14,*,end=33)tmp11,tmp22,tmpint
          this%ZminRfgt = tmp11
          this%ZmaxRfgt = tmp22
          this%NzIntvRfgt = tmpint
          !print*,"Nz: ",NzIntvRft(ib),ZminRft(ib),ZmaxRft(ib)
          !print*,"Nr: ",NrIntvRft(ib),RminRft(ib),RmaxRft(ib)
        endif
33      continue
        call MPI_BCAST(this%NxIntvRfgt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%NyIntvRfgt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%NzIntvRfgt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%XminRfgt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%XmaxRfgt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%YminRfgt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%YmaxRfgt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ZminRfgt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%ZmaxRfgt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        deallocate(this%Exgridt)
        deallocate(this%Eygridt)
        deallocate(this%Ezgridt)
        deallocate(this%Bxgridt)
        deallocate(this%Bygridt)
        deallocate(this%Bzgridt)
        allocate(this%Exgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Eygridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Ezgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bxgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bygridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bzgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        this%Exgridt = cmplx(0.0d0,0.0d0)
        this%Eygridt = cmplx(0.0d0,0.0d0)
        this%Ezgridt = cmplx(0.0d0,0.0d0)
        this%Bxgridt = cmplx(0.0d0,0.0d0)
        this%Bygridt = cmplx(0.0d0,0.0d0)
        this%Bzgridt = cmplx(0.0d0,0.0d0)

        if(myrank.eq.0) then
          n = 0
50        continue
              read(14,*,end=77)tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
              n = n+1
              k = (n-1)/((this%NxIntvRfgt+1)*(this%NyIntvRfgt+1))+1
              j = (n-(k-1)*(this%NxIntvRfgt+1)*(this%NyIntvRfgt+1)-1)/(this%NxIntvRfgt+1) + 1
              i = n - (k-1)*(this%NxIntvRfgt+1)*(this%NyIntvRfgt+1) - &
                  (j-1)*(this%NxIntvRfgt+1)
              this%Exgridt(i,j,k) = tmp1
              this%Eygridt(i,j,k) = tmp2
              this%Ezgridt(i,j,k) = tmp3
              this%Bxgridt(i,j,k) = tmp4
              this%Bygridt(i,j,k) = tmp5
              this%Bzgridt(i,j,k) = tmp6
          goto 50
77        continue
          close(14)
          Ndatalc = n
          !print*,"Ndata in 0: ",Ndatalc
        endif
100     format(2(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%Exgridt(1,1,1),Ndatalc,MPI_DOUBLE_COMPLEX,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%Eygridt(1,1,1),Ndatalc,MPI_DOUBLE_COMPLEX,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%Ezgridt(1,1,1),Ndatalc,MPI_DOUBLE_COMPLEX,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%Bxgridt(1,1,1),Ndatalc,MPI_DOUBLE_COMPLEX,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%Bygridt(1,1,1),Ndatalc,MPI_DOUBLE_COMPLEX,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(this%Bzgridt(1,1,1),Ndatalc,MPI_DOUBLE_COMPLEX,0,&
             MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndata
        end subroutine read3t_Data

      end module Dataclass
