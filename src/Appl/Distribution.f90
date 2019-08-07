!----------------------------------------------------------------
! (c) Copyright, 2013 by the Regents of the University of California.
! Distributionclass: Initial distribution of charged beam bunch class in 
!                    Beam module of APPLICATION layer.
! MODULE  : ... Distributionclass
! VERSION : ... 1.7
! DATE  : ... 01/07/2013
!> @author 
!> Ji Qiang  
! DESCRIPTION:
!> This class defines initial distributions for the charged particle beam bunch information in the accelerator.
!        
! Comments:
!----------------------------------------------------------------
      module Distributionclass
        use Pgrid2dclass
        use CompDomclass
        use Timerclass
        use NumConstclass
        use PhysConstclass
        use BeamBunchclass
      contains
        
        !--------------------------------------------------------------------------------------
        !> @brief
        !> sample the particles with intial distribution.
        !--------------------------------------------------------------------------------------
        subroutine sample_Dist(this,distparam,nparam,flagdist,geom,grid,Flagbc,ib,Nb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (BeamBunch), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: flagdist,ib,Nb
        integer :: myid, myidx, myidy,seedsize,i,isize
        !integer seedarray(1)
        !integer*8, allocatable, dimension(:) :: seedarray
        integer, allocatable, dimension(:) :: seedarray
        integer :: totnp,npx,npy,meanpts20
        real rancheck

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        meanpts20 = (this%Npt/totnp)*20
!        seedarray(1)=(100001+myid)*(myid+7)
!        call random_seed(put=seedarray(1:1))
!        write(6,*)'seedarray=',seedarray

!this is done in AccSimulator
!        call random_seed(SIZE=seedsize)
!        allocate(seedarray(seedsize))
!        do i = 1, seedsize
!          !seedarray(i) = (1000+5*myid)*(myid+7)+i-1 !seed 1
!          !seedarray(i) = (2000+5*myid)*(myid+7)+i-1  !seed 2
!          !seedarray(i) = (3000+5*myid)*(myid+7)+i-1  !seed 3
!          !seedarray(i) = ib*400000+(3000+5*myid)*(myid+7)+i-1  !seed 3
!          !seedarray(i) = ib*(300+200*myid)*(myid+7)+i-1  !seed 6
!          seedarray(i) = 10.0d0 + myid*1.0d0*meanpts20+i*1.0d0*myid
!        enddo
!        call random_seed(PUT=seedarray)
!        call random_number(rancheck)
!        !randomnize the random number generator
!        do i = 1, 1000
!          call random_number(rancheck)
!        enddo
!!        write(6,*)'myid,rancheck=',seedarray,myid,rancheck

        !from the SI (m) unit to dimensionless unit
        distparam(1) = distparam(1)/Scxlt
        distparam(6) = distparam(6)/Scxlt
        distparam(8) = distparam(8)/Scxlt
        distparam(13) = distparam(13)/Scxlt
        distparam(15) = distparam(15)/Scxlt
        distparam(20) = distparam(20)/Scxlt
        if(flagdist.eq.1) then
          !6d uniform distribution
          call Uniform_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.2) then
          !6d Gaussian distribution
          call Gauss3_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.3) then
          !6d Waterbag distribution
          call Waterbag_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.4) then
          !3d Waterbag distribution in spatial and 3d Gaussian distribution in
          !momentum space
          call Semigauss_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.5) then
          !transverse KV distribution and longitudinal uniform distribution
          call KV3d_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.16) then
          !read in an initial distribution with format from IMPACT-T
          call read_Dist(this,nparam,distparam,ib)
        else if(flagdist.eq.24) then
          call readParmela_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.25) then
          call readElegant_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.27) then
          call CylcoldZSob_Dist(this,nparam,distparam,grid)
        else
          call Combine_Dist(this,nparam,distparam,grid,ib,Nb,flagdist)
        endif

!        deallocate(seedarray)

        end subroutine sample_Dist
        
        !--------------------------------------------------------------------------------------
        !> @brief
        !> 6d uniform distribution
        !--------------------------------------------------------------------------------------
        subroutine Uniform_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6,1) :: a
        double precision, dimension(2) :: x1, x2
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj
        double precision :: t0,gamma,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call random_number(x11)
        print*,"x11: ",x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = avgpts

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(6,avgpts))
        this%Pts1 = 0.0
!        print*,"avgpts: ",avgpts
    
        do ii = 1, avgpts
          call random_number(r1)
          r1 = (2*r1 - 1.0d0)*sqrt(3.0d0)
          call random_number(r2)
          r2 = (2*r2 - 1.0d0)*sqrt(3.0d0)
          this%Pts1(1,ii) = xmu1 + sig1*r1/sq12
          this%Pts1(2,ii) = xmu2 + sig2*(-muxpx*r2/sq12+r2)
          call random_number(r1)
          r1 = (2*r1 - 1.0d0)*sqrt(3.0d0)
          call random_number(r2)
          r2 = (2*r2 - 1.0d0)*sqrt(3.0d0)
          this%Pts1(3,ii) = xmu3 + sig3*r1/sq34
          this%Pts1(4,ii) = xmu4 + sig4*(-muypy*r2/sq34+r2)
          call random_number(r1)
          r1 = (2*r1 - 1.0d0)*sqrt(3.0d0)
          call random_number(r2)
          r2 = (2*r2 - 1.0d0)*sqrt(3.0d0)
          this%Pts1(5,ii) = xmu5 + sig5*r1/sq56
          this%Pts1(6,ii) = xmu6 + sig6*(-muzpz*r2/sq56+r2)
        enddo

        this%Nptlocal = numpts0

!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Uniform_Dist

        !--------------------------------------------------------------------------------------
        !> @brief
        !> 6d Gaussian distribution
        !--------------------------------------------------------------------------------------
        subroutine Gauss3_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp
!        integer seedarray(1)
        double precision :: t0,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        if(mod(avgpts,10).ne.0) then
!          print*,"The number of particles has to be an integer multiple of 10Nprocs"
!          stop
        endif
        
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(6,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
        intvsamp = 1
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo
          
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss3_Dist

        subroutine normdv(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(2), intent(out) :: y
        double precision :: twopi,x1,x2,epsilon

        epsilon = 1.0e-18

        twopi = 4.0d0*asin(1.0d0)
        call random_number(x2)
10      call random_number(x1)
!        x1 = 0.5
!10      x2 = 0.6
        if(x1.eq.0.0d0) goto 10
!        if(x1.eq.0.0d0) x1 = epsilon
        y(1) = sqrt(-2.0d0*log(x1))*cos(twopi*x2)
        y(2) = sqrt(-2.0d0*log(x1))*sin(twopi*x2)

        end subroutine normdv

        subroutine normVec(y,num)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: num
        double precision, dimension(2,num), intent(out) :: y
        double precision :: twopi,epsilon
        double precision, dimension(num) :: x1,x2
        integer :: i

        epsilon = 1.0d-16

        twopi = 4.0d0*asin(1.0d0)
        call random_number(x2)
        call random_number(x1)
        do i = 1, num
          if(x1(i).eq.0.0d0) x1(i) = epsilon
          y(1,i) = sqrt(-2.0d0*log(x1(i)))*cos(twopi*x2(i))
          y(2,i) = sqrt(-2.0d0*log(x1(i)))*sin(twopi*x2(i))
        enddo

        end subroutine normVec

        !--------------------------------------------------------------------------------------
        !> @brief
        !> 6d Waterbag distribution.
        !> sample the particles with intial distribution using rejection method. 
        !--------------------------------------------------------------------------------------
        subroutine Waterbag_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy,iran,intvsamp
!        integer seedarray(2)
        double precision :: t0,x11
        double precision, allocatable, dimension(:) :: ranum6

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        seedarray(2)=(101+2*myid)*(myid+4)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray)
        call random_number(x11)
        print*,"x11: ",myid,x11

        avgpts = this%Npt/(npx*npy)
        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs" 
        !  stop
        !endif
 
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(6,avgpts))
          this%Pts1 = 0.0
        endif
        numpts = 0
        isamz = 0
        isamy = 0
        intvsamp = avgpts
        !intvsamp = 10
        allocate(ranum6(6*intvsamp))

        do 
          ! rejection sample.
10        continue 
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call random_number(ranum6)
          endif
          iran = 6*mod(isamz-1,intvsamp)
          r1 = 2.0d0*ranum6(iran+1)-1.0
          r2 = 2.0d0*ranum6(iran+2)-1.0
          r3 = 2.0d0*ranum6(iran+3)-1.0
          r4 = 2.0d0*ranum6(iran+4)-1.0
          r5 = 2.0d0*ranum6(iran+5)-1.0
          r6 = 2.0d0*ranum6(iran+6)-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0d0) goto 10
          isamy = isamy + 1
          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!x-px:
          x1 = r1*sqrt(8.0d0)
          x2 = r2*sqrt(8.0d0)
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          x3 = r3*sqrt(8.0d0)
          x4 = r4*sqrt(8.0d0)
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          x5 = r5*sqrt(8.0d0)
          x6 = r6*sqrt(8.0d0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo

        deallocate(ranum6)
          
        this%Nptlocal = avgpts
        print*,"avgpts: ",avgpts
       
!        print*,avgpts,isamz,isamy

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Waterbag_Dist

        !--------------------------------------------------------------------------------------
        !> @brief
        !> transverse KV distribution and longitudinal uniform distribution
        !--------------------------------------------------------------------------------------
        subroutine KV3d_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy
!        integer seedarray(1)
        double precision :: t0,x11,twopi

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(6,avgpts))
        this%Pts1 = 0.0
        twopi = 4*asin(1.0d0)

        do numpts = 1, avgpts
          call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          r4 = sqrt(r1)
          r5 = sqrt(1.0-r1)
          r2 = r2*twopi
          r3 = r3*twopi
          x1 = 2*r4*cos(r2)
          x2 = 2*r4*sin(r2)
          x3 = 2*r5*cos(r3)
          x4 = 2*r5*sin(r3)
!x-px:
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation.
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          call random_number(r5)
          r5 = 2*r5 - 1.0
          call random_number(r6)
          r6 = 2*r6 - 1.0
          x5 = r5*sqrt(3.0d0)
          x6 = r6*sqrt(3.0d0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
        enddo
          
        this%Nptlocal = avgpts
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine KV3d_Dist

        !--------------------------------------------------------------------------------------
        !> @brief
        !> 3d Waterbag distribution in spatial and 3d Gaussian distribution in
        !> momentum space
        !--------------------------------------------------------------------------------------
        subroutine Semigauss_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6,2) :: a
        double precision, dimension(3) :: x1, x2
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2, r, r3
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj
        double precision :: t0

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(6,avgpts))
        this%Pts1 = 0.0
        do ii = 1, avgpts
          ! rejection sample.
10        call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          r1 = 2.0d0*r1-1.0
          r2 = 2.0d0*r2-1.0
          r3 = 2.0d0*r3-1.0
          if(r1**2+r2**2+r3**2.gt.1.0d0) goto 10
          x2(1) = r1
          x2(2) = r2
          x2(3) = r3
          call normdv2(x1)

          !x-px:
!         Correct Gaussian distribution.
          a(1,1) = xmu1 + sig1*x2(1)/sq12*sqrt(5.0d0)
          a(2,1) = xmu2 + sig2*(-muxpx*x2(1)/sq12+x1(1))
!         Rob's Gaussian distribution.
          !a(1,1) = xmu1 + sig1*x2(1)*sqrt(5.0d0)
          !a(2,1) = xmu2 + sig2*(muxpx*x2(1)+sq12*x1(1))
          !y-py
!         Correct Gaussian distribution.
          a(3,1) = xmu3 + sig3*x2(2)/sq34*sqrt(5.0d0)
          a(4,1) = xmu4 + sig4*(-muypy*x2(2)/sq34+x1(2))
!         Rob's Gaussian distribution.
          !a(3,1) = xmu3 + sig3*x2(2)*sqrt(5.0d0)
          !a(4,1) = xmu4 + sig4*(muypy*x2(2)+sq34*x1(2))
          !z-pz
!         Correct Gaussian distribution.
          a(5,1) = xmu5 + sig5*x2(3)/sq56*sqrt(5.0d0)
          a(6,1) = xmu6 + sig6*(-muzpz*x2(3)/sq56+x1(3))
!         Rob's Gaussian distribution.
          !a(5,1) = xmu5 + sig5*x2(3)*sqrt(5.0d0)
          !a(6,1) = xmu6 + sig6*(muzpz*x2(3)+sq56*x1(3))

          do j = 1, 6
             this%Pts1(j,ii) = a(j,1)
          enddo

        enddo

        numpts0 = avgpts

        this%Nptlocal = numpts0
!        print*,"numpts0: ",numpts0

!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Semigauss_Dist

        subroutine normdv2(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(3), intent(out) :: y
        double precision :: sumtmp,x
        integer :: i

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(1) = sumtmp - 6.0

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(2) = sumtmp - 6.0

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(3) = sumtmp - 6.0

        end subroutine normdv2

        !--------------------------------------------------------------------------------------
        !> @brief
        !> read in an initial distribution with format (x(m), px/mc, y(m), py/mc,...)
        !--------------------------------------------------------------------------------------
        subroutine read_Dist(this,nparam,distparam,ib)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,ib
        double precision, dimension(nparam) :: distparam
        integer :: i,j,jlow,jhigh,avgpts,myid,nproc,ierr,nptot,nleft
        double precision, dimension(6) :: tmptcl
        double precision :: sum1,sum2
        character*12 name1
        character*13 name2
        character*14 name3
        integer :: k,l
 
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

        name1 = 'partclx.data'
        name2 = 'partclxx.data'
        name3 = 'partclxxx.data'

         if(ib.eq.1) then
            open(unit=12,file='partcl.data',status='old')
         else if(ib.lt.10) then
            name1(7:7) = char(ib+48)
            open(unit=12,file=name1,status='old')
         else if(ib.lt.100) then
            i = ib/10
            j = ib - 10*i
            name2(7:7) = char(i+48)
            name2(8:8) = char(j+48)
            open(unit=12,file=name2,status='old')
         else if(ib.lt.1000) then
            i = ib/100
            j = ib - 100*i
            k = j/10
            l = j - 10*k
            name3(7:7) = char(i+48)
            name3(8:8) = char(k+48)
            name3(9:9) = char(l+48)
            open(unit=12,file=name3,status='old')
          else
            print*,"over maximum # of input particle files:...."
            stop
          endif
 
        sum1 = 0.0
        sum2 = 0.0
 
          read(12,*)nptot
          avgpts = nptot/nproc
          nleft = nptot - avgpts*nproc
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          allocate(this%Pts1(6,avgpts))
          this%Pts1 = 0.0
          !jlow = myid*avgpts + 1
          !jhigh = (myid+1)*avgpts
          print*,"avgpts, jlow, and jhigh: ",avgpts,jlow,jhigh
          do j = 1, nptot
            read(12,*)tmptcl(1:6)
            sum1 = sum1 + tmptcl(1)
            sum2 = sum2 + tmptcl(3)
            if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              this%Pts1(1:6,i) = tmptcl(1:6)
            endif
!            if(myid.eq.0) print*,i,sum1,sum2
          enddo
          print*,"sumx1,sumy1: ",sum1/nptot,sum2/nptot
 
          close(12)
 
          this%Nptlocal = avgpts
          !change length to the dimensionless unit
          do i = 1, avgpts
            this%Pts1(1,i) = this%Pts1(1,i)/Scxlt + distparam(6)
            this%Pts1(2,i) = this%Pts1(2,i) + distparam(7)
            this%Pts1(3,i) = this%Pts1(3,i)/Scxlt + distparam(13)
            this%Pts1(4,i) = this%Pts1(4,i) + distparam(14)
            this%Pts1(5,i) = this%Pts1(5,i)/Scxlt + distparam(20)
            this%Pts1(6,i) = this%Pts1(6,i)  + distparam(21)
          enddo
 
        end subroutine read_Dist

        subroutine readParmela_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        real*8 :: beta0,sumeng
        integer :: jlow,jhigh,nleft,avgpts

        pi = 2*asin(1.0d0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)*Scxlt
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)*Scxlt
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)*Scxlt
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(7,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        sumeng = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:7,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
          sumeng = sumeng + Ptcl(6,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        sumeng = sumeng/inipts
        print*,"sumx2: ",sumx2,sumy2,sumeng
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxlt
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        gamma0 = 1.0+sumeng*1.0e6/mccq  
        beta0 = sqrt(1.0-1./gamma0**2)
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0

        avgpts = inipts/totnp
        nleft = nptot - avgpts*totnp
        if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
        else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
        endif

        allocate(this%Pts1(6,avgpts))
        this%Pts1 = 0.0

        do j = 1, inipts
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
            i = j - jlow + 1
            this%Pts1(1,i) = (Ptcl(1,j)/100.0d0)*xscale + xmu1
            this%Pts1(3,i) = (Ptcl(3,j)/100.0d0)*yscale + xmu3
            gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
            betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
            this%Pts1(2,i) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
            this%Pts1(4,i) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
            this%Pts1(5,i) = -Ptcl(5,j)/Rad2deg*beta0*Clight/(2*Pi*Scfreq)*zscale + xmu5
            this%Pts1(6,i) = gamma*betaz
          endif
        enddo

        this%Nptlocal = avgpts

        do i = 1, this%Nptlocal
          this%Pts1(1,i) = this%Pts1(1,i)/Scxlt
          this%Pts1(3,i) = this%Pts1(3,i)/Scxlt
          this%Pts1(5,i) = this%Pts1(5,i)/Scxlt
        enddo

        deallocate(Ptcl)

        end subroutine readParmela_Dist

        subroutine readElegant_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        real*8 :: beta0,sumeng
        integer :: jlow,jhigh,nleft,avgpts

        pi = 2*asin(1.0d0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)*Scxlt
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)*Scxlt
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)*Scxlt
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(7,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        sumeng = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:7,i)
          sumx = sumx + Ptcl(2,i)
          sumx2 = sumx2 + Ptcl(2,i)*Ptcl(2,i)
          sumy = sumy + Ptcl(4,i)
          sumy2 = sumy2 + Ptcl(4,i)*Ptcl(4,i)
          sumeng = sumeng + Ptcl(7,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        sumeng = sumeng/inipts
        print*,"sumx2: ",sumx2,sumy2,sumeng
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxlt
        mccq = this%Mass
        gamma0 = 1.0+sumeng*1.0e6/mccq  
        beta0 = sqrt(1.0-1./gamma0**2)
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0

        avgpts = inipts/totnp
        nleft = nptot - avgpts*totnp
        if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
        else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
        endif

        allocate(this%Pts1(6,avgpts))
        this%Pts1 = 0.0

        do j = 1, inipts
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
            i = j - jlow + 1
            this%Pts1(1,i) = Ptcl(2,j)*xscale + xmu1
            this%Pts1(3,i) = Ptcl(4,j)*yscale + xmu3
            this%Pts1(5,i) = -Ptcl(6,j)*Clight*zscale + xmu5
            gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
            this%Pts1(2,i) = Ptcl(3,j)*gammabet*pxscale + xmu2
            this%Pts1(4,i) = Ptcl(5,j)*gammabet*pyscale + xmu4
            this%Pts1(6,i) = gammabet*pzscale + xmu6
          endif
        enddo

        this%Nptlocal = avgpts

        do i = 1, this%Nptlocal
          this%Pts1(1,i) = this%Pts1(1,i)/Scxlt
          !this%Pts1(2,i) = 0.0
          this%Pts1(3,i) = this%Pts1(3,i)/Scxlt
          !this%Pts1(4,i) = 0.0
          this%Pts1(5,i) = this%Pts1(5,i)/Scxlt
          !this%Pts1(6,i) = xmu6
        enddo

        deallocate(Ptcl)

        end subroutine readElegant_Dist

        subroutine normdv1d(y)
        implicit none
        include 'mpif.h'
        double precision, intent(out) :: y
        double precision :: sumtmp,x
        integer :: i
 
        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y = sumtmp - 6.0
 
        end subroutine normdv1d

        subroutine CylcoldZSob_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,x1,x3,cs,ss
        double precision, allocatable, dimension(:) :: r1,r2,r3,r4 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,ierr
!        integer seedarray(1)
        double precision :: t0,x11,twopi,tmpmax,tmpmaxgl,shiftz
        real*8 :: vx,vy,r44,vr1,vr2,vzmax,r,fvalue
        integer :: isamz
!for quiet start of the modulated uniform current profile in z
        integer :: nptsob,k,ilow,ihigh,Nmax,iseed,ndim
        real*8 :: eps,epsilon,xz,xmod,rk,psi
        real*8, dimension(6) :: xtmp

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17) !for modulization amplitude
        zscale = distparam(18) !for wavelength (m)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        do i = 1, 3000
          call random_number(x11)
        enddo
!        print*,myid,x11

        if(this%Npt<=100000000) then
          nptsob = this%Npt
        else
          nptsob = 100000000
          print*,"maximum number of total particle is ",nptsob
        endif

        !avgpts = this%Npt/(npx*npy)
        avgpts = nptsob/(npx*npy)
        nptsob = avgpts*npx*npy
        this%Npt = nptsob

        ilow = myid*avgpts
        ihigh = (myid+1)*avgpts

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(6,avgpts))
        twopi = 4*asin(1.0d0)

        xmod = muzpz
        rk = twopi/(zscale/Scxlt)
        Nmax = 500
        eps = 1.0d-8
        epsilon = 1.0e-18

        ndim = 6
        iseed = 100
        call sobseq(-2,xtmp,iseed)
        do k = 1, nptsob
          call sobseq(ndim,xtmp,iseed)
          if(k.gt.ilow .and. k.le.ihigh) then
            i = k - ilow
            x1 = sig1*sqrt(xtmp(1))
            x3 = sig3*sqrt(xtmp(1))
            cs = cos(xtmp(2)*twopi)
            ss = sin(xtmp(2)*twopi)
            this%Pts1(1,i) = xmu1 + x1*cs
            this%Pts1(3,i) = xmu3 + x3*ss
            if(xtmp(3).eq.0.0d0) xtmp(3) = epsilon
            this%Pts1(2,i) = xmu2 + sig2* &
                     sqrt(-2.0d0*log(xtmp(3)))*cos(twopi*xtmp(4))
            this%Pts1(4,i) = xmu4 + sig4* &
                     sqrt(-2.0d0*log(xtmp(3)))*sin(twopi*xtmp(4))
            xz = (2*xtmp(5)-1.0d0)*sigz
            call psiroot(rk,xz,xmod,psi,eps,Nmax)
            this%Pts1(5,i) = xmu5 + psi
            if(xtmp(6).eq.0.0d0) xtmp(6) = epsilon
            this%Pts1(6,i) = xmu6 + sigpz*sqrt(-2.0d0*log(xtmp(6)))
          endif
        enddo
        
        this%Nptlocal = avgpts
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine CylcoldZSob_Dist

      subroutine psiroot(rk,xx,xmod,psi,eps,Nmax)
      implicit none
      integer :: Nmax
      real*8 :: rk,xx,xmod,eps,psi
      integer :: i
      real*8 :: ps0,ps1,fps0,dfps0
 
      ps0 = xx
 
      do i = 1, Nmax
        fps0 = ps0+xmod*sin(rk*ps0)/rk-xx
        !print*,"ps0: ",i,ps0,fps0
        if(abs(fps0).le.eps) then
          psi = ps0
          return
        else
          dfps0 = 1.0d0+xmod*cos(rk*ps0)
          ps1 = ps0 - fps0/dfps0
          ps0 = ps1
        endif
      enddo
 
      if(i.ge.Nmax) then
        print*,"Not converged in psiroot"
        stop
      endif
 
      end subroutine psiroot

      SUBROUTINE sobseq(n,x,iseed)
      INTEGER n,MAXBIT,MAXDIM,iseed
      REAL*8 x(*)
      PARAMETER (MAXBIT=30,MAXDIM=6)
      INTEGER i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT),&
      iv(MAXBIT*MAXDIM),ix(MAXDIM),mdeg(MAXDIM)
      REAL*8 fac
      SAVE ip,mdeg,ix,iv,in,fac
      EQUIVALENCE (iv,iu)
      DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
      DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
      if (n.lt.0) then
        do 14 k=1,MAXDIM
          do 11 j=1,mdeg(k)
            iu(k,j)=iu(k,j)*2**(MAXBIT-j)
11        continue
          do 13 j=mdeg(k)+1,MAXBIT
            ipp=ip(k)
            i=iu(k,j-mdeg(k))
            i=ieor(i,i/2**mdeg(k))
            do 12 l=mdeg(k)-1,1,-1
              if(iand(ipp,1).ne.0)i=ieor(i,iu(k,j-l))
              ipp=ipp/2
12          continue
            iu(k,j)=i
13        continue
14      continue
        fac=1./2.**MAXBIT
        in=iseed
      else
        im=in
        do 15 j=1,MAXBIT
          if(iand(im,1).eq.0)goto 1
          im=im/2
15      continue
        pause 'MAXBIT too small in sobseq'
1       im=(j-1)*MAXDIM
        do 16 k=1,min(n,MAXDIM)
          ix(k)=ieor(ix(k),iv(im+k))
          x(k)=ix(k)*fac
16      continue
        in=in+1
      endif
      return
      END subroutine sobseq

      !-----------------------------------------------------------------------------------------------
      !> @brief
      !> generating initial particle distribution based on the combination of
      !> transverse spatial distribution (2 types), longitudinal spatial distribution (3 types)
      !> and 3D momentum distribution (4 types). For example, flagdist = 111 denotes
      !> type 1 from transverse spatial distribution, type 1 from longitudinal spatial distribution,
      !> and type 1 from 3D momentum distribution. In this case, it denotes a transverse uniform
      !> elliptical, longitudinal flat-top with linear ramping in spatial, and 3d full Gaussian 
      !> distribution in momentum space.
      !> flagdist = ijk, where 
      !> i = 1 (transverse uniform ellipse) and 2 (Gaussian with cut-off);
      !> j = 1 (flat-top with linear ramping), 2 (flat-top with 2sigma Gaussain ramping), and
      !>     3 (Gaussian with cut-off);
      !> k = 1, (3D Gaussian momentum), 2 (transverse Gaussian momentum, longitudinal semi-Gaussian), 
      !>     3 (streak camera model), and 4 (3 step model).
      !------------------------------------------------------------------------------------------------
        subroutine Combine_Dist(this,nparam,distparam,grid,ib,Nb,flagdist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,ib,Nb,flagdist
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        integer :: totnp,npy,npx
        integer :: avgpts,avgpts0
        integer :: myid,myidx,myidy,i,ierr
        double precision :: t0,x11
        real*8, allocatable, dimension(:) :: tmptcl1
        real*8, allocatable, dimension(:,:) :: tmptcl2
        real*8, allocatable, dimension(:,:) :: tmptcl3
        real*8 :: currtmp,zflat,zrise,zrisehalf,cutx,cuty,cutz
        integer :: iptgl,iptlc,npt0,ii,jj,kk,ierrdist
        real*8 :: emax,emass,wkf,Eph,Ef,Tem,Ewk

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        do i = 1, 3000
          call random_number(x11)
        enddo
        print*,myid,x11
 
        avgpts = this%Npt/(npx*npy)
        avgpts0 = this%Npt/(npx*npy)*Nb

        ii = flagdist/100
        jj = (flagdist-ii*100)/10
        kk = flagdist-ii*100-jj*10
    
        ierrdist = 0

        allocate(tmptcl1(avgpts0))
        !longitudinal spatial distribution (behind cathode z=0)
        if(jj.eq.1) then
          zflat = distparam(15)
          zrise = distparam(18)/Scxlt
          call LongFlattoplinz_Dist(tmptcl1,avgpts,avgpts0,zrise,zflat,xmu5,ib,Nb)
          iptlc = avgpts
          npt0 = this%Npt
          call MPI_ALLREDUCE(iptlc,iptgl,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
          this%Npt = iptgl
          currtmp = this%Current
          this%Current = currtmp*this%Npt/dble(npt0)
        else if(jj.eq.2) then
          zflat = distparam(15)
          zrisehalf = distparam(18)/Scxlt
          call LongFlattopGasz_Dist(tmptcl1,avgpts,avgpts0,zrisehalf,zflat,xmu5,ib,Nb)
          iptlc = avgpts
          npt0 = this%Npt
          call MPI_ALLREDUCE(iptlc,iptgl,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
          this%Npt = iptgl
          currtmp = this%Current
          this%Current = currtmp*this%Npt/dble(npt0)
        else if(jj.eq.3) then
          cutz = distparam(18)
          call LongGaussz_Dist(tmptcl1,avgpts,sigz,xmu5,cutz)
        else
          print*,"wrong initial distribution code in jj!"
          ierrdist = 1
        endif
        !transverse spatial distribution
        allocate(tmptcl2(2,avgpts))
        if(ii.eq.1) then
          call TranUnifxy_Dist(tmptcl2,avgpts,sigx,sigy,xmu1,xmu3)
        else if(ii.eq.2) then
          cutx = distparam(4)
          cuty = distparam(11)
          call TranGaussxy_Dist(tmptcl2,avgpts,sigx,sigy,xmu1,xmu3,cutx,cuty)
        else
          print*,"wrong initial distribution code in ii!"
          ierrdist = 1
        endif
        !3D momentum distribution
        allocate(tmptcl3(3,avgpts))
        if(kk.eq.1) then
          call PxPyPzGauss3d_Dist(tmptcl3,avgpts,sigpx,sigpy,sigpz,xmu2,xmu4,xmu6)
        else if(kk.eq.2) then
          call PxPyGaussPzSg_Dist(tmptcl3,avgpts,sigpx,sigpy,sigpz,xmu2,xmu4,xmu6)
        else if(kk.eq.3) then
          emax = distparam(16) !maximum energy (eV)
          wkf = distparam(19) !work function (eV)
          emass = this%Mass ! (eV)
          call PxPyPzStreakCam_Dist(tmptcl3,avgpts,emax,wkf,emass,xmu2,xmu4,xmu6)
        else if(kk.eq.4) then
          Eph = distparam(16) !photon energy (eV)
          Ef = distparam(17)  !Fermi energy level (eV)
          !Tem = distparam(18) !cathode temperature (eV)
          Tem = 0.026d0 !room temperature (eV)
          Ewk = distparam(19) !effective work function (eV)
          emass = this%Mass
          call PxPyPz3step_Dist(tmptcl3,avgpts,Ef,Eph,Tem,Ewk,emass,xmu2,xmu4,xmu6)
        else
          print*,"wrong initial distribution code in kk!"
          ierrdist = 1
        endif

        this%Nptlocal = avgpts
        if(ierrdist.ne.1) then
          allocate(this%Pts1(6,avgpts))
          do i = 1, avgpts
            this%Pts1(1,i) = tmptcl2(1,i)
            this%Pts1(2,i) = tmptcl3(1,i)
            this%Pts1(3,i) = tmptcl2(2,i)
            this%Pts1(4,i) = tmptcl3(2,i)
            this%Pts1(5,i) = tmptcl1(i)
            this%Pts1(6,i) = tmptcl3(3,i)
          enddo
        else
          print*,"wrong initial distribution code!"
          stop
        endif

        deallocate(tmptcl1)
        deallocate(tmptcl2)
        deallocate(tmptcl3)

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Combine_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> transverse spatial uniform elleptical distribution. (id=1)
        !------------------------------------------------------------------------------------------------
        subroutine TranUnifxy_Dist(ptxy,avgpts,sigx,sigy,xmu1,xmu3)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts
        real*8, intent(out), dimension(2,avgpts) :: ptxy
        real*8, intent(in) :: sigx,sigy,xmu1,xmu3
        double precision, dimension(avgpts) :: r1,r2
        real*8 :: x1,x3,twopi,cs,ss
        integer :: numpts
  
        twopi = 4*asin(1.0d0)
        call random_number(r1)
        call random_number(r2)

        do numpts = 1, avgpts
          x1 = sigx*sqrt(r1(numpts))
          x3 = sigy*sqrt(r1(numpts))
          cs = cos(r2(numpts)*twopi)
          ss = sin(r2(numpts)*twopi)
          ptxy(1,numpts) = xmu1 + x1*cs
          ptxy(2,numpts) = xmu3 + x3*ss
        enddo

        end subroutine TranUnifxy_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> transverse spatial Gaussian distribution with cut-off. (id=2)
        !------------------------------------------------------------------------------------------------
        subroutine TranGaussxy_Dist(ptxy,avgpts,sigx,sigy,xmu1,xmu3,cutx,cuty)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts
        real*8, intent(out), dimension(2,avgpts) :: ptxy
        real*8, intent(in) :: sigx,sigy,xmu1,xmu3,cutx,cuty
        integer :: numpts
  
        call normVecCut(ptxy,avgpts,cutx,cuty)

        do numpts = 1, avgpts
          ptxy(1,numpts) = xmu1 + sigx*ptxy(1,numpts)
          ptxy(2,numpts) = xmu3 + sigy*ptxy(2,numpts)
        enddo

        end subroutine TranGaussxy_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> longitudinal spatial flat-top distribution with linear ramp (id=1).
        !------------------------------------------------------------------------------------------------
        subroutine LongFlattoplinz_Dist(ptz,avgpts,avgpts0,zscale,pzscale,xmu5,ib,Nb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts0,ib,Nb
        integer, intent(inout) :: avgpts
        real*8, intent(out), dimension(avgpts0) :: ptz
        real*8, intent(in) :: xmu5,zscale,pzscale
        real*8, dimension(avgpts0) :: tmptcl
        integer :: numpts,iptlc
        real*8 :: zrise,zflat,totleng,deltaz,zminbin,zmaxbin,zmingl,z1,z2,z3
        real*8 :: rr,zz,rr2,fvalue
  
        zrise = zscale !linear rise part
        zflat = pzscale !flat-top part

        totleng = zflat + 2*zrise
        deltaz = totleng/Nb
        zminbin = -ib*deltaz
        zmaxbin = -(ib-1)*deltaz
        zmingl = -totleng
        z1 = zmingl + zrise
        z2 = z1 + zflat
        z3 = 0.0d0
        numpts = 0
        iptlc = 0
        tmptcl = 0.0d0

        do
10        call random_number(rr)
          zz = -totleng*rr
          if(zz.ge.zmingl .and. zz.lt.z1 ) then
            fvalue = (zz-zmingl)/(z1-zmingl)
          else if( zz.ge.z1 .and. zz.lt.z2) then
            fvalue = 1.0d0
          else if( zz.ge.z2 .and. zz.le.z3) then
            fvalue = 1.0d0 - (zz-z2)/(z3-z2)
          else
          endif
          call random_number(rr2)
          if(rr2.gt.fvalue) goto 10
          numpts = numpts + 1
          if( (zz.le.zmaxbin) .and. (zz.gt.zminbin) ) then
            iptlc = iptlc + 1
            tmptcl(iptlc) = zz
          endif
          if(numpts.ge.avgpts0) exit
        enddo

        avgpts = iptlc

        do numpts = 1, avgpts
         ptz(numpts) = xmu5 + tmptcl(numpts)
        enddo

        end subroutine LongFlattoplinz_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> longitudinal spatial flat-top distribution with 2sigma Gaussian ramp. (id=2)
        !------------------------------------------------------------------------------------------------
        subroutine LongFlattopGasz_Dist(ptz,avgpts,avgpts0,zscale,pzscale,xmu5,ib,Nb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts0,ib,Nb
        integer, intent(inout) :: avgpts
        real*8, intent(out), dimension(avgpts0) :: ptz
        real*8, intent(in) :: xmu5,zscale,pzscale
        real*8, dimension(avgpts0) :: tmptcl
        integer :: numpts,iptlc
        real*8 :: zrise,zflat,totleng,deltaz,zminbin,zmaxbin,zmingl,z1,z2,z3
        real*8 :: rr,zz,rr2,fvalue
  

        zrise = 2*zscale !2sigma Gaussian rise part
        zflat = pzscale !flat-top part

        totleng = zflat + 2*zrise
        deltaz = totleng/Nb
        zminbin = -ib*deltaz
        zmaxbin = -(ib-1)*deltaz
        zmingl = -totleng
        z1 = zmingl + zrise
        z2 = z1 + zflat
        z3 = 0.0d0
        numpts = 0
        iptlc = 0
        tmptcl = 0.0d0

        do
10        call random_number(rr)
          zz = -totleng*rr
          if(zz.ge.zmingl .and. zz.lt.z1 ) then
            fvalue = exp(-((zz-z1)/zscale)**2)
          else if( zz.ge.z1 .and. zz.lt.z2) then
            fvalue = 1.0d0
          else if( zz.ge.z2 .and. zz.le.z3) then
            fvalue = exp(-((zz-z2)/zscale)**2)
          else
          endif
          call random_number(rr2)
          if(rr2.gt.fvalue) goto 10
          numpts = numpts + 1
          if( (zz.le.zmaxbin) .and. (zz.gt.zminbin) ) then
            iptlc = iptlc + 1
            tmptcl(iptlc) = zz
          endif
          if(numpts.ge.avgpts0) exit
        enddo

        avgpts = iptlc

        do numpts = 1, avgpts
         ptz(numpts) = xmu5 + tmptcl(numpts)
        enddo

        end subroutine LongFlattopGasz_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> longitudinal spatial Gaussian distribution with cut-off. (id=3)
        !------------------------------------------------------------------------------------------------
        subroutine LongGaussz_Dist(ptz,avgpts,sigz,xmu5,cutz)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts
        real*8, intent(out), dimension(avgpts) :: ptz
        real*8, intent(in) :: sigz,xmu5,cutz
        real*8, dimension(2,avgpts) :: pt
        integer :: numpts
  
        call normVecCut(pt,avgpts,cutz,cutz)

        do numpts = 1, avgpts
          ptz(numpts) = xmu5 + sigz*pt(1,numpts)-cutz*sigz
          !ptz(numpts) = xmu5 + sigz*pt(1,numpts)
        enddo

        end subroutine LongGaussz_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> 3D spatial Waterbag distribution. (id=34)
        !------------------------------------------------------------------------------------------------
        subroutine Waterxyz_Dist(ptxyz,avgpts,sigx,sigy,sigz,xmu1,xmu3,xmu5)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts
        real*8, intent(out), dimension(3,avgpts) :: ptxyz
        real*8, intent(in) :: sigx,sigy,sigz,xmu1,xmu3,xmu5
        integer :: numpts
        real*8 :: r1,r2,r3
  
        do numpts = 1, avgpts
          ! rejection sample.
10        call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          r1 = 2.0d0*r1-1.0d0
          r2 = 2.0d0*r2-1.0d0
          r3 = 2.0d0*r3-1.0d0
          if(r1**2+r2**2+r3**2.gt.1.0d0) goto 10
          ptxyz(1,numpts) = xmu1 + sigx*r1
          ptxyz(2,numpts) = xmu3 + sigy*r2
          ptxyz(3,numpts) = xmu5 + sigy*r3
        enddo

        end subroutine Waterxyz_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> 3D Momentum Gaussian distribution. (id=1)
        !------------------------------------------------------------------------------------------------
        subroutine PxPyPzGauss3d_Dist(pmxyz,avgpts,sigpx,sigpy,sigpz,xmu2,xmu4,xmu6)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts
        real*8, intent(out), dimension(3,avgpts) :: pmxyz
        real*8, intent(in) :: sigpx,sigpy,sigpz,xmu2,xmu4,xmu6
        integer :: numpts
        real*8, dimension(3) :: yy
 
        !3D Gaussian distribution
        do numpts = 1, avgpts
          call normdv2(yy)
          pmxyz(1,numpts) = xmu2 + sigpx*yy(1)
          pmxyz(2,numpts) = xmu4 + sigpy*yy(2)
          pmxyz(3,numpts) = xmu6 + sigpz*yy(3)
        enddo
 
        end subroutine PxPyPzGauss3d_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> Momentum Gaussian distribution in transverse while longitudinal xexp(-x^2/sigz^2). (id=2)
        !> The longitudinal distribution follows Bird's book, p.129.
        !------------------------------------------------------------------------------------------------
        subroutine PxPyGaussPzSg_Dist(pmxyz,avgpts,sigpx,sigpy,sigpz,xmu2,xmu4,xmu6)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts
        real*8, intent(out), dimension(3,avgpts) :: pmxyz
        real*8, intent(in) :: sigpx,sigpy,sigpz,xmu2,xmu4,xmu6
        integer :: numpts 
        real*8, dimension(2,avgpts) :: yy
        real*8 :: epsilon,rr

        epsilon = 1.0d-18
         
        call normVec(yy,avgpts)

        !3D Gaussian distribution
        do numpts = 1, avgpts
          pmxyz(1,numpts) = xmu2 + sigpx*yy(1,numpts)
          pmxyz(2,numpts) = xmu4 + sigpy*yy(2,numpts)
          call random_number(rr)
          if(rr.eq.0.0d0) rr = epsilon
          pmxyz(3,numpts) = xmu6 + sigpz*sqrt(-log(rr))
        enddo
 
        end subroutine PxPyGaussPzSg_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> Momentum distribution following the streak camera model. (id=3)
        !> both emax and wkf are in the units of eV.
        !> emass is also in eV.
        !------------------------------------------------------------------------------------------------
        subroutine PxPyPzStreakCam_Dist(pmxyz,avgpts,emax,wkf,emass,xmu2,xmu4,xmu6)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts
        real*8, intent(out), dimension(3,avgpts) :: pmxyz
        real*8, intent(in) :: emax,wkf,xmu2,xmu4,xmu6,emass
        integer :: numpts,ii
        real*8, dimension(3) :: yy
        real*8 :: epoint,fmax,hpi,r,rr1,rr2,fvalue,e1,gam,gambet,theta,phi

          epoint = wkf/3.0d0
          fmax = epoint/(epoint+wkf)**4
          hpi = dasin(1.0d0)*4
          numpts = 0
          ii = 0
          do
            ! rejection sample
10          call random_number(r)
            rr1 = r*emax
            fvalue = rr1/(rr1+wkf)**4/fmax
            call random_number(rr2)
            if(rr2.gt.fvalue) goto 10
            e1 = rr1
            numpts = numpts + 1
            if(numpts.gt.avgpts) exit
            gam = e1/eMass + 1.0d0
            gambet =  sqrt(gam*gam - 1.0d0)
            call random_number(r)
            rr1 = 2*r-1.0d0
            theta = acos(-rr1)/2
            call random_number(r)
            rr1 = r
            phi = hpi*rr1
            ii = ii + 1
            pmxyz(1,ii) = xmu2 + gambet*sin(theta)*cos(phi)
            pmxyz(2,ii) = xmu4 + gambet*sin(theta)*sin(phi)
            pmxyz(3,ii) = xmu6 + gambet*cos(theta)
          enddo

        end subroutine PxPyPzStreakCam_Dist

        !-----------------------------------------------------------------------------------------------
        !> @brief
        !> Momentum distribution following the streak camera model. (id=4)
        !> momentum: following the 3 step model in David Dowell's PRSTAB paper
        !> v.12, 074201, 2009.
        !> momentum input parameters are:
        !> laser photon energy (eV), cathode temperature (eV),
        !> cathode effective work function (eV) (work function - surface field),
        !> Fermi energy of cathode material, emass (eV).
        !------------------------------------------------------------------------------------------------
        subroutine PxPyPz3step_Dist(pmxyz,avgpts,Ef,Eph,Tem,Ewk,emass,xmu2,xmu4,xmu6)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: avgpts
        real*8, intent(out), dimension(3,avgpts) :: pmxyz
        real*8, intent(in) :: Ef,Eph,Tem,emass,Ewk,xmu2,xmu4,xmu6
        integer :: numpts,ii,isamz
        real*8 :: Emax,Emin,fmax,twopi,r,rr,fvalue,theta,phi
        real*8 :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,vr1,vr2

        twopi = 4*asin(1.0d0)
        Emax = Ef+Eph+10*Tem
        Emin = Ef+Ewk
        if(Emin.ge.Emax) then
          print*,"The input photon enery is too small:"
          print*,"photon energy should be greater than the sum of Fermi energy and the effective work function"
          stop 
        endif
        tmp1 = Eph/Tem/2
        tmp2 = exp(-tmp1)
        fmax = (1.0d0/(1.0d0+tmp2))**2
        isamz = 0
        do  
          ! rejection sample.
10        call random_number(r) 
          vr1 = Emin+r*(Emax-Emin)
          tmp3 = (vr1-Ef)/Tem
          if(tmp3.gt.0.0d0) then
            tmp4 = exp(-tmp3)
            tmp5 = 1.0d0/(1.0d0+tmp4)
          else
            tmp4 = exp(tmp3)
            tmp5 = tmp4/(1.0d0+tmp4)
          endif
          tmp6 = (vr1-Ef-Eph)/Tem
          if(tmp6.gt.0.0d0) then
            tmp7 = exp(-tmp6)
            tmp8 = tmp7/(1.0d0+tmp7)
          else
            tmp7 = exp(tmp6)
            tmp8 = 1.0d0/(1.0d0+tmp7)
          endif
          fvalue = tmp5*tmp8/fmax
          call random_number(vr2)
          if(vr2.gt.fvalue) goto 10
          call random_number(r)
          tmp9 = vr1*r*r
          if(tmp9.lt.(Ef+Ewk)) goto 10
          isamz = isamz + 1
          if(isamz.le.avgpts) then
            theta = acos(r)
            call random_number(rr)
            phi = twopi*rr
            pmxyz(1,isamz) = xmu2+sqrt(2*vr1/emass)*sin(theta)*cos(phi)  
            pmxyz(2,isamz) = xmu4+sqrt(2*vr1/emass)*sin(theta)*sin(phi) 
            pmxyz(3,isamz) = xmu6+sqrt(2*(tmp9-Ef-Ewk)/emass)  
          else
            exit
          endif
        enddo

        end subroutine PxPyPz3step_Dist
 
        subroutine normVecCut(y,num,cutx,cuty)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: num
        double precision, dimension(2,num), intent(out) :: y
        real*8, intent(in) :: cutx,cuty
        double precision :: twopi,epsilon
        double precision, dimension(10*num) :: x1,x2
        integer :: i,ipt1,ipt2
        real*8 :: tmp1
 
        epsilon = 1.0d-18
 
        twopi = 4.0*dasin(1.0d0)
        call random_number(x2)
        call random_number(x1)

        ipt1 = 0
        do i = 1, 10*num
          if(x1(i).eq.0.0d0) x1(i) = epsilon
          tmp1 = sqrt(-2.0*log(x1(i)))
          if(abs(tmp1).le.cutx) then
            ipt1 = ipt1 + 1
            y(1,ipt1) = tmp1*cos(twopi*x2(i))
            y(2,ipt1) = tmp1*sin(twopi*x2(i))
          else
          endif
          if(ipt1.ge.num) exit
        enddo

!        ipt2 = 0
!        do i = 1, 4*num
!          if(x1(i).eq.0.0d0) x1(i) = epsilon
!          tmp1 = sqrt(-2.0*log(x1(i)))*sin(twopi*x2(i))
!          if(abs(tmp1).le.cuty) then
!            ipt2 = ipt2 + 1
!            y(2,ipt2) = tmp1
!          else
!          endif
!          if(ipt2.ge.num) exit
!        enddo

!        if(ipt1.lt.num .or. ipt2.lt.num) then
        if(ipt1.lt.num) then
          print*,"input parameter wrong in normvecCut!"
          stop
        endif
 
        end subroutine normVecCut

      end module Distributionclass
