!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! SolRFclass: Solenoid with imbeded RF field beam line element class
!             in Lattice module of APPLICATION layer.
! MODULE  : ... SolRFclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang 
! DESCRIPTION:
!> This class defines the linear transfer map and RF field
!> for the Sol-RF beam line elment.
! Comments:
!----------------------------------------------------------------
      module SolRFclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 12
        type SolRF
          !Itype = 105
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : scale
          !      (3) : RF frequency
          !      (4) : theta0
          !      (5) : file ID
          !      (6) : radius
          !      (7) : x misalignment error
          !      (8) : y misalignment error
          !      (9) : rotation error x
          !      (10) : rotation error y
          !      (11) : rotation error z
          !      (12) : Bz0
        end type SolRF
        interface getparam_SolRF
          module procedure getparam1_SolRF,  &
                          getparam2_SolRF,   &
                          getparam3_SolRF
        end interface
        interface setparam_SolRF
          module procedure setparam1_SolRF,  &
                          setparam2_SolRF, setparam3_SolRF
        end interface
      contains
        subroutine construct_SolRF(this,numseg,nmpstp,type,blength)
        implicit none
        type (SolRF), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_SolRF
   
        subroutine setparam1_SolRF(this,i,value)
        implicit none
        type (SolRF), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_SolRF

        subroutine setparam2_SolRF(this,values)
        implicit none
        type (SolRF), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_SolRF

        subroutine setparam3_SolRF(this,numseg,nmpstp,type,blength)
        implicit none
        type (SolRF), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_SolRF
   
        subroutine getparam1_SolRF(this,i,blparam) 
        implicit none 
        type (SolRF), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_SolRF
  
        subroutine getparam2_SolRF(this,blparams)
        implicit none
        type (SolRF), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_SolRF

        subroutine getparam3_SolRF(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (SolRF), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_SolRF
       
        !--------------------------------------------------------------------------------------
        !> @brief
        !> interpolate the field from the SolRF rf cavity onto bunch location.
        !--------------------------------------------------------------------------------------
        subroutine getaxfldE_SolRF(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        type (SolRF), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1
        double precision:: zz,hstep,slope,zedge,escale,zlen
        integer :: klo,khi,k
        integer :: my_rank,ierr

        zedge = this%Param(1)
        escale = this%Param(2)
        zz=z-zedge

!-------------------------------------------------------------------
! on axis field from interpolation
!        klo=1
!        khi=Ndata
!1       if(khi-klo.gt.1) then
!          k=(khi+klo)/2
!          if(zdat(k).gt.zz)then
!             khi=k
!          else
!             klo=k
!          endif
!          goto 1
!        endif
!        hstep=zdat(khi)-zdat(klo)
!        slope=(edat(khi)-edat(klo))/hstep
!        ez1 =edat(klo)+slope*(zz-zdat(klo))
!        slope=(epdat(khi)-epdat(klo))/hstep
!        ezp1=epdat(klo)+slope*(zz-zdat(klo))
!        ezpp1 = 0.0
!-------------------------------------------------------------------
! on axis field from analytical function
        zlen = this%Length
        pi = 2*asin(1.0)
        ez1 = 30.0*sin(26.72*2*pi*zz/zlen)
        ezp1 = 30.0*26.72*2*pi/zlen*cos(26.72*2*pi*zz/zlen)
        ezpp1 = 0.0

        ez1 = ez1*escale
        ezp1 = ezp1*escale
        ezpp1 = ezpp1*escale

        !print*,"ez1: ",ez1,escale,zedge,zz
!        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
!        if(my_rank.eq.1) then
!         write(19,101)z,ez1/1.0e6,ezp1/1.0e6,zz,float(m),float(init)
!        endif
!        if(my_rank.eq.1) then
!         write(19,101)z,gt
!        endif
  101   format(6(1x,1pe13.6))

        end subroutine getaxfldE_SolRF

        subroutine getBgradfld_SolRF(z,this,b0,bgrad)
        implicit none
        include 'mpif.h'
        type (SolRF), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: b0,bgrad

        !uniform bz field.
        b0 = this%Param(12)
        bgrad = 0.0

        end subroutine getBgradfld_SolRF

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getflderr_SolRF(pos,extfld,this,dx,dy,anglex,angley,&
                                    anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (SolRF), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid
        double precision, dimension(3) :: temp,tmp
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq 
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        b0 = this%Param(12)
        !move into the tank local coordinate.
!        zz=pos(3)-zedge
        !uniform bz field.

!        dx = this%Param(5)
!        dy = this%Param(6)
!        anglex = this%Param(7)
!        angley = this%Param(8)
!        anglez = this%Param(9)

        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = pos(3) - zedge
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

        bgrad = 0.0
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = -bgrad/2*tmp(1)
        extfld(5) = -bgrad/2*tmp(2)
        extfld(6) = b0

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        !move into the tank local coordinate.
!        zz=pos(3)-zmid
        zz=pos(3)-zedge
        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = zz
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

! get external RF field on axis from analytical function
!-----------------------------------------------------------------
!        zz=pos(3)-zmid
!        ez1 = Fcoef(1)/2
!        ezp1 = 0.0
!        ezpp1 = 0.0
!        do i = 2, (Ndata-1)/2+1
!          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
!                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
!          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
!                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
!          ezpp1 = ezp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
!                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
!        enddo
        ez1 = 30.0*sin(26.72*2*pi*zz/len)
        ezp1 = 30.0*26.72*2*pi/len*cos(26.72*2*pi*zz/len)
        ezpp1 = 0.0

        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
!-----------------------------------------------------------------
        tt = pos(4) 
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        tmpex = -tmp(1)*ezp1*tmpcos/2
        tmpey = -tmp(2)*ezp1*tmpcos/2
        tmpez = (ez1+f1*(tmp(1)**2+tmp(2)**2))*tmpcos
        tmpbx = tmp(2)/(xl*clite)*(ez1*tmpsin/2)
        tmpby = -tmp(1)/(xl*clite)*(ez1*tmpsin/2)
        tmpbz = 0.0

        !for E field
        tmp(1) = tmpex
        tmp(2) = tmpey*cos(anglex)-tmpez*sin(anglex)
        tmp(3) = tmpey*sin(anglex)+tmpez*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1) = extfld(1)+temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(2) = extfld(2)+temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(3) = extfld(3)+temp(3)

        !for B field
        tmp(1) = tmpbx
        tmp(2) = tmpby*cos(anglex)-tmpbz*sin(anglex)
        tmp(3) = tmpby*sin(anglex)+tmpbz*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = extfld(4)+temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = extfld(5)+temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = extfld(6)+temp(3)

        end subroutine getflderr_SolRF
        
        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getfld_SolRF(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (SolRF), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        !move into the tank local coordinate.
        zz=pos(3)-zedge
        !uniform bz field.
        b0 = this%Param(12)
        bgrad = 0.0
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = -bgrad/2*pos(1)
        extfld(5) = -bgrad/2*pos(2)
        extfld(6) = b0
        !This search is based on the assumption that the data is uniformly
        !distributed along z.
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
!        zz=pos(3)-zedge
!        ez1 = Fcoef(1)/2
!        ezp1 = 0.0
!        ezpp1 = 0.0
!        do i = 2, (Ndata-1)/2+1
!          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
!                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
!          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
!                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
!          ezpp1 = ezp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
!                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
!        enddo
        ez1 = 30.0*sin(26.72*2*pi*zz/len)
        ezp1 = 30.0*26.72*2*pi/len*cos(26.72*2*pi*zz/len)
        ezpp1 = 0.0

        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
!-----------------------------------------------------------------
        tt = pos(4)
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        tmpex = -pos(1)*ezp1*tmpcos/2
        tmpey = -pos(2)*ezp1*tmpcos/2
        f1 = -(ezpp1+ez1/xl/xl)/4
        tmpez = (ez1+f1*(pos(1)**2+pos(2)**2))*tmpcos
        tmpbx = pos(2)/(xl*clite)*(ez1*tmpsin/2)
        tmpby = -pos(1)/(xl*clite)*(ez1*tmpsin/2)
        tmpbz = 0.0
        extfld(1) = extfld(1) + tmpex
        extfld(2) = extfld(2) + tmpey
        extfld(3) = extfld(3) + tmpez
        extfld(4) = extfld(4) + tmpbx
        extfld(5) = extfld(5) + tmpby
        extfld(6) = extfld(6) + tmpbz

        end subroutine getfld_SolRF

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getfldt_SolRF(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (SolRF), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl,xlrep
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid,rr,zstart1,zend1,zlength1,&
                   zstart2,zend2,zlength2,zlen,f1p,r2,zlc,ezppp,bscale
        integer :: i,ntmp,numpar1,numpar2

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
!        print*,"zedge: ",zedge,len,pos(3)
    
        ez1 = 0.0
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        if((pos(3).gt.zedge).and.(pos(3).lt.(zedge+len))) then
          ww = this%Param(3)*2*pi
!          xl = clite/ww  !real frequency has to be used here
          xlrep = ww/clite  !real frequency has to be used here
          theta0 = this%Param(4)*asin(1.0)/90
          !move into the tank local coordinate.
          zlc=pos(3)-zedge
          numpar1 = fldata%Fcoeft(1)+0.1
          !//Here, zstart1 is the starting RF field location with 
          !//respect to the "zedge " as origin.
          !zstart1 can be negative
          zstart1 = fldata%Fcoeft(2)
          zend1 = fldata%Fcoeft(3)
          zlength1 = fldata%Fcoeft(4)
          extfld = 0.0
!          print*,"zstart1: ",zstart1,zend1,zlc,numpar1
          if( (zlc.ge.zstart1).and.(zlc.le.zend1)) then
            zmid = zlength1/2
            zz = pos(3)-zedge-zstart1-zmid
            zlen = zlength1
            !first 4 parameters in Fcoeft are not Forier coefficients.
            ntmp = 4
            !//find the field on the axis and its 1rst, 2nd, and 3rd
            !//derivatives which will be used to calculate the off-axis
            !//field distribution. 
            ez1 = fldata%Fcoeft(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar1-1)/2+1
              ez1 = ez1 + fldata%Fcoeft(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
!               print*,fldata%Fcoeft(2*(i-ntmp)-2)
!               print*,fldata%Fcoeft(2*(i-ntmp)-1)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-fldata%Fcoeft(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                fldata%Fcoeft(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-fldata%Fcoeft(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
              ezppp = ezppp+((i-1)*2*pi/zlen)**3*&
                (fldata%Fcoeft(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
            enddo
            ez1=ez1*escale
            ezp1=ezp1*escale
            ezpp1=ezpp1*escale
            ezppp = ezppp*escale
            tt = pos(4)
            tmpcos = cos(ww*tt+theta0)
            tmpsin = sin(ww*tt+theta0)
!            f1 = -(ezpp1+ez1/xl/xl)/4
!            f1p = -(ezppp+ezp1/xl/xl)/4
            f1 = -(ezpp1+ez1*xlrep*xlrep)/4
            f1p = -(ezppp+ezp1*xlrep*xlrep)/4
!            f1 = 0.0
!            f1p = 0.0
            r2 = pos(1)**2+pos(2)**2
            tmpex = -pos(1)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpey = -pos(2)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpez = (ez1+f1*r2)*tmpcos
            tmpbx = pos(2)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpby = -pos(1)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
!            tmpex = -pos(1)*(ezp1/2+f1p*r2/4)
!            tmpey = -pos(2)*(ezp1/2+f1p*r2/4)
!            tmpez = (ez1+f1*r2)
!            tmpbx = pos(2)/(xl*clite)*(ez1/2+f1*r2/4)
!            tmpby = -pos(1)/(xl*clite)*(ez1/2+f1*r2/4)
            tmpbz = 0.0
            extfld(1) = extfld(1) + tmpex
            extfld(2) = extfld(2) + tmpey
            extfld(3) = extfld(3) + tmpez
            extfld(4) = extfld(4) + tmpbx
            extfld(5) = extfld(5) + tmpby
            extfld(6) = extfld(6) + tmpbz
          else
            extfld = 0.0
          endif

!          print*,"strange",zlc,ez1,ezp1,ezpp1,tmpez,tmpcos,tmpsin,ww,tt,theta0
!          write(10,101)pos(3),zlc,ez1/escale,ezp1,ezpp1,tmpez,tmpcos,tmpsin,tt,theta0
!101       format(10(1x,e14.5))
!          call flush_(10)

          !//# of parameters for solenoid B fields.
          numpar2 = fldata%Fcoeft(5+numpar1) + 0.1
          !zstart2 can be negative
          zstart2 = fldata%Fcoeft(6+numpar1) 
          zend2 = fldata%Fcoeft(7+numpar1) 
          zlength2 = fldata%Fcoeft(8+numpar1) 
          bscale = this%Param(12)
          if( (zlc.ge.zstart2) .and. (zlc.le.zend2)) then
            zmid = zlength2/2
            zz = pos(3) - zedge - zstart2 - zmid
            zlen = zlength2
            ntmp = numpar1+8
            !//find the B field on axis and its 1st,2nd,3rd derivatives
            !//which will be used to calculate the offaxis B field
            ez1 = fldata%Fcoeft(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar2-1)/2+1
              ez1 = ez1 + fldata%Fcoeft(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-fldata%Fcoeft(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                fldata%Fcoeft(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-fldata%Fcoeft(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
              ezppp = ezppp+((i-1)*2*pi/zlen)**3*&
                (fldata%Fcoeft(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
            enddo
            r2 = pos(1)*pos(1)+pos(2)*pos(2)
            extfld(4) = extfld(4)-bscale*ezp1/2*pos(1)+bscale*ezppp*pos(1)*r2/16  
            extfld(5) = extfld(5)-bscale*ezp1/2*pos(2)+bscale*ezppp*pos(2)*r2/16  
            extfld(6) = extfld(6)+bscale*ez1-bscale*ezpp1*r2/4
!            write(11,102)pos(3),zlc,ez1*bscale,ezp1*bscale,ezpp1,ezppp
!!            write(11,102)pos(3),zlc,ez1,ezp1*bscale,ezpp1,ezppp
!102         format(6(1x,e14.5))
!            call flush_(11)
          else
!            extfld = 0.0
          endif
!          print*,zlc,extfld
        else
          extfld = 0.0
        endif

        end subroutine getfldt_SolRF

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getflderrt_SolRF(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (SolRF), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl,xlrep
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid,rr,zstart1,zend1,zlength1,&
                   zstart2,zend2,zlength2,zlen,f1p,r2,zlc,ezppp,bscale
        integer :: i,ntmp,numpar1,numpar2
        double precision :: dx,dy,anglex,angley,anglez
        double precision, dimension(3) :: temp,tmp

        clite = 299792458.e0
        pi = 2*asin(1.0)

        len = this%Length
        zedge = this%Param(1)
        escale = this%Param(2)
        dx = this%Param(7)
        dy = this%Param(8)
        anglex = this%Param(9)
        angley = this%Param(10)
        anglez = this%Param(11)
!        print*,"zedge: ",zedge,len,pos(3)

        ez1 = 0.0
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        if((pos(3).gt.zedge).and.(pos(3).lt.(zedge+len))) then

          ww = this%Param(3)*2*pi
!          xl = clite/ww  !real frequency has to be used here
          xlrep = ww/clite  !real frequency has to be used here
          theta0 = this%Param(4)*asin(1.0)/90
          !move into the tank local coordinate.
          !zlc=pos(3)-zedge
          temp(1) = pos(1) - dx
          temp(2) = pos(2) - dy
          tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
          tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
          tmp(3) = pos(3) - zedge
          temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
          tmp(1) = temp(1)
          tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
          tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
          zlc = tmp(3)

          numpar1 = fldata%Fcoeft(1)+0.1
          !//Here, zstart1 is the starting RF field location with 
          !//respect to the "zedge " as origin.
          !zstart1 can be negative
          zstart1 = fldata%Fcoeft(2)
          zend1 = fldata%Fcoeft(3)
          zlength1 = fldata%Fcoeft(4)
          extfld = 0.0
!          print*,"zstart1: ",zstart1,zend1,zlc,numpar1
          if( (zlc.ge.zstart1).and.(zlc.le.zend1)) then
            zmid = zlength1/2
            !zz = pos(3)-zedge-zstart1-zmid
            zz = zlc-zstart1-zmid
            zlen = zlength1
            !first 4 parameters in Fcoeft are not Forier coefficients.
            ntmp = 4
            !//find the field on the axis and its 1rst, 2nd, and 3rd
            !//derivatives which will be used to calculate the off-axis
            !//field distribution. 
            ez1 = fldata%Fcoeft(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar1-1)/2+1
              ez1 = ez1 + fldata%Fcoeft(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
!               print*,fldata%Fcoeft(2*(i-ntmp)-2)
!               print*,fldata%Fcoeft(2*(i-ntmp)-1)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-fldata%Fcoeft(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                fldata%Fcoeft(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-fldata%Fcoeft(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
              ezppp = ezppp+((i-1)*2*pi/zlen)**3*&
                (fldata%Fcoeft(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
            enddo
            ez1=ez1*escale
            ezp1=ezp1*escale
            ezpp1=ezpp1*escale
            ezppp = ezppp*escale
            tt = pos(4)
            tmpcos = cos(ww*tt+theta0)
            tmpsin = sin(ww*tt+theta0)
!            f1 = -(ezpp1+ez1/xl/xl)/4
!            f1p = -(ezppp+ezp1/xl/xl)/4
            f1 = -(ezpp1+ez1*xlrep*xlrep)/4
            f1p = -(ezppp+ezp1*xlrep*xlrep)/4
!            f1 = 0.0
!            f1p = 0.0
            !r2 = pos(1)**2+pos(2)**2
            r2 = tmp(1)**2+tmp(2)**2
            tmpex = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpey = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpez = (ez1+f1*r2)*tmpcos
            tmpbx = tmp(2)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpby = -tmp(1)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
!            tmpex = -pos(1)*(ezp1/2+f1p*r2/4)
!            tmpey = -pos(2)*(ezp1/2+f1p*r2/4)
!            tmpez = (ez1+f1*r2)
!            tmpbx = pos(2)/(xl*clite)*(ez1/2+f1*r2/4)
!            tmpby = -pos(1)/(xl*clite)*(ez1/2+f1*r2/4)
            tmpbz = 0.0
            extfld(1) = extfld(1) + tmpex
            extfld(2) = extfld(2) + tmpey
            extfld(3) = extfld(3) + tmpez
            extfld(4) = extfld(4) + tmpbx
            extfld(5) = extfld(5) + tmpby
            extfld(6) = extfld(6) + tmpbz
          else
            extfld = 0.0
          endif

!          print*, zlc,ez1,ezp1,ezpp1,tmpez,tmpcos,tmpsin,ww,tt,theta0
!          write(10,101)pos(3),zlc,ez1,ezp1,ezpp1,tmpez,tmpcos,tmpsin,tt,theta0
!101       format(10(1x,e14.5))
!          call flush(10)

          !//# of parameters for solenoid B fields.
          numpar2 = fldata%Fcoeft(5+numpar1) + 0.1
          !zstart2 can be negative
          zstart2 = fldata%Fcoeft(6+numpar1) 
          zend2 = fldata%Fcoeft(7+numpar1) 
          zlength2 = fldata%Fcoeft(8+numpar1) 
          bscale = this%Param(12)
          if( (zlc.ge.zstart2) .and. (zlc.le.zend2)) then
            zmid = zlength2/2
            !zz = pos(3) - zedge - zstart2 - zmid
            zz = zlc - zstart2 - zmid
            zlen = zlength2
            ntmp = numpar1+8
            !//find the B field on axis and its 1st,2nd,3rd derivatives
            !//which will be used to calculate the offaxis B field
            ez1 = fldata%Fcoeft(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar2-1)/2+1
              ez1 = ez1 + fldata%Fcoeft(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-fldata%Fcoeft(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                fldata%Fcoeft(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-fldata%Fcoeft(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
              ezppp = ezppp+((i-1)*2*pi/zlen)**3*&
                (fldata%Fcoeft(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
            enddo
            r2 = tmp(1)*tmp(1)+tmp(2)*tmp(2)
            extfld(4) = extfld(4)-bscale*ezp1/2*tmp(1)+bscale*ezppp*tmp(1)*r2/16  
            extfld(5) = extfld(5)-bscale*ezp1/2*tmp(2)+bscale*ezppp*tmp(2)*r2/16  
            extfld(6) = extfld(6)+bscale*ez1-bscale*ezpp1*r2/4
!            write(11,102)pos(3),zlc,ez1*bscale,ezp1*bscale,ezpp1,ezppp
!102         format(6(1x,e14.5))
!            call flush_(11)
          else
!            extfld = 0.0
          endif
!          print*,zlc,extfld
          !for E field
          tmp(1) = extfld(1)
          tmp(2) = extfld(2)*cos(anglex)-extfld(3)*sin(anglex)
          tmp(3) = extfld(2)*sin(anglex)+extfld(3)*cos(anglex)
          temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
          extfld(1) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
          extfld(2) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
          extfld(3) = temp(3)
 
          !for B field
          tmp(1) = extfld(4)
          tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
          tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
          temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
          extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
          extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
          extfld(6) = temp(3)
        else
          extfld = 0.0
        endif

        end subroutine getflderrt_SolRF

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getvecAt_SolRF(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (SolRF), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(3), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl,xlrep
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid,rr,zstart1,zend1,zlength1,&
                   zstart2,zend2,zlength2,zlen,f1p,r2,zlc,ezppp,bscale
        integer :: i,ntmp,numpar1,numpar2

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
!        print*,"zedge: ",zedge,len,pos(3)
    
        ez1 = 0.0
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        if((pos(3).gt.zedge).and.(pos(3).lt.(zedge+len))) then
          ww = this%Param(3)*2*pi
!          xl = clite/ww  !real frequency has to be used here
          xlrep = ww/clite  !real frequency has to be used here
          theta0 = this%Param(4)*asin(1.0)/90
          !move into the tank local coordinate.
          zlc=pos(3)-zedge
          numpar1 = fldata%Fcoeft(1)+0.1
          !//Here, zstart1 is the starting RF field location with 
          !//respect to the "zedge " as origin.
          !zstart1 can be negative
          zstart1 = fldata%Fcoeft(2)
          zend1 = fldata%Fcoeft(3)
          zlength1 = fldata%Fcoeft(4)
          extfld = 0.0
!          print*,"zstart1: ",zstart1,zend1,zlc,numpar1

          !//# of parameters for solenoid B fields.
          numpar2 = fldata%Fcoeft(5+numpar1) + 0.1
          !zstart2 can be negative
          zstart2 = fldata%Fcoeft(6+numpar1) 
          zend2 = fldata%Fcoeft(7+numpar1) 
          zlength2 = fldata%Fcoeft(8+numpar1) 
          bscale = this%Param(12)
          if( (zlc.ge.zstart2) .and. (zlc.le.zend2)) then
            zmid = zlength2/2
            zz = pos(3) - zedge - zstart2 - zmid
            zlen = zlength2
            ntmp = numpar1+8
            !//find the B field on axis and its 1st,2nd,3rd derivatives
            !//which will be used to calculate the offaxis B field
            ez1 = fldata%Fcoeft(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar2-1)/2+1
              ez1 = ez1 + fldata%Fcoeft(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-fldata%Fcoeft(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - fldata%Fcoeft(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
            enddo
            r2 = pos(1)*pos(1)+pos(2)*pos(2)
            extfld(1) = -bscale*ez1/2*pos(2)+bscale*ezpp1*pos(2)*r2/16  
            extfld(2) = bscale*ez1/2*pos(1)-bscale*ezpp1*pos(1)*r2/16  
            extfld(3) = 0.0d0
          else
!            extfld = 0.0
          endif
!          print*,zlc,extfld
        else
          extfld = 0.0
        endif

        end subroutine getvecAt_SolRF

      end module SolRFclass
