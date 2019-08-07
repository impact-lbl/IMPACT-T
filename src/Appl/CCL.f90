!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! CCLclass: Coupled-cavity-linac beam line element class
!             in Lattice module of APPLICATION layer.
! 
! MODULE        : ... CCLclass
! VERSION       : ... 1.0
!> @author
!> Ji Qiang 
! DESCRIPTION: 
!> This class defines the linear transfer map and RF field for the CCL beam line elment.
! Comments:
!----------------------------------------------------------------
      module CCLclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 11
        type CCL
          !Itype = 103
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
        end type CCL
        interface getparam_CCL
          module procedure getparam1_CCL,  &
                          getparam2_CCL,   &
                          getparam3_CCL
        end interface
        interface setparam_CCL
          module procedure setparam1_CCL,  &
                           setparam2_CCL, setparam3_CCL
        end interface
      contains
        subroutine construct_CCL(this,numseg,nmpstp,type,blength)
        implicit none
        type (CCL), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_CCL
   
        subroutine setparam1_CCL(this,i,value)
        implicit none
        type (CCL), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_CCL

        subroutine setparam2_CCL(this,values)
        implicit none
        type (CCL), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_CCL

        subroutine setparam3_CCL(this,numseg,nmpstp,type,blength)
        implicit none
        type (CCL), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_CCL

        subroutine getparam1_CCL(this,i,blparam) 
        implicit none 
        type (CCL), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_CCL
  
        subroutine getparam2_CCL(this,blparams)
        implicit none
        type (CCL), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_CCL

        subroutine getparam3_CCL(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (CCL), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_CCL
       
        !--------------------------------------------------------------------------------------
        !> @brief
        !> interpolate the field from the CCL rf cavity onto bunch location.
        !--------------------------------------------------------------------------------------
        subroutine getaxfldE_CCL(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        type (CCL), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1
        double precision:: zz,hstep,slope,zedge,escale
        integer :: klo,khi,k
        integer :: my_rank,ierr

        zedge = this%Param(1)
        escale = this%Param(2)
        zz=z-zedge
        klo=1
        khi=Ndata
1       if(khi-klo.gt.1) then
          k=(khi+klo)/2
          if(zdat(k).gt.zz)then
             khi=k
          else
             klo=k
          endif
          goto 1
        endif
        hstep=zdat(khi)-zdat(klo)
        slope=(edat(khi)-edat(klo))/hstep
        ez1 =edat(klo)+slope*(zz-zdat(klo))
        slope=(epdat(khi)-epdat(klo))/hstep
        ezp1=epdat(klo)+slope*(zz-zdat(klo))
        slope=(eppdat(khi)-eppdat(klo))/hstep
        ezpp1=eppdat(klo)+slope*(zz-zdat(klo))
        ez1 = ez1*escale
        ezp1 = ezp1*escale
        ezpp1 = ezpp1*escale

!        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
!        if(my_rank.eq.1) then
!         write(19,101)z,ez1/1.0e6,ezp1/1.0e6,zz,float(m),float(init)
!        endif
!        if(my_rank.eq.1) then
!         write(19,101)z,gt
!        endif
  101   format(6(1x,1pe13.6))

        end subroutine getaxfldE_CCL

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external RF field on axis from analytical Fourier coefficients
        !--------------------------------------------------------------------------------------
        subroutine  getaxfldEfc_CCL(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (CCL), intent(in) :: this
        double precision, intent(out) :: ez1,ezp1,ezpp1
        double precision :: zedge,escale,len,zz,zmid
        double precision :: pi
        integer :: i

        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2.0

        !move into the tank local coordinate.
        !zz=z-zedge
        zz=z-zmid
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        do i = 2, (Ndata-1)/2+1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))            
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))            
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale

        end subroutine getaxfldEfc_CCL
        
        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getflderr_CCL(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (CCL), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: r2
        double precision :: zmid
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

!        dx = this%Param(7)
!        dy = this%Param(8)
!        anglex = this%Param(9)
!        angley = this%Param(10)
!        anglez = this%Param(11)

        !move into the tank local coordinate.
        !zz=pos(3)-zmid
        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = pos(3) - zmid
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        !Ndata has to be odd number
        do i = 2, (Ndata-1)/2+1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
        tt = pos(4) 
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
        
        r2 = tmp(1)**2+tmp(2)**2
        extfld(1) = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
        extfld(2) = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
        extfld(3) = (ez1+f1*r2)*tmpcos
        extfld(4) = tmp(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        extfld(5) = -tmp(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        extfld(6) = 0.0

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

        end subroutine getflderr_CCL
       
        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getfld_CCL(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (CCL), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz,r2
        double precision :: c5,c6,s5,s6,zmid
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        extfld = 0.0
        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        !move into the tank local coordinate.
        !zz=pos(3)-zedge
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
!        c5 = 0.0
!        c6 = 0.0
!        s5 = Fcoef(11)
!        s6 = 0.0
!        !print*,"s5: ",s5
!        ez1 = c5*cos(4*2*pi*zz/len)+s5*sin(4*2*pi*zz/len)+&
!            c6*cos(5*2*pi*zz/len)+s6*sin(5*2*pi*zz/len)
!        ezp1 = 8*pi/len*(-c5*sin(4*2*pi*zz/len)+s5*cos(4*2*pi*zz/len))+&
!            10*pi/len*(-c6*sin(5*2*pi*zz/len)+s6*cos(5*2*pi*zz/len))
!        ezpp1=(8*pi/len)**2*(-c5*cos(4*2*pi*zz/len)-s5*sin(4*2*pi*zz/len))+&
!            (10*pi/len)**2*(-c6*cos(5*2*pi*zz/len)-s6*sin(5*2*pi*zz/len))
        zz=pos(3)-zmid
        !ez1 = Fcoef(1)
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        !Ndata has to be odd number
        do i = 2, (Ndata-1)/2+1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
        tt = pos(4)
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
!        f1 = 0.0
!        f1p = 0.0
        r2 = pos(1)**2+pos(2)**2
        tmpex = -pos(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -pos(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = pos(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -pos(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpbz = 0.0
        extfld(1) = extfld(1) + tmpex
        extfld(2) = extfld(2) + tmpey
        extfld(3) = extfld(3) + tmpez
        extfld(4) = extfld(4) + tmpbx
        extfld(5) = extfld(5) + tmpby
        extfld(6) = extfld(6) + tmpbz

        end subroutine getfld_CCL

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getfldt_CCL(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (CCL), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz,r2
        double precision :: c5,c6,s5,s6,zmid,xlrep
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        extfld = 0.0
        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)*2*pi
!        xl = clite/ww  !real frequency has to be used here
        xlrep = ww/clite  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        !move into the tank local coordinate.
        !zz=pos(3)-zedge
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
!        c5 = 0.0
!        c6 = 0.0
!        s5 = Fcoef(11)
!        s6 = 0.0
!        !print*,"s5: ",s5
!        ez1 = c5*cos(4*2*pi*zz/len)+s5*sin(4*2*pi*zz/len)+&
!            c6*cos(5*2*pi*zz/len)+s6*sin(5*2*pi*zz/len)
!        ezp1 = 8*pi/len*(-c5*sin(4*2*pi*zz/len)+s5*cos(4*2*pi*zz/len))+&
!            10*pi/len*(-c6*sin(5*2*pi*zz/len)+s6*cos(5*2*pi*zz/len))
!        ezpp1=(8*pi/len)**2*(-c5*cos(4*2*pi*zz/len)-s5*sin(4*2*pi*zz/len))+&
!            (10*pi/len)**2*(-c6*cos(5*2*pi*zz/len)-s6*sin(5*2*pi*zz/len))
        zz=pos(3)-zmid
        !ez1 = Fcoef(1)
        ez1 = fldata%Fcoeft(1)/2
!        print*,"Fcoeft1: ",fldata%Fcoeft(1)
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        !Ndatat has to be odd number
        do i = 2, (fldata%Ndatat-1)/2+1
!          print*,"Fcoeft: ",2*i-2,fldata%Fcoeft(2*i-2)
!          print*,"Fcoeft: ",2*i-1,fldata%Fcoeft(2*i-1)
          ez1 = ez1 + fldata%Fcoeft(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                fldata%Fcoeft(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*&
                (-fldata%Fcoeft(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                fldata%Fcoeft(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*&
                  (-fldata%Fcoeft(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - fldata%Fcoeft(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*&
                  (fldata%Fcoeft(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - fldata%Fcoeft(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
        tt = pos(4)
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
!        f1 = -(ezpp1+ez1/xl/xl)/4
!        f1p = -(ezppp+ezp1/xl/xl)/4
        f1 = -(ezpp1+ez1*xlrep*xlrep)/4
        f1p = -(ezppp+ezp1*xlrep*xlrep)/4
!        f1 = 0.0
!        f1p = 0.0
        r2 = pos(1)**2+pos(2)**2
        tmpex = -pos(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -pos(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = pos(2)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -pos(1)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
        tmpbz = 0.0
        extfld(1) = extfld(1) + tmpex
        extfld(2) = extfld(2) + tmpey
        extfld(3) = extfld(3) + tmpez
        extfld(4) = extfld(4) + tmpbx
        extfld(5) = extfld(5) + tmpby
        extfld(6) = extfld(6) + tmpbz

!        write(12,101)pos(3),zz,ez1/escale,ezp1/escale,ezpp1/escale
!101     format(5(1x,e14.5))
!        call flush_(12)
        end subroutine getfldt_CCL
      end module CCLclass
