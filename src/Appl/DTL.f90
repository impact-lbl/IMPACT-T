!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! DTLclass: Drift-tube-linac beam line element class
!             in Lattice module of APPLICATION layer.
! MODULE  : ... DTLclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang
! DESCRIPTION: 
!> This class defines the linear transfer map and RF field for the DTL beam line elment.
!              
! Comments:
!----------------------------------------------------------------
      module DTLclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 25
        type DTL
          !Itype = 101
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : scale
          !      (3) : RF frequency
          !      (4) : theta0
          !      (5) : file ID
          !      (6) : radius
          !      (7) : quad 1 length
          !      (8) : quad 1 gradient
          !      (9) : quad 2 length
          !      (10) : quad 2 gradient
          !      (11) : x misalignment error for Quad 1
          !      (12) : y misalignment error for Quad 1
          !      (13) : rotation error x for Quad 1
          !      (14) : rotation error y for Quad 1
          !      (15) : rotation error z for Quad 1
          !      (16) : x misalignment error for Quad 2
          !      (17) : x misalignment error for Quad 2
          !      (18) : rotation error x for Quad 2
          !      (19) : rotation error y for Quad 2
          !      (20) : rotation error z for Quad 2
          !      (21) : x misalignment error for RF cavity
          !      (22) : y misalignment error for RF cavity
          !      (23) : rotation error x for RF cavity
          !      (24) : rotation error y for RF cavity
          !      (25) : rotation error z for RF cavity
        end type DTL
        interface getparam_DTL
          module procedure getparam1_DTL,  &
                          getparam2_DTL,   &
                          getparam3_DTL
        end interface
        interface setparam_DTL
          module procedure setparam1_DTL,  &
                          setparam2_DTL, setparam3_DTL
        end interface
      contains
        subroutine construct_DTL(this,numseg,nmpstp,type,blength)
        implicit none
        type (DTL), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_DTL
   
        subroutine setparam1_DTL(this,i,value)
        implicit none
        type (DTL), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_DTL

        subroutine setparam2_DTL(this,values)
        implicit none
        type (DTL), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_DTL

        subroutine setparam3_DTL(this,numseg,nmpstp,type,blength)
        implicit none
        type (DTL), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_DTL
   
        subroutine getparam1_DTL(this,i,blparam) 
        implicit none 
        type (DTL), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_DTL
  
        subroutine getparam2_DTL(this,blparams)
        implicit none
        type (DTL), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_DTL

        subroutine getparam3_DTL(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (DTL), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_DTL
       
        !--------------------------------------------------------------------------------------
        !> @brief
        !> interpolate the field from the DTL rf cavity onto bunch location.
        !--------------------------------------------------------------------------------------
        subroutine getaxfldE_DTL(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        type (DTL), intent(in) :: this
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
        ezpp1 = 0.0
        ez1 = ez1*escale
        ezp1 = ezp1*escale

        !print*,"ez1: ",ez1,escale,zedge,zz
!        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
!        if(my_rank.eq.1) then
!         write(19,101)z,ez1/1.0e6,ezp1/1.0e6,zz,float(m),float(init)
!        endif
!        if(my_rank.eq.1) then
!         write(19,101)z,gt
!        endif
  101   format(6(1x,1pe13.6))

        end subroutine getaxfldE_DTL

        subroutine getBgradfld_DTL(z,this,bgrad)
        implicit none
        include 'mpif.h'
        type (DTL), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: bgrad
        double precision:: zz,bgrad1,bgrad2,len1,len2,len,zedge

        len = this%Length
        zedge = this%Param(1)
        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)
        zz=z-zedge
        if(zz.lt.len1) then
          bgrad = bgrad1
        else if(zz.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif

        end subroutine getBgradfld_DTL

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getflderrold_DTL(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (DTL), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision :: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
        double precision :: dx1,dy1,anglex1,angley1,anglez1
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
        !move into the tank local coordinate.
        !zz=pos(3)-zedge

        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)

!        dx = this%Param(11)
!        dy = this%Param(12)
!        anglex = this%Param(13)
!        angley = this%Param(14)
!        anglez = this%Param(15)
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

        if(zz.lt.len1) then
          bgrad = bgrad1
        else if(zz.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = bgrad*tmp(2)
        extfld(5) = bgrad*tmp(1)
        extfld(6) = 0.0

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        dx1 = this%Param(16)
        dy1 = this%Param(17)
        anglex1 = this%Param(18)
        angley1 = this%Param(19)
        anglez1 = this%Param(20)

        zz=pos(3)-zmid
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

        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (Ndata-1)/2 + 1
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
        tmpex = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = tmp(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -tmp(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
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

        end subroutine getflderrold_DTL
        
        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with different field and Quad displacement and rotation errors.
        !> Here, there 2 displacement errors and rotation errors for Quad 
        !--------------------------------------------------------------------------------------
        subroutine  getflderr_DTL(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (DTL), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision :: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
        double precision :: dx1,dy1,anglex1,angley1,anglez1
        double precision :: dx0,dy0,anglex0,angley0,anglez0,zz0
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
        !move into the tank local coordinate.
        zz0=pos(3)-zedge

        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)

        if(zz0.lt.len1) then
          bgrad = bgrad1
        else if(zz0.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif

        if(zz0.lt.len1) then
          bgrad = bgrad1
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

          extfld(1) = 0.0
          extfld(2) = 0.0
          extfld(3) = 0.0
          extfld(4) = bgrad*tmp(2)
          extfld(5) = bgrad*tmp(1)
          extfld(6) = 0.0

          tmp(1) = extfld(4)
          tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
          tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
          temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
          extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
          extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
          extfld(6) = temp(3)
        else if(zz0.lt.len2) then
          bgrad = 0.0
          extfld = 0.0
        else
          dx0 = this%Param(16)
          dy0 = this%Param(17)
          anglex0 = this%Param(18)
          angley0 = this%Param(19)
          anglez0 = this%Param(20)
          temp(1) = pos(1) - dx0
          temp(2) = pos(2) - dy0
          tmp(1) = temp(1)*cos(anglez0) + temp(2)*sin(anglez0)
          tmp(2) = -temp(1)*sin(anglez0) + temp(2)*cos(anglez0)
          tmp(3) = pos(3) - zedge
          temp(1) = tmp(1)*cos(angley0)+tmp(3)*sin(angley0)
          temp(2) = tmp(2)
          temp(3) = -tmp(1)*sin(angley0)+tmp(3)*cos(angley0)
          tmp(1) = temp(1)
          tmp(2) = temp(2)*cos(anglex0)+temp(3)*sin(anglex0)
          tmp(3) = -temp(2)*sin(anglex0)+temp(3)*cos(anglex0)
          zz = tmp(3)

          extfld(1) = 0.0
          extfld(2) = 0.0
          extfld(3) = 0.0
          extfld(4) = bgrad*tmp(2)
          extfld(5) = bgrad*tmp(1)
          extfld(6) = 0.0

          tmp(1) = extfld(4)
          tmp(2) = extfld(5)*cos(anglex0)-extfld(6)*sin(anglex0)
          tmp(3) = extfld(5)*sin(anglex0)+extfld(6)*cos(anglex0)
          temp(1) = tmp(1)*cos(angley0)-tmp(3)*sin(angley0)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley0)+tmp(3)*cos(angley0)
          extfld(4) = temp(1)*cos(anglez0) - temp(2)*sin(anglez0)
          extfld(5) = temp(1)*sin(anglez0) + temp(2)*cos(anglez0)
          extfld(6) = temp(3)
        endif
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        dx1 = this%Param(21)
        dy1 = this%Param(22)
        anglex1 = this%Param(23)
        angley1 = this%Param(24)
        anglez1 = this%Param(25)

        zz=pos(3)-zmid
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

        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (Ndata-1)/2 + 1
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
        tmpex = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = tmp(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -tmp(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
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

        end subroutine getflderr_DTL

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getfld_DTL(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (DTL), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
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
        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)
        if(zz.lt.len1) then
          bgrad = bgrad1
        else if(zz.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = bgrad*pos(2)
        extfld(5) = bgrad*pos(1)
        extfld(6) = 0.0
        !This search is based on the assumption that the data is uniformly
        !distributed along z.
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        zz=pos(3)-zmid
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (Ndata-1)/2 + 1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
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

        end subroutine getfld_DTL

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getfldt_DTL(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (DTL), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl,xlrep
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)*2*pi
!        xl = clite/ww  !real frequency has to be used here
        xlrep = ww/clite  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        !move into the tank local coordinate.
        zz=pos(3)-zedge
        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)
        if(zz.lt.len1) then
          bgrad = bgrad1
        else if(zz.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = bgrad*pos(2)
        extfld(5) = bgrad*pos(1)
        extfld(6) = 0.0
        !This search is based on the assumption that the data is uniformly
        !distributed along z.
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        zz=pos(3)-zmid
        ez1 = fldata%Fcoeft(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (fldata%Ndatat-1)/2 + 1
          ez1 = ez1 + fldata%Fcoeft(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                fldata%Fcoeft(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-fldata%Fcoeft(2*i-2)*&
                 sin((i-1)*2*pi*zz/len)+&
                fldata%Fcoeft(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*&
                  (-fldata%Fcoeft(2*i-2)*cos((i-1)*2*pi*zz/len)&
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
        f1 = -(ezpp1+ez1*xlrep*xlrep)/4
        f1p = -(ezppp+ezp1*xlrep*xlrep)/4
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

        end subroutine getfldt_DTL
      end module DTLclass
