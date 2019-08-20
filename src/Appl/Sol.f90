!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! Solclass: Solenoid beam line element class
!             in Lattice module of APPLICATION layer.
! MODULE  : ... Solclass
! VERSION : ... 1.0
!> @authors
!> Ji Qiang
! DESCRIPTION:
!> This class defines the linear transfer map and field
!> for the Solenoid beam line elment.
! Comments:
!----------------------------------------------------------------
      module Solclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 9
        type Sol
          !Itype = 3
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : Bz0
          !      (3) : file ID
          !      (4) : radius
          !      (5) : x misalignment error
          !      (6) : y misalignment error
          !      (7) : rotation error x
          !      (8) : rotation error y
          !      (9) : rotation error z
        end type Sol
        interface getparam_Sol
          module procedure getparam1_Sol,  &
                          getparam2_Sol,   &
                          getparam3_Sol
        end interface
        interface setparam_Sol
          module procedure setparam1_Sol,  &
                          setparam2_Sol, setparam3_Sol
        end interface
      contains
        subroutine construct_Sol(this,numseg,nmpstp,type,blength)
        implicit none
        type (Sol), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Sol
   
        subroutine setparam1_Sol(this,i,value)
        implicit none
        type (Sol), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Sol

        subroutine setparam2_Sol(this,values)
        implicit none
        type (Sol), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_Sol

        subroutine setparam3_Sol(this,numseg,nmpstp,type,blength)
        implicit none
        type (Sol), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Sol
   
        subroutine getparam1_Sol(this,i,blparam) 
        implicit none 
        type (Sol), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Sol
  
        subroutine getparam2_Sol(this,blparams)
        implicit none
        type (Sol), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Sol

        subroutine getparam3_Sol(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Sol), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Sol
       

        subroutine getBgradfld_Sol(z,this,b0,bgrad)
        implicit none
        include 'mpif.h'
        type (Sol), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: b0,bgrad

        !uniform bz field.
        b0 = this%Param(2)
        bgrad = 0.0

        end subroutine getBgradfld_Sol

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getflderr_Sol(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Sol), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge
        double precision:: zz,bgrad,b0
        double precision, dimension(3) :: temp,tmp

        zedge = this%Param(1)
        b0 = this%Param(2)
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

        end subroutine getflderr_Sol
        
        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors and
        !> with fringe field of Solenoid. (f(z) = b0 + bb*z^2 + cc*z^4,f(x0)=0,f'(x0)=0)
        !--------------------------------------------------------------------------------------
        subroutine  getfld_Sol(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Sol), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge
        double precision:: zz,bgrad,b0,bb,cc,zshift,x0

        zedge = this%Param(1)
        x0 = 2*this%Param(4)
        !move into the tank local coordinate.
        zz=pos(3)-zedge
        !uniform bz field.
        b0 = this%Param(2)
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0

        if((zz.ge.0.0).and.(zz.lt.x0)) then
          bb = -2*b0/(x0*x0)
          cc = b0/(x0*x0*x0*x0)
          bgrad = 2*bb*(zz-x0)+4*cc*(zz-x0)*(zz-x0)*(zz-x0)
          extfld(4) = -bgrad/2*pos(1)
          extfld(5) = -bgrad/2*pos(2)
          extfld(6) = b0 + bb*(zz-x0)*(zz-x0) + cc*(zz-x0)* &
                              (zz-x0)*(zz-x0)*(zz-x0)
        else if((zz.ge.x0).and.(zz.lt.this%Length-x0)) then
          extfld(4) = 0.0
          extfld(5) = 0.0
          extfld(6) = b0
        else if((zz.ge.this%Length-x0).and.(zz.le.this%Length)) then
          bb = -2*b0/(x0*x0)
          cc = b0/(x0*x0*x0*x0)
          zshift = this%Length - x0
          bgrad = 2*bb*(zz-zshift)+4*cc*(zz-zshift)*(zz-zshift)*(zz-zshift)
          extfld(4) = -bgrad/2*pos(1)
          extfld(5) = -bgrad/2*pos(2)
          extfld(6) = b0 + bb*(zz-zshift)*(zz-zshift) + cc*(zz-zshift)* &
                             (zz-zshift)*(zz-zshift)*(zz-zshift)
        else
          !stop
          extfld(4) = 0.0
          extfld(5) = 0.0
          extfld(6) = 0.0
        endif

        end subroutine getfld_Sol

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get the discrete Br, Bz as a function or 
        !> "r" at given "z".
        !--------------------------------------------------------------------------------------
        subroutine  getfldt_Sol(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Sol), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld
        double precision :: hz,hr,zz,efz,efr,br,bz,rr,tmpcos,tmpsin
        double precision :: scale,tt,z
        integer :: iz,iz1,ir,ir1

        z = pos(3)
        if((z.gt.this%Param(1)).and.(z.lt.(this%Param(1)+this%Length))) then
          scale = this%Param(2)
          zz = z-this%Param(1)
          hz = (fldata%ZmaxRft-fldata%ZminRft)/fldata%NzIntvRft
          hr = (fldata%RmaxRft-fldata%RminRft)/fldata%NrIntvRft
          if(zz.le.0) then
            iz = 1
            efz = zz/hz
          else if(zz.ge.(fldata%ZmaxRft-fldata%ZminRft)) then
            iz = fldata%NzIntvRft
            efz = (zz-(iz-1)*hz)/hz
          else
            iz = zz/hz + 1
            efz = (zz-(iz-1)*hz)/hz
          endif
          iz1 = iz+1
 
          rr = sqrt(pos(1)*pos(1)+pos(2)*pos(2))
          ir = rr/hr + 1
          if(rr.eq.0.0) then
            rr = 1.0e-12
          endif
          if(ir.gt.fldata%NrIntvRft) then
             print*,"ir: ",ir,rr,pos(1),pos(2),fldata%NzIntvRft,&
             fldata%NrIntvRft,fldata%ZmaxRft,fldata%ZminRft,fldata%RmaxRft,&
             fldata%RminRft
!             print*,"ir: ",ir,rr,pos(1),pos(2)
             stop
          endif
          efr = (rr-(ir-1)*hr)/hr
          ir1 = ir+1
 
          ! get the external field along "r" at given "z" by linear interpolation
          ! Here, we have only Br, Bz.
          br = (fldata%brdatat(ir,iz)*(1.0-efz)*(1.0-efr)+&
               fldata%brdatat(ir,iz1)*efz*(1.0-efr)+&
               fldata%brdatat(ir1,iz1)*efz*efr+fldata%brdatat(ir1,iz)*&
               (1.0-efz)*efr)
          bz = (fldata%bzdatat(ir,iz)*(1.0-efz)*(1.0-efr)+fldata%bzdatat(ir,iz1)*&
               efz*(1.0-efr)+&
               fldata%bzdatat(ir1,iz1)*efz*efr+fldata%bzdatat(ir1,iz)*&
               (1.0-efz)*efr)
 
          extfld(1) = 0.0
          extfld(2) = 0.0
          extfld(3) = 0.0
!          extfld(4) = 0.0
!          extfld(5) = 0.0
!          extfld(6) = 0.0
          extfld(4) = scale*br*pos(1)/rr
          extfld(5) = scale*br*pos(2)/rr
          extfld(6) = scale*bz
        else
          extfld = 0.0
        endif
!        write(14,100)z,1.0*ir,1.0*iz,br,bz,&
!         efz,efr,fldata%brdatat(ir,iz),fldata%brdatat(ir1,iz),&
!         fldata%bzdatat(ir,iz),fldata%bzdatat(ir1,iz),extfld(6)

100     format(12(1x,e13.6))

        end subroutine getfldt_Sol

      end module Solclass
