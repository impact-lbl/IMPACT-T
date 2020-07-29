!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! Quadrupoleclass: Quadrupole beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang
! Description: This class defines the linear transfer map and field
!              for the quadrupole beam line elment.
! Comments:
!----------------------------------------------------------------
      module Quadrupoleclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 11
        type Quadrupole
          !Itype = 1
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : quad gradient
          !      (3) : file ID
          !      (4) : radius
          !      (5) : x misalignment error
          !      (6) : y misalignment error
          !      (7) : rotation error x
          !      (8) : rotation error y
          !      (9) : rotation error z
          !      (10) : frequency of rf quadrupole (Hz)
          !      (11) : phase of rf quadrupole (degree)
        end type Quadrupole
        interface getparam_Quadrupole
          module procedure getparam1_Quadrupole,  &
                          getparam2_Quadrupole,   &
                          getparam3_Quadrupole
        end interface
        interface setparam_Quadrupole
          module procedure setparam1_Quadrupole,  &
                          setparam2_Quadrupole, setparam3_Quadrupole
        end interface
      contains
        subroutine construct_Quadrupole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Quadrupole), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Quadrupole
   
        subroutine setparam1_Quadrupole(this,i,value)
        implicit none
        type (Quadrupole), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Quadrupole

        subroutine setparam2_Quadrupole(this,values)
        implicit none
        type (Quadrupole), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_Quadrupole

        subroutine setparam3_Quadrupole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Quadrupole), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Quadrupole
   
        subroutine getparam1_Quadrupole(this,i,blparam) 
        implicit none 
        type (Quadrupole), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Quadrupole
  
        subroutine getparam2_Quadrupole(this,blparams)
        implicit none
        type (Quadrupole), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Quadrupole

        subroutine getparam3_Quadrupole(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Quadrupole), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Quadrupole
       

        !get external field with displacement and rotation errors.
        subroutine  getflderr_Quadrupole(pos,extfld,this,dx,dy,anglex,&
                                         angley,anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge
        double precision, dimension(3) :: temp,tmp
        real*8 :: bgradp,bgradpp
        real*8 :: eps,twopi,phs

        zedge = this%Param(1)
        !if(this%Param(3).gt.1.0e-5) then
        if(this%Param(3).gt.0.0) then
          zz = pos(3)-zedge
          !call getfldfrg_Quadrupole(zz,this,bgrad)
          !call getfldfrgAna_Quadrupole(zz,this,bgrad)
          call getfldfrgAna2_Quadrupole(zz,this,bgrad,bgradp,bgradpp)
        else
          bgrad = this%Param(2)
        endif

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

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
!        extfld(4) = bgrad*tmp(2)
!        extfld(5) = bgrad*tmp(1)
!        extfld(6) = 0.0
        if(this%Param(3).gt.0.0) then
          extfld(4) = bgrad*tmp(2) - &
                      bgradpp*(tmp(2)**3+3*tmp(1)**2*tmp(2))/12
          extfld(5) = bgrad*tmp(1) - &
                      bgradpp*(tmp(1)**3+3*tmp(1)*tmp(2)**2)/12
          extfld(6) = bgradp*tmp(1)*tmp(2)
        else
          extfld(4) = bgrad*tmp(2)
          extfld(5) = bgrad*tmp(1)
          extfld(6) = 0.0
        endif

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        eps = 1.0d-5
        twopi = 4*asin(1.0d0)
        if(this%Param(10).gt.eps) then
          phs = twopi*this%Param(10)*pos(4)+this%Param(11)/360.0d0*twopi
          extfld(4) = extfld(4)*cos(phs)
          extfld(5) = extfld(5)*cos(phs)
          extfld(6) = extfld(6)*cos(phs)
        endif

        end subroutine getflderr_Quadrupole
        
        !get external field without displacement and rotation errors.
        !here, the skew quad can can be modeled with nonzero anglez
        subroutine getfld_Quadrupole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Quadrupole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge,bgradp,bgradpp,anglez
        real*8, dimension(2) :: tmp
        real*8, dimension(3) :: temp
        real*8 :: eps,twopi,phs

        zedge = this%Param(1)
        zz=pos(3)-zedge
        if(this%Param(3).gt.0.0) then
          !call getfldfrg_Quadrupole(zz,this,bgrad)
          !call getfldfrgAna_Quadrupole(zz,this,bgrad)
          call getfldfrgAna2_Quadrupole(zz,this,bgrad,bgradp,bgradpp)
          !bgradp = 0.0
          !bgradpp = 0.0
        else
          bgrad = this%Param(2)
        endif

        !j.q. 08/27/08 added skew capability to quad.
        anglez = this%Param(9)
        tmp(1) = pos(1)*cos(anglez) + pos(2)*sin(anglez)
        tmp(2) = -pos(1)*sin(anglez) + pos(2)*cos(anglez)

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        if(this%Param(3).gt.0.0) then
          temp(1) = bgrad*tmp(2) - &
                      bgradpp*(tmp(2)**3+3*tmp(1)**2*tmp(2))/12
          temp(2) = bgrad*tmp(1) - &
                      bgradpp*(tmp(1)**3+3*tmp(1)*tmp(2)**2)/12

          temp(3) = bgradp*tmp(1)*tmp(2)
        else
          !extfld(4) = bgrad*pos(2)
          !extfld(5) = bgrad*pos(1)
          !extfld(6) = 0.0
          temp(1) = bgrad*tmp(2)
          temp(2) = bgrad*tmp(1)
          temp(3) = 0.0
        endif

        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        eps = 1.0d-5
        twopi = 4*asin(1.0d0)
        if(this%Param(10).gt.eps) then
          phs = twopi*this%Param(10)*pos(4)+this%Param(11)/360.0d0*twopi
          extfld(4) = extfld(4)*cos(phs)
          extfld(5) = extfld(5)*cos(phs)
          extfld(6) = extfld(6)*cos(phs)
        endif

        end subroutine getfld_Quadrupole

        !interpolate the field from the SC rf cavity onto bunch location.
        subroutine getfldfrg_Quadrupole(zz,this,bgrad)
        implicit none
        include 'mpif.h'
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: zz
        double precision, intent(out) :: bgrad
        double precision:: hstep,slope
        integer :: klo,khi,k
        integer :: my_rank,ierr

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
        bgrad =edat(klo)+slope*(zz-zdat(klo))

        end subroutine getfldfrg_Quadrupole

        subroutine getfldfrgAna_Quadrupole(zz,this,bgrad)
        implicit none
        include 'mpif.h'
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: zz
        double precision, intent(out) :: bgrad
        double precision :: bb,dd,z1,z2

        bb = this%Param(2)
        z1 = 2*(this%Length - this%Param(3))/2
        z2 = this%Length - z1
        dd = 2*this%Param(4) 
        if(zz.lt.z1) then
          bgrad = bb/(1.0+exp(-0.00004+4.518219*(-(zz-z1)/dd-1.5)))
        else if(zz.gt.z2) then
          bgrad = bb/(1.0+exp(-0.00004+4.518219*((zz-z2)/dd-1.5)))
        else 
          bgrad = bb
        endif

        !write(1,*)zz,bgrad

        end subroutine getfldfrgAna_Quadrupole

        subroutine getfldfrgAna2_Quadrupole(zz,this,bgrad,bgradp,bgradpp)
        implicit none
        include 'mpif.h'
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: zz
        double precision, intent(out) :: bgrad
        double precision :: bgradp,bgradpp
        double precision :: bb,dd,z10,z20,z3,z4,tmpz
        double precision :: c1,c2,s1,s2,s3,s4
        double precision :: tmp1
 
        c1 = -0.00004d0
        c2 = 4.518219d0
        bb = this%Param(2)
        dd = 2*this%Param(4)
        !coordinate 0 points for Enge function.
        z10 = (this%Length - this%Param(3))/2 !for entrance
        z20 = this%Length - z10               !for exit
        !fringe field range inside coordinate 0 point.
        !this could be different for different Enge coefficients.
        z3 = z10 + 1.5*dd  !for entrance
        z4 = z20 - 1.5*dd  !for exit
        if(zz.lt.z3) then
          tmpz = -(zz - z10)/dd
          s1 = c1 + c2*tmpz
          bgrad = bb/(1.0+exp(s1))
          s2 = -c2/dd !first derivative
          tmp1 = exp(s1)*s2
          bgradp = -bb*tmp1/(1.0+exp(s1))**2
          s3 = 0.0 !second derivative
          s4 = exp(s1)*(s3+s2*s2)
          bgradpp = bb*(-s4/(1.0+exp(s1))**2+2*tmp1*tmp1/(1.0+exp(s1))**3)
        else if(zz.gt.z4) then
          tmpz = (zz - z20)/dd
          s1 = c1 + c2*tmpz
          bgrad = bb/(1.0+exp(s1))
          s2 = c2/dd
          tmp1 = exp(s1)*s2
          bgradp = -bb*tmp1/(1.0+exp(s1))**2
          s3 = 0.0
          s4 = exp(s1)*(s3 + s2*s2)
          bgradpp = bb*(-s4/(1.0+exp(s1))**2+2*tmp1*tmp1/(1.0+exp(s1))**3)
        else
          bgrad = bb
          bgradp = 0.0
          bgradpp = 0.0
        endif

        !write(1,*)zz,bgrad,bgradp,bgradpp
 
        end subroutine getfldfrgAna2_Quadrupole

        !get external field with displacement and rotation errors.
        subroutine  getflderrt_Quadrupole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Quadrupole), intent(in) :: this
        double precision :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge,bgradp,bgradpp
        double precision, dimension(3) :: temp,tmp
        real*8 :: eps,twopi,phs

        zedge = this%Param(1)
        if(this%Param(3).gt.0.0) then
          zz = pos(3)-zedge
          !call getfldfrg_Quadrupole(zz,this,bgrad)
          call getfldfrgAna2_Quadrupole(zz,this,bgrad,bgradp,bgradpp)
        else
          bgrad = this%Param(2)
        endif

        dx = this%Param(5)
        dy = this%Param(6)
        anglex = this%Param(7)
        angley = this%Param(8)
        anglez = this%Param(9)

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
!        extfld(4) = bgrad*tmp(2)
!        extfld(5) = bgrad*tmp(1)
!        extfld(6) = 0.0
        if(this%Param(3).gt.0.0) then
          extfld(4) = bgrad*tmp(2) - &
                      bgradpp*(tmp(2)**3+3*tmp(1)**2*tmp(2))/12
          extfld(5) = bgrad*tmp(1) - &
                      bgradpp*(tmp(1)**3+3*tmp(1)*tmp(2)**2)/12

          extfld(6) = bgradp*tmp(1)*tmp(2)
        else
          extfld(4) = bgrad*tmp(2)
          extfld(5) = bgrad*tmp(1)
          extfld(6) = 0.0
        endif

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        eps = 1.0d-5
        twopi = 4*asin(1.0d0)
        if(this%Param(10).gt.eps) then
          phs = twopi*this%Param(10)*pos(4)+this%Param(11)/360.0d0*twopi
          extfld(4) = extfld(4)*cos(phs)
          extfld(5) = extfld(5)*cos(phs)
          extfld(6) = extfld(6)*cos(phs)
        endif

        end subroutine getflderrt_Quadrupole
      end module Quadrupoleclass
