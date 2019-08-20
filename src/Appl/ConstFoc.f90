!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! ConstFocclass: 3D constant focusing beam line element class
!             in Lattice module of APPLICATION layer.
! MODULE  : ... ConstFocclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang
! DESCRIPTION:
!> This class defines the linear transfer map and field for the 3d constant focusing beam line elment.
!
! Comments:
!----------------------------------------------------------------
      module ConstFocclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 5
        type ConstFoc
          !Itype = 2
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : x focusing gradient: kx0^2
          !      (3) : y focusing gradient: ky0^2 
          !      (4) : z focusing gradient: kz0^2 
          !      (5) : radius 
        end type ConstFoc
        interface getparam_ConstFoc
          module procedure getparam1_ConstFoc,  &
                          getparam2_ConstFoc,   &
                          getparam3_ConstFoc
        end interface
        interface setparam_ConstFoc
          module procedure setparam1_ConstFoc,  &
                          setparam2_ConstFoc, setparam3_ConstFoc
        end interface
      contains
        subroutine construct_ConstFoc(this,numseg,nmpstp,type,blength)
        implicit none
        type (ConstFoc), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_ConstFoc
   
        subroutine setparam1_ConstFoc(this,i,value)
        implicit none
        type (ConstFoc), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_ConstFoc

        subroutine setparam2_ConstFoc(this,values)
        implicit none
        type (ConstFoc), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_ConstFoc

        subroutine setparam3_ConstFoc(this,numseg,nmpstp,type,blength)
        implicit none
        type (ConstFoc), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_ConstFoc
   
        subroutine getparam1_ConstFoc(this,i,blparam) 
        implicit none 
        type (ConstFoc), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_ConstFoc
  
        subroutine getparam2_ConstFoc(this,blparams)
        implicit none
        type (ConstFoc), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_ConstFoc

        subroutine getparam3_ConstFoc(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (ConstFoc), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_ConstFoc
       

        subroutine  getfldZ_ConstFoc(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (ConstFoc), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,zedge
        double precision :: xkperp02x,xkperp02y,xklong02,qmcc,gambetz,pz

        zedge = this%Param(1)
        xkperp02x = this%Param(2)
        xkperp02y = this%Param(3)
        xklong02 = this%Param(4)
        !needs to be set
        qmcc = 1.0
        gambetz = 1.0
        pz = gambetz 

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = xklong02*pos(4)*pz*pz*pz*Scxl/qmcc
        extfld(4) = -xkperp02x*pos(2)*pz/qmcc/Clight
        extfld(5) = xkperp02y*pos(1)*pz/qmcc/Clight
        extfld(6) = 0.0

        end subroutine getfldZ_ConstFoc
        
        subroutine  getfld_ConstFoc(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (ConstFoc), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,zedge
        double precision :: xkperp02x,xkperp02y,xklong02,qmcc,gambetz,pz

        zedge = this%Param(1)
        xkperp02x = this%Param(2)
        xkperp02y = this%Param(3)
        xklong02 = this%Param(4)
        !needs to be set
        qmcc = 1.0
        gambetz = 1.0
        pz = gambetz 

        extfld(1) = -xkperp02x*pos(1)
        extfld(2) = -xkperp02y*pos(2)
        extfld(3) = -xklong02*pos(3) 
        extfld(4) = 0.0
        extfld(5) = 0.0
        extfld(6) = 0.0

        end subroutine getfld_ConstFoc
      end module ConstFocclass
