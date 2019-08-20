!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! DriftTubeclass: Drift space beam line element class
!             in Lattice module of APPLICATION layer.
! MODULE  : ... DriftTubeclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang
! DESCRIPTION: 
!> This class defines the linear transfer map and field for the drift space beam line elment.
!    
! Comments:
!----------------------------------------------------------------
      module DriftTubeclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 2
        type DriftTube
          !Itype = 0
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : radius
        end type DriftTube
        interface getparam_DriftTube
          module procedure getparam1_DriftTube,  &
                          getparam2_DriftTube,   &
                          getparam3_DriftTube
        end interface
        interface setparam_DriftTube
          module procedure setparam1_DriftTube,  &
                          setparam2_DriftTube, setparam3_DriftTube
        end interface
      contains
        subroutine construct_DriftTube(this,numseg,nmpstp,type,blength)
        implicit none
        type (DriftTube), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_DriftTube
   
        subroutine setparam1_DriftTube(this,i,value)
        implicit none
        type (DriftTube), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_DriftTube

        subroutine setparam2_DriftTube(this,values)
        implicit none
        type (DriftTube), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_DriftTube

        subroutine setparam3_DriftTube(this,numseg,nmpstp,type,blength)
        implicit none
        type (DriftTube), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_DriftTube
   
        subroutine getparam1_DriftTube(this,i,blparam) 
        implicit none 
        type (DriftTube), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_DriftTube
  
        subroutine getparam2_DriftTube(this,blparams)
        implicit none
        type (DriftTube), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_DriftTube

        subroutine getparam3_DriftTube(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (DriftTube), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_DriftTube
       

        subroutine  getfld_DriftTube(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (DriftTube), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = 0.0
        extfld(5) = 0.0
        extfld(6) = 0.0

        end subroutine getfld_DriftTube
        
      end module DriftTubeclass
