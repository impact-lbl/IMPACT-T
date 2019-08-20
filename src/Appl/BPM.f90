!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! BPMclass: Beam position monitor class in Lattice module of APPLICATION layer.
!
! MODULE  : ... BPMclass
! VERSION : ... 1.0
!
!> @author
!> Ji Qiang
!
! DESCRIPTION: 
!> This class defines the different beam diagnostics at given beam position.        
! Comments:
!----------------------------------------------------------------
      module BPMclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 9
        type BPM
          !Itype < 0
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : radius
          !      (3) : xmax
          !      (4) : pxmax
          !      (5) : ymax
          !      (6) : pymax
          !      (7) : zmax
          !      (8) : pzmax
        end type BPM
        interface getparam_BPM
          module procedure getparam1_BPM,  &
                          getparam2_BPM,   &
                          getparam3_BPM
        end interface
        interface setparam_BPM
          module procedure setparam1_BPM,  &
                           setparam2_BPM, setparam3_BPM
        end interface
      contains
        subroutine construct_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_BPM
   
        subroutine setparam1_BPM(this,i,value)
        implicit none
        type (BPM), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_BPM

        subroutine setparam2_BPM(this,values)
        implicit none
        type (BPM), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_BPM

        subroutine setparam3_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_BPM
   
        subroutine getparam1_BPM(this,i,blparam) 
        implicit none 
        type (BPM), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_BPM
  
        subroutine getparam2_BPM(this,blparams)
        implicit none
        type (BPM), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_BPM

        subroutine getparam3_BPM(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (BPM), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_BPM

        subroutine shift_BPM(Pts1,itype,innp,nptot)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: itype,innp,nptot
        double precision, pointer, dimension(:,:) :: Pts1
        double precision:: x0lc,px0lc,y0lc,py0lc
        double precision, dimension(4) :: tmplc,tmpgl
        integer :: i,j,ierr

        tmplc = 0.0
        tmpgl = 0.0
        if(itype.eq.(-2)) then
          x0lc = 0.0
          px0lc = 0.0
          y0lc = 0.0
          py0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(1,i)
            px0lc = px0lc + Pts1(2,i)
            y0lc = y0lc + Pts1(3,i)
            py0lc = py0lc + Pts1(4,i)
          enddo

          tmplc(1) = x0lc
          tmplc(2) = px0lc
          tmplc(3) = y0lc
          tmplc(4) = py0lc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(2) = tmpgl(2)/nptot
          tmpgl(3) = tmpgl(3)/nptot
          tmpgl(4) = tmpgl(4)/nptot

          do i = 1, innp
            Pts1(1,i) = Pts1(1,i) - tmpgl(1)
            Pts1(2,i) = Pts1(2,i) - tmpgl(2)
            Pts1(3,i) = Pts1(3,i) - tmpgl(3)
            Pts1(4,i) = Pts1(4,i) - tmpgl(4)
          enddo
        else if(itype.eq.(-1)) then
          x0lc = 0.0
          y0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(1,i)
            y0lc = y0lc + Pts1(3,i)
          enddo

          tmplc(1) = x0lc
          tmplc(3) = y0lc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(3) = tmpgl(3)/nptot

          do i = 1, innp
            Pts1(1,i) = Pts1(1,i) - tmpgl(1)
            Pts1(3,i) = Pts1(3,i) - tmpgl(3)
          enddo
        else
        endif

        end subroutine shift_BPM
       
      end module BPMclass
