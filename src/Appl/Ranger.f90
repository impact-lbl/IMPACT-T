!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! Rangerclass: particle range class of APPLICATION                        
!                 layer.
! MODULE  : ... Rangerclass 
! VERSION : ... 1.0
!> @author
!> Ji Qiang
! DESCRIPTION:
!> find the global range of computation domain class implementation.
! Comments:  
!----------------------------------------------------------------
      module Rangerclass
        use Pgrid2dclass
        interface globalrange
          module procedure globalrange1,globalrange2,globalrange3,globalrange4
        end interface
      contains
      !--------------------------------------------------------------------------------------
      !> @brief
      !> find the range (minimum and maxium of a single bunch) and the centroid
      !> of a single beam bunch.
      !--------------------------------------------------------------------------------------
      subroutine singlerange(partcls,nplc,nptot,range,center)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nplc, nptot
      double precision, dimension(:,:) :: partcls
      double precision, dimension(6), intent(out) :: range,center
      double precision, dimension(3) :: localrange1,localrange2,temp1,temp2
      double precision, dimension(6) :: centerlc
      integer :: i, j, ierr

      do i = 1, 3
        localrange1(i) = 1000000.0
        localrange2(i) = -1000000.0
        centerlc(2*i-1) = 0.0
        centerlc(2*i) = 0.0
      enddo

!      print*,"nplc, partcls: ",nplc,partcls
      do i = 1, nplc
        do j = 1, 3
          if(localrange1(j) > partcls(2*j-1,i)) then
            localrange1(j) = partcls(2*j-1,i)
          endif
          if(localrange2(j) < partcls(2*j-1,i)) then
            localrange2(j) = partcls(2*j-1,i)
          endif
          centerlc(2*j-1) = centerlc(2*j-1) + partcls(2*j-1,i)
          centerlc(2*j) = centerlc(2*j) + partcls(2*j,i)
        enddo
      enddo

!      print*,"localrange1: ",localrange1
!      print*,"localrange2: ",localrange1
!      print*,"centerlc: ",centerlc
      
      call MPI_ALLREDUCE(localrange1,temp1,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,&
                         ierr);
      call MPI_ALLREDUCE(localrange2,temp2,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,&
                         ierr);
      call MPI_ALLREDUCE(centerlc,center,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,&
                         ierr);

      if(nptot.gt.0) then
        do i = 1, 3
          range(2*i-1) = temp1(i)
          range(2*i) = temp2(i)
          center(2*i-1) = center(2*i-1)/nptot
          center(2*i) = center(2*i)/nptot
        enddo
      else
        do i = 1, 3
          range(2*i-1) = temp1(i)
          range(2*i) = temp2(i)
          center(2*i-1) = center(2*i-1)
          center(2*i) = center(2*i)
        enddo
      endif

      end subroutine singlerange

      !--------------------------------------------------------------------------------------
      !> @brief
      !> find the global range of "nbunch" of beam bunch with global periodic
      !> length "perdlen".
      !--------------------------------------------------------------------------------------
      subroutine globalrange1(Brange,perdlen,grange,gammaz,zcent,nptot,nbunch)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nbunch
      double precision, dimension(:,:) :: Brange
      double precision, intent(in) :: perdlen
      double precision, dimension(6), intent(out):: grange
      double precision, intent(out) :: gammaz,zcent
      integer, dimension(:) :: nptot
      double precision :: sumz, sumpz, eps,xmin,xmax,ymin,ymax
      integer :: i, totptcl

!      eps = 0.002
      eps = 0.0
    
      totptcl = 0
      do i = 1, nbunch
        totptcl = totptcl + nptot(i)
      enddo

      grange(1) = Brange(1,1)
      grange(2) = Brange(2,1)
      grange(3) = Brange(3,1)
      grange(4) = Brange(4,1)
      do i = 1, nbunch
        if(grange(1) > Brange(1,i)) then
          grange(1) = Brange(1,i)
        endif
        if(grange(2) < Brange(2,i)) then
          grange(2) = Brange(2,i)
        endif
        if(grange(3) > Brange(3,i)) then
          grange(3) = Brange(3,i)
        endif
        if(grange(4) < Brange(4,i)) then
          grange(4) = Brange(4,i)
        endif
      enddo

      xmin = grange(1) - eps*(grange(2)-grange(1))
      xmax = grange(2) + eps*(grange(2)-grange(1))
      ymin = grange(3) - eps*(grange(4)-grange(3))
      ymax = grange(4) + eps*(grange(4)-grange(3))
      grange(1) = xmin
      grange(2) = xmax
      grange(3) = ymin
      grange(4) = ymax

      sumz = 0.0
      sumpz = 0.0
      do i = 1, nbunch
        sumz = sumz + Brange(11,i)*nptot(i)
        sumpz = sumpz + Brange(12,i)*nptot(i)
      enddo
      sumz = sumz/totptcl
      sumpz = sumpz/totptcl
      zcent = sumz
      gammaz = sqrt(1.0+sumpz*sumpz)

      grange(5) = sumz - 0.5*perdlen
      grange(6) = sumz + 0.5*perdlen

      end subroutine globalrange1

      !--------------------------------------------------------------------------------------
      !> @brief
      !> find the global range of "nbunch" of beam bunch.
      !--------------------------------------------------------------------------------------
      subroutine globalrange2(Brange,grange,gammaz,zcent,nptot,nbunch)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nbunch
      double precision, dimension(:,:) :: Brange
      double precision, dimension(6), intent(out):: grange
      double precision, intent(out) :: gammaz,zcent
      integer, dimension(:) :: nptot
      double precision :: sumz, sumpz, eps,xmin,xmax,ymin,ymax,zmin,zmax
      integer :: i, totptcl

!      eps = 0.002
      eps = 0.0
   
      totptcl = 0
      do i = 1, nbunch
        totptcl = totptcl + nptot(i)
      enddo

      grange(1) = Brange(1,1)
      grange(2) = Brange(2,1)
      grange(3) = Brange(3,1)
      grange(4) = Brange(4,1)
      grange(5) = Brange(5,1)
      grange(6) = Brange(6,1)
      do i = 1, nbunch
        if(grange(1) > Brange(1,i)) &
          grange(1) = Brange(1,i)
        if(grange(2) < Brange(2,i)) &
          grange(2) = Brange(2,i)
        if(grange(3) > Brange(3,i)) &
          grange(3) = Brange(3,i)
        if(grange(4) < Brange(4,i)) &
          grange(4) = Brange(4,i)
        if(grange(5) > Brange(5,i)) &
          grange(5) = Brange(5,i)
        if(grange(6) < Brange(6,i)) &
          grange(6) = Brange(6,i)
      enddo

      xmin = grange(1) - eps*(grange(2)-grange(1))
      xmax = grange(2) + eps*(grange(2)-grange(1))
      ymin = grange(3) - eps*(grange(4)-grange(3))
      ymax = grange(4) + eps*(grange(4)-grange(3))
      zmin = grange(5) - eps*(grange(6)-grange(5))
      zmax = grange(6) + eps*(grange(6)-grange(5))
      grange(1) = xmin
      grange(2) = xmax
      grange(3) = ymin
      grange(4) = ymax
      grange(5) = zmin
      grange(6) = zmax

      sumz = 0.0
      sumpz = 0.0
      do i = 1, nbunch
        sumz = sumz + Brange(11,i)*nptot(i)
        sumpz = sumpz + Brange(12,i)*nptot(i)
      enddo
      if(totptcl.gt.0) then
        sumz = sumz/totptcl
        sumpz = sumpz/totptcl
      else
        sumz = -1.0e6
      endif
      zcent = sumz
      gammaz = sqrt(1.0+sumpz*sumpz)

      end subroutine globalrange2

      !--------------------------------------------------------------------------------------
      !> @brief
      !> find the global range of "nbunch" of beam bunch.
      !> find the global range of "nbunch" of beam bunch with global periodic
      !> length "perdlen".
      !--------------------------------------------------------------------------------------
      subroutine globalrange3(Brange,grange,gammaz,nptot,nbunch,nx,ny,perdlen,zcent)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nbunch,nx,ny
      double precision, dimension(:,:) :: Brange
      double precision, dimension(6), intent(out):: grange
      double precision, intent(out) :: gammaz,zcent
      double precision, intent(in) :: perdlen
      integer, dimension(:) :: nptot
      double precision :: sumz, sumpz, epsx,epsy,xmin,xmax,ymin,ymax
      integer :: i, totptcl

      epsx = 1.0/(nx-3.0)
      epsy = 1.0/(ny-3.0)
     
      totptcl = 0
      do i = 1, nbunch
        totptcl = totptcl + nptot(i)
      enddo

      grange(1) = Brange(1,1)
      grange(2) = Brange(2,1)
      grange(3) = Brange(3,1)
      grange(4) = Brange(4,1)
      do i = 1, nbunch
        if(grange(1) > Brange(1,i)) &
          grange(1) = Brange(1,i)
        if(grange(2) < Brange(2,i)) &
          grange(2) = Brange(2,i)
        if(grange(3) > Brange(3,i)) &
          grange(3) = Brange(3,i)
        if(grange(4) < Brange(4,i)) &
          grange(4) = Brange(4,i)
      enddo
      xmin = grange(1) - epsx*(grange(2)-grange(1))
      xmax = grange(2) + epsx*(grange(2)-grange(1))
      ymin = grange(3) - epsy*(grange(4)-grange(3))
      ymax = grange(4) + epsy*(grange(4)-grange(3))
      grange(1) = xmin
      grange(2) = xmax
      grange(3) = ymin
      grange(4) = ymax
  
      sumz = 0.0
      sumpz = 0.0
      do i = 1, nbunch
        sumz = sumz + Brange(11,i)*nptot(i)
        sumpz = sumpz + Brange(12,i)*nptot(i)
      enddo
      sumz = sumz/totptcl
      sumpz = sumpz/totptcl
      zcent = sumz
      gammaz = sqrt(1.0+sumpz*sumpz)

      grange(5) = sumz - 0.5*perdlen
      grange(6) = sumz + 0.5*perdlen
  
      end subroutine globalrange3

      !--------------------------------------------------------------------------------------
      !> @brief
      !> find the global range of "nbunch" of beam bunch.
      !> find the global range of "nbunch" of beam bunch with global periodic
      !> length "perdlen".
      !--------------------------------------------------------------------------------------
      subroutine globalrange4(Brange,grange,gammaz,nptot,nbunch,nx,ny,nz,zcent)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nbunch,nx,ny,nz
      double precision, dimension(:,:) :: Brange
      double precision, dimension(6), intent(out):: grange
      double precision, intent(out) :: gammaz,zcent
      integer, dimension(:) :: nptot
      double precision :: sumz, sumpz, epsx,epsy,epsz,xmin,xmax,ymin,ymax,&
                          zmin,zmax
      integer :: i, totptcl

      epsx = 1.0/(nx-3.0)
      epsy = 1.0/(ny-3.0)
      epsz = 1.0/(nz-3.0)

      totptcl = 0
      do i = 1, nbunch
        totptcl = totptcl + nptot(i)
      enddo

      grange(1) = Brange(1,1)
      grange(2) = Brange(2,1)
      grange(3) = Brange(3,1)
      grange(4) = Brange(4,1)
      grange(5) = Brange(5,1)
      grange(6) = Brange(6,1)
      do i = 1, nbunch
        if(grange(1) > Brange(1,i)) &
          grange(1) = Brange(1,i)
        if(grange(2) < Brange(2,i)) &
          grange(2) = Brange(2,i)
        if(grange(3) > Brange(3,i)) &
          grange(3) = Brange(3,i)
        if(grange(4) < Brange(4,i)) &
          grange(4) = Brange(4,i)
        if(grange(5) > Brange(5,i)) &
          grange(5) = Brange(5,i)
        if(grange(6) < Brange(6,i)) &
          grange(6) = Brange(6,i)
      enddo

      xmin = grange(1) - epsx*(grange(2)-grange(1))
      xmax = grange(2) + epsx*(grange(2)-grange(1))
      ymin = grange(3) - epsy*(grange(4)-grange(3))
      ymax = grange(4) + epsy*(grange(4)-grange(3))
      zmin = grange(5) - epsz*(grange(6)-grange(5))
      zmax = grange(6) + epsz*(grange(6)-grange(5))
      grange(1) = xmin
      grange(2) = xmax
      grange(3) = ymin
      grange(4) = ymax
      grange(5) = zmin
      grange(6) = zmax

      sumz = 0.0
      sumpz = 0.0
      do i = 1, nbunch
        sumz = sumz + Brange(11,i)*nptot(i)
        sumpz = sumpz + Brange(12,i)*nptot(i)
      enddo
      sumz = sumz/totptcl
      sumpz = sumpz/totptcl
      zcent = sumz
      gammaz = sqrt(1.0+sumpz*sumpz)

      end subroutine globalrange4

      end module Rangerclass
