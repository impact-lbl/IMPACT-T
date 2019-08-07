!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! FieldQuantclass: 3D field quantity class in Field module of APPLICATION
!                 layer.
! MODULE  : ... FieldQuantclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang
! DESCRIPTION:
!> This class defines a 3-D field quantity in the accelerator.
!> The field quantity can be updated at each step. 
! Comments:
!----------------------------------------------------------------
      module FieldQuantclass
        use Timerclass
        use CompDomclass
        use Pgrid2dclass
        use FFTclass
        use Transposeclass
        use PhysConstclass
        use Fldmgerclass
        type FieldQuant
!          private
          !> num of mesh points in x and y directions.
          integer :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal
          !> num field quantity array.
          double precision, pointer, dimension(:,:,:) :: FieldQ
        end type FieldQuant
      contains
        !--------------------------------------------------------------------------------------
        !> @brief
        !> Initialize field class.
        !--------------------------------------------------------------------------------------
        subroutine construct_FieldQuant(this,innx,inny,innz,geom,grid) 
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz 
        integer :: myid, myidx, myidy, nptot,nproccol,nprocrow
        integer, allocatable, dimension(:,:,:) :: LocalTable

        call getsize_Pgrid2d(grid,nptot,nproccol,nprocrow) 
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        allocate(LocalTable(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(nproccol.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprocrow.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
  
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal))
        this%FieldQ = 0.0d0

        deallocate(LocalTable)

        end subroutine construct_FieldQuant
   
        !--------------------------------------------------------------------------------------
        !> @brief
        !> set field quantity.
        !--------------------------------------------------------------------------------------
        subroutine set_FieldQuant(this,innx,inny,innz,geom,grid,&
                                  nprx,npry) 
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz, nprx, npry 
        integer :: myid, myidx, myidy
        integer, dimension(2,0:nprx-1,0:npry-1)::LocalTable

        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(npry.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprx.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
        this%Nxlocal = innx
  
        deallocate(this%FieldQ) 
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal)) 
        this%FieldQ = 0.0d0

        end subroutine set_FieldQuant

        !--------------------------------------------------------------------------------------
        !> @brief
        !> find the E and B fields in the lab frame
        !> from the potential on the grid of beam frame.
        !--------------------------------------------------------------------------------------
        subroutine gradEB_FieldQuant(innx,inny,innz,temppotent,ptsgeom,&
              grid,Flagbc,gammaz,FlagImage,egxout,egyout,egzout,bgxout,&
              bgyout,bgzout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz, Flagbc,FlagImage
        type (CompDom), intent(in) :: ptsgeom
        double precision,dimension(innx,inny,innz),intent(inout) :: temppotent
        double precision,dimension(innx,inny,innz),intent(inout) :: egxout,&
               egyout,egzout
        double precision,dimension(innx,inny,innz),intent(inout) :: bgxout,&
               bgyout,bgzout
        type (Pgrid2d), intent(in) :: grid
        double precision, intent(in) :: gammaz
        double precision, dimension(3) :: msize
        double precision :: hxi, hyi, hzi
        double precision :: betC
        double precision,dimension(innx,inny,innz) :: egx,&
               egy,egz
        integer :: totnp,nproccol,nprocrow,myid,myidx,myidy
        integer :: i, j, k, yadd, zadd, innp
!        integer :: comm2d,commcol,commrow

        call getsize_Pgrid2d(grid,totnp,nproccol,nprocrow)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        if(nproccol.gt.1) then
          yadd = 1
        else
          yadd = 0
        endif
        if(nprocrow.gt.1) then
          zadd = 1
        else
          zadd = 0
        endif

        call getmsize_CompDom(ptsgeom,msize)
        !get the real length for calculation Ex,Ey,Ez
        hxi = 1.0d0/(msize(1)*Scxlt)
        hyi = 1.0d0/(msize(2)*Scxlt)
        !gammaz is due to relativistic factor
        hzi = 1.0d0/(msize(3)*Scxlt*gammaz)

        !call MPI_BARRIER(comm2d,ierr)
        !if(myid.eq.16) then
        !  print*,"before guardexch:"
        !endif

        ! Exchange potential information on guard grid cells.
        if((Flagbc.eq.1).or.(Flagbc.eq.5)) then ! 3D open 
          if(totnp.ne.1) then
            call guardexch1_Fldmger(temppotent,innx,inny,innz,grid)
          endif
        else if((Flagbc.eq.2).or.(Flagbc.eq.6)) then ! 2D open, 1D periodic
        else if(Flagbc.eq.3) then
        else if(Flagbc.eq.4) then
        else
          print*,"no such boundary condition!!!!"
          stop
        endif

! cache optimization here.
        !print*,"yadd: ",yadd,zadd,innx,inny,innz,1.0d0/hxi,1.0d0/hyi,1.0d0/hzi
        !Ex
        !egx = 0.0
        do k = zadd+1, innz-zadd
          do j = yadd+1, inny-yadd
            !egx(1,j,k) = hxi*(temppotent(1,j,k)-temppotent(2,j,k))
            egx(1,j,k) = hxi*(1.5d0*temppotent(1,j,k)-2.0d0*temppotent(2,j,k)+&
                              0.5d0*temppotent(3,j,k))
            do i = 2, innx-1
              egx(i,j,k) = 0.5d0*hxi*(temppotent(i-1,j,k)- &
                           temppotent(i+1,j,k))
            enddo
            !egx(innx,j,k) = hxi*(temppotent(innx-1,j,k)- &
            !                     temppotent(innx,j,k))
            egx(innx,j,k) = hxi*(-0.5d0*temppotent(innx-2,j,k)+&
                                  2.0d0*temppotent(innx-1,j,k)- &
                                  1.5d0*temppotent(innx,j,k))
          enddo
        enddo

        !Ey
        !egy = 0.0
        if((Flagbc.eq.3).or.(Flagbc.eq.4)) then  ! periodic in theta

        if(nproccol.gt.1) then
          do k = 1, innz
            do j = 2, inny-1
              do i = 2, innx
                egy(i,j,k) = 0.5d0*hyi*(temppotent(i,j-1,k) -  &
                             temppotent(i,j+1,k))*hxi/(i-1)
              enddo
              egy(1,j,k) = 0.0
            enddo
          enddo
        else
          do k = 1, innz
            do j = 2, inny-1
              do i = 2, innx
                egy(i,j,k) = 0.5d0*hyi*(temppotent(i,j-1,k) -  &
                             temppotent(i,j+1,k))*hxi/(i-1)
              enddo
                egy(1,j,k) = 0.0
            enddo
          enddo
          do k = 1, innz
            do i = 2, innx
              egy(i,1,k) = 0.5d0*hyi*(temppotent(i,inny-1,k) -  &
                           temppotent(i,2,k))*hxi/(i-1)
              egy(i,inny,k) = egy(i,1,k)
            enddo
            egy(1,1,k) = 0.0
            egy(1,inny,k) = 0.0
          enddo
        endif

        else  ! open in Y

        if(nproccol.gt.1) then 
          if(myidy.eq.0) then
            do k = zadd+1, innz-zadd
              do i = 1, innx
                !egy(i,yadd+1,k) = hyi*(temppotent(i,yadd+1,k)- &
                !                       temppotent(i,yadd+2,k))
                egy(i,yadd+1,k) = hyi*(1.5d0*temppotent(i,yadd+1,k)- &
                                       2*temppotent(i,yadd+2,k)+ &
                                       0.5d0*temppotent(i,yadd+3,k) )
              enddo
            enddo
            do k = zadd+1, innz-zadd
              do j = yadd+2, inny-yadd
                do i = 1, innx
                  egy(i,j,k) = 0.5d0*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          else if(myidy.eq.(nproccol-1)) then
            do k = zadd+1, innz-zadd
              do i = 1, innx
                !egy(i,inny-yadd,k) = hyi*(temppotent(i,inny-yadd-1,k)- &
                !                          temppotent(i,inny-yadd,k))
                egy(i,inny-yadd,k) = hyi*(-0.5d0*temppotent(i,inny-yadd-2,k)+ &
                                          2*temppotent(i,inny-yadd-1,k)- &
                                          1.5d0*temppotent(i,inny-yadd,k))
              enddo
            enddo
            do k = zadd+1, innz-zadd
              do j = yadd+1, inny-yadd-1
                do i = 1, innx
                  egy(i,j,k) = 0.5d0*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          else
            do k = zadd+1, innz-zadd
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egy(i,j,k) = 0.5d0*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          endif
        else
          do k = zadd+1, innz-zadd
            do i = 1, innx
              !egy(i,1,k) = hyi*(temppotent(i,1,k)-temppotent(i,2,k))
              egy(i,1,k) = hyi*(1.5d0*temppotent(i,1,k)-2*temppotent(i,2,k)+&
                                0.5d0*temppotent(i,3,k) )
            enddo
            do j = 2, inny-1
              do i = 1, innx
                egy(i,j,k) = 0.5d0*hyi*(temppotent(i,j-1,k)- &
                             temppotent(i,j+1,k))
              enddo
            enddo
            do i = 1, innx
              !egy(i,inny,k) = hyi*(temppotent(i,inny-1,k)- &
              !                     temppotent(i,inny,k))
              egy(i,inny,k) = hyi*(-0.5d0*temppotent(i,inny-2,k)+&
                                    2*temppotent(i,inny-1,k)- &
                                   1.5d0*temppotent(i,inny,k))
            enddo
          enddo
        endif

        endif

        !Ez
        !egz = 0.0
        if(nprocrow.gt.1) then 
          if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then ! 3D open
            if(myidx.eq.0) then
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  !egz(i,j,zadd+1) = hzi*(temppotent(i,j,zadd+1)- &
                  !                       temppotent(i,j,zadd+2))
                  egz(i,j,zadd+1) = hzi*(1.5d0*temppotent(i,j,zadd+1)- &
                                         2*temppotent(i,j,zadd+2)+ &
                                         0.5d0*temppotent(i,j,zadd+3) )
                enddo
              enddo
              do k = zadd+2, innz-zadd
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5d0*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
            else if(myidx.eq.(nprocrow-1)) then
              do k = zadd+1, innz-zadd-1
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5d0*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  !egz(i,j,innz-zadd) = hzi*(temppotent(i,j,innz-zadd-1)- &
                  !                          temppotent(i,j,innz-zadd))
                  egz(i,j,innz-zadd) = hzi*(-0.5d0*temppotent(i,j,innz-zadd-2)+&
                                             2*temppotent(i,j,innz-zadd-1)- &
                                             1.5d0*temppotent(i,j,innz-zadd))
                enddo
              enddo
            else
              do k = zadd+1, innz-zadd
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5d0*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
            endif
          else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then 
          ! 2D open, 1D periodic
            do k = zadd+1, innz-zadd
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egz(i,j,k) = 0.5d0*hzi*(temppotent(i,j,k-1)- &
                               temppotent(i,j,k+1))
                enddo
              enddo
            enddo
          else
            print*,"no such boundary condition!!!"
            stop
          endif
        else
          if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then ! 3D open
            do j = yadd+1, inny-yadd
              do i = 1, innx
                !egz(i,j,1) = hzi*(temppotent(i,j,1)-temppotent(i,j,2))
                egz(i,j,1) = hzi*(1.5d0*temppotent(i,j,1)-2*temppotent(i,j,2)+&
                                  0.5d0*temppotent(i,j,3))
              enddo
            enddo
          else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then 
          ! 2D open, 1D periodic
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,1)=0.5d0*hzi*(temppotent(i,j,innz-1)-temppotent(i,j,2))
              enddo
            enddo
          else
          endif
          do k = 2, innz-1
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,k) = 0.5d0*hzi*(temppotent(i,j,k-1)- &
                                      temppotent(i,j,k+1))
              enddo
            enddo
          enddo
          if((Flagbc.eq.1).or.(Flagbc.eq.3).or.(Flagbc.eq.5)) then
            do j = yadd+1, inny-yadd
              do i = 1, innx
                !egz(i,j,innz) = hzi*(temppotent(i,j,innz-1)- &
                !                     temppotent(i,j,innz))
                egz(i,j,innz) = hzi*(-0.5d0*temppotent(i,j,innz-2)+&
                                      2*temppotent(i,j,innz-1)- &
                                      1.5d0*temppotent(i,j,innz))
              enddo
            enddo
          else if((Flagbc.eq.2).or.(Flagbc.eq.4).or.(Flagbc.eq.6)) then
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,innz) = 0.5d0*hzi*(temppotent(i,j,innz-1)- &
                                     temppotent(i,j,2))
              enddo
            enddo
          else
          endif
        endif

        ! find field from potential.
!        egx = 0.5d0*hxi*(cshift(temppotent,-1,1) -  &
!                       cshift(temppotent,1,1))
!        egy = 0.5d0*hyi*(cshift(temppotent,-1,2) -  &
!                       cshift(temppotent,1,2))
!        egz = 0.5d0*hzi*(cshift(temppotent,-1,3) -  &
!                      cshift(temppotent,1,3))

!        if(totnp.eq.1) then
!          egz(:,:,2) = 0.5d0*hzi*(temppotent(:,:,innz-1)-&
!                                temppotent(:,:,3))
!          egz(:,:,innz-1) = 0.5d0*hzi*(temppotent(:,:,innz-2)-&
!                                temppotent(:,:,2))
!          egy(:,2,:) = 0.5d0*hyi*(temppotent(:,inny-1,:)-&
!                                temppotent(:,3,:))
!          egy(:,inny-1,:) = 0.5d0*hyi*(temppotent(:,inny-2,:)-&
!                                temppotent(:,2,:))
!        endif

        !Send the E field to the neibhoring guard grid to do the CIC
        !interpolation.
        if(totnp.ne.1) then
          call boundint4_Fldmger(egx,egy,egz,innx,inny,innz,grid)
        endif

        if(FlagImage.ne.1) then
          betC = sqrt(gammaz**2-1.0)/gammaz/Clight
        else
          betC = -sqrt(gammaz**2-1.0)/gammaz/Clight
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              egxout(i,j,k) = egxout(i,j,k) + gammaz*egx(i,j,k)
              bgyout(i,j,k) = bgyout(i,j,k) + gammaz*egx(i,j,k)*betC
            enddo
          enddo
        enddo
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              egyout(i,j,k) = egyout(i,j,k) + gammaz*egy(i,j,k)
              bgxout(i,j,k) = bgxout(i,j,k) - gammaz*egy(i,j,k)*betC
            enddo
          enddo
        enddo
        egzout = egzout + egz
        bgzout = 0.0

        !print*,"grad: ",sum(egx),sum(egy),sum(egz),sum(temppotent),&
        !                hxi,hyi,hzi,innx,inny,innz,myidx,myidy

        end subroutine gradEB_FieldQuant

        !--------------------------------------------------------------------------------------
        !> @brief
        !> update potential (solving Possion's equation) with 3D isolated 
        !> boundary conditions. Here, image charge potential along z also calculated
        !--------------------------------------------------------------------------------------
        subroutine update3Otnew_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr,gammaz,tmppot,flagImage,&
          zshift)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr,flagImage
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        double precision, dimension(nxlc,nylc,nzlc), intent(out):: tmppot
        type (Pgrid2d), intent(in) :: grid
        type (FieldQuant), intent(inout) :: this
        double precision, intent(in) :: gammaz,zshift
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho,tmprho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc2,nypzlc2
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              rho(i,j,k) = source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo
        
        ! +1 is from the real to complex fft.
        nsxy1 = (inxglb+1)/nproccol
        nsxy2 = (inxglb+1) - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo

        nsyz1 = 2*inyglb/nprocrow
        nsyz2 = 2*inyglb - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo

        nxpylc2 = xpystable(myidy)
        nypzlc2 = ypzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)*Scxlt
        hy = msize(2)*Scxlt
        !gammaz is due to the relativistic effect
        hz = msize(3)*gammaz*Scxlt

        ! Open boundary conditions!
        call openBC3Dnew(innx,inny,innz,rho,hx,hy,hz, &
        nxpylc2,nypzlc2,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,&
        inxglb,inyglb,inzglb,tmprho,flagImage,zshift)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*hx*hy*hz
            enddo
          enddo
        enddo

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              tmppot(i,j+jadd,k+kadd) = tmprho(i,j,k)*hx*hy*hz
            enddo
          enddo
        enddo

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update3Otnew_FieldQuant
 
        !--------------------------------------------------------------------------------------
        !> @brief
        !> Solving Poisson's equation with open BCs.
        !--------------------------------------------------------------------------------------
        subroutine openBC3Dnew(innx,inny,innz,rho,hx,hy,hz,&
           nxpylc2,nypzlc2,myidz,myidy,npz,npy,commrow,commcol,comm2d, &
           pztable,pytable,ypzstable,xpystable,inxglb,&
           inyglb,inzglb,rhoImg,flagImage,zshift)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc2,nypzlc2
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d,flagImage
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
        double precision, dimension(innx,inny,innz), intent(inout):: rho
        double precision, dimension(innx,inny,innz), intent(inout):: rhoImg
        double precision, intent(in) :: hx, hy, hz, zshift
        double precision :: scalex,scaley,scalez
        integer :: i,j,k,n1,n2,n3,nylc22,nzlc22
        double complex, allocatable, dimension(:,:,:) :: rho2out,rho3out
        double complex, allocatable, dimension(:,:,:) :: grn
        integer :: ginny,ginnz,ierr

        n1 = 2*inxglb
        n2 = 2*inyglb
        n3 = 2*inzglb

        nylc22 = nxpylc2
        nzlc22 = nypzlc2

        scalex = 1.0d0
        scaley = 1.0d0
        scalez = 1.0d0

        allocate(rho2out(n3,nylc22,nzlc22))

        call fft3d1_FFT(n1,n2,n3,innz,inny,nylc22,nzlc22,&
                1,scalex,scaley,scalez,rho,ypzstable,pztable,&
        xpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy, &
        rho2out)

        !c compute FFT of the Green function on the grid:
        ! here the +1 is from the unsymmetry of green function
        ! on double-sized grid.
        if(myidz.eq.(npz-1)) then
           ginnz = innz + 1
        else
           ginnz = innz
        endif
        if(myidy.eq.(npy-1)) then
           ginny = inny + 1
        else
           ginny = inny
        endif
        allocate(grn(n3,nylc22,nzlc22))
        !call greenf1t(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
        !       hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
        !          ypzstable,pztable,xpystable,pytable,grn)
        call greenf1tIntnew2(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn)

        ! multiply transformed charge density and transformed Green 
        ! function:
        allocate(rho3out(n3,nylc22,nzlc22))
        do k = 1, nzlc22
          do j = 1, nylc22
            do i = 1, n3
              rho3out(i,j,k) = rho2out(i,j,k)*grn(i,j,k)
            enddo
          enddo
        enddo

        deallocate(grn)

        ! inverse FFT:
        scalex = 1.0d0/dble(n1)
        scaley = 1.0d0/dble(n2)
        scalez = 1.0d0/dble(n3)
        call invfft3d1_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho3out,pztable,ypzstable,&
        pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

        deallocate(rho3out)

        if(flagImage.eq.1) then
          !c compute FFT of the shifted (along z) Green function on the grid:
          ! here the +1 is from the unsymmetry of green function
          ! on double-sized grid.
          ginnz = innz
          if(myidy.eq.(npy-1)) then
             ginny = inny + 1
          else
             ginny = inny
          endif
          allocate(grn(n3,nylc22,nzlc22))
          call greenf1tIntshift2(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn,zshift)
          !call greenf1tshift(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
          !     hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
          !        ypzstable,pztable,xpystable,pytable,grn,zshift)

          ! multiply transformed charge density and transformed Green
          ! function: - sign is due to image charge
          do k = 1, nzlc22
            do j = 1, nylc22
              do i = 1, n3
                rho2out(i,j,k) = -rho2out(i,j,k)*grn(i,j,k)
              enddo
            enddo
          enddo
 
          deallocate(grn)

          ! inverse FFT:
          scalex = 1.0d0/dble(n1)
          scaley = 1.0d0/dble(n2)
          scalez = 1.0d0/dble(n3)
          call invfft3d1Img_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho2out,pztable,ypzstable,&
          pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
          rhoImg)
        else
          rhoImg = 0.0
        endif

        deallocate(rho2out)

        end subroutine openBC3Dnew
 
        !--------------------------------------------------------------------------------------
        !> @brief
        !> green function for extended array.
        !--------------------------------------------------------------------------------------
        subroutine greenf1tIntnew2(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision, dimension (nx+2,nsizey+1,nsizez+1) :: grntmp
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0
        double precision :: rr,aa,bb,cc,dd,ee,ff,ss
        double complex :: gg,gg2
        double complex :: ggrr
        double precision, dimension(2) :: xx,yy,zz
        double precision, dimension(3) :: vv
        integer :: n,i0,j0,k0
        double precision :: recfourpi

        recfourpi = 1.0d0/(8.0d0*asin(1.0d0))

        call starttime_Timer(t0)

!        if(myidx.eq.0 .and. myidy.eq.0) then
!          print*,"into integrated Green function....."
!        endif
        gxrtable = xrtable
        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

        do k0 = 1, nsizez+1
          do j0 = 1, nsizey+1
            do i0 = 1, nx+2
              jj = j0 + nblocky
              kk = k0 + nblockx
              kkk = kk - 1
              jjj = jj - 1
              iii = i0 - 1

              vv(1) = iii*hx-hx/2
              vv(2) = jjj*hy-hy/2
              vv(3) = kkk*hz-hz/2
            
                rr = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
                aa = -vv(3)**2*atan(vv(1)*vv(2)/(vv(3)*rr))/2
                bb = -vv(2)**2*atan(vv(1)*vv(3)/(vv(2)*rr))/2
                cc = -vv(1)**2*atan(vv(2)*vv(3)/(vv(1)*rr))/2
                dd = vv(2)*vv(3)*log(vv(1)+rr)
                ee = vv(1)*vv(3)*log(vv(2)+rr)
                ff = vv(1)*vv(2)*log(vv(3)+rr)
                !aa = vv(1)**2*vv(3)+(vv(2)**2+vv(3)**2)*vv(3) + &
                !            vv(1)*vv(3)*rr
                !bb = vv(1)**2*vv(2) + vv(1)*vv(2)*rr
                !cc = vv(1)*(vv(1)**2+vv(2)**2)+vv(3)*vv(1)*(vv(3)+rr)
                !dd = vv(3)*vv(2)*(vv(3) + rr)
                !ee = vv(1)**2*vv(2)+vv(2)*(vv(2)**2+vv(3)*(vv(3)+rr))
                !ff = vv(1)*vv(3)*(vv(3)+rr)
                !ss = 4*vv(2)*vv(3)*log(vv(1)+rr) + &
                !        4*vv(1)*vv(3)*log(vv(2)+rr) + &
                !        4*vv(1)*vv(2)*log(vv(3)+rr)

                !gg2 = cmplx(0.0d0,vv(3)**2)*log(cmplx(aa**2-bb**2,2*aa*bb)/ &
                !        (aa**2+bb**2) ) + cmplx(0.0d0,vv(1)**2)*&
                !        log(cmplx(cc**2-dd**2,2*cc*dd )/(cc**2+dd**2))+&
                !        cmplx(0.0d0,vv(2)**2)*log(cmplx(ee**2-ff**2,2*ee*ff)/ &
                !        (ee**2+ff**2) ) + ss

                gg2 = aa + bb + cc + dd + ee + ff
               
              !grntmp(i0,j0,k0) = recfourpi*real(gg2)/(4*hx*hy*hz)
              grntmp(i0,j0,k0) = recfourpi*gg2/(hx*hy*hz)

            enddo
          enddo
        enddo

        do k0 = 1, nsizez
          do j0 = 1, nsizey
            do i0 = 1, nx+1
              grn(i0,j0,k0) = grntmp(i0+1,j0+1,k0+1)-grntmp(i0,j0+1,k0+1)-grntmp(i0+1,j0,k0+1)+&
                              grntmp(i0,j0,k0+1)-grntmp(i0+1,j0+1,k0)+grntmp(i0,j0+1,k0)+&
                              grntmp(i0+1,j0,k0)-grntmp(i0,j0,k0)
            enddo
          enddo
        enddo

!        print*,"inside green: ",myidx,myidy,sum(grntmp),sum(grn),hx,hy,hz

!              xx(1) = iii*hx-hx/2
!              xx(2) = iii*hx+hx/2
!              yy(1) = jjj*hy-hy/2
!              yy(2) = jjj*hy+hy/2
!              zz(1) = kkk*hz-hz/2
!              zz(2) = kkk*hz+hz/2
!       
!              !find the integrated Green function.
!              n = 0
!              do k = 1, 2
!                do j = 1, 2
!                  do i = 1, 2
!                    n = n+1
!                    rr(n) = sqrt(xx(i)**2+yy(j)**2+zz(k)**2)
!                    aa(n) = xx(i)**2*zz(k)+(yy(j)**2+zz(k)**2)*zz(k) + &
!                            xx(i)*zz(k)*rr(n)
!                    bb(n) = xx(i)**2*yy(j) + xx(i)*yy(j)*rr(n)
!                    cc(n) = xx(i)*(xx(i)**2+yy(j)**2)+zz(k)*xx(i)*(zz(k)+rr(n))
!                    dd(n) = zz(k)*yy(j)*(zz(k) + rr(n))
!                    ee(n) = xx(i)**2*yy(j)+yy(j)*(yy(j)**2+zz(k)*(zz(k)+rr(n)))
!                    ff(n) = xx(i)*zz(k)*(zz(k)+rr(n))
!                    ss(n) = 4*yy(j)*zz(k)*log(xx(i)+rr(n)) + &
!                            4*xx(i)*zz(k)*log(yy(j)+rr(n)) + &
!                            4*xx(i)*yy(j)*log(zz(k)+rr(n))
!                    gg(n) = cmplx(0.0d0,zz(k)**2)*log(cmplx(aa(n)**2-bb(n)**2,2*aa(n)*bb(n))/ &
!                            (aa(n)**2+bb(n)**2) ) + cmplx(0.0d0,xx(i)**2)*&
!                            log(cmplx(cc(n)**2-dd(n)**2,2*cc(n)*dd(n) )/(cc(n)**2+dd(n)**2))+&
!                            cmplx(0.0d0,yy(j)**2)*log(cmplx(ee(n)**2-ff(n)**2,2*ee(n)*ff(n))/ &
!                            (ee(n)**2+ff(n)**2) )
!                    gg2(n) = ss(n) +  gg(n)
!                  enddo
!                enddo
!              enddo
!              ggrr = (-gg2(1)+gg2(2)+gg2(3)-gg2(4)+gg2(5)-gg2(6)-gg2(7)+gg2(8))/4
!              grn(i0,j0,k0) = real(ggrr)/(hx*hy*hz)
!            enddo
!          enddo
!        enddo


        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0d0
        scaley = 1.0d0
        scalez = 1.0d0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along y and z dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2+1,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2+1)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2+1
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = n3/2+2,n3
              tmp3(i,j) = x1(n3-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf1tIntnew2

        !--------------------------------------------------------------------------------------
        !> @brief
        !> green function for extended array.
        !--------------------------------------------------------------------------------------
        subroutine greenf1tIntshift2(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout,zshift)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz, zshift
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision, dimension (nx+2,nsizey+1,nsizez+1) :: grntmp
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1,x2
        double complex, allocatable, dimension(:,:,:) :: x0
        double precision :: rr,aa,bb,cc,dd,ee,ff,ss
        double complex :: gg,gg2
        double complex :: ggrr
        double precision, dimension(2) :: xx,yy,zz
        double precision, dimension(3) :: vv
        integer :: n,i0,j0,k0
        double precision :: recfourpi

        call starttime_Timer(t0)

        recfourpi = 1.0d0/(8.0d0*asin(1.0d0))
!        if(myidx.eq.0 .and. myidy.eq.0) then
!          print*,"into integrated Green function....."
!        endif
        gxrtable = xrtable
!        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

        do k0 = 1, nsizez+1
          do j0 = 1, nsizey+1
            do i0 = 1, nx+2
              jj = j0 + nblocky
              kk = k0 + nblockx
              kkk = kk - 1
              jjj = jj - 1
              iii = i0 - 1

              vv(1) = iii*hx-hx/2
              vv(2) = jjj*hy-hy/2
              vv(3) = kkk*hz-hz/2 + zshift

                rr = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
                aa = vv(1)**2*vv(3)+(vv(2)**2+vv(3)**2)*vv(3) + &
                            vv(1)*vv(3)*rr
                bb = vv(1)**2*vv(2) + vv(1)*vv(2)*rr
                cc = vv(1)*(vv(1)**2+vv(2)**2)+vv(3)*vv(1)*(vv(3)+rr)
                dd = vv(3)*vv(2)*(vv(3) + rr)
                ee = vv(1)**2*vv(2)+vv(2)*(vv(2)**2+vv(3)*(vv(3)+rr))
                ff = vv(1)*vv(3)*(vv(3)+rr)
                ss = 4*vv(2)*vv(3)*log(vv(1)+rr) + &
                        4*vv(1)*vv(3)*log(vv(2)+rr) + &
                        4*vv(1)*vv(2)*log(vv(3)+rr)

                gg2 = cmplx(0.0d0,vv(3)**2)*log(cmplx(aa**2-bb**2,2*aa*bb)/ &
                        (aa**2+bb**2) ) + cmplx(0.0d0,vv(1)**2)*&
                        log(cmplx(cc**2-dd**2,2*cc*dd )/(cc**2+dd**2))+&
                        cmplx(0.0d0,vv(2)**2)*log(cmplx(ee**2-ff**2,2*ee*ff)/ &
                        (ee**2+ff**2) ) + ss
     
              grntmp(i0,j0,k0) = recfourpi*real(gg2)/(4*hx*hy*hz)
            enddo
          enddo
        enddo

        do k0 = 1, nsizez
          do j0 = 1, nsizey
            do i0 = 1, nx+1
              grn(i0,j0,k0) = grntmp(i0+1,j0+1,k0+1)-grntmp(i0,j0+1,k0+1)-grntmp(i0+1,j0,k0+1)+&
                              grntmp(i0,j0,k0+1)-grntmp(i0+1,j0+1,k0)+grntmp(i0,j0+1,k0)+&
                              grntmp(i0+1,j0,k0)-grntmp(i0,j0,k0)
            enddo
          enddo
        enddo
!        print*,"inside image shifted Green: ",myidx,myidy,sum(grntmp),sum(grn),&
!                                               hx,hy,hz

        scalex = 1.0d0
        scaley = 1.0d0
        scalez = 1.0d0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along x and y dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along x dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2)
        deallocate(x0)

        !shifted green function in the doubled half space
        do k0 = 1, nsizez+1
          do j0 = 1, nsizey+1
            do i0 = 1, nx+2
              jj = j0 + nblocky
              kk = k0 + nblockx
              kkk = kk - 1
              jjj = jj - 1
              iii = i0 - 1
              vv(1) = iii*hx-hx/2
              vv(2) = jjj*hy-hy/2
              vv(3) = zshift - hz*(nz-kkk) - hz/2

                rr = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
                aa = vv(1)**2*vv(3)+(vv(2)**2+vv(3)**2)*vv(3) + &
                            vv(1)*vv(3)*rr
                bb = vv(1)**2*vv(2) + vv(1)*vv(2)*rr
                cc = vv(1)*(vv(1)**2+vv(2)**2)+vv(3)*vv(1)*(vv(3)+rr)
                dd = vv(3)*vv(2)*(vv(3) + rr)
                ee = vv(1)**2*vv(2)+vv(2)*(vv(2)**2+vv(3)*(vv(3)+rr))
                ff = vv(1)*vv(3)*(vv(3)+rr)
                ss = 4*vv(2)*vv(3)*log(vv(1)+rr) + &
                        4*vv(1)*vv(3)*log(vv(2)+rr) + &
                        4*vv(1)*vv(2)*log(vv(3)+rr)

                gg2 = cmplx(0.0d0,vv(3)**2)*log(cmplx(aa**2-bb**2,2*aa*bb)/ &
                        (aa**2+bb**2) ) + cmplx(0.0d0,vv(1)**2)*&
                        log(cmplx(cc**2-dd**2,2*cc*dd )/(cc**2+dd**2))+&
                        cmplx(0.0d0,vv(2)**2)*log(cmplx(ee**2-ff**2,2*ee*ff)/ &
                        (ee**2+ff**2) ) + ss

              grntmp(i0,j0,k0) = recfourpi*real(gg2)/(4*hx*hy*hz)
            enddo
          enddo
        enddo
        do k0 = 1, nsizez
          do j0 = 1, nsizey
            do i0 = 1, nx+1
              grn(i0,j0,k0) = grntmp(i0+1,j0+1,k0+1)-grntmp(i0,j0+1,k0+1)-grntmp(i0+1,j0,k0+1)+&
                              grntmp(i0,j0,k0+1)-grntmp(i0+1,j0+1,k0)+grntmp(i0,j0+1,k0)+&
                              grntmp(i0+1,j0,k0)-grntmp(i0,j0,k0)
            enddo
          enddo
        enddo

        scalex = 1.0d0
        scaley = 1.0d0
        scalez = 1.0d0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along x and y dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along x dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x2(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x2,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x2(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x2(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x2)
        allocate(x2(n3/2,nsizexy,nsizeyz))
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x2,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2)
        deallocate(x0)


        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = n3/2+1,n3
              tmp3(i,j) = x2(i-nz,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)
        deallocate(x2)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf1tIntshift2

       !--------------------------------------------------------------------------------------
       !> @brief
       !> longitudinal and transverse wakefield
       !--------------------------------------------------------------------------------------
       subroutine wakefield_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                                       hz,aa,gg,leng,flagbtw)
       implicit none
       include 'mpif.h'
       integer, intent(in) :: Nz,flagbtw
       double precision, intent(in) :: hz, aa, gg, leng
       double precision, dimension(Nz,2), intent(in) :: recvdensz
       double precision, dimension(Nz), intent(in) :: xwakez,ywakez
       double precision, dimension(Nz), intent(out) :: exwake,ezwake,eywake
       double precision, dimension(2*Nz,1) :: densz2n,densz2nout,&
                         greenwake,greenwakeout
       integer :: kz,twonz,one,ksign,kkzz,i
       double precision :: scale,zz00,pilc,Z0,alpha1,alpha,zz
       double precision :: coef1,coef2,densconst,offset,zzmax
       real*8 :: tmptmp,pbtw1,pbtw2,pbtw3,pbtw4
       real*8 :: leng1,leng2

       pilc = 2*asin(1.0d0)
  
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1) 
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       !longitudinal wakefield function
       !the following longitudinal wakefield function is from
       !"Short-range dipole wakefields in accelerating structures for the NLC",
       !K. L. F. Bane, SLAC-PUB-9663, 2003.
       pilc = 2*asin(1.0d0)
       alpha1 = 0.4648
       alpha = 1.-alpha1*sqrt(gg/leng)-(1.-2*alpha1)*(gg/leng)
       !slac wake
       !zz00 = gg*(aa/(alpha*leng))**2/8 
       !fermi wake
       zz00 = 0.41*(aa/leng)**0.8*(gg/leng)**1.6*aa
       Z0 = 120*pilc
       !greenwake(1,1) = 0.0
       zz = 0.0d0
       !parameter for Elettra BTW
       pbtw1 = 1226.0d0
       pbtw2 = 3.0d-4
       pbtw3 = 0.494
       pbtw4 = 494.0

       !here, leng1 and leng2 are the length factors.
       !Elegant uses RF module length instead of cavity length.  
       !leng1 = 1.0500439034821d0
       !leng2 = 1.33408842738323d0
       !here, the 1st factor is from Macro's correction. The
       !2nd factor is from ratio of MAD cavity length to the IMPACT cavity length
       !leng1 = 0.749566724436742*0.787077969247848
       !the following modification based on Marco's NGLS talk 1.03774/1.3188
       leng1 = 0.78688d0
       !leng2 = 1.33408842738323d0
       !This factor is from ratio of MAD cavity length to the IMPACT cavity length
       !leng2 = 0.769215375902761
       !the following modification based on Marco's NGLS talk 0.346/0.4497453
       leng2 = 0.769324d0

       if(flagbtw.eq.1) then !for Backward TWS
         !The 0.5 is from S+
         greenwake(1,1) = 0.5d12*(pbtw1*exp(-sqrt(zz/pbtw2))+&
            pbtw3*2*(sqrt(zz+0.5d0*hz*Scxlt))/(hz*Scxlt)+pbtw4*sqrt(zz))
         do kz = 2, Nz+1
           zz = (kz-1)*hz*Scxlt
           greenwake(kz,1) = 1.0d12*(pbtw1*exp(-sqrt(zz/pbtw2))+&
             pbtw3*2*(sqrt(zz+0.5d0*hz*Scxlt)-sqrt(zz-0.5d0*hz*Scxlt))/(hz*Scxlt)+&
             pbtw4*sqrt(zz))
         enddo
       else if(flagbtw.eq.2) then !for Tesla standing wave structure
         !The 0.5 is from S+
         greenwake(1,1) = leng1*0.5d12*38.1*(1.165d0*exp(-sqrt(zz/3.65d-3))-&
            0.165d0)
         do kz = 2, Nz+1
           zz = (kz-1)*hz*Scxlt
           greenwake(kz,1) = leng1*1.0d12*38.1*(1.165d0*exp(-sqrt(zz/3.65d-3))-&
             0.165d0)
         enddo
       else if(flagbtw.eq.3) then !for Tesla 3rd harm. standing wave structure
         !The 0.5 is from S+
         !greenwake(1,1) = leng2*0.5d12*130.0d0*(1.075d0*exp(-sqrt(zz/2.25d-3))-&
         !   0.075d0)
         !3/18/08 changed following Sasha's new formulae
         !greenwake(1,1) = leng2*0.5d12*130.0d0*(1.02*exp(-sqrt(zz/2.0d-3))-&
         !12/15/2011 this follows Macro's new formulae
         greenwake(1,1) = leng2*0.5d12*229.8d0*exp(-sqrt(zz/0.84d-3))
         do kz = 2, Nz+1
           zz = (kz-1)*hz*Scxlt
!           greenwake(kz,1) = leng2*1.0d12*130.0d0*(1.075d0*exp(-sqrt(zz/2.25d-3))-&
!             0.075d0)
!           greenwake(kz,1) = leng2*1.0d12*130.0d0*(1.02*exp(-sqrt(zz/2.0d-3))-&
!             0.02)
           greenwake(kz,1) = leng2*1.0d12*229.8d0*exp(-sqrt(zz/0.84d-3))
         enddo
       else !for TWS
         !The 0.5 is from S+
         greenwake(1,1) = 0.5d0*Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)
         !do kz = 1, Nz+1
         do kz = 2, Nz+1
           zz = (kz-1)*hz*Scxlt
           greenwake(kz,1) = Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)
         enddo
       endif
       do kz = Nz+2, twonz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 2, Nz
         greenwake(kz,1) = 0.0
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       !do kz = 1, twonz
       !  greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       !enddo
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

       !the "-" sign is from definition of longitudinal wake function, see Wangler's book
       do kz = 1, Nz
         ezwake(kz) = -greenwakeout(kz,1)*hz*Scxlt
       enddo
       
       !The following is for comparison purpose
       !calculate the short range longitudinal wakefield by direct summation
!       do kz = 1, Nz
!         densz2n(kz,1) = 0.5d0*Z0*Clight/(pilc*aa*aa)*recvdensz(kz,1)*hz*Scxlt
!         do kkzz = kz+1, Nz
!           zz = (kkzz-kz)*hz*Scxlt
!       !    !densz2n(kz,1) = densz2n(kz,1) + 4*Z0*Clight*zz00*(1.0-(1.+sqrt(zz/zz00))*&
!       !    !                 exp(-sqrt(zz/zz00)))/(pilc*aa*aa*aa*aa)*recvdensz(kkzz,1)*&
!       !    !                 recvdensz(kkzz,2)*hz*Scxlt
!       !    !densz2n(kz,1) = densz2n(kz,1) + coef1*(zz/zz00)*recvdensz(kkzz,1)* &
!       !    !                                       recvdensz(kkzz,2)*hz*Scxlt
!       !    !densz2n(kz,1) = densz2n(kz,1) + coef1*(1.0-(1.+sqrt(zz/zz00))*&
!       !    !                 exp(-sqrt(zz/zz00)))*recvdensz(kkzz,1)*&
!       !    !                 recvdensz(kkzz,2)*hz*Scxlt
!           densz2n(kz,1) = densz2n(kz,1) + Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)*& 
!                            recvdensz(kkzz,1)*hz*Scxlt
!       
!         enddo
!       enddo
    
       !write(17,*)erwake
       !write(17,*)densz2n(1:Nz,1)
       !call flush_(17)
!       do i = 1, Nz
!         zz = (i-1)*hz*Scxlt
!         tmptmp = Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)
!         write(19,*)zz,-densz2n(i,1),ezwake(i),recvdensz(i,1),tmptmp
!       enddo
       !exwake(:) = densz2n(1:Nz,1)

!------------------------------------------------------------------------------
       !if((flagbtw.eq.2) .or. (flagbtw.eq.3)) then !no transverse wake for standing wave cavity so far
       if(flagbtw.eq.3) then !no transverse wake for standing wave cavity so far
         exwake = 0.0d0
         eywake = 0.0d0
         goto 100
       endif

       !calculate the transverse wakefield effects
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*xwakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)
 
       !tranverse wakefield
       !the following transverse wakefield function is from
       !"Short-range dipole wakefields in accelerating structures for the NLC",
       !K. L. F. Bane, SLAC-PUB-9663, 2003. here, 1.1 is an average as suggested       !on page 11 for cell 45.
       !for LCLS slac
       !zz00 = 1.1*0.169*aa*(aa/leng)**1.17*(gg/aa)**0.38
       !for Fermi Elettra
       zz00 = 1.0d0*0.169*aa*(aa/leng)**1.17*(gg/aa)**0.38
       coef1 = 4*Z0*Clight*zz00/(pilc*aa*aa*aa*aa)
       !densconst = 0.5e-6
       !offset = 0.001
       !zzmax = 0.002
       !coef2 = coef1/zz00*densconst*offset
       !print*,"aa wake: ",aa,leng,gg,zz00,Z0,coef1,coef2,offset,densconst
       ! parameter for Elettra BTW linac  
       pbtw1 = 2.8d4
       pbtw2 = 1.2d-4
       pbtw3 = 1.2d-4
       pbtw4 = 1.4d4
       if(flagbtw.eq.1) then !for Backward TWS
         greenwake(1,1) = 0.0
         do kz = 1, Nz+1
           zz = (kz-1)*hz*Scxlt
           greenwake(kz,1) = 1.0d12*(pbtw1*(1.0-(1.0+sqrt(zz/pbtw2))*exp(-sqrt(zz/pbtw3))) + &
                             pbtw4*sqrt(zz))
         enddo
       else if(flagbtw.eq.2) then !for Tesla standing wave structure
         greenwake(1,1) = 0.0
         pbtw1 = 121.d0
         pbtw2 = 0.92d-3
         pbtw3 = 0.92d-3
         pbtw4 = 0.0d0
         do kz = 1, Nz+1
           zz = (kz-1)*hz
           greenwake(kz,1) = leng1*1.0d12*(pbtw1*(1.0-(1.0+sqrt(zz/pbtw2))*exp(-sqrt(zz/pbtw3))) + &
                             pbtw4*sqrt(zz))
         enddo

       else ! for TWS
         greenwake(1,1) = 0.0
         do kz = 1, Nz+1
           zz = (kz-1)*hz*Scxlt
           greenwake(kz,1) = coef1*(1.0-(1.+sqrt(zz/zz00))*exp(-sqrt(zz/zz00)))
         enddo
       endif

       do kz = Nz+2, 2*Nz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 1, Nz
         greenwake(kz,1) = 0.0
       enddo

       ksign = 1
       !scale = 1./(twonz)
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
                greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   densz2nout)

       do kz = 1, Nz
         exwake(kz) = densz2nout(kz,1)*hz*Scxlt
       enddo

       !The following is for comparison purpose
       !calculate the short range transverse wakefield by direct summation
       !do kz = 1, Nz
       !  densz2n(kz,1) = 0.0
       !  do kkzz = kz+1, Nz
       !    zz = (kkzz-kz)*hz*Scxlt
       !    !densz2n(kz,1) = densz2n(kz,1) + 4*Z0*Clight*zz00*(1.0-(1.+sqrt(zz/zz00))*&
       !    !                 exp(-sqrt(zz/zz00)))/(pilc*aa*aa*aa*aa)*recvdensz(kkzz,1)*&
       !    !                 recvdensz(kkzz,2)*hz*Scxlt
       !    !densz2n(kz,1) = densz2n(kz,1) + coef1*(zz/zz00)*recvdensz(kkzz,1)* &
       !    !                                       recvdensz(kkzz,2)*hz*Scxlt
       !    densz2n(kz,1) = densz2n(kz,1) + coef1*(1.0-(1.+sqrt(zz/zz00))*&
       !                     exp(-sqrt(zz/zz00)))*recvdensz(kkzz,1)*&
       !                     recvdensz(kkzz,2)*hz*Scxlt
       !
       !  enddo
       !enddo
     
       !write(17,*)erwake
       !write(17,*)densz2n(1:Nz,1)
       !call flush_(17)
       !do i = 1, Nz
       !  zz = (i-1)*hz*Scxlt
       !  write(19,*)zz,densz2n(i,1),coef2*(0.5d0*zz**2-zzmax*zz+0.5d0*zzmax**2),exwake(i)
       !enddo
       !exwake(:) = densz2n(1:Nz,1)

       !vertical "Y" wakefields
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*ywakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)
 
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   greenwakeout)

       do kz = 1, Nz
         eywake(kz) = greenwakeout(kz,1)*hz*Scxlt
       enddo

100    continue

       end subroutine wakefield_FieldQuant

       !--------------------------------------------------------------------------------------
       !> @brief
       !> longitudinal and transverse wakefield on a beam using the
       !> readin longitudinal and transverse wake functions
       !--------------------------------------------------------------------------------------
       subroutine wakereadin_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                          hz,leng,ndatawk,wklong,wktran)
       implicit none
       include 'mpif.h'
       integer, intent(in) :: Nz,ndatawk
       double precision, intent(in) :: hz, leng
       double precision, dimension(Nz,2), intent(in) :: recvdensz
       double precision, dimension(Nz), intent(in) :: xwakez,ywakez
       double precision, dimension(ndatawk), intent(in) :: wklong,wktran
       double precision, dimension(Nz), intent(out) :: exwake,ezwake,eywake
       double precision, dimension(2*Nz,1) :: densz2n,densz2nout,&
                         greenwake,greenwakeout
       integer :: kz,twonz,one,ksign,kkzz,i,iz,iz1
       double precision :: scale,zz,zziz
       real*8 :: hzwake

       hzwake = leng/(ndatawk-1)
  
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1) 
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       if(Nz*hz*Scxlt.gt.leng) then
         print*,"warning: total bunch length is greater than the readin wake field range!!"
       endif

       !longitudinal wakefield function
       greenwake(1,1) = 0.5d0*wklong(1)
       do kz = 2, Nz+1
         zz = (kz-1)*hz*Scxlt
         iz = zz/hzwake + 1
         iz1 = iz + 1
         if(iz1.gt.ndatawk) then
           iz = ndatawk - 1 
           iz1 = ndatawk  
         endif
         zziz = (iz-1)*hzwake
         greenwake(kz,1) = wklong(iz)+(wklong(iz1)-wklong(iz))*(zz-zziz)/hzwake
       enddo

       do kz = Nz+2, twonz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 2, Nz
         greenwake(kz,1) = 0.0
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

       !the "-" sign is from definition of longitudinal wake function, see Wangler's book
       do kz = 1, Nz
         ezwake(kz) = -greenwakeout(kz,1)*hz*Scxlt
       enddo
       
!------------------------------------------------------------------------------
       !calculate the transverse wakefield effects
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*xwakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)
 
       !tranverse wakefield
       do kz = 1, Nz+1
         zz = (kz-1)*hz*Scxlt
         iz = zz/hzwake + 1
         iz1 = iz + 1
         if(iz1.gt.ndatawk) then
           iz = ndatawk - 1
           iz1 = ndatawk
         endif
         zziz = (iz-1)*hzwake
         greenwake(kz,1) = wktran(iz)+(wktran(iz1)-wktran(iz))*(zz-zziz)/hzwake
       enddo
       do kz = Nz+2, 2*Nz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 1, Nz
         greenwake(kz,1) = 0.0
       enddo

       ksign = 1
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
                greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   densz2nout)
       do kz = 1, Nz
         exwake(kz) = densz2nout(kz,1)*hz*Scxlt
       enddo

       !vertical "Y" wakefields
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*ywakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)
 
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   greenwakeout)

       do kz = 1, Nz
         eywake(kz) = greenwakeout(kz,1)*hz*Scxlt
       enddo

       end subroutine wakereadin_FieldQuant

        !--------------------------------------------------------------------------------------
        !> @author J. Q. 11/3/08
        !> @brief
        !> new version of the csrwake calculation. 
        !> this subroutine calculate the 1d csr wakefield including
        !> entrance, stead-state, and transitions effects
        !> here, hx, rho,...are in real unit. the return ezwake is V/m
        !--------------------------------------------------------------------------------------
        subroutine csrwakeTr_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,ezwake)
        implicit none
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp
        integer :: Ni,il,i,Nmax,j,Nsup,islp,Nisup,jj0,jj,jj1,isup
        real*8 :: tcoef,hxsup,csrss,xxp0,xxsup,rhoxxp,xblg,deltasmax2

        epstol = 1.0d-10 !tolerance for root finding
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817e-12)
        xk = -2.0d0/(3*r0**2)**(1.0d0/3.0d0)*xconst

        Nsup = 10

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24
           Ni = i
           ezwake(i) = 0.0
           xblg = (i-1)*hx

           !print*,"xxl: ",i,xx,blength,xxl,xblg
           if(xx.le.0) then
           else if(xx.le.blength) then
             if(xxl.ge.xblg) then !S-S
               do j = 1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo

               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
 
               Nisup = isup*Nsup
 
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
 
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0 
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo
             else if(xxl.lt.hx) then !Transient
               if(4*xxl.ge.(xx-ptmin)) then
               else 
                 il = (xx-4*xxl-ptmin)/hx + 1  
                 dx = xx-4*xxl-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)   
               endif

               !contribution from case B
               il = Ni-1
               dx = xx-xxl-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) + tmprho/xxl**(1.0d0/3.0d0)

               hxsup = xxl/Nsup
               jj0 = Ni-1
               xxp0 = xx - xxl
               Nisup = Nsup
               jj = jj0
               jj1 = jj0+1
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif

                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif

                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo

             else
               islp = (xx-xxl-ptmin)/hx + 1

               !transient contribution from 0 to xbl-xxl 
               if(4*xxl.ge.(xx-ptmin)) then
               else
                 il = (xx-4*xxl-ptmin)/hx + 1                 
                 dx = xx-4*xxl-(ptmin+(il-1)*hx) 
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/xxl**(1.0d0/3.0d0)
               endif 

!               print*,"islp: ",islp,xxl,xx,hx

               !case B

               !from islp to islp + 1 
               isup = xxl/hx
               hxsup = (xxl-isup*hx)/Nsup
               jj0 = islp
               xxp0 = xx - xxl
               Nisup = Nsup
               dx = xxp0 - (ptmin+(jj0-1)*hx)
               tmprho = rhonew(jj0) + dx/hx*(rhonew(jj0+1)-rhonew(jj0))
               ezwake(i) = ezwake(i) + tmprho/xxl**(1.0d0/3.0d0)
               jj = jj0
               jj1 = jj0+1
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
 
                 xx2 = xxp0 + (j-1)*hxsup
                 csrss = (xx-xx2)**(-1.0d0/3.0d0)
 
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo

               !ezwake(i) = ezwake(i) + rhonew(islp)/xxl**(1.0d0/3.0d0)
               !do j = islp, Ni-1
               do j = islp+1, Ni-1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni-1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 !csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 ezwake(i) = ezwake(i)+tcoef*&
                            (xx-xx2)**(-1.0d0/3.0d0)*rhonewp(j)*hx
               enddo
               !subgrid integration from Ni-1 to Ni
               hxsup = hx/Nsup
               isup = 1
               jj0 = Ni-isup
               xxp0 = ptmin+(jj0-1)*hx
               Nisup = isup*Nsup
               do j = 1, Nisup+1
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Nisup+1) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = xxp0 + (j-1)*hxsup
                 if(j.ne.Nisup+1) then
                   csrss = (xx-xx2)**(-1.0d0/3.0d0)
                 else
                   csrss = 2*(hxsup)**(-1.0d0/3.0d0)-&
                           (2*hxsup)**(-1.0d0/3.0d0)
                 endif
                 jj = jj0
                 jj1 = jj0+1
                 xxsup = xx2 - (ptmin+(jj-1)*hx)
                 rhoxxp = rhonewp(jj)+(rhonewp(jj1)-rhonewp(jj))*xxsup/hx
                 ezwake(i) = ezwake(i)+tcoef*csrss*rhoxxp*hxsup
               enddo
             endif
             ezwake(i) = ezwake(i)*xk
           else
             xxbar = (xx - blength)/r0
             deltasmax = (r0*phim**3/24.0d0)*(phim+4*xxbar)/(phim+xxbar)
             !print*,"xxbar: ",i,xx,xxbar,deltasmax,xblg
             if(deltasmax.ge.xblg) then !case D
               !neglect the contributions from s to s
               !do j = 1, Ni-1
               do j = 1, Ni
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i) + tcoef*rhonewp(j)/(psi+2*xxbar)*hx
               enddo
             else if(deltasmax.lt.hx) then !case C
               !print*,"deltasmaxC: ",deltasmax,hx
               deltasmax2 = (r0*phim**2/6.0d0)*(phim+3*xxbar)
               if(deltasmax2.ge.xx-ptmin) then
               else
                 il = (xx-deltasmax2-ptmin)/hx + 1
                 dx = xx-deltasmax2-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
               endif

               deltasmax2 = deltasmax
               il = (xx-deltasmax2-ptmin)/hx + 1
               dx = xx-deltasmax2-(ptmin+(il-1)*hx)
               tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
               ezwake(i) = ezwake(i) + tmprho/(phim+2*xxbar)

               !neglect the integral from xx-deltasmax to xx within hx
             else
               islp = (xx-deltasmax-ptmin)/hx + 1

               !case C from 0 to xblg-deltasmax
               deltasmax2 = (r0*phim**2/6.0d0)*(phim+3*xxbar)
               if(deltasmax2.ge.xx-ptmin) then
               else
                 il = (xx-deltasmax2-ptmin)/hx + 1
                 dx = xx-deltasmax2-(ptmin+(il-1)*hx)
                 tmprho = rhonew(il) + dx/hx*(rhonew(il+1)-rhonew(il))
                 ezwake(i) = ezwake(i) - tmprho/(phim+2*xxbar)
               endif 

               !case D
               ezwake(i) = ezwake(i) + rhonew(islp)/(phim+2*xxbar)
               !do j = islp, Ni-1
               do j = islp, Ni
                 if(j==1) then
                   tcoef = 0.5d0
                 else if (j==Ni) then
                   tcoef = 0.5d0
                 else
                   tcoef = 1.0d0
                 endif
                 xx2 = ptmin + (j-1)*hx
                 deltas = xx - xx2
                 call psiroot(r0,xxbar,deltas,psi,epstol,Nmax)
                 ezwake(i) = ezwake(i) + tcoef*rhonewp(j)/(psi+2*xxbar)*hx
               enddo
             endif
             ezwake(i) = -ezwake(i)*4/r0*xconst
           endif
        enddo

        end subroutine csrwakeTr_FieldQuant

        subroutine setval_FieldQuant(this,i,j,k,value)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        integer, intent(in) :: i, j, k
        double precision, intent(in) :: value

        this%FieldQ(i,j,k) = value

        end subroutine setval_FieldQuant

        double precision function get_FieldQuant(this,i,j,k)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(in) :: i, j, k

        get_FieldQuant = this%FieldQ(i,j,k)

        end function get_FieldQuant

        subroutine getglb_FieldQuant(this,temp)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        type (FieldQuant), intent(out) :: temp
        integer :: i, j, k, lcnz,lcny,lcnx
        double precision :: value

        lcnx = this%Nxlocal
        lcny = this%Nylocal
        lcnz = this%Nzlocal
    
        do k = 1, lcnz
          do j = 1, lcny
            do i = 1, lcnx
              value = get_FieldQuant(this,i,j,k)
              call setval_FieldQuant(temp,i,j,k,value)
            enddo
          enddo 
        enddo

        end subroutine getglb_FieldQuant

        subroutine getlcgrid_FieldQuant(this,nxlc,nylc,nzlc)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(out) :: nxlc,nylc,nzlc

        nxlc = this%Nxlocal
        nylc = this%Nylocal
        nzlc = this%Nzlocal

        end subroutine getlcgrid_FieldQuant

      !--------------------------------------------------------------------------------------
      !> @brief
      !> find the psi in equation 12 of Stupakov and Emma's paper
      !--------------------------------------------------------------------------------------
      subroutine psiroot(r0,xx,deltas,psi,eps,Nmax)
      implicit none
      integer :: Nmax
      real*8 :: r0, xx, deltas,eps,psi
      integer :: i
      real*8 :: ps0,ps1,fps0,dfps0


      ps0 = (24*deltas/r0)**(1.0d0/3.0d0)

      do i = 1, Nmax
        fps0 = r0*ps0**4+4*xx*r0*ps0**3-24*deltas*ps0-24*deltas*xx
        !print*,"ps0: ",i,ps0,fps0
        if(abs(fps0).le.eps) then
          psi = ps0
          return
        else
          dfps0 = 4*r0*ps0**3+12*xx*r0*ps0**2-24*deltas
          ps1 = ps0 - fps0/dfps0
          ps0 = ps1 
        endif
      enddo

      print*,"out of maximum interation: ",i

      end subroutine psiroot

        subroutine destruct_FieldQuant(this)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this

        deallocate(this%FieldQ) 

        end subroutine destruct_FieldQuant

        !--------------------------------------------------------------------------------------
        !> @brief
        !> This subroutine calculates the 1d csr wakefield including
        !> entrance, stead-state, and transitions effects.
        !> Here, hx, rho,...are in real units. The return ezwake is in V/m.
        !> The current version uses IGF corresponding to the four cases of Saldin et al.
        !--------------------------------------------------------------------------------------
        subroutine csrwakeTrIGF_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,gam,ezwake)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength,gam
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp,xslpN
        integer :: Ni,il,i,Nmax,j,islpN,islpN1
        integer :: myidlc,ierr,islp,islp1,islp0,jstart
        real*8 :: aa,uuh,uuh24,xslp,csrss,xxp,ssh,xk2,tcoef
        real*8 :: bb,cc,yh,phih,phpy2,csrtr,xblg,xxbarh,phimh,psimax
        real*8 :: psip2x2,psipx2,csrdrm1,csrdr1,phpxy2,phpy,csrss1,csrss2,&
                  csrtr1,csrtr2,csrdrm2,csrdr2
!        real*8 :: IcsrCaseA,IcsrCaseB,IcsrCaseC,IcsrCaseD

!        call MPI_COMM_RANK(MPI_COMM_WORLD, myidlc, ierr)

!        epstol = 1.0d-10 !tolerance for root finding
        epstol = 2.0d-9
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817d-12)
        xk2 = gam/r0

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24

           Ni = i
           ezwake(i) = 0.0

           xblg = (i-1)*hx

           !if(myidlc.eq.0) print*,"ii:",i,xxl,xx

           if(xx.le.0) then
           else if(xx.le.blength) then

             xslp = (xx/r0)*r0/2/gam/gam + xxl
             islp0 = (xx-xslp-ptmin)/hx
             islp1 = islp0 + 1

             !Case A 
             phih = xx/r0*gam
             ! IGF integral over the interior sample points for Case A:
             do j = 2, islp1-1
                 !write(*,*) 'Inside Case A!'
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrtr1-csrtr2)
             enddo
             if(islp1.ge.2) then
             !Add the upper end IGF integral for case A
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrtr1-csrtr2)
             !Add the lower end IGF integral for case A
                 xxp = ptmin + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(1)*(csrtr1-csrtr2)
             ! Special case
             else if(islp1.eq.1) then
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrtr1-csrtr2)
             endif

             ! Case B (steady-state regime)
             jstart = max(islp1,1)
             ! IGF integral over the interior sample points for Case B:
             do j = jstart+2,Ni-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
             enddo
             !add end integrals
             if(jstart.le.Ni-2) then
             !Add the upper end IGF integral for case B
                 xxp = ptmin + (Ni-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = ptmin + (Ni-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(Ni)*(csrss1-csrss2)
             !Add the lower end IGF integral for case B
                 if(islp1.lt.1) then
                   j = 1
                   xxp = ptmin + (j-1)*hx + hx/2
                   ssh = (xx-xxp)*gam**3/r0
                   csrss2 = IcsrCaseB(ssh,xk2)  
                   xxp = ptmin + (j-1)*hx 
                   ssh = (xx-xxp)*gam**3/r0
                   csrss1 = IcsrCaseB(ssh,xk2)
                   ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
                 else
                   j = islp1+1
                   xxp = ptmin + (j-1)*hx + hx/2
                   ssh = (xx-xxp)*gam**3/r0
                   csrss2 = IcsrCaseB(ssh,xk2)
                   xxp = xx-xslp
                   ssh = (xx-xxp)*gam**3/r0
                   csrss1 = IcsrCaseB(ssh,xk2)
                   ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
                 endif
             else if(jstart.eq.Ni-1) then
                 j = Ni
                 xxp = ptmin + (j-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = xx-xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
             endif

           !  ezwake(i) = ezwake(i)*xconst

             !if(myidlc.eq.0) print*,"xxl: ",xxl,xx,ezwake(i),r0,ptmin,phim
             !if(myidlc.eq.0) print*,"xxl: ",xx,ezwake(i),xxl,xx-ptmin,xslp,hx

           else

             !if(myidlc.eq.0) print*,"CD:",i,islp,xslp,hx,xx,ezwake(i)

             xxbar = (xx - blength)
             xslp= (r0*phim+xxbar)/2/gam/gam + r0*phim**3/24*&
                   (r0*phim+4*xxbar)/(r0*phim+xxbar)
             islp0 = (xx-xslp-ptmin)/hx
             islp1 = islp0 + 1
             xslpN = xxbar/2.d0/gam/gam
             islpN = (xx-xslpN-ptmin)/hx
             islpN1 = islpN + 1

             ! Case C
             xxbarh = xxbar*gam/r0
             phimh = phim*gam
             ! IGF integral over the interior sample points for Case C:
             do j = 2, islp1-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrdrm1-csrdrm2)
             enddo
             if(islp1.ge.2) then
             !Add the upper end IGF integral for case C
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrdrm1-csrdrm2)
             !Add the lower end IGF integral for case C
                 xxp = ptmin + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(1)*(csrdrm1-csrdrm2)
             ! Special case
             else if(islp1.eq.1) then
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrdrm1-csrdrm2)
             endif


             !Case D 
             psimax = phim*gam
             xxbarh = xxbar*gam/r0
             jstart = max(islp1,0)
             !write(13,*) 'xslp and xslpN',xslp,xslpN
             ! IGF integral over the interior sample points for Case D:
             do j = jstart+2,islpN1-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrdr1-csrdr2)
             enddo
             if(islpN1.ge.(jstart+2)) then
             !Add the upper end IGF integral for case D
                 ssh = xslpN*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 xxp = ptmin + (islpN1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(islpN1)*(csrdr1-csrdr2)
             !Add the lower end IGF integral for case D
                 xxp = ptmin + (jstart)*hx + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ssh = xslp*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(jstart+1)*(csrdr1-csrdr2)
             else if(islpN1.eq.(jstart+1)) then
                 ssh = xslpN*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ssh = xslp*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(jstart+1)*(csrdr1-csrdr2)
             endif


           endif
           ezwake(i) = ezwake(i)*xconst
        enddo

        end subroutine csrwakeTrIGF_FieldQuant

           function IcsrCaseA(phih,ssh,xk2)
           implicit none
           double precision:: phih,ssh,xk2,IcsrCaseA
           double precision:: bb,cc,yh,phpy
           bb = 2*phih+phih**3/3-2*ssh
           cc = phih**2+phih**4/12-2*ssh*phih
           yh = (-bb+sqrt(bb*bb-4*cc))/2
           phpy = (phih+yh)
           IcsrCaseA = xk2*(-(2*phpy+phih**3)/&
                    (phpy**2+phih**4/4)+1.0d0/ssh)
           end function IcsrCaseA


           function IcsrCaseB(ssh,xk2)
           implicit none
           double precision:: ssh,xk2,IcsrCaseB
           double precision:: aa,uuh
           aa = sqrt(64.0+144.0d0*ssh**2)
           uuh = (aa+12*ssh)**(1.0d0/3.0d0)-(aa-12*ssh)**(1.0d0/3.0d0)
           IcsrCaseB = xk2*(-4*uuh*(uuh**2+8)/ &
                    ((uuh**2+4.0d0)*(uuh**2+12.0d0)))
           end function IcsrCaseB


           function IcsrCaseC(phimh,xxbarh,ssh,xk2)
           implicit none
           double precision:: phimh,xxbarh,ssh,xk2,IcsrCaseC
           double precision:: bb,cc,yh,phpxya,phpxyb
           bb = 2*(phimh+xxbarh)-2*ssh+phimh**3/3+phimh**2*xxbarh
           cc = (phimh+xxbarh)**2+phimh**2*(phimh**2+4*phimh*xxbarh)/12-&
                      2*ssh*(phimh+xxbarh)
           yh = (-bb+sqrt(bb*bb-4*cc))/2
           phpxya = (phimh+xxbarh+yh)
           phpxyb = phimh*xxbarh+phimh*phimh/2.d0
           IcsrCaseC = xk2*(-2.d0*(phpxya+phimh*phpxyb)/ &
                          (phpxya**2+phpxyb**2)+1.0d0/ssh)
           end function IcsrCaseC

           function IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
           implicit none
           double precision:: xxbarh,ssh,xk2,epstol,IcsrCaseD
           double precision:: psi,psipxa,psipxb,psimax,x1,x2
           integer:: Nmax
           x1 = -epstol
           x2 = psimax*(1.d0+epstol)
           call root(x1,x2,epstol,Nmax,ssh,xxbarh,psi)
           psipxa = (psi+xxbarh)
           psipxb = psi*(xxbarh+psi/2.d0)
           IcsrCaseD = xk2*(-2.d0*(psipxa+psi*psipxb)/&
                    (psipxa**2+psipxb**2)+1.0d0/ssh)
           end function IcsrCaseD

    !--------------------------------------------------------------------------------------
    !> @brief
    !> This routine computes the root of the function that is evaluated in 
    !> the subroutine 'funcd'. It is based on the subroutine 'root' of 
    !> Numerical Recipes 9.4, which makes use of a Newton-Raphson method 
    !> with root bracketing.  It has been modified to handle the two bracket 
    !> endpoints carefully. The routine searches for a root in the interval
    !> [x1,x2] with a tolerance given by 'xacc', and returns this value
    !> as 'rtsafe'.  The maximum number of iterations allowed is 'maxit'.
    !> C.E.M.
    !--------------------------------------------------------------------------------------
     subroutine root(x1,x2,xacc,maxit,zeta,xxh,rtsafe)
     implicit none
     double precision:: rtsafe,x1,x2,xacc
     double precision:: xxh,zeta
     integer:: j,maxit
     double precision:: df,dx,dxold,f,fh,fl,temp,xh,xl
     call funcd(x1,xxh,zeta,fl,df)
     call funcd(x2,xxh,zeta,fh,df)
     if((fl>0.d0.and.fh>0.d0).or.(fl<0.d0.and.fh<0.d0)) then
           pause 'root must be bracketed in rtsafe'
           write(*,*) 'psimax,fl,fh = ',x2,fl,fh
     endif
     if(dabs(fl)< xacc) then
       rtsafe=x1
       return
     else if(dabs(fh)< xacc) then
       rtsafe=x2
       return
     else if(fl<0.d0) then
       xl=x1
       xh=x2
     else
       xh=x1
       xl=x2
     endif
     rtsafe=0.5d0*(x1+x2)
     dxold=dabs(x2-x1)
     dx=dxold
     call funcd(rtsafe,xxh,zeta,f,df)
     do j=1,maxit
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f)>0.d0.or. &
&         dabs(2.d0*f)>dabs(dxold*df)) then
          dxold=dx
          dx=0.5d0*(xh-xl)
          rtsafe=xl+dx
          if(xl==rtsafe) return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp==rtsafe) return
        endif
        if(abs(dx)<xacc) return
        call funcd(rtsafe,xxh,zeta,f,df)
        if(f<0.d0) then
           xl=rtsafe
        else
           xh=rtsafe
        endif
     enddo
     pause 'root finding exceeding maximum iterations'
     return
     end subroutine

    !--------------------------------------------------------------------------------------
    !> @brief
    !> This routine evaluates the function whose root produces
    !> the retarded angle psi that is required for evaluating
    !> the CSR kernel in Case D of Saldin et al.  The value
    !> of the function is output as 'f', and its derivative
    !> is output as 'derivf'.  C.E.M.
    !--------------------------------------------------------------------------------------
     subroutine funcd(psi,xxh,deltas,f,derivf)
     implicit none
     double precision:: deltas,psi,term1,xxh,gamma,f,derivf
     double precision:: alpha,kappa,tau,theta
     f = psi**4/12+xxh*psi**3/3+psi**2+(2*xxh-2*deltas)*psi-&
               2*deltas*xxh+xxh*xxh
     derivf = psi**3/3+xxh*psi**2+2*psi+2*xxh-2*deltas
     end subroutine

!!! DWA
       subroutine dwakefield_FieldQuant(Nx,Ny,Nz,hx,hy,hz,chgdens,dexwake,&
           deywake,dezwake,dbxwake,dbywake,dbzwake,dtyp,deps,dnkx,dnky,&
           dkx,dky,amp,gammaz,xmin,ymin)
       implicit none
       include 'mpif.h'
       integer, intent(in) :: Nx,Ny,Nz
       double precision, intent(in) :: dtyp,dnkx,dnky,hx,hy,hz,deps,gammaz
       double precision, intent(in) :: xmin,ymin
       double precision, dimension(:,:,:), intent(in) :: chgdens
       double precision, dimension(:), intent(in) :: dkx
       double precision, dimension(:,:,:), intent(in) :: dky,amp
       double precision, dimension(Nx,Ny,Nz), intent(out) :: dexwake,&
                        dezwake,deywake,dbxwake,dbywake,dbzwake
       integer :: i,j,k,l,m,n,p,q,e,par
       double precision :: x,y,z,x0,y0,z0,kx,ky,kz,dv,lx,conv2V,conv2T
       double precision, dimension(:) :: fs(Nz),fa(Nz),lam(Nz)
       double precision, dimension(:) :: lams(Nz),lama(Nz),laa(Nz)
       double precision, dimension(:) :: laas(Nz),laaa(Nz)

       dexwake=0.0
       deywake=0.0
       dezwake=0.0
       dbxwake=0.0
       dbywake=0.0
       dbzwake=0.0
       dv=hx*hy*hz

       if(dtyp.eq.0)then
        lx=PI/dkx(1)
        do e=1,2
         do l=1,nint(dnkx)
          kx=dkx(l)
          do m=1,nint(dnky)
           if(2*((m-1)/2).eq.m-1)then
            par=0
           else
            par=1
           endif
           ky=dky(l,m,e)
           kz=sqrt((kx*kx+ky*ky)/(deps-1)) ! omega/c
           do i=1,Nx
            x=xmin+(i-1)*hx
            do j=1,Ny
             y=ymin+(j-1)*hy
             do k=1,Nz
              z=(k-1)*hz
              fs(k)=amp(l,m,e)*cos(kz*z)
              fa(k)=amp(l,m,e)*sin(kz*z)
              lam(k)=0
              laa(k)=0
              do n=1,Nx
               x0=xmin+(n-1)*hx
               do p=1,Ny
                y0=ymin+(p-1)*hy
                lam(k)=lam(k)+chgdens(n,p,k)*cosh(kx*y0)*cos(kx*x0)
                laa(k)=laa(k)+chgdens(n,p,k)*sinh(kx*y0)*cos(kx*x0)
               enddo ! p
              enddo ! n
             enddo ! k
             if(e.eq.1)then
              call dwakeConv_FieldQuant(Nz,lam,hz,fs,lams)
              call dwakeConv_FieldQuant(Nz,lam,hz,fa,lama)
             else
              call dwakeConv_FieldQuant(Nz,laa,hz,fs,laas)
              call dwakeConv_FieldQuant(Nz,laa,hz,fa,laaa)
             endif
             do k=1,Nz
              if(par.eq.0.and.e.eq.1)then ! LSM-SYM
               dexwake(i,j,k)=dexwake(i,j,k)-lama(k)*sin(kx*x)*cosh(kx*y)&
               *kx/kz
               deywake(i,j,k)=deywake(i,j,k)+lama(k)*cos(kx*x)*sinh(kx*y)&
               *(kx*kx+kz*kz)/kx/kz
               dezwake(i,j,k)=dezwake(i,j,k)+lams(k)*cos(kx*x)*cosh(kx*y)
               dbxwake(i,j,k)=dbxwake(i,j,k)-lama(k)*cos(kx*x)*sinh(kx*y)&
               *kz/kx
               dbzwake(i,j,k)=dbzwake(i,j,k)+lams(k)*sin(kx*x)*sinh(kx*y)
              elseif(par.eq.1.and.e.eq.1)then ! LSE-SYM
               dexwake(i,j,k)=dexwake(i,j,k)+lama(k)*sin(kx*x)*cosh(kx*y)&
               *kz/kx
               dezwake(i,j,k)=dezwake(i,j,k)+lams(k)*cos(kx*x)*cosh(kx*y)
               dbxwake(i,j,k)=dbxwake(i,j,k)+lama(k)*cos(kx*x)*sinh(kx*y)&
               *kx/kz
               dbywake(i,j,k)=dbywake(i,j,k)+lama(k)*sin(kx*x)*cosh(kx*y)&
               *(kx*kx+kz*kz)/kx/kz
               dbzwake(i,j,k)=dbzwake(i,j,k)+lams(k)*sin(kx*x)*sinh(kx*y)
              elseif(par.eq.0.and.e.eq.2)then ! LSM-ASYM
               dexwake(i,j,k)=dexwake(i,j,k)-laaa(k)*sin(kx*x)*sinh(kx*y)&
               *kx/kz
               deywake(i,j,k)=deywake(i,j,k)+laaa(k)*cos(kx*x)*cosh(kx*y)&
               *(kx*kx+kz*kz)/kx/kz
               dezwake(i,j,k)=dezwake(i,j,k)+laas(k)*cos(kx*x)*sinh(kx*y)
               dbxwake(i,j,k)=dbxwake(i,j,k)-laaa(k)*cos(kx*x)*cosh(kx*y)&
               *kz/kx
               dbzwake(i,j,k)=dbzwake(i,j,k)+laas(k)*sin(kx*x)*cosh(kx*y)
              elseif(par.eq.1.and.e.eq.2)then ! LSE-ASYM
               dexwake(i,j,k)=dexwake(i,j,k)+laaa(k)*sin(kx*x)*sinh(kx*y)&
               *kz/kx
               dezwake(i,j,k)=dezwake(i,j,k)+laas(k)*cos(kx*x)*sinh(kx*y)
               dbxwake(i,j,k)=dbxwake(i,j,k)+laaa(k)*cos(kx*x)*cosh(kx*y)&
               *kx/kz
               dbywake(i,j,k)=dbywake(i,j,k)+laaa(k)*sin(kx*x)*sinh(kx*y)&
               *(kx*kx+kz*kz)/kx/kz
               dbzwake(i,j,k)=dbzwake(i,j,k)+laas(k)*sin(kx*x)*cosh(kx*y)
              endif
             enddo ! k (z)
            enddo ! j (y)
           enddo ! i (x)
          enddo ! m
         enddo ! l
        enddo ! e
       
       ! take care of units: V/m and T
        conv2V=dv*gammaz/lx/hz
        dexwake=conv2V*dexwake
        deywake=conv2V*deywake
        dezwake=-conv2V*dezwake
        conv2T=conv2V/Clight
        dbxwake=conv2T*dbxwake
        dbywake=conv2T*dbywake
        dbzwake=-conv2T*dbzwake
       else
        lama=0.0
        lam=0.0
        do n=1,Nx
         do p=1,Ny
          do k=1,Nz
           lam(k)=lam(k)+chgdens(n,p,k)
          enddo
         enddo ! p
        enddo ! n

        do l=1,nint(dnkx)+1
         do m=1,nint(dnky)
          ky=dky(l,m,1)
          kz=ky/sqrt(deps-1) ! omega/c
          do k=1,Nz
           z=(k-1)*hz
           fs(k)=amp(l,m,1)*cos(kz*z)
          enddo ! k
          call dwakeConv_FieldQuant(Nz,lam,hz,fs,lams)
          do k=1,Nz
           lama(k)=lama(k)+lams(k)
          enddo
         enddo ! m
        enddo ! l
 
        do i=1,Nx
         do j=1,Ny
          do k=1,Nz
           dezwake(i,j,k)=lama(k)
          enddo
         enddo
        enddo
       
        conv2V=dv*gammaz/deps/PI/hz
        dezwake=-conv2V*dezwake
       endif ! dtype

       end subroutine dwakefield_FieldQuant


       ! make convolution for d-wakefields
       subroutine dwakeConv_FieldQuant(Nz,weightz,&
                          hz,kern,lam)
       implicit none
       include 'mpif.h'
       integer, intent(in) :: Nz
       double precision, intent(in) :: hz
       double precision, dimension(Nz), intent(in) :: weightz,kern
       double precision, dimension(Nz), intent(out) :: lam
       double precision, dimension(2*Nz,1) :: densz2n,densz2nout,&
                         greenwake,greenwakeout
       integer :: kz,twonz,one,ksign,kkzz,i,iz,iz1
       double precision :: scale,zz,zziz

       do kz = 1, Nz
          densz2n(kz,1) = weightz(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       !longitudinal wakefield function
       greenwake(1,1) = 0.5*kern(1)
       do kz = 2, Nz+1
         zz = (kz-1)*hz
         iz = zz/hz + 1
         iz1 = iz + 1
         if(iz1.gt.Nz) then
           iz = Nz - 1
           iz1 = Nz
         endif
         zziz = (iz-1)*hz
         greenwake(kz,1) = kern(iz)+(kern(iz1)-kern(iz))*(zz-zziz)/hz
       enddo

       do kz = Nz+2, twonz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 2, Nz
         greenwake(kz,1) = 0.0
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

       !the "-" sign is from definition of longitudinal wake function, see Wangler's book
       do kz = 1, Nz
         lam(kz) = greenwakeout(kz,1)*hz
       enddo

       end subroutine dwakeConv_FieldQuant

      subroutine dispEqSlab_FieldQuant(a,b,eps,dlx,nkx,nky,dkx,dky)
      implicit none
      include 'mpif.h'
      double precision, intent(in) :: a,b,eps,dlx,nkx,nky
      double precision, dimension(:), intent(out) :: dkx
      double precision, dimension(:,:,:), intent(out) :: dky
      integer :: i,j,k,nx,ny,parity
      double precision :: dx,c1,c2,c3,x1,x2,sm,cx,f1,cf

      nx=nkx
      ny=nky
      dx=PI/2.0
      sm=1.0e-6
      do i=1,nx
        dkx(i)=(2*i-1)*PI/dlx
!! solve eqn.: c1*cot(x)=c2*x (LSM)
!! solve eqn.: c1*cot(x)=-c3/x (LSE)
        do k=1,2
          if(k.eq.1)then
            c1=1.0/tanh(dkx(i)*a); ! monopoles
          else
            c1=tanh(dkx(i)*a); ! dipoles
          endif
          c2=1.0/dkx(i)/(b-a)/eps;
          c3=dkx(i)*(b-a);
          do j=1,ny
            x1=(j-1)*dx
            x2=j*dx
            x1=x1+1e-10
            x2=x2-1e-10
            parity=int((j-1)/2);
            if(parity*2.eq.(j-1))then
              parity=0 ! LSM
            else
              parity=1 ! LSE
            endif
21          if(abs(x2-x1).gt.sm)then
              cx=(x1+x2)/2.
              if(parity.eq.0)then
                f1=func1(x1,c1,c2)
                cf=func1(cx,c1,c2)
              else
                f1=func2(x1,c1,c3)
                cf=func2(cx,c1,c3)     
              endif
              if(f1*cf.lt.0)then
               x2=cx
              else
               x1=cx
              endif
              goto 21
            endif ! abs...
            dky(i,j,k)=cx/(b-a) ! ky
!            dky(i,j,k)=1.0/sqrt(eps-1)*sqrt(cx*cx/(b-a)/(b-a)+dkx(i)*dkx(i))
          enddo ! j
        enddo ! k
      enddo ! i

      end subroutine dispEqSlab_FieldQuant

      subroutine ampSlab_FieldQuant(a,b,eps,nkx,nky,dkx,dky,amp)
      implicit none
      include 'mpif.h'
      double precision, intent(in) :: a,b,eps,nkx,nky
      double precision, dimension(:), intent(in) :: dkx
      double precision, dimension(:,:,:), intent(in) :: dky
      double precision, dimension(:,:,:), intent(out) :: amp
      integer :: i,j,k,nx,ny,parity
!      double precision :: ampLSM_SYM,ampLSM_ASYM,ampLSE_SYM,ampLSE_ASYM
      
      nx=nkx
      ny=nky
      do i=1,nx
        do k=1,2
          do j=1,ny
            parity=int((j-1)/2);
            if(parity*2-(j-1).eq.0)then
              ! LSM
              if(k.eq.1)then
              ! monopole
                amp(i,j,k)=ampLSM_SYM(a,b,eps,dkx(i),dky(i,j,k))
              else
              ! dipole
                amp(i,j,k)=ampLSM_ASYM(a,b,eps,dkx(i),dky(i,j,k))
              endif
            else
              ! LSE
              if(k.eq.1)then
              ! monopole
                amp(i,j,k)=ampLSE_SYM(a,b,eps,dkx(i),dky(i,j,k))
              else
              ! dipole
                amp(i,j,k)=ampLSE_ASYM(a,b,eps,dkx(i),dky(i,j,k))
              endif
            endif
          enddo
        enddo
      enddo

      end subroutine ampSlab_FieldQuant

      subroutine dispEqCyl_FieldQuant(a,b,eps,nkx,nky,dky)
      implicit none
      include 'mpif.h'
      double precision, intent(in) :: a,b,eps,nkx,nky
      double precision, dimension(:,:,:), intent(out) :: dky
      integer :: i,j,found
      double precision :: nbig,nbbig,s_start,s_step,ds,actual_s
      !double precision :: x1,x2,x,sg,calcLFS,BESJ0
      double precision :: x1,x2,x,sg

      sg=2*PI/a
      nbig=1e3
      nbbig=1e10
      s_start=sg/nbig
      s_step=sg/nbig
      ds=s_step/nbbig
      actual_s=s_start
      found=0
      do i=1,nint(nkx)+1
       if(i.eq.1)then ! monopole case (m=0)
        do while(found.lt.nky)
         if(calcLFS(actual_s,a,b,eps)*calcLFS(actual_s+s_step,a,b,eps).lt.0)then
          x1=actual_s
          x2=actual_s+s_step
          do while(x2-x1.gt.ds)
           x=x1+(x2-x1)/2.
           if(calcLFS(x1,a,b,eps)*calcLFS(x,a,b,eps).lt.0)then
            x2=x
           else
            x1=x
           endif
          enddo ! while "x2-x1"
          found=found+1
          dky(i,found,1)=x1+(x2-x1)/2
         endif
         actual_s=actual_s+s_step
        enddo ! while "found"
       endif ! i=1
      enddo ! i

      end subroutine dispEqCyl_FieldQuant

      subroutine ampCyl_FieldQuant(a,b,eps,nkx,nky,dky,amp)
      implicit none
      include 'mpif.h'
      double precision, intent(in) :: a,b,eps,nkx,nky
      double precision, dimension(:,:,:), intent(in) :: dky
      double precision, dimension(:,:,:), intent(out) :: amp
      integer :: i,j,k,nx,ny
      !double precision :: x,dx,deriv1,deriv2,deriv,nbig,calcLFS
      double precision :: x,dx,deriv1,deriv2,deriv,nbig
!      double precision :: calcR0

      nbig=1e3
      do i=1,nint(nkx)+1
       if(i.eq.1)then
        do j=1,nint(nky)
         x=dky(i,j,1)
         dx=x/nbig
         deriv1=(calcLFS(x+dx,a,b,eps)-calcLFS(x,a,b,eps))/dx
         deriv2=-(calcLFS(x-dx,a,b,eps)-calcLFS(x,a,b,eps))/dx
         deriv=(deriv1+deriv2)/2.0
         amp(i,j,1)=calcR0(x,a,b)/deriv/a
        enddo
       endif
      enddo
      end subroutine ampCyl_FieldQuant

      double precision function func1(x,c1,c2)
      double precision x,c1,c2
        func1=c1/tan(x)-c2*x
      return
      end function func1

      double precision function func2(x,c1,c3)
      double precision x,c1,c3
        func2=c1/tan(x)+c3/x
      return
      end function func2

      double precision function ampLSM_SYM(a,b,eps,kx,ky)
      double precision a,b,eps,kx,ky,t1,t2,t3
        t1=sinh(2*kx*a)/2/kx
        t2=eps*cosh(kx*a)*cosh(kx*a)/sin(ky*(b-a))/sin(ky*(b-a))
        t3=(b-a)*(1.0+eps*kx*kx/ky/ky)/2.
        t3=t3-sin(2*ky*(b-a))/4./ky*(1.0-eps*kx*kx/ky/ky)
        t2=t2*t3
        ampLSM_SYM=1.0/(t1+t2)
      return
      end function ampLSM_SYM

      double precision function ampLSM_ASYM(a,b,eps,kx,ky)
      double precision a,b,eps,kx,ky,t1,t2,t3
        t1=cosh(2*kx*a)/2/kx
        t2=eps*sinh(kx*a)*sinh(kx*a)/sin(ky*(b-a))/sin(ky*(b-a))
        t3=(b-a)*(1.0+eps*kx*kx/ky/ky)/2.
        t3=t3-sin(2*ky*(b-a))/4./ky*(1.0-eps*kx*kx/ky/ky)
        t2=t2*t3
        ampLSM_ASYM=1.0/(t1+t2)
      return
      end function ampLSM_ASYM

      double precision function ampLSE_SYM(a,b,eps,kx,ky)
      double precision a,b,eps,kx,ky,t1,t2,t3
        t1=sinh(2*kx*a)/2/kx
        t2=cosh(kx*a)*cosh(kx*a)/sin(ky*(b-a))/sin(ky*(b-a))
        t3=(b-a)/2.*(eps+ky*ky/kx/kx)
        t3=t3-sin(2*ky*(b-a))/4./ky*(eps-ky*ky/kx/kx)
        t2=t2*t3
        ampLSE_SYM=1.0/(t1+t2)
      return
      end function ampLSE_SYM

      double precision function ampLSE_ASYM(a,b,eps,kx,ky)
      double precision a,b,eps,kx,ky,t1,t2,t3
        t1=cosh(2*kx*a)/2/kx
        t2=sinh(kx*a)*sinh(kx*a)/sin(ky*(b-a))/sin(ky*(b-a))
        t3=(b-a)/2.*(eps+ky*ky/kx/kx)
        t3=t3-sin(2*ky*(b-a))/4./ky*(eps-ky*ky/kx/kx)
        t2=t2*t3
        ampLSE_ASYM=1.0/(t1+t2)
      return
      end function ampLSE_ASYM

      double precision function calcLFS(s,a,b,eps)
      double precision s,a,b,eps,calcR0P,calcR0
      real sb,da,jderiv1,yderiv1,jderiv2
      real yderiv2,jderiv,yderiv,t1,t2
      sb=1000.
      t1=s*b
      t2=s*a
      da=s*a/sb
      jderiv1=(BESJ0(t2+da)-BESJ0(t2))/da
      yderiv1=(BESY0(t2+da)-BESY0(t2))/da
      jderiv2=-(BESJ0(t2-da)-BESJ0(t2))/da
      yderiv2=-(BESY0(t2-da)-BESY0(t2))/da
      jderiv=(jderiv1+jderiv2)/2
      yderiv=(yderiv1+yderiv2)/2
      calcR0P=BESY0(t1)*jderiv-BESJ0(t1)*yderiv
      calcR0=BESY0(t1)*BESJ0(t2)-BESJ0(t1)*BESY0(t2)
      calcLFS=calcR0P+s*a*calcR0/eps/2.0
      return
      end function calcLFS

      double precision function calcR0(s,a,b)
      double precision s,a,b
      real t1,t2
      t1=s*b
      t2=s*a
      calcR0=BESY0(t1)*BESJ0(t2)-BESJ0(t1)*BESY0(t2)
      return
      end function calcR0

!!! end DWA

      end module FieldQuantclass