!--------------------------------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! BeamBunchclass: Charged beam bunch class in Beam module of APPLICATION layer.
! 
! MODULE    : ... BeamBunchclass
! VERSION   : ... 1.0
!> @author
!> Ji Qiang
!
! DESCRIPTION: 
!> This class defines the charged particle beam bunch information in the accelerator.
! Comments: 1) 
!--------------------------------------------------------------------------------------
      module BeamBunchclass
        use Timerclass
        use PhysConstclass
        use BeamLineElemclass
        use CompDomclass
        use Pgrid2dclass
        use Fldmgerclass
        use Dataclass
        type BeamBunch
          !private
          double precision :: Current !< beam current
          double precision :: Mass !< part. mass
          double precision :: Charge !< charge
          integer :: Npt !< num of total global macroparticles
          integer :: Nptlocal !< num of total local particles
          double precision, pointer, dimension(:,:) :: Pts1 !< particles type one
          double precision, dimension(6) :: refptcl !< reference particle
        end type BeamBunch
      contains
        !--------------------------------------------------------------------------------------
        !> @author Ji Qiang
        !> @brief
        !> Initialize Beambunch class.
        !> @param[in] incurr, inkin, inmass, incharge, innp, phasini 
        !--------------------------------------------------------------------------------------
        subroutine construct_BeamBunch(this,incurr,inkin,inmass,incharge,innp,&
                                       phasini)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: incurr,inkin,inmass,&
                                        incharge,phasini
        integer, intent(in) :: innp
        integer :: myid, myidx, myidy,comm2d,commrow,commcol,ierr
        integer :: nptot,nprocrow,nproccol
   
        this%Current = incurr
        this%Mass = inmass
        this%Charge = incharge
        this%Npt = innp

        this%refptcl = 0.0d0
        this%refptcl(5) = phasini
        this%refptcl(6) = -(inkin/this%Mass + 1.0d0)

        end subroutine construct_BeamBunch

        !--------------------------------------------------------------------------------------
        !> @author Ji Qiang
        !> @brief
        !> Set local # of particles.
        !> @param[in] innpt 
        !--------------------------------------------------------------------------------------
        subroutine setnpt_BeamBunch(this,innpt)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innpt
        type (BeamBunch), intent(inout) :: this

        this%Nptlocal = innpt

        end subroutine setnpt_BeamBunch

        !--------------------------------------------------------------------------------------
        !> @author Ji Qiang
        !> @brief
        !> Get local # of particles.
        !> @param[out] outnpt 
        !--------------------------------------------------------------------------------------
        subroutine getnpt_BeamBunch(this,outnpt)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        integer, intent(out) :: outnpt

        outnpt = this%Nptlocal

        end subroutine getnpt_BeamBunch

        !--------------------------------------------------------------------------------------
        !> @author Ji Qiang
        !> @brief
        !> Drift half step in positions.
        !> Here, x, y, z are normalized by C * Dt
        !> tau - normalized step size (by Dt).
        !> Only particle with z > 0 is drifted.
        !> @param[in] tau, betazini 
        !> @param[inout] t
        !--------------------------------------------------------------------------------------
        subroutine drifthalf_BeamBunch(this,t,tau,betazini)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: t
        double precision, intent (in) :: tau,betazini
        double precision :: t0,recpgam
        integer :: i

        call starttime_Timer(t0)

        do i = 1, this%Nptlocal
          ! get 1.0d0/gamma of each particle
          recpgam = 1.0d0/sqrt(1.0d0+this%Pts1(2,i)**2+this%Pts1(4,i)**2+&
                                 this%Pts1(6,i)**2)
          if(this%Pts1(5,i).gt.0.0) then
            this%Pts1(1,i) = this%Pts1(1,i)+0.5d0*tau*this%Pts1(2,i)*recpgam
            this%Pts1(3,i) = this%Pts1(3,i)+0.5d0*tau*this%Pts1(4,i)*recpgam
            this%Pts1(5,i) = this%Pts1(5,i)+0.5d0*tau*this%Pts1(6,i)*recpgam
          endif
        enddo

        t_map1 = t_map1 + elapsedtime_Timer(t0)

        end subroutine drifthalf_BeamBunch

        !--------------------------------------------------------------------------------------
        !> @author Ji Qiang
        !> @brief
        !> Particle emission
        !> For particle with z < 0, they are just shifted long z
        !> This is used to simulate the process of emission from photocathod
        !> @param[inout] t
        !> @param[in] tau, betazini 
        !-------------------------------------------------------------------------------------- 
        subroutine driftemission_BeamBunch(this,t,tau,betazini)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: t
        double precision, intent (in) :: tau,betazini
        double precision :: t0,recpgam
        integer :: i

        call starttime_Timer(t0)

        do i = 1, this%Nptlocal
          if(this%Pts1(5,i).le.0.0d0 .and. this%Pts1(6,i).ge.0.0d0) then
            this%Pts1(5,i) = this%Pts1(5,i)+tau*betazini
          endif
        enddo

        t_map1 = t_map1 + elapsedtime_Timer(t0)

        end subroutine driftemission_BeamBunch

        !--------------------------------------------------------------------------------------
        !> @author Ji Qiang
        !> @brief
        !> Drift half step in positions.
        !> Here, x, y, z are normalized by C * Dt
        !> tau - normalized step size (by Dt).
        !> @param[inout] t
        !> @param[in] tau
        !-------------------------------------------------------------------------------------- 
        subroutine drifthalforg_BeamBunch(this,t,tau)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: t
        double precision, intent (in) :: tau
        double precision :: t0,recpgam
        integer :: i

        call starttime_Timer(t0)

        do i = 1, this%Nptlocal
          !//get 1.0d0/gamma of each particle
          recpgam = 1.0d0/sqrt(1.0d0+this%Pts1(2,i)**2+this%Pts1(4,i)**2+&
                                 this%Pts1(6,i)**2)
          this%Pts1(1,i) = this%Pts1(1,i)+0.5d0*tau*this%Pts1(2,i)*recpgam
          this%Pts1(3,i) = this%Pts1(3,i)+0.5d0*tau*this%Pts1(4,i)*recpgam
          this%Pts1(5,i) = this%Pts1(5,i)+0.5d0*tau*this%Pts1(6,i)*recpgam
        enddo

        t_map1 = t_map1 + elapsedtime_Timer(t0)

        end subroutine drifthalforg_BeamBunch

        !drift all particles along dz
        subroutine driftz_BeamBunch(this,dz)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent (in) :: dz
        integer :: i

        do i = 1, this%Nptlocal
          this%Pts1(5,i) = this%Pts1(5,i) + dz
        enddo
        
        end subroutine driftz_BeamBunch

        ! Here, all indices of potential are local to processor.
        ! Advance the particles in the velocity space using the force
        ! from the external field and the self space charge force
        ! interpolated from the grid to particles. (Lorentz force)
        subroutine kick1t_BeamBunch(this,beamelem,zbeamelem,idrfile,nbeamln,t,tau,innx,&
              inny,innz,temppotent,ptsgeom,grid,Flagbc,flagerr,gammaz,&
              ibinit,ibend)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (BeamLineElem), dimension(:), intent(in) :: beamelem
        double precision, dimension(:,:), intent(in) :: zbeamelem
        integer, dimension(:,:), intent(in) :: idrfile
        integer, intent(in) :: innx, inny, innz, Flagbc,flagerr,nbeamln,&
                               ibinit,ibend
        type (CompDom), intent(in) :: ptsgeom
        double precision,dimension(innx,inny,innz),intent(inout) :: temppotent
        type (Pgrid2d), intent(in) :: grid
        double precision, intent(in) :: t,tau,gammaz
        double precision, dimension(innx,inny,innz) :: egx,egy,egz
        double precision, dimension(3) :: msize
        double precision :: hxi, hyi, hzi
        double precision :: t0,tg,chge,mass,curr,gam
        integer :: totnp,nproccol,nprocrow,myid,myidx,myidy
        integer :: i, j, k, yadd, zadd, innp
!        integer :: comm2d,commcol,commrow

        call starttime_Timer(t0)

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

        curr = this%Current
        innp = this%Nptlocal
        tg = this%refptcl(5)
        gam = -this%refptcl(6)
        mass = this%Mass
        chge = this%Charge

        if(curr.gt.0.0) then
!------------------------------------------------------------------
! current greater than 0
    
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
            egx(1,j,k) = hxi*(temppotent(1,j,k)-temppotent(2,j,k))
            do i = 2, innx-1
              egx(i,j,k) = 0.5d0*hxi*(temppotent(i-1,j,k)- &
                           temppotent(i+1,j,k))
            enddo
            egx(innx,j,k) = hxi*(temppotent(innx-1,j,k)- &
                                 temppotent(innx,j,k))
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
                egy(i,yadd+1,k) = hyi*(temppotent(i,yadd+1,k)- &
                                       temppotent(i,yadd+2,k))
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
                egy(i,inny-yadd,k) = hyi*(temppotent(i,inny-yadd-1,k)- &
                                          temppotent(i,inny-yadd,k))
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
              egy(i,1,k) = hyi*(temppotent(i,1,k)-temppotent(i,2,k))
            enddo
            do j = 2, inny-1
              do i = 1, innx
                egy(i,j,k) = 0.5d0*hyi*(temppotent(i,j-1,k)- &
                             temppotent(i,j+1,k))
              enddo
            enddo
            do i = 1, innx
              egy(i,inny,k) = hyi*(temppotent(i,inny-1,k)- &
                                   temppotent(i,inny,k))
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
                  egz(i,j,zadd+1) = hzi*(temppotent(i,j,zadd+1)- &
                                         temppotent(i,j,zadd+2))
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
                  egz(i,j,innz-zadd) = hzi*(temppotent(i,j,innz-zadd-1)- &
                                            temppotent(i,j,innz-zadd))
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
                egz(i,j,1) = hzi*(temppotent(i,j,1)-temppotent(i,j,2))
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
                egz(i,j,innz) = hzi*(temppotent(i,j,innz-1)- &
                                     temppotent(i,j,innz))
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

        if((Flagbc.eq.3).or.(Flagbc.eq.4)) then ! round pipe
!          if(flagerr.eq.1) then
!            call scatter2rerr_BeamBunch(innp,innx,inny,innz,this%Pts1,egx,egy,egz,&
!            ptsgeom,nprocrow,nproccol,myidx,myidy,tg,gam,curr,chge,mass,tau,z,&
!            beamelem)
!          else
!            call scatter2r_BeamBunch(innp,innx,inny,innz,this%Pts1,egx,egy,egz,&
!            ptsgeom,nprocrow,nproccol,myidx,myidy,tg,gam,curr,chge,mass,tau,z,&
!            beamelem)
!          endif
        else
!          if(flagerr.eq.1) then
!            call scatter2err_BeamBunch(innp,innx,inny,innz,this%Pts1,egx,egy,egz,&
!            ptsgeom,nprocrow,nproccol,myidx,myidy,tg,gam,curr,chge,mass,tau,z,&
!            beamelem)
!          else
            call scatter2t_BeamBunch(innp,innx,inny,innz,this%Pts1,egx,egy,egz,&
            ptsgeom,nprocrow,nproccol,myidx,myidy,t,gammaz,chge,mass,tau,&
            beamelem,zbeamelem,idrfile,nbeamln,ibinit,ibend)
!          endif
        endif

        else
!------------------------------------------------------------------
! current is 0
          if(flagerr.eq.1) then
!            call scatter20err_BeamBunch(innp,this%Pts1,tg,gam,chge,mass,tau,z,&
!                                     beamelem)
          else
!            call scatter20t_BeamBunch(innp,this%Pts1,t,chge,mass,tau,&
!                              beamelem,zbeamelem,idrfile,nbeamln,ibinit,ibend)
          endif
        endif

        this%refptcl(6) = -gam

        !call MPI_BARRIER(comm2d,ierr)
        !if(myid.eq.0) then
        !  print*,"after kick:"
        !endif

        t_force = t_force + elapsedtime_Timer(t0)

        end subroutine kick1t_BeamBunch

        subroutine scatter2t_BeamBunch(innp,innx,inny,innz,rays,exg,&
        eyg,ezg,ptsgeom,npx,npy,myidx,myidy,tg,gammaz,chge,mass,dt,&
        beamelem,zbeamelem,idrfile,nbeamln,ibinit,ibend)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy,nbeamln
        double precision, intent (inout), dimension(6,innp) :: rays
        double precision, intent (in), dimension(innx,inny,innz) :: exg,eyg,ezg 
        double precision, intent (in) :: dt,mass,chge,tg
        double precision, intent (in) :: gammaz
        type (BeamLineElem), dimension(:), intent(in) :: beamelem
        double precision, dimension(:,:), intent(in) :: zbeamelem
        integer, dimension(:,:), intent(in) :: idrfile
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab, cd, ef
        integer :: n,i,ii,ibinit,ibend
        double precision :: t0
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin,zmax
        double precision, dimension(3) :: msize
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: range, extfld,tmpfld
        type (CompDom) :: ptsgeom
        integer :: ix,jx,kx,ix1,jx1,kx1,ierr,kadd,jadd
        double precision :: beta0,qmcc,recpgamma,betC,coefLz,&
                            umx,umy,umz,upx,upy,upz,tmp,a1,a2,a3,a4,s1,s2,&
                            s3,exn,eyn,ezn,ex,ey,ez,bx,by,bz,zz
        integer, dimension(Maxoverlap) :: idbeamln
        integer :: noverlap,nnn

        call starttime_Timer( t0 )

        beta0 = sqrt(gammaz**2-1.0d0)/gammaz
        betC = beta0/Clight
        qmcc = chge/mass
        coefLz = qmcc*Scxlt
    
        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0d0/msize(1)
        hyi = 1.0d0/msize(2)
        hzi = 1.0d0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)
        zmax = range(6)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        idbeamln = 1

        do n = 1, innp
          if(rays(5,n)>0) then !//only particles outside the emission plane get accelerated

          ix=(rays(1,n)-xmin)*hxi + 1
          ab=((xmin-rays(1,n))+ix*hx)*hxi
          jx=(rays(3,n)-ymin)*hyi + 1 + jadd
          cd=((ymin-rays(3,n))+(jx-jadd)*hy)*hyi
          kx=(rays(5,n)-zmin)*hzi + 1 + kadd
          ef=((zmin-rays(5,n))+(kx-kadd)*hz)*hzi

          ix1 = ix + 1
          jx1 = jx + 1
          kx1 = kx + 1

          exn = (exg(ix,jx,kx)*ab*cd*ef  &
                  +exg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +exg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +exg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +exg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +exg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +exg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +exg(ix1,jx,kx)*(1.0d0-ab)*cd*ef)

          eyn = (eyg(ix,jx,kx)*ab*cd*ef  &
                  +eyg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +eyg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +eyg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +eyg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +eyg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +eyg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +eyg(ix1,jx,kx)*(1.0d0-ab)*cd*ef) 

          ezn = ezg(ix,jx,kx)*ab*cd*ef  &
                  +ezg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +ezg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +ezg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +ezg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +ezg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +ezg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +ezg(ix1,jx,kx)*(1.0d0-ab)*cd*ef

          !find which element the particle belongs to.
          zz = rays(5,n)*Scxlt
          !counter how many elements overlaped at given location
          noverlap = 0
          do i = ibinit, ibend
            if( (zz.ge.zbeamelem(1,i)) .and. &
                (zz.le.zbeamelem(2,i)) ) then
               noverlap = noverlap + 1
               idbeamln(noverlap) = i
!               exit
            endif
          enddo
 
          pos(1) = rays(1,n)*Scxlt
          pos(2) = rays(3,n)*Scxlt
          pos(3) = zz
          pos(4) = tg
          !get external field from all overlaped fields at one location
          extfld = 0.0
          !print*,"noverlap: ",noverlap
          tmpfld = 0.0
          do nnn = 1, noverlap
            !call getfldt_BeamLineElem(beamelem(idbeamln(nnn)),pos,&
            !tmpfld,idrfile(3,idbeamln(nnn)))
            extfld = extfld + tmpfld
          enddo

          ex = gammaz*exn+extfld(1)
          ey = gammaz*eyn+extfld(2)
          ez = ezn+extfld(3)
          bx = -gammaz*eyn*betC+extfld(4)
          by = gammaz*exn*betC+extfld(5)
          bz = extfld(6)

!The following is commented out because the exn, eyn might include
!transverse wakefield.
          !here, the space-charge B field is lumped into the E field
!          ex = exn/gammaz+extfld(1)
!          ey = eyn/gammaz+extfld(2)
!          ez = ezn+extfld(3)
!          bx = extfld(4)
!          by = extfld(5)
!          bz = extfld(6)

          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefLz*ex*0.5d0*dt
          umy = rays(4,n) + coefLz*ey*0.5d0*dt
          umz = rays(6,n) + coefLz*ez*0.5d0*dt
          recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
          !recpgamma = 1.0d0/gammaz
          tmp = 0.5d0*Clight*dt*recpgamma*coefLz
          a1 = tmp*bx
          a2 = tmp*by
          a3 = tmp*bz
          a4 = 1.0d0+a1*a1+a2*a2+a3*a3
          s1 = umx + tmp*(umy*bz-umz*by)
          s2 = umy - tmp*(umx*bz-umz*bx)
          s3 = umz + tmp*(umx*by-umy*bx)
          upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
          upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
          rays(2,n) = upx + coefLz*ex*0.5d0*dt
          rays(4,n) = upy + coefLz*ey*0.5d0*dt
          rays(6,n) = upz + coefLz*ez*0.5d0*dt

          endif
        enddo

        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine scatter2t_BeamBunch

        subroutine scatter20t_BeamBunch(innp,rays,&
        tg,chge,mass,dt,beamelem,zbeamelem,idrfile,nbeamln,ibinit,&
        ibend,fldmap,flagerr)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,nbeamln,ibinit,ibend,flagerr
        double precision, intent (inout), dimension(6,innp) :: rays
        double precision, intent (in) :: dt,mass,chge,tg
        type (BeamLineElem), dimension(:), intent(in) :: beamelem
        double precision, dimension(:,:), intent(in) :: zbeamelem
        type (fielddata), dimension(:), intent(in) :: fldmap
        integer, dimension(:,:), intent(in) :: idrfile
        integer :: n,i
        double precision :: t0
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: range, extfld,tmpfld
        double precision :: qmcc,recpgamma,coefLz,&
                            umx,umy,umz,upx,upy,upz,tmp,a1,a2,a3,a4,s1,s2,&
                            s3,ex,ey,ez,bx,by,bz,zz
        integer, dimension(Maxoverlap) :: idbeamln
        integer :: noverlap,nnn

        call starttime_Timer( t0 )

        qmcc = chge/mass
        coefLz = qmcc*Scxlt
    
!        print*,"coefLz: ",coefLz,qmcc,Scxlt,dt,tg
        idbeamln = 1
        do n = 1, innp
          !find which element the particle belongs to.
          zz = rays(5,n)*Scxlt
          !counter how many elements overlaped at given location
          noverlap = 0
          do i = ibinit, ibend
            if( (zz.ge.zbeamelem(1,i)) .and. &
                (zz.le.zbeamelem(2,i)) ) then
               noverlap = noverlap + 1
               idbeamln(noverlap) = i
!               exit
            endif
          enddo
 
          pos(1) = rays(1,n)*Scxlt
          pos(2) = rays(3,n)*Scxlt
          pos(3) = zz
          pos(4) = tg
          !get external field from all overlaped fields at one location
          extfld = 0.0
!          print*,"noverlap: ",noverlap
          if(flagerr.eq.1) then
            do nnn = 1, noverlap
              call getflderrt_BeamLineElem(beamelem(idbeamln(nnn)),pos,&
              tmpfld,fldmap(idrfile(3,idbeamln(nnn))))
              extfld = extfld + tmpfld
            enddo
          else
            do nnn = 1, noverlap
              call getfldt_BeamLineElem(beamelem(idbeamln(nnn)),pos,&
              tmpfld,fldmap(idrfile(3,idbeamln(nnn))))
              extfld = extfld + tmpfld
!            if(nnn.eq.2) then
!              write(16,*)zz,tmpfld(4),tmpfld(5),tmpfld(6)
!            endif
            enddo
          endif

!          print*,"pos...",idbeamln,idrfile(3,idbeamln),pos

!          extfld = 0.0
          ex = extfld(1)
          ey = extfld(2)
          ez = extfld(3)
          bx = extfld(4)
          by = extfld(5)
          bz = extfld(6)
!          ez = 1.0e6

!          print*,pos(3),ex,ey,ez,bx,by,bz,pos(4)
!          write(12,*)pos(3),ex,ey,ez,bx,by,bz,pos(4)
!          call flush_(12)
          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefLz*ex*0.5d0*dt
          umy = rays(4,n) + coefLz*ey*0.5d0*dt
          umz = rays(6,n) + coefLz*ez*0.5d0*dt
          recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
          tmp = 0.5d0*Clight*dt*recpgamma*coefLz
          a1 = tmp*bx
          a2 = tmp*by
          a3 = tmp*bz
          a4 = 1.0d0+a1*a1+a2*a2+a3*a3
          s1 = umx + tmp*(umy*bz-umz*by)
          s2 = umy - tmp*(umx*bz-umz*bx)
          s3 = umz + tmp*(umx*by-umy*bx)
          upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
          upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
          rays(2,n) = upx + coefLz*ex*0.5d0*dt
          rays(4,n) = upy + coefLz*ey*0.5d0*dt
          rays(6,n) = upz + coefLz*ez*0.5d0*dt
        enddo
!        print*,"rays6: ",rays(6,1)

        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine scatter20t_BeamBunch

        !implicit Boris type integrator
        subroutine kick2t_BeamBunch(innp,innx,inny,innz,rays,exg,&
        eyg,ezg,bxg,byg,bzg,ptsgeom,npx,npy,myidx,myidy,tg,&
        chge,mass,dt,beamelem,zbeamelem,idrfile,nbeamln,ibinit,ibend,&
        fldmap,flagerr)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy,&
                              nbeamln,flagerr
        double precision, intent (inout), dimension(6,innp) :: rays
        double precision, intent (in), dimension(innx,inny,innz) :: exg,eyg,ezg 
        double precision, intent (in), dimension(innx,inny,innz) :: bxg,byg,bzg 
        double precision, intent (in) :: dt,mass,chge,tg
        type (BeamLineElem), dimension(:), intent(in) :: beamelem
        double precision, dimension(:,:), intent(in) :: zbeamelem
        integer, dimension(:,:), intent(in) :: idrfile
        type (fielddata), dimension(:), intent(in) :: fldmap
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab, cd, ef
        integer :: n,i,ii,ibinit,ibend
        double precision :: t0
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin,zmax
        double precision, dimension(3) :: msize
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: range, extfld
        type (CompDom) :: ptsgeom
        integer :: ix,jx,kx,ix1,jx1,kx1,ierr,kadd,jadd
        double precision :: qmcc,recpgamma,coefLz,&
                            umx,umy,umz,upx,upy,upz,tmp,a1,a2,a3,a4,s1,s2,&
                            s3,exn,eyn,ezn,ex,ey,ez,bx,by,bz,zz,bxn,byn
        !for residence correction purpose
        double precision :: dev,frac,frac1,stepsize
        double precision, dimension(6) :: extfld1,extfld2,tmpfld
        double precision, dimension(4) :: pos2
        integer, dimension(Maxoverlap) :: idbeamln 
        integer :: noverlap,nnn
        integer :: ntmp1,ntmp2

        call starttime_Timer( t0 )

        qmcc = chge/mass
        coefLz = qmcc*Scxlt
    
        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0d0/msize(1)
        hyi = 1.0d0/msize(2)
        hzi = 1.0d0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)
        zmax = range(6)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        idbeamln = 1

        ntmp1 = 0
        ntmp2 = 0

        do n = 1, innp
          if(rays(5,n)>0) then !//only particles outside the emission plane get accelerated

          ix=(rays(1,n)-xmin)*hxi + 1
          ab=((xmin-rays(1,n))+ix*hx)*hxi
          jx=(rays(3,n)-ymin)*hyi + 1 + jadd
          cd=((ymin-rays(3,n))+(jx-jadd)*hy)*hyi
          kx=(rays(5,n)-zmin)*hzi + 1 + kadd
          ef=((zmin-rays(5,n))+(kx-kadd)*hz)*hzi
        

          ix1 = ix + 1
          if( (ix1.gt.innx).or.(ix.lt.1)) then
            print*,"ix: ",ix,rays(1,n),xmin
            stop
          endif
          jx1 = jx + 1
          if( (jx1.gt.inny).or.(jx.lt.1)) then
            print*,"jx: ",jx,rays(3,n),ymin
            stop
          endif
          kx1 = kx + 1
          if( (kx1.gt.innz).or.(kx.lt.1)) then
            print*,"kx: ",kx,rays(5,n),zmin
            stop
          endif

          exn = (exg(ix,jx,kx)*ab*cd*ef  &
                  +exg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +exg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +exg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +exg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +exg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +exg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +exg(ix1,jx,kx)*(1.0d0-ab)*cd*ef)

          eyn = (eyg(ix,jx,kx)*ab*cd*ef  &
                  +eyg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +eyg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +eyg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +eyg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +eyg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +eyg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +eyg(ix1,jx,kx)*(1.0d0-ab)*cd*ef) 

          ezn = ezg(ix,jx,kx)*ab*cd*ef  &
                  +ezg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +ezg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +ezg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +ezg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +ezg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +ezg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +ezg(ix1,jx,kx)*(1.0d0-ab)*cd*ef

          bxn = bxg(ix,jx,kx)*ab*cd*ef  &
                  +bxg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +bxg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +bxg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +bxg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +bxg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +bxg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +bxg(ix1,jx,kx)*(1.0d0-ab)*cd*ef

          byn = byg(ix,jx,kx)*ab*cd*ef  &
                  +byg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +byg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +byg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +byg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +byg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +byg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +byg(ix1,jx,kx)*(1.0d0-ab)*cd*ef

          !find which element the particle belongs to.
          zz = rays(5,n)*Scxlt
          !counter how many elements overlaped at given location
          noverlap = 0
          do i = ibinit, ibend
            if( (zz.ge.zbeamelem(1,i)) .and. &
                (zz.le.zbeamelem(2,i)) ) then
               noverlap = noverlap + 1 
               idbeamln(noverlap) = i
!               exit
            endif
          enddo

          pos(1) = rays(1,n)*Scxlt
          pos(2) = rays(3,n)*Scxlt
          pos(3) = zz
          pos(4) = tg
          !get external field from all overlaped fields at one location
          extfld = 0.0
          if(flagerr.eq.1) then 
            do nnn = 1, noverlap
              call getflderrt_BeamLineElem(beamelem(idbeamln(nnn)),pos,&
              tmpfld,fldmap(idrfile(3,idbeamln(nnn))))
              extfld = extfld + tmpfld
            enddo
          else
            do nnn = 1, noverlap
              call getfldt_BeamLineElem(beamelem(idbeamln(nnn)),pos,&
              tmpfld,fldmap(idrfile(3,idbeamln(nnn))))
              extfld = extfld + tmpfld
            enddo
          endif
!          !for residence correction
! wait until the discussion with Warp people
!          dev = zbeamelem(2,idbeamln)/Scxlt - rays(5,n)
!          recpgamma = 1.0d0/sqrt(1.0d0+rays(2,n)**2+rays(4,n)**2+rays(6,n)**2)
!          stepsize = rays(6,n)*recpgamma*dt
!          frac = 2*dev/stepsize !half step
!          if(frac.gt.1.0d0) then !no need for correction
!            call getfldt_BeamLineElem(beamelem(idbeamln),pos,extfld,fldmap(idrfile(3,idbeamln)))
!          else
!            frac1 = dev/stepsize+0.5
!            call getfldt_BeamLineElem(beamelem(idbeamln),pos,extfld1,fldmap(idrfile(3,idbeamln)))
!            pos2(1) = pos(1)
!            pos2(2) = pos(2)
!            pos2(3) = pos(3) + 0.5d0*stepsize*Scxlt
!            pos2(4) = pos(4)
!            if(idbeamln.lt.nbeamln) then
!               call getfldt_BeamLineElem(beamelem(idbeamln+1),pos2,&
!               extfld2,fldmap(idrfile(3,idbeamln+1)))
!            else
!               extfld2 = 0.0
!            endif
!            extfld = frac1*extfld1+extfld2*(1.0d0-frac1)
!          endif

          ex = exn+extfld(1)
          ey = eyn+extfld(2)
          ez = ezn+extfld(3)
          bx = bxn+extfld(4)
          by = byn+extfld(5)
          bz = extfld(6)
! test with no B field from space-charge
!          bx = extfld(4)
!          by = extfld(5)
!          bz = extfld(6)

          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefLz*ex*0.5d0*dt
          umy = rays(4,n) + coefLz*ey*0.5d0*dt
          umz = rays(6,n) + coefLz*ez*0.5d0*dt
          recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
          tmp = 0.5d0*Clight*dt*recpgamma*coefLz
          a1 = tmp*bx
          a2 = tmp*by
          a3 = tmp*bz
          a4 = 1.0d0+a1*a1+a2*a2+a3*a3
          s1 = umx + tmp*(umy*bz-umz*by)
          s2 = umy - tmp*(umx*bz-umz*bx)
          s3 = umz + tmp*(umx*by-umy*bx)
          upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
          upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
          rays(2,n) = upx + coefLz*ex*0.5d0*dt
          rays(4,n) = upy + coefLz*ey*0.5d0*dt
          rays(6,n) = upz + coefLz*ez*0.5d0*dt

!          ntmp1 = ntmp1 + 1
!          if(rays(6,n).lt.0.0) then
!            !print*,"rays: ",n,rays(5,n),rays(6,n),upz,ezn,extfld(3)
!            ntmp2 = ntmp2 + 1
!          endif

          endif
        enddo

!        print*,"ntmp1, ",ntmp1,ntmp2,myidx,myidy
        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine kick2t_BeamBunch

        !check the particles outside the computational domain for each bunch/bin
        subroutine lost_BeamBunch(this,xrad,yrad,zleng,zcent,nplc,nptot)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: xrad,yrad,zleng,zcent
        integer, intent(out) :: nplc,nptot
        integer :: i
        integer :: ilost,i0,ierr
        double precision, dimension(6) :: ptrange
        double precision :: rr,rrap
 
        ptrange(1) = -xrad/Scxlt
        ptrange(2) = xrad/Scxlt
        ptrange(3) = -yrad/Scxlt
        ptrange(4) = yrad/Scxlt
        !ptrange(5) = zcent-0.5d0*zleng/Scxlt
        !ptrange(6) = zcent+0.5d0*zleng/Scxlt
        ptrange(5) = 0.0
        ptrange(6) = zleng/Scxlt
        rrap = ptrange(1)**2 
        ilost = 0
        do i0 = 1, this%Nptlocal
          i = i0 - ilost
          this%Pts1(1,i) = this%Pts1(1,i0)
          this%Pts1(2,i) = this%Pts1(2,i0)
          this%Pts1(3,i) = this%Pts1(3,i0)
          this%Pts1(4,i) = this%Pts1(4,i0)
          this%Pts1(5,i) = this%Pts1(5,i0)
          this%Pts1(6,i) = this%Pts1(6,i0)
          rr = this%Pts1(1,i0)**2+this%Pts1(3,i0)**2
          if(rr.ge.rrap) then
            ilost = ilost + 1
          else if(this%Pts1(1,i0).le.ptrange(1)) then
            ilost = ilost + 1
          else if(this%Pts1(1,i0).ge.ptrange(2)) then
            ilost = ilost + 1
          else if(this%Pts1(3,i0).le.ptrange(3)) then
            ilost = ilost + 1
          else if(this%Pts1(3,i0).ge.ptrange(4)) then
            ilost = ilost + 1
          else if(this%Pts1(5,i0).le.ptrange(5) .and. &
                  this%Pts1(6,i0).lt.0.0) then
            ilost = ilost + 1
          else if(this%Pts1(5,i0).ge.ptrange(6)) then
            ilost = ilost + 1
!          else if(this%Pts1(6,i0).lt.0.0) then !this does not allow particles move in negative direction
!            ilost = ilost + 1
          else
          endif
        enddo
        do i = this%Nptlocal - ilost + 1, this%Nptlocal
          this%Pts1(1,i) = 0.0
          this%Pts1(2,i) = 0.0
          this%Pts1(3,i) = 0.0
          this%Pts1(4,i) = 0.0
          this%Pts1(5,i) = -1.0e5
          this%Pts1(6,i) = 0.0
        enddo
        this%Nptlocal = this%Nptlocal - ilost
        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Current = this%Current*nptot/this%Npt
        this%Npt = nptot
 
        end subroutine lost_BeamBunch


        !check the particles outside the computational domain for each bunch/bin
        subroutine lostXY_BeamBunch(this,xradmin,xradmax,yradmin,yradmax,nplc,nptot)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: xradmin,xradmax,yradmin,yradmax
        integer, intent(out) :: nplc,nptot
        integer :: i
        integer :: ilost,i0,ierr
        double precision, dimension(6) :: ptrange
        double precision :: rr,rrap
 
        ptrange(1) = xradmin/Scxlt
        ptrange(2) = xradmax/Scxlt
        ptrange(3) = yradmin/Scxlt
        ptrange(4) = yradmax/Scxlt
        !ptrange(5) = zcent-0.5d0*zleng/Scxlt
        !ptrange(6) = zcent+0.5d0*zleng/Scxlt
        !ptrange(5) = 0.0
        !ptrange(6) = zleng/Scxlt
        rrap = ptrange(1)**2 
        ilost = 0
        do i0 = 1, this%Nptlocal
          i = i0 - ilost
          this%Pts1(1,i) = this%Pts1(1,i0)
          this%Pts1(2,i) = this%Pts1(2,i0)
          this%Pts1(3,i) = this%Pts1(3,i0)
          this%Pts1(4,i) = this%Pts1(4,i0)
          this%Pts1(5,i) = this%Pts1(5,i0)
          this%Pts1(6,i) = this%Pts1(6,i0)
          if(this%Pts1(1,i0).le.ptrange(1)) then
            ilost = ilost + 1
          else if(this%Pts1(1,i0).ge.ptrange(2)) then
            ilost = ilost + 1
          else if(this%Pts1(3,i0).le.ptrange(3)) then
            ilost = ilost + 1
          else if(this%Pts1(3,i0).ge.ptrange(4)) then
            ilost = ilost + 1
          else
          endif
        enddo
        do i = this%Nptlocal - ilost + 1, this%Nptlocal
          this%Pts1(1,i) = 0.0
          this%Pts1(2,i) = 0.0
          this%Pts1(3,i) = 0.0
          this%Pts1(4,i) = 0.0
          this%Pts1(5,i) = -1.0e15
          this%Pts1(6,i) = 0.0
        enddo
        this%Nptlocal = this%Nptlocal - ilost
        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Current = this%Current*nptot/this%Npt
        this%Npt = nptot
 
        end subroutine lostXY_BeamBunch

        !//Calculate the space-charge E and B forces from point-to-point
        !//summation including relativistic effects, and update the particle
        !//momenta using an implicit second order leap-frog scheme.
        !//here, t, dt are in dimensionless unit of (Dt),
        !//totchrg is the total charge in the bunch in unit of Coulomb,
        !//r0 is the spherical ball radius (dimensionless c Dt),
        !//nptlc is the local initial # of particles,
        !//npttot is the total # of particles.
        !//nptrue is the actual total # of particles
        !//the partcile which really participate in the Nbody calculation is controlled by
        !//rays(5,n) > 0 (i.e. z>0).
        subroutine kickpt2pt_BeamBunch(nptlc,rays,&
        tg,chge,mass,dt,beamelem,zbeamelem,idrfile,nbeamln,ibinit,&
        ibend,fldmap,totchrg,r0,npttot,nptrue)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nptlc,nbeamln,ibinit,ibend,npttot,nptrue
        double precision, intent (inout), dimension(6,nptlc) :: rays
        double precision, intent (in) :: dt,mass,chge,tg,totchrg,r0
        type (BeamLineElem), dimension(:), intent(in) :: beamelem
        double precision, dimension(:,:), intent(in) :: zbeamelem
        type (fielddata), dimension(:), intent(in) :: fldmap
        integer, dimension(:,:), intent(in) :: idrfile
        double precision, dimension(6,npttot) :: ptstot
        integer :: n,i
        double precision :: t0
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: range, extfld,tmpfld
        double precision :: qmcc,recpgamma,coefLz,&
                            umx,umy,umz,upx,upy,upz,tmp,a1,a2,a3,a4,s1,s2,&
                            s3,ex,ey,ez,bx,by,bz,zz
        integer, dimension(Maxoverlap) :: idbeamln
        integer :: noverlap,nnn
        double precision :: betajx,betajy,betajz,gammaj,rijx,rijy,&
                            rijz,rij,rijx2,rijy2,rijz2,rij2
        double precision :: tmp1,tmp2,tmp3,coefE
        integer :: j,ierr,sixnpt,jstart,jend,myid,sgn,sixnpttot
        double precision :: dev,frac,frac1,stepsize
        integer :: it
        double precision :: gambetz,betz,gamz

        call starttime_Timer( t0 )

        qmcc = chge/mass
        coefLz = qmcc*Scxlt
        !//coeficients used in the calculation of E field
        !//coeficients of B field is coefE/c
        if(chge.gt.0.0) then
          sgn = 1
        else
          sgn = -1
        endif
        coefE = 1.0d0/(4*Pi*Epsilon0*Scxlt*Scxlt)*sgn*totchrg/nptrue
    
        !collect all particles from the other processors to the local processor.
        sixnpt = nptlc*6
        sixnpttot = nptlc*6
        call MPI_ALLGATHER(rays,sixnpt,MPI_DOUBLE_PRECISION,ptstot,&
             sixnpttot,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

        !get my processor ID.
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        !print*,"coefLz: ",coefLz,qmcc,Scxlt,dt,tg


! for the purpose of comparison with mesh poisson solver
!        gambetz = 0.0
!        do n = 1, npttot 
!          gambetz = gambetz + ptstot(6,n)
!        enddo
!        gambetz = gambetz/npttot
!        gamz = sqrt(1.0d0+gambetz**2)
!        betz = gambetz/gamz

        idbeamln = 1
        do n = 1, nptlc
          if(rays(5,n).gt.0.0) then
 
          !find which element the particle belongs to.
          zz = rays(5,n)*Scxlt
          !counter how many elements overlaped at given location
          noverlap = 0
          do i = ibinit, ibend
            if( (zz.ge.zbeamelem(1,i)) .and. &
                (zz.le.zbeamelem(2,i)) ) then
               noverlap = noverlap + 1
               idbeamln(noverlap) = i
!               exit
            endif
          enddo
 
          pos(1) = rays(1,n)*Scxlt
          pos(2) = rays(3,n)*Scxlt
          pos(3) = zz
          pos(4) = tg
          !get external field from all overlaped fields at one location
          extfld = 0.0
!          print*,"noverlap: ",noverlap
          do nnn = 1, noverlap
            call getfldt_BeamLineElem(beamelem(idbeamln(nnn)),pos,&
            tmpfld,fldmap(idrfile(3,idbeamln(nnn))))
            extfld = extfld + tmpfld
!            if(nnn.eq.2) then
!              write(16,*)zz,tmpfld(4),tmpfld(5),tmpfld(6)
!            endif
          enddo
!           !uniform Ez within 1 mm
!           if( (pos(3).gt.0.0).and.(pos(3).le.1.0e-3)) then
!             extfld(3) = -10.0e6 
!           else
!!             extfld(3) = -10.0e6*exp(-0.5d0*((pos(3)-1.0e-3)/1.0e-3)**2)
!           endif
!           dev = 1.0e-3/Scxlt - rays(5,n)
!           recpgamma = 1.0d0/sqrt(1.0d0+rays(2,n)**2+rays(4,n)**2+rays(6,n)**2)
!           stepsize = rays(6,n)*recpgamma*dt
!           frac = 2*dev/stepsize !half step
!           if(frac.gt.1.0d0) then !no need for correction
!             extfld(3) = -10.0e6 
!           else if(frac.le.0.0) then
!             extfld(3) = 0.0
!           else
!             frac1 = dev/stepsize+0.5
!             extfld(3) = -10.0e6*frac1 
!           endif

!          print*,"pos...",idbeamln,idrfile(3,idbeamln),pos
!          extfld = 0.0
!          ez = 1.0e6
!          print*,pos(3),ex,ey,ez,bx,by,bz,pos(4)

          ex = 0.0 
          ey = 0.0
          ez = 0.0
          bx = 0.0
          by = 0.0
          bz = 0.0
          !//we need some coefficients for Q/(4 pi epsilon0) and 1/C
          !//find the contributions from 1 to i-1 particles of local ptcs.
          it = 0 !count # of particles within the cut-off radius
          do j = 1, n-1
            if(rays(5,j).gt.0.0) then

            rijx = rays(1,n)-rays(1,j)
            rijy = rays(3,n)-rays(3,j)
            rijz = rays(5,n)-rays(5,j)
            gammaj = sqrt(1+rays(2,j)**2+rays(4,j)**2+rays(6,j)**2) 
            betajx = rays(2,j)/gammaj
            betajy = rays(4,j)/gammaj
            betajz = rays(6,j)/gammaj
            ! for the purpose of comparison with mesh poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1) 
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3 
            else
              tmp3 = r0**3
!              it = it+1
            endif
              ex = ex + gammaj*(rijx2-tmp2*betajx)/tmp3 
              ey = ey + gammaj*(rijy2-tmp2*betajy)/tmp3 
              ez = ez + gammaj*(rijz2-tmp2*betajz)/tmp3 
              bx = bx + gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by + gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz + gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
!            else
!              tmp3 = r0**3
!              ex = ex + rijx2/tmp3
!              ey = ey + rijy2/tmp3
!              ex = ez + rijz2/tmp3
!            endif

            endif
          enddo

          !//find the contributions from i+1 to nplocal particles.
          do j = n+1,nptlc
            if(rays(5,j).gt.0.0) then

            rijx = rays(1,n)-rays(1,j)
            rijy = rays(3,n)-rays(3,j)
            rijz = rays(5,n)-rays(5,j)
            gammaj = sqrt(1+rays(2,j)**2+rays(4,j)**2+rays(6,j)**2)
            betajx = rays(2,j)/gammaj
            betajy = rays(4,j)/gammaj
            betajz = rays(6,j)/gammaj
            ! for the purpose of comparison with mesh poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1)
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3
            else
              tmp3 = r0**3
!              it = it+1
            endif
              ex = ex + gammaj*(rijx2-tmp2*betajx)/tmp3
              ey = ey + gammaj*(rijy2-tmp2*betajy)/tmp3
              ez = ez + gammaj*(rijz2-tmp2*betajz)/tmp3
              bx = bx + gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by + gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz + gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
!            else
!              tmp3 = r0**3
!              ex = ex + rijx2/tmp3
!              ey = ey + rijy2/tmp3
!              ex = ez + rijz2/tmp3
!            endif

            endif
          enddo

          !//find the contributions from the other processor particles.
          jstart = 1
          jend = myid*nptlc
          do j = jstart,jend
            if(ptstot(5,j).gt.0.0) then

            rijx = rays(1,n)-ptstot(1,j)
            rijy = rays(3,n)-ptstot(3,j)
            rijz = rays(5,n)-ptstot(5,j)
            gammaj = sqrt(1+ptstot(2,j)**2+ptstot(4,j)**2+ptstot(6,j)**2)
            betajx = ptstot(2,j)/gammaj
            betajy = ptstot(4,j)/gammaj
            betajz = ptstot(6,j)/gammaj
            ! for the purpose of comparison with mesh poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1)
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3
            else
              tmp3 = r0**3
!              it = it+1
            endif
              ex = ex + gammaj*(rijx2-tmp2*betajx)/tmp3
              ey = ey + gammaj*(rijy2-tmp2*betajy)/tmp3
              ez = ez + gammaj*(rijz2-tmp2*betajz)/tmp3
              bx = bx + gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by + gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz + gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
!            else
!              tmp3 = r0**3
!              ex = ex + rijx2/tmp3
!              ey = ey + rijy2/tmp3
!              ex = ez + rijz2/tmp3
!            endif

            endif
          enddo

          jstart =(myid+1)*nptlc+1 
          jend = npttot
          do j = jstart,jend
            if(ptstot(5,j).gt.0.0) then

            rijx = rays(1,n)-ptstot(1,j)
            rijy = rays(3,n)-ptstot(3,j)
            rijz = rays(5,n)-ptstot(5,j)
            gammaj = sqrt(1+ptstot(2,j)**2+ptstot(4,j)**2+ptstot(6,j)**2)
            betajx = ptstot(2,j)/gammaj
            betajy = ptstot(4,j)/gammaj
            betajz = ptstot(6,j)/gammaj
            ! for the purpose of comparison with mesh poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1)
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3
            else
              tmp3 = r0**3
!              it = it+1
            endif
              ex = ex + gammaj*(rijx2-tmp2*betajx)/tmp3
              ey = ey + gammaj*(rijy2-tmp2*betajy)/tmp3
              ez = ez + gammaj*(rijz2-tmp2*betajz)/tmp3
              bx = bx + gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by + gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz + gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
!            else
!              tmp3 = r0**3
!              ex = ex + rijx2/tmp3
!              ey = ey + rijy2/tmp3
!              ex = ez + rijz2/tmp3
!            endif

            endif
          enddo

!          ex = extfld(1)
!          ey = extfld(2)
!          ez = extfld(3)
!          bx = extfld(4)
!          by = extfld(5)
!          bz = extfld(6)
          ex = coefE*ex+extfld(1)
          ey = coefE*ey+extfld(2)
          ez = coefE*ez+extfld(3)
! This is a test of the point-2-point solver, no B field from space-charge
!          bx = extfld(4)
!          by = extfld(5)
!          bz = extfld(6)
          bx = coefE*bx/Clight+extfld(4)
          by = coefE*by/Clight+extfld(5)
          bz = coefE*bz/Clight+extfld(6)

          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefLz*ex*0.5d0*dt
          umy = rays(4,n) + coefLz*ey*0.5d0*dt
          umz = rays(6,n) + coefLz*ez*0.5d0*dt
          recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
          tmp = 0.5d0*Clight*dt*recpgamma*coefLz
          a1 = tmp*bx
          a2 = tmp*by
          a3 = tmp*bz
          a4 = 1.0d0+a1*a1+a2*a2+a3*a3
          s1 = umx + tmp*(umy*bz-umz*by)
          s2 = umy - tmp*(umx*bz-umz*bx)
          s3 = umz + tmp*(umx*by-umy*bx)
          upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
          upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
          rays(2,n) = upx + coefLz*ex*0.5d0*dt
          rays(4,n) = upy + coefLz*ey*0.5d0*dt
          rays(6,n) = upz + coefLz*ez*0.5d0*dt

          endif
        enddo
!        print*,"rays6: ",rays(6,1)
        print*,"# of partcle within radius: ",it,r0
!bx,by,bz,extfld(4),&
!                                              extfld(5),extfld(6)
!        print*,"gambetz: ",gambetz,gamz,betz,npttot

        !//Here, Charge is in unit of electron integer, 
        !//Mass is in the unit of eV. 
        !//coeficients for 1/E0, the coeficients for 1/B0 is just coefE0/c
        !//coeficients for the Lorents force 

        !print*,myid,nptlc,npttot

        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine kickpt2pt_BeamBunch

        !//Include image charge from the cathode
        !//Calculate the space-charge E and B forces from point-to-point
        !//summation including relativistic effects, and update the particle
        !//momenta using an implicit second order leap-frog scheme.
        !//here, t, dt are in dimensionless unit of (Dt),
        !//totchrg is the total charge in the bunch in unit of Coulomb,
        !//r0 is the spherical ball radius (dimensionless c Dt),
        !//nptlc is the initial local # of particles,
        !//npttot is the total # of particles.
        !//nptrue is the actual total # of particles
        !//the partcile which really participate in the Nbody calculation is controlled by
        !//rays(5,n) > 0 (i.e. z>0).
        subroutine kickpt2ptImg_BeamBunch(nptlc,rays,&
        tg,chge,mass,dt,beamelem,zbeamelem,idrfile,nbeamln,ibinit,&
        ibend,fldmap,totchrg,r0,npttot,nptrue)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nptlc,nbeamln,ibinit,ibend,npttot,nptrue
        double precision, intent (inout), dimension(6,nptlc) :: rays
        double precision, intent (in) :: dt,mass,chge,tg,totchrg,r0
        type (BeamLineElem), dimension(:), intent(in) :: beamelem
        double precision, dimension(:,:), intent(in) :: zbeamelem
        type (fielddata), dimension(:), intent(in) :: fldmap
        integer, dimension(:,:), intent(in) :: idrfile
        double precision, dimension(6,npttot) :: ptstot
        integer :: n,i
        double precision :: t0
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: range, extfld,tmpfld
        double precision :: qmcc,recpgamma,coefLz,&
                            umx,umy,umz,upx,upy,upz,tmp,a1,a2,a3,a4,s1,s2,&
                            s3,ex,ey,ez,bx,by,bz,zz
        integer, dimension(Maxoverlap) :: idbeamln
        integer :: noverlap,nnn
        double precision :: betajx,betajy,betajz,gammaj,rijx,rijy,&
                            rijz,rij,rijx2,rijy2,rijz2,rij2
        double precision :: tmp1,tmp2,tmp3,coefE
        integer :: j,ierr,sixnpt,jstart,jend,myid,sgn,sixnpttot
        double precision :: dev,frac,frac1,stepsize
        double precision :: gambetz,betz,gamz

        call starttime_Timer( t0 )

        qmcc = chge/mass
        coefLz = qmcc*Scxlt
        !//coeficients used in the calculation of E field
        !//coeficients of B field is coefE/c
        if(chge.gt.0.0) then
          sgn = 1
        else
          sgn = -1
        endif
        coefE = 1.0d0/(4*Pi*Epsilon0*Scxlt*Scxlt)*sgn*totchrg/nptrue
    
        !collect all particles from the other processors to the local processor.
        sixnpt = nptlc*6
        !sixnpttot = npttot*6
        sixnpttot = nptlc*6
        !print*,"before all gather in kick: ",sixnpt,nptlc,npttot
        call MPI_ALLGATHER(rays,sixnpt,MPI_DOUBLE_PRECISION,ptstot,&
             sixnpttot,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        !print*,"after all gather in kick: ",sixnpt,nptlc,npttot

        !get my processor ID.
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
!        print*,"coefLz: ",coefLz,qmcc,Scxlt,dt,tg

! for the purpose of comparison with mesh Poisson solver
!        gambetz = 0.0
!        do n = 1, npttot
!          gambetz = gambetz + ptstot(6,n)
!        enddo
!        gambetz = gambetz/npttot
!        gamz = sqrt(1.0d0+gambetz**2)
!        betz = gambetz/gamz

        idbeamln = 1
        do n = 1, nptlc
          if(rays(5,n).gt.0.0) then
 
          !find which element the particle belongs to.
          zz = rays(5,n)*Scxlt
          !counter how many elements overlaped at given location
          noverlap = 0
          do i = ibinit, ibend
            if( (zz.ge.zbeamelem(1,i)) .and. &
                (zz.le.zbeamelem(2,i)) ) then
               noverlap = noverlap + 1
               idbeamln(noverlap) = i
!               exit
            endif
          enddo
 
          pos(1) = rays(1,n)*Scxlt
          pos(2) = rays(3,n)*Scxlt
          pos(3) = zz
          pos(4) = tg
          !get external field from all overlaped fields at one location
          extfld = 0.0
!          print*,"noverlap: ",noverlap
          do nnn = 1, noverlap
            call getfldt_BeamLineElem(beamelem(idbeamln(nnn)),pos,&
            tmpfld,fldmap(idrfile(3,idbeamln(nnn))))
            extfld = extfld + tmpfld
!            if(nnn.eq.2) then
!              write(16,*)zz,tmpfld(4),tmpfld(5),tmpfld(6)
!            endif
          enddo
!           !uniform Ez within 1 mm
!           if( (pos(3).gt.0.0).and.(pos(3).le.1.0e-3)) then
!             extfld(3) = -10.0e6 
!           else
!!             extfld(3) = -10.0e6*exp(-0.5d0*((pos(3)-1.0e-3)/1.0e-3)**2)
!           endif
!           dev = 1.0e-3/Scxlt - rays(5,n)
!           recpgamma = 1.0d0/sqrt(1.0d0+rays(2,n)**2+rays(4,n)**2+rays(6,n)**2)
!           stepsize = rays(6,n)*recpgamma*dt
!           frac = 2*dev/stepsize !half step
!           if(frac.gt.1.0d0) then !no need for correction
!             extfld(3) = -10.0e6 
!           else if(frac.le.0.0) then
!             extfld(3) = 0.0
!           else
!             frac1 = dev/stepsize+0.5
!             extfld(3) = -10.0e6*frac1 
!           endif

!          print*,"pos...",idbeamln,idrfile(3,idbeamln),pos
!          extfld = 0.0
!          ez = 1.0e6
!          print*,pos(3),ex,ey,ez,bx,by,bz,pos(4)

          ex = 0.0 
          ey = 0.0
          ez = 0.0
          bx = 0.0
          by = 0.0
          bz = 0.0
          !//we need some coefficients for Q/(4 pi epsilon0) and 1/C
          !//find the contributions from 1 to j-1 particles of local ptcs.
          !//space-charge of the beam itself
          do j = 1, n-1
            if(rays(5,j).gt.0.0) then

            rijx = rays(1,n)-rays(1,j)
            rijy = rays(3,n)-rays(3,j)
            rijz = rays(5,n)-rays(5,j)
            gammaj = sqrt(1+rays(2,j)**2+rays(4,j)**2+rays(6,j)**2) 
            betajx = rays(2,j)/gammaj
            betajy = rays(4,j)/gammaj
            betajz = rays(6,j)/gammaj
            ! for the purpose of comparison with mesh Poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1) 
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3 
            else
              tmp3 = r0**3
            endif
              ex = ex + gammaj*(rijx2-tmp2*betajx)/tmp3 
              ey = ey + gammaj*(rijy2-tmp2*betajy)/tmp3 
              ez = ez + gammaj*(rijz2-tmp2*betajz)/tmp3 
              bx = bx + gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by + gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz + gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
            !  ex = ex + rijx2/tmp3
            !  ey = ey + rijy2/tmp3
            !  ex = ez + rijz2/tmp3

            endif
          enddo

          !//find the contributions from i+1 to nplocal particles.
          do j = n+1,nptlc
            if(rays(5,j).gt.0.0) then

            rijx = rays(1,n)-rays(1,j)
            rijy = rays(3,n)-rays(3,j)
            rijz = rays(5,n)-rays(5,j)
            gammaj = sqrt(1+rays(2,j)**2+rays(4,j)**2+rays(6,j)**2)
            betajx = rays(2,j)/gammaj
            betajy = rays(4,j)/gammaj
            betajz = rays(6,j)/gammaj
            ! for the purpose of comparison with mesh Poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1)
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3
            else
              tmp3 = r0**3
            endif
              ex = ex + gammaj*(rijx2-tmp2*betajx)/tmp3
              ey = ey + gammaj*(rijy2-tmp2*betajy)/tmp3
              ez = ez + gammaj*(rijz2-tmp2*betajz)/tmp3
              bx = bx + gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by + gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz + gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
            !  ex = ex + rijx2/tmp3
            !  ey = ey + rijy2/tmp3
            !  ex = ez + rijz2/tmp3

            endif
          enddo

          !//Image space-charge of the beam 
          !//find the contributions from 1 to j-1 particles of local ptcs.
          do j = 1, nptlc
            if(rays(5,j).gt.0.0) then
            if(j.ne.n) then

            rijx = rays(1,n)-rays(1,j)
            rijy = rays(3,n)-rays(3,j)
            !due to image charge
            rijz = rays(5,n)+rays(5,j)
            gammaj = sqrt(1+rays(2,j)**2+rays(4,j)**2+rays(6,j)**2) 
            betajx = rays(2,j)/gammaj
            betajy = rays(4,j)/gammaj
            !due to image charge
            betajz = -rays(6,j)/gammaj
            ! for the purpose of comparison with mesh Poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = -betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1) 
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3 
            else
              tmp3 = r0**3
            endif
              !due to image charge
              ex = ex - gammaj*(rijx2-tmp2*betajx)/tmp3 
              ey = ey - gammaj*(rijy2-tmp2*betajy)/tmp3 
              ez = ez - gammaj*(rijz2-tmp2*betajz)/tmp3 
              bx = bx - gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by - gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz - gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
              !due to image charge
            !  ex = ex - rijx2/tmp3
            !  ey = ey - rijy2/tmp3
            !  ex = ez - rijz2/tmp3

            endif
            endif
          enddo

          !//space charge from beam itself
          !//find the contributions from the other processor particles.
          jstart = 1
          jend = myid*nptlc
          do j = jstart,jend
            if(ptstot(5,j).gt.0.0) then

            rijx = rays(1,n)-ptstot(1,j)
            rijy = rays(3,n)-ptstot(3,j)
            rijz = rays(5,n)-ptstot(5,j)
            gammaj = sqrt(1+ptstot(2,j)**2+ptstot(4,j)**2+ptstot(6,j)**2)
            betajx = ptstot(2,j)/gammaj
            betajy = ptstot(4,j)/gammaj
            betajz = ptstot(6,j)/gammaj
            ! for the purpose of comparison with mesh Poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1)
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3
            else
              tmp3 = r0**3
            endif
              ex = ex + gammaj*(rijx2-tmp2*betajx)/tmp3
              ey = ey + gammaj*(rijy2-tmp2*betajy)/tmp3
              ez = ez + gammaj*(rijz2-tmp2*betajz)/tmp3
              bx = bx + gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by + gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz + gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
            !  ex = ex + rijx2/tmp3
            !  ey = ey + rijy2/tmp3
            !  ex = ez + rijz2/tmp3

            endif
          enddo

          !//image space charge from beam
          !//find the contributions from the other processor particles.
          jstart = 1
          jend = myid*nptlc
          do j = jstart,jend
            if(ptstot(5,j).gt.0.0) then

            rijx = rays(1,n)-ptstot(1,j)
            rijy = rays(3,n)-ptstot(3,j)
            rijz = rays(5,n)+ptstot(5,j) !due to image charge
            gammaj = sqrt(1+ptstot(2,j)**2+ptstot(4,j)**2+ptstot(6,j)**2)
            betajx = ptstot(2,j)/gammaj
            betajy = ptstot(4,j)/gammaj
            betajz = -ptstot(6,j)/gammaj !image charge
            ! for the purpose of comparison with mesh Poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = -betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1)
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3
            else
              tmp3 = r0**3
            endif
              !image charge
              ex = ex - gammaj*(rijx2-tmp2*betajx)/tmp3
              ey = ey - gammaj*(rijy2-tmp2*betajy)/tmp3
              ez = ez - gammaj*(rijz2-tmp2*betajz)/tmp3
              bx = bx - gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by - gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz - gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
              !image charge
            !  ex = ex - rijx2/tmp3
            !  ey = ey - rijy2/tmp3
            !  ex = ez - rijz2/tmp3

            endif
          enddo

          !//space charge from beam itself
          jstart =(myid+1)*nptlc+1 
          jend = npttot
          do j = jstart,jend
            if(ptstot(5,j).gt.0.0) then

            rijx = rays(1,n)-ptstot(1,j)
            rijy = rays(3,n)-ptstot(3,j)
            rijz = rays(5,n)-ptstot(5,j)
            gammaj = sqrt(1+ptstot(2,j)**2+ptstot(4,j)**2+ptstot(6,j)**2)
            betajx = ptstot(2,j)/gammaj
            betajy = ptstot(4,j)/gammaj
            betajz = ptstot(6,j)/gammaj
            ! for the purpose of comparison with mesh Poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1)
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3
            else
              tmp3 = r0**3
            endif
              ex = ex + gammaj*(rijx2-tmp2*betajx)/tmp3
              ey = ey + gammaj*(rijy2-tmp2*betajy)/tmp3
              ez = ez + gammaj*(rijz2-tmp2*betajz)/tmp3
              bx = bx + gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by + gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz + gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
            !  ex = ex + rijx2/tmp3
            !  ey = ey + rijy2/tmp3
            !  ex = ez + rijz2/tmp3

            endif
          enddo

          !//image space charge from beam
          jstart =(myid+1)*nptlc+1 
          jend = npttot
          do j = jstart,jend
            if(ptstot(5,j).gt.0.0) then

            rijx = rays(1,n)-ptstot(1,j)
            rijy = rays(3,n)-ptstot(3,j)
            rijz = rays(5,n)+ptstot(5,j) !image charge
            gammaj = sqrt(1+ptstot(2,j)**2+ptstot(4,j)**2+ptstot(6,j)**2)
            betajx = ptstot(2,j)/gammaj
            betajy = ptstot(4,j)/gammaj
            betajz = -ptstot(6,j)/gammaj !image charge
            ! for the purpose of comparison with mesh Poisson solver
            !gammaj = gamz
            !betajx = 0.0
            !betajy = 0.0
            !betajz = -betz
            tmp1 = gammaj**2*(rijx*betajx+rijy*betajy+rijz*betajz)/(gammaj+1)
            rijx2 = rijx + tmp1*betajx
            rijy2 = rijy + tmp1*betajy
            rijz2 = rijz + tmp1*betajz
            rij2 = sqrt(rijx2**2+rijy2**2+rijz2**2)
              tmp2 = gammaj*(betajx*rijx2+betajy*rijy2+betajz*rijz2)/(gammaj+1)
            if(rij2.gt.r0) then
              tmp3 = rij2**3
            else
              tmp3 = r0**3
            endif
              ex = ex - gammaj*(rijx2-tmp2*betajx)/tmp3
              ey = ey - gammaj*(rijy2-tmp2*betajy)/tmp3
              ez = ez - gammaj*(rijz2-tmp2*betajz)/tmp3
              bx = bx - gammaj*(betajy*rijz2-betajz*rijy2)/tmp3
              by = by - gammaj*(betajz*rijx2-betajx*rijz2)/tmp3
              bz = bz - gammaj*(betajx*rijy2-betajy*rijx2)/tmp3
            !  ex = ex - rijx2/tmp3
            !  ey = ey - rijy2/tmp3
            !  ex = ez - rijz2/tmp3

            endif
          enddo

!          ex = extfld(1)
!          ey = extfld(2)
!          ez = extfld(3)
!          bx = extfld(4)
!          by = extfld(5)
!          bz = extfld(6)
          ex = coefE*ex+extfld(1)
          ey = coefE*ey+extfld(2)
          ez = coefE*ez+extfld(3)
          bx = coefE*bx/Clight+extfld(4)
          by = coefE*by/Clight+extfld(5)
          bz = coefE*bz/Clight+extfld(6)
!          bx = extfld(4)
!          by = extfld(5)
!          bz = extfld(6)

          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefLz*ex*0.5d0*dt
          umy = rays(4,n) + coefLz*ey*0.5d0*dt
          umz = rays(6,n) + coefLz*ez*0.5d0*dt
          recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
          tmp = 0.5d0*Clight*dt*recpgamma*coefLz
          a1 = tmp*bx
          a2 = tmp*by
          a3 = tmp*bz
          a4 = 1.0d0+a1*a1+a2*a2+a3*a3
          s1 = umx + tmp*(umy*bz-umz*by)
          s2 = umy - tmp*(umx*bz-umz*bx)
          s3 = umz + tmp*(umx*by-umy*bx)
          upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
          upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
          rays(2,n) = upx + coefLz*ex*0.5d0*dt
          rays(4,n) = upy + coefLz*ey*0.5d0*dt
          rays(6,n) = upz + coefLz*ez*0.5d0*dt

          endif
        enddo
!        print*,"rays6: ",rays(6,1)

       
        !//Here, Charge is in unit of electron integer, 
        !//Mass is in the unit of eV. 
        !//coeficients for 1/E0, the coeficients for 1/B0 is just coefE0/c
        !//coeficients for the Lorents force 

        !print*,myid,nptlc,npttot

        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine kickpt2ptImg_BeamBunch

        !rotate to the particle coordinates to local beam coordinate of "ptref".
        !Here, "ptref" coordinate has a rotation "theta" in x-z plane with
        !respect to the orginal coordinate.
        subroutine rottoT_BeamBunch(this,ptref,ptrange,poscent)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, dimension(6) :: ptref
        double precision, dimension(6), intent(out) :: ptrange
        double precision  :: poslc,poscent
        double precision  :: cs,ss,gamma
        double precision, dimension(6) :: temp
        integer :: i,ierr
 
        cs = ptref(6)/sqrt(ptref(6)**2 + ptref(2)**2)
        ss = ptref(2)/sqrt(ptref(6)**2 + ptref(2)**2)
        do i = 1, 3
          ptrange(2*i-1) = 1.0e20
          ptrange(2*i) = -1.0e20
        enddo

        do i = 1, this%Nptlocal
          temp(1) = this%Pts1(1,i)
          temp(2) = this%Pts1(2,i)
          temp(3) = this%Pts1(3,i) 
          temp(4) = this%Pts1(4,i)
          temp(5) = this%Pts1(5,i) 
          temp(6) = this%Pts1(6,i)
          this%Pts1(1,i) = temp(1)*cs - temp(5)*ss
          this%Pts1(2,i) = temp(2)*cs - temp(6)*ss
          this%Pts1(3,i) = temp(3)
          this%Pts1(4,i) = temp(4)
          this%Pts1(5,i) = temp(1)*ss + temp(5)*cs
          this%Pts1(6,i) = temp(2)*ss + temp(6)*cs
        enddo
        do i = 1, this%Nptlocal
          if(ptrange(1).gt.this%Pts1(1,i)) then
            ptrange(1) = this%Pts1(1,i) 
          endif
          if(ptrange(2).lt.this%Pts1(1,i)) then
            ptrange(2) = this%Pts1(1,i) 
          endif
          if(ptrange(3).gt.this%Pts1(3,i)) then
            ptrange(3) = this%Pts1(3,i) 
          endif
          if(ptrange(4).lt.this%Pts1(3,i)) then
            ptrange(4) = this%Pts1(3,i) 
          endif
          if(ptrange(5).gt.this%Pts1(5,i)) then
            ptrange(5) = this%Pts1(5,i) 
          endif
          if(ptrange(6).lt.this%Pts1(5,i)) then
            ptrange(6) = this%Pts1(5,i) 
          endif
        enddo

        !find the centroid z location
        poslc = 0.0d0
        do i = 1, this%Nptlocal
          poslc = poslc + this%Pts1(5,i)
        enddo
        call MPI_ALLREDUCE(poslc,poscent,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        poscent = poscent/this%Npt

        end subroutine rottoT_BeamBunch

        subroutine rotbackT_BeamBunch(this,ptref)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, dimension(6) :: ptref
        double precision, dimension(6) :: tmp
        double precision  :: cs,ss
        integer :: i
 
        cs = ptref(6)/sqrt(ptref(6)**2 + ptref(2)**2)
        ss = ptref(2)/sqrt(ptref(6)**2 + ptref(2)**2)
        do i = 1, this%Nptlocal
          tmp(1) = this%Pts1(1,i)*cs + this%Pts1(5,i)*ss
          tmp(2) = this%Pts1(2,i)*cs + this%Pts1(6,i)*ss
          tmp(3) = this%Pts1(3,i)
          tmp(4) = this%Pts1(4,i)
          tmp(5) = -this%Pts1(1,i)*ss + this%Pts1(5,i)*cs
          tmp(6) = -this%Pts1(2,i)*ss + this%Pts1(6,i)*cs
          this%Pts1(1,i) = tmp(1)
          this%Pts1(2,i) = tmp(2)
          this%Pts1(3,i) = tmp(3) 
          this%Pts1(4,i) = tmp(4) 
          this%Pts1(5,i) = tmp(5)
          this%Pts1(6,i) = tmp(6)
        enddo
 
        end subroutine rotbackT_BeamBunch

        !convert to the local Cartesian coordinate at the entry of bend
        subroutine convEntr_BeamBunch(this,zorgin,gamin)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: zorgin,gamin
        integer :: i
        double precision :: gam
 
        gam = -this%refptcl(6)
 
        !go to the local Cartesian coordinates at the beginning of bend
        do i = 1, this%Nptlocal
          this%Pts1(5,i) = this%Pts1(5,i) - zorgin/Scxlt 
        enddo

        this%refptcl(1) = 0.0
        this%refptcl(2) = 0.0
        this%refptcl(3) = 0.0
        this%refptcl(4) = 0.0
        this%refptcl(5) = this%refptcl(5) - zorgin/Scxlt
        this%refptcl(6) = sqrt(gamin**2  - 1.0d0) !//gamma beta_z
          
        end subroutine convEntr_BeamBunch

        !The Cartesian coordinate has been rotated following the exit direction 
        !of the reference particle. The exit angle of the reference particle
        !should correspond to the bend angle.
        subroutine convExit_BeamBunch(this,zorgin)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        real*8 :: zorgin
        integer :: i
        double precision :: gamma0
        double precision :: cs,ss
        double precision, dimension(6) :: temp
 
        !rotate and shift to the local coordinates of the reference particle.
        !However, there is no shift of the momentum
        cs = this%refptcl(6)/sqrt(this%refptcl(6)**2 + this%refptcl(2)**2)
        ss = this%refptcl(2)/sqrt(this%refptcl(6)**2 + this%refptcl(2)**2)
        do i = 1, this%Nptlocal
          temp(1) = this%Pts1(1,i) - this%refptcl(1)
          temp(2) = this%Pts1(2,i) 
          temp(3) = this%Pts1(3,i) - this%refptcl(3)
          temp(4) = this%Pts1(4,i)
          temp(5) = this%Pts1(5,i) - this%refptcl(5)
          temp(6) = this%Pts1(6,i)
          this%Pts1(1,i) = temp(1)*cs - temp(5)*ss
          this%Pts1(2,i) = temp(2)*cs - temp(6)*ss
          this%Pts1(3,i) = temp(3)
          this%Pts1(4,i) = temp(4)
          this%Pts1(5,i) = temp(1)*ss + temp(5)*cs + zorgin/Scxlt
          this%Pts1(6,i) = temp(2)*ss + temp(6)*cs
        enddo
        temp(2) = this%refptcl(2)
        temp(6) = this%refptcl(6)
        this%refptcl(1) = 0.0 
        this%refptcl(2) = 0.0
        this%refptcl(3) = 0.0 
        this%refptcl(4) = 0.0 
        this%refptcl(5) = 0.0 + zorgin/Scxlt
        this%refptcl(6) = sqrt(temp(2)**2 + temp(6)**2)

        gamma0 = sqrt(1.0d0+this%refptcl(6)**2)
        this%refptcl(6) = -gamma0
 
        end subroutine convExit_BeamBunch

        !The Cartesian coordinate has been rotated following the exit direction 
        !of the reference particle. The exit angle of the reference particle
        !should correspond to the bend angle.
        subroutine convExitold_BeamBunch(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer :: i
        double precision :: gamma0
        double precision :: cs,ss
        double precision, dimension(6) :: temp
 
        !rotate and shift to the local coordinates of the reference particle.
        !However, there is no shift of the momentum
        cs = this%refptcl(6)/sqrt(this%refptcl(6)**2 + this%refptcl(2)**2)
        ss = this%refptcl(2)/sqrt(this%refptcl(6)**2 + this%refptcl(2)**2)
        do i = 1, this%Nptlocal
          temp(1) = this%Pts1(1,i) - this%refptcl(1)
          temp(2) = this%Pts1(2,i) 
          temp(3) = this%Pts1(3,i) - this%refptcl(3)
          temp(4) = this%Pts1(4,i)
          temp(5) = this%Pts1(5,i) - this%refptcl(5)
          temp(6) = this%Pts1(6,i)
          this%Pts1(1,i) = temp(1)*cs - temp(5)*ss
          this%Pts1(2,i) = temp(2)*cs - temp(6)*ss
          this%Pts1(3,i) = temp(3)
          this%Pts1(4,i) = temp(4)
          this%Pts1(5,i) = temp(1)*ss + temp(5)*cs
          this%Pts1(6,i) = temp(2)*ss + temp(6)*cs
        enddo
        temp(2) = this%refptcl(2)
        temp(6) = this%refptcl(6)
        this%refptcl(1) = 0.0 
        this%refptcl(2) = 0.0
        this%refptcl(3) = 0.0 
        this%refptcl(4) = 0.0 
        this%refptcl(5) = 0.0
        this%refptcl(6) = sqrt(temp(2)**2 + temp(6)**2)

        gamma0 = sqrt(1.0d0+this%refptcl(6)**2)
        this%refptcl(6) = -gamma0
 
        end subroutine convExitold_BeamBunch

        !//drift half step in positions.
        !//Here, x, y, z are normalized by C * Dt
        !//tau - normalized step size (by Dt).
        !//the time "t" is normalized by the scaling frequency.
        subroutine driftbackhalf_BeamBunch(this,t,tau)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: t
        double precision, intent (in) :: tau
        double precision :: xl,pz
        double precision :: t0,recpgam
        integer :: i
 
        call starttime_Timer(t0)
 
        do i = 1, this%Nptlocal
          !//get 1.0d0/gamma of each particle
          recpgam = 1.0d0/sqrt(1.0d0+this%Pts1(2,i)**2+this%Pts1(4,i)**2+&
                                 this%Pts1(6,i)**2)
          this%Pts1(1,i) = this%Pts1(1,i)-0.5d0*tau*this%Pts1(2,i)*recpgam
          this%Pts1(3,i) = this%Pts1(3,i)-0.5d0*tau*this%Pts1(4,i)*recpgam
          this%Pts1(5,i) = this%Pts1(5,i)-0.5d0*tau*this%Pts1(6,i)*recpgam
        enddo
        recpgam = 1.0d0/sqrt(1.0d0+this%refptcl(2)**2+this%refptcl(4)**2+&
                           this%refptcl(6)**2)
        this%refptcl(1) = this%refptcl(1)-0.5d0*tau*this%refptcl(2)*recpgam
        this%refptcl(3) = this%refptcl(3)-0.5d0*tau*this%refptcl(4)*recpgam
        this%refptcl(5) = this%refptcl(5)-0.5d0*tau*this%refptcl(6)*recpgam
 
        t_map1 = t_map1 + elapsedtime_Timer(t0)
 
        end subroutine driftbackhalf_BeamBunch

        !//drift half step in positions.
        !//Here, x, y, z are normalized by C * Dt
        !//tau - normalized step size (by Dt).
        !//the time "t" is normalized by the scaling frequency.
        subroutine drifthalfBd_BeamBunch(this,t,tau)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: t
        double precision, intent (in) :: tau
        double precision :: xl,pz
        double precision :: t0,recpgam
        integer :: i
 
        call starttime_Timer(t0)
 
        do i = 1, this%Nptlocal
          !//get 1.0d0/gamma of each particle
          recpgam = 1.0d0/sqrt(1.0d0+this%Pts1(2,i)**2+this%Pts1(4,i)**2+&
                                 this%Pts1(6,i)**2)
          this%Pts1(1,i) = this%Pts1(1,i)+0.5d0*tau*this%Pts1(2,i)*recpgam
          this%Pts1(3,i) = this%Pts1(3,i)+0.5d0*tau*this%Pts1(4,i)*recpgam
          this%Pts1(5,i) = this%Pts1(5,i)+0.5d0*tau*this%Pts1(6,i)*recpgam
        enddo
        recpgam = 1.0d0/sqrt(1.0d0+this%refptcl(2)**2+this%refptcl(4)**2+&
                           this%refptcl(6)**2)
        this%refptcl(1) = this%refptcl(1)+0.5d0*tau*this%refptcl(2)*recpgam
        this%refptcl(3) = this%refptcl(3)+0.5d0*tau*this%refptcl(4)*recpgam
        this%refptcl(5) = this%refptcl(5)+0.5d0*tau*this%refptcl(6)*recpgam

        t_map1 = t_map1 + elapsedtime_Timer(t0)
 
        end subroutine drifthalfBd_BeamBunch

        !interpolate the space-charge fields E and B in the lab frame
        !to individual particle + external fields.
        subroutine kick2tBd_BeamBunch(innp,innx,inny,innz,rays,exg,&
        eyg,ezg,bxg,byg,bzg,ptsgeom,npx,npy,myidx,myidy,tg,&
        chge,mass,dt,beamelem,zbeamelem,idrfile,nbeamln,idbd,&
        fldmap,refpt)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,&
                               myidy,nbeamln,idbd
        double precision, intent (inout), dimension(6,innp) :: rays
        double precision, intent (in), dimension(innx,inny,innz) :: exg,eyg,ezg 
        double precision, intent (in), dimension(innx,inny,innz) :: bxg,byg,bzg 
        double precision, intent (in) :: dt,mass,chge,tg
        type (BeamLineElem), dimension(:), intent(in) :: beamelem
        double precision, dimension(:,:), intent(in) :: zbeamelem
        integer, dimension(:,:), intent(in) :: idrfile
        type (fielddata), dimension(:), intent(in) :: fldmap
        double precision, dimension(6) :: refpt
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab, cd, ef
        integer :: n,i,ii
        double precision :: t0
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin,zmax
        double precision, dimension(3) :: msize
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: range, extfld
        type (CompDom) :: ptsgeom
        integer :: ix,jx,kx,ix1,jx1,kx1,ierr,kadd,jadd
        double precision :: qmcc,recpgamma,coefLz,&
                            umx,umy,umz,upx,upy,upz,tmp,a1,a2,a3,a4,s1,s2,&
                            s3,exn,eyn,ezn,ex,ey,ez,bx,by,bz,zz,bxn,byn
        !for residence correction purpose
        double precision :: dev,frac,frac1,stepsize
        double precision, dimension(6) :: tmpfld
        integer :: ntmp1,ntmp2
        double precision :: cs,ss

        call starttime_Timer( t0 )

        qmcc = chge/mass
        coefLz = qmcc*Scxlt
        !cos and sin of the rotation angle. In a straight machine,
        !the rotation angle is 0.
        cs = refpt(6)/sqrt(refpt(2)**2+refpt(6)**2)
        ss = refpt(2)/sqrt(refpt(2)**2+refpt(6)**2)

    
        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0d0/msize(1)
        hyi = 1.0d0/msize(2)
        hzi = 1.0d0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)
        zmax = range(6)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        ntmp1 = 0
        ntmp2 = 0

        pos(1) = refpt(1)*Scxlt
        pos(2) = refpt(3)*Scxlt
        pos(3) = refpt(5)*Scxlt
        pos(4) = tg
        !get external field from all overlaped fields at one location
        extfld = 0.0
        call getfldt_BeamLineElem(beamelem(idbd),pos,&
        extfld,fldmap(idrfile(3,idbd)))

        ex = extfld(1)
        ey = extfld(2)
        ez = extfld(3)
        bx = extfld(4)
        by = extfld(5)
        bz = extfld(6)

        !//advance the momenta of particles using implicit central
        !//difference scheme from Birdall and Longdon's book.
        umx = refpt(2) + coefLz*ex*0.5d0*dt
        umy = refpt(4) + coefLz*ey*0.5d0*dt
        umz = refpt(6) + coefLz*ez*0.5d0*dt
        recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
        tmp = 0.5d0*Clight*dt*recpgamma*coefLz
        a1 = tmp*bx
        a2 = tmp*by
        a3 = tmp*bz
        a4 = 1.0d0+a1*a1+a2*a2+a3*a3
        s1 = umx + tmp*(umy*bz-umz*by)
        s2 = umy - tmp*(umx*bz-umz*bx)
        s3 = umz + tmp*(umx*by-umy*bx)
        upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
        upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
        upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
        refpt(2) = upx + coefLz*ex*0.5d0*dt
        refpt(4) = upy + coefLz*ey*0.5d0*dt
        refpt(6) = upz + coefLz*ez*0.5d0*dt
 
        do n = 1, innp
          ix=(rays(1,n)-xmin)*hxi + 1
          ab=((xmin-rays(1,n))+ix*hx)*hxi
          jx=(rays(3,n)-ymin)*hyi + 1 + jadd
          cd=((ymin-rays(3,n))+(jx-jadd)*hy)*hyi
          kx=(rays(5,n)-zmin)*hzi + 1 + kadd
          ef=((zmin-rays(5,n))+(kx-kadd)*hz)*hzi
        

          ix1 = ix + 1
          if( (ix1.gt.innx).or.(ix.lt.1)) then
            print*,"ix: ",ix,rays(1,n),xmin
            stop
          endif
          jx1 = jx + 1
          if( (jx1.gt.inny).or.(jx.lt.1)) then
            print*,"jx: ",jx,rays(3,n),ymin
            stop
          endif
          kx1 = kx + 1
          if( (kx1.gt.innz).or.(kx.lt.1)) then
            print*,"kx: ",kx,rays(5,n),zmin
            stop
          endif

          exn = (exg(ix,jx,kx)*ab*cd*ef  &
                  +exg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +exg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +exg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +exg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +exg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +exg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +exg(ix1,jx,kx)*(1.0d0-ab)*cd*ef)

          eyn = (eyg(ix,jx,kx)*ab*cd*ef  &
                  +eyg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +eyg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +eyg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +eyg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +eyg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +eyg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +eyg(ix1,jx,kx)*(1.0d0-ab)*cd*ef) 

          ezn = ezg(ix,jx,kx)*ab*cd*ef  &
                  +ezg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +ezg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +ezg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +ezg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +ezg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +ezg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +ezg(ix1,jx,kx)*(1.0d0-ab)*cd*ef

          bxn = bxg(ix,jx,kx)*ab*cd*ef  &
                  +bxg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +bxg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +bxg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +bxg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +bxg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +bxg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +bxg(ix1,jx,kx)*(1.0d0-ab)*cd*ef

          byn = byg(ix,jx,kx)*ab*cd*ef  &
                  +byg(ix,jx1,kx)*ab*(1.0d0-cd)*ef &
                  +byg(ix,jx1,kx1)*ab*(1.0d0-cd)*(1.0d0-ef) &
                  +byg(ix,jx,kx1)*ab*cd*(1.0d0-ef) &
                  +byg(ix1,jx,kx1)*(1.0d0-ab)*cd*(1.0d0-ef) &
                  +byg(ix1,jx1,kx1)*(1.0d0-ab)*(1.0d0-cd)*(1.0d0-ef)&
                  +byg(ix1,jx1,kx)*(1.0d0-ab)*(1.0d0-cd)*ef &
                  +byg(ix1,jx,kx)*(1.0d0-ab)*cd*ef

          !get field in Cartesian coordinate from analytical function.
          !we need to rotate the coordinate back to the orginal coordinate.
          pos(1) = (cs*rays(1,n)+ss*rays(5,n))*Scxlt
          pos(2) = rays(3,n)*Scxlt
          pos(3) = (-ss*rays(1,n) + cs*rays(5,n))*Scxlt
          pos(4) = tg

          !get external field from all overlaped fields at one location
          call getfldt_BeamLineElem(beamelem(idbd),pos,&
          extfld,fldmap(idrfile(3,idbd)))

          ex = exn+extfld(1)*cs - extfld(3)*ss
          ey = eyn+extfld(2)
          ez = ezn+extfld(1)*ss + extfld(3)*cs
          bx = bxn+extfld(4)*cs - extfld(6)*ss
          by = byn+extfld(5)
          bz = extfld(4)*ss + extfld(6)*cs

          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefLz*ex*0.5d0*dt
          umy = rays(4,n) + coefLz*ey*0.5d0*dt
          umz = rays(6,n) + coefLz*ez*0.5d0*dt
          recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
          tmp = 0.5d0*Clight*dt*recpgamma*coefLz
          a1 = tmp*bx
          a2 = tmp*by
          a3 = tmp*bz
          a4 = 1.0d0+a1*a1+a2*a2+a3*a3
          s1 = umx + tmp*(umy*bz-umz*by)
          s2 = umy - tmp*(umx*bz-umz*bx)
          s3 = umz + tmp*(umx*by-umy*bx)
          upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
          upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
          rays(2,n) = upx + coefLz*ex*0.5d0*dt
          rays(4,n) = upy + coefLz*ey*0.5d0*dt
          rays(6,n) = upz + coefLz*ez*0.5d0*dt

        enddo

        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine kick2tBd_BeamBunch

        !interpolate the E and B in the lab frame
        !to individual particle, external fields.
        subroutine kick2tBd0_BeamBunch(innp,rays,tg,&
        chge,mass,dt,beamelem,zbeamelem,idrfile,nbeamln,idbd,&
        fldmap,refpt)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,nbeamln,idbd
        double precision, intent (inout), dimension(6,innp) :: rays
        double precision, intent (in) :: dt,mass,chge,tg
        type (BeamLineElem), dimension(:), intent(in) :: beamelem
        double precision, dimension(:,:), intent(in) :: zbeamelem
        integer, dimension(:,:), intent(in) :: idrfile
        type (fielddata), dimension(:), intent(in) :: fldmap
        double precision, dimension(6) :: refpt
        integer :: n,i
        double precision :: t0
        double precision, dimension(4) :: pos
        double precision, dimension(6) :: extfld
        integer :: ierr
        double precision :: qmcc,recpgamma,coefLz,&
                            umx,umy,umz,upx,upy,upz,tmp,a1,a2,a3,a4,s1,s2,&
                            s3,ex,ey,ez,bx,by,bz,zz
        !for residence correction purpose
        double precision :: dev,frac,frac1,stepsize
        double precision :: cs,ss

        call starttime_Timer( t0 )

        qmcc = chge/mass
        coefLz = qmcc*Scxlt
        !cos and sin of the rotation angle. In a straight machine,
        !the rotation angle is 0.
        cs = refpt(6)/sqrt(refpt(2)**2+refpt(6)**2)
        ss = refpt(2)/sqrt(refpt(2)**2+refpt(6)**2)

        pos(1) = refpt(1)*Scxlt
        pos(2) = refpt(3)*Scxlt
        pos(3) = refpt(5)*Scxlt
        pos(4) = tg
        !get external field from all overlaped fields at one location
        extfld = 0.0
        call getfldt_BeamLineElem(beamelem(idbd),pos,&
        extfld,fldmap(idrfile(3,idbd)))

        ex = extfld(1)
        ey = extfld(2)
        ez = extfld(3)
        bx = extfld(4)
        by = extfld(5)
        bz = extfld(6)

        !//advance the momenta of particles using implicit central
        !//difference scheme from Birdall and Longdon's book.
        umx = refpt(2) + coefLz*ex*0.5d0*dt
        umy = refpt(4) + coefLz*ey*0.5d0*dt
        umz = refpt(6) + coefLz*ez*0.5d0*dt
        recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
        tmp = 0.5d0*Clight*dt*recpgamma*coefLz
        a1 = tmp*bx
        a2 = tmp*by
        a3 = tmp*bz
        a4 = 1.0d0+a1*a1+a2*a2+a3*a3
        s1 = umx + tmp*(umy*bz-umz*by)
        s2 = umy - tmp*(umx*bz-umz*bx)
        s3 = umz + tmp*(umx*by-umy*bx)
        upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
        upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
        upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
        refpt(2) = upx + coefLz*ex*0.5d0*dt
        refpt(4) = upy + coefLz*ey*0.5d0*dt
        refpt(6) = upz + coefLz*ez*0.5d0*dt
 
        do n = 1, innp
          !get field in Cartesian coordinate from analytical function.
          !we need to rotate the coordinate back to the orginal coordinate.
          pos(1) = (cs*rays(1,n)+ss*rays(5,n))*Scxlt
          pos(2) = rays(3,n)*Scxlt
          pos(3) = (-ss*rays(1,n) + cs*rays(5,n))*Scxlt
          pos(4) = tg

          !get external field from all overlaped fields at one location
          call getfldt_BeamLineElem(beamelem(idbd),pos,&
          extfld,fldmap(idrfile(3,idbd)))

          ex = extfld(1)*cs - extfld(3)*ss
          ey = extfld(2)
          ez = extfld(1)*ss + extfld(3)*cs
          bx = extfld(4)*cs - extfld(6)*ss
          by = extfld(5)
          bz = extfld(4)*ss + extfld(6)*cs

          !//advance the momenta of particles using implicit central
          !//difference scheme from Birdall and Longdon's book.
          umx = rays(2,n) + coefLz*ex*0.5d0*dt
          umy = rays(4,n) + coefLz*ey*0.5d0*dt
          umz = rays(6,n) + coefLz*ez*0.5d0*dt
          recpgamma = 1.0d0/sqrt(1.0d0+umx*umx+umy*umy+umz*umz)
          tmp = 0.5d0*Clight*dt*recpgamma*coefLz
          a1 = tmp*bx
          a2 = tmp*by
          a3 = tmp*bz
          a4 = 1.0d0+a1*a1+a2*a2+a3*a3
          s1 = umx + tmp*(umy*bz-umz*by)
          s2 = umy - tmp*(umx*bz-umz*bx)
          s3 = umz + tmp*(umx*by-umy*bx)
          upx = ((1.0d0+a1*a1)*s1+(a1*a2+a3)*s2+(a1*a3-a2)*s3)/a4
          upy = ((a1*a2-a3)*s1+(1.0d0+a2*a2)*s2+(a2*a3+a1)*s3)/a4
          upz = ((a1*a3+a2)*s1+(a2*a3-a1)*s2+(1.0d0+a3*a3)*s3)/a4
          rays(2,n) = upx + coefLz*ex*0.5d0*dt
          rays(4,n) = upy + coefLz*ey*0.5d0*dt
          rays(6,n) = upz + coefLz*ez*0.5d0*dt
        enddo

        t_ntrslo = t_ntrslo + elapsedtime_Timer( t0 )

        end subroutine kick2tBd0_BeamBunch

        !//release all the memory held by BeamBunch.
        subroutine destruct_BeamBunch(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(out) :: this

        deallocate(this%Pts1)

        end subroutine destruct_BeamBunch

      end module BeamBunchclass
