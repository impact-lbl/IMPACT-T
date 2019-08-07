!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! Outputclass: Output class in I/O module of CONTROL layer.
!
! MODULE  : ... Outputclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang
!
! DESCRIPTION: 
!> This class defines functions to print out the charged
!> particle beam information in the accelerator.
! Comments:
!----------------------------------------------------------------
      module Outputclass
        use Timerclass
        use BeamBunchclass
        use PhysConstclass

      contains
        !> calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        !> <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
        subroutine diagnostic1_Output(z,this)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), intent(inout) :: this
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet,gambet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(27) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax,z0gl,z0avg
        integer :: npctmin,npctmax

        call starttime_Timer(t0)

        qmc = this%Mass/1.0e6
        xl = Scxlt
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt
        !print*,"pts1: ",this%Pts1(1,1),this%Pts1(2,1),my_rank,innp

        den1 = 1.0/dble(nptot)
        den2 = den1*den1

        z0lc = 0.0
        do i = 1, innp
          z0lc = z0lc + this%Pts1(5,i)
        enddo
        call MPI_ALLREDUCE(z0lc,z0gl,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)
        z0avg = z0gl*den1

        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0

        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          x0lc3 = x0lc3 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)
          x0lc4 = x0lc4 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)*&
                  this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          px0lc3 = px0lc3 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)
          px0lc4 = px0lc4 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)*&
                   this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          y0lc3 = y0lc3 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)
          y0lc4 = y0lc4 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)*&
                  this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          py0lc3 = py0lc3 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)
          py0lc4 = py0lc4 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)*&
                   this%Pts1(4,i)
          sqsum5local = sqsum5local + (this%Pts1(5,i)-z0avg)**2
          z0lc3 = z0lc3 + abs((this%Pts1(5,i)-z0avg)**3)
          z0lc4 = z0lc4 + (this%Pts1(5,i)-z0avg)**4
          zpzlocal = zpzlocal + (this%Pts1(5,i)-z0avg)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc3 = pz0lc3 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc4 = pz0lc4 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)*&
                            this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
            if(localmax(5).lt.abs(this%Pts1(5,i)-z0avg)) then
               localmax(5) = abs(this%Pts1(5,i)-z0avg)
            endif
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        
        call MPI_REDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        gambet = sqrt(gam**2-1.0)
        bet = sqrt(gam**2-1.0)/gam
        energy = (gam-1.)*qmc

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
!          sqsum5 = sqz - z0*z0
          sqsum5 = sqz 
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
!          zpz = tmpgl(15)*den1 - z0*pz0
          zpz = tmpgl(15)*den1 
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
!          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z03 = cubz**(1.0/3.0)
!          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          z04 = sqrt(sqrt(fthz))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(abs(sqsum1))
          pxrms = sqrt(abs(sqsum2))
          yrms = sqrt(abs(sqsum3))
          pyrms = sqrt(abs(sqsum4))
          zrms = sqrt(abs(sqsum5))
          pzrms = sqrt(abs(sqsum6))
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          gam = sqrt(1.0+px0**2+py0**2+pz0**2)
          energy = qmc*(gam-1.0)
          bet = sqrt(1.0-(1.0/gam)**2)
          write(18,99)z,z0avg*xl,gam,energy,bet,sqrt(glrmax)*xl
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl
!          write(26,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz/epz,epz*xl
          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx,epx*xl
          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy,epy*xl
          write(26,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz,epz*xl

          write(27,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
          write(28,101)z,npctmin,npctmax,nptot
          write(29,100)z,x03*xl,px03,y03*xl,py03,z03*xl,&
                       pz03
          write(30,100)z,x04*xl,px04,y04*xl,py04,z04*xl,&
                       pz04

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(30)
        endif

99      format(6(1x,e13.6))
100      format(7(1x,e13.6))
101     format(1x,e13.6,3I10)

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1_Output

        !> calculate averaged <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        !> <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance from
        !> multiple bunch/bin.
        subroutine diagnostic1avg_Output(z,this,Nbunch)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), dimension(:), intent(inout) :: this
        integer, intent(in) :: Nbunch
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(29) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax,z0gl,z0avg,testmax,pz0avg
        double precision :: gamlc,gam2lc,gam2avg,tmpgam,gamdel
        integer :: npctmin,npctmax,ib,innpmb,i1,i2
        real*8 :: gamavg,gamgl
        real*8, dimension(3) :: tmp3lc,tmp3gl

        call starttime_Timer(t0)

        qmc = this(1)%Mass/1.0e6
        xl = Scxlt
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        nptot = 0
        innpmb = 0
        z0lc = 0.0
        gamlc = 0.0d0
        pz0lc = 0.0
        tmp3lc = 0.0d0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          innpmb = innpmb + innp
          nptot = nptot + this(ib)%Npt
          do i = 1, innp
            tmp3lc(1) = tmp3lc(1) + this(ib)%Pts1(5,i)
            tmp3lc(2) = tmp3lc(2) + this(ib)%Pts1(6,i)
            tmp3lc(3) = tmp3lc(3) + sqrt(1.0+this(ib)%Pts1(2,i)**2+&
               this(ib)%Pts1(4,i)**2+this(ib)%Pts1(6,i)**2)
          enddo
        enddo

        call MPI_ALLREDUCE(tmp3lc,tmp3gl,3,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)

        den1 = 1.0d0/nptot
        den2 = den1*den1
        z0avg = tmp3gl(1)/nptot
        pz0avg = tmp3gl(2)/nptot
        gamavg = tmp3gl(3)/nptot
        !print*,"z0avg: ",z0avg,pz0avg,gamavg


        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0
        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            !localmax(i) = abs(this(1)%Pts1(i,1))
            localmax(i) = 0.0
          enddo
          !lcrmax = this(1)%Pts1(1,1)**2+this(1)%Pts1(3,1)**2
          lcrmax = 0.0
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
!        testmax = 0.0
        gamlc = 0.0
        gam2lc = 0.0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          do i = 1, innp
            tmpgam = sqrt(1.0+this(ib)%Pts1(2,i)**2+&
               this(ib)%Pts1(4,i)**2+this(ib)%Pts1(6,i)**2)
            gamlc = gamlc + tmpgam
            !gam2lc = gam2lc + tmpgam**2
            gam2lc = gam2lc + (tmpgam-gamavg)**2
            x0lc = x0lc + this(ib)%Pts1(1,i)
            sqsum1local = sqsum1local + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
            x0lc3 = x0lc3 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
            x0lc4 = x0lc4 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*&
                    this(ib)%Pts1(1,i)
            xpxlocal = xpxlocal + this(ib)%Pts1(1,i)*this(ib)%Pts1(2,i)
            px0lc = px0lc + this(ib)%Pts1(2,i)
            sqsum2local = sqsum2local + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
            px0lc3 = px0lc3 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
            px0lc4 = px0lc4 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*&
                     this(ib)%Pts1(2,i)
            y0lc = y0lc + this(ib)%Pts1(3,i)
            sqsum3local = sqsum3local + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
            y0lc3 = y0lc3 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
            y0lc4 = y0lc4 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*&
                    this(ib)%Pts1(3,i)
            ypylocal = ypylocal + this(ib)%Pts1(3,i)*this(ib)%Pts1(4,i)
            py0lc = py0lc + this(ib)%Pts1(4,i)
            sqsum4local = sqsum4local + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc3 = py0lc3 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc4 = py0lc4 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*&
                     this(ib)%Pts1(4,i)
            sqsum5local = sqsum5local + (this(ib)%Pts1(5,i)-z0avg)**2
            z0lc3 = z0lc3 + abs((this(ib)%Pts1(5,i)-z0avg)**3)
            z0lc4 = z0lc4 + (this(ib)%Pts1(5,i)-z0avg)**4
            zpzlocal = zpzlocal + (this(ib)%Pts1(5,i)-z0avg)*(this(ib)%Pts1(6,i)-pz0avg)
            pz0lc = pz0lc + this(ib)%Pts1(6,i)
            sqsum6local = sqsum6local + (this(ib)%Pts1(6,i)-pz0avg)**2 
            pz0lc3 = pz0lc3 + (abs(this(ib)%Pts1(6,i)-pz0avg)**3) 
            pz0lc4 = pz0lc4 + (this(ib)%Pts1(6,i)-pz0avg)**4 
            do j = 1, 4
              if(localmax(j).lt.abs(this(ib)%Pts1(j,i))) then
                 localmax(j) = abs(this(ib)%Pts1(j,i))
              endif
            enddo
            if(localmax(5).lt.abs(this(ib)%Pts1(5,i)-z0avg)) then
               localmax(5) = abs(this(ib)%Pts1(5,i)-z0avg)
               i1 = i
            endif
!            if(testmax .lt. abs(this(ib)%Pts1(5,i)) ) then
!              testmax = this(ib)%Pts1(5,i)
!              i2 = i
!              print*,"testmax: ",i,testmax,this(ib)%Pts1(5,i),&
!                        abs(this(ib)%Pts1(5,i))
!            endif
            if(localmax(6).lt.abs(this(ib)%Pts1(6,i)-pz0avg)) then
                 localmax(6) = abs(this(ib)%Pts1(6,i)-pz0avg)
            endif
            if(lcrmax.lt.(this(ib)%Pts1(1,i)**2+this(ib)%Pts1(3,i)**2)) then
              lcrmax = this(ib)%Pts1(1,i)**2 + this(ib)%Pts1(3,i)**2
            endif
          enddo
        enddo

        !print*,"z0avg: ",z0avg*xl,localmax(5)*xl,&
        !       this(1)%Pts1(5,i1)*xl
        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        tmplc(28) = gamlc
        tmplc(29) = gam2lc
        
        call MPI_REDUCE(tmplc,tmpgl,29,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
!          sqsum5 = sqz - z0*z0
          sqsum5 = sqz 
          sqpz = tmpgl(12)*den1
!          sqsum6 = sqpz - pz0*pz0
          sqsum6 = sqpz
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
!          zpz = tmpgl(15)*den1 - z0*pz0
          zpz = tmpgl(15)*den1 
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
!          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z03 = cubz**(1.0/3.0)
!          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          z04 = sqrt(sqrt(fthz))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = cubpz**(1.0/3.0)
          pz04 = sqrt(sqrt(fthpz))
          !pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          !pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
          !       3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(abs(sqsum1))
          pxrms = sqrt(abs(sqsum2))
          yrms = sqrt(abs(sqsum3))
          pyrms = sqrt(abs(sqsum4))
          zrms = sqrt(abs(sqsum5))
          pzrms = sqrt(abs(sqsum6))
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          gam = tmpgl(28)*den1
!          gam = sqrt(1.0+px0**2+py0**2+pz0**2)
          energy = qmc*(gam-1.0)
          bet = sqrt(1.0-(1.0/gam)**2)
          gam2avg = tmpgl(29)*den1
          !gamdel = sqrt(abs(gam2avg - gam**2))
          gamdel = sqrt(gam2avg)
          write(18,100)z,z0avg*xl,gam,energy,bet,sqrt(glrmax)*xl,gamdel
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl
!          write(26,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz/epz,epz*xl
          write(24,102)z,z0avg*xl,x0*xl,xrms*xl,px0,pxrms,-xpx*xl,epx*xl
          write(25,102)z,z0avg*xl,y0*xl,yrms*xl,py0,pyrms,-ypy*xl,epy*xl
          write(26,100)z,z0avg*xl,zrms*xl,pz0,pzrms,-zpz*xl,epz*xl

          write(27,102)z,z0avg*xl,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
          write(28,101)z,z0avg*xl,npctmin,npctmax,nptot
          write(29,102)z,z0avg*xl,x03*xl,px03,y03*xl,py03,z03*xl,&
                       pz03
          write(30,102)z,z0avg*xl,x04*xl,px04,y04*xl,py04,z04*xl,&
                       pz04

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(30)
        endif

99      format(6(1x,e16.8))
100      format(7(1x,e18.10))
101     format(1x,e16.8,e16.8,3I10)
102      format(8(1x,e16.8))

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1avg_Output

        !> calculate averaged <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        !> <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance from
        !> multiple bunch/bin.
        subroutine diagnostic1avgB_Output(z,this,Nbunch)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), dimension(:), intent(inout) :: this
        integer, intent(in) :: Nbunch
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(29) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax,z0gl,z0avg,testmax,pz0avg
        double precision :: gamlc,gam2lc,gam2avg,tmpgam,gamdel
        integer :: npctmin,npctmax,ib,innpmb,i1,i2
        real*8 :: gamavg,gamgl
        real*8, dimension(3) :: tmp3lc,tmp3gl

        call starttime_Timer(t0)

        qmc = this(1)%Mass/1.0e6
        xl = Scxlt
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        nptot = 0
        innpmb = 0
        z0lc = 0.0
        gamlc = 0.0d0
        pz0lc = 0.0
        tmp3lc = 0.0d0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          innpmb = innpmb + innp
          nptot = nptot + this(ib)%Npt
          do i = 1, innp
            tmp3lc(1) = tmp3lc(1) + this(ib)%Pts1(5,i)
            tmp3lc(2) = tmp3lc(2) + this(ib)%Pts1(6,i)
            tmp3lc(3) = tmp3lc(3) + sqrt(1.0+this(ib)%Pts1(2,i)**2+&
               this(ib)%Pts1(4,i)**2+this(ib)%Pts1(6,i)**2)
          enddo
        enddo

        call MPI_ALLREDUCE(tmp3lc,tmp3gl,3,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)

        den1 = 1.0d0/nptot
        den2 = den1*den1
        z0avg = tmp3gl(1)/nptot
        pz0avg = tmp3gl(2)/nptot
        gamavg = tmp3gl(3)/nptot
        !print*,"z0avg: ",z0avg,pz0avg,gamavg


        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0
        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            !localmax(i) = abs(this(1)%Pts1(i,1))
            localmax(i) = 0.0
          enddo
          !lcrmax = this(1)%Pts1(1,1)**2+this(1)%Pts1(3,1)**2
          lcrmax = 0.0
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
!        testmax = 0.0
        gamlc = 0.0
        gam2lc = 0.0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          do i = 1, innp
            tmpgam = sqrt(1.0+this(ib)%Pts1(2,i)**2+&
               this(ib)%Pts1(4,i)**2+this(ib)%Pts1(6,i)**2)
            gamlc = gamlc + tmpgam
            !gam2lc = gam2lc + tmpgam**2
            gam2lc = gam2lc + (tmpgam-gamavg)**2
            x0lc = x0lc + this(ib)%Pts1(1,i)
            sqsum1local = sqsum1local + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
            x0lc3 = x0lc3 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
            x0lc4 = x0lc4 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*&
                    this(ib)%Pts1(1,i)
            xpxlocal = xpxlocal + this(ib)%Pts1(1,i)*this(ib)%Pts1(2,i)
            px0lc = px0lc + this(ib)%Pts1(2,i)
            sqsum2local = sqsum2local + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
            px0lc3 = px0lc3 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
            px0lc4 = px0lc4 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*&
                     this(ib)%Pts1(2,i)
            y0lc = y0lc + this(ib)%Pts1(3,i)
            sqsum3local = sqsum3local + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
            y0lc3 = y0lc3 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
            y0lc4 = y0lc4 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*&
                    this(ib)%Pts1(3,i)
            ypylocal = ypylocal + this(ib)%Pts1(3,i)*this(ib)%Pts1(4,i)
            py0lc = py0lc + this(ib)%Pts1(4,i)
            sqsum4local = sqsum4local + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc3 = py0lc3 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc4 = py0lc4 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*&
                     this(ib)%Pts1(4,i)
            sqsum5local = sqsum5local + (this(ib)%Pts1(5,i)-z0avg)**2
            z0lc3 = z0lc3 + abs((this(ib)%Pts1(5,i)-z0avg)**3)
            z0lc4 = z0lc4 + (this(ib)%Pts1(5,i)-z0avg)**4
            zpzlocal = zpzlocal + (this(ib)%Pts1(5,i)-z0avg)*(this(ib)%Pts1(6,i)-pz0avg)
            pz0lc = pz0lc + this(ib)%Pts1(6,i)
            sqsum6local = sqsum6local + (this(ib)%Pts1(6,i)-pz0avg)**2 
            pz0lc3 = pz0lc3 + (abs(this(ib)%Pts1(6,i)-pz0avg)**3) 
            pz0lc4 = pz0lc4 + (this(ib)%Pts1(6,i)-pz0avg)**4 
            do j = 1, 4
              if(localmax(j).lt.abs(this(ib)%Pts1(j,i))) then
                 localmax(j) = abs(this(ib)%Pts1(j,i))
              endif
            enddo
            if(localmax(5).lt.abs(this(ib)%Pts1(5,i)-z0avg)) then
               localmax(5) = abs(this(ib)%Pts1(5,i)-z0avg)
               i1 = i
            endif
!            if(testmax .lt. abs(this(ib)%Pts1(5,i)) ) then
!              testmax = this(ib)%Pts1(5,i)
!              i2 = i
!              print*,"testmax: ",i,testmax,this(ib)%Pts1(5,i),&
!                        abs(this(ib)%Pts1(5,i))
!            endif
            if(localmax(6).lt.abs(this(ib)%Pts1(6,i)-pz0avg)) then
                 localmax(6) = abs(this(ib)%Pts1(6,i)-pz0avg)
            endif
            if(lcrmax.lt.(this(ib)%Pts1(1,i)**2+this(ib)%Pts1(3,i)**2)) then
              lcrmax = this(ib)%Pts1(1,i)**2 + this(ib)%Pts1(3,i)**2
            endif
          enddo
        enddo

        !print*,"z0avg: ",z0avg*xl,localmax(5)*xl,&
        !       this(1)%Pts1(5,i1)*xl
        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        tmplc(28) = gamlc
        tmplc(29) = gam2lc
        
        call MPI_REDUCE(tmplc,tmpgl,29,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
!          sqsum5 = sqz - z0*z0
          sqsum5 = sqz 
          sqpz = tmpgl(12)*den1
!          sqsum6 = sqpz - pz0*pz0
          sqsum6 = sqpz
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
!          zpz = tmpgl(15)*den1 - z0*pz0
          zpz = tmpgl(15)*den1 
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
!          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z03 = cubz**(1.0/3.0)
!          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          z04 = sqrt(sqrt(fthz))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = cubpz**(1.0/3.0)
          pz04 = sqrt(sqrt(fthpz))
          !pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          !pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
          !       3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(abs(sqsum1))
          pxrms = sqrt(abs(sqsum2))
          yrms = sqrt(abs(sqsum3))
          pyrms = sqrt(abs(sqsum4))
          zrms = sqrt(abs(sqsum5))
          pzrms = sqrt(abs(sqsum6))
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          gam = tmpgl(28)*den1
!          gam = sqrt(1.0+px0**2+py0**2+pz0**2)
          energy = qmc*(gam-1.0)
          bet = sqrt(1.0-(1.0/gam)**2)
          gam2avg = tmpgl(29)*den1
          !gamdel = sqrt(abs(gam2avg - gam**2))
          gamdel = sqrt(gam2avg)


!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl
!          write(26,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz/epz,epz*xl
          write(34,102)z,z0avg*xl,x0*xl,xrms*xl,px0,pxrms,-xpx*xl,epx*xl
          write(35,102)z,z0avg*xl,y0*xl,yrms*xl,py0,pyrms,-ypy*xl,epy*xl
          write(36,100)z,z0avg*xl,zrms*xl,pz0,pzrms,-zpz*xl,epz*xl

          write(37,102)z,z0avg*xl,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
!          write(38,100)z,z0avg*xl,gam,energy,bet,sqrt(glrmax)*xl,gamdel
          write(38,100)z,this(1)%refptcl(1)*xl,this(1)%refptcl(2),this(1)%refptcl(3)*xl,&
                      this(1)%refptcl(4),this(1)%refptcl(5)*xl,this(1)%refptcl(6)
!          write(38,101)z,z0avg*xl,npctmin,npctmax,nptot
!          write(29,102)z,z0avg*xl,x03*xl,px03,y03*xl,py03,z03*xl,&
!                       pz03

!          call flush(18)
          call flush(34)
          call flush(35)
          call flush(36)
          call flush(37)
          call flush(38)
!          call flush(48)
!          call flush(29)
        endif

99      format(6(1x,e16.8))
100      format(7(1x,e16.8))
101     format(1x,e16.8,e16.8,3I10)
102      format(8(1x,e16.8))

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1avgB_Output

        !> calculate averaged <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        !> <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance from
        !> multiple bunch/bin at fixed z (i.e. bunch center).
        subroutine diagnostic1avgZ_Output(z,this,Nbunch)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), dimension(:), intent(inout) :: this
        integer, intent(in) :: Nbunch
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(29) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax,z0gl,z0avg,testmax
        double precision :: gamlc,gam2lc,gam2avg,tmpgam,gamdel
        double precision :: tmpx,tmpy,deltat,recpgamma
        integer :: npctmin,npctmax,ib,innpmb,i1,i2

        call starttime_Timer(t0)

        qmc = this(1)%Mass/1.0e6
        xl = Scxlt
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        nptot = 0
        innpmb = 0
        z0lc = 0.0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          innpmb = innpmb + innp
          nptot = nptot + this(ib)%Npt
          do i = 1, innp
            z0lc = z0lc + this(ib)%Pts1(5,i)
          enddo
        enddo
        call MPI_ALLREDUCE(z0lc,z0gl,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)
        den1 = 1.0/dble(nptot)
        den2 = den1*den1
        z0avg = z0gl*den1

        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0
        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            !localmax(i) = abs(this(1)%Pts1(i,1))
            localmax(i) = 0.0
          enddo
          !lcrmax = this(1)%Pts1(1,1)**2+this(1)%Pts1(3,1)**2
          lcrmax = 0.0
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
!        testmax = 0.0
        gamlc = 0.0
        gam2lc = 0.0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          do i = 1, innp
            tmpgam = sqrt(1.0+this(ib)%Pts1(2,i)**2+&
               this(ib)%Pts1(4,i)**2+this(ib)%Pts1(6,i)**2)
            gamlc = gamlc + tmpgam  
            gam2lc = gam2lc + tmpgam**2
            recpgamma = 1.0/tmpgam
            deltat = (this(ib)%Pts1(5,i)-z0avg)&
                     /(recpgamma*this(ib)%Pts1(6,i))
            !drift to the "fixed z" location
            tmpx = this(ib)%Pts1(1,i)-recpgamma*this(ib)%Pts1(2,i)*deltat
            x0lc = x0lc + tmpx
            sqsum1local = sqsum1local + tmpx*tmpx
            x0lc3 = x0lc3 + tmpx*tmpx*tmpx
            x0lc4 = x0lc4 + tmpx*tmpx*tmpx*tmpx
!            xpxlocal = xpxlocal + tmpx*this(ib)%Pts1(2,i)/this(ib)%Pts1(6,i)
            xpxlocal = xpxlocal + tmpx*this(ib)%Pts1(2,i)
!            x0lc = x0lc + this(ib)%Pts1(1,i)
!            sqsum1local = sqsum1local + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
!            x0lc3 = x0lc3 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
!            x0lc4 = x0lc4 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*&
!                    this(ib)%Pts1(1,i)
!            xpxlocal = xpxlocal + this(ib)%Pts1(1,i)*this(ib)%Pts1(2,i)
            px0lc = px0lc + this(ib)%Pts1(2,i)
!            px0lc = px0lc + this(ib)%Pts1(2,i)/this(ib)%Pts1(6,i)
            sqsum2local = sqsum2local + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
!            sqsum2local = sqsum2local + (this(ib)%Pts1(2,i)/this(ib)%Pts1(6,i))**2
            px0lc3 = px0lc3 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
            px0lc4 = px0lc4 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*&
                     this(ib)%Pts1(2,i)
            !drift to the "fixed z" location
            tmpy = this(ib)%Pts1(3,i)-recpgamma*this(ib)%Pts1(4,i)*deltat
            y0lc = y0lc + tmpy
            sqsum3local = sqsum3local + tmpy*tmpy
            y0lc3 = y0lc3 + tmpy*tmpy*tmpy
            y0lc4 = y0lc4 + tmpy*tmpy*tmpy*tmpy
            ypylocal = ypylocal + tmpy*this(ib)%Pts1(4,i)
!            y0lc = y0lc + this(ib)%Pts1(3,i)
!            sqsum3local = sqsum3local + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
!            y0lc3 = y0lc3 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
!            y0lc4 = y0lc4 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*&
!                    this(ib)%Pts1(3,i)
!            ypylocal = ypylocal + this(ib)%Pts1(3,i)*this(ib)%Pts1(4,i)
            py0lc = py0lc + this(ib)%Pts1(4,i)
            sqsum4local = sqsum4local + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc3 = py0lc3 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc4 = py0lc4 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*&
                     this(ib)%Pts1(4,i)
            sqsum5local = sqsum5local + (this(ib)%Pts1(5,i)-z0avg)**2
            z0lc3 = z0lc3 + abs((this(ib)%Pts1(5,i)-z0avg)**3)
            z0lc4 = z0lc4 + (this(ib)%Pts1(5,i)-z0avg)**4
            zpzlocal = zpzlocal + (this(ib)%Pts1(5,i)-z0avg)*this(ib)%Pts1(6,i)
            pz0lc = pz0lc + this(ib)%Pts1(6,i)
            sqsum6local = sqsum6local + this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)
            pz0lc3 = pz0lc3 + this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)
            pz0lc4 = pz0lc4 + this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)*&
                              this(ib)%Pts1(6,i)
            do j = 1, 4
              if(localmax(j).lt.abs(this(ib)%Pts1(j,i))) then
                 localmax(j) = abs(this(ib)%Pts1(j,i))
              endif
            enddo
            if(localmax(5).lt.abs(this(ib)%Pts1(5,i)-z0avg)) then
               localmax(5) = abs(this(ib)%Pts1(5,i)-z0avg)
               i1 = i
            endif
!            if(testmax .gt. this(ib)%Pts1(5,i) ) then
!              testmax = this(ib)%Pts1(5,i)
!              i2 = i
!            endif
            if(localmax(6).lt.abs(this(ib)%Pts1(6,i))) then
                 localmax(6) = abs(this(ib)%Pts1(6,i))
            endif
            if(lcrmax.lt.(this(ib)%Pts1(1,i)**2+this(ib)%Pts1(3,i)**2)) then
              lcrmax = this(ib)%Pts1(1,i)**2 + this(ib)%Pts1(3,i)**2
            endif
          enddo
        enddo

!        print*,"z0avg: ",z0avg*xl,localmax(5)*xl,testmax*xl,xl,&
!               this(1)%Pts1(5,i1)*xl,this(1)%Pts1(5,i2)*xl
        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        tmplc(28) = gamlc
        tmplc(29) = gam2lc
        
        call MPI_REDUCE(tmplc,tmpgl,29,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
!          sqsum5 = sqz - z0*z0
          sqsum5 = sqz 
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
!          zpz = tmpgl(15)*den1 - z0*pz0
          zpz = tmpgl(15)*den1 
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
!          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z03 = cubz**(1.0/3.0)
!          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          z04 = sqrt(sqrt(fthz))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(abs(sqsum1))
          pxrms = sqrt(abs(sqsum2))
          yrms = sqrt(abs(sqsum3))
          pyrms = sqrt(abs(sqsum4))
          zrms = sqrt(abs(sqsum5))
          pzrms = sqrt(abs(sqsum6))
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          gam = tmpgl(28)*den1
!          gam = sqrt(1.0+px0**2+py0**2+pz0**2)
          energy = qmc*(gam-1.0)
          bet = sqrt(1.0-(1.0/gam)**2)
          gam2avg = tmpgl(29)*den1
          gamdel = sqrt(abs(gam2avg - gam**2))
          write(18,100)z,z0avg*xl,gam,energy,bet,sqrt(glrmax)*xl,gamdel
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl
!          write(26,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz/epz,epz*xl
          write(24,102)z,z0*xl,x0*xl,xrms*xl,px0,pxrms,-xpx*xl,epx*xl
          write(25,102)z,z0*xl,y0*xl,yrms*xl,py0,pyrms,-ypy*xl,epy*xl
          write(26,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz*xl,epz*xl

          write(27,102)z,z0*xl,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
          write(28,101)z,z0*xl,npctmin,npctmax,nptot
          write(29,102)z,z0*xl,x03*xl,px03,y03*xl,py03,z03*xl,&
                       pz03
          write(30,102)z,z0*xl,x04*xl,px04,y04*xl,py04,z04*xl,&
                       pz04

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(30)
        endif

99      format(6(1x,e16.8))
100      format(7(1x,e16.8))
101     format(1x,e16.8,e16.8,3I10)
102      format(8(1x,e16.8))

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1avgZ_Output

        !> the 6D phase space output has (x(m), px/mc, y(m), py/mc, z(m), pz/mc).
        subroutine phase_Output(nfile,this,samplePeriod)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: samplePeriod
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf

        if (samplePeriod .eq. 0) then
           samplePeriod = 1
        endif

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(6,mnpt))
        sixnpt = 6*this%Nptlocal

        call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 6*nptlist

        if(my_rank.eq.0) then
          open(nfile,status='unknown')
          do i = 1, this%Nptlocal,samplePeriod
!            write(nfile,100)this%Pts1(1,i),this%Pts1(2,i),this%Pts1(3,i),&
!                            this%Pts1(4,i),this%Pts1(5,i),this%Pts1(6,i)
            write(nfile,100)this%Pts1(1,i)*Scxlt,this%Pts1(2,i),&
                            this%Pts1(3,i)*Scxlt,&
                            this%Pts1(4,i),this%Pts1(5,i)*Scxlt,&
                            this%Pts1(6,i)
          enddo
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
        
            do j = 1, nptlist(i)/6,samplePeriod
!              write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
!                              recvbuf(4,j),recvbuf(5,j),recvbuf(6,j)
              write(nfile,100)recvbuf(1,j)*Scxlt,recvbuf(2,j),&
                              recvbuf(3,j)*Scxlt,&
                              recvbuf(4,j),recvbuf(5,j)*Scxlt,recvbuf(6,j)
            enddo
          enddo
          close(nfile)
          call flush(nfile)
        else
          call MPI_SEND(this%Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

!100     format(6(1x,e17.9))
100     format(6(1x,e20.12))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phase_Output

        !> output 2D particle number density.
        subroutine dens2d_Output(nstep,nfile,this,totnptcls,xmnin,xmxin,&
        pxmnin,pxmxin,ymnin,ymxin,pymnin,pymxin,zmnin,zmxin,pzmnin,pzmxin)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: totnptcls
        double precision, intent(in) :: xmnin,xmxin,pxmnin,pxmxin,ymnin,&
        ymxin,pymnin,pymxin,zmnin,zmxin,pzmnin,pzmxin
        integer, parameter :: nxcell = 128
        integer, parameter :: npxcell = 128
        integer, parameter :: nycell = 128
        integer, parameter :: npycell = 128
        integer, parameter :: nzcell = 128
        integer, parameter :: npzcell = 128
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,k,ii,req,jj,kk,icc
        integer :: ix,ipx,iy,ipy,iz,ipz,ic,nsend,nfile1,nfile2
        integer :: nfile3,nfile4,nfile5
        integer, dimension(nxcell,npxcell) :: count1
        integer, dimension(nxcell,npxcell) :: count2
        integer, dimension(nxcell,npxcell) :: count3
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(nxcell,nzcell) :: count5
        integer, dimension(nycell,nzcell) :: count6
        integer, dimension(3,nxcell*npxcell) :: sendbuf1,recvbuf1
        integer, dimension(3,nycell*npycell) :: sendbuf2,recvbuf2
        integer, dimension(3,nzcell*npzcell) :: sendbuf3,recvbuf3
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, dimension(3,nxcell*nzcell) :: sendbuf5,recvbuf5
        integer, dimension(3,nycell*nzcell) :: sendbuf6,recvbuf6
        integer, allocatable, dimension(:) :: nptlist
        double precision, dimension(6) :: range,prange
        double precision, dimension(6) :: lcrange1,lcrange2,temp1
        double precision :: hxx,hyy,hzz,hpxx,hpyy,hpzz,invvol
        double precision ::eps,tmp1,tmp2,tmp3,xmax,xmin,pxmax,pxmin,&
        ymax,ymin,pymax,pymin,zmin,zmax,pzmin,pzmax
        character*3 name
        integer :: ioerr
        integer gethostname

!        ioerr = gethostname(name,3)
        eps = 1.0e-4
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange1(i+3) = this%Pts1(2*i,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
            lcrange2(i+3) = this%Pts1(2*i,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange1(i+3) = 5000.0
            lcrange2(i) = -5000.0
            lcrange2(i+3) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
            ii = i + 3
            if(lcrange1(ii).gt.this%Pts1(2*i,j)) then
              lcrange1(ii) = this%Pts1(2*i,j)
            endif
            if(lcrange2(ii).lt.this%Pts1(2*i,j)) then
              lcrange2(ii) = this%Pts1(2*i,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,6,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)

        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo 
        do i = 1, 3
          prange(2*i-1) = temp1(i+3) + temp1(i+3)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,6,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo 
        do i = 1, 3
          prange(2*i) = temp1(i+3) + temp1(i+3)*eps
        enddo 

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        zmin = min(zmnin,range(5))
        zmax = max(zmxin,range(6))
        pxmin = min(pxmnin,prange(1))
        pxmax = max(pxmxin,prange(2))
        pymin = min(pymnin,prange(3))
        pymax = max(pymxin,prange(4))
        pzmin = min(pzmnin,prange(5))
        pzmax = max(pzmxin,prange(6))

        hxx = (xmax-xmin)/nxcell
        hpxx = (pxmax-pxmin)/npxcell
        hyy = (ymax-ymin)/nycell
        hpyy = (pymax-pymin)/npycell
        hzz = (zmax-zmin)/nzcell
        hpzz = (pzmax-pzmin)/npzcell

        count1 = 0 
        count2 = 0 
        count3 = 0 
        count4 = 0
        count5 = 0
        count6 = 0
        do i = 1, this%Nptlocal
          ix = int((this%Pts1(1,i)-range(1))/hxx) + 1
          ipx = int((this%Pts1(2,i)-prange(1))/hpxx) + 1
          count1(ix,ipx) = count1(ix,ipx) + 1
          iy = int((this%Pts1(3,i)-range(3))/hyy) + 1
          ipy = int((this%Pts1(4,i)-prange(3))/hpyy) + 1
          count2(iy,ipy) = count2(iy,ipy) + 1
          iz = int((this%Pts1(5,i)-range(5))/hzz) + 1
          ipz = int((this%Pts1(6,i)-prange(5))/hpzz) + 1
          count3(iz,ipz) = count3(iz,ipz) + 1
          count4(ix,iy) = count4(ix,iy) + 1
          count5(ix,iz) = count5(ix,iz) + 1
          count6(iy,iz) = count6(iy,iz) + 1
        enddo

        invvol = 1.0/(hxx*hpxx)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npxcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf1(3,ic) = count1(i,j)
            enddo
          enddo
        else
          do j = 1, npxcell
            do i = 1, nxcell
              if(count1(i,j).gt.0) then
                icc = icc + 1
                sendbuf1(1,icc) = i
                sendbuf1(2,icc) = j
                sendbuf1(3,icc) = count1(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        allocate(nptlist(0:np-1))
        nptlist = 0

        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf1(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf1(1,j)
              jj = recvbuf1(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf1(3,icc) = sendbuf1(3,icc) + recvbuf1(3,j)
            enddo
          enddo

!          open(unit=nfile,file="/n/"//name//&
!               "/scratch/jiqiang/XPx.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile,file='XPx.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,prange(1)+0.5*hpxx,prange(1)+npxcell*hpxx-&
           0.5*hpxx
          do j = 1, npxcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = prange(1) + j*hpxx - 0.5*hpxx
!              write(nfile)sendbuf1(3,ic)
!              write(nfile)tmp1,tmp2,sendbuf1(3,ic)*tmp3
              write(nfile,100)tmp1,tmp2,sendbuf1(3,ic)*tmp3
            enddo
          enddo
          close(nfile)
        else
          call MPI_ISEND(sendbuf1(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hyy*hpyy)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npycell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              sendbuf2(3,ic) = count2(i,j) 
            enddo
          enddo
        else
          do j = 1, npycell
            do i = 1, nycell
              if(count2(i,j).gt.0) then
                icc = icc + 1
                sendbuf2(1,icc) = i
                sendbuf2(2,icc) = j
                sendbuf2(3,icc) = count2(i,j) 
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf2(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf2(1,j)
              jj = recvbuf2(2,j)
              icc = (jj-1)*nycell+ii
              sendbuf2(3,icc) = sendbuf2(3,icc) + recvbuf2(3,j)
            enddo
          enddo

          nfile1 = nfile + 1
!          open(unit=nfile1,file="/n/"//name//&
!               "/scratch/jiqiang/YPy.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile1,file='YPy.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile1,*)nstep,range(3)+0.5*hyy,range(3)+nycell*hyy &
           -0.5*hyy,prange(3)+0.5*hpyy,prange(3)+npycell*hpyy-&
           0.5*hpyy
          do j = 1, npycell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              tmp1 = range(3) + i*hyy - 0.5*hyy
              tmp2 = prange(3) + j*hpyy - 0.5*hpyy
!              write(nfile1)sendbuf2(3,ic)
!              write(nfile1)tmp1,tmp2,sendbuf2(3,ic)*tmp3
              write(nfile1,100)tmp1,tmp2,sendbuf2(3,ic)*tmp3
            enddo
          enddo
          close(nfile1)
        else
          call MPI_ISEND(sendbuf2(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hzz*hpzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npzcell
            do i = 1, nzcell
              ic = (j-1)*nzcell + i
              sendbuf3(3,ic) = count3(i,j) 
            enddo
          enddo
        else
          do j = 1, npzcell
            do i = 1, nzcell
              if(count3(i,j).gt.0) then
                icc = icc + 1
                sendbuf3(1,icc) = i
                sendbuf3(2,icc) = j
                sendbuf3(3,icc) = count3(i,j) 
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf3(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf3(1,j)
              jj = recvbuf3(2,j)
              icc = (jj-1)*nzcell + ii
              sendbuf3(3,icc) = sendbuf3(3,icc) + recvbuf3(3,j)
            enddo
          enddo

          nfile2 = nfile + 2
!          open(unit=nfile2,file="/n/"//name//&
!               "/scratch/jiqiang/ZPz.data",status="unknown",&
!               position="append",form="unformatted")
         open(unit=nfile2,file='ZPz.data',status='unknown',&
              position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile2,*)nstep,range(5)+0.5*hzz,range(5)+nzcell*hzz &
           -0.5*hzz,prange(5)+0.5*hpzz,prange(5)+npzcell*hpzz-&
           0.5*hpzz
          do j = 1, npzcell
            do i = 1, nzcell
              ic = (j-1)*nzcell + i
              tmp1 = range(5) + i*hzz - 0.5*hzz
              tmp2 = prange(5) + j*hpzz - 0.5*hpzz
!              write(nfile2)sendbuf3(3,ic)
!              write(nfile2)tmp1,tmp2,sendbuf3(3,ic)*tmp3
              write(nfile2,100)tmp1,tmp2,sendbuf3(3,ic)*tmp3
            enddo
          enddo
          close(nfile2)
        else
          call MPI_ISEND(sendbuf3(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hxx*hyy)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
!          open(unit=nfile3,file="/n/"//name//&
!                "/scratch/jiqiang/XY.data",status="unknown",&
!                position="append",form="unformatted")
          open(unit=nfile3,file='XY.data',status='unknown',&
                position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile3,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,range(3)+0.5*hyy,range(3)+nycell*hyy-&
           0.5*hyy
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = range(3) + j*hyy - 0.5*hyy
!              write(nfile3)sendbuf4(3,ic)
!              write(nfile3)tmp1,tmp2,sendbuf4(3,ic)*tmp3
              write(nfile3,100)tmp1,tmp2,sendbuf4(3,ic)*tmp3
            enddo
          enddo
          close(nfile3)
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hxx*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nzcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf5(3,ic) = count5(i,j)
            enddo
          enddo
        else
          do j = 1, nzcell
            do i = 1, nxcell
              if(count5(i,j).gt.0) then
                icc = icc + 1
                sendbuf5(1,icc) = i
                sendbuf5(2,icc) = j
                sendbuf5(3,icc) = count5(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf5(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf5(1,j)
              jj = recvbuf5(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf5(3,icc) = sendbuf5(3,icc) + recvbuf5(3,j)
            enddo
          enddo

          nfile4 = nfile + 4
!          open(unit=nfile4,file="/n/"//name//&
!               "/scratch/jiqiang/XZ.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile4,file='XZ.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile4,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,range(5)+0.5*hzz,range(5)+nzcell*hzz-&
           0.5*hzz
          do j = 1, nzcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = range(5) + j*hzz - 0.5*hzz
!              write(nfile4)sendbuf5(3,ic)
!              write(nfile4)tmp1,tmp2,sendbuf5(3,ic)*tmp3
              write(nfile4,100)tmp1,tmp2,sendbuf5(3,ic)*tmp3
            enddo
          enddo
          close(nfile4)
        else
          call MPI_ISEND(sendbuf5(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hyy*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nzcell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              sendbuf6(3,ic) = count6(i,j)
            enddo
          enddo
        else
          do j = 1, nzcell
            do i = 1, nycell
              if(count6(i,j).gt.0) then
                icc = icc + 1
                sendbuf6(1,icc) = i
                sendbuf6(2,icc) = j
                sendbuf6(3,icc) = count6(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf6(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf6(1,j)
              jj = recvbuf6(2,j)
              icc = (jj-1)*nycell + ii
              sendbuf6(3,icc) = sendbuf6(3,icc) + recvbuf6(3,j)
            enddo
          enddo

          nfile5 = nfile + 5
!          open(unit=nfile5,file="/n/"//name//&
!               "/scratch/jiqiang/YZ.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile5,file='YZ.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile5,*)nstep,range(3)+0.5*hyy,range(3)+nycell*hyy &
           -0.5*hyy,range(5)+0.5*hzz,range(5)+nzcell*hzz-&
           0.5*hzz
          do j = 1, nzcell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              tmp1 = range(3) + i*hyy - 0.5*hyy
              tmp2 = range(5) + j*hzz - 0.5*hzz
!              write(nfile5)sendbuf6(3,ic)
!              write(nfile5)tmp1,tmp2,sendbuf6(3,ic)*tmp3
              write(nfile5,100)tmp1,tmp2,sendbuf6(3,ic)*tmp3
            enddo
          enddo
          close(nfile5)
        else
          call MPI_ISEND(sendbuf6(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        deallocate(nptlist)

100     format(3(1x,e14.7))

        end subroutine dens2d_Output

        ! Terminate MPI
        subroutine end_Output(time)
        implicit none
        include 'mpif.h'
        double precision, intent(inout) :: time
        double precision :: endtime, mtime
        integer :: my_rank,ierr

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        endtime = MPI_WTIME()
        time = endtime - time
        call MPI_REDUCE(time,mtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        !for measurement of memory
        !call system_stats()

        if(my_rank.eq.0) then
          print*,"time: ",mtime
        endif

!        call MPI_Finalize(ierr)

        end subroutine end_Output

        subroutine inpoint_Output(nfile,this,z,inb,jstp,nprocrow,nproccol,&
               geom,nx,ny,nz,myidx,myidy,nptot,iout,itsz,isteer,islout,dtless)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: z,dtless
        type(CompDom), intent(inout) :: geom
        integer, intent(inout) :: inb,jstp
        integer, intent(inout) :: nptot,iout,itsz,isteer,islout
        integer, intent(in) :: nprocrow,nproccol,nx,ny,nz,myidx,myidy
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        integer :: nplc,ierr
        character*6 name1
        character*7 name2
        character*8 name3
        character*9 name4
        integer :: i,j,k,l,m,n

        name1 = 'fort.x'
        name2 = 'fort.xx'
        name3 = 'fort.xxx'
        name4 = 'fort.xxxx'

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))

!        open(nfile,status="old",form="unformatted")
        if(nfile < 10) then
            name1(6:6) = char(nfile+48)
            open(9,file=name1,status="unknown",form="unformatted")
        else if((nfile.ge.10).and.(nfile.lt.100)) then
            i = nfile/10
            j = nfile - 10*i
            name2(6:6) = char(i+48)
            name2(7:7) = char(j+48)
            open(9,file=name2,status="unknown",form="unformatted")
        else if((nfile.ge.100).and.(nfile.lt.1000)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name3(6:6) = char(i+48)
            name3(7:7) = char(k+48)
            name3(8:8) = char(l+48)
            open(9,file=name3,status="unknown",form="unformatted")
        else
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(6:6) = char(i+48)
            name4(7:7) = char(k+48)
            name4(8:8) = char(m+48)
            name4(9:9) = char(n+48)
            open(9,file=name4,status="unknown",form="unformatted")
        endif
 
!        print*,"nfile: ",nfile
        read(9)z,dtless
        read(9)inb,jstp,iout,itsz,isteer,islout
        read(9)msize(1:3)
        read(9)range(1:6)
!        print*,"zz: ",z,inb,jstp
        read(9)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        read(9)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        read(9)this%Nptlocal
        read(9)this%refptcl(1:6)
        allocate(this%Pts1(6,this%Nptlocal))
        read(9)this%Pts1(1:6,1:this%Nptlocal)

        close(9)

        geom%Meshsize = msize
        geom%SpatRange = range
!        allocate(geom%LcTabrg(4,0:nprocrow-1,0:nproccol-1))
!        allocate(geom%LcTabnm(2,0:nprocrow-1,0:nproccol-1))
        geom%lcTabnm = Localnum
        geom%LcTabrg = Localrange

        deallocate(Localnum)
        deallocate(Localrange)

        geom%Meshnum(1) = nx
        geom%Meshnum(2) = ny
        geom%Meshnum(3) = nz

        geom%Mshlocal(3) = geom%LcTabnm(1,myidx,myidy)
        geom%Mshlocal(2) = geom%LcTabnm(2,myidx,myidy)
        geom%Mshlocal(1) = geom%Meshnum(1)

        geom%Sptrnglocal(5) = geom%LcTabrg(1,myidx,myidy)
        geom%Sptrnglocal(6) = geom%LcTabrg(2,myidx,myidy)
        geom%Sptrnglocal(3) = geom%LcTabrg(3,myidx,myidy)
        geom%Sptrnglocal(4) = geom%LcTabrg(4,myidx,myidy)
        geom%Sptrnglocal(1) = geom%SpatRange(1)
        geom%Sptrnglocal(2) = geom%SpatRange(2)

        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Npt = nptot

        end subroutine inpoint_Output

        subroutine outpoint_Output(nfile,this,z,inb,jstp,nprocrow,nproccol,&
                   geom,iout,itsz,isteer,islout,dtless)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: z,dtless
        type(CompDom), intent(in) :: geom
        integer, intent(in) :: inb,jstp,nprocrow,nproccol,&
                               iout,itsz,isteer,islout
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,Localnum)

        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))
        call getlctabrg_CompDom(geom,Localrange)
  
        call getmsize_CompDom(geom,msize)
        call getrange_CompDom(geom,range)

        open(nfile,status="unknown",form="unformatted")

        write(nfile)z,dtless
        write(nfile)inb,jstp,iout,itsz,isteer,islout
        write(nfile)msize(1:3)
        write(nfile)range(1:6)
        write(nfile)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        write(nfile)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        deallocate(Localnum)
        deallocate(Localrange)

        write(nfile)this%Nptlocal
        write(nfile)this%refptcl(1:6)
        write(nfile)this%Pts1(1:6,1:this%Nptlocal)

        close(nfile)

        end subroutine outpoint_Output

      !--------------------------------------------------------------------------------------
      !> @author Ji Qiang
      !> @date November 6, 2008
      !> @brief
      !> This program calculate the current profile, slice emittance
      !> and uncorrelated energy spread using linear deposition.
      !--------------------------------------------------------------------------------------
      subroutine sliceprocdep_Output(pts,innp,npt,nslice,qchg,pmass,nfile)
      implicit none
      include 'mpif.h'
      integer :: nfile,innp,nslice,npt
      real*8, pointer, dimension(:,:) :: pts
      real*8, dimension(nslice) :: epx,epy,gam,xx,xx2,&
      px,px2,yy,yy2,py,py2,count,xpx,ypy,gam2uncor,gam2uncor2
      real*8, dimension(nslice,12) :: tmparry,tmpglb
      real*8 :: qchg,pmass
      integer :: i,iz,iz1,nsend,ierr,my_rank
      real*8 :: zmin,zmax,hz,zz,zavg,detaabb,&
                deltagamma,gammaz,ab,clite,sclcur,zmingl,zmaxgl

      clite = 2.99792e8
      sclcur = clite*qchg/npt

      zmin = 1.e12
      zmax = -1.e12
      do i = 1, innp
        if(zmin.ge.pts(5,i)) then
          zmin = pts(5,i)
        endif
        if(zmax.le.pts(5,i)) then
          zmax = pts(5,i)
        endif
      enddo
      call MPI_ALLREDUCE(zmin,zmingl,1,MPI_DOUBLE_PRECISION,&
                            MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(zmax,zmaxgl,1,MPI_DOUBLE_PRECISION,&
                            MPI_MAX,MPI_COMM_WORLD,ierr)

      !avoid index overflow
      zmin = zmingl - (zmaxgl-zmingl)*0.5e-5
      zmax = zmaxgl + (zmaxgl-zmingl)*0.5e-5
      hz = (zmax-zmin)/(nslice-1)
      !print*,"zmin, zmax and hz: ",zmin*scxlt,zmax*scxlt,hz*scxlt
      count = 0.0
      xx = 0.0
      xx2 = 0.0
      px = 0.0
      px2 = 0.0
      xpx = 0.0
      yy = 0.0
      yy2 = 0.0
      py = 0.0 
      py2 = 0.0 
      ypy = 0.0
      gam = 0.0
      zavg = 0.0
      zavg = (zmin+zmax)/2
      !calculate the mean of each slice by linear deposition
      do i = 1, innp
        iz = (pts(5,i)-zmin)/hz + 1
        if(iz.lt.1) iz = 1
        if(iz.ge.nslice) iz = nslice-1
        iz1 = iz + 1
        ab = ((zmin-pts(5,i))+iz*hz)/hz
        count(iz) = count(iz) + ab
        count(iz1) = count(iz1) + 1.0d0 - ab
        xx(iz) = xx(iz) + pts(1,i)*ab
        xx(iz1) = xx(iz1) + pts(1,i)*(1.0d0-ab) 
        xx2(iz) = xx2(iz) + pts(1,i)**2*ab
        xx2(iz1) = xx2(iz1) + pts(1,i)**2*(1.0d0-ab)
        px(iz) = px(iz) + pts(2,i)*ab
        px(iz1) = px(iz1) + pts(2,i)*(1.0d0-ab)
        px2(iz) = px2(iz) + pts(2,i)**2*ab
        px2(iz1) = px2(iz1) + pts(2,i)**2*(1.0d0-ab)
        xpx(iz) = xpx(iz) + pts(1,i)*pts(2,i)*ab
        xpx(iz1) = xpx(iz1) + pts(1,i)*pts(2,i)*(1.0d0-ab)
        yy(iz) = yy(iz) + pts(3,i)*ab
        yy(iz1) = yy(iz1) + pts(3,i)*(1.0d0-ab)
        yy2(iz) = yy2(iz) + pts(3,i)**2*ab
        yy2(iz1) = yy2(iz1) + pts(3,i)**2*(1.0d0-ab)
        py(iz) = py(iz) + pts(4,i)*ab
        py(iz1) = py(iz1) + pts(4,i)*(1.0d0-ab)
        py2(iz) = py2(iz) + pts(4,i)**2*ab
        py2(iz1) = py2(iz1) + pts(4,i)**2*(1.0d0-ab)
        ypy(iz) = ypy(iz) + pts(3,i)*pts(4,i)*ab
        ypy(iz1) = ypy(iz1) + pts(3,i)*pts(4,i)*(1.0d0-ab)
        gam(iz) = gam(iz) + sqrt(pts(2,i)**2+pts(4,i)**2+pts(6,i)**2+1)*ab
        gam(iz1) = gam(iz1) + &
                   sqrt(pts(2,i)**2+pts(4,i)**2+pts(6,i)**2+1)*(1.0d0-ab)
      enddo 
      tmparry(:,1) = count
      tmparry(:,2) = xx
      tmparry(:,3) = xx2
      tmparry(:,4) = px
      tmparry(:,5) = px2
      tmparry(:,6) = xpx
      tmparry(:,7) = yy
      tmparry(:,8) = yy2
      tmparry(:,9) = py
      tmparry(:,10) = py2
      tmparry(:,11) = ypy
      tmparry(:,12) = gam

      nsend = nslice*12

      tmpglb = 0.0d0
      call MPI_ALLREDUCE(tmparry,tmpglb,nsend,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)
      count = tmpglb(:,1)
      xx = tmpglb(:,2)
      xx2 = tmpglb(:,3)
      px = tmpglb(:,4)
      px2 = tmpglb(:,5)
      xpx = tmpglb(:,6)
      yy = tmpglb(:,7)
      yy2 = tmpglb(:,8)
      py = tmpglb(:,9)
      py2 = tmpglb(:,10)
      ypy = tmpglb(:,11)
      gam = tmpglb(:,12)

      !print*,"sum count: ",sum(count)
      do i = 1, nslice
        if(count(i).gt.0) then
          xx(i) = xx(i)/count(i)
          xx2(i) = xx2(i)/count(i)
          px(i) = px(i)/count(i)
          px2(i) = px2(i)/count(i)
          xpx(i) = xpx(i)/count(i)
          yy(i) = yy(i)/count(i)
          yy2(i) = yy2(i)/count(i)
          py(i) = py(i)/count(i)
          py2(i) = py2(i)/count(i)
          ypy(i) = ypy(i)/count(i)
          gam(i) = gam(i)/count(i)
        else
          xx(i) = 0.0d0
          xx2(i) = 0.0d0
          px(i) = 0.0d0
          px2(i) = 0.0d0
          yy(i) = 0.0d0
          yy2(i) = 0.0d0
          py(i) = 0.0d0
          py2(i) = 0.0d0
          xpx(i) = 0.0d0
          ypy(i) = 0.0d0
          gam(i) = 0.0d0
        endif
      enddo

      do i = 1, nslice
        epx(i) = sqrt((xx2(i)-xx(i)**2)*(px2(i)-px(i)**2)-&
                      (xpx(i)-xx(i)*px(i))**2)
        epy(i) = sqrt((yy2(i)-yy(i)**2)*(py2(i)-py(i)**2)-&
                      (ypy(i)-yy(i)*py(i))**2)
      enddo

      gam2uncor = 0.0
      do i = 1, innp
        iz = (pts(5,i)-zmin)/hz + 1
        if(iz.lt.1) iz = 1
        if(iz.ge.nslice) iz = nslice-1
        iz1 = iz + 1
        ab = ((zmin-pts(5,i))+iz*hz)/hz
        gammaz = gam(iz)*ab+ gam(iz1)*(1.0d0-ab)
        deltagamma = sqrt(pts(2,i)**2+pts(4,i)**2+pts(6,i)**2+1.0)-gammaz
        !write(11,111)pts(5,i),deltagamma,gammaz
        gam2uncor(iz) = gam2uncor(iz) + deltagamma**2*ab
        gam2uncor(iz1) = gam2uncor(iz1) + deltagamma**2*(1.0d0-ab)
      enddo
!111   format(3(1x,e15.7))
      gam2uncor2 = 0.0d0
      call MPI_ALLREDUCE(gam2uncor,gam2uncor2,nslice,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

      do i = 1, nslice
        gam2uncor2(i) = sqrt(gam2uncor2(i)/count(i))
      enddo

      call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

      if(my_rank.eq.0) then
        do i = 1, nslice
          zz = zmin + (i-1)*hz  - zavg
          write(nfile,777)zz*scxlt,count(i),count(i)/(hz*scxlt)*sclcur,epx(i)*scxlt,&
                  epy(i)*scxlt,gam(i)*pmass,gam2uncor2(i)*pmass
        enddo
      endif

      call flush(nfile)

777   format(7(1x,e15.7))

      end subroutine sliceprocdep_Output

        !> for 2 bunch diagnostic, only work on 1 PE.
        !> calculate averaged <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        !> <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance from
        !> multiple bunch/bin at fixed z (i.e. bunch center).
        subroutine diagnostic1avgZtest_Output(z,this,Nbunch)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), dimension(:), intent(inout) :: this
        integer, intent(in) :: Nbunch
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(29) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax,z0gl,z0avg,testmax
        double precision :: gamlc,gam2lc,gam2avg,tmpgam,gamdel
        double precision :: tmpx,tmpy,deltat,recpgamma
        integer :: npctmin,npctmax,ib,innpmb,i1,i2

        call starttime_Timer(t0)

        qmc = this(1)%Mass/1.0e6
        xl = Scxlt
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        do ib = 1, Nbunch

          nptot = 0
          innpmb = 0
          z0lc = 0.0
          innp = this(ib)%Nptlocal
          innpmb = innpmb + innp
          !nptot = nptot + this(ib)%Npt
          nptot = this(ib)%Npt
          do i = 1, innp
            z0lc = z0lc + this(ib)%Pts1(5,i)
          enddo
          !enddo

        call MPI_ALLREDUCE(z0lc,z0gl,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)

        den1 = 1.0/dble(nptot)
        den2 = den1*den1
        z0avg = z0gl*den1

        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0

        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            !localmax(i) = abs(this(1)%Pts1(i,1))
            localmax(i) = 0.0
          enddo
          !lcrmax = this(1)%Pts1(1,1)**2+this(1)%Pts1(3,1)**2
          lcrmax = 0.0
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
!        testmax = 0.0
        gamlc = 0.0
        gam2lc = 0.0
!        do ib = 1, Nbunch

          innp = this(ib)%Nptlocal
          do i = 1, innp
            tmpgam = sqrt(1.0+this(ib)%Pts1(2,i)**2+&
               this(ib)%Pts1(4,i)**2+this(ib)%Pts1(6,i)**2)
            gamlc = gamlc + tmpgam  
            gam2lc = gam2lc + tmpgam**2
            recpgamma = 1.0/tmpgam
            deltat = (this(ib)%Pts1(5,i)-z0avg)&
                     /(recpgamma*this(ib)%Pts1(6,i))
            !drift to the "fixed z" location
            tmpx = this(ib)%Pts1(1,i)-recpgamma*this(ib)%Pts1(2,i)*deltat
            x0lc = x0lc + tmpx
            sqsum1local = sqsum1local + tmpx*tmpx
            x0lc3 = x0lc3 + tmpx*tmpx*tmpx
            x0lc4 = x0lc4 + tmpx*tmpx*tmpx*tmpx
!            xpxlocal = xpxlocal + tmpx*this(ib)%Pts1(2,i)/this(ib)%Pts1(6,i)
            xpxlocal = xpxlocal + tmpx*this(ib)%Pts1(2,i)
!            x0lc = x0lc + this(ib)%Pts1(1,i)
!            sqsum1local = sqsum1local + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
!            x0lc3 = x0lc3 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
!            x0lc4 = x0lc4 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*&
!                    this(ib)%Pts1(1,i)
!            xpxlocal = xpxlocal + this(ib)%Pts1(1,i)*this(ib)%Pts1(2,i)
            px0lc = px0lc + this(ib)%Pts1(2,i)
!            px0lc = px0lc + this(ib)%Pts1(2,i)/this(ib)%Pts1(6,i)
            sqsum2local = sqsum2local + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
!            sqsum2local = sqsum2local + (this(ib)%Pts1(2,i)/this(ib)%Pts1(6,i))**2
            px0lc3 = px0lc3 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
            px0lc4 = px0lc4 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*&
                     this(ib)%Pts1(2,i)
            !drift to the "fixed z" location
            tmpy = this(ib)%Pts1(3,i)-recpgamma*this(ib)%Pts1(4,i)*deltat
            y0lc = y0lc + tmpy
            sqsum3local = sqsum3local + tmpy*tmpy
            y0lc3 = y0lc3 + tmpy*tmpy*tmpy
            y0lc4 = y0lc4 + tmpy*tmpy*tmpy*tmpy
            ypylocal = ypylocal + tmpy*this(ib)%Pts1(4,i)
!            y0lc = y0lc + this(ib)%Pts1(3,i)
!            sqsum3local = sqsum3local + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
!            y0lc3 = y0lc3 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
!            y0lc4 = y0lc4 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*&
!                    this(ib)%Pts1(3,i)
!            ypylocal = ypylocal + this(ib)%Pts1(3,i)*this(ib)%Pts1(4,i)
            py0lc = py0lc + this(ib)%Pts1(4,i)
            sqsum4local = sqsum4local + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc3 = py0lc3 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc4 = py0lc4 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*&
                     this(ib)%Pts1(4,i)
            sqsum5local = sqsum5local + (this(ib)%Pts1(5,i)-z0avg)**2
            z0lc3 = z0lc3 + abs((this(ib)%Pts1(5,i)-z0avg)**3)
            z0lc4 = z0lc4 + (this(ib)%Pts1(5,i)-z0avg)**4
            zpzlocal = zpzlocal + (this(ib)%Pts1(5,i)-z0avg)*this(ib)%Pts1(6,i)
            pz0lc = pz0lc + this(ib)%Pts1(6,i)
            sqsum6local = sqsum6local + this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)
            pz0lc3 = pz0lc3 + this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)
            pz0lc4 = pz0lc4 + this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)*this(ib)%Pts1(6,i)*&
                              this(ib)%Pts1(6,i)
            do j = 1, 4
              if(localmax(j).lt.abs(this(ib)%Pts1(j,i))) then
                 localmax(j) = abs(this(ib)%Pts1(j,i))
              endif
            enddo
            if(localmax(5).lt.abs(this(ib)%Pts1(5,i)-z0avg)) then
               localmax(5) = abs(this(ib)%Pts1(5,i)-z0avg)
               i1 = i
            endif
!            if(testmax .gt. this(ib)%Pts1(5,i) ) then
!              testmax = this(ib)%Pts1(5,i)
!              i2 = i
!            endif
            if(localmax(6).lt.abs(this(ib)%Pts1(6,i))) then
                 localmax(6) = abs(this(ib)%Pts1(6,i))
            endif
            if(lcrmax.lt.(this(ib)%Pts1(1,i)**2+this(ib)%Pts1(3,i)**2)) then
              lcrmax = this(ib)%Pts1(1,i)**2 + this(ib)%Pts1(3,i)**2
            endif
          enddo
!        enddo

!        print*,"z0avg: ",z0avg*xl,localmax(5)*xl,testmax*xl,xl,&
!               this(1)%Pts1(5,i1)*xl,this(1)%Pts1(5,i2)*xl
        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        tmplc(28) = gamlc
        tmplc(29) = gam2lc
        
        call MPI_REDUCE(tmplc,tmpgl,29,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
!          sqsum5 = sqz - z0*z0
          sqsum5 = sqz 
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
!          zpz = tmpgl(15)*den1 - z0*pz0
          zpz = tmpgl(15)*den1 
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
!          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z03 = cubz**(1.0/3.0)
!          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          z04 = sqrt(sqrt(fthz))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(abs(sqsum1))
          pxrms = sqrt(abs(sqsum2))
          yrms = sqrt(abs(sqsum3))
          pyrms = sqrt(abs(sqsum4))
          zrms = sqrt(abs(sqsum5))
          pzrms = sqrt(abs(sqsum6))
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          gam = tmpgl(28)*den1
!          gam = sqrt(1.0+px0**2+py0**2+pz0**2)
          energy = qmc*(gam-1.0)
          bet = sqrt(1.0-(1.0/gam)**2)
          gam2avg = tmpgl(29)*den1
          gamdel = sqrt(abs(gam2avg - gam**2))

          if(ib.eq.1) then
            write(18,100)z,z0avg*xl,gam,energy,bet,sqrt(glrmax)*xl,gamdel
            write(24,102)z,z0*xl,x0*xl,xrms*xl,px0,pxrms,-xpx*xl,epx*xl
            write(25,102)z,z0*xl,y0*xl,yrms*xl,py0,pyrms,-ypy*xl,epy*xl
            write(26,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz*xl,epz*xl
            write(27,102)z,z0*xl,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
            write(28,101)z,z0*xl,npctmin,npctmax,nptot
          else
            write(188,100)z,z0avg*xl,gam,energy,bet,sqrt(glrmax)*xl,gamdel
            write(244,102)z,z0*xl,x0*xl,xrms*xl,px0,pxrms,-xpx*xl,epx*xl
            write(255,102)z,z0*xl,y0*xl,yrms*xl,py0,pyrms,-ypy*xl,epy*xl
            write(266,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz*xl,epz*xl
            write(277,102)z,z0*xl,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
            write(288,101)z,z0*xl,npctmin,npctmax,nptot
          endif

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
        endif

        enddo

99      format(6(1x,e16.8))
100      format(7(1x,e16.8))
101     format(1x,e16.8,e16.8,3I10)
102      format(8(1x,e16.8))

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1avgZtest_Output
      end module Outputclass

