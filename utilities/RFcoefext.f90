!----------------------------------------------------------------
! (c) Copyright, 2006 by the Regents of the University of California.
! Version: 1.2
! Author: Ji Qiang
! Description: find the Fourier coefficients of RF data on axis.
! The input data has been extended to make a symmetric data for Fourier expansion
! Input: rfdata.in: Ez vs. z (m)
! Output: rfdata.out: reconstruct the field from the Fourier coefficients.
!         rfdatax: is the one to be used for the IMPACT-T
!         simulation.
      program rfcoef
      integer, parameter :: Ndata = 50000
      integer, parameter :: Ncoef = 200
      double precision, dimension(Ndata) :: zdata,edata
      double precision, dimension(Ncoef) :: Fcoef,Fcoef2

      print*,"Note: in rfdata.in file, z must be in unit of meter!!!!"
      open(3,file="rfdata.in",status="old")
      n = 0
10    continue
!        read(3,*,end=100)tmp1,tmp2,tmp3,tmp4
        read(3,*,end=100)tmp1,tmp2
        n = n+1
        !zdata(n) = tmp1/1000.
        !edata(n) = tmp2/10000.
        zdata(n) = tmp1
        edata(n) = tmp2
      goto 10
100   continue
      close(3)

      emax = 0.0
      do i = 1, n
!        edata(i) = -(edata(i)-edata(n))
        if(abs(emax).le.abs(edata(i))) then
          emax = edata(i)
        endif
!        zdata(i) = zdata(i)/100
      enddo

      do i = 1, n
        edata(i) = edata(i)/emax
        write(44,*)zdata(i),edata(i)
      enddo
      close(44)
      
      zst = zdata(1)

      zhalf = zdata(n)-zst
      print*,"zhalf: ",zhalf
      do i = n+1, 2*n
        zdata(i) = zdata(i-n)-zst + zhalf
        edata(i) = edata(i-n)
      enddo
      do i = 1, n
        zdata(i) = -(zdata(2*n+1-i)) + 2*zhalf
        edata(i) = edata(2*n+1-i)
      enddo
      ndatareal = 2*n

      print*,"How many Fourier coeficients you want?"
      read(*,*)ncoefreal

      zlen = zdata(ndatareal)-zdata(1)
      zhalf = zlen/2.0d0
      zmid = (zdata(ndatareal)+zdata(1))/2
      h = zlen/(ndatareal-1)
      pi = 2*asin(1.0d0)
      print*,"The RF data number is: ",ndatareal,zlen,zmid,h

      do j = 1, ncoefreal
        zz = zdata(1) - zmid
        Fcoef(j) = (-0.5*edata(1)*cos((j-1)*2*pi*zz/zlen)*h)/zhalf
        Fcoef2(j) = (-0.5*edata(1)*sin((j-1)*2*pi*zz/zlen)*h)/zhalf
        zz = zdata(ndatareal) - zmid
        Fcoef(j) = Fcoef(j)-(0.5*edata(ndatareal)*cos((j-1)*2*pi*zz/zlen)*h)&
                            /zhalf
        Fcoef2(j) = Fcoef2(j)-(0.5*edata(ndatareal)*sin((j-1)*2*pi*zz/zlen)*h)&
                            /zhalf
      enddo

      do i = 1, ndatareal
        zz = (i-1)*h+zdata(1)
        klo=1
        khi=ndatareal
1       if(khi-klo.gt.1) then
          k=(khi+klo)/2
          if(zdata(k).gt.zz)then
             khi=k
          else
             klo=k
          endif
          goto 1
        endif
        hstep=zdata(khi)-zdata(klo)
        slope=(edata(khi)-edata(klo))/hstep
        ez1 =edata(klo)+slope*(zz-zdata(klo))

        zz = zdata(1)+(i-1)*h - zmid
        do j = 1, ncoefreal
          Fcoef(j) = Fcoef(j) + (ez1*cos((j-1)*2*pi*zz/zlen)*h)/zhalf
          Fcoef2(j) = Fcoef2(j) + (ez1*sin((j-1)*2*pi*zz/zlen)*h)/zhalf
        enddo
      enddo

      open(7,file="rfcoef.out",status="unknown")
      do j = 1, ncoefreal
        write(7,*)j,Fcoef(j),Fcoef2(j)
      enddo
      close(7)

      open(8,file="rfdatax",status="unknown")
      write(8,*)Fcoef(1)
      do j = 2, ncoefreal
        write(8,*)Fcoef(j)
        write(8,*)Fcoef2(j)
      enddo
      close(8)

      open(8,file="rfdata.out",status="unknown")
      do i = 1, ndatareal
        zz = zdata(i) - zmid
        tmpsum = 0.5*Fcoef(1)
        tmpsump = 0.0
        tmpsumpp = 0.0
        do j = 2,ncoefreal
         tmpsum = tmpsum + Fcoef(j)*cos((j-1)*2*pi*zz/zlen) + &
                  Fcoef2(j)*sin((j-1)*2*pi*zz/zlen)
         tmpsump = tmpsump-(j-1)*2*pi*Fcoef(j)*sin((j-1)*2*pi*zz/zlen)/zlen +&
                  (j-1)*2*pi*Fcoef2(j)*cos((j-1)*2*pi*zz/zlen)/zlen
         tmpsumpp = tmpsumpp-((j-1)*2*pi/zlen)**2*&
                  (Fcoef(j)*cos((j-1)*2*pi*zz/zlen) + &
                   Fcoef2(j)*sin((j-1)*2*pi*zz/zlen))
        enddo
        write(8,*)zdata(i),tmpsum,tmpsump,tmpsumpp
      enddo
      close(8)

      print*,"input # of data points:"
      read(*,*)nout
      hz = (zdata(ndatareal)-zdata(1))/(nout-1)
      zmid = (zdata(ndatareal)+zdata(1))/2
      vtmp = 0.0

      open(8,file="rfdatax2",status="unknown")

      write(8,77)nout,zdata(1),zdata(ndatareal),vtmp
77    format(I10,3(1x,e15.7))
      do i = 1, nout
        zz = (i-1)*hz+zdata(1) - zmid
        tmpsum = 0.5*Fcoef(1)
        tmpsump = 0.0
        tmpsumpp = 0.0
        tmpsump3 = 0.0
        do j = 2,ncoefreal
         tmpsum = tmpsum + Fcoef(j)*cos((j-1)*2*pi*zz/zlen) + &
                  Fcoef2(j)*sin((j-1)*2*pi*zz/zlen)
         tmpsump = tmpsump-(j-1)*2*pi*Fcoef(j)*sin((j-1)*2*pi*zz/zlen)/zlen +&
                  (j-1)*2*pi*Fcoef2(j)*cos((j-1)*2*pi*zz/zlen)/zlen
         tmpsumpp = tmpsumpp-((j-1)*2*pi/zlen)**2*&
                    (Fcoef(j)*cos((j-1)*2*pi*zz/zlen) + &
                     Fcoef2(j)*sin((j-1)*2*pi*zz/zlen))
         tmpsump3 = tmpsump3+((j-1)*2*pi/zlen)**3*&
                    (Fcoef(j)*sin((j-1)*2*pi*zz/zlen) - &
                     Fcoef2(j)*cos((j-1)*2*pi*zz/zlen))
        enddo
        write(8,101)tmpsum,tmpsump,tmpsumpp,tmpsump3
      enddo
      close(8)

101   format(4(1x,e17.9))
      stop
      end program rfcoef
