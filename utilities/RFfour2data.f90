!------------------------------------------------------------------------------
! (c) Copyright, 2024 by the Regents of the University of California.
! Version: 1.2
! Author: Ji Qiang
! Description: Given the Fourier coefficients of RF field on axis and produce a
!              field data and its derivatives on axis. 
! Input: rfdata.in:   
! Output: rfdata.out: reconstruct the field from the Fourier coefficients.
!         rfdatax: is the one to be used for the IMPACT-T
!         simulation.
!         rfdata.tmp: shifted Ez vs. z. This file can be used as rfdata.in to
!                     generate the Fourier coefficients of shifted RF data. This
!                     coefficients together with the original rf data Fourier 
!                     coefficients will be used to model the fields in TWS.
      program rfcoef
      implicit double precision (a-h,o-z)
      integer, parameter :: Ndata = 50000
      integer, parameter :: Ncoef = 200
      double precision, dimension(Ndata) :: zdata,edata
      double precision, dimension(Ncoef) :: Fcoef,Fcoef2

      emax = 0.0d0
      open(3,file="rfdata.in",status="old")
      read(3,*)ntmp
      ncoefreal = (ntmp-1)/2 + 1
      print*, "ncoefreal: ",ncoefreal
      read(3,*)zst
      read(3,*)zend
      read(3,*)tmp1
      read(3,*)Fcoef(1)
      do i = 2, ncoefreal
        read(3,*)Fcoef(i)
        read(3,*)Fcoef2(i)
      enddo
      close(3)
     
      zdata(1) = zst
      zdata(ndatareal) = zend

      print*,"How many data points you want?"
      read(*,*)ndatareal

      print*,"zdata1: ",zdata(1),zdata(ndatareal)
      zlen = zdata(ndatareal)-zdata(1)
      zhalf = zlen/2.0
      zmid = (zdata(ndatareal)+zdata(1))/2
      h = zlen/(ndatareal-1)
      pi = 2*asin(1.0)
      print*,"The RF data number is: ",ndatareal,zlen,zmid,h

      open(8,file="rfdata.out",status="unknown")
      do i = 1, ndatareal
        !zz = zdata(i) - zmid
        zz = (i-1)*h+zdata(1) - zmid
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
        !write(8,*)zz,tmpsum,tmpsump,tmpsumpp
        write(8,111)zz,tmpsum,tmpsump,tmpsumpp,tmpsump3
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

111   format(5(1x,e17.9))

      dd = 0.0349753333333333
      print*,"input shift length:"
      read(*,*)dd
      open(9,file="rfdata.tmp",status="unknown")
      do i = 1, ndatareal
        zz = zdata(i) - zmid + dd
        tmpsum = 0.5*Fcoef(1)
        do j = 2,ncoefreal
         tmpsum = tmpsum + Fcoef(j)*cos((j-1)*2*pi*zz/zlen) + &
                  Fcoef2(j)*sin((j-1)*2*pi*zz/zlen)
        enddo
        write(9,*)zdata(i),tmpsum
      enddo
      close(9)

      stop
      end program rfcoef
