!----------------------------------------------------------------
! (c) Copyright, 2012 by the Regents of the University of California.
! Version: 1.0
! Author: Ji Qiang
!prepare the 4 input files for a chicane used in the Impact-T code
!The 4 output files are fort.1, fort.2, fort.3, and fort.4
!
        program chicaneImpt
        implicit double precision (a-h,o-z)
        real*8, dimension(4) :: xk1,dz1,xk2,dz2,xk3,dz3,xk4,dz4
        real*8 :: gambet

        print*,"input gamma, B field, and rBend length:"
        read(*,*)gamma,bfld,dleng

        print*,"first dipole bending upward (1) or downward (-1)?"
        iflagbend = -1
        read(*,*)iflagbend
        if(iflagbend.ne.1 .and. iflagbend.ne.-1) then
          print*,"wrong input! (has to be 1 or -1)"
          stop
        endif

        emass = 0.511005d6
        clite = 2.99792458d8
        gambet = sqrt(gamma**2-1.0d0)
        rho = emass*gambet/(clite*abs(bfld))
        !bending angle
        theta = asin(dleng/rho)
        !arc length
        arcleng = rho*theta
        print*,"bending angle:",theta
        print*,"bending radius:",rho
        print*,"arc length:",arcleng
 
        theta = theta*iflagbend
        print*,"entrance fringe field length, flat field lengh, and &
                &exit fringe field length:"
        dl1 = 0.21d0
        dl2 = 0.0d0
        dl3 = 0.21d0
        read(*,*)dl1,dl2,dl3
        !1st rectanuglar bend
        xk1(1) = 0.0d0
        dz1(1) = 0.0d0
        xk2(1) = 0.0d0
        dz2(1) = dz1(1)+dl1
        xk3(1) = 0.0d0
        dz3(1) = dz2(1)+dl2
        xk4(1) = 0.0d0
        dz4(1) = dz3(1)+dl3
        !2nd rectanuglar bend
        xk1(2) = tan(theta)
        dz1(2) = 0.0d0
        xk2(2) = tan(theta)
        dz2(2) = dz1(2)+dl1/cos(theta)
        xk3(2) = tan(theta)
        dz3(2) = dz2(2)+dl2/cos(theta)
        xk4(2) = tan(theta)
        dz4(2) = dz3(2)+dl3/cos(theta)
        !3rd rectanuglar bend
        xk1(3) = 0.0d0
        dz1(3) = 0.0d0
        xk2(3) = 0.0d0
        dz2(3) = dz1(3)+dl1
        xk3(3) = 0.0d0
        dz3(3) = dz2(3)+dl2
        xk4(3) = 0.0d0
        dz4(3) = dz3(3)+dl3
        !4th rectanuglar bend
        xk1(4) = -tan(theta)
        dz1(4) = 0.0d0
        xk2(4) = -tan(theta)
        dz2(4) = dz1(4)+dl1/cos(theta)
        xk3(4) = -tan(theta)
        dz3(4) = dz2(4)+dl2/cos(theta)
        xk4(4) = -tan(theta)
        dz4(4) = dz3(4)+dl3/cos(theta)
        print*,"input the shift z01 and z02 in Enge function:"
        z01 = 0.21d0
        z02 = 0.21d0
        read(*,*)z01,z02
        print*,"input the Enge function coefficients, c1-c8:"
        c1 = -0.714651d0
        c2 = 3.37476d0
        c3 = -1.58194d0
        c4 = 0.447227d0
        c5 = 0.0d0
        c6 = 0.0d0
        c7 = 0.0d0
        c8 = 0.0d0
        read(*,*)c1,c2,c3,c4,c5,c6,c7,c8
        print*,"input CSR flag (0 - no CSR, 1 - CSR):"
        flagcsr = 0.0d0
        read(*,*)flagcsr
        print*,"input the effective starting and ending location of &
               &the bend (for CSR only):"
        zcsr1 = 0.1d0
        zcsr2 = 0.5d0
        read(*,*)zcsr1,zcsr2
        !write output files 1-4: those files should be renamed as rfdata...
        !in the Impact-T simulation
        write(1,*)flagcsr
        write(1,*)gamma
        write(1,*)xk1(1)
        write(1,*)dz1(1)
        write(1,*)xk2(1)
        write(1,*)dz2(1)
        write(1,*)xk3(1)
        write(1,*)dz3(1)
        write(1,*)xk4(1)
        write(1,*)dz4(1)
        write(1,*)z01
        write(1,*)z02
        write(1,*)c1
        write(1,*)c2
        write(1,*)c3
        write(1,*)c4
        write(1,*)c5
        write(1,*)c6
        write(1,*)c7
        write(1,*)c8
        write(1,*)zcsr1
        write(1,*)zcsr2
        close(1)
        write(2,*)flagcsr
        write(2,*)gamma
        write(2,*)xk1(2)
        write(2,*)dz1(2)
        write(2,*)xk2(2)
        write(2,*)dz2(2)
        write(2,*)xk3(2)
        write(2,*)dz3(2)
        write(2,*)xk4(2)
        write(2,*)dz4(2)
        write(2,*)z01
        write(2,*)z02
        write(2,*)c1
        write(2,*)c2
        write(2,*)c3
        write(2,*)c4
        write(2,*)c5
        write(2,*)c6
        write(2,*)c7
        write(2,*)c8
        write(2,*)zcsr1
        write(2,*)zcsr2
        close(2)
        write(3,*)flagcsr
        write(3,*)gamma
        write(3,*)xk1(3)
        write(3,*)dz1(3)
        write(3,*)xk2(3)
        write(3,*)dz2(3)
        write(3,*)xk3(3)
        write(3,*)dz3(3)
        write(3,*)xk4(3)
        write(3,*)dz4(3)
        write(3,*)z01
        write(3,*)z02
        write(3,*)c1
        write(3,*)c2
        write(3,*)c3
        write(3,*)c4
        write(3,*)c5
        write(3,*)c6
        write(3,*)c7
        write(3,*)c8
        write(3,*)zcsr1
        write(3,*)zcsr2
        close(3)
        write(4,*)flagcsr
        write(4,*)gamma
        write(4,*)xk1(4)
        write(4,*)dz1(4)
        write(4,*)xk2(4)
        write(4,*)dz2(4)
        write(4,*)xk3(4)
        write(4,*)dz3(4)
        write(4,*)xk4(4)
        write(4,*)dz4(4)
        write(4,*)z01
        write(4,*)z02
        write(4,*)c1
        write(4,*)c2
        write(4,*)c3
        write(4,*)c4
        write(4,*)c5
        write(4,*)c6
        write(4,*)c7
        write(4,*)c8
        write(4,*)zcsr1
        write(4,*)zcsr2
        close(4)

        end program chicaneImpt
