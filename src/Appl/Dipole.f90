!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! Dipoleclass: Dipole beam line element class
!             in Lattice module of APPLICATION layer.
! MODULE  : ... Dipoleclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang
! DESCRIPTION: 
!> This class defines the linear transfer map and field for the Dipole beam line elment.
!              
! Comments:
!----------------------------------------------------------------
      module Dipoleclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 10
        type Dipole
          !Itype = 4
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : x field strength
          !      (3) : y field strength
          !      (4) : file ID
          !      (5) : radius
          !      (6) : x misalignment error
          !      (7) : y misalignment error
          !      (8) : rotation error x
          !      (9) : rotation error y
          !      (10) : rotation error z
        end type Dipole
        interface getparam_Dipole
          module procedure getparam1_Dipole,  &
                          getparam2_Dipole,   &
                          getparam3_Dipole
        end interface
        interface setparam_Dipole
          module procedure setparam1_Dipole,  &
                          setparam2_Dipole, setparam3_Dipole
        end interface
      contains
        subroutine construct_Dipole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Dipole), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Dipole
   
        subroutine setparam1_Dipole(this,i,value)
        implicit none
        type (Dipole), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Dipole

        subroutine setparam2_Dipole(this,values)
        implicit none
        type (Dipole), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_Dipole

        subroutine setparam3_Dipole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Dipole), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Dipole
   
        subroutine getparam1_Dipole(this,i,blparam) 
        implicit none 
        type (Dipole), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Dipole
  
        subroutine getparam2_Dipole(this,blparams)
        implicit none
        type (Dipole), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Dipole

        subroutine getparam3_Dipole(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Dipole), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Dipole
       
        !--------------------------------------------------------------------------------------
        !> @brief
        !> The linear map calculation for the dipole is not correct
        !--------------------------------------------------------------------------------------
        subroutine maplinear_Dipole(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (Dipole), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,escale,tfin
        integer  :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        escale = this%Param(2)
        mpstp = this%Mapstp

        y(1)=0.0
        y(2)=0.0
        y(3)=0.0
        y(4)=0.0
        y(5)=refpt(5)
        y(6)=refpt(6)
        y(7)=1.0
        y(8)=0.0
        y(9)=0.0
        y(10)=1.0
        y(11)=1.0
        y(12)=0.0
        y(13)=0.0
        y(14)=1.0
        y(15)=1.0
        y(16)=0.0
        y(17)=0.0
        y(18)=1.0
        h=tau/mpstp
        t0=t
        gi=-refpt(6)
        betai=sqrt(gi**2-1.)/gi
        ui=gi*betai
        squi=sqrt(ui)
        sqi3=sqrt(ui)**3

        dlti=0.0
        thli=0.0

        call rk6i_Dipole(h,mpstp,t0,y,18,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3

        tfin = t + tau

        dltf = 0.0
        thlf = 0.0

        xm(1,1)= y(7)*squi/squf+y(9)*squi/squf*dlti
        xm(2,1)=(y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*&
                 dlti
        xm(1,2)= y(9)/(squi*squf)
        xm(2,2)=(y(10)-y(9)*dltf)*squf/squi
        xm(3,3)= y(11)*squi/squf+y(13)*squi/squf*dlti
        xm(4,3)= &
        (y(12)-y(11)*dltf)*squi*squf+(y(14)-y(13)*dltf)*squi*squf*dlti
        xm(3,4)= y(13)/(squi*squf)
        xm(4,4)=(y(14)-y(13)*dltf)*squf/squi
        xm(5,5)= y(15)*sqi3/sqf3+y(17)*sqi3/sqf3*thli
        xm(6,5)= &
        (y(16)-y(15)*thlf)*sqi3*sqf3+(y(18)-y(17)*thlf)*sqi3*sqf3*thli
        xm(5,6)= y(17)/(sqi3*sqf3)
        xm(6,6)=(y(18)-y(17)*thlf)*sqf3/sqi3
        !print*,"xm: ",xm(1,1),xm(2,1),xm(1,2),xm(2,2),xm(3,3),xm(4,3),xm(4,4)

        end subroutine maplinear_Dipole

        subroutine rk6i_Dipole(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (Dipole), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_Dipole(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_Dipole(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_Dipole

        subroutine intfunc1_Dipole(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (Dipole), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55
        double precision :: zedge,qmcc,brho,bgrad

        zedge = this%Param(1)
        bgrad = this%Param(2)
        qmcc = Bchg/Bmass

        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        f(4) = 0.0

        ! synchronous particle:
        gamma0=-y(6)
        beta0=sqrt(gamma0**2-1.)/gamma0
        gbet=sqrt((gamma0-1.0)*(gamma0+1.0))

        f(5)=1.0/(beta0*Scxl)
        f(6)=0.0

        ! matrix elements
        brho=gbet/Clight/qmcc
        s11=bgrad/brho*Scxl
        s33=-bgrad/brho*Scxl
        s55=  0.0

        f(7)=y(8)/Scxl
        f(8)=-s11*y(7)
        f(9)=y(10)/Scxl
        f(10)=-s11*y(9)
        f(11)=y(12)/Scxl
        f(12)=-s33*y(11)
        f(13)=y(14)/Scxl
        f(14)=-s33*y(13)
        f(15)=y(16)/Scxl
        f(16)=-s55*y(15)
        f(17)=y(18)/Scxl
        f(18)=-s55*y(17)

        end subroutine intfunc1_Dipole

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getflderr_Dipole(pos,extfld,this,dx,dy,anglex,angley,&
                                     anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        type (Dipole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bx,by,zedge
        double precision, dimension(3) :: temp,tmp

        zedge = this%Param(1)
!        zz = pos(3)-zedge

!        dx = this%Param(6)
!        dy = this%Param(7)
!        anglex = this%Param(8)
!        angley = this%Param(9)
!        anglez = this%Param(10)

        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = pos(3) - zedge
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = this%Param(2)
        extfld(5) = this%Param(3)
        extfld(6) = 0.0

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getflderr_Dipole
        
        !get external field without displacement and rotation errors.
        subroutine  getfld_Dipole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Dipole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bx,by,zedge

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = this%Param(2)
        extfld(5) = this%Param(3)
        extfld(6) = 0.0

        end subroutine getfld_Dipole

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field without displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getfldtlinear_Dipole(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Dipole), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,zedge,z1,z2,z3,z4,xx1,zz1,ss,k1,lengfrg,alpha
        integer :: faceid
 
        !zz = pos(3) - this%Param(1)
        zz = pos(3)
        !print*,"inside Dipole: ",pos(1),pos(2),pos(3),Fcoef(1),&
        !          Fcoef(2),Fcoef(3),Fcoef(4),Fcoef(5),Fcoef(6),&
        !          zz
        !print*,"inside Dipole: ",pos,fldata%Fcoeft(1:10)
 
        faceid = fldata%Fcoeft(1)+0.01 !switch id of pole face
        !The pole face is characterized by z = a x + b
        !Fcoef(2) is the gamma of the reference particle
        !Fcoef(3)-(10) are the coeficients for the geometry of bend
        extfld = 0.0
        z1 = fldata%Fcoeft(3)*pos(1) + fldata%Fcoeft(4)
        z2 = fldata%Fcoeft(5)*pos(1) + fldata%Fcoeft(6)
        z3 = fldata%Fcoeft(7)*pos(1) + fldata%Fcoeft(8)
        z4 = fldata%Fcoeft(9)*pos(1) + fldata%Fcoeft(10)
        if((zz.ge.z1).and.(zz.le.z2)) then
        !if((zz.gt.z1).and.(zz.lt.z2)) then
          !inside fringe field region of entrance
          k1 = fldata%Fcoeft(3)
          !this angle may have opposite sign as the conventional used one.
          alpha = atan(k1)
          xx1 = (pos(1)+pos(3)*k1-fldata%Fcoeft(4)*k1)/(1.0+k1*k1)
          zz1 = k1*xx1+fldata%Fcoeft(4)
          ss = sqrt((xx1-pos(1))**2+(zz1-pos(3))**2)
          lengfrg = (fldata%Fcoeft(6)-fldata%Fcoeft(4))*cos(alpha)
          extfld(4) = -this%Param(3)*sin(alpha)*pos(2)/lengfrg
          extfld(5) = this%Param(3)*ss/lengfrg
          extfld(6) = this%Param(3)*cos(alpha)*pos(2)/lengfrg
        else if((zz.gt.z2).and.(zz.lt.z3)) then
        !else if((zz.ge.z2).and.(zz.le.z3)) then
          !inside the constant field region
          extfld(5) = this%Param(3)
        else if((zz.ge.z3).and.(zz.le.z4)) then
        !else if((zz.gt.z3).and.(zz.lt.z4)) then
          !inside fringe field region of exit
          k1 = fldata%Fcoeft(7)
          !this angle may have opposite sign as the conventional used one.
          alpha = atan(k1)
          xx1 = (pos(1)+pos(3)*k1-fldata%Fcoeft(8)*k1)/(1.0+k1*k1)
          zz1 = k1*xx1+fldata%Fcoeft(8)
          ss = sqrt((xx1-pos(1))**2+(zz1-pos(3))**2)
          lengfrg = (fldata%Fcoeft(10)-fldata%Fcoeft(8))*cos(alpha)
          extfld(4) = this%Param(3)*sin(alpha)*pos(2)/lengfrg
          extfld(5) = this%Param(3)-this%Param(3)*ss/lengfrg
          extfld(6) = -this%Param(3)*cos(alpha)*pos(2)/lengfrg
        endif
 
        !print*,"dipole: ",pos,extfld(5)
        end subroutine getfldtlinear_Dipole

        subroutine  getfldtold_Dipole(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Dipole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        type (fielddata), intent(in) :: fldata
        double precision:: zz,bx,by,zedge,z1,z2
        integer :: faceid
 
        faceid = fldata%Fcoeft(1)+0.01 !switch id of pole face
        !Fcoeft(2) is the input gamma for the reference particle
        !Fcoeft(3)-(6) are the coeficients for the geometry of bend
        extfld = 0.0
        if(faceid.eq.1) then !parallel face
          z1 = fldata%Fcoeft(3)*pos(3) + fldata%Fcoeft(4)
          z2 = fldata%Fcoeft(5)*pos(3) + fldata%Fcoeft(6)
          if((pos(1).gt.z2).and.(pos(1).lt.z1)) then
            extfld(5) = this%Param(3)
          endif
        else if(faceid.eq.2) then !non-parallel face, bending toward x < 0
          z1 = fldata%Fcoeft(3)*pos(3) + fldata%Fcoeft(4)
          z2 = fldata%Fcoeft(5)*pos(3) + fldata%Fcoeft(6)
          if((pos(1).gt.z2).and.(pos(1).gt.z1)) then
            extfld(5) = this%Param(3)
          endif
        else if(faceid.eq.3) then !non-parallel face, bending toward x > 0
          z1 = fldata%Fcoeft(3)*pos(3) + fldata%Fcoeft(4)
          z2 = fldata%Fcoeft(5)*pos(3) + fldata%Fcoeft(6)
          if((pos(1).lt.z2).and.(pos(1).lt.z1)) then
            extfld(5) = this%Param(3)
          endif
        else
          print*,"wrong input option!!!"
          stop
        endif
 
        end subroutine getfldtold_Dipole

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external dipole field without displacement and rotation errors.
        !> here the fringe field of dipole is approximated using an Enge function.
        !--------------------------------------------------------------------------------------
        subroutine  getfldt_Dipole(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Dipole), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,zedge,z1,z2,z3,z4,xx1,zz1,ss,k1,lengfrg,alpha
        integer :: faceid
        real*8 :: dd,s1,s2,s3,s4,tmpz,tmp1,bb,bbp,bbpp,dd1,dd2
        real*8, dimension(8) :: cc
 
        !half width of the fringe field region
        dd1 = fldata%Fcoeft(11)/2
        dd2 = fldata%Fcoeft(12)/2
        !here, c1 and c2 are default values from Track manual.
        !c1 = -0.003183 
        !c2 = 2.32
        cc(1) = fldata%Fcoeft(13) 
        cc(2) = fldata%Fcoeft(14) 
        cc(3) = fldata%Fcoeft(15) 
        cc(4) = fldata%Fcoeft(16) 
        cc(5) = fldata%Fcoeft(17) 
        cc(6) = fldata%Fcoeft(18) 
        cc(7) = fldata%Fcoeft(19) 
        cc(8) = fldata%Fcoeft(20) 
        !given these c1 and c2, the fringe field region has to -2.5dd to 2.5dd to 
        !make continuous change of field.
        !zz = pos(3) - this%Param(1)
        zz = pos(3)
!        print*,"inside Dipole: ",pos(1),pos(2),pos(3),Fcoef(1),&
!                  Fcoef(2),Fcoef(3),Fcoef(4),Fcoef(5),Fcoef(6),zz
 
        faceid = fldata%Fcoeft(1)+0.01 !switch id of pole face
        !The pole face is characterized by z = a x + b
        !Fcoef(2) is the gamma of the reference particle
        !Fcoef(3)-(10) are the coeficients for the geometry of bend
        extfld = 0.0
        z1 = fldata%Fcoeft(3)*pos(1) + fldata%Fcoeft(4)
        z2 = fldata%Fcoeft(5)*pos(1) + fldata%Fcoeft(6)
        z3 = fldata%Fcoeft(7)*pos(1) + fldata%Fcoeft(8)
        z4 = fldata%Fcoeft(9)*pos(1) + fldata%Fcoeft(10)
        dd = 2*this%Param(5)
        if((zz.ge.z1).and.(zz.le.z2)) then
          !inside fringe field region of entrance
          k1 = fldata%Fcoeft(3)
          !this angle may have opposite sign as the conventional used one.
          alpha = atan(k1)
          xx1 = (pos(1)+pos(3)*k1-fldata%Fcoeft(4)*k1)/(1.0+k1*k1)
          zz1 = k1*xx1+fldata%Fcoeft(4)
          !distance from point to L1
          ss = sqrt((xx1-pos(1))**2+(zz1-pos(3))**2)
          !tmpz = -(ss-2.5*dd)/dd !shift 2.5 is due to c0 and c1 used
          tmpz = -(ss-dd1)/dd !shift 2.5 is due to c0 and c1 used
          s1 = cc(1)+cc(2)*tmpz+cc(3)*tmpz**2+cc(4)*tmpz**3+cc(5)*tmpz**4+&
               cc(6)*tmpz**5+cc(7)*tmpz**6+cc(8)*tmpz**7
          bb = 1./(1.0+exp(s1))
          !first derivative of s1
          s2 = -(cc(2)+2*cc(3)*tmpz+3*cc(4)*tmpz**2+4*cc(5)*tmpz**3+ &
                 5*cc(6)*tmpz**4+6*cc(7)*tmpz**5+7*cc(8)*tmpz**6)/dd 
          tmp1 = exp(s1)*s2
          bbp = -tmp1/(1.0+exp(s1))**2
          !second derivative of s1
          s3 = (2*cc(3)+6*cc(4)*tmpz+12*cc(5)*tmpz**2+&
                20*cc(6)*tmpz**3+30*cc(7)*tmpz**4+42*cc(8)*tmpz**5)/dd/dd 
          s4 = exp(s1)*(s3+s2*s2)
          bbpp = (-s4/(1.0+exp(s1))**2+2*tmp1*tmp1/(1.0+exp(s1))**3)

          extfld(4) = -this%Param(3)*(bbp*pos(2))*sin(alpha)
          extfld(5) = this%Param(3)*(bb-bbpp/2*pos(2)*pos(2))
          extfld(6) = this%Param(3)*(bbp*pos(2))*cos(alpha)
        else if((zz.gt.z2).and.(zz.lt.z3)) then
          !inside the constant field region
          extfld(5) = this%Param(3)
        else if((zz.ge.z3).and.(zz.le.z4)) then
          !inside fringe field region of exit
          k1 = fldata%Fcoeft(7)
          !this angle may have opposite sign as the conventional used one.
          alpha = atan(k1)
          xx1 = (pos(1)+pos(3)*k1-fldata%Fcoeft(8)*k1)/(1.0+k1*k1)
          zz1 = k1*xx1+fldata%Fcoeft(8)
          !distance from point to L3
          ss = sqrt((xx1-pos(1))**2+(zz1-pos(3))**2)
          !tmpz = (ss-2.5*dd)/dd !shift 2.5 is due to c0 and c1
          tmpz = (ss-dd2)/dd !shift 2.5 is due to c0 and c1
          s1 = cc(1)+cc(2)*tmpz+cc(3)*tmpz**2+cc(4)*tmpz**3+cc(5)*tmpz**4+&
               cc(6)*tmpz**5+cc(7)*tmpz**6+cc(8)*tmpz**7
          bb = 1.0d0/(1.0+exp(s1))
          s2 = (cc(2)+2*cc(3)*tmpz+3*cc(4)*tmpz**2+4*cc(5)*tmpz**3+ &
                 5*cc(6)*tmpz**4+6*cc(7)*tmpz**5+7*cc(8)*tmpz**6)/dd 
          tmp1 = exp(s1)*s2
          bbp = -tmp1/(1.0+exp(s1))**2
          s3 = (2*cc(3)+6*cc(4)*tmpz+12*cc(5)*tmpz**2+&
                20*cc(6)*tmpz**3+30*cc(7)*tmpz**4+42*cc(8)*tmpz**5)/dd/dd 
          s4 = exp(s1)*(s3 + s2*s2)
          bbpp = (-s4/(1.0+exp(s1))**2+2*tmp1*tmp1/(1.0+exp(s1))**3)
 
          extfld(4) = -this%Param(3)*(bbp*pos(2))*sin(alpha)
          extfld(5) = this%Param(3)*(bb-bbpp/2*pos(2)*pos(2))
          extfld(6) = this%Param(3)*(bbp*pos(2))*cos(alpha)
        endif
 
!        print*,"dipole: ",pos,extfld(5)
        end subroutine getfldt_Dipole

      end module Dipoleclass
