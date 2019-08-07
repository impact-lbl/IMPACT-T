!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! EMfldCartclass: ElectroMagnetic field data container class
!             in Lattice module of APPLICATION layer.
! MODULE  : ... EMfldCartclass
! VERSION : ... 1.0
!> @author  
!> Ji Qiang
! DESCRIPTION: 
!> This class contains discrete EM field data (as a function of
!> x,y,z) or (r,z) and analytical representation of EM field data (user can
!> supply the function form). The linear transfer map is also
!> computed base on the field on the axis. 
! Comments:
!----------------------------------------------------------------
      module EMfldCartclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 11
        type EMfldCart
          !Itype = 111
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : scale
          !      (3) : RF frequency
          !      (4) : theta0
          !      (5) : file ID
          !      (6) : radius
          !      (7) : x misalignment error
          !      (8) : y misalignment error
          !      (9) : rotation error x
          !      (10) : rotation error y
          !      (11) : rotation error z
        end type EMfldCart
        interface getparam_EMfldCart
          module procedure getparam1_EMfldCart,  &
                          getparam2_EMfldCart,   &
                          getparam3_EMfldCart
        end interface
        interface setparam_EMfldCart
          module procedure setparam1_EMfldCart,  &
                           setparam2_EMfldCart, setparam3_EMfldCart
        end interface
      contains
        subroutine construct_EMfldCart(this,numseg,nmpstp,type,blength)
        implicit none
        type (EMfldCart), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_EMfldCart
   
        subroutine setparam1_EMfldCart(this,i,value)
        implicit none
        type (EMfldCart), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_EMfldCart

        subroutine setparam2_EMfldCart(this,values)
        implicit none
        type (EMfldCart), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_EMfldCart

        subroutine setparam3_EMfldCart(this,numseg,nmpstp,type,blength)
        implicit none
        type (EMfldCart), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_EMfldCart

        subroutine getparam1_EMfldCart(this,i,blparam) 
        implicit none 
        type (EMfldCart), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_EMfldCart
  
        subroutine getparam2_EMfldCart(this,blparams)
        implicit none
        type (EMfldCart), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_EMfldCart

        subroutine getparam3_EMfldCart(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (EMfldCart), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_EMfldCart
       
        subroutine maplinear_EMfldCart(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (EMfldCart), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,escale,ww,theta0,qmcc
        double precision :: ein,e1in,e2in,ef,e1f,e2f
        double precision :: sinthi,costhi,uprimi,qpwi,tfin,sinthf,costhf
        double precision :: uprimf,qpwf
        integer :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        escale = this%Param(2)
        ww = this%Param(3)/Scfreq
        theta0 = this%Param(4)*asin(1.0)/90
        qmcc = Bchg/Bmass
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

!        call getaxfldE_EMfldCart(t,this,ein,e1in,e2in)
        call getaxfldEfc_EMfldCart(t,this,ein,e1in,e2in)
        sinthi=sin(ww*refpt(5)+theta0)
        costhi=cos(ww*refpt(5)+theta0)
        uprimi=qmcc*ein/betai*costhi
        qpwi=0.5*qmcc*Scxl/(ui*ww)
        dlti=Scxl*(0.5*uprimi/ui-qpwi*e1in*sinthi)
        thli=1.5*Scxl*(uprimi/ui)

        call rk6i_EMfldCart(h,mpstp,t0,y,18,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3
        tfin = t + tau
!        call getaxfldE_EMfldCart(tfin,this,ef,e1f,e2f)
        call getaxfldEfc_EMfldCart(tfin,this,ef,e1f,e2f)
        sinthf=sin(ww*refpt(5)+theta0)
        costhf=cos(ww*refpt(5)+theta0)
        uprimf=qmcc*ef/betaf*costhf
        qpwf=0.5*qmcc*Scxl/(uf*ww)
        dltf=Scxl*(0.5*uprimf/uf-qpwf*e1f*sinthf)
        thlf=1.5*Scxl*(uprimf/uf)

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

        end subroutine maplinear_EMfldCart

        subroutine rk6i_EMfldCart(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (EMfldCart), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_EMfldCart(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_EMfldCart(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_EMfldCart(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_EMfldCart(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_EMfldCart(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_EMfldCart(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_EMfldCart(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_EMfldCart(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_EMfldCart

        subroutine intfunc1_EMfldCart(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (EMfldCart), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55,s11tmp
        double precision :: zedge,escale,ez1,ezp1,ezpp1,ww,theta0,qmcc,&
                            sinphi,cosphi,rfdsgn
        integer :: my_rank, ierr

        zedge = this%Param(1)
!        call getaxfldE_EMfldCart(t,this,ez1,ezp1,ezpp1)
        call getaxfldEfc_EMfldCart(t,this,ez1,ezp1,ezpp1)
        ww = this%Param(3)/Scfreq 
        theta0 = this%Param(4)*asin(1.0)/90 
        qmcc = Bchg/Bmass

        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        f(4) = 0.0

        ! synchronous particle:
        gamma0=-y(6)
        beta0=sqrt(gamma0**2-1.)/gamma0
        gbet=sqrt((gamma0-1.0)*(gamma0+1.0))
        sinphi=sin(ww*y(5)+theta0)
        cosphi=cos(ww*y(5)+theta0)
        rfdsgn=ez1*cosphi
        f(5)=1.0/(beta0*Scxl)
        f(6)=-qmcc*rfdsgn

        ! matrix elements
        s11tmp=0.5e0*(1.0+0.5*gamma0**2)* &
         (qmcc*rfdsgn/beta0**2/gamma0**2)**2 + &
          qmcc*0.5/beta0**3/gamma0**3*ez1*sinphi*ww/Scxl
        s11=s11tmp*Scxl
        s33=s11tmp*Scxl
        s55=  &
         -1.5e0*qmcc/beta0**2/gamma0*ezp1*cosphi &
         +(beta0**2+0.5)/beta0**3/gamma0*qmcc*ez1*sinphi*ww/Scxl &
         +1.5e0*(1.d0-0.5*gamma0**2)* &
         (qmcc*rfdsgn/beta0**2/gamma0**2)**2
        s55=Scxl*s55

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

        end subroutine intfunc1_EMfldCart

        !--------------------------------------------------------------------------------------
        !> @brief
        !> interpolate the field from the EMfldCart rf cavity onto bunch location.
        !--------------------------------------------------------------------------------------
        subroutine getaxfldE_EMfldCart(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        type (EMfldCart), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1
        double precision:: zz,hstep,slope,zedge,escale
        integer :: klo,khi,k
        integer :: my_rank,ierr

        zedge = this%Param(1)
        escale = this%Param(2)
        zz=z-zedge
        klo=1
        khi=Ndata
1       if(khi-klo.gt.1) then
          k=(khi+klo)/2
          if(zdat(k).gt.zz)then
             khi=k
          else
             klo=k
          endif
          goto 1
        endif
        hstep=zdat(khi)-zdat(klo)
        slope=(edat(khi)-edat(klo))/hstep
        ez1 =edat(klo)+slope*(zz-zdat(klo))
        slope=(epdat(khi)-epdat(klo))/hstep
        ezp1=epdat(klo)+slope*(zz-zdat(klo))
        slope=(eppdat(khi)-eppdat(klo))/hstep
        ezpp1=eppdat(klo)+slope*(zz-zdat(klo))
        ez1 = ez1*escale
        ezp1 = ezp1*escale
        ezpp1 = ezpp1*escale

!        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
!        if(my_rank.eq.1) then
!         write(19,101)z,ez1/1.0e6,ezp1/1.0e6,zz,float(m),float(init)
!        endif
!        if(my_rank.eq.1) then
!         write(19,101)z,gt
!        endif
  101   format(6(1x,1pe13.6))

        end subroutine getaxfldE_EMfldCart

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external RF field on axis from analytical function
        !> Here, we have used a Fouier expansion representation of external field.
        !> Users should supply the field function as they want.
        !--------------------------------------------------------------------------------------
        subroutine  getaxfldEfc_EMfldCart(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (EMfldCart), intent(in) :: this
        double precision, intent(out) :: ez1,ezp1,ezpp1
        double precision :: zedge,escale,len,zz
        double precision :: pi,zmid
        integer :: i

        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2

        zz=z-zmid
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        !Ndata has to be odd number
        do i = 2, (Ndata-1)/2+1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale

        end subroutine getaxfldEfc_EMfldCart
        
        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field Ex, Ey, Ez, Bx, Bx, Bz at given position x, y, z, t from
        !> analytical function. Here we have used Fourier expansion of function. The
        !> user should supply his own analytical function if needed.
        !--------------------------------------------------------------------------------------
        subroutine  getfld_EMfldCart(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (EMfldCart), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz,r2
        double precision :: c5,c6,s5,s6,zmid
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        extfld = 0.0
        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq 
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        !move into the tank local coordinate.
        !zz=pos(3)-zedge
!-----------------------------------------------------------------
! get external RF field on axis from analytical function
! Here, we have used a Fouier expansion representation of external field.
! Users should supply the field function as they want.
        zz=pos(3)-zmid
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        !Ndata has to be odd number
        do i = 2, (Ndata-1)/2+1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
! calculate the E and B field using the analytical ez0 on the axis.
! This is based on the assumption of azmuthal symmetry and TM mode.
        tt = pos(4) 
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
        r2 = pos(1)**2+pos(2)**2
        tmpex = -pos(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -pos(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = pos(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -pos(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpbz = 0.0

        extfld(1) = extfld(1) + tmpex
        extfld(2) = extfld(2) + tmpey
        extfld(3) = extfld(3) + tmpez
        extfld(4) = extfld(4) + tmpbx
        extfld(5) = extfld(5) + tmpby
        extfld(6) = extfld(6) + tmpbz

        end subroutine getfld_EMfldCart

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getflderrold_EMfldCart(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (EMfldCart), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz,r2
        double precision :: zmid,dx,dy,anglex,angley,anglez
        double precision, dimension(3) :: temp,tmp
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq 
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90

        dx = this%Param(7)
        dy = this%Param(8)
        anglex = this%Param(9)
        angley = this%Param(10)
        anglez = this%Param(11)

        !move into the tank local coordinate.
        !zz=pos(3)-zedge
        !zz=pos(3)-zmid

        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = pos(3)-zmid
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

!-----------------------------------------------------------------
! get external RF field on axis from analytical function
! Here, we have used a Fouier expansion representation of external field.
! Users should supply the field function as they want.
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        !Ndata has to be odd number
        do i = 2, (Ndata-1)/2+1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
! calculate the E and B field using the analytical ez0 on the axis.
! This is based on the assumption of azmuthal symmetry and TM mode.
        tt = pos(4) 
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
        r2 = tmp(1)**2+tmp(2)**2
        extfld(1) = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
        extfld(2) = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
        extfld(3) = (ez1+f1*r2)*tmpcos
        extfld(4) = tmp(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        extfld(5) = -tmp(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        extfld(6) = 0.0

        !for E field
        tmp(1) = extfld(1)
        tmp(2) = extfld(2)*cos(anglex)-extfld(3)*sin(anglex)
        tmp(3) = extfld(2)*sin(anglex)+extfld(3)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(2) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(3) = temp(3)

        !for B field
        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getflderrold_EMfldCart

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine  getflderr_EMfldCart(pos,extfld,this,dx,dy,anglex,angley,anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (EMfldCart), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz,r2
        double precision :: zmid
        double precision, dimension(3) :: temp,tmp
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq 
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90

        !move into the tank local coordinate.
        !zz=pos(3)-zedge
        !zz=pos(3)-zmid

        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = pos(3) - zmid
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

!-----------------------------------------------------------------
! get external RF field on axis from analytical function
! Here, we have used a Fouier expansion representation of external field.
! Users should supply the field function as they want.
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        !Ndata has to be odd number
        do i = 2, (Ndata-1)/2+1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
! calculate the E and B field using the analytical ez0 on the axis.
! This is based on the assumption of azmuthal symmetry and TM mode.
        tt = pos(4) 
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
        r2 = tmp(1)**2+tmp(2)**2
        extfld(1) = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
        extfld(2) = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
        extfld(3) = (ez1+f1*r2)*tmpcos
        extfld(4) = tmp(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        extfld(5) = -tmp(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        extfld(6) = 0.0

        !for E field
        tmp(1) = extfld(1)
        tmp(2) = extfld(2)*cos(anglex)-extfld(3)*sin(anglex)
        tmp(3) = extfld(2)*sin(anglex)+extfld(3)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(2) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(3) = temp(3)

        !for B field
        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getflderr_EMfldCart

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get the discrete Er,Etheta,Ez, Br, Btheta, Bz as a function or "r" at given "z".
        !--------------------------------------------------------------------------------------
        subroutine getfld6_EMfldCart(this,z,extfld6)
        type (EMfldCart), intent(in) :: this
        double precision, intent(in) :: z
        double precision, dimension(6,NrIntvRf+1), intent(out) :: extfld6
        double precision :: zz,hz,ef

        hz = (ZmaxRf-ZminRf)/NzIntvRf
        zz = z-this%Param(1)
        iz = zz/hz + 1 
        if(iz.eq.(NzIntvRf+1)) then
           iz = iz - 1
        endif
        iz1 = iz+1
        ef = (iz*hz - zz)/hz
        
! get the external field along "r" at given "z" by linear interpolation
! Here, we have only Er, Ez, and Btheta for TM mode.
        do j = 1, NrIntvRf+1
          extfld6(1,j) = erdata(iz,j)*ef + erdata(iz1,j)*(1.d0-ef)  
          extfld6(2,j) = 0.0 
          extfld6(3,j) = ezdata(iz,j)*ef + ezdata(iz1,j)*(1.d0-ef)  
          extfld6(4,j) = 0.0
          extfld6(5,j) = btdata(iz,j)*ef + btdata(iz1,j)*(1.d0-ef)
          extfld6(6,j) = 0.0
        enddo

        end subroutine getfld6_EMfldCart
        
        subroutine getfld6err_EMfldCart(this,pos,extfld,dx,dy,anglex,angley,&
                                    anglez)
        type (EMfldCart), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos 
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zz,hz,efz,hr,efr,len,zedge,escale,ww,theta0,&
                            tt,er,ez,bt,tmpcos,tmpsin
        double precision, dimension(3) :: tmp,temp
        integer :: iz,iz1,ir,ir1

        len = this%Length
        zedge = this%Param(1)
        escale = this%Param(2)
        ww = this%Param(3)/Scfreq
        theta0 = this%Param(4)*asin(1.0)/90

        hz = (ZmaxRf-ZminRf)/NzIntvRf
        hr = (RmaxRf-RminRf)/NrIntvRf

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
        tt = pos(4)
        tmpcos = cos(ww*tt+theta0)*escale
        tmpsin = sin(ww*tt+theta0)*escale

        if(zz.le.0) then
          iz = 1
          efz = zz/hz
        else if(zz.ge.(ZmaxRf-ZminRf)) then
          iz = NzIntvRf
          efz = (zz-(iz-1)*hz)/hz
        else
          iz = zz/hz + 1 
          efz = (zz-(iz-1)*hz)/hz
        endif
        iz1 = iz+1
        
        rr = sqrt(tmp(1)*tmp(1)+tmp(2)*tmp(2))
        ir = rr/hr + 1
        if(rr.eq.0.0) then
          rr = 1.0e-10
        endif
        if(ir.gt.NrIntvRf) then
           print*,"ir: ",ir,rr,tmp(1),tmp(2)
           stop
        endif
        efr = (rr-(ir-1)*hr)/hr
        ir1 = ir+1

        ! get the external field along "r" at given "z" by linear interpolation
        ! Here, we have only Er, Ez, and Btheta for TM mode.
        er = (erdata(iz,ir)*(1.d0-efz)*(1.d0-efr)+erdata(iz1,ir)*efz*(1.d0-efr)+&
             erdata(iz1,ir1)*efz*efr+erdata(iz,ir1)*(1.d0-efz)*efr)*tmpcos
        ez = (ezdata(iz,ir)*(1.d0-efz)*(1.d0-efr)+ezdata(iz1,ir)*efz*(1.d0-efr)+&
             ezdata(iz1,ir1)*efz*efr+ezdata(iz,ir1)*(1.d0-efz)*efr)*tmpcos
        bt = (btdata(iz,ir)*(1.d0-efz)*(1.d0-efr)+btdata(iz1,ir)*efz*(1.d0-efr)+&
             btdata(iz1,ir1)*efz*efr+btdata(iz,ir1)*(1.d0-efz)*efr)*tmpsin

        extfld(1) = er*tmp(1)/rr
        extfld(2) = er*tmp(2)/rr
        extfld(3) = ez
        extfld(4) = -bt*tmp(2)/rr
        extfld(5) = bt*tmp(1)/rr
        extfld(6) = 0.0
       
        !for E field
        tmp(1) = extfld(1)
        tmp(2) = extfld(2)*cos(anglex)-extfld(3)*sin(anglex)
        tmp(3) = extfld(2)*sin(anglex)+extfld(3)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(2) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(3) = temp(3)

        !for B field
        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getfld6err_EMfldCart

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get the discrete Ex,Ey,Ez, Bx, By, Bz as a function of x and y at given "z".
        !--------------------------------------------------------------------------------------
        subroutine getfld6xyz_EMfldCart(this,z,extfld6xyz)
        type (EMfldCart), intent(in) :: this
        double precision, intent(in) :: z
        double precision, dimension(6,NxIntvRfg+1,NyIntvRfg+1),&
                          intent(out) :: extfld6xyz
        double precision :: zz,hz,ef

        hz = (ZmaxRfg-ZminRfg)/NzIntvRfg
        zz = z-this%Param(1)
        iz = zz/hz + 1 
        if(iz.eq.(NzIntvRfg+1)) then
           iz = iz - 1
        endif
        iz1 = iz+1
        ef = (iz*hz - zz)/hz
        
! get the external field along "r" at given "z" by linear interpolation
! Here, we have only Er, Ez, and Btheta for TM mode.
        do j = 1, NyIntvRfg+1
          do i = 1, NxIntvRfg+1
            extfld6xyz(1,i,j) = Exgrid(i,j,iz)*ef + Exgrid(i,j,iz1)*(1.d0-ef)  
            extfld6xyz(2,i,j) = Eygrid(i,j,iz)*ef + Eygrid(i,j,iz1)*(1.d0-ef)  
            extfld6xyz(3,i,j) = Ezgrid(i,j,iz)*ef + Ezgrid(i,j,iz1)*(1.d0-ef)  
            extfld6xyz(4,i,j) = Bxgrid(i,j,iz)*ef + Bxgrid(i,j,iz1)*(1.d0-ef)  
            extfld6xyz(5,i,j) = Bygrid(i,j,iz)*ef + Bygrid(i,j,iz1)*(1.d0-ef)  
            extfld6xyz(6,i,j) = Bzgrid(i,j,iz)*ef + Bzgrid(i,j,iz1)*(1.d0-ef)  
          enddo
        enddo

        end subroutine getfld6xyz_EMfldCart
        
        subroutine getfld6xyzerr_EMfldCart(this,pos,extfld,dx,dy,anglex,angley,&
                                    anglez)
        type (EMfldCart), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos 
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zz,hz,efz,hx,efx,hy,efy,len,zedge,escale,ww,theta0,&
                        tt,exn,eyn,ezn,bxn,byn,bzn,tmpcos,tmpsin
        double precision, dimension(3) :: tmp,temp
        integer :: iz,iz1,ix,ix1,iy,iy1

        len = this%Length
        zedge = this%Param(1)
        escale = this%Param(2)
        ww = this%Param(3)/Scfreq
        theta0 = this%Param(4)*asin(1.0)/90

        hx = (XmaxRfg-XminRfg)/NxIntvRfg
        hy = (YmaxRfg-YminRfg)/NyIntvRfg
        hz = (ZmaxRfg-ZminRfg)/NzIntvRfg

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
        tt = pos(4)
        tmpcos = cos(ww*tt+theta0)*escale
        tmpsin = sin(ww*tt+theta0)*escale

        if(tmp(1).le.XminRfg) then
          ix = 1
          efx = (XminRfg-tmp(1)+ix*hx)/hx
        else if(tmp(1).ge.XmaxRfg) then
          ix = NxIntvRfg
          efx = (XminRfg-tmp(1)+ix*hx)/hx
        else
          ix = (tmp(1)-XminRfg)/hx + 1 
          efx = (XminRfg-tmp(1)+ix*hx)/hx
        endif
        ix1 = ix+1
        
        if(tmp(2).le.YminRfg) then
          iy = 1
          efy = (YminRfg-tmp(2)+iy*hy)/hy
        else if(tmp(2).ge.YmaxRfg) then
          iy = NyIntvRfg
          efy = (YminRfg-tmp(2)+iy*hy)/hy
        else
          iy = (tmp(2)-YminRfg)/hy + 1 
          efy = (YminRfg-tmp(2)+iy*hy)/hy
        endif
        iy1 = iy+1

        if(zz.le.0.0) then
          iz = 1
          efz = 1.0 - zz/hz
        else if(zz.ge.(ZmaxRfg-ZminRfg)) then
          iz = NzIntvRfg
          efz = 1.0 - (zz-(iz-1)*hz)/hz
        else
          iz = zz/hz + 1 
          efz = 1.0 - (zz-(iz-1)*hz)/hz
        endif
        iz1 = iz+1

        exn = Exgrid(ix,iy,iz)*efx*efy*efz + &
          Exgrid(ix,iy1,iz)*efx*(1-efy)*efz  + &
          Exgrid(ix,iy1,iz1)*efx*(1-efy)*(1-efz) + &
          Exgrid(ix,iy,iz1)*efx*efy*(1-efz) + &
          Exgrid(ix1,iy,iz1)*(1-efx)*efy*(1-efz) + &
          Exgrid(ix1,iy1,iz1)*(1-efx)*(1-efy)*(1-efz) + &
          Exgrid(ix1,iy1,iz)*(1-efx)*(1-efy)*efz + &
          Exgrid(ix1,iy,iz)*(1-efx)*efy*efz 
        eyn = Eygrid(ix,iy,iz)*efx*efy*efz + &
          Eygrid(ix,iy1,iz)*efx*(1-efy)*efz  + &
          Eygrid(ix,iy1,iz1)*efx*(1-efy)*(1-efz) + &
          Eygrid(ix,iy,iz1)*efx*efy*(1-efz) + &
          Eygrid(ix1,iy,iz1)*(1-efx)*efy*(1-efz) + &
          Eygrid(ix1,iy1,iz1)*(1-efx)*(1-efy)*(1-efz) + &
          Eygrid(ix1,iy1,iz)*(1-efx)*(1-efy)*efz + &
          Eygrid(ix1,iy,iz)*(1-efx)*efy*efz 
        ezn = Ezgrid(ix,iy,iz)*efx*efy*efz + &
          Ezgrid(ix,iy1,iz)*efx*(1-efy)*efz  + &
          Ezgrid(ix,iy1,iz1)*efx*(1-efy)*(1-efz) + &
          Ezgrid(ix,iy,iz1)*efx*efy*(1-efz) + &
          Ezgrid(ix1,iy,iz1)*(1-efx)*efy*(1-efz) + &
          Ezgrid(ix1,iy1,iz1)*(1-efx)*(1-efy)*(1-efz) + &
          Ezgrid(ix1,iy1,iz)*(1-efx)*(1-efy)*efz + &
          Ezgrid(ix1,iy,iz)*(1-efx)*efy*efz 
        bxn = Bxgrid(ix,iy,iz)*efx*efy*efz + &
          Bxgrid(ix,iy1,iz)*efx*(1-efy)*efz  + &
          Bxgrid(ix,iy1,iz1)*efx*(1-efy)*(1-efz) + &
          Bxgrid(ix,iy,iz1)*efx*efy*(1-efz) + &
          Bxgrid(ix1,iy,iz1)*(1-efx)*efy*(1-efz) + &
          Bxgrid(ix1,iy1,iz1)*(1-efx)*(1-efy)*(1-efz) + &
          Bxgrid(ix1,iy1,iz)*(1-efx)*(1-efy)*efz + &
          Bxgrid(ix1,iy,iz)*(1-efx)*efy*efz 
        byn = Bygrid(ix,iy,iz)*efx*efy*efz + &
          Bygrid(ix,iy1,iz)*efx*(1-efy)*efz  + &
          Bygrid(ix,iy1,iz1)*efx*(1-efy)*(1-efz) + &
          Bygrid(ix,iy,iz1)*efx*efy*(1-efz) + &
          Bygrid(ix1,iy,iz1)*(1-efx)*efy*(1-efz) + &
          Bygrid(ix1,iy1,iz1)*(1-efx)*(1-efy)*(1-efz) + &
          Bygrid(ix1,iy1,iz)*(1-efx)*(1-efy)*efz + &
          Bygrid(ix1,iy,iz)*(1-efx)*efy*efz 
        bzn = Bzgrid(ix,iy,iz)*efx*efy*efz + &
          Bzgrid(ix,iy1,iz)*efx*(1-efy)*efz  + &
          Bzgrid(ix,iy1,iz1)*efx*(1-efy)*(1-efz) + &
          Bzgrid(ix,iy,iz1)*efx*efy*(1-efz) + &
          Bzgrid(ix1,iy,iz1)*(1-efx)*efy*(1-efz) + &
          Bzgrid(ix1,iy1,iz1)*(1-efx)*(1-efy)*(1-efz) + &
          Bzgrid(ix1,iy1,iz)*(1-efx)*(1-efy)*efz + &
          Bzgrid(ix1,iy,iz)*(1-efx)*efy*efz 

        extfld(1) = exn*tmpcos
        extfld(2) = eyn*tmpcos
        extfld(3) = ezn*tmpcos
        extfld(4) = bxn*tmpsin
        extfld(5) = byn*tmpsin
        extfld(6) = bzn*tmpsin
       
        !transform 
        !for E field
        tmp(1) = extfld(1)
        tmp(2) = extfld(2)*cos(anglex)-extfld(3)*sin(anglex)
        tmp(3) = extfld(2)*sin(anglex)+extfld(3)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(2) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(3) = temp(3)

        !for B field
        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getfld6xyzerr_EMfldCart

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get the discrete Ex,Ey,Ez, Bx, By, Bz at given x, y, z, t.
        !--------------------------------------------------------------------------------------
        subroutine  getfldt_EMfldCart(pos,extfld,this,fldata)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (EMfldCart), intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld
        double precision :: hx,hy,hz,zz,ab,cd,ef
        integer :: ix,jx,kx,ix1,jx1,kx1
        double precision :: escale,ww,theta0,tt,tmpcos,tmpsin,z
        double complex :: tmptime

        z = pos(3)
        if((z.gt.this%Param(1)).and.(z.lt.(this%Param(1)+this%Length))) then 
          escale = this%Param(2)
          ww = this%Param(3)*4*asin(1.0)
          theta0 = this%Param(4)*asin(1.0)/90
          tt = pos(4)
!          tmpcos = cos(ww*tt+theta0)*escale
          !tmpsin = sin(ww*tt+theta0)*escale
!          tmpsin = -sin(ww*tt+theta0)*escale
          tmptime = dcmplx(cos(ww*tt+theta0),sin(ww*tt+theta0))*escale
          hx = (fldata%XmaxRfgt-fldata%XminRfgt)/fldata%NxIntvRfgt
          ix = (pos(1)-fldata%XminRfgt)/hx + 1
          ab = ((fldata%XminRfgt-pos(1))+ix*hx)/hx
          hy = (fldata%YmaxRfgt-fldata%YminRfgt)/fldata%NyIntvRfgt
          jx = (pos(2)-fldata%YminRfgt)/hy + 1
          cd = ((fldata%YminRfgt-pos(2))+jx*hy)/hy
          hz = (fldata%ZmaxRfgt-fldata%ZminRfgt)/fldata%NzIntvRfgt
          zz = pos(3)-this%Param(1)+fldata%ZminRfgt
          kx = (zz-fldata%ZminRfgt)/hz + 1
          ef = ((fldata%ZminRfgt-zz)+kx*hz)/hz
          ix1 = ix + 1
          jx1 = jx + 1
          kx1 = kx + 1

!          print*,"EMCART, ix,jx,kx: ",ix,jx,kx

          extfld(1) = dreal((fldata%Exgridt(ix,jx,kx)*ab*cd*ef  &
                  +fldata%Exgridt(ix,jx1,kx)*ab*(1.d0-cd)*ef &
                  +fldata%Exgridt(ix,jx1,kx1)*ab*(1.d0-cd)*(1.d0-ef) &
                  +fldata%Exgridt(ix,jx,kx1)*ab*cd*(1.d0-ef) &
                  +fldata%Exgridt(ix1,jx,kx1)*(1.d0-ab)*cd*(1.d0-ef) &
                  +fldata%Exgridt(ix1,jx1,kx1)*(1.d0-ab)*(1.d0-cd)*(1.d0-ef)&
                  +fldata%Exgridt(ix1,jx1,kx)*(1.d0-ab)*(1.d0-cd)*ef &
                  +fldata%Exgridt(ix1,jx,kx)*(1.d0-ab)*cd*ef)*tmptime)
 
          extfld(2) = dreal((fldata%Eygridt(ix,jx,kx)*ab*cd*ef  &
                  +fldata%Eygridt(ix,jx1,kx)*ab*(1.d0-cd)*ef &
                  +fldata%Eygridt(ix,jx1,kx1)*ab*(1.d0-cd)*(1.d0-ef) &
                  +fldata%Eygridt(ix,jx,kx1)*ab*cd*(1.d0-ef) &
                  +fldata%Eygridt(ix1,jx,kx1)*(1.d0-ab)*cd*(1.d0-ef) &
                  +fldata%Eygridt(ix1,jx1,kx1)*(1.d0-ab)*(1.d0-cd)*(1.d0-ef)&
                  +fldata%Eygridt(ix1,jx1,kx)*(1.d0-ab)*(1.d0-cd)*ef &
                  +fldata%Eygridt(ix1,jx,kx)*(1.d0-ab)*cd*ef)*tmptime)
 
          extfld(3) = dreal((fldata%Ezgridt(ix,jx,kx)*ab*cd*ef  &
                  +fldata%Ezgridt(ix,jx1,kx)*ab*(1.d0-cd)*ef &
                  +fldata%Ezgridt(ix,jx1,kx1)*ab*(1.d0-cd)*(1.d0-ef) &
                  +fldata%Ezgridt(ix,jx,kx1)*ab*cd*(1.d0-ef) &
                  +fldata%Ezgridt(ix1,jx,kx1)*(1.d0-ab)*cd*(1.d0-ef) &
                  +fldata%Ezgridt(ix1,jx1,kx1)*(1.d0-ab)*(1.d0-cd)*(1.d0-ef)&
                  +fldata%Ezgridt(ix1,jx1,kx)*(1.d0-ab)*(1.d0-cd)*ef &
                  +fldata%Ezgridt(ix1,jx,kx)*(1.d0-ab)*cd*ef)*tmptime)

          extfld(4) = dreal((fldata%Bxgridt(ix,jx,kx)*ab*cd*ef  &
                  +fldata%Bxgridt(ix,jx1,kx)*ab*(1.d0-cd)*ef &
                  +fldata%Bxgridt(ix,jx1,kx1)*ab*(1.d0-cd)*(1.d0-ef) &
                  +fldata%Bxgridt(ix,jx,kx1)*ab*cd*(1.d0-ef) &
                  +fldata%Bxgridt(ix1,jx,kx1)*(1.d0-ab)*cd*(1.d0-ef) &
                  +fldata%Bxgridt(ix1,jx1,kx1)*(1.d0-ab)*(1.d0-cd)*(1.d0-ef)&
                  +fldata%Bxgridt(ix1,jx1,kx)*(1.d0-ab)*(1.d0-cd)*ef &
                  +fldata%Bxgridt(ix1,jx,kx)*(1.d0-ab)*cd*ef)*tmptime)
 
          extfld(5) = dreal((fldata%Bygridt(ix,jx,kx)*ab*cd*ef  &
                  +fldata%Bygridt(ix,jx1,kx)*ab*(1.d0-cd)*ef &
                  +fldata%Bygridt(ix,jx1,kx1)*ab*(1.d0-cd)*(1.d0-ef) &
                  +fldata%Bygridt(ix,jx,kx1)*ab*cd*(1.d0-ef) &
                  +fldata%Bygridt(ix1,jx,kx1)*(1.d0-ab)*cd*(1.d0-ef) &
                  +fldata%Bygridt(ix1,jx1,kx1)*(1.d0-ab)*(1.d0-cd)*(1.d0-ef)&
                  +fldata%Bygridt(ix1,jx1,kx)*(1.d0-ab)*(1.d0-cd)*ef &
                  +fldata%Bygridt(ix1,jx,kx)*(1.d0-ab)*cd*ef)*tmptime)
 
          extfld(6) = dreal((fldata%Bzgridt(ix,jx,kx)*ab*cd*ef  &
                  +fldata%Bzgridt(ix,jx1,kx)*ab*(1.d0-cd)*ef &
                  +fldata%Bzgridt(ix,jx1,kx1)*ab*(1.d0-cd)*(1.d0-ef) &
                  +fldata%Bzgridt(ix,jx,kx1)*ab*cd*(1.d0-ef) &
                  +fldata%Bzgridt(ix1,jx,kx1)*(1.d0-ab)*cd*(1.d0-ef) &
                  +fldata%Bzgridt(ix1,jx1,kx1)*(1.d0-ab)*(1.d0-cd)*(1.d0-ef)&
                  +fldata%Bzgridt(ix1,jx1,kx)*(1.d0-ab)*(1.d0-cd)*ef &
                  +fldata%Bzgridt(ix1,jx,kx)*(1.d0-ab)*cd*ef)*tmptime)
        else
          extfld = 0.0
        endif
!        write(9,100)pos(1),pos(2),pos(3),ix*1.0,jx*1.0,kx*1.0,&
!                    extfld(1:6),fldata%Exgridt(ix,jx,kx),fldata%Eygridt(ix,jx,kx),&
!                    fldata%Ezgridt(ix,jx,kx)
!100     format(15(1x,e14.6))
!        call flush_(9)

        end subroutine getfldt_EMfldCart
      end module EMfldCartclass
