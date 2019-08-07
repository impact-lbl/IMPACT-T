!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! PhysConstclass: physical constants class in CONSTANTS module of 
!                 DATA STRUCTURE layer.
! MODULE  : ... PhysConstclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang 
!
! DESCRIPTION: 
!> This class defines the physical constant parameters used
!> in the simulation.
! Comments:
!----------------------------------------------------------------
      module PhysConstclass
        use mpistub
        implicit none
      
        !physical parameters and constants ---------------------
        double precision :: Pi
        double precision :: Clight !< speed of light in vaccum
        double precision :: Scxl !< length scale in z domain code
        double precision :: Rad2deg !< conversion factor from radian to degree
        double precision :: Epsilon0 !< permittivity of vacuum
        double precision :: Scfreq !< time scale
        double precision :: Scxlt !< length scale in time domain code
        double precision :: Scq0 !< elementary charge
        double precision :: Scphi0 !< electical potential
        double precision :: Sct0  !< time scale
        double precision :: Scm0 !< mass scale
        double precision :: Sce0 !< electrical field scale
        double precision :: Scb0 !< magnetic field scale
        double precision :: Scrho0 !< charge density scale
      contains
        subroutine constructT_PhysConst(dt,reffreq)
        double precision, intent(in) :: dt,reffreq

        Clight = 299792458.0d0 
        Pi = 2*asin(1.0d0)
        Scxlt = Clight*dt
        Rad2deg = 180.0d0/Pi
        Epsilon0 = 8.854187817d-12
        Scq0 = 1.60217733d-19
        Scphi0 = 1.0d0
        Scxl = 0.0
        Scfreq = reffreq
        Sct0 = dt
        Scm0 = Scq0/(Clight*Clight)
        Sce0 = Scphi0/Scxlt
        Scb0 = Scm0/(Scxlt*Scq0)
        Scrho0 = Epsilon0/(Scxlt*Scxlt)

        end subroutine constructT_PhysConst
 
      end module PhysConstclass
