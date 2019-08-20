!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! BeamLineElemclass: Beam line element base class in Lattice module of APPLICATION layer.
!
! MODULE  : ... BeamLineElemclass
! VERSION : ... 1.0
! 
!> @author
!> Ji Qiang
!
! DESCRIPTION: 
!> This class defines the base beam line element class for different lattice element class.
! Comments:
!----------------------------------------------------------------
      module BeamLineElemclass
        use Quadrupoleclass 
        use Multipoleclass 
        use DriftTubeclass
        use CCLclass
        use CCDTLclass
        use DTLclass
        use SCclass
        use BPMclass
        use ConstFocclass
        use SolRFclass
        use Solclass
        use Dipoleclass
        use EMfldclass
        use EMfldAnaclass
        use EMfldCylclass
        use EMfldCartclass
        type BeamLineElem
!          private
          type (BPM), pointer :: pbpm
          type (CCL), pointer :: pccl
          type (CCDTL), pointer :: pccdtl
          type (DTL), pointer :: pdtl
          type (DriftTube), pointer :: pdrift
          type (Quadrupole), pointer :: pquad
          type (Multipole), pointer :: pmult
          type (SC), pointer :: psc
          type (ConstFoc), pointer :: pcf
          type (SolRF), pointer :: pslrf
          type (Sol), pointer :: psl
          type (Dipole), pointer :: pdipole
          type (EMfld), pointer :: pemfld
          type (EMfldAna), pointer :: pemfldana
          type (EMfldCart), pointer :: pemfldcart
          type (EMfldCyl), pointer :: pemfldcyl
        end type BeamLineElem
        interface assign_BeamLineElem
          module procedure assign_ccl,assign_ccdtl,assign_dtl,assign_quad,&
          assign_drift,assign_sc,assign_bpm,assign_cf,assign_slrf,assign_sl,&
          assign_dipole,assign_emfld,assign_emfldana,assign_emfldcart,&
          assign_emfldcyl,assign_mult
        end interface
        interface getparam_BeamLineElem
          module procedure getparam1_BeamLineElem, &
                           getparam2_BeamLineElem, &
                           getparam3_BeamLineElem
        end interface
        interface setparam_BeamLineElem
          module procedure setparam1_BeamLineElem, &
                           setparam2_BeamLineElem, &
                           setparam3_BeamLineElem
        end interface
      contains
        function assign_quad(tquad) result(ppquad)
        type (BeamLineElem) :: ppquad
        type (Quadrupole), target, intent(in) :: tquad

        ppquad%pquad => tquad
        nullify(ppquad%pdrift)
        nullify(ppquad%pccl)
        nullify(ppquad%pccdtl)
        nullify(ppquad%pdtl)
        nullify(ppquad%psc)
        nullify(ppquad%pbpm)
        nullify(ppquad%pcf)
        nullify(ppquad%pslrf)
        nullify(ppquad%psl)
        nullify(ppquad%pdipole)
        nullify(ppquad%pemfld)
        nullify(ppquad%pemfldana)
        nullify(ppquad%pemfldcart)
        nullify(ppquad%pemfldcyl)
        nullify(ppquad%pmult)

        end function assign_quad

        function assign_drift(tdrift) result(ppdrift)
        type (BeamLineElem) :: ppdrift
        type (DriftTube), target, intent(in) :: tdrift

        ppdrift%pdrift => tdrift
        nullify(ppdrift%pquad)
        nullify(ppdrift%pccl)
        nullify(ppdrift%pccdtl)
        nullify(ppdrift%pdtl)
        nullify(ppdrift%psc)
        nullify(ppdrift%pbpm)
        nullify(ppdrift%pcf)
        nullify(ppdrift%pslrf)
        nullify(ppdrift%psl)
        nullify(ppdrift%pdipole)
        nullify(ppdrift%pemfld)
        nullify(ppdrift%pemfldana)
        nullify(ppdrift%pemfldcart)
        nullify(ppdrift%pemfldcyl)
        nullify(ppdrift%pmult)

        end function assign_drift
         
        function assign_ccl(tccl) result(ppccl)
        type (BeamLineElem) :: ppccl
        type (CCL), target, intent(in) :: tccl

        ppccl%pccl => tccl
        nullify(ppccl%pquad)
        nullify(ppccl%pdrift)
        nullify(ppccl%pbpm)
        nullify(ppccl%pccdtl)
        nullify(ppccl%pdtl)
        nullify(ppccl%psc)
        nullify(ppccl%pcf)
        nullify(ppccl%pslrf)
        nullify(ppccl%psl)
        nullify(ppccl%pdipole)
        nullify(ppccl%pemfld)
        nullify(ppccl%pemfldana)
        nullify(ppccl%pemfldcart)
        nullify(ppccl%pemfldcyl)
        nullify(ppccl%pmult)

        end function assign_ccl
         
        function assign_ccdtl(tccdtl) result(ppccdtl)
        type (BeamLineElem) :: ppccdtl
        type (CCDTL), target, intent(in) :: tccdtl

        ppccdtl%pccdtl => tccdtl
        nullify(ppccdtl%pquad)
        nullify(ppccdtl%pdrift)
        nullify(ppccdtl%pbpm)
        nullify(ppccdtl%pccl)
        nullify(ppccdtl%pdtl)
        nullify(ppccdtl%psc)
        nullify(ppccdtl%pcf)
        nullify(ppccdtl%pslrf)
        nullify(ppccdtl%psl)
        nullify(ppccdtl%pdipole)
        nullify(ppccdtl%pemfld)
        nullify(ppccdtl%pemfldana)
        nullify(ppccdtl%pemfldcart)
        nullify(ppccdtl%pemfldcyl)
        nullify(ppccdtl%pmult)

        end function assign_ccdtl

        function assign_dtl(tdtl) result(ppdtl)
        type (BeamLineElem) :: ppdtl
        type (DTL), target, intent(in) :: tdtl

        ppdtl%pdtl => tdtl
        nullify(ppdtl%pquad)
        nullify(ppdtl%pdrift)
        nullify(ppdtl%pbpm)
        nullify(ppdtl%pccl)
        nullify(ppdtl%pccdtl)
        nullify(ppdtl%psc)
        nullify(ppdtl%pcf)
        nullify(ppdtl%pslrf)
        nullify(ppdtl%psl)
        nullify(ppdtl%pdipole)
        nullify(ppdtl%pemfld)
        nullify(ppdtl%pemfldana)
        nullify(ppdtl%pemfldcart)
        nullify(ppdtl%pemfldcyl)
        nullify(ppdtl%pmult)

        end function assign_dtl

        function assign_sc(tsc) result(ppsc)
        type (BeamLineElem) :: ppsc
        type (SC), target, intent(in) :: tsc

        ppsc%psc => tsc
        nullify(ppsc%pquad)
        nullify(ppsc%pdrift)
        nullify(ppsc%pbpm)
        nullify(ppsc%pccl)
        nullify(ppsc%pccdtl)
        nullify(ppsc%pdtl)
        nullify(ppsc%pcf)
        nullify(ppsc%pslrf)
        nullify(ppsc%psl)
        nullify(ppsc%pdipole)
        nullify(ppsc%pemfld)
        nullify(ppsc%pemfldana)
        nullify(ppsc%pemfldcart)
        nullify(ppsc%pemfldcyl)
        nullify(ppsc%pmult)

        end function assign_sc

        function assign_bpm(tbpm) result(ppbpm)
        type (BeamLineElem) :: ppbpm
        type (BPM), target, intent(in) :: tbpm

        ppbpm%pbpm => tbpm
        nullify(ppbpm%pquad)
        nullify(ppbpm%pdrift)
        nullify(ppbpm%pccl)
        nullify(ppbpm%pccdtl)
        nullify(ppbpm%pdtl)
        nullify(ppbpm%psc)
        nullify(ppbpm%pcf)
        nullify(ppbpm%pslrf)
        nullify(ppbpm%psl)
        nullify(ppbpm%pdipole)
        nullify(ppbpm%pemfld)
        nullify(ppbpm%pemfldana)
        nullify(ppbpm%pemfldcart)
        nullify(ppbpm%pemfldcyl)
        nullify(ppbpm%pmult)

        end function assign_bpm

        function assign_cf(tcf) result(ppcf)
        type (BeamLineElem) :: ppcf
        type (ConstFoc), target, intent(in) :: tcf

        ppcf%pcf => tcf
        nullify(ppcf%pbpm)
        nullify(ppcf%pquad)
        nullify(ppcf%pdrift)
        nullify(ppcf%pccl)
        nullify(ppcf%pccdtl)
        nullify(ppcf%pdtl)
        nullify(ppcf%psc)
        nullify(ppcf%pslrf)
        nullify(ppcf%psl)
        nullify(ppcf%pdipole)
        nullify(ppcf%pemfld)
        nullify(ppcf%pemfldana)
        nullify(ppcf%pemfldcart)
        nullify(ppcf%pemfldcyl)
        nullify(ppcf%pmult)

        end function assign_cf

        function assign_slrf(tslrf) result(ppslrf)
        type (BeamLineElem) :: ppslrf
        type (SolRF), target, intent(in) :: tslrf

        ppslrf%pslrf => tslrf
        nullify(ppslrf%pbpm)
        nullify(ppslrf%pquad)
        nullify(ppslrf%pdrift)
        nullify(ppslrf%pccl)
        nullify(ppslrf%pccdtl)
        nullify(ppslrf%pdtl)
        nullify(ppslrf%psc)
        nullify(ppslrf%pcf)
        nullify(ppslrf%psl)
        nullify(ppslrf%pdipole)
        nullify(ppslrf%pemfld)
        nullify(ppslrf%pemfldana)
        nullify(ppslrf%pemfldcart)
        nullify(ppslrf%pemfldcyl)
        nullify(ppslrf%pmult)

        end function assign_slrf

        function assign_sl(tsl) result(ppsl)
        type (BeamLineElem) :: ppsl
        type (Sol), target, intent(in) :: tsl

        ppsl%psl => tsl
        nullify(ppsl%pbpm)
        nullify(ppsl%pquad)
        nullify(ppsl%pdrift)
        nullify(ppsl%pccl)
        nullify(ppsl%pccdtl)
        nullify(ppsl%pdtl)
        nullify(ppsl%psc)
        nullify(ppsl%pcf)
        nullify(ppsl%pslrf)
        nullify(ppsl%pdipole)
        nullify(ppsl%pemfld)
        nullify(ppsl%pemfldana)
        nullify(ppsl%pemfldcart)
        nullify(ppsl%pemfldcyl)
        nullify(ppsl%pmult)

        end function assign_sl

        function assign_dipole(tdipole) result(ppdipole)
        type (BeamLineElem) :: ppdipole
        type (Dipole), target, intent(in) :: tdipole

        ppdipole%pdipole => tdipole
        nullify(ppdipole%psl)
        nullify(ppdipole%pbpm)
        nullify(ppdipole%pquad)
        nullify(ppdipole%pdrift)
        nullify(ppdipole%pccl)
        nullify(ppdipole%pccdtl)
        nullify(ppdipole%pdtl)
        nullify(ppdipole%psc)
        nullify(ppdipole%pcf)
        nullify(ppdipole%pslrf)
        nullify(ppdipole%pemfld)
        nullify(ppdipole%pemfldana)
        nullify(ppdipole%pemfldcart)
        nullify(ppdipole%pemfldcyl)
        nullify(ppdipole%pmult)

        end function assign_dipole

        function assign_emfld(temfld) result(ppemfld)
        type (BeamLineElem) :: ppemfld
        type (EMfld), target, intent(in) :: temfld

        ppemfld%pemfld => temfld
        nullify(ppemfld%psl)
        nullify(ppemfld%pbpm)
        nullify(ppemfld%pquad)
        nullify(ppemfld%pdrift)
        nullify(ppemfld%pccl)
        nullify(ppemfld%pccdtl)
        nullify(ppemfld%pdtl)
        nullify(ppemfld%psc)
        nullify(ppemfld%pcf)
        nullify(ppemfld%pslrf)
        nullify(ppemfld%pdipole)
        nullify(ppemfld%pemfldana)
        nullify(ppemfld%pemfldcart)
        nullify(ppemfld%pemfldcyl)
        nullify(ppemfld%pmult)

        end function assign_emfld

        function assign_emfldana(temfldana) result(ppemfldana)
        type (BeamLineElem) :: ppemfldana
        type (EMfldana), target, intent(in) :: temfldana

        ppemfldana%pemfldana => temfldana
        nullify(ppemfldana%psl)
        nullify(ppemfldana%pbpm)
        nullify(ppemfldana%pquad)
        nullify(ppemfldana%pdrift)
        nullify(ppemfldana%pccl)
        nullify(ppemfldana%pccdtl)
        nullify(ppemfldana%pdtl)
        nullify(ppemfldana%psc)
        nullify(ppemfldana%pcf)
        nullify(ppemfldana%pslrf)
        nullify(ppemfldana%pdipole)
        nullify(ppemfldana%pemfld)
        nullify(ppemfldana%pemfldcart)
        nullify(ppemfldana%pemfldcyl)
        nullify(ppemfldana%pmult)

        end function assign_emfldana

        function assign_emfldcart(temfldcart) result(ppemfldcart)
        type (BeamLineElem) :: ppemfldcart
        type (EMfldcart), target, intent(in) :: temfldcart

        ppemfldcart%pemfldcart => temfldcart
        nullify(ppemfldcart%psl)
        nullify(ppemfldcart%pbpm)
        nullify(ppemfldcart%pquad)
        nullify(ppemfldcart%pdrift)
        nullify(ppemfldcart%pccl)
        nullify(ppemfldcart%pccdtl)
        nullify(ppemfldcart%pdtl)
        nullify(ppemfldcart%psc)
        nullify(ppemfldcart%pcf)
        nullify(ppemfldcart%pslrf)
        nullify(ppemfldcart%pdipole)
        nullify(ppemfldcart%pemfld)
        nullify(ppemfldcart%pemfldana)
        nullify(ppemfldcart%pemfldcyl)
        nullify(ppemfldcart%pmult)

        end function assign_emfldcart

        function assign_emfldcyl(temfldcyl) result(ppemfldcyl)
        type (BeamLineElem) :: ppemfldcyl
        type (EMfldcyl), target, intent(in) :: temfldcyl

        ppemfldcyl%pemfldcyl => temfldcyl
        nullify(ppemfldcyl%psl)
        nullify(ppemfldcyl%pbpm)
        nullify(ppemfldcyl%pquad)
        nullify(ppemfldcyl%pdrift)
        nullify(ppemfldcyl%pccl)
        nullify(ppemfldcyl%pccdtl)
        nullify(ppemfldcyl%pdtl)
        nullify(ppemfldcyl%psc)
        nullify(ppemfldcyl%pcf)
        nullify(ppemfldcyl%pslrf)
        nullify(ppemfldcyl%pdipole)
        nullify(ppemfldcyl%pemfld)
        nullify(ppemfldcyl%pemfldana)
        nullify(ppemfldcyl%pemfldcart)
        nullify(ppemfldcyl%pmult)

        end function assign_emfldcyl

        function assign_mult(tmult) result(ppmult)
        type (BeamLineElem) :: ppmult
        type (Multipole), target, intent(in) :: tmult

        ppmult%pmult => tmult
        nullify(ppmult%psl)
        nullify(ppmult%pbpm)
        nullify(ppmult%pquad)
        nullify(ppmult%pdrift)
        nullify(ppmult%pccl)
        nullify(ppmult%pccdtl)
        nullify(ppmult%pdtl)
        nullify(ppmult%psc)
        nullify(ppmult%pcf)
        nullify(ppmult%pslrf)
        nullify(ppmult%pdipole)
        nullify(ppmult%pemfld)
        nullify(ppmult%pemfldana)
        nullify(ppmult%pemfldcart)
        nullify(ppmult%pemfldcyl)

        end function assign_mult

        subroutine getparam1_BeamLineElem(this,i,blparam)
        implicit none 
        type (BeamLineElem), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        if(associated(this%pquad)) then
          call getparam_Quadrupole(this%pquad,i,blparam)
        elseif(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,i,blparam)
        elseif(associated(this%pccl)) then
          call getparam_CCL(this%pccl,i,blparam)
        elseif(associated(this%pccdtl)) then
          call getparam_CCDTL(this%pccdtl,i,blparam)
        elseif(associated(this%pdtl)) then
          call getparam_DTL(this%pdtl,i,blparam)
        elseif(associated(this%psc)) then
          call getparam_SC(this%psc,i,blparam)
        elseif(associated(this%pbpm)) then
          call getparam_BPM(this%pbpm,i,blparam)
        elseif(associated(this%pcf)) then
          call getparam_ConstFoc(this%pcf,i,blparam)
        elseif(associated(this%pslrf)) then
          call getparam_SolRF(this%pslrf,i,blparam)
        elseif(associated(this%psl)) then
          call getparam_Sol(this%psl,i,blparam)
        elseif(associated(this%pdipole)) then
          call getparam_Dipole(this%pdipole,i,blparam)
        elseif(associated(this%pemfld)) then
          call getparam_EMfld(this%pemfld,i,blparam)
        elseif(associated(this%pemfldana)) then
          call getparam_EMfldAna(this%pemfldana,i,blparam)
        elseif(associated(this%pemfldcart)) then
          call getparam_EMfldCart(this%pemfldcart,i,blparam)
        elseif(associated(this%pemfldcyl)) then
          call getparam_EMfldCyl(this%pemfldcyl,i,blparam)
        elseif(associated(this%pmult)) then
          call getparam_Multipole(this%pmult,i,blparam)
        endif

        end subroutine getparam1_BeamLineElem
  
        subroutine getparam2_BeamLineElem(this,blparams)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        if(associated(this%pquad)) then
          call getparam_Quadrupole(this%pquad,blparams)
        elseif(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,blparams)
        elseif(associated(this%pccl)) then
          call getparam_CCL(this%pccl,blparams)
        elseif(associated(this%pccdtl)) then
          call getparam_CCDTL(this%pccdtl,blparams)
        elseif(associated(this%pdtl)) then
          call getparam_DTL(this%pdtl,blparams)
        elseif(associated(this%psc)) then
          call getparam_SC(this%psc,blparams)
        elseif(associated(this%pbpm)) then
          call getparam_BPM(this%pbpm,blparams)
        elseif(associated(this%pcf)) then
          call getparam_ConstFoc(this%pcf,blparams)
        elseif(associated(this%pslrf)) then
          call getparam_SolRF(this%pslrf,blparams)
        elseif(associated(this%psl)) then
          call getparam_Sol(this%psl,blparams)
        elseif(associated(this%pdipole)) then
          call getparam_Dipole(this%pdipole,blparams)
        elseif(associated(this%pemfld)) then
          call getparam_EMfld(this%pemfld,blparams)
        elseif(associated(this%pemfldana)) then
          call getparam_EMfldAna(this%pemfldana,blparams)
        elseif(associated(this%pemfldcart)) then
          call getparam_EMfldCart(this%pemfldcart,blparams)
        elseif(associated(this%pemfldcyl)) then
          call getparam_EMfldCyl(this%pemfldcyl,blparams)
        elseif(associated(this%pmult)) then
          call getparam_Multipole(this%pmult,blparams)
        endif

        end subroutine getparam2_BeamLineElem

        subroutine getparam3_BeamLineElem(this,blength,bnseg,bmapstp,&
                                          btype)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        if(associated(this%pquad)) then
          call getparam_Quadrupole(this%pquad,blength,bnseg,bmapstp,&
                                    btype)
        elseif(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,blength,bnseg,bmapstp,&
                                  btype)
        elseif(associated(this%pccl)) then
          call getparam_CCL(this%pccl,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pccdtl)) then
          call getparam_CCDTL(this%pccdtl,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pdtl)) then
          call getparam_DTL(this%pdtl,blength,bnseg,bmapstp,btype)
        elseif(associated(this%psc)) then
          call getparam_SC(this%psc,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pbpm)) then
          call getparam_BPM(this%pbpm,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pcf)) then
          call getparam_ConstFoc(this%pcf,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pslrf)) then
          call getparam_SolRF(this%pslrf,blength,bnseg,bmapstp,btype)
        elseif(associated(this%psl)) then
          call getparam_Sol(this%psl,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pdipole)) then
          call getparam_Dipole(this%pdipole,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pemfld)) then
          call getparam_EMfld(this%pemfld,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pemfldana)) then
          call getparam_EMfldAna(this%pemfldana,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pemfldcart)) then
          call getparam_EMfldCart(this%pemfldcart,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pemfldcyl)) then
          call getparam_EMfldCyl(this%pemfldcyl,blength,bnseg,bmapstp,btype)
        elseif(associated(this%pmult)) then
          call getparam_Multipole(this%pmult,blength,bnseg,bmapstp,btype)
        endif

        end subroutine getparam3_BeamLineElem
       
        subroutine getradius_BeamLineElem(this,piperadius)
        implicit none 
        type (BeamLineElem), intent(in) :: this
        double precision, intent(out) :: piperadius

        if(associated(this%pquad)) then
          call getparam_Quadrupole(this%pquad,4,piperadius)
        elseif(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,2,piperadius)
        elseif(associated(this%pccl)) then
          call getparam_CCL(this%pccl,6,piperadius)
        elseif(associated(this%pccdtl)) then
          call getparam_CCDTL(this%pccdtl,6,piperadius)
        elseif(associated(this%pdtl)) then
          call getparam_DTL(this%pdtl,6,piperadius)
        elseif(associated(this%psc)) then
          call getparam_SC(this%psc,6,piperadius)
        elseif(associated(this%pbpm)) then
          call getparam_BPM(this%pbpm,2,piperadius)
        elseif(associated(this%pcf)) then
          call getparam_ConstFoc(this%pcf,5,piperadius)
        elseif(associated(this%pslrf)) then
          call getparam_SolRF(this%pslrf,6,piperadius)
        elseif(associated(this%psl)) then
          call getparam_Sol(this%psl,4,piperadius)
        elseif(associated(this%pdipole)) then
          call getparam_Dipole(this%pdipole,5,piperadius)
        elseif(associated(this%pemfld)) then
          call getparam_EMfld(this%pemfld,6,piperadius)
        elseif(associated(this%pemfldana)) then
          call getparam_EMfldAna(this%pemfldana,6,piperadius)
        elseif(associated(this%pemfldcart)) then
          call getparam_EMfldCart(this%pemfldcart,6,piperadius)
        elseif(associated(this%pemfldcyl)) then
          call getparam_EMfldCyl(this%pemfldcyl,6,piperadius)
        elseif(associated(this%pmult)) then
          call getparam_Multipole(this%pmult,5,piperadius)
        endif

        end subroutine getradius_BeamLineElem
  
        subroutine geterr_BeamLineElem(this,xerr,yerr,anglerrx,anglerry,&
                                       anglerrz)
        implicit none 
        type (BeamLineElem), intent(in) :: this
        double precision, intent(out) :: xerr,yerr,anglerrx,anglerry,anglerrz

        if(associated(this%pquad)) then
          xerr = this%pquad%Param(5)
          yerr = this%pquad%Param(6)
          anglerrx = this%pquad%Param(7)
          anglerry = this%pquad%Param(8)
          anglerrz = this%pquad%Param(9)
        elseif(associated(this%pdrift)) then
          xerr = 0.0
          yerr = 0.0
          anglerrx = 0.0
          anglerry = 0.0
          anglerrz = 0.0
        elseif(associated(this%pccl)) then
          xerr = this%pccl%Param(7)
          yerr = this%pccl%Param(8)
          anglerrx = this%pccl%Param(9)
          anglerry = this%pccl%Param(10)
          anglerrz = this%pccl%Param(11)
        elseif(associated(this%pccdtl)) then
          xerr = this%pccdtl%Param(7)
          yerr = this%pccdtl%Param(8)
          anglerrx = this%pccdtl%Param(9)
          anglerry = this%pccdtl%Param(10)
          anglerrz = this%pccdtl%Param(11)
        elseif(associated(this%pdtl)) then
          xerr = this%pdtl%Param(11)
          yerr = this%pdtl%Param(12)
          anglerrx = this%pdtl%Param(13)
          anglerry = this%pdtl%Param(14)
          anglerrz = this%pdtl%Param(15)
        elseif(associated(this%psc)) then
          xerr = this%psc%Param(7)
          yerr = this%psc%Param(8)
          anglerrx = this%psc%Param(9)
          anglerry = this%psc%Param(10)
          anglerrz = this%psc%Param(11)
        elseif(associated(this%pbpm)) then
          xerr = 0.0
          yerr = 0.0
          anglerrx = 0.0
          anglerry = 0.0
          anglerrz = 0.0
        elseif(associated(this%pcf)) then
          xerr = 0.0
          yerr = 0.0
          anglerrx = 0.0
          anglerry = 0.0
          anglerrz = 0.0
        elseif(associated(this%pslrf)) then
          xerr = this%pslrf%Param(7)
          yerr = this%pslrf%Param(8)
          anglerrx = this%pslrf%Param(9)
          anglerry = this%pslrf%Param(10)
          anglerrz = this%pslrf%Param(11)
        elseif(associated(this%psl)) then
          xerr = this%psl%Param(5)
          yerr = this%psl%Param(6)
          anglerrx = this%psl%Param(7)
          anglerry = this%psl%Param(8)
          anglerrz = this%psl%Param(9)
        elseif(associated(this%pdipole)) then
          xerr = this%pdipole%Param(6)
          yerr = this%pdipole%Param(7)
          anglerrx = this%pdipole%Param(8)
          anglerry = this%pdipole%Param(9)
          anglerrz = this%pdipole%Param(10)
        elseif(associated(this%pemfld)) then
          xerr = this%pemfld%Param(7)
          yerr = this%pemfld%Param(8)
          anglerrx = this%pemfld%Param(9)
          anglerry = this%pemfld%Param(10)
          anglerrz = this%pemfld%Param(11)
        elseif(associated(this%pemfldana)) then
          xerr = this%pemfldana%Param(7)
          yerr = this%pemfldana%Param(8)
          anglerrx = this%pemfldana%Param(9)
          anglerry = this%pemfldana%Param(10)
          anglerrz = this%pemfldana%Param(11)
        elseif(associated(this%pemfldcart)) then
          xerr = this%pemfldcart%Param(7)
          yerr = this%pemfldcart%Param(8)
          anglerrx = this%pemfldcart%Param(9)
          anglerry = this%pemfldcart%Param(10)
          anglerrz = this%pemfldcart%Param(11)
        elseif(associated(this%pemfldcyl)) then
          xerr = this%pemfldcyl%Param(7)
          yerr = this%pemfldcyl%Param(8)
          anglerrx = this%pemfldcyl%Param(9)
          anglerry = this%pemfldcyl%Param(10)
          anglerrz = this%pemfldcyl%Param(11)
        elseif(associated(this%pmult)) then
          xerr = this%pmult%Param(6)
          yerr = this%pmult%Param(7)
          anglerrx = this%pmult%Param(8)
          anglerry = this%pmult%Param(9)
          anglerrz = this%pmult%Param(10)
        endif

        end subroutine geterr_BeamLineElem
  
        subroutine setparam1_BeamLineElem(this,i,blparam)
        implicit none 
        type (BeamLineElem), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: blparam

        if(associated(this%pquad)) then
          call setparam_Quadrupole(this%pquad,i,blparam)
        elseif(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,i,blparam)
        elseif(associated(this%pccl)) then
          call setparam_CCL(this%pccl,i,blparam)
        elseif(associated(this%pccdtl)) then
          call setparam_CCDTL(this%pccdtl,i,blparam)
        elseif(associated(this%pdtl)) then
          call setparam_DTL(this%pdtl,i,blparam)
        elseif(associated(this%psc)) then
          call setparam_SC(this%psc,i,blparam)
        elseif(associated(this%pbpm)) then
          call setparam_BPM(this%pbpm,i,blparam)
        elseif(associated(this%pcf)) then
          call setparam_ConstFoc(this%pcf,i,blparam)
        elseif(associated(this%pslrf)) then
          call setparam_SolRF(this%pslrf,i,blparam)
        elseif(associated(this%psl)) then
          call setparam_Sol(this%psl,i,blparam)
        elseif(associated(this%pdipole)) then
          call setparam_Dipole(this%pdipole,i,blparam)
        elseif(associated(this%pemfld)) then
          call setparam_EMfld(this%pemfld,i,blparam)
        elseif(associated(this%pemfldana)) then
          call setparam_EMfldAna(this%pemfldana,i,blparam)
        elseif(associated(this%pemfldcart)) then
          call setparam_EMfldCart(this%pemfldcart,i,blparam)
        elseif(associated(this%pemfldcyl)) then
          call setparam_EMfldCyl(this%pemfldcyl,i,blparam)
        elseif(associated(this%pmult)) then
          call setparam_Multipole(this%pmult,i,blparam)
        endif

        end subroutine setparam1_BeamLineElem
  
        subroutine setparam2_BeamLineElem(this,blparams)
        implicit none
        type (BeamLineElem), intent(inout) :: this
        double precision, dimension(:), intent(in) :: blparams

        if(associated(this%pquad)) then
          call setparam_Quadrupole(this%pquad,blparams)
        elseif(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,blparams)
        elseif(associated(this%pccl)) then
          call setparam_CCL(this%pccl,blparams)
        elseif(associated(this%pccdtl)) then
          call setparam_CCDTL(this%pccdtl,blparams)
        elseif(associated(this%pdtl)) then
          call setparam_DTL(this%pdtl,blparams)
        elseif(associated(this%psc)) then
          call setparam_SC(this%psc,blparams)
        elseif(associated(this%pbpm)) then
          call setparam_BPM(this%pbpm,blparams)
        elseif(associated(this%pcf)) then
          call setparam_ConstFoc(this%pcf,blparams)
        elseif(associated(this%pslrf)) then
          call setparam_SolRF(this%pslrf,blparams)
        elseif(associated(this%psl)) then
          call setparam_Sol(this%psl,blparams)
        elseif(associated(this%pdipole)) then
          call setparam_Dipole(this%pdipole,blparams)
        elseif(associated(this%pemfld)) then
          call setparam_EMfld(this%pemfld,blparams)
        elseif(associated(this%pemfldana)) then
          call setparam_EMfldAna(this%pemfldana,blparams)
        elseif(associated(this%pemfldcart)) then
          call setparam_EMfldCart(this%pemfldcart,blparams)
        elseif(associated(this%pemfldcyl)) then
          call setparam_EMfldCyl(this%pemfldcyl,blparams)
        elseif(associated(this%pmult)) then
          call setparam_Multipole(this%pmult,blparams)
        endif

        end subroutine setparam2_BeamLineElem

        subroutine setparam3_BeamLineElem(this,blength,bnseg,bmapstp,&
                                          btype)
        implicit none
        type (BeamLineElem), intent(inout) :: this
        double precision, intent(in) :: blength
        integer, intent(in) :: bnseg,bmapstp,btype

        if(associated(this%pquad)) then
          call setparam_Quadrupole(this%pquad,bnseg,bmapstp,&
                                    btype,blength)
        elseif(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,bnseg,bmapstp,&
                                  btype,blength)
        elseif(associated(this%pccl)) then
          call setparam_CCL(this%pccl,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pccdtl)) then
          call setparam_CCDTL(this%pccdtl,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pdtl)) then
          call setparam_DTL(this%pdtl,bnseg,bmapstp,btype,blength)
        elseif(associated(this%psc)) then
          call setparam_SC(this%psc,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pbpm)) then
          call setparam_BPM(this%pbpm,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pcf)) then
          call setparam_ConstFoc(this%pcf,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pslrf)) then
          call setparam_SolRF(this%pslrf,bnseg,bmapstp,btype,blength)
        elseif(associated(this%psl)) then
          call setparam_Sol(this%psl,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pdipole)) then
          call setparam_Dipole(this%pdipole,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pemfld)) then
          call setparam_EMfld(this%pemfld,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pemfldana)) then
          call setparam_EMfldAna(this%pemfldana,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pemfldcart)) then
          call setparam_EMfldCart(this%pemfldcart,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pemfldcyl)) then
          call setparam_EMfldCyl(this%pemfldcyl,bnseg,bmapstp,btype,blength)
        elseif(associated(this%pmult)) then
          call setparam_Multipole(this%pmult,bnseg,bmapstp,&
                                    btype,blength)
        endif

        end subroutine setparam3_BeamLineElem
       
        subroutine getfld_BeamLineElem(this,pos,extfld)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(6), intent(out) :: extfld

        if(associated(this%pquad)) then
          call getfld_Quadrupole(pos,extfld,this%pquad)
        elseif(associated(this%pdrift)) then
          call getfld_DriftTube(pos,extfld,this%pdrift)
        elseif(associated(this%pccl)) then
          call getfld_CCL(pos,extfld,this%pccl)
        elseif(associated(this%pccdtl)) then
          call getfld_CCDTL(pos,extfld,this%pccdtl)
        elseif(associated(this%pdtl)) then
          call getfld_DTL(pos,extfld,this%pdtl)
        elseif(associated(this%psc)) then
          call getfld_SC(pos,extfld,this%psc)
        elseif(associated(this%pbpm)) then
          !call getfld_BPM(pos,extfld,this%pbpm)
          print*,"no field for BPM!!"
          extfld = 0.0
        elseif(associated(this%pcf)) then
          call getfld_ConstFoc(pos,extfld,this%pcf)
        elseif(associated(this%pslrf)) then
          call getfld_SolRF(pos,extfld,this%pslrf)
        elseif(associated(this%psl)) then
          call getfld_Sol(pos,extfld,this%psl)
        elseif(associated(this%pdipole)) then
          call getfld_Dipole(pos,extfld,this%pdipole)
        elseif(associated(this%pemfld)) then
          call getfld_EMfld(pos,extfld,this%pemfld)
        elseif(associated(this%pemfldana)) then
          call getfld_EMfldAna(pos,extfld,this%pemfldana)
        elseif(associated(this%pemfldcart)) then
          call getfld_EMfldCart(pos,extfld,this%pemfldcart)
        elseif(associated(this%pemfldcyl)) then
          call getfld_EMfldCyl(pos,extfld,this%pemfldcyl)
        elseif(associated(this%pmult)) then
          call getfld_Multipole(pos,extfld,this%pmult)
        endif

        end subroutine getfld_BeamLineElem

        !--------------------------------------------------------------------------------------
        !> @brief
        !> get external field with displacement and rotation errors.
        !--------------------------------------------------------------------------------------
        subroutine getflderr_BeamLineElem(this,pos,extfld,dx,dy,anglex,&
                                          angley,anglez)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(6), intent(out) :: extfld

        if(associated(this%pquad)) then
          call getflderr_Quadrupole(pos,extfld,this%pquad,dx,dy,anglex,&
                                    angley,anglez)
        elseif(associated(this%pdrift)) then
          call getfld_DriftTube(pos,extfld,this%pdrift)
        elseif(associated(this%pccl)) then
          call getflderr_CCL(pos,extfld,this%pccl,dx,dy,anglex,&
                                    angley,anglez)
        elseif(associated(this%pccdtl)) then
          call getflderr_CCDTL(pos,extfld,this%pccdtl,dx,dy,anglex,&
                                    angley,anglez)
        elseif(associated(this%pdtl)) then
          call getflderr_DTL(pos,extfld,this%pdtl,dx,dy,anglex,&
                                    angley,anglez)
        elseif(associated(this%psc)) then
          call getflderr_SC(pos,extfld,this%psc,dx,dy,anglex,&
                                    angley,anglez)
        elseif(associated(this%pbpm)) then
          !call getfld_BPM(pos,extfld,this%pbpm)
          print*,"no field for BPM!!"
          extfld = 0.0
        elseif(associated(this%pcf)) then
          call getfld_ConstFoc(pos,extfld,this%pcf)
        elseif(associated(this%pslrf)) then
          call getflderr_SolRF(pos,extfld,this%pslrf,dx,dy,anglex,&
                               angley,anglez)
        elseif(associated(this%psl)) then
          call getflderr_Sol(pos,extfld,this%psl,dx,dy,anglex,&
                               angley,anglez)
        elseif(associated(this%pdipole)) then
          call getflderr_Dipole(pos,extfld,this%pdipole,dx,dy,anglex,&
                               angley,anglez)
        elseif(associated(this%pemfld)) then
          call getflderr_EMfld(pos,extfld,this%pemfld,dx,dy,anglex,&
                               angley,anglez)
        elseif(associated(this%pemfldana)) then
          call getflderr_EMfldAna(pos,extfld,this%pemfldana,dx,dy,anglex,&
                               angley,anglez)
        elseif(associated(this%pemfldcart)) then
          call getflderr_EMfldCart(pos,extfld,this%pemfldcart,dx,dy,anglex,&
                               angley,anglez)
        elseif(associated(this%pemfldcyl)) then
          call getflderr_EMfldCyl(pos,extfld,this%pemfldcyl,dx,dy,anglex,&
                               angley,anglez)
        elseif(associated(this%pmult)) then
          call getflderr_Multipole(pos,extfld,this%pmult,dx,dy,anglex,&
                                    angley,anglez)
        endif

        end subroutine getflderr_BeamLineElem

        subroutine getaxfldE_BeamLineElem(this,z,ez1,ezp1,ezpp1)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1

        if(associated(this%pquad)) then
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        elseif(associated(this%pdrift)) then
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        elseif(associated(this%pccl)) then
          call getaxfldE_CCL(z,this%pccl,ez1,ezp1,ezpp1)
        elseif(associated(this%pccdtl)) then
          call getaxfldE_CCDTL(z,this%pccdtl,ez1,ezp1,ezpp1)
        elseif(associated(this%pdtl)) then
          call getaxfldE_DTL(z,this%pdtl,ez1,ezp1,ezpp1)
        elseif(associated(this%psc)) then
          call getaxfldE_SC(z,this%psc,ez1,ezp1,ezpp1)
        elseif(associated(this%pbpm)) then
          print*,"no field in BPM!!"
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        elseif(associated(this%pcf)) then
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        elseif(associated(this%pslrf)) then
          call getaxfldE_SolRF(z,this%pslrf,ez1,ezp1,ezpp1)
        elseif(associated(this%psl)) then
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        elseif(associated(this%pdipole)) then
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        elseif(associated(this%pemfld)) then
          call getaxfldE_EMfld(z,this%pemfld,ez1,ezp1,ezpp1)
        elseif(associated(this%pmult)) then
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        endif

        end subroutine getaxfldE_BeamLineElem

        subroutine getfldt_BeamLineElem(this,pos,extfld,fldata)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(6), intent(out) :: extfld
        type (fielddata), intent(in) :: fldata

        if(associated(this%pquad)) then
          call getfld_Quadrupole(pos,extfld,this%pquad)
        elseif(associated(this%pdrift)) then
          call getfld_DriftTube(pos,extfld,this%pdrift)
        elseif(associated(this%pccl)) then
          call getfldt_CCL(pos,extfld,this%pccl,fldata)
        elseif(associated(this%pccdtl)) then
          call getfldt_CCDTL(pos,extfld,this%pccdtl,fldata)
        elseif(associated(this%pdtl)) then
          call getfldt_DTL(pos,extfld,this%pdtl,fldata)
        elseif(associated(this%psc)) then
          call getfldt_SC(pos,extfld,this%psc,fldata)
        elseif(associated(this%pbpm)) then
          !call getfld_BPM(pos,extfld,this%pbpm)
          print*,"no field for BPM!!"
          extfld = 0.0
        elseif(associated(this%pcf)) then
          call getfld_ConstFoc(pos,extfld,this%pcf)
        elseif(associated(this%pslrf)) then
          call getfldt_SolRF(pos,extfld,this%pslrf,fldata)
        elseif(associated(this%psl)) then
          call getfldt_Sol(pos,extfld,this%psl,fldata)
        elseif(associated(this%pdipole)) then
          call getfldt_Dipole(pos,extfld,this%pdipole,fldata)
        elseif(associated(this%pemfld)) then
          print*,"wrong element type..."
          stop
        elseif(associated(this%pemfldana)) then
          call getfldt_EMfldAna(pos,extfld,this%pemfldana,fldata)
        elseif(associated(this%pemfldcart)) then
          call getfldt_EMfldCart(pos,extfld,this%pemfldcart,fldata)
        elseif(associated(this%pemfldcyl)) then
          call getfldt_EMfldCyl(pos,extfld,this%pemfldcyl,fldata)
        elseif(associated(this%pmult)) then
          call getfld_Multipole(pos,extfld,this%pmult)
        endif

        end subroutine getfldt_BeamLineElem

        subroutine getflderrt_BeamLineElem(this,pos,extfld,fldata)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(6), intent(out) :: extfld
        type (fielddata), intent(in) :: fldata

        if(associated(this%pquad)) then
          call getflderrt_Quadrupole(pos,extfld,this%pquad)
        elseif(associated(this%pdrift)) then
          call getfld_DriftTube(pos,extfld,this%pdrift)
        elseif(associated(this%pccl)) then
          call getfldt_CCL(pos,extfld,this%pccl,fldata)
        elseif(associated(this%pccdtl)) then
          call getfldt_CCDTL(pos,extfld,this%pccdtl,fldata)
        elseif(associated(this%pdtl)) then
          call getfldt_DTL(pos,extfld,this%pdtl,fldata)
        elseif(associated(this%psc)) then
          call getfldt_SC(pos,extfld,this%psc,fldata)
        elseif(associated(this%pbpm)) then
          !call getfld_BPM(pos,extfld,this%pbpm)
          print*,"no field for BPM!!"
          extfld = 0.0
        elseif(associated(this%pcf)) then
          call getfld_ConstFoc(pos,extfld,this%pcf)
        elseif(associated(this%pslrf)) then
          call getflderrt_SolRF(pos,extfld,this%pslrf,fldata)
        elseif(associated(this%psl)) then
          call getfldt_Sol(pos,extfld,this%psl,fldata)
        elseif(associated(this%pdipole)) then
          call getfldt_Dipole(pos,extfld,this%pdipole,fldata)
        elseif(associated(this%pemfld)) then
          print*,"wrong element type..."
          stop
        elseif(associated(this%pemfldana)) then
          call getfldt_EMfldAna(pos,extfld,this%pemfldana,fldata)
        elseif(associated(this%pemfldcart)) then
          call getfldt_EMfldCart(pos,extfld,this%pemfldcart,fldata)
        elseif(associated(this%pemfldcyl)) then
          call getfldt_EMfldCyl(pos,extfld,this%pemfldcyl,fldata)
        elseif(associated(this%pmult)) then
          call getflderrt_Multipole(pos,extfld,this%pmult)
        endif

        end subroutine getflderrt_BeamLineElem
      end module BeamLineElemclass
