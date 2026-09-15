!
!      Open Boundary Condition Module
!
module openbc_mod

#include <def-undef.h>

use precision_mod
use pconst_mod, only : IMT, JMT, KM, imt_global, jmt_global

!
!       The Given Value ( 1D Array ) Outside of the Open Boundary For 2-D External Barotropic Mode
!       They are given in inirun.F90 by using the given initial values
!

        REAL(r8), DIMENSION(2,JMT) :: UBTRE, UBTRW, VBTRE, VBTRW
        REAL(r8), DIMENSION(IMT,2) :: VBTRS, VBTRN, UBTRS, UBTRN
        REAL(r8), DIMENSION(2,JMT) :: H0BTRE, H0BTRW
        REAL(r8), DIMENSION(IMT,2) :: H0BTRS, H0BTRN
!
!
!       The Given Value ( 2D Array ) Outside of the Open Boundary For 3-D internal Baroclinic Mode
!       They are given in inirun.F90 by using the given intial values
!

        REAL(r8), DIMENSION(JMT,KM) :: TBCLE, SBCLE
        REAL(r8), DIMENSION(JMT,KM) :: TBCLW, SBCLW
        REAL(r8), DIMENSION(IMT,KM) :: TBCLS, SBCLS
        REAL(r8), DIMENSION(IMT,KM) :: TBCLN, SBCLN

        REAL(r8), DIMENSION(2,JMT,KM) :: UBCLE, VBCLE
        REAL(r8), DIMENSION(2,JMT,KM) :: UBCLW, VBCLW
        REAL(r8), DIMENSION(IMT,2,KM) :: UBCLS, VBCLS
        REAL(r8), DIMENSION(IMT,2,KM) :: UBCLN, VBCLN

        INTEGER, PARAMETER :: imt_global_gl=3602
        INTEGER, PARAMETER :: jmt_global_gl=1683

        REAL(r8), DIMENSION(IMT,JMT) :: H0_MON_LOC
        REAL(r8), DIMENSION(IMT,JMT) :: UB_MON_LOC
        REAL(r8), DIMENSION(IMT,JMT) :: VB_MON_LOC

        REAL(r8), DIMENSION(IMT,JMT,KM) :: TT_MON_LOC, SS_MON_LOC
        REAL(r8), DIMENSION(IMT,JMT,KM) :: UU_MON_LOC, VV_MON_LOC

        REAL(r8) :: STARTLON=95.0D0
        REAL(r8) :: STARTLAT=48.0D0

        INTEGER ::  BC_COUNT
        REAL(r8) :: BC_CRSS
        REAL(r8) :: BC_FLUX
        REAL(r8) :: UBAR_XS

        INTEGER, PARAMETER :: ID_OBC_E = 1
        INTEGER, PARAMETER :: ID_OBC_W = 2
        INTEGER, PARAMETER :: ID_OBC_S = 3
        INTEGER, PARAMETER :: ID_OBC_N = 4

!        REAL(r8) :: TAUOUT = 86400.0D0*365.0D0
!        REAL(r8) :: TAUIN  = 86400.0D0*5.0D0

!        REAL(r8) :: TAUOUT_TS = 86400.0D0*365.0D0
!        REAL(r8) :: TAUIN_TS = 86400.0D0*5.0D0

!        REAL(r8) :: TAUOUT_H0 = 86400.0D0*365.0D0
!        REAL(r8) :: TAUIN_H0 = 86400.0D0*5.0D0
        
        REAL(r8) :: TAUOUT = 86400.0D0*30.0D0
        REAL(r8) :: TAUIN  = 86400.0D0*1.0D0

        REAL(r8) :: TAUOUT_TS = 86400.0D0*15.0D0
        REAL(r8) :: TAUIN_TS = 86400.0D0*1.0D0

        REAL(r8) :: TAUOUT_H0 = 86400.0D0*15.0D0
        REAL(r8) :: TAUIN_H0 = 86400.0D0*1.0D0

!--------------------------------------------------------------------------------------------------
!add by zhangjk 20241127
        real(r8) :: tau_in_c5  = 1.0d0 / (5.0d0 * 86400.0d0)
        real(r8) :: tau_out_c5 = 1.0d0 / (25.0d0 * 86400.0d0)
        REAL(R8) :: TAU_2DUV = 86400.0D0*1.0D0
        REAL(R8) :: TAU_2DH0 = 86400.0D0*(1.0d0/24.0D0)

!add by zhangjk 20250311
        real(r8) :: ub_correct

!add by zhangjk 0250314
        real(r8) :: dtb_openbc, dtc_openbc
        real(r8), dimension(imt, jmt) :: ubpp, vbpp, h0pp
        real(r8), dimension(imt, jmt, km) :: unow, vnow 
!--------------------------------------------------------------------------------------------------

#if (defined OBCSPONGE )
        REAL(r8) :: width_sponge=15.0D0
#endif

#if (defined UV_NUDGE )
        REAL(r8) :: width_uvnudge=30.0D0
        REAL(r8), DIMENSION(IMT,JMT) :: UVnudgcof
#else 
        REAL(r8) :: width_uvnudge=30.0D0
        REAL(r8), DIMENSION(IMT,JMT) :: UVnudgcof
#endif

#if (defined TS_NUDGE )
        REAL(r8) :: width_tsnudge=30.0D0
        REAL(r8), DIMENSION(IMT,JMT) :: TSnudgcof 
#else
        REAL(r8) :: width_tsnudge=30.0D0
        REAL(r8), DIMENSION(IMT,JMT) :: TSnudgcof 
#endif

!
!     ID_OBC_H0_X = 1     Zero Gradient BC
!                 = 2     Clamped BC
!                 = 3     Implicit Chapman BC
!                 = 4     Explicit Chapman BC
!                 = 5     Implicit Upstream Radiation BC
!                 = 6     Flather BC
!                 = 7     recode Implicit Chapman BC (add by zhangjk 20250221)
!                 = 8     recode Explicit Chapman BC (add by zhangjk 20250311)
!                 = 9     use work in the Implicit Chapman BC (add by zhangjk 20250313)
!                 = 0     No BC
!
        INTEGER, PARAMETER :: ID_OBC_H0_E = 1
        INTEGER, PARAMETER :: ID_OBC_H0_W = 1
        INTEGER, PARAMETER :: ID_OBC_H0_S = 1
        INTEGER, PARAMETER :: ID_OBC_H0_N = 1

!
!     ID_OBC_UB_E/W or ID_OBC_VB_S/N  = 1     Zero Gradient BC 
!                                     = 2     Clamped BC 
!                                     = 3     Upstream Advection BC 
!                                     = 4     Implicit Chapman BC
!                                     = 5     Explicit Chapman BC
!                                     = 6     Implicit Upstream Radiation BC
!                                     = 7     Flather BC
!                                     = 8     The Modified Flather BC by Shchepetkin (Maison et al., 2010)
!                                     = 9     recode The Modified Flather BC (add by zhangjk 20250216)
!                                     = 10    recode Flather BC (add by zhangjk 20250223)
!                                     = 11    use wka in the Modified Flather BC (add by zhangjk 20250313)
!                                     = 0     No BC
!
        INTEGER, PARAMETER :: ID_OBC_UB_E = 7
        INTEGER, PARAMETER :: ID_OBC_UB_W = 7

        INTEGER, PARAMETER :: ID_OBC_VB_S = 7
        INTEGER, PARAMETER :: ID_OBC_VB_N = 7

!
!     ID_OBC_UB_S/N or ID_OBC_VB_E/W  = 1     Zero Gradient BC 
!                                     = 2     Clamped BC 
!                                     = 3     Upstream Advection BC 
!                                     = 4     Implicit Chapman BC
!                                     = 5     Explicit Chapman BC
!                                     = 6     Implicit Upstream Radiation BC
!                                     = 7     recode Implicit Chapman BC (add by zhangjk 20250216)
!                                     = 8     recode Explicit Chapman BC (add by zhangjk 20250311)
!                                     = 9     use wka in the Implicit Chapman BC (add by zhangjk 20250313)
!                                     = 0     No BC
!
        INTEGER, PARAMETER :: ID_OBC_UB_S = 3
        INTEGER, PARAMETER :: ID_OBC_UB_N = 3

        INTEGER, PARAMETER :: ID_OBC_VB_E = 3
        INTEGER, PARAMETER :: ID_OBC_VB_W = 3
!
!
!     ID_OBC_U/V_X = 1     Zero Gradient BC 
!                  = 2     Clamped BC 
!                  = 3     Upstream Advection BC 
!                  = 4     Implicit Upstream Radiation BC
!                  = 5     recode Implicit Upstream Radiation BC (add by zhangjk 20241227)
!                  = 0     No BC
!
!
        INTEGER, PARAMETER :: ID_OBC_U_E  = 5
        INTEGER, PARAMETER :: ID_OBC_U_W  = 5

        INTEGER, PARAMETER :: ID_OBC_V_S  = 5
        INTEGER, PARAMETER :: ID_OBC_V_N  = 5

        INTEGER, PARAMETER :: ID_OBC_U_S  = 3
        INTEGER, PARAMETER :: ID_OBC_U_N  = 3

        INTEGER, PARAMETER :: ID_OBC_V_E  = 3
        INTEGER, PARAMETER :: ID_OBC_V_W  = 3

!
!     ID_OBC_TS_X = 1     Zero Gradient BC 
!                 = 2     Clamped BC 
!                 = 3     Upstream Advection BC 
!                 = 4     Implicit Upstream Radiation BC
!                 = 5     recode Implicit Upstream Radiation BC (add by zhangjk 20250316)
!                 = 0     No BC
!
!
        INTEGER, PARAMETER :: ID_OBC_TS_E = 3
        INTEGER, PARAMETER :: ID_OBC_TS_W = 3
        INTEGER, PARAMETER :: ID_OBC_TS_S = 3
        INTEGER, PARAMETER :: ID_OBC_TS_N = 3
!
!
!       ID_OBC_2DRad_XX = 1  2D Implicit Upstream Radiation
!                         0  1D Radiation
!
        INTEGER, PARAMETER :: ID_OBC_2DRad_H0 = 1
        INTEGER, PARAMETER :: ID_OBC_2DRad_UB = 1
        INTEGER, PARAMETER :: ID_OBC_2DRad_VB = 1
        INTEGER, PARAMETER :: ID_OBC_2DRad_U = 1
        INTEGER, PARAMETER :: ID_OBC_2DRad_V = 1
        INTEGER, PARAMETER :: ID_OBC_2DRad_TS = 1
!
!       ID_OBC_2DRad_NPO_XX = 0 2D Implicit Upstream Radiation
!                           = 1  NPO 
!
        INTEGER, PARAMETER :: ID_OBC_2DRad_NPO_H0 = 0 
        INTEGER, PARAMETER :: ID_OBC_2DRad_NPO_UB = 0
        INTEGER, PARAMETER :: ID_OBC_2DRad_NPO_VB = 0
        INTEGER, PARAMETER :: ID_OBC_2DRad_NPO_U = 1 
        INTEGER, PARAMETER :: ID_OBC_2DRad_NPO_V = 1 
        INTEGER, PARAMETER :: ID_OBC_2DRad_NPO_TS = 1 
        integer, parameter :: ID_OBC_2DRad_NPO_UV5 = 1  ! 0 : don't use NPO in U/V openbc 5 ;
        integer, parameter :: ID_OBC_2DRad_NPO_TS5 = 1  ! 0 : don't use NPO in T/S openbc 5 ;

!
!       ID_OBC_VCON_UVB = 1  Volume constrain
!                         2  add by zhangjk 20250104
!                         0  No Volume constrain
!
        INTEGER, PARAMETER :: ID_OBC_VCON_UVB = 0    

!
!       ID_OBC_SET_MASS = 1 set         add by zhangjk 20250311
!                       = 0 no set
!
        INTEGER, PARAMETER :: ID_OBC_SET_MASS = 0

        INTEGER, PARAMETER :: ID_H0_ZEROGRADI_NUDGING = 1
        INTEGER, PARAMETER :: ID_UBVB_CHAPMAN_NUDGING = 0
        
        INTEGER, PARAMETER :: ID_OBC_SMO_3DUV_NM = 0
        INTEGER, PARAMETER :: ID_OBC_SMO_3DUV_TG = 0

end module openbc_mod
