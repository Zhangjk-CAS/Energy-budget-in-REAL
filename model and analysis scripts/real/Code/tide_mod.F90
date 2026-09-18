!
!      tide Condition Module
!
module tide_mod

#include <def-undef.h>

use precision_mod
use pconst_mod, only : IMT, JMT, imt_global, jmt_global

#ifdef OBCDT
use openbc_mod, only : STARTLAT, STARTLON
#endif

        INTEGER, PARAMETER :: Ntide=10
  
        REAL(r8) :: Tintg_tide
!
        REAL(r8), DIMENSION(Ntide) :: OMGATIDE

        REAL(r8), DIMENSION(Ntide) :: FNOD
        REAL(r8), DIMENSION(Ntide) :: UNOD
        REAL(r8), DIMENSION(Ntide) :: V0NOD

        REAL(r8), DIMENSION(imt_global,jmt_global) :: H0TIDE_TPXO7
        REAL(r8), DIMENSION(imt_global,jmt_global) :: THTIDE_TPXO7

        REAL(r8), DIMENSION(imt,jmt,Ntide) :: H0TIDE_TPXO7_LOC
        REAL(r8), DIMENSION(imt,jmt,Ntide) :: THTIDE_TPXO7_LOC

        REAL(r8), DIMENSION(IMT,Ntide) :: H0TIDES, H0TIDEN
        REAL(r8), DIMENSION(IMT,Ntide) :: THTIDES, THTIDEN

        REAL(r8), DIMENSION(JMT,Ntide) :: H0TIDEE, H0TIDEW
        REAL(r8), DIMENSION(JMT,Ntide) :: THTIDEE, THTIDEW

        REAL(r8), DIMENSION(IMT,JMT) :: ptide_elevat

        REAL(r8), PARAMETER :: BETA_TIDE = 1.052D0

        REAL(r8), PARAMETER :: CD_COASTAL = 0.001 
        REAL(r8), PARAMETER :: CD_INDIAN = 0.1
        REAL(r8), PARAMETER :: CD_ATLANTIC = 0.4
        REAL(r8), PARAMETER :: CD_PACIFIC = 0.7

        REAL(r8), DIMENSION(IMT,JMT) :: CD_TIDE

end module tide_mod
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
