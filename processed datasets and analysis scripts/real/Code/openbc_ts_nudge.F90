!
!

      SUBROUTINE openbc_ts_nudge

#include <def-undef.h>
use   precision_mod

use   param_mod,  only: IMT, JMT, KM, I, J, K, nx_proc, ny_proc
use   pconst_mod, only: VIT, ix, iy, DTS   
use   tracer_mod, only: AT, ATB

#if ( defined SPMD )
use msg_mod
#endif

#ifdef OBCDT
!use   openbc_mod, only: TT_MON_LOC, SS_MON_LOC, TAUOUT_TS, width_tsnudge, TSnudgcof
use   openbc_mod, only: TT_MON_LOC, SS_MON_LOC, TAUIN_TS, width_tsnudge, TSnudgcof
#endif

#ifdef TS_NUDGE

      IMPLICIT NONE

!      INTEGER :: Iwest
!      INTEGER :: Ieast
!      INTEGER :: Jnorth
!      INTEGER :: Jsouth

      INTEGER :: ibnd
      REAL(r8) :: taunudg
      REAL(r8) :: WRK
      REAL(r8) :: pitmp

!      Iwest = INT(width_tsnudge)
!      Ieast = IMT+1-INT(width_tsnudge)
!
!      Jnorth = INT(width_tsnudge)
!      Jsouth = JMT+1-INT(width_tsnudge)
!
!      taunudg=TAUOUT_TS
      taunudg=TAUIN_TS

      pitmp = 3.1415926

      WRK = 0.0D0

      ibnd = INT(width_tsnudge)

      IF ( ix.eq.0 ) THEN

!$OMP PARALLEL DO PRIVATE (I,J)
           DO I = 1, IMT
              ibnd = min(I,ibnd)
              DO J = 1, JMT
                 WRK = 0.5D0*(cos(pitmp*float(ibnd-1)/(width_tsnudge-1.0))+1.0)
                 TSnudgcof(I,J) = 1.0D0/taunudg*WRK
              ENDDO
           ENDDO

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (I,J)
              DO I = 1, IMT
                 DO J = 1, JMT
                    AT(I,J,K,1) = AT(I,J,K,1)-DTS*TSnudgcof(I,J)*(ATB(I,J,K,1)-TT_MON_LOC(I,J,K)*VIT(I,J,K)) 
                    AT(I,J,K,2) = AT(I,J,K,2)-DTS*TSnudgcof(I,J)*(ATB(I,J,K,2)-SS_MON_LOC(I,J,K)*VIT(I,J,K))                
                 ENDDO
              ENDDO
           ENDDO
                    
      ENDIF

      IF ( ix.eq.(nx_proc-1) ) THEN

!$OMP PARALLEL DO PRIVATE (I,J)
           DO I = 1, IMT
              ibnd = min(IMT+1-I,ibnd)
              DO J = 1, JMT
                 WRK = 0.5D0*(cos(pitmp*float(ibnd-1)/(width_tsnudge-1.0))+1.0)
                 TSnudgcof(I,J) = 1.0D0/taunudg*WRK
              ENDDO
           ENDDO

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (I,J)
              DO I = 1, IMT
                 DO J = 1, JMT
                    AT(I,J,K,1) = AT(I,J,K,1)-DTS*TSnudgcof(I,J)*(ATB(I,J,K,1)-TT_MON_LOC(I,J,K)*VIT(I,J,K)) 
                    AT(I,J,K,2) = AT(I,J,K,2)-DTS*TSnudgcof(I,J)*(ATB(I,J,K,2)-SS_MON_LOC(I,J,K)*VIT(I,J,K))               
                 ENDDO
              ENDDO
           ENDDO
                    
      ENDIF


      IF ( iy.eq.0 ) THEN

!$OMP PARALLEL DO PRIVATE (I,J)
           DO I = 1, IMT
              DO J = 1, JMT
                 ibnd = min(J,ibnd)
                 WRK = 0.5D0*(cos(pitmp*float(ibnd-1)/(width_tsnudge-1.0))+1.0)
                 TSnudgcof(I,J) = 1.0D0/taunudg*WRK
              ENDDO
           ENDDO

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (J,I)
              DO J = 1, JMT
                 DO I = 1, IMT
                    AT(I,J,K,1) = AT(I,J,K,1)-DTS*TSnudgcof(I,J)*(ATB(I,J,K,1)-TT_MON_LOC(I,J,K)*VIT(I,J,K)) 
                    AT(I,J,K,2) = AT(I,J,K,2)-DTS*TSnudgcof(I,J)*(ATB(I,J,K,2)-SS_MON_LOC(I,J,K)*VIT(I,J,K))                
                 ENDDO
              ENDDO
           ENDDO

      ENDIF

      IF ( iy.eq.(ny_proc-1) ) THEN

!$OMP PARALLEL DO PRIVATE (I,J)
           DO I = 1, IMT
              DO J = 1, JMT
                 ibnd = min(JMT+1-J,ibnd)
                 WRK = 0.5D0*(cos(pitmp*float(ibnd-1)/(width_tsnudge-1.0))+1.0)
                 TSnudgcof(I,J) = 1.0D0/taunudg*WRK
              ENDDO
           ENDDO

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (J,I)
              DO J = 1, JMT
                 DO I = 1, IMT
                    AT(I,J,K,1) = AT(I,J,K,1)-DTS*TSnudgcof(I,J)*(ATB(I,J,K,1)-TT_MON_LOC(I,J,K)*VIT(I,J,K)) 
                    AT(I,J,K,2) = AT(I,J,K,2)-DTS*TSnudgcof(I,J)*(ATB(I,J,K,2)-SS_MON_LOC(I,J,K)*VIT(I,J,K))                
                 ENDDO
              ENDDO
           ENDDO

      ENDIF

#endif

      RETURN

      END SUBROUTINE openbc_ts_nudge
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
