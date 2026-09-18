!
!

      SUBROUTINE openbc_uv_nudge

#include <def-undef.h>
use   precision_mod

use   param_mod,  only: IMT, JMT, KM, I, J, K, nx_proc, ny_proc
use   pconst_mod, only: VIT, VIV, ix, iy, DTC   
use   dyn_mod,    only: U, UP, V, VP

#if ( defined SPMD )
use msg_mod
#endif

#ifdef OBCDT
!use   openbc_mod, only: UU_MON_LOC, VV_MON_LOC, TAUOUT, width_uvnudge, UVnudgcof
use   openbc_mod, only: UU_MON_LOC, VV_MON_LOC, TAUIN, width_uvnudge, UVnudgcof
#endif

#ifdef UV_NUDGE

      IMPLICIT NONE

!      INTEGER :: Iwest
!      INTEGER :: Ieast
!      INTEGER :: Jnorth
!      INTEGER :: Jsouth
!
      INTEGER :: ibnd

      REAL(r8) :: taunudg
      REAL(r8) :: pitmp
      REAL(r8) :: WRK

!      Iwest = INT(width_uvnudge)
!      Ieast = IMT+1-INT(width_uvnudge)
!
!      Jnorth = INT(width_uvnudge)
!      Jsouth = JMT+1-INT(width_uvnudge)
!
!      taunudg=TAUOUT
      taunudg=TAUIN

      pitmp = 3.1415926

      WRK = 0.0D0

      ibnd = INT(width_uvnudge)

      IF ( ix.eq.0 ) THEN

!$OMP PARALLEL DO PRIVATE (I,J)
           DO I = 1, IMT
              ibnd = min(I,ibnd)
              DO J = 1, JMT
                 WRK = 0.5D0*(cos(pitmp*float(ibnd-1)/(width_uvnudge-1.0))+1.0)
                 UVnudgcof(I,J) = 1.0D0/taunudg*WRK
              ENDDO
           ENDDO

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (I,J)
              DO I = 1, IMT
                 DO J = 1, JMT
                    U(I,J,K) = U(I,J,K)-DTC*UVnudgcof(I,J)*(UP(I,J,K)-UU_MON_LOC(I,J,K)*VIV(I,J,K)) 
                    V(I,J,K) = V(I,J,K)-DTC*UVnudgcof(I,J)*(VP(I,J,K)-VV_MON_LOC(I,J,K)*VIV(I,J,K))                      
                 ENDDO
              ENDDO
           ENDDO
                    
      ENDIF

      IF ( ix.eq.(nx_proc-1) ) THEN

!$OMP PARALLEL DO PRIVATE (I,J)
           DO I = 1, IMT
              ibnd = min(IMT+1-I,ibnd)
              DO J = 1, JMT
                 WRK = 0.5D0*(cos(pitmp*float(ibnd-1)/(width_uvnudge-1.0))+1.0)
                 UVnudgcof(I,J) = 1.0D0/taunudg*WRK
              ENDDO
           ENDDO

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (I,J)
              DO I = 1, IMT
                 DO J = 1, JMT
                    U(I,J,K) = U(I,J,K)-DTC*UVnudgcof(I,J)*(UP(I,J,K)-UU_MON_LOC(I,J,K)*VIV(I,J,K)) 
                    V(I,J,K) = V(I,J,K)-DTC*UVnudgcof(I,J)*(VP(I,J,K)-VV_MON_LOC(I,J,K)*VIV(I,J,K))                        
                 ENDDO
              ENDDO
           ENDDO
                    
      ENDIF


      IF ( iy.eq.0 ) THEN

!$OMP PARALLEL DO PRIVATE (I,J)
           DO I = 1, IMT
              DO J = 1, JMT
                 ibnd = min(J,ibnd)
                 WRK = 0.5D0*(cos(pitmp*float(ibnd-1)/(width_uvnudge-1.0))+1.0)
                 UVnudgcof(I,J) = 1.0D0/taunudg*WRK
              ENDDO
           ENDDO

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (J,I)
              DO J = 1, JMT
                 DO I = 1, IMT
                    U(I,J,K) = U(I,J,K)-DTC*UVnudgcof(I,J)*(UP(I,J,K)-UU_MON_LOC(I,J,K)*VIV(I,J,K)) 
                    V(I,J,K) = V(I,J,K)-DTC*UVnudgcof(I,J)*(VP(I,J,K)-VV_MON_LOC(I,J,K)*VIV(I,J,K))                        
                 ENDDO
              ENDDO
           ENDDO

      ENDIF

      IF ( iy.eq.(ny_proc-1) ) THEN

!$OMP PARALLEL DO PRIVATE (I,J)
           DO I = 1, IMT
              DO J = 1, JMT
                 ibnd = min(JMT+1-J,ibnd)
                 WRK = 0.5D0*(cos(pitmp*float(ibnd-1)/(width_uvnudge-1.0))+1.0)
                 UVnudgcof(I,J) = 1.0D0/taunudg*WRK
              ENDDO
           ENDDO

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (J,I)
              DO J = 1, JMT
                 DO I = 1, IMT
                    U(I,J,K) = U(I,J,K)-DTC*UVnudgcof(I,J)*(UP(I,J,K)-UU_MON_LOC(I,J,K)*VIV(I,J,K)) 
                    V(I,J,K) = V(I,J,K)-DTC*UVnudgcof(I,J)*(VP(I,J,K)-VV_MON_LOC(I,J,K)*VIV(I,J,K))                        
                 ENDDO
              ENDDO
           ENDDO

      ENDIF

#endif

      RETURN

      END SUBROUTINE openbc_uv_nudge
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
