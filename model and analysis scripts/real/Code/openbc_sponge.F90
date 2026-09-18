!
!

      SUBROUTINE openbc_sponge

#include <def-undef.h>
use   precision_mod

use   param_mod,  only: IMT, JMT, KM, I, J, K, nx_proc, ny_proc
use   pconst_mod, only: VIT, VIV, ix, iy, AM3, AH3, AM !,am3, ah3, AM, AH   

#if ( defined SPMD )
use msg_mod
#endif

#ifdef OBCDT
use   openbc_mod
#endif

#ifdef OBCSPONGE

      IMPLICIT NONE

      REAL(r8) :: outerF_UV
      REAL(r8) :: outerF_TS
      REAL(r8) :: innerF
      REAL(r8) :: val

      REAL(r8) :: fac_uv
      REAL(r8) :: fac_ts

      INTEGER :: Iwest
      INTEGER :: Ieast
      INTEGER :: Jnorth
      INTEGER :: Jsouth

      fac_uv = 20.0D+0
      fac_ts = 20.0D+0

      Iwest = INT(width_sponge)
      Ieast = IMT+1-INT(width_sponge)

      Jnorth = INT(width_sponge)
      Jsouth = JMT+1-INT(width_sponge)

      IF ( ix.eq.0 ) THEN

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (I,J)
              DO I = 1, Iwest
                 DO J = 1, JMT
                    innerF = AM(J)
                    outerF_UV = innerF*fac_uv
                    val = innerF+(outerF_UV-innerF)*(width_sponge-1.0D0*I)/width_sponge
#ifdef BIHAR
                    AM(J) = MIN(MAX(val,outerF_UV),innerF)

#else
                    AM(J) = MAX(MIN(val,outerF_UV),innerF)
#endif

                    innerF = AM3(I,J,K)*VIV(I,J,K)
                    outerF_UV = innerF*fac_uv
                    val = innerF+(outerF_UV-innerF)*(width_sponge-1.0D0*I)/width_sponge
#ifdef BIHAR
                    AM3(I,J,K) = MIN(MAX(val,outerF_UV),innerF)

#else
                    AM3(I,J,K) = MAX(MIN(val,outerF_UV),innerF)
#endif

                    innerF = AH3(I,J,K)*VIT(I,J,K)
                    outerF_TS = innerF*fac_ts
                    val = innerF+(outerF_TS-innerF)*(width_sponge-1.0D0*I)/width_sponge
#ifdef BIHAR

                    AH3(I,J,K) = MIN(MAX(val,outerF_TS),innerF)
#else
                    AH3(I,J,K) = MAX(MIN(val,outerF_TS),innerF)
#endif
                 ENDDO
              ENDDO
           ENDDO
                    
      ENDIF

      IF ( ix.eq.(nx_proc-1) ) THEN

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (I,J)
              DO I = Ieast, IMT
                 DO J = 1, JMT
                    innerF = AM(J)
                    outerF_UV = innerF*fac_uv
                    val = outerF_UV+(innerF-outerF_UV)*(1.0D0*(IMT+1-I))/width_sponge
#ifdef BIHAR
                    AM(J) = MIN(MAX(val,outerF_UV),innerF)
#else
                    AM(J) = MAX(MIN(val,outerF_UV),innerF)
#endif

                    innerF = AM3(I,J,K)*VIV(I,J,K)
                    outerF_UV = innerF*fac_uv
                    val = outerF_UV+(innerF-outerF_UV)*(1.0D0*(IMT+1-I))/width_sponge
#ifdef BIHAR
                    AM3(I,J,K) = MIN(MAX(val,outerF_UV),innerF)
#else
                    AM3(I,J,K) = MAX(MIN(val,outerF_UV),innerF)
#endif

                    innerF = AH3(I,J,K)*VIT(I,J,K)
                    outerF_TS = innerF*fac_ts
                    val = outerF_TS+(innerF-outerF_TS)*(1.0D0*(IMT+1-I))/width_sponge
#ifdef BIHAR
                    AH3(I,J,K) = MIN(MAX(val,outerF_TS),innerF)
#else
                    AH3(I,J,K) = MAX(MIN(val,outerF_TS),innerF)
#endif
                 ENDDO
              ENDDO
           ENDDO
                    
      ENDIF


      IF ( iy.eq.0 ) THEN

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (J,I)
              DO J = 1, Jnorth
                 DO I = 1, IMT
                    innerF = AM(J)
                    outerF_UV = innerF*fac_uv
                    val = innerF+(outerF_UV-innerF)*(width_sponge-1.0D0*J)/width_sponge
#ifdef BIHAR
                    AM(J) = MIN(MAX(val,outerF_UV),innerF)
#else
                    AM(J) = MAX(MIN(val,outerF_UV),innerF)
#endif

                    innerF = AM3(I,J,K)*VIV(I,J,K)
                    outerF_UV = innerF*fac_uv
                    val = innerF+(outerF_UV-innerF)*(width_sponge-1.0D0*J)/width_sponge
#ifdef BIHAR
                    AM(J) = MIN(MAX(val,outerF_UV),innerF)
                    AM3(I,J,K) = MIN(MAX(val,outerF_UV),innerF)

#else
                    AM(J) = MAX(MIN(val,outerF_UV),innerF)
                    AM3(I,J,K) = MAX(MIN(val,outerF_UV),innerF)
#endif

                    innerF = AH3(I,J,K)*VIT(I,J,K)
                    outerF_TS = innerF*fac_ts
                    val = innerF+(outerF_TS-innerF)*(width_sponge-1.0D0*J)/width_sponge
#ifdef BIHAR
                    AH3(I,J,K) = MIN(MAX(val,outerF_TS),innerF)
#else
                    AH3(I,J,K) = MAX(MIN(val,outerF_TS),innerF)
#endif
                 ENDDO
              ENDDO
           ENDDO

      ENDIF

      IF ( iy.eq.(ny_proc-1) ) THEN

           DO K = 1, KM
!$OMP PARALLEL DO PRIVATE (J,I)
              DO J = Jsouth, JMT
                 DO I = 1, IMT
                    innerF = AM(J)
                    outerF_UV = innerF*fac_uv
                    val = outerF_UV+(innerF-outerF_UV)*(1.0D0*(JMT+1-J))/width_sponge
#ifdef BIHAR
                    AM(J) = MIN(MAX(val,outerF_UV),innerF)
#else
                    AM(J) = MAX(MIN(val,outerF_UV),innerF)
#endif

                    innerF = AM3(I,J,K)*VIV(I,J,K)
                    outerF_UV = innerF*fac_uv
                    val = outerF_UV+(innerF-outerF_UV)*(1.0D0*(JMT+1-J))/width_sponge
#ifdef BIHAR
                    AM(J) = MIN(MAX(val,outerF_UV),innerF)
                    AM3(I,J,K) = MIN(MAX(val,outerF_UV),innerF)
#else
                    AM(J) = MAX(MIN(val,outerF_UV),innerF)
                    AM3(I,J,K) = MAX(MIN(val,outerF_UV),innerF)
#endif

                    innerF = AH3(I,J,K)*VIT(I,J,K)
                    outerF_TS = innerF*fac_ts
                    val = outerF_TS+(innerF-outerF_TS)*(1.0D0*(JMT+1-J))/width_sponge
#ifdef BIHAR
                    AH3(I,J,K) = MIN(MAX(val,outerF_TS),innerF)
#else
                    AH3(I,J,K) = MAX(MIN(val,outerF_TS),innerF)
#endif
                 ENDDO
              ENDDO
           ENDDO

      ENDIF

#endif

      RETURN

      END SUBROUTINE openbc_sponge
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
