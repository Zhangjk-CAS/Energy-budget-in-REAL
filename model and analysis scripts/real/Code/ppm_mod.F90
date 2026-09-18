!
!      Open Boundary Condition Module
!
module ppm_mod

      use precision_mod
      use param_mod
      use pconst_mod


      contains

      SUBROUTINE ADVEC_PPM_UV( UPPM, VPPM, WPPM, VARTMP, DLTMP )

#include <def-undef.h>
      use   precision_mod
      use   param_mod,  only: IMT, JMT, KM, I, J, K
      use   pconst_mod, only: OUX, OUY, ZKT, ZKP, DTC

      IMPLICIT NONE

      REAL(r8) :: UPPM(IMT,JMT,KM)
      REAL(r8) :: VPPM(IMT,JMT,KM)
      REAL(r8) :: WPPM(IMT,JMT,KM)

      REAL(r8) :: VARTMP(IMT,JMT,KM)
      REAL(r8) :: DLTMP(IMT,JMT,KM)

      integer :: NST
      REAL(r8), ALLOCATABLE :: VQC(:)
      REAL(r8), ALLOCATABLE :: DXYZ(:)
      INTEGER :: INDXYZ

      REAL(r8), ALLOCATABLE :: UVWhaf(:)
      REAL(r8), ALLOCATABLE :: flxqR(:), flxqL(:)

      REAL(r8) :: CMIU

      REAL(r8) :: uhafR, uhafL
!
      REAL(r8), ALLOCATABLE :: FPAVE(:,:,:), FCAVE(:,:,:)
      REAL(r8), ALLOCATABLE :: FDAVE(:,:,:), FEAVE(:,:,:), FNP1AVE(:,:,:)
!
!
      NST = 0

      DO K = 1, KM
         DO J = 1, JMT
            DO I = 1, IMT
               DLTMP(I,J,K) = 0.0D0
            ENDDO
        ENDDO
      ENDDO

      IF ( ALLOCATED( FPAVE ) ) DEALLOCATE( FPAVE )
      ALLOCATE( FPAVE(IMT,JMT,KM) )

      IF ( ALLOCATED( FCAVE ) ) DEALLOCATE( FCAVE )
      ALLOCATE( FCAVE(IMT,JMT,KM) )

      IF ( ALLOCATED( FDAVE ) ) DEALLOCATE( FDAVE )
      ALLOCATE( FDAVE(IMT,JMT,KM) )

      IF ( ALLOCATED( FEAVE ) ) DEALLOCATE( FEAVE )
      ALLOCATE( FEAVE(IMT,JMT,KM) )

      IF ( ALLOCATED( FNP1AVE ) ) DEALLOCATE( FNP1AVE )
      ALLOCATE( FNP1AVE(IMT,JMT,KM) )

      DO J = 1, JMT

         DO K = 1, KM

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(IMT) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(IMT) )

            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(IMT) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(IMT) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(IMT) )

            DO I = 1, IMT
               VQC(I) = VARTMP(I,J,K)
            ENDDO

            DO I = 1, IMT
               DXYZ(I) = 1.0D0/OUX(J)
            ENDDO
            DO I = 2, IMT
               UVWhaf(I) = 0.50D0*(UPPM(I-1,J,K)+UPPM(I,J,K))
            ENDDO
            UVWhaf(1) = 0.0D0

            INDXYZ = 1
!
!            To Compute the Flux Through the Cell

            CALL QRLFLX_PPM_NQGRID( INDXYZ, IMT, VQC, DXYZ, UVWhaf, DTC, NST, flxqR, flxqL )
!
            DO I = 1, IMT

               IF ( I.EQ.1 ) THEN

                  uhafR = UVWhaf(I+1)
                  uhafL = 0.0D0

                  CMIU = DTC*OUX(J)

                  FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                  FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUX(J)

               ELSE IF ( I.GE.2 .AND. I.LE.(IMT-1) ) THEN

                  uhafR = UVWhaf(I+1)
                  uhafL = UVWhaf(I)

                  CMIU = DTC*OUX(J)
                  FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                  FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUX(J)

               ELSE IF ( I.EQ.IMT ) THEN

                  uhafR = 0.0D0
                  uhafL = UVWhaf(I)

                  CMIU = DTC*OUX(J)
                  FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                  FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUX(J)

               ENDIF

            ENDDO   ! end I

         ENDDO   ! end K

         DO I = 1, IMT

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(KM) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(KM) )

            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(KM) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(KM) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(KM) )

            DO K = 1, KM
               VQC(K) = FCAVE(I,J,K)
            ENDDO

            DO K = 1, KM
               DXYZ(K) = ZKP(K)-ZKP(K+1)
            ENDDO

            DO K = 1, KM-1
               UVWhaf(K) = -WPPM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0

            INDXYZ = 3
!
            CALL QRLFLX_PPM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTC, NST, flxqR, flxqL )
!
!            To Compute the Flux Through the Cell
!

            DO K = 1, KM

               CMIU = DTC/DXYZ(K)
               FDAVE(I,J,K) = FCAVE(I,J,K) + CMIU*(flxqL(K)-flxqR(K))

            ENDDO ! end K

         ENDDO  ! end I

      ENDDO  ! end J
      DO K = 1, KM
      DO J = 1, JMT
      DO I = 1, IMT
         FCAVE(I,J,K) = FDAVE(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      DO I = 1, IMT

         DO K = 1, KM

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(JMT) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(JMT) )
!
            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(JMT) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(JMT) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(JMT) )
            DO J = 1, JMT
               VQC(J) = FCAVE(I,J,K)
            ENDDO

            DO J = 1, JMT
               DXYZ(J) = 1.0D0/OUY(J)
            ENDDO

            DO J = 1, JMT-1
               UVWhaf(J) = 0.50D0*(VPPM(I,J,K)+VPPM(I,J+1,K))
            ENDDO
            UVWhaf(JMT) = 0.0D0

            INDXYZ = 2
!

            CALL QRLFLX_PPM_NQGRID( INDXYZ, JMT, VQC, DXYZ, UVWhaf, DTC, NST, flxqR, flxqL )
!
!            To Compute the Flux Through the Cell
            DO J = 1, JMT

               IF ( J.EQ.1 ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = 0.0D0

                    CMIU = DTC*OUY(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))
                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUY(J)

               ELSE IF ( J.GE.2 .AND. J.LE.(JMT-1) ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = UVWhaf(J-1)

                    CMIU = DTC*OUY(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))

                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUY(J)
               ELSE IF ( J.EQ.JMT ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = UVWhaf(J-1)

                    CMIU = DTC*OUY(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))

                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUY(J)

               ENDIF

            ENDDO  ! end J

         ENDDO  ! end K
         DO J = 1, JMT

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(KM) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(KM) )
!
            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(KM) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(KM) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(KM) )

            DO K = 1, KM
               VQC(K) = FEAVE(I,J,K)
            ENDDO
!
            DO K = 1, KM
               DXYZ(K) = ZKP(K)-ZKP(K+1)
            ENDDO

            DO K = 1, KM-1
               UVWhaf(K) = -WPPM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0
            INDXYZ = 3
!
            CALL QRLFLX_PPM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTC, NST, flxqR, flxqL )
!
!            To Compute the Flux Through the Cell
!
            DO K = 1, KM

               CMIU = DTC/DXYZ(K)
               FNP1AVE(I,J,K) = FDAVE(I,J,K) + CMIU*(flxqL(K)-flxqR(K))

            ENDDO ! end K

         ENDDO  ! end j

      ENDDO  ! end I
      DO K = 1, KM
         DO J = 1, JMT
            DO I = 1, IMT

               DLTMP(I,J,K) = ( FNP1AVE(I,J,K)-VARTMP(I,J,K) )/DTC

            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE( VQC )
      DEALLOCATE( DXYZ )

      DEALLOCATE( UVWhaf )
      DEALLOCATE( flxqL )
      DEALLOCATE( flxqR )

      DEALLOCATE( FPAVE )
      DEALLOCATE( FCAVE )
      DEALLOCATE( FDAVE )
      DEALLOCATE( FEAVE )
      DEALLOCATE( FNP1AVE )

      RETURN

      END SUBROUTINE ADVEC_PPM_UV

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE ADVEC_PPM_TS( UPPM, VPPM, WPPM, VARTMP, NST, DLTMP, AXTMP, AYTMP, AZTMP )

#include <def-undef.h>
use   precision_mod
use   param_mod,  only: IMT, JMT, KM, KMP1, I, J, K
use   pconst_mod, only: OTX, DYT, ZKT, ZKP, DTS

      IMPLICIT NONE

      REAL(r8) :: UPPM(IMT,JMT,KM)
      REAL(r8) :: VPPM(IMT,JMT,KM)
      REAL(r8) :: WPPM(IMT,JMT,KMP1)
      REAL(r8) :: VARTMP(IMT,JMT,KM)
      REAL(r8) :: DLTMP(IMT,JMT,KM)

      integer, intent(in) :: NST

      REAL(r8) :: AXTMP(IMT,JMT,KM)
      REAL(r8) :: AYTMP(IMT,JMT,KM)
      REAL(r8) :: AZTMP(IMT,JMT,KM)

      REAL(r8), ALLOCATABLE :: VQC(:)
      REAL(r8), ALLOCATABLE :: DXYZ(:)
      INTEGER :: INDXYZ
!
      REAL(r8), ALLOCATABLE :: flxqR(:), flxqL(:)

      REAL(r8), ALLOCATABLE :: UVWhaf(:)

      REAL(r8) :: CMIU

      REAL(r8) :: uhafR, uhafL
!
      REAL(r8), ALLOCATABLE :: FPAVE(:,:,:), FCAVE(:,:,:)
      REAL(r8), ALLOCATABLE :: FDAVE(:,:,:), FEAVE(:,:,:), FNP1AVE(:,:,:)
!
      DO K = 1, KM
         DO J = 1, JMT
            DO I = 1, IMT
               DLTMP(I,J,K) = 0.0D0
               AXTMP(I,J,K) = 0.0D0
               AYTMP(I,J,K) = 0.0D0
               AZTMP(I,J,K) = 0.0D0
            ENDDO
        ENDDO
      ENDDO


      IF ( ALLOCATED( FPAVE ) ) DEALLOCATE( FPAVE )
      ALLOCATE( FPAVE(IMT,JMT,KM) )

      IF ( ALLOCATED( FCAVE ) ) DEALLOCATE( FCAVE )
      ALLOCATE( FCAVE(IMT,JMT,KM) )

      IF ( ALLOCATED( FDAVE ) ) DEALLOCATE( FDAVE )
      ALLOCATE( FDAVE(IMT,JMT,KM) )

      IF ( ALLOCATED( FEAVE ) ) DEALLOCATE( FEAVE )
      ALLOCATE( FEAVE(IMT,JMT,KM) )

      IF ( ALLOCATED( FNP1AVE ) ) DEALLOCATE( FNP1AVE )
      ALLOCATE( FNP1AVE(IMT,JMT,KM) )

      DO J = 2, JMT

         DO K = 1, KM

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(IMT) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(IMT) )
!
            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(IMT) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(IMT) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(IMT) )

            DO I = 1, IMT
               VQC(I) = VARTMP(I,J,K)
            ENDDO

            DO I = 1, IMT
               DXYZ(I) = 1.0D0/OTX(J)
            ENDDO
            DO I = 1, IMT
               UVWhaf(I) = 0.50D0*(UPPM(I,J-1,K)+UPPM(I,J,K))
            ENDDO

            INDXYZ = 1

            CALL QRLFLX_PPM_NQGRID( INDXYZ, IMT, VQC, DXYZ, UVWhaf, DTS, NST, flxqR, flxqL )
! 
!
!            To Compute the Flux Through the Cell
!
            DO I = 1, IMT

              IF ( I.EQ.1 ) THEN

                   uhafR = UVWhaf(I+1)
                   uhafL = UVWhaf(I)

                   CMIU = DTS*OTX(J)
                   FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                   FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)*OTX(J)

!                  AXTMP(I,J,K) = ( FCAVE(I,J,K)-VARTMP(I,J,K) )/DTS
                   AXTMP(I,J,K) = ( FPAVE(I,J,K)-VARTMP(I,J,K) )/DTS

              ELSE IF ( I.GE.2 .AND. I.LE.(IMT-1) ) THEN

                   uhafR = UVWhaf(I+1)
                   uhafL = UVWhaf(I)

                   CMIU = DTS*OTX(J)
                   FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                   FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)*OTX(J)
!
                   AXTMP(I,J,K) = ( FPAVE(I,J,K)-VARTMP(I,J,K) )/DTS

               ELSE IF ( I.EQ.IMT ) THEN

                   uhafR = 0.0D0
                   uhafL = UVWhaf(I)

                   CMIU = DTS*OTX(J)
                   FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                   FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)*OTX(J)
!
!                   AXTMP(I,J,K) = ( FCAVE(I,J,K)-VARTMP(I,J,K) )/DTS
                   AXTMP(I,J,K) = ( FPAVE(I,J,K)-VARTMP(I,J,K) )/DTS

               ENDIF

            ENDDO    ! end I

         ENDDO    ! end K

         DO I = 1, IMT

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(KM) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(KM) )
!
            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(KM) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(KM) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(KM) )

            DO K = 1, KM
               VQC(K) = FCAVE(I,J,K)
            ENDDO
!
            DO K = 1, KM
               DXYZ(K) = ZKP(K)-ZKP(K+1)
            ENDDO

            DO K = 1, KM-1
               UVWhaf(K) = -WPPM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0

            INDXYZ = 3

            CALL QRLFLX_PPM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTS, NST, flxqR, flxqL )
! 
!
!            To Compute the Flux Through the Cell
!
            DO K = 1, KM

                  CMIU = DTS/DXYZ(K)
                  FDAVE(I,J,K) = FCAVE(I,J,K) + CMIU*(flxqL(K)-flxqR(K))
                  AZTMP(I,J,K) = ( FDAVE(I,J,K)-FCAVE(I,J,K) )/DTS

            ENDDO ! end K

         ENDDO  ! end I

      ENDDO  !end J

      DO K = 1, KM
      DO I = 1, IMT
         FPAVE(I,1,K) = FPAVE(I,2,K)
         FCAVE(I,1,K) = FCAVE(I,2,K)
         FDAVE(I,1,K) = FDAVE(I,2,K)
         AXTMP(I,1,K) = AXTMP(I,2,K)
         AZTMP(I,1,K) = AZTMP(I,2,K)
      ENDDO
      ENDDO

      DO K = 1, KM
      DO J = 1, JMT
      DO I = 1, IMT
         FCAVE(I,J,K) = FDAVE(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      DO I = 1, IMT-1

         DO K = 1, KM

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(JMT) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(JMT) )
!
            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(JMT) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(JMT) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(JMT) )

            DO J = 1, JMT
               VQC(J) = FCAVE(I,J,K)
            ENDDO

            DO J = 1, JMT
               DXYZ(J) = DYT(J)
            ENDDO

            DO J = 1, JMT
               UVWhaf(J) = 0.50D0*(VPPM(I,J,K)+VPPM(I+1,J,K))
            ENDDO

            INDXYZ = 2

            CALL QRLFLX_PPM_NQGRID( INDXYZ, JMT, VQC, DXYZ, UVWhaf, DTS, NST, flxqR, flxqL )
!
!
!            To Compute the Flux Through the Cell
!
            DO J = 1, JMT

               IF ( J.EQ.1 ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = 0.0D0

                    CMIU = DTS/DYT(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))

                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)/DYT(J)
!
!                   AYTMP(I,J,K) = ( FEAVE(I,J,K) - FPAVE(I,J,K) )/DTS
                    AYTMP(I,J,K) = ( FDAVE(I,J,K) - FPAVE(I,J,K) )/DTS

               ELSE IF ( J.GE.2 .AND. J.LE.(JMT-1) ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = UVWhaf(J-1)
!
                    CMIU = DTS/DYT(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))

                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)/DYT(J)
!
!                   AYTMP(I,J,K) = ( FEAVE(I,J,K) - FPAVE(I,J,K) )/DTS
                    AYTMP(I,J,K) = ( FDAVE(I,J,K) - FPAVE(I,J,K) )/DTS

               ELSE IF ( J.EQ.JMT ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = UVWhaf(J-1)
!
                    CMIU = DTS/DYT(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))

                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)/DYT(J)
!
!                   AYTMP(I,J,K) = ( FEAVE(I,J,K) - FPAVE(I,J,K) )/DTS
                    AYTMP(I,J,K) = ( FDAVE(I,J,K) - FPAVE(I,J,K) )/DTS

               ENDIF

            ENDDO  ! end J

         ENDDO  ! end K
         DO J = 1, JMT

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(KM) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(KM) )
!
            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(KM) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(KM) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(KM) )

            DO K = 1, KM
               VQC(K) = FEAVE(I,J,K)
            ENDDO
!
            DO K = 1, KM
               DXYZ(K) = ZKP(K)-ZKP(K+1)
            ENDDO
            DO K = 1, KM-1
              UVWhaf(K) = -WPPM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0

            INDXYZ = 3
!
            CALL QRLFLX_PPM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTS, NST, flxqR, flxqL )
! 
!
!            To Compute the Flux Through the Cell
!
            DO K = 1, KM

               CMIU = DTS/DXYZ(K)
               FNP1AVE(I,J,K) = FDAVE(I,J,K) + CMIU*(flxqL(K)-flxqR(K))
               AZTMP(I,J,K) = ( FNP1AVE(I,J,K)-FDAVE(I,J,K) )/DTS

            ENDDO ! end K

         ENDDO  ! end j

      ENDDO  ! end I

      DO K = 1, KM
         DO J = 1, JMT
            FDAVE(IMT,J,K) = FDAVE(IMT-1,J,K)
            FEAVE(IMT,J,K) = FEAVE(IMT-1,J,K)
            FNP1AVE(IMT,J,K) = FNP1AVE(IMT-1,J,K)
            AYTMP(IMT,J,K) = AYTMP(IMT-1,J,K)
            AZTMP(IMT,J,K) = AZTMP(IMT-1,J,K)
         ENDDO
      ENDDO
!

      DO K = 1, KM
         DO J = 1, JMT
            DO I = 1, IMT

               DLTMP(I,J,K) = ( FNP1AVE(I,J,K)-VARTMP(I,J,K) )/DTS

            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE( VQC )
      DEALLOCATE( DXYZ )
      DEALLOCATE( UVWhaf )
      DEALLOCATE( flxqL )
      DEALLOCATE( flxqR )

      DEALLOCATE( FPAVE )
      DEALLOCATE( FCAVE )
      DEALLOCATE( FDAVE )
      DEALLOCATE( FEAVE )
      DEALLOCATE( FNP1AVE )

      RETURN

      END SUBROUTINE ADVEC_PPM_TS
!
!-----------------------------------------------------------------------------------------------------------
!
      SUBROUTINE QRLFLX_PPM_NQGRID( KXYZ, IXYZ, VQ, DltXYZ, Suvw, Dtsc, IndST, flxR, flxL)

       use   precision_mod

       IMPLICIT NONE

       INTEGER :: KXYZ
       INTEGER :: IXYZ
       REAL(r8) :: VQ(IXYZ), DltXYZ(IXYZ)
       REAL(r8) :: Suvw(IXYZ)
       REAL(r8) :: flxR(IXYZ), flxL(IXYZ)

       INTEGER ::  IndST
       REAL(r8) :: Dtsc

       REAL(r8), ALLOCATABLE :: dltQ(:)
       REAL(r8), ALLOCATABLE :: qR(:), qL(:)
       REAL(r8), ALLOCATABLE :: alphaR(:), alphaL(:)


       REAL(r8) :: COF

       REAL(r8) :: deltX1, deltX2, deltX3, deltX4
       REAL(r8) :: deltX5, deltX6, deltX7
       REAL(r8) :: dltc, dltm, dltp
       REAL(r8) :: dltca, dltma, dltpa

       REAL(r8) :: VQmin, VQmax

       INTEGER :: N

       REAL(r8) :: tmpa, tmpb, tmpc

       REAL(r8) :: VQstr, VQlac

       REAL(r8) :: uhafR, uhafL
       REAL(r8) :: FMIU
!
       INTEGER :: Istpn

       REAL(r8) :: VQminR, VQmaxR
       REAL(r8) :: VQminL, VQmaxL
       REAL(r8) :: Rjm1, Rjp1
       REAL(r8) :: VQvlin, VQvlac
       REAL(r8) :: VQvmin, VQvmax
       REAL(r8) :: VQvhfjL, VQvhfjR
       REAL(r8) :: pjp2, pjmp2, pjpp2
       REAL(r8) :: aj, bj, abjmax
       REAL(r8) :: EPSM
       REAL(r8) :: ck
       REAL(r8) :: cp2
       REAL(r8) :: cstpn
       IF ( ALLOCATED( dltQ ) ) DEALLOCATE( dltQ )
       ALLOCATE( dltQ(IXYZ) )

       IF ( ALLOCATED( qR ) ) DEALLOCATE( qR )
       ALLOCATE( qR(IXYZ) )

       IF ( ALLOCATED( qL )) DEALLOCATE( qL )
       ALLOCATE( qL(IXYZ) )

       IF ( ALLOCATED( alphaR ) ) DEALLOCATE( alphaR )
       ALLOCATE( alphaR(IXYZ) )

       IF ( ALLOCATED( alphaL ) ) DEALLOCATE( alphaL )
       ALLOCATE( alphaL(IXYZ) )


       COF = 2.0D0

       IF ( IndST.EQ.2 ) THEN
            EPSM = 1.0D-12
       ELSE
            EPSM = 1.0D-6
       ENDIF

       ck = 2.0D0

       cp2 = 0.9D0
!       cp2 = 1.0D0
       Istpn = 1

       cstpn = 3.0D0
!
!            Apply the conventional van Leer limiter
!
       DO N = 2, IXYZ-1

          deltX1 = DltXYZ(N)/( DltXYZ(N-1)+DltXYZ(N)+DltXYZ(N+1) )
          deltX2 = ( 2.0D0*DltXYZ(N-1)+DltXYZ(N) )/( DltXYZ(N+1)+DltXYZ(N) )
          deltX3 = ( DltXYZ(N)+2.0D0*DltXYZ(N+1) )/( DltXYZ(N-1)+DltXYZ(N) )
!
!          dltc = deltX1*(deltX2*(VQ(N+1)-VQ(N))+deltX3*(VQ(N)-VQ(N-1)))
!          dltm = VQ(N)-VQ(N-1)
!          dltp = VQ(N+1)-VQ(N)
!
!          IF ( (dltm*dltp) .LE. 0.0 ) THEN
!                dltQ(N) = 0.0
!          ELSE
!                dltma = COF*DABS(dltm)
!                dltpa = COF*DABS(dltp)
!                dltca = DABS(dltc)
!                dltQ(N) = DMIN1( dltca, dltma, dltpa ) * DSIGN(1.0D0,dltc)
!          ENDIF
!
          dltc = deltX1*(deltX2*(VQ(N+1)-VQ(N))+deltX3*(VQ(N)-VQ(N-1)))
          VQmin = DMIN1( VQ(N-1), VQ(N), VQ(N+1) )
          VQmax = DMAX1( VQ(N-1), VQ(N), VQ(N+1) )
          dltm = VQ(N) - VQmin
          dltp = VQmax -VQ(N)

          dltma = COF*dltm
          dltpa = COF*dltp
          dltca = DABS(dltc)
          dltQ(N) = DMIN1( dltca, dltma, dltpa ) * DSIGN(1.0D0,dltc)

       ENDDO
       dltQ(1) = 0.0D0
       dltQ(IXYZ) = 0.0D0
!
!             Compute the values at the left and right edges of cells
!
       DO N = 2, IXYZ-2
          deltX1 = DltXYZ(N)/( DltXYZ(N)+DltXYZ(N+1) )
          deltX2 = 1.0D0/( DltXYZ(N-1)+DltXYZ(N)+DltXYZ(N+1)+DltXYZ(N+2) )
          deltX3 = ( 2.0D0*DltXYZ(N+1)*DltXYZ(N) )/( DltXYZ(N)+DltXYZ(N+1) )
          deltX4 = ( DltXYZ(N-1)+DltXYZ(N) )/( 2.0D0*DltXYZ(N)+DltXYZ(N+1) )
          deltX5 = ( DltXYZ(N+2)+DltXYZ(N+1) )/( 2.0D0*DltXYZ(N+1)+DltXYZ(N) )
          deltX6 = DltXYZ(N)*( DltXYZ(N-1)+DltXYZ(N) )/(2.0D0*DltXYZ(N)+DltXYZ(N+1) )
          deltX7 = DltXYZ(N+1)*( DltXYZ(N+1)+DltXYZ(N+2) )/(DltXYZ(N)+2.0D0*DltXYZ(N+1) )

          qR(N) = VQ(N)+deltX1*(VQ(N+1)-VQ(N))+deltX2*(deltX3*(deltX4-deltX5)*(VQ(N+1)-VQ(N))- &
                  deltX6*dltQ(N+1)+deltX7*dltQ(N) )
          VQmin = DMIN1( VQ(N-1), VQ(N), VQ(N+1) )
          VQmax = DMAX1( VQ(N-1), VQ(N), VQ(N+1) )
          qR(N) = DMAX1( VQmin, qR(N) )
          qR(N) = DMIN1( VQmax, qR(N) )
       ENDDO
       qR(1) = 0.5D0*(VQ(1)+VQ(2))
!       qR(1) = VQ(1)
       qR(IXYZ) = VQ(IXYZ)
       qR(IXYZ-1) = 0.5D0*(VQ(IXYZ-1)+VQ(IXYZ))
!       qR(IXYZ-1) = VQ(IXYZ)

       DO N = 2, IXYZ
          qL(N) = qR(N-1)
       ENDDO
       qL(1) = VQ(1)
!
!       IF ( IndST.EQ.2 ) THEN
!            qL(1) = DMAX1( EPSM, qL(1) )
!            qR(1) = DMAX1( EPSM, qR(1) )
!       ENDIF
!
!          Dai Fushan, Based on Huynh H. T., 1995 (A Piecwise-Parabolic
!          Dual-Mesh Methold for the  
!                                                  Eulaer Equations.
!                                                  AIAA-95-1739-CP)
!
       DO N = 3, IXYZ-2

          pjp2 = (6.0D0*(qR(N)+qL(N))-12.0D0*VQ(N))/DltXYZ(N)**2

          tmpa = 6.0D0*VQ(N)/(DltXYZ(N)**2 *( DltXYZ(N-1)+DltXYZ(N) ) )
          tmpb = 6.0D0*VQ(N-1)/(DltXYZ(N-1)*DltXYZ(N-1)*(DltXYZ(N-1)+DltXYZ(N)))
          tmpc = 6.0D0*qL(N)/(DltXYZ(N-1)*DltXYZ(N))
          pjmp2 = tmpa+tmpb+tmpc

          tmpa = 6.0D0*VQ(N)/( DltXYZ(N)**2 *( DltXYZ(N+1)+DltXYZ(N) ) )
          tmpb = 6.0D0*VQ(N+1)/(DltXYZ(N+1)*DltXYZ(N+1)*(DltXYZ(N+1)+DltXYZ(N)))
          tmpc = 6.0D0*qR(N)/(DltXYZ(N+1)*DltXYZ(N))
          pjpp2 = tmpa+tmpb+tmpc

          aj = (pjmp2-pjp2/ck)*(pjmp2-ck*pjp2)
          bj = (pjpp2-pjp2/ck)*(pjpp2-ck*pjp2)

          abjmax = DMAX1(aj,bj)

          IF ( abjmax .GT. EPSM ) THEN

             VQstr = VQ(N) + 2.0D0*(VQ(N)-qL(N))

             tmpa = 2.0D0+cp2/(DltXYZ(N-1)+DltXYZ(N))
             tmpb = cp2*DltXYZ(N)*DltXYZ(N)/(DltXYZ(N-1)*DltXYZ(N-1)*(DltXYZ(N-1)+DltXYZ(N)))
             tmpc = 1.0D0+cp2*DltXYZ(N)/DltXYZ(N-1)

             VQlac = tmpa*VQ(N)+tmpb*VQ(N-1)-tmpc*qL(N)
!
!        Lin Shian-Jian, 2004: A "Vertical Lagrangian" Finite-Volume Dynamical
!        Core for Global Models. Monthly Weather Review, Vol.132, 2293-2307
!
!               VQstr = VQ(N) + dltQ(N)
!               VQlac = VQ(N) + 0.75D0*(dltQ(N)-dltQ(N-2)) + 0.50D0*dltQ(N)
!
             VQminR = DMIN1( VQ(N), VQstr, VQlac )
             VQmaxR = DMAX1( VQ(N), VQstr, VQlac )

             VQstr = VQ(N) + 2.0D0*(VQ(N)-qR(N))

             tmpa = 2.0D0+cp2/(DltXYZ(N+1)+DltXYZ(N))
             tmpb = cp2*DltXYZ(N)*DltXYZ(N)/(DltXYZ(N+1)*DltXYZ(N+1)*(DltXYZ(N+1)+DltXYZ(N)))
             tmpc = 1.0D0+cp2*DltXYZ(N)/DltXYZ(N+1)

             VQlac = tmpa*VQ(N)+tmpb*VQ(N+1)-tmpc*qR(N)
!
!        Lin Shian-Jian, 2004: A "Vertical Lagrangian" Finite-Volume Dynamical
!        Core for Global Models. Monthly Weather Review, Vol.132, 2293-2307
!
!               VQstr = VQ(N) - dltQ(N)
!               VQlac = VQ(N) + 0.75D0*(dltQ(N+2)-dltQ(N)) - 0.50D0*dltQ(N)
!
             VQminL = DMIN1( VQ(N), VQstr, VQlac )
             VQmaxL = DMAX1( VQ(N), VQstr, VQlac )
             IF ( Istpn.EQ.1 ) THEN
                IF ( (VQ(N-1)-VQ(N-2)).NE.0.0D0 ) THEN
                   Rjm1 = DABS( (VQ(N)-VQ(N-1))/(VQ(N-1)-VQ(N-2)) )
                   IF ( ((Rjm1-cstpn)*(Rjm1-1.0D0/cstpn)).GT.0.0D0 ) THEN
                        VQvlin = VQ(N+1)+( VQ(N+1)-qR(N+1) )
                        tmpa = 2.0D0+2.0D0/(DltXYZ(N+2)+DltXYZ(N+1))
                        tmpb = 2.0D0*DltXYZ(N+1)*DltXYZ(N+1)/(DltXYZ(N+2)*DltXYZ(N+2)*(DltXYZ(N+2)+DltXYZ(N+1)))
                        tmpc = 1.0D0+2.0D0*DltXYZ(N+1)/DltXYZ(N+2)
!
!                        tmpa = 2.0D0+cp2/(DltXYZ(N+2)+DltXYZ(N+1))
!                        tmpb =
!                        cp2*DltXYZ(N+1)*DltXYZ(N+1)/(DltXYZ(N+2)*DltXYZ(N+2)*(DltXYZ(N+2)+DltXYZ(N+1)))
!                        tmpc = 1.0D0+cp2*DltXYZ(N+1)/DltXYZ(N+2)

                        VQvlac = tmpa*VQ(N+1)+tmpb*VQ(N+2)-tmpc*qR(N+1)
                        VQvmin = DMIN1( VQ(N+1), VQvlin, VQvlac )
                        VQvmax = DMAX1( VQ(N+1), VQvlin, VQvlac )
!
                        VQvhfjR = DMIN1( DMAX1( qR(N),VQvmin ), VQvmax )
                    ELSE
                        VQvhfjR = qR(N)
                    ENDIF
                ELSE
                    VQvhfjR = qR(N)
                ENDIF
                qR(N) = DMIN1( DMAX1( VQvhfjR,VQminR), VQmaxR )

             ELSE
                qR(N) = DMIN1( DMAX1( qR(N),VQminR), VQmaxR )
             ENDIF
!
!             VQstr = VQ(N) + 2.0D0*(VQ(N)-qR(N))
!             tmpa = 2.0D0+2.0D0/(DltXYZ(N+1)+DltXYZ(N))
!             tmpb =
!             2.0D0*DltXYZ(N)*DltXYZ(N)/(DltXYZ(N+1)*DltXYZ(N+1)*(DltXYZ(N+1)+DltXYZ(N)))
!             tmpc = 1.0D0+2.0D0*DltXYZ(N)/DltXYZ(N+1)
!             VQlac = tmpa*VQ(N)+tmpb*VQ(N+1)-tmpc*qR(N)
!             VQminL = DMIN1( VQ(N), VQstr, VQlac )
!             VQmaxL = DMAX1( VQ(N), VQstr, VQlac )
!
             IF ( Istpn.EQ.1 ) THEN
                IF ( (VQ(N+1)-VQ(N)).NE.0.0D0 ) THEN
                   Rjp1 = DABS( (VQ(N+2)-VQ(N+1))/(VQ(N+1)-VQ(N)) )

                   IF ( ((Rjp1-cstpn)*(Rjp1-1.0D0/cstpn)).GT.0.0D0 ) THEN
                        VQvlin = VQ(N-1)+( VQ(N-1)-qL(N-1) )
                        tmpa = 2.0D0+2.0D0/(DltXYZ(N-2)+DltXYZ(N-1))
                        tmpb = 2.0D0*DltXYZ(N-1)*DltXYZ(N-1)/(DltXYZ(N-2)*DltXYZ(N-2)*(DltXYZ(N-2)+DltXYZ(N-1)))
                        tmpc = 1.0D0+2.0D0*DltXYZ(N-1)/DltXYZ(N-2)
!
!                        tmpa = 2.0D0+cp2/(DltXYZ(N-2)+DltXYZ(N-1))
!                        tmpb =
!                        cp2*DltXYZ(N-1)*DltXYZ(N-1)/(DltXYZ(N-2)*DltXYZ(N-2)*(DltXYZ(N-2)+DltXYZ(N-1)))
!                        tmpc = 1.0D0+cp2*DltXYZ(N-1)/DltXYZ(N-2)
                        VQvlac = tmpa*VQ(N-1)+tmpb*VQ(N-2)-tmpc*qL(N-1)
                        VQvmin = DMIN1( VQ(N-1), VQvlin, VQvlac )
                        VQvmax = DMAX1( VQ(N-1), VQvlin, VQvlac )

                        VQvhfjL = DMIN1( DMAX1( qL(N),VQvmin ), VQvmax )
                   ELSE
                        VQvhfjL = qL(N)
                   ENDIF
                ELSE
                   VQvhfjL = qL(N)
                ENDIF
                qL(N) = DMIN1( DMAX1( VQvhfjL,VQminL), VQmaxL )

             ELSE
                qL(N) = DMIN1( DMAX1( qL(N),VQminL), VQmaxL )
             ENDIF

          ENDIF
!
       ENDDO
       pjp2 = (6.0D0*(qR(2)+qL(2))-12.0D0*VQ(2))/DltXYZ(2)**2
       tmpa = 6.0D0*VQ(2)/(DltXYZ(2)**2 *( DltXYZ(1)+DltXYZ(2) ) )
       tmpb = 6.0D0*VQ(1)/(DltXYZ(1)*DltXYZ(1)*(DltXYZ(1)+DltXYZ(2)))
       tmpc = 6.0D0*qL(2)/(DltXYZ(1)*DltXYZ(2))
       pjmp2 = tmpa+tmpb+tmpc
       tmpa = 6.0D0*VQ(2)/( DltXYZ(2)**2 *( DltXYZ(3)+DltXYZ(2) ) )
       tmpb = 6.0D0*VQ(3)/(DltXYZ(3)*DltXYZ(3)*(DltXYZ(3)+DltXYZ(2)))
       tmpc = 6.0D0*qR(2)/(DltXYZ(3)*DltXYZ(2))
       pjpp2 = tmpa+tmpb+tmpc
       aj = (pjmp2-pjp2/ck)*(pjmp2-ck*pjp2)
       bj = (pjpp2-pjp2/ck)*(pjpp2-ck*pjp2)

       abjmax = DMAX1(aj,bj)

       IF ( abjmax .GT. EPSM ) THEN

            VQstr = VQ(2) + 2.0D0*(VQ(2)-qL(2))

            tmpa = 2.0D0+cp2/(DltXYZ(1)+DltXYZ(2))
            tmpb = cp2*DltXYZ(2)*DltXYZ(2)/(DltXYZ(1)*DltXYZ(1)*(DltXYZ(1)+DltXYZ(2)))
            tmpc = 1.0D0+cp2*DltXYZ(2)/DltXYZ(1)

            VQlac = tmpa*VQ(2)+tmpb*VQ(1)-tmpc*qL(2)

            VQminR = DMIN1( VQ(2), VQstr, VQlac )
            VQmaxR = DMAX1( VQ(2), VQstr, VQlac )
!!
!!            qR(2) = DMIN1( DMAX1( qR(2),VQminR), VQmaxR )
!!
            VQstr = VQ(2) + 2.0D0*(VQ(2)-qR(2))

            tmpa = 2.0D0+cp2/(DltXYZ(3)+DltXYZ(2))
            tmpb = cp2*DltXYZ(2)*DltXYZ(2)/(DltXYZ(3)*DltXYZ(3)*(DltXYZ(3)+DltXYZ(2)))
            tmpc = 1.0D0+cp2*DltXYZ(2)/DltXYZ(3)

            VQlac = tmpa*VQ(2)+tmpb*VQ(3)-tmpc*qR(2)

            VQminL = DMIN1( VQ(2), VQstr, VQlac )
            VQmaxL = DMAX1( VQ(2), VQstr, VQlac )

            qL(2) = DMIN1( DMAX1( qL(2),VQminL), VQmaxL )

            qR(2) = DMIN1( DMAX1( qR(2),VQminR), VQmaxR )

       ENDIF
!!       qR(1) = qL(2)
!!!       qL(1) = qL(2)
!!
!!       VQstr = VQ(1) + 2.0D0*(VQ(1)-qR(1))
!!!       tmpa = 2.0D0+2.0D0/(DltXYZ(2)+DltXYZ(1))
!!       tmpb =
!2.0D0*DltXYZ(1)*DltXYZ(1)/(DltXYZ(2)*DltXYZ(2)*(DltXYZ(2)+DltXYZ(1)))
!!       tmpc = 1.0D0+2.0D0*DltXYZ(2)/DltXYZ(3)
!!       tmpa = 2.0D0+1.0D0/(DltXYZ(2)+DltXYZ(1))
!!       tmpb =
!1.0D0*DltXYZ(1)*DltXYZ(1)/(DltXYZ(2)*DltXYZ(2)*(DltXYZ(2)+DltXYZ(1)))
!!       tmpc = 1.0D0+1.0D0*DltXYZ(1)/DltXYZ(2)
!!       VQlac = tmpa*VQ(1)+tmpb*VQ(2)-tmpc*qR(1)
!!       VQminL = DMIN1( VQ(1), VQstr, VQlac )
!!       VQmaxL = DMAX1( VQ(1), VQstr, VQlac )
!!
!!       qL(1) = DMIN1( DMAX1( qL(1),VQminL), VQmaxL )                
!!
       pjp2 = (6.0D0*(qR(IXYZ-1)+qL(IXYZ-1))-12.0D0*VQ(IXYZ-1))/DltXYZ(IXYZ-1)**2
       tmpa = 6.0D0*VQ(IXYZ-1)/(DltXYZ(IXYZ-1)**2 *( DltXYZ(IXYZ-2)+DltXYZ(IXYZ-1) ) )
       tmpb = 6.0D0*VQ(IXYZ-2)/(DltXYZ(IXYZ-2)*DltXYZ(IXYZ-2)*(DltXYZ(IXYZ-2)+DltXYZ(IXYZ-1)))
       tmpc = 6.0D0*qL(IXYZ-1)/(DltXYZ(IXYZ-2)*DltXYZ(IXYZ-1))
       pjmp2 = tmpa+tmpb+tmpc
       tmpa = 6.0D0*VQ(IXYZ-1)/( DltXYZ(IXYZ-1)**2 *(DltXYZ(IXYZ)+DltXYZ(IXYZ-1) ) )
       tmpb = 6.0D0*VQ(IXYZ)/(DltXYZ(IXYZ)*DltXYZ(IXYZ)*(DltXYZ(IXYZ)+DltXYZ(IXYZ-1)))
       tmpc = 6.0D0*qR(IXYZ-1)/(DltXYZ(IXYZ)*DltXYZ(IXYZ-1))
       pjpp2 = tmpa+tmpb+tmpc
       aj = (pjmp2-pjp2/ck)*(pjmp2-ck*pjp2)
       bj = (pjpp2-pjp2/ck)*(pjpp2-ck*pjp2)
       abjmax = DMAX1(aj,bj)

       IF ( abjmax .GT. EPSM ) THEN

            VQstr = VQ(IXYZ-1) + 2.0D0*(VQ(IXYZ-1)-qL(IXYZ-1))

            tmpa = 2.0D0+cp2/(DltXYZ(IXYZ-2)+DltXYZ(IXYZ-1))
            tmpb = cp2*DltXYZ(IXYZ-1)*DltXYZ(IXYZ-1)/(DltXYZ(IXYZ-2)*DltXYZ(IXYZ-2)*(DltXYZ(IXYZ-2)+DltXYZ(IXYZ-1)))
            tmpc = 1.0D0+cp2*DltXYZ(IXYZ-1)/DltXYZ(IXYZ-2)

            VQlac = tmpa*VQ(IXYZ-1)+tmpb*VQ(IXYZ-2)-tmpc*qL(IXYZ-1)

            VQminR = DMIN1( VQ(IXYZ-1), VQstr, VQlac )
            VQmaxR = DMAX1( VQ(IXYZ-1), VQstr, VQlac )
!!
!!            qR(IXYZ-1) = DMIN1( DMAX1( qR(IXYZ-1),VQminR), VQmaxR )
!!
            VQstr = VQ(IXYZ-1) + 2.0D0*(VQ(IXYZ-1)-qR(IXYZ-1))

            tmpa = 2.0D0+cp2/(DltXYZ(IXYZ)+DltXYZ(IXYZ-1))
            tmpb = cp2*DltXYZ(IXYZ-1)*DltXYZ(IXYZ-1)/(DltXYZ(IXYZ)*DltXYZ(IXYZ)*(DltXYZ(IXYZ)+DltXYZ(IXYZ-1)))
            tmpc = 1.0D0+cp2*DltXYZ(IXYZ-1)/DltXYZ(IXYZ)

            VQlac = tmpa*VQ(IXYZ-1)+tmpb*VQ(IXYZ)-tmpc*qR(IXYZ-1)

            VQminL = DMIN1( VQ(IXYZ-1), VQstr, VQlac )
            VQmaxL = DMAX1( VQ(IXYZ-1), VQstr, VQlac )
!!
            qL(IXYZ-1) = DMIN1( DMAX1( qL(IXYZ-1),VQminL), VQmaxL )

            qR(IXYZ-1) = DMIN1( DMAX1( qR(IXYZ-1),VQminR), VQmaxR )

       ENDIF
!!       qR(IXYZ) = qR(IXYZ-1)
!!       qL(IXYZ) = qR(IXYZ-1)
!
!!       VQstr = VQ(IXYZ) + 2.0D0*(VQ(IXYZ)-qL(IXYZ))
!!       tmpa = 2.0D0+2.0D0/(DltXYZ(IXYZ-1)+DltXYZ(IXYZ))
!!       tmpb = 2.0D0*DltXYZ(IXYZ)*DltXYZ(IXYZ)/(DltXYZ(IXYZ-1)*DltXYZ(IXYZ-1)*(DltXYZ(IXYZ-1)+DltXYZ(IXYZ)))
!!       tmpc = 1.0D0+2.0D0*DltXYZ(IXYZ)/DltXYZ(IXYZ-1)
!!       tmpa = 2.0D0+1.0D0/(DltXYZ(IXYZ-1)+DltXYZ(IXYZ))
!!       tmpb = 1.0D0*DltXYZ(IXYZ)*DltXYZ(IXYZ)/(DltXYZ(IXYZ-1)*DltXYZ(IXYZ-1)*(DltXYZ(IXYZ-1)+DltXYZ(IXYZ)))
!!       tmpc = 1.0D0+1.0D0*DltXYZ(IXYZ)/DltXYZ(IXYZ-1)
!!       VQlac = tmpa*VQ(IXYZ)+tmpb*VQ(IXYZ-1)-tmpc*qL(IXYZ)
!!       VQminR = DMIN1( VQ(IXYZ), VQstr, VQlac )
!!       VQminR = DMAX1( EPSM, VQminR )
!!
!!       VQmaxR = DMAX1( VQ(IXYZ), VQstr, VQlac )
!!       qR(IXYZ) = DMIN1( DMAX1( qR(IXYZ),VQminR), VQmaxR )
!! 
       DO N = 1, IXYZ
          alphaR(N) = qR(N)-VQ(N)
          alphaL(N) = qL(N)-VQ(N)
       ENDDO

!
!            To Compute the Flux Through the Cell
!
       DO N = 1, IXYZ

          IF ( N.EQ.1 ) THEN

               IF ( KXYZ.EQ.1 ) THEN
                    uhafR = Suvw(N+1)
                    uhafL = Suvw(N)
               ELSE IF ( KXYZ.EQ.2 ) THEN
                    uhafR = Suvw(N)
                    uhafL = 0.0D0
               ELSE IF ( KXYZ.EQ.3 ) THEN
                    uhafR = Suvw(N)
                    uhafL = 0.0D0
               ENDIF
               IF( uhafR .GT. 0.0D0 ) THEN

                   FMIU = DMIN1( 1.0D0, DABS(uhafR*Dtsc/DltXYZ(N)) )

                   flxR(N) = uhafR*( VQ(N)+alphaR(N)+0.5D0*FMIU*(alphaL(N)-alphaR(N)- &
                              (3.0D0-2.0D0*FMIU)*(alphaL(N)+alphaR(N))) )

               ELSE IF ( uhafR.EQ.0.0D0 ) THEN

                   flxR(N) = 0.0D0
               ELSE

                   FMIU = DMIN1( 1.0D0, DABS(uhafR*Dtsc/DltXYZ(N+1)) )

                   flxR(N) = uhafR*(VQ(N)+alphaR(N)+0.5D0*FMIU*(alphaR(N+1)-alphaL(N+1)- &
                              (3.0D0-2.0D0*FMIU)*(alphaR(N+1)+alphaL(N+1))))

               ENDIF
               IF ( uhafL .GE. 0.0D0 ) THEN

                    flxL(N) = 0.0D0
               ELSE

                    FMIU = DMIN1( 1.0D0, DABS(uhafL*Dtsc/DltXYZ(N) ) )

                    flxL(N) = uhafL*(VQ(N)+alphaL(N)+0.5D0*FMIU*(alphaR(N)-alphaL(N)- &
                              (3.0D0-2.0D0*FMIU)*(alphaR(N)+alphaL(N))))

               ENDIF

          ELSE IF ( N.GE.2 .AND. N.LE.(IXYZ-1) ) THEN

               IF ( KXYZ.EQ.1 ) THEN
                    uhafR = Suvw(N+1)
                    uhafL = Suvw(N)
               ELSE IF ( KXYZ.EQ.2 ) THEN
                    uhafR = Suvw(N)
                    uhafL = Suvw(N-1)
               ELSE IF ( KXYZ.EQ.3 ) THEN
                    uhafR = Suvw(N)
                    uhafL = Suvw(N-1)
               ENDIF
               IF( uhafR .GT. 0.0D0 ) THEN

                   FMIU = DMIN1( 1.0D0, DABS(uhafR*Dtsc/DltXYZ(N)) )

                   flxR(N) = uhafR*( VQ(N)+alphaR(N)+0.5D0*FMIU*(alphaL(N)-alphaR(N)- &
                              (3.0D0-2.0D0*FMIU)*(alphaL(N)+alphaR(N))) )

               ELSE IF ( uhafR.EQ.0.0D0 ) THEN

                   flxR(N) = 0.0D0
               ELSE

                   FMIU = DMIN1( 1.0D0, DABS(uhafR*Dtsc/DltXYZ(N+1)) )

                   flxR(N) = uhafR*( VQ(N)+alphaR(N)+0.5D0*FMIU*(alphaR(N+1)-alphaL(N+1)- &
                              (3.0D0-2.0D0*FMIU)*(alphaR(N+1)+alphaL(N+1))) )

               ENDIF

               IF( uhafL .GT. 0.0D0 ) THEN

                   FMIU = DMIN1( 1.0D0, DABS(uhafL*Dtsc/DltXYZ(N-1)) )

                   flxL(N) = uhafL*( VQ(N-1)+alphaR(N-1)+0.5D0*FMIU*(alphaL(N-1)-alphaR(N-1)- &
                             (3.0D0-2.0D0*FMIU)*(alphaL(N-1)+alphaR(N-1))))

               ELSE IF ( uhafL.EQ.0.0D0 ) THEN
                   flxL(N) = 0.0D0

               ELSE
                   FMIU = DMIN1( 1.0D0, DABS(uhafL*Dtsc/DltXYZ(N)) )

                   flxL(N) = uhafL*( VQ(N)+alphaL(N)+0.5D0*FMIU*(alphaR(N)-alphaL(N)- &
                              (3.0D0-2.0D0*FMIU)*(alphaR(N)+alphaL(N))))

               ENDIF

          ELSE IF ( N.EQ.IXYZ ) THEN

               IF ( KXYZ.EQ.1 ) THEN
                    uhafR = 0.0D0
                    uhafL = Suvw(N)
               ELSE IF ( KXYZ.EQ.2 ) THEN
                    uhafR = Suvw(N)
                    uhafL = Suvw(N-1)
               ELSE IF ( KXYZ.EQ.3 ) THEN
                    uhafR = 0.0D0
                    uhafL = Suvw(N-1)
               ENDIF
               IF ( uhafR.GT.0.0D0 ) THEN

                   FMIU = DMIN1( 1.0D0, DABS(uhafR*Dtsc/DltXYZ(N)) )

                   flxR(N) = uhafR*( VQ(N)+alphaR(N)+0.5D0*FMIU*(alphaL(N)-alphaR(N)- &
                              (3.0D0-2.0D0*FMIU)*(alphaL(N)+alphaR(N))) )

               ELSE
                   flxR(N) = 0.0D0
               ENDIF

               IF( uhafL .GT. 0.0D0 ) THEN

                   FMIU = DMIN1( 1.0D0, DABS(uhafL*Dtsc/DltXYZ(N-1)) )

                   flxL(N) = uhafL*( VQ(N-1)+alphaR(N-1)+0.5D0*FMIU*(alphaL(N-1)-alphaR(N-1)- &
                              (3.0D0-2.0D0*FMIU)*(alphaL(N-1)+alphaR(N-1))) )

               ELSE IF ( uhafL.EQ.0.0D0 ) THEN

                   flxL(N) = 0.0D0
               ELSE

                   FMIU = DMIN1( 1.0D0, DABS(uhafL*Dtsc/DltXYZ(N)) )

                   flxL(N) = uhafL*( VQ(N)+alphaL(N)+0.5D0*FMIU*(alphaR(N)-alphaL(N)- &
                              (3.0D0-2.0D0*FMIU)*(alphaR(N)+alphaL(N))) )

               ENDIF

          ENDIF

       ENDDO

       DEALLOCATE( dltQ )
       DEALLOCATE( qR )
       DEALLOCATE( qL )
       DEALLOCATE( alphaR )
       DEALLOCATE( alphaL )

       RETURN

      END SUBROUTINE QRLFLX_PPM_NQGRID


end module ppm_mod
