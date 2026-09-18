!
!      Open Boundary Condition Module
!
module prm_mod

      use precision_mod
      use param_mod
      use pconst_mod


      contains

      SUBROUTINE ADVEC_PRM_UV( UPRM, VPRM, WPRM, VARTMP, DLTMP )

#include <def-undef.h>
      use   precision_mod
      use   param_mod,  only: IMT, JMT, KM, I, J, K
      use   pconst_mod, only: OUX, OUY, ZKT, ZKP, DTC

      IMPLICIT NONE

      REAL(r8) :: UPRM(IMT,JMT,KM)
      REAL(r8) :: VPRM(IMT,JMT,KM)
      REAL(r8) :: WPRM(IMT,JMT,KM)

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
               UVWhaf(I) = 0.50D0*(UPRM(I-1,J,K)+UPRM(I,J,K))
            ENDDO
            UVWhaf(1) = 0.0D0

            INDXYZ = 1
!
!            To Compute the Flux Through the Cell

            CALL QRLFLX_PRM_NQGRID( INDXYZ, IMT, VQC, DXYZ, UVWhaf, DTC, NST, flxqR, flxqL )
!
            DO I = 1, IMT

               IF ( I.EQ.1 ) THEN

                  uhafR = UVWhaf(I+1)
                  uhafL = UVWhaf(I)

                  CMIU = OUX(J)

                  FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                  FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUX(J)

               ELSE IF ( I.GE.2 .AND. I.LE.(IMT-1) ) THEN

                  uhafR = UVWhaf(I+1)
                  uhafL = UVWhaf(I)

                  CMIU = OUX(J)
                  FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                  FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUX(J)

               ELSE IF ( I.EQ.IMT ) THEN

                  uhafR = 0.0D0
                  uhafL = UVWhaf(I)

                  CMIU = OUX(J)
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
               UVWhaf(K) = -WPRM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0

            INDXYZ = 3
!
            CALL QRLFLX_PRM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTC, NST, flxqR, flxqL )
!
!            To Compute the Flux Through the Cell
!

            DO K = 1, KM

               CMIU = 1.0D0/DXYZ(K)
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
               UVWhaf(J) = 0.50D0*(VPRM(I,J,K)+VPRM(I,J+1,K))
            ENDDO
            UVWhaf(JMT) = 0.0D0

            INDXYZ = 2
!

            CALL QRLFLX_PRM_NQGRID( INDXYZ, JMT, VQC, DXYZ, UVWhaf, DTC, NST, flxqR, flxqL )
!
!            To Compute the Flux Through the Cell
            DO J = 1, JMT

               IF ( J.EQ.1 ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = 0.0D0

                    CMIU = OUY(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))
                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUY(J)

               ELSE IF ( J.GE.2 .AND. J.LE.(JMT-1) ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = UVWhaf(J-1)

                    CMIU = OUY(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))

                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTC*(uhafR-uhafL)*OUY(J)
               ELSE IF ( J.EQ.JMT ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = UVWhaf(J-1)

                    CMIU = OUY(J)

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
               UVWhaf(K) = -WPRM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0
            INDXYZ = 3
!
            CALL QRLFLX_PRM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTC, NST, flxqR, flxqL )
!
!            To Compute the Flux Through the Cell
!
            DO K = 1, KM

               CMIU = 1.0D0/DXYZ(K)
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

      END SUBROUTINE ADVEC_PRM_UV

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE ADVEC_PRM_TS( UPRM, VPRM, WPRM, VARTMP, NST, DLTMP, AXTMP, AYTMP, AZTMP )

#include <def-undef.h>
use   precision_mod
use   param_mod,  only: IMT, JMT, KM, KMP1, I, J, K
use   pconst_mod, only: OTX, DYT, ZKT, ZKP, DTS

      IMPLICIT NONE

      REAL(r8) :: UPRM(IMT,JMT,KM)
      REAL(r8) :: VPRM(IMT,JMT,KM)
      REAL(r8) :: WPRM(IMT,JMT,KMP1)
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
               UVWhaf(I) = 0.50D0*(UPRM(I,J-1,K)+UPRM(I,J,K))
            ENDDO

            INDXYZ = 1

            CALL QRLFLX_PRM_NQGRID( INDXYZ, IMT, VQC, DXYZ, UVWhaf, DTS, NST, flxqR, flxqL )
! 
!
!            To Compute the Flux Through the Cell
!
            DO I = 1, IMT

              IF ( I.EQ.1 ) THEN

                   uhafR = UVWhaf(I+1)
                   uhafL = UVWhaf(I)

                   CMIU = OTX(J)
                   FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                   FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)*OTX(J)

!                  AXTMP(I,J,K) = ( FCAVE(I,J,K)-VARTMP(I,J,K) )/DTS
                   AXTMP(I,J,K) = ( FPAVE(I,J,K)-VARTMP(I,J,K) )/DTS

              ELSE IF ( I.GE.2 .AND. I.LE.(IMT-1) ) THEN

                   uhafR = UVWhaf(I+1)
                   uhafL = UVWhaf(I)

                   CMIU = OTX(J)
                   FPAVE(I,J,K) = VARTMP(I,J,K) + CMIU*(flxqL(I)-flxqR(I))

                   FCAVE(I,J,K) = FPAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)*OTX(J)
!
                   AXTMP(I,J,K) = ( FPAVE(I,J,K)-VARTMP(I,J,K) )/DTS

               ELSE IF ( I.EQ.IMT ) THEN

                   uhafR = 0.0D0
                   uhafL = UVWhaf(I)

                   CMIU = OTX(J)
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
               UVWhaf(K) = -WPRM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0

            INDXYZ = 3

            CALL QRLFLX_PRM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTS, NST, flxqR, flxqL )
! 
!
!            To Compute the Flux Through the Cell
!
            DO K = 1, KM

                  CMIU = 1.0D0/DXYZ(K)
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
               UVWhaf(J) = 0.50D0*(VPRM(I,J,K)+VPRM(I+1,J,K))
            ENDDO

            INDXYZ = 2

            CALL QRLFLX_PRM_NQGRID( INDXYZ, JMT, VQC, DXYZ, UVWhaf, DTS, NST, flxqR, flxqL )
!
!
!            To Compute the Flux Through the Cell
!
            DO J = 1, JMT

               IF ( J.EQ.1 ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = 0.0D0

                    CMIU = 1.0D0/DYT(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))

                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)/DYT(J)
!
!                   AYTMP(I,J,K) = ( FEAVE(I,J,K) - FPAVE(I,J,K) )/DTS
                    AYTMP(I,J,K) = ( FDAVE(I,J,K) - FPAVE(I,J,K) )/DTS

               ELSE IF ( J.GE.2 .AND. J.LE.(JMT-1) ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = UVWhaf(J-1)
!
                    CMIU = 1.0D0/DYT(J)

                    FDAVE(I,J,K) = FPAVE(I,J,K) + CMIU*(flxqL(J)-flxqR(J))

                    FEAVE(I,J,K) = FDAVE(I,J,K) + VARTMP(I,J,K)*DTS*(uhafR-uhafL)/DYT(J)
!
!                   AYTMP(I,J,K) = ( FEAVE(I,J,K) - FPAVE(I,J,K) )/DTS
                    AYTMP(I,J,K) = ( FDAVE(I,J,K) - FPAVE(I,J,K) )/DTS

               ELSE IF ( J.EQ.JMT ) THEN

                    uhafR = UVWhaf(J)
                    uhafL = UVWhaf(J-1)
!
                    CMIU = 1.0D0/DYT(J)

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
              UVWhaf(K) = -WPRM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0

            INDXYZ = 3
!
            CALL QRLFLX_PRM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTS, NST, flxqR, flxqL )
! 
!
!            To Compute the Flux Through the Cell
!
            DO K = 1, KM

               CMIU = 1.0D0/DXYZ(K)
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

      END SUBROUTINE ADVEC_PRM_TS
!
!-----------------------------------------------------------------------------------------------------------
!
      SUBROUTINE QRLFLX_PRM_NQGRID( KXYZ, IXYZ, VQ, DltXYZ, Suvw, Dtsc, IndST, flxR, flxL)

       use   precision_mod

       IMPLICIT NONE

       INTEGER :: KXYZ
       INTEGER :: IXYZ
       REAL(r8) :: VQ(IXYZ), DltXYZ(IXYZ)
       REAL(r8) :: Suvw(IXYZ)
       REAL(r8) :: Dtsc
       INTEGER :: IndST
       REAL(r8) :: flxR(IXYZ), flxL(IXYZ)

       REAL(r8), ALLOCATABLE :: dltQ(:)
       REAL(r8), ALLOCATABLE :: qR(:), qL(:)

       REAL(r8) :: COF

       REAL(r8) :: deltX1, deltX2, deltX3, deltX4
       REAL(r8) :: deltX5, deltX6, deltX7

       REAL(r8) :: dltc, dltm, dltp
       REAL(r8) :: dltca, dltma, dltpa
       REAL(r8) :: VQmin, VQmax

       INTEGER ::  N

       REAL(r8) :: uhafR, uhafL

       REAL(r8) :: FMIU
       REAL(r8) :: Atmp, Btmp
       REAL(r8) :: Ctmp, Dtmp
       REAL(r8) :: Gama, Beta

       REAL(r8) :: Fmin, Fmax, Fave

       REAL(r8) :: EPStmp

       IF ( ALLOCATED( dltQ ) ) DEALLOCATE( dltQ )
       ALLOCATE( dltQ(IXYZ) )

       IF ( ALLOCATED( qR ) ) DEALLOCATE( qR )
       ALLOCATE( qR(IXYZ) )

       IF ( ALLOCATED( qL )) DEALLOCATE( qL )
       ALLOCATE( qL(IXYZ) )

       COF = 3.0D0
       EPStmp = 1.0D-20
!
       DO N = 2, IXYZ-1

          deltX1 = DltXYZ(N)/( DltXYZ(N-1)+DltXYZ(N)+DltXYZ(N+1) )
          deltX2 = ( 2.0D0*DltXYZ(N-1)+DltXYZ(N) )/( DltXYZ(N+1)+DltXYZ(N) )
          deltX3 = ( DltXYZ(N)+2.0D0*DltXYZ(N+1) )/( DltXYZ(N-1)+DltXYZ(N) )

          dltc = deltX1*(deltX2*(VQ(N+1)-VQ(N))+deltX3*(VQ(N)-VQ(N-1)))
          dltm = VQ(N)-VQ(N-1)
          dltp = VQ(N+1)-VQ(N)

          IF ( (dltm*dltp) .LE. 0.0 ) THEN
                dltQ(N) = 0.0
          ELSE
                dltma = COF*DABS(dltm)
                dltpa = COF*DABS(dltp)
                dltca = DABS(dltc)
                dltQ(N) = DMIN1( dltca, dltma, dltpa ) * DSIGN(1.0D0,dltc)
          ENDIF

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

          qR(N) =VQ(N)+deltX1*(VQ(N+1)-VQ(N))+deltX2*(deltX3*(deltX4-deltX5)*(VQ(N+1)-VQ(N))- &
                  deltX6*dltQ(N+1)+deltX7*dltQ(N) )

          VQmin = DMIN1( VQ(N-1), VQ(N), VQ(N+1) )
          VQmax = DMAX1( VQ(N-1), VQ(N), VQ(N+1) )
          qR(N) = DMAX1( VQmin, qR(N) )
          qR(N) = DMIN1( VQmax, qR(N) )
       ENDDO
       qR(1) = 0.5D0*(VQ(1)+VQ(2))
       qR(IXYZ) = VQ(IXYZ)
       qR(IXYZ-1) = 0.5D0*(VQ(IXYZ-1)+VQ(IXYZ))
       DO N = 2, IXYZ
          qL(N) = qR(N-1)
       ENDDO
       qL(1) = VQ(1)
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

!                   FMIU = DMIN1( 1.0D0, DABS(uhafR*Dtsc/DltXYZ(N)) )

                   FMIU = uhafR*Dtsc

                   Atmp = qR(N)

                   Fmin = DMIN1( qR(N), qL(N) )
                   Fmax = DMAX1( qR(N), qL(N) )
                   Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                   Ctmp = DABS( qR(N)-Fave )
                   Dtmp = DABS( Fave-qL(N) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0 )/(-DltXYZ(N))
                   Btmp = ( Gama*Fave-qR(N) )/(-DltXYZ(N))

                   flxR(N) = ( Atmp*FMIU - Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU )

               ELSE IF ( uhafR.EQ.0.0D0 ) THEN

                   flxR(N) = 0.0D0
               ELSE

!                   FMIU = DMIN1( 1.0D0, DABS(-uhafR*Dtsc/DltXYZ(N+1)) )

                   FMIU = -uhafR*Dtsc

                   Atmp = qR(N)

                   Fmin = DMIN1( qR(N+1), qL(N+1) )
                   Fmax = DMAX1( qR(N+1), qL(N+1) )
                   Fave = DMAX1( DMIN1( VQ(N+1), Fmax ), Fmin )
                   Ctmp = DABS( qR(N)-Fave )
                   Dtmp = DABS( Fave-qR(N+1) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0 )/DltXYZ(N+1)
                   Btmp = ( Gama*Fave-qR(N) )/DltXYZ(N+1)

                   flxR(N) = -1.0D0*( Atmp*FMIU + Btmp*FMIU*FMIU )/ ( 1.0D0+Beta*FMIU )
               ENDIF

               IF ( uhafL .GE. 0.0D0 ) THEN

                    flxL(N) = 0.0D0
               ELSE

!                    FMIU = DMIN1( 1.0D0, DABS(-uhafL*Dtsc/DltXYZ(N) ) )
                    FMIU = -uhafL*Dtsc

                    Atmp = qL(N)

                    Fmin = DMIN1( qL(N), qR(N) )
                    Fmax = DMAX1( qL(N), qR(N) )
                    Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                    Ctmp = DABS( qL(N)-Fave )
                    Dtmp = DABS( Fave-qR(N) )
                    Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                    Beta = ( Gama-1.0D0 )/DltXYZ(N)
                    Btmp = ( Gama*Fave-qL(N) )/DltXYZ(N)

                    flxL(N) = -1.0D0*( Atmp*FMIU+Btmp*FMIU*FMIU )/( 1.0D0+Beta*FMIU )

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

!                   FMIU = DMIN1( 1.0D0, DABS(uhafR*Dtsc/DltXYZ(N)) )
                   FMIU = uhafR*Dtsc

                   Atmp = qR(N)

                   Fmin = DMIN1( qR(N), qL(N) )
                   Fmax = DMAX1( qR(N), qL(N) )
                   Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                   Ctmp = DABS( qR(N)-Fave )
                   Dtmp = DABS( Fave-qL(N) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0)/(-DltXYZ(N))
                   Btmp = ( Gama*Fave-qR(N) )/(-DltXYZ(N))

                   flxR(N) = ( Atmp*FMIU - Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU)

               ELSE IF ( uhafR.EQ.0.0D0 ) THEN

                   flxR(N) = 0.0D0
               ELSE
                   FMIU = -uhafR*Dtsc

                   Atmp = qR(N)

                   Fmin = DMIN1( qR(N+1), qL(N+1) )
                   Fmax = DMAX1( qR(N+1), qL(N+1) )
                   Fave = DMAX1( DMIN1( VQ(N+1), Fmax ), Fmin )
                   Ctmp = DABS( qR(N)-Fave )
                   Dtmp = DABS( Fave-qR(N+1) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0 )/DltXYZ(N+1)
                   Btmp = ( Gama*Fave-qR(N) )/DltXYZ(N+1)

                   flxR(N) = -1.0D0*( Atmp*FMIU + Btmp*FMIU*FMIU )/ ( 1.0D0+Beta*FMIU )

               ENDIF

               IF( uhafL .GT. 0.0D0 ) THEN

!                   FMIU = DMIN1( 1.0D0, DABS(uhafL*Dtsc/DltXYZ(N-1)) )
                   FMIU = uhafL*Dtsc

                   Atmp = qL(N)

                   Fmin = DMIN1( qL(N-1), qR(N-1) )
                   Fmax = DMAX1( qL(N-1), qR(N-1) )
                   Fave = DMAX1( DMIN1( VQ(N-1), Fmax ), Fmin )
                   Ctmp = DABS( qL(N)-Fave )
                   Dtmp = DABS( Fave-qL(N-1) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0 )/(-DltXYZ(N-1))
                   Btmp = ( Gama*Fave-qL(N) )/(-DltXYZ(N-1))

                   flxL(N) = ( Atmp*FMIU-Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU )

               ELSE IF ( uhafL.EQ.0.0D0 ) THEN

                   flxL(N) = 0.0D0
               ELSE

!                   FMIU = DMIN1( 1.0D0, DABS(-uhafL*Dtsc/DltXYZ(N)) )
                   FMIU = -uhafL*Dtsc

                   Atmp = qL(N)

                   Fmin = DMIN1( qL(N), qR(N) )
                   Fmax = DMAX1( qL(N), qR(N) )
                   Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                   Ctmp = DABS( qL(N)-Fave )
                   Dtmp = DABS( Fave-qR(N) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0 )/DltXYZ(N)
                   Btmp = ( Gama*Fave-qL(N) )/DltXYZ(N)

                   flxL(N) = -1.0D0*( Atmp*FMIU+Btmp*FMIU*FMIU )/( 1.0D0+Beta*FMIU )

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

!                   FMIU = DMIN1( 1.0D0, DABS(uhafR*Dtsc/DltXYZ(N)) )
                   FMIU = uhafR*Dtsc

                   Atmp = qR(N)

                   Fmin = DMIN1( qR(N), qL(N) )
                   Fmax = DMAX1( qR(N), qL(N) )
                   Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                   Ctmp = DABS( qR(N)-Fave )
                   Dtmp = DABS( Fave-qL(N) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0 )/(-DltXYZ(N))
                   Btmp = ( Gama*Fave-qR(N) )/(-DltXYZ(N))

                   flxR(N) = ( Atmp*FMIU - Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU)

               ELSE
                   flxR(N) = 0.0D0
               ENDIF

               IF( uhafL .GT. 0.0D0 ) THEN

!                   FMIU = DMIN1( 1.0D0, DABS(uhafL*Dtsc/DltXYZ(N-1)) )
                   FMIU = uhafL*Dtsc

                   Atmp = qL(N)

                   Fmin = DMIN1( qL(N-1), qR(N-1) )
                   Fmax = DMAX1( qL(N-1), qR(N-1) )
                   Fave = DMAX1( DMIN1( VQ(N-1), Fmax ), Fmin )
                   Ctmp = DABS( qL(N)-Fave )
                   Dtmp = DABS( Fave-qL(N-1) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0 )/(-DltXYZ(N-1))
                   Btmp = ( Gama*Fave-qL(N) )/(-DltXYZ(N-1))

                   flxL(N) = ( Atmp*FMIU-Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU )

               ELSE IF ( uhafL.EQ.0.0D0 ) THEN

                   flxL(N) = 0.0D0
               ELSE

!                   FMIU = DMIN1( 1.0D0, DABS(-uhafL*Dtsc/DltXYZ(N)) )
                   FMIU = -uhafL*Dtsc

                   Atmp = qL(N)

                   Fmin = DMIN1( qL(N), qR(N) )
                   Fmax = DMAX1( qL(N), qR(N) )
                   Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                   Ctmp = DABS( qL(N)-Fave )
                   Dtmp = DABS( Fave-qR(N) )
                   Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                   Beta = ( Gama-1.0D0 )/DltXYZ(N)
                   Btmp = ( Gama*Fave-qL(N) )/DltXYZ(N)

                   flxL(N) = -1.0D0*( Atmp*FMIU+Btmp*FMIU*FMIU )/( 1.0D0+Beta*FMIU )

               ENDIF

          ENDIF

       ENDDO
       DEALLOCATE( dltQ )
       DEALLOCATE( qR )
       DEALLOCATE( qL )

       RETURN

      END SUBROUTINE QRLFLX_PRM_NQGRID


end module prm_mod
