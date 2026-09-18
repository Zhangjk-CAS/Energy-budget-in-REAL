SUBROUTINE OBC_H0

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, I, J, nx_proc, ny_proc
  use   pconst_mod, only: VIT, ix, iy
  use   dyn_mod,    only: H0
  use   work_mod

#if (defined OBCDT)
  use openbc_mod
#endif

#if (defined OBTIDE)
  use tide_mod, only: H0TIDE_TPXO7_LOC, THTIDE_TPXO7_LOC, &
      H0TIDEN, THTIDEN, H0TIDES, THTIDES, &
      H0TIDEE, THTIDEE, H0TIDEW, THTIDEW
#endif

#if (defined OBCDT)

  IF ( iy .EQ. 0 ) THEN
#if (defined OBTIDE)

    DO K = 1, Ntide
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, imt
        H0TIDEN(I,K) = H0TIDE_TPXO7_LOC(I,3,K)*VIT(I,3,1)
        THTIDEN(I,K) = THTIDE_TPXO7_LOC(I,3,K)*VIT(I,3,1)
      ENDDO
    ENDDO
#endif

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, imt
      H0BTRN(I,1) = H0_MON_LOC(I,3)*VIT(I,3,1)
      H0BTRN(I,2) = H0_MON_LOC(I,4)*VIT(I,4,1)
    ENDDO

    CALL OBC_2DH0( ID_OBC_N )
  ELSE IF ( iy .EQ. (ny_proc-1) ) THEN
#if (defined OBTIDE)

    DO K = 1, Ntide
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, imt
        H0TIDES(I,K) = H0TIDE_TPXO7_LOC(I,jmt-2,K)*VIT(I,jmt-2,1)
        THTIDES(I,K) = THTIDE_TPXO7_LOC(I,jmt-2,K)*VIT(I,jmt-2,1)
      ENDDO
    ENDDO
#endif

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      H0BTRS(I,1) = H0_MON_LOC(I,JMT-3)*VIT(I,JMT-3,1)
      H0BTRS(I,2) = H0_MON_LOC(I,JMT-2)*VIT(I,JMT-2,1)
    ENDDO

    CALL OBC_2DH0( ID_OBC_S )
  endif

  IF ( ix .EQ. 0 ) THEN
#if (defined OBTIDE)

    DO K = 1, Ntide
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, jmt
        H0TIDEW(J,K) = H0TIDE_TPXO7_LOC(3,J,K)*VIT(3,J,1)
        THTIDEW(J,K) = THTIDE_TPXO7_LOC(3,J,K)*VIT(3,J,1)
      ENDDO
    ENDDO
#endif

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      H0BTRW(1,J) = H0_MON_LOC(3,J)*VIT(3,J,1)
      H0BTRW(2,J) = H0_MON_LOC(4,J)*VIT(4,J,1)
    ENDDO

    CALL OBC_2DH0( ID_OBC_W )
  ELSE IF ( ix .EQ. ( nx_proc-1 ) ) THEN
#if (defined OBTIDE)

    DO K = 1, Ntide
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, jmt
        H0TIDEE(J,K) = H0TIDE_TPXO7_LOC(imt-2,J,K)*VIT(imt-2,J,1)
        THTIDEE(J,K) = THTIDE_TPXO7_LOC(imt-2,J,K)*VIT(imt-2,J,1)
      ENDDO
    ENDDO
#endif

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      H0BTRE(1,J) = H0_MON_LOC(IMT-3,J)*VIT(IMT-3,J,1)
      H0BTRE(2,J) = H0_MON_LOC(IMT-2,J)*VIT(IMT-2,J,1)
    ENDDO

    CALL OBC_2DH0( ID_OBC_E )
  ENDIF

#endif

  RETURN
END SUBROUTINE OBC_H0

SUBROUTINE OBC_UB

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, I, J, nx_proc, ny_proc
  use   pconst_mod, only: VIV, ix, iy
  use   dyn_mod,    only: UB

#if (defined OBCDT)
  use openbc_mod
#endif

#if (defined OBCDT)

  IF ( iy .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, imt
      UBTRN(I,1) = UB_MON_LOC(I,3)*VIV(I,3,1)
      UBTRN(I,2) = UB_MON_LOC(I,4)*VIV(I,4,1)
    ENDDO

    CALL OBC_2DUB( ID_OBC_N )
  ELSE IF ( iy .EQ. (ny_proc-1) ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      UBTRS(I,1) = UB_MON_LOC(I,JMT-4)*VIV(I,JMT-4,1)
      UBTRS(I,2) = UB_MON_LOC(I,JMT-3)*VIV(I,JMT-3,1)
    ENDDO

    CALL OBC_2DUB( ID_OBC_S )
  endif

  IF ( ix .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      UBTRW(1,J) = UB_MON_LOC(4,J)*VIV(4,J,1)
      UBTRW(2,J) = UB_MON_LOC(5,J)*VIV(5,J,1)
    ENDDO

    CALL OBC_2DUB( ID_OBC_W )
  ELSE IF ( ix .EQ. ( nx_proc-1 ) ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      UBTRE(1,J) = UB_MON_LOC(IMT-3,J)*VIV(IMT-3,J,1)
      UBTRE(2,J) = UB_MON_LOC(IMT-2,J)*VIV(IMT-2,J,1)
    ENDDO

    CALL OBC_2DUB( ID_OBC_E )
  ENDIF

#endif

  RETURN
END SUBROUTINE OBC_UB

SUBROUTINE OBC_VB

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, I, J, nx_proc, ny_proc
  use   pconst_mod, only: VIV, ix, iy
  use   dyn_mod,    only: VB

#if (defined OBCDT)
  use openbc_mod
#endif

#if (defined OBCDT)

  IF ( iy .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, imt
      VBTRN(I,1) = VB_MON_LOC(I,3)*VIV(I,3,1)
      VBTRN(I,2) = VB_MON_LOC(I,4)*VIV(I,4,1)
    ENDDO

    CALL OBC_2DVB( ID_OBC_N )
  ELSE IF ( iy .EQ. (ny_proc-1) ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      VBTRS(I,1) = VB_MON_LOC(I,JMT-4)*VIV(I,JMT-4,1)
      VBTRS(I,2) = VB_MON_LOC(I,JMT-3)*VIV(I,JMT-3,1)
    ENDDO

    CALL OBC_2DVB( ID_OBC_S )
  endif

  IF ( ix .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      VBTRW(1,J) = VB_MON_LOC(4,J)*VIV(4,J,1)
      VBTRW(2,J) = VB_MON_LOC(5,J)*VIV(5,J,1)
    ENDDO

    CALL OBC_2DVB( ID_OBC_W )
  ELSE IF ( ix .EQ. ( nx_proc-1 ) ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      VBTRE(1,J) = VB_MON_LOC(IMT-3,J)*VIV(IMT-3,J,1)
      VBTRE(2,J) = VB_MON_LOC(IMT-2,J)*VIV(IMT-2,J,1)
    ENDDO

    CALL OBC_2DVB( ID_OBC_E )
  ENDIF

#endif

  RETURN
END SUBROUTINE OBC_VB

SUBROUTINE OBC_UU

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, KM, I, J, K, nx_proc, ny_proc
  use   pconst_mod, only: VIV, ix, iy
  use   dyn_mod,    only: U

#if (defined OBCDT)
  use openbc_mod
#endif

#if (defined OBCDT)

  IF( iy .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (K,I)
    DO K = 1, KM
      DO I = 1, IMT
        UBCLN(I,1,K) = UU_MON_LOC(I,2,K)*VIV(I,3,K)
        UBCLN(I,2,K) = UU_MON_LOC(I,3,K)*VIV(I,4,K)
      ENDDO
    ENDDO

    CALL OBC_3DU( ID_OBC_N )
  ELSE IF ( iy .eq. (ny_proc-1) ) THEN
    !$OMP PARALLEL DO PRIVATE (K,I)
    DO K = 1, KM
      DO I = 1, IMT
        UBCLS(I,1,K) = UU_MON_LOC(I,JMT-3,K)*VIV(I,JMT-4,K)
        UBCLS(I,2,K) = UU_MON_LOC(I,JMT-2,K)*VIV(I,JMT-3,K)
      ENDDO
    ENDDO

    CALL OBC_3DU( ID_OBC_S )
  endif

  IF( ix .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (K,J)
    DO K = 1, KM
      DO J = 1, JMT
        UBCLW(1,J,K) = UU_MON_LOC(3,J,K)*VIV(4,J,K)
        UBCLW(2,J,K) = UU_MON_LOC(4,J,K)*VIV(5,J,K)
      ENDDO
    ENDDO

    CALL OBC_3DU( ID_OBC_W )
  ELSE IF ( ix .EQ. ( nx_proc-1 ) ) THEN
    !$OMP PARALLEL DO PRIVATE (K,J)
    DO K = 1, KM
      DO J = 1, JMT
        UBCLE(1,J,K) = UU_MON_LOC(IMT-2,J,K)*VIV(IMT-3,J,K)
        UBCLE(2,J,K) = UU_MON_LOC(IMT-1,J,K)*VIV(IMT-2,J,K)
      ENDDO
    ENDDO

    CALL OBC_3DU( ID_OBC_E )
  ENDIF

#endif

  RETURN
END SUBROUTINE OBC_UU

SUBROUTINE OBC_VV

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, KM, I, J, K, nx_proc, ny_proc
  use   pconst_mod, only: VIV, ix, iy
  use   dyn_mod,    only: V

#if (defined OBCDT)
  use openbc_mod
#endif

#if (defined OBCDT)

  IF( iy .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (K,I)
    DO K = 1, KM
      DO I = 1, IMT
        VBCLN(I,1,K) = VV_MON_LOC(I,2,K)*VIV(I,3,K)
        VBCLN(I,2,K) = VV_MON_LOC(I,3,K)*VIV(I,4,K)
      ENDDO
    ENDDO

    CALL OBC_3DV( ID_OBC_N )
  ELSE IF ( iy .eq. (ny_proc-1) ) THEN
    !$OMP PARALLEL DO PRIVATE (K,I)
    DO K = 1, KM
      DO I = 1, IMT
        VBCLS(I,1,K) = VV_MON_LOC(I,JMT-3,K)*VIV(I,JMT-4,K)
        VBCLS(I,2,K) = VV_MON_LOC(I,JMT-2,K)*VIV(I,JMT-3,K)
      ENDDO
    ENDDO

    CALL OBC_3DV( ID_OBC_S )
  endif
  IF( ix .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (K,J)
    DO K = 1, KM
      DO J = 1, JMT
        VBCLW(1,J,K) = VV_MON_LOC(3,J,K)*VIV(4,J,K)
        VBCLW(2,J,K) = VV_MON_LOC(4,J,K)*VIV(5,J,K)
      ENDDO
    ENDDO

    CALL OBC_3DV( ID_OBC_W )
  ELSE IF ( ix .EQ. ( nx_proc-1 ) ) THEN
    !$OMP PARALLEL DO PRIVATE (K,J)
    DO K = 1, KM
      DO J = 1, JMT
        VBCLE(1,J,K) = VV_MON_LOC(IMT-2,J,K)*VIV(IMT-3,J,K)
        VBCLE(2,J,K) = VV_MON_LOC(IMT-1,J,K)*VIV(IMT-2,J,K)
      ENDDO
    ENDDO

    CALL OBC_3DV( ID_OBC_E )
  ENDIF

#endif

  RETURN
END SUBROUTINE OBC_VV

SUBROUTINE OBC_TS

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, KM, I, J, K, nx_proc, ny_proc
  use   pconst_mod, only: VIT, ix, iy
  use   tracer_mod, only: AT

#if ( defined SPMD )
  use msg_mod
#endif

#if (defined OBCDT)
  use   openbc_mod
#endif

  IMPLICIT NONE

#if (defined OBCDT)

  IF ( iy .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (K,I)
    DO K = 1, KM
      DO I = 1, IMT
        TBCLN(I,K) = TT_MON_LOC(I,2,K)*VIT(I,3,K)
        SBCLN(I,K) = SS_MON_LOC(I,2,K)*VIT(I,3,K)
      ENDDO
    ENDDO

    CALL OBC_3DTS( ID_OBC_N )
  ELSE IF ( iy .eq. (ny_proc-1) ) then
    !$OMP PARALLEL DO PRIVATE (K,I)
    DO K = 1, KM
      DO I = 1, IMT
        TBCLS(I,K) = TT_MON_LOC(I,JMT-1,K)*VIT(I,JMT-2,K)
        SBCLS(I,K) = SS_MON_LOC(I,JMT-1,K)*VIT(I,JMT-2,K)
      ENDDO
    ENDDO

    CALL OBC_3DTS( ID_OBC_S )
  endif

  IF ( ix .EQ. 0 ) THEN
    !$OMP PARALLEL DO PRIVATE (K,J)
    DO K = 1, KM
      DO J = 1, JMT
        TBCLW(J,K) = TT_MON_LOC(2,J,K)*VIT(3,J,K)
        SBCLW(J,K) = SS_MON_LOC(2,J,K)*VIT(3,J,K)
      ENDDO
    ENDDO

    CALL OBC_3DTS( ID_OBC_W )
  ELSE IF ( ix .EQ. ( nx_proc-1 ) ) THEN
    !$OMP PARALLEL DO PRIVATE (K,J)
    DO K = 1, KM
      DO J = 1, JMT
        TBCLE(J,K) = TT_MON_LOC(IMT-1,J,K)*VIT(IMT-2,J,K)
        SBCLE(J,K) = SS_MON_LOC(IMT-1,J,K)*VIT(IMT-2,J,K)
      ENDDO
    ENDDO

    CALL OBC_3DTS( ID_OBC_E )
  ENDIF

#endif

  RETURN

END SUBROUTINE OBC_TS

SUBROUTINE OBC_3DU(JEWSN)

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, JET, KM, I, J, K
  use   pconst_mod, only: VIV, DYT, OTX, DTC, DTC2, IY
  use   dyn_mod,    only: U, V, UP, VP

#if ( defined SPMD )
  use msg_mod
#endif

#if (defined OBCDT)
  use   openbc_mod
#endif

  IMPLICIT NONE
  INTEGER :: JEWSN
  REAL(r8) :: Tau_obc
  REAL(r8) :: EPSCF
  REAL(r8) :: CFF1, CFF2
  REAL(r8) :: UDTDX
  REAL(r8) :: DVDTTMP, DVDNMTMP, DVDTGTMP, DVDTG2TMP, CFTG, CFNM
  REAL(r8) :: CMIU
  EPSCF = 1.0D-20

  IF ( JEWSN .EQ. 1 ) THEN
    ! Eastern Boundary Condition
    IF ( ID_OBC_U_E .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For U at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          U(IMT-2,J,K) = U(IMT-3,J,K) * VIV(IMT-3,J,K)
          U(IMT-1,J,K) = U(IMT-2,J,K)
          U(IMT,J,K) = U(IMT-2,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_E .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For U at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          U(IMT-2,J,K) = UBCLE(2,J,K) * VIV(IMT-2,J,K)
          U(IMT-1,J,K) = U(IMT-2,J,K)
          U(IMT,J,K) = U(IMT-2,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_E .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For U at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          UDTDX = U(IMT-2,J,K) * DTC * OTX(J)

          IF ( UDTDX.LE.0.0D0 ) THEN
            U(IMT-2,J,K) = UP(IMT-2,J,K) - UDTDX * ( UBCLE(2,J,K) - UP(IMT-2,J,K) )
          ELSE
            U(IMT-2,J,K) = UP(IMT-2,J,K) - UDTDX * ( UP(IMT-2,J,K) - UP(IMT-3,J,K) )
          ENDIF
          U(IMT-1,J,K) = U(IMT-2,J,K)
          U(IMT,J,K) = U(IMT-2,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_E .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For U at Eastern Edge
      ! (Marchesiello P., 2001, Ocean Modelling 3 (2001) 1-20 )

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 2, JMT-1
          DVDTTMP = UP(IMT-3,J,K) - Unow(IMT-3,J,K)
          DVDNMTMP = Unow(IMT-3,J,K) - Unow(IMT-4,J,K)
          DVDTG2TMP = UP(IMT-3,J+1,K) - UP(IMT-3,J-1,K)

          IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
            DVDTTMP = 0.0D0
            Tau_obc = TAUIN
          ELSE
            Tau_obc = TAUOUT
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            DVDTGTMP = UP(IMT-3,J,K) - UP(IMT-3,J-1,K)
          ELSE
            DVDTGTMP = UP(IMT-3,J+1,K) - UP(IMT-3,J,K)
          ENDIF

          IF ( ID_OBC_2DRad_U .EQ. 0 ) THEN
            DVDTGTMP = 0.0D0
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            CFF1 = DVDNMTMP*DVDNMTMP + &
                DVDTGTMP*DVDTGTMP / ( DYT(J)*DYT(J) ) * &
                ( 0.50D0/OTX(J)+0.50D0*OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP * DYT(J)*DYT(J) / &
                (  ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )*( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
          ELSE
            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP / ( DYT(J+1)*DYT(J+1) ) * &
                (  0.5D0/OTX(J)+0.50D0/OTX(J+1) )*( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP * DYT(J+1)*DYT(J+1) / &
                ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
          ENDIF

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = MAX(CFF2,EPSCF)
          CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
          CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

          IF ( ID_OBC_2DRad_NPO_U .EQ. 1 ) THEN
            CFTG = 0.0D0
          ENDIF

          CMIU = CFF1 + CFNM
          Unow(IMT-2,J,K) = ( CFF1*UP(IMT-2,J,K) + &
              CFNM*Unow(IMT-3,J,K) -  &
              MAX(CFTG,0.0D0)*(CFF1/CFF2)*(UP(IMT-2,J,K)-UP(IMT-2,J-1,K)) - &
              MIN(CFTG,0.0D0)*(CFF1/CFF2)*(UP(IMT-2,J+1,K)-UP(IMT-2,J,K)) )/CMIU

          Unow(IMT-2,J,K) = Unow(IMT-2,J,K) + &
              DTC_OPENBC/Tau_obc * ( ( 0.2*UBCLE(1,J,K)+0.4*UBCLE(2,J,K)+ &
              0.2*UBCLE(2,J-1,K)+0.2*UBCLE(2,J+1,K) )-UP(IMT-2,J,K) ) * CFF1/CMIU
          Unow(IMT-1,J,K) = Unow(IMT-2,J,K)
          Unow(IMT,J,K) = Unow(IMT-2,J,K)
        ENDDO
        do j = 2,jmt-1
          unow(imt-2,j,k) = 0.25d0*unow(imt-2,j-1,k) + 0.50d0*unow(imt-2,j,k) + 0.25d0*unow(imt-2,j+1,k)
        enddo
        Unow(IMT-2,1,K) = Unow(IMT-2,2,K)
        Unow(IMT-1,1,K) = Unow(IMT-1,2,K)
        Unow(IMT,1,K) = Unow(IMT,2,K)
        Unow(IMT-2,JMT,K) = Unow(IMT-2,JMT-1,K)
        Unow(IMT-1,JMT,K) = Unow(IMT-1,JMT-1,K)
        Unow(IMT,JMT,K) = Unow(IMT,JMT-1,K)
      ENDDO
    else if ( ID_OBC_U_E .eq. 5 ) then
      do k = 1, km
        do j = 2, jmt-1
          dvdttmp =  up(imt-3,j  ,k) - unow(imt-3,j  ,k)
          dvdnmtmp = unow (imt-3,j  ,k) - unow(imt-4,j  ,k)

          if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
            tau_obc = tau_in_c5
            dvdttmp = 0.0d0
          else
            tau_obc = tau_out_c5
          endif
          tau_obc = tau_obc*dtc_openbc

          if (dvdttmp*(up(imt-3,j+1,k) - up(imt-3,j-1,k)) .gt. 0.0d0) then
            dvdtgtmp = up(imt-3,j  ,k) - up(imt-3,j-1,k)
          else
            dvdtgtmp = up(imt-3,j+1,k) - up(imt-3,j  ,k)
          end if

          cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp, 3.0d-7)
          cfnm = dvdttmp*dvdnmtmp
          if (ID_OBC_2DRad_NPO_UV5 == 0 ) then
            cftg = dvdttmp*dvdtgtmp
          else
            cftg = 0.0d0
          end if

          unow(imt-2,j  ,k) = (cff1*up(imt-2,j  ,k)+cfnm*unow(imt-3,j  ,k)-         &
              &       max(cftg,0.0d0)*(up(imt-2,j  ,k)-up(imt-2,j-1,k))-              &
              &       min(cftg,0.0d0)*(up(imt-2,j+1,k)-up(imt-2,j  ,k)))/max(cff1+cfnm,1.0d-20)

          unow(imt-2,j  ,k) = unow(imt-2,j  ,k) + tau_obc * ((UBCLE(2,J,K))-unow(imt-2,j,k))
          unow(imt-1,j,k) = unow(imt-2,j,k) * VIV(IMT-1,J,K)
          unow(imt  ,j,k) = unow(imt-2,j,k) * VIV(IMT  ,J,K)
        enddo
        if (ID_OBC_SMO_3DUV_NM .eq. 1 ) then
          do j = 2,jmt-1
            unow(imt-2,j,k) = 0.50d0*unow(imt-2,j,k) + 0.50d0*unow(imt-2,j-1,k)
          enddo
        endif
        unow(imt-2,1  ,k) = unow(imt-2,2    ,k) * VIV(IMT-2,1,K)
        unow(imt-1,1  ,k) = unow(imt-1,2    ,k) * VIV(IMT-1,1,K)
        unow(imt  ,1  ,k) = unow(imt  ,2    ,k) * VIV(IMT  ,1,K)
        unow(imt-2,jmt,k) = unow(imt-2,jmt-1,k) * VIV(IMT-2,JMT,K)
        unow(imt-1,jmt,k) = unow(imt-1,jmt-1,k) * VIV(IMT-1,JMT,K)
        unow(imt  ,jmt,k) = unow(imt  ,jmt-1,k) * VIV(IMT  ,JMT,K)
      enddo
    ENDIF
  ELSE IF ( JEWSN .EQ. 2 ) THEN
    ! Western Boundary Condition

    IF ( ID_OBC_U_W .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For U at Western Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          U(4,J,K) = U(5,J,K) * VIV(5,J,K)
          U(3,J,K) = U(4,J,K) * VIV(3,J,K)
          U(2,J,K) = U(4,J,K) * VIV(2,J,K)
          U(1,J,K) = U(4,J,K) * VIV(1,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_W .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For U at Western Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          U(4,J,K) = UBCLW(1,J,K) * VIV(4,J,K)
          U(3,J,K) = U(4,J,K) * VIV(3,J,K)
          U(2,J,K) = U(4,J,K) * VIV(2,J,K)
          U(1,J,K) = U(4,J,K) * VIV(1,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_W .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For U at Western Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          UDTDX = U(4,J,K) * DTC * OTX(J)

          IF ( UDTDX.GT.0.0D0 ) THEN
            U(4,J,K) = UP(4,J,K) - UDTDX * ( UP(4,J,K) - UBCLW(1,J,K) )
          ELSE
            U(4,J,K) = UP(4,J,K) - UDTDX * ( UP(5,J,K) - UP(4,J,K) )
          ENDIF
          U(3,J,K) = U(4,J,K) * VIV(3,J,K)
          U(2,J,K) = U(4,J,K) * VIV(2,J,K)
          U(1,J,K) = U(4,J,K) * VIV(1,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_W .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For U at Western Edge
      ! (Marchesiello P., 2001, Ocean Modelling 3 (2001) 1-20 )

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 2, JMT-1
          DVDTTMP = UP(5,J,K) - Unow(5,J,K)
          DVDNMTMP = Unow(5,J,K) - Unow(6,J,K)
          DVDTG2TMP = UP(5,J+1,K) - UP(5,J-1,K)

          IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
            DVDTTMP = 0.0D0
            Tau_obc = TAUIN
          ELSE
            Tau_obc = TAUOUT
          ENDIF

          IF ( ( DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            DVDTGTMP = UP(5,J,K) - UP(5,J-1,K)
          ELSE
            DVDTGTMP = UP(5,J+1,K) - UP(5,J,K)
          ENDIF

          IF ( ID_OBC_2DRad_U .EQ. 0 ) THEN
            DVDTGTMP = 0.0D0
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP / (  DYT(J)*DYT(J) ) * &
                (  0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP * DYT(J) * DYT(J) / &
                ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
          ELSE
            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP / ( DYT(J+1)*DYT(J+1) ) * &
                ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP * DYT(J+1) * DYT(J+1) / &
                ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
          ENDIF

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = MAX(CFF2,EPSCF)
          CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
          CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

          IF ( ID_OBC_2DRad_NPO_U .EQ. 1 ) THEN
            CFTG = 0.0D0
          ENDIF

          CMIU = CFF1 + CFNM
          Unow(4,J,K) = ( CFF1*UP(4,J,K) + &
              CFNM*Unow(5,J,K) -  &
              MAX(CFTG,0.0D0)*(CFF1/CFF2)*(UP(4,J,K)-UP(4,J-1,K)) - &
              MIN(CFTG,0.0D0)*(CFF1/CFF2)*(UP(4,J+1,K)-UP(4,J,K)) ) / CMIU

          Unow(4,J,K) = Unow(4,J,K) + &
              DTC_OPENBC/Tau_obc * ( ( 0.2*UBCLW(1,J-1,K)+0.2*UBCLW(1,J+1,K)+ &
              0.4*UBCLW(1,J,K)+0.2*UBCLE(2,J,K) )-UP(4,J,K) ) * CFF1/CMIU

          Unow(3,J,K) = Unow(4,J,K) * VIV(3,J,K)
          Unow(2,J,K) = Unow(4,J,K) * VIV(2,J,K)
          Unow(1,J,K) = Unow(4,J,K) * VIV(1,J,K)
        ENDDO
        do j = 2,jmt-1
          unow(4,j,k) = 0.25d0*unow(4,j-1,k) + 0.50d0*unow(4,j,k) + 0.25d0*unow(4,j+1,k)
        enddo
        Unow(4,1,K) = Unow(4,2,K) * VIV(4,1,K)
        Unow(3,1,K) = Unow(3,2,K) * VIV(3,1,K)
        Unow(2,1,K) = Unow(2,2,K) * VIV(2,1,K)
        Unow(1,1,K) = Unow(1,2,K) * VIV(1,1,K)
        Unow(4,JMT,K) = Unow(4,JMT-1,K) * VIV(4,JMT,K)
        Unow(3,JMT,K) = Unow(3,JMT-1,K) * VIV(3,JMT,K)
        Unow(2,JMT,K) = Unow(2,JMT-1,K) * VIV(2,JMT,K)
        Unow(1,JMT,K) = Unow(1,JMT-1,K) * VIV(1,JMT,K)
      ENDDO
    else if ( ID_OBC_U_W .eq. 5 ) then
      do k = 1, km
        do j = 2, jmt-1
          dvdttmp  = up(5,j,k) - unow(5,j,k)
          dvdnmtmp = unow (5,j,k) - unow(6,j,k)

          if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
            tau_obc = tau_in_c5
            dvdttmp = 0.0d0
          else
            tau_obc = tau_out_c5
          endif
          tau_obc = tau_obc*dtc_openbc

          if (dvdttmp*(up(5,j+1,k)-up(5,j-1,k)) .gt. 0.0d0) then
            dvdtgtmp = up(5,j  ,k) - up(5,j-1,k)
          else
            dvdtgtmp = up(5,j+1,k) - up(5,j  ,k)
          endif

          cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
          cfnm = dvdttmp*dvdnmtmp  ! here the phase v is actually -Cnm
          if (ID_OBC_2DRad_NPO_UV5 == 0 ) then
            cftg = dvdttmp*dvdtgtmp
          else
            cftg = 0.0d0
          end if

          unow(4,j,k) = (cff1*up(4,j,k)+cfnm*unow(5,j,k)-                   &
              &       max(cftg,0.0d0)*(up(4,j  ,k)-up(4,j-1,k))-              &
              &       min(cftg,0.0d0)*(up(4,j+1,k)-up(4,j  ,k)))/max(cff1+cfnm,1.0d-20)

          unow(4,j,k) = unow(4,j,k) + tau_obc * ((UBCLW(1,J,K))-unow(4,j,k))
          Unow(3,J,K) = Unow(4,J,K) * VIV(3,J,K)
          Unow(2,J,K) = Unow(4,J,K) * VIV(2,J,K)
          Unow(1,J,K) = Unow(4,J,K) * VIV(1,J,K)
        enddo
        if (ID_OBC_SMO_3DUV_NM .eq. 1 ) then
          do j = 2,jmt-1
            unow(4,j,k) = 0.50d0*unow(4,j,k) + 0.50d0*unow(4,j-1,k)
          enddo
        endif
        Unow(4,1,K) = Unow(4,2,K) * VIV(4,1,K)
        Unow(3,1,K) = Unow(3,2,K) * VIV(3,1,K)
        Unow(2,1,K) = Unow(2,2,K) * VIV(2,1,K)
        Unow(1,1,K) = Unow(1,2,K) * VIV(1,1,K)
        Unow(4,JMT,K) = Unow(4,JMT-1,K) * VIV(4,JMT,K)
        Unow(3,JMT,K) = Unow(3,JMT-1,K) * VIV(3,JMT,K)
        Unow(2,JMT,K) = Unow(2,JMT-1,K) * VIV(2,JMT,K)
        Unow(1,JMT,K) = Unow(1,JMT-1,K) * VIV(1,JMT,K)
      enddo
    ENDIF
  ELSE IF ( JEWSN .EQ. 3 ) THEN
    ! Southern Boundary Condition
    IF ( ID_OBC_U_S .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For U at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          U(I,JET-3,K) = U(I,JET-4,K) * VIV(I,JET-4,K)
          U(I,JET-2,K) = U(I,JET-3,K) * VIV(I,JET-2,K)
          U(I,JET-1,K) = U(I,JET-3,K) * VIV(I,JET-1,K)
          U(I,JET,K) = U(I,JET-3,K) * VIV(I,JET,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_S .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For U at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          U(I,JET-3,K) = UBCLS(I,2,K) * VIV(I,JET-3,K)
          U(I,JET-2,K) = U(I,JET-3,K) * VIV(I,JET-2,K)
          U(I,JET-1,K) = U(I,JET-3,K) * VIV(I,JET-1,K)
          U(I,JET,K) = U(I,JET-3,K) * VIV(I,JET,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_S .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For U at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          UDTDX = Vp(I,JET-3,K) * DTC_OPENBC / DYT(JET-3)

          IF ( UDTDX.LE.0.0D0 ) THEN
            Unow(I,JET-3,K) = Up(I,JET-3,K) - UDTDX * ( UBCLS(I,2,K) - Up(I,JET-3,K) )
          ELSE
            Unow(I,JET-3,K) = Up(I,JET-3,K) - UDTDX * ( Up(I,JET-3,K) - Up(I,JET-4,K) )
          ENDIF
          Unow(I,JET-2,K) = Unow(I,JET-3,K) * VIV(I,JET-2,K)
          Unow(I,JET-1,K) = Unow(I,JET-3,K) * VIV(I,JET-1,K)
          Unow(I,JET,K)   = Unow(I,JET-3,K) * VIV(I,JET,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_S .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For 3D U at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 2, IMT-1
          DVDTTMP = UP(I,JET-4,K) - Unow(I,JET-4,K)
          DVDNMTMP = Unow(I,JET-4,K) - Unow(I,JET-5,K)
          DVDTG2TMP = UP(I+1,JET-4,K) - UP(I-1,JET-4,K)

          IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
            DVDTTMP = 0.0D0
            Tau_obc = TAUIN
          ELSE
            Tau_obc = TAUOUT
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            DVDTGTMP = UP(I,JET-4,K) - UP(I-1,JET-4,K)
          ELSE
            DVDTGTMP = UP(I+1,JET-4,K) - UP(I,JET-4,K)
          ENDIF

          IF ( ID_OBC_2DRad_U .EQ. 0 ) THEN
            DVDTGTMP = 0.0D0
          ENDIF

          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP * DYT(JET-4) * DYT(JET-4) / &
              ( ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) * ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP / ( DYT(JET-4)*DYT(JET-4) ) * &
              ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) * ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) )

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = MAX(CFF2,EPSCF)
          CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
          CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

          IF ( ID_OBC_2DRad_NPO_U .EQ. 1 ) THEN
            CFTG = 0.0D0
          ENDIF

          CMIU = CFF1 + CFNM
          Unow(I,JET-3,K) = ( CFF1*UP(I,JET-3,K) + &
              CFNM*Unow(I,JET-4,K) -  &
              MAX(CFTG,0.0D0)*(CFF1/CFF2)*(UP(I,JET-3,K)-UP(I-1,JET-3,K)) - &
              MIN(CFTG,0.0D0)*(CFF1/CFF2)*(UP(I+1,JET-3,K)-UP(I,JET-3,K)) ) / CMIU

          Unow(I,JET-3,K) = Unow(I,JET-3,K) + &
              DTC_OPENBC/Tau_obc * ( ( 0.2*UBCLS(I,1,K)+0.4*UBCLS(I,2,K)+0.2*UBCLS(I-1,2,K)+ &
              0.2*UBCLS(I+1,2,K) )-UP(I,JET-3,K) ) * CFF1/CMIU

          Unow(I,JET-2,K) = Unow(I,JET-3,K) * VIV(I,JET-2,K)
          Unow(I,JET-1,K) = Unow(I,JET-3,K) * VIV(I,JET-1,K)
          Unow(I,JET,K)   = Unow(I,JET-3,K) * VIV(I,JET,K)
        ENDDO
        do i = 2,imt-1
          unow(i,JET-3,k) = 0.25d0*unow(i-1,JET-3,k) + 0.50d0*unow(i,JET-3,k) + 0.25d0*unow(i+1,JET-3,k)
        enddo
        Unow(1,JET-3,K) = Unow(2,JET-3,K) * VIV(1,JET-3,K)
        Unow(1,JET-2,K) = Unow(2,JET-2,K) * VIV(1,JET-2,K)
        Unow(1,JET-1,K) = Unow(2,JET-1,K) * VIV(1,JET-1,K)
        Unow(1,JET,K)   = Unow(2,JET,K) * VIV(1,JET,K)
        Unow(IMT,JET-3,K) = Unow(IMT-1,JET-3,K) * VIV(IMT,JET-3,K)
        Unow(IMT,JET-2,K) = Unow(IMT-1,JET-2,K) * VIV(IMT,JET-2,K)
        Unow(IMT,JET-1,K) = Unow(IMT-1,JET-1,K) * VIV(IMT,JET-1,K)
        Unow(IMT,JET,K)   = Unow(IMT-1,JET,K) * VIV(IMT,JET,K)
      ENDDO
    else if ( ID_OBC_U_S .eq. 5 ) then
      do k = 1, km
        do i = 2,imt-1
          dvdttmp  = up(i  ,jet-4,k)-unow(i  ,jet-4,k)
          dvdnmtmp = unow (i  ,jet-4,k)-unow(i  ,jet-5,k)

          if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
            tau_obc = tau_in_c5
            dvdttmp = 0.0d0
          else
            tau_obc = tau_out_c5
          endif
          tau_obc  = tau_obc*dtc_openbc

          if (dvdttmp*(up(i+1,jet-4,k)-up(i-1,jet-4,k)) .gt. 0.0d0) then
            dvdtgtmp = up(i  ,jet-4,k)-up(i-1,jet-4,k)
          else
            dvdtgtmp = up(i+1,jet-4,k)-up(i  ,jet-4,k)
          endif

          cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
          cfnm = dvdttmp*dvdnmtmp
          if (ID_OBC_2DRad_NPO_UV5 == 0 ) then
            cftg = dvdttmp*dvdtgtmp
          else
            cftg = 0.0d0
          end if

          unow(i  ,jet-3,k) = (cff1*up(i  ,jet-3,k)+cfnm*unow(i  ,jet-4,k)-         &
              &       max(cftg,0.0d0)*(up(i  ,jet-4,k)-up(i-1,jet-4,k))-              &
              &       min(cftg,0.0d0)*(up(i+1,jet-4,k)-up(i  ,jet-4,k)))/max(cff1+cfnm,1.0d-20)

          unow(i  ,jet-3,k) = unow(i  ,jet-3,k) + tau_obc * (( UBCLS(I,2,K) )-unow(i,jet-3,k))
          Unow(I,JET-2,K) = Unow(I,JET-3,K) * VIV(I,JET-2,K)
          Unow(I,JET-1,K) = Unow(I,JET-3,K) * VIV(I,JET-1,K)
          Unow(I,JET,K)   = Unow(I,JET-3,K) * VIV(I,JET,K)
        enddo
        if (ID_OBC_SMO_3DUV_TG .eq. 1 ) then
          do i = 2,imt-1
            unow(i,JET-3,k) = 0.50d0*unow(i,JET-3,k) + 0.50d0*unow(i-1,JET-3,k)
          enddo
        endif
        Unow(1,JET-3,K) = Unow(2,JET-3,K) * VIV(1,JET-3,K)
        Unow(1,JET-2,K) = Unow(2,JET-2,K) * VIV(1,JET-2,K)
        Unow(1,JET-1,K) = Unow(2,JET-1,K) * VIV(1,JET-1,K)
        Unow(1,JET,K)   = Unow(2,JET,K) * VIV(1,JET,K)
        Unow(IMT,JET-3,K) = Unow(IMT-1,JET-3,K) * VIV(IMT,JET-3,K)
        Unow(IMT,JET-2,K) = Unow(IMT-1,JET-2,K) * VIV(IMT,JET-2,K)
        Unow(IMT,JET-1,K) = Unow(IMT-1,JET-1,K) * VIV(IMT,JET-1,K)
        Unow(IMT,JET,K)   = Unow(IMT-1,JET,K) * VIV(IMT,JET,K)
      enddo
    ENDIF
  ELSE IF ( JEWSN .EQ. 4 ) THEN
    ! Northern Boundary Condition

    IF ( ID_OBC_U_N .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For U at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          U(I,3,K) = U(I,4,K) * VIV(I,4,K)
          U(I,2,K) = U(I,3,K) * VIV(I,2,K)
          U(I,1,K) = U(I,3,K) * VIV(I,2,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_N .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For U at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          U(I,3,K) = UBCLN(I,1,K) * VIV(I,3,K)
          U(I,2,K) = U(I,3,K) * VIV(I,2,K)
          U(I,1,K) = U(I,3,K) * VIV(I,2,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_N .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For U at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          UDTDX = Vp(I,3,K) * DTC_OPENBC / DYT(3)

          IF ( UDTDX.GT.0.0D0 ) THEN
            Unow(I,3,K) = Up(I,3,K) - UDTDX * ( Up(I,3,K) - UBCLN(I,1,K) )
          ELSE
            Unow(I,3,K) = Up(I,3,K) - UDTDX * ( Up(I,4,K) - Up(I,3,K) )
          ENDIF

          Unow(I,2,K) = Unow(I,3,K) * VIV(I,2,K)
          Unow(I,1,K) = Unow(I,3,K) * VIV(I,1,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_U_N .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For U at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 2, IMT-1
          DVDTTMP = UP(I,4,K) - Unow(I,4,K)
          DVDNMTMP = Unow(I,4,K) - Unow(I,5,K)
          DVDTG2TMP = UP(I+1,4,K) - UP(I-1,4,K)

          IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
            DVDTTMP = 0.0D0
            Tau_obc = TAUIN
          ELSE
            Tau_obc = TAUOUT
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            DVDTGTMP = UP(I,4,K) - UP(I-1,4,K)
          ELSE
            DVDTGTMP = UP(I+1,4,K) - UP(I,4,K)
          ENDIF

          IF ( ID_OBC_2DRad_U .EQ. 0 ) THEN
            DVDTGTMP = 0.0D0
          ENDIF

          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP * DYT(5)*DYT(5) / &
              ( ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) * ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP / ( DYT(5)*DYT(5) ) * &
              ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) * ( 0.50D0/OTX(4)+0.50D0/OTX(5) )

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = MAX(CFF2,EPSCF)
          CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
          CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

          IF ( ID_OBC_2DRad_NPO_U .EQ. 1 ) THEN
            CFTG = 0.0D0
          ENDIF

          CMIU = CFF1 + CFNM
          Unow(I,3,K) = ( CFF1*UP(I,3,K) + &
              CFNM*Unow(I,4,K) -  &
              MAX(CFTG,0.0D0)*(CFF1/CFF2)*(UP(I,3,K)-UP(I-1,3,K)) - &
              MIN(CFTG,0.0D0)*(CFF1/CFF2)*(UP(I+1,3,K)-UP(I,3,K)) ) / CMIU

          Unow(I,3,K) = Unow(I,3,K) + &
              DTC_OPENBC/Tau_obc * ( ( 0.2*UBCLN(I-1,1,K)+0.2*UBCLN(I+1,1,K)+ &
              0.4*UBCLN(I,1,K)+0.2*UBCLN(I,2,K) )-UP(I,3,K) ) * CFF1/CMIU

          Unow(I,2,K) = Unow(I,3,K) * VIV(I,2,K)
          Unow(I,1,K) = Unow(I,3,K) * VIV(I,1,K)
        ENDDO
        do i = 2,imt-1
          unow(i,3,k) = 0.25d0*unow(i-1,3,k) + 0.50d0*unow(i,3,k) + 0.25d0*unow(i+1,3,k)
        enddo
        Unow(1,3,K) = Unow(2,3,k) * VIV(1,3,K)
        Unow(1,2,K) = Unow(2,2,K) * VIV(1,2,K)
        Unow(1,1,K) = Unow(2,1,K) * VIV(1,1,K)
        Unow(IMT,3,K) = Unow(IMT-1,3,K) * VIV(IMT,3,K)
        Unow(IMT,2,K) = Unow(IMT-1,2,K) * VIV(IMT,2,K)
        Unow(IMT,1,K) = Unow(IMT-1,1,K) * VIV(IMT,1,K)
      ENDDO
    else if ( ID_OBC_U_N .eq. 5 ) then
      do k = 1, km
        do i = 2, imt-1
          dvdttmp  = up(i,4,k)-unow(i,4,k)
          dvdnmtmp = unow (i,4,k)-unow(i,5,k)

          if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
            tau_obc = tau_in_c5
            dvdttmp = 0.0d0
          else
            tau_obc = tau_out_c5
          endif
          tau_obc = tau_obc*dtc_openbc

          if (dvdttmp*(up(i+1,4,k)-up(i-1,4,k)) .gt. 0.0d0) then
            dvdtgtmp = up(i  ,4,k) - up(i-1,4,k)
          else
            dvdtgtmp = up(i+1,4,k) - up(i  ,4,k)
          endif

          cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
          cfnm = dvdttmp*dvdnmtmp
          if (ID_OBC_2DRad_NPO_UV5 == 0 ) then
            cftg = dvdttmp*dvdtgtmp
          else
            cftg = 0.0d0
          end if

          unow(i,3,k) = (cff1*up(i,3,k)+cfnm*unow(i,4,k)-                   &
              &       max(cftg,0.0d0)*(up(i  ,3,k)-up(i-1,3,k))-              &
              &       min(cftg,0.0d0)*(up(i+1,3,k)-up(i  ,3,k)))/max(cff1+cfnm,1.0d-20)

          unow(i,3,k) = unow(i,3,k) + tau_obc * (( UBCLN(I,1,K) )-unow(i,1,k))
          Unow(I,2,K) = Unow(I,3,K) * VIV(I,2,K)
          Unow(I,1,K) = Unow(I,3,K) * VIV(I,1,K)
        enddo
        if (ID_OBC_SMO_3DUV_TG .eq. 1 ) then
          do i = 2,imt-1
            unow(i,3,k) = 0.50d0*unow(i,3,k) + 0.50d0*unow(i-1,3,k)
          enddo
        endif
        Unow(1,3,K) = Unow(2,3,k) * VIV(1,3,K)
        Unow(1,2,K) = Unow(2,2,K) * VIV(1,2,K)
        Unow(1,1,K) = Unow(2,1,K) * VIV(1,1,K)
        Unow(IMT,3,K) = Unow(IMT-1,3,K) * VIV(IMT,3,K)
        Unow(IMT,2,K) = Unow(IMT-1,2,K) * VIV(IMT,2,K)
        Unow(IMT,1,K) = Unow(IMT-1,1,K) * VIV(IMT,1,K)
      enddo
    ENDIF
  ENDIF

  !$OMP PARALLEL DO PRIVATE (K,J,I)
  DO K=1, KM
    DO J=1, JMT
      DO I=1, IMT
        Unow(I,J,K) = Unow(I,J,K) * VIV(I,J,K)
      ENDDO
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE OBC_3DU

SUBROUTINE OBC_3DV(JEWSN)

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, JET, KM, I, J, K
  use   pconst_mod, only: VIV, DTC, DYT, OTX, DTC2
  use   dyn_mod,    only: V, U, VP, UP

#if ( defined SPMD )
  use msg_mod
#endif

#ifdef OBCDT
  use   openbc_mod
#endif

  IMPLICIT NONE
  INTEGER :: JEWSN
  REAL(r8) :: Tau_obc
  REAL(r8) :: EPSCF
  REAL(r8) :: CFF1, CFF2
  REAL(r8) :: UDTDX
  REAL(r8) :: DVDTTMP, DVDNMTMP, DVDTGTMP, DVDTG2TMP, CFTG, CFNM
  REAL(r8) :: CMIU
  EPSCF = 1.0D-20

  IF ( JEWSN .EQ. 1 ) THEN
    ! Eastern Boundary Condition
    IF ( ID_OBC_V_E .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For V at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          V(IMT-2,J,K) = V(IMT-3,J,K) * VIV(IMT-3,J,K)
          V(IMT-1,J,K) = V(IMT-2,J,K) * VIV(IMT-1,J,K)
          V(IMT,J,K) = V(IMT-2,J,K)  * VIV(IMT,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_E .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For V at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          V(IMT-2,J,K) = VBCLE(2,J,K) * VIV(IMT-2,J,K)
          V(IMT-1,J,K) = V(IMT-2,J,K) * VIV(IMT-1,J,K)
          V(IMT,J,K) = V(IMT-2,J,K) * VIV(IMT,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_E .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For V at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          UDTDX = Up(IMT-2,J,K) * DTC_OPENBC * OTX(J)

          IF ( UDTDX.LE.0.0D0 ) THEN
            Vnow(IMT-2,J,K) = Vp(IMT-2,J,K) - UDTDX * ( VBCLE(2,J,K) - Vp(IMT-2,J,K) )
          ELSE
            Vnow(IMT-2,J,K) = Vp(IMT-2,J,K) - UDTDX * ( Vp(IMT-2,J,K) - Vp(IMT-3,J,K) )
          ENDIF
          Vnow(IMT-1,J,K) = Vnow(IMT-2,J,K) * VIV(IMT-1,J,K)
          Vnow(IMT,J,K) = Vnow(IMT-2,J,K) * VIV(IMT,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_E .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For V at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 2, JMT-1
          DVDTTMP = VP(IMT-3,J,K) - Vnow(IMT-3,J,K)
          DVDNMTMP = Vnow(IMT-3,J,K) - Vnow(IMT-4,J,K)
          DVDTG2TMP = VP(IMT-3,J+1,K) - VP(IMT-3,J-1,K)

          IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
            DVDTTMP = 0.0D0
            Tau_obc = TAUIN
          ELSE
            Tau_obc = TAUOUT
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            DVDTGTMP = VP(IMT-3,J,K) - VP(IMT-3,J-1,K)
          ELSE
            DVDTGTMP = VP(IMT-3,J+1,K) - VP(IMT-3,J,K)
          ENDIF

          IF ( ID_OBC_2DRad_V .EQ. 0 ) THEN
            DVDTGTMP = 0.0D0
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP / ( DYT(J)*DYT(J) ) * &
                ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP * DYT(J) * DYT(J) / &
                ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
          ELSE
            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP / ( DYT(J+1)*DYT(J+1) ) * &
                ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP*DVDNMTMP * DYT(J+1)*DYT(J+1) / &
                ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) *( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
          ENDIF

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = MAX(CFF2,EPSCF)
          CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
          CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

          IF ( ID_OBC_2DRad_NPO_V .EQ. 1 ) THEN
            CFTG = 0.0D0
          ENDIF

          CMIU = CFF1 + CFNM
          Vnow(IMT-2,J,K) = ( CFF1*VP(IMT-2,J,K) + &
              CFNM*Vnow(IMT-3,J,K) -  &
              MAX(CFTG,0.0D0)*(CFF1/CFF2)*(VP(IMT-2,J,K)-VP(IMT-2,J-1,K)) - &
              MIN(CFTG,0.0D0)*(CFF1/CFF2)*(VP(IMT-2,J+1,K)-VP(IMT-2,J,K)) )/CMIU

          Vnow(IMT-2,J,K) = Vnow(IMT-2,J,K) + &
              DTC_OPENBC/Tau_obc * ( ( 0.2*VBCLE(1,J,K)+0.4*VBCLE(2,J,K)+ &
              0.2*VBCLE(2,J-1,K)+0.2*VBCLE(2,J+1,K) )-VP(IMT-2,J,K) ) * CFF1/CMIU
          Vnow(IMT-1,J,K) = Vnow(IMT-2,J,K) * VIV(IMT-1,J,K)
          Vnow(IMT,J,K) = Vnow(IMT-2,J,K) * VIV(IMT,J,K)
        ENDDO
        do j = 2,jmt-1
          vnow(imt-2,j,k) = 0.25d0*vnow(imt-2,j-1,k) + 0.50d0*vnow(imt-2,j,k) + 0.25d0*vnow(imt-2,j+1,k)
        enddo
        Vnow(IMT-2,1,K) = Vnow(IMT-2,2,K) * VIV(IMT-2,1,K)
        Vnow(IMT-1,1,K) = Vnow(IMT-1,2,K) * VIV(IMT-1,1,K)
        Vnow(IMT,1,K) = Vnow(IMT,2,K) * VIV(IMT,1,K)
        Vnow(IMT-2,JMT,K) = Vnow(IMT-2,JMT-1,K) * VIV(IMT-2,JMT,K)
        Vnow(IMT-1,JMT,K) = Vnow(IMT-1,JMT-1,K) * VIV(IMT-1,JMT,K)
        Vnow(IMT,JMT,K) = Vnow(IMT,JMT-1,K) * VIV(IMT,JMT,K)
      ENDDO
    else if ( ID_OBC_V_E .eq. 5 ) then
      do k = 1, km
        do j = 2, jmt-1
          dvdttmp =  vp(imt-3,j  ,k) - vnow(imt-3,j  ,k)
          dvdnmtmp = vnow (imt-3,j  ,k) - vnow(imt-4,j  ,k)
          if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
            tau_obc = tau_in_c5
            dvdttmp = 0.0d0
          else
            tau_obc = tau_out_c5
          endif
          tau_obc = tau_obc*dtc_openbc

          if (dvdttmp*(vp(imt-3,j+1,k) - vp(imt-3,j-1,k)) .gt. 0.0d0) then
            dvdtgtmp = vp(imt-3,j  ,k) - vp(imt-3,j-1,k)
          else
            dvdtgtmp = vp(imt-3,j+1,k) - vp(imt-3,j  ,k)
          endif

          cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp, 3.0d-7)
          cfnm = dvdttmp*dvdnmtmp
          if (ID_OBC_2DRad_NPO_UV5 == 0 ) then
            cftg = dvdttmp*dvdtgtmp
          else
            cftg = 0.0d0
          end if

          vnow(imt-2,j  ,k) = (cff1*vp(imt-2,j  ,k)+cfnm*vnow(imt-3,j  ,k)-         &
              &       max(cftg,0.0d0)*(vp(imt-2,j  ,k)-vp(imt-2,j-1,k))-              &
              &       min(cftg,0.0d0)*(vp(imt-2,j+1,k)-vp(imt-2,j  ,k)))/max(cff1+cfnm,1.0d-20)

          vnow(imt-2,j  ,k) = vnow(imt-2,j  ,k) + tau_obc * (( VBCLE(2,J,K) )-vnow(imt-2,j,k))
          vnow(imt-1,j,k) = vnow(imt-2,j,k)
          vnow(imt  ,j,k) = vnow(imt-2,j,k)
        enddo
        if (ID_OBC_SMO_3DUV_TG .eq. 1 ) then
          do j = 2,jmt-1
            vnow(imt-2,j,k) = 0.50d0*vnow(imt-2,j,k) + 0.50d0*vnow(imt-2,j-1,k)
          enddo
        endif
        vnow(imt-2,1  ,k) = vnow(imt-2,2    ,k)
        vnow(imt-1,1  ,k) = vnow(imt-1,2    ,k)
        vnow(imt  ,1  ,k) = vnow(imt  ,2    ,k)
        vnow(imt-2,jmt,k) = vnow(imt-2,jmt-1,k)
        vnow(imt-1,jmt,k) = vnow(imt-1,jmt-1,k)
        vnow(imt  ,jmt,k) = vnow(imt  ,jmt-1,k)
      enddo
    ENDIF
  ELSE IF ( JEWSN .EQ. 2 ) THEN
    ! Western Boundary Condition
    IF ( ID_OBC_V_W .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For V at Western Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          V(4,J,K) = V(5,J,K) * VIV(5,J,K)
          V(3,J,K) = V(4,J,K) * VIV(3,J,K)
          V(2,J,K) = V(4,J,K) * VIV(2,J,K)
          V(1,J,K) = V(4,J,K) * VIV(1,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_W .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For V at Western Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          V(4,J,K) = VBCLW(1,J,K) * VIV(4,J,K)
          V(3,J,K) = V(4,J,K) * VIV(3,J,K)
          V(2,J,K) = V(4,J,K) * VIV(2,J,K)
          V(1,J,K) = V(4,J,K) * VIV(1,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_W .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For V at Western Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 1, JMT
          UDTDX = Up(4,J,K) * DTC_OPENBC * OTX(J)
          IF ( UDTDX.GT.0.0D0 ) THEN
            Vp(4,J,K) = VP(4,J,K) - UDTDX * ( VP(4,J,K) - VBCLW(1,J,K) )
          ELSE
            Vp(4,J,K) = VP(4,J,K) - UDTDX * ( VP(5,J,K) - VP(4,J,K) )
          ENDIF
          Vnow(3,J,K) = Vnow(4,J,K) * VIV(3,J,K)
          Vnow(2,J,K) = Vnow(4,J,K) * VIV(2,J,K)
          Vnow(1,J,K) = Vnow(4,J,K) * VIV(1,J,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_W .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For V at Western Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 2, JMT-1
          DVDTTMP = VP(5,J,K) - Vnow(5,J,K)
          DVDNMTMP = Vnow(5,J,K) - Vnow(6,J,K)
          DVDTG2TMP = VP(5,J+1,K) - VP(5,J-1,K)

          IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
            DVDTTMP = 0.0D0
            Tau_obc = TAUIN
          ELSE
            Tau_obc = TAUOUT
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            DVDTGTMP = VP(5,J,K) - VP(5,J-1,K)
          ELSE
            DVDTGTMP = VP(5,J+1,K) - VP(5,J,K)
          ENDIF

          IF ( ID_OBC_2DRad_V .EQ. 0 ) THEN
            DVDTGTMP = 0.0D0
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP / ( DYT(J)*DYT(J) ) * &
                ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP * DYT(J) * DYT(J) / &
                ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
          ELSE
            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP / ( DYT(J+1)*DYT(J+1) ) * &
                ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP * DYT(J+1) * DYT(J+1) / &
                ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
          ENDIF

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = MAX(CFF2,EPSCF)
          CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
          CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

          IF ( ID_OBC_2DRad_NPO_V .EQ. 1 ) THEN
            CFTG = 0.0D0
          ENDIF

          CMIU = CFF1 + CFNM
          Vnow(4,J,K) = ( CFF1*VP(4,J,K) + &
              CFNM*Vnow(5,J,K) -  &
              MAX(CFTG,0.0D0)*(CFF1/CFF2)*(VP(4,J,K)-VP(4,J-1,K)) - &
              MIN(CFTG,0.0D0)*(CFF1/CFF2)*(VP(4,J+1,K)-VP(4,J,K)) ) / CMIU

          Vnow(4,J,K) = Vnow(4,J,K) + &
              DTC_OPENBC/Tau_obc * ( ( 0.2*VBCLW(1,J-1,K)+0.2*VBCLW(1,J+1,K)+ &
              0.4*VBCLW(1,J,K)+0.2*VBCLW(2,J,K) )-VP(4,J,K) ) * CFF1/CMIU

          Vnow(3,J,K) = Vnow(4,J,K) * VIV(3,J,K)
          Vnow(2,J,K) = Vnow(4,J,K) * VIV(2,J,K)
          Vnow(1,J,K) = Vnow(4,J,K) * VIV(1,J,K)
        ENDDO
        do j = 2,jmt-1
          vnow(4,j,k) = 0.25d0*vnow(4,j-1,k) + 0.50d0*vnow(4,j,k) + 0.25d0*vnow(4,j+1,k)
        enddo
        Vnow(4,1,K) = Vnow(4,2,K) * VIV(4,1,K)
        Vnow(3,1,K) = Vnow(3,2,K) * VIV(3,1,K)
        Vnow(2,1,K) = Vnow(2,2,K) * VIV(2,1,K)
        Vnow(1,1,K) = Vnow(1,2,K) * VIV(1,1,K)
        Vnow(4,JMT,K) = Vnow(4,JMT-1,K) * VIV(4,JMT,K)
        Vnow(3,JMT,K) = Vnow(3,JMT-1,K) * VIV(3,JMT,K)
        Vnow(2,JMT,K) = Vnow(2,JMT-1,K) * VIV(2,JMT,K)
        Vnow(1,JMT,K) = Vnow(1,JMT-1,K) * VIV(1,JMT,K)
      ENDDO
    else if ( ID_OBC_V_W .eq. 5 ) then
      do k = 1, km
        do j = 2, jmt-1
          dvdttmp  = vp(5,j,k) - vnow(5,j,k)
          dvdnmtmp = vnow (5,j,k) - vnow(6,j,k)

          if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
            tau_obc = tau_in_c5
            dvdttmp = 0.0d0
          else
            tau_obc = tau_out_c5
          endif
          tau_obc = tau_obc*dtc_openbc

          if (dvdttmp*(vp(5,j+1,k)-vp(5,j-1,k)) .gt. 0.0d0) then
            dvdtgtmp = vp(5,j  ,k) - vp(5,j-1,k)
          else
            dvdtgtmp = vp(5,j+1,k) - vp(5,j  ,k)
          endif

          cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
          cfnm = dvdttmp*dvdnmtmp  ! here the phase v is actually -Cnm
          if (ID_OBC_2DRad_NPO_UV5 == 0 ) then
            cftg = dvdttmp*dvdtgtmp
          else
            cftg = 0.0d0
          end if

          vnow(4,j,k) = (cff1*vp(4,j,k)+cfnm*vnow(5,j,k)-                   &
              &       max(cftg,0.0d0)*(vp(4,j  ,k)-vp(4,j-1,k))-              &
              &       min(cftg,0.0d0)*(vp(4,j+1,k)-vp(4,j  ,k)))/max(cff1+cfnm,1.0d-20)

          vnow(4,j,k) = vnow(4,j,k) + tau_obc * (( VBCLW(1,J,K) )-vnow(4,j,k))
          vnow(3,J,K) = vnow(4,J,K) * VIV(3,J,K)
          vnow(2,J,K) = vnow(4,J,K) * VIV(2,J,K)
          vnow(1,J,K) = vnow(4,J,K) * VIV(1,J,K)
        enddo
        if (ID_OBC_SMO_3DUV_TG .eq. 1 ) then
          do j = 2,jmt-1
            vnow(4,j,k) = 0.50d0*vnow(4,j,k) + 0.50d0*vnow(4,j-1,k)
          enddo
        endif
        vnow(4,1,K) = vnow(4,2,K) * VIV(4,1,K)
        vnow(3,1,K) = vnow(3,2,K) * VIV(3,1,K)
        vnow(2,1,K) = vnow(2,2,K) * VIV(2,1,K)
        vnow(1,1,K) = vnow(1,2,K) * VIV(1,1,K)
        vnow(4,JMT,K) = vnow(4,JMT-1,K) * VIV(4,JMT,K)
        vnow(3,JMT,K) = vnow(3,JMT-1,K) * VIV(3,JMT,K)
        vnow(2,JMT,K) = vnow(2,JMT-1,K) * VIV(2,JMT,K)
        vnow(1,JMT,K) = vnow(1,JMT-1,K) * VIV(1,JMT,K)
      enddo
    ENDIF
  ELSE IF ( JEWSN .EQ. 3 ) THEN
    ! Southern Boundary Condition

    IF ( ID_OBC_V_S .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For V at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          V(I,JET-3,K) = V(I,JET-4,K) * VIV(I,JET-4,K)
          V(I,JET-2,K) = V(I,JET-3,K) * VIV(I,JET-2,K)
          V(I,JET-1,K) = V(I,JET-3,K) * VIV(I,JET-1,K)
          V(I,JET,K) = V(I,JET-3,K) * VIV(I,JET,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_S .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For V at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          V(I,JET-3,K) = VBCLW(I,2,K) * VIV(I,JET-3,K)
          V(I,JET-2,K) = V(I,JET-3,K) * VIV(I,JET-2,K)
          V(I,JET-1,K) = V(I,JET-3,K) * VIV(I,JET-1,K)
          V(I,JET,K) = V(I,JET-3,K) * VIV(I,JET,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_S .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For V at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          UDTDX = V(I,JET-3,K) * DTC / DYT(JET-3)

          IF ( UDTDX.LE.0.0D0 ) THEN
            V(I,JET-3,K) = VP(I,JET-3,K) - UDTDX * ( VBCLS(I,2,K) - VP(I,JET-3,K) )
          ELSE
            V(I,JET-3,K) = VP(I,JET-3,K) - UDTDX * ( VP(I,JET-3,K) - VP(I,JET-4,K) )
          ENDIF

          V(I,JET-2,K) = V(I,JET-3,K) * VIV(I,JET-2,K)
          V(I,JET-1,K) = V(I,JET-3,K) * VIV(I,JET-1,K)
          V(I,JET,K)   = V(I,JET-3,K) * VIV(I,JET,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_S .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition (Marchesiello P., 2001, Ocean Modelling 3 (2001) 1-20 )
      ! 2-D Radiation Boundary Condition Or Normal Projection of Oblique Radiation (NPO)
      ! The Implicit Upstream Radiation Condition For 3D V at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 2, IMT-1
          DVDTTMP = VP(I,JET-4,K) - Vnow(I,JET-4,K)
          DVDNMTMP = Vnow(I,JET-4,K) - Vnow(I,JET-5,K)
          DVDTG2TMP = VP(I+1,JET-4,K) - VP(I-1,JET-4,K)

          IF ( ( DVDTTMP*DVDNMTMP ) .LT. 0.0D0 ) THEN
            DVDTTMP = 0.0D0
            Tau_obc = TAUIN
          ELSE
            Tau_obc = TAUOUT
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            DVDTGTMP = VP(I,JET-4,K) - VP(I-1,JET-4,K)
          ELSE
            DVDTGTMP = VP(I+1,JET-4,K) - VP(I,JET-4,K)
          ENDIF

          IF ( ID_OBC_2DRad_V .EQ. 0 ) THEN
            DVDTGTMP = 0.0D0
          ENDIF

          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP * DYT(JET-4) * DYT(JET-4) / &
              ( ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) * ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP / ( DYT(JET-4)*DYT(JET-4) ) * &
              ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) * ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) )

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = MAX(CFF2,EPSCF)
          CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
          CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

          IF ( ID_OBC_2DRad_NPO_V .EQ. 1 ) THEN
            CFTG = 0.0D0
          ENDIF

          CMIU = CFF1 + CFNM
          Vnow(I,JET-3,K) = ( CFF1*VP(I,JET-3,K) + &
              CFNM*Vnow(I,JET-4,K) -  &
              MAX(CFTG,0.0D0)*(CFF1/CFF2)*(VP(I,JET-3,K)-VP(I-1,JET-3,K)) - &
              MIN(CFTG,0.0D0)*(CFF1/CFF2)*(VP(I+1,JET-3,K)-VP(I,JET-3,K)) ) / CMIU

          Vnow(I,JET-3,K) = Vnow(I,JET-3,K) + &
              DTC_OPENBC/Tau_obc * ( ( 0.2*VBCLS(I,1,K)+0.4*VBCLS(I,2,K)+ &
              0.2*VBCLS(I-1,2,K)+0.2*VBCLS(I+1,2,K) )-VP(I,JET-3,K) ) * CFF1/CMIU

          Vnow(I,JET-2,K) = Vnow(I,JET-3,K) * VIV(I,JET-2,K)
          Vnow(I,JET-1,K) = Vnow(I,JET-3,K) * VIV(I,JET-1,K)
          Vnow(I,JET,K)   = Vnow(I,JET-3,K) * VIV(I,JET,K)
        ENDDO
        do i = 2,imt-1
          vnow(i,JET-3,k) = 0.25d0*vnow(i-1,JET-3,k) + 0.50d0*vnow(i,JET-3,k) + 0.25d0*vnow(i+1,JET-3,k)
        enddo
        Vnow(1,JET-3,K) = Vnow(2,JET-3,K) * VIV(1,JET-3,K)
        Vnow(1,JET-2,K) = Vnow(2,JET-2,K) * VIV(1,JET-2,K)
        Vnow(1,JET-1,K) = Vnow(2,JET-1,K) * VIV(1,JET-1,K)
        Vnow(1,JET,K)   = Vnow(2,JET,K) * VIV(1,JET,K)
        Vnow(IMT,JET-3,K) = Vnow(IMT-1,JET-3,K) * VIV(IMT,JET-3,K)
        Vnow(IMT,JET-2,K) = Vnow(IMT-1,JET-2,K) * VIV(IMT,JET-2,K)
        Vnow(IMT,JET-1,K) = Vnow(IMT-1,JET-1,K) * VIV(IMT,JET-1,K)
        Vnow(IMT,JET,K)   = Vnow(IMT-1,JET,K) * VIV(IMT,JMT,K)
      ENDDO
    else if ( ID_OBC_V_S .eq. 5 ) then
      do k = 1, km
        do i = 2,imt-1
          dvdttmp  = vp(i  ,jet-4,k)-vnow(i  ,jet-4,k)
          dvdnmtmp = vnow (i  ,jet-4,k)-vnow(i  ,jet-5,k)

          if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
            tau_obc = tau_in_c5
            dvdttmp = 0.0d0
          else
            tau_obc = tau_out_c5
          endif
          tau_obc  = tau_obc*dtc_openbc

          if (dvdttmp*(vp(i+1,jet-4,k)-vp(i-1,jet-4,k)) .gt. 0.0d0) then
            dvdtgtmp = vp(i  ,jet-4,k)-vp(i-1,jet-4,k)
          else
            dvdtgtmp = vp(i+1,jet-4,k)-vp(i  ,jet-4,k)
          endif

          cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
          cfnm = dvdttmp*dvdnmtmp
          if (ID_OBC_2DRad_NPO_UV5 == 0 ) then
            cftg = dvdttmp*dvdtgtmp
          else
            cftg = 0.0d0
          end if

          vnow(i  ,jet-3,k) = (cff1*vp(i  ,jet-3,k)+cfnm*vnow(i  ,jet-4,k)-         &
              &       max(cftg,0.0d0)*(vp(i  ,jet-4,k)-vp(i-1,jet-4,k))-              &
              &       min(cftg,0.0d0)*(vp(i+1,jet-4,k)-vp(i  ,jet-4,k)))/max(cff1+cfnm,1.0d-20)

          vnow(i  ,jet-3,k) = vnow(i  ,jet-3,k) + tau_obc * (( VBCLS(I,2,K) )-vnow(i,jet-3,k))
          vnow(I,JET-2,K) = vnow(I,JET-3,K) * VIV(I,JET-2,K)
          vnow(I,JET-1,K) = vnow(I,JET-3,K) * VIV(I,JET-1,K)
          vnow(I,JET,K)   = vnow(I,JET-3,K) * VIV(I,JET,K)
        enddo
        if (ID_OBC_SMO_3DUV_NM .eq. 1 ) then
          do i = 2,imt-1
            vnow(i,JET-3,k) = 0.50d0*vnow(i,JET-3,k) + 0.50d0*vnow(i-1,JET-3,k)
          enddo
        endif
        vnow(1,JET-3,K) = vnow(2,JET-3,K) * VIV(1,JET-3,K)
        vnow(1,JET-2,K) = vnow(2,JET-2,K) * VIV(1,JET-2,K)
        vnow(1,JET-1,K) = vnow(2,JET-1,K) * VIV(1,JET-1,K)
        vnow(1,JET,K)   = vnow(2,JET,K) * VIV(1,JET,K)
        vnow(IMT,JET-3,K) = vnow(IMT-1,JET-3,K) * VIV(IMT,JET-3,K)
        vnow(IMT,JET-2,K) = vnow(IMT-1,JET-2,K) * VIV(IMT,JET-2,K)
        vnow(IMT,JET-1,K) = vnow(IMT-1,JET-1,K) * VIV(IMT,JET-1,K)
        vnow(IMT,JET,K)   = vnow(IMT-1,JET,K) * VIV(IMT,JET,K)
      enddo
    ENDIF
  ELSE IF ( JEWSN .EQ. 4 ) THEN
    ! Northern Boundary Condition

    IF ( ID_OBC_V_N .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For V at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          V(I,3,K) = V(I,4,K) * VIV(I,4,K)
          V(I,2,K) = V(I,3,K) * VIV(I,2,K)
          V(I,1,K) = V(I,3,K) * VIV(I,1,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_N .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For V at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          V(I,3,K) = VBCLN(I,1,K) * VIV(I,3,K)
          V(I,2,K) = V(I,3,K) * VIV(I,2,K)
          V(I,1,K) = V(I,3,K) * VIV(I,1,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_N .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For V at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 1, IMT
          UDTDX = V(I,3,K) * DTC / DYT(3)

          IF ( UDTDX.GT.0.0D0 ) THEN
            V(I,3,K) = VP(I,3,K) - UDTDX * ( VP(I,3,K) - VBCLN(I,1,K) )
          ELSE
            V(I,3,K) = VP(I,3,K) - UDTDX * ( VP(I,4,K) - VP(I,3,K) )
          ENDIF

          V(I,2,K) = V(I,3,K) * VIV(I,2,K)
          V(I,1,K) = V(I,3,K) * VIV(I,1,K)
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_V_N .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For V at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 2, IMT-1
          DVDTTMP = VP(I,4,K) - Vnow(I,4,K)
          DVDNMTMP = Vnow(I,4,K) - Vnow(I,5,K)
          DVDTG2TMP = VP(I+1,4,K) - VP(I-1,4,K)

          IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
            DVDTTMP = 0.0D0
            Tau_obc = TAUIN
          ELSE
            Tau_obc = TAUOUT
          ENDIF

          IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
            DVDTGTMP = VP(I,4,K) - VP(I-1,4,K)
          ELSE
            DVDTGTMP = VP(I+1,4,K) - VP(I,4,K)
          ENDIF

          IF ( ID_OBC_2DRad_V .EQ. 0 ) THEN
            DVDTGTMP = 0.0D0
          ENDIF

          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP * DYT(5) * DYT(5) / &
              ( ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) * ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP / ( DYT(5)*DYT(5) ) * &
              ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) * ( 0.50D0/OTX(4)+0.50D0/OTX(5) )

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = MAX(CFF2,EPSCF)
          CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
          CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

          IF ( ID_OBC_2DRad_NPO_V .EQ. 1 ) THEN
            CFTG = 0.0D0
          ENDIF

          CMIU = CFF1 + CFNM
          Vnow(I,3,K) = ( CFF1*VP(I,3,K) + &
              CFNM*Vnow(I,4,K) -  &
              MAX(CFTG,0.0D0)*(CFF1/CFF2)*(VP(I,3,K)-VP(I-1,3,K)) - &
              MIN(CFTG,0.0D0)*(CFF1/CFF2)*(VP(I+1,3,K)-VP(I,3,K)) ) / CMIU

          Vnow(I,3,K) = Vnow(I,3,K) + &
              DTC_OPENBC/Tau_obc * ( ( 0.2*VBCLN(I-1,1,K)+0.2*VBCLN(I+1,1,K)+ &
              0.4*VBCLN(I,1,K)+0.2*VBCLN(I,2,K) )-VP(I,3,K) ) * CFF1/CMIU

          Vnow(I,2,K) = Vnow(I,3,K) * VIV(I,2,K)
          Vnow(I,1,K) = Vnow(I,3,K) * VIV(I,1,K)
        ENDDO
        do i = 2,imt-1
          vnow(i,3,k) = 0.25d0*vnow(i-1,3,k) + 0.50d0*vnow(i,3,k) + 0.25d0*vnow(i+1,3,k)
        enddo
        Vnow(1,3,K) = Vnow(2,3,k) * VIV(1,3,K)
        Vnow(1,2,K) = Vnow(2,2,K) * VIV(1,2,K)
        Vnow(1,1,K) = Vnow(2,1,K) * VIV(1,1,K)
        Vnow(IMT,3,K) = Vnow(IMT-1,3,K) * VIV(IMT,3,K)
        Vnow(IMT,2,K) = Vnow(IMT-1,2,K) * VIV(IMT,2,K)
        Vnow(IMT,1,K) = Vnow(IMT-1,1,K) * VIV(IMT,1,K)
      ENDDO
    else if ( ID_OBC_V_N .eq. 5 ) then
      do k = 1, km
        do i = 2, imt-1
          dvdttmp  = vp(i,4,k)-vnow(i,4,k)
          dvdnmtmp = vnow (i,4,k)-vnow(i,5,k)

          if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
            tau_obc = tau_in_c5
            dvdttmp = 0.0d0
          else
            tau_obc = tau_out_c5
          endif
          tau_obc = tau_obc*dtc_openbc

          if (dvdttmp*(vp(i+1,4,k)-vp(i-1,4,k)) .gt. 0.0d0) then
            dvdtgtmp = vp(i  ,4,k) - vp(i-1,4,k)
          else
            dvdtgtmp = vp(i+1,4,k) - vp(i  ,4,k)
          endif

          cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
          cfnm = dvdttmp*dvdnmtmp
          if (ID_OBC_2DRad_NPO_UV5 == 0 ) then
            cftg = dvdttmp*dvdtgtmp
          else
            cftg = 0.0d0
          end if

          vnow(i,3,k) = (cff1*vp(i,3,k)+cfnm*vnow(i,4,k)-                   &
              &       max(cftg,0.0d0)*(vp(i  ,3,k)-vp(i-1,3,k))-              &
              &       min(cftg,0.0d0)*(vp(i+1,3,k)-vp(i  ,3,k)))/max(cff1+cfnm,1.0d-20)

          vnow(i,3,k) = vnow(i,3,k) + tau_obc * (( VBCLN(I,1,K) )-vnow(i,1,k))
          vnow(I,2,K) = vnow(I,3,K) * VIV(I,2,K)
          vnow(I,1,K) = vnow(I,3,K) * VIV(I,1,K)
        enddo
        if (ID_OBC_SMO_3DUV_NM .eq. 1 ) then
          do i = 2,imt-1
            vnow(i,3,k) = 0.50d0*vnow(i,3,k) + 0.50d0*vnow(i-1,3,k)
          enddo
        endif
        vnow(1,3,K) = vnow(2,3,k) * VIV(1,3,K)
        vnow(1,2,K) = vnow(2,2,K) * VIV(1,2,K)
        vnow(1,1,K) = vnow(2,1,K) * VIV(1,1,K)
        vnow(IMT,3,K) = vnow(IMT-1,3,K) * VIV(IMT,3,K)
        vnow(IMT,2,K) = vnow(IMT-1,2,K) * VIV(IMT,2,K)
        vnow(IMT,1,K) = vnow(IMT-1,1,K) * VIV(IMT,1,K)
      enddo
    ENDIF
  ENDIF

  !$OMP PARALLEL DO PRIVATE (K,J,I)
  DO K=1, KM
    DO J=1, JMT
      DO I=1, IMT
        Vnow(I,J,K) = Vnow(I,J,K) * VIV(I,J,K)
      ENDDO
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE OBC_3DV

SUBROUTINE OBC_3DTS(JEWSN)

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, IMM, JET, KM, I, J, K
  use   pconst_mod, only: VIT, DTS, DYR, OUX, OUY, ODZT
  use   dyn_mod,    only: V, U, WS
  use   tracer_mod, only: AT, ATB

#if ( defined SPMD )
  use msg_mod
#endif

#ifdef OBCDT
  use   openbc_mod
#endif

  IMPLICIT NONE
  INTEGER :: JEWSN
  REAL(r8) :: Tau_obc
  REAL(r8) :: CFF1, CFF2
  REAL(r8) :: UDTDX
  REAL(r8) :: DVDTTMP, DVDNMTMP, DVDTGTMP, DVDTG2TMP, CFTG, CFNM
  REAL(r8) :: EPSCF
  REAL(r8) :: CMIU
  INTEGER :: NOST
  EPSCF = 1.0D-20

  IF ( JEWSN .EQ. 1 ) THEN
    ! Eastern Boundary Condition

    IF ( ID_OBC_TS_E .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For T and S at Eastern Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,J)
        DO K = 1, KM
          DO J = 1, JMT
            AT(IMT-2,J,K,NOST) = AT(IMT-3,J,K,NOST) * VIT(IMT-2,J,K)
            AT(IMT-1,J,K,NOST) = AT(IMT-2,J,K,NOST) * VIT(IMT-1,J,K)
            AT(IMT,J,K,NOST) = AT(IMT-2,J,K,NOST) * VIT(IMT,J,K)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_TS_E .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For T and S at Eastern Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,J)
        DO K = 1, KM
          DO J = 1, JMT
            IF ( NOST.EQ.1 ) THEN
              AT(IMT-2,J,K,NOST) = TBCLE(J,K) * VIT(IMT-2,J,K)
            ELSE IF ( NOST.EQ.2 ) THEN
              AT(IMT-2,J,K,NOST) = SBCLE(J,K) * VIT(IMT-2,J,K)
            ENDIF

            AT(IMT-1,J,K,NOST) = AT(IMT-2,J,K,NOST) * VIT(IMT-1,J,K)
            AT(IMT,J,K,NOST) = AT(IMT-2,J,K,NOST) * VIT(IMT,J,K)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_TS_E .EQ. 3 ) THEN
      ! The upstream advection scheme for T and S at the eastern edge
      ! This scheme follows formula (B-3) in the POM2K manual:
      ! Tt + U Tx = 0
      ! See Mellor (2003), Users Guide for a Three-Dimensional Primitive
      ! Equation Numerical Ocean Model, and Chen et al. (2013), equation (8).

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 2, JMT-1
          UDTDX = 0.5D0*( U(IMT-2,J-1,K) + U(IMT-2,J,K) ) * DTS * OUX(J)
          IF ( UDTDX .LE. 0.0D0 ) THEN
            AT(IMT-2,J,K,1) = ATB(IMT-2,J,K,1) - UDTDX * ( TBCLE(J,K) -  &
                ATB(IMT-2,J,K,1) )
            AT(IMT-2,J,K,2) = ATB(IMT-2,J,K,2) - UDTDX * ( SBCLE(J,K) -  &
                ATB(IMT-2,J,K,2) )
          ELSE
            AT(IMT-2,J,K,1) = ATB(IMT-2,J,K,1) - UDTDX * ( ATB(IMT-2,J,K,1) &
                - ATB(IMT-3,J,K,1) )
            AT(IMT-2,J,K,2) = ATB(IMT-2,J,K,2) - UDTDX * ( ATB(IMT-2,J,K,2) &
                - ATB(IMT-3,J,K,2) )

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(IMT-2,J,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(IMT-2,J,K-1,1) = AT(IMT-2,J,K-1,1) - CMIU * ( &
                    ATB(IMT-2,J,K-1,1) - ATB(IMT-2,J,K,1) )
                AT(IMT-2,J,K-1,2) = AT(IMT-2,J,K-1,2) - CMIU * ( &
                    ATB(IMT-2,J,K-1,2) - ATB(IMT-2,J,K,2) )
              ELSE
                AT(IMT-2,J,K,1) = AT(IMT-2,J,K,1) - CMIU * ( &
                    ATB(IMT-2,J,K-1,1) - ATB(IMT-2,J,K,1) )
                AT(IMT-2,J,K,2) = AT(IMT-2,J,K,2) - CMIU * ( &
                    ATB(IMT-2,J,K-1,2) - ATB(IMT-2,J,K,2) )
              ENDIF
            ENDIF
          ENDIF
          AT(IMT-1,J,K,1) = AT(IMT-2,J,K,1) * VIT(IMT-1,J,K)
          AT(IMT,J,K,1) = AT(IMT-2,J,K,1) * VIT(IMT,J,K)
          AT(IMT-1,J,K,2) = AT(IMT-2,J,K,2) * VIT(IMT-1,J,K)
          AT(IMT,J,K,2) = AT(IMT-2,J,K,2) * VIT(IMT,J,K)
        ENDDO
        AT(IMT-2,1,K,1) = AT(IMT-2,2,K,1) * VIT(IMT-2,1,K)
        AT(IMT-1,1,K,1) = AT(IMT-1,2,K,1) * VIT(IMT-1,1,K)
        AT(IMT,1,K,1) = AT(IMT,2,K,1) * VIT(IMT,1,K)
        AT(IMT-2,1,K,2) = AT(IMT-2,2,K,2) * VIT(IMT-2,1,K)
        AT(IMT-1,1,K,2) = AT(IMT-1,2,K,2) * VIT(IMT-1,1,K)
        AT(IMT,1,K,2) = AT(IMT,2,K,2) * VIT(IMT,1,K)
        AT(IMT-2,JMT,K,1) = AT(IMT-2,JMT-1,K,1) * VIT(IMT-2,JMT,K)
        AT(IMT-1,JMT,K,1) = AT(IMT-1,JMT-1,K,1) * VIT(IMT-1,JMT,K)
        AT(IMT,JMT,K,1) = AT(IMT,JMT-1,K,1) * VIT(IMT,JMT,K)
        AT(IMT-2,JMT,K,2) = AT(IMT-2,JMT-1,K,2) * VIT(IMT-2,JMT,K)
        AT(IMT-1,JMT,K,2) = AT(IMT-1,JMT-1,K,2) * VIT(IMT-1,JMT,K)
        AT(IMT,JMT,K,2) = AT(IMT,JMT-1,K,2) * VIT(IMT,JMT,K)
      ENDDO
    ELSE IF ( ID_OBC_TS_E .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For T and S at Eastern Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,J)
        DO K = 1, KM
          DO J =2, JMT-1
            DVDTTMP = ATB(IMT-3,J,K,NOST) - AT(IMT-3,J,K,NOST)
            DVDNMTMP = AT(IMT-3,J,K,NOST) - AT(IMT-4,J,K,NOST)
            DVDTG2TMP = ATB(IMT-3,J+1,K,NOST) - ATB(IMT-3,J-1,K,NOST)

            IF ( DVDTTMP*DVDNMTMP .LT. 0.0D0 ) THEN
              DVDTTMP = 0.0D0
              Tau_obc = TAUIN_TS
            ELSE
              Tau_obc = TAUOUT_TS
            ENDIF

            IF( ( DVDTTMP*DVDTG2TMP).GT.0.0) THEN
              DVDTGTMP=ATB(IMT-3,J,K,NOST) - ATB(IMT-3,J-1,K,NOST)
            ELSE
              DVDTGTMP=ATB(IMT-3,J+1,K,NOST) - ATB(IMT-3,J,K,NOST)
            ENDIF

            IF ( ID_OBC_2DRad_TS .EQ. 0 ) THEN
              DVDTGTMP = 0.0D0
            ENDIF

            IF( (DVDTTMP*DVDTG2TMP).GT.0.0) THEN
              CFF1 = DVDNMTMP * DVDNMTMP + &
                  DVDTGTMP * DVDTGTMP / ( DYR(J-1)*DYR(J-1) ) * &
                  ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) * ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) )

              CFF2 = DVDTGTMP * DVDTGTMP + &
                  DVDNMTMP * DVDNMTMP * DYR(J-1) * DYR(J-1) / &
                  ( ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) *( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) )
            ELSE
              CFF1 = DVDNMTMP * DVDNMTMP + &
                  DVDTGTMP * DVDTGTMP / ( DYR(J)*DYR(J) ) * &
                  ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) * ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) )
              CFF2 = DVDTGTMP * DVDTGTMP + &
                  DVDNMTMP * DVDNMTMP * DYR(J) * DYR(J) / &
                  ( ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) *( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) )
            ENDIF

            CFF1 = MAX(CFF1,EPSCF)
            CFF2 = MAX(CFF2,EPSCF)
            CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
            CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

            IF ( ID_OBC_2DRad_NPO_TS .EQ. 1 ) THEN
              CFTG = 0.0D0
            ENDIF

            CMIU = CFF1 + CFNM
            AT(IMT-2,J,K,NOST) = ( CFF1*ATB(IMT-2,J,K,NOST) + &
                CFNM*AT(IMT-3,J,K,NOST) -  &
                MAX(CFTG,0.0D0)*(CFF1/CFF2)*(ATB(IMT-2,J,K,NOST)-ATB(IMT-2,J-1,K,NOST)) - &
                MIN(CFTG,0.0D0)*(CFF1/CFF2)*(ATB(IMT-2,J+1,K,NOST)-ATB(IMT-2,J,K,NOST)) )/CMIU

            IF (NOST.EQ.1) THEN
              AT(IMT-2,J,K,NOST) = AT(IMT-2,J,K,NOST) + &
                  DTS/Tau_obc * ( TBCLE(J,K)-ATB(IMT-2,J,K,NOST) ) * CFF1/CMIU
            ELSE IF (NOST.EQ.2) THEN
              AT(IMT-2,J,K,NOST) = AT(IMT-2,J,K,NOST) + &
                  DTS/Tau_obc * ( SBCLE(J,K)-ATB(IMT-2,J,K,NOST) ) * CFF1/CMIU
            ENDIF

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(IMT-2,J,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(IMT-2,J,K-1,NOST) = AT(IMT-2,J,K-1,NOST) - &
                    CMIU * ( ATB(IMT-2,J,K-1,NOST) - ATB(IMT-2,J,K,NOST) )
              ELSE
                AT(IMT-2,J,K,NOST) = AT(IMT-2,J,K,NOST) - &
                    CMIU * ( ATB(IMT-2,J,K-1,NOST) - ATB(IMT-2,J,K,NOST) )
              ENDIF
            ENDIF

            AT(IMT-1,J,K,NOST) = AT(IMT-2,J,K,NOST) * VIT(IMT-1,J,K)
            AT(IMT,J,K,NOST) = AT(IMT-2,J,K,NOST) * VIT(IMT-1,J,K)
          ENDDO

          AT(IMT-2,1,K,NOST) = AT(IMT-2,2,K,NOST) * VIT(IMT-2,1,K)
          AT(IMT-1,1,K,NOST) = AT(IMT-1,2,K,NOST) * VIT(IMT-1,1,K)
          AT(IMT,1,K,NOST) = AT(IMT,2,K,NOST) * VIT(IMT,1,K)
          AT(IMT-2,JMT,K,NOST) = AT(IMT-2,JMT-1,K,NOST) * VIT(IMT-2,JMT,K)
          AT(IMT-1,JMT,K,NOST) = AT(IMT-1,JMT-1,K,NOST) * VIT(IMT-1,JMT,K)
          AT(IMT,JMT,K,NOST) = AT(IMT,JMT-1,K,NOST) * VIT(IMT,JMT,K)
        ENDDO
      ENDDO
    else if ( ID_OBC_TS_E .eq. 5 ) then
      do nost = 1, 2
        do k = 1, km
          do j = 2, jmt-1
            dvdttmp =  atb(imt-3,j  ,k,nost) - at(imt-3,j  ,k,nost)
            dvdnmtmp = at (imt-3,j  ,k,nost) - at(imt-4,j  ,k,nost)

            if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
              tau_obc = tau_in_c5
              dvdttmp = 0.0d0
            else
              tau_obc = tau_out_c5
            endif
            tau_obc = tau_obc*dtc_openbc

            if (dvdttmp*(atb(imt-3,j+1,k,nost) - atb(imt-3,j-1,k,nost)) .gt. 0.0d0) then
              dvdtgtmp = atb(imt-3,j  ,k,nost) - atb(imt-3,j-1,k,nost)
            else
              dvdtgtmp = atb(imt-3,j+1,k,nost) - atb(imt-3,j  ,k,nost)
            end if

            cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp, 3.0d-7)
            cfnm = dvdttmp*dvdnmtmp
            if (ID_OBC_2DRad_NPO_TS5 == 0 ) then
              cftg = dvdttmp*dvdtgtmp
            else
              cftg = 0.0d0
            end if

            at(imt-2,j  ,k,nost) = (cff1*atb(imt-2,j  ,k,nost)+cfnm*at(imt-3,j  ,k,nost)-         &
                &       max(cftg,0.0d0)*(atb(imt-2,j  ,k,nost)-atb(imt-2,j-1,k,nost))-              &
                &       min(cftg,0.0d0)*(atb(imt-2,j+1,k,nost)-atb(imt-2,j  ,k,nost)))/max(cff1+cfnm,1.0d-20)

            if (nost .eq. 1) then
              at(imt-2,j  ,k,nost) = at(imt-2,j  ,k,nost) + tau_obc * ((TBCLE(J,K))-at(imt-2,j,k,nost))
            elseif (nost .eq. 2) then
              at(imt-2,j  ,k,nost) = at(imt-2,j  ,k,nost) + tau_obc * ((SBCLE(J,K))-at(imt-2,j,k,nost))
            endif

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(IMT-2,J,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(IMT-2,J,K-1,NOST) = AT(IMT-2,J,K-1,NOST) - &
                    CMIU * ( ATB(IMT-2,J,K-1,NOST) - ATB(IMT-2,J,K,NOST) )
              ELSE
                AT(IMT-2,J,K,NOST) = AT(IMT-2,J,K,NOST) - &
                    CMIU * ( ATB(IMT-2,J,K-1,NOST) - ATB(IMT-2,J,K,NOST) )
              ENDIF
            ENDIF

            at(imt-1,j,k,nost) = at(imt-2,j,k,nost) * VIT(IMT-1,J,K)
            at(imt  ,j,k,nost) = at(imt-2,j,k,nost) * VIT(IMT  ,J,K)
          enddo

          AT(IMT-2,1,K,NOST) = AT(IMT-2,2,K,NOST) * VIT(IMT-2,1,K)
          AT(IMT-1,1,K,NOST) = AT(IMT-1,2,K,NOST) * VIT(IMT-1,1,K)
          AT(IMT,1,K,NOST) = AT(IMT,2,K,NOST) * VIT(IMT,1,K)
          AT(IMT-2,JMT,K,NOST) = AT(IMT-2,JMT-1,K,NOST) * VIT(IMT-2,JMT,K)
          AT(IMT-1,JMT,K,NOST) = AT(IMT-1,JMT-1,K,NOST) * VIT(IMT-1,JMT,K)
          AT(IMT,JMT,K,NOST) = AT(IMT,JMT-1,K,NOST) * VIT(IMT,JMT,K)
        enddo
      enddo
    ENDIF
  ELSE IF ( JEWSN .EQ. 2 ) THEN
    ! Western Boundary Conditions

    IF ( ID_OBC_TS_W .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For T and S at Western Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,J)
        DO K = 1, KM
          DO J = 1, JMT
            AT(3,J,K,NOST) = AT(4,J,K,NOST) * VIT(3,J,K)
            AT(2,J,K,NOST) = AT(3,J,K,NOST) * VIT(2,J,K)
            AT(1,J,K,NOST) = AT(3,J,K,NOST) * VIT(1,J,K)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_TS_W .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For T and S at Western Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,J)
        DO K = 1, KM
          DO J = 1, JMT
            IF ( NOST.EQ.1 ) THEN
              AT(3,J,K,NOST) = TBCLW(J,K) * VIT(3,J,K)
            ELSE IF ( NOST.EQ.2 ) THEN
              AT(3,J,K,NOST) = SBCLW(J,K) * VIT(3,J,K)
            ENDIF

            AT(2,J,K,NOST) = AT(3,J,K,NOST) * VIT(2,J,K)
            AT(1,J,K,NOST) = AT(3,J,K,NOST) * VIT(1,J,K)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_TS_W .EQ.3 ) THEN
      ! The Upstream Advection Scheme For T and S at Western Edge

      !$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1, KM
        DO J = 2, JMT-1
          UDTDX = 0.5D0*( U(4,J-1,K) + U(4,J,K) ) * DTS * OUX(J)
          IF ( UDTDX .GE. 0.0D0 ) THEN
            AT(3,J,K,1) = ATB(3,J,K,1) - UDTDX * ( ATB(3,J,K,1) - TBCLW(J,K) )
            AT(3,J,K,2) = ATB(3,J,K,2) - UDTDX * ( ATB(3,J,K,2) - SBCLW(J,K) )
          ELSE
            AT(3,J,K,1) = ATB(3,J,K,1) - UDTDX * ( ATB(4,J,K,1) - ATB(3,J,K,1) )
            AT(3,J,K,2) = ATB(3,J,K,2) - UDTDX * ( ATB(4,J,K,2) - ATB(3,J,K,2) )

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(3,J,K) * DTS * ODZT(K-1)
              IF ( CMIU.GE.0.0 ) THEN
                AT(3,J,K-1,1) = AT(3,J,K-1,1) - CMIU * ( ATB(3,J,K-1,1) - ATB(3,J,K,1) )
                AT(3,J,K-1,2) = AT(3,J,K-1,2) - CMIU * ( ATB(3,J,K-1,2) - ATB(3,J,K,2) )
              ELSE
                AT(3,J,K,1) = AT(3,J,K,1) - CMIU * ( ATB(3,J,K-1,1) - ATB(3,J,K,1) )
                AT(3,J,K,2) = AT(3,J,K,2) - CMIU * ( ATB(3,J,K-1,2) - ATB(3,J,K,2) )
              ENDIF
            ENDIF
          ENDIF
          AT(2,J,K,1) = AT(3,J,K,1)*VIT(2,J,K)
          AT(1,J,K,1) = AT(3,J,K,1)*VIT(1,J,K)
          AT(2,J,K,2) = AT(3,J,K,2)*VIT(2,J,K)
          AT(1,J,K,2) = AT(3,J,K,2)*VIT(1,J,K)
        ENDDO
        AT(3,1,K,1) = AT(3,2,K,1)*VIT(3,1,K)
        AT(2,1,K,1) = AT(2,2,K,1)*VIT(2,1,K)
        AT(1,1,K,1) = AT(1,2,K,1)*VIT(1,1,K)
        AT(3,1,K,2) = AT(3,2,K,2)*VIT(3,1,K)
        AT(2,1,K,2) = AT(2,2,K,2)*VIT(2,1,K)
        AT(1,1,K,2) = AT(1,2,K,2)*VIT(1,1,K)
        AT(3,JMT,K,1) = AT(3,JMT-1,K,1)*VIT(3,JMT,K)
        AT(2,JMT,K,1) = AT(2,JMT-1,K,1)*VIT(2,JMT,K)
        AT(1,JMT,K,1) = AT(1,JMT-1,K,1)*VIT(1,JMT,K)
        AT(3,JMT,K,2) = AT(3,JMT-1,K,2)*VIT(3,JMT,K)
        AT(2,JMT,K,2) = AT(2,JMT-1,K,2)*VIT(2,JMT,K)
        AT(1,JMT,K,2) = AT(1,JMT-1,K,2)*VIT(1,JMT,K)
      ENDDO
    ELSE IF ( ID_OBC_TS_W .EQ.4 ) THEN
      ! The Implicit Upstream Radiation Condition For T and S at Western Edge

      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,J)
        DO K = 1, KM
          DO J =2, JMT-1
            DVDTTMP = ATB(4,J,K,NOST) - AT(4,J,K,NOST)
            DVDNMTMP = AT(4,J,K,NOST) - AT(5,J,K,NOST)
            DVDTG2TMP = ATB(4,J+1,K,NOST) - ATB(4,J-1,K,NOST)

            IF ( ( DVDTTMP*DVDNMTMP ) .LT. 0.0D0 ) THEN
              DVDTTMP = 0.0D0
              Tau_obc = TAUIN_TS
            ELSE
              Tau_obc = TAUOUT_TS
            ENDIF

            IF( ( DVDTTMP*DVDTG2TMP ).GT.0.0) THEN
              DVDTGTMP=ATB(4,J,K,NOST) - ATB(4,J-1,K,NOST)
            ELSE
              DVDTGTMP=ATB(4,J+1,K,NOST) - ATB(4,J,K,NOST)
            ENDIF

            IF ( ID_OBC_2DRad_TS .EQ. 0 ) THEN
              DVDTGTMP = 0.0D0
            ENDIF

            IF( (DVDTTMP*DVDTG2TMP).GT.0.0) THEN
              CFF1 = DVDNMTMP * DVDNMTMP + &
                  DVDTGTMP * DVDTGTMP / ( DYR(J-1)*DYR(J-1) ) * &
                  ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) * ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) )

              CFF2 = DVDTGTMP * DVDTGTMP + &
                  DVDNMTMP * DVDNMTMP * DYR(J-1) * DYR(J-1) / &
                  ( ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) * ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) )
            ELSE
              CFF1 = DVDNMTMP * DVDNMTMP + &
                  DVDTGTMP * DVDTGTMP / ( DYR(J)*DYR(J) ) * &
                  ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) * ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) )

              CFF2 = DVDTGTMP * DVDTGTMP + &
                  DVDNMTMP * DVDNMTMP * DYR(J) * DYR(J) / &
                  ( ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) * ( 0.50D0/OUX(J-1)+0.50D0/OUX(J) ) )
            ENDIF

            CFF1 = MAX(CFF1,EPSCF)
            CFF2 = MAX(CFF2,EPSCF)
            CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
            CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

            IF ( ID_OBC_2DRad_NPO_TS .EQ. 1 ) then
              CFTG = 0.0D0
            ENDIF

            CMIU = CFF1 + CFNM
            AT(3,J,K,NOST) = ( CFF1*ATB(3,J,K,NOST) + &
                CFNM*AT(4,J,K,NOST) -  &
                MAX(CFTG,0.0D0)*(CFF1/CFF2)*(ATB(3,J,K,NOST)-ATB(3,J-1,K,NOST)) - &
                MIN(CFTG,0.0D0)*(CFF1/CFF2)*(ATB(3,J+1,K,NOST)-ATB(3,J,K,NOST)) )/CMIU

            IF ( NOST.EQ.1) THEN
              AT(3,J,K,NOST) = AT(3,J,K,NOST) + &
                  DTS/Tau_obc*CFF1*(TBCLW(J,K)-ATB(3,J,K,NOST))/CMIU
            ELSE IF (NOST.EQ.2) THEN
              AT(3,J,K,NOST) = AT(3,J,K,NOST) + &
                  DTS/Tau_obc * CFF1 * ( SBCLW(J,K)-ATB(3,J,K,NOST) )/CMIU
            ENDIF

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(3,J,K) * DTS * ODZT(K-1)
              IF ( CMIU.GE.0.0 ) THEN
                AT(3,J,K-1,1) = AT(3,J,K-1,1) - CMIU * ( ATB(3,J,K-1,1) - ATB(3,J,K,1) )
                AT(3,J,K-1,2) = AT(3,J,K-1,2) - CMIU * ( ATB(3,J,K-1,2) - ATB(3,J,K,2) )
              ELSE
                AT(3,J,K,1) = AT(3,J,K,1) - CMIU * ( ATB(3,J,K-1,1) - ATB(3,J,K,1) )
                AT(3,J,K,2) = AT(3,J,K,2) - CMIU * ( ATB(3,J,K-1,2) - ATB(3,J,K,2) )
              ENDIF
            ENDIF

            AT(2,J,K,NOST) = AT(3,J,K,NOST) * VIT(2,J,K)
            AT(1,J,K,NOST) = AT(3,J,K,NOST) * VIT(1,J,K)
          ENDDO
          AT(3,1,K,NOST) = 0.5*AT(3,2,K,NOST)*VIT(3,1,K) + 0.5*ATB(3,1,K,NOST)
          AT(2,1,K,NOST) = 0.5*AT(2,2,K,NOST)*VIT(2,1,K) + 0.5*ATB(2,1,K,NOST)
          AT(1,1,K,NOST) = 0.5*AT(1,2,K,NOST)*VIT(1,1,K) + 0.5*ATB(1,1,K,NOST)
          AT(3,JMT,K,NOST) = 0.5*AT(3,JMT-1,K,NOST)*VIT(3,JMT,K) + 0.5*ATB(3,JMT,K,NOST)
          AT(2,JMT,K,NOST) = 0.5*AT(2,JMT-1,K,NOST)*VIT(2,JMT,K) + 0.5*ATB(2,JMT,K,NOST)
          AT(1,JMT,K,NOST) = 0.5*AT(1,JMT-1,K,NOST)*VIT(1,JMT,K) + 0.5*ATB(1,JMT,K,NOST)
        ENDDO
      ENDDO
    else if ( ID_OBC_TS_W .eq. 5 ) then
      do nost = 1, 2
        do k = 1, km
          do j = 2, jmt-1
            dvdttmp  = atb(4,j,k,nost) - at(4,j,k,nost)
            dvdnmtmp = at (4,j,k,nost) - at(5,j,k,nost)

            if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
              tau_obc = tau_in_c5
              dvdttmp = 0.0d0
            else
              tau_obc = tau_out_c5
            endif
            tau_obc = tau_obc*dtc_openbc

            if (dvdttmp*(atb(4,j+1,k,nost)-atb(4,j-1,k,nost)) .gt. 0.0d0) then
              dvdtgtmp = atb(4,j  ,k,nost) - atb(4,j-1,k,nost)
            else
              dvdtgtmp = atb(4,j+1,k,nost) - atb(4,j  ,k,nost)
            endif

            cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
            cfnm = dvdttmp*dvdnmtmp  ! here the phase v is actually -Cnm
            if (ID_OBC_2DRad_NPO_TS5 == 0 ) then
              cftg = dvdttmp*dvdtgtmp
            else
              cftg = 0.0d0
            end if

            at(3,j,k,nost) = (cff1*atb(3,j,k,nost)+cfnm*at(4,j,k,nost)-                   &
                &       max(cftg,0.0d0)*(atb(3,j  ,k,nost)-atb(3,j-1,k,nost))-              &
                &       min(cftg,0.0d0)*(atb(3,j+1,k,nost)-atb(3,j  ,k,nost)))/max(cff1+cfnm,1.0d-20)

            if (nost .eq. 1) then
              at(3,j,k,nost) = at(3,j,k,nost) + tau_obc * ((TBCLW(J,K))-at(3,j,k,nost))
            elseif (nost .eq. 2) then
              at(3,j,k,nost) = at(3,j,k,nost) + tau_obc * ((SBCLW(J,K))-at(3,j,k,nost))
            endif

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(3,J,K) * DTS * ODZT(K-1)
              IF ( CMIU.GE.0.0 ) THEN
                AT(3,J,K-1,1) = AT(3,J,K-1,1) - CMIU * ( ATB(3,J,K-1,1) - ATB(3,J,K,1) )
                AT(3,J,K-1,2) = AT(3,J,K-1,2) - CMIU * ( ATB(3,J,K-1,2) - ATB(3,J,K,2) )
              ELSE
                AT(3,J,K,1) = AT(3,J,K,1) - CMIU * ( ATB(3,J,K-1,1) - ATB(3,J,K,1) )
                AT(3,J,K,2) = AT(3,J,K,2) - CMIU * ( ATB(3,J,K-1,2) - ATB(3,J,K,2) )
              ENDIF
            ENDIF

            at(2,J,K,nost) = at(3,J,K,nost) * VIT(2,J,K)
            at(1,J,K,nost) = at(3,J,K,nost) * VIT(1,J,K)
          enddo

          AT(3,1,K,NOST) = AT(3,2,K,NOST)*VIT(3,1,K)
          AT(2,1,K,NOST) = AT(2,2,K,NOST)*VIT(2,1,K)
          AT(1,1,K,NOST) = AT(1,2,K,NOST)*VIT(1,1,K)
          AT(3,JMT,K,NOST) = AT(3,JMT-1,K,NOST)*VIT(3,JMT,K)
          AT(2,JMT,K,NOST) = AT(2,JMT-1,K,NOST)*VIT(2,JMT,K)
          AT(1,JMT,K,NOST) = AT(1,JMT-1,K,NOST)*VIT(1,JMT,K)
        enddo
      enddo
    ENDIF
  ELSE IF ( JEWSN .EQ. 3 ) THEN
    ! Southern Boundary Conditions
    IF ( ID_OBC_TS_S .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For T and S at Southern Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,I)
        DO K = 1, KM
          DO I = 1, IMT
            AT(I,JET-2,K,NOST) = AT(I,JET-3,K,NOST) * VIT(I,JET-2,K)
            AT(I,JET-1,K,NOST) = AT(I,JET-2,K,NOST) * VIT(I,JET-1,K)
            AT(I,JET,K,NOST) = AT(I,JET-2,K,NOST) * VIT(I,JET,K)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_TS_S .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For T and S at Southern Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,I)
        DO K = 1, KM
          DO I = 1, IMT
            IF ( NOST.EQ.1 ) THEN
              AT(I,JET-2,K,NOST) = TBCLS(I,K) * VIT(I,JET-2,K)
            ELSE IF ( NOST.EQ.2 ) THEN
              AT(I,JET-2,K,NOST) = SBCLS(I,K) * VIT(I,JET-2,K)
            ENDIF

            AT(I,JET-1,K,NOST) = AT(I,JET-2,K,NOST) * VIT(I,JET-1,K)
            AT(I,JET,K,NOST) = AT(I,JET-2,K,NOST) * VIT(I,JET,K)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_TS_S .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For T and S at Southern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 2, IMM
          UDTDX = 0.5D0 * ( V(I,JET-3,K)+V(I+1,JET-3,K) ) * DTS * OUY(JET-3)

          IF ( UDTDX.LE.0.0D0 ) THEN
            AT(I,JET-2,K,1) = ATB(I,JET-2,K,1) - UDTDX * ( TBCLS(I,K) - ATB(I,JET-2,K,1) )
            AT(I,JET-2,K,2) = ATB(I,JET-2,K,2) - UDTDX * ( SBCLS(I,K) - ATB(I,JET-2,K,2) )
          ELSE
            AT(I,JET-2,K,1) = ATB(I,JET-2,K,1) - UDTDX * ( ATB(I,JET-2,K,1) - ATB(I,JET-3,K,1) )
            AT(I,JET-2,K,2) = ATB(I,JET-2,K,2) - UDTDX * ( ATB(I,JET-2,K,2) - ATB(I,JET-3,K,2) )

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(I,JET-2,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(I,JET-2,K-1,1) = AT(I,JET-2,K-1,1) - CMIU * ( ATB(I,JET-2,K-1,1) - ATB(I,JET-2,K,1) )
                AT(I,JET-2,K-1,2) = AT(I,JET-2,K-1,2) - CMIU * ( ATB(I,JET-2,K-1,2) - ATB(I,JET-2,K,2) )
              ELSE
                AT(I,JET-2,K,1) = AT(I,JET-2,K,1) - CMIU * ( ATB(I,JET-2,K-1,1) - ATB(I,JET-2,K,1) )
                AT(I,JET-2,K,2) = AT(I,JET-2,K,2) - CMIU * ( ATB(I,JET-2,K-1,2) - ATB(I,JET-2,K,2) )
              ENDIF
            ENDIF
          ENDIF
          AT(I,JET-1,K,1) = AT(I,JET-2,K,1) * VIT(I,JET-1,K)
          AT(I,JET,K,1) = AT(I,JET-2,K,1) * VIT(I,JET,K)
          AT(I,JET-1,K,2) = AT(I,JET-2,K,2) * VIT(I,JET-1,K)
          AT(I,JET,K,2) = AT(I,JET-2,K,2) * VIT(I,JET,K)
        ENDDO
        AT(1,JET-2,K,1) = AT(2,JET-2,K,1) * VIT(1,JET-2,K)
        AT(1,JET-1,K,1) = AT(2,JET-1,K,1) * VIT(1,JET-1,K)
        AT(1,JET,K,1) = AT(2,JET,K,1) * VIT(1,JET,K)
        AT(IMT,JET-2,K,1) = AT(IMM,JET-2,K,1) * VIT(IMT,JET-2,K)
        AT(IMT,JET-1,K,1) = AT(IMM,JET-1,K,1) * VIT(IMT,JET-1,K)
        AT(IMT,JET,K,1) = AT(IMM,JET,K,1) * VIT(IMT,JET,K)
        AT(1,JET-2,K,2) = AT(2,JET-2,K,2) * VIT(1,JET-2,K)
        AT(1,JET-1,K,2) = AT(2,JET-1,K,2) * VIT(1,JET-1,K)
        AT(1,JET,K,2) = AT(2,JET,K,2) * VIT(1,JET,K)
        AT(IMT,JET-2,K,2) = AT(IMM,JET-2,K,2) * VIT(IMT,JET-2,K)
        AT(IMT,JET-1,K,2) = AT(IMM,JET-1,K,2) * VIT(IMT,JET-1,K)
        AT(IMT,JET,K,2) = AT(IMM,JET,K,2) * VIT(IMT,JET,K)
      ENDDO
    ELSE IF ( ID_OBC_TS_S .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For T and S at Southern Edge

      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,I)
        DO K = 1, KM
          DO I =2, IMT-1
            DVDTTMP = ATB(I,JET-3,K,NOST) - AT(I,JET-3,K,NOST)
            DVDNMTMP = AT(I,JET-3,K,NOST) - AT(I,JET-4,K,NOST)
            DVDTG2TMP = ATB(I+1,JET-3,K,NOST) - ATB(I-1,JET-3,K,NOST)

            IF ( ( DVDTTMP*DVDNMTMP ) .LT. 0.0 ) THEN
              DVDTTMP = 0.0D0
              Tau_obc = TAUIN_TS
            ELSE
              Tau_obc = TAUOUT_TS
            ENDIF

            IF( ( DVDTTMP*DVDTG2TMP ).GT.0.0) THEN
              DVDTGTMP = ATB(I,JET-3,K,NOST) - ATB(I-1,JET-3,K,NOST)
            ELSE
              DVDTGTMP = ATB(I+1,JET-3,K,NOST) - ATB(I,JET-3,K,NOST)
            ENDIF

            IF ( ID_OBC_2DRad_TS .EQ. 0 ) THEN
              DVDTGTMP = 0.0D0
            ENDIF

            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP * DYR(JET-4) * DYR(JET-4) / &
                ( ( 0.50D0/OUX(JET-4)+0.50D0/OUX(JET-3) ) * ( 0.50D0/OUX(JET-4)+0.50D0/OUX(JET-3) ) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP / ( DYR(JET-4)*DYR(JET-4) ) * &
                ( 0.50D0/OUX(JET-4)+0.50D0/OUX(JET-3) ) * ( 0.50D0/OUX(JET-4)+0.50D0/OUX(JET-3) )

            CFF1 = MAX(CFF1,EPSCF)
            CFF2 = MAX(CFF2,EPSCF)
            CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
            CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

            IF ( ID_OBC_2DRad_NPO_TS .EQ. 1 ) THEN
              CFTG = 0.0D0
            ENDIF

            CMIU = CFF1 + CFNM
            AT(I,JET-2,K,NOST) = ( CFF1*ATB(I,JET-2,K,NOST) + &
                CFNM*AT(I,JET-3,K,NOST) -  &
                MAX(CFTG,0.0D0)*(CFF1/CFF2)*(ATB(I,JET-2,K,NOST)-ATB(I-1,JET-2,K,NOST)) - &
                MIN(CFTG,0.0D0)*(CFF1/CFF2)*(ATB(I+1,JET-2,K,NOST)-ATB(I,JET-2,K,NOST)) )/CMIU

            IF (NOST.EQ.1) THEN
              AT(I,JET-2,K,NOST) = AT(I,JET-2,K,NOST) + &
                  DTS/Tau_obc*CFF1*(TBCLS(I,K)-ATB(I,JET-2,K,NOST))/CMIU
            ELSE IF (NOST.EQ.2) THEN
              AT(I,JET-2,K,NOST) = AT(I,JET-2,K,NOST) + &
                  DTS/Tau_obc*CFF1*(SBCLS(I,K)-ATB(I,JET-2,K,NOST))/CMIU
            ENDIF

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(I,JET-2,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(I,JET-2,K-1,1) = AT(I,JET-2,K-1,1) - CMIU * ( ATB(I,JET-2,K-1,1) - ATB(I,JET-2,K,1) )
                AT(I,JET-2,K-1,2) = AT(I,JET-2,K-1,2) - CMIU * ( ATB(I,JET-2,K-1,2) - ATB(I,JET-2,K,2) )
              ELSE
                AT(I,JET-2,K,1) = AT(I,JET-2,K,1) - CMIU * ( ATB(I,JET-2,K-1,1) - ATB(I,JET-2,K,1) )
                AT(I,JET-2,K,2) = AT(I,JET-2,K,2) - CMIU * ( ATB(I,JET-2,K-1,2) - ATB(I,JET-2,K,2) )
              ENDIF
            ENDIF

            AT(I,JET-1,K,NOST) = AT(I,JET-2,K,NOST) * VIT(I,JET-1,K)
            AT(I,JET,K,NOST) = AT(I,JET-2,K,NOST) * VIT(I,JET,K)
          ENDDO
          AT(1,JET-2,K,NOST) = 0.5*AT(2,JET-2,K,NOST)*VIT(1,JET-2,K) + 0.5*ATB(1,JET-2,K,NOST)
          AT(1,JET-1,K,NOST) = 0.5*AT(2,JET-1,K,NOST)*VIT(1,JET-1,K) + 0.5*ATB(1,JET-1,K,NOST)
          AT(1,JET,K,NOST) = 0.5*AT(2,JET,K,NOST)*VIT(1,JET,K) + 0.5*ATB(1,JET,K,NOST)
          AT(IMT,JET-2,K,NOST) = 0.5*AT(IMT-1,JET-2,K,NOST)*VIT(IMT,JET-2,K) + 0.5*ATB(IMT,JET-2,K,NOST)
          AT(IMT,JET-1,K,NOST) = 0.5*AT(IMT-1,JET-1,K,NOST)*VIT(IMT,JET-1,K) + 0.5*ATB(IMT,JET-1,K,NOST)
          AT(IMT,JET,K,NOST) = 0.5*AT(IMT-1,JET,K,NOST)*VIT(IMT,JET,K) + 0.5*ATB(IMT,JET,K,NOST)
        ENDDO
      ENDDO
    else if ( ID_OBC_TS_S .eq. 5 ) then
      do nost = 1, 2
        do k = 1, km
          do i = 2,imt-1
            dvdttmp  = atb(i  ,jet-3,k,nost)-at(i  ,jet-3,k,nost)
            dvdnmtmp = at (i  ,jet-3,k,nost)-at(i  ,jet-4,k,nost)

            if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
              tau_obc = tau_in_c5
              dvdttmp = 0.0d0
            else
              tau_obc = tau_out_c5
            endif
            tau_obc  = tau_obc*dtc_openbc

            if (dvdttmp*(atb(i+1,jet-3,k,nost)-atb(i-1,jet-3,k,nost)) .gt. 0.0d0) then
              dvdtgtmp = atb(i  ,jet-3,k,nost)-atb(i-1,jet-3,k,nost)
            else
              dvdtgtmp = atb(i+1,jet-3,k,nost)-atb(i  ,jet-3,k,nost)
            endif

            cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
            cfnm = dvdttmp*dvdnmtmp
            if (ID_OBC_2DRad_NPO_TS5 == 0 ) then
              cftg = dvdttmp*dvdtgtmp
            else
              cftg = 0.0d0
            end if

            at(i  ,jet-2,k,nost) = (cff1*atb(i  ,jet-2,k,nost)+cfnm*at(i  ,jet-3,k,nost)-         &
                &       max(cftg,0.0d0)*(atb(i  ,jet-3,k,nost)-atb(i-1,jet-3,k,nost))-              &
                &       min(cftg,0.0d0)*(atb(i+1,jet-3,k,nost)-atb(i  ,jet-3,k,nost)))/max(cff1+cfnm,1.0d-20)

            if (nost .eq. 1) then
              at(i  ,jet-2,k,nost) = at(i  ,jet-2,k,nost) + tau_obc * (( TBCLS(I,K) )-at(i,jet-2,k,nost))
            elseif (nost .eq. 2) then
              at(i  ,jet-2,k,nost) = at(i  ,jet-2,k,nost) + tau_obc * (( SBCLS(I,K) )-at(i,jet-2,k,nost))
            endif

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(I,JET-2,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(I,JET-2,K-1,1) = AT(I,JET-2,K-1,1) - CMIU * ( ATB(I,JET-2,K-1,1) - ATB(I,JET-2,K,1) )
                AT(I,JET-2,K-1,2) = AT(I,JET-2,K-1,2) - CMIU * ( ATB(I,JET-2,K-1,2) - ATB(I,JET-2,K,2) )
              ELSE
                AT(I,JET-2,K,1) = AT(I,JET-2,K,1) - CMIU * ( ATB(I,JET-2,K-1,1) - ATB(I,JET-2,K,1) )
                AT(I,JET-2,K,2) = AT(I,JET-2,K,2) - CMIU * ( ATB(I,JET-2,K-1,2) - ATB(I,JET-2,K,2) )
              ENDIF
            ENDIF

            at(I,JET-1,K,nost) = at(I,JET-2,K,nost) * VIT(I,JET-1,K)
            at(I,JET  ,K,nost) = at(I,JET-2,K,nost) * VIT(I,JET,K)
          enddo

          AT(1,JET-2,K,NOST) = AT(2,JET-2,K,NOST)*VIT(1,JET-2,K)
          AT(1,JET-1,K,NOST) = AT(2,JET-1,K,NOST)*VIT(1,JET-1,K)
          AT(1,JET  ,K,NOST) = AT(2,JET,K,NOST)*VIT(1,JET,K)
          AT(IMT,JET-2,K,NOST) = AT(IMT-1,JET-2,K,NOST)*VIT(IMT,JET-2,K)
          AT(IMT,JET-1,K,NOST) = AT(IMT-1,JET-1,K,NOST)*VIT(IMT,JET-1,K)
          AT(IMT,JET  ,K,NOST) = AT(IMT-1,JET,K,NOST)*VIT(IMT,JET,K)
        enddo
      enddo
    ENDIF
  ELSE IF (JEWSN .EQ. 4 ) THEN
    ! Northern Boundary Condition

    IF ( ID_OBC_TS_N .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For T and S at Northern Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,I)
        DO K = 1, KM
          DO I = 1, IMT
            AT(I,3,K,NOST) = AT(I,4,K,NOST) * VIT(I,3,K)
            AT(I,2,K,NOST) = AT(I,3,K,NOST) * VIT(I,2,K)
            AT(I,1,K,NOST) = AT(I,3,K,NOST) * VIT(I,1,K)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_TS_N .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For T and S at Northern Edge
      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,I)
        DO K = 1, KM
          DO I = 1, IMT
            IF ( NOST.EQ.1 ) THEN
              AT(I,3,K,NOST) = TBCLN(I,K) * VIT(I,3,K)
            ELSE IF ( NOST.EQ.2 ) THEN
              AT(I,3,K,NOST) = SBCLN(I,K) * VIT(I,3,K)
            ENDIF

            AT(I,2,K,NOST) = AT(I,3,K,NOST) * VIT(I,2,K)
            AT(I,1,K,NOST) = AT(I,3,K,NOST) * VIT(I,1,K)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( ID_OBC_TS_N .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For T and S at Northern Edge

      !$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 1, KM
        DO I = 2, IMM
          UDTDX = 0.5D0 * ( V(I,3,K) + V(I+1,3,K) ) * DTS * OUY(3)
          IF ( UDTDX.GT.0.0D0 ) THEN
            AT(I,3,K,1) = ATB(I,3,K,1) - UDTDX * ( ATB(I,3,K,1) - TBCLN(I,K) )
            AT(I,3,K,2) = ATB(I,3,K,2) - UDTDX * ( ATB(I,3,K,2) - SBCLN(I,K) )
          ELSE
            AT(I,3,K,1) = ATB(I,3,K,1) - UDTDX * ( ATB(I,4,K,1) - ATB(I,3,K,1) )
            AT(I,3,K,2) = ATB(I,3,K,2) - UDTDX * ( ATB(I,4,K,2) - ATB(I,3,K,2) )

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(I,3,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(I,3,K-1,1) = AT(I,3,K-1,1) - CMIU * ( ATB(I,3,K-1,1) - ATB(I,3,K,1) )
                AT(I,3,K-1,2) = AT(I,3,K-1,2) - CMIU * ( ATB(I,3,K-1,2) - ATB(I,3,K,2) )
              ELSE
                AT(I,3,K,1) = AT(I,3,K,1) - CMIU * ( ATB(I,3,K-1,1) - ATB(I,3,K,1) )
                AT(I,3,K,2) = AT(I,3,K,2) - CMIU * ( ATB(I,3,K-1,2) - ATB(I,3,K,2) )
              ENDIF
            ENDIF
          ENDIF
          AT(I,2,K,1) = AT(I,3,K,1) * VIT(I,2,K)
          AT(I,1,K,1) = AT(I,3,K,1) * VIT(I,1,K)
          AT(I,2,K,2) = AT(I,3,K,2) * VIT(I,2,K)
          AT(I,1,K,2) = AT(I,3,K,2) * VIT(I,1,K)
        ENDDO
        AT(1,3,K,1) = AT(2,3,K,1) * VIT(1,3,K)
        AT(1,2,K,1) = AT(2,2,K,1) * VIT(1,2,K)
        AT(1,1,K,1) = AT(2,1,K,1) * VIT(1,1,K)
        AT(IMT,3,K,1) = AT(IMT-1,3,K,1) * VIT(IMT,3,K)
        AT(IMT,2,K,1) = AT(IMT-1,2,K,1) * VIT(IMT,2,K)
        AT(IMT,1,K,1) = AT(IMT-1,1,K,1) * VIT(IMT,1,K)
        AT(1,3,K,2) = AT(2,3,K,2) * VIT(1,3,K)
        AT(1,2,K,2) = AT(2,2,K,2) * VIT(1,2,K)
        AT(1,1,K,2) = AT(2,1,K,2) * VIT(1,1,K)
        AT(IMT,3,K,2) = AT(IMT-1,3,K,2) * VIT(IMT,3,K)
        AT(IMT,2,K,2) = AT(IMT-1,2,K,2) * VIT(IMT,2,K)
        AT(IMT,1,K,2) = AT(IMT-1,1,K,2) * VIT(IMT,1,K)
      ENDDO
    ELSE IF ( ID_OBC_TS_N .EQ. 4 ) THEN
      ! The Implicit Upstream Radiation Condition For T and S at Northern Edge

      DO NOST = 1, 2
        !$OMP PARALLEL DO PRIVATE (K,I)
        DO K = 1, KM
          DO I =2, IMT-1
            DVDTTMP = ATB(I,4,K,NOST) - AT(I,4,K,NOST)
            DVDNMTMP = AT(I,4,K,NOST) - AT(I,5,K,NOST)
            DVDTG2TMP = ATB(I+1,4,K,NOST) - ATB(I-1,4,K,NOST)

            IF ( ( DVDTTMP*DVDNMTMP ) .LT. 0.0D0 ) THEN
              DVDTTMP = 0.0D0
              Tau_obc = TAUIN_TS
            ELSE
              Tau_obc = TAUOUT_TS
            ENDIF

            IF( ( DVDTTMP*DVDTG2TMP ) .GT. 0.0) THEN
              DVDTGTMP = ATB(I,4,K,NOST) - ATB(I-1,4,K,NOST)
            ELSE
              DVDTGTMP = ATB(I+1,4,K,NOST) - ATB(I,4,K,NOST)
            ENDIF

            IF ( ID_OBC_2DRad_TS .EQ. 0 ) THEN
              DVDTGTMP = 0.0D0
            ENDIF

            CFF1 = DVDNMTMP * DVDNMTMP + &
                DVDTGTMP * DVDTGTMP * DYR(4) * DYR(4) / &
                ( ( 0.50D0/OUX(3)+0.50D0/OUX(4) ) * ( 0.50D0/OUX(3)+0.50D0/OUX(4) ) )

            CFF2 = DVDTGTMP * DVDTGTMP + &
                DVDNMTMP * DVDNMTMP / ( DYR(4)*DYR(4) ) * &
                ( 0.50D0/OUX(3)+0.50D0/OUX(4) ) * ( 0.50D0/OUX(3)+0.50D0/OUX(4) )

            CFF1 = MAX(CFF1,EPSCF)
            CFF2 = MAX(CFF2,EPSCF)
            CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
            CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

            IF ( ID_OBC_2DRad_NPO_TS .EQ. 1 ) THEN
              CFTG = 0.0D0
            ENDIF

            CMIU = CFF1 + CFNM
            AT(I,3,K,NOST) = ( CFF1*ATB(I,3,K,NOST) + &
                CFNM*AT(I,4,K,NOST) -  &
                MAX(CFTG,0.0D0)*(CFF1/CFF2)*(ATB(I,3,K,NOST)-ATB(I-1,3,K,NOST)) - &
                MIN(CFTG,0.0D0)*(CFF1/CFF2)*(ATB(I+1,3,K,NOST)-ATB(I,3,K,NOST)) )/CMIU

            IF (NOST.EQ.1) THEN
              AT(I,3,K,NOST) = AT(I,3,K,NOST) + &
                  DTS/Tau_obc*CFF1*(TBCLN(I,K)-ATB(I,3,K,NOST))/CMIU
            ELSE IF (NOST.EQ.2) THEN
              AT(I,3,K,NOST) = AT(I,3,K,NOST) + &
                  DTS/Tau_obc*CFF1*(SBCLN(I,K)-ATB(I,3,K,NOST))/CMIU
            ENDIF

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(I,3,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(I,3,K-1,1) = AT(I,3,K-1,1) - CMIU * ( ATB(I,3,K-1,1) - ATB(I,3,K,1) )
                AT(I,3,K-1,2) = AT(I,3,K-1,2) - CMIU * ( ATB(I,3,K-1,2) - ATB(I,3,K,2) )
              ELSE
                AT(I,3,K,1) = AT(I,3,K,1) - CMIU * ( ATB(I,3,K-1,1) - ATB(I,3,K,1) )
                AT(I,3,K,2) = AT(I,3,K,2) - CMIU * ( ATB(I,3,K-1,2) - ATB(I,3,K,2) )
              ENDIF
            ENDIF

            AT(I,2,K,NOST) = AT(I,3,K,NOST) * VIT(I,2,K)
            AT(I,1,K,NOST) = AT(I,3,K,NOST) * VIT(I,1,K)
          ENDDO
          AT(1,3,K,NOST) = 0.5*AT(2,3,K,NOST) * VIT(1,3,K) + 0.5*ATB(1,3,K,NOST)
          AT(1,2,K,NOST) = 0.5*AT(2,2,K,NOST) * VIT(1,2,K) + 0.5*ATB(1,2,K,NOST)
          AT(1,1,K,NOST) = 0.5*AT(2,1,K,NOST) * VIT(1,1,K) + 0.5*ATB(1,1,K,NOST)
          AT(IMT,3,K,NOST) = 0.5*AT(IMT-1,3,K,NOST) * VIT(IMT,3,K) + 0.5*ATB(IMT,3,K,NOST)
          AT(IMT,2,K,NOST) = 0.5*AT(IMT-1,2,K,NOST) * VIT(IMT,2,K) + 0.5*ATB(IMT,2,K,NOST)
          AT(IMT,1,K,NOST) = 0.5*AT(IMT-1,1,K,NOST) * VIT(IMT,1,K) + 0.5*ATB(IMT,1,K,NOST)
        ENDDO
      ENDDO
    else if ( ID_OBC_TS_N .eq. 5 ) then
      do nost = 1, 2
        do k = 1, km
          do i = 2, imt-1
            dvdttmp  = atb(i,4,k,nost)-at(i,4,k,nost)
            dvdnmtmp = at (i,4,k,nost)-at(i,5,k,nost)

            if (dvdttmp*dvdnmtmp .lt. 0.0d0) then
              tau_obc = tau_in_c5
              dvdttmp = 0.0d0
            else
              tau_obc = tau_out_c5
            endif
            tau_obc = tau_obc*dtc_openbc

            if (dvdttmp*(atb(i+1,4,k,nost)-atb(i-1,4,k,nost)) .gt. 0.0d0) then
              dvdtgtmp = atb(i  ,4,k,nost) - atb(i-1,4,k,nost)
            else
              dvdtgtmp = atb(i+1,4,k,nost) - atb(i  ,4,k,nost)
            endif

            cff1 = max(dvdnmtmp*dvdnmtmp+dvdtgtmp*dvdtgtmp,3.0d-7)
            cfnm = dvdttmp*dvdnmtmp
            if (ID_OBC_2DRad_NPO_TS5 == 0 ) then
              cftg = dvdttmp*dvdtgtmp
            else
              cftg = 0.0d0
            end if

            at(i,3,k,nost) = (cff1*atb(i,3,k,nost)+cfnm*at(i,4,k,nost)-                   &
                &       max(cftg,0.0d0)*(atb(i  ,3,k,nost)-atb(i-1,3,k,nost))-              &
                &       min(cftg,0.0d0)*(atb(i+1,3,k,nost)-atb(i  ,3,k,nost)))/max(cff1+cfnm,1.0d-20)

            if (nost .eq. 1) then
              at(i,3,k,nost) = at(i,3,k,nost) + tau_obc * (( TBCLN(I,K) )-at(i,1,k,nost))
            elseif (nost .eq. 2) then
              at(i,3,k,nost) = at(i,3,k,nost) + tau_obc * (( SBCLN(I,K) )-at(i,1,k,nost))
            endif

            IF ( K.GE.2 .AND. K.LE.KM ) THEN
              CMIU = WS(I,3,K) * DTS * ODZT(K-1)
              IF ( CMIU .GE. 0.0 ) THEN
                AT(I,3,K-1,1) = AT(I,3,K-1,1) - CMIU * ( ATB(I,3,K-1,1) - ATB(I,3,K,1) )
                AT(I,3,K-1,2) = AT(I,3,K-1,2) - CMIU * ( ATB(I,3,K-1,2) - ATB(I,3,K,2) )
              ELSE
                AT(I,3,K,1) = AT(I,3,K,1) - CMIU * ( ATB(I,3,K-1,1) - ATB(I,3,K,1) )
                AT(I,3,K,2) = AT(I,3,K,2) - CMIU * ( ATB(I,3,K-1,2) - ATB(I,3,K,2) )
              ENDIF
            ENDIF

            at(I,2,K,nost) = at(I,3,K,nost) * VIT(I,2,K)
            at(I,1,K,nost) = at(I,3,K,nost) * VIT(I,1,K)
          enddo

          AT(1,3,K,NOST) = AT(2,3,K,NOST) * VIT(1,3,K)
          AT(1,2,K,NOST) = AT(2,2,K,NOST) * VIT(1,2,K)
          AT(1,1,K,NOST) = AT(2,1,K,NOST) * VIT(1,1,K)
          AT(IMT,3,K,NOST) = AT(IMT-1,3,K,NOST) * VIT(IMT,3,K)
          AT(IMT,2,K,NOST) = AT(IMT-1,2,K,NOST) * VIT(IMT,2,K)
          AT(IMT,1,K,NOST) = AT(IMT-1,1,K,NOST) * VIT(IMT,1,K)
        enddo
      enddo
    ENDIF
  ENDIF

  !$OMP PARALLEL DO PRIVATE (K,J,I)
  DO K = 1, KM
    DO J = 1, JMT
      DO I = 1, IMT
        AT(I,J,K,1) = AT(I,J,K,1) * VIT(I,J,K)
        AT(I,J,K,2) = AT(I,J,K,2) * VIT(I,J,K)
      ENDDO
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE OBC_3DTS

SUBROUTINE OBC_2DUB(JEWSN)

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, KM, JST, JET, JMM, I, J, K
  use   pconst_mod, only: VIT, VIV, DTB, DTB2, DYT, DYR, OTX, OHBU, DZP, G, OUX, OUY, OHBT, iy
  use   dyn_mod,    only: UB, VB, H0, H0P, UBP, VBP
  use   work_mod

#if ( defined SPMD )
  use msg_mod
#endif

#ifdef OBCDT
  use   openbc_mod
#endif

  IMPLICIT NONE
  INTEGER :: JEWSN
  REAL(r8) :: EPSCF
  REAL(r8) :: Tau_obc
  REAL(r8) :: UDTDX
  REAL(r8) :: DVDTTMP, DVDNMTMP, DVDTGTMP, DVDTG2TMP, CFTG, CFNM
  REAL(r8) :: CMIU
  REAL(r8) :: CFF, CFF1, CFF2, CFF3
  REAL(r8) :: Chalf, COUR
  REAL(r8) :: ZSTAR
  REAL(r8) :: HTMP
  EPSCF = 1.0D-20

  IF ( JEWSN .EQ. 1 ) THEN
    ! Eastern Boundary Condition
    IF ( ID_OBC_UB_E .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For UB at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        UB(IMT-2,J) = UB(IMT-3,J) * VIV(IMT-2,J,1)
        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT,J) = UB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_E .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For UB at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        UB(IMT-2,J) = UBTRE(2,J)*VIV(IMT-2,J,1)
        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT,J) = UB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_E .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For UB at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        UDTDX = UB(IMT-2,J) * DTB * OTX(J)

        IF ( UDTDX.LT.0.0D0 ) THEN
          UB(IMT-2,J) = UBP(IMT-2,J) - UDTDX * ( UBTRE(2,J) - UBP(IMT-2,J) )
        ELSE
          UB(IMT-2,J) = UBP(IMT-2,J) - UDTDX * ( UBP(IMT-2,J) - UBP(IMT-3,J) )
        ENDIF

        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT,J) = UB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_E .EQ. 4 ) THEN
      ! The Implicit Chapman Boundary Condition For UB at Eastern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        IF ( OHBU(IMT-3,J).GT.0.0D0 .AND. OHBU(IMT-3,J).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(IMT-3,J) + &
              0.25D0*( H0P(IMT-4,J) * VIT(IMT-4,J,1) + H0P(IMT-4,J+1) * VIT(IMT-4,J+1,1)  + &
              H0P(IMT-3,J) * VIT(IMT-3,J,1) + H0P(IMT-3,J+1) * VIT(IMT-3,J+1,1) )

          CMIU = SQRT( G*HTMP )*DTB / (0.5D0/OTX(J)+0.5D0/OTX(J+1))
          UB(IMT-2,J) = ( 0.25D0*UBP(IMT-2,J-1) + 0.50D0*UBP(IMT-2,J) + &
              0.25D0*UBP(IMT-2,J+1) ) / ( 1.0D0+CMIU ) + &
              ( 0.25D0*UB(IMT-3,J-1)  + 0.50D0*UB(IMT-3,J)  + &
              0.25D0*UB(IMT-3,J+1)  ) * CMIU / ( 1.0D0+CMIU )
        ELSE
          UB(IMT-2,J) = 0.0D0
        ENDIF

        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT,J)   = UB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO

      UB(IMT-2,1) = UB(IMT-2,2) * VIV(IMT-2,1,1)
      UB(IMT-1,1) = UB(IMT-1,2) * VIV(IMT-1,1,1)
      UB(IMT,1) = UB(IMT,2) * VIV(IMT,1,1)
      UB(IMT-2,JMT) = UB(IMT-2,JMT-1) * VIV(IMT-2,JMT,1)
      UB(IMT-1,JMT) = UB(IMT-1,JMT-1) * VIV(IMT-1,JMT,1)
      UB(IMT,JMT) = UB(IMT,JMT-1) * VIV(IMT,JMT,1)
    ELSE IF ( ID_OBC_UB_E .EQ. 5 ) THEN
      ! The Explicit Chapman Boundary Condition For UB at Eastern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        IF ( OHBU(IMT-3,J).GT.0.0D0 .AND. OHBU(IMT-3,J).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(IMT-3,J) + &
              0.25D0*( H0P(IMT-4,J) * VIT(IMT-4,J,1) + H0P(IMT-4,J+1) * VIT(IMT-4,J+1,1)  + &
              H0P(IMT-3,J) * VIT(IMT-3,J,1) + H0P(IMT-3,J+1) * VIT(IMT-3,J+1,1) )

          CMIU = SQRT( G*HTMP )*DTB / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          CMIU = MIN( CMIU, 0.9999 )
          UB(IMT-2,J) = ( 0.25D0*UBP(IMT-2,J-1) + 0.50D0*UBP(IMT-2,J) + &
              0.25D0*UBP(IMT-2,J+1) ) * ( 1.0D0-CMIU ) + &
              ( 0.25D0*UBP(IMT-3,J-1) + 0.50D0*UBP(IMT-3,J) + &
              0.25D0*UBP(IMT-3,J+1) ) * CMIU
        ELSE
          UB(IMT-2,J) = 0.0D0
        ENDIF

        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT,J)   = UB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO

      UB(IMT-2,1) = UB(IMT-2,2) * VIV(IMT-2,1,1)
      UB(IMT-1,1) = UB(IMT-1,2) * VIV(IMT-1,1,1)
      UB(IMT,1) = UB(IMT,2) * VIV(IMT,1,1)
      UB(IMT-2,JMT) = UB(IMT-2,JMT-1) * VIV(IMT-2,JMT,1)
      UB(IMT-1,JMT) = UB(IMT-1,JMT-1) * VIV(IMT-1,JMT,1)
      UB(IMT,JMT) = UB(IMT,JMT-1) * VIV(IMT,JMT,1)
    ELSE IF ( ID_OBC_UB_E .EQ. 6 ) THEN
      ! The Implicit Upstream Radiation Condition For UB at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        DVDTTMP = UBP(IMT-3,J) - UB(IMT-3,J)
        DVDNMTMP = UB(IMT-3,J) - UB(IMT-4,J)
        DVDTG2TMP = UBP(IMT-3,J+1) - UBP(IMT-3,J-1)

        IF ( (DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
          DVDTTMP = 0.0D0
          Tau_obc = TAUIN
        ELSE
          Tau_obc = TAUOUT
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          DVDTGTMP = UBP(IMT-3,J) - UBP(IMT-3,J-1)
        ELSE
          DVDTGTMP = UBP(IMT-3,J+1) - UBP(IMT-3,J)
        ENDIF

        IF ( ID_OBC_2DRad_UB .EQ. 0 ) THEN
          DVDTGTMP = 0.0D0
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP / ( DYT(J)*DYT(J) ) * &
              ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP * DYT(J) * DYT(J) / &
              ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
        ELSE
          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP / ( DYT(J+1)*DYT(J+1) ) * &
              ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP * DYT(J+1) * DYT(J+1) / &
              ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
        ENDIF

        CFF1 = MAX(CFF1,EPSCF)
        CFF2 = MAX(CFF2,EPSCF)
        CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
        CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

        IF ( ID_OBC_2DRad_NPO_UB .EQ. 1 ) THEN
          CFTG = 0.0D0
        ENDIF

        CMIU = CFF1 + CFNM
        UB(IMT-2,J) = ( CFF1*UBP(IMT-2,J) + &
            CFNM*UB(IMT-3,J) -  &
            MAX(CFTG,0.0D0)*(CFF1/CFF2)*(UBP(IMT-2,J)-UBP(IMT-2,J-1))- &
            MIN(CFTG,0.0D0)*(CFF1/CFF2)*(UBP(IMT-2,J+1)-UBP(IMT-2,J)) )/CMIU

        UB(IMT-2,J) = UB(IMT-2,J) + DTB/Tau_obc * ( (0.2*UBTRE(1,J)+0.4*UBTRE(2,J)+ &
            0.2*UBTRE(2,J-1)+0.2*UBTRE(2,J+1))-UBP(IMT-2,J) )*CFF1/CMIU

        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT,J) = UB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
      UB(IMT-2,1) = UB(IMT-2,2) * VIV(IMT-2,1,1)
      UB(IMT-1,1) = UB(IMT-1,2) * VIV(IMT-1,1,1)
      UB(IMT,1) = UB(IMT,2) * VIV(IMT,1,1)
      UB(IMT-2,JMT) = UB(IMT-2,JMT-1) * VIV(IMT-2,JMT,1)
      UB(IMT-1,JMT) = UB(IMT-1,JMT-1) * VIV(IMT-1,JMT,1)
      UB(IMT,JMT) = UB(IMT,JMT-1) * VIV(IMT,JMT,1)
    ELSE IF ( ID_OBC_UB_E .EQ. 7 ) THEN
      ! Flather Boundary Condition For UB at Eastern Edge
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO
      CMIU = 1.0d-10

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMM
        IF ( OHBU(IMT-2,J).GT.0.0D0 .AND. OHBU(IMT-2,J).LE.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(IMT-2,J) + &
              0.25D0 * ( H0PP(IMT-3,J) * VIT(IMT-3,J,1) + H0PP(IMT-3,J+1) * VIT(IMT-3,J+1,1) +  &
              H0PP(IMT-2,J) * VIT(IMT-2,J,1) + H0PP(IMT-2,J+1) * VIT(IMT-2,J+1,1) )

          UB(IMT-2,J) = UBTRE(2,J) + SQRT( G/HTMP ) *       &
              ( 0.25D0 * ( H0PP(IMT-3,J) * VIT(IMT-3,J,1) + H0PP(IMT-3,J+1) * VIT(IMT-3,J+1,1) +   &
              H0PP(IMT-2,J) * VIT(IMT-2,J,1) + H0PP(IMT-2,J+1) * VIT(IMT-2,J+1,1) ) - &
              0.25D0 * ( H0BTRE(1,J) + H0BTRE(1,J+1) + H0BTRE(2,J) + H0BTRE(2,J+1) ) )
        ELSE
          UB(IMT-2,J ) = 0.0D0
        ENDIF

        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT,J) = UB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
      UB(IMT-2,JMT) = UB(IMT-2,JMM) * VIT(IMT-2,JMT,1)
      UB(IMT-1,JMT) = UB(IMT-1,JMM) * VIT(IMT-1,JMT,1)
      UB(IMT,JMT) = UB(IMT,JMM) * VIT(IMT,JMT,1)
    ELSE IF ( ID_OBC_UB_E .EQ. 8 ) THEN
      ! The Modified Flather Boundary Condition For UB at Eastern Edge by
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO
      CMIU = 1.0d-10

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMM
        IF ( OHBU(IMT-2,J).GT.0.0D0 .AND. OHBU(IMT-2,J).LE.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(IMT-2,J)
          CFF1 = SQRT(G/HTMP)
          COUR = DTB_OPENBC*CFF1*HTMP/( 0.5D0*(1.0D0/OTX(J)+1.0D0/OTX(J+1)) )
          ZSTAR = ( 0.5D0 + COUR ) * &
              0.5D0 * ( H0Pp(IMT-3,J)*VIT(IMT-3,J,1) + H0Pp(IMT-3,J+1)*VIT(IMT-3,J+1,1) ) + &
              ( 0.5D0 - COUR ) * &
              0.5D0 * ( H0Pp(IMT-2,J)*VIT(IMT-2,J,1) + H0Pp(IMT-2,J+1)*VIT(IMT-2,J+1,1) )

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0 - Chalf/COUR )**2
            CFF3 = 0.5D0 * ( H0(IMT-3,J)*VIT(IMT-3,J,1) + H0(IMT-3,J+1)*VIT(IMT-3,J+1,1) ) + &
                COUR * &
                0.5D0 * ( H0Pp(IMT-2,J)*VIT(IMT-2,J,1) + H0Pp(IMT-2,J+1)*VIT(IMT-2,J+1,1) ) - &
                ( 1.0D0 + COUR ) * &
                0.5D0 * ( H0Pp(IMT-3,J)*VIT(IMT-3,J,1) + H0Pp(IMT-3,J+1)*VIT(IMT-3,J+1,1) )
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF
          UB(IMT-2,J) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * UBPp(IMT-2,J) + COUR * UBPp(IMT-3,J) + &
              UBTRE(2,J) + &
              CFF1 * ( ZSTAR - 0.5D0 * ( H0BTRE(2,J) + H0BTRE(2,J+1) ) ) )
        ELSE
          UB(IMT-2,J ) = 0.0D0
        ENDIF

        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT,J) = UB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
      UB(IMT-2,JMT) = UB(IMT-2,JMM) * VIT(IMT-2,JMT,1)
      UB(IMT-1,JMT) = UB(IMT-1,JMM) * VIT(IMT-1,JMT,1)
      UB(IMT,JMT) = UB(IMT,JMM) * VIT(IMT,JMT,1)
    ELSE IF ( ID_OBC_UB_E .EQ. 9 ) THEN
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT-1
        IF (OHBU(IMT-3,J) .GT. 1.0D-20) THEN
          HTMP = 1.0D0/OHBU(IMT-3,J)
          CFF1 = SQRT(G/HTMP)
          COUR = DTB_OPENBC*CFF1*HTMP*OUX(J)
          ZSTAR = ( 0.5D0 + COUR ) * (H0PP(IMT-3,J+1) + H0PP(IMT-3,J))*0.5D0 +     &
              &       ( 0.5D0 - COUR ) * (H0PP(IMT-2,J+1) + H0PP(IMT-2,J))*0.5D0

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0 - Chalf/COUR )**2
            CFF3 = 0.50D0*(H0(IMT-3,J+1) + H0(IMT-3,J)) + &
                COUR * &
                0.50D0*(H0PP(IMT-2,J+1) + H0PP(IMT-2,J)) - &
                ( 1.0D0 + COUR ) * &
                0.50D0*(H0PP(IMT-3,J+1) + H0PP(IMT-3,J))
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF
          UB(IMT-2,J) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * UBPP(IMT-2,J) + COUR * UBPP(IMT-3,J) + &
              UBTRE(2,J) + &
              CFF1 * ( ZSTAR - 0.50D0*( H0BTRE(2,J) + H0BTRE(2,J+1)) ) )
        ENDIF
        UB(IMT-1,J) = UB(IMT-2,J) * VIV(IMT-1,J,1)
        UB(IMT  ,J) = UB(IMT-2,J) * VIV(IMT  ,J,1)
      ENDDO
      UB(IMT-1,JMT) = UB(IMT-1,JMT-1) * VIT(IMT-1,JMT,1)
      UB(IMT  ,JMT) = UB(IMT  ,JMT-1) * VIT(IMT  ,JMT,1)
    ELSE IF ( ID_OBC_UB_E .EQ. 10 ) THEN
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        HTMP = 1.0D0/OHBU(IMT-2,J) + 0.5D0 * (H0Pp(IMT-3,J) + H0Pp(IMT-2,J))
        IF (HTMP .LT. 1.0D+5 .AND. HTMP .GT. 1.0D-5) THEN
          UB(IMT-2,J) = UBTRE(2,J) + SQRT( G/HTMP ) *       &
              ( 0.50D0 * ( H0Pp(IMT-3,J) + H0Pp(IMT-2,J) ) -  H0BTRE(1,J))
        ELSE
          UB(IMT-2,J ) = 0
        ENDIF

        UB(IMT-2,J) = UB(IMT-2,J) * VIV(IMT-2,J,1)
        UB(IMT-1,J) = UB(IMT-2,J)
        UB(IMT,J) = UB(IMT-2,J)
      ENDDO
    ELSE IF ( ID_OBC_UB_E .EQ. 11 ) THEN
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT-1
        HTMP = 1.0D0/OHBU(IMT-2,J)

        IF (HTMP .LT. 1.0D+5 .AND. HTMP .GT. 1.0D0) THEN
          CFF1 = SQRT(G/HTMP)
          COUR = DTB2*CFF1*HTMP*OUX(J)
          ZSTAR = ( 0.5D0 + COUR ) * (H0P(IMT-3,J+1) + H0P(IMT-3,J))*0.5D0 +     &
              &       ( 0.5D0 - COUR ) * (H0P(IMT-2,J+1) + H0P(IMT-2,J))*0.5D0

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0 - Chalf/COUR )**2
            CFF3 = 0.50D0*(WORK(IMT-3,J+1) + WORK(IMT-3,J)) + &
                COUR * &
                0.50D0*(H0P(IMT-2,J+1) + H0P(IMT-2,J)) - &
                ( 1.0D0 + COUR ) * &
                0.50D0*(H0P(IMT-3,J+1) + H0P(IMT-3,J))
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF
          WKA(IMT-2,J,1) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * UBP(IMT-2,J) + COUR * UBP(IMT-3,J) + &
              UBTRE(2,J) + &
              CFF1 * ( ZSTAR - 0.50D0*( H0BTRE(2,J) + H0BTRE(2,J+1)) ) )
        ELSE
          WKA(IMT-2,J,1) = WKA(IMT-3,J,1)
        ENDIF

        WKA(IMT-1,J,1) = WKA(IMT-2,J,1) * VIV(IMT-1,J,1)
        WKA(IMT,J,1) = WKA(IMT-2,J,1) * VIV(IMT,J,1)
      ENDDO
      WKA(IMT-2,JMT,1) = WKA(IMT-2,JMT-1,1) * VIT(IMT-2,JMT,1)
      WKA(IMT-1,JMT,1) = WKA(IMT-1,JMT-1,1) * VIT(IMT-1,JMT,1)
      WKA(IMT  ,JMT,1) = WKA(IMT  ,JMT-1,1) * VIT(IMT  ,JMT,1)
    ENDIF
  ELSE IF ( JEWSN .EQ. 2 ) THEN
    ! Western Boundary Condition
    IF ( ID_OBC_UB_W .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For UB at Western Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        UB(4,J) = UB(5,J) * VIV(4,J,1)
        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_W .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For UB at Western Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        UB(4,J) = UBTRW(1,J) * VIV(4,J,1)
        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_W .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For UB at Western Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        UDTDX = UB(4,J) * DTB * OTX(J)

        IF ( UDTDX.LT.0.0D0 ) THEN
          UB(4,J) = UBP(4,J) - UDTDX * ( UBTRW(1,J) - UBP(4,J) )
        ELSE
          UB(4,J) = UBP(4,J) - UDTDX * ( UBP(4,J) - UBP(5,J) )
        ENDIF

        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_W .EQ. 4 ) THEN
      ! The Implicit Chapman Boundary Condition For UB at Western Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        IF ( OHBU(5,J).GT.0.0D0 .AND. OHBU(5,J).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(5,J) + &
              0.25D0*( H0P(4,J) * VIT(4,J,1) + H0P(4,J+1) * VIT(4,J+1,1) + &
              H0P(5,J) * VIT(5,J,1) + H0P(5,J+1) * VIT(5,J+1,1) )

          CMIU = SQRT( G*HTMP )*DTB / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          UB(4,J) = ( 0.25D0*UBP(4,J-1) + 0.50D0*UBP(4,J) + 0.25D0*UBP(4,J+1) ) / ( 1.0D0+CMIU ) + &
              ( 0.25D0*UB(5,J-1)  + 0.50D0*UB(5,J)  + 0.25D0*UB(5,J+1)  ) * CMIU / ( 1.0D0+CMIU )
        ELSE
          UB(4,J) = 0.0D0
        ENDIF

        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO

      UB(4,1) = UB(4,2) * VIV(4,1,1)
      UB(3,1) = UB(3,2) * VIV(3,1,1)
      UB(2,1) = UB(2,2) * VIV(2,1,1)
      UB(1,1) = UB(1,2) * VIV(1,1,1)
      UB(4,JMT) = UB(4,JMT-1) * VIV(4,JMT,1)
      UB(3,JMT) = UB(3,JMT-1) * VIV(3,JMT,1)
      UB(2,JMT) = UB(2,JMT-1) * VIV(2,JMT,1)
      UB(1,JMT) = UB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_UB_W .EQ. 5 ) THEN
      ! The Explicit Chapman Boundary Condition For UB at Western Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        IF ( OHBU(5,J).GT.0.0D0 .AND. OHBU(5,J).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(5,J) + &
              0.25D0*( H0P(4,J) * VIT(4,J,1) + H0P(4,J+1) * VIT(4,J+1,1) + &
              H0P(5,J) * VIT(5,J,1) + H0P(5,J+1) * VIT(5,J+1,1) )

          CMIU = SQRT( G*HTMP )*DTB / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          CMIU = MIN( CMIU, 0.9999 )
          UB(4,J) = ( 0.25D0*UBP(4,J-1) + 0.50D0*UBP(4,J) + 0.25D0*UBP(4,J+1) ) * ( 1.0D0 - CMIU ) +  &
              ( 0.25D0*UBP(5,J-1) + 0.50D0*UBP(5,J) + 0.25D0*UBP(5,J+1) ) * CMIU
        ELSE
          UB(4,J) = 0.0D0
        ENDIF

        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO

      UB(4,1) = UB(4,2) * VIV(4,1,1)
      UB(3,1) = UB(3,2) * VIV(3,1,1)
      UB(2,1) = UB(2,2) * VIV(2,1,1)
      UB(1,1) = UB(1,2) * VIV(1,1,1)
      UB(4,JMT) = UB(4,JMT-1) * VIV(4,JMT,1)
      UB(3,JMT) = UB(3,JMT-1) * VIV(3,JMT,1)
      UB(2,JMT) = UB(2,JMT-1) * VIV(2,JMT,1)
      UB(1,JMT) = UB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_UB_W .EQ. 6 ) THEN
      ! The Implicit Upstream Radiation Condition For UB at Western Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        DVDTTMP = UBP(5,J) - UB(5,J)
        DVDNMTMP = UB(5,J) - UB(6,J)
        DVDTG2TMP = UBP(5,J+1) - UBP(5,J-1)

        IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
          DVDTTMP = 0.0D0
          Tau_obc = TAUIN
        ELSE
          Tau_obc = TAUOUT
        ENDIF

        IF ( ( DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          DVDTGTMP = UBP(5,J) - UBP(5,J-1)
        ELSE
          DVDTGTMP = UBP(5,J+1) - UBP(5,J)
        ENDIF

        IF ( ID_OBC_2DRad_UB .EQ. 0 ) THEN
          DVDTGTMP = 0.0D0
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP / ( DYT(J)*DYT(J) ) * &
              ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP * DYT(J) * DYT(J) / &
              ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
        ELSE
          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP / ( DYT(J+1)*DYT(J+1) ) * &
              ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP * DYT(J+1) * DYT(J+1) / &
              ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
        ENDIF

        CFF1 = MAX(CFF1,EPSCF)
        CFF2 = MAX(CFF2,EPSCF)
        CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
        CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

        IF ( ID_OBC_2DRad_NPO_UB .EQ. 1 ) THEN
          CFTG = 0.0D0
        ENDIF

        CMIU = CFF1 + CFNM
        UB(4,J) = ( CFF1*UBP(4,J) + &
            CFNM*UB(5,J)  - &
            MAX(CFTG,0.0D0)*(CFF1/CFF2)*(UBP(4,J)-UBP(4,J-1)) - &
            MIN(CFTG,0.0D0)*(CFF1/CFF2)*(UBP(4,J+1)-UBP(4,J)) ) / CMIU

        UB(4,J) = UB(4,J) + DTB/Tau_obc * ( (UBTRW(1,J-1)*0.2+UBTRW(1,J+1)*0.2+ &
            UBTRW(1,J)*0.4+UBTRW(2,J)*0.2)-UBP(4,J) ) * CFF1/CMIU

        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO
      UB(4,1) = UB(4,2) * VIV(4,1,1)
      UB(3,1) = UB(3,2) * VIV(3,1,1)
      UB(2,1) = UB(2,2) * VIV(2,1,1)
      UB(1,1) = UB(1,2) * VIV(1,1,1)
      UB(4,JMT) = UB(4,JMT-1) * VIV(4,JMT,1)
      UB(3,JMT) = UB(3,JMT-1) * VIV(3,JMT,1)
      UB(2,JMT) = UB(2,JMT-1) * VIV(2,JMT,1)
      UB(1,JMT) = UB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_UB_W .EQ. 7 ) THEN
      ! Flather Boundary Condition For UB at Western Edge
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO
      CMIU = 1.0d-10

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMM
        IF ( OHBU(4,J).GT.0.0D0 .AND. OHBU(4,J).LE.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(4,J) + &
              0.25D0 * ( H0PP(3,J) * VIT(3,J,1) + H0PP(3,J+1) * VIT(3,J+1,1) +  &
              H0PP(4,J) * VIT(4,J,1) + H0PP(4,J+1) * VIT(4,J+1,1) )

          UB(4,J) = UBTRW(1,J) + SQRT( G/HTMP ) *          &
              ( 0.25D0 * ( H0PP(3,J) * VIT(3,J,1) + H0PP(3,J+1) * VIT(3,J+1,1) +   &
              H0PP(4,J) * VIT(4,J,1) + H0PP(4,J+1) * VIT(4,J+1,1) ) - &
              0.25D0 * ( H0BTRW(1,J) + H0BTRW(1,J+1) + H0BTRW(2,J) + H0BTRW(2,J+1) ) )
        ELSE
          UB(4,J) = 0.0D0
        ENDIF
        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO
      UB(4,JMT) = UB(4,JMT-1) * VIV(4,JMT,1)
      UB(3,JMT) = UB(3,JMT-1) * VIV(3,JMT,1)
      UB(2,JMT) = UB(2,JMT-1) * VIV(2,JMT,1)
      UB(1,JMT) = UB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_UB_W .EQ. 8 ) THEN
      ! The Modified Flather Boundary Condition For UB at Western Edge
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO
      CMIU = 1.0d-10

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMM
        IF ( OHBU(4,J).GT.0.0D0 .AND. OHBU(4,J).LE.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(4,J)
          CFF1 = SQRT(G/HTMP)
          COUR = DTB_OPENBC*CFF1*HTMP / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          ZSTAR = ( 0.5D0 + COUR ) * &
              0.5D0 * ( H0Pp(4,J)*VIT(4,J,1) + H0Pp(4,J+1)*VIT(4,J+1,1) ) + &
              ( 0.5D0 - COUR ) * &
              0.5D0 * ( H0Pp(3,J)*VIT(3,J,1) + H0Pp(3,J+1)*VIT(3,J+1,1) )

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0 - Chalf/COUR )**2
            CFF3 = 0.5D0 * ( H0(4,J)*VIT(4,J,1) + H0(4,J+1)*VIT(4,J+1,1) ) + &
                COUR * &
                0.5D0 * ( H0Pp(3,J)*VIT(3,J,1) + H0Pp(3,J+1)*VIT(3,J+1,1) ) - &
                ( 1.0D0 + COUR ) * &
                0.5D0 * ( H0Pp(4,J)*VIT(4,J,1) + H0Pp(4,J+1)*VIT(4,J+1,1) )
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF
          UB(4,J) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * UBPp(4,J) + COUR * UBPp(5,J) + &
              UBTRW(1,J) - &
              CFF1 * ( ZSTAR - 0.5D0 * ( H0BTRW(1,J) + H0BTRW(1,J+1) ) ) )
        ELSE
          UB(4,J) = 0.0D0
        ENDIF
        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO
      UB(4,JMT) = UB(4,JMT-1) * VIV(4,JMT,1)
      UB(3,JMT) = UB(3,JMT-1) * VIV(3,JMT,1)
      UB(2,JMT) = UB(2,JMT-1) * VIV(2,JMT,1)
      UB(1,JMT) = UB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_UB_W .EQ. 9 ) THEN
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT-1
        IF (OHBU(5,J) .GT. 1.0D-20) THEN
          HTMP = 1.0D0/OHBU(5,J)
          CFF1 = SQRT(G/HTMP)
          COUR = DTB_OPENBC*CFF1*HTMP*OUX(J)
          ZSTAR = ( 0.5D0 + COUR ) * (H0PP(4,J+1) + H0PP(4,J))*0.5D0 +       &
              ( 0.5D0 - COUR ) * (H0PP(3,J+1) + H0PP(3,J))*0.5D0

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0 - Chalf/COUR )**2
            CFF3 = 0.5D0*(H0(4,J+1) + H0(4,J)) + &
                COUR * &
                0.5D0*(H0PP(3,J+1) + H0PP(3,J)) - &
                ( 1.0D0 + COUR ) * &
                0.50D0*(H0PP(4,J+1) + H0PP(4,J))
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF
          UB(4,J) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * UBPP(4,J) + COUR * UBPP(5,J) + &
              UBTRW(1,J) - &
              CFF1 * ( ZSTAR - 0.50D0*(H0BTRW(1,J+1) + H0BTRW(1,J))) )
        ENDIF
        UB(3,J) = UB(4,J) * VIV(3,J,1)
        UB(2,J) = UB(4,J) * VIV(2,J,1)
        UB(1,J) = UB(4,J) * VIV(1,J,1)
      ENDDO
      UB(3,JMT) = UB(3,JMT-1) * VIV(3,JMT,1)
      UB(2,JMT) = UB(2,JMT-1) * VIV(2,JMT,1)
      UB(1,JMT) = UB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_UB_W .EQ. 10 ) THEN
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        HTMP = 1.0D0/OHBU(4,J) + 0.50D0 * (H0Pp(3,J) + H0Pp(4,J))
        IF (HTMP .LT. 1.0D+5 .AND. HTMP .GT. 1.0D-5) THEN
          UB(4,J) = UBTRW(1,J) - SQRT( G/HTMP ) *          &
              ( 0.50D0 * ( H0Pp(3,J)  + H0Pp(4,J)) -  H0BTRW(2,J))
        ELSE
          UB(4,J) = 0
        ENDIF
        UB(4,J) = UB(4,J) * VIV(4,J,1)
        UB(3,J) = UB(4,J)
        UB(2,J) = UB(4,J)
        UB(1,J) = UB(4,J)
      ENDDO
    ELSE IF ( ID_OBC_UB_W .EQ. 11 ) THEN
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT-1
        HTMP = 1.0D0/OHBU(4,J)

        IF (HTMP .LT. 1.0D+5 .AND. HTMP .GT. 1.0D0) THEN
          CFF1 = SQRT(G/HTMP)
          COUR = DTB2*CFF1*HTMP*OUX(J)
          ZSTAR = ( 0.5D0 + COUR ) * (H0P(4,J+1) + H0P(4,J))*0.5D0 +       &
              ( 0.5D0 - COUR ) * (H0P(3,J+1) + H0P(3,J))*0.5D0

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0 - Chalf/COUR )**2
            CFF3 = 0.5D0*(WORK(4,J+1) + WORK(4,J)) + &
                COUR * &
                0.5D0*(H0P(3,J+1) + H0P(3,J)) - &
                ( 1.0D0 + COUR ) * &
                0.50D0*(H0P(4,J+1) + H0P(4,J))
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF
          WKA(4,J,1) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * UBP(4,J) + COUR * UBP(5,J) + &
              UBTRW(1,J) - &
              CFF1 * ( ZSTAR - 0.50D0*(H0BTRW(1,J+1) + H0BTRW(1,J))) )
        ELSE
          WKA(4,J,1) = WKA(5,J,1)
        ENDIF
        WKA(3,J,1) = WKA(4,J,1) * VIV(3,J,1)
        WKA(2,J,1) = WKA(4,J,1) * VIV(2,J,1)
        WKA(1,J,1) = WKA(4,J,1) * VIV(1,J,1)
      ENDDO
      WKA(4,JMT,1) = WKA(4,JMT-1,1) * VIV(4,JMT,1)
      WKA(3,JMT,1) = WKA(3,JMT-1,1) * VIV(3,JMT,1)
      WKA(2,JMT,1) = WKA(2,JMT-1,1) * VIV(2,JMT,1)
      WKA(1,JMT,1) = WKA(1,JMT-1,1) * VIV(1,JMT,1)
    ENDIF
  ELSE IF ( JEWSN .EQ. 3 ) THEN
    ! Southern Boundary Condition
    IF ( ID_OBC_UB_S .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For UB at Southern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I =1, IMT
        UB(I,JET-3) = UB(I,JET-4) * VIV(I,JET-3,1)
        UB(I,JET-2) = UB(I,JET-3) * VIV(I,JET-2,1)
        UB(I,JET-1) = UB(I,JET-3) * VIV(I,JET-1,1)
        UB(I,JET) = UB(I,JET-3) * VIV(I,JET,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_S .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For UB at Southern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I =1, IMT
        UB(I,JET-3) = UBTRS(I,2)*VIV(I,JET-3,1)
        UB(I,JET-2) = UB(I,JET-3) * VIV(I,JET-2,1)
        UB(I,JET-1) = UB(I,JET-3) * VIV(I,JET-1,1)
        UB(I,JET) = UB(I,JET-3) * VIV(I,JET,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_S .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For UB at Southern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        UDTDX = VBpp(I,JET-3) * DTB_OPENBC / DYT(JET-3)

        IF ( UDTDX.LT.0.0D0 ) THEN
          UB(I,JET-3) = UBPp(I,JET-3) - UDTDX * ( UBTRS(I,2) - UBPp(I,JET-3) )
        ELSE
          UB(I,JET-3) = UBPp(I,JET-3) - UDTDX * ( UBPp(I,JET-3) - UBPp(I,JET-4) )
        ENDIF
        UB(I,JET-2) = UB(I,JET-3) * VIV(I,JET-2,1)
        UB(I,JET-1) = UB(I,JET-3) * VIV(I,JET-1,1)
        UB(I,JET)   = UB(I,JET-3) * VIV(I,JET,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_S .EQ. 4 ) THEN
      ! The Implicit Chapman Boundary Condition For UB at Southern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        IF ( OHBU(I,JET-4).GT. 0.0D0 .AND. OHBU(I,JET-4).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,JET-4) + &
              0.25D0*( H0Pp(I-1,JET-4) * VIT(I-1,JET-4,1) + H0Pp(I,JET-4) * VIT(I,JET-4,1) +&
              H0Pp(I-1,JET-3) * VIT(I-1,JET-3,1) + H0Pp(I,JET-3) * VIT(I,JET-3,1) )

          CMIU = SQRT( G*HTMP )*DTB_OPENBC / DYT(JET-3)
          UB(I,JET-3) = ( 0.25D0*UBPp(I-1,JET-3) + 0.50D0*UBPp(I,JET-3) + &
              0.25D0*UBPp(I+1,JET-3) ) / ( 1.0D0+CMIU ) + &
              ( 0.25D0*UB(I-1,JET-4)  + 0.50D0*UB(I,JET-4)  + &
              0.25D0*UB(I+1,JET-4)  ) * CMIU/( 1.0D0+CMIU )
        ELSE
          UB(I,JET-3) = 0.0D0
        ENDIF
        UB(I,JET-2) = UB(I,JET-3) * VIV(I,JET-2,1)
        UB(I,JET-1) = UB(I,JET-3) * VIV(I,JET-1,1)
        UB(I,JET) = UB(I,JET-3) * VIV(I,JET,1)
      ENDDO

      UB(1,JET-3) = UB(2,JET-3) * VIV(1,JET-3,1)
      UB(1,JET-2) = UB(2,JET-2) * VIV(1,JET-2,1)
      UB(1,JET-1) = UB(2,JET-1) * VIV(1,JET-1,1)
      UB(1,JET) = UB(2,JET) * VIV(1,JET,1)
      UB(IMT,JET-3) = UB(IMT-1,JET-3) * VIV(IMT,JET-3,1)
      UB(IMT,JET-2) = UB(IMT-1,JET-2) * VIV(IMT,JET-2,1)
      UB(IMT,JET-1) = UB(IMT-1,JET-1) * VIV(IMT,JET-1,1)
      UB(IMT,JET) = UB(IMT-1,JET) * VIV(IMT,JET,1)
    ELSE IF ( ID_OBC_UB_S .EQ. 5 ) THEN
      ! The Explicit Chapman Boundary Condition For UB at Southern Edge
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        IF ( OHBU(I,JET-4).GT.0.0D0 .AND. OHBU(I,JET-4).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,JET-4) + &
              0.25D0*( H0P(I-1,JET-4) * VIT(I-1,JET-4,1) + H0P(I,JET-4) * VIT(I,JET-4,1) + &
              H0P(I-1,JET-3) * VIT(I-1,JET-3,1) + H0P(I,JET-3) * VIT(I,JET-3,1) )

          CMIU = SQRT( G*HTMP )*DTB/DYT(JET-3)
          CMIU = MIN( CMIU, 0.9999 )
          UB(I,JET-3) = ( 0.25D0*UBP(I-1,JET-3) + 0.50D0*UBP(I,JET-3) + &
              0.25D0*UBP(I+1,JET-3) ) * ( 1.0D0-CMIU ) + &
              ( 0.25D0*UBP(I-1,JET-4) + 0.50D0*UBP(I,JET-4) + &
              0.25D0*UBP(I+1,JET-4) ) * CMIU
        ELSE
          UB(I,JET-3) = 0.0D0
        ENDIF
        UB(I,JET-2) = UB(I,JET-3) * VIV(I,JET-2,1)
        UB(I,JET-1) = UB(I,JET-3) * VIV(I,JET-1,1)
        UB(I,JET) = UB(I,JET-3) * VIV(I,JET,1)
      ENDDO
      UB(1,JET-3) = UB(2,JET-3) * VIV(1,JET-3,1)
      UB(1,JET-2) = UB(2,JET-2) * VIV(1,JET-2,1)
      UB(1,JET-1) = UB(2,JET-1) * VIV(1,JET-1,1)
      UB(1,JET) = UB(2,JET) * VIV(1,JET,1)
      UB(IMT,JET-3) = UB(IMT-1,JET-3) * VIV(IMT,JET-3,1)
      UB(IMT,JET-2) = UB(IMT-1,JET-2) * VIV(IMT,JET-2,1)
      UB(IMT,JET-1) = UB(IMT-1,JET-1) * VIV(IMT,JET-1,1)
      UB(IMT,JET) = UB(IMT-1,JET) * VIV(IMT,JET,1)
    ELSE IF ( ID_OBC_UB_S .EQ. 6 ) THEN
      ! The Implicit Upstream Radiation Condition For 2D UB at Southern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        DVDTTMP = UBPP(I,JET-4) - UB(I,JET-4)
        DVDNMTMP = UB(I,JET-4) - UB(I,JET-5)
        DVDTG2TMP = UBPP(I+1,JET-4) - UBPP(I-1,JET-4)

        IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
          DVDTTMP = 0.0D0
          Tau_obc = TAUIN
        ELSE
          Tau_obc = TAUOUT
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          DVDTGTMP = UBPP(I,JET-4) - UBPP(I-1,JET-4)
        ELSE
          DVDTGTMP = UBPP(I+1,JET-4) - UBPP(I,JET-4)
        ENDIF

        IF ( ID_OBC_2DRad_UB .EQ. 0 ) THEN
          DVDTGTMP = 0.0D0
        ENDIF

        CFF1 = DVDNMTMP * DVDNMTMP + &
            DVDTGTMP * DVDTGTMP * DYT(JET-4) * DYT(JET-4) / &
            ( ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) * ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) )

        CFF2 = DVDTGTMP * DVDTGTMP + &
            DVDNMTMP * DVDNMTMP / ( DYT(JET-4)*DYT(JET-4) ) * &
            ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) * ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) )

        CFF1 = MAX(CFF1,EPSCF)
        CFF2 = MAX(CFF2,EPSCF)
        CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
        CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

        IF ( ID_OBC_2DRad_NPO_UB .EQ. 1 ) THEN
          CFTG = 0.0D0
        ENDIF

        CMIU = CFF1 + CFNM
        UB(I,JET-3) = ( CFF1*UBPP(I,JET-3) + &
            CFNM*UB(I,JET-4) -  &
            MAX(CFTG,0.0D0)*(CFF1/CFF2)*(UBPP(I,JET-3)-UBPP(I-1,JET-3))- &
            MIN(CFTG,0.0D0)*(CFF1/CFF2)*(UBPP(I+1,JET-3)-UBPP(I,JET-3)) ) / CMIU

        UB(I,JET-3) = UB(I,JET-3) + DTB_OPENBC/Tau_obc * ( ( 0.2*UBTRS(I-1,2) + 0.2*UBTRS(I+1,2) + &
            0.4*UBTRS(I,2)+0.2*UBTRS(I,1) )- UBPP(I,JET-3) ) * CFF1/CMIU

        UB(I,JET-2) = UB(I,JET-3) * VIV(I,JET-2,1)
        UB(I,JET-1) = UB(I,JET-3) * VIV(I,JET-1,1)
        UB(I,JET)   = UB(I,JET-3) * VIV(I,JET,1)
      ENDDO
      UB(1,JET-3) = UB(2,JET-3) * VIV(1,JET-3,1)
      UB(1,JET-2) = UB(2,JET-2) * VIV(1,JET-2,1)
      UB(1,JET-1) = UB(2,JET-1) * VIV(1,JET-1,1)
      UB(1,JET)   = UB(2,JET) * VIV(1,JET,1)
      UB(IMT,JET-3) = UB(IMT-1,JET-3) * VIV(IMT,JET-3,1)
      UB(IMT,JET-2) = UB(IMT-1,JET-2) * VIV(IMT,JET-2,1)
      UB(IMT,JET-1) = UB(IMT-1,JET-1) * VIV(IMT,JET-1,1)
      UB(IMT,JET)   = UB(IMT-1,JET) * VIV(IMT,JMT,1)
    ELSE IF ( ID_OBC_UB_S .EQ. 7 ) THEN
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        IF (OHBU(I,JET-4) .GT. 1.0D-20) THEN
          HTMP = 1.0D0/OHBU(I,JET-4)
          CMIU = SQRT( G*HTMP ) * DTB_OPENBC * OUY(JET-3)
          UB(I,JET-3) = UBPP(I,JET-3) / ( 1.0D0+CMIU ) + &
              UB(I,JET-4) * CMIU/( 1.0D0+CMIU )
        ELSE
          UB(I,JET-3) = 0
        ENDIF
        IF (ID_UBVB_CHAPMAN_NUDGING .EQ. 1) THEN
          UB(I,JET-3) = UB(I,JET-3) + DTB_OPENBC/TAU_2DUV*(UBTRS(I,2) - UBPP(I,JET-3))
        ENDIF

        UB(I,JET-2) = UB(I,JET-3) * VIV(I,JET-2,1)
        UB(I,JET-1) = UB(I,JET-3) * VIV(I,JET-1,1)
        UB(I,JET  ) = UB(I,JET-3) * VIV(I,JET  ,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_S .EQ. 8 ) THEN
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        IF ( OHBU(I,JET-3).GT.0.0D0  ) THEN
          HTMP = 1.0D0/OHBU(I,JET-3)  ! +
          CMIU = SQRT( G*HTMP )*DTB/DYT(JET-3)
          CMIU = MIN( CMIU, 0.9999 )
          UB(I,JET-3) = UB(I,JET-3) - CMIU * ( UB(I,JET-3) - UB(I,JET-4) )
        ENDIF
        UB(I,JET-2) = UB(I,JET-3) * VIV(I,JET-2,1)
        UB(I,JET-1) = UB(I,JET-3) * VIV(I,JET-1,1)
        UB(I,JET) = UB(I,JET-3) * VIV(I,JET,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_S .EQ. 9 ) THEN
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        HTMP = 1.0D0/OHBU(I,JET-4)

        IF (HTMP .LT. 1.0D+5) THEN
          CMIU = SQRT( G*HTMP ) * DTB2 * OUY(JET-3)
          WKA(I,JET-3,1) = UBP(I,JET-3) / ( 1.0D0+CMIU ) + &
              WKA(I,JET-4,1) * CMIU/( 1.0D0+CMIU )
        ELSE
          WKA(I,JET-3,1) = WKA(I,JET-4,1)
        ENDIF
        WKA(I,JET-2,1) = WKA(I,JET-3,1) * VIV(I,JET-2,1)
        WKA(I,JET-1,1) = WKA(I,JET-3,1) * VIV(I,JET-1,1)
        WKA(I,JET  ,1) = WKA(I,JET-3,1) * VIV(I,JET  ,1)
      ENDDO

      WKA(1,JET-3,1) = WKA(2,JET-3,1) * VIV(1,JET-3,1)
      WKA(1,JET-2,1) = WKA(2,JET-2,1) * VIV(1,JET-2,1)
      WKA(1,JET-1,1) = WKA(2,JET-1,1) * VIV(1,JET-1,1)
      WKA(1,JET  ,1) = WKA(2,JET  ,1) * VIV(1,JET  ,1)
    ENDIF
  ELSE IF ( JEWSN .EQ. 4 ) THEN
    ! Northern Boundary Condition
    IF ( ID_OBC_UB_N .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For UB at Northern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        UB(I,3) = UB(I,4) * VIV(I,3,1)
        UB(I,2) = UB(I,3) * VIV(I,2,1)
        UB(I,1) = UB(I,3) * VIV(I,1,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_N .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For UB at Northern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        UB(I,3) = UBTRN(I,1)*VIV(I,3,1)
        UB(I,2) = UB(I,3) * VIV(I,2,1)
        UB(I,1) = UB(I,3) * VIV(I,1,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_N .EQ.3 ) THEN
      ! The Upstream Advection Scheme For UB at Northern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        UDTDX = VBpp(I,3) * DTB_OPENBC / DYT(3)

        IF ( UDTDX.GT.0.0D0 ) THEN
          UB(I,3) = UBPp(I,3) - UDTDX * ( UBPp(I,3) - UBTRN(I,1)  )
        ELSE
          UB(I,3) = UBPp(I,3) - UDTDX * ( UBPp(I,4) - UBPp(I,3) )
        ENDIF
        UB(I,2) = UB(I,3) * VIV(I,2,1)
        UB(I,1) = UB(I,3) * VIV(I,1,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_N .EQ.4 ) THEN
      ! The Implicit Chapman Boundary Condition For UB at Northern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        IF ( OHBU(I,4).GT.0.0D0 .AND. OHBU(I,4).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,4) + &
              0.25D0*( H0Pp(I-1,4) * VIT(I-1,4,1) + H0Pp(I-1,3) * VIT(I-1,3,1) + &
              H0Pp(I,4) * VIT(I,4,1) + H0Pp(I,3) * VIT(I,3,1) )

          CMIU = SQRT( G*HTMP )*DTB_OPENBC / DYT(4)
          UB(I,3) = ( 0.25D0*UBPp(I-1,3) + 0.50D0*UBPp(I,3) + 0.25D0*UBPp(I+1,3) ) / ( 1.0D0+CMIU ) + &
              ( 0.25D0*UB(I-1,4)  + 0.50D0*UB(I,4)  + 0.25D0*UB(I+1,4) ) * CMIU/( 1.0D0+CMIU )
        ELSE
          UB(I,3) = 0.0D0
        ENDIF
        UB(I,2) = UB(I,3) * VIV(I,2,1)
        UB(I,1) = UB(I,3) * VIV(I,1,1)
      ENDDO

      UB(1,3) = UB(2,3) * VIV(1,3,1)
      UB(1,2) = UB(2,2) * VIV(1,2,1)
      UB(1,1) = UB(2,1) * VIV(1,1,1)
      UB(IMT,3) = UB(IMT-1,3) * VIV(IMT,3,1)
      UB(IMT,2) = UB(IMT-1,2) * VIV(IMT,2,1)
      UB(IMT,1) = UB(IMT-1,1) * VIV(IMT,1,1)
    ELSE IF ( ID_OBC_UB_N .EQ.5 ) THEN
      ! The Explicit Chapman Boundary Condition For UB at Northern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        IF ( OHBU(I,4).GT.0.0D0 .AND. OHBU(I,4).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,4) + &
              0.25D0*( H0P(I-1,4) * VIT(I-1,4,1) + H0P(I,4) * VIT(I,4,1) +  &
              H0P(I-1,3) * VIT(I-1,3,1) + H0P(I,3) * VIT(I,3,1) )

          CMIU = SQRT( G*HTMP )*DTB / DYT(4)
          CMIU = MIN( CMIU, 0.9999 )
          UB(I,3) = ( 0.25D0*UBP(I-1,3) + 0.50D0*UBP(I,3) + 0.25D0*UBP(I+1,3) ) * ( 1.0D0-CMIU ) + &
              ( 0.25D0*UBP(I-1,4) + 0.50D0*UBP(I,4) + 0.25D0*UBP(I+1,4) ) * CMIU
        ELSE
          UB(I,3) = 0.0D0
        ENDIF

        UB(I,2) = UB(I,3) * VIV(I,2,1)
        UB(I,1) = UB(I,3) * VIV(I,1,1)
      ENDDO

      UB(1,3) = UB(2,3) * VIV(1,3,1)
      UB(1,2) = UB(2,2) * VIV(1,2,1)
      UB(1,1) = UB(2,1) * VIV(1,1,1)
      UB(IMT,3) = UB(IMT-1,3) * VIV(IMT,3,1)
      UB(IMT,2) = UB(IMT-1,2) * VIV(IMT,2,1)
      UB(IMT,1) = UB(IMT-1,1) * VIV(IMT,1,1)
    ELSE IF ( ID_OBC_UB_N .EQ.6 ) THEN
      ! The Implicit Upstream Radiation Condition For UB at Northern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        DVDTTMP = UBPP(I,4) - UB(I,4)
        DVDNMTMP = UB(I,4) - UB(I,5)
        DVDTG2TMP = UBPP(I+1,4) - UBPP(I-1,4)

        IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
          DVDTTMP = 0.0D0
          Tau_obc = TAUIN
        ELSE
          Tau_obc = TAUOUT
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          DVDTGTMP = UBPP(I,4) - UBPP(I-1,4)
        ELSE
          DVDTGTMP = UBPP(I+1,4) - UBPP(I,4)
        ENDIF

        IF ( ID_OBC_2DRad_UB .EQ. 0 ) THEN
          DVDTGTMP = 0.0D0
        ENDIF

        CFF1 = DVDNMTMP * DVDNMTMP + &
            DVDTGTMP * DVDTGTMP * DYT(5)*DYT(5) / &
            ( ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) * ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) )

        CFF2 = DVDTGTMP * DVDTGTMP + &
            DVDNMTMP * DVDNMTMP / ( DYT(5)*DYT(5) ) * &
            ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) * ( 0.50D0/OTX(4)+0.50D0/OTX(5) )

        CFF1 = MAX(CFF1,EPSCF)
        CFF2 = MAX(CFF2,EPSCF)
        CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
        CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

        IF ( ID_OBC_2DRad_NPO_UB .EQ. 1 ) THEN
          CFTG = 0.0D0
        ENDIF

        CMIU = CFF1 + CFNM
        UB(I,3) = ( CFF1*UBPP(I,3) + &
            CFNM*UB(I,4) -  &
            MAX(CFTG,0.0D0)*(CFF1/CFF2)*(UBPP(I,3)-UBPP(I-1,3)) - &
            MIN(CFTG,0.0D0)*(CFF1/CFF2)*(UBPP(I+1,3)-UBPP(I,3)) ) / CMIU

        UB(I,3) = UB(I,3) + DTB_OPENBC/Tau_obc * ( ( 0.2*UBTRN(I-1,1)+0.2*UBTRN(I+1,1)+ &
            0.4*UBTRN(I,1)+0.2*UBTRN(I,2) )- UBPP(I,3) ) * CFF1/CMIU

        UB(I,2) = UB(I,3) * VIV(I,2,1)
        UB(I,1) = UB(I,3) * VIV(I,1,1)
      ENDDO
      UB(1,3) = UB(2,3) * VIV(1,3,1)
      UB(1,2) = UB(2,2) * VIV(1,2,1)
      UB(1,1) = UB(2,1) * VIV(1,1,1)
      UB(IMT,3) = UB(IMT-1,3) * VIV(IMT,3,1)
      UB(IMT,2) = UB(IMT-1,2) * VIV(IMT,2,1)
      UB(IMT,1) = UB(IMT-1,1) * VIV(IMT,1,1)
    ELSE IF ( ID_OBC_UB_N .EQ. 7 ) THEN
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        IF (OHBU(I,4) .GT. 1.0D-20) THEN
          HTMP = 1.0D0/OHBU(I,4)
          CMIU = SQRT( G*HTMP )* DTB_OPENBC * OUY(4)
          UB(I,3) = UBPP(I,3) / ( 1.0D0+CMIU ) + &
              UB(I,4) * CMIU/( 1.0D0+CMIU )
        ELSE
          UB(I,3) = 0
        ENDIF
        IF (ID_UBVB_CHAPMAN_NUDGING .EQ. 1) THEN
          UB(I,3) = UB(I,3) + DTB_OPENBC/TAU_2DUV*(UBTRN(I,1) - UBPP(I,3))
        ENDIF
        UB(I,2) = UB(I,3) * VIV(I,2,1)
        UB(I,1) = UB(I,3) * VIV(I,1,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_N .EQ. 8 ) THEN
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        IF ( OHBU(I,3).GT.0.0D0 ) THEN
          HTMP = 1.0D0/OHBU(I,3)
          CMIU = SQRT( G*HTMP )*DTB / DYT(4)
          CMIU = MIN( CMIU, 0.9999 )
          UB(I,3) = UB(I,3) - CMIU * ( UB(I,3) - UB(I,4) )
        ELSE
          UB(I,3) = 0.0D0
        ENDIF

        UB(I,2) = UB(I,3) * VIV(I,2,1)
        UB(I,1) = UB(I,3) * VIV(I,1,1)
      ENDDO
    ELSE IF ( ID_OBC_UB_N .EQ. 9 ) THEN
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        HTMP = 1.0D0/OHBU(I,4)

        IF (HTMP .LT. 1.0D+5) THEN
          CMIU = SQRT( G*HTMP )*DTB2 * OUY(4)
          WKA(I,3,1) = UBP(I,3) / ( 1.0D0+CMIU ) + &
              WKA(I,4,1) * CMIU/( 1.0D0+CMIU )
        ELSE
          WKA(I,3,1) = WKA(I,4,1)
        ENDIF
        WKA(I,2,1) = WKA(I,3,1) * VIV(I,2,1)
        WKA(I,1,1) = WKA(I,3,1) * VIV(I,1,1)
      ENDDO

      WKA(1,3,1) = WKA(2,3,1) * VIV(1,3,1)
      WKA(1,2,1) = WKA(2,2,1) * VIV(1,2,1)
      WKA(1,1,1) = WKA(2,1,1) * VIV(1,1,1)
    ENDIF
  ENDIF

  !$OMP PARALLEL DO PRIVATE (J,I)
  DO J=1, JMT
    DO I=1, IMT
      UB(I,J) = UB(I,J) * VIV(I,J,1)
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE OBC_2DUB

SUBROUTINE OBC_2DVB(JEWSN)

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, KM, JET,I, J, K
  use   pconst_mod, only: VIT, VIV, G, OHBU, OHBT, DTB, DTB2, DYT, DYR, OTX, DZP, OUX, iy
  use   dyn_mod,    only: UB, VB, H0, H0P, UBP, VBP
  use   work_mod

#if ( defined SPMD )
  use msg_mod
#endif

#ifdef OBCDT
  use   openbc_mod
#endif

  IMPLICIT NONE
  INTEGER :: JEWSN
  REAL(r8) :: EPSCF
  REAL(r8) :: Tau_obc
  REAL(r8) :: HTMP
  REAL(r8) :: UDTDX
  REAL(r8) :: DVDTTMP, DVDNMTMP, DVDTGTMP, DVDTG2TMP, CFTG, CFNM
  REAL(r8) :: CMIU
  REAL(r8) :: CFF, CFF1, CFF2, CFF3
  REAL(r8) :: Chalf, COUR
  REAL(r8) :: ZSTAR

  EPSCF = 1.0D-20

  IF ( JEWSN .EQ. 1 ) THEN
    ! Eastern Boundary Condition
    IF ( ID_OBC_VB_E .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For VB at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        VB(IMT-2,J) = VB(IMT-3,J)*VIV(IMT-2,J,1)
        VB(IMT-1,J) = VB(IMT-2,J) * VIV(IMT-1,J,1)
        VB(IMT,J) = VB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_E .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For VB at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        VB(IMT-2,J) = VBTRE(2,J)*VIV(IMT-2,J,1)
        VB(IMT-1,J) = VB(IMT-2,J) * VIV(IMT-1,J,1)
        VB(IMT,J) = VB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_E .EQ. 3) THEN
      ! The Upstream Advection Scheme For VB at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        UDTDX = UBpp(IMT-2,J) * DTB_OPENBC * OTX(J)

        IF ( UDTDX.LT.0.0D0 ) THEN
          VB(IMT-2,J) = VBPp(IMT-2,J) - UDTDX * ( VBTRE(2,J) - VBPp(IMT-2,J) )
        ELSE
          VB(IMT-2,J) = VBPp(IMT-2,J) - UDTDX * ( VBPp(IMT-2,J) - VBPp(IMT-3,J) )
        ENDIF

        VB(IMT-1,J) = VB(IMT-2,J) * VIV(IMT-1,J,1)
        VB(IMT,J) = VB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_E .EQ. 4 ) THEN
      ! The Implicit Chapman Boundary Condition For VB at Eastern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        IF ( OHBU(IMT-3,J).GT.0.0D0 .AND. OHBU(IMT-3,J).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(IMT-3,J) + &
              0.25D0*( H0Pp(IMT-4,J) * VIT(IMT-4,J,1) + H0Pp(IMT-4,J+1) * VIT(IMT-4,J+1,1) + &
              H0Pp(IMT-3,J) * VIT(IMT-3,J,1) + H0Pp(IMT-3,J+1) * VIT(IMT-3,J+1,1) )

          CMIU = SQRT( G*HTMP )*DTB_OPENBC / ( 0.5D0/OTX(J) +0.5D0/OTX(J+1) )
          VB(IMT-2,J) = ( 0.25D0*VBPp(IMT-2,J-1) + 0.50D0*VBPp(IMT-2,J) + &
              0.25D0*VBPp(IMT-2,J+1) ) / (1.0D0+CMIU ) + &
              ( 0.25D0*VB(IMT-3,J-1)  + 0.50D0*VB(IMT-3,J)  + &
              0.25D0*VB(IMT-3,J+1)  ) * CMIU / ( 1.0D0 + CMIU )
        ELSE
          VB(IMT-2,J) = 0.0D0
        ENDIF

        VB(IMT-1,J) = VB(IMT-2,J) * VIV(IMT-1,J,1)
        VB(IMT,J)   = VB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO

      VB(IMT-2,1) = VB(IMT-2,2) * VIV(IMT-2,1,1)
      VB(IMT-1,1) = VB(IMT-1,2) * VIV(IMT-1,1,1)
      VB(IMT,1) = VB(IMT,2) * VIV(IMT,1,1)
      VB(IMT-2,JMT) = VB(IMT-2,JMT-1) * VIV(IMT-2,JMT,1)
      VB(IMT-1,JMT) = VB(IMT-1,JMT-1) * VIV(IMT-1,JMT,1)
      VB(IMT,JMT) = VB(IMT,JMT-1) * VIV(IMT,JMT,1)
    ELSE IF ( ID_OBC_VB_E .EQ. 5 ) THEN
      ! The Explicit Chapman Boundary Condition For VB at Eastern Edge
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        IF ( OHBU(IMT-3,J).GT.0.0D0 .AND. OHBU(IMT-3,J).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(IMT-3,J) + &
              0.25D0*( H0P(IMT-4,J) * VIT(IMT-4,J,1) + H0P(IMT-4,J+1) * VIT(IMT-4,J+1,1)  + &
              H0P(IMT-3,J) * VIT(IMT-3,J,1) + H0P(IMT-3,J+1) * VIT(IMT-3,J+1,1) )

          CMIU = SQRT( G*HTMP )*DTB / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          CMIU = MIN( CMIU, 0.9999 )
          VB(IMT-2,J) = ( 0.25D0*VBP(IMT-2,J-1) + 0.50D0*VBP(IMT-2,J) + &
              0.25D0*VBP(IMT-2,J+1) ) * ( 1.0D0-CMIU ) + &
              ( 0.25D0*VBP(IMT-3,J-1) + 0.50D0*VBP(IMT-3,J) + &
              0.25D0*VBP(IMT-3,J+1) ) * CMIU
        ELSE
          VB(IMT-2,J) = 0.0D0
        ENDIF

        VB(IMT-1,J) = VB(IMT-2,J) * VIV(IMT-1,J,1)
        VB(IMT,J)   = VB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO

      VB(IMT-2,1) = VB(IMT-2,2) * VIV(IMT-2,1,1)
      VB(IMT-1,1) = VB(IMT-1,2) * VIV(IMT-1,1,1)
      VB(IMT,1) = VB(IMT,2) * VIV(IMT,1,1)
      VB(IMT-2,JMT) = VB(IMT-2,JMT-1) * VIV(IMT-2,JMT,1)
      VB(IMT-1,JMT) = VB(IMT-1,JMT-1) * VIV(IMT-1,JMT,1)
      VB(IMT,JMT) = VB(IMT,JMT-1) * VIV(IMT,JMT,1)
    ELSE IF ( ID_OBC_VB_E .EQ. 6 ) THEN
      ! The Implicit Upstream Radiation Condition For VB at Eastern Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        DVDTTMP = VBPP(IMT-3,J) - VB(IMT-3,J)
        DVDNMTMP = VB(IMT-3,J) - VB(IMT-4,J)
        DVDTG2TMP = VBPP(IMT-3,J+1) - VBPP(IMT-3,J-1)

        IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
          DVDTTMP = 0.0D0
          Tau_obc = TAUIN
        ELSE
          Tau_obc = TAUOUT
        ENDIF

        IF ( ( DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          DVDTGTMP = VBPP(IMT-3,J) - VBPP(IMT-3,J-1)
        ELSE
          DVDTGTMP = VBPP(IMT-3,J+1) - VBPP(IMT-3,J)
        ENDIF

        IF ( ID_OBC_2DRad_VB .EQ. 0 ) THEN
          DVDTGTMP = 0.0D0
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP / ( DYT(J)*DYT(J) ) * &
              ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP * DYT(J) * DYT(J) / &
              ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
        ELSE
          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP / ( DYT(J+1)*DYT(J+1) ) * &
              ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

          CFF1 = MAX(CFF1,EPSCF)
          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP * DYT(J+1) * DYT(J+1) / &
              ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
        ENDIF
        CFF1 = MAX(CFF1,EPSCF)
        CFF2 = MAX(CFF2,EPSCF)
        CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
        CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

        IF ( ID_OBC_2DRad_NPO_VB .EQ. 1 ) THEN
          CFTG = 0.0D0
        ENDIF

        CMIU = CFF1 + CFNM
        VB(IMT-2,J) = ( CFF1*VBPP(IMT-2,J) + CFNM*VB(IMT-3,J) - &
            MAX(CFTG,0.0D0)*(CFF1/CFF2)*(VBPP(IMT-2,J)-VBPP(IMT-2,J-1))- &
            MIN(CFTG,0.0D0)*(CFF1/CFF2)*(VBPP(IMT-2,J+1)-VBPP(IMT-2,J)))/CMIU

        VB(IMT-2,J) = VB(IMT-2,J) + DTB_OPENBC/Tau_obc * ( (VBTRE(1,J)*0.2+VBTRE(2,J)*0.4+ &
            VBTRE(2,J-1)*0.2+VBTRE(2,J+1)*0.2)-VBPP(IMT-2,J) ) * CFF1/CMIU

        VB(IMT-1,J) = VB(IMT-2,J) * VIV(IMT-1,J,1)
        VB(IMT,J) = VB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
      VB(IMT-2,1) = VB(IMT-2,2) * VIV(IMT-2,1,1)
      VB(IMT-1,1) = VB(IMT-1,2) * VIV(IMT-1,1,1)
      VB(IMT,1) = VB(IMT,2) * VIV(IMT,1,1)
      VB(IMT-2,JMT) = VB(IMT-2,JMT-1) * VIV(IMT-2,JMT,1)
      VB(IMT-1,JMT) = VB(IMT-1,JMT-1) * VIV(IMT-1,JMT,1)
      VB(IMT,JMT) = VB(IMT,JMT-1) * VIV(IMT,JMT,1)
    ELSE IF ( ID_OBC_VB_E .EQ. 7 ) THEN
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        IF (OHBU(IMT-3,J) .GT. 1.0D-20) THEN
          HTMP = 1.0D0/OHBU(IMT-3,J)
          CMIU = SQRT( G*HTMP )* DTB_OPENBC * OUX(J)
          VB(IMT-2,J) = VBPP(IMT-2,J) / (1.0D0+CMIU ) + &
              VB(IMT-3,J) * CMIU / ( 1.0D0 + CMIU )
        ELSE
          VB(IMT-2,J) = 0
        ENDIF
        IF (ID_UBVB_CHAPMAN_NUDGING .EQ. 1) THEN
          VB(IMT-2,J) = VB(IMT-2,J) + DTB_OPENBC/TAU_2DUV*(VBTRE(2,J) - VBPP(IMT-2,J))
        ENDIF
        VB(IMT-1,J) = VB(IMT-2,J) * VIV(IMT-1,J,1)
        VB(IMT,J)   = VB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_E .EQ. 8 ) THEN
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        IF ( OHBU(IMT-2,J).GT.0.0D0 ) THEN
          HTMP = 1.0D0/OHBU(IMT-2,J)
          CMIU = SQRT( G*HTMP )*DTB / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          CMIU = MIN( CMIU, 0.9999 )
          VB(IMT-2,J) = ( 1.0D0 - CMIU ) * VB(IMT-2,J) + CMIU * VB(IMT-3,J)
        ELSE
          VB(IMT-2,J) = 0.0D0
        ENDIF

        VB(IMT-1,J) = VB(IMT-2,J) * VIV(IMT-1,J,1)
        VB(IMT,J)   = VB(IMT-2,J) * VIV(IMT,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_E .EQ. 9 ) THEN
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        HTMP = 1.0D0/OHBU(IMT-3,J)

        IF (HTMP .LT. 1.0D+5) THEN
          CMIU = SQRT( G*HTMP )*DTB2 * OUX(J)
          WKA(IMT-2,J,2) = VBP(IMT-2,J) / (1.0D0+CMIU ) + &
              WKA(IMT-3,J,2) * CMIU / ( 1.0D0 + CMIU )
        ELSE
          WKA(IMT-2,J,2) = WKA(IMT-3,J,2)
        ENDIF

        WKA(IMT-1,J,2) = WKA(IMT-2,J,2) * VIV(IMT-1,J,1)
        WKA(IMT,J,2)   = WKA(IMT-2,J,2) * VIV(IMT,J,1)
      ENDDO
      WKA(IMT-2,1,2) = WKA(IMT-2,2,2) * VIV(IMT-2,1,1)
      WKA(IMT-1,1,2) = WKA(IMT-1,2,2) * VIV(IMT-1,1,1)
      WKA(IMT  ,1,2) = WKA(IMT  ,2,2) * VIV(IMT,1,1)
    ENDIF
  ELSE IF ( JEWSN .EQ. 2 ) THEN
    ! Western Boundary Condition
    IF ( ID_OBC_VB_W .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For VB at Western Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        VB(4,J) = VB(5,J) * VIV(4,J,1)
        VB(3,J) = VB(4,J) * VIV(3,J,1)
        VB(2,J) = VB(4,J) * VIV(2,J,1)
        VB(1,J) = VB(4,J) * VIV(1,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_W .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For VB at Western Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        VB(4,J) = VBTRW(1,J)*VIV(4,J,1)
        VB(3,J) = VB(4,J) * VIV(3,J,1)
        VB(2,J) = VB(4,J) * VIV(2,J,1)
        VB(1,J) = VB(4,J) * VIV(1,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_W .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For VB at Western Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        UDTDX = UBpp(4,J) * DTB_OPENBC * OTX(J)
        IF ( UDTDX.GT.0.0D0 ) THEN
          VB(4,J) = VBPp(4,J) - UDTDX * ( VBPp(4,J) - VBTRW(1,J) )
        ELSE
          VB(4,J) = VBPp(4,J) - UDTDX * ( VBPp(5,J) - VBPp(4,J) )
        ENDIF
        VB(3,J) = VB(4,J) * VIV(3,J,1)
        VB(2,J) = VB(4,J) * VIV(2,J,1)
        VB(1,J) = VB(4,J) * VIV(1,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_W .EQ. 4 ) THEN
      ! The Implicit Chapman Boundary Condition For VB at Western Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        IF ( OHBU(5,J).GT.0.0D0 .AND. OHBU(5,J).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(5,J) + &
              0.25D0*( H0P(4,J) * VIT(4,J,1) + H0P(4,J+1) * VIT(4,J+1,1) +  &
              H0P(5,J) * VIT(5,J,1) + H0P(5,J+1) * VIT(5,J+1,1) )

          CMIU = SQRT( G*HTMP )*DTB / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          VB(4,J) = ( 0.25D0*VBP(4,J-1) + 0.50D0*VBP(4,J) + 0.25D0*VBP(4,J+1) ) / ( 1.0D0+CMIU ) + &
              ( 0.25D0*VB(5,J-1)  + 0.50D0*VB(5,J)  + 0.25D0*VB(5,J+1)  ) *  CMIU / ( 1.0D0+CMIU )
        ELSE
          VB(4,J) = 0.0D0
        ENDIF
        VB(3,J) = VB(4,J) * VIV(3,J,1)
        VB(2,J) = VB(4,J) * VIV(2,J,1)
        VB(1,J) = VB(4,J) * VIV(1,J,1)
      ENDDO

      VB(4,1) = VB(4,2) * VIV(4,1,1)
      VB(3,1) = VB(3,2) * VIV(3,1,1)
      VB(2,1) = VB(2,2) * VIV(2,1,1)
      VB(1,1) = VB(1,2) * VIV(1,1,1)
      VB(4,JMT) = VB(4,JMT-1) * VIV(4,JMT,1)
      VB(3,JMT) = VB(3,JMT-1) * VIV(3,JMT,1)
      VB(2,JMT) = VB(2,JMT-1) * VIV(2,JMT,1)
      VB(1,JMT) = VB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_VB_W .EQ. 5 ) THEN
      ! The Explicit Chapman Boundary Condition For VB at Western Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        IF ( OHBU(5,J).GT.0.0D0 .AND. OHBU(5,J).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(5,J) + &
              0.25D0*( H0P(4,J) * VIT(4,J,1) + H0P(4,J+1) * VIT(4,J+1,1) + &
              H0P(5,J) * VIT(5,J,1) + H0P(5,J+1) * VIT(5,J+1,1) )

          CMIU = SQRT( G*HTMP )*DTB / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          CMIU = MIN( CMIU, 0.9999 )
          VB(4,J) = ( 0.25D0*VBP(4,J-1) + 0.50D0*VBP(4,J) + 0.25D0*VBP(4,J+1) ) * ( 1.0D0-CMIU ) + &
              ( 0.25D0*VBP(5,J-1) + 0.50D0*VBP(5,J) + 0.25D0*VBP(5,J+1) ) *  CMIU
        ELSE
          VB(4,J) = 0.0D0
        ENDIF
        VB(3,J) = VB(4,J) * VIV(3,J,1)
        VB(2,J) = VB(4,J) * VIV(2,J,1)
        VB(1,J) = VB(4,J) * VIV(1,J,1)
      ENDDO

      VB(4,1) = VB(4,2) * VIV(4,1,1)
      VB(3,1) = VB(3,2) * VIV(3,1,1)
      VB(2,1) = VB(2,2) * VIV(2,1,1)
      VB(1,1) = VB(1,2) * VIV(1,1,1)
      VB(4,JMT) = VB(4,JMT-1) * VIV(4,JMT,1)
      VB(3,JMT) = VB(3,JMT-1) * VIV(3,JMT,1)
      VB(2,JMT) = VB(2,JMT-1) * VIV(2,JMT,1)
      VB(1,JMT) = VB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_VB_W .EQ. 6 ) THEN
      ! The Implicit Upstream Radiation Condition For VB at Western Edge

      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 2, JMT-1
        DVDTTMP = VBPP(5,J) - VB(5,J)
        DVDNMTMP = VB(5,J) - VB(6,J)
        DVDTG2TMP = VBPP(5,J+1) - VBPP(5,J-1)

        IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
          DVDTTMP = 0.0D0
          Tau_obc = TAUIN
        ELSE
          Tau_obc = TAUOUT
        ENDIF

        IF ( ( DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          DVDTGTMP = VBPP(5,J) - VBPP(5,J-1)
        ELSE
          DVDTGTMP = VBPP(5,J+1) - VBPP(5,J)
        ENDIF

        IF ( ID_OBC_2DRad_VB .EQ. 0 ) THEN
          DVDTGTMP = 0.0D0
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP / ( DYT(J)*DYT(J) ) * &
              ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP * DYT(J) * DYT(J) / &
              ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
        ELSE
          CFF1 = DVDNMTMP * DVDNMTMP + &
              DVDTGTMP * DVDTGTMP / ( DYT(J+1)*DYT(J+1) ) * &
              ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) )

          CFF2 = DVDTGTMP * DVDTGTMP + &
              DVDNMTMP * DVDNMTMP * DYT(J+1) * DYT(J+1) / &
              ( ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) * ( 0.50D0/OTX(J)+0.50D0/OTX(J+1) ) )
        ENDIF

        CFF1 = MAX(CFF1,EPSCF)
        CFF2 = MAX(CFF2,EPSCF)
        CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
        CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

        IF ( ID_OBC_2DRad_NPO_VB .EQ. 1 ) THEN
          CFTG = 0.0D0
        ENDIF

        CMIU = CFF1 + CFNM
        VB(4,J) = ( CFF1*VBPP(4,J) + &
            CFNM*VB(5,J) -  &
            MAX(CFTG,0.0D0)*(CFF1/CFF2)*(VBPP(4,J)-VBPP(4,J-1)) - &
            MIN(CFTG,0.0D0)*(CFF1/CFF2)*(VBPP(4,J+1)-VBPP(4,J)) ) / CMIU

        VB(4,J) = VB(4,J) + DTB_OPENBC/Tau_obc * ( (VBTRW(1,J-1)*0.2+VBTRW(1,J+1)*0.2+ &
            VBTRW(1,J)*0.4+VBTRW(2,J)*0.2)-VBPP(4,J) ) * CFF1/CMIU

        VB(3,J) = VB(4,J) * VIV(3,J,1)
        VB(2,J) = VB(4,J) * VIV(2,J,1)
        VB(1,J) = VB(4,J) * VIV(1,J,1)
      ENDDO
      VB(4,1) = VB(4,2) * VIV(4,1,1)
      VB(3,1) = VB(3,2) * VIV(3,1,1)
      VB(2,1) = VB(2,2) * VIV(2,1,1)
      VB(1,1) = VB(1,2) * VIV(1,1,1)
      VB(4,JMT) = VB(4,JMT-1) * VIV(4,JMT,1)
      VB(3,JMT) = VB(3,JMT-1) * VIV(3,JMT,1)
      VB(2,JMT) = VB(2,JMT-1) * VIV(2,JMT,1)
      VB(1,JMT) = VB(1,JMT-1) * VIV(1,JMT,1)
    ELSE IF ( ID_OBC_VB_W .EQ. 7 ) THEN
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        IF (OHBU(5,J) .GT. 1.0D-20) THEN
          HTMP = 1.0D0/OHBU(5,J)
          CMIU = SQRT( G*HTMP )* DTB_OPENBC * OUX(J)
          VB(4,J) = VBPP(4,J) / ( 1.0D0+CMIU ) + &
              VB(5,J) *  CMIU / ( 1.0D0+CMIU )
        ELSE
          VB(4,J) = 0
        ENDIF
        IF (ID_UBVB_CHAPMAN_NUDGING .EQ. 1) THEN
          VB(4,J) = VB(4,J) + DTB_OPENBC/TAU_2DUV*(VBTRW(1,J) - VBPP(4,J))
        ENDIF
        VB(3,J) = VB(4,J) * VIV(3,J,1)
        VB(2,J) = VB(4,J) * VIV(2,J,1)
        VB(1,J) = VB(4,J) * VIV(1,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_W .EQ. 8 ) THEN
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        IF ( OHBU(4,J).GT.0.0D0  ) THEN
          HTMP = 1.0D0/OHBU(4,J)
          CMIU = SQRT( G*HTMP )*DTB / ( 0.5D0/OTX(J) + 0.5D0/OTX(J+1) )
          CMIU = MIN( CMIU, 0.9999 )
          VB(4,J) = ( 1.0D0 - CMIU ) * VB(4,J) + CMIU * VB(5,J)
        ELSE
          VB(4,J) = 0.0D0
        ENDIF
        VB(3,J) = VB(4,J) * VIV(3,J,1)
        VB(2,J) = VB(4,J) * VIV(2,J,1)
        VB(1,J) = VB(4,J) * VIV(1,J,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_W .EQ. 9 ) THEN
      !$OMP PARALLEL DO PRIVATE (J)
      DO J = 1, JMT
        HTMP = 1.0D0/OHBU(5,J)

        IF (HTMP .LT. 1.0D+5) THEN
          CMIU = SQRT( G*HTMP )*DTB2 * OUX(J)
          WKA(4,J,2) = VBP(4,J) / ( 1.0D0+CMIU ) + &
              WKA(5,J,2) *  CMIU / ( 1.0D0+CMIU )
        ELSE
          WKA(4,J,2) = WKA(5,J,2)
        ENDIF
        WKA(3,J,2) = WKA(4,J,2) * VIV(3,J,1)
        WKA(2,J,2) = WKA(4,J,2) * VIV(2,J,1)
        WKA(1,J,2) = WKA(4,J,2) * VIV(1,J,1)
      ENDDO

      WKA(4,1,2) = WKA(4,2,2) * VIV(4,1,1)
      WKA(3,1,2) = WKA(3,2,2) * VIV(3,1,1)
      WKA(2,1,2) = WKA(2,2,2) * VIV(2,1,1)
      WKA(1,1,2) = WKA(1,2,2) * VIV(1,1,1)
    ENDIF
  ELSE IF ( JEWSN .EQ. 3 ) THEN
    ! Southern Boundary Condition

    IF ( ID_OBC_VB_S .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For VB at Southern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I =1, IMT
        VB(I,JET-3) = VB(I,JET-4) * VIV(I,JET-3,1)
        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET) = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_S .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For VB at Southern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I =1, IMT
        VB(I,JET-3) = VBTRS(I,2)*VIV(I,JET-3,1)
        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET) = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_S .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For VB at Southern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        UDTDX = VB(I,JET-3) * DTB / DYT(JET-3)

        IF ( UDTDX.LT.0.0D0 ) THEN
          VB(I,JET-3) = VBP(I,JET-3) - UDTDX * ( VBTRS(I,2) -VBP(I,JET-3) )
        ELSE
          VB(I,JET-3) = VBP(I,JET-3) - UDTDX * ( VBP(I,JET-3) - VBP(I,JET-4) )
        ENDIF
        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET)   = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_S .EQ. 4 ) THEN
      ! The Implicit Chapman Boundary Condition For VB at Southern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        IF ( OHBU(I,JET-4).GT.0.0D0 .AND. OHBU(I,JET-4).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,JET-4) + &
              0.25D0*( H0P(I-1,JET-4) * VIT(I-1,JET-4,1) + H0P(I-1,JET-3) * VIT(I-1,JET-3,1) +  &
              H0P(I,JET-4) * VIT(I,JET-4,1) + H0P(I,JET-3) * VIT(I,JET-3,1) )

          CMIU = SQRT( G*HTMP )*DTB / DYT(JET-3)
          VB(I,JET-3) = (0.25D0*VBP(I-1,JET-3)+0.50D0*VBP(I,JET-3)+0.25D0*VBP(I+1,JET-3))/(1.0D0+CMIU) + &
              (0.25D0*VB(I-1,JET-4)+0.50D0*VB(I,JET-4)+0.25D0*VB(I+1,JET-4))*CMIU/(1.0D0+CMIU)
        ELSE
          VB(I,JET-3) = 0.0D0
        ENDIF
        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET) = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO
      VB(1,JET-3) = VB(2,JET-3) * VIV(1,JET-3,1)
      VB(1,JET-2) = VB(2,JET-2) * VIV(1,JET-2,1)
      VB(1,JET-1) = VB(2,JET-1) * VIV(1,JET-1,1)
      VB(1,JET) = VB(2,JET) * VIV(1,JET,1)
      VB(IMT,JET-3) = VB(IMT-1,JET-3) * VIV(IMT,JET-3,1)
      VB(IMT,JET-2) = VB(IMT-1,JET-2) * VIV(IMT,JET-2,1)
      VB(IMT,JET-1) = VB(IMT-1,JET-1) * VIV(IMT,JET-1,1)
      VB(IMT,JET) = VB(IMT-1,JET) * VIV(IMT,JET,1)
    ELSE IF ( ID_OBC_VB_S .EQ. 5 ) THEN
      ! The Explicit Chapman Boundary Condition For VB at Southern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        IF ( OHBU(I,JET-4).GT.0.0D0 .AND. OHBU(I,JET-4).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,JET-4) + &
              0.25D0*( H0P(I-1,JET-4) * VIT(I-1,JET-4,1) + H0P(I,JET-4) * VIT(I,JET-4,1) + &
              H0P(I-1,JET-3) * VIT(I-1,JET-3,1) + H0P(I,JET-3) * VIT(I,JET-3,1) )

          CMIU = SQRT( G*HTMP )*DTB / DYT(JET-3)
          CMIU = MIN( CMIU, 0.9999 )
          VB(I,JET-3) = (0.25D0*VBP(I-1,JET-3)+0.50D0*VBP(I,JET-3)+0.25D0*VBP(I+1,JET-3))*(1.0D0-CMIU)+ &
              ( 0.25D0*VBP(I-1,JET-4) + 0.50D0*VBP(I,JET-4) + 0.25D0*VBP(I+1,JET-4) ) * CMIU
        ELSE
          VB(I,JET-3) = 0.0D0
        ENDIF
        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET) = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO

      VB(1,JET-3) = VB(2,JET-3) * VIV(1,JET-3,1)
      VB(1,JET-2) = VB(2,JET-2) * VIV(1,JET-2,1)
      VB(1,JET-1) = VB(2,JET-1) * VIV(1,JET-1,1)
      VB(1,JET) = VB(2,JET) * VIV(1,JET,1)
      VB(IMT,JET-3) = VB(IMT-1,JET-3) * VIV(IMT,JET-3,1)
      VB(IMT,JET-2) = VB(IMT-1,JET-2) * VIV(IMT,JET-2,1)
      VB(IMT,JET-1) = VB(IMT-1,JET-1) * VIV(IMT,JET-1,1)
      VB(IMT,JET) = VB(IMT-1,JET) * VIV(IMT,JET,1)
    ELSE IF ( ID_OBC_VB_S .EQ. 6 ) THEN
      ! The Implicit Upstream Radiation Condition For 2D VB at Southern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        DVDTTMP = VBP(I,JET-4) - VB(I,JET-4)
        DVDNMTMP = VB(I,JET-4) - VB(I,JET-5)
        DVDTG2TMP = VBP(I+1,JET-4) - VBP(I-1,JET-4)

        IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
          DVDTTMP = 0.0D0
          Tau_obc = TAUIN
        ELSE
          Tau_obc = TAUOUT
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          DVDTGTMP = VBP(I,JET-4) - VBP(I-1,JET-4)
        ELSE
          DVDTGTMP = VBP(I+1,JET-4) - VBP(I,JET-4)
        ENDIF

        IF ( ID_OBC_2DRad_VB .EQ. 0 ) THEN
          DVDTGTMP = 0.0D0
        ENDIF

        CFF1 = DVDNMTMP * DVDNMTMP + &
            DVDTGTMP * DVDTGTMP * DYT(JET-4) * DYT(JET-4) / &
            ( ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) * ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) )

        CFF2 = DVDTGTMP * DVDTGTMP + &
            DVDNMTMP * DVDNMTMP / ( DYT(JET-4)*DYT(JET-4) ) * &
            ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) ) * ( 0.50D0/OTX(JET-4)+0.50D0/OTX(JET-3) )

        CFF1 = MAX(CFF1,EPSCF)
        CFF2 = MAX(CFF2,EPSCF)
        CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
        CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

        IF ( ID_OBC_2DRad_NPO_VB .EQ. 1 ) THEN
          CFTG = 0.0D0
        ENDIF

        CMIU = CFF1 + CFNM
        VB(I,JET-3) = ( CFF1*VBP(I,JET-3) + &
            CFNM*VB(I,JET-4) -  &
            MAX(CFTG,0.0D0)*(CFF1/CFF2)*(VBP(I,JET-3)-VBP(I-1,JET-3)) - &
            MIN(CFTG,0.0D0)*(CFF1/CFF2)*(VBP(I+1,JET-3)-VBP(I,JET-3)) ) / CMIU

        VB(I,JET-3) = VB(I,JET-3) + DTB/Tau_obc * ( ( 0.2*VBTRS(I-1,2)+0.2*VBTRS(I+1,2)+ &
            0.4*VBTRS(I,2)+0.2*VBTRS(I,1) )- VBP(I,JET-3) ) * CFF1/CMIU

        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET)   = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO
      VB(1,JET-3) = VB(2,JET-3) * VIV(1,JET-3,1)
      VB(1,JET-2) = VB(2,JET-2) * VIV(1,JET-2,1)
      VB(1,JET-1) = VB(2,JET-1) * VIV(1,JET-1,1)
      VB(1,JET)   = VB(2,JET) * VIV(1,JET,1)
      VB(IMT,JET-3) = VB(IMT-1,JET-3) * VIV(IMT,JET-3,1)
      VB(IMT,JET-2) = VB(IMT-1,JET-2) * VIV(IMT,JET-2,1)
      VB(IMT,JET-1) = VB(IMT-1,JET-1) * VIV(IMT,JET-1,1)
      VB(IMT,JET)   = VB(IMT-1,JET) * VIV(IMT,JET,1)
    ELSE IF ( ID_OBC_VB_S .EQ. 7 ) THEN
      ! Flather Boundary Condition For VB at Southern Edge
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO
      CMIU = 1.0d-10

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT
        IF ( OHBU(I,JET-3).GT.0.0D0 .AND. OHBU(I,JET-3).LE.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,JET-3) + &
              0.25D0 * ( H0PP(I-1,JET-2) * VIT(I-1,JET-2,1) + H0PP(I,JET-2) * VIT(I,JET-2,1) + &
              H0PP(I-1,JET-3) * VIT(I-1,JET-3,1) + H0PP(I,JET-3) * VIT(I,JET-3,1) )

          VB(I,JET-3) = VBTRS(I,2) + SQRT( G / HTMP ) *       &
              (0.25D0*( H0PP(I-1,JET-2) * VIT(I-1,JET-2,1) + H0PP(I,JET-2) * VIT(I,JET-2,1) +  &
              H0PP(I-1,JET-3) * VIT(I-1,JET-3,1) + H0PP(I,JET-3) * VIT(I,JET-3,1) )- &
              0.25D0*( H0BTRS(I-1,1) + H0BTRS(I,1) +  H0BTRS(I-1,2) + H0BTRS(I,2) ) )
        ELSE
          VB(I,JET-3) = 0.0D0
        ENDIF
        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET) = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO
      VB(1,JET-3) = VB(2,JET-3) * VIV(1,JET-3,1)
      VB(1,JET-2) = VB(2,JET-2) * VIV(1,JET-2,1)
      VB(1,JET-1) = VB(2,JET-1) * VIV(1,JET-1,1)
      VB(1,JET) = VB(2,JET) * VIV(1,JET,1)
    ELSE IF ( ID_OBC_VB_S .EQ. 8 ) THEN
      ! The Modified Flather Boundary Condition For VB at Southern Edge by
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO
      CMIU = 1.0d-10

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT
        IF ( OHBU(I,JET-3).GT.0.0D0 .AND. OHBU(I,JET-3).LE.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,JET-3)
          CFF1 = SQRT(G/HTMP)
          COUR = DTB_OPENBC*CFF1*HTMP/DYR(JET-3)
          ZSTAR = ( 0.5D0+COUR ) * &
              0.5D0 * ( H0Pp(I-1,JET-3)*VIT(I-1,JET-3,1) + H0Pp(I,JET-3)*VIT(I,JET-3,1) ) + &
              ( 0.5D0-COUR ) * &
              0.5D0 * ( H0Pp(I-1,JET-2)*VIT(I-1,JET-2,1) + H0Pp(I,JET-2)*VIT(I,JET-2,1) )

          IF ( COUR.GT.Chalf ) THEN
            CFF2 = ( 1.0D0-Chalf/COUR )**2
            CFF3 = 0.5D0 * ( H0(I-1,JET-3)*VIT(I-1,JET-3,1) + H0(I,JET-3)*VIT(I,JET-3,1) ) + &
                COUR * &
                0.5D0 * ( H0Pp(I-1,JET-2)*VIT(I-1,JET-2,1) + H0Pp(I,JET-2)*VIT(I,JET-2,1) ) - &
                ( 1.0D0 + COUR ) * &
                0.5D0 * ( H0Pp(I-1,JET-3)*VIT(I-1,JET-3,1) + H0Pp(I,JET-3)*VIT(I,JET-3,1) )
            ZSTAR = ZSTAR + CFF2*CFF3
          ENDIF

          VB(I,JET-3) = 0.5D0 * &
              ( ( 1.0D0-COUR ) * VBPp(I,JET-3) + COUR * VBPp(I,JET-4) + &
              VBTRS(I,2) + &
              CFF1 * ( ZSTAR - 0.5D0 * ( H0BTRS(I-1,2) + H0BTRS(I,2) ) ) )
        ELSE
          VB(I,JET-3) = 0.0D0
        ENDIF
        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET) = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO
      VB(1,JET-3) = VB(2,JET-3) * VIV(1,JET-3,1)
      VB(1,JET-2) = VB(2,JET-2) * VIV(1,JET-2,1)
      VB(1,JET-1) = VB(2,JET-1) * VIV(1,JET-1,1)
      VB(1,JET) = VB(2,JET) * VIV(1,JET,1)
    ELSE IF ( ID_OBC_VB_S .EQ. 9 ) THEN
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT
        IF (OHBU(I,JET-4) .GT. 1.0D-20) THEN
          HTMP = 1.0D0/OHBU(I,JET-4)
          CFF1 = SQRT(G/HTMP)
          COUR = DTB_OPENBC*CFF1*HTMP/DYR(JET-3)
          ZSTAR = ( 0.5D0+COUR ) * (H0PP(I-1,JET-3) + H0PP(I,JET-3))*0.50D0 +         &
              ( 0.5D0-COUR ) * (H0PP(I-1,JET-2) + H0PP(I,JET-2))*0.50D0

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0-Chalf/COUR )**2
            CFF3 = 0.50D0*(H0(I-1,JET-3) + H0(I,JET-3))+ &
                COUR * &
                0.50D0*(H0PP(I-1,JET-2) + H0PP(I,JET-2))- &
                ( 1.0D0 + COUR ) * &
                0.50D0*(H0PP(I-1,JET-3) + H0PP(I,JET-3))
            ZSTAR = ZSTAR + CFF2*CFF3
          ENDIF

          VB(I,JET-3) = 0.5D0 * &
              ( ( 1.0D0-COUR ) * VBPP(I,JET-3) + COUR * VBPP(I,JET-4) + &
              VBTRS(I,2) + &
              CFF1 * ( ZSTAR - 0.50D0*(H0BTRS(I-1,2) + H0BTRS(I,2))) )
        ENDIF
        VB(I,JET-2) = VB(I,JET-3) * VIV(I,JET-2,1)
        VB(I,JET-1) = VB(I,JET-3) * VIV(I,JET-1,1)
        VB(I,JET) = VB(I,JET-3) * VIV(I,JET,1)
      ENDDO
      VB(1,JET-2) = VB(2,JET-2) * VIV(1,JET-2,1)
      VB(1,JET-1) = VB(2,JET-1) * VIV(1,JET-1,1)
      VB(1,JET  ) = VB(2,JET  ) * VIV(1,JET  ,1)
    ELSE IF ( ID_OBC_VB_S .EQ. 10 ) THEN
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        HTMP = 1.0D0/OHBU(I,JET-3) + &
            0.50D0 * ( H0Pp(I,JET-2) + H0Pp(I,JET-3) )
        IF (HTMP .LT. 1.0D+5 .AND. HTMP .GT. 1.0D-5) THEN
          VB(I,JET-3) = VBTRS(I,2) + SQRT( G / HTMP ) *       &
              (0.50D0*( H0Pp(I,JET-2) + H0Pp(I,JET-3)) -  H0BTRS(I,1))
        ELSE
          VB(I,JET-3) = 0
        ENDIF
        VB(I,JET-3) = VB(I,JET-3) * VIV(I,JET-3,1)
        VB(I,JET-2) = VB(I,JET-3)
        VB(I,JET-1) = VB(I,JET-3)
        VB(I,JET) = VB(I,JET-3)
      ENDDO
    ELSE IF ( ID_OBC_VB_S .EQ. 11 ) THEN
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT
        HTMP = 1.0D0/OHBU(I,JET-3)

        IF (HTMP .LT. 1.0D+5 .AND. HTMP .GT. 1.0D0) THEN
          CFF1 = SQRT(G/HTMP)
          COUR = DTB2*CFF1*HTMP/DYR(JET-3)
          ZSTAR = ( 0.5D0+COUR ) * (H0P(I-1,JET-3) + H0P(I,JET-3))*0.50D0 +         &
              ( 0.5D0-COUR ) * (H0P(I-1,JET-2) + H0P(I,JET-2))*0.50D0

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0-Chalf/COUR )**2
            CFF3 = 0.50D0*(WORK(I-1,JET-3) + WORK(I,JET-3))+ &
                COUR * &
                0.50D0*(H0P(I-1,JET-2) + H0P(I,JET-2))- &
                ( 1.0D0 + COUR ) * &
                0.50D0*(H0P(I-1,JET-3) + H0P(I,JET-3))
            ZSTAR = ZSTAR + CFF2*CFF3
          ENDIF

          WKA(I,JET-3,2) = 0.5D0 * &
              ( ( 1.0D0-COUR ) * VBP(I,JET-3) + COUR * VBP(I,JET-4) + &
              VBTRS(I,2) + &
              CFF1 * ( ZSTAR - 0.50D0*(H0BTRS(I-1,2) + H0BTRS(I,2))) )
        ELSE
          WKA(I,JET-3,2) = WKA(I,JET-4,2)
        ENDIF
        WKA(I,JET-2,2) = WKA(I,JET-3,2) * VIV(I,JET-2,1)
        WKA(I,JET-1,2) = WKA(I,JET-3,2) * VIV(I,JET-1,1)
        WKA(I,JET  ,2) = WKA(I,JET-3,2) * VIV(I,JET,1)
      ENDDO
      WKA(1,JET-3,2) = WKA(2,JET-3,2) * VIV(1,JET-3,1)
      WKA(1,JET-2,2) = WKA(2,JET-2,2) * VIV(1,JET-2,1)
      WKA(1,JET-1,2) = WKA(2,JET-1,2) * VIV(1,JET-1,1)
      WKA(1,JET  ,2) = WKA(2,JET  ,2) * VIV(1,JET  ,1)
    ENDIF
  ELSE IF ( JEWSN .EQ. 4 ) THEN
    ! Northern Boundary Condition

    IF ( ID_OBC_VB_N .EQ. 1 ) THEN
      ! The Zero Gradient Boundary Condition For VB at Northern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        VB(I,3) = VB(I,4) * VIV(I,3,1)
        VB(I,2) = VB(I,3) * VIV(I,2,1)
        VB(I,1) = VB(I,3) * VIV(I,1,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_N .EQ. 2 ) THEN
      ! The Clamped Boundary Condition For VB at Northern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        VB(I,3) = VBTRN(I,1)*VIV(I,3,1)
        VB(I,2) = VB(I,3) * VIV(I,2,1)
        VB(I,1) = VB(I,3) * VIV(I,1,1)
      ENDDO
    ELSE IF ( ID_OBC_VB_N .EQ. 3 ) THEN
      ! The Upstream Advection Scheme For VB at Northern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        UDTDX = VB(I,3) * DTB / DYT(3)

        IF ( UDTDX.GT.0.0D0 ) THEN
          VB(I,3) = VBP(I,3) - UDTDX * ( VBP(I,3) - VBTRN(I,1)  )
        ELSE
          VB(I,3) = VBP(I,3) - UDTDX * ( VBP(I,4) - VBP(I,3) )
        ENDIF
        VB(I,2) = VB(I,3)
        VB(I,1) = VB(I,3)
      ENDDO
    ELSE IF ( ID_OBC_VB_N .EQ. 4 ) THEN
      ! The Implicit Chapman Boundary Condition For VB at Northern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        IF ( OHBU(I,4).GT.0.0D0 .AND. OHBU(I,4).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,4) + &
              0.25D0*( H0P(I-1,4) * VIT(I-1,4,1) + H0P(I-1,3) * VIT(I-1,3,1) + &
              H0P(I,4) * VIT(I,4,1) + H0P(I,3) * VIT(I,3,1) )

          CMIU = SQRT( G*HTMP )*DTB / DYT(4)
          VB(I,3) = ( 0.25D0*VBP(I-1,3) + 0.50D0*VBP(I,3) + 0.25D0*VBP(I+1,3) ) / ( 1.0D0+CMIU ) + &
              ( 0.25D0*VB(I-1,4)  + 0.50D0*VB(I,4)  + 0.25D0*VB(I+1,4)  ) * CMIU/( 1.0D0+CMIU )
        ELSE
          VB(I,3) = 0.0D0
        ENDIF
        VB(I,2) = VB(I,3) * VIV(I,2,1)
        VB(I,1) = VB(I,3) * VIV(I,1,1)
      ENDDO
      VB(1,3) = VB(2,3) * VIV(1,3,1)
      VB(1,2) = VB(2,2) * VIV(1,2,1)
      VB(1,1) = VB(2,1) * VIV(1,1,1)
      VB(IMT,3) = VB(IMT-1,3) * VIV(IMT,3,1)
      VB(IMT,2) = VB(IMT-1,2) * VIV(IMT,2,1)
      VB(IMT,1) = VB(IMT-1,1) * VIV(IMT,1,1)
    ELSE IF ( ID_OBC_VB_N .EQ. 5 ) THEN
      ! The Explicit Chapman Boundary Condition For VB at Northern Edge

      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        IF ( OHBU(I,4).GT.0.0D0 .AND. OHBU(I,4).LT.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0/OHBU(I,4) + &
              0.25D0*( H0P(I-1,4) * VIT(I-1,4,1) + H0P(I,4) * VIT(I,4,1) + &
              H0P(I-1,3) * VIT(I-1,3,1) + H0P(I,3) * VIT(I,3,1) )

          CMIU = SQRT( G*HTMP )*DTB/DYT(4)
          CMIU = MIN( CMIU, 0.9999 )
          VB(I,3) = ( 0.25D0*VBP(I-1,3) + 0.50D0*VBP(I,3) + 0.25D0*VBP(I+1,3) ) * (1.0D0-CMIU ) + &
              ( 0.25D0*VBP(I-1,4) + 0.50D0*VBP(I,4) + 0.25D0*VBP(I+1,4) ) * CMIU
        ELSE
          VB(I,3) = 0.0D0
        ENDIF

        VB(I,2) = VB(I,3) * VIV(I,2,1)
        VB(I,1) = VB(I,3) * VIV(I,1,1)
      ENDDO

      VB(1,3) = VB(2,3) * VIV(1,3,1)
      VB(1,2) = VB(2,2) * VIV(1,2,1)
      VB(1,1) = VB(2,1) * VIV(1,1,1)
      VB(IMT,3) = VB(IMT-1,3) * VIV(IMT,3,1)
      VB(IMT,2) = VB(IMT-1,2) * VIV(IMT,2,1)
      VB(IMT,1) = VB(IMT-1,1) * VIV(IMT,1,1)
    ELSE IF ( ID_OBC_VB_N .EQ. 6 ) THEN
      ! The Implicit Upstream Radiation Condition For VB at Northern Edge

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT-1
        DVDTTMP = VBP(I,4) - VB(I,4)
        DVDNMTMP = VB(I,4) - VB(I,5)
        DVDTG2TMP = VBP(I+1,4) - VBP(I-1,4)

        IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
          DVDTTMP = 0.0D0
          Tau_obc = TAUIN
        ELSE
          Tau_obc = TAUOUT
        ENDIF

        IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
          DVDTGTMP = VBP(I,4) - VBP(I-1,4)
        ELSE
          DVDTGTMP = VBP(I+1,4) - VBP(I,4)
        ENDIF

        IF ( ID_OBC_2DRad_VB .EQ. 0 ) THEN
          DVDTGTMP = 0.0D0
        ENDIF

        CFF1 = DVDNMTMP * DVDNMTMP + &
            DVDTGTMP * DVDTGTMP * DYT(5)*DYT(5) / &
            ( ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) * ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) )

        CFF2 = DVDTGTMP * DVDTGTMP + &
            DVDNMTMP * DVDNMTMP / ( DYT(5)*DYT(5) ) * &
            ( 0.50D0/OTX(4)+0.50D0/OTX(5) ) * ( 0.50D0/OTX(4)+0.50D0/OTX(5) )

        CFF1 = MAX(CFF1,EPSCF)
        CFF2 = MAX(CFF2,EPSCF)
        CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
        CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

        IF ( ID_OBC_2DRad_NPO_VB .EQ. 1 ) THEN
          CFTG = 0.0D0
        ENDIF

        CMIU = CFF1 + CFNM
        VB(I,3) = ( CFF1*VBP(I,3) + &
            CFNM*VB(I,4) -  &
            MAX(CFTG,0.0D0)*(CFF1/CFF2)*(VBP(I,3)-VBP(I-1,3)) - &
            MIN(CFTG,0.0D0)*(CFF1/CFF2)*(VBP(I+1,3)-VBP(I,3)) ) / CMIU

        VB(I,3) = VB(I,3) + DTB/Tau_obc * ( ( 0.2*VBTRN(I-1,1)+0.2*VBTRN(I+1,1)+ &
            0.4*VBTRN(I,1)+0.2*VBTRN(I,2) )- VBP(I,3) ) * CFF1/CMIU

        VB(I,2) = VB(I,3) * VIV(I,2,1)
        VB(I,1) = VB(I,3) * VIV(I,1,1)
      ENDDO
      VB(1,3) = VB(2,3) * VIV(1,3,1)
      VB(1,2) = VB(2,2) * VIV(1,2,1)
      VB(1,1) = VB(2,1) * VIV(1,1,1)
      VB(IMT,3) = VB(IMT-1,3) * VIV(IMT,3,1)
      VB(IMT,2) = VB(IMT-1,2) * VIV(IMT,2,1)
      VB(IMT,1) = VB(IMT-1,1) * VIV(IMT,1,1)
    ELSE IF ( ID_OBC_VB_N .EQ. 7 ) THEN
      ! Flather Boundary Condition For VB at Northern Edge
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO
      CMIU = 1.0d-10

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT
        IF ( OHBU(I,3).GT.0.0D0 .AND. OHBU(I,3).LE.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0 / OHBU(I,3) + &
              0.25D0 * ( H0PP(I-1,3) * VIT(I-1,3,1) + H0PP(I,3) * VIT(I,3,1) +  &
              H0PP(I-1,4) * VIT(I-1,4,1) + H0PP(I,4) * VIT(I,4,1) )

          VB(I,3) = VBTRN(I,1) + SQRT( G / HTMP ) *           &
              ( 0.25D0 * ( H0PP(I-1,3) * VIT(I-1,3,1) + H0PP(I,3) * VIT(I,3,1) +   &
              H0PP(I-1,4) * VIT(I-1,4,1) + H0PP(I,4) * VIT(I,4,1) ) - &
              0.25D0 * ( H0BTRN(I-1,1) + H0BTRN(I,1) + H0BTRN(I-1,2) + H0BTRN(I,2) ) )
        ELSE
          VB(I,3) = 0.0D0
        ENDIF
        VB(I,2) = VB(I,3) * VIV(I,2,1)
        VB(I,1) = VB(I,3) * VIV(I,1,1)
      ENDDO
      VB(1,3) = VB(2,3) * VIV(1,3,1)
      VB(1,2) = VB(2,2) * VIV(1,2,1)
      VB(1,1) = VB(2,1) * VIV(1,1,1)
    ELSE IF ( ID_OBC_VB_N .EQ. 8 ) THEN
      ! The Modified Flather Boundary Condition For VB at Northern Edge by

      Chalf = 1.0D0/(2.0D0+SQRT(2.0))
      CMIU = 0.0D0
      DO K = 1, 5
        CMIU = CMIU + DZP(K)
      ENDDO
      CMIU = 1.0d-10

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT
        IF ( OHBU(I,3).GT.0.0D0 .AND. OHBU(I,3).LE.(1.0D0/CMIU) ) THEN
          HTMP = 1.0D0 / OHBU(I,3)
          CFF1 = SQRT(G/HTMP)
          COUR = DTB*CFF1*HTMP/DYR(4)
          ZSTAR = ( 0.5D0+COUR ) * &
              0.5D0 * ( H0Pp(I-1,4)*VIT(I-1,4,1) + H0Pp(I,4)*VIT(I,4,1) ) + &
              ( 0.5D0-COUR ) * &
              0.5D0 * ( H0Pp(I-1,3)*VIT(I-1,3,1) + H0Pp(I,3)*VIT(I,3,1) )
          IF ( COUR.GT.Chalf ) THEN
            CFF2 = ( 1.0D0-Chalf/COUR )**2
            CFF3 = 0.5D0 * ( H0(I-1,4)*VIT(I-1,4,1) + H0(I,4)*VIT(I,4,1) ) + &
                COUR * &
                0.5D0 * ( H0Pp(I-1,3)*VIT(I-1,3,1) + H0Pp(I,3)*VIT(I,3,1) ) - &
                ( 1.0D0 + COUR ) * &
                0.5D0 * ( H0Pp(I-1,4)*VIT(I-1,4,1) + H0Pp(I,4)*VIT(I,4,1) )
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF

          VB(I,3) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * VBPp(I,3) + COUR * VBPp(I,4) + &
              VBTRN(I,1) - &
              CFF1 * ( ZSTAR - 0.5D0 * ( H0BTRN(I-1,1) + H0BTRN(I,1) ) ) )
        ELSE
          VB(I,3) = 0.0D0
        ENDIF
        VB(I,2) = VB(I,3) * VIV(I,2,1)
        VB(I,1) = VB(I,3) * VIV(I,1,1)
      ENDDO
      VB(1,3) = VB(2,3) * VIV(1,3,1)
      VB(1,2) = VB(2,2) * VIV(1,2,1)
      VB(1,1) = VB(2,1) * VIV(1,1,1)
    ELSE IF ( ID_OBC_VB_N .EQ. 9 ) THEN
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT
        IF (OHBU(I,4) .GT. 1.0D-20) THEN
          HTMP = 1.0D0 / OHBU(I,4)
          CFF1 = SQRT(G/HTMP)
          COUR = DTB_OPENBC*CFF1*HTMP/DYR(4)
          ZSTAR = ( 0.5D0+COUR ) * (H0PP(I-1,4) + H0PP(I,4))*0.50D0 +          &
              ( 0.5D0-COUR ) * (H0PP(I-1,3) + H0PP(I,3))*0.50D0

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0-Chalf/COUR )**2
            CFF3 = 0.50D0*(H0(I-1,4) + H0(I,4)) + &
                COUR * &
                0.50D0*(H0PP(I-1,3) + H0PP(I,3)) - &
                ( 1.0D0 + COUR ) * &
                0.50D0*(H0PP(I-1,4) + H0PP(I,4))
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF

          VB(I,3) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * VBPP(I,3) + COUR * VBPP(I,4) + &
              VBTRN(I,1) - &
              CFF1 * ( ZSTAR - 0.50D0*(H0BTRN(I-1,1) + H0BTRN(I,1)) ) )
        ENDIF
        VB(I,2) = VB(I,3) * VIV(I,2,1)
        VB(I,1) = VB(I,3) * VIV(I,1,1)
      ENDDO
      VB(1,2) = VB(2,2) * VIV(1,2,1)
      VB(1,1) = VB(2,1) * VIV(1,1,1)
    ELSE IF ( ID_OBC_VB_N .EQ. 10 ) THEN
      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 1, IMT
        HTMP = 1.0D0 / OHBU(I,3) + &
            0.50D0 * ( H0Pp(I,3) + H0Pp(I,4) )
        IF (HTMP .LT. 1.0D+5 .AND. HTMP .GT. 1.0D-5) THEN
          VB(I,3) = VBTRN(I,1) - SQRT( G / HTMP ) *           &
              ( 0.50D0 * ( H0Pp(I,3) + H0Pp(I,4)) - H0BTRN(I,2))
        ELSE
          VB(I,3) = 0
        ENDIF
        VB(I,3) = VB(I,3) * VIV(I,3,1)
        VB(I,2) = VB(I,3)
        VB(I,1) = VB(I,3)
      ENDDO
    ELSE IF ( ID_OBC_VB_N .EQ. 11 ) THEN
      Chalf = 1.0D0/(2.0D0+SQRT(2.0))

      !$OMP PARALLEL DO PRIVATE (I)
      DO I = 2, IMT
        HTMP = 1.0D0 / OHBU(I,3)

        IF (HTMP .LT. 1.0D+5 .AND. HTMP .GT. 1.0D0) THEN
          CFF1 = SQRT(G/HTMP)
          COUR = DTB2*CFF1*HTMP/DYR(4)
          ZSTAR = ( 0.5D0+COUR ) * (H0P(I-1,4) + H0P(I,4))*0.50D0 +          &
              ( 0.5D0-COUR ) * (H0P(I-1,3) + H0P(I,3))*0.50D0

          IF ( COUR .GT. Chalf ) THEN
            CFF2 = ( 1.0D0-Chalf/COUR )**2
            CFF3 = 0.50D0*(WORK(I-1,4) + WORK(I,4)) + &
                COUR * &
                0.50D0*(H0P(I-1,3) + H0P(I,3)) - &
                ( 1.0D0 + COUR ) * &
                0.50D0*(H0P(I-1,4) + H0P(I,4))
            ZSTAR = ZSTAR + CFF2 * CFF3
          ENDIF

          WKA(I,3,2) = 0.5D0 * &
              ( ( 1.0D0 - COUR ) * VBP(I,3) + COUR * VBP(I,4) + &
              VBTRN(I,1) - &
              CFF1 * ( ZSTAR - 0.50D0*(H0BTRN(I-1,1) + H0BTRN(I,1)) ) )
        ELSE
          WKA(I,3,2) = WKA(I,4,2)
        ENDIF
        WKA(I,2,2) = WKA(I,3,2) * VIV(I,2,1)
        WKA(I,1,2) = WKA(I,3,2) * VIV(I,1,1)
      ENDDO
      WKA(1,3,2) = WKA(2,3,2) * VIV(1,3,1)
      WKA(1,2,2) = WKA(2,2,2) * VIV(1,2,1)
      WKA(1,1,2) = WKA(2,1,2) * VIV(1,1,1)
    ENDIF
  ENDIF

  !$OMP PARALLEL DO PRIVATE (J,I)
  DO J=1, JMT
    DO I=1, IMT
      VB(I,J) = VB(I,J) * VIV(I,J,1)
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE OBC_2DVB

SUBROUTINE OBC_2DH0(JEWSN)

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, KM, JET, I, J, K
  use   pconst_mod, only: VIT, VIV, G, OHBT, DYR, OUX, DTB, DTB2, DZP, OTX, DYT, ix, iy
  use   dyn_mod,    only: UB, VB, H0, H0P, UBP, VBP
  use   work_mod

#if ( defined SPMD )
  use msg_mod
#endif

#if (defined OBCDT)
  use   openbc_mod
#endif

#if (defined OBTIDE)
  use   tide_mod, only: Ntide, OMGATIDE, Tintg_tide, V0NOD, UNOD, FNOD, &
      H0TIDEE, THTIDEE, H0TIDEW, THTIDEW, &
      H0TIDES, THTIDES, H0TIDEN, THTIDEN
#endif

  IMPLICIT NONE
  INTEGER :: JEWSN
  REAL(r8) :: EPSCF
  REAL(r8) :: Tau_obc
  REAL(r8) :: CFF1, CFF2
  REAL(r8) :: HTMP
  REAL(r8) :: UDTDX
  REAL(r8) :: DVDTTMP, DVDNMTMP, DVDTGTMP, DVDTG2TMP, CFTG, CFNM
  REAL(r8) :: CMIU

#if (defined OBTIDE)
  REAL(r8) :: H0tide
  REAL(r8) :: var_undef
  INTEGER :: KP
  REAL(r8) :: lon_tide

  var_undef = 1.0e+35
#endif

  EPSCF = 1.0D-20

  IF ( JEWSN .EQ. 1 ) THEN
    ! Eastern Open Boundary Condition
#if (defined OBTIDE)

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      lon_tide = 0.0D0
      H0tide = 0.0D0
      DO KP = 1, Ntide
        IF ( KP.LE.4 ) THEN
          IF ( H0TIDEE(J,KP) .NE. var_undef .and. &
              THTIDEE(J,KP) .NE. var_undef ) THEN
          H0tide = H0tide + H0TIDEE(J,KP)*FNOD(KP+4)* &
              DCOS( OMGATIDE(KP+4)*Tintg_tide*DTB + &
              lon_tide + V0NOD(KP+4)+UNOD(KP+4) + THTIDEE(J,KP) )
        ENDIF
      ELSE IF ( KP.LE.8 ) THEN
        IF ( H0TIDEE(J,KP) .NE. var_undef .and. &
            THTIDEE(J,KP) .NE. var_undef ) THEN
        H0tide = H0tide + H0TIDEE(J,KP)*FNOD(KP-4)* &
            DCOS( OMGATIDE(KP-4)*Tintg_tide*DTB + &
            lon_tide + V0NOD(KP-4)+UNOD(KP-4) + THTIDEE(J,KP) )
      ENDIF
    ELSE
      IF ( H0TIDEE(J,KP) .NE. var_undef .and. &
          THTIDEE(J,KP) .NE. var_undef ) THEN
      H0tide = H0tide + H0TIDEE(J,KP)*FNOD(KP)* &
          DCOS( OMGATIDE(KP)*Tintg_tide*DTB + &
          lon_tide + V0NOD(KP)+UNOD(KP) + THTIDEE(J,KP) )
    ENDIF
  ENDIF
  ENDDO

  H0(IMT-2,J) = ( H0BTRE(2,J)+H0tide )*VIT(IMT-2,J,1)
  H0(IMT-1,J)=H0(IMT-2,J)
  H0(IMT,J)=H0(IMT-2,J)
  ENDDO

#else

  IF ( ID_OBC_H0_E .EQ. 1 ) THEN
    ! The Zero Gradient Boundary Condition For H0 at Eastern Edge

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      H0(IMT-2,J) = H0(IMT-3,J) * VIT(IMT-2,J,1)
      H0(IMT-1,J) = H0(IMT-2,J) * VIT(IMT-1,J,1)
      H0(IMT,J) = H0(IMT-2,J) * VIT(IMT,J,1)

      IF (ID_H0_ZEROGRADI_NUDGING .EQ. 1) THEN
        H0(IMT-2,J) = H0(IMT-2,J) + DTB_OPENBC/TAU_2DH0*(H0BTRE(2,J) - H0PP(IMT-2,J))
      ENDIF
    ENDDO
  ELSE IF ( ID_OBC_H0_E .EQ. 2 ) THEN
    ! The Clamped Boundary Condition For H0 at Eastern Edge

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      H0(IMT-2,J) = H0BTRE(2,J) * VIT(IMT-2,J,1)
      H0(IMT-1,J) = H0(IMT-2,J) * VIT(IMT-1,J,1)
      H0(IMT,J) = H0(IMT-2,J) * VIT(IMT,J,1)
    ENDDO
  ELSE IF ( ID_OBC_H0_E .EQ. 3 ) THEN
    ! The Implicit Chapman Boundary Condition ( i.e. Gravity Wave  Radiation Condition ) For H0 at Eastern Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 2, JMT-1
      IF ( OHBT(IMT-3,J) .GT. 0.0D0 .AND. OHBT(IMT-3,J).LT.(1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(IMT-3,J) + H0P(IMT-3,J)
        CMIU = SQRT( G*HTMP )*DTB / ( 0.50D0/OUX(J-1) + 0.50D0/OUX(J) )
        H0(IMT-2,J) = ( 0.25D0*H0P(IMT-2,J-1)+0.50D0*H0P(IMT-2,J)+0.25D0*H0P(IMT-2,J+1) )/(1.0D0+CMIU) + &
            ( 0.25D0*H0(IMT-3,J-1)+0.50D0*H0(IMT-3,J)+0.25D0*H0(IMT-3,J+1))*CMIU/(1.0D0+CMIU)
      ELSE
        H0(IMT-2,J) = 0.0D0
      ENDIF

      H0(IMT-1,J) = H0(IMT-2,J) * VIT(IMT-1,J,1)
      H0(IMT,J)   = H0(IMT-2,J) * VIT(IMT,J,1)
    ENDDO
    H0(IMT-2,1) = H0(IMT-2,2) * VIT(IMT-2,1,1)
    H0(IMT-1,1) = H0(IMT-1,2) * VIT(IMT-1,1,1)
    H0(IMT,1) = H0(IMT,2) * VIT(IMT,1,1)
    H0(IMT-2,JMT) = H0(IMT-2,JMT-1) * VIT(IMT-2,JMT,1)
    H0(IMT-1,JMT) = H0(IMT-1,JMT-1) * VIT(IMT-1,JMT,1)
    H0(IMT,JMT) = H0(IMT,JMT-1) * VIT(IMT,JMT,1)
  ELSE IF ( ID_OBC_H0_E .EQ. 4 ) THEN
    ! The Explicit Chapman Boundary Condition For H0 at Eastern Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 2, JMT-1
      IF ( OHBT(IMT-3,J) .GT. 0.0D0 .AND. OHBT(IMT-3,J).LT.(1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(IMT-3,J) + H0P(IMT-3,J)
        CMIU = SQRT( G*HTMP )*DTB / (0.50D0/OUX(J-1) + 0.50D0/OUX(J) )
        CMIU = MIN(0.999D0, CMIU)
        H0(IMT-2,J) = (1.0D0-CMIU)*( 0.25D0*H0P(IMT-2,J-1)+0.50D0*H0P(IMT-2,J)+0.25D0*H0P(IMT-2,J+1) ) + &
            CMIU *( 0.25D0*H0P(IMT-3,J-1)+0.50D0*H0P(IMT-3,J)+0.25D0*H0P(IMT-3,J+1) )
      ELSE
        H0(IMT-2,J) = 0.0D0
      ENDIF

      H0(IMT-1,J) = H0(IMT-2,J) * VIT(IMT-1,J,1)
      H0(IMT,J)   = H0(IMT-2,J) * VIT(IMT,J,1)
    ENDDO
    H0(IMT-2,1) = H0(IMT-2,2) * VIT(IMT-2,1,1)
    H0(IMT-1,1) = H0(IMT-1,2) * VIT(IMT-1,1,1)
    H0(IMT,1) = H0(IMT,2) * VIT(IMT,1,1)
    H0(IMT-2,JMT) = H0(IMT-2,JMT-1) * VIT(IMT-2,JMT,1)
    H0(IMT-1,JMT) = H0(IMT-1,JMT-1) * VIT(IMT-1,JMT,1)
    H0(IMT,JMT) = H0(IMT,JMT-1) * VIT(IMT,JMT,1)
  ELSE IF ( ID_OBC_H0_E .EQ. 5 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 2, JMT-1
      DVDTTMP = H0P(IMT-3,J) - H0(IMT-3,J)
      DVDNMTMP = H0(IMT-3,J) - H0(IMT-4,J)
      DVDTG2TMP = H0P(IMT-3,J+1) - H0P(IMT-3,J-1)

      IF ( (DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
        DVDTTMP = 0.0D0
        Tau_obc = TAUIN_H0
      ELSE
        Tau_obc = TAUOUT_H0
      ENDIF

      IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
        DVDTGTMP = H0P(IMT-3,J) - H0P(IMT-3,J-1)
      ELSE
        DVDTGTMP = H0P(IMT-3,J+1) - H0P(IMT-3,J)
      ENDIF

      IF ( ID_OBC_2DRad_H0 .EQ. 0 ) THEN
        DVDTGTMP = 0.0D0
      ENDIF

      IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
        CFF1 = DVDNMTMP * DVDNMTMP + &
            DVDTGTMP * DVDTGTMP / ( DYR(J-1)*DYR(J-1) ) * &
            ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) * ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) )

        CFF2 = DVDTGTMP * DVDTGTMP + &
            DVDNMTMP * DVDNMTMP * DYR(J-1) * DYR(J-1) / &
            ( ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) * ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) )
      ELSE
        CFF1 = DVDNMTMP * DVDNMTMP + &
            DVDTGTMP * DVDTGTMP / ( DYR(J)*DYR(J) ) * &
            ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) * ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) )

        CFF2 = DVDTGTMP * DVDTGTMP + &
            DVDNMTMP * DVDNMTMP * DYR(J) * DYR(J) / &
            ( ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) * ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) )
      ENDIF

      CFF1 = MAX(CFF1,EPSCF)
      CFF2 = MAX(CFF2,EPSCF)
      CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
      CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

      IF ( ID_OBC_2DRad_NPO_H0 .EQ. 1 ) THEN
        CFTG = 0.0D0
      ENDIF

      CMIU = CFF1 + CFNM
      H0(IMT-2,J) = ( CFF1*H0P(IMT-2,J) + &
          CFNM*H0(IMT-3,J) -  &
          MAX(CFTG,0.0D0)*(CFF1/CFF2)*(H0P(IMT-2,J)-H0P(IMT-2,J-1))- &
          MIN(CFTG,0.0D0)*(CFF1/CFF2)*(H0P(IMT-2,J+1)-H0P(IMT-2,J)) )/CMIU

      H0(IMT-2,J) = H0(IMT-2,J) + DTB/Tau_obc * ( H0BTRE(2,J)-H0P(IMT-2,J) )*CFF1/CMIU
      H0(IMT-1,J) = H0(IMT-2,J) * VIT(IMT-1,J,1)
      H0(IMT,J) = H0(IMT-2,J) * VIT(IMT,J,1)
    ENDDO
    H0(IMT-2,1) = H0(IMT-2,2) * VIT(IMT-2,1,1)
    H0(IMT-1,1) = H0(IMT-1,2) * VIT(IMT-1,1,1)
    H0(IMT,1) = H0(IMT,2) * VIT(IMT,1,1)
    H0(IMT-2,JMT) = H0(IMT-2,JMT-1) * VIT(IMT-2,JMT,1)
    H0(IMT-1,JMT) = H0(IMT-1,JMT-1) * VIT(IMT-1,JMT,1)
    H0(IMT,JMT) = H0(IMT,JMT-1) * VIT(IMT,JMT,1)
  ELSE IF ( ID_OBC_H0_E .EQ. 6 ) THEN
    ! The Flather Boundary Condition For H0 at Eastern Edge
    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 2, JMT
      IF ( OHBT(IMT-3,J) .GT. 0.0D0 .AND. OHBT(IMT-3,J) .LE. (1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(IMT-3,J) + H0P(IMT-3,J)
        H0(IMT-3,J) = H0BTRE(1,J) + &
            ( 0.25D0*( UBP(IMT-3,J-1) * VIV(IMT-3,J-1,1) + UBP(IMT-3,J) * VIV(IMT-3,J,1) + &
            UBP(IMT-2,J-1) * VIV(IMT-2,J-1,1) + UBP(IMT-2,J) * VIV(IMT-2,J,1) ) - &
            0.25D0*( UBTRE(1,J-1) + UBTRE(1,J) + UBTRE(2,J-1) + UBTRE(2,J) ) ) / SQRT( G/HTMP )
      ELSE
        H0(IMT-3,J) = 0.0D0
      ENDIF

      H0(IMT-2,J) = H0(IMT-3,J) * VIT(IMT-2,J,1)
      H0(IMT-1,J) = H0(IMT-3,J) * VIT(IMT-1,J,1)
      H0(IMT,J) = H0(IMT-3,J) * VIT(IMT,J,1)
    ENDDO
    H0(IMT-3,1) = H0(IMT-3,2) * VIT(IMT-3,1,1)
    H0(IMT-2,1) = H0(IMT-2,2) * VIT(IMT-2,1,1)
    H0(IMT-1,1) = H0(IMT-1,2) * VIT(IMT-1,1,1)
    H0(IMT,1) = H0(IMT,2) * VIT(IMT,1,1)
  ELSE IF ( ID_OBC_H0_E .EQ. 7 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      IF (OHBT(IMT-3,J) .GT. 1.0d-10) THEN
        HTMP = 1.0D0/OHBT(IMT-3,J)
        CMIU = SQRT( G*HTMP )* DTB_OPENBC *OTX(J)
        H0(IMT-2,J) = ( H0PP(IMT-2,J) + CMIU*H0(IMT-3,J) )/( 1.0D0 + CMIU )
        ! CMIU*(H0(IMT-3,J+1) + H0(IMT-3,J)))/( 1.0D0 + CMIU )
        H0(IMT-1,J) = H0(IMT-2,J) * VIT(IMT-1,J,1)
        H0(IMT,J)   = H0(IMT-2,J) * VIT(IMT,J,1)
      ENDIF
    ENDDO
  ELSE IF ( ID_OBC_H0_E .EQ. 8 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT-1
      IF ( OHBT(IMT-2,J) .GT. 0.0D0) THEN
        HTMP = 1.0D0/OHBT(IMT-2,J)
        CMIU = SQRT( G*HTMP )*DTB / (0.50D0/OUX(J-1) + 0.50D0/OUX(J) )
        CMIU = MIN(0.999D0, CMIU)
        H0(IMT-2,J) = ( 1.0D0-CMIU ) * (H0(IMT-2,J) + H0(IMT-2,J+1)) * 0.5d0 &
            + CMIU*(H0(IMT-3,J) + H0(IMT-3,J+1)) * 0.5d0
      ENDIF

      H0(IMT-1,J) = H0(IMT-2,J) * VIT(IMT-1,J,1)
      H0(IMT,J)   = H0(IMT-2,J) * VIT(IMT,J,1)
    ENDDO

    H0(IMT-2,JMT) = H0(IMT-2,JMT-1) * VIT(IMT-2,JMT,1)
    H0(IMT-1,JMT) = H0(IMT-1,JMT-1) * VIT(IMT-1,JMT,1)
    H0(IMT,JMT) = H0(IMT,JMT-1) * VIT(IMT,JMT,1)
  ELSE IF ( ID_OBC_H0_E .EQ. 9 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT-1
      IF (OHBT(IMT-2,J) .GT. 1.0d-10) THEN
        HTMP = 1.0D0/OHBT(IMT-2,J)
        CMIU = SQRT( G*HTMP )*DTB2 *OTX(J)
        WORK(IMT-2,J) = 0.5d0 * ( H0P(IMT-2,J+1) + H0P(IMT-2,J) + &
            CMIU*(WORK(IMT-3,J+1) + WORK(IMT-3,J)))/( 1.0D0 + CMIU )
        ! CMIU*(WORK(IMT-3,J)))/( 1.0D0 + CMIU )

        WORK(IMT-1,J) = WORK(IMT-2,J) * VIT(IMT-1,J,1)
        WORK(IMT,J)   = WORK(IMT-2,J) * VIT(IMT,J,1)
      ELSE
        WORK(IMT-2,J) = WORK(IMT-2,J)
        WORK(IMT-1,J) = WORK(IMT-2,J) * VIT(IMT-1,J,1)
        WORK(IMT,J)   = WORK(IMT-2,J) * VIT(IMT,J,1)
      ENDIF
    ENDDO
  ENDIF

#endif
  ELSE IF ( JEWSN .EQ. 2 ) THEN
    ! Western Open Boundary Condition
#if (defined OBTIDE)

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      lon_tide = 0.0D0
      H0tide = 0.0D0
      DO KP = 1, Ntide
        IF ( KP.LE.4 ) THEN
          IF ( H0TIDEW(J,KP) .NE. var_undef .AND. &
              THTIDEW(J,KP) .NE. var_undef ) THEN

          H0tide = H0tide + H0TIDEW(J,KP)*FNOD(KP+4) * &
              DCOS( OMGATIDE(KP+4)*Tintg_tide*DTB + &
              lon_tide + V0NOD(KP+4)+UNOD(KP+4) + THTIDEW(J,KP) )
        ENDIF
      ELSE IF ( KP.LE.8 ) THEN
        IF ( H0TIDEW(J,KP) .NE. var_undef .AND. &
            THTIDEW(J,KP) .NE. var_undef ) THEN

        H0tide = H0tide + H0TIDEW(J,KP)*FNOD(KP-4) * &
            DCOS( OMGATIDE(KP-4)*Tintg_tide*DTB + &
            lon_tide + V0NOD(KP-4)+UNOD(KP-4) + THTIDEW(J,KP) )
      ENDIF
    ELSE
      IF ( H0TIDEW(J,KP) .NE. var_undef .AND. &
          THTIDEW(J,KP) .NE. var_undef ) THEN

      H0tide = H0tide + H0TIDEW(J,KP)*FNOD(KP) * &
          DCOS( OMGATIDE(KP)*Tintg_tide*DTB + &
          lon_tide + V0NOD(KP)+UNOD(KP) + THTIDEW(J,KP) )
    ENDIF
  ENDIF
  ENDDO

  H0(3,J) = ( H0BTRW(1,J)+H0tide )*VIT(3,J,1)
  H0(2,J) = H0(3,J)
  H0(1,J) = H0(3,J)
  ENDDO

#else

  IF ( ID_OBC_H0_W .EQ. 1 ) THEN
    ! The Zero Gradient Boundary Condition For H0 at Western Edge

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      H0(3,J) = H0(4,J) * VIT(3,J,1)
      H0(2,J) = H0(3,J) * VIT(2,J,1)
      H0(1,J) = H0(3,J) * VIT(1,J,1)

      IF (ID_H0_ZEROGRADI_NUDGING .EQ. 1) THEN
        H0(3,J) = H0(3,J) + DTB_OPENBC/TAU_2DH0*(H0BTRW(1,J) - H0PP(3,J))
      ENDIF
    ENDDO
  ELSE IF ( ID_OBC_H0_W .EQ. 2 ) THEN
    ! The Clamped Boundary Condition For H0 at Western Edge

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT
      H0(3,J) = H0BTRW(1,J) * VIT(3,J,1)
      H0(2,J) = H0(3,J) * VIT(2,J,1)
      H0(1,J) = H0(3,J) * VIT(1,J,1)
    ENDDO
  ELSE IF ( ID_OBC_H0_W .EQ. 3 ) THEN
    ! The Implicit Chapman Boundary Condition For H0 at Western Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 2, JMT-1
      IF ( OHBT(4,J) .GT. 0.0D0 .AND. OHBT(4,J).LT.(1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(4,J) + H0P(4,J)
        CMIU = SQRT( G*HTMP )*DTB / ( 0.50D0/OUX(J-1) + 0.50D0/OUX(J) )
        H0(3,J) = ( 0.25D0*H0P(3,J-1)+0.50D0*H0P(3,J)+0.25D0*H0P(3,J+1) )/(1.0D0+CMIU) + &
            ( 0.25D0*H0(4,J-1)+0.50D0*H0(4,J)+0.25D0*H0(4,J+1) )*CMIU/(1.0D0+CMIU)
      ELSE
        H0(3,J) = 0.0D0
      ENDIF

      H0(2,J) = H0(3,J) * VIT(2,J,1)
      H0(1,J) = H0(3,J) * VIT(1,J,1)
    ENDDO
    H0(3,1) = H0(3,2) * VIT(3,1,1)
    H0(2,1) = H0(2,2) * VIT(2,1,1)
    H0(1,1) = H0(1,2) * VIT(1,1,1)
    H0(3,JMT) = H0(3,JMT-1) * VIT(3,JMT,1)
    H0(2,JMT) = H0(2,JMT-1) * VIT(2,JMT,1)
    H0(1,JMT) = H0(1,JMT-1) * VIT(1,JMT,1)
  ELSE IF ( ID_OBC_H0_W .EQ. 4 ) THEN
    ! The Explicit Chapman Boundary Condition For H0 at Western Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 2, JMT-1
      IF ( OHBT(4,J) .GT. 0.0D0 .AND. OHBT(4,J).LT.(1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(4,J) + H0P(4,J)
        CMIU = SQRT( G*HTMP )*DTB / ( 0.50D0/OUX(J-1) + 0.50D0/OUX(J) )
        CMIU = MIN(0.999D0, CMIU)
        H0(3,J) = ( 0.25D0*H0P(3,J-1)+0.50D0*H0P(3,J)+0.25D0*H0P(3,J+1) )*(1.0D0-CMIU) + &
            ( 0.25D0*H0P(4,J-1)+0.50D0*H0P(4,J)+0.25D0*H0P(4,J+1) )*CMIU
      ELSE
        H0(3,J) = 0.0D0
      ENDIF
      H0(2,J) = H0(3,J) * VIT(2,J,1)
      H0(1,J) = H0(3,J) * VIT(1,J,1)
    ENDDO

    H0(3,1) = H0(3,2) * VIT(3,1,1)
    H0(2,1) = H0(2,2) * VIT(2,1,1)
    H0(1,1) = H0(1,2) * VIT(1,1,1)
    H0(3,JMT) = H0(3,JMT-1) * VIT(3,JMT,1)
    H0(2,JMT) = H0(2,JMT-1) * VIT(2,JMT,1)
    H0(1,JMT) = H0(1,JMT-1) * VIT(1,JMT,1)
  ELSE IF ( ID_OBC_H0_W .EQ. 5 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 2, JMT-1
      DVDTTMP = H0P(4,J) - H0(4,J)
      DVDNMTMP = H0(4,J) - H0(5,J)
      DVDTG2TMP = H0P(4,J+1) - H0P(4,J-1)

      IF ( (DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
        DVDTTMP = 0.0D0
        Tau_obc = TAUIN_H0
      ELSE
        Tau_obc = TAUOUT_H0
      ENDIF

      IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
        DVDTGTMP = H0P(4,J) - H0P(4,J-1)
      ELSE
        DVDTGTMP = H0P(4,J+1) - H0P(4,J)
      ENDIF

      IF ( ID_OBC_2DRad_H0 .EQ. 0 ) THEN
        DVDTGTMP = 0.0D0
      ENDIF

      IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
        CFF1 = DVDNMTMP * DVDNMTMP + &
            DVDTGTMP * DVDTGTMP / ( DYR(J-1)*DYR(J-1) ) * &
            ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) * ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) )

        CFF2 = DVDTGTMP * DVDTGTMP + &
            DVDNMTMP * DVDNMTMP * DYR(J-1) * DYR(J-1) / &
            ( ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) * ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) )
      ELSE
        CFF1 = DVDNMTMP * DVDNMTMP + &
            DVDTGTMP * DVDTGTMP / ( DYR(J)*DYR(J) ) * &
            ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) * ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) )

        CFF2 = DVDTGTMP * DVDTGTMP + &
            DVDNMTMP * DVDNMTMP * DYR(J) * DYR(J) / &
            ( ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) * ( 0.50D0/OUX(J)+0.50D0/OUX(J-1) ) )
      ENDIF

      CFF1 = MAX(CFF1,EPSCF)
      CFF2 = MAX(CFF2,EPSCF)
      CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
      CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

      IF ( ID_OBC_2DRad_NPO_H0 .EQ. 1 ) THEN
        CFTG = 0.0D0
      ENDIF

      CMIU = CFF1 + CFNM
      H0(3,J) = ( CFF1*H0P(3,J) + &
          CFNM*H0(4,J) -  &
          MAX(CFTG,0.0D0)*(CFF1/CFF2)*(H0P(3,J)-H0P(3,J-1))- &
          MIN(CFTG,0.0D0)*(CFF1/CFF2)*(H0P(3,J+1)-H0P(3,J)) )/CMIU

      H0(3,J) = H0(3,J) + DTB/Tau_obc * ( H0BTRW(1,J)-H0P(3,J) )*CFF1/CMIU
      H0(2,J) = H0(3,J) * VIT(2,J,1)
      H0(1,J) = H0(3,J) * VIT(1,J,1)
    ENDDO
    H0(3,1) = H0(3,2) * VIT(3,1,1)
    H0(2,1) = H0(2,2) * VIT(2,1,1)
    H0(1,1) = H0(1,2) * VIT(1,1,1)
    H0(3,JMT) = H0(3,JMT-1) * VIT(3,JMT,1)
    H0(2,JMT) = H0(2,JMT-1) * VIT(2,JMT,1)
    H0(1,JMT) = H0(1,JMT-1) * VIT(1,JMT,1)
  ELSE IF ( ID_OBC_H0_W .EQ. 6 ) THEN
    ! The Flather Boundary Condition For H0 at Western Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU+DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 2, JMT
      IF ( OHBT(4,J) .GT. 0.0D0 .AND. OHBT(4,J) .LE. (1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(4,J) + H0P(4,J)
        H0(4,J) = H0BTRW(2,J) - &
            ( 0.25D0*( UBP(4,J-1) * VIV(4,J-1,1) + UBP(4,J) * VIV(4,J,1) +  &
            UBP(5,J-1) * VIV(5,J-1,1) + UBP(5,J) * VIV(5,J,1) )- &
            0.25D0*( UBTRW(1,J-1)+UBTRW(1,J)+UBTRW(2,J-1)+UBTRW(2,J) ) )/SQRT(G/HTMP)
      ELSE
        H0(4,J) = 0.0D0
      ENDIF

      H0(3,J) = H0(4,J) * VIT(3,J,1)
      H0(2,J) = H0(3,J) * VIT(2,J,1)
      H0(1,J) = H0(3,J) * VIT(1,J,1)
    ENDDO
    H0(4,1) = H0(4,2) * VIT(4,1,1)
    H0(3,1) = H0(3,2) * VIT(3,1,1)
    H0(2,1) = H0(2,2) * VIT(2,1,1)
    H0(1,1) = H0(1,2) * VIT(1,1,1)
  ELSE IF ( ID_OBC_H0_W .EQ. 7 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1 , JMT
      IF (OHBT(4,J) .GT. 1.0d-10) THEN
        HTMP = 1.0D0/OHBT(4,J)
        CMIU = SQRT( G*HTMP )* DTB_OPENBC * OTX(J)
        H0(3,J) = ( H0PP(3,J) + CMIU*H0(4,J) )/( 1.0D0 + CMIU )
        ! CMIU*(H0(4,J+1) + H0(4,J)) )/( 1.0D0 + CMIU )
        H0(2,J) = H0(3,J) * VIT(2,J,1)
        H0(1,J) = H0(3,J) * VIT(1,J,1)
      ENDIF
    ENDDO
  ELSE IF ( ID_OBC_H0_W .EQ. 8 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1, JMT-1
      IF ( OHBT(3,J) .GT. 0.0D0  ) THEN
        HTMP = 1.0D0/OHBT(3,J)
        CMIU = SQRT( G*HTMP )*DTB / ( 0.50D0/OUX(J-1) + 0.50D0/OUX(J) )
        CMIU = MIN(0.999D0, CMIU)
        H0(3,J) = ( 1.0D0-CMIU ) * (H0(3,J) + H0(3,J+1)) * 0.5D0 &
            + CMIU*(H0(4,J) + H0(4,J+1)) *0.5D0
      ENDIF
      H0(2,J) = H0(3,J) * VIT(2,J,1)
      H0(1,J) = H0(3,J) * VIT(1,J,1)
    ENDDO

    H0(3,JMT) = H0(3,JMT-1) * VIT(3,JMT,1)
    H0(2,JMT) = H0(2,JMT-1) * VIT(2,JMT,1)
    H0(1,JMT) = H0(1,JMT-1) * VIT(1,JMT,1)
  ELSE IF ( ID_OBC_H0_W .EQ. 9 ) THEN
    !$OMP PARALLEL DO PRIVATE (J)
    DO J = 1 , JMT-1
      IF (OHBT(3,J) .GT. 1.0d-10) THEN
        HTMP = 1.0D0/OHBT(3,J)
        CMIU = SQRT( G*HTMP )*DTB2 * OTX(J)
        WORK(3,J) = 0.50d0 * ( H0P(3,J+1) + H0P(3,J) + &
            CMIU*(WORK(4,J+1) + WORK(4,J)) )/( 1.0D0 + CMIU )
        ! CMIU*(WORK(4,J)) )/( 1.0D0 + CMIU )

        WORK(2,J) = WORK(3,J) * VIT(2,J,1)
        WORK(1,J) = WORK(3,J) * VIT(1,J,1)
      ELSE
        WORK(3,J) = WORK(3,J)
        WORK(2,J) = WORK(3,J) * VIT(2,J,1)
        WORK(1,J) = WORK(3,J) * VIT(1,J,1)
      ENDIF
    ENDDO
  ENDIF

#endif
  ELSE IF ( JEWSN .EQ. 3 ) THEN
    ! Southern Open Boundary Condition
#if (defined OBTIDE)
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      lon_tide = 0.0D0
      H0tide = 0.0D0
      DO KP = 1, Ntide
        IF ( KP.LE.4 ) THEN
          IF ( H0TIDES(I,KP) .NE. var_undef .AND. &
              THTIDES(I,KP) .NE. var_undef ) THEN

          H0tide = H0tide + H0TIDES(I,KP) * FNOD(KP+4) *  &
              DCOS( OMGATIDE(KP+4)*Tintg_tide*DTB + &
              lon_tide + V0NOD(KP+4)+UNOD(KP+4) + THTIDES(I,KP) )
        ENDIF
      ELSE IF ( KP.LE.8 ) THEN
        IF ( H0TIDES(I,KP) .NE. var_undef .AND. &
            THTIDES(I,KP) .NE. var_undef ) THEN

        H0tide = H0tide + H0TIDES(I,KP) * FNOD(KP-4) *  &
            DCOS( OMGATIDE(KP-4)*Tintg_tide*DTB + &
            lon_tide + V0NOD(KP-4)+UNOD(KP-4) + THTIDES(I,KP) )
      ENDIF
    ELSE
      IF ( H0TIDES(I,KP) .NE. var_undef .AND. &
          THTIDES(I,KP) .NE. var_undef ) THEN

      H0tide = H0tide + H0TIDES(I,KP) * FNOD(KP) *  &
          DCOS( OMGATIDE(KP)*Tintg_tide*DTB + &
          lon_tide + V0NOD(KP)+UNOD(KP) + THTIDES(I,KP) )
    ENDIF
  ENDIF
  ENDDO

  H0(I,JET-2) = ( H0BTRS(I,2)+H0tide )*VIT(I,JET-2,1)
  H0(I,JET-1)=H0(I,JET-2)
  H0(I,JET)=H0(I,JET-2)
  ENDDO

#else

  IF ( ID_OBC_H0_S .EQ. 1 ) THEN
    ! The Zero Gradient Boundary Condition For H0 at Southern Edge

    !$OMP PARALLEL DO PRIVATE (I)
    DO I =1, IMT
      H0(I,JET-2) = H0(I,JET-3) * VIT(I,JET-2,1)
      H0(I,JET-1) = H0(I,JET-2) * VIT(I,JET-1,1)
      H0(I,JET) = H0(I,JET-2) * VIT(I,JET,1)

      IF (ID_H0_ZEROGRADI_NUDGING .EQ. 1) THEN
        H0(I,JET-2) = H0(I,JET-2) + DTB_OPENBC/TAU_2DH0*(H0BTRS(I,2) - H0PP(I,JET-2))
      ENDIF
    ENDDO
  ELSE IF ( ID_OBC_H0_S .EQ. 2 ) THEN
    ! The Clamped Gradient Boundary Condition For H0 at Southern Edge

    !$OMP PARALLEL DO PRIVATE (I)
    DO I =1, IMT
      H0(I,JET-2) = H0BTRS(I,2) * VIT(I,JET-2,1)
      H0(I,JET-1) = H0(I,JET-2) * VIT(I,JET-1,1)
      H0(I,JET) = H0(I,JET-2) * VIT(I,JET,1)
    ENDDO
  ELSE IF ( ID_OBC_H0_S .EQ. 3 ) THEN
    ! The Implicit Chapman Boundary Condition For H0 at Southern Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT-1
      IF ( OHBT(I,JET-3) .GT. 0.0D0 .AND. OHBT(I,JET-3).LT.(1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(I,JET-3) + H0P(I,JET-3)
        CMIU = SQRT( G*HTMP )*DTB / DYR(JET-3)
        H0(I,JET-2) = ( 0.25D0*H0P(I-1,JET-2)+0.50D0*H0P(I,JET-2)+0.25D0*H0P(I+1,JET-2) )/(1.0D0+CMIU) + &
            ( 0.25D0*H0(I-1,JET-3)+0.50D0*H0(I,JET-3)+0.25D0*H0(I+1,JET-3) )*CMIU/(1.0D0+CMIU)
      ELSE
        H0(I,JET-2) = 0.0D0
      ENDIF

      H0(I,JET-1) = H0(I,JET-2) * VIT(I,JET-1,1)
      H0(I,JET)   = H0(I,JET-2) * VIT(I,JET,1)
    ENDDO
    H0(IMT,JET-2) = H0(IMT-1,JET-2) * VIT(IMT,JET-2,1)
    H0(IMT,JET-1) = H0(IMT-1,JET-1) * VIT(IMT,JET-1,1)
    H0(IMT,JET) = H0(IMT-1,JET) * VIT(IMT,JET,1)
    H0(1,JET-2) = H0(2,JET-2) * VIT(1,JET-2,1)
    H0(1,JET-1) = H0(2,JET-1) * VIT(1,JET-1,1)
    H0(1,JET) = H0(2,JET) * VIT(1,JET,1)
  ELSE IF ( ID_OBC_H0_S .EQ. 4 ) THEN
    ! The Explicit Chapman Boundary Condition For H0 at Southern Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT-1
      IF ( OHBT(I,JET-3) .GT. 0.0D0 .AND. OHBT(I,JET-3).LT.(1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(I,JET-3) + H0P(I,JET-3)
        CMIU = SQRT( G*HTMP )*DTB / DYR(JET-3)
        CMIU = MIN(0.999D0, CMIU)
        H0(I,JET-2) = ( 0.25D0*H0P(I-1,JET-2)+0.50D0*H0P(I,JET-2)+0.25D0*H0P(I+1,JET-2) )*(1.0D0-CMIU) + &
            ( 0.25D0*H0P(I-1,JET-3)+0.50D0*H0P(I,JET-3)+0.25D0*H0P(I+1,JET-3) )*CMIU
      ELSE
        H0(I,JET-2) = 0.0D0
      ENDIF

      H0(I,JET-1) = H0(I,JET-2) * VIT(I,JET-1,1)
      H0(I,JET)   = H0(I,JET-2) * VIT(I,JET,1)
    ENDDO

    H0(IMT,JET-2) = H0(IMT-1,JET-2) * VIT(IMT,JET-2,1)
    H0(IMT,JET-1) = H0(IMT-1,JET-1) * VIT(IMT,JET-1,1)
    H0(IMT,JET) = H0(IMT-1,JET) * VIT(IMT,JET,1)
    H0(1,JET-2) = H0(2,JET-2) * VIT(1,JET-2,1)
    H0(1,JET-1) = H0(2,JET-1) * VIT(1,JET-1,1)
    H0(1,JET) = H0(2,JET) * VIT(1,JET,1)
  ELSE IF ( ID_OBC_H0_S .EQ. 5 ) THEN
    ! The Implicit Upstream Radiation Condition For H0 at Southern Edge

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT-1
      DVDTTMP = H0P(I,JET-3) - H0(I,JET-3)
      DVDNMTMP = H0(I,JET-3) - H0(I,JET-4)
      DVDTG2TMP = H0P(I+1,JET-3) - H0P(I-1,JET-3)

      IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
        DVDTTMP = 0.0D0
        Tau_obc = TAUIN_H0
      ELSE
        Tau_obc = TAUOUT_H0
      ENDIF

      IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
        DVDTGTMP = H0P(I,JET-3) - H0P(I-1,JET-3)
      ELSE
        DVDTGTMP = H0P(I+1,JET-3) - H0P(I,JET-3)
      ENDIF

      IF ( ID_OBC_2DRad_H0 .EQ. 0 ) THEN
        DVDTGTMP = 0.0D0
      ENDIF

      CFF1 = DVDNMTMP * DVDNMTMP + &
          DVDTGTMP * DVDTGTMP * DYR(JET-4) * DYR(JET-4) / &
          ( ( 0.50D0/OUX(JET-4)+0.50D0/OUX(JET-3) ) * ( 0.50D0/OUX(JET-4)+0.50D0/OUX(JET-3) ) )

      CFF2 = DVDTGTMP * DVDTGTMP + &
          DVDNMTMP * DVDNMTMP / ( DYR(JET-4)*DYR(JET-4) ) * &
          ( 0.50D0/OUX(JET-4)+0.50D0/OUX(JET-3) ) * ( 0.50D0/OUX(JET-4)+0.50D0/OUX(JET-3) )

      CFF1 = MAX(CFF1,EPSCF)
      CFF2 = MAX(CFF2,EPSCF)
      CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
      CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

      IF ( ID_OBC_2DRad_NPO_H0 .EQ. 1 ) THEN
        CFTG = 0.0D0
      ENDIF

      CMIU = CFF1 + CFNM
      H0(I,JET-2) = ( CFF1*H0P(I,JET-2) + &
          CFNM*H0(I,JET-3) -  &
          MAX(CFTG,0.0D0)*(CFF1/CFF2)*(H0P(I,JET-2)-H0P(I-1,JET-2))- &
          MIN(CFTG,0.0D0)*(CFF1/CFF2)*(H0P(I+1,JET-2)-H0P(I,JET-2)) ) / CMIU

      H0(I,JET-2) = H0(I,JET-2) + DTB/Tau_obc * ( H0BTRS(I,2)-H0P(I,JET-2) ) * CFF1/CMIU
      H0(I,JET-1) = H0(I,JET-2) * VIT(I,JET-1,1)
      H0(I,JET)   = H0(I,JET-2) * VIT(I,JET,1)
    ENDDO
    H0(1,JET-2) = H0(2,JET-2) * VIT(1,JET-2,1)
    H0(1,JET-1) = H0(2,JET-1) * VIT(1,JET-1,1)
    H0(1,JET)   = H0(2,JET) * VIT(1,JET,1)
    H0(IMT,JET-2) = H0(IMT-1,JET-2) * VIT(IMT,JET-2,1)
    H0(IMT,JET-1) = H0(IMT-1,JET-1) * VIT(IMT,JET-1,1)
    H0(IMT,JET)   = H0(IMT-1,JET) * VIT(IMT,JMT,1)
  ELSE IF ( ID_OBC_H0_S .EQ. 6 ) THEN
    ! The Flather Boundary Condition For H0 at Southern Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT-1
      IF ( OHBT(I,JET-3) .GT. 0.0D0 .AND. OHBT(I,JET-3) .LE. (1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(I,JET-3) + H0P(I,JET-3)
        H0(I,JET-3) = H0BTRS(I,1) + &
            ( 0.25D0 * ( UBP(I,JET-4)*VIV(I,JET-4,1)+UBP(I+1,JET-4)*VIV(I+1,JET-4,1) + &
            UBP(I,JET-3)*VIV(I,JET-3,1)+UBP(I+1,JET-3)*VIV(I+1,JET-3,1) ) - &
            0.25D0 * ( UBTRS(I,1)+UBTRS(I+1,1)+UBTRS(I,2)+UBTRS(I+1,2) ) ) / SQRT(G/HTMP)
      ELSE
        H0(I,JET-3) = 0.0D0
      ENDIF

      H0(I,JET-2) = H0(I,JET-3) * VIT(I,JET-2,1)
      H0(I,JET-1) = H0(I,JET-2) * VIT(I,JET-1,1)
      H0(I,JET) = H0(I,JET-2) * VIT(I,JET,1)
    ENDDO

    H0(IMT,JET-3) = H0(IMT-1,JET-3) * VIT(IMT,JET-3,1)
    H0(IMT,JET-2) = H0(IMT-1,JET-2) * VIT(IMT,JET-2,1)
    H0(IMT,JET-1) = H0(IMT-1,JET-1) * VIT(IMT,JET-1,1)
    H0(IMT,JET) = H0(IMT-1,JET) * VIT(IMT,JMT,1)
  ELSE IF ( ID_OBC_H0_S .EQ. 7 ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      IF (OHBT(I,JET-3) .GT. 1.0d-10) THEN
        HTMP = 1.0D0/OHBT(I,JET-3)
        CMIU = SQRT( G*HTMP )* DTB_OPENBC / DYT(JET-3)
        H0(I,JET-2) = ( H0PP(I,JET-2) + CMIU*H0(I,JET-3) )/( 1.0D0 + CMIU )
        ! CMIU*(H0(I-1,JET-3) + H0(I,JET-3)) )/( 1.0D0 + CMIU )
        H0(I,JET-1) = H0(I,JET-2) * VIT(I,JET-1,1)
        H0(I,JET)   = H0(I,JET-2) * VIT(I,JET,1)
      ENDIF
    ENDDO
  ELSE IF ( ID_OBC_H0_S .EQ. 8 ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT
      IF ( OHBT(I,JET-2) .GT. 0.0D0  ) THEN
        HTMP = 1.0D0/OHBT(I,JET-2)
        CMIU = SQRT( G*HTMP )*DTB / DYR(JET-3)
        CMIU = MIN(0.999D0, CMIU)
        H0(I,JET-2) = ( 1.0D0-CMIU ) * (H0(I,JET-2) + H0(I-1,JET-2)) * 0.5D0  &
            + CMIU * (H0(I,JET-3) + H0(I-1,JET-3)) * 0.5D0
      ENDIF

      H0(I,JET-1) = H0(I,JET-2) * VIT(I,JET-1,1)
      H0(I,JET)   = H0(I,JET-2) * VIT(I,JET,1)
    ENDDO

    H0(1,JET-2) = H0(2,JET-2) * VIT(1,JET-2,1)
    H0(1,JET-1) = H0(2,JET-1) * VIT(1,JET-1,1)
    H0(1,JET) = H0(2,JET) * VIT(1,JET,1)
  ELSE IF ( ID_OBC_H0_S .EQ. 9 ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT
      IF (OHBT(I,JET-2) .GT. 1.0d-10) THEN
        HTMP = 1.0D0/OHBT(I,JET-2)
        CMIU = SQRT( G*HTMP )*DTB2 / DYT(JET-3)
        WORK(I,JET-2) = 0.50d0*( H0P(I-1,JET-2) + H0P(I,JET-2) + &
            CMIU*(WORK(I-1,JET-3) + WORK(I,JET-3)) )/( 1.0D0 + CMIU )
        ! CMIU*(WORK(I,JET-3)) )/( 1.0D0 + CMIU )

        WORK(I,JET-1) = WORK(I,JET-2) * VIT(I,JET-1,1)
        WORK(I,JET)   = WORK(I,JET-2) * VIT(I,JET,1)
      ELSE
        WORK(I,JET-2) = WORK(I,JET-2)
        WORK(I,JET-1) = WORK(I,JET-2) * VIT(I,JET-1,1)
        WORK(I,JET)   = WORK(I,JET-2) * VIT(I,JET,1)
      ENDIF
    ENDDO
  ENDIF

#endif
  ELSE IF ( JEWSN .EQ. 4 ) THEN
    ! Northern Open Boundary Condition
#if (defined OBTIDE)
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      lon_tide = 0.0D0
      H0tide = 0.0D0

      DO KP = 1, Ntide
        IF ( KP.LE.4 ) THEN
          IF ( H0TIDEN(I,KP) .NE. var_undef .AND. &
              THTIDEN(I,KP) .NE. var_undef ) THEN

          H0tide = H0tide + H0TIDEN(I,KP) * FNOD(KP+4) *  &
              DCOS( OMGATIDE(KP+4)*Tintg_tide*DTB + &
              lon_tide + V0NOD(KP+4) + UNOD(KP+4) + THTIDEN(I,KP) )
        ENDIF
      ELSE IF ( KP.LE.8 ) THEN
        IF ( H0TIDEN(I,KP) .NE. var_undef .AND. &
            THTIDEN(I,KP) .NE. var_undef ) THEN

        H0tide = H0tide + H0TIDEN(I,KP) * FNOD(KP-4) *  &
            DCOS( OMGATIDE(KP-4)*Tintg_tide*DTB + &
            lon_tide + V0NOD(KP-4) + UNOD(KP-4) + THTIDEN(I,KP) )
      ENDIF
    ELSE
      IF ( H0TIDEN(I,KP) .NE. var_undef .AND. &
          THTIDEN(I,KP) .NE. var_undef ) THEN

      H0tide = H0tide + H0TIDEN(I,KP) * FNOD(KP) *  &
          DCOS( OMGATIDE(KP)*Tintg_tide*DTB + &
          lon_tide + V0NOD(KP) + UNOD(KP) + THTIDEN(I,KP) )
    ENDIF
  ENDIF
  ENDDO

  H0(I,3) = ( H0BTRN(I,1)+H0tide )*VIT(I,3,1)
  H0(I,2) = H0(I,3)
  H0(I,1) = H0(I,3)
  ENDDO

#else

  IF ( ID_OBC_H0_N .EQ. 1 ) THEN
    ! The Zero Gradient Boundary Condition For H0 at Northern Edge

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      H0(I,3) = H0(I,4) * VIT(I,3,1)
      H0(I,2) = H0(I,3) * VIT(I,2,1)
      H0(I,1) = H0(I,3) * VIT(I,1,1)

      IF (ID_H0_ZEROGRADI_NUDGING .EQ. 1) THEN
        H0(I,3) = H0(I,3) + DTB_OPENBC/TAU_2DH0*(H0BTRN(I,1) - H0PP(I,3))
      ENDIF
    ENDDO
  ELSE IF ( ID_OBC_H0_N .EQ. 2 ) THEN
    ! The Clamped Boundary Condition For H0 at Northern Edge

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      H0(I,3) = H0BTRN(I,1) * VIT(I,3,1)
      H0(I,2) = H0(I,3) * VIT(I,2,1)
      H0(I,1) = H0(I,3) * VIT(I,1,1)
    ENDDO
  ELSE IF ( ID_OBC_H0_N .EQ. 3 ) THEN
    ! The Implicit Chapman Boundary Condition For H0 at Northern Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT-1
      IF ( OHBT(I,4) .GT. 0.0D0 .AND. OHBT(I,4).LT.(1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(I,4) + H0P(I,4)
        CMIU = SQRT( G*HTMP )*DTB / DYR(3)
        H0(I,3) = ( 0.25D0*H0P(I-1,3)+0.50D0*H0P(I,3)+0.25D0*H0P(I+1,3) )/(1.0D0+CMIU) + &
            ( 0.25D0*H0(I-1,4)+0.50D0*H0(I,4)+0.25D0*H0(I+1,4) )*CMIU/(1.0D0+CMIU)
      ELSE
        H0(I,3) = 0.0D0
      ENDIF

      H0(I,2) = H0(I,3) * VIT(I,2,1)
      H0(I,1) = H0(I,3) * VIT(I,1,1)
    ENDDO

    H0(IMT,3) = H0(IMT-1,3) * VIT(IMT,3,1)
    H0(IMT,2) = H0(IMT-1,2) * VIT(IMT,2,1)
    H0(IMT,1) = H0(IMT-1,1) * VIT(IMT,1,1)
    H0(1,3) = H0(2,3) * VIT(1,3,1)
    H0(1,2) = H0(2,2) * VIT(1,2,1)
    H0(1,1) = H0(2,1) * VIT(1,1,1)
  ELSE IF ( ID_OBC_H0_N .EQ. 4 ) THEN
    ! The EXplicit Chapman Boundary Condition For H0 at Northern Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT-1
      IF ( OHBT(I,4) .GT. 0.0D0 .AND. OHBT(I,4).LT.(1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(I,4) + H0P(I,4)
        CMIU = SQRT( G*HTMP )*DTB / DYR(3)
        CMIU = MIN(0.999D0, CMIU)
        H0(I,3) = ( 0.25D0*H0P(I-1,3)+0.50D0*H0P(I,3)+0.25D0*H0P(I+1,3) )*(1.0D0-CMIU) + &
            ( 0.25D0*H0P(I-1,4)+0.50D0*H0P(I,4)+0.25D0*H0P(I+1,4) )*CMIU
      ELSE
        H0(I,3) = 0.0D0
      ENDIF

      H0(I,2) = H0(I,3) * VIT(I,2,1)
      H0(I,1) = H0(I,3) * VIT(I,1,1)
    ENDDO

    H0(IMT,3) = H0(IMT-1,3) * VIT(IMT,3,1)
    H0(IMT,2) = H0(IMT-1,2) * VIT(IMT,2,1)
    H0(IMT,1) = H0(IMT-1,1) * VIT(IMT,1,1)
    H0(1,3) = H0(2,3) * VIT(1,3,1)
    H0(1,2) = H0(2,2) * VIT(1,2,1)
    H0(1,1) = H0(2,1) * VIT(1,1,1)
  ELSE IF ( ID_OBC_H0_N .EQ. 5 ) THEN
    ! The Implicit Upstream Radiation Condition For H0 at Northern Edge

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT-1
      DVDTTMP = H0P(I,4) - H0(I,4)
      DVDNMTMP = H0(I,4) - H0(I,5)
      DVDTG2TMP = H0P(I+1,4) - H0P(I-1,4)

      IF ( ( DVDTTMP*DVDNMTMP) .LT. 0.0D0 ) THEN
        DVDTTMP = 0.0D0
        Tau_obc = TAUIN_H0
      ELSE
        Tau_obc = TAUOUT_H0
      ENDIF

      IF ( (DVDTTMP*DVDTG2TMP) .GT. 0.0 ) THEN
        DVDTGTMP = H0P(I,4) - H0P(I-1,4)
      ELSE
        DVDTGTMP = H0P(I+1,4) - H0P(I,4)
      ENDIF

      IF ( ID_OBC_2DRad_H0 .EQ. 0 ) THEN
        DVDTGTMP = 0.0D0
      ENDIF

      CFF1 = DVDNMTMP * DVDNMTMP + &
          DVDTGTMP * DVDTGTMP * DYR(4) * DYR(4) / &
          ( ( 0.50D0/OUX(4)+0.50D0/OUX(3) ) * ( 0.50D0/OUX(4)+0.50D0/OUX(3) ) )

      CFF2 = DVDTGTMP * DVDTGTMP + &
          DVDNMTMP * DVDNMTMP / ( DYR(4)*DYR(4) ) * &
          ( 0.50D0/OUX(4)+0.50D0/OUX(3) ) * ( 0.50D0/OUX(4)+0.50D0/OUX(3) )

      CFF1 = MAX(CFF1,EPSCF)
      CFF2 = MAX(CFF2,EPSCF)
      CFNM = MIN(CFF1,MAX(DVDTTMP*DVDNMTMP,-CFF1))
      CFTG = MIN(CFF2,MAX(DVDTTMP*DVDTGTMP,-CFF2))

      IF ( ID_OBC_2DRad_NPO_H0 .EQ. 1 ) THEN
        CFTG = 0.0D0
      ENDIF

      CMIU = CFF1 + CFNM
      H0(I,3) = ( CFF1*H0P(I,3) + &
          CFNM*H0(I,4) -  &
          MAX(CFTG,0.0D0)*(CFF1/CFF2)*(H0P(I,3)-H0P(I-1,3))- &
          MIN(CFTG,0.0D0)*(CFF1/CFF2)*(H0P(I+1,3)-H0P(I,3)) ) / CMIU

      H0(I,3) = H0(I,3) + DTB/Tau_obc * ( H0BTRN(I,1)-H0P(I,3) ) * CFF1/CMIU
      H0(I,2) = H0(I,3) * VIT(I,2,1)
      H0(I,1) = H0(I,3) * VIT(I,1,1)
    ENDDO

    H0(1,3) = H0(2,3) * VIT(1,3,1)
    H0(1,2) = H0(2,2) * VIT(1,2,1)
    H0(1,1) = H0(2,1) * VIT(1,1,1)
    H0(IMT,3) = H0(IMT-1,3) * VIT(IMT,3,1)
    H0(IMT,2) = H0(IMT-1,2) * VIT(IMT,2,1)
    H0(IMT,1) = H0(IMT-1,1) * VIT(IMT,1,1)
  ELSE IF ( ID_OBC_H0_N .EQ. 6 ) THEN
    ! The Flather Boundary Condition For H0 at Northern Edge

    CMIU = 0.0D0
    DO K = 1, 5
      CMIU = CMIU + DZP(K)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT-1
      IF ( OHBT(I,4) .GT. 0.0D0 .AND. OHBT(I,4) .LE. (1.0D0/CMIU) ) THEN
        HTMP = 1.0D0/OHBT(I,4) + H0P(I,4)
        H0(I,4) = H0BTRN(I,2) - &
            ( 0.25D0 *( UBP(I,3)*VIV(I,3,1) + UBP(I+1,3)*VIV(I+1,3,1) + &
            UBP(I,4)*VIV(I,4,1) + UBP(I+1,4)*VIV(I+1,4,1) ) - &
            0.25D0 *( UBTRN(I,1)+UBTRN(I+1,1)+UBTRN(I,2)+UBTRN(I+1,2) ) )/SQRT(G/HTMP)
      ELSE
        H0(I,4) = 0.0D0
      ENDIF

      H0(I,3) = H0(I,4) * VIT(I,3,1)
      H0(I,2) = H0(I,3) * VIT(I,2,1)
      H0(I,1) = H0(I,3) * VIT(I,1,1)
    ENDDO

    H0(IMT,4) = H0(IMT-1,4) * VIT(IMT,4,1)
    H0(IMT,3) = H0(IMT-1,3) * VIT(IMT,3,1)
    H0(IMT,2) = H0(IMT-1,2) * VIT(IMT,2,1)
    H0(IMT,1) = H0(IMT-1,1) * VIT(IMT,1,1)
  ELSE IF ( ID_OBC_H0_N .EQ. 7 ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 1, IMT
      IF (OHBT(I,4) .GT. 1.0d-10) THEN
        HTMP = 1.0D0/OHBT(I,4)
        CMIU = SQRT( G*HTMP )* DTB_OPENBC / DYT(3)
        H0(I,3) = ( H0PP(I,3) + CMIU*H0(I,4) )/( 1.0D0 + CMIU )
        ! CMIU*(H0(I-1,4) + H0(I,4)) )/( 1.0D0 + CMIU )
        H0(I,2) = H0(I,3) * VIT(I,2,1)
        H0(I,1) = H0(I,3) * VIT(I,1,1)
      ENDIF
    ENDDO
  ELSE IF ( ID_OBC_H0_N .EQ. 8 ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT
      IF ( OHBT(I,3) .GT. 0.0D0) THEN
        HTMP = 1.0D0/OHBT(I,3)
        CMIU = SQRT( G*HTMP )*DTB / DYR(3)
        CMIU = MIN(0.999D0, CMIU)
        H0(I,3) = ( 1.0D0-CMIU ) * (H0(I,3) + H0(I-1,3)) * 0.5D0  &
            + CMIU * (H0(I,4) + H0(I-1,4)) * 0.5D0
      ENDIF

      H0(I,2) = H0(I,3) * VIT(I,2,1)
      H0(I,1) = H0(I,3) * VIT(I,1,1)
    ENDDO

    H0(1,3) = H0(2,3) * VIT(1,3,1)
    H0(1,2) = H0(2,2) * VIT(1,2,1)
    H0(1,1) = H0(2,1) * VIT(1,1,1)
  ELSE IF ( ID_OBC_H0_N .EQ. 9 ) THEN
    !$OMP PARALLEL DO PRIVATE (I)
    DO I = 2, IMT
      IF (OHBT(I,3) .GT. 1.0d-10) THEN
        HTMP = 1.0D0/OHBT(I,3)
        CMIU = SQRT( G*HTMP )*DTB2 / DYT(3)
        WORK(I,3) = 0.50d0*( H0P(I-1,3) + H0P(I,3) + &
            CMIU*(WORK(I-1,4) + WORK(I,4)) )/( 1.0D0 + CMIU )
        ! CMIU*(WORK(I,4)) )/( 1.0D0 + CMIU )

        WORK(I,2) = WORK(I,3) * VIT(I,2,1)
        WORK(I,1) = WORK(I,3) * VIT(I,1,1)
      ELSE
        WORK(I,3) = WORK(I,3)
        WORK(I,2) = WORK(I,3) * VIT(I,2,1)
        WORK(I,1) = WORK(I,3) * VIT(I,1,1)
      ENDIF
    ENDDO
  ENDIF

#endif
  ENDIF

  DO J=1, JMT
    DO I=1, IMT
      H0(I,J)=H0(I,J)*VIT(I,J,1)
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE OBC_2DH0

subroutine corners_2d()

  ! handle points at corners of 2d variables

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, KM, I, J, K
  use   dyn_mod,    only: H0, UB, VB
  use   pconst_mod, only: VIT, VIV, ny_proc, nx_proc, ix, iy
  implicit none

  ! northeast
  if (ix == nx_proc-1 .and. iy == 0) then
    h0(imt-2,3) = 0.5d0*(h0(imt-3,3)+h0(imt-2,4))*vit(imt-2,3,1)
    ub(imt-2,3) = 0.5d0*(ub(imt-3,3)+ub(imt-2,4))*viv(imt-2,3,1)
    vb(imt-2,3) = 0.5d0*(vb(imt-3,3)+vb(imt-2,4))*viv(imt-2,3,1)
  endif

  ! northwest
  if (ix == 0 .and. iy == 0) then
    h0(3,3) = 0.5d0*(h0(4,3)+h0(3,4))*vit(3,3,1)
    ub(4,3) = 0.5d0*(ub(5,3)+ub(4,4))*viv(4,3,1)
    vb(4,3) = 0.5d0*(vb(5,3)+vb(4,4))*viv(4,3,1)
  endif

  ! southeast
  if (ix == nx_proc-1 .and. iy == ny_proc-1) then
    h0(imt-2,jmt-2) = 0.5d0*(h0(imt-3,jmt-2)+h0(imt-2,jmt-3))*vit(imt-2,jmt-2,1)
    ub(imt-2,jmt-3) = 0.5d0*(ub(imt-3,jmt-3)+ub(imt-2,jmt-4))*viv(imt-2,jmt-3,1)
    vb(imt-2,jmt-3) = 0.5d0*(vb(imt-3,jmt-3)+vb(imt-2,jmt-4))*viv(imt-2,jmt-3,1)
  endif

  ! southwest
  if (ix == 0 .and. iy == ny_proc-1) then
    h0(3,jmt-2) = 0.5d0*(h0(4,jmt-2)+h0(3,jmt-3))*vit(3,jmt-2,1)
    ub(4,jmt-3) = 0.5d0*(ub(5,jmt-3)+ub(4,jmt-4))*viv(4,jmt-3,1)
    vb(4,jmt-3) = 0.5d0*(vb(5,jmt-3)+vb(4,jmt-4))*viv(4,jmt-3,1)
  endif

end subroutine corners_2d

subroutine corners_3d()

  ! handle points at corners of 3d variables

#include <def-undef.h>
  use   precision_mod
  use   param_mod,  only: IMT, JMT, KM, I, J, K
  use   dyn_mod,    only: U, V
  use   pconst_mod, only: VIT, VIV, ny_proc, nx_proc, ix, iy
  implicit none

  ! northeast
  if (ix == nx_proc-1 .and. iy == 0) then
    do k = 1,km
      u(imt-2,3,k) = 0.5d0*(u(imt-3,3,k)+u(imt-2,4,k))*viv(imt-2,3,1)
      v(imt-2,3,k) = 0.5d0*(v(imt-3,3,k)+v(imt-2,4,k))*viv(imt-2,3,1)
    enddo
  endif

  ! northwest
  if (ix == 0 .and. iy == 0) then
    do k = 1,km
      u(4,3,k) = 0.5d0*(u(5,3,k)+u(4,4,k))*viv(4,3,1)
      v(4,3,k) = 0.5d0*(v(5,3,k)+v(4,4,k))*viv(4,3,1)
    enddo
  endif

  ! southeast
  if (ix == nx_proc-1 .and. iy == ny_proc-1) then
    do k = 1,km
      u(imt-2,jmt-3,k) = 0.5d0*(u(imt-3,jmt-3,k)+u(imt-2,jmt-4,k))*viv(imt-2,jmt-3,1)
      v(imt-2,jmt-3,k) = 0.5d0*(v(imt-3,jmt-3,k)+v(imt-2,jmt-4,k))*viv(imt-2,jmt-3,1)
    enddo
  endif

  ! southwest
  if (ix == 0 .and. iy == ny_proc-1) then
    do k = 1,km
      u(4,jmt-3,k) = 0.5d0*(u(5,jmt-3,k)+u(4,jmt-4,k))*viv(4,jmt-3,1)
      v(4,jmt-3,k) = 0.5d0*(v(5,jmt-3,k)+v(4,jmt-4,k))*viv(4,jmt-3,1)
    enddo
  endif

end subroutine corners_3d

SUBROUTINE output_east_obc_dat

#include <def-undef.h>
#if (defined OBCDT)
  use   param_mod,  only: JMT, KM, nx_proc, ny_proc, ierr
  use   pconst_mod, only: ix, iy, j_global
  use   openbc_mod, only: H0BTRE, UBTRE, VBTRE, UBCLE, VBCLE
  use   msg_mod
#endif
  IMPLICIT NONE

#if (defined OBCDT)
  integer, parameter :: unit_h0_obc  = 71
  integer, parameter :: unit_ub_obc  = 72
  integer, parameter :: unit_vb_obc  = 73
  integer, parameter :: unit_ubc_obc = 74
  integer, parameter :: unit_vbc_obc = 75
  integer :: pe_obc, jb_obc, kb_obc
  logical, save :: dumped_obc = .false.
  character(len=8) :: status_obc

  if (dumped_obc) return

  do pe_obc = 0, ny_proc - 1
    if (ix .eq. nx_proc - 1 .and. iy .eq. pe_obc) then
      if (pe_obc .eq. 0) then
        status_obc = 'replace'
      else
        status_obc = 'old'
      endif

      open(unit_h0_obc, file='h0btre_east.dat', form='formatted', &
          status=status_obc, position='append')
      do jb_obc = 1, JMT
        write(unit_h0_obc, '(E24.16,1X,E24.16,1X,I8,1X,I8)') &
            H0BTRE(1,jb_obc), H0BTRE(2,jb_obc), iy, j_global(jb_obc)
      enddo
      close(unit_h0_obc)

      open(unit_ub_obc, file='ubtre_east.dat', form='formatted', &
          status=status_obc, position='append')
      do jb_obc = 1, JMT
        write(unit_ub_obc, '(E24.16,1X,E24.16,1X,I8,1X,I8)') &
            UBTRE(1,jb_obc), UBTRE(2,jb_obc), iy, j_global(jb_obc)
      enddo
      close(unit_ub_obc)

      open(unit_vb_obc, file='vbtre_east.dat', form='formatted', &
          status=status_obc, position='append')
      do jb_obc = 1, JMT
        write(unit_vb_obc, '(E24.16,1X,E24.16,1X,I8,1X,I8)') &
            VBTRE(1,jb_obc), VBTRE(2,jb_obc), iy, j_global(jb_obc)
      enddo
      close(unit_vb_obc)

      open(unit_ubc_obc, file='ubcle_east.dat', form='formatted', &
          status=status_obc, position='append')
      do jb_obc = 1, JMT
        do kb_obc = 1, KM
          write(unit_ubc_obc, '(E24.16,1X,E24.16,1X,I8,1X,I8)') &
              UBCLE(1,jb_obc,kb_obc), UBCLE(2,jb_obc,kb_obc), &
              iy, j_global(jb_obc)
        enddo
      enddo
      close(unit_ubc_obc)

      open(unit_vbc_obc, file='vbcle_east.dat', form='formatted', &
          status=status_obc, position='append')
      do jb_obc = 1, JMT
        do kb_obc = 1, KM
          write(unit_vbc_obc, '(E24.16,1X,E24.16,1X,I8,1X,I8)') &
              VBCLE(1,jb_obc,kb_obc), VBCLE(2,jb_obc,kb_obc), &
              iy, j_global(jb_obc)
        enddo
      enddo
      close(unit_vbc_obc)
    endif

    call mpi_barrier(mpi_comm_ocn, ierr)
  enddo

  dumped_obc = .true.
#endif

  RETURN
END SUBROUTINE output_east_obc_dat
