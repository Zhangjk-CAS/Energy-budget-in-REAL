!
!

      SUBROUTINE OBC_FLUX

#include <def-undef.h>
use   precision_mod

use   param_mod,  only: IMT, JMT, JET, I, J, mytid, n_proc, nx_proc, ny_proc
use   pconst_mod, only: VIT, VIV, OHBU, ix, iy, DYR, OUX !, OHBT, DYT, OTX   
use   dyn_mod,    only: UB, VB, H0

#if ( defined SPMD )
use msg_mod
#endif

#ifdef OBCDT
use   openbc_mod

      IMPLICIT NONE

      REAL(r8) :: MY_CRSS
      REAL(r8) :: MY_FLUX
      REAL(r8) :: CFF

      INTEGER :: NSUB
      INTEGER :: SIZE, STEP, IERR
      REAL(r8) :: BUFF(2)

      MY_CRSS = 0.0D0
      MY_FLUX = 0.0D0

      IF ( ix.eq.0 ) THEN

!$OMP PARALLEL DO PRIVATE (J)

           DO J = 1, JMT-1
!           DO J = 1, JMT
!
!!               IF( VIV(4,J,1) .GT. 0 ) THEN
!!                   CFF = ( H0(4,J) + DZPH(4,J) )*DYR(J)
!!                   MY_CRSS = MY_CRSS + CFF
!!                   MY_FLUX = MY_FLUX + CFF * 0.5D0*( UB(4,J-1)+UB(4,J) ) 
!!               ENDIF
!!
!              IF( VIT(4,J,1) .GT. 0 ) THEN
!                  CFF = ( H0(4,J) + 1.0D0/OHBT(4,J) ) * DYT(J)
!                  MY_CRSS = MY_CRSS + CFF
!                  MY_FLUX = MY_FLUX + CFF * &
!                            0.25 * ( UB(4,J-1)*VIV(4,J-1,1) + UB(5,J-1)*VIV(5,J-1,1) + &
!                                     UB(4,J) * VIV(4,J,1) + UB(5,J) * VIV(5,J,1) )
!              ENDIF
!
              IF( VIV(4,J,1) .GT. 0 ) THEN
                  CFF = ( ( H0(3,J)*VIT(3,J,1) + H0(3,J+1)*VIT(3,J+1,1) +  &
                            H0(4,J)*VIT(4,J,1) + H0(4,J+1)*VIT(4,J+1,1) )/ &
                          ( VIT(3,J,1)+VIT(3,J+1,1)+VIT(4,J,1)+VIT(4,J+1,1) )+ &
                          1.0D0/OHBU(4,J) ) * DYR(J)

                  MY_CRSS = MY_CRSS + CFF
                  MY_FLUX = MY_FLUX + CFF * UB(4,J)*VIV(4,J,1)
              ENDIF

           ENDDO
      ENDIF

      IF ( ix.eq.(nx_proc-1) ) THEN

!$OMP PARALLEL DO PRIVATE (J)

           DO J = 1, JMT-1
!           DO J = 1, JMT
!
!!             IF( VIV(IMT-2,J,1) .GT. 0 ) THEN
!!                 CFF = ( H0(IMT-2,J)+DZPH(IMT-2,J) )*DYR(J)
!!                 MY_CRSS = MY_CRSS + CFF                
!!                 MY_FLUX = MY_FLUX + CFF * (0.5D0*(UB(IMT-2,J-1)+UB(IMT-2,J)))
!!             ENDIF
!!
!              IF( VIT(IMT-2,J,1) .GT. 0 ) THEN
!                  CFF = ( H0(IMT-2,J) + 1.0D0/OHBT(IMT-2,J) ) *DYT(J)
!                  MY_CRSS = MY_CRSS + CFF
!                  MY_FLUX = MY_FLUX + CFF * 0.5D0 * ( UB(IMT-2,J-1) * VIV(IMT-2,J-1,1) + &
!                                                      UB(IMT-2,J) * VIV(IMT-2,J,1) )
!              ENDIF
!
              IF( VIV(IMT-2,J,1) .GT. 0 ) THEN
                  CFF = ( ( H0(IMT-3,J)*VIT(IMT-3,J,1)+H0(IMT-2,J)*VIT(IMT-2,J,1)+ &
                            H0(IMT-3,J+1)*VIT(IMT-3,J+1,1)+H0(IMT-2,J+1)*VIT(IMT-2,J+1,1) )/ &
                          ( VIT(IMT-3,J,1)+VIT(IMT-2,J,1)+VIT(IMT-3,J+1,1)+VIT(IMT-2,J+1,1) )+ &
                          1.0D0/OHBU(IMT-2,J) ) *DYR(J)

                  MY_CRSS = MY_CRSS + CFF
                  MY_FLUX = MY_FLUX - CFF* UB(IMT-2,J)
              ENDIF

           ENDDO

      ENDIF

      IF ( iy.eq.0 ) THEN

!$OMP PARALLEL DO PRIVATE (I)

           DO I = 2, IMT
!           DO I = 1, IMT
!
!!              IF ( VIV(I,3,1) .GT. 0 ) THEN
!!                   CFF = (H0(I,3)+DZPH(I,3))/OUX(3)
!!                   MY_CRSS = MY_CRSS + CFF
!!                   MY_FLUX = MY_FLUX + CFF * 0.5D0*(VB(I,3)+VB(I+1,3))
!!              ENDIF
!!
!              IF ( VIT(I,3,1) .GT. 0 ) THEN
!                   CFF = ( H0(I,3) + 1.0D0/OHBT(I,3) )/OTX(3)
!                   MY_CRSS = MY_CRSS + CFF
!                   MY_FLUX = MY_FLUX + CFF * 0.5D0*( VB(I,3)*VIV(I,3,1)+VB(I+1,3)*VIV(I+1,3,1) )
!              ENDIF
!
              IF ( VIV(I,3,1) .GT. 0 ) THEN
                   CFF = ( ( H0(I-1,3)*VIT(I-1,3,1)+H0(I,3)*VIT(I,3,1) + &
                             H0(I-1,4)*VIT(I-1,4,1)+H0(I,4)*VIT(I,4,1) )/ &
                           ( VIT(I-1,3,1)+VIT(I,3,1)+VIT(I-1,4,1)+VIT(I,4,1) ) + &
                           1.0D0/OHBU(I,3) )/OUX(3)

                   MY_CRSS = MY_CRSS + CFF
                   MY_FLUX = MY_FLUX + CFF * VB(I,3)*VIV(I,3,1)
              ENDIF

           ENDDO

      ENDIF

      IF ( iy.eq.(ny_proc-1) ) THEN

!$OMP PARALLEL DO PRIVATE (I)

           DO I = 2, IMT
!           DO I = 1, IMT
!!
!!              IF ( VIV(I,JET-3,1) .GT. 0 ) THEN
!!                   CFF = (H0(I,JET-3)+DZPH(I,JET-3))/OUX(JET-3)
!!                   MY_CRSS = MY_CRSS + CFF
!!                   MY_FLUX = MY_FLUX + CFF * 0.5D0*(VBP(I,JET-3)+VBP(I+1,JET-3))
!!              ENDIF
!!
!              IF ( VIT(I,JET-3,1) .GT. 0 ) THEN
!                   CFF = ( H0(I,JET-3) + 1.0D0/OHBT(I,JET-3) )/OTX(JET-3)
!                   MY_CRSS = MY_CRSS + CFF
!                   MY_FLUX = MY_FLUX + CFF * &
!                             0.25*( VB(I,JET-4)*VIV(I,JET-4,1) + VB(I+1,JET-4)*VIV(I+1,JET-4,1) + &
!                                    VB(I,JET-3)*VIV(I,JET-3,1) + VB(I+1,JET-3)*VIV(I+1,JET-3,1) )
!              ENDIF
!
              IF ( VIV(I,JET-3,1) .GT. 0 ) THEN
                   CFF = ( ( H0(I-1,JET-3)*VIT(I-1,JET-3,1)+H0(I,JET-3)*VIT(I,JET-3,1) + &
                             H0(I-1,JET-2)*VIT(I-1,JET-2,1)+H0(I,JET-2)*VIT(I,JET-2,1) )/ &
                           ( VIT(I-1,JET-3,1)+VIT(I,JET-3,1)+VIT(I-1,JET-2,1)+VIT(I,JET-2,1) ) + &
                           1.0D0/OHBU(I,JET-3) )/OUX(JET-3)

                   MY_CRSS = MY_CRSS + CFF
                   MY_FLUX = MY_FLUX - CFF * VB(I,JET-3)*VIV(I,JET-3,1)
              ENDIF

           ENDDO
      ENDIF

      IF ( n_proc == 1 ) THEN
           NSUB = 1
      ELSE
           NSUB = nx_proc*ny_proc
      ENDIF

!$OMP CRITICAL (obc_flux)

      IF (BC_COUNT.EQ.0) THEN
          BC_CRSS = 0.0D0
          BC_FLUX = 0.0D0
      ENDIF

      BC_COUNT = BC_COUNT + 1
      BC_CRSS = BC_CRSS + MY_CRSS
      BC_FLUX = BC_FLUX + MY_FLUX

      IF ( BC_COUNT.EQ.NSUB ) THEN
           BC_COUNT = 0

#ifdef SPMD
           SIZE = n_proc
10         STEP = (SIZE+1)/2
           IF ( mytid.GE.STEP .AND. mytid.LT.SIZE ) THEN
                BUFF(1) = BC_CRSS
                BUFF(2) = BC_FLUX
                CALL mpi_send( BUFF, 2, MPI_PR, mytid-STEP, tag_2d, mpi_comm_ocn,IERR)

           ELSE IF ( mytid.LT.(SIZE-STEP) ) THEN             
                CALL mpi_recv( BUFF, 2, MPI_PR, mytid+STEP, tag_2d,mpi_comm_ocn,status,IERR)
                BC_CRSS = BC_CRSS + BUFF(1)
                BC_FLUX = BC_FLUX + BUFF(2)
           ENDIF

           SIZE = STEP
           IF ( SIZE.GT.1 ) GO TO 10

           BUFF(1) = BC_CRSS
           BUFF(2) = BC_FLUX

           CALL mpi_bcast( BUFF, 2, MPI_PR, 0, mpi_comm_ocn, IERR)
           BC_CRSS = BUFF(1)
           BC_FLUX = BUFF(2)
#endif
           
           UBAR_XS = BC_FLUX/BC_CRSS
      ENDIF

!$OMP END CRITICAL (obc_flux)  

#endif

      RETURN

      END SUBROUTINE OBC_FLUX
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
     
      SUBROUTINE CONSERVE_MASS

#include <def-undef.h>
use   precision_mod

use   param_mod,  only: IMT, JMT, I, J, nx_proc, ny_proc
use   pconst_mod, only: VIV, ix, iy
use   dyn_mod,    only: UB, VB

#ifdef OBCDT

use   openbc_mod, only: UBAR_XS


      IMPLICIT NONE

      IF ( ix.eq.(nx_proc-1) ) THEN

           DO J = 1, JMT
              UB(IMT-2,J) = ( UB(IMT-2,J) + UBAR_XS ) * VIV(IMT-2,J,1)
              UB(IMT-1,J) = ( UB(IMT-1,J) + UBAR_XS ) * VIV(IMT-1,J,1)
              UB(IMT,J) = ( UB(IMT,J) + UBAR_XS ) * VIV(IMT,J,1)
           ENDDO
      ENDIF

      IF ( ix.eq.0 ) THEN
           DO J = 1, JMT
              UB(4,J) = ( UB(4,J) - UBAR_XS ) * VIV(4,J,1)
              UB(3,J) = ( UB(3,J) - UBAR_XS ) * VIV(3,J,1)
              UB(2,J) = ( UB(2,J) - UBAR_XS ) * VIV(2,J,1)
              UB(1,J) = ( UB(1,J) - UBAR_XS ) * VIV(1,J,1)
           ENDDO
      ENDIF 

      IF ( iy.eq.(ny_proc-1) ) THEN
           DO I = 1, IMT
              VB(I,JMT-3) = ( VB(I,JMT-3) + UBAR_XS )*VIV(I,JMT-3,1)
              VB(I,JMT-2) = ( VB(I,JMT-2) + UBAR_XS )*VIV(I,JMT-2,1)
              VB(I,JMT-1) = ( VB(I,JMT-1) + UBAR_XS )*VIV(I,JMT-1,1)
              VB(I,JMT) = ( VB(I,JMT) + UBAR_XS )*VIV(I,JMT,1)
           ENDDO
      ENDIF

      IF ( iy.eq.0 ) THEN
           DO I = 1, IMT
              VB(I,3) = ( VB(I,3) - UBAR_XS )*VIV(I,3,1)
              VB(I,2) = ( VB(I,2) - UBAR_XS )*VIV(I,2,1)
              VB(I,1) = ( VB(I,1) - UBAR_XS )*VIV(I,1,1)
           ENDDO 
      ENDIF

#endif

      RETURN
      END SUBROUTINE CONSERVE_MASS
