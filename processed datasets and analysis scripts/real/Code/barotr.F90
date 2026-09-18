!  CVS: $Id: barotr.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE BAROTR
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
#ifdef SPMD
use msg_mod
#endif
!
! Dfs Added For Use In Open Boundary Condition
#if (defined OBCDT)
use openbc_mod
#endif
! Dfs End
!
! DFS 2021-05
!
#if (defined OBTIDE)
use tide_mod, only: Tintg_tide
#endif

#if (defined POTTIDE)
use tide_mod, only: Tintg_tide, BETA_TIDE, ptide_elevat, &
                    CD_COASTAL, CD_PACIFIC, CD_TIDE
#endif
!
! DFS
!
      IMPLICIT NONE
!
!#if (defined OBCDT)
!      INTEGER :: LMN
!#endif
!
!#if(defined OBTIDE)
!      REAL(r8) :: VIV_tide(IMT,JMT)
!      REAL(r8) :: VIT_tide(IMT,JMT)
!#endif
!
!
      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2
      real(r8), dimension(imt,jmt) :: hdtk
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
!

       fil_lat1=65.0D0
       fil_lat2=65.0D0

!#endif
!
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
      wka=0
      work=0
!---------------------------------------------------------------------
!     EULER BACKWARD SCHEME IS USED FOR THE FIRST STEP OF EVERY MONTH
!     IEB=0: LEAP-FROG SCHEME; IEB=1: EULER BACKWARD SCHEME
!---------------------------------------------------------------------
      IEB = 0 ; IEB_LOOP=0

      IF (ISB == 0)  THEN
         IEB = 1 ; IEB_LOOP=1
      END IF
!
!
       baro_loop : DO NC = 1,NBB+IEB_LOOP

!      if (mytid.eq.0) print*,nc
!      call energy
!
! DFS 2021-12
!
#if ( defined POTTIDE )

             call POT_TIDE

!#ifdef SPMD
!       call exchange_2d(elevat)
!#endif

#endif
!
! DFS 2021-12
!

      if (IEB==1.or.ISB>1) then

!---------------------------------------------------------------------
!     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
!---------------------------------------------------------------------

#if ( defined SMAG1)
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,11)= UBP (I,J)
               WKA (I,J,12)= VBP (I,J)
            END DO
         END DO
         CALL SMAG2 (1)
#if (defined SMAG_FZ)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,5)= VIV (I,J,1)* (0.5* OUX (J)* (WKA (I +1,J,7) &
                            - WKA (I -1,J,7)) &
               - R2E (J)* WKA (I,J +1,8) + R2F (J)* WKA (I,J -1,8))
               WKA (I,J,6)= VIV (I,J,1)* (0.5* OUX (J)* (WKA (I +1,J,9) &
                            - WKA (I -1,J,9)) &
               - R3E (J)* WKA (I,J +1,10) + R3F (J)* WKA (I,J -1,10)    &
                            + R4E (J)* WKA (I,J,7))
            END DO
         END DO

!-new
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,5)= VIV (I,J,1)* (0.5* OUX (J)* (WKA (I +1,J,7) &
                            - WKA (I -1,J,7)) &
               - R2E (J)* WKA (I,J +1,8) + R2F (J)* WKA (I,J -1,8))
               WKA (I,J,6)= VIV (I,J,1)* (0.5* OUX (J)* (WKA (I +1,J,9) &
                            - WKA (I -1,J,9)) &
               - R3E (J)* WKA (I,J +1,10) + R3F (J)* WKA (I,J -1,10)    &
                            + R4E (J)* WKA (I,J,7))
            END DO
         END DO
#endif
#else
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,7)= AM3(I,J,1)*VIV(I,J,1)*(R1D(J)*(UBP(I,J+1)-UBP(I,J))&
                           -R1C (J)*(UBP(I,J)-UBP(I,J-1)) &
                           +SOUX (J)* (UBP(I+1,J)-2.0*UBP (I,J)+UBP(I-1,J)) &
                           +CV1(J)*UBP(I,J)+CV2(J)*(VBP(I+1,J)-VBP(I-1,J)))
               WKA (I,J,8)= AM3 (I,J,1)* VIV(I,J,1)*(R1D(J)*(VBP(I,J+1)-VBP (I,J)) &
                           -R1C(J)*(VBP(I,J)-VBP(I,J-1))&
                           +SOUX (J)* (VBP(I +1,J) -2.0* VBP (I,J)+VBP (I -1,J)) &
                           + CV1 (J)*VBP(I,J)-CV2(J)*(UBP(I+1,J)-UBP(I-1,J)))
            END DO
         END DO
!
!  DFS 2024-02-29
!
!       The original code as following
!
!!
!!
!      if (nx_proc==1) then
!         do j=1,jsm,jem
!            wka(1,j,7)=wka(imm,j,7)
!            wka(imt,j,7)=wka(2,j,7)
!            wka(1,j,8)=wka(imm,j,8)
!            wka(imt,j,8)=wka(2,j,8)
!         end do
!     end if
!!
!!
!
#if ( defined OBCDT )
      if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (J)
         do J = JSM, JEM
            wka(1,J,7) = wka(2,J,7)
            wka(IMT,J,7) = wka(IMM,J,7)

            wka(1,J,8) = wka(2,J,8)
            wka(IMT,J,8) = wka(IMM,J,8)
         end do

!$OMP PARALLEL DO PRIVATE (I)
         do I = 1, IMT
            wka(I,1,7) = wka(I,2,7)
            wka(I,JMT,7) = wka(I,JEM,7)

            wka(I,1,8) = wka(I,2,8)
            wka(I,JMT,8) = wka(I,JEM,8)
         end do
      end if
#else
      if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (J)
         do j=1,jsm,jem
            wka(1,j,7)=wka(imm,j,7)
            wka(imt,j,7)=wka(2,j,7)
            wka(1,j,8)=wka(imm,j,8)
            wka(imt,j,8)=wka(2,j,8)
         end do
      end if
#endif
!
!
!

#ifdef SPMD
         call exchange_2d(wka(1,1,7),1,1)
         call exchange_2d(wka(1,1,8),1,1)
#endif

!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,5)= VIV (I,J,1)* ( R1D (J)* (WKA (I,J +1,7)     &
                           -WKA(I,J,7))-R1C (J)*(WKA(I,J,7)-WKA(I,J-1,7))+  &
                           SOUX (J)* (WKA (I+1,J,7) -2.0*               &
               WKA (I,J,7) + WKA (I -1,J,7)) + CV1 (J)* WKA (I,J,7)     &
                            + CV2 (J)*(WKA (I+1,J,8)-WKA (I-1,J,8)))
               WKA (I,J,6)= VIV (I,J,1)* ( R1D (J)* (WKA (I,J +1,8)     &
                            - WKA (I,J,8)) &
               - R1C (J)* (WKA (I,J,8) - WKA (I,J -1,8)) + SOUX (J)* (  &
                           WKA (I +1,J,8) -2.0*&
               WKA (I,J,8) + WKA (I -1,J,8)) + CV1 (J)* WKA (I,J,8)     &
                            - CV2 (J)* (WKA (I +1,J,7)-WKA (I-1,J,7)))
            END DO
         END DO
!
!
!  DFS 2024-02-29
!
!     The original code as following
!
!!
!!
!     if (nx_proc == 1) then
!         do j=jsm,jem
!            wka(1,j,5)=wka(imm,j,5)
!            wka(imt,j,5)=wka(2,j,5)
!            wka(1,j,6)=wka(imm,j,6)
!            wka(imt,j,6)=wka(2,j,6)
!         end do
!     end if
!!
!!
!
#if ( defined OBCDT ) 
     if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (J)
         do J = JSM, JEM
            wka(1,J,5) = wka(2,J,5)
            wka(IMT,J,5) = wka(IMM,J,5)
            wka(1,J,6) = wka(2,J,6)
            wka(IMT,J,6) = wka(IMM,J,6)
         end do

!$OMP PARALLEL DO PRIVATE (I)
         do I = 1, IMT
            wka(I,1,5) = wka(I,2,5)
            wka(I,JMT,5) = wka(I,JEM,5)
            wka(I,1,6) = wka(I,2,6)
            wka(I,JMT,6) = wka(I,JEM,6)
         end do
     end if
#else

     if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (J)
         do j=jsm,jem
            wka(1,j,5)=wka(imm,j,5)
            wka(imt,j,5)=wka(2,j,5)
            wka(1,j,6)=wka(imm,j,6)
            wka(imt,j,6)=wka(2,j,6)
         end do
     end if
#endif
!
!  DFS 2024-02-29
!

!!!!!!
#else
!$OMP PARALLEL DO PRIVATE (J,I,am_viv)
         DO J = JSM,JEM
            DO I = 2,IMM
               am_viv=AM3 (I,J,1)* VIV (I,J,1)
               WKA (I,J,5)=am_viv*(R1D(J)*(UBP(I,J+1)-UBP(I,J))-R1C(J)*  &
                           (UBP(I,J)-UBP(I,J-1))+SOUX(J)*(UBP(I+1,J)-    &
                           2.0*UBP(I,J)+UBP(I-1,J))+CV1(J)*UBP(I,J)+     &
                           CV2(J)*(VBP(I+1,J)-VBP(I-1,J)))
               WKA (I,J,6)=am_viv*(R1D(J)*(VBP(I,J+1)-VBP(I,J))-R1C(J)*  &
                           (VBP(I,J)-VBP(I,J-1))+SOUX(J)*(VBP(I+1,J)-    &
                           2.0* VBP (I,J)+VBP(I-1,J))+CV1(J)*VBP (I,J)-  &
                           CV2(J)*(UBP(I+1,J)-UBP(I-1,J)))
            END DO
         END DO

#endif
#endif

            IF (mod(isb,36)  == 1 ) THEN
!$OMP PARALLEL DO PRIVATE (J,I)
               DO J = JSM,JEM
               DO I = 2,IMM
                  DLUB (I,J)= DLUB (I,J) + WKA (I,J,5)
                  DLVB (I,J)= DLVB (I,J) + WKA (I,J,6)
               END DO
               END DO
            END IF
!
         END IF

!---------------------------------------------------------------------
!     + (g'-1)g*dH/dr
!---------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (J,I,gstar)
         DO J = JSM,JEM
            DO I = 2,IMM
               gstar=(WGP (I,J) -1.0)*G *0.5
!
! DFS 2021-12
!
!lhl20100801
!
#if ( defined POTTIDE )

               WKA (I,J,1) = WKA (I,J,5) &
               + gstar*OUX(J)*(H0(I,J) - H0(I-1,J) + H0(I,J+1) - H0(I-1,J+1))*BETA_TIDE &
!               + gstar*OUX(J)*(H0(I,J) - H0(I-1,J) + H0(I,J+1) - H0(I-1,J+1))*1.052 &
!               + gstar*OUX(J)*(H0(I,J) - H0(I-1,J) + H0(I,J+1) - H0(I-1,J+1))*0.052 &
               - gstar*OUX(J)*(ptide_elevat(I,J) - ptide_elevat(I-1,J) + ptide_elevat(I,J+1) &
               - ptide_elevat(I-1,J+1) )

               WKA (I,J,2) = WKA (I,J,6) &
               + gstar*OUY(J)*(H0(I,J+1) - H0(I,J) + H0(I-1,J+1) - H0(I-1,J))*BETA_TIDE &
!               + gstar*OUY(J)*(H0(I,J+1) - H0(I,J) + H0(I-1,J+1) - H0(I-1,J))*1.052 &
!               + gstar*OUY(J)*(H0(I,J+1) - H0(I,J) + H0(I-1,J+1) - H0(I-1,J))*0.052 &
               - gstar*OUY(J)*( ptide_elevat(I,J+1) - ptide_elevat(I,J) + ptide_elevat(I-1,J+1) &
               - ptide_elevat(I-1,J) )

#else

               WKA (I,J,1) = WKA (I,J,5) &
               + gstar*OUX(J)*(H0(I,J) - H0(I-1,J) + H0(I,J+1) - H0(I-1,J+1))

               WKA (I,J,2) = WKA (I,J,6) &
               + gstar*OUY(J)*(H0(I,J+1) - H0(I,J) + H0(I-1,J+1) - H0(I-1,J))

#endif
!
!!lhl20100801
!
! DFS 2021-12
!
            END DO
         END DO


!---------------------------------------------------------------------
!     COMPUTING H0 AT U/V POINTS
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WORK (I,J) = 0.25* (H0 (I,J) + H0 (I -1,J) + H0 (I,J +1) &
                            + H0 (I -1,J +1))
            END DO

!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!     DO J=2,JMM
         END DO

!
!  DFS 2024-02-29
!
!     The oringinal code as following
!
!!
!!
!     if (nx_proc == 1) then
!         do j=jsm,jem
!            WORK (1,J) = WORK (IMM,J)
!            WORK (IMT,J) = WORK (2,J)
!         end do
!      end if
!!
!!
!
#if ( defined OBCDT ) 
     if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (J)
         do J = JSM, JEM
            WORK (1,J) = WORK (2,J)
            WORK (IMT,J) = WORK (IMM,J)
         end do

!$OMP PARALLEL DO PRIVATE (I)
         do I = 1, IMT
            WORK (I,1) = WORK (I,2)
            WORK (I,JMT) = WORK (I,JEM)
         end do
     end if
#else
     if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (J)
         do j=jsm,jem
            WORK (1,J) = WORK (IMM,J)
            WORK (IMT,J) = WORK (2,J)
         end do
     end if
#endif
!
!  DFS 2024-02-29
!

#ifdef SPMD
       call exchange_2d(work,1,1)
#endif

!Yu
!---------------------------------------------------------------------
!     COMPUTING DU & DV
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,1)= VIV (I,J,1)* ( WKA (I,J,1) + DLUB (I,J)     &
                              - FF (J)* VBP (I,J) + &
               PAX (I,J) + PXB (I,J) - WORK (I,J)* WHX (I,J) )
               WKA (I,J,2)= VIV (I,J,1)* ( WKA (I,J,2) + DLVB (I,J)     &
                              + FF (J)* UBP (I,J) + &
               PAY (I,J) + PYB (I,J) - WORK (I,J)* WHY (I,J) )
            END DO
         END DO

!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------

         IF (ISB == 0) THEN
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = JSM,JEM
               DO I = 2,IMM
                  WKA (I,J,3)= EBEA (J)* WKA (I,J,1) - EBEB (J)* WKA (I,J,2)
                  WKA (I,J,4)= EBEA (J)* WKA (I,J,2) + EBEB (J)* WKA (I,J,1)
               END DO
            END DO
         ELSE
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = JSM,JEM
               DO I = 2,IMM
                  WKA (I,J,3)= EBLA (J)* WKA (I,J,1) - EBLB (J)* WKA (I,J,2)
                  WKA (I,J,4)= EBLA (J)* WKA (I,J,2) + EBLB (J)* WKA (I,J,1)
               END DO
            END DO
         END IF


!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
        if (nx_proc /= 1) then
#ifdef SPMD
        call exchange_2d(wka(1,1,3),1,0)
        call exchange_2d(wka(1,1,4),1,0)
#endif
         else
#if ( defined OBCDT ) 
!$OMP PARALLEL DO PRIVATE (J)
         DO J = JSM,JEM
            WKA (1,J,3) = WKA (2,J,3)
            WKA (IMT,J,3) = WKA (IMM,J,3)
            WKA (1,J,4) = WKA (2,J,4)
            WKA (IMT,J,4) = WKA (IMM,J,4)
         END DO
#else
!$OMP PARALLEL DO PRIVATE (J)
         DO J = JSM,JEM
            WKA (1,J,3) = WKA (IMM,J,3)
            WKA (IMT,J,3) = WKA (2,J,3)
            WKA (1,J,4) = WKA (IMM,J,4)
            WKA (IMT,J,4) = WKA (2,J,4)
         END DO
#endif
         end if

!---------------------------------------------------------------------
!     COMPUTING DH0
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,1)= UB (I,J)* (DZPH (I,J) + WORK (I,J))
               WKA (I,J,2)= VB (I,J)* (DZPH (I,J) + WORK (I,J))
            END DO
         END DO

! !ADD BY ZHANGJK 20250331
! IF (ID_OBC_SET_MASS .EQ. 1) THEN 
! CALL set_mass()
! ENDIF 

!-----------------------------------------------------------------
! calculate hdtk
! add by zhangjk 20250113

!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               hdtk (i,j) = AM3(I,J,1)*VIT(I,J,1)*(R1D(J)*(H0P(I,J+1)-H0P(I,J))&
                           -R1C (J)*(H0P(I,J)-H0P(I,J-1)) &
                           +SOTX (J)* (H0P(I+1,J)-2.0*H0P (I,J)+H0P(I-1,J)) &
                           +CV1(J)*H0P(I,J)+CV2(J)*(H0P(I+1,J)-H0P(I-1,J)))
               ! WKA (I,J,7)= AM3(I,J,1)*VIV(I,J,1)*(R1D(J)*(UBP(I,J+1)-UBP(I,J))&
               !             -R1C (J)*(UBP(I,J)-UBP(I,J-1)) &
               !             +SOUX (J)* (UBP(I+1,J)-2.0*UBP (I,J)+UBP(I-1,J)) &
               !             +CV1(J)*UBP(I,J)+CV2(J)*(VBP(I+1,J)-VBP(I-1,J)))
            END DO
         END DO
! call exchange_2d(hdtk,1,1)
#if (defined OBCDT)
!---------------------------------------
!ADD BY ZHANGJK 20250329
IF (IX == NX_PROC - 1) THEN
DO J = 2, JMT-1   
hdtk(IMT-2,J) = hdtk(IMT-3,J)
ENDDO 
ENDIF 

IF (IX == 0) THEN 
DO J = 2,JMT-1 
hdtk(3,J) = hdtk(4,J)
ENDDO 
ENDIF 

IF (IY == 0) THEN 
DO I = 2, IMT-1 
hdtk(I,3) = hdtk(I,4)
ENDDO 
ENDIF 

IF (IY == NY_PROC - 1) THEN 
DO I = 2, IMT-1 
hdtk(I,JET-2) = hdtk(I,JET-3)
ENDDO 
ENDIF  
!-------------------------------------
#endif 


!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               ! WORK (I,J)=VIT(I,J,1)*(-1)*(0.5*OTX(J)*(WKA(I+1,J,1)+WKA(I+1,J-1,1) &
               !            -WKA(I,J,1)-WKA(I,J-1,1))+2.0*R2A(J)*(WKA(I,J,2) + &
               !            WKA(I+1,J,2))-2.0*R2B(J)*(WKA (I,J-1,2)+WKA (I+1,J-1,2)))
! !add dissipation 20250224
               WORK (I,J)=VIT(I,J,1)*(-1)*(0.5*OTX(J)*(WKA(I+1,J,1)+WKA(I+1,J-1,1) &
                          -WKA(I,J,1)-WKA(I,J-1,1))+2.0*R2A(J)*(WKA(I,J,2) + &
                          WKA(I+1,J,2))-2.0*R2B(J)*(WKA (I,J-1,2)+WKA (I+1,J-1,2))) - VIT(I,J,1)*hdtk(i,j)*1.0d-6
                        !   print*,work(i,j),hdtk(i,j),vit(i,j,1)
            END DO
!     ENDDO

!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
         END DO
         
! #if (defined OBCDT)
! !---------------------------------------
! !ADD BY ZHANGJK 20250317
! IF (IX == NX_PROC - 1) THEN
! DO J = 2, JMT-1   
! CALL GHOST_POINTS(WORK(IMT-5,J), WORK(IMT-4,J), WORK(IMT-3,J), WORK(IMT-2,J))
! WORK(IMT-2,J) = WORK(IMT-2,J) * VIT(IMT-2,J,1) * VIT(IMT-3,J,1) * VIT(IMT-4,J,1) * VIT(IMT-5,J,1)
! ENDDO 
! ENDIF 

! IF (IX == 0) THEN 
! DO J = 2,JMT-1 
! CALL GHOST_POINTS(WORK(6,J), WORK(5,J), WORK(4,J), WORK(3,J))
! WORK(3,J) = WORK(3,J) * VIT(3,J,1) * VIT(4,J,1) * VIT(5,J,1) * VIT(6,J,1)
! ENDDO 
! ENDIF 

! IF (IY == 0) THEN 
! DO I = 2, IMT-1 
! CALL GHOST_POINTS(WORK(I,6), WORK(I,5), WORK(I,4), WORK(I,3))
! WORK(I,3) = WORK(I,3) * VIT(I,3,1) * VIT(I,4,1) * VIT(I,5,1) * VIT(I,6,1)
! ENDDO 
! ENDIF 

! IF (IY == NY_PROC - 1) THEN 
! DO I = 2, IMT-1 
! CALL GHOST_POINTS(WORK(I,JET-5), WORK(I,JET-4), WORK(I,JET-3), WORK(I,JET-2))
! WORK(I,JET-2) = WORK(I,JET-2) * VIT(I,JET-2,1) * VIT(I,JET-3,1) * VIT(I,JET-4,1) * VIT(I,JET-5,1)
! ENDDO 
! ENDIF  
! !-------------------------------------
! #endif 

!
!  DFS 2024-02-29
!
!     The oringinal code as following
!
!!
!!
!       if (nx_proc == 1) then
!         do j=jsm,jem
!            WORK (1,J) = WORK (IMM,J)
!            WORK (IMT,J) = WORK (2,J)
!         end do
!       end if 
!!
!!
!
#if ( defined OBCDT )
       if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (J)
         do J = JSM, JEM
            WORK (1,J) = WORK (2,J)
            WORK (IMT,J) = WORK (IMM,J)
         end do

!$OMP PARALLEL DO PRIVATE (I)
         do I = 1, IMT
            WORK (I,1) = WORK (I,2)
            WORK (I,JMT) = WORK (I,JEM)
         end do
       end if 
#else
       if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (J)
         do j=jsm,jem
            WORK (1,J) = WORK (IMM,J)
            WORK (IMT,J) = WORK (2,J)
         end do
       end if 
#endif
!
!  DFS 2024-02-29
!

#ifdef SPMD
       call exchange_2d(work,1,0)
#endif

!---------------------------------------------------------------------
!     PREDICTING VB , UB & H0
!---------------------------------------------------------------------
!
! DFS 2021-12
!
#if ( defined POTTIDE )
!$OMP PARALLEL DO PRIVATE (J,I)
     DO J = 1,JMT
        DO I = 1,IMT

           IF ( 1/OHBU(I,J).LT.1000 ) THEN

               WKA(I,J,3)=WKA(I,J,3)-(UBP(I,J)**2+VBP(I,J)**2)**0.5 &
                          *UBP(I,J)/(1/OHBU(I,J)+H0(I,J))*CD_COASTAL

               WKA(I,J,4)=WKA(I,J,4)-(UBP(I,J)**2+VBP(I,J)**2)**0.5 &
                          *VBP(I,J)/(1/OHBU(I,J)+H0(I,J))*CD_COASTAL

           ELSE
!
!               WKA(I,J,3)=WKA(I,J,3)-(UBP(I,J)**2+VBP(I,J)**2)**0.5 &
!                          *UBP(I,J)/(1/OHBU(I,J)+H0(I,J))*CD_PACIFIC
!
!               WKA(I,J,4)=WKA(I,J,4)-(UBP(I,J)**2+VBP(I,J)**2)**0.5 &
!                          *VBP(I,J)/(1/OHBU(I,J)+H0(I,J))*CD_PACIFIC
!
               WKA(I,J,3)=WKA(I,J,3)-(UBP(I,J)**2+VBP(I,J)**2)**0.5 &
                          *UBP(I,J)/(1/OHBU(I,J)+H0(I,J))*CD_TIDE(I,J)

               WKA(I,J,4)=WKA(I,J,4)-(UBP(I,J)**2+VBP(I,J)**2)**0.5 &
                          *VBP(I,J)/(1/OHBU(I,J)+H0(I,J))*CD_TIDE(I,J)
!
           ENDIF
!
!       WKA(I,J,3)=WKA(I,J,3)-(UBP(I,J)**2+VBP(I,J)**2)**0.5&
!                 *UBP(I,J)/(1/OHBU(I,J)+H0(I,J))*CD_TIDE(I,J)
!       WKA(I,J,4)=WKA(I,J,4)-(UBP(I,J)**2+VBP(I,J)**2)**0.5&
!                *VBP(I,J)/(1/OHBU(I,J)+H0(I,J))*CD_TIDE(I,J)
!
        END DO
     ENDDO
#endif
!
! DFS 2021-12
!
!YU  Oct. 24,2005
         CALL SMUV (WKA(1,1,3) ,VIV,1,fil_lat1)
         CALL SMUV (WKA(1,1,4) ,VIV,1,fil_lat1)
         CALL SMZ0 (WORK,VIT,fil_lat1)
!        call FILTER_TRACER(work,vit_1d,FIL_LAT1,1)
!YU  Oct. 24,2005
!
         IF (ISB < 1) THEN

!
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 1,IMT
               UB (I,J)= UBP (I,J) + WKA (I,J,3)* DTB
               VB (I,J)= VBP (I,J) + WKA (I,J,4)* DTB
               H0 (I,J)= H0P (I,J) + WORK (I,J) * DTB
            END DO
         END DO

!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------

            CALL SMUV (UB ,VIV,1,fil_lat1)
            CALL SMUV (VB ,VIV,1,fil_lat1)
            CALL SMZ0 (H0,VIT,fil_lat1)
!           call FILTER_TRACER(h0,vit_1d,FIL_LAT1,1)
#ifdef SPMD
       call exchange_2d(ub,1,1)
       call exchange_2d(vb,1,1)
       call exchange_2d(h0,1,1)         
!       call exchange_pack(ub,vb,h0)
#endif

#if (defined OBCDT)
!---------------------------------------------------------
! open boundary condition when ISB < 1
! add by zhangjk 20250314
!---------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
               UBPP (I,J) = UBP (I,J)
               VBPP (I,J) = VBP (I,J)
               H0PP (I,J) = H0P (I,J)
            END DO
         END DO

      DTB_OPENBC = DTB 

! zhangjk 20260613 
! add exchange 
! #ifdef SPMD
!    call exchange_2d_obc(h0pp,3,4,imt-3,imt-2,jmt-3,jmt-2,3,4)
!    call exchange_2d_obc(ubpp,4,5,imt-3,imt-2,1,0,1,0)
!    call exchange_2d_obc(vbpp,1,0,1,0,jmt-4,jmt-3,3,4)
! #endif

      CALL OBC_H0

! zhangjk 20260613 
! add exchange 
#ifdef SPMD
   call exchange_2d_obc(h0,3,3,imt-2,imt-2,jmt-2,jmt-2,3,3)
#endif

      CALL OBC_UB
      CALL OBC_VB

! zhangjk 20260613 
! add exchange 
#ifdef SPMD
   call exchange_2d_obc(ub,4,4,imt-2,imt-2,jmt-3,jmt-3,3,3)
   call exchange_2d_obc(vb,4,4,imt-2,imt-2,jmt-3,jmt-3,3,3)
   ! call exchange_2d_obc(h0,3,3,imt-2,imt-2,jmt-2,jmt-2,3,3)
#endif

#endif
! #ifdef SPMD
!       call exchange_2d(ub,1,1)
!       call exchange_2d(vb,1,1)
!       call exchange_2d(h0,1,1)         
! #endif
! ---------------------------------------------------------

         IF (IEB == 0) THEN
            ISB = ISB +1
!$OMP PARALLEL DO PRIVATE (J,I)
!Yu         DO J = JSM,JEM
            DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
               DO I = 1,IMT
                  H0F (I,J) = H0F (I,J) + H0 (I,J)
                  H0BF (I,J) = H0BF (I,J) + H0 (I,J)
               END DO
            END DO
            cycle baro_loop
         END IF

         IEB = 0

         cycle baro_loop

      ELSE


!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 1,IMT
               WKA (I,J,1) = UBP (I,J) + WKA (I,J,3)* DTB2
               WKA (I,J,2) = VBP (I,J) + WKA (I,J,4)* DTB2
               WORK(I,J)   = H0P (I,J) + WORK (I,J)* DTB2
            END DO
         END DO

!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------


#ifdef SPMD
!         call exchange_pack(wka(1,1,1),wka(1,1,2),work)
       call exchange_2d(wka(1,1,1),1,1)
       call exchange_2d(wka(1,1,2),1,1)
       call exchange_2d(work,1,1)
#endif

#if (defined OBCDT)
!----------------------------------------------------------------------------------
!ADD BY ZHANGJK 20250314
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
               UBPP (I,J) = UBP (I,J)
               VBPP (I,J) = VBP (I,J)
               H0PP (I,J) = H0P (I,J)
            END DO
         END DO
!----------------------------------------------------------------------------------
#endif


!
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
               UBP (I,J) = AFB2* UB (I,J) + AFB1* (UBP (I,J) + WKA (I,J,1))
               UB (I,J) = WKA (I,J,1)*VIV(I,J,1)
               VBP (I,J) = AFB2* VB (I,J) + AFB1* (VBP (I,J) + WKA (I,J,2))
               VB (I,J) = WKA (I,J,2)*VIV(I,J,1)
               H0P (I,J) = AFB2* H0 (I,J) + AFB1* (H0P (I,J) + WORK(I,J))
               H0 (I,J) = WORK (I,J)
            END DO
         END DO

!---------------------------------------------------------
! open boundary condition when ISB >= 1
! add by zhangjk 20250314
!---------------------------------------------------------
#if (defined OBCDT)
   DTB_OPENBC = DTB2 

! zhangjk 20260613 
! add exchange 
! #ifdef SPMD
!    call exchange_2d_obc(ubpp,4,4,imt-2,imt-2,jmt-3,jmt-3,3,3)
!    call exchange_2d_obc(vbpp,4,4,imt-2,imt-2,jmt-3,jmt-3,3,3)
!    call exchange_2d_obc(h0pp,3,3,imt-2,imt-2,jmt-2,jmt-2,3,3)
! #endif

      CALL OBC_H0

! zhangjk 20260613 
! add exchange 
#ifdef SPMD
   call exchange_2d_obc(h0,3,3,imt-2,imt-2,jmt-2,jmt-2,3,3)
#endif

      CALL OBC_UB
      CALL OBC_VB
! zhangjk 20260613 
! add exchange 
#ifdef SPMD
   call exchange_2d_obc(ub,4,4,imt-2,imt-2,jmt-3,jmt-3,3,3)
   call exchange_2d_obc(vb,4,4,imt-2,imt-2,jmt-3,jmt-3,3,3)
   ! call exchange_2d_obc(h0,3,3,imt-2,imt-2,jmt-2,jmt-2,3,3)
#endif

#endif
! #ifdef SPMD
!       call exchange_2d(ub,1,1)
!       call exchange_2d(vb,1,1)
!       call exchange_2d(h0,1,1)         
! #endif
!---------------------------------------------------------
!
!
!YU  Oct. 24,2005
!lhl0711         IF (MOD(ISB,1200)==0) THEN
         IF (MOD(ISB,1440)==1) THEN
            CALL SMUV (UB ,VIV,1,fil_lat2)
            CALL SMUV (VB ,VIV,1,fil_lat2)
            CALL SMZ0 (H0 ,VIT,fil_lat2)
!           call FILTER_TRACER(h0,vit_1d,FIL_LAT2,1)
            CALL SMUV (UBP,VIV,1,fil_lat2)
            CALL SMUV (VBP,VIV,1,fil_lat2)
            CALL SMZ0 (H0P,VIT,fil_lat2)
!           call FILTER_TRACER(h0p,vit_1d,FIL_LAT2,1)
         END IF
!
!        IF ( isb == 1 .and. mod(month,2) ==0 ) THEN
!           CALL SMOOTH (UB ,VIV(1,1,1),0.5,0)
!           CALL SMOOTH (VB ,VIV(1,1,1),0.5,0)
!           CALL SMOOTH (UBP,VIV(1,1,1),0.5,0)
!           CALL SMOOTH (VBP,VIV(1,1,1),0.5,0)
!           CALL SMOOTH (H0 ,VIT(1,1,1),0.5,1)
!           CALL SMOOTH (H0P,VIT(1,1,1),0.5,1)
!
!           call exchange_2d(ub)
!           call exchange_2d(vb)
!           call exchange_2d(ubp)
!           call exchange_2d(vbp)
!           call exchange_2d(h0)
!           call exchange_2d(h0p)
!        END IF

!YU  Oct. 24,2005
!

         ISB = ISB +1
      END IF
!
! DFS 2021-10
!
#if (defined OBTIDE) || (defined POTTIDE)
!     
       IF ( ISB >= 1 ) Tintg_tide = Tintg_tide + 1.0D0
!
!       IF ( ISB>=1 ) Tintg_tide = Tintg_tide + DTB
!
#endif
!
! DFS
!
! DFS 2021-10
!
!#ifdef TIDE
!lhl20100801
!      IF (ISB >= 1) time_tidal=time_tidal+1.0
!!      IF (ISB >= 1) time_tidal=time_tidal+IDTB
!!      if (mytid.eq.0) print*,time_tidal
!!lhl20100801
!#endif
!
!
! Added by Dai in 2024.11
!
! #if (defined OBCDT)
!         CALL OBC_H0
!         CALL OBC_UB
!         CALL OBC_VB
! #endif

        IF ( ID_OBC_VCON_UVB .EQ. 1) THEN
             CALL OBC_FLUX
             CALL CONSERVE_MASS
        ENDIF

      !   IF ( ID_OBC_VCON_UVB .EQ. 2) THEN
      !        CALL CORRECT_UBVB()
      !   ENDIF

! End by Dai in 2024.11
!
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
            DO I = 1,IMT
               H0F (I,J) = H0F (I,J) + H0 (I,J)
               H0BF (I,J) = H0BF (I,J) + H0 (I,J)
            END DO
         END DO

      END DO baro_loop

      deallocate(dlub,dlvb)
      RETURN
      END SUBROUTINE BAROTR


