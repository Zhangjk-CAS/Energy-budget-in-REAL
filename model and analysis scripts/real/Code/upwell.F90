!  CVS: $Id: upwell.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ===============================
      SUBROUTINE UPWELL (UWK,VWK,H0WK)
!     ===============================
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
!
!
 
      IMPLICIT NONE
      REAL(r8):: UWK (IMT,JMT,KM),VWK (IMT,JMT,KM),H0WK (IMT,JMT)
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
      allocate(uk(imt,jmt,km),vk(imt,jmt,km))

!M
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JET
         DO I = 1,IMT
            WORK (I,J) = 0.0D0
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JET
            DO I = 1,IMT
               WKA (I,J,K)= 0.0D0
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JSM,JEM
         DO I = 2,IMM
            WORK (I,J)= 0.25D0* (H0WK (I,J) + H0WK (I -1,J) + H0WK (I,    &
                        J +1) + H0WK (I -1,J +1))
         END DO
      END DO
!
!
!   DFS 2024-03-02
!
!      The original code as following
!
!!
!!
!         if (nx_proc==1 ) then
!         do j=jsm,jem
!         WORK (1,J)= WORK (IMM,J)
!         WORK (IMT,J)= WORK (2,J)
!         end do
!         end if
!!
!!
!
#if ( defined OBCDT )
         if (nx_proc==1 ) then
!$OMP PARALLEL DO PRIVATE (J)
         do j = jsm,jem
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
         if (nx_proc==1 ) then
!$OMP PARALLEL DO PRIVATE (J)
         do j=jsm,jem
         WORK (1,J)= WORK (IMM,J)
         WORK (IMT,J)= WORK (2,J)
         end do
         end if
#endif
!
!   DFS 2024-03-02
!

#ifdef SPMD
      call exchange_2d(work,1,1)
#endif
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UK (I,J,K)= (1.0D0+ WORK (I,J)* OHBU (I,J))* UWK (I,J,K)
               VK (I,J,K)= (1.0D0+ WORK (I,J)* OHBU (I,J))* VWK (I,J,K)
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,K)= 0.5D0* OTX (J)* ( (UK (I +1,J,K) + UK (I +1, J -1,K)) - (UK (I,J,K) + UK (I,J -1,K))) 
               WKA (I,J,K)= WKA(I,J,K)+ R2A (J)*2.0D0* (VK (I,J,K) + VK (I+1,J,K))-R2B (J)*2.0D0* (VK (I,J -1,K) + VK (I +1,J -1,K))      
            END DO
         END DO
      END DO
!
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            WORK (I,J)= 0.0D0
         END DO
      END DO
 
 
      DO K = 1,KM
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WORK (I,J)= WORK(I,J) - DZP (K)* WKA (I,J,K)* VIT (I,J,K)
            END DO
         END DO
      END DO
 
 
      DO K = 2,KM
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WS (I,J,K)= VIT (I,J,K)* (WS (I,J,K -1) + &
               DZP (K -1)* (WORK (I,J)* OHBT (I,J) + WKA (I,J,K -1)))
            END DO
         END DO
      END DO
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            WORK (I,J)= 1.0D0/ (1.0D0+ H0WK (I,J)* OHBT (I,J))
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 2,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               WS (I,J,K)= WS (I,J,K)* WORK (I,J)
            END DO
         END DO
      END DO
 
!

!!
!!    The original code
!! 
!  if (nx_proc ==1 ) then
!!$OMP PARALLEL DO PRIVATE (K,J)
!      DO K = 2,KM
!         DO J = JSM,JEM
!            WS (1,J,K) = WS (IMM,J,K)
!            WS (IMT,J,K) = WS (2,J,K)
!         END DO
!      END DO
!   end if
!
!!
!!
!
#if ( defined OBCDT )
      if (nx_proc ==1 ) then
!$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 2,KM
         DO J = JSM,JEM
            WS (1,J,K) = WS (2,J,K) * VIT(1,J,K)
            WS (IMT,J,K) = WS (IMM,J,K) * VIT(IMT,J,K)
         END DO
      END DO

!$OMP PARALLEL DO PRIVATE (K,I)
      DO K = 2, KM
         DO I = 1, IMT
            WS (I,1,K) = WS (I,2,K) * VIT(I,1,K)
            WS (I,JMT,K) = WS (I,JEM,K) * VIT(I,JMT,K)
         END DO
      END DO
      end if
#else

   if (nx_proc ==1 ) then
!$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 2,KM
         DO J = JSM,JEM
            WS (1,J,K) = WS (IMM,J,K)
            WS (IMT,J,K) = WS (2,J,K)
         END DO
      END DO
   end if

#endif
!
!    DFS 2024-02-29
!

#ifdef SPMD
      call exchange_3d(ws,kmp1,1,1)
#endif
      deallocate(uk,vk)

      RETURN
      END SUBROUTINE UPWELL
