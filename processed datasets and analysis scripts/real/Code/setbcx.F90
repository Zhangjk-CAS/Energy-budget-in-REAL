!  CVS: $Id: setbcx.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     ===================================
      SUBROUTINE setbcx (a, imttmp, jmtorkmtmp)
!     ===================================
use precision_mod
use param_mod
      IMPLICIT NONE
      INTEGER :: imttmp,jmtorkmtmp,ktmp
      REAL(r8) :: a(imttmp,jmtorkmtmp)
!
!
!   DFS 2024-03-02
!
!      The original code as following
!
!!
!!
!      if (nx_proc==1) then
!      DO ktmp = 1,jmtorkmtmp
!         a (1,k) = a (imttmp -1,k)
!         a (imttmp,k) = a (2,k)
!      END DO
!      else
!      DO ktmp=1,jmtorkmtmp
!      call exchange_boundary(a(1,k),1)
!      END DO
!      end if
!!
!!
!
      if (nx_proc==1) then

#if ( defined OBCDT )

      DO ktmp = 1,jmtorkmtmp
         a (1,k) = a (2,k)
         a (imttmp,k) = a (imttmp-1,k)
      END DO

#else
      DO ktmp = 1,jmtorkmtmp
         a (1,k) = a (imttmp -1,k)
         a (imttmp,k) = a (2,k)
      END DO
#endif

      else
      DO ktmp=1,jmtorkmtmp
      call exchange_boundary(a(1,k),1)
      END DO
      end if
!
!  DFS 2024-03-02
!
 
      RETURN
      END SUBROUTINE setbcx
 
#else
      SUBROUTINE setbcx ()
      RETURN
      END SUBROUTINE setbcx
#endif 
 
