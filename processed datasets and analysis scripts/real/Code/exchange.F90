!  CVS: $Id: exchange.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      subroutine exchange_2d(aa,ixx,iyy)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: n1,n2,ixx,iyy
      real(r8) :: aa(imt,jmt),bb(imt),cc(imt),dd(imt),ee(imt),bbb(jmt),ccc(jmt),ddd(jmt),eee(jmt)
!
       n2=mod(mytid,nx_proc)
       n1=(mytid-n2)/nx_proc
!       write(*,*) "n2=",n2,"n1=",n1, "myid=", mytid
       if (ny_proc.ne.1) then
       if (iyy==1) then
!$OMP PARALLEL DO PRIVATE (I)
         do i=1,imt
            bb(i)=aa(i,2)
            cc(i)=aa(i,jem)
          end do

        if (n1< ny_proc-1) then
         call mpi_send(cc,imt,MPI_PR,mytid+nx_proc,tag_2d,mpi_comm_ocn,ierr)
        end if

        if (n1>0) then
         call mpi_recv(dd,imt,MPI_PR,mytid-nx_proc,tag_2d,mpi_comm_ocn,status,ierr)
        end if

        if (n1> 0) then
         call mpi_send(bb,imt,MPI_PR,mytid-nx_proc,tag_2d,mpi_comm_ocn,ierr)
        end if

        if (n1< ny_proc-1) then
         call mpi_recv(ee,imt,MPI_PR,mytid+nx_proc,tag_2d,mpi_comm_ocn,status,ierr)
        end if

        if (n1==0) then
!$OMP PARALLEL DO PRIVATE (I)
         do i=1,imt
            aa(i,jmt)=ee(i)
         end do
     else if (n1==ny_proc-1) then
!$OMP PARALLEL DO PRIVATE (I)
         do i=1,imt
            aa(i,1)=dd(i)
         end do
     else
!$OMP PARALLEL DO PRIVATE (I)
         do i=1,imt
            aa(i,1)=dd(i)
            aa(i,jmt)=ee(i)
         end do
      end if

       end if  !iyy
      else
      end if
       
      if (nx_proc.ne.1) then
       if (ixx==1) then
!$OMP PARALLEL DO PRIVATE (J)
          do j=1,jmt
            bbb(j)=aa(2,j)
            ccc(j)=aa(imt-1,j)
          end do  
         if (n2< nx_proc-1) then
         call mpi_send(ccc,jmt,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,ierr)
         end if

         if (n2>0) then
         call mpi_recv(ddd,jmt,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
         end if

         if (n2>0) then
         call mpi_send(bbb,jmt,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,ierr)
         end if

         if (n2< nx_proc-1) then
         call mpi_recv(eee,jmt,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,status,ierr)
         end if         

        if (n2==0) then

         call mpi_send(bbb,jmt,MPI_PR,mytid+nx_proc-1,tag_3d,mpi_comm_ocn,ierr)
         call mpi_recv(ddd,jmt,MPI_PR,mytid+nx_proc-1,tag_3d,mpi_comm_ocn,status,ierr)

!$OMP PARALLEL DO PRIVATE (J)
         do j=1,jmt
            aa(imt,j)=eee(j)
            aa(1,j)=ddd(j)
         end do
          
     else if (n2==nx_proc-1) then
         

         call mpi_recv(eee,jmt,MPI_PR,mytid-nx_proc+1,tag_3d,mpi_comm_ocn,status,ierr)
         call mpi_send(ccc,jmt,MPI_PR,mytid-nx_proc+1,tag_3d,mpi_comm_ocn,ierr)
          
         

!$OMP PARALLEL DO PRIVATE (J)
         do j=1,jmt
            aa(1,j)=ddd(j)
            aa(imt,j)=eee(j)
         end do
     else
!$OMP PARALLEL DO PRIVATE (J)
         do j=1,jmt
            aa(1,j)=ddd(j)
            aa(imt,j)=eee(j)
         end do
      end if
               
       end if  !ixx
      else
      end if

#endif
     end subroutine exchange_2d
!
!
      subroutine exchange_3d(aa,kk,ixx,iyy)
!     ================
!     To compute bounary of subdomain for the each processor.
 
#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
use dyn_mod
      IMPLICIT NONE

      integer :: kk,n1,n2,ixx,iyy
      real(r8)    :: aa(imt,jmt,kk),bb(imt,kk),cc(imt,kk),dd(imt,kk),ee(imt,kk)
      real(r8)    :: bbb(jmt,kk),ccc(jmt,kk),ddd(jmt,kk),eee(jmt,kk)

       n2=mod(mytid,nx_proc)
       n1=(mytid-n2)/nx_proc
       if (ny_proc.ne.0) then
       if(iyy==1) then
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            bb(i,k)=aa(i,2,k)
            cc(i,k)=aa(i,jem,k)
         end do
         end do


         if (n1< ny_proc-1) then
           call mpi_send(cc,imt*kk,MPI_PR,mytid+nx_proc,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n1>0) then
           call mpi_recv(dd,imt*kk,MPI_PR,mytid-nx_proc,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n1> 0) then
           call mpi_send(bb,imt*kk,MPI_PR,mytid-nx_proc,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n1< ny_proc-1) then
           call mpi_recv(ee,imt*kk,MPI_PR,mytid+nx_proc,tag_3d,mpi_comm_ocn,status,ierr)
         end if
   
     if (n1==0) then
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,jmt,k)=ee(i,k)
         end do
         end do
     else if (n1==ny_proc-1) then
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,1,k)=dd(i,k)
         end do
         end do
     else
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,1,k)=dd(i,k)
            aa(i,jmt,k)=ee(i,k)
         end do
         end do
      end if
     
     end if  !iyy
     else 
     end if
!PY     
     call mpi_barrier(mpi_comm_ocn,ierr)

     if (nx_proc.ne.1) then
     if (ixx==1) then
!M
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,kk
         do j=1,jmt
            bbb(j,k)=aa(2,j,k)
            ccc(j,k)=aa(imt-1,j,k)
         end do
         end do


         if (n2< nx_proc-1) then
           call mpi_send(ccc,jmt*kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2>0) then
           call mpi_recv(ddd,jmt*kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n2> 0) then
           call mpi_send(bbb,jmt*kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2< nx_proc-1) then
           call mpi_recv(eee,jmt*kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,status,ierr)
         end if
   
        if (n2==0) then
         call mpi_send(bbb,jmt*kk,MPI_PR,mytid+nx_proc-1,tag_3d,mpi_comm_ocn,ierr)
         call mpi_recv(ddd,jmt*kk,MPI_PR,mytid+nx_proc-1,tag_3d,mpi_comm_ocn,status,ierr)
         
        end if
!PY
        if (n2==nx_proc-1) then

!PY         write(*,*) mytid
         call mpi_recv(eee,jmt*kk,MPI_PR,mytid-nx_proc+1,tag_3d,mpi_comm_ocn,status,ierr)
         call mpi_send(ccc,jmt*kk,MPI_PR,mytid-nx_proc+1,tag_3d,mpi_comm_ocn,ierr)
         
!!$OMP PARALLEL DO PRIVATE (K,I)
!         do k=1,kk
!         do j=1,jmt
!            aa(1,j,k)=ddd(j,k)
!            aa(imt,j,k)=eee(j,k)
!         end do
!         end do
        end if 
!     else
!M
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,kk
         do j=1,jmt
            aa(1,j,k)=ddd(j,k)
            aa(imt,j,k)=eee(j,k)
         end do
         end do
!      end if     
     end if  !ixx
     else 
     end if
     

#endif
      return
      end subroutine exchange_3d


      subroutine exchange_3d_iso(aa,kk,ixx,iyy)
!     ================
!     To compute bounary of subdomain for the each processor.
 
#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: kk,n1,n2,ixx,iyy
      real(r8)    :: aa(imt,kk,jmt),bb(imt,kk),cc(imt,kk),dd(imt,kk),ee(imt,kk)
      real(r8)    :: bbb(kk,jmt),ccc(kk,jmt),ddd(kk,jmt),eee(kk,jmt)

       n2=mod(mytid,nx_proc)
       n1=(mytid-n2)/nx_proc
      
       if (ny_proc.ne.1) then
       if(iyy==1) then
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            bb(i,k)=aa(i,k,2)
            cc(i,k)=aa(i,k,jem)
         end do
         end do


         if (n1< ny_proc-1) then
           call mpi_send(cc,imt*kk,MPI_PR,mytid+nx_proc,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n1>0) then
           call mpi_recv(dd,imt*kk,MPI_PR,mytid-nx_proc,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n1> 0) then
           call mpi_send(bb,imt*kk,MPI_PR,mytid-nx_proc,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n1< ny_proc-1) then
           call mpi_recv(ee,imt*kk,MPI_PR,mytid+nx_proc,tag_3d,mpi_comm_ocn,status,ierr)
         end if
   
      if (n1==0) then
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,k,jmt)=ee(i,k)
         end do
         end do
     else if (n1==ny_proc-1) then
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,k,1)=dd(i,k)
         end do
         end do
     else
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,k,1)=dd(i,k)
            aa(i,k,jmt)=ee(i,k)
         end do
         end do
      end if
     
     end if  !iyy
     else
     end if
     
     if (nx_proc.ne.1) then
     if (ixx==1) then
!M
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,kk
         do j=1,jmt
            bbb(j,k)=aa(2,j,k)
            ccc(j,k)=aa(imt-1,j,k)
         end do
         end do

         if (n2< nx_proc-1) then
           call mpi_send(ccc,jmt*kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2>0) then
           call mpi_recv(ddd,jmt*kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n2> 0) then
           call mpi_send(bbb,jmt*kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2< nx_proc-1) then
           call mpi_recv(eee,jmt*kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,status,ierr)
         end if
   
        if (n2==0) then
         call mpi_send(bbb,jmt*kk,MPI_PR,mytid+nx_proc-1,tag_3d,mpi_comm_ocn,ierr)
         call mpi_recv(ddd,jmt*kk,MPI_PR,mytid+nx_proc-1,tag_3d,mpi_comm_ocn,status,ierr)
!M
!$OMP PARALLEL DO PRIVATE (K,J)         
         do k=1,kk
         do j=1,jmt
            aa(imt,j,k)=eee(j,k)
            aa(1,j,k)=ddd(j,k)
         end do
         end do
     else if (n2==nx_proc-1) then

!PY        write(*,*) mytid
         call mpi_recv(eee,jmt*kk,MPI_PR,mytid-nx_proc+1,tag_3d,mpi_comm_ocn,status,ierr)
         call mpi_send(ccc,jmt*kk,MPI_PR,mytid-nx_proc+1,tag_3d,mpi_comm_ocn,ierr)

!M         
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,kk
         do j=1,jmt
            aa(1,j,k)=ddd(j,k)
            aa(imt,j,k)=eee(j,k)
         end do
         end do
     else
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,kk
         do j=1,jmt
            aa(1,j,k)=ddd(j,k)
            aa(imt,j,k)=eee(j,k)
         end do
         end do
      end if     
     end if  !ixx
     else
     end if 
     

#endif
      return
      end subroutine exchange_3d_iso

!      ================
      subroutine exchange_pack(aa,bb,cc,ixx,iyy)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      real(r8):: aa(imt,jmt),bb(imt,jmt),cc(imt,jmt)
      real(r8):: outbuf_1(imt*3),outbuf_2(imt*3),outbuf_3(jmt*3),outbuf_4(jmt*3)
      integer:: bsize,pos,n1,n2,ixx,iyy
      bsize=imt*3
!
       n2=mod(mytid,nx_proc)
       n1=(mytid-n2)/nx_proc
      
       if (ny_proc.ne.1) then
       if(iyy==1) then
!$omp parallel do private(i)
      do i=1,imt
         outbuf_2(i)=aa(i,jem)
         outbuf_2(imt+i)=bb(i,jem)
         outbuf_2(2*imt+i)=cc(i,jem)
      end do

      if (n1< ny_proc-1) then
         call mpi_send(outbuf_2,imt*3,MPI_PR,mytid+nx_proc,tag_2d,mpi_comm_ocn,ierr)
      end if

      if (n1>0) then
         call mpi_recv(outbuf_1,imt*3,MPI_PR,mytid-nx_proc,tag_2d,mpi_comm_ocn,status,ierr)
!$omp parallel do private(i)
        do i=1,imt
           aa(i,1)=outbuf_1(i)
           bb(i,1)=outbuf_1(imt+i)
           cc(i,1)=outbuf_1(2*imt+i)
        end do
      end if

!$omp parallel do private(i)
      do i=1,imt
         outbuf_1(i)=aa(i,2)
         outbuf_1(imt+i)=bb(i,2)
         outbuf_1(2*imt+i)=cc(i,2)
      end do

      if (n1> 0) then
         call mpi_send(outbuf_1,imt*3,MPI_PR,mytid-nx_proc,tag_2d,mpi_comm_ocn,ierr)
      end if

      if (n1< ny_proc-1) then
         call mpi_recv(outbuf_2,imt*3,MPI_PR,mytid+nx_proc,tag_2d,mpi_comm_ocn,status,ierr)
!$omp parallel do private(i)
        do i=1,imt
           aa(i,jmt)=outbuf_2(i)
           bb(i,jmt)=outbuf_2(imt+i)
           cc(i,jmt)=outbuf_2(2*imt+i)
        end do
      end if
      
      end if  !iyy
      else
      end if
      
      if (nx_proc.ne.1) then
      if(ixx==1) then
!M
!$OMP PARALLEL DO PRIVATE (J)
      do j=1,jmt
         outbuf_4(j)=aa(imt-1,j)
         outbuf_4(jmt+j)=bb(imt-1,j)
         outbuf_4(2*jmt+j)=cc(imt-1,j)
      end do

      if (n2< nx_proc-1) then
         call mpi_send(outbuf_4,jmt*3,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,ierr)
      else if(n2==nx_proc-1) then
         call mpi_send(outbuf_4,jmt*3,MPI_PR,mytid-nx_proc+1,tag_2d,mpi_comm_ocn,ierr)
      end if

      if (n2>0) then
         call mpi_recv(outbuf_3,jmt*3,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
!$omp parallel do private(J)
        do j=1,jmt
           aa(1,j)=outbuf_3(j)
           bb(1,j)=outbuf_3(jmt+j)
           cc(1,j)=outbuf_3(2*jmt+j)
        end do
      else if(n2==0) then
          call mpi_recv(outbuf_3,jmt*3,MPI_PR,mytid+nx_proc-1,tag_2d,mpi_comm_ocn,status,ierr)  
!$omp parallel do private(J)
        do j=1,jmt
           aa(1,j)=outbuf_3(j)
           bb(1,j)=outbuf_3(jmt+j)
           cc(1,j)=outbuf_3(2*jmt+j)    
        end do  
      end if

!$omp parallel do private(J)
      do j=1,jmt
         outbuf_3(j)=aa(2,j)
         outbuf_3(jmt+j)=bb(2,j)
         outbuf_3(2*jmt+j)=cc(2,j)
      end do

      if (n2> 0) then
         call mpi_send(outbuf_3,jmt*3,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,ierr)
      else if(n2==0) then
         call mpi_send(outbuf_3,jmt*3,MPI_PR,mytid+nx_proc-1,tag_2d,mpi_comm_ocn,ierr)
      end if

      if (n2< nx_proc-1) then
         call mpi_recv(outbuf_4,jmt*3,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,status,ierr)
!$omp parallel do private(J)
        do j=1,jmt
           aa(imt,j)=outbuf_4(j)
           bb(imt,j)=outbuf_4(jmt+j)
           cc(imt,j)=outbuf_4(2*jmt+j)
        end do
      else if(n2==nx_proc-1) then
          call mpi_recv(outbuf_4,jmt*3,MPI_PR,mytid-nx_proc+1,tag_2d,mpi_comm_ocn,status,ierr)
!$omp parallel do private(j)
        do j=1,jmt
           aa(imt,j)=outbuf_4(j)
           bb(imt,j)=outbuf_4(jmt+j)
           cc(imt,j)=outbuf_4(2*jmt+j)
        end do
      end if       
      end if  !ixx
      else
      end if
     
#endif
     end subroutine exchange_pack


!     ===================
      subroutine exchange_boundary(aa,iyy)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: n1,n2,iyy
      real(r8) :: aa(imt)
!
       n2=mod(mytid,nx_proc)
       n1=(mytid-n2)/nx_proc
!PY       write(*,*) "n2=",n2,"n1=",n1 

!PY       if (n1==iyy) then
         if (iyy.eq.0) then
         if (n2< nx_proc-1) then
         call mpi_send(aa(imt-1),1,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,ierr)
         end if

         if (n2>0) then
         call mpi_recv(aa(1),1,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
         end if

         if (n2>0) then
         call mpi_send(aa(2),1,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2< nx_proc-1) then
         call mpi_recv(aa(imt),1,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,status,ierr)
         end if
         end if

         if (iyy.eq.1) then
         if (n2==0) then
         call mpi_send(aa(2),1,MPI_PR,mytid+nx_proc-1,tag_1d,mpi_comm_ocn,ierr)
         call mpi_recv(aa(1),1,MPI_PR,mytid+nx_proc-1,tag_1d,mpi_comm_ocn,status,ierr)
         end if
         
         if (n2==nx_proc-1) then
         call mpi_send(aa(imt-1),1,MPI_PR,mytid-nx_proc+1,tag_1d,mpi_comm_ocn,ierr)
         call mpi_recv(aa(imt),1,MPI_PR,mytid-nx_proc+1,tag_1d,mpi_comm_ocn,status,ierr)
         end if
         end if
          
!PY       end if


#endif
     end subroutine exchange_boundary

!     ===================
      subroutine exchange_1D_boundary(aa,lj,kk,iyy)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: lj,kk,lk,n1,n2,iyy
      real(r8) :: aa(IMT,JMT,KK)
      real(r8) :: send_buf(KK), recv_buf(KK)
!
      n2=mod(mytid,nx_proc)
      n1=(mytid-n2)/nx_proc
!PY       write(*,*) "n2=",n2,"n1=",n1 

!PY       if (n1==iyy) then
      if (iyy.eq.0) then
         if (n2<nx_proc-1) then
            do k = 1, KK
               send_buf(k) = aa(imt-1,lj,k)
            enddo
            call mpi_send(send_buf,kk,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,ierr)
         end if

         if (n2>0) then
            call mpi_recv(recv_buf,kk,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
            do k = 1, KK
               aa(1,lj,k) = recv_buf(k)
            enddo
         end if

         if (n2>0) then
            do k = 1, KK
               send_buf(k) = aa(2,lj,k)
            enddo
            call mpi_send(send_buf,kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2< nx_proc-1) then
            call mpi_recv(recv_buf,kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,status,ierr)
            do k = 1, KK
               aa(imt,lj,k) = recv_buf(k)
            enddo
         end if
      end if

      if (iyy.eq.1) then
         if (n2==0) then
            do k = 1, KK
               send_buf(k) = aa(2,lj,k)
            enddo
            call mpi_send(send_buf,kk,MPI_PR,mytid+nx_proc-1,tag_1d,mpi_comm_ocn,ierr)
            call mpi_recv(recv_buf,kk,MPI_PR,mytid+nx_proc-1,tag_1d,mpi_comm_ocn,status,ierr)
            do k = 1, KK
               aa(1,lj,k) = recv_buf(k)
            enddo
         end if
         
         if (n2==nx_proc-1) then
            do k = 1, KK
               send_buf(k) = aa(imt-1,lj,k)
            enddo
            call mpi_send(send_buf,kk,MPI_PR,mytid-nx_proc+1,tag_1d,mpi_comm_ocn,ierr)
            call mpi_recv(recv_buf,kk,MPI_PR,mytid-nx_proc+1,tag_1d,mpi_comm_ocn,status,ierr)
            do k = 1, KK
               aa(imt,lj,k) = recv_buf(k)
            enddo
         end if
      end if
#endif
     end subroutine exchange_1D_boundary

!     ===================
      subroutine exchange_2D_boundary(aa,kk,NN,JNY,iyy)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: li,lj,kk,lk,n1,n2,iyy
      real(r8) :: aa(IMT,JMT,KK)
      real(r8) :: send_buf(KK*JMT), recv_buf(KK*JMT)
      integer :: NN(JMT), JNY, transfer_num
!
      n2=mod(mytid,nx_proc)
      n1=(mytid-n2)/nx_proc
!PY       write(*,*) "n2=",n2,"n1=",n1 

!PY       if (n1==iyy) then
      if (iyy.eq.0) then
         if (n2<nx_proc-1) then
            transfer_num = 0
            do k = 1, KK
               do lj = jst, jmt
                  if (NN(lj) .ge. JNY) then
                     transfer_num = transfer_num + 1
                     send_buf(transfer_num) = aa(imt-1,lj,k)
                  endif
               enddo
            enddo
            call mpi_send(send_buf,transfer_num,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,ierr)
         end if

         if (n2>0) then
            transfer_num = 0
            do lj = jst, jmt
               if (NN(lj) .ge. JNY) then
                  transfer_num = transfer_num + kk
               endif
            enddo 
            call mpi_recv(recv_buf,transfer_num,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
            transfer_num = 0
            do k = 1, KK
               do lj = jst, jmt
                  if (NN(lj) .ge. JNY) then
                     transfer_num = transfer_num + 1
                     aa(1,lj,k) = recv_buf(transfer_num)
                  endif
               enddo
            enddo
         end if

         if (n2>0) then
            transfer_num = 0
            do k = 1, KK
               do lj = jst, jmt
                  if (NN(lj) .ge. JNY) then
                     transfer_num = transfer_num + 1
                     send_buf(transfer_num) = aa(2,lj,k)
                  endif
               enddo
            enddo
            call mpi_send(send_buf,transfer_num,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2< nx_proc-1) then
            transfer_num = 0
            do lj = jst, jmt
               if (NN(lj) .ge. JNY) then
                  transfer_num = transfer_num + kk
               endif
            enddo 
            call mpi_recv(recv_buf,transfer_num,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,status,ierr)
            transfer_num = 0
            do k = 1, KK
               do lj = jst, jmt
                  if (NN(lj) .ge. JNY) then
                     transfer_num = transfer_num + 1
                     aa(imt,lj,k) = recv_buf(transfer_num)
                  endif
               enddo
            enddo
         end if
      end if

      if (iyy.eq.1) then
         if (n2==0) then
            transfer_num = 0
            do k = 1, KK
               do lj = jst, jmt
                  if (NN(lj) .ge. JNY) then
                     transfer_num = transfer_num + 1
                     send_buf(transfer_num) = aa(2,lj,k)
                  endif
               enddo
            enddo
            call mpi_send(send_buf,transfer_num,MPI_PR,mytid+nx_proc-1,tag_1d,mpi_comm_ocn,ierr)
            call mpi_recv(recv_buf,transfer_num,MPI_PR,mytid+nx_proc-1,tag_1d,mpi_comm_ocn,status,ierr)
            transfer_num = 0
            do k = 1, KK
               do lj = jst, jmt
                  if (NN(lj) .ge. JNY) then
                     transfer_num = transfer_num + 1
                     aa(1,lj,k) = recv_buf(transfer_num)
                  endif
               enddo
            enddo
         end if
         
         if (n2==nx_proc-1) then
            transfer_num = 0
            do k = 1, KK
               do lj = jst, jmt
                  if (NN(lj) .ge. JNY) then
                     transfer_num = transfer_num + 1
                     send_buf(transfer_num) = aa(imt-1,lj,k)
                  endif
               enddo
            enddo
            call mpi_send(send_buf,transfer_num,MPI_PR,mytid-nx_proc+1,tag_1d,mpi_comm_ocn,ierr)
            call mpi_recv(recv_buf,transfer_num,MPI_PR,mytid-nx_proc+1,tag_1d,mpi_comm_ocn,status,ierr)
            transfer_num = 0
            do k = 1, KK
               do lj = jst, jmt
                  if (NN(lj) .ge. JNY) then
                     transfer_num = transfer_num + 1
                     aa(imt,lj,k) = recv_buf(transfer_num)
                  endif
               enddo
            enddo
         end if
      end if
#endif
     end subroutine exchange_2D_boundary



!     ================
      subroutine exchange_2d_real(aa,ixx,iyy)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: n1,n2,ixx,iyy
      real(r4) :: aa(imt,jmt),bb(imt),cc(imt),dd(imt),ee(imt),bbb(jmt),ccc(jmt),ddd(jmt),eee(jmt)
!
       n2=mod(mytid,nx_proc)
       n1=(mytid-n2)/nx_proc
!       write(*,*) "n2=",n2,"n1=",n1, "myid=", mytid
       if (ny_proc.ne.1) then
       if (iyy==1) then
!$OMP PARALLEL DO PRIVATE (I)
         do i=1,imt
            bb(i)=aa(i,2)
            cc(i)=aa(i,jem)
          end do

        if (n1< ny_proc-1) then
         call mpi_send(cc,imt,MPI_PR1,mytid+nx_proc,tag_2d,mpi_comm_ocn,ierr)
        end if

        if (n1>0) then
         call mpi_recv(dd,imt,MPI_PR1,mytid-nx_proc,tag_2d,mpi_comm_ocn,status,ierr)
        end if

        if (n1> 0) then
         call mpi_send(bb,imt,MPI_PR1,mytid-nx_proc,tag_2d,mpi_comm_ocn,ierr)
        end if

        if (n1< ny_proc-1) then
         call mpi_recv(ee,imt,MPI_PR1,mytid+nx_proc,tag_2d,mpi_comm_ocn,status,ierr)
        end if

        if (n1==0) then
!$OMP PARALLEL DO PRIVATE (I)
         do i=1,imt
            aa(i,jmt)=ee(i)
         end do
     else if (n1==ny_proc-1) then
!$OMP PARALLEL DO PRIVATE (I)
         do i=1,imt
            aa(i,1)=dd(i)
         end do
     else
!$OMP PARALLEL DO PRIVATE (I)
         do i=1,imt
            aa(i,1)=dd(i)
            aa(i,jmt)=ee(i)
         end do
      end if

       end if  !iyy
      else
      end if
       
      if (nx_proc.ne.1) then
       if (ixx==1) then
!$OMP PARALLEL DO PRIVATE (J)
          do j=1,jmt
            bbb(j)=aa(2,j)
            ccc(j)=aa(imt-1,j)
          end do  
         if (n2< nx_proc-1) then
         call mpi_send(ccc,jmt,MPI_PR1,mytid+1,tag_2d,mpi_comm_ocn,ierr)
         end if

         if (n2>0) then
         call mpi_recv(ddd,jmt,MPI_PR1,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
         end if

         if (n2>0) then
         call mpi_send(bbb,jmt,MPI_PR1,mytid-1,tag_2d,mpi_comm_ocn,ierr)
         end if

         if (n2< nx_proc-1) then
         call mpi_recv(eee,jmt,MPI_PR1,mytid+1,tag_2d,mpi_comm_ocn,status,ierr)
         end if         

        if (n2==0) then

         call mpi_send(bbb,jmt,MPI_PR1,mytid+nx_proc-1,tag_3d,mpi_comm_ocn,ierr)
         call mpi_recv(ddd,jmt,MPI_PR1,mytid+nx_proc-1,tag_3d,mpi_comm_ocn,status,ierr)

!$OMP PARALLEL DO PRIVATE (J)
         do j=1,jmt
            aa(imt,j)=eee(j)
            aa(1,j)=ddd(j)
         end do
          
     else if (n2==nx_proc-1) then
         

         call mpi_recv(eee,jmt,MPI_PR1,mytid-nx_proc+1,tag_3d,mpi_comm_ocn,status,ierr)
         call mpi_send(ccc,jmt,MPI_PR1,mytid-nx_proc+1,tag_3d,mpi_comm_ocn,ierr)
          
         

!$OMP PARALLEL DO PRIVATE (J)
         do j=1,jmt
            aa(1,j)=ddd(j)
            aa(imt,j)=eee(j)
         end do
     else
!$OMP PARALLEL DO PRIVATE (J)
         do j=1,jmt
            aa(1,j)=ddd(j)
            aa(imt,j)=eee(j)
         end do
      end if
               
       end if  !ixx
      else
      end if

#endif
    end subroutine exchange_2d_real



! zhangjk 20260612 
!     ================
      subroutine exchange_2d_obc(aa,iw1,iw2,ie1,ie2,js1,js2,jn1,jn2)
!     ================

!     Exchange only OBC points located on processor interfaces.   

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: iw1,iw2,ie1,ie2,js1,js2,jn1,jn2
      integer :: n1,n2,side_obc,iloc_obc,jloc_obc
      integer :: ibeg_obc,iend_obc,jbeg_obc,jend_obc
      real(r8) :: aa(imt,jmt),bb,cc,dd,ee

       n2=mod(mytid,nx_proc)
       n1=(mytid-n2)/nx_proc

       if (ny_proc.ne.1) then
         do side_obc=1,2
            if (side_obc.eq.1) then
               if (ix.ne.0) cycle
               ibeg_obc=iw1
               iend_obc=iw2
            else
               if (ix.ne.nx_proc-1) cycle
               ibeg_obc=ie1
               iend_obc=ie2
            end if

            do iloc_obc=ibeg_obc,iend_obc
            bb=aa(iloc_obc,2)
            cc=aa(iloc_obc,jem)

            if (n1<ny_proc-1) then
               call mpi_send(cc,1,MPI_PR,mytid+nx_proc,tag_2d,mpi_comm_ocn,ierr)
            end if

            if (n1>0) then
               call mpi_recv(dd,1,MPI_PR,mytid-nx_proc,tag_2d,mpi_comm_ocn,status,ierr)
            end if

            if (n1>0) then
               call mpi_send(bb,1,MPI_PR,mytid-nx_proc,tag_2d,mpi_comm_ocn,ierr)
            end if

            if (n1<ny_proc-1) then
               call mpi_recv(ee,1,MPI_PR,mytid+nx_proc,tag_2d,mpi_comm_ocn,status,ierr)
            end if

            if (n1==0) then
               aa(iloc_obc,jmt)=ee
            else if (n1==ny_proc-1) then
               aa(iloc_obc,1)=dd
            else
               aa(iloc_obc,1)=dd
               aa(iloc_obc,jmt)=ee
            end if
            end do
         end do
       end if

       if (nx_proc.ne.1) then
         do side_obc=1,2
            if (side_obc.eq.1) then
               if (iy.ne.ny_proc-1) cycle
               jbeg_obc=js1
               jend_obc=js2
            else
               if (iy.ne.0) cycle
               jbeg_obc=jn1
               jend_obc=jn2
            end if

            do jloc_obc=jbeg_obc,jend_obc
            bb=aa(2,jloc_obc)
            cc=aa(imt-1,jloc_obc)

            if (n2<nx_proc-1) then
               call mpi_send(cc,1,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,ierr)
            end if

            if (n2>0) then
               call mpi_recv(dd,1,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
            end if

            if (n2>0) then
               call mpi_send(bb,1,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,ierr)
            end if

            if (n2<nx_proc-1) then
               call mpi_recv(ee,1,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,status,ierr)
            end if

            if (n2==0) then
               aa(imt,jloc_obc)=ee
            else if (n2==nx_proc-1) then
               aa(1,jloc_obc)=dd
            else
               aa(1,jloc_obc)=dd
               aa(imt,jloc_obc)=ee
            end if
            end do
         end do
       end if

#endif
      return
      end subroutine exchange_2d_obc



! zhangjk 20260612 
!     ================
      subroutine exchange_3d_obc(aa,kk,iw1,iw2,ie1,ie2,js1,js2,jn1,jn2)
!     ================
!     Exchange only OBC columns located on processor interfaces.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: kk,iw1,iw2,ie1,ie2,js1,js2,jn1,jn2
      integer :: k_obc,n1,n2,side_obc,iloc_obc,jloc_obc
      integer :: ibeg_obc,iend_obc,jbeg_obc,jend_obc
      real(r8) :: aa(imt,jmt,kk),bb(kk),cc(kk),dd(kk),ee(kk)

       n2=mod(mytid,nx_proc)
       n1=(mytid-n2)/nx_proc

       if (ny_proc.ne.1) then
         do side_obc=1,2
            if (side_obc.eq.1) then
               if (ix.ne.0) cycle
               ibeg_obc=iw1
               iend_obc=iw2
            else
               if (ix.ne.nx_proc-1) cycle
               ibeg_obc=ie1
               iend_obc=ie2
            end if

         do iloc_obc=ibeg_obc,iend_obc
!$OMP PARALLEL DO PRIVATE (K_OBC)
         do k_obc=1,kk
            bb(k_obc)=aa(iloc_obc,2,k_obc)
            cc(k_obc)=aa(iloc_obc,jem,k_obc)
         end do

         if (n1<ny_proc-1) then
            call mpi_send(cc,kk,MPI_PR,mytid+nx_proc,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n1>0) then
            call mpi_recv(dd,kk,MPI_PR,mytid-nx_proc,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n1>0) then
            call mpi_send(bb,kk,MPI_PR,mytid-nx_proc,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n1<ny_proc-1) then
            call mpi_recv(ee,kk,MPI_PR,mytid+nx_proc,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n1==0) then
!$OMP PARALLEL DO PRIVATE (K_OBC)
            do k_obc=1,kk
               aa(iloc_obc,jmt,k_obc)=ee(k_obc)
            end do
         else if (n1==ny_proc-1) then
!$OMP PARALLEL DO PRIVATE (K_OBC)
            do k_obc=1,kk
               aa(iloc_obc,1,k_obc)=dd(k_obc)
            end do
         else
!$OMP PARALLEL DO PRIVATE (K_OBC)
            do k_obc=1,kk
               aa(iloc_obc,1,k_obc)=dd(k_obc)
               aa(iloc_obc,jmt,k_obc)=ee(k_obc)
            end do
         end if
         end do
         end do
       end if

       if (nx_proc.ne.1) then
         do side_obc=1,2
            if (side_obc.eq.1) then
               if (iy.ne.ny_proc-1) cycle
               jbeg_obc=js1
               jend_obc=js2
            else
               if (iy.ne.0) cycle
               jbeg_obc=jn1
               jend_obc=jn2
            end if

         do jloc_obc=jbeg_obc,jend_obc
!$OMP PARALLEL DO PRIVATE (K_OBC)
         do k_obc=1,kk
            bb(k_obc)=aa(2,jloc_obc,k_obc)
            cc(k_obc)=aa(imt-1,jloc_obc,k_obc)
         end do

         if (n2<nx_proc-1) then
            call mpi_send(cc,kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2>0) then
            call mpi_recv(dd,kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n2>0) then
            call mpi_send(bb,kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,ierr)
         end if

         if (n2<nx_proc-1) then
            call mpi_recv(ee,kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,status,ierr)
         end if

         if (n2==0) then
!$OMP PARALLEL DO PRIVATE (K_OBC)
            do k_obc=1,kk
               aa(imt,jloc_obc,k_obc)=ee(k_obc)
            end do
         else if (n2==nx_proc-1) then
!$OMP PARALLEL DO PRIVATE (K_OBC)
            do k_obc=1,kk
               aa(1,jloc_obc,k_obc)=dd(k_obc)
            end do
         else
!$OMP PARALLEL DO PRIVATE (K_OBC)
            do k_obc=1,kk
               aa(1,jloc_obc,k_obc)=dd(k_obc)
               aa(imt,jloc_obc,k_obc)=ee(k_obc)
            end do
         end if
         end do
         end do
       end if

#endif
      return
      end subroutine exchange_3d_obc