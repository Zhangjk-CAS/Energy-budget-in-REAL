!     =================
      SUBROUTINE CORE_DAILY(TNUM)
!     =================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use forc_mod
use dyn_mod, only: u,v,buffer
use tracer_mod, only: at
#if ( defined SPMD )
use msg_mod
#endif
!
!
! DFS 2024-07
!
#if (defined SWF_DIURNAL)
    use output_mod, only: spval
#endif
!
! DFS 2024-07
!
      IMPLICIT NONE
#include <netcdf.inc>

      integer :: mon_day,irec,tnum
      real(r8),dimension(s_imt,s_jmt) :: t10,u10,v10,slp,q10,swhf,lwhf
      real(r8),dimension(s_imt,s_jmt) :: precr,precs
!zhangjk 20250428,add input runoff
      real(r8),dimension(s_imt_ro,s_jmt_ro) :: friver

      real(r8),dimension(imt,jmt) :: model_sst,es,qs,zz,uu,vv,windx,windy,theta
      real(r8),dimension(imt,jmt) :: core_sensible,core_latent,core_tau
!
      real(r8), parameter :: tok=273.15
      real(r8), parameter :: epsln=1e-25
      real(r8),dimension(imt,jmt) :: tmp1,tmp2
      real(r8),dimension(imt_global,jmt_global) :: tmp3
!
!  DFS 2024-06
!
real(r8) :: Rlwv, Rswv, Rrain
integer  :: Itnum
#if ( defined TOGA_COARE_FLUX )

      real(r8) :: Zum, Ztm, Zqm
      real(r8) :: tsdepth

      real(r8) :: wsx, wsy
      real(r8) :: WWs, TTs, TTa, QQa, PPsf
      real(r8) :: Rlwv, Rswv, Rrain
      real(r8) :: ftime
      integer  :: Itnum

      real(r8) :: qcol_ac, tau_ac
      real(r8) :: hf_old, ef_old, tau_old
      real(r8) :: rf_old
      real(r8) :: tskin

      real(r8) :: fsens, flaten, evpr
      real(r8) :: tauxy, usr
#endif
#if (defined SWF_DIURNAL)
      real(r8) :: fnumdy
      real(r8) :: fstep
      integer  :: nstep
      real(r8) :: swv_rec(NSS)
      integer  :: ktm
#endif
#if (defined SWF_DIURNAL || defined TOGA_COARE_FLUX )
      INTEGER, DIMENSION(imt) :: iglobali
      INTEGER, DIMENSION(jmt) :: jglobalj
      real(r8) :: glat, glon
#endif
!
!
! DFS 2024-06
!

      if ( TNUM.gt.1 ) goto 111

! decide recode number
!      IYFM,MON0,IDAY
!
!      mon_day=0
!      do i=1,mon0-1
!         mon_day=mon_day+nmonth(i)
!      enddo
!
!  DFS 2024-05
!
      MON_DAY = 0

      DO I =1, MON0-1
         IF ( MOD(IYRUN,4).EQ.0 ) THEN
              MON_DAY = MON_DAY + N2MONTH(I)
         ELSE
              MON_DAY = MON_DAY + NMONTH(I)
         ENDIF
      ENDDO
!
! start from a 1000-year spinup
!lhl20150914 irec=(iyfm-1)*365+mon_day+iday-365*(13-1)
!      irec=(iyfm-1)*365+mon_day+iday-365*(11-1)+365*42
!
!      irec=(iyfm-1)*365+mon_day+iday-365*28+365*42
!
!lhl20150914
!      irec=(iyfm-1)*365+mon_day+iday+1
!climatology forcing
!      irec=mon_day+iday+1
!      irec=mon_day+iday
!
! DFS 2024-05
!
      ! irec = ( iyfm-1 )*365 + ( 42-28 )*365 + MON_DAY + IDAY + 11
!      irec = ( iyfm-1 )*365 + ( 42-28 )*365 + MON_DAY + IDAY 
!      irec = ( IYRUN-1948 )*365 + 11 + MON_DAY + IDAY

! zhangjk 20250410
      irec = ( iyfm-1 )*365 + MON_DAY + IDAY
!
! DFS 2024
!
      if (mytid.eq.0) then
!      write(*,*) "iyfm=",iyfm
!      write(*,*) "mon0=",mon0
!      write(*,*) "iday=",iday
!      write(*,*) "mon_day=",mon_day
      write(*,*) "irec=",irec
!      write(*,*) s_imt,s_jmt
      endif

! read in core data
! note the dimensions are s_imt, s_jmt
#ifdef SPMD
      if (mytid.eq.0) then

      ! call read_core(irec,"t_10.db.1948-2007.daymean.05APR2010.nc",t10)
      ! call read_core(irec,"u_10.db.1948-2007.daymean.05APR2010.nc",u10)
      ! call read_core(irec,"v_10.db.1948-2007.daymean.05APR2010.nc",v10)
      ! call read_core(irec,"slp.db.1948-2007.daymean.05APR2010.nc",slp)
      ! call read_core(irec,"q_10.db.1948-2007.daymean.05APR2010.nc",q10)
      ! call read_core(irec,"swdn.db.1948-2007.daymean.05APR2010.nc",swhf)
      ! call read_core(irec,"lwdn.db.1948-2007.daymean.05APR2010.nc",lwhf)
      ! call read_core(irec,"rain.db.1948-2007.daymean.05APR2010.nc",precr)
      ! call read_core(irec,"snow.db.1948-2007.daymean.05APR2010.nc",precs)

! zhangjk 20250410
      call read_core5(irec,"../force-era/t02m-2013-2024.nc",t10)
      call read_core5(irec,"../force-era/u10m-2013-2024.nc",u10)
      call read_core5(irec,"../force-era/v10m-2013-2024.nc",v10)
      call read_core5(irec,"../force-era/slpp-2013-2024.nc",slp)
      call read_coreQ(irec,"../force-era/huss-2013-2024.nc",q10)
      call read_core5(irec,"../force-era/dwsw-2013-2024.nc",swhf)
      call read_core5(irec,"../force-era/dwlw-2013-2024.nc",lwhf)
      call read_core5(irec,"../force-era/prec-2013-2024.nc",precr)
      call read_core5(irec,"../force-era/snow-2013-2024.nc",precs)
! zhangjk 20250428
      call read_core_ro(irec,"../force-era/rivr_climatological.nc",friver)

      precr = precr - precs 
!      call read_core(irec,"t_10.db.clim.daymean.15JUNE2009.nc",t10)
!      call read_core(irec,"u_10.db.clim.daymean.15JUNE2009.nc",u10)
!      call read_core(irec,"v_10.db.clim.daymean.15JUNE2009.nc",v10)
!      call read_core(irec,"slp.db.clim.daymean.15JUNE2009.nc",slp)
!      call read_core(irec,"q_10.db.clim.daymean.15JUNE2009.nc",q10)
!      call read_core(irec,"swdn.db.clim.daymean.15JUNE2009.nc",swhf)
!      call read_core(irec,"lwdn.db.clim.daymean.15JUNE2009.nc",lwhf)
!      call read_core(irec,"rain.db.clim.daymean.15JUNE2009.nc",precr)
!      call read_core(irec,"snow.db.clim.daymean.15JUNE2009.nc",precs)

      endif
!
! DFS 2024-09
!
      if ( allocated(tsa3) ) deallocate(tsa3)
      allocate( tsa3(imt,jmt,12) )

      if ( allocated(wspdu3) ) deallocate(wspdu3)
      allocate( wspdu3(imt,jmt,12) )

      if ( allocated(wspdv3) ) deallocate(wspdv3)
      allocate( wspdv3(imt,jmt,12) )

      if ( allocated(psa3) ) deallocate(psa3)
      allocate( psa3(imt,jmt,12) )

      if ( allocated(qar3) ) deallocate(qar3)
      allocate( qar3(imt,jmt,12) )

      if ( allocated(swv3) ) deallocate(swv3)
      allocate( swv3(imt,jmt,12) )

      if ( allocated(lwv3) ) deallocate(lwv3)
      allocate( lwv3(imt,jmt,12) )

      if ( allocated(rain3) ) deallocate(rain3)
      allocate( rain3(imt,jmt,12) )

      if ( allocated(snow3) ) deallocate(snow3)
      allocate( snow3(imt,jmt,12) )

! zhangjk 20250428, add input runoff
      if (mytid.eq.0) then
      call interplation_ro(friver,tmp3)
      endif
      call global_distribute(tmp3,runoff3(1,1,1))

!
! DFS 2024-09
!
!
! interplate to T grid
! the dimensions are imt and jmt_global
      if (mytid.eq.0) then
      call interplation(t10,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,tsa3(1,1,1))

      if (mytid.eq.0) then
      call interplation(u10,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,wspdu3(1,1,1))

      if (mytid.eq.0) then
      call interplation(v10,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,wspdv3(1,1,1))

      if (mytid.eq.0) then
      call interplation(slp,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,psa3(1,1,1))

      if (mytid.eq.0) then
      call interplation(q10,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,qar3(1,1,1))

      if (mytid.eq.0) then
      call interplation(swhf,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,swv3(1,1,1))

      if (mytid.eq.0) then
      call interplation(lwhf,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,lwv3(1,1,1))

      if (mytid.eq.0) then
      call interplation(precr,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,rain3(1,1,1))

      if (mytid.eq.0) then
      call interplation(precs,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,snow3(1,1,1))

!      call interplation(t10,tsa3_io(1,1,1))
!      call interplation(u10,wspdu3_io(1,1,1))
!      call interplation(v10,wspdv3_io(1,1,1))
!      call interplation(slp,psa3_io(1,1,1))
!      call interplation(q10,qar3_io(1,1,1))
!      call interplation(swhf,swv3_io(1,1,1))
!      call interplation(lwhf,lwv3_io(1,1,1))
!      call interplation(precr,rain3_io(1,1,1))
!      call interplation(precs,snow3_io(1,1,1))
!
!      do j=1,jmt_global
!         do i=1,imt
!            tsa3_io(i,j,1)=tsa3_io(i,j,1)*vit_global(i,j,1)
!            wspdu3_io(i,j,1)=wspdu3_io(i,j,1)*vit_global(i,j,1)
!            wspdv3_io(i,j,1)=wspdv3_io(i,j,1)*vit_global(i,j,1)
!            psa3_io(i,j,1)=psa3_io(i,j,1)*vit_global(i,j,1)
!            qar3_io(i,j,1)=qar3_io(i,j,1)*vit_global(i,j,1)
!            swv3_io(i,j,1)=swv3_io(i,j,1)*vit_global(i,j,1)
!            lwv3_io(i,j,1)=lwv3_io(i,j,1)*vit_global(i,j,1)
!            rain3_io(i,j,1)=rain3_io(i,j,1)*vit_global(i,j,1)
!            snow3_io(i,j,1)=snow3_io(i,j,1)*vit_global(i,j,1)
!         end do
!      end do
!
!      endif
!
!       call global_to_local_4d(tsa3_io(1,1,1),tsa3(1,1,1),1,1)
!       call global_to_local_4d(wspdu3_io(1,1,1),wspdu3(1,1,1),1,1)
!       call global_to_local_4d(wspdv3_io(1,1,1),wspdv3(1,1,1),1,1)
!       call global_to_local_4d(psa3_io(1,1,1),psa3(1,1,1),1,1)
!       call global_to_local_4d(qar3_io(1,1,1),qar3(1,1,1),1,1)
!       call global_to_local_4d(swv3_io(1,1,1),swv3(1,1,1),1,1)
!       call global_to_local_4d(lwv3_io(1,1,1),lwv3(1,1,1),1,1)
!       call global_to_local_4d(rain3_io(1,1,1),rain3(1,1,1),1,1)
!       call global_to_local_4d(snow3_io(1,1,1),snow3(1,1,1),1,1)
!
!
!  DFS 2024-07
!
#if ( defined TOGA_COARE_FLUX )

      if ( idx_start.eq.1 ) then
           do i = 1, imt
           do j = 1, jmt
              sst_skin(i,j) = 0.0
              wrm_qcol_ac(i,j) = 0.0
              wrm_tauw_ac(i,j) = 0.0
              wrm_sflx_old(i,j) = 0.0
              wrm_eflx_old(i,j) = 0.0
              wrm_tauw_old(i,j) = 0.0
              wrm_rf_old(i,j) = 0.0
           enddo
           enddo
      endif
#endif

#if ( defined TOGA_COARE_FLUX || defined SWF_DIURNAL )

      do i = 1, imt
         iglobali(i) = ix*(imt-num_overlap)+i
      enddo

      do i = 1, imt
         if ( iglobali(i)<=imt_global ) then
               fflon_loc(i) = lon( iglobali(i) )
         else
               fflon_loc(i) = lon( imt_global )
         endif
      enddo

      do j = 1, jmt
         if ( iy==0 ) then
              jglobalj(j) = j
         else
              jglobalj(j) = iy*(jmt-num_overlap)+j
         endif
      enddo

      do j = 1, jmt
         if ( jglobalj(j)<=jmt_global ) then
               fflat_loc(j) = lat( jglobalj(j) )
         else
               fflat_loc(j) = lat( jmt_global )
         endif
      enddo

#endif

#if (defined SWF_DIURNAL)
      if ( allocated(swv_new) ) deallocate(swv_new)
      allocate( swv_new(imt,jmt,NSS) )

      do i=1, imt
      do j=1, jmt
         do ktm =1, NSS
            swv_new(i,j,ktm) = 0.0
         enddo
      enddo
      enddo

      do ktm =1, NSS
         swv_rec(ktm) = 0.0
      enddo

      fnumdy = MON_DAY+IDAY

      fstep = DTS
      nstep = NSS

!      do i = 1, imt
!         iglobali(i) = ix*(imt-num_overlap)+i
!      enddo
!
!      do i = 1, imt
!         if ( iglobali(i)<=imt_global ) then
!               fflon_loc(i) = lon( iglobali(i) )
!         else
!               fflon_loc(i) = lon( imt_global )
!         endif
!      enddo
!
!      do j = 1, jmt
!         if ( iy==0 ) then
!              jglobalj(j) = j
!         else
!              jglobalj(j) = iy*(jmt-num_overlap)+j
!         endif
!      enddo
!
!      do j = 1, jmt
!         if ( jglobalj(j)<=jmt_global ) then
!               fflat_loc(j) = lat( jglobalj(j) )
!         else
!               fflat_loc(j) = lat( jmt_global )
!         endif
!      enddo
!
      do j=1,jmt
         glat = fflat_loc(j)/180.0*pi
 
         do i=1,imt
            glon = fflon_loc(i)/180.0*pi

            if ( swv3(i,j,1).ne.spval ) then

                 rswv = swv3(i,j,1)
                 CALL SWF_RECONSTRUCTION(glat, glon, fnumdy, rswv, fstep, nstep, swv_rec)

                 do ktm = 1, nstep
                    swv_new(i,j,ktm) = swv_rec(ktm)
                 enddo
            else
                 do ktm = 1, nstep
                    swv_new(i,j,ktm) = spval
                 enddo
            endif
         enddo
      enddo

!      open(101, file='swv_new.dat')
!           do i = 1, imt
!              if ( fflon_loc(i).eq.120 ) then
!                 do j =1, 10
!                    write(101,*) j
!                    write(101,*) (swv_new(i,j,ktm),ktm=1,nstep)
!                 enddo
!              endif
!          enddo
      close(101) 
#endif
!
!  DFS 2024-07
!
#endif

111    continue

         do j = jsm,jem
            do i = 2,imm
              uu(i,j)=vit(i,j,1)*(u(i,j,1)+u(i+1,j,1)+u(i,j-1,1)+u(i+1,j-1,1)) &
              /(viv(i,j,1)+viv(i+1,j,1)+viv(i,j-1,1)+viv(i+1,j-1,1)+epsln)
              vv(i,j)=vit(i,j,1)*(v(i,j,1)+v(i+1,j,1)+v(i,j-1,1)+v(i+1,j-1,1)) &
              /(viv(i,j,1)+viv(i+1,j,1)+viv(i,j-1,1)+viv(i+1,j-1,1)+epsln)
            end do
         end do
!
!   DFS 2024-03-02
!
!      The original code as following
!
!!
!!
!     if (nx_proc == 1) then
!            do j=jst,jem
!              uu(1,j)=uu(imm,j)
!              uu(imt,j)=uu(2,j)
!              vv(1,j)=vv(imm,j)
!              vv(imt,j)=vv(2,j)
!            end do
!     endif
!!
!!
!
#if ( defined OBCDT )
     if (nx_proc == 1) then
         do j = jst, jem
            uu(1,j) = uu(2,j)
            uu(imt,j) = uu(imm,j)
            vv(1,j) = vv(2,j)
            vv(imt,j) = vv(imm,j)
         end do

         do i = 1, imt
            uu(i,1) = uu(i,2)
            uu(i,jmt) = uu(i,jem)
            vv(i,1) = vv(i,2)
            vv(i,jmt) = vv(i,jem)
         end do
     endif
#else
     if (nx_proc == 1) then
            do j=jst,jem
              uu(1,j)=uu(imm,j)
              uu(imt,j)=uu(2,j)
              vv(1,j)=vv(imm,j)
              vv(imt,j)=vv(2,j)
            end do
     endif
#endif

#ifdef SPMD
       call exchange_2d(uu,1,1)
       call exchange_2d(vv,1,1)
#endif
!

! transfer core data to what the subroutine need

!
! DFS 2024-06
!
#if ( defined TOGA_COARE_FLUX )

      Zum = 10.0
      Ztm = 10.0
      Zqm = 10.0

      tsdepth = 2.5

      Itnum = TNUM

      fstep = DTS

      do j=1,jmt
         do i=1,imt

            glat = fflat_loc(j)
            glon = fflon_loc(i)

            ftime = float(TNUM)*DTS

            if ( vit(i,j,1).gt.0.5 ) then            
                 
                 wsx = ( wspdu3(i,j,1)-uu(i,j))*vit(i,j,1)
!                 wsy = ( wspdv3(i,j,1)+vv(i,j))*vit(i,j,1)
                 wsy = ( wspdv3(i,j,1)-vv(i,j))*vit(i,j,1)

                 WWs = sqrt( wsx*wsx+wsy*wsy+1.0 )
!                 WWs = sqrt( (wspdu3(i,j,1)*wspdu3(i,j,1)+wspdv3(i,j,1)*wspdv3(i,j,1))*vit(i,j,1) )- &
!                       sqrt( (uu(i,j)*uu(i,j)+vv(i,j)*vv(i,j))*vit(i,j,1) )
!
!                 if ( WWs.le.0.0 ) WWs = 0.5
!
                 TTs  = at(i,j,1,1)*vit(i,j,1)
                 TTa  = ( tsa3(i,j,1)-tok )*vit(i,j,1)
                 QQa  = qar3(i,j,1)*vit(i,j,1)
                 PPsf = ( psa3(i,j,1)/100.0 )*vit(i,j,1)

                 Rlwv = lwv3(i,j,1) * vit(i,j,1) * (1.0d0-seaice(i,j))

#if (defined SWF_DIURNAL)
                 Rswv = swv_new(i,j,Itnum) * vit(i,j,1) * (1.0d0-seaice(i,j) )
#else
                 Rswv = swv3(i,j,1) * vit(i,j,1) * (1.0d0-seaice(i,j))
#endif

                 Rrain = rain3(i,j,1)

                 tskin = sst_skin(i,j)
                 qcol_ac = wrm_qcol_ac(i,j)
                 tau_ac =  wrm_tauw_ac(i,j)
                 hf_old = wrm_sflx_old(i,j)
                 ef_old = wrm_eflx_old(i,j)
                 tau_old = wrm_tauw_old(i,j)
                 rf_old = wrm_rf_old(i,j)

                 call toga_coare_ocean_fluxes( Zum, Ztm, Zqm,                   &
                                               tsdepth,                         &
                                               WWs, TTs, TTa, QQa, PPsf,        &
                                               Rlwv, Rswv,                      &
                                               Rrain,                           &
                                               glat, glon,                      &
                                               Itnum, ftime,              &
!                                               jwarm, jcool, jwave,             &
                                               qcol_ac, tau_ac,                 &
!                                               hf_old, ef_old, tau_old,         &
                                               hf_old, ef_old, tau_old, rf_old, &
                                               tskin,                           &
                                               fsens, flaten, evpr, tauxy, usr )
                 sst_skin(i,j) = tskin
                 wrm_qcol_ac(i,j) = qcol_ac
                 wrm_tauw_ac(i,j) = tau_ac
                 wrm_sflx_old(i,j) = hf_old
                 wrm_eflx_old(i,j) = ef_old
                 wrm_tauw_old(i,j) = tau_old
                 wrm_rf_old(i,j) = rf_old

                 sshf(i,j) = fsens * ( 1.0d0-seaice(i,j) )
                 lthf(i,j) = ( flaten - snow3(i,j,1)*3.335e+5 ) * vit(i,j,1) * ( 1.0d0-seaice(i,j) )

!                 lwv(i,j) = ( 0.95*lwv3(i,j,1) - 0.95*5.67E-8*model_sst(i,j)**4)*vit(i,j,1)*(1.0d0-seaice(i,j))
!                 lwv(i,j) = ( 0.97*lwv3(i,j,1) - 0.97*5.67E-8*(TTs+tok)**4 ) * vit(i,j,1) * ( 1.0d0-seaice(i,j) )
                 lwv(i,j) = ( 0.97*lwv3(i,j,1) - 0.97*5.67E-8*(tskin+tok)**4 ) * vit(i,j,1) * ( 1.0d0-seaice(i,j) )

                 tmp1(i,j) =  tauxy * wsx
                 tmp2(i,j) = -tauxy * wsy

                 nswv(i,j) = ( lwv(i,j) + sshf(i,j) + lthf(i,j) )

#if (defined SWF_DIURNAL)
                 swv(i,j) = (1-0.066)*swv_new(i,j,Itnum)*vit(i,j,1)*(1.0d0-seaice(i,j))
#else
                 swv(i,j) = (1-0.066)*swv3(i,j,1)*vit(i,j,1)*(1.0d0-seaice(i,j))
#endif

!                 fresh(i,j) = -( flaten/(2.5e+6) + rain3(i,j,1) + snow3(i,j,1) + runoff(i,j) ) * vit(i,j,1) * &
                 fresh(i,j) = -( evpr + rain3(i,j,1) + snow3(i,j,1) + runoff3(i,j,1) ) * vit(i,j,1) * &
                               ( 1.0d0-seaice(i,j) )

                 ustar(i,j) = usr
            else
                 sshf(i,j) = 0.0D0
                 lthf(i,j) = 0.0D0
                 lwv(i,j) = 0.0D0
                 tmp1(i,j) = 0.0D0
                 tmp2(i,j) = 0.0D0
                 nswv(i,j) = 0.0D0
                 fresh(i,j) = 0.0D0
                 ustar(i,j) = 0.0D0
            endif

         enddo
      enddo
!   
!   tau to U/V grid
!
      do j = jst,jem
         do i = 2,imm
            su(i,j)= 0.25*(tmp1(i,j)+tmp1(i-1,j)+tmp1(i,j+1)+tmp1(i-1,j+1))*viv(i,j,1)
            sv(i,j)= 0.25*(tmp2(i,j)+tmp2(i-1,j)+tmp2(i,j+1)+tmp2(i-1,j+1))*viv(i,j,1)
         end do
      end do

#else

!add by zhangjk 20250703
Itnum = TNUM !when undef TOGA_COARE_FLUX, update the itnum 

      do j=1,jmt
         do i=1,imt

! relative speed to surface currents
         windx(i,j)=(wspdu3(i,j,1)-uu(i,j))*vit(i,j,1)
         windy(i,j)=(wspdv3(i,j,1)-vv(i,j))*vit(i,j,1)
!
!  1.0 is from mom4
         wspd3(i,j,1)=sqrt(windx(i,j)**2+windy(i,j)**2+1.0)*vit(i,j,1)
! using a transient temperature, not daily mean
         model_sst(i,j)=(at(i,j,1,1)+tok)*vit(i,j,1)
         zz(i,j)=10.
         qs(i,j)=0.98*640380*exp(-5107.4/model_sst(i,j))/1.22*vit(i,j,1)
! temperature to potential temperature
         theta(i,j)=tsa3(i,j,1)*(100000.0/psa3(i,j,1))**0.286*vit(i,j,1)
         end do
      end do

! compute heat flux
       call ncar_ocean_fluxes(wspd3(1,1,1),theta(1,1),model_sst(1,1),qar3(1,1,1),qs(1,1),zz(1,1),vit(1,1,1),&
           core_sensible(1,1),core_latent(1,1),core_tau,ustar(1,1))

       do j=1,jmt
          do i=1,imt
            sshf(i,j)=core_sensible(i,j)*vit(i,j,1)*(1.0d0-seaice(i,j))
            lthf(i,j)=(core_latent(i,j)-snow3(i,j,1)*3.335e+5)*vit(i,j,1)*(1.0d0-seaice(i,j))
            lwv(i,j)=(0.95*lwv3(i,j,1)-0.95*5.67E-8*model_sst(i,j)**4)*vit(i,j,1)*(1.0d0-seaice(i,j))
            tmp1(i,j)=core_tau(i,j)*windx(i,j)*vit(i,j,1)
            tmp2(i,j)=-core_tau(i,j)*windy(i,j)*vit(i,j,1)
            nswv(i,j)=(lwv(i,j)+sshf(i,j)+lthf(i,j))
!
! DFS 2024-07
!
#if (defined SWF_DIURNAL)
            swv(i,j) = (1-0.066)*swv_new(i,j,Itnum)*vit(i,j,1)*(1.0d0-seaice(i,j))
#else
            swv(i,j)=(1-0.066)*swv3(i,j,1)*vit(i,j,1)*(1.0d0-seaice(i,j))
#endif
!
! DFS 2024-07
!
            fresh(i,j)=-(core_latent(i,j)/(2.5e+6)+rain3(i,j,1)+snow3(i,j,1)+runoff(i,j))*vit(i,j,1)*(1.0d0-seaice(i,j))
            ustar(i,j)=ustar(i,j)*vit(i,j,1)
         end do
      end do

! tau to U/V grid
        do j = jst,jem
           do i = 2,imm
             su(i,j)= 0.25*(tmp1(i,j)+tmp1(i-1,j)+tmp1(i,j+1)+tmp1(i-1,j+1))*viv(i,j,1)
             sv(i,j)= 0.25*(tmp2(i,j)+tmp2(i-1,j)+tmp2(i,j+1)+tmp2(i-1,j+1))*viv(i,j,1)
           end do
        end do
!
#endif
!
!
!     DFS 2024-03-02
!
!        The original code as following
!!
!!
!     if (nx_proc == 1) then
!            do j=jst,jem
!           su(1,j)=su(imm,j)
!           su(imt,j)=su(2,j)
!           sv(1,j)=sv(imm,j)
!           sv(imt,j)=sv(2,j)
!            end do
!     end if
!
#if ( defined OBCDT )
     if (nx_proc == 1) then
            do j = jst, jem
               su(1,j) = su(2,j)
               su(imt,j) = su(imm,j)
               sv(1,j) = sv(2,j)
               sv(imt,j) = sv(imm,j)
            end do

            do i = 1, imt
               su(i,1) = su(i,2)
               su(i,jmt) = su(i,jem)
               sv(i,1) = sv(i,2)
               sv(i,jmt) = sv(i,jem)
            end do
     end if
#else
     if (nx_proc == 1) then
            do j=jst,jem
               su(1,j)=su(imm,j)
               su(imt,j)=su(2,j)
               sv(1,j)=sv(imm,j)
               sv(imt,j)=sv(2,j)
            end do
     end if
#endif
!
!     DFS 2024-03-02
!

#ifdef SPMD
       call exchange_2d(su,1,1)
       call exchange_2d(sv,1,1)
#endif

!       allocate(buffer(imt_global,jmt_global))
!       call local_to_global_4d_double(su,buffer,1,1)
!       if (mytid.eq.0) then
!       write(133,*) buffer
!       close(133)
!       endif
!       call local_to_global_4d_double(sv,buffer,1,1)
!       if (mytid.eq.0) then
!       write(134,*) buffer
!       close(134)
!       endif
!       deallocate(buffer)
!       stop


      return
      end subroutine CORE_DAILY


!---------------------------------------------
      subroutine read_core5(nnn,fname,var)
!---------------------------------------------
use precision_mod
      use param_mod, only: s_imt,s_jmt
      implicit none
#include <netcdf.inc>

      integer :: start(3),count(3)
      integer :: ncid,iret,nnn
      real(r8) :: var(s_imt,s_jmt)
      character (len=180) :: fname

      start(1)=1;count(1)=s_imt
      start(2)=1;count(2)=s_jmt
      start(3)=nnn;count(3)=1

      iret=nf_open(fname,nf_nowrite,ncid)
      call check_err (iret)

      iret=nf_get_vara_double(ncid,5,start,count,var)
      call check_err (iret)

      iret=nf_close(ncid)
      call check_err (iret)

      return
      end subroutine read_core5

!---------------------------------------------
      subroutine read_coreQ(nnn,fname,var)
!---------------------------------------------
use precision_mod
      use param_mod, only: s_imt,s_jmt
      implicit none
#include <netcdf.inc>

      integer :: start(4),count(4)
      integer :: ncid,iret,nnn
      real(r8) :: var(s_imt,s_jmt)
      character (len=180) :: fname

      start(1)=1;count(1)=s_imt
      start(2)=1;count(2)=s_jmt
      start(4)=nnn;count(4)=1
      start(3)=1;count(3)=1

      iret=nf_open(fname,nf_nowrite,ncid)
      call check_err (iret)

      iret=nf_get_vara_double(ncid,6,start,count,var)
      call check_err (iret)

      iret=nf_close(ncid)
      call check_err (iret)

      return
      end subroutine read_coreQ

!---------------------------------------------
      subroutine read_core_ro(nnn,fname,var)
!---------------------------------------------
use precision_mod
      use param_mod, only: s_imt_ro,s_jmt_ro
      implicit none
#include <netcdf.inc>

      integer :: start(2),count(2)
      integer :: ncid,iret,nnn
      real(r8) :: var(s_imt_ro,s_jmt_ro)
      character (len=180) :: fname

      start(1)=1;count(1)=s_imt_ro
      start(2)=1;count(2)=s_jmt_ro
      ! start(3)=nnn;count(3)=1

      iret=nf_open(fname,nf_nowrite,ncid)
      call check_err (iret)

      iret=nf_get_vara_double(ncid,3,start,count,var)
      call check_err (iret)

      iret=nf_close(ncid)
      call check_err (iret)

      return
      end subroutine read_core_ro

!---------------------------------------------
      subroutine interplation(source,object)
!---------------------------------------------
!
!  
!  DFS 2024-05
!
#include <def-undef.h>
!
!  DFS
!
use precision_mod
      use param_mod, only: imt,jmt,imt_global,jmt_global,s_imt,s_jmt,mytid
      use pconst_mod, only: s_lon,s_lat
      use output_mod, only: spval
      implicit none
!
!  DFS 2024-05
!
#if ( defined OBCDT )
      integer, parameter :: iwk=s_imt,jwk=s_jmt
#else
      integer, parameter :: iwk=s_imt+2,jwk=s_jmt+2
#endif
!
!  DFS 2024-05
!
!      integer, parameter :: iwk=s_imt+2,jwk=s_jmt+2
!
      real(r8) :: source(s_imt,s_jmt)
      real(r8) :: object(imt_global,jmt_global)
      real(r8) :: s_work(iwk,jwk),s_wx(iwk),s_wy(jwk)
      integer :: i,j

       do j=1,s_jmt
       do i=1,s_imt
       if(source(i,j).lt. -1.0e+30)then
        source(i,j)=spval
       endif
       enddo
       enddo
!
!  DFS 2024-05
!
#if ( defined OBCDT )

       do i=1,iwk
          s_wx(i)=s_lon(i)
       enddo

! s_lat is already from north to south!
       do j=1,jwk
          s_wy(j)=s_lat(j)
       enddo
!
! source from north to south!
      do j=1,s_jmt
      do i=1,s_imt
         s_work(i,j)=source(i,s_jmt-j+1)
      enddo
      enddo

#else

       do i=2,iwk-1
       s_wx(i)=s_lon(i-1)
       enddo
       s_wx(1)=s_wx(2)-1.875
       s_wx(iwk)=s_wx(iwk-1)+1.875

! s_lat is already from north to south!
       do j=2,jwk-1
       s_wy(j)=s_lat(j-1)
       enddo
       s_wy(1)=s_wy(2)+1.875
       s_wy(jwk)=s_wy(jwk-1)-1.875
!	if(mytid.eq.0) write(*,*)s_wx
!	if(mytid.eq.0) write(*,*)s_wy
!	stop
!
! source from north to south!
      do j=1,s_jmt
      do i=1,s_imt
      s_work(i+1,j+1)=source(i,s_jmt-j+1)
      enddo
!wrong      s_work(1,j+1)=source(imt,s_jmt-j+1)
      s_work(1,j+1)=source(s_imt,s_jmt-j+1)
      s_work(iwk,j+1)=source(1,s_jmt-j+1)
      enddo
!
      do i=1,s_imt
      s_work(i+1,1)=source(i,2)
      s_work(i+1,jwk)=source(i,1)
      enddo
      s_work(  1,  1)=source(s_imt,  2)
      s_work(iwk,jwk)=source(  1,1)
!
#endif
!
! DFS 2024-05
!
      call hsetup(s_work,iwk,jwk,s_wx,s_wy,object(1,1))

      return
      end subroutine interplation

!zhangjk 20250528
!interplate runoff
!---------------------------------------------
      subroutine interplation_ro(source,object)
!---------------------------------------------
!
!  
!  DFS 2024-05
!
#include <def-undef.h>
!
!  DFS
!
use precision_mod
      use param_mod, only: imt,jmt,imt_global,jmt_global,s_imt_ro,s_jmt_ro,mytid
      use pconst_mod, only: s_lon_ro,s_lat_ro
      use output_mod, only: spval
      implicit none
!
!  DFS 2024-05
!
#if ( defined OBCDT )
      integer, parameter :: iwk=s_imt_ro,jwk=s_jmt_ro
#else
      integer, parameter :: iwk=s_imt_ro+2,jwk=s_jmt_ro+2
#endif
!
!  DFS 2024-05
!
!      integer, parameter :: iwk=s_imt+2,jwk=s_jmt+2
!
      real(r8) :: source(s_imt_ro,s_jmt_ro)
      real(r8) :: object(imt_global,jmt_global)
      real(r8) :: s_work(iwk,jwk),s_wx(iwk),s_wy(jwk)
      integer :: i,j

       do j=1,s_jmt_ro
       do i=1,s_imt_ro
       if(source(i,j).lt. -1.0e+30)then
        source(i,j)=0.0d0 !spval
       endif
       enddo
       enddo
!
!  DFS 2024-05
!
#if ( defined OBCDT )

       do i=1,iwk
          s_wx(i)=s_lon_ro(i)
       enddo

! s_lat is already from north to south!
       do j=1,jwk
          s_wy(j)=s_lat_ro(j)
       enddo
!
! source from north to south!
      do j=1,s_jmt_ro
      do i=1,s_imt_ro
         s_work(i,j)=source(i,s_jmt_ro-j+1)
      enddo
      enddo

#else

       do i=2,iwk-1
       s_wx(i)=s_lon_ro(i-1)
       enddo
       s_wx(1)=s_wx(2)-1.875
       s_wx(iwk)=s_wx(iwk-1)+1.875

! s_lat is already from north to south!
       do j=2,jwk-1
       s_wy(j)=s_lat_ro(j-1)
       enddo
       s_wy(1)=s_wy(2)+1.875
       s_wy(jwk)=s_wy(jwk-1)-1.875
!	if(mytid.eq.0) write(*,*)s_wx
!	if(mytid.eq.0) write(*,*)s_wy
!	stop
!
! source from north to south!
      do j=1,s_jmt_ro
      do i=1,s_imt_ro
      s_work(i+1,j+1)=source(i,s_jmt_ro-j+1)
      enddo
!wrong      s_work(1,j+1)=source(imt,s_jmt-j+1)
      s_work(1,j+1)=source(s_imt_ro,s_jmt_ro-j+1)
      s_work(iwk,j+1)=source(1,s_jmt_ro-j+1)
      enddo
!
      do i=1,s_imt_ro
      s_work(i+1,1)=source(i,2)
      s_work(i+1,jwk)=source(i,1)
      enddo
      s_work(  1,  1)=source(s_imt_ro,  2)
      s_work(iwk,jwk)=source(  1,1)
!
#endif
!
! DFS 2024-05
!
      call hsetup(s_work,iwk,jwk,s_wx,s_wy,object(1,1))

      return
      end subroutine interplation_ro

!---------------------------------------------
      subroutine hsetup(a,mx,my,alon,alat,b)
!---------------------------------------------
!     input : a
!     output: b
!     A bi-linear interpolation will be used for the initial guess, then
!     a refilling procedure be called to redefine the missing data.
!     mx,my            x and y grids number of a
!     alon,alat        lontitude and latitude of a
!     spval            missing flag
!     nf               1 for T grid; 0 for U grid
!

!
!  DFS 2024-05
!
#include <def-undef.h>
!
!  DFS
!

use precision_mod
      use param_mod, only: imt_global,jmt_global,s_imt,s_jmt
      use pconst_mod, only: lon,lat,vit_global_surface
      use output_mod, only: spval
      implicit none
!
      integer :: mx,my,ic,jc,ip,jp,i,j
!      integer :: ind(imt_global,jmt_global)
      real(r8):: b(imt_global,jmt_global)
      real(r8):: a(mx,my),alat(my),alon(mx)
      real(r8), parameter :: isp=99999.0d0
      real(r8):: r1,r2,s1,s2,b1,b2
!
!   Initiallize
!
!  ---nf = 1 on T girds

!  ---b results
      do j=1,jmt_global
      do i=1,imt_global
      b(i,j)=spval
      enddo
      enddo
!

!
!  DFS 2024-05
!
#if ( defined OBCDT )

      do 100 j=1,jmt_global
      do 100 i=1,imt_global
!
      if(vit_global_surface(i,j).lt.0.01) goto 100
!
!  ---find adjacent two grids on x direction (ic)
      ic=isp
      jc=isp
      do 45 ip=2,mx
         if(alon(ip-1).le.lon(i).and.alon(ip).ge.lon(i))then
            ic=ip
            goto 33
         endif
   45 continue
   33 continue

!
!  ---find adjacent two grids on y direction (jc)
      do 50 jp=2,my
         if(alat(jp).le.lat(j).and.alat(jp-1).ge.lat(j))then
            jc=jp
            goto 44
         endif
   50 continue
   44 continue

!
!  ---break if adjacent grids has no found
      if(ic.eq.isp.or.jc.eq.isp) then
      write(*,*)ic,jc
      write(*,*)i,j,alon(i),alat(j)
      write(*,*)lat
      stop 999
      endif
!
!  Bilinear interpolater
!
      r1=(lon(i)-alon(ic-1))/(alon(ic)-alon(ic-1))
      r2=(alon(ic)-lon(i))/(alon(ic)-alon(ic-1))
      s1=(lat(j)-alat(jc-1))/(alat(jc)-alat(jc-1))
      s2=(alat(jc)-lat(j))/(alat(jc)-alat(jc-1))
!


      if(a(ic-1,jc-1).eq.spval.or.a(ic,jc-1).eq.spval.or.&
         a(ic-1,jc).eq.spval.or.a(ic,jc).eq.spval)then

         b1=spval
         b2=spval
         b(i,j)=spval
      else

         b1=a(ic-1,jc-1)*r2+a(ic,jc-1)*r1
         b2=a(ic-1,jc)*r2+a(ic,jc)*r1
         b(i,j)=s1*b2+s2*b1
      endif

  100 continue

#else

      do 100 j=1,jmt_global
      do 100 i=1,imt_global-2
!
      if(vit_global_surface(i,j).lt.0.01) goto 100
!
!  ---find adjacent two grids on x direction (ic)
      ic=isp
      jc=isp
      do 45 ip=2,mx
      if(alon(ip-1).le.lon(i).and.alon(ip).ge.lon(i))then
      ic=ip
      goto 33
      endif
   45 continue
   33 continue

!
!  ---find adjacent two grids on y direction (jc)
      do 50 jp=2,my
      if(alat(jp).le.lat(j).and.alat(jp-1).ge.lat(j))then
      jc=jp
      goto 44
      endif
   50 continue
   44 continue

!
!  ---break if adjacent grids has no found
      if(ic.eq.isp.or.jc.eq.isp) then
      write(*,*)ic,jc
      write(*,*)i,j,alon(i),alat(j)
      write(*,*)lat
      stop 999
      endif
!
!  Bilinear interpolater
!
      r1=(lon(i)-alon(ic-1))/(alon(ic)-alon(ic-1))
      r2=(alon(ic)-lon(i))/(alon(ic)-alon(ic-1))
      s1=(lat(j)-alat(jc-1))/(alat(jc)-alat(jc-1))
      s2=(alat(jc)-lat(j))/(alat(jc)-alat(jc-1))
!


      if(a(ic-1,jc-1).eq.spval.or.a(ic,jc-1).eq.spval.or.&
       a(ic-1,jc).eq.spval.or.a(ic,jc).eq.spval)then

      b1=spval
      b2=spval
      b(i,j)=spval
      else

      b1=a(ic-1,jc-1)*r2+a(ic,jc-1)*r1
      b2=a(ic-1,jc)*r2+a(ic,jc)*r1

      b(i,j)=s1*b2+s2*b1

      endif
!
  100 continue
!
!    set cyclic b. c.
!
      do j=1,jmt_global
      b(imt_global-1,j)     = b(1,j)
      b(imt_global,j)       = b(2,j)
      enddo
!
#endif
!
!  DFS
!
      return
      end subroutine hsetup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
!
! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
! Stephen.Griffies@noaa.gov updated the code with the bug fix. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes (u_del, t, ts, q, qs, z, avail, &
                              sh,lh,tau,ustar)
!                              cd, ch, ce, ustar, bstar       )
use precision_mod
      use param_mod, only: imt,jmt,imt_global,jmt_global,s_imt,s_jmt,mytid
      implicit none
    real(r8)   , intent(in)   , dimension(imt,jmt) :: u_del, t, ts, q, qs, z
    real(r8)   , intent(in)   , dimension(imt,jmt) :: avail
    real(r8)   , intent(inout), dimension(imt,jmt) :: lh,sh,tau
    real(r8)   , dimension(imt,jmt) :: cd, ch, ce, ustar, bstar
!    real   , intent(inout), dimension(imt,jmt) :: cd, ch, ce, ustar, bstar

  real(r8) :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real(r8) :: cd_rt                                ! full drag coefficients @ z
  real(r8) :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real(r8) :: u, u10, tv, tstar, qstar, z0, xx, stab

  integer, parameter :: n_itts = 2
  real(r8), parameter :: grav = 9.80, vonkarm = 0.40,  L=2.5e6, cp=1000.5, r0=1.22
  integer               i, j, jj


  do j=1,jmt
  do i=1,imt
!  do i=1,size(u_del)
    if (avail(i,j) > 0.5 ) then
      tv = t(i,j)*(1+0.608*q(i,j));
      u = max(u_del(i,j), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
      u10 = u;                                                ! first guess 10m wind
    
      cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
      cd_n10_rt = sqrt(cd_n10);
      ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
      stab = 0.5 + sign(0.5,t(i,j)-ts(i,j))
      ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c
  
      cd(i,j) = cd_n10;                                         ! first guess for exchange coeff's at z
      ch(i,j) = ch_n10;
      ce(i,j) = ce_n10;
      do jj=1,n_itts                                           ! Monin-Obukhov iteration
        cd_rt = sqrt(cd(i,j));
        ustar(i,j) = cd_rt*u;                                   ! L-Y eqn. 7a
        tstar    = (ch(i,j)/cd_rt)*(t(i,j)-ts(i,j));                ! L-Y eqn. 7b
        qstar    = (ce(i,j)/cd_rt)*(q(i,j)-qs(i,j));                ! L-Y eqn. 7c
        bstar(i,j) = grav*(tstar/tv+qstar/(q(i,j)+1/0.608));

        zeta     = vonkarm*bstar(i,j)*z(i,j)/(ustar(i,j)*ustar(i,j)); ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
        x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
        x2 = max(x2, 1.0);                                    ! undocumented NCAR
        x = sqrt(x2);
    
        if (zeta > 0) then
          psi_m = -5*zeta;                                    ! L-Y eqn. 8c
          psi_h = -5*zeta;                                    ! L-Y eqn. 8c
        else
          psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
          psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
        end if
    
        u10 = u/(1+cd_n10_rt*(log(z(i,j)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9


        cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10);
        ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
        stab = 0.5 + sign(0.5,zeta)
        ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
        z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic
    
        xx = (log(z(i,j)/10)-psi_m)/vonkarm;
        cd(i,j) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
        xx = (log(z(i,j)/10)-psi_h)/vonkarm;
!!$        ch(i,j) = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2;                !       b (bug)
!!$        ce(i,j) = ce_n10/(1+ce_n10*xx/cd_n10_rt)**2;                !       c (bug)
        ch(i,j) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i,j)/cd_n10) ! 10b (corrected code aug2007)
        ce(i,j) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i,j)/cd_n10) ! 10c (corrected code aug2007)
      end do
    end if
  end do
  end do

    sh=r0*cp*ch*(t-ts)*u_del
    lh=r0*ce*l*(q-qs)*u_del
    tau=r0*cd*u_del

    return
end subroutine ncar_ocean_fluxes
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
#if ( defined TOGA_COARE_FLUX )
      subroutine toga_coare_ocean_fluxes( ZU, ZT, ZQ,                      &
                                          tsdepth,                         &
                                          Wsd, Ts, Ta, Qa, Psf,             &
                                          rlwv, rswv,                      &
                                          rainx,                           &
                                          glat,                            &
                                          glon, Itnum, ftime,              &
!                                          jwarm, jcool, jwave,             &
                                          qcol_ac, tau_ac,                 &
!                                          hf_old, ef_old, tau_old,         &
                                          hf_old, ef_old, tau_old, rf_old, &
                                          tskin,                           &
                                          sflux, eflux, evp, tau, ustr )

      use precision_mod

      use forc_mod, only :  jwarm, jcool, jwave, jamset, jday1, idx_start, &
                            fxp, tk_pwp

      implicit none

      real(r8) :: ZU, ZT, ZQ
      real(r8) :: tsdepth

      real(r8) :: Wsd, Ts, Ta, Qa, Psf
      real(r8) :: rlwv, rswv
      real(r8) :: rainx
      real(r8) :: glat, glon
      real(r8) :: ftime
      integer  :: Itnum
!      integer  :: jwarm, jcool, jwave

      real(r8) :: qcol_ac, tau_ac
      real(r8) :: hf_old, ef_old, tau_old
      real(r8) :: rf_old

      real(r8) :: tskin
      real(r8) :: sflux, eflux, evp, tau, ustr

      real(r8) :: tstr, qstr

      real(r8) :: tauw
      real(r8) :: rnl, rns

      real(r8) :: zibld

      real(r8) :: tctotk
      real(r8) :: glocal
      real(r8) :: Beta, vonK, fdg
      real(r8) :: Rgas, xlv, Cpa, Cpv, rhoa, visa
      real(r8) :: Cpw, rhow, visw

      real(r8) :: soltime      
      real(r8) :: rich, ctd1, ctd2, dtime
      real(r8) :: qr_out, q_pwp, qjoule

      real(r8) :: dt_wrm
      real(r8) :: dsea
      real(r8) :: time_old
      
      real(r8) :: Qas
      real(r8) :: Qs

      real(r8) :: tkt
      real(r8) :: al, be, tcw, bigc, wetc 
      real(r8) :: hsb, hlb, qout, dels, qcol
      real(r8) :: alq, xlamx

      real(r8) :: Du, Dt, Dq
      real(r8) :: Tak, TV, TVstr
      real(r8) :: OBHL
      real(r8) :: ZUL, ZTL, ZQL
      real(r8) :: PUZ, PTZ, PQZ
      real(r8) :: z0, z0t, z0q
      real(r8) :: RR

      real(r8) :: RT, RQ

      real(r8) :: S

      real(r8) :: dter, dqer
      real(r8) :: Bf
      real(r8) :: Wg, DU_Wg

      real(r8) :: dwat, dtmp
!      real(r8) :: dqs_dt
      real(r8) :: alfac
!      real(r8) :: Wb
      real(r8) :: rainflux
!      real(r8) :: tau_rain
!
      real(r8) :: tsw

      real(r8) :: hwave, twave
      real(r8) :: cwave, lwave
      real(r8) :: twopi

      real(r8) :: u10, z010, Cd10, Ch10, Ct10, z0t10
      real(r8) :: Cd, Ct, CC
      real(r8) :: Ribcu, Ribu
      real(r8) :: zetu
      real(r8) :: charn

      integer :: nits
      integer :: idex

!
!       zibld: atmospheric boundary layer depth
!
      zibld = 600.0   
!
!        Factors
!
      Beta = 1.2      ! evaluated from Fairall's low windspeed turbulence data
      vonK = 0.4      ! von Karman's "constant"
      fdg = 1.00      !based on results from Flux workshop August 1995

      tctotk = 273.16  ! Celsius to Kelvin
! 
!        Air constants and coefficients
!
      Rgas = 287.1                     !J/kg/K     gas const. dry air

      Cpa = 1004.67                    !J/kg/K specific heat of dry air (Businger 1982)

      Cpv = Cpa * ( 1+0.84*Qa )      !Moist air - currently not used (Businger 1982)

      rhoa = (Psf*100.0)/( Rgas*(Ta+tctotk)*(1.0+0.61*Qa) ) !kg/m3  Moist air density ( " )

      visa = 1.326e-5*( 1.0 + 6.542e-3*Ta + 8.301e-6 *Ta*Ta - 4.84e-9*Ta*Ta*Ta)   !m2/s
          !Kinematic viscosity of dry air - Andreas (1989) CRREL Rep. 89-11

! 
!        Water constants and coefficients
!
      Cpw = 4000.0                  !J/kg/K specific heat water
      rhow = 1022.0                 !kg/m3  density water
      visw = 1.e-6                  !m2/s kinematic viscosity water

      al = 2.1e-5*(Ts+3.2)**0.79                  ! water thermal expansion coefft.
!
!         Cool skin constants
!
      be = 0.026                                  !salinity expansion coefft.
      tcw = 0.6                                   !W/m/K   Thermal conductivity water
      bigc = 16.0*glocal*Cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)

      CALL gravity( glat, glocal)

      if ( idx_start .eq. 1 ) then
           tsw = Ts                      !henceforth redefined after warm/cool calculations
           dter = 0.3*jcool              !or in "if" block below when jwarm=0.
           tskin = Ts-dter                 !for initial Rnl.  sst: skin temperature C (sst = Ts - dter + dsea)
      endif

!
!       Compute net radiation, updated in flux loop
!           oceanic emissivity 0.97 broadband; albedo 0.055 daily average
!
      if ( (jwarm*jcool).eq.0 ) then
            rnl = 0.97 * ( 5.67E-8*(Ts+tctotk)**4 - rlwv )
      else
            rnl = 0.97 * ( 5.67E-8*(tskin+tctotk)**4 - rlwv )
      endif
!      rnl = 0.97 * ( 5.67E-8*(tskin+tctotk)**4 - rlwv )  ! Net longwave (up=+)

      rns = 0.945 * rswv                               ! Net shortwave (into water)
!
!         START Warm Layer - check switch
!            
      if ( jwarm.eq.0 ) then  ! jump over warm layer calculation
           tsw = Ts
           go to 15           ! convert humidities and calculate fluxes
      endif
     
      soltime = mod( (glon/15.0+ftime/3600.0+24.0),24.0)*3600.0

      if ( idx_start.eq.1) then
           jday1 = 1
           go to 16
      endif

      if ( soltime .lt. time_old ) then
!
!      reset all variables at local midnight
!     
           jday1 = 0
           jamset =0

           tau_ac  = 0.0
           qcol_ac = 0.0
           dt_wrm = 0.0
!
!      initial guess at warm layer parameters expected in early morning
!      fxp=0.5 implies a shallow heating layer to start the integration;
!      tk_pwp=19.0 implies warm layer thickness is a maximum from the day
!      before and is not meant to match this timestep's fxp.     
!
           fxp = 0.5
           tk_pwp = 19.0
           tsw = Ts

           go to 16

      else if ( soltime.gt.21600.0 .and. jday1.eq.1 ) then  
!
!         6 am too late to start on first day
!
           dt_wrm = 0.0
           tsw = Ts

           go to 16

      else
!
!          compute warm layer. Rnl and "_old"s from previous timestep
!
           rich = 0.65                                    !critical Rich. No.
           ctd1 = sqrt(2.0*rich*Cpw/(al*glocal*rhow))        !u*^2 integrated so
           ctd2 = sqrt(2.0*al*glocal/(rich*rhow))/(Cpw**1.5) !has /rhow in both

           dtime = soltime - time_old                      !delta time

           qr_out = rnl + hf_old + ef_old + rf_old         !flux out from previous pass
!           qr_out = rnl + hf_old + ef_old                   !flux out from previous pass
           q_pwp = fxp*rns - qr_out                         !effective net warming

           if ( q_pwp .lt. 50.0 .and. jamset.eq.0 ) then   !integration threshold
                tsw = Ts
                go to 16
           endif

           jamset = 1                                   ! indicates integration has started

           tau_ac = tau_ac + max(0.002,tau_old)*dtime   ! momentum integral

           if ( (qcol_ac+q_pwp*dtime) .gt. 0.0 ) then

                do idex = 1, 5                           !iterate for warm layer thickness
                   fxp = 1.0 - ( 0.28*0.014*(1.0-dexp(-tk_pwp/0.014)) + &
                                 0.27*0.357*(1.0-dexp(-tk_pwp/0.357)) + &
                                 0.45*12.82*(1.0-dexp(-tk_pwp/12.82)) )/tk_pwp  !solar absorb. prof

                   qjoule = ( fxp*rns - qr_out ) * dtime

                   if ( (qcol_ac+qjoule) .gt. 0.0 ) then
                        tk_pwp = min( 19.0, ctd1*tau_ac/sqrt(qcol_ac+qjoule) ) ! warm layer thickness
                   endif
                enddo
           else
                fxp = 0.75
                tk_pwp = 19.0
                qjoule = (fxp*rns-qr_out)*dtime  
           endif

           qcol_ac = qcol_ac + qjoule               !integrate heat input
           if ( qcol_ac.gt.0.0 ) then
                dt_wrm = ctd2*(qcol_ac)**1.5/tau_ac  !pwp model warming
           else
                dt_wrm = 0.0
           endif
        
      endif
         
      if ( tk_pwp.lt.tsdepth ) then            !sensor deeper than pwp layer
           dsea = dt_wrm                          !all warming must be added to ts
      else                                      !warming deeper than sensor
           dsea = dt_wrm*tsdepth/tk_pwp          !assume linear temperature profile
      endif

      tsw = Ts + dsea                             !add warming above sensor for new ts

16    time_old = soltime
!
!     End of Warm Layer 
!
!  
15    CALL Teten_Humidity( Ta, Psf, Qas )         !Teten's formula returns sat. air in mb

      Qas = 0.62197 * ( Qas/( Psf-0.378*Qas ) )   !convert from mb to spec. humidity  kg/kg

      CALL Teten_Humidity( Tsw, Psf, Qs)           !sea QS returned in mb      
      Qs = Qs * 0.98                              !reduced for salinity Kraus 1972 p. 46
      Qs = 0.62197 * ( Qs/( Psf-0.378*Qs ) )      !convert from mb to spec. humidity  kg/kg

      xlv = ( 2.501-0.00237*tsw )*1e+6              !J/kg  latent heat of vaporization at TS
      wetc = 0.622*xlv*Qs/(Rgas*(tsw+tctotk)**2)    !correction for dq;slope of sat. vap.
!
      Tak = Ta + tctotk
!
!       Wave parameters
!
      twopi = 6.2831852
!
!      Default values for wave height and period for equilibrium sea
!
      hwave = 0.018*Wsd*Wsd*(1.0+0.015*Wsd)
      twave = 0.729*Wsd
      cwave = glocal*twave/twopi
      lwave = cwave*twave
!
!         Initial guesses
!
      dter = 0.3*jcool                    ! cool skin Dt
      dqer = wetc*dter                    ! cool skin Dq

      Wg = 0.5                            !Gustiness factor initial guess
      z0 = 0.0001                         !roughness initial guess

      tkt = 0.001*jcool                   !guess sublayer thickness
!
!       Air-sea differences - includes warm layer in Dt and Dq
!
      Du = ( Wsd**2.0 + Wg**2.0 )**0.5  !include gustiness in wind spd. difference

      Dt = tsw-Ta-0.0098*ZT              !potential temperature diff        
      Dq = Qs-Qa
!
!   **************** neutral coefficients ******************
!
      u10 = Du*dlog(10.0/z0)/dlog(ZU/z0)
      ustr = 0.035*u10
      z010 = 0.011*ustr*ustr/glocal + 0.11*visa/ustr

      Cd10 = (vonK/dlog(10.0/z010))**2
      Ch10 = 0.00115
      Ct10 = Ch10/sqrt(Cd10)
      z0t10 = 10.0/dexp(vonK/Ct10)
      Cd = (vonK/dlog(ZU/z010))**2
!
!     
! ************* Grachev and Fairall (JAM, 1997) **********
!
      Ct = vonK/dlog(ZT/z0t10)                  ! Temperature transfer coefficient
      CC = vonK*Ct/Cd                           ! z/L vs Rib linear coefficient
      Ribcu = -ZU/(zibld*0.004*Beta**3)         ! Saturation or plateau Rib
      Ribu = -glocal*ZU*((Dt-dter)+0.61*Tak*Dq)/(Tak*Du**2)

      if ( Ribu.lt.0.0 ) then
          zetu = CC*Ribu/(1.0+Ribu/Ribcu)         ! Unstable G and F
      else
          zetu = CC*Ribu*(1.0+27.0/9.0*Ribu/CC)   ! Stable
      endif

      OBHL = ZU/zetu                              ! MO length

      if ( zetu.gt.50.0 ) then
           nits=1
      else
           nits=3   ! number of iterations
!           nits = 10   ! number of iterations
      endif
!
!     First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate zo and z/L
!

      ZUL = ZU/OBHL
      CALL PSI_COARE3(1, ZUL, PUZ)
!      CALL PSI_COARE_SHEBA(1, ZUL, PUZ)

      ustr = Du*vonK/(dlog(ZU/z010)-PUZ)

      ZTL = ZT/OBHL
      CALL PSI_COARE3(2, ZTL, PTZ) 
!      CALL PSI_COARE_SHEBA(2, ZTL, PTZ) 

      tstr = -(Dt-dter)*vonK/(dlog(ZT/z0t10)-PTZ)

      ZQL = ZQ/OBHL
      CALL PSI_COARE3(2, ZQL, PQZ) 
!      CALL PSI_COARE_SHEBA(2, ZQL, PQZ) 

      qstr = -(Dq-dqer)*vonK/(dlog(ZQ/z0t10)-PQZ)
      
      charn = 0.011                  !then modify Charnock for high wind speeds Chris' data
      if ( Du.gt.10.0 ) charn = 0.011+(0.018-0.011)*(Du-10.0)/(18.0-10.0)
      if ( Du.gt.18.0 ) charn=0.018

!      
! **** Iterate across u*(t*,q*),z0(z0t,z0q) and z/L including cool skin ****
!
      do idex = 1, nits

         if ( Jwave.eq.0 ) then

             z0 = charn*ustr*ustr/glocal + 0.11*visa/ustr               !after Smith 1988

         else if ( Jwave.eq.1 ) then

             z0 = (50.0/twopi)*lwave*(ustr/cwave)**4.5 + 0.11*visa/ustr  !Oost et al.

         else if ( Jwave.eq.2 ) then

             z0 = 1200.0*hwave*(hwave/lwave)**4.5 + 0.11*visa/ustr      !Taylor and Yelland 

         endif

         RR = z0*ustr/visa

!!
!!        *** z0q and z0t fitted to results from several ETL cruises ************
!!
         z0q = min(1.15e-4,5.5e-5/RR**0.6)
         z0t = z0q
!
!            The other methold to calculate z0t, z0q, but not coare3 way
!
!         CALL LKB_COAMPS( RR, RT, RQ )
!         z0t = RT*visa/ustr
!         z0q = RQ*visa/ustr
!
         TV = Tak*(1.0+0.61*Qa)
         TVstr = tstr*(1.0+0.61*Qa) + 0.61*Tak*qstr
         OBHL = TV*ustr*ustr/(glocal*vonK*TVstr )

         ZUL = ZU/OBHL
         CALL PSI_COARE3(1, ZUL, PUZ)
!         CALL PSI_COARE_SHEBA(1, ZUL, PUZ)

         ZTL = ZT/OBHL
         CALL PSI_COARE3(2, ZTL, PTZ) 
!         CALL PSI_COARE_SHEBA(2, ZTL, PTZ) 

         ZQL = ZQ/OBHL
         CALL PSI_COARE3(2, ZQL, PQZ) 
!         CALL PSI_COARE_SHEBA(2, ZQL, PQZ) 

         dqer = wetc*dter*jcool

         ustr = Du*vonK/( dlog(ZU/z0)-PUZ )
         tstr = -(Dt-dter)*vonK/( dlog(ZT/z0t)-PTZ )
         qstr = -(Dq-dqer)*vonK/( dlog(ZQ/z0q)-PQZ )

         TVstr = tstr*(1.0+0.61*Qa)+(0.61*Tak*qstr)
         Bf = -glocal/Tak*ustr*TVstr

         if ( Bf.gt.0.0 ) then
              Wg = Beta*(Bf*zibld)**0.333  
         else
              Wg = 0.2
         endif

         Du = sqrt(Wsd**2.0+Wg**2.0)       !include gustiness in wind spd.
     
!         rnl = 0.97*(5.67e-8*(tsw-dter+tctotk)**4-rlwv)  !Recompute net longwave; cool skin=-2W/m2
!
!           Cool skin
!
         if ( Jcool.ne.0 ) then                    !Cool skin

              rnl = 0.97*(5.67e-8*(tsw-dter+tctotk)**4-rlwv)  !Recompute net longwave; cool skin=-2W/m2

              hsb = -rhoa * Cpa * ustr * tstr
              hlb = -rhoa * xlv * ustr * qstr
              qout = rnl + hsb + hlb

              dels = rns*( 0.065+11.0*tkt-6.6e-5/tkt*(1.0-dexp(-tkt/8.0e-4)) )      ! Eq.16 Ohlmann
              qcol = qout - dels
              
              alq = al*qcol + be*hlb*Cpw/xlv                    ! Eq. 7 Buoy flux water

              if ( alq.gt.0.0 ) then                                      !originally (qcol.gt.0)
                   xlamx = 6.0/(1.0+(bigc*alq/ustr**4)**0.75)**0.333      !Eq 13 Saunders coeff.
                   tkt = xlamx*visw/(sqrt(rhoa/rhow)*ustr)                !Eq.11 Sublayer thickness
              else
                   xlamx = 6.0                                            !prevent excessive warm skins
                   tkt = min(0.01,xlamx*visw/(sqrt(rhoa/rhow)*ustr))      !Limit tkt
              endif

              dter = qcol*tkt/tcw                                 ! Eq.12 Cool skin
              dqer = wetc*dter

         endif          !end cool skin

      enddo   !end iterations                 
!
!
!         Compute surface fluxes and other parameters
!
      tskin = tsw-dter*jcool                  ! final skin temperature this timestep

      S = sqrt( Wsd*Wsd + wg*wg )            ! velocity incl. gustiness param.
      tau = rhoa*ustr*ustr/S          
      tauw = rhoa*ustr*ustr*Wsd/S           ! kinematic units,  stress N/m2

      sflux = Cpa * rhoa * ustr * tstr    ! sensible fluxes W/m2,  !!!downward
      eflux = xlv * rhoa * ustr * qstr    ! latent fluxes W/m2
      evp = rhoa * ustr * qstr
!!
!!           compute heat flux due to rainfall
!!
      dwat = 2.11e-5 * ( ( Ta+tctotk )/tctotk )**1.94                ! water vapour diffusivity
      dtmp = ( 1.0+3.309e-3*Ta-1.44e-6*Ta*Ta ) * 0.02411/(rhoa*Cpa)  !heat diffusivity
!      dqs_dt = Qas*xlv/( Rgas*(Ta+tctotk)**2 )                       !Clausius-Clapeyron
      alfac = 1.0/( 1.0+(wetc*xlv*dwat)/(Cpa*dtmp) )                 !wet bulb factor
!
!!      rainflux = rainx*alfac*Cpw*( (sst-Ta)+(Qs-Qa-dqer)*xlv/Cpa )/3600.0
      rainflux = rainx*alfac*Cpw*( (tskin-Ta)+(Qs-Qa-dqer)*xlv/Cpa )
!!
!!
!!       Compute momentum flux due to rainfall
!!
!!       tau_rain = 0.85*rainx/3600*Ws
!       tau_rain = 0.85*rainx*Ws
!!  
!!           Webb correction to latent heat flux already in ef via zoq/rr function s0 return Wbar
!!
!       Wb = -1.61*ustr*qstr/(1.0+1.61*Qa)-ustr*tstr/(Ta+tctotK)
!!
!
!          Save fluxes for next timestep warm layer integrals
!
!       if ( jwarm.ne.0 ) then 
            tau_old = tauw 
            ef_old = -eflux
            hf_old = -sflux
            rf_old = rainflux
!       else
!            tau_old = 0.0
!            ef_old = 0.0
!            hf_old = 0.0
!!            rf_old = 0.0
!       endif

      return

      end subroutine toga_coare_ocean_fluxes
#endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      Subroutine gravity( flat, gloc )
!
!       calculates g as a funciton of latitude using the 1980 IUGG formula
!         
!       Bulletin Geodesique, Vol 62, No 3, 1988 (Geodesist's Handbook)
!       p 356, 1980 Gravity Formula (IUGG, H. Moritz)
!       units are in m/sec^2 and have a relative precision of 1 part
!       in 10^10 (0.1 microGal)
!       code by M. Zumberge.
!
!       check values are:
!
!        g = 9.780326772 at latitude  0.0
!        g = 9.806199203 at latitude 45.0
!        g = 9.832186368 at latitude 90.0
!
      use precision_mod
      implicit none

      real(r8) :: flat, gloc
      real(r8) :: gam, c1, c2, c3, c4, phi
 
      gam = 9.7803267715
      c1 = 0.0052790414
      c2 = 0.0000232718
      c3 = 0.0000001262
      c4 = 0.0000000007
      phi = flat * 3.14159265358979 / 180.0
      gloc = gam * ( 1.0 + c1 * ((sin(phi))**2) + &
                           c2 * ((sin(phi))**4) + &
                           c3 * ((sin(phi))**6) + &
                           c4 * ((sin(phi))**8) )

      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SWF_RECONSTRUCTION( flat, flon, fnumdy, swv, tstep, nstp, swv_rec)
!
!        diurnal short wave reconstrucktion from daily mean 
!        based on 
!           (1) D. J. Bernie, et al., 2007, Impact of resolving the diurnal
!        cycle in an ocean-atmosphere GCM. Part 1: a diurnal forced OGCM, Clim
!        Dyn (2007) 29: 575-590
!
!           (2)  Henry Rachele and Arnold Tunick 1994,
!        Journal of Applied Meteorology,1994,33,964-976.
! 
      use output_mod, only: spval

      use precision_mod
      implicit none

      real(r8) :: flat, flon
      real(r8) :: fnumdy
      real(r8) :: swv
      real(r8) :: tstep
      integer  :: nstp
      real(r8) :: swv_rec(nstp)

      real(r8) :: sigma
      real(r8) :: pi
      real(r8) :: biga, bigb, bigc, bigd
      real(r8) :: tdawn, tduck
      real(r8) :: tx

      real(r8) :: d
      real(r8) :: gama
      real(r8) :: beta
      real(r8) :: epsim
      real(r8) :: sinsigma, cossigma
      real(r8) :: bigM

      real(r8) :: bigf
      real(r8) :: ftmp

      real(r8) :: tmid

      real(r8) :: sstar
      real(r8) :: fintg
      real(r8) :: swvi

      real(r8) :: ti
      real(r8) :: tind
      real(r8) :: tloc
      real(r8) :: tdawn_loc
      real(r8) :: tduck_loc
      integer  :: ktm

      real(r8) :: fftmp
      
      pi = 3.1415926

      d = (fnumdy-1.0)*2.0*pi/365.242
      gama = 279.9348*pi/180.0
      beta = gama + &
             ( 0.4087*sin(gama) + 1.8724*cos(gama) - 0.0182*sin(2.0*gama) + &
               0.0083*cos(2.0*gama) )*pi/180.0

      epsim = 23.4438*pi/180.0

      sinsigma = sin(epsim)*sin(beta)
      sigma = asin(sinsigma)
      cossigma = cos(sigma)

      bigm = ( 12.0 + 0.12357*sin(d) - 0.004289*cos(d) + &
               0.153809*sin(2.0*d) + 0.060783*cos(2.0*d) )/24.0
 
      biga = sin(flat)*sinsigma
      bigb = cos(flat)*cossigma

      bigc = flon - 2.0*bigm*pi

      bigd = 2.0*pi

      tmid = bigm - flon/(2.0*pi)

      bigf = abs(biga)/abs(bigb)

      if ( bigf.le.1.0 ) then

           tx = ( acos(-biga/bigb)-bigc )/bigd

           ftmp = -bigb*bigd*sin(bigc+bigd*tx)

           if ( ftmp.gt.0.0 ) then
                tdawn = tx
                tduck = 2.0*tmid-tdawn
           else
                tduck = tx
                tdawn = 2.0*tmid-tduck
           endif
      else
           if ( biga.lt.bigb ) then
                tdawn = 0.0
                tduck = 0.0
           else
                tdawn = 0.0
                tduck = 1.0
           endif

      endif

      fintg = biga*(tduck-tdawn) + bigb/bigd*(sin(bigd*tduck+bigc)-sin(bigd*tdawn+bigc))
      sstar = swv*86400.0/fintg

      do ktm = 1, nstp
         ti = tstep*(ktm-0.5)
         tind = ti/86400.0
!         tloc = tind + flon/pi*180.0/15.0/24.0
!
         tloc = tind + flon/(2.0*pi)
         if ( tloc.gt.1.0 ) then
              tloc = tloc-1.0
         else if ( tloc.lt.0.0 ) then
              tloc = tloc+1.0
         endif

         tdawn_loc = tdawn + flon/(2.0*pi)
         tduck_loc = tduck + flon/(2.0*pi)

         if ( tloc.ge.tdawn_loc .and. tloc.le.tduck_loc ) then      

             fftmp = biga+bigb*cos(bigc+bigd*tind)
             if ( fftmp.gt.0.0 ) then
                  swvi = biga*tstep/86400.0 + &
                         bigb/bigd * ( sin(bigd*(ti+0.5*tstep)/86400.0+bigc) - &
                                 sin(bigd*(ti-0.5*tstep)/86400.0+bigc) )
 
                  swv_rec(ktm) = sstar*swvi/tstep
             else
                  swv_rec(ktm) = 0.0
             endif
         else
             swv_rec(ktm) = 0.0
         endif
      enddo

      return

      end subroutine SWF_RECONSTRUCTION

!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Teten_Humidity( T, P, Qsat)                                 
!
! Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532 
!
      use precision_mod
      implicit none

      real(r8) :: T, P, Qsat
     
      Qsat = (1.0007+3.46e-6*P)*6.1121*dexp(17.502*T/(240.97+T))

      return

      end
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

      SUBROUTINE LKB_COAMPS( rr, rt, rq)

      use precision_mod
      implicit none

      real(r8) :: rr, rt, rq      
      
      if (rr.le.0.11) then
          rt=0.177
          rq=0.292
      else if (rr.le.0.8) then
!         rt=1.376*rr**.929
!         rq=1.808*rr**.826
          rt = 1.376 * exp(0.929*log(rr))
          rq = 1.808 * exp(0.826*log(rr))
      else if (rr.le.3.0) then
!         rt=1.026*rr**(-.599) 
!         rq=1.393*rr**(-.528)
          rt = 1.026 * exp(-0.599*log(rr))
          rq = 1.393 * exp(-0.528*log(rr))
      else if (rr.le.10.0) then
!         rt=1.625*rr**(-1.018)
!         rq=1.956*rr**(-.870)
          rt = 1.625 * exp(-1.018*log(rr))
          rq = 1.956 * exp(-0.870*log(rr))
      else if (rr.le.30.0) then
!         rt=4.661*rr**(-1.475)
!         rq=4.994*rr**(-1.297)
          rt = 4.661 * exp(-1.475*log(rr))
          rq = 4.994 * exp(-1.297*log(rr))
      else if (rr.le.100.0) then
!         rt=34.904*rr**(-2.067)
!         rq=30.709*rr**(-1.845)
          rt = 34.904 * exp(-2.067*log(rr))
          rq = 30.709 * exp(-1.845*log(rr))
      else if (rr.le.300.0) then
!         rt=1667.19*rr**(-2.907)
!         rq=1448.68*rr**(-2.682)
          rt = 1667.19 * exp(-2.907*log(rr))
          rq = 1448.68 * exp(-2.682*log(rr))
      else if (rr.le.1000.0) then
!         rt=5.88e5*rr**(-3.935)
!         rq=2.98e5*rr**(-3.616)
          rt = 5.88e5 * exp(-3.935*log(rr))
          rq = 2.98e5 * exp(-3.616*log(rr))
      endif

      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE PSI_COARE_SHEBA( id, zll, psi)
      
      use precision_mod
      implicit none
      
      integer :: id
      real(r8) :: zll, psi

      real(r8) :: x, y
      real(r8) :: psik, psic, f
 
      real(r8) :: am, bm, dm
      real(r8) :: ah, bh, ch, dh

      if (zll.lt.0.) then

            if(id.eq.1) then

               x = ( 1.0-16.0*zll )**0.25 
               y = ( 1.0-12.87*zll )**0.33333 

               psik = 2.d0*log((1.d0+x)/2.d0) + log((1.d0+x*x)/2.d0) - &
                      2.d0*atan(x) + 2.d0*atan(1.D0)

               psic = 1.5*log( (1.d0 + y + y*y )/3.d0 ) - &
                      sqrt(3.d0) * atan( ( 1.d0 + 2.d0*y )/sqrt(3.d0) ) + &
                      4.d0 * atan(1.d0)/sqrt(3.d0)

               f = 1.d0/( 1.d0 + zll*zll )
               psi = f * psik + (1.d0-f) * psic

            else             
               x = ( 1.d0-15.d0*zll )**.5
               y = ( 1.d0-34.15*zll )**0.33333  
!              
!               x = ( 1.d0-16.d0*zll )**.5
!               y = ( 1.d0-12.87*zll )**0.33333
!
               psik = 2.d0*log( (1.d0+x)/2.d0 )
               psic = 1.5*log( (1.d0 + y + y*y )/3.d0 ) - &
                      sqrt(3.d0) * atan( (1.d0+2.d0*y )/sqrt(3.d0) ) + &
                     4.d0 * atan(1.d0)/sqrt(3.d0)

                f = 1.d0/(1.d0+zll*zll)
                psi = f*psik + (1.0-f)*psic

            endif

      else if (zll.eq.0.) then

            psi=0.

      else if (zll.gt.0.) then
!
!              SHEBA
!
            am = 5.d0
            bm = am/6.5
            dm = ((1.d0-bm)/bm)**0.33333
            ah = 5.d0
            bh = 5.d0
            ch = 3.d0
            dh = sqrt(5.d0)
            
            if ( id.EQ.1 ) then

               x = (1.d0+zll)**0.33333
               psi = -3.d0*am/bm*(x-1.d0) + am*dm/(2.0d0*bm) * &
                     ( 2.d0*log((x+dm)/(1.d0+dm) ) - &
                     log( (x*x - x*dm + dm*dm )/( 1.0d0 - dm + dm*dm ) ) + &
                     2.d0 * sqrt(3.d0) * ( atan((2.d0*x-dm)/(sqrt(3.d0)*dm)) - &
                     atan( (2.d0-dm)/(sqrt(3.d0)*dm) ) ) )
            else
               psi = -0.5 * bh * log( 1.d0 + ch*zll + zll*zll ) + &
                     ( -ah/dh + 0.5*bh*ch/dh ) * &
                     ( log( (2.d0*zll + ch - dh )/(2.d0*zll + ch + dh ) ) - &
                       log( (ch-dh)/(ch+bh) ) )
            endif
      endif
      return
      end
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
      SUBROUTINE PSI_COARE25( ID, ZLL, PSI )

      use precision_mod
      implicit none
     
      integer :: ID
      real(r8) :: ZLL, PSI
      real(r8) :: CHIK, PSIK, F, CHIC, PSIC

      if ( ZLL.lt.0.0 ) then
           F = 1.D0/( 1.0D0 + ZLL*ZLL )
           CHIK = ( 1.D0-16.0D0*ZLL )**0.25

           if ( ID.eq.1 ) then
                PSIK = 2.0D0 * dlog( (1.D0+CHIK)/2.D0 ) + dlog( (1.D0 + CHIK*CHIK)/2.D0 ) - &
                       2.0D0 * ATAN(CHIK) + 2.0D0*ATAN(1.D0)
           else
                PSIK = 2.0D0 * dlog( (1.0D0 + CHIK*CHIK )/2.D0 )
           endif

           CHIC = ( 1.0D0-12.87*ZLL )**0.333    !for very unstable conditions

           PSIC = 1.5 * dlog( ( CHIC*CHIC + CHIC + 1.0D0 )/3.0D0 ) - &
                  ( 3.0D0**0.5 ) * ATAN( ( 2*CHIC + 1.0D0 )/(3.0D0**0.5) ) + &
                  4.0D0 * ATAN(1.0D0)/(3.0D0**0.5)
!
!         match Kansas and free-conv. forms with weighting F
!
           PSI = F*PSIK + ( 1.0D0-F )*PSIC

      else if ( ZLL.eq.0.0 ) then

           PSI = 0.0D0
      else
           PSI = -4.7 * ZLL

      endif
      
      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE PSI_COARE3( ID, ZLL, PSI )

      use precision_mod
      implicit none

      integer :: ID
      real(r8) :: ZLL, PSI
      real(r8) :: C
      real(r8) :: CHIK, PSIK, F, CHIC, PSIC


      IF ( ZLL.LT.0.0 ) THEN

           F = ZLL*ZLL/(1.0D0+ZLL*ZLL)

           IF ( ID.EQ.1 ) THEN
                CHIK = (1.D0-15.0D0*ZLL )**0.25
                PSIK = 2.0D0 * DLOG( (1.D0+CHIK)/2.D0 ) + DLOG( (1.D0+CHIK*CHIK)/2.D0 ) - &
                       2.0D0 * ATAN(CHIK) + 2.0D0*ATAN(1.D0)

                CHIC = ( 1.0D0-10.15*ZLL )**0.333                                              ! Convective
                PSIC = 1.5*DLOG( (CHIC*CHIC+CHIC+1.0D0)/3.0D0 ) - &
                       SQRT(3.0)*ATAN((2*CHIC+1.0D0)/SQRT(3.0)) + 4.0D0*ATAN(1.0D0)/SQRT(3.0)
           ELSE
                CHIK = (1.D0-15.0D0*ZLL )**0.5
                PSIK = 2.0D0 * DLOG( (1.0D0+CHIK)/2.D0 )

                CHIC = ( 1.0D0-34.15*ZLL )**0.333                                          !for very unstable conditions
                PSIC = 1.5*DLOG((CHIC*CHIC+CHIC+1.0D0)/3.0D0 ) - &
                       SQRT(3.0)*ATAN((2.0*CHIC+1.0D0)/SQRT(3.0)) + 4.0D0*ATAN(1.0D0)/SQRT(3.0)
           ENDIF

!
!         match Kansas and free-conv. forms with weighting F
!
           PSI = (1.0-F)*PSIK + F*PSIC

      ELSE IF ( ZLL.EQ.0.0 ) then

           PSI = 0.0D0
      ELSE
           IF ( ID.EQ.1 ) THEN
                C = MIN(50.0,0.35*ZLL)
                PSI = -( (1.0+1.0*ZLL)**1.0+0.6667*(ZLL-14.28)/DEXP(C)+8.525 )
           ELSE
                C = MIN(50.0,0.35*ZLL)
                PSI = -((1.0+2.0*ZLL/3.0)**1.5+0.6667*(ZLL-14.28)/DEXP(C)+8.525 )
           ENDIF

      ENDIF

      return
      end

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
