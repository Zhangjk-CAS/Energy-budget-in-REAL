#define LOGMSG()
!write(mytid+600,'(a,i4)')"LICOM",__LINE__
#define LOGMSGCarbon()
!write(600+mytid,'(a,i4)')"Licom CARBON",__LINE__

!  CVS: $Id: licom.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
 program licom
!    *.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*
!    +                                                               +
!    + ============================================================= +
!    +         LASG/IAP  CLIMATE OCEAN MODEL Version2.0 (LICOM2.0)
!    + ============================================================= +
!
!             A primitive-equation ocean model upgraded by the
!
!    +    State Key Laboratory of Numerical Modelling for            +
!    +    Atmospheric Sciences and Geophysical Fluid Dynamics (LASG) +
!    +    INSTITUTE OF ATMOSPHERIC PHYSICS (IAP)                     +
!    +    CHINESE ACADEMY OF SCIENCE (CAS)                           +
!    +    P.O. Box 9804, Beijing 100029, P.R.China                   +
!    +                                                               +
!    +                   December, 2002                              +
!    +                                                               +
!    +   AUTHORS: Hailong LIU  (lhl@lasg.iap.ac.cn)                  +
!    +            Xuehong ZHANG (zxh@lasg.iap.ac.cn)                  +
!    +            Yongqiang YU (yyq@lasg.iap.ac.cn)                  +
!    +                                                               +
!    +                                                               +
!
!
!       This model is based upon, but differs substantially from the
!                                     ~~~~~~~~~~~~~~~~~~~~~
!
!                        LASG/IAP L30T63 OGCM
!                   ==================================
!
!                        which was implemented by
!                      Xiangze Jin and Xuehong ZHANG
!
!
!             +------------------------------------------+
!             | DISTRIBUTION TERMS AND CONDITIONS NOTICE |
!             +------------------------------------------+
!
! (c) Copyright 2002 State Key Laboratory of Numerical Modelling for
! Atmospheric Sciences and Geophysical Fluid Dynamics (LASG) /
! Institute of Atmospheric Physics  (IAP)
!
! This software, the LASG/IAP Climate Ocean Model (LICOM), version 2.0 , was
! upgraded by State Key Laboratory of Numerical Modelling for
! Atmospheric Sciences and Geophysical Fluid Dynamics (LASG), Institute
! of Atmospheric Physics  (IAP), Chinese Academy of Sciences (CAS)
!
! Access and use of this software shall impose the following obligations
! and understandings on the user.  The user is granted the right,
! without any fee or cost, to use, copy, modify, alter, enhance and
! distribute this software, and any derivative works thereof, and its
! supporting documentation for any purpose whatsoever, except commercial
! sales, provided that this entire notice appears in all copies of the
! software, derivative works and supporting documentation.  Further, the
! user agrees to credit LASG/IAP in any publications that result
! from the use of this software or in any software package that includes
! this software.  The names LASG/IAP, however, may not be used in
! any advertising or publicity to endorse or promote any products or
! commercial entity unless specific written permission is obtained from
! LASG/IAP.
!
! THE LICOM MATERIALS ARE MADE AVAILABLE WITH THE UNDERSTANDING THAT
! UCAR/NCAR/CGD IS NOT OBLIGATED TO PROVIDE (AND WILL NOT PROVIDE) THE
! USER WITH ANY SUPPORT, CONSULTING, TRAINING, OR ASSISTANCE OF ANY KIND
! WITH REGARD TO THE USE, OPERATION AND PERFORMANCE OF THIS SOFTWARE, NOR
! TO PROVIDE THE USER WITH ANY UPDATES, REVISIONS, NEW VERSIONS, OR "BUG
! FIXES."
!
! THIS SOFTWARE IS PROVIDED BY LASG/IAP "AS IS" AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED.  IN NO EVENT SHALL UCAR/NCAR/CGD BE LIABLE FOR ANY
! SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER,
! INCLUDING BUT NOT LIMITED TO CLAIMS ASSOCIATED WITH THE LOSS OF DATA
! OR PROFITS, WHICH MAY RESULT FROM AN ACTION IN CONTRACT, NEGLIGENCE OR
! OTHER TORTIOUS CLAIM THAT ARISES OUT OF OR IN CONNECTION WITH THE
! ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
!
!-----------------------------------------------------------------------
!
! Purpose: Entry point for LICOM
!
!-----------------------------NOTICE------------------------------------
!
!
! Method: Call appropriate initialization, time-stepping, and finalization routines.
!
! Author: Hailong Liu and Yongqiang Yu, Dec. 31, 2002
!
!-----------------------------------------------------------------------

#include <def-undef.h>
use param_mod
use pconst_mod
#ifdef COUP
use shr_msg_mod
use shr_sys_mod
use control_mod
#endif
#if ( defined SPMD ) || ( defined COUP)
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
use tracer_mod
use pmix_mod
use forc_mod
#ifdef USE_OCN_CARBON
use carbon_mod
use cforce_mod
#endif
!
! Added by Dai in 2017.07
#if (defined OBCDT)
use openbc_mod
#endif
! by Dai in 2017.07
!
! DFS 2021-05
!
#if (defined OBTIDE) || (defined POTTIDE)
use tide_mod
#endif 
!
! DFS
!

      IMPLICIT NONE
#include <netcdf.inc>
     integer :: start_day
!
! DFS 2021-05
!
#if (defined OBTIDE) || (defined POTTIDE)
     integer :: iy0_tide
     integer :: imon_tide
     integer :: idays_tide
     integer :: i_tide

     real :: Tyls
     real :: Tday  
#endif
!
! DFS
!
! DFS 2024-07
!
#if ( defined TOGA_COARE_FLUX )
     integer :: idx_day
#endif
!
! DFS 2024-07     
!
!---------------------------------------------------------------------
!     Initilizing Message Passing
!---------------------------------------------------------------------
#ifdef SPMD
mpi_comm_ocn=0
#endif
#ifdef COUP
     call shr_msg_stdio('ocn')
     call msg_pass('connect')
#endif
!
      mytid=0
!
	LOGMSG()
#ifdef SPMD
#if ( ! defined COUP )
      call mpi_init(ierr)
      call mpi_comm_dup(mpi_comm_world, mpi_comm_ocn, ierr)
!      write(*,*)"COMM",mpi_comm_world, mpi_comm_ocn, ierr
#endif
	LOGMSG()
      call mpi_comm_rank (mpi_comm_ocn, mytid, ierr)
!      write(6,*)"LLLLLLLLLLLL"
      call mpi_comm_size (mpi_comm_ocn, nproc, ierr)
!      write(6,*) "MYTID=",mytid,"Number of Processors is",nproc
#endif
	LOGMSG()

!---------------------------------------------------------------------
!     SET THE CONSTANTS USED IN THE MODEL
!---------------------------------------------------------------------
#ifdef COUP
      call shr_sys_flush(6)
#endif
	LOGMSG()
      CALL CONST
	LOGMSG()
#ifdef COUP
      call shr_sys_flush(6)
#endif
#ifdef SHOW_TIME
      call run_time('CONST')
#endif

!***********************************************************************
!          SET SOME CONSTANTS FOR THE BOGCM
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL CTRLC
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('CTRLC')
#endif
#endif

!---------------------------------------------------------------------
!     SET MODEL'S RESOLUTION,TOPOGRAPHY AND THE CONSTANT
!     PARAMETERS RELATED TO LATITUDES (J)
!---------------------------------------------------------------------
	LOGMSG()
      CALL GRIDS
#ifdef SHOW_TIME
      call run_time('GRIDS')
#endif

!---------------------------------------------------------------------
!     SET SURFACE FORCING FIELDS (1: Annual mean; 0: Seasonal cycle)
!---------------------------------------------------------------------
	LOGMSG()
      CALL RDRIVER
#ifdef SHOW_TIME
      call run_time('RDRIVER')
#endif

!***********************************************************************
!      SET FORCING DATA USED IN CARBON CYCLYE
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
       CALL CFORCE
      LOGMSGCarbon()
#ifdef SHOW_TIME
       call run_time('CFORCE')
#endif
#endif

!---------------------------------------------------------------------
!     INITIALIZATION
!---------------------------------------------------------------------
	LOGMSG()
      CALL INIRUN
#ifdef SHOW_TIME
      call run_time('INRUN')
#endif

!
! DFS 2021-05
!
#if (defined OBTIDE)
     
      CALL Read_TPXO7

#endif
!
! DFS
!**********************************************************************
!      INITIALIZATION CARBON MODEL
!**********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL INIRUN_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INIRUN_PT')
#endif
#endif
#ifdef CANUTO
      call turb_ini
#endif
!---------------------------------------------------------------------
!     INITIALIZATION OF ISOPYCNAL MIXING
!---------------------------------------------------------------------
	LOGMSG()
#ifdef ISO
      CALL ISOPYI
#ifdef SHOW_TIME
      call run_time('ISOPYI')
#endif
#endif
!

      MEND = MONTH + NUMBER

      if (mytid==0) then
         WRITE (6,FMT='(A,I5,A,I5)') 'NUMBER=',NUMBER,',MONTH=',MONTH
!*************************************************************************
      IF(NSTART==1) THEN
        WRITE (6,FMT='(A)') 'The physical model is initial run!'
      ELSE
        WRITE (6,FMT='(A)') 'The physical model is continuous run!'
      ENDIF
#ifdef USE_OCN_CARBON
      WRITE (6,FMT='(A,I5,A,I5)') 'MONTHR=',monthR, &
                 ',yearR=', yearStart + (monthR-1)/12

      IF(NSTARTC==1) THEN
        WRITE (6,FMT='(A)') 'The passive-tracer model is initial run!'
      ELSE
        WRITE (6,FMT='(A)') 'The passive-tracer model is continuous run!'
      ENDIF

#endif
!*****************************************************************************

#ifdef COUP
         call shr_sys_flush(6)
#endif
      endif
      start_day = number_day
!
!
! DFS 2021-10
!
!#if (defined OBTIDE) || (defined POTTIDE)
!!
!      iy0_tide = ( MONTH-1 )/12
!      imon_tide = MONTH - iy0_tide*12
!
!!      idays_tide = 0
!      idays_tide = number_day
!      if ( imon_tide.gt.1 ) then
!           do i_tide = 1, imon_tide-1
!              idays_tide = idays_tide + NMONTH(i_tide)
!           enddo
!      endif
!       
!!      Tyls = 1979.0 + ( iy0_tide+1.0 )
!      Tyls = 1961.0 + ( iy0_tide+1.0 )
!      Tday = idays_tide*1.0
!
!      CALL ASTRONOMIC_NOD( Tyls, Tday ) 
!
!      Tintg_tide = 0.0D0
!#endif
!
! DFS
!
! DFS 2024-07
!
#if ( defined TOGA_COARE_FLUX )
      idx_day = 0

      idx_start = 1

      jwarm = 1
      jcool = 1
      jwave = 0
!
!         initialize warm layer variables and switches
!
      fxp = 0.5
      tk_pwp = 19.0
      jamset = 0
      jday1 = 1      
#endif
!
! DFS 2024-07
!
      loop1 : do

      IY0 = (MONTH -1)/12
      IYFM = IY0+1
!     IYFM is the number of the current year

      MON0 = MONTH - IY0*12
!      IMD = NMONTH (MON0)
!     IMD is the total day number of the current month
!

!
!  DFS 2024-05
!
      ! IYRUN = 1948 + ( IYFM - 1 ) + ( 42-28 )
      IYRUN = the_first_year + ( IYFM - 1 )   !zhangjk 20250417
      IF ( MOD(IYRUN,4).EQ.0 ) THEN
           IMD = N2MONTH(MON0)
      ELSE
           IMD = NMONTH (MON0)
      ENDIF
!
!  DFS 2024
!

      ISB = 0
      ISC = 0
      IST = 0
!     ISB/C/T : SWITCH ON EULER FORWARD OR LEAP-FROG SCHEME FOR
!     BAROTROPIC, BAROCLINIC, AND THERMOHALINE PROCESSES RESPECTIVELY

!**********************************************************************
#ifdef USE_OCN_CARBON
      yearR=(monthR-1)/12+1
      isp=0
#endif

!**********************************************************************
!
! Added by Dai in 2017.07
!#ifdef OBCDT
!      call openbc_pre
!#endif
! by Dai in 2017.07
!
!---------------------------------------------------------------------
!     THE CYCLE OF THE CURRENT MONTH
!---------------------------------------------------------------------
      loop2 : DO IDAY = start_day, IMD
!
! DFS 2021-10
!
#if (defined OBTIDE) || (defined POTTIDE)
!
      iy0_tide = ( MONTH-1 )/12
      imon_tide = MONTH - iy0_tide*12

!      idays_tide = number_day
      idays_tide = IDAY
      if ( imon_tide.gt.1 ) then
           do i_tide = 1, imon_tide-1
              idays_tide = idays_tide + NMONTH(i_tide)
           enddo
      endif
       
!      Tyls = 1961.0 + ( iy0_tide+1.0 )
      Tyls = IYRUN
      Tday = idays_tide*1.0

      CALL ASTRONOMIC_NOD( Tyls, Tday ) 

!      Tintg_tide = 0.0D0
#endif
!
! DFS
!
! Added by Dai in 2024.10
#ifdef OBCDT
           IF ( IDAY.EQ.start_day ) THEN
                call openbc_pre
           ENDIF
#endif
! by Dai in 2024.10
!
!
! DFS 2024-07
! 
#if ( defined TOGA_COARE_FLUX )
      idx_day = idx_day + 1
      if ( idx_day.gt.365 ) idx_day=3
#endif
!
! DFS 2024-07
!

#ifdef COUP
      if (nstart==2.and.iday==1) then
         nstart = 0
         cycle loop2
      end if
      cdate=((month-mon0)/12+1)*10000+mon0*100+iday
      sec=0
!
	LOGMSG()
      call msg_pass('send')
	LOGMSG()
      call msg_pass('recv')
	LOGMSG()
!
      call post_cpl
	LOGMSG()
#ifdef SPMD
      call mpi_bcast(stop_now,1,mpi_integer,0,mpi_comm_ocn,ierr)
#endif
      if (stop_now==1) exit loop1
#else
      if (month==mend) exit loop1
#endif
	LOGMSG()
      if (mytid==0) then
         OPEN (56,FILE ='modeltime',FORM ='FORMATTED',STATUS ='UNKNOWN')
         WRITE (56,FMT='(A,A,I4,I8)') ' It is running on ', ABMON (MON0),IDAY,1+ IY0
!**********************************************************************
#ifdef USE_OCN_CARBON
         WRITE (56,FMT='(A,I4,A,I8)') ' The ptracer model is running on ', &
                                      IDAY, ABMON (monthR-(yearR-1)*12),yearStart+yearR-1
#endif
!*********************************************************************
         CLOSE (56)
      endif
	LOGMSG()
!
! DFS 2024-08
!
#ifndef OBCDT

!---------------------------------------------------------------------
!     INTERPOLATE THE OBSERVED MONTHLY MEAN DATA TO CERTAIN DAY
!---------------------------------------------------------------------
	LOGMSG()
#if (!defined COUP)
      CALL INTFOR
!         CALL INTFOR (IYFM,MON0)
#endif
!
!*********************************************************************
!   INTERPOLATE THE FORCING DATA TO CERTAIN DAY
!*********************************************************************
#ifdef USE_OCN_CARBON
#if (!defined COUP)
      LOGMSGCarbon()
      CALL INTFORW
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INTFORW')
#endif
#endif
      LOGMSGCarbon()
      CALL INTFOR_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INTFOR_PT')
#endif
#endif

	LOGMSG()
!
!
#endif
!
! DFS 2024-08
!
#if (defined SOLARCHLORO)
       call sw_absor
! compute the short wave pentration dependent on chlorophy
#endif

#ifdef  SHOW_TIME
         call run_time('INTFOR')
#endif
!
!---------------------------------------------------------------------
!     THERMAL CYCLE
!---------------------------------------------------------------------

!zhangjk 20250417
#ifdef DAILYMEAN_OP
call mm00(km)
#endif 

         DO II = 1,NSS
!
!
! DFS 2024-07
!
#if ( defined TOGA_COARE_FLUX )
           IF ( idx_day.EQ.1 ) THEN
                IF ( II.EQ.1 ) THEN
                     idx_start = 1
                ELSE
                     idx_start = 2
                ENDIF
           ELSE 
                idx_start = 2
           ENDIF
#endif
!
! DFS 2024-07
!
!
!lhl20110731
           CALL CORE_DAILY(II)
!lhl20110731
!
!     COMPUTE DENSITY, BAROCLINIC PRESSURE AND THE RELAVANT VARIABLES
	LOGMSG()
            CALL READYT
#ifdef  SHOW_TIME
            call run_time('READYT')
#endif
!
!
!---------------------------------------------------------------------
!     BAROCLINIC & BAROTROPIC CYCLE
!---------------------------------------------------------------------
            DO JJ = 1,NCC

!     COMPUTE MOMENTUM ADVECTION, DIFFUSION & THEIR VERTICAL INTEGRALS
	LOGMSG()
               CALL READYC
 	LOGMSG()
#ifdef  SHOW_TIME
               call run_time('READYC')
#endif
!
!
!     PREDICTION OF BAROTROPIC MODE
	LOGMSG()
               CALL BAROTR
#ifdef  SHOW_TIME
               call run_time('BAROTR')
#endif

!     PREDICTION OF BAROCLINIC MODE
	LOGMSG()
               CALL BCLINC
#ifdef  SHOW_TIME
               call run_time('BCLINC')
#endif
!lhl
#ifdef  DEBUG
	LOGMSG()
#ifdef  SHOW_TIME
#endif
#endif
!lhl

            END DO

	LOGMSG()
!*******************************************************************
!     PREDICTION OF PASSIVE TRACER
!*******************************************************************
#ifdef USE_OCN_CARBON
#ifdef carbonDebug
      print*, "II in licom======================================================================",II
#endif
#ifdef printcall
#ifdef SPMD
      print*,"call ptracer in licom, mytid=",mytid
#else
      print*,"call ptracer in licom"
#endif
#endif
      LOGMSGCarbon()
      CALL NEWPTRACER
      LOGMSGCarbon()

#ifdef SHOW_TIME
      call run_time('PTRACER')
#endif
#endif
!

!      CALL ENERGY
	LOGMSG()
            CALL TRACER
	LOGMSG()
#ifdef  SHOW_TIME
            call run_time('TRACER')
#endif
      CALL ENERGY
	LOGMSG()
            CALL ICESNOW
#ifdef  SHOW_TIME
            call run_time('ICESNOW')
#endif

!***********************************************************************
!     PERFORM CONVECTIVE ADJUSTMENT IF UNSTABLE STRATIFICATION OCCURS
!************************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
            CALL CONVADJ_PT
      LOGMSGCarbon()
#else
!lhl1204
!#if (!defined CANUTO)
            CALL CONVADJ
!#endif
#endif

#ifdef  SHOW_TIME
            call run_time('CONVADJ')
#endif

!zhangjk 20250417
#ifdef DAILYMEAN_OP
!     ACCUMULATE SOME VARIABLES FOR !!daily!! OUTPUT
	LOGMSG()
         CALL ACCUMM
#ifdef  SHOW_TIME
         call run_time('ACCUMM')
#endif
#endif 

!add hourly output by zhangjk 20250610 
#ifdef HOURLY_OP 
call sSAVE_hourly(ii)
#endif 

!*************************************************************************
!
         END DO
      
         CALL ENERGY

!     COMPENSATE THE LOSS OF GROSS MASS
	LOGMSG()
!
!    The Following line is canceled by DFS in 2017-08
!         CALL ADDPS
!    DFS in2017-08
!
#ifndef OBCDT
        CALL ADDPS
#endif 
!
#ifdef  SHOW_TIME
         call run_time('ADDPS')
#endif
!
	LOGMSG()
#ifdef  SHOW_TIME
         call run_time('ENERGY')
#endif
!
!     MONITOR THE MODEL INTEGRATION

!the accumm canceled by zhangjk 20250417
#ifndef DAILYMEAN_OP
!     ACCUMULATE SOME VARIABLES FOR MONTHLY OUTPUT
	LOGMSG()
         CALL ACCUMM
#ifdef  SHOW_TIME
         call run_time('ACCUMM')
#endif
#endif
!*********************************************************************
!      ACCUMULATE SOME VARIABLES IN CARBON MODEL
!*********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL ACCUMM_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('ACCUMM_PT')
#endif
#endif

#ifdef COUP
         call shr_sys_flush(6)
#endif
!
	LOGMSG()
#ifdef COUP
         call flux_cpl
#endif

      CALL SSAVECDF
      END DO loop2

	LOGMSG()
      start_day = 1
      MONTH = MONTH +1

!---------------------------------------------------------------------
!     SAVE, SMOOTH & RESET SOME ARRAYS
!---------------------------------------------------------------------
	LOGMSG()
!     CALL SSAVECDF
#ifdef  SHOW_TIME
      call run_time('SSAVECDF')
#endif
!*********************************************************************
!     SAVE, SMOOTH & RESET SOME ARRAYS IN CARBON MODEL
!*********************************************************************
#ifdef USE_OCN_CARBON
#ifdef printcall
#ifdef SPMD
      print*,"call ssave_pt in licom, mytid=",mytid
#else
      print*,"call ssave_pt in licom"
#endif
#endif
      CALL SSAVE_PT
#ifdef SHOW_TIME
      call run_time('SSAVE_PT')
#endif
#endif
      end do loop1

#ifdef COUP
      call msg_pass('disconnect')
#else
#ifdef SPMD
      call mpi_finalize(ierr)
#endif
#endif
      close(6)
!
    deallocate(su3,sv3,psa3,tsa3,qar3,uva3,swv3,cld3,sss3,sst3 ,nswv3,dqdt3,chloro3)
    deallocate(seaice3,runoff3)
    deallocate(wspd3,wspdu3,wspdv3,lwv3,rain3,snow3) 

!
      STOP
      END  program LICOM
