
      SUBROUTINE openbc_pre

#include <def-undef.h>

!      use param_mod, only: imt_global, jmt_global, &
!                           imt, jmt, km, i, j, k, mytid
!      USE pconst_mod, only: IYFM, MON0, dzp 

      use param_mod
      USE pconst_mod
 
      use openbc_mod
      use PRECISION_MOD
      use forc_mod, only: sst, sss

#if ( defined SPMD )
      use msg_mod
#endif

      IMPLICIT NONE
 
#include <netcdf.inc>

!      INTEGER :: TNUM

      INTEGER ::   IDYEAR, IDMON

      INTEGER ::   ID_CORE_UBVB
      INTEGER ::   ID_UBVB_OUT

      REAL(r4) :: BUFF_GL(imt_global_gl)
      REAL(r4) :: BUFF_GL_SN(jmt_global_gl)
      REAL(r4) :: spvalue

      CHARACTER :: CHMON*2
      CHARACTER :: CHYEAR*4
      CHARACTER :: FNAMEZ*13
      CHARACTER :: FNAMET*13
      CHARACTER :: FNAMES*13
      CHARACTER :: FNAMEU*13
      CHARACTER :: FNAMEV*13

      CHARACTER :: FNAMEOUT*28

      INTEGER :: iret, ncid, varid
      integer*4,  dimension(4) :: startd(4)
      integer*4,  dimension(4) :: countd(4)

      INTEGER :: JGLSTART, JGLEND
      INTEGER :: IGLSTART, IGLEND

      INTEGER :: IL, JL

      INTEGER :: ncid_rg
      INTEGER :: ind_id_rg
!
!      INTEGER, DIMENSION(imt_global) :: ind_rg
!      INTEGER, DIMENSION(jmt_global) :: ind_rg_sn
!
!      INTEGER :: ind_id_rg_2d
!      INTEGER, DIMENSION(imt_global,2) :: ind_rg_2d
!      INTEGER, DIMENSION(2,jmt_global) :: ind_rg_sn_2d
!
!      REAL(r4) :: BUFF_GL_2D(imt_global_gl,2)
!      REAL(r4) :: BUFF_GL_SN_2D(2,jmt_global_gl)
!
!      real, DIMENSION(imt_global,2,km) :: find_rg
!      real, DIMENSION(2,jmt_global,km) :: find_rg_sn
!
      integer*4,  dimension(3) :: start2(3)
      integer*4,  dimension(3) :: count2(3)

      REAL(r4), allocatable, dimension(:,:) :: BUFF_GL_GL
      INTEGER, DIMENSION(imt_global,jmt_global) :: ind_rggl_2d
      INTEGER :: ind_id_rggl

      REAL(r8), allocatable, dimension(:,:) :: UV_IJMT_GL
      REAL(r8), allocatable, dimension(:,:) :: UB_IJMT_GL
      REAL(r8), allocatable, dimension(:,:) :: VB_IJMT_GL

      REAL(r8) :: TMPH_GL(imt_global,jmt_global)
      REAL(r8) :: TMPVAR_GL(imt_global,jmt_global)
      REAL(r8) :: TMPH
 
      REAL(r8), allocatable, dimension(:,:) :: TT_IJMT_GL
      REAL(r8), allocatable, dimension(:,:) :: SS_IJMT_GL

      REAL(r8), allocatable, dimension(:,:) :: H0_IJMT_GL
!
!      REAL(r8), dimension(imt,jmt) :: H0OBCELOC
!      REAL(r8), dimension(imt,jmt) :: H0OBCWLOC
!      REAL(r8), dimension(imt,jmt) :: H0OBCSLOC
!      REAL(r8), dimension(imt,jmt) :: H0OBCNLOC
!
!      REAL(r8), dimension(imt,jmt,km) :: TTOBCELOC
!      REAL(r8), dimension(imt,jmt,km) :: TTOBCWLOC
!      REAL(r8), dimension(imt,jmt,km) :: TTOBCSLOC
!      REAL(r8), dimension(imt,jmt,km) :: TTOBCNLOC
!
!      REAL(r8), dimension(imt,jmt,km) :: SSOBCELOC
!      REAL(r8), dimension(imt,jmt,km) :: SSOBCWLOC
!      REAL(r8), dimension(imt,jmt,km) :: SSOBCSLOC
!      REAL(r8), dimension(imt,jmt,km) :: SSOBCNLOC
!
!      REAL(r8), dimension(imt,jmt,km) :: UUOBCELOC
!      REAL(r8), dimension(imt,jmt,km) :: UUOBCWLOC
!      REAL(r8), dimension(imt,jmt,km) :: UUOBCSLOC
!      REAL(r8), dimension(imt,jmt,km) :: UUOBCNLOC
!
!      REAL(r8), dimension(imt,jmt,km) :: VVOBCELOC
!      REAL(r8), dimension(imt,jmt,km) :: VVOBCWLOC
!      REAL(r8), dimension(imt,jmt,km) :: VVOBCSLOC
!      REAL(r8), dimension(imt,jmt,km) :: VVOBCNLOC
!
!      REAL(r8), dimension(imt,jmt) :: UBOBCELOC
!      REAL(r8), dimension(imt,jmt) :: UBOBCWLOC
!      REAL(r8), dimension(imt,jmt) :: UBOBCSLOC
!      REAL(r8), dimension(imt,jmt) :: UBOBCNLOC
!
!      REAL(r8), dimension(imt,jmt) :: VBOBCELOC
!      REAL(r8), dimension(imt,jmt) :: VBOBCWLOC
!      REAL(r8), dimension(imt,jmt) :: VBOBCSLOC
!      REAL(r8), dimension(imt,jmt) :: VBOBCNLOC
!
!
!
      spvalue = 1.e+35

      JGLSTART = int((90-STARTLAT)*10)+1
      JGLEND = JGLSTART + jmt_global - 1

      IGLSTART = int(STARTLON*10) + 1
      IGLEND = IGLSTART + imt_global - 1


      
      IDYEAR = IYFM-1 + the_first_year
      WRITE(CHYEAR,'(I4.4)') IDYEAR

      IDMON = MON0
      WRITE(CHMON,'(I2.2)') IDMON
!
!      DO K = 1, km
!         DO IL = 1, imt_global
!            DO JL = 1, jmt
!               TTOBCSGL(IL,JL,K) = 0.0D0
!               SSOBCSGL(IL,JL,K) = 0.0D0
!               TTOBCNGL(IL,JL,K) = 0.0D0
!               SSOBCNGL(IL,JL,K) = 0.0D0
!               UUOBCSGL(IL,JL,K) = 0.0D0
!               VVOBCSGL(IL,JL,K) = 0.0D0
!               UUOBCNGL(IL,JL,K) = 0.0D0
!               VVOBCNGL(IL,JL,K) = 0.0D0
!            ENDDO 
!         ENDDO
!         DO JL = 1, jmt_global
!            DO IL = 1, imt
!               TTOBCEGL(IL,JL,K) = 0.0D0
!               SSOBCEGL(IL,JL,K) = 0.0D0
!               TTOBCWGL(IL,JL,K) = 0.0D0
!               SSOBCWGL(IL,JL,K) = 0.0D0
!               UUOBCEGL(IL,JL,K) = 0.0D0
!               VVOBCEGL(IL,JL,K) = 0.0D0
!               UUOBCWGL(IL,JL,K) = 0.0D0
!               VVOBCWGL(IL,JL,K) = 0.0D0
!            ENDDO
!         ENDDO
!         
!      ENDDO
!
!      DO IL = 1, imt_global
!         DO JL = 1, jmt
!            UBOBCSGL(IL,JL) = 0.0D0
!            UBOBCNGL(IL,JL) = 0.0D0
!            VBOBCSGL(IL,JL) = 0.0D0
!            VBOBCNGL(IL,JL) = 0.0D0
!         
!            H0OBCSGL(IL,JL) = 0.0D0
!            H0OBCNGL(IL,JL) = 0.0D0
!         ENDDO
!      ENDDO
!
!      DO JL = 1, jmt_global
!         DO IL = 1, imt
!            UBOBCEGL(IL,JL) = 0.0D0
!            UBOBCWGL(IL,JL) = 0.0D0
!            VBOBCEGL(IL,JL) = 0.0D0
!            VBOBCWGL(IL,JL) = 0.0D0
!
!            H0OBCEGL(IL,JL) = 0.0D0
!            H0OBCWGL(IL,JL) = 0.0D0
!         ENDDO
!      ENDDO     
!
!!$OMP PARALLEL DO PRIVATE (K,J,I)
!      DO K = 1, KM
!         DO J = 1, jmt
!         DO I = 1, imt
!            TTOBCSLOC(I,J,K) = 0.D0
!            SSOBCSLOC(I,J,K) = 0.D0
!            TTOBCNLOC(I,J,K) = 0.D0
!            SSOBCNLOC(I,J,K) = 0.D0
!         ENDDO
!         ENDDO
!      ENDDO
!
!!$OMP PARALLEL DO PRIVATE (K,J,I)
!      DO K =1, KM
!         DO J = 1, jmt
!         DO I = 1, imt
!            TTOBCELOC(I,J,K) = 0.D0
!            SSOBCELOC(I,J,K) = 0.D0
!            TTOBCWLOC(I,J,K) = 0.D0
!            SSOBCWLOC(I,J,K) = 0.D0
!         ENDDO
!         ENDDO
!      ENDDO
!
!!$OMP PARALLEL DO PRIVATE (K,J,I)
!      DO K = 1, KM
!         DO J = 1, jmt
!            DO I = 1, imt
!               UUOBCSLOC(I,J,K) = 0.D0
!               VVOBCSLOC(I,J,K) = 0.D0
!               UUOBCNLOC(I,J,K) = 0.D0
!               VVOBCNLOC(I,J,K) = 0.D0
!            ENDDO
!         ENDDO
!      ENDDO
!
!!$OMP PARALLEL DO PRIVATE (K,J,I)
!      DO K = 1, KM
!         DO J = 1, jmt
!            DO I = 1, imt
!               UUOBCELOC(I,J,K) = 0.D0
!               VVOBCELOC(I,J,K) = 0.D0
!               UUOBCWLOC(I,J,K) = 0.D0
!               VVOBCWLOC(I,J,K) = 0.D0
!            ENDDO
!         ENDDO
!      ENDDO
!
!!$OMP PARALLEL DO PRIVATE (J,I)
!      DO J = 1, jmt
!         DO I = 1, imt
!            H0OBCSLOC(I,J) = 0.D0
!            H0OBCNLOC(I,J) = 0.D0
!            UBOBCSLOC(I,J) = 0.D0
!            UBOBCNLOC(I,J) = 0.D0
!            VBOBCSLOC(I,J) = 0.D0
!            VBOBCNLOC(I,J) = 0.D0
!         ENDDO
!      ENDDO
!
!!$OMP PARALLEL DO PRIVATE (I,J) 
!      DO I = 1, imt
!         DO J = 1, jmt
!            H0OBCELOC(I,J) = 0.D0
!            H0OBCWLOC(I,J) = 0.D0
!            UBOBCELOC(I,J) = 0.D0
!            UBOBCWLOC(I,J) = 0.D0
!            VBOBCELOC(I,J) = 0.D0
!            VBOBCWLOC(I,J) = 0.D0
!         ENDDO
!      ENDDO
!

!$OMP PARALLEL DO PRIVATE (I,J) 
      DO I = 1, IMT
      DO J = 1, JMT
         H0_MON_LOC(I,J) = 0.0D0
         UB_MON_LOC(I,J) = 0.0D0
         VB_MON_LOC(I,J) = 0.0D0
      ENDDO
      ENDDO

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1, KM 
         DO J = 1, JMT
         DO I = 1, IMT
            UU_MON_LOC(I,J,K) = 0.0D0
            VV_MON_LOC(I,J,K) = 0.0D0
            TT_MON_LOC(I,J,K) = 0.0D0
            SS_MON_LOC(I,J,K) = 0.0D0
         ENDDO
         ENDDO
      ENDDO
!     
         if ( allocated(BUFF_GL_GL) ) deallocate(BUFF_GL_GL)
         allocate( BUFF_GL_GL(imt_global_gl,jmt_global_gl) )

         if ( allocated( H0_IJMT_GL ) ) deallocate(H0_IJMT_GL)
         allocate( H0_IJMT_GL(imt_global,jmt_global) )

         FNAMET="z0-"//CHYEAR//"-"//CHMON//".nc"
!
!             Read Temperature From LICOM Global Products For Southern BC 
!
         if ( mytid==0 ) then

         iret = nf_open('./boundary-data-glorys/'//FNAMET, nf_nowrite, ncid)
         call check_err(iret)

         iret = nf_open('ind_final_eas.nc', nf_nowrite, ncid_rg)
         call check_err(iret)

            startd(1)=1; countd(1)=imt_global_gl
            startd(2)=1; countd(2)=jmt_global_gl
            startd(3)=1; countd(3)=1
            startd(4)=1; countd(4)=1
            ! start2(1)=1; count2(1)=imt_global_gl
            ! start2(2)=1; count2(2)=jmt_global_gl
            ! start2(3)=1; count2(3)=1
            iret = nf_inq_varid(ncid, 'ssh', varid)
            call check_err(iret)
            iret = nf_get_vara_real(ncid, varid, startd, countd, &
                                    BUFF_GL_GL)
            call check_err(iret)

            start2(1)=1; count2(1)=imt_global
            start2(2)=1; count2(2)=jmt_global
            start2(3)=1; count2(3)=1
            iret = nf_inq_varid(ncid_rg, 'ind', ind_id_rggl)
            call check_err(iret)
            iret = nf_get_vara_int( ncid_rg, ind_id_rggl, start2, &
                                    count2,ind_rggl_2d)
            call check_err(iret)

            DO JL = 1, jmt_global
               DO IL=1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                          spvalue ) THEN
                       H0_IJMT_GL(IL,JL) = 1.D0*BUFF_GL_GL(IL+ &
                        (IGLSTART-1),JL+(JGLSTART-1))*ind_rggl_2d(IL,JL)
                  ELSE
                       H0_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

            DO JL = 1, 2
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                          spvalue ) THEN
                       H0_IJMT_GL(IL,JL) = 1.0D0* &
                           BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       H0_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO JL = jmt_global-1, jmt_global
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                          spvalue ) THEN
                       H0_IJMT_GL(IL,JL) = 1.0D0* &
                          BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       H0_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = 1, 2
               DO JL = 3, jmt_global-2
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                          spvalue ) THEN
                       H0_IJMT_GL(IL,JL) = 1.0D0* &
                         BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       H0_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = imt_global-1, imt_global
               DO JL = 3, jmt_global-2
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                          spvalue ) THEN
                       H0_IJMT_GL(IL,JL) = 1.0D0* &
                          BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       H0_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

!
!            DO J = 1, jmt
!               DO I = 1, imt_global
!                  H0OBCSGL(I,J) = H0_IJMT_GL(I,jmt_global-jmt+J)
!                  H0OBCNGL(I,J) = H0_IJMT_GL(I,J)
!               ENDDO
!            ENDDO
!
!            DO I = 1, imt
!               DO J = 1, jmt_global
!                  H0OBCEGL(I,J) = H0_IJMT_GL(imt_global-imt+I,J)
!                  H0OBCWGL(I,J) = H0_IJMT_GL(I,J)
!               ENDDO
!            ENDDO
!
         iret = nf_close(ncid)
         call check_err(iret)
         iret = nf_close(ncid_rg)
         call check_err(iret)

         endif

         CALL global_distribute( H0_IJMT_GL, H0_MON_LOC )

!         CALL global_to_local_2d( H0OBCEGL, H0OBCELOC )
!         CALL global_to_local_2d( H0OBCWGL, H0OBCWLOC )
!
!         CALL distribute_2d_sn( H0OBCSGL, H0OBCSLOC )
!         CALL distribute_2d_sn( H0OBCNGL, H0OBCNLOC )
!
!!$OMP PARALLEL DO PRIVATE (J)
!         DO J = 1, jmt
!            H0BTRE(1,J) = H0OBCELOC(IMT-3,J)
!            H0BTRE(2,J) = H0OBCELOC(IMT-2,J)
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (J)
!         DO J = 1, jmt
!            H0BTRW(1,J) = H0OBCWLOC(3,J)
!            H0BTRW(2,J) = H0OBCWLOC(4,J)
!         ENDDO
!
!!$OMP PARALLEL DO PRIVATE (I)
!         DO I = 1, imt
!            H0BTRS(I,1) = H0OBCSLOC(I,JMT-3)
!            H0BTRS(I,2) = H0OBCSLOC(I,JMT-2)
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (I)
!         DO I = 1, imt
!            H0BTRN(1,I) = H0OBCNLOC(I,3)
!            H0BTRN(2,I) = H0OBCNLOC(I,4)
!         ENDDO
!        
         deallocate( H0_IJMT_GL)
         deallocate( BUFF_GL_GL)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
         if ( allocated(BUFF_GL_GL) ) deallocate(BUFF_GL_GL)
         allocate( BUFF_GL_GL(imt_global_gl,jmt_global_gl) )

         if ( allocated( TT_IJMT_GL ) ) deallocate(TT_IJMT_GL)
         allocate( TT_IJMT_GL(imt_global,jmt_global) )

         FNAMET="tt-"//CHYEAR//"-"//CHMON//".nc"
!
!             Read Temperature From LICOM Global Products For Southern BC 
!

         if ( mytid==0 ) then
         iret = nf_open('./boundary-data-glorys/'//FNAMET, nf_nowrite, ncid)
         call check_err(iret)
         iret = nf_open('ind_final_eas.nc', nf_nowrite, ncid_rg)
         call check_err(iret)
         endif

         DO K =1, km

            if ( mytid == 0 ) then

            startd(1)=1; countd(1)=imt_global_gl
            startd(2)=1; countd(2)=jmt_global_gl
            startd(3)=k; countd(3)=1
            startd(4)=1; countd(4)=1

            iret = nf_inq_varid(ncid, 'tt', varid)
            call check_err(iret)
            iret = nf_get_vara_real(ncid, varid, startd, countd, &
                                    BUFF_GL_GL)
            call check_err(iret)

            start2(1)=1; count2(1)=imt_global
            start2(2)=1; count2(2)=jmt_global
            start2(3)=k; count2(3)=1
            iret = nf_inq_varid(ncid_rg, 'ind', ind_id_rggl)
            call check_err(iret)
            iret = nf_get_vara_int( ncid_rg, ind_id_rggl, start2, &
                                    count2,ind_rggl_2d)
            call check_err(iret)

            DO JL = 1, jmt_global
               DO IL=1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                              spvalue ) THEN
                       TT_IJMT_GL(IL,JL) = 1.0D0*BUFF_GL_GL(IL+ &
                        (IGLSTART-1),JL+(JGLSTART-1))*ind_rggl_2d(IL,JL)
                  ELSE
                       TT_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

            DO JL = 1, 2
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                              spvalue ) THEN
                       TT_IJMT_GL(IL,JL) = 1.0D0* &
                         BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       TT_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO JL = jmt_global-1, jmt_global
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                              spvalue ) THEN
                       TT_IJMT_GL(IL,JL) = 1.0D0* &
                          BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       TT_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = 1, 2
               DO JL = 3, jmt_global-2
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                              spvalue ) THEN
                       TT_IJMT_GL(IL,JL) = 1.0D0* &
                          BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       TT_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = imt_global-1, imt_global
               DO JL = 3, jmt_global-2
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                              spvalue ) THEN
                       TT_IJMT_GL(IL,JL) = 1.0D0* &
                          BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       TT_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

!            DO I = 1, imt_global
!            DO J = 1, jmt
!               TTOBCSGL(I,J,K) = TT_IJMT_GL(I,jmt_global-jmt+J)
!               TTOBCNGL(I,J,K) = TT_IJMT_GL(I,J)
!            ENDDO
!            ENDDO
!            
!            DO I = 1, imt
!            DO J = 1, jmt_global
!               TTOBCEGL(I,J,K) = TT_IJMT_GL(imt_global-imt+I,j)
!               TTOBCWGL(I,J,K) = TT_IJMT_GL(I,J)
!            ENDDO
!            ENDDO
!
            endif

            CALL global_distribute( TT_IJMT_GL, TT_MON_LOC(:,:,K))

         ENDDO
!zhangjk 250510
         SST(:,:) = TT_MON_LOC(:,:,1)

!         CALL global_to_local_3d( TTOBCEGL, TTOBCELOC, KM)
!         CALL global_to_local_3d( TTOBCWGL, TTOBCWLOC, KM)
!
!         CALL distribute_3d_sn( TTOBCSGL, TTOBCSLOC, KM)
!         CALL distribute_3d_sn( TTOBCNGL, TTOBCNLOC, KM)
!
!!$OMP PARALLEL DO PRIVATE (K,J)
!         DO K = 1, KM
!         DO J = 1, jmt
!            TBCLE(J,K) = TTOBCELOC(IMT-2,J,K)
!            TBCLW(J,K) = TTOBCWLOC(3,J,K)
!         ENDDO
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (K,I)
!         DO K = 1, KM
!         DO I = 1, imt
!            TBCLS(I,K) = TTOBCSLOC(I,JMT-2,K)
!            TBCLN(I,K) = TTOBCNLOC(I,3,K)
!         ENDDO
!         ENDDO
!      
         if ( mytid == 0 ) then

         iret = nf_close(ncid)
         call check_err(iret)
         iret = nf_close(ncid_rg)
         call check_err(iret)

         endif

         deallocate( TT_IJMT_GL )
!
!
!
         if ( allocated( SS_IJMT_GL ) ) deallocate(SS_IJMT_GL)
         allocate( SS_IJMT_GL(imt_global,jmt_global) )

         FNAMET="ss-"//CHYEAR//"-"//CHMON//".nc"
!
!             Read Temperature From LICOM Global Products For Southern BC 
!

         if ( mytid == 0 ) then

         iret = nf_open('./boundary-data-glorys/'//FNAMET, nf_nowrite, ncid)
         call check_err(iret)

         iret = nf_open('ind_final_eas.nc', nf_nowrite, ncid_rg)
         call check_err(iret)

         endif

         DO K =1, km

            if ( mytid == 0 ) then

            startd(1)=1; countd(1)=imt_global_gl
            startd(2)=1; countd(2)=jmt_global_gl
            startd(3)=k; countd(3)=1
            startd(4)=1; countd(4)=1

            iret = nf_inq_varid(ncid, 'ss', varid)
            call check_err(iret)
            iret = nf_get_vara_real(ncid, varid, startd, countd, &
                                    BUFF_GL_GL)
            call check_err(iret)

            start2(1)=1; count2(1)=imt_global
            start2(2)=1; count2(2)=jmt_global
            start2(3)=k; count2(3)=1
            iret = nf_inq_varid(ncid_rg, 'ind', ind_id_rggl)
            call check_err(iret)
            iret = nf_get_vara_int( ncid_rg, ind_id_rggl, start2, &
                                    count2,ind_rggl_2d)
            call check_err(iret)

            DO JL = 1, jmt_global
               DO IL=1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                      SS_IJMT_GL(IL,JL) = 1.0D0* &
                       ( (BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))- &
                                   35.0D0)*0.001D0 )*ind_rggl_2d(IL,JL)
                  ELSE
                       SS_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

            DO JL = 1, 2
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       SS_IJMT_GL(IL,JL) = 1.0D0* &
                        ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))- &
                          35.0D0 )*0.001D0
                  ELSE
                       SS_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO JL = jmt_global-1, jmt_global
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                      SS_IJMT_GL(IL,JL) = 1.0D0* &
                       ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))- &
                         35.0D0 )*0.001D0
                  ELSE
                       SS_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = 1, 2
               DO JL = 3, jmt_global-2
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       SS_IJMT_GL(IL,JL) = 1.0D0* &
                        ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))- &
                          35.0D0 )*0.001D0
                  ELSE
                       SS_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = imt_global-1, imt_global
               DO JL = 3, jmt_global-2
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       SS_IJMT_GL(IL,JL) = 1.0D0* &
                        ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))- &
                          35.0D0 )*0.001D0
                  ELSE
                       SS_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

!            DO I = 1, imt_global
!            DO J = 1, jmt
!               SSOBCSGL(I,J,K) = SS_IJMT_GL(I,jmt_global-jmt+J)
!               SSOBCNGL(I,J,K) = SS_IJMT_GL(I,J)
!            ENDDO
!            ENDDO
!            
!            DO I = 1, imt
!            DO j = 1, jmt_global
!               SSOBCEGL(I,J,K) = SS_IJMT_GL(imt_global-imt+I,j)
!               SSOBCWGL(I,J,K) = SS_IJMT_GL(I,J)
!            ENDDO
!            ENDDO
!
            endif

            CALL global_distribute( SS_IJMT_GL, SS_MON_LOC(1,1,K) )

         ENDDO
!zhangjk 250510
         SSS(:,:) = SS_MON_LOC(:,:,1)

!         CALL global_to_local_3d( SSOBCEGL, SSOBCELOC, KM)
!         CALL global_to_local_3d( SSOBCWGL, SSOBCWLOC, KM)
!
!         CALL distribute_3d_sn( SSOBCSGL, SSOBCSLOC, KM)
!         CALL distribute_3d_sn( SSOBCNGL, SSOBCNLOC, KM)
!
!!$OMP PARALLEL DO PRIVATE (K,J)
!         DO K = 1, KM
!         DO J = 1, jmt
!            SBCLE(J,K) = SSOBCELOC(IMT-2,J,K)
!            SBCLW(J,K) = SSOBCWLOC(3,J,K)
!         ENDDO
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (K,I)
!         DO K = 1, KM
!         DO I = 1, imt
!            SBCLS(I,K) = SSOBCSLOC(I,JMT-2,K)
!            SBCLN(I,K) = SSOBCNLOC(I,3,K)
!         ENDDO
!         ENDDO
!     

         if ( mytid == 0 ) then

         iret = nf_close(ncid)
         call check_err(iret)
         iret = nf_close(ncid_rg)
         call check_err(iret)

         endif

         deallocate(SS_IJMT_GL)
         deallocate(BUFF_GL_GL)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
         if ( allocated(BUFF_GL_GL) ) deallocate(BUFF_GL_GL)
         allocate( BUFF_GL_GL(imt_global_gl,jmt_global_gl) )

         if ( allocated( UV_IJMT_GL ) ) deallocate(UV_IJMT_GL)
         allocate( UV_IJMT_GL(imt_global,jmt_global) )
        
         if ( allocated( UB_IJMT_GL ) ) deallocate(UB_IJMT_GL)
         allocate( UB_IJMT_GL(imt_global,jmt_global) )
        
         DO JL = 1, jmt_global
         DO IL = 1, imt_global
            TMPH_GL(IL,JL) = 0.0D0
            TMPVAR_GL(IL,JL) = 0.0D0
         ENDDO
         ENDDO

         FNAMEU="uu-"//CHYEAR//"-"//CHMON//".nc"
!
!             Read U From LICOM Global Products For Nudging 
!
         if ( mytid==0 ) then

         iret = nf_open('./boundary-data-glorys/'//FNAMEU, nf_nowrite, ncid)
         call check_err(iret)

         iret = nf_open('ind_final_eas_uv.nc', nf_nowrite, ncid_rg)
         call check_err(iret)

         endif

         DO K =1, km

            if ( mytid==0 ) then

            startd(1)=1; countd(1)=imt_global_gl
            startd(2)=1; countd(2)=jmt_global_gl
            startd(3)=k; countd(3)=1
            startd(4)=1; countd(4)=1

            iret = nf_inq_varid(ncid, 'uu', varid)
            call check_err(iret)
            iret = nf_get_vara_real(ncid, varid, startd, countd, &
                                    BUFF_GL_GL)
            call check_err(iret)

            start2(1)=1; count2(1)=imt_global
            start2(2)=1; count2(2)=jmt_global
            start2(3)=k; count2(3)=1
            iret = nf_inq_varid(ncid_rg, 'ind', ind_id_rggl)
            call check_err(iret)
            iret = nf_get_vara_int( ncid_rg, ind_id_rggl, start2, &
                                    count2,ind_rggl_2d)
            call check_err(iret)

            DO JL = 1, jmt_global
               DO IL=1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0*BUFF_GL_GL(IL+ &
                       (IGLSTART-1),JL+(JGLSTART-1))*ind_rggl_2d(IL,JL)
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

            DO JL = 1, 2
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0* &
                         BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO JL = jmt_global-2, jmt_global
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0* &
                           BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = 1, 3
               DO JL = 3, jmt_global-3
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0* &
                           BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = imt_global-1, imt_global
               DO JL = 3, jmt_global-3
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0* &
                          BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

!            DO JL = 1, jmt
!               DO IL = 1, imt_global
!                  UUOBCSGL(IL,JL,K) = UV_IJMT_GL(IL,jmt_global-jmt+JL)
!
!                  UUOBCNGL(IL,JL,K) = UV_IJMT_GL(IL,JL)
!               ENDDO
!            ENDDO
!
!            DO IL = 1, imt
!               DO JL = 1, jmt_global
!                  UUOBCEGL(IL,JL,K) = UV_IJMT_GL(imt_global-imt+IL,JL)
!
!                  UUOBCWGL(IL,JL,K) = UV_IJMT_GL(IL,JL)
!               ENDDO
!            ENDDO
!
            DO JL = 1, jmt_global
            DO IL = 1, imt_global
               TMPH_GL(IL,JL) = TMPH_GL(IL,JL) + 1.0D0* &
                                ind_rggl_2d(IL,JL)*dzp(k)
               TMPVAR_GL(IL,JL) = TMPVAR_GL(IL,JL) + &
                       1.0D0*ind_rggl_2d(IL,JL)*dzp(k)*UV_IJMT_GL(IL,JL)
            ENDDO
            ENDDO


            endif

            CALL global_distribute( UV_IJMT_GL, UU_MON_LOC(1,1,K) )

         ENDDO

         if ( mytid==0 ) then
         iret = nf_close(ncid)
         call check_err(iret)
         iret = nf_close(ncid_rg)
         call check_err(iret)
         endif


         if ( mytid==0 ) then
            DO JL = 1, jmt_global
            DO IL = 1, imt_global
               TMPH = TMPH_GL(IL,JL)
               IF ( TMPH.GT.0.0D0 ) THEN
                    UB_IJMT_GL(IL,JL) = TMPVAR_GL(IL,JL)/TMPH
               ELSE
                    UB_IJMT_GL(IL,JL) = 0.0D0
               ENDIF
            ENDDO
            ENDDO

!            DO J = 1, jmt
!               DO I = 1, imt_global
!                  UBOBCSGL(I,J) = UB_IJMT_GL(I,jmt_global-jmt+J)
!                  UBOBCNGL(I,J) = UB_IJMT_GL(I,J)
!               ENDDO
!            ENDDO
!            DO I = 1, imt
!               DO J = 1, jmt_global
!                  UBOBCEGL(I,J) = UB_IJMT_GL(imt_global-imt+I,J)
!                  UBOBCWGL(I,J) = UB_IJMT_GL(I,J)
!               ENDDO
!            ENDDO

         endif

         CALL global_distribute( UB_IJMT_GL, UB_MON_LOC )

!         CALL global_to_local_3d( UUOBCEGL, UUOBCELOC, KM)
!         CALL global_to_local_3d( UUOBCWGL, UUOBCWLOC, KM)
!
!         CALL distribute_3d_sn( UUOBCSGL, UUOBCSLOC, KM)
!         CALL distribute_3d_sn( UUOBCNGL, UUOBCNLOC, KM)
!
!!$OMP PARALLEL DO PRIVATE (K,J)
!         DO K = 1, KM
!         DO J = 1, jmt
!            UBCLE(1,J,K) = UUOBCELOC(IMT-3,J,K)
!            UBCLE(2,J,K) = UUOBCELOC(IMT-2,J,K)
!            UBCLW(1,J,K) = UUOBCWLOC(4,J,K)
!            UBCLW(2,J,K) = UUOBCWLOC(5,J,K)
!         ENDDO
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (K,I)
!         DO K = 1, KM
!         DO I = 1, imt
!            UBCLS(I,1,K) = UUOBCSLOC(I,JMT-4,K)
!            UBCLS(I,2,K) = UUOBCSLOC(I,JMT-3,K)
!            UBCLN(I,1,K) = UUOBCNLOC(I,3,K)
!            UBCLN(I,2,K) = UUOBCNLOC(I,4,K)
!         ENDDO
!         ENDDO
!
!         CALL global_to_local_2d( UBOBCEGL, UBOBCELOC )
!         CALL global_to_local_2d( UBOBCWGL, UBOBCWLOC )
!
!         CALL distribute_2d_sn( UBOBCSGL, UBOBCSLOC )
!         CALL distribute_2d_sn( UBOBCNGL, UBOBCNLOC )
!
!
!!$OMP PARALLEL DO PRIVATE (J)
!         DO J = 1, jmt
!            UBTRE(1,J) = UBOBCELOC(IMT-3,J)
!            UBTRE(2,J) = UBOBCELOC(IMT-2,J)
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (J)
!         DO J = 1, jmt
!            UBTRW(1,J) = UBOBCWLOC(4,J)
!            UBTRW(2,J) = UBOBCWLOC(5,J)
!         ENDDO
!
!!$OMP PARALLEL DO PRIVATE (I)
!         DO I = 1, imt
!            UBTRS(I,1) = UBOBCSLOC(I,JMT-4)
!            UBTRS(I,2) = UBOBCSLOC(I,JMT-3)
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (I)
!         DO I = 1, imt
!            UBTRN(I,1) = UBOBCNLOC(I,3)
!            UBTRN(I,2) = UBOBCNLOC(I,4)
!         ENDDO
!        
         deallocate( UB_IJMT_GL)
!
!
!
         if ( allocated(BUFF_GL_GL) ) deallocate(BUFF_GL_GL)
         allocate( BUFF_GL_GL(imt_global_gl,jmt_global_gl) )

         if ( allocated( UV_IJMT_GL ) ) deallocate(UV_IJMT_GL)
         allocate( UV_IJMT_GL(imt_global,jmt_global) )
        
         if ( allocated( VB_IJMT_GL ) ) deallocate(VB_IJMT_GL)
         allocate( VB_IJMT_GL(imt_global,jmt_global) )
        
         DO JL = 1, jmt_global
         DO IL = 1, imt_global
            TMPH_GL(IL,JL) = 0.0D0
            TMPVAR_GL(IL,JL) = 0.0D0
         ENDDO
         ENDDO

         FNAMEU="vv-"//CHYEAR//"-"//CHMON//".nc"
!
!             Read V From LICOM Global Products For Nudging 
!
         if ( mytid==0 ) then

         iret = nf_open('./boundary-data-glorys/'//FNAMEU, nf_nowrite, ncid)
         call check_err(iret)

         iret = nf_open('ind_final_eas_uv.nc', nf_nowrite, ncid_rg)
         call check_err(iret)

         endif

         DO K =1, km

            if ( mytid==0 ) then

            startd(1)=1; countd(1)=imt_global_gl
            startd(2)=1; countd(2)=jmt_global_gl
            startd(3)=k; countd(3)=1
            startd(4)=1; countd(4)=1

            iret = nf_inq_varid(ncid, 'vv', varid)
            call check_err(iret)
            iret = nf_get_vara_real(ncid, varid, startd, countd, &
                                    BUFF_GL_GL)
            call check_err(iret)
            
            BUFF_GL_GL = BUFF_GL_GL * (1.0d0)

            start2(1)=1; count2(1)=imt_global
            start2(2)=1; count2(2)=jmt_global
            start2(3)=k; count2(3)=1
            iret = nf_inq_varid(ncid_rg, 'ind', ind_id_rggl)
            call check_err(iret)
            iret = nf_get_vara_int( ncid_rg, ind_id_rggl, start2, &
                                    count2,ind_rggl_2d)
            call check_err(iret)

            DO JL = 1, jmt_global
               DO IL=1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                           spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0*BUFF_GL_GL(IL+ &
                        (IGLSTART-1),JL+(JGLSTART-1))*ind_rggl_2d(IL,JL)
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO

            DO JL = 1, 2
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0* &
                         BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO JL = jmt_global-2, jmt_global
               DO IL = 1, imt_global
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0* &
                           BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = 1, 3
               DO JL = 3, jmt_global-3
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0* &
                           BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
            DO IL = imt_global-1, imt_global
               DO JL = 3, jmt_global-3
                  IF ( BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1)).NE. &
                                                         spvalue ) THEN
                       UV_IJMT_GL(IL,JL) = 1.0D0* &
                          BUFF_GL_GL(IL+(IGLSTART-1),JL+(JGLSTART-1))
                  ELSE
                       UV_IJMT_GL(IL,JL) = 0.0D0
                  ENDIF
               ENDDO
            ENDDO
!            DO JL = 1, jmt
!               DO IL = 1, imt_global
!                  VVOBCSGL(IL,JL,K) = UV_IJMT_GL(IL,jmt_global-jmt+JL)
!                  VVOBCNGL(IL,JL,K) = UV_IJMT_GL(IL,JL)
!               ENDDO
!            ENDDO
!
!            DO IL = 1, imt
!               DO JL = 1, jmt_global
!                  VVOBCEGL(IL,JL,K) = UV_IJMT_GL(imt_global-imt+IL,JL)
!                  VVOBCWGL(IL,JL,K) = UV_IJMT_GL(IL,JL)
!               ENDDO
!            ENDDO
!
            DO JL = 1, jmt_global
            DO IL = 1, imt_global
               TMPH_GL(IL,JL) = TMPH_GL(IL,JL) + 1.0D0* &
                                ind_rggl_2d(IL,JL)*dzp(k)
               TMPVAR_GL(IL,JL) = TMPVAR_GL(IL,JL) + 1.0D0* &
                             ind_rggl_2d(IL,JL)*dzp(k)*UV_IJMT_GL(IL,JL)
            ENDDO
            ENDDO

            endif

            CALL global_distribute( UV_IJMT_GL, VV_MON_LOC(1,1,K) )

         ENDDO

         if ( mytid==0 ) then
         iret = nf_close(ncid)
         call check_err(iret)
         iret = nf_close(ncid_rg)
         call check_err(iret)
         endif

         if ( mytid==0 ) then
              DO JL = 1, jmt_global
              DO IL = 1, imt_global

                 TMPH = TMPH_GL(IL,JL)
                 IF ( TMPH.GT.0.0D0 ) THEN
                      VB_IJMT_GL(IL,JL) = TMPVAR_GL(IL,JL)/TMPH
                 ELSE
                      VB_IJMT_GL(IL,JL) = 0.0D0
                 ENDIF
              ENDDO
              ENDDO
!
!              DO J = 1, jmt
!              DO I = 1, imt_global
!                 VBOBCSGL(I,J) = VB_IJMT_GL(I,jmt_global-jmt+J)
!                 VBOBCNGL(I,J) = VB_IJMT_GL(I,J)
!              ENDDO
!              ENDDO
!              DO I = 1, imt
!              DO J = 1, jmt_global
!                 VBOBCEGL(I,J) = VB_IJMT_GL(imt_global-imt+I,J)
!                 VBOBCWGL(I,J) = VB_IJMT_GL(I,J)
!              ENDDO
!              ENDDO

         endif

         CALL global_distribute( VB_IJMT_GL, VB_MON_LOC )


!         CALL global_to_local_3d( VVOBCEGL, VVOBCELOC, KM)
!         CALL global_to_local_3d( VVOBCWGL, VVOBCWLOC, KM)
!
!         CALL distribute_3d_sn( VVOBCSGL, VVOBCSLOC, KM)
!         CALL distribute_3d_sn( VVOBCNGL, VVOBCNLOC, KM)
!
!
!!$OMP PARALLEL DO PRIVATE (K,J)
!         DO K = 1, KM
!         DO J = 1, jmt
!            VBCLE(1,J,K) = VVOBCELOC(IMT-3,J,K)
!            VBCLE(2,J,K) = VVOBCELOC(IMT-2,J,K)
!            VBCLW(1,J,K) = VVOBCWLOC(4,J,K)
!            VBCLW(2,J,K) = VVOBCWLOC(5,J,K)
!         ENDDO
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (K,I)
!         DO K = 1, KM
!         DO I = 1, imt
!            VBCLS(I,1,K) = VVOBCSLOC(I,JMT-4,K)
!            VBCLS(I,2,K) = VVOBCSLOC(I,JMT-3,K)
!            VBCLN(I,1,K) = VVOBCNLOC(I,3,K)
!            VBCLN(I,2,K) = VVOBCNLOC(I,4,K)
!         ENDDO
!         ENDDO
!
!         CALL global_to_local_2d( VBOBCEGL, VBOBCELOC )
!         CALL global_to_local_2d( VBOBCWGL, VBOBCWLOC )
!
!         CALL distribute_2d_sn( VBOBCSGL, VBOBCSLOC )
!         CALL distribute_2d_sn( VBOBCNGL, VBOBCNLOC )
!
!
!!$OMP PARALLEL DO PRIVATE (J)
!         DO J = 1, jmt
!            VBTRE(1,J) = VBOBCELOC(IMT-3,J)
!            VBTRE(2,J) = VBOBCELOC(IMT-2,J)
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (J)
!         DO J = 1, jmt
!            VBTRW(1,J) = VBOBCWLOC(4,J)
!            VBTRW(2,J) = VBOBCWLOC(5,J)
!         ENDDO
!
!!$OMP PARALLEL DO PRIVATE (I)
!         DO I = 1, imt
!            VBTRS(I,1) = VBOBCSLOC(I,JMT-4)
!            VBTRS(I,2) = VBOBCSLOC(I,JMT-3)
!         ENDDO
!        
!!$OMP PARALLEL DO PRIVATE (I)
!         DO I = 1, imt
!            VBTRN(I,1) = VBOBCNLOC(I,3)
!            VBTRN(I,2) = VBOBCNLOC(I,4)
!         ENDDO
!
         deallocate( VB_IJMT_GL )
         deallocate( UV_IJMT_GL )
         deallocate( BUFF_GL_GL )

      RETURN

      END SUBROUTINE openbc_pre
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
      SUBROUTINE UBVB_VINTEG_IJMT_GL( find_gl, VARGL, VARBGL )

      use precision_mod

      use param_mod, only: imt_global, jmt_global, km, i, j, k

      use pconst_mod, only: dzp

      IMPLICIT NONE

      REAL(r8) :: find_gl(imt_global,jmt_global,km)
      REAL(r8) :: VARGL(imt_global,jmt_global,km)

      REAL(r8) :: VARBGL(imt_global,jmt_global)

      REAL(r8) :: TMPH, TMPHOB
      REAL(r8) :: TMPVAR
      
      DO i = 1, imt_global
         DO j = 1, jmt_global          
            VARBGL(i,j) = 0.0D0
         ENDDO
      ENDDO

      DO j = 1, jmt_global
         DO i = 1, imt_global
         
            TMPH = 0.0D0
            TMPHOB = 0.0D0
            TMPVAR = 0.0D0

            DO k = 1, km
               TMPH = TMPH + find_gl(i,j,k)*dzp(k)
            ENDDO

            IF ( TMPH .GT. 0.0D0 ) THEN
                 TMPHOB = 1.0D0/TMPH
            ELSE
                 TMPHOB = 0.0D0
            ENDIF

            DO k = 1, km
               TMPVAR = TMPVAR+VARGL(i,j,k)*dzp(k)*TMPHOB*find_gl(i,j,k)
            ENDDO

            VARBGL(i,j) = TMPVAR

         ENDDO
      ENDDO

      RETURN

      END SUBROUTINE UBVB_VINTEG_IJMT_GL

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
