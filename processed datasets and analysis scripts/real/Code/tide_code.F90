!
!      tidie codee 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE ASTRONOMIC_NOD( Tyear, Tdays )

#include <def-undef.h>
      use precision_mod

#if (defined OBTIDE || defined POTTIDE)
      use tide_mod, only: Ntide, FNOD, UNOD, V0NOD 
!#endif

      implicit none

      REAL :: Tyear, Tdays

      INTEGER, PARAMETER :: Ntide3=30
      INTEGER, PARAMETER :: Ntide4=40
      INTEGER, PARAMETER :: Ntide5=50

      REAL(r8), PARAMETER, DIMENSION(Ntide4) :: &
           FNODCO = (/1.0060_r8,  0.1150_r8, -0.0088_r8,  0.0006_r8, &  ! K1
                      1.0089_r8,  0.1871_r8, -0.0147_r8,  0.0014_r8, &  ! O1
                      1.0000_r8,  0.0000_r8,  0.0000_r8,  0.0000_r8, &  ! P1
                      1.0089_r8,  0.1871_r8, -0.0147_r8,  0.0014_r8, &  ! Q1
                      1.0004_r8, -0.0373_r8,  0.0002_r8,  0.0000_r8, &  ! M2
                      1.0000_r8,  0.0000_r8,  0.0000_r8,  0.0000_r8, &  ! S2
                      1.0004_r8, -0.0373_r8,  0.0002_r8,  0.0000_r8, &  ! N2
                      1.0241_r8,  0.2861_r8,  0.0083_r8, -0.0015_r8, &  ! K2
                      1.0429_r8,  0.4135_r8, -0.0040_r8,  0.0000_r8, &  ! Mf
                      1.0000_r8, -0.1300_r8,  0.0013_r8,  0.0000_r8/)  ! Mm

      REAL(r8), PARAMETER, DIMENSION(Ntide3) :: &
                UNODCO = (/ -8.86_r8,  0.68_r8, -0.07_r8, &
                            10.80_r8, -1.34_r8,  0.19_r8, & 
                             0.00_r8,  0.00_r8,  0.00_r8, &
                            10.80_r8, -1.34_r8,  0.19_r8, &
                            -2.14_r8,  0.00_r8,  0.00_r8, &
                             0.00_r8,  0.00_r8,  0.00_r8, &
                            -2.14_r8,  0.00_r8,  0.00_r8, &
                           -17.74_r8,  0.68_r8, -0.04_r8, &
                           -23.74_r8,  2.68_r8, -0.38_r8, &
                             0.00_r8,  0.00_r8,  0.00_r8/) 

      REAL(r8), PARAMETER, DIMENSION(Ntide5) :: &
                V0CO = (/0.00_r8,  1.00_r8,  0.00_r8,  0.00_r8,   90.0_r8, &
                        -2.00_r8,  1.00_r8,  0.00_r8,  0.00_r8,  270.0_r8, &
                         0.00_r8, -1.00_r8,  0.00_r8,  0.00_r8,  270.0_r8, &
                        -3.00_r8,  1.00_r8,  1.00_r8,  0.00_r8,  270.0_r8, &
                        -2.00_r8,  2.00_r8,  0.00_r8,  0.00_r8,    0.0_r8, &
                         0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8,    0.0_r8, &
                        -3.00_r8,  2.00_r8,  1.00_r8,  0.00_r8,    0.0_r8, &
                         0.00_r8,  2.00_r8,  0.00_r8,  0.00_r8,    0.0_r8, &
                         2.00_r8,  0.00_r8,  0.00_r8,  0.00_r8,    0.0_r8, &
                         1.00_r8,  0.00_r8, -1.00_r8,  0.00_r8,    0.0_r8/)
 
      REAL(r8) :: Astro_S
      REAL(r8) :: Astro_H
      REAL(r8) :: Astro_P
      REAL(r8) :: Astro_Ps
      REAL(r8) :: Astro_N

      INTEGER :: Lmon
      REAL(r8) :: Tts
      REAL(r8) :: DtoR

      INTEGER :: KT
      INTEGER :: KP

      Lmon = int( (Tyear-1901)/4 )
      Tts  = ( (Tyear-1900.0D0)*365.0D0 + 1.0D0*Lmon + 1.0D0*Tdays )/36525.0D0
      
      DtoR = 3.1415926D0/180.0D0

      Astro_S  = ( 277.0248 + 481267.8906*Tts + 0.0020*Tts*Tts ) * DtoR
      Astro_H  = ( 280.1895 + 36000.7689*Tts + 0.0003*Tts*Tts ) * DtoR
      Astro_P  = ( 334.3853 + 4069.0340*Tts - 0.0103*Tts*TtS ) * DtoR
      Astro_Ps = ( 281.2209 + 1.7192*Tts + 0.0005*Tts*Tts ) * DtoR
      Astro_N  = ( -100.8432 - 1934.1420*Tts + 0.0021*Tts*Tts ) * DtoR

      DO KT = 1, Ntide
         KP = 4*KT
         FNOD(KT) = FNODCO(KP-3) + FNODCO(KP-2)*DCOS(Astro_N) + &
                                 FNODCO(KP-1)*DCOS(2.0D0*Astro_N) + &
                                 FNODCO(KP)*DCOS(3.0D0*Astro_N)
      ENDDO 


      DO KT = 1, Ntide
         KP = 3*KT
         UNOD(KT) = UNODCO(KP-2)*DSIN(Astro_N) + UNODCO(KP-1)* &
                    DSIN(2.0D0*Astro_N) + UNODCO(KP)*DSIN(3.0D0*Astro_N)
      ENDDO 


      DO KT = 1, Ntide
         KP = 5*KT
         V0NOD(KT) = V0CO(KP-4)*Astro_S + V0CO(KP-3)*Astro_H + &
                     V0CO(KP-2)*Astro_P + &
                     V0CO(KP-1)*Astro_Ps + V0CO(KP)*DtoR
      ENDDO 

#endif

      RETURN 

      ENDSUBROUTINE ASTRONOMIC_NOD
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE POT_TIDE

#include <def-undef.h>
      use precision_mod

      use param_mod, only :  imt_global, jmt_global, imt, jmt, i, j,&
                              num_overlap, mytid, nx_proc, ny_proc

      use pconst_mod, only : vit, viv, lat, lon, ix, iy, DTB, PI, &
                             basin, viv !viv_global

#if ( defined POTTIDE )

      use tide_mod, only: Ntide, FNOD, UNOD, V0NOD, &
                          Tintg_tide, OMGATIDE, ptide_elevat, &
                          CD_COASTAL, CD_INDIAN, CD_ATLANTIC, &
                          CD_PACIFIC, CD_TIDE  
!#endif

      implicit none
!
!                               w(1/day) 1 beta h a (m) 
! 1 K1 Luni-solar declinational 0.7292117 0.736 0.141565 
! 2 O1 Principal lunar declinational 0.6759774 0.695 0.100661 
! 3 P1 Principal solar declinational 0.7252295 0.706 0.046848 
! 4 Q1 Larger lunar elliptic 0.6495854 0.695 0.019273 
! 5 M2 Principal lunar 1.405189 0.693 0.242334 
! 6 S2 Principal solar 1.454441 0.693 0.112743 
! 7 N2 Largerl lunar elliptic 1.378797 0.693 0.046397 
! 8 K2 Luni-solar declinational 1.458423 0.693 0.030684 
!
      INTEGER, PARAMETER :: TN = 8

      REAL(r8), PARAMETER,DIMENSION(TN) :: &
           amp = (/0.141565_r8, 0.100661_r8, 0.046848_r8, 0.019273_r8,&
                   0.242334_r8, 0.112743_r8,0.046397_r8,0.030684_r8/)

      REAL(r8), PARAMETER,DIMENSION(TN) :: &
           freq = (/0.7292117_r8,0.6759774_r8,0.7252295_r8,0.6495854_r8,&
                    1.405189_r8,1.454441_r8,1.378797_r8,1.458423_r8/)

      REAL(r8), PARAMETER,DIMENSION(TN) :: &
           beta = (/0.736_r8,0.695_r8,0.706_r8,0.695_r8,&
                    0.693_r8,0.693_r8,0.693_r8,0.693_r8/)

      REAL(r8) :: ptlat
      REAL(r8) :: Ptlon

      REAL(r8) :: pot_lon(imt)
      REAL(r8) :: pot_lat(jmt)

      INTEGER, DIMENSION(imt) :: iglobali
      INTEGER, DIMENSION(jmt) :: jglobalj

      INTEGER :: KP

!$OMP PARALLEL DO PRIVATE (i,j)
      do i= 1, imt
      do j= 1, jmt

            ptide_elevat(i,j) = 0.0D0

            CD_TIDE(i,j) = 0.0D0
      enddo
      enddo

!$OMP PARALLEL DO PRIVATE (i)
      do i= 1, imt
            pot_lon(i) = 0.0d0
      enddo

!$OMP PARALLEL DO PRIVATE (j)
      do j= 1, jmt
            pot_lat(j) = 0.0d0
      enddo

!$OMP PARALLEL DO PRIVATE (i)
      do i = 1, imt
         iglobali(i) = ix*(imt-num_overlap)+i
      enddo

!$OMP PARALLEL DO PRIVATE (i)
      do i = 1, imt
         if ( iglobali(i)<=imt_global ) then
               pot_lon(i) = lon( iglobali(i) )
         else
               pot_lon(i) = lon( imt_global )
         endif
      enddo

!$OMP PARALLEL DO PRIVATE (j)
      do j = 1, jmt
         if ( iy==0 ) then
              jglobalj(j) = j
         else
              jglobalj(j) = iy*(jmt-num_overlap)+j
         endif
      enddo

!$OMP PARALLEL DO PRIVATE (j)
      do j = 1, jmt
         if ( jglobalj(j)<=jmt_global ) then
               pot_lat(j) = lat( jglobalj(j) )
         else
               pot_lat(j) = lat( jmt_global )
         endif
      enddo

!$OMP PARALLEL DO PRIVATE (j,i)
      do j = 1, jmt
      do i = 1, imt

         if ( (iglobali(i)<=imt_global) ) then
            if ( jglobalj(j)<=jmt_global ) then
               if ( viv(i,j,1).gt.0.5 ) then

                  if ( (basin(iglobali(i),jglobalj(j)).eq.1) .or. &
                       (basin(iglobali(i),jglobalj(j)).eq.2) .or. &
                       (basin(iglobali(i),jglobalj(j)).eq.5) ) then

                       CD_TIDE(i,j)=CD_ATLANTIC

                  else if ( (basin(iglobali(i),jglobalj(j)).eq.3) .or. &
                            (basin(iglobali(i),jglobalj(j)).eq.6) ) then

                       CD_TIDE(i,j) = CD_INDIAN

                  else

                       CD_TIDE(i,j) = CD_PACIFIC

                  endif
               else
                  CD_TIDE(i,j) = 0.0
               endif
            else
               if ( viv(i,j,1).gt.0.5 ) then

                  if ( ( basin(iglobali(i),jmt_global).eq.1) .or. &
                       ( basin(iglobali(i),jmt_global).eq.2) .or. &
                       ( basin(iglobali(i),jmt_global).eq.5) ) then

                       CD_TIDE(i,j)=CD_ATLANTIC

                  else if ( ( basin(iglobali(i),jmt_global).eq.3) .or. &
                            ( basin(iglobali(i),jmt_global).eq.6) ) then

                       CD_TIDE(i,j) = CD_INDIAN

                  else

                       CD_TIDE(i,j) = CD_PACIFIC

                  endif
               else
                  CD_TIDE(i,j) = 0.0
               endif
            endif

         else
            if ( jglobalj(j)<=jmt_global ) then
               if ( viv(i,j,1).gt.0.5 ) then

                  if ( (basin(imt_global,jglobalj(j)).eq.1) .or. &
                       (basin(imt_global,jglobalj(j)).eq.2) .or. &
                       (basin(imt_global,jglobalj(j)).eq.5) ) then

                       CD_TIDE(i,j)=CD_ATLANTIC

                  else if ( (basin(imt_global,jglobalj(j)).eq.3) .or. &
                            (basin(imt_global,jglobalj(j)).eq.6) ) then

                       CD_TIDE(i,j) = CD_INDIAN

                  else

                       CD_TIDE(i,j) = CD_PACIFIC

                  endif
               else
                  CD_TIDE(i,j) = 0.0
               endif
            else
               if ( viv(i,j,1).gt.0.5 ) then

                  if ( (basin(imt_global,jmt_global).eq.1) .or. &
                       (basin(imt_global,jmt_global).eq.2) .or. &
                       (basin(imt_global,jmt_global).eq.5) ) then

                       CD_TIDE(i,j)=CD_ATLANTIC

                  else if ( (basin(imt_global,jmt_global).eq.3) .or. &
                            (basin(imt_global,jmt_global).eq.6) ) then

                       CD_TIDE(i,j) = CD_INDIAN

                  else

                       CD_TIDE(i,j) = CD_PACIFIC

                  endif
               else
                  CD_TIDE(i,j) = 0.0
               endif
            endif
         endif
      enddo
      enddo


      DO KP = 1, TN

!$OMP PARALLEL DO PRIVATE (I,J)
         DO I = 1, IMT

            ptlon = pot_lon(I)*2.*PI/360.0

            DO J= 1, JMT

               ptlat = pot_lat(J)*2.*PI/360.0

!            DO KP = 1, TN

               IF( KP > 4 )THEN
                  ptide_elevat(I,J) = ptide_elevat(I,J) + &
                                      vit(I,J,1)*beta(KP)*amp(KP)* &
                                      FNOD(KP)*COS(ptlat)*COS(ptlat)* &
                                      COS(freq(KP)*1e-4*Tintg_tide*DTB+&
                                          2.0*ptlon +UNOD(KP)+V0NOD(KP))
!                                      vit(I,J,1)*beta(KP)*amp(KP)*COS(ptlat)*COS(ptlat)  &
!                                      *COS( freq(KP)*1e-4*Tintg_tide*DTB + 2.0*ptlon )
               ELSE
                  ptide_elevat(I,J) = ptide_elevat(I,J) + &
                                      vit(I,J,1)*beta(KP)*amp(KP)* &
                                      FNOD(KP)*SIN(2.0*ptlat)* &
                                      COS(freq(KP)*1e-4*Tintg_tide*DTB+&
                                           ptlon+V0NOD(KP)+UNOD(KP))
!                                      vit(I,J,1)*beta(KP)*amp(KP)*SIN(2.0*ptlat) &
!                                      *COS( freq(KP)*1e-4*Tintg_tide*DTB + ptlon ) 
               END IF
            END DO

         END DO

      END DO

#endif

      RETURN

      END SUBROUTINE POT_TIDE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

      SUBROUTINE Read_TPXO7

#include <def-undef.h>

      use precision_mod

      use pconst_mod, only : imt_global, jmt_global, imt, jmt, i, j, k
      
      use param_mod, only : mytid

#ifdef OBCDT
      use openbc_mod, only : STARTLAT, STARTLON
#endif

#if (defined OBTIDE )
      use tide_mod, only: Ntide, OMGATIDE, H0TIDE_TPXO7, THTIDE_TPXO7, &
                          H0TIDE_TPXO7_LOC, THTIDE_TPXO7_LOC
!#endif
      implicit none

#include <netcdf.inc>

!
!      include'/THL7/home/liuhailong/netcdf-3.6.0-x64/include/netcdf.inc'
!
 
      integer, parameter :: imt_gl=1440
      integer, parameter :: jmt_gl=721

      integer, parameter :: kperi = Ntide
      integer, parameter :: imt_rg = imt_global
      integer, parameter :: jmt_rg = jmt_global

      real(r4), dimension(imt_gl,jmt_gl) :: ssh_r_gl
      real(r4), dimension(imt_gl,jmt_gl) :: ssh_i_gl

      real(r4), dimension(imt_gl,jmt_gl) :: ssh_amp_gl
      real(r4), dimension(imt_gl,jmt_gl) :: ssh_pha_gl
      real(r4), dimension(imt_gl,jmt_gl) :: h_topo_gl

      real(r4), dimension(imt_rg,jmt_rg) :: ssh_r_rg
      real(r4), dimension(imt_rg,jmt_rg) :: ssh_i_rg

      real(r4), dimension(imt_rg,jmt_rg) :: ssh_amp_rg
      real(r4), dimension(imt_rg,jmt_rg) :: ssh_pha_rg
      
      real(r4), dimension(imt_gl) :: lon_ssh_gl
      real(r4), dimension(jmt_gl) :: lat_ssh_gl
      real(r4), dimension(kperi) :: peri_ssh

      real(r4), dimension(imt_rg) :: lon_ssh_rg
      real(r4), dimension(jmt_rg) :: lat_ssh_rg
      
      integer :: lat_ssh_id
      integer :: lon_ssh_id

      integer :: lev_id
      integer :: lev_id1

!      integer :: lev_peri_id
!      integer :: lev_rg_id

      integer :: lat_ssh_rg_id
      integer :: lon_ssh_rg_id
!
      integer :: lat_ssh_dimid
      integer :: lon_ssh_dimid

      integer :: lev_dimid
!      integer :: lev_rg_dimid

      integer :: lat_ssh_rg_dimid
      integer :: lon_ssh_rg_dimid

      integer :: ssh_r_id
      integer :: ssh_i_id
!
      integer :: h_topo_id
!      integer :: h_topo_rg_id

      integer :: ssh_r_rg_id
      integer :: ssh_i_rg_id

      integer :: var_id
      integer :: var1_id
      integer :: var2_id
!        
      integer, parameter :: NDIMS=3


!      integer :: var_rg_id
!      integer :: var1_rg_id
!      integer :: var2_rg_id
!        
      integer :: dimids(NDIMS)

      integer :: start1(NDIMS), count1(NDIMS)

      integer :: start2(NDIMS), count2(NDIMS)

      integer :: start3(2), count3(2)

      integer :: dimids_rg(NDIMS)

      integer, parameter :: NLATS=jmt_gl
      integer, parameter :: NLONS=imt_gl

      integer, parameter :: NLEV=kperi
!      integer, parameter :: NLEV_RG=kperi
!
      integer, parameter :: NLATS_RG=jmt_rg
      integer, parameter :: NLONS_RG=imt_rg

      character*7, parameter :: LAT_NAME = "lat_ssh"

      character*7, parameter :: LON_NAME ="lon_ssh"

      character*7, parameter :: LEV_NAME="periods"
!
      character*7, parameter :: VAR1_NAME="ssh_amp"

      character*7, parameter :: VAR2_NAME="ssh_pha"

!
      integer :: iret
      integer :: in_ncid
 
      integer :: out_ncid

      integer :: out_rg_ncid

      character :: fnmout*12

      character :: fnmout_rg*16

!      character :: fnmin*8

      real(r4) :: stlat, stlon
!
      real(r4) :: rpi
      real(r4) :: radtod

      integer :: ixf, jyf
      integer :: ixt, jyt
      real(r4) :: dxf, dyf

      real(r4), dimension(imt_rg,jmt_rg) :: latt, lont

      real(r4) :: var_undef
      real(r8) :: varundef

      DO I = 1, imt_global
      DO J = 1, jmt_global
         H0TIDE_TPXO7(I,J) = 0.0D0
         THTIDE_TPXO7(I,J) = 0.0D0
      ENDDO
      ENDDO

      DO K = 1, kperi
!$OMP PARALLEL DO PRIVATE (J,I)
         DO I = 1, IMT
         DO J = 1, JMT
            H0TIDE_TPXO7_LOC(I,J,K) = 0.0D0
            THTIDE_TPXO7_LOC(I,J,K) = 0.0D0
         ENDDO
         ENDDO
      ENDDO

!
!      fnmin="TPXO7.nc"
!
      fnmout="TPXO7_SSH.nc"
!
      fnmout_rg="TPXO7_SSH_EAS.nc"

      if (mytid==0)then
          write(6,*)"BEGINNING-----Read_TPXO7 !"
      endif

      stlat = 1.0*STARTLAT-(jmt_global-1)*0.1
      stlon = 1.0*STARTLON
!
      dxf = 0.25
      dyf = 0.25

      rpi = 3.1415926
      radtod = 180.0/rpi

      var_undef = 1.e+35
      varundef = 1.e+35

      if ( mytid==0 ) then

           iret = nf_open('./tide-data/TPXO7.nc', nf_nowrite, in_ncid)
!      iret = nf_open(fnmin, nf_nowrite, in_ncid)
!      endif

      call check(nf_create(fnmout,OR(NF_CLOBBER,NF_64BIT_OFFSET), &
                      out_ncid))

      call check(nf_def_dim(out_ncid,LON_NAME,NLONS,lon_ssh_dimid))

      call check(nf_def_dim(out_ncid,LAT_NAME,NLATS,lat_ssh_dimid))

      call check(nf_def_dim(out_ncid,LEV_NAME,NLEV,lev_dimid))
!!
!!      call check(nf_def_dim(out_ncid,TIME_NAME,NTIME,time_ssh_dimid))
!!
      dimids = (/lon_ssh_dimid,lat_ssh_dimid,lev_dimid/)

      call check(nf_def_var(out_ncid,LON_NAME,NF_REAL,1, &
                   lon_ssh_dimid, lon_ssh_id))

      call check(nf_put_att_text(out_ncid,lon_ssh_id,'long_name',&
                                 17,'longitude for ssh'))
      call check(nf_put_att_text(out_ncid,lon_ssh_id,'units',&
                                 12,'degrees east'))

      call check(nf_def_var(out_ncid,LAT_NAME,NF_REAL,1, &
                 lat_ssh_dimid, lat_ssh_id))

      call check(nf_put_att_text(out_ncid,lat_ssh_id,'long_name',&
                                 16,'latitude for ssh'))
      call check(nf_put_att_text(out_ncid,lat_ssh_id,'units',&
                                 13,'degrees north'))

      call check(nf_def_var(out_ncid,LEV_NAME,NF_REAL,1,&
                  lev_dimid,lev_id))

      call check(nf_put_att_text(out_ncid,lev_id,'long_name',&
                                 15,'periods for ssh'))
      call check(nf_put_att_text(out_ncid,lev_id,'units',&
                                 5,'hours'))
      call check(nf_def_var(out_ncid,VAR1_NAME,NF_REAL, &
                   3, dimids,var1_id))
      call check(nf_put_att_text(out_ncid,var1_id,'long_name',&
                                 19,'The Ampitude of SSH'))
      call check(nf_put_att_text(out_ncid,var1_id,'units',&
                                 1,'m'))
      call check(nf_put_att_text(out_ncid,var1_id,'FillValue',&
                                 4,'NaNf'))
!!                                 7,'1.e+35f'))

      call check(nf_def_var(out_ncid,VAR2_NAME,NF_REAL, &
                   3, dimids,var2_id))
      call check(nf_put_att_text(out_ncid,var2_id,'long_name',&
                                 24,'The Delayed Angle of SSH'))
      call check(nf_put_att_text(out_ncid,var2_id,'units',&
                                 7,'degrees'))
      call check(nf_put_att_text(out_ncid,var2_id,'FillValue',&
                                 4,'NaNf'))
!!                                 7,'1.e+35f'))
      call check(nf_enddef(out_ncid))

      print *, 'pass enddef out_ncid'

      call check(nf_create(fnmout_rg,OR(NF_CLOBBER,NF_64BIT_OFFSET), &
                      out_rg_ncid))

      call check(nf_def_dim(out_rg_ncid,LON_NAME,NLONS_RG,&
                            lon_ssh_rg_dimid))

      call check(nf_def_dim(out_rg_ncid,LAT_NAME,NLATS_RG,&
                            lat_ssh_rg_dimid))

      call check(nf_def_dim(out_rg_ncid,LEV_NAME,NLEV,lev_dimid))

      dimids_rg = (/lon_ssh_rg_dimid,lat_ssh_rg_dimid,lev_dimid/)

      call check(nf_def_var(out_rg_ncid,LON_NAME,NF_REAL,1, &
                   lon_ssh_rg_dimid, lon_ssh_rg_id))

      call check(nf_put_att_text(out_rg_ncid,lon_ssh_rg_id,'long_name',&
                                 17,'longitude for ssh'))
      call check(nf_put_att_text(out_rg_ncid,lon_ssh_rg_id,'units',&
                                 12,'degrees east'))

      call check(nf_def_var(out_rg_ncid,LAT_NAME,NF_REAL,1, &
                 lat_ssh_rg_dimid, lat_ssh_rg_id))

      call check(nf_put_att_text(out_rg_ncid,lat_ssh_rg_id,'long_name',&
                                 16,'latitude for ssh'))
      call check(nf_put_att_text(out_rg_ncid,lat_ssh_rg_id,'units',&
                                 13,'degrees north'))

      call check(nf_def_var(out_rg_ncid,LEV_NAME,NF_REAL,1,&
                   lev_dimid,lev_id))

      call check(nf_put_att_text(out_rg_ncid,lev_id,'long_name',&
                                 15,'periods for ssh'))
      call check(nf_put_att_text(out_rg_ncid,lev_id,'units',&
                                 5,'hours'))

      call check(nf_def_var(out_rg_ncid,VAR1_NAME,NF_REAL, &
                   3, dimids_rg,var1_id))
      call check(nf_put_att_text(out_rg_ncid,var1_id,'long_name',&
                                 19,'The Ampitude of SSH'))
      call check(nf_put_att_text(out_rg_ncid,var1_id,'units',&
                                 1,'m'))
      call check(nf_put_att_text(out_rg_ncid,var1_id,'FillValue',&
!                                 4,'NaNf'))
                                 6,'1.e+35'))
      call check(nf_def_var(out_rg_ncid,VAR2_NAME,NF_REAL, &
                   3, dimids_rg,var2_id))
      call check(nf_put_att_text(out_rg_ncid,var2_id,'long_name',&
                                 24,'The Delayed Angle of SSH'))
      call check(nf_put_att_text(out_rg_ncid,var2_id,'units',&
                                 7,'degrees'))
      call check(nf_put_att_text(out_rg_ncid,var2_id,'FillValue',&
!                                 4,'NaNf'))
                                 6,'1.e+35'))

      call check(nf_enddef(out_rg_ncid))

      print *, 'pass enddef out_rg'

      iret = nf_inq_varid(in_ncid, 'lat_r', lat_ssh_id)
      call check(iret)
      iret = nf_get_vara_real(in_ncid, lat_ssh_id,1, jmt_gl, lat_ssh_gl)
      print *, 'pass lat_ssh_id'

      iret = nf_inq_varid(in_ncid, 'lon_r', lon_ssh_id)
      call check(iret)
      iret = nf_get_vara_real(in_ncid, lon_ssh_id,1, imt_gl,lon_ssh_gl)
      print *, 'pass lon_ssh_id'

      iret = nf_inq_varid(in_ncid, 'periods', lev_id1)
      call check(iret)
      iret = nf_get_vara_real(in_ncid, lev_id1,1, NLEV, peri_ssh)

       print *, 'passs lev_peri_id'

      do i = 1, imt_rg
         lon_ssh_rg(i) = stlon + (i-1)*0.1
      enddo
      do j = 1, jmt_rg
         lat_ssh_rg(j) = stlat + (j-1)*0.1
      enddo
!        
      ixf = imt_gl
      jyf = jmt_gl

      ixt = imt_rg
      jyt = jmt_rg

      do j = 1, jyt
      do i = 1, ixt
         lont(i,j) = lon_ssh_rg(i)
         latt(i,j) = lat_ssh_rg(j)
      enddo
      enddo
!
      call check(nf_put_vara_real(out_ncid,lon_ssh_id,1,NLONS,&
                 lon_ssh_gl))

      call check(nf_put_vara_real(out_ncid,lat_ssh_id,1,NLATS,&
                 lat_ssh_gl))

      call check(nf_put_vara_real(out_ncid,lev_id,1,NLEV,&
                 peri_ssh))

      print *, 'pass put gl'

      call check(nf_put_vara_real(out_rg_ncid,lon_ssh_rg_id,1,NLONS_RG,&
                 lon_ssh_rg))
!
      call check(nf_put_vara_real(out_rg_ncid,lat_ssh_rg_id,1,NLATS_RG,&
                 lat_ssh_rg))
!
      call check(nf_put_vara_real(out_rg_ncid,lev_id,1,NLEV,&
                 peri_ssh))
!
!      do k = 1, NLEV
!         OMGATIDE(k) = 2.0*rpi/(peri_ssh(k)*3600.0)
!      enddo
!
!
      print *, 'pass put rg'

      do k = 1, 4
         OMGATIDE(k) = 2.0*rpi/(peri_ssh(k+4)*3600.0)
      enddo
      do k = 1, 4
         OMGATIDE(k+4) = 2.0*rpi/(peri_ssh(k)*3600.0)
      enddo
      do k = 9, NLEV
         OMGATIDE(k) = 2.0*rpi/(peri_ssh(k)*3600.0)
      enddo
!
!      do k = 1, NLEV
!         PeriTIDE(k) = peri_ssh(k)
!      enddo
!
      start3(1)=1; count3(1)=imt_gl
      start3(2)=1; count3(2)=jmt_gl

      iret = nf_inq_varid(in_ncid, 'h', h_topo_id)
      call check(iret)
      iret = nf_get_vara_real(in_ncid, h_topo_id, start3, count3,&
                                h_topo_gl)
      endif

      do k = 1, NLEV
!
        if ( mytid==0 ) then

        start1(1)=1; count1(1)=imt_gl
        start1(2)=1; count1(2)=jmt_gl
        start1(3)=k; count1(3)=1

        iret = nf_inq_varid(in_ncid, 'ssh_r', var_id)
        call check(iret)
        iret = nf_get_vara_real(in_ncid, var_id, start1, count1,&
                                ssh_r_gl)

        iret = nf_inq_varid(in_ncid, 'ssh_i', var_id)
        call check(iret)
        iret = nf_get_vara_real(in_ncid, var_id, start1, count1,&
                                ssh_i_gl)

        do j = 1, jmt_gl
        do i = 1, imt_gl
!
            if ( isnan(ssh_r_gl(i,j)) ) ssh_r_gl(i,j) = var_undef
            if ( isnan(ssh_i_gl(i,j)) ) ssh_i_gl(i,j) = var_undef
        enddo
        enddo

        do i=1, imt_gl
           do j=1, jmt_gl

              if ( ssh_r_gl(i,j) .ne. var_undef .and. &
                   ssh_r_gl(i,j) .ne. var_undef ) then
 
                   ssh_amp_gl(i,j)=sqrt(ssh_r_gl(i,j)**2 + &
                                     ssh_i_gl(i,j)**2)
                   ssh_pha_gl(i,j)=atan2((-1.0*ssh_i_gl(i,j)),&
                                     ssh_r_gl(i,j))*radtod
                   if ( ssh_pha_gl(i,j).lt.0.0 ) then
                        ssh_pha_gl(i,j)=ssh_pha_gl(i,j)+360.0
                   endif
              else
                  ssh_amp_gl(i,j) = var_undef
                  ssh_pha_gl(i,j) = var_undef
              endif
           enddo
        enddo
!
        start2(1)=1; count2(1)=imt_gl
        start2(2)=1; count2(2)=jmt_gl
        start2(3)=k; count2(3)=1

        call check(nf_put_vara_real(out_ncid,var1_id,start2, &
                                                    count2, ssh_amp_gl))

        call check(nf_put_vara_real(out_ncid,var2_id,start2, &
                                                    count2, ssh_pha_gl))
!
!        CALL BILINE_INTERP( ixf, jyf, ixt, jyt, dxf, dyf, &
!                            lon_ssh_gl, lat_ssh_gl, ssh_r_gl, &
!                            latt, lont, ssh_r_rg )
! 
!        CALL BILINE_INTERP( ixf, jyf, ixt, jyt, dxf, dyf, &
!                            lon_ssh_gl, lat_ssh_gl, ssh_i_gl, &
!                            lat2d, lon2d, ssh_i_rg )
!!
        CALL BILINE_INTERP_DAI( ixf, jyf, ixt, jyt, dxf, dyf, &
                                lon_ssh_gl, lat_ssh_gl, ssh_r_gl, &
                                var_undef, latt, lont, ssh_r_rg )

        CALL BILINE_INTERP_DAI( ixf, jyf, ixt, jyt, dxf, dyf, &
                                lon_ssh_gl, lat_ssh_gl, ssh_i_gl,&
                                var_undef, latt, lont, ssh_i_rg )
        do i =1, imt_rg
           do j=1, jmt_rg

              if ( ssh_r_rg(i,j).ne.var_undef .and. &
                   ssh_i_rg(i,j).ne.var_undef ) then

                   if ( ssh_r_rg(i,j).eq.0.0 .and. &
                        ssh_i_rg(i,j).eq.0.0 ) then

                        ssh_amp_rg(i,j) = 0.0
                        ssh_pha_rg(i,j) = 0.0
                   else
                        ssh_amp_rg(i,j) = sqrt(ssh_r_rg(i,j)**2 + &
                                          ssh_i_rg(i,j)**2)

                        ssh_pha_rg(i,j) = (-1.0)*atan2(ssh_i_rg(i,j),&
                                           ssh_r_rg(i,j))*radtod

                        if ( ssh_pha_rg(i,j).lt.0.0 ) then
                           ssh_pha_rg(i,j)=ssh_pha_rg(i,j)+360.0
                        endif

                        ssh_pha_rg(i,j) = ssh_pha_rg(i,j)/radtod
                   endif

              else
                 ssh_amp_rg(i,j) = var_undef
                 ssh_pha_rg(i,j) = var_undef
              endif

           enddo
        enddo

        start2(1)=1; count2(1)=imt_rg
        start2(2)=1; count2(2)=jmt_rg
        start2(3)=k; count2(3)=1

        call check(nf_put_vara_real(out_rg_ncid,var1_id,start2, &
                                                    count2, ssh_amp_rg))

        call check(nf_put_vara_real(out_rg_ncid,var2_id,start2, &
                                                    count2, ssh_pha_rg))

        do i = 1, imt_rg
           do j = 1, jmt_rg

              if ( ssh_amp_rg(i,j) .ne. var_undef .and. &
                   ssh_pha_rg(i,j) .ne. var_undef ) then
 
                   H0TIDE_TPXO7(i,j) = 1.0D0*ssh_amp_rg(i,jmt_rg-j+1)
                   THTIDE_TPXO7(i,j) = 1.0D0*ssh_pha_rg(i,jmt_rg-j+1) 
              else
                   H0TIDE_TPXO7(i,j) = varundef
                   THTIDE_TPXO7(i,j) = varundef
              endif

           enddo
        enddo

        endif

        CALL global_distribute ( H0TIDE_TPXO7, H0TIDE_TPXO7_LOC(:,:,K) )

        CALL global_distribute ( THTIDE_TPXO7, THTIDE_TPXO7_LOC(:,:,K)) 

!!
!!            if ( ssh_amp_rg(i,3) .ne. var_undef .and. &
!!                 ssh_pha_rg(i,3) .ne. var_undef ) then
!!
!!                 H0TIDESGL(i,k) = ssh_amp_rg(i,3)
!!                 THTIDESGL(i,k) = ssh_pha_rg(i,3)
!!
!            else
!                 H0TIDESGL(i,k) = varundef
!                 THTIDESGL(i,k) = varundef
!            endif

!
!         do i = 1, imt_rg
!!
!            if ( ssh_amp_rg(i,2) .ne. var_undef .and. &
!                 ssh_pha_rg(i,2) .ne. var_undef ) then
!
!                 H0TIDESGL(i,k) = ssh_amp_rg(i,2)
!                 THTIDESGL(i,k) = ssh_pha_rg(i,2) 
!!
!!            if ( ssh_amp_rg(i,3) .ne. var_undef .and. &
!!                 ssh_pha_rg(i,3) .ne. var_undef ) then
!!
!!                 H0TIDESGL(i,k) = ssh_amp_rg(i,3)
!!                 THTIDESGL(i,k) = ssh_pha_rg(i,3)
!!
!            else
!                 H0TIDESGL(i,k) = varundef
!                 THTIDESGL(i,k) = varundef
!            endif
!
!            if ( ssh_amp_rg(i,jmt_rg-1) .ne. var_undef .and. &
!                 ssh_pha_rg(i,jmt_rg-1) .ne. var_undef ) then
!
!                 H0TIDENGL(i,k) = ssh_amp_rg(i,jmt_rg-1)            
!                 THTIDENGL(i,k) = ssh_pha_rg(i,jmt_rg-1)
!!
!!            if ( ssh_amp_rg(i,jmt_rg-2) .ne. var_undef .and. &
!!                 ssh_pha_rg(i,jmt_rg-2) .ne. var_undef ) then
!!
!!                 H0TIDENGL(i,k) = ssh_amp_rg(i,jmt_rg-2)            
!!                 THTIDENGL(i,k) = ssh_pha_rg(i,jmt_rg-2)
!!
!            else
!         
!                 H0TIDENGL(i,k) = varundef            
!                 THTIDENGL(i,k) = varundef
!            endif
!
!         enddo
!
!         do j = 1, jmt_rg
!
!            if ( ssh_amp_rg(imt_rg-1,j) .ne. var_undef .and. &
!                 ssh_pha_rg(imt_rg-1,j) .ne. var_undef ) then
!
!                H0TIDEEGL(jmt_rg-j+1,k) = ssh_amp_rg(imt_rg-1,j)
!                THTIDEEGL(jmt_rg-j+1,k) = ssh_pha_rg(imt_rg-1,j)
!
!!            if ( ssh_amp_rg(imt_rg-2,j) .ne. var_undef .and. &
!!                 ssh_pha_rg(imt_rg-2,j) .ne. var_undef ) then
!!
!!                H0TIDEEGL(jmt_rg-j+1,k) = ssh_amp_rg(imt_rg-2,j)
!!                THTIDEEGL(jmt_rg-j+1,k) = ssh_pha_rg(imt_rg-2,j)
!!
!!
!            else
!
!                H0TIDEEGL(jmt_rg-j+1,k) = varundef
!                THTIDEEGL(jmt_rg-j+1,k) = varundef
!
!            endif
!
!            if ( ssh_amp_rg(2,j) .ne. var_undef .and. &
!                 ssh_pha_rg(2,j) .ne. var_undef ) then
!
!                 H0TIDEWGL(jmt_rg-j+1,k) = ssh_amp_rg(2,j)
!                 THTIDEWGL(jmt_rg-j+1,k) = ssh_pha_rg(2,j)
!
!!            if ( ssh_amp_rg(3,j) .ne. var_undef .and. &
!!                 ssh_pha_rg(3,j) .ne. var_undef ) then
!!
!!                 H0TIDEWGL(jmt_rg-j+1,k) = ssh_amp_rg(3,j)
!!                 THTIDEWGL(jmt_rg-j+1,k) = ssh_pha_rg(3,j)
!!
!            else
!
!                 H0TIDEWGL(jmt_rg-j+1,k) = varundef
!                 THTIDEWGL(jmt_rg-j+1,k) = varundef
!
!            endif           
!
!         enddo
!

      enddo
!
      if ( mytid==0 ) then
             
      call check(nf_close(in_ncid))
! 
      call check(nf_close(out_ncid))
!!
      call check(nf_close(out_rg_ncid))
!
      endif

      if (mytid==0)then
          write(6,*)"END-----------Read_TPXO7 !"
      endif 

#endif

      return

      end SUBROUTINE Read_TPXO7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      SUBROUTINE check (iret)

#include <netcdf.inc>
!
!      include'/THL7/home/liuhailong/netcdf-3.6.0-x64/include/netcdf.inc'
!
      INTEGER :: iret
      IF (iret /= NF_NOERR) THEN
         PRINT *, nf_strerror (iret)
         STOP
      END IF
      END SUBROUTINE check
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!      subroutine BILINE_INTERP( ixf,jyf,ixt,jyt,dxf,dyf,&
!                                lonf,latf,varf,&
!                                latt,lont,vart )
!
!#include <netcdf.inc>
!
!      use precision_mod
!
!      implicit none
!
!      integer :: ixf,jyf
!      integer :: ixt,jyt
!      real(r4) :: dxf, dyf
!
!      real(r4),dimension(ixf,jyf):: varf
!      real(r4),dimension(ixf):: lonf
!      real(r4),dimension(jyf):: latf
!
!      real(r4),dimension(ixt,jyt):: latt,lont,vart
!
!      integer,dimension(ixt,jyt):: isw,ise,inw,ine
!      integer,dimension(ixt,jyt):: jsw,jse,jnw,jne
!      real(r4),dimension(ixt,jyt):: ww,we
!      real(r4),dimension(ixt,jyt):: ws,wn
!
!      integer i,j,k,ie,iw,jn,js
!!  
!!find the weight & indexes  
!!  
!      do j = 1, jyt
!      do i = 1, ixt
!!  
!!x direction  
!!  
!         iw = int( (lont(i,j)-lonf(1))/dxf )+1
!         ie = iw+1
!         if ( ie.gt.ixf ) then
!             ww(i,j) = (360.-lont(i,j))/dxf
!             we(i,j)=1-ww(i,j)
!             ie=1
!         else
!             ww(i,j)=(lonf(ie)-lont(i,j))/dxf
!             we(i,j)=1-ww(i,j)
!         endif
!
!         ise(i,j)=ie
!         ine(i,j)=ie
!         isw(i,j)=iw
!         inw(i,j)=iw
!      enddo
!      enddo
!!  
!!y direction  
!!  
!      do j=1,jyt
!      do i=1,ixt
!
!         js = int( (latt(i,j)-latf(1))/dyf )+1
!
!!         if(js.eq.jyf) js=js-1
!         if(js.ge.jyf) js=js-1
!
!         jn = js+1
!         ws(i,j) = (latf(jn)-latt(i,j))/(latf(jn)-latf(js))
!         wn(i,j) = 1.0-ws(i,j)
!
!         jnw(i,j)=jn
!         jne(i,j)=jn
!         jsw(i,j)=js
!         jse(i,j)=js
!      enddo
!      enddo
!!  
!!interpolate bilinearly  
!!  
!!  
!      do j = 1, jyt
!      do i = 1, ixt
!
!         vart(i,j) = ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j)) + &
!                     ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j)) + &
!                     we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j)) + &
!                     we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))
!      enddo
!      enddo
!
!      return
!
!      end subroutine BILINE_INTERP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine BILINE_INTERP_DAI( ixf,jyf,ixt,jyt,dxf,dyf,&
                                    lonf,latf,varf,var_undef, &
                                    latt,lont,vart )

      use precision_mod

      implicit none

#include <netcdf.inc>
!
!      include'/THL7/home/liuhailong/netcdf-3.6.0-x64/include/netcdf.inc'
!
      integer :: ixf,jyf
      integer :: ixt,jyt
      real(r4) :: dxf, dyf

      real(r4),dimension(ixf,jyf):: varf
      real(r4),dimension(ixf):: lonf
      real(r4),dimension(jyf):: latf

      real(r4),dimension(ixt,jyt):: latt,lont,vart

      real(r4) :: var_undef

      integer,dimension(ixt,jyt):: isw,ise,inw,ine
      integer,dimension(ixt,jyt):: jsw,jse,jnw,jne
      real(r4),dimension(ixt,jyt):: ww,we
      real(r4),dimension(ixt,jyt):: ws,wn
      integer :: i,j,k,ie,iw,jn,js

      integer :: ncid_rg_ind, ind_id
      integer, dimension(ixt,jyt) :: ind_rg
      integer :: start1(3), count1(3)

      integer :: iret

!      real(r4) :: var_undef
!
!      var_undef = 1.e+35
!
!      iret=nf_open('ind_final_eas_2.nc',nf_nowrite,ncid_rg_ind)
!
      iret=nf_open('ind_final_eas.nc',nf_nowrite,ncid_rg_ind)
!
      call check(iret)
      iret=nf_inq_varid(ncid_rg_ind,'ind',ind_id)
      call check(iret)
      start1(1)=1; count1(1)=ixt
      start1(2)=1; count1(2)=jyt
      start1(3)=1; count1(3)=1
      iret=nf_get_vara_int(ncid_rg_ind,ind_id,start1,count1,ind_rg)
      call check(nf_close(ncid_rg_ind))
!  
!find the weight & indexes  
!  
      do j = 1, jyt
      do i = 1, ixt
!  
!x direction  
!  
         iw = int( (lont(i,j)-lonf(1))/dxf )+1
         ie = iw+1

         if ( ie.gt.ixf ) then
             ww(i,j) = (360.-lont(i,j))/dxf
             we(i,j)=1-ww(i,j)
             ie=1
         else
             ww(i,j)=(lonf(ie)-lont(i,j))/dxf
             we(i,j)=1-ww(i,j)
         endif

         ise(i,j)=ie
         ine(i,j)=ie
         isw(i,j)=iw
         inw(i,j)=iw
      enddo
      enddo
!  
!y direction  
!  
      do j=1,jyt
      do i=1,ixt

         js = int( (latt(i,j)-latf(1))/dyf )+1

!         if(js.eq.jyf) js=js-1
         if(js.ge.jyf) js=js-1

         jn = js+1
         ws(i,j) = (latf(jn)-latt(i,j))/dyf
         wn(i,j) = 1.0-ws(i,j)

         jnw(i,j)=jn
         jne(i,j)=jn
         jsw(i,j)=js
         jse(i,j)=js
      enddo
      enddo
!  
!interpolate bilinearly  
!  
!  
      do j = 1, jyt
      do i = 1, ixt

         if ( varf(inw(i,j),jnw(i,j)).ne.var_undef .and. &
              varf(isw(i,j),jsw(i,j)).ne.var_undef .and. &
              varf(ise(i,j),jse(i,j)).ne.var_undef .and. &
              varf(ine(i,j),jne(i,j)).ne.var_undef         ) then

              vart(i,j) = ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j)) + &
                          ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j)) + &
                          we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j)) + &
                          we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))             

!              if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).ne.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).ne.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).ne.var_undef         ) then

                   vart(i,j)=ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j))+ &
                             we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j))+ &
                             we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).ne.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).ne.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).ne.var_undef         ) then

                   vart(i,j)=ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j))+ &
                             we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j))+ &
                             we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).ne.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).ne.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).ne.var_undef         ) then

                   vart(i,j)=ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j))+ &
                             ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j))+ &
                             we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).ne.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).ne.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).ne.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

                   vart(i,j)=ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j))+&
                          ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j)) + &
                          we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).ne.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).ne.var_undef         ) then

                   vart(i,j)=we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j))+ &
                             we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))

!                   vart(i,j)=ws(i,j)*varf(ise(i,j),jse(i,j))+ &
!                             wn(i,j)*varf(ine(i,j),jne(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).ne.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).ne.var_undef         ) then

                   vart(i,j)=ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j))+ &
                             we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).ne.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).ne.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

                   vart(i,j)=ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j))+ &
                             we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j))

!                   vart(i,j)=ww(i,j)*varf(isw(i,j),jsw(i,j))+ &
!                             we(i,j)*varf(ise(i,j),jse(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).ne.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).ne.var_undef         ) then

                   vart(i,j)=ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j))+ &
                             we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))

!                   vart(i,j)=ww(i,j)*varf(inw(i,j),jnw(i,j))+ &
!                             we(i,j)*varf(ine(i,j),jne(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef
!
         else if ( varf(inw(i,j),jnw(i,j)).ne.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).ne.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

                   vart(i,j)=ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j))+ &
                             we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).ne.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).ne.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

                   vart(i,j)=ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j))+ &
                             ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).ne.var_undef         ) then

                   vart(i,j)= we(i,j)*wn(i,j)*varf(ine(i,j),jne(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).ne.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

                   vart(i,j)=we(i,j)*ws(i,j)*varf(ise(i,j),jse(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).ne.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

                   vart(i,j)=ww(i,j)*ws(i,j)*varf(isw(i,j),jsw(i,j))

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).ne.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

                   vart(i,j) =ww(i,j)*wn(i,j)*varf(inw(i,j),jnw(i,j)) 

!                   if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

         else if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
                   varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
                   varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
                   varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

                   vart(i,j) = var_undef

         endif

         if ( ind_rg(i,jyt-j+1).eq.0 ) vart(i,j)=var_undef

      enddo
      enddo

      do j = 2, jyt-1
      do i = 2, ixt-1

         if ( varf(inw(i,j),jnw(i,j)).eq.var_undef .and. &
              varf(isw(i,j),jsw(i,j)).eq.var_undef .and. &
              varf(ise(i,j),jse(i,j)).eq.var_undef .and. &
              varf(ine(i,j),jne(i,j)).eq.var_undef         ) then

              if ( ind_rg(i,jyt-j+1).ne.0 ) then

                   if ( vart(i-1,j) .ne. var_undef ) then
                        vart(i,j) = vart(i-1,j)
                   else if ( vart(i,j-1) .ne. var_undef ) then
                        vart(i,j) = vart(i,j-1)
                   else if ( vart(i+1,j) .ne. var_undef ) then
                        vart(i,j) = vart(i+1,j)
                   else if ( vart(i,j+1) .ne. var_undef ) then
                        vart(i,j) = vart(i,j+1)
                   else
                        vart(i,j) = 0.0
                   endif
              endif

         endif

      enddo
      enddo

      return

      end subroutine BILINE_INTERP_DAI
!
