!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine sSAVEcdf 
#include <def-undef.h>
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use output_mod
use cdf_mod
       
#include <netcdf.inc> 
!
      character (len=18) :: fname
      character (len=22) :: pname 
      character (len=4 ) :: stryr
      integer :: klevel
      integer :: vartime
      real(r4), dimension(lon_len,lat_len,1,1) :: varvar2d
      real(r4), dimension(lon_len,lat_len,lev_len,1) :: varvar3d 
      real(r8), dimension(lon_len,lat_len,lev_len  ) :: vit_output_global, viv_output_global
      real(r8), dimension(    imt,    jmt,     km  ) :: vit_output       , viv_output
      
      pname='../output-era-glorys2/'
      write(stryr, '(i4)') the_first_year
      

!
         number_day = iday +1
         number_month = month
         if ( number_day  > imd ) then
             number_day = 1
             number_month = month + 1
         end if
         nwmf= iyfm -1 + the_first_year
!
     if ( mod(iday,rest_freq) == 0 .or. iday == imd) then
       if ( mytid == 0 ) then
         fname(1:8)='fort.22.'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon0
         write(fname(9:12),'(i4.4)')nwmf
         write(fname(17:18),'(i2.2)') number_day
!
         if ( number_day ==1 .and. mon0 < 12) then
             write(fname(14:15),'(i2.2)')mon0+1
         end if
!
         if ( number_day ==1 .and. mon0 == 12) then
             write(fname(14:15),'(i2.2)')mon0-11
             write(fname(9:12),'(i4.4)')nwmf+1
         end if

         open (17, file="rpointer.ocn", form='formatted')
         write(17,'(a18)') fname
         close(17)
if (number_day > -1 ) then
         open(22,file=trim(pname)//fname,form='unformatted')
endif 
       end if
#ifdef SPMD
!
if (number_day > -1 ) then
         allocate(buffer(imt_global,jmt_global))

         call local_to_global_4d_double(h0,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!
         do klevel=1,km
         call local_to_global_4d_double(u(1,1,klevel),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!
         do klevel=1,km
         call local_to_global_4d_double(v(1,1,klevel),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!
         do klevel=1,km
         call local_to_global_4d_double(at(1,1,klevel,1),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!
         do klevel=1,km
         call local_to_global_4d_double(at(1,1,klevel,2),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!lhl20110728
         do klevel=1,km
         call local_to_global_4d_double(ws(1,1,klevel),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do

         call local_to_global_4d_double(su,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         call local_to_global_4d_double(sv,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         call local_to_global_4d_double(swv,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         call local_to_global_4d_double(lwv,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         call local_to_global_4d_double(sshf,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         call local_to_global_4d_double(lthf,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         call local_to_global_4d_double(fresh,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!lhl20110728
         if (mytid==0) then
             write(22) number_month, number_day
             close(22)
         end if
!
         deallocate( buffer)
endif 

#ifdef DAILYMEAN_OP
!------------------------output dailymean-----------------------------------------------------------------------
!add by zhangjk 20250414

!prepare ----------------------------
         vartime = (iyfm -1 + the_first_year) * 10000 + month * 100 + iday 
         vit_output = real(vit,kind=8)
         call local_to_global_4d_double(vit_output,vit_output_global,km,1)
         viv_output = real(viv,kind=8)
         call local_to_global_4d_double(viv_output,viv_output_global,km,1)
!output z0 ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(1:7)='daym-z0'
         fname(16:18)='.nc'
         write(fname(12:13),'(i2.2)')mon0 
         write(fname(8:11),'(i4.4)')iyfm -1 + the_first_year
         write(fname(14:15),'(i2.2)')iday 
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', 1, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(6:7), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 18, 'sea surface height')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, 0)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

         call local_to_global_4d(z0mon,varvar2d,1,1)
         varvar2d = varvar2d / nss 
         if (mytid==0) then
         where (vit_output_global(:,:,1) < 0.5D0) varvar2d(:,:,1,1) = spval 
         start4(1)=1 ; count4(1)=imt_global
         start4(2)=1 ; count4(2)=jmt_global
         start4(3)=1 ; count4(3)=1
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar2d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif


!output uu ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(6:7)='uu'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(6:7), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

         call local_to_global_4d(usmon,varvar3d,km,1)
         varvar3d = varvar3d / nss 
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif

!output vv ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(6:7)='vv'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(6:7), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

         call local_to_global_4d(vsmon,varvar3d,km,1)
      !redefine the direction of v
         varvar3d = varvar3d / nss  * (-1)
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif

! !output ww ----------------------------------------------------------------------
! if (mytid==0) then
!    !get file name ----------------------
!          fname(6:7)='ww'
!    !create file ------------------------
!          iert = nf_create(pname//fname, NF_CLOBBER, ncid)
!          call check_err(iert)
!    !define dimension -------------------
!          iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
!          call check_err(iert)
!          iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
!          call check_err(iert)
!          iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
!          call check_err(iert)
!          iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
!          call check_err(iert)         
!    !define variables -------------------
!          iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
!          CALL check_err (iret)
!          iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
!          CALL check_err (iret)
!          iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
!          CALL check_err (iret)
!          iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
!          CALL check_err (iret)
!          iert = nf_def_var(ncid, fname(6:7), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
!          call check_err(iert)      
!    !assign attributes -------------------
!          iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
!          CALL check_err (iret)
!          iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
!          CALL check_err (iret)
!          iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
!          CALL check_err (iret)
!          iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
!          CALL check_err (iret)

!          iert = nf_enddef(ncid)
!          call check_err(iert)
!    !value to variables ------------------
!          iret = nf_put_var_real (ncid, lon_id, lon)
!          CALL check_err (iret)
!          iret = nf_put_var_real (ncid, lat_id, lat)
!          CALL check_err (iret)
!          iret = nf_put_var_real (ncid, lev_id, lev)
!          CALL check_err (iret)
!          iret = nf_put_var_int (ncid, itime_id, vartime)
!          CALL check_err (iret)
! endif

!          call local_to_global_4d(wsmon,varvar3d,km,1)
!          varvar3d = varvar3d / nss 
!          if (mytid==0) then
!          where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
!          start4(1)=1 ; count4(1)=lon_len
!          start4(2)=1 ; count4(2)=lat_len
!          start4(3)=1 ; count4(3)=lev_len 
!          start4(4)=1 ; count4(4)=1
!          iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
!          call check_err(iert)
!    ! close file ------------------------
!          iert = nf_close(ncid)
!          call check_err(iert)
!          endif

!output tt ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(6:7)='tt'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(6:7), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

         call local_to_global_4d(tsmon,varvar3d,km,1)
         varvar3d = varvar3d / nss 
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif

!output ss ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(6:7)='ss'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(6:7), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

         call local_to_global_4d(ssmon,varvar3d,km,1)
      !alternate unit of ss
         varvar3d = varvar3d / nss *1000 + 35 
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif
!---------------------------end output dailymean----------------------------------------------------------------
#endif 

#else
        write(22) h0
!
        do k =1, km
           write(22) ((u(i,j,k),i=1,imt_global),j=1,jmt_global)
        end do
!
        do k =1, km
           write(22) ((v(i,j,k),i=1,imt_global),j=1,jmt_global)
        end do
!
        do k =1, km
           write(22) ((at(i,j,k,1),i=1,imt_global),j=1,jmt_global)
        end do
!
        do k =1, km
           write(22) ((at(i,j,k,2),i=1,imt_global),j=1,jmt_global)
        end do
!lhl20110728
        do k =1, km
           write(22) ((ws(i,j,k),i=1,imt_global),j=1,jmt_global)
        end do
!lhl20110728
!
        write(22) number_month, number_day
        close(22)
!
#endif
      end if
!

!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
            DO i = 1,imt
               up (i,j,k) = u (i,j,k)
               vp (i,j,k) = v (i,j,k)
               utf (i,j,k) = u (i,j,k)
               vtf (i,j,k) = v (i,j,k)
               atb (i,j,k,1) = at (i,j,k,1)
               atb (i,j,k,2) = at (i,j,k,2)
            END DO
         END DO
      END DO

!$OMP PARALLEL DO PRIVATE (j,i)
!      DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
!         DO i = 1,imt
!            h0p (i,j)= h0 (i,j)
!            ubp (i,j)= ub (i,j)
!            vbp (i,j)= vb (i,j)
!            h0f (i,j)= h0 (i,j)
!            h0bf (i,j)= h0 (i,j)
!         END DO
!      END DO
!      ISB = 0
!      ISC = 0
!      IST = 0
!
!

      return
      end subroutine sSAVEcdf 


subroutine sSAVE_hourly(iii)

#include <def-undef.h>

use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use output_mod
use cdf_mod
       
#include <netcdf.inc> 
!
      character (len=18) :: fname
      character (len=17) :: pname 
      character (len=4 ) :: stryr
      integer :: klevel
      integer :: vartime
      integer :: iii 
      real(r4), dimension(imt    ,jmt    ,1      ,1) :: varvar2d4 
      real(r4), dimension(imt    ,jmt    ,km     ,1) :: varvar3d4 
      real(r4), dimension(lon_len,lat_len,1      ,1) :: varvar2d
      real(r4), dimension(lon_len,lat_len,lev_len,1) :: varvar3d 
      real(r8), dimension(lon_len,lat_len,lev_len  ) :: vit_output_global, viv_output_global
      real(r8), dimension(    imt,    jmt,     km  ) :: vit_output       , viv_output
      
      pname='../output-pttide/'
      write(stryr, '(i4)') the_first_year
      

!
         number_day = iday +1
         number_month = month
         if ( number_day  > imd ) then
             number_day = 1
             number_month = month + 1
         end if
         nwmf= iyfm -1 + the_first_year
!


#ifdef SPMD

if (mod(int(iii*dts) ,3600) == 0) then 

!------------------------output hourly-----------------------------------------------------------------------
!add by zhangjk 20250610

!prepare ----------------------------
         vartime = (iyfm -1 + the_first_year) * 10000 + month * 100 + iday 
         vit_output = real(vit,kind=8)
         call local_to_global_4d_double(vit_output,vit_output_global,km,1)
         viv_output = real(viv,kind=8)
         call local_to_global_4d_double(viv_output,viv_output_global,km,1)
!output z0 ----------------------------------------------------------------------
if (mytid==0) then
print*, iii 
   !get file name ----------------------
         fname(1:3)='z0-'
         fname(16:18)='.nc'
         write(fname(8:9),'(i2.2)')mon0 
         write(fname(4:7),'(i4.4)')iyfm -1 + the_first_year
         write(fname(10:11),'(i2.2)')iday
         write(fname(12:15),'(i4.4)') int(iii*dts/3600*100)
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', 1, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(1:2), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 18, 'sea surface height')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, 0)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

do i = 1,imt 
do j = 1,jmt 
varvar2d4(i,j,1,1) = h0 (i,j) 
enddo 
enddo 
         call local_to_global_4d(varvar2d4,varvar2d,1,1)
         if (mytid==0) then
         where (vit_output_global(:,:,1) < 0.5D0) varvar2d(:,:,1,1) = spval 
         start4(1)=1 ; count4(1)=imt_global
         start4(2)=1 ; count4(2)=jmt_global
         start4(3)=1 ; count4(3)=1
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar2d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif

!output uu ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(1:2)='uu'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(1:2), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

do i = 1,imt 
do j = 1,jmt 
do k = 1,km 
varvar3d4(i,j,k,1) = u (i,j,k) 
enddo 
enddo 
enddo 
         call local_to_global_4d(varvar3d4,varvar3d,km,1)
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif

!output vv ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(1:2)='vv'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(1:2), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

do i = 1,imt 
do j = 1,jmt 
do k = 1,km 
varvar3d4(i,j,k,1) = v (i,j,k) 
enddo 
enddo 
enddo 
         call local_to_global_4d(varvar3d4,varvar3d,km,1)
      !redefine the direction of v
         varvar3d = varvar3d * (-1)
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif

!output ww ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(1:2)='ww'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(1:2), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

do i = 1,imt 
do j = 1,jmt 
do k = 1,km 
varvar3d4(i,j,k,1) = ws (i,j,k) 
enddo 
enddo 
enddo 
         call local_to_global_4d(varvar3d4,varvar3d,km,1)
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif

!output tt ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(1:2)='tt'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(1:2), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

do i = 1,imt 
do j = 1,jmt 
do k = 1,km 
varvar3d4(i,j,k,1) = at (i,j,k,1) 
enddo 
enddo 
enddo 
         call local_to_global_4d(varvar3d4,varvar3d,km,1)
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif

!output ss ----------------------------------------------------------------------
if (mytid==0) then
   !get file name ----------------------
         fname(1:2)='ss'
   !create file ------------------------
         iert = nf_create(pname//fname, NF_CLOBBER, ncid)
         call check_err(iert)
   !define dimension -------------------
         iert = nf_def_dim(ncid, 'lon', lon_len, imt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lat', lat_len, jmt_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'lev', lev_len, km_dimid)
         call check_err(iert)
         iert = nf_def_dim(ncid, 'time', NF_UNLIMITED, ktime_dimid)
         call check_err(iert)         
   !define variables -------------------
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, (/imt_dimid/), lon_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, (/jmt_dimid/), lat_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, (/km_dimid/), lev_id)
         CALL check_err (iret)
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, (/ktime_dimid/), itime_id)
         CALL check_err (iret)
         iert = nf_def_var(ncid, fname(1:2), NF_REAL, 4, (/imt_dimid, jmt_dimid, km_dimid, ktime_dimid/), i_varid)
         call check_err(iert)      
   !assign attributes -------------------
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 5, 'depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'units', 21, 'days since '//stryr//'-01-01')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, itime_id, 'calendar', 8, 'standard')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, i_varid, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, i_varid, 'missing_value', NF_REAL, 1, spval)
         CALL check_err (iret)

         iert = nf_enddef(ncid)
         call check_err(iert)
   !value to variables ------------------
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_int (ncid, itime_id, vartime)
         CALL check_err (iret)
endif

do i = 1,imt 
do j = 1,jmt 
do k = 1,km 
varvar3d4(i,j,k,1) = at (i,j,k,2) 
enddo 
enddo 
enddo 
         call local_to_global_4d(varvar3d4,varvar3d,km,1)
      !alternate unit of ss
         varvar3d = varvar3d *1000 + 35 
         if (mytid==0) then
         where (viv_output_global < 0.5D0) varvar3d(:,:,:,1) = spval 
         start4(1)=1 ; count4(1)=lon_len
         start4(2)=1 ; count4(2)=lat_len
         start4(3)=1 ; count4(3)=lev_len 
         start4(4)=1 ; count4(4)=1
         iert = nf_put_vara_real(ncid, i_varid, start4, count4, varvar3d)
         call check_err(iert)
   ! close file ------------------------
         iert = nf_close(ncid)
         call check_err(iert)
         endif
!---------------------------end output dailymean----------------------------------------------------------------
endif 

#endif 

end subroutine sSAVE_hourly 



#if (defined NETCDF) || (defined ALL)
      SUBROUTINE check_err (iret)
#include <netcdf.inc>
      INTEGER :: iret
      IF (iret /= NF_NOERR) THEN
         PRINT *, nf_strerror (iret)
         STOP
      END IF
      END SUBROUTINE check_err
#endif

