      module data_licom
!
      include '/public/software/mathlib/libs-intel/netcdf/4.4.1/&
                include/netcdf.inc'
      !include'/THL7/home/liuhailong/netcdf-3.6.0-x64/include/netcdf.inc'
!
      integer,parameter :: imt=551,jmt=641,km=55,ntra=2,tmt=1
!
      real(kind=8),dimension(imt,jmt):: data_in
      real,dimension(imt,jmt):: data_out
      real lat(jmt),lon(imt),lev(km)
      integer,dimension(imt,jmt)::ind,indv
      real, parameter :: spval =1.0e+35

      character*2 cid,monid
      character*4 yearid
      integer dyear,dmonth
      integer nday,mday(12),mmday(12)
      data mday/0,31,59,90,120,151,181,212,243,273,304,334/
      data mmday/31,28,31,30,31,30,31,31,30,31,30,31/
      integer iid,i,j,k,ii
      
      integer ncid,iret,ncid1,ncid2,ncid3,ncid4,ncid5,ncid6,&
              ncid7,ncid8,ncid9,ncid10,ncid11,ncid12,ncid13
      character fname1*16,fname2*16,fname3*16,&
                  fname4*16,fname5*16,fname6*16,&
                  fname7*16,fname8*16,fname9*16,&
                  fname10*16,fname11*16,fname12*16,&
                  fname13*16
      character tt*8,dd*10
      integer lon_dim,lat_dim,lev_dim,time_dim,lev1_dim
      integer lon_len,lat_len,lev_len,time_len,lev1_len
      parameter (lon_len=551,lat_len=641,lev_len=1,&
                   lev1_len=km,time_len=1) 

      integer lat_id,lon_id,lev_id,time_id,ssh_id,uu_id,&
              tt_id,ss_id,vv_id,lev1_id,ww_id,&
              tx_id,ty_id,sw_id,lw_id,sh_id,lh_id,fr_id
      integer  lat_rank,lon_rank,lev_rank,time_rank,&
              ssh_rank,lev1_rank,uu_rank,tt_rank,ss_rank,vv_rank,ww_rank,&
              tx_rank,ty_rank,sw_rank,lw_rank,sh_rank,lh_rank,fr_rank

      parameter (lat_rank=1,lon_rank=1,lev_rank=1,time_rank=1,&
                 lev1_rank=1,uu_rank=4,tt_rank=4,ss_rank=4,vv_rank=4,&
                 ssh_rank=4,ww_rank=4, tx_rank=4,ty_rank=4,&
                 sw_rank=4,lw_rank=4, sh_rank=4,lh_rank=4,fr_rank=4)
      integer lat_dims(lat_rank),lon_dims(lon_rank), &
              lev_dims(lev_rank),time_dims(time_rank),&
             ssh_dims(ssh_rank), lev1_dims(lev1_rank),uu_dims(uu_rank),&
              vv_dims(vv_rank), tt_dims(tt_rank), ss_dims(ss_rank),&
              ww_dims(ww_rank),&
              tx_dims(tx_rank), ty_dims(ty_rank), sw_dims(sw_rank),&
              lw_dims(lw_rank), sh_dims(sh_rank), lh_dims(ss_rank),&
              fr_dims(fr_rank)

      integer start1(1),count1(1)
      integer start3(3),count3(3)
      integer start4(4),count4(4)
!      real xx(lon_len),yy(lat_len)

      end module data_licom
!
      program fort22
      use data_licom
      IMPLICIT NONE

! read index
      iret=nf_open('ind_final_eas_uv.nc',nf_nowrite,ncid)
      call check_err (iret)
      write(*,*)'ready to read INDEX.DATA_3602X1683X55.nc!'
      start1(1)=1 ; count1(1)=imt
      iret=nf_get_vara_real(ncid,   1,start1,count1,lon)
      call check_err (iret)
      write(*,*)'read lon ok' 

      start1(1)=1 ; count1(1)=jmt
      iret=nf_get_vara_real(ncid,   2,start1,count1,lat)
      call check_err (iret)
       write(*,*)'read lat ok'
 
      start1(1)=1 ; count1(1)=km
      iret=nf_get_vara_real(ncid,   3,start1,count1,lev)
      call check_err (iret)
      write(*,*)'read lev ok'

!set year, month
      do dyear =2016,2016
!      dyear =33
!
!      do dmonth=1,12
      do dmonth=2,2

!      do iid=1,mmday(dmonth)
      do iid=29,29!mmday(dmonth)
        write(cid,'(I2.2)') iid
        write(monid,'(I2.2)') dmonth
        write(yearid,'(I4.4)') dyear
!
      fname1="z0-"//yearid//"-"//monid//"-"//cid//".nc"
      fname2="uu-"//yearid//"-"//monid//"-"//cid//".nc"
      fname3="vv-"//yearid//"-"//monid//"-"//cid//".nc"
      fname4="tt-"//yearid//"-"//monid//"-"//cid//".nc"
      fname5="ss-"//yearid//"-"//monid//"-"//cid//".nc"
      fname6="ww-"//yearid//"-"//monid//"-"//cid//".nc"
      fname7="tx-"//yearid//"-"//monid//"-"//cid//".nc"
      fname8="ty-"//yearid//"-"//monid//"-"//cid//".nc"
      fname9="sw-"//yearid//"-"//monid//"-"//cid//".nc"
      fname10="lw-"//yearid//"-"//monid//"-"//cid//".nc"
      fname11="sh-"//yearid//"-"//monid//"-"//cid//".nc"
      fname12="lh-"//yearid//"-"//monid//"-"//cid//".nc"
      fname13="fr-"//yearid//"-"//monid//"-"//cid//".nc"
!
! ssh
      write(*,*) fname1
      iret = nf_create(fname1, NF_CLOBBER, ncid1)
      call check_err(iret)
 
      iret = nf_def_dim(ncid1,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid1,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid1,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim (ncid1, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid1,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid1,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(ncid1,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid1, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)

      ssh_dims (4) = time_dim
      ssh_dims (3) = lev_dim
      ssh_dims (2) = lat_dim
      ssh_dims (1) = lon_dim
      iret = nf_def_var (ncid1, 'ssh', NF_REAL, ssh_rank, ssh_dims, ssh_id)
      CALL check_err (iret)


      iret = nf_put_att_text(ncid1, lat_id, 'long_name', 21, 'latitude on T grids' )  
      call check_err(iret)
      iret = nf_put_att_text(ncid1, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid1, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid1, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid1, lev_id, 'long_name', 18, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid1, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid1, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid1, time_id, 'units', 21, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(ncid1, time_id, 'calendar', 7, '365days')
      call check_err(iret)

      iret = nf_put_att_text (ncid1, ssh_id, 'long_name', 18, 'sea surface height')
      call check_err (iret)
      iret = nf_put_att_text (ncid1, ssh_id, 'units', 5, 'meter')
      call check_err (iret)
      iret = nf_put_att_real (ncid1, ssh_id, 'missing_value', NF_REAL, 1, spval) !M
      call check_err (iret)
      write(*,*)'assign attributes OK'

      iret = NF_PUT_ATT_TEXT (NCID1, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID1, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)
      iret = nf_enddef(ncid1)
      call check_err(iret)

!    write dimensions
      iret = nf_put_var_real(ncid1, lon_id, lon)
      call check_err(iret)
       write(*,*)'prepare data for storing  lon ok'
 
      iret = nf_put_var_real(ncid1, lat_id, lat)
      call check_err(iret)
      write(*,*)'prepare data for storing  lat ok'

      iret = nf_put_var_real(ncid1, lev_id, lev)
      call check_err(iret)
      write(*,*)'prepare data for storing  lev ok'

      start1(1)=1
      count1(1)=time_len
      nday= (dyear-1)*365+ mday(dmonth)+iid-1
      iret = nf_put_vara_int(ncid1, time_id,start1,count1,nday)
      call check_err(iret)
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=1 ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)

      call index_deal(ind,1,indv)


!
! Open fort.22
      print*,'./fort.22.'//yearid//'-'//monid//'-'//cid//''
      open(22,file='./fort.22.'//yearid//'-'//monid//'-'//cid//'',form='unformatted',status='unknown')
      write(*,*)cid
      write(*,*)'Start to Read'
      read(22)data_in
      write(*,*)'Read h0 OK'
      
      do j=1,jmt
         do i=1,imt
         data_out(i,j)=data_in(i,j)
         end do
      end do

!      do j=1,240
!         do i=1,imt
!         data_out(i,j)=spval
!         end do 
!      end do 
!

      do j=1,jmt
         do i=1,imt
         if(ind(i,j)==0) data_out(i,j)=spval
         end do 
      end do 


         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4) =time_len
      iret = nf_put_vara_real(ncid1,  ssh_id, start4, count4,data_out)
      call check_err(iret)
      write(*,*)'prepare data for storing  ssh ok'

      iret=nf_close(ncid1)
      call check_err(iret)
!
!
      iret = nf_create(fname2, NF_CLOBBER, ncid2)
      call check_err(iret)
!
      call nc_defdim(ncid2)
!
         uu_dims (4) = time_dim
         uu_dims (3) = lev1_dim
         uu_dims (2) = lat_dim
         uu_dims (1) = lon_dim
         iret = nf_def_var (ncid2, 'uu', NF_REAL, uu_rank, uu_dims, uu_id)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid2, uu_id, 'long_name', 14, 'zonal velocity')
 
      CALL check_err (iret)
      iret = nf_put_att_text (ncid2, uu_id, 'units', 5, ' m/s ')
      CALL check_err (iret)
      iret = nf_put_att_real (ncid2, uu_id, '_FillValue', NF_REAL, 1, spval)
      CALL check_err (iret)

      iret = NF_PUT_ATT_TEXT (NCID2, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID2, NF_GLOBAL, 'history', 20, tt//'  '//dd)
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID2, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)

      iret = nf_enddef(ncid2)
      call check_err(iret)

      iret = nf_put_var_real(ncid2, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid2, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid2, lev1_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid2, time_id,start1,count1,nday)
      call check_err(iret)

      do k=1 ,km
      read(22) data_in
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=k ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,k,indv)
!
         do j=1,jmt
         do i=1,imt
              data_out(i,j)= data_in(i,j)
         end do
         end do
!
!           do j=1,240
!               do i=1,imt
!               data_out(i,j)=spval
!              end do 
!           end do 

           do j=1,jmt
               do i=1,imt
               if(indv(i,j)==0)then
                 data_out(i,j)=spval
               else  
               end if
              end do 
           end do 


         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= k
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= 1
         count4 (4) =time_len
         iret = nf_put_vara_real(ncid2,  uu_id, start4, count4,data_out)
         call check_err(iret)

      end do
      iret = nf_close (ncid2)
      call check_err (iret)
      write(*,*)'uu ok'
!
!
      iret = nf_create(fname3, NF_CLOBBER, ncid3)
      call check_err(iret)
!
      call nc_defdim(ncid3)
!
         vv_dims (4) = time_dim
         vv_dims (3) = lev1_dim
         vv_dims (2) = lat_dim
         vv_dims (1) = lon_dim
         iret = nf_def_var (ncid3, 'vv', NF_REAL, vv_rank, vv_dims, vv_id)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid3, vv_id, 'long_name', 19, 'meridional velocity')
!
         CALL check_err (iret)
         iret = nf_put_att_text (ncid3, vv_id, 'units', 5, ' m/s ')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid3, vv_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
      iret = NF_PUT_ATT_TEXT (NCID3, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID3, NF_GLOBAL, 'history', 20, tt//'  '//dd)
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID3, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)

      iret = nf_enddef(ncid3)
      call check_err(iret)


!----------------------------------------------------------
!     prepare data for storing
!----------------------------------------------------------

      iret = nf_put_var_real(ncid3, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid3, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid3, lev1_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid3, time_id,start1,count1,nday)
      call check_err(iret)
!
      do k=1 ,km
         read(22) data_in
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=k ; count3(3)=1
 
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,k,indv)
!
         do j=1,jmt
         do i=1,imt
              data_out(i,j)= data_in(i,j)
         end do
         end do
!
!           do j=1,240
!               do i=1,imt
!               data_out(i,j)=spval
!              end do 
!           end do 

           do j=1,jmt
               do i=1,imt
               if(indv(i,j)==0)then
                 data_out(i,j)=spval
               else  
               end if
              end do 
           end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= k
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= 1
         count4 (4) =time_len
         iret = nf_put_vara_real(ncid3,  vv_id, start4, count4,data_out)
         call check_err(iret)

      end do
      iret = nf_close (ncid3)
      call check_err (iret)
      write(*,*)'vv ok'
!
      iret = nf_create(fname4, NF_CLOBBER, ncid4)
      call check_err(iret)
!
      call nc_defdim(ncid4)
!
         tt_dims (4) = time_dim
         tt_dims (3) = lev1_dim
         tt_dims (2) = lat_dim
         tt_dims (1) = lon_dim
         iret = nf_def_var (ncid4, 'tt', NF_REAL, tt_rank, tt_dims, tt_id)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid4, tt_id, 'long_name', 11, 'temperature')
!
         CALL check_err (iret)
         iret = nf_put_att_text (ncid4, tt_id, 'units', 5, ' C  ')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid4, tt_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
      iret = NF_PUT_ATT_TEXT (NCID4, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID4, NF_GLOBAL, 'history', 20, tt//'  '//dd)
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID4, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)

      iret = nf_enddef(ncid4)
      call check_err(iret)


!----------------------------------------------------------
!     prepare data for storing
!----------------------------------------------------------

      iret = nf_put_var_real(ncid4, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid4, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid4, lev1_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid4, time_id,start1,count1,nday)
      call check_err(iret)
!
      do k=1 ,km
         read(22) data_in
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=k ; count3(3)=1
 
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,k,indv)
!
         do j=1,jmt
         do i=1,imt
              data_out(i,j)= data_in(i,j)
         end do
         end do
!
!           do j=1,240
!               do i=1,imt
!               data_out(i,j)=spval
!              end do 
!           end do 

           do j=1,jmt
               do i=1,imt
               if(ind(i,j)==0)then
                 data_out(i,j)=spval
               else  
               end if
              end do 
           end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= k
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= 1
         count4 (4) =time_len
         iret = nf_put_vara_real(ncid4,  tt_id, start4, count4,data_out)
         call check_err(iret)

      end do
      iret = nf_close (ncid4)
      call check_err (iret)
      write(*,*)'tt ok'

! ss
      iret = nf_create(fname5, NF_CLOBBER, ncid5)
      call check_err(iret)
!
      call nc_defdim(ncid5)
!
         ss_dims (4) = time_dim
         ss_dims (3) = lev1_dim
         ss_dims (2) = lat_dim
         ss_dims (1) = lon_dim
         iret = nf_def_var (ncid5, 'ss', NF_REAL, ss_rank, ss_dims, ss_id)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid5, ss_id, 'long_name', 8, 'salinity')
!
         CALL check_err (iret)
         iret = nf_put_att_text (ncid5, ss_id, 'units', 5, ' psu ')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid5, ss_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
      iret = NF_PUT_ATT_TEXT (NCID5, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID5, NF_GLOBAL, 'history', 20, tt//'  '//dd)
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID5, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)

      iret = nf_enddef(ncid5)
      call check_err(iret)

      iret = nf_put_var_real(ncid5, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid5, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid5, lev1_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid5, time_id,start1,count1,nday)
      call check_err(iret)

!
      do k=1 ,km
         read(22) data_in
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=k ; count3(3)=1
 
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,k,indv)
!
         do j=1,jmt
         do i=1,imt
              data_out(i,j)= data_in(i,j)*1000.D0+35.D0
         end do
         end do
!
!           do j=1,240
!               do i=1,imt
!               data_out(i,j)=spval
!              end do 
!           end do 

           do j=1,jmt
               do i=1,imt
               if(ind(i,j)==0)then
                 data_out(i,j)=spval
               else  
               end if
              end do 
           end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= k
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= 1
         count4 (4) =time_len
         iret = nf_put_vara_real(ncid5,  ss_id, start4, count4,data_out)
         call check_err(iret)

      end do
      iret = nf_close (ncid5)
      call check_err (iret)
      write(*,*)'ss ok'

! ww
      write(*,*) "begin creat fname6 OK"
      iret = nf_create(fname6, NF_CLOBBER, ncid6)
      call check_err(iret)
      write(*,*) "creat fname6 OK"
!
      call nc_defdim(ncid6)
!
      ww_dims (4) = time_dim
      ww_dims (3) = lev1_dim
      ww_dims (2) = lat_dim
      ww_dims (1) = lon_dim
      iret = nf_def_var (ncid6, 'ww', NF_REAL, ww_rank, ww_dims, ww_id)
      CALL check_err (iret)
      iret = nf_put_att_text (ncid6, ww_id, 'long_name', 17, 'Vertical Velocity')
      CALL check_err (iret)
      iret = nf_put_att_text (ncid6, ww_id, 'units', 5, ' m/s ')
      CALL check_err (iret)
      iret = nf_put_att_real (ncid6, ww_id, '_FillValue', NF_REAL, 1, spval)
      CALL check_err (iret)
!
      iret = NF_PUT_ATT_TEXT (NCID6, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID6, NF_GLOBAL, 'history', 20, tt//'  '//dd)
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID6, NF_GLOBAL, 'source', 35, &
               'LASG/IAP Climate system Ocean Model')
      call check_err(iret)

      iret = nf_enddef(ncid6)
      call check_err(iret)
      write(*,*) "define fname6 OK"

      iret = nf_put_var_real(ncid6, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid6, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid6, lev1_id, lev)
      call check_err(iret)
      write(*,*) "put dimension fname6 OK"

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid6, time_id,start1,count1,nday)
      call check_err(iret)
      write(*,*) "put time fname6 OK"
!
      do k=1 ,km
      write(*,*) "read 22 fname6 OK k=", k
         read(22) data_in
      write(*,*) "read 22 fname6 OK"
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=k ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      write(*,*) "read index fname6 OK"
!
         do j=1,jmt
         do i=1,imt
              data_out(i,j)= data_in(i,j)
         end do
         end do
!
!           do j=1,240
!               do i=1,imt
!               data_out(i,j)=spval
!              end do 
!           end do 

           do j=1,jmt
               do i=1,imt
               if(ind(i,j)==0)then
                 data_out(i,j)=spval
               else  
               end if
              end do 
           end do 
      write(*,*) "deal data fname6 OK"
!
         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= k
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= 1
         count4 (4) =time_len
         iret = nf_put_vara_real(ncid6,  ww_id, start4, count4,data_out)
         call check_err(iret)
      write(*,*) "put data fname6 OK"

      end do
      iret = nf_close (ncid6)
      call check_err (iret)
      write(*,*) "close fname6 OK"

!---------------------------tx
      write(*,*) fname7
      iret = nf_create(fname7, NF_CLOBBER, ncid7)
      call check_err(iret)
 
      iret = nf_def_dim(ncid7,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid7,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid7,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim (ncid7, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid7,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid7,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(ncid7,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid7, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)

      tx_dims (4) = time_dim
      tx_dims (3) = lev_dim
      tx_dims (2) = lat_dim
      tx_dims (1) = lon_dim
      iret = nf_def_var (ncid7, 'tx', NF_REAL, tx_rank, tx_dims, tx_id)
      CALL check_err (iret)


      iret = nf_put_att_text(ncid7, lat_id, 'long_name', 21, 'latitude on T grids' )  
      call check_err(iret)
      iret = nf_put_att_text(ncid7, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid7, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid7, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid7, lev_id, 'long_name', 18, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid7, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid7, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid7, time_id, 'units', 21, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(ncid7, time_id, 'calendar', 7, '365days')
      call check_err(iret)

      iret = nf_put_att_text (ncid7, tx_id, 'long_name', 17, 'zonal wind stress')
      call check_err (iret)
      iret = nf_put_att_text (ncid7, tx_id, 'units', 2, 'Pa')
      call check_err (iret)
      iret = nf_put_att_real (ncid7, tx_id, 'missing_value', NF_REAL, 1, spval) !M
      call check_err (iret)

      iret = NF_PUT_ATT_TEXT (NCID7, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID7, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)
      iret = nf_enddef(ncid7)
      call check_err(iret)

      iret = nf_put_var_real(ncid7, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid7, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid7, lev_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid7, time_id,start1,count1,nday)
      call check_err(iret)
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=1 ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,1,indv)
!
      read(22)data_in
      
      do j=1,jmt
         do i=1,imt
         data_out(i,j)=data_in(i,j)
         end do
      end do

!      do j=1,240
!         do i=1,imt
!         data_out(i,j)=spval
!         end do 
!      end do 

      do j=1,jmt
         do i=1,imt
         if(indv(i,j)==0) data_out(i,j)=spval
         end do 
      end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4) =time_len
      iret = nf_put_vara_real(ncid7,  tx_id, start4, count4,data_out)
      call check_err(iret)

      iret=nf_close(ncid7)
      call check_err(iret)

!---------------------------ty
      write(*,*) fname8
      iret = nf_create(fname8, NF_CLOBBER, ncid8)
      call check_err(iret)
      write(*,*) fname8
 
      iret = nf_def_dim(ncid8,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid8,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid8,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim (ncid8, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)
      write(*,*) fname8

      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid8,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid8,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(ncid8,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid8, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)
      write(*,*) fname8

      ty_dims (4) = time_dim
      ty_dims (3) = lev_dim
      ty_dims (2) = lat_dim
      ty_dims (1) = lon_dim
      iret = nf_def_var (ncid8, 'ty', NF_REAL, ty_rank, ty_dims, ty_id)
      CALL check_err (iret)
      write(*,*) fname8


      iret = nf_put_att_text(ncid8, lat_id, 'long_name', 21, 'latitude on T grids' )  
      call check_err(iret)
      iret = nf_put_att_text(ncid8, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid8, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid8, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid8, lev_id, 'long_name', 18, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid8, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid8, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid8, time_id, 'units', 21, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(ncid8, time_id, 'calendar', 7, '365days')
      call check_err(iret)

      iret = nf_put_att_text (ncid8, ty_id, 'long_name', 22, 'meridional wind stress')
      call check_err (iret)
      iret = nf_put_att_text (ncid8, ty_id, 'units', 2, 'Pa')
      call check_err (iret)
      iret = nf_put_att_real (ncid8, ty_id, 'missing_value', NF_REAL, 1, spval) !M
      call check_err (iret)

      iret = NF_PUT_ATT_TEXT (NCID8, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID8, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)
      iret = nf_enddef(ncid8)
      call check_err(iret)

      write(*,*) fname8
      iret = nf_put_var_real(ncid8, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid8, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid8, lev_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid8, time_id,start1,count1,nday)
      call check_err(iret)
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=1 ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,1,indv)
      write(*,*) fname8
!
      read(22)data_in
      
      do j=1,jmt
         do i=1,imt
         data_out(i,j)=data_in(i,j)*-1
         end do
      end do

!      do j=1,240
!         do i=1,imt
!         data_out(i,j)=spval
!         end do 
!      end do 

      do j=1,jmt
         do i=1,imt
         if(indv(i,j)==0) data_out(i,j)=spval
         end do 
      end do 
      write(*,*) fname8
!
         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4) =time_len
      iret = nf_put_vara_real(ncid8,  ty_id, start4, count4,data_out)
      call check_err(iret)
      write(*,*) fname8

      iret=nf_close(ncid8)
      call check_err(iret)

!---------------------------sw
      write(*,*) fname9
      iret = nf_create(fname9, NF_CLOBBER, ncid9)
      call check_err(iret)
 
      iret = nf_def_dim(ncid9,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid9,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid9,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim (ncid9, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid9,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid9,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(ncid9,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid9, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)

      sw_dims (4) = time_dim
      sw_dims (3) = lev_dim
      sw_dims (2) = lat_dim
      sw_dims (1) = lon_dim
      iret = nf_def_var (ncid9, 'sw', NF_REAL, sw_rank, sw_dims, sw_id)
      CALL check_err (iret)


      iret = nf_put_att_text(ncid9, lat_id, 'long_name', 21, 'latitude on T grids' )  
      call check_err(iret)
      iret = nf_put_att_text(ncid9, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid9, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid9, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid9, lev_id, 'long_name', 19, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid9, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid9, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid9, time_id, 'units', 21, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(ncid9, time_id, 'calendar', 7, '365days')
      call check_err(iret)

      iret = nf_put_att_text (ncid9, sw_id, 'long_name', 15, 'Solar Radiation')
      call check_err (iret)
      iret = nf_put_att_text (ncid9, sw_id, 'units', 5, 'W/m^2')
      call check_err (iret)
      iret = nf_put_att_real (ncid9, sw_id, 'missing_value', NF_REAL, 1, spval) !M
      call check_err (iret)

      iret = NF_PUT_ATT_TEXT (NCID9, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID9, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)
      iret = nf_enddef(ncid9)
      call check_err(iret)

      iret = nf_put_var_real(ncid9, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid9, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid9, lev_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid9, time_id,start1,count1,nday)
      call check_err(iret)
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=1 ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,1,indv)
!
      read(22)data_in
      
      do j=1,jmt
         do i=1,imt
         data_out(i,j)=data_in(i,j)
         end do
      end do

!      do j=1,240
!         do i=1,imt
!         data_out(i,j)=spval
!         end do 
!      end do 

      do j=1,jmt
         do i=1,imt
         if(ind(i,j)==0) data_out(i,j)=spval
         end do 
      end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4) =time_len
      iret = nf_put_vara_real(ncid9,  sw_id, start4, count4,data_out)
      call check_err(iret)

      iret=nf_close(ncid9)
      call check_err(iret)


!---------------------------lw
      write(*,*) fname10
      iret = nf_create(fname10, NF_CLOBBER, ncid10)
      call check_err(iret)
 
      iret = nf_def_dim(ncid10,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid10,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid10,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim (ncid10, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid10,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid10,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(ncid10,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid10, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)

      lw_dims (4) = time_dim
      lw_dims (3) = lev_dim
      lw_dims (2) = lat_dim
      lw_dims (1) = lon_dim
      iret = nf_def_var (ncid10, 'lw', NF_REAL, lw_rank, lw_dims, lw_id)
      CALL check_err (iret)


      iret = nf_put_att_text(ncid10, lat_id, 'long_name', 21, 'latitude on T grids' )  
      call check_err(iret)
      iret = nf_put_att_text(ncid10, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid10, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid10, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid10, lev_id, 'long_name', 19, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid10, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid10, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid10, time_id, 'units', 21, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(ncid10, time_id, 'calendar', 7, '365days')
      call check_err(iret)

      iret = nf_put_att_text (ncid10, lw_id, 'long_name', 18, 'Longwave Radiation')
      call check_err (iret)
      iret = nf_put_att_text (ncid10, lw_id, 'units', 5, 'W/m^2')
      call check_err (iret)
      iret = nf_put_att_real (ncid10, lw_id, 'missing_value', NF_REAL, 1, spval) !M
      call check_err (iret)

      iret = NF_PUT_ATT_TEXT (NCID10, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID10, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)
      iret = nf_enddef(ncid10)
      call check_err(iret)

      iret = nf_put_var_real(ncid10, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid10, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid10, lev_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid10, time_id,start1,count1,nday)
      call check_err(iret)
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=1 ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,1,indv)
!
      read(22)data_in
      
      do j=1,jmt
         do i=1,imt
         data_out(i,j)=data_in(i,j)
         end do
      end do

!      do j=1,240
!         do i=1,imt
!         data_out(i,j)=spval
!         end do 
!      end do 

      do j=1,jmt
         do i=1,imt
         if(ind(i,j)==0) data_out(i,j)=spval
         end do 
      end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4) =time_len
      iret = nf_put_vara_real(ncid10,  lw_id, start4, count4,data_out)
      call check_err(iret)

      iret=nf_close(ncid10)
      call check_err(iret)


!---------------------------sh
      write(*,*) fname11
      iret = nf_create(fname11, NF_CLOBBER, ncid11)
      call check_err(iret)
 
      iret = nf_def_dim(ncid11,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid11,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid11,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim (ncid11, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid11,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid11,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(ncid11,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid11, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)

      sh_dims (4) = time_dim
      sh_dims (3) = lev_dim
      sh_dims (2) = lat_dim
      sh_dims (1) = lon_dim
      iret = nf_def_var (ncid11, 'sh', NF_REAL, sh_rank, sh_dims, sh_id)
      CALL check_err (iret)


      iret = nf_put_att_text(ncid11, lat_id, 'long_name', 21, 'latitude on T grids' )  
      call check_err(iret)
      iret = nf_put_att_text(ncid11, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid11, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid11, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid11, lev_id, 'long_name', 19, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid11, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid11, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid11, time_id, 'units', 21, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(ncid11, time_id, 'calendar', 7, '365days')
      call check_err(iret)

      iret = nf_put_att_text (ncid11, sh_id, 'long_name', 18, 'Sensible Heat Flux')
      call check_err (iret)
      iret = nf_put_att_text (ncid11, sh_id, 'units', 5, 'W/m^2')
      call check_err (iret)
      iret = nf_put_att_real (ncid11, sh_id, 'missing_value', NF_REAL, 1, spval) !M
      call check_err (iret)

      iret = NF_PUT_ATT_TEXT (NCID11, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID11, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)
      iret = nf_enddef(ncid11)
      call check_err(iret)

      iret = nf_put_var_real(ncid11, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid11, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid11, lev_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid11, time_id,start1,count1,nday)
      call check_err(iret)
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=1 ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,1,indv)
!
      read(22)data_in
      
      do j=1,jmt
         do i=1,imt
         data_out(i,j)=data_in(i,j)
         end do
      end do

!      do j=1,240
!         do i=1,imt
!         data_out(i,j)=spval
!         end do 
!      end do 

      do j=1,jmt
         do i=1,imt
         if(ind(i,j)==0) data_out(i,j)=spval
         end do 
      end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4) =time_len
      iret = nf_put_vara_real(ncid11,  sh_id, start4, count4,data_out)
      call check_err(iret)

      iret=nf_close(ncid11)
      call check_err(iret)


!---------------------------lh
      write(*,*) fname12
      iret = nf_create(fname12, NF_CLOBBER, ncid12)
      call check_err(iret)
 
      iret = nf_def_dim(ncid12,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid12,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid12,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim (ncid12, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid12,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid12,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(ncid12,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid12, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)

      lh_dims (4) = time_dim
      lh_dims (3) = lev_dim
      lh_dims (2) = lat_dim
      lh_dims (1) = lon_dim
      iret = nf_def_var (ncid12, 'lh', NF_REAL, lh_rank, lh_dims, lh_id)
      CALL check_err (iret)


      iret = nf_put_att_text(ncid12, lat_id, 'long_name', 21, 'latitude on T grids' )  
      call check_err(iret)
      iret = nf_put_att_text(ncid12, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid12, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid12, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid12, lev_id, 'long_name', 19, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid12, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid12, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid12, time_id, 'units', 21, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(ncid12, time_id, 'calendar', 7, '365days')
      call check_err(iret)

      iret = nf_put_att_text (ncid12, lh_id, 'long_name', 16, 'Latent Heat Flux')
      call check_err (iret)
      iret = nf_put_att_text (ncid12, lh_id, 'units', 5, 'W/m^2')
      call check_err (iret)
      iret = nf_put_att_real (ncid12, lh_id, 'missing_value', NF_REAL, 1, spval) !M
      call check_err (iret)

      iret = NF_PUT_ATT_TEXT (NCID12, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID12, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)
      iret = nf_enddef(ncid12)
      call check_err(iret)

      iret = nf_put_var_real(ncid12, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid12, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid12, lev_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid12, time_id,start1,count1,nday)
      call check_err(iret)
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=1 ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,1,indv)
!
      read(22)data_in
      
      do j=1,jmt
         do i=1,imt
         data_out(i,j)=data_in(i,j)
         end do
      end do

!      do j=1,240
!         do i=1,imt
!         data_out(i,j)=spval
!         end do 
!      end do 

      do j=1,jmt
         do i=1,imt
         if(ind(i,j)==0) data_out(i,j)=spval
         end do 
      end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4) =time_len
      iret = nf_put_vara_real(ncid12,  lh_id, start4, count4,data_out)
      call check_err(iret)

      iret=nf_close(ncid12)
      call check_err(iret)

!---------------------------fr
      write(*,*) fname13
      iret = nf_create(fname13, NF_CLOBBER, ncid13)
      call check_err(iret)
 
      iret = nf_def_dim(ncid13,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid13,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid13,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim (ncid13, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)

      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid13,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid13,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(ncid13,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid13, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)

      fr_dims (4) = time_dim
      fr_dims (3) = lev_dim
      fr_dims (2) = lat_dim
      fr_dims (1) = lon_dim
      iret = nf_def_var (ncid13, 'fr', NF_REAL, fr_rank, fr_dims, fr_id)
      CALL check_err (iret)


      iret = nf_put_att_text(ncid13, lat_id, 'long_name', 21, 'latitude on T grids' )  
      call check_err(iret)
      iret = nf_put_att_text(ncid13, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid13, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid13, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid13, lev_id, 'long_name', 19, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(ncid13, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(ncid13, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid13, time_id, 'units', 21, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(ncid13, time_id, 'calendar', 7, '365days')
      call check_err(iret)

      iret = nf_put_att_text (ncid13, fr_id, 'long_name', 11, 'Fresh Water')
      call check_err (iret)
      iret = nf_put_att_text (ncid13, fr_id, 'units', 3, 'm/s')
      call check_err (iret)
      iret = nf_put_att_real (ncid13, fr_id, 'missing_value', NF_REAL, 1, spval) !M
      call check_err (iret)

      iret = NF_PUT_ATT_TEXT (NCID13, NF_GLOBAL, 'title', 4, 'test')
      call check_err(iret)
      iret = NF_PUT_ATT_TEXT (NCID13, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
      call check_err(iret)
      iret = nf_enddef(ncid13)
      call check_err(iret)

      iret = nf_put_var_real(ncid13, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid13, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid13, lev_id, lev)
      call check_err(iret)

      start1(1)=1
      count1(1)=time_len
      iret = nf_put_vara_int(ncid13, time_id,start1,count1,nday)
      call check_err(iret)
!
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt
      start3(3)=1 ; count3(3)=1
      iret=nf_get_vara_int(ncid, 4,start3,count3,ind)
      call check_err (iret)
      call index_deal(ind,1,indv)
!
      read(22)data_in
      
      do j=1,jmt
         do i=1,imt
         data_out(i,j)=data_in(i,j)
         end do
      end do

!      do j=1,240
!         do i=1,imt
!         data_out(i,j)=spval
!         end do 
!      end do 

      do j=1,jmt
         do i=1,imt
         if(ind(i,j)==0) data_out(i,j)=spval
         end do 
      end do 

         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4) =time_len
      iret = nf_put_vara_real(ncid13, fr_id, start4, count4,data_out)
      call check_err(iret)

      iret=nf_close(ncid13)
      call check_err(iret)

      close(22)
      end do

      end do
end do
      iret = nf_close (ncid)
      call check_err (iret)
!
      end 

!-----------------------------------------------------------------------
      subroutine check_err(iret)
!-----------------------------------------------------------------------
      include '/public/software/mathlib/libs-intel/netcdf/4.4.1/&
                include/netcdf.inc'
      !include '/THL7/home/liuhailong/netcdf-3.6.0-x64/include/netcdf.inc'
      integer iret
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end subroutine check_err

!
      subroutine nc_defdim(nncid)
      use data_licom
      integer :: nncid
!
      iret = nf_def_dim(nncid ,  'lat',  lat_len,  lat_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid ,  'lon',  lon_len,  lon_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid,  'lev1',  lev1_len,  lev1_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid, 'time', time_len, time_dim)
      call check_err(iret)
!
      lon_dims(1) = lon_dim
      iret = nf_def_var(nncid,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(nncid,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
      call check_err(iret)
      lev1_dims(1) = lev1_dim
      iret = nf_def_var(nncid,  'lev1', NF_REAL,  lev1_rank,  lev1_dims,  lev1_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(nncid, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)


      iret = nf_put_att_text(nncid, lat_id, 'long_name', 21, "latitude (on T grids)" )
      call check_err(iret)
      iret = nf_put_att_text(nncid, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lev1_id, 'long_name', 18, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lev1_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'units', 23, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'calendar', 7, '365days')
      call check_err(iret)
!
      return
      end
!
      subroutine index_deal(ind,k,indv)
      integer,parameter :: imt=551,jmt=641
!        
!          2018-09
!      integer,parameter :: imt=521,jmt=551
!
      integer,dimension(imt,jmt)::ind,indv
      integer i,j,k

!      if ( k < 37) then
!      do j=1512,1524
!      do i=3047,3049
!      ind(i,j)=1.0
!      end do
!      end do
!      end if
!
!      do j=956,964
!      do i=1461,1496
!      ind(i,j)=0.0
!      end do
!      end do
!
!      do j=827,835
!      do i=1251,1261
!      ind(i,j)=0.0
!      end do
!      end do
!
!      do j=850,859
!      do i=1175,1192
!      ind(i,j)=0.0
!      end do
!      end do
!
!      do j=935,942
!      do i=1516,1528
!      ind(i,j)=0.0
!      end do
!      end do
!
!      if (k < 42) then
!      do j=921,929
!      do i=1498,1510
!      ind(i,j)=1.0
!      end do
!      end do
!      end if
!
      do j = 1,jmt-1
         do i = 2,imt
         indv(I,J)=ind (i-1,j)*ind(i-1,j+1)*ind(i,j)*IND(i,j+1)
         end do
         indv(1,j)=indv(imt-1,j)
      END DO
      do i = 1,imt
         indv (i,jmt)= 0.0
      end do

      return
      end
