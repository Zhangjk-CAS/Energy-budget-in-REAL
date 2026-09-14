import numpy as np
import xarray as xr
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.dates as mdates
from datetime import datetime, timedelta

# lon_range=[95.3, 149.8];lat_range=[-15.7,47.8]
lon_range=[100, 150];lat_range=[-16,48]
# lon_range=[105, 125];lat_range=[0, 25]
# lon_range=[105, 140];lat_range=[0,40]

# key_strs = ['uu2014','uu2015','uu2016','uu2017','uu2018']
key_strs = ['uu-2014','uu-2015','uu-2016','uu-2017','uu-2018']
vn = 'uu'
vn2= 'vv'
vn3= 'tx'
vn4= 'ty'

# pnstrs = [f'edgeuaemke-jl']
# paths = ["../output-jra-licom3/fort.22/"]
pnstrs = [f'edgeuaemke-glorys-all']
paths = ["../exe/boundary-data-glorys-daily/"]
for nforce in range(1):
    pnstr = pnstrs[nforce]
    path = paths[nforce]
    
    nowdir = os.path.abspath(path)
    number_day = 0 
    regionalmean = np.zeros((2400,1))
    weighted_averages = np.zeros((55))
    meansst = np.zeros((2400,1))
    dz = np.zeros((55))

    # dsxm = xr.open_dataset(path + 'tx-2014-2018.nc', decode_times=False)
    # dsym = xr.open_dataset(path + 'ty-2014-2018.nc', decode_times=False)
    dsum = xr.open_dataset(path + 'ave/uu-2014-2018.nc', decode_times=False)
    dsvm = xr.open_dataset(path + 'ave/vv-2014-2018.nc', decode_times=False)
    # xm = dsxm['tx'].values[0,0,:,:]
    # ym = dsym['ty'].values[0,0,:,:]
    um = dsum['uu'].values[0,1:,:,:]
    vm = dsvm['vv'].values[0,1:,:,:]
    depth = dsum['depth'].values[1:]
    ndep = len(depth)
    dep3d = np.broadcast_to(depth.reshape(-1, 1, 1), um.shape)
    dz_mid = (dep3d[1:, :, :] + dep3d[:-1, :, :]) / 2
    dz = np.zeros_like(dep3d)
    dz[0, :, :] = depth[0] / 2
    dz[1:-1, :, :] = dz_mid[1:, :, :] - dz_mid[:-1, :, :]
    dz[-1, :, :] = (depth[-1] - dz_mid[-1, :, :]) * 2


    for item in sorted(os.listdir(nowdir)):
        fn = []
        if any(key in item for key in key_strs):
            fn = path + item
            fn2= fn.replace(vn,vn2)
            fn3= fn.replace(vn,vn3)
            fn4= fn.replace(vn,vn4)

            print('reading netcdf '+fn+' -'+vn)
            print(fn3)

            ds = xr.open_dataset(fn, decode_times=False)
            var0 = ds[vn].values[0,1:,:,:] 
            lon0 = ds['lon'].values
            lat0 = ds['lat'].values
            depth = ds['depth'].values[1:]
            
            ds2 = xr.open_dataset(fn2, decode_times=False)
            var2 = ds2[vn2].values[0,1:,:,:] 

            lon_loc = [0, 0];lat_loc = [0, 0]
            lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
            lon_loc[1] = np.argmin(abs(lon0-lon_range[1])) + 1 
            lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
            lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) + 1  
            lon = lon0[lon_loc[0]:lon_loc[1]]
            lat = lat0[lat_loc[0]:lat_loc[1]]
            nlon = len(lon) 
            nlat = len(lat) 
            if number_day == 0 :
                deudx_e = np.zeros((2400,ndep,nlat))
                deudx_w = np.zeros((2400,ndep,nlat))
                devdy_s = np.zeros((2400,ndep,nlon))
                devdy_n = np.zeros((2400,ndep,nlon))

            R = 6371000 
            dlon = np.deg2rad(np.diff(lon)[0]) 
            dlat = -np.deg2rad(np.diff(lat)[0])

            lon2d, lat2d = np.meshgrid(lon, lat)

            # Calculate grid cell areas
            dy = R * dlat
            dx2 = abs(R * np.cos(np.deg2rad(lat2d)) * dlon)
            dx = np.repeat(dx2[np.newaxis, :, :], ndep, axis=0)
            areas = dx * dy  

            # varx = np.squeeze(var3[lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
            # vary = np.squeeze(var4[lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
            varu = np.squeeze(var0[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            varv = np.squeeze(var2[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            varum= np.squeeze(um[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            varvm= np.squeeze(vm[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            # ind = np.squeeze(ind0[lev-1, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])

            vare = 1025*((varu-varum)*(varu) + (varv-varvm)*(varv)) 
            areas = areas * varu/varu  

            # vareu = -vare * (varu-varum) 
            # varev = -vare * (varv-varvm)
            if number_day == 0:
                varpu = np.zeros_like(varu)
                varpv = np.zeros_like(varu)
            varpu += -vare * (varu-varum) 
            varpv += -vare * (varv-varvm)
            print(np.nanmean(varpu,axis=(1,2)))
            print(np.nanmean(varpu))

            # deudx_e[number_day, :, :] = (vareu[:,:,-1] - vareu[:,:,-2]) / ((dx[:,:,-1]+dx[:,:,-2])/2) 
            # deudx_w[number_day, :, :] = (vareu[:,:, 1] - vareu[:,:, 0]) / ((dx[:,:, 1]+dx[:,:, 0])/2) 
            # devdy_s[number_day, :, :] = (varev[:,-1,:] - varev[:,-2,:]) / dy 
            # devdy_n[number_day, :, :] = (varev[:, 1,:] - varev[:, 0,:]) / dy 
            # deudx_e[number_day, :, :] = (vareu[:,:,-1]) 
            # deudx_w[number_day, :, :] = (vareu[:,:, 0]) 
            # devdy_s[number_day, :, :] = (varev[:,-1,:]) 
            # devdy_n[number_day, :, :] = (varev[:, 0,:]) 
            # print(vareu.shape)

            print([np.nansum(deudx_e[number_day, :]),np.nansum(deudx_w[number_day, :]),
                np.nansum(devdy_s[number_day, :]),np.nansum(devdy_n[number_day, :])])
            number_day = number_day + 1

            ds.close()
            ds2.close()

    varpu = varpu / number_day 
    varpv = varpv / number_day
    print(np.nanmean(varpu))

    def save_to_netcdf(filename, varpu, varpv, lon, lat, depth, time_index=None):

        if time_index is None:
            ds = xr.Dataset(
                {
                    "varpu": (("depth", "lat", "lon"), varpu),
                    "varpv": (("depth", "lat", "lon"), varpv),
                },
                coords={
                    "depth": depth,
                    "lat": lat,
                    "lon": lon,
                },
            )
        else:
            ds = xr.Dataset(
                {
                    "varpu": (("time", "depth", "lat", "lon"), varpu[np.newaxis, ...]),
                    "varpv": (("time", "depth", "lat", "lon"), varpv[np.newaxis, ...]),
                },
                coords={
                    "time": [time_index],
                    "depth": depth,
                    "lat": lat,
                    "lon": lon,
                },
            )

        ds.to_netcdf(filename)
        ds.close()


    fnot = 'data/' + pnstr + '.nc'
    save_to_netcdf(fnot, varpu, varpv, lon, lat, depth)
