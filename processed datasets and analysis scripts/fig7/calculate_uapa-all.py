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
vn3= 'rho'

# pnstrs = [f'edgeuapa-eg-all-TH']
# paths = ["../output-era-glorys/fort.22/"]
pnstrs = [f'TH/edgeuapa-eg-all', f'TH/edgeuapa-jg-all', f'TH/edgeuapa-jl-all']
paths = ["../output-era-glorys/fort.22/TH/", "../output-jra-glorys/fort.22/TH/", "../output-jra-licom3/fort.22/TH/"]
for nforce in range(3):
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
    dsrm = xr.open_dataset(path + 'ave/rho-2014-2018.nc', decode_times=False)
    dshm = xr.open_dataset(path + 'ave/z0-2014-2018.nc', decode_times=False)
    # xm = dsxm['tx'].values[0,0,:,:]
    # ym = dsym['ty'].values[0,0,:,:]
    um = dsum['uu'].values[0,:,:,:]
    vm = dsvm['vv'].values[0,:,:,:]
    rhom = dsrm['rho'].values[0,:,:,:]
    sshm = dshm['ssh'].values[0,:,:,:]
    depth = dsrm['lev1'].values[:]
    ndep = len(depth)
    dep3d = np.broadcast_to(depth.reshape(-1, 1, 1), rhom.shape)
    dz_mid = (dep3d[1:, :, :] + dep3d[:-1, :, :]) / 2
    dz = np.zeros_like(dep3d)
    dz[0, :, :] = depth[0] / 2
    dz[1:-1, :, :] = dz_mid[1:, :, :] - dz_mid[:-1, :, :]
    dz[-1, :, :] = (depth[-1] - dz_mid[-1, :, :]) * 2
    dzm = dz.copy()
    dzm[0,:,:] = dzm[0,:,:] + sshm 

    g = 9.8 
    pm= np.cumsum(np.squeeze(rhom) * g * dzm ,axis=0) + sshm * g * 1025

    for item in sorted(os.listdir(nowdir)):
        fn = []
        if any(key in item for key in key_strs):
            fn = path + item
            fn2= fn.replace(vn,vn2)
            fn3= fn.replace(vn,'rho')
            fn4= fn.replace(vn,'z0')

            print('reading netcdf '+fn+' -'+vn)
            print(fn3)

            ds = xr.open_dataset(fn, decode_times=False)
            var0 = ds[vn].values[0,:,:,:] - um 
            lon0 = ds['lon'].values
            lat0 = ds['lat'].values
            
            ds2 = xr.open_dataset(fn2, decode_times=False)
            var2 = ds2[vn2].values[0,:,:,:] - vm 

            ds4 = xr.open_dataset(fn4, decode_times=False)
            var4 = ds4['ssh'].values[0,:,:,:] 

            ds3 = xr.open_dataset(fn3, decode_times=False)
            var3 = ds3[vn3].values[0,:,:,:] 
            dzday = dz.copy()
            # dzday[0,:,:] = dzday[0,:,:] + var4
            p = np.cumsum(np.squeeze(var3) * g * dzday ,axis=0) + var4 * g * 1025
            pa = p - pm 

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
            varua = np.squeeze(var0[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            varva = np.squeeze(var2[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            varpa = np.squeeze(pa[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 

            areas = areas * varua/varua  

            if number_day == 0:
                varpu = np.zeros_like(varua)
                varpv = np.zeros_like(varua)
            varpu += varpa * varua
            varpv += varpa * varva 

            print(np.nanmean(varpu))
            print(np.nanmean(varpa))
            print(np.nanmean(varua))

            deudx_e[number_day, :, :] = (varpu[:,:,-8]) 
            print([np.nansum(deudx_e[number_day, 0, :]),np.nansum(deudx_w[number_day, 0, :]),
                np.nansum(devdy_s[number_day, 0, :]),np.nansum(devdy_n[number_day, 0, :])])
            number_day = number_day + 1

            ds.close()
            ds2.close()
            ds3.close()

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
