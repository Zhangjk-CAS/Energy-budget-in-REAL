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
lon_range=[100, 148];lat_range=[-15,45]
# lon_range=[105, 125];lat_range=[0, 25]
# lon_range=[105, 140];lat_range=[0,40]

# key_strs = ['uu2014','uu2015','uu2016','uu2017','uu2018']
key_strs = ['ww-2014','ww-2015','ww-2016','ww-2017','ww-2018']
vn = 'ww'
vn2= 'rho'
# vn3= 'rho'

pnstrs = [f'epe2eke-eg-dailyBC']
paths = ["../output-era2/fort.22/"]
for nforce in range(1):
    pnstr = pnstrs[nforce]
    path = paths[nforce]

    nowdir = os.path.abspath(path)
    number_day = 0 
    regionalmean = np.zeros((2400,1))
    weighted_averages = np.zeros((55))
    meansst = np.zeros((2400,1))
    dz = np.zeros((55))

    # dsind = xr.open_dataset('ind_final_eas_uv.nc', decode_times=False)
    # ind0 = dsind['ind'].values[:,:,:] + 0.0 
    # ind0[ind0<0.5] = np.nan 

    # dsdep = xr.open_dataset('./optimized_ocean_bathymetry_10_RL.nc', decode_times=False)
    # oceandepth0 = dsdep['depth'].values

    # dsxm = xr.open_dataset(path + 'tx-2014-2018.nc', decode_times=False)
    # dsym = xr.open_dataset(path + 'ty-2014-2018.nc', decode_times=False)
    # dsum = xr.open_dataset(path + 'ave/uu-2014-2018.nc', decode_times=False)
    # dsvm = xr.open_dataset(path + 'ave/vv-2014-2018.nc', decode_times=False)
    dsrm = xr.open_dataset(path + 'ave/rho-2014-2018.nc', decode_times=False)
    dswm = xr.open_dataset(path + 'ave/ww-2014-2018.nc', decode_times=False)
    # xm = dsxm['tx'].values[0,0,:,:]
    # ym = dsym['ty'].values[0,0,:,:]
    # um = dsum['uu'].values[0,:,:,:]
    # vm = dsvm['vv'].values[0,:,:,:]
    rhom = dsrm['rho'].values[0,:,:,:]
    wm = dswm['ww'].values[0,:,:,:]
    depth = dsrm['lev1'].values[:]
    ndep = len(depth)
    dep3d = np.broadcast_to(depth.reshape(-1, 1, 1), rhom.shape)
    dz_mid = (dep3d[1:, :, :] + dep3d[:-1, :, :]) / 2
    dz = np.zeros_like(dep3d)
    dz[0, :, :] = depth[0] / 2
    dz[1:-1, :, :] = dz_mid[1:, :, :] - dz_mid[:-1, :, :]
    dz[-1, :, :] = (depth[-1] - dz_mid[-1, :, :]) * 2
    dzm = dz.copy()
    dzm[0,:,:] = dzm[0,:,:] 

    g = 9.8 
    pm= np.cumsum(np.squeeze(rhom) * g * dzm ,axis=0)

    for item in sorted(os.listdir(nowdir)):
        fn = []
        if any(key in item for key in key_strs):
            fn = path + item
            fn2= fn.replace(vn,'rho')
            # fn3= fn.replace(vn,vn3)
            # fn4= fn.replace(vn,'ww')

            print('reading netcdf '+fn+' -'+vn)
            print(fn2)

            ds = xr.open_dataset(fn, decode_times=False)
            var0 = ds[vn].values[0,:,:,:] - wm 
            lon0 = ds['lon'].values
            lat0 = ds['lat'].values
            
            ds2 = xr.open_dataset(fn2, decode_times=False)
            var2 = ds2[vn2].values[0,:,:,:] - rhom  

            lon_loc = [0, 0];lat_loc = [0, 0]
            lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
            lon_loc[1] = np.argmin(abs(lon0-lon_range[1])) + 1 
            lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
            lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) + 1 
            lon = lon0[lon_loc[0]:lon_loc[1]]
            lat = lat0[lat_loc[0]:lat_loc[1]]
            nlon = len(lon) 
            nlat = len(lat) 

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
            varw = np.squeeze(var0[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            varrho = np.squeeze(var2[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 

            areas[np.isnan(varrho)] = np.nan   

            if number_day == 0:
                p2e = - g * varrho * varw 
            else:
                p2e = p2e - g * varrho * varw 

            print(p2e.shape)

            print([np.nansum(p2e)])
            number_day = number_day + 1

            ds.close()
            ds2.close()

    p2e = p2e/number_day 

    # fnote = 'data/'+pnstr+'-e.txt'
    # np.savetxt(fnote, meansst, fmt="%.6f")
    # fnotw = 'data/'+pnstr+'-w.txt'
    # np.savetxt(fnote, meansst, fmt="%.6f")
    # fnots = 'data/'+pnstr+'-s.txt'
    # np.savetxt(fnote, meansst, fmt="%.6f")
    # fnotn = 'data/'+pnstr+'-n.txt'
    # np.savetxt(fnote, meansst, fmt="%.6f")

    def save_to_netcdf(p2e, output_filename):

        depth_coords = depth 
        lat_coords = lat 
        lon_coords = lon 
        
        ds = xr.Dataset(
            {
                'p2e': (['depth', 'lat', 'lon'], p2e)
            },
            coords={
                'depth': depth_coords,
                'lat': lat_coords,
                'lon': lon_coords
            }
        )
        
        ds.to_netcdf(output_filename)
        print(f"Data saved to {output_filename}")
        
        return ds

    fnot = 'data/' + pnstr + '.nc'
    save_to_netcdf(p2e, fnot)