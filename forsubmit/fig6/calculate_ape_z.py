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


inputstr = 'nwp'
npath = 0

if inputstr == 'nwp':
    lon_range=[100, 148];lat_range=[-15,45]
elif inputstr == 'scs':
    lon_range=[105, 125];lat_range=[0, 25]
elif inputstr == 'ke':
    lon_range=[130, 140];lat_range=[25, 35]
elif inputstr == 'stcc':
    lon_range=[125, 145];lat_range=[18, 24]
# lon_range=[95, 150];lat_range=[-16,47]
# lon_range=[105, 140];lat_range=[0,40] 
# lon_range=[105, 125];lat_range=[0, 25]
# lon_range=[130, 140];lat_range=[25, 35]
# lon_range=[125, 145];lat_range=[18, 24]

paths = ["../exe/boundary-data-glorys-daily/ave/"]#"../output-era-glorys/fort.22/TH/ave/", "../output-jra-glorys/fort.22/TH/ave/", "../output-jra-licom3/fort.22/TH/ave/"]
pnstrs= ["-glorys-"]#'-era-glorys-', '-jra-glorys-', '-jra-licom3-']
path = paths[npath] 
pnstr = 'ape' + pnstrs[npath] + inputstr 

print(pnstr)

rho = 1025

key_strs = ['rho-2014-2018.nc']
# key_strs = ['uu-2013','uu-2014','uu-2015','uu-2016','uu-2017','uu-2018']
vn = 'rho'
vn2= 'vv'
vn3= 'rho'
# pnstr = 'eke-jra-licom3-nwp'
# path = "../output-jra-glorys/"
# path = "../output-FaD-era/"
# path = "../output-FaD-jra/"
# path = "../../mercatordata/mercator_compare/RegionalData/"
nowdir = os.path.abspath(path)
number_day = 0
regionalmean = np.zeros((2400,55))
weighted_averages = np.zeros((55))
meansst = np.zeros((2400,1))
dz = np.zeros((55))

# dsum = xr.open_dataset(path + 'uu-2014-2018.nc', decode_times=False)
# varum = dsum['uu'].values[:,:,:,:]
# dsvm = xr.open_dataset(path + 'vv-2014-2018.nc', decode_times=False)
# varvm = dsvm['vv'].values[:,:,:,:]

for item in sorted(os.listdir(nowdir)):
    fn = []
    if any(key in item for key in key_strs):
        fn = path + item
        # fn2= fn.replace(vn,vn2)

        # fn3= path + 'monthly-5year-uu-' + item[11:13] + '.nc' #'annual-ave/uu-2014-2018.nc'
        # fn4= path + 'monthly-5year-vv-' + item[11:13] + '.nc'

        fn3= path + 'uu-2014-2018.nc'
        fn4= path + 'vv-2014-2018.nc'

        print(fn)
        print('reading netcdf '+fn+' -'+vn)

        # ds3 = xr.open_dataset(fn3, decode_times=False)
        # varum = ds3[vn].values[:,:,:,:]
        # ds4 = xr.open_dataset(fn4, decode_times=False)
        # varvm = ds4[vn2].values[:,:,:,:]

        ds = xr.open_dataset(fn, decode_times=False)
        var0 = ds[vn].values[:,:,:,:] #- varum 
        lon0 = ds['lon'].values
        lat0 = ds['lat'].values
        depth = ds['lev'].values[:]

        # ds2 = xr.open_dataset(fn2, decode_times=False)
        # var2 = ds2[vn2].values[:,:,:,:] #- varvm 

        lon_loc = [0, 0];lat_loc = [0, 0]
        lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
        lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
        lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
        lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
        lon = lon0[lon_loc[0]:lon_loc[1]]
        lat = lat0[lat_loc[0]:lat_loc[1]]
        nlon = len(lon) 
        nlat = len(lat) 
        ndep = len(depth)

        if number_day == 0 :
            eke = np.zeros((len(depth),len(lat),len(lon)))

        g = 9.8 
        rho0 = 1025 
        R = 6371000 
        dlon = np.deg2rad(np.diff(lon)[0]) 
        dlat = -np.deg2rad(np.diff(lat)[0])

        lon2d, lat2d = np.meshgrid(lon, lat)

        # Calculate grid cell areas
        dy = R * dlat
        dx2 = abs(R * np.cos(np.deg2rad(lat2d)) * dlon)
        dx = np.repeat(dx2[np.newaxis, :, :], ndep, axis=0)
        areas = dx * dy  
        dz_mid = (depth[1:] + depth[:-1]) / 2
        dz = np.zeros_like(depth)
        dz[0] = depth[0] / 2
        dz[1:-1] = dz_mid[1:] - dz_mid[:-1]
        dz[-1] = (depth[-1] - dz_mid[-1]) * 2
        print(areas.shape)
        print(var0.shape)
            
        varrho = np.squeeze(var0[0, :, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
        areas = areas * varrho/varrho   

        # get reference N^2 
        rhoref = np.nansum(areas * varrho, axis=(1,2))/np.nansum(areas, axis=(1,2))
        rhoa = varrho - rhoref[:, np.newaxis, np.newaxis] 
        N2 = g / rho0 * np.gradient(np.squeeze(rhoref), axis=0) / dz

        # Mask unstable or excessively weak stratification
        N2_min = 1.0e-8  # s^-2
        invalid_n2 = (~np.isfinite(N2)) | (N2 <= N2_min)

        N2_safe = N2.astype(np.float64, copy=True)
        N2_safe[invalid_n2] = np.nan

        print(
            f"Invalid N2 levels: {invalid_n2.sum()}/{N2.size}; "
            f"valid range: {np.nanmin(N2_safe):.3e}–"
            f"{np.nanmax(N2_safe):.3e} s^-2"
        )
        #----------------------------------------------------------
        # get clean rhoa**2
        valid_rho = (
            np.isfinite(varrho)
            & (np.abs(varrho) < 1.0e10)
        )

        bottom_mask = np.zeros_like(valid_rho, dtype=bool)

        for j in range(nlat):
            for i in range(nlon):

                wet_k = np.flatnonzero(valid_rho[:, j, i])

                if wet_k.size > 20:
                    remove_k = wet_k[-5:]
                    bottom_mask[remove_k, j, i] = True

        rhoa_clean = np.where(bottom_mask, np.nan, rhoa)

        rhoa2 = np.square(rhoa_clean.astype(np.float64))
        #-----------------------------------------------------------

        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            ape = (
                g**2 * rhoa2
                / (2.0 * rho0 * N2_safe[:, np.newaxis, np.newaxis])
            )

        ape[~np.isfinite(ape)] = np.nan

        eke = eke + ape 

        number_day = number_day + 1
        # regionalmean[number_day-1,0] = np.sum(weighted_averages * dz )
        # meansst[number_day-1,0] = weighted_averages[kk+1] * dz[kk]
        print(regionalmean[number_day-1,0])

        ds.close()

ape [abs(varrho)<1e-10] = np.nan 
# fnot1 = 'data/average-ske-'+pnstr+'.txt'
# fnot1 = f"data/average-ske-{pnstr}-{kk}.txt"
# fnot2 = 'data/'+pnstr+'.txt'

# # np.savetxt(fnot1, meansst, fmt="%.6f")
# np.savetxt(fnot2, regionalmean, fmt="%.6f")

fnot3 = 'data/'+pnstr+'.nc'

eke_da = xr.DataArray(
    eke,
    coords={
        "depth": depth,
        "lat": lat,
        "lon": lon
    },
    dims=["depth", "lat", "lon"],
    name="eke"
)

ds = xr.Dataset({"eke": eke_da})
ds.to_netcdf(fnot3)