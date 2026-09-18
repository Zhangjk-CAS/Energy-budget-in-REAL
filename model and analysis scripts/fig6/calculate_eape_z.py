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

paths = ["../output-era2/fort.22/"]
pnstrs= ['-licom3Global']
path = paths[npath] 
pnstr = 'eape' + pnstrs[npath] + inputstr 

print(pnstr)

key_strs = ['uu-2014','uu-2015','uu-2016','uu-2017','uu-2018']
# key_strs = ['rho-2014','rho-2015','rho-2016','rho-2017','rho-2018']
vn = 'rho'
vn2= 'uu'

nowdir = os.path.abspath(path)
number_day = 0
regionalmean = np.zeros((2400,55))
weighted_averages = np.zeros((55))
meansst = np.zeros((2400,1))
dz = np.zeros((55))

dsrm = xr.open_dataset(path + 'ave/rho-2014-2018.nc', decode_times=False)
rhom = dsrm['rho'].values[0,:,:,:]

for item in sorted(os.listdir(nowdir)):
    fn = []
    if any(key in item for key in key_strs):
        fnu= path + item
        fn = fnu.replace(vn2,'rho')

        print(fn)
        print('reading netcdf '+fn+' -'+vn)

        ds = xr.open_dataset(fn, decode_times=False)
        var0 = ds[vn].values[:,:,:,:] - rhom  
        lon0 = ds['lon'].values
        lat0 = ds['lat'].values
        depth = ds['lev1'].values[:]

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
        rhomin = np.squeeze(rhom[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
        areas = areas * varrho/varrho  
        rhomean = np.nansum(areas * varrho, axis=(1,2))/np.nansum(areas, axis=(1,2))
        print(np.nanmean(varrho**2, axis=(1,2)))

        # get reference N^2 
        rhoref = np.nansum(areas * rhomin, axis=(1,2))/np.nansum(areas, axis=(1,2))
        N2 = g / rho0 * np.gradient(np.squeeze(rhoref), axis=0) / dz 
        eape = g**2 * varrho**2 / (2 * rho0 * N2[:, np.newaxis, np.newaxis])

        eke = eke + eape 

        number_day = number_day + 1

        print(np.nanmean(eape))

        ds.close()

fnot3 = 'data/'+pnstr+'.nc'
print('-------------')
print(np.nanmean(eke))
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
