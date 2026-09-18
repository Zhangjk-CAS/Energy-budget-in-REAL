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

# paths = ["../output-jra-glorys/fort.22/ave/", "../output-FaD-era/fort.22/ave/", "../output-FaD-jra/fort.22/ave/"]
paths = ["../output-era2/fort.22/ave/"]
pnstrs= ['-eg-dailyBC']
path = paths[npath] 
pnstr = 'mke' + pnstrs[npath] + inputstr 

print(pnstr)

rho = 1025

key_strs = ['uu-2014-2018']
vn = 'uu'
vn2= 'vv'
nowdir = os.path.abspath(path)
number_day = 0
regionalmean = np.zeros((2400,55))
weighted_averages = np.zeros((55))
meansst = np.zeros((2400,1))
dz = np.zeros((55))

for item in sorted(os.listdir(nowdir)):
    fn = []
    if any(key in item for key in key_strs):
        fn = path + item
        fn2= fn.replace(vn,vn2)

        print(fn)
        print('reading netcdf '+fn+' -'+vn)

        ds = xr.open_dataset(fn, decode_times=False)
        var0 = ds[vn].values[:,:,:,:] #- varum 
        lon0 = ds['lon'].values
        lat0 = ds['lat'].values
        depth = ds['lev1'].values[:]

        ds2 = xr.open_dataset(fn2, decode_times=False)
        var2 = ds2[vn2].values[:,:,:,:] #- varvm 

        lon_loc = [0, 0];lat_loc = [0, 0]
        lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
        lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
        lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
        lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
        lon = lon0[lon_loc[0]:lon_loc[1]]
        lat = lat0[lat_loc[0]:lat_loc[1]]

        if number_day == 0 :
            eke = np.zeros((len(depth),len(lat),len(lon)))

        dz[0] = depth[0] - 0
        dz[1:] = np.diff(depth)
        dz = abs(dz)
        R = 6371000 
        dlon = np.deg2rad(np.diff(lon)[0]) 
        dlat = -np.deg2rad(np.diff(lat)[0])

        lon2d, lat2d = np.meshgrid(lon, lat)

        # Calculate grid cell areas
        dy = R * dlat
        dx = abs(R * np.cos(np.deg2rad(lat2d)) * dlon)
        areas = dx * dy  

        if number_day == 0 :
            ind = np.squeeze(var0[0, :, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])

        for lev in range(55):
            
            varu = np.squeeze(var0[0,lev, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
            varv = np.squeeze(var2[0,lev, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
            var = 0.5 * (varu*varu + varv*varv) * rho
            areas = areas * var / var
            if number_day == 0 :
                ind[lev, :, :] = var

            weighted_sum = np.nansum(var * areas)
            total_area = np.nansum(areas[~np.isnan(var)])
            weighted_averages[lev] = weighted_sum / total_area 

            eke[lev,:,:] = eke[lev,:,:] + var
            
            print([weighted_averages[lev], total_area])

        regionalmean[number_day,:] = weighted_averages 
        number_day = number_day + 1
        # regionalmean[number_day-1,0] = np.sum(weighted_averages * dz )
        # meansst[number_day-1,0] = weighted_averages[kk+1] * dz[kk]
        print(regionalmean[number_day-1,0])

        ds.close()
        ds2.close()


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