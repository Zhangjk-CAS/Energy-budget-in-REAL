import numpy as np
import netCDF4 as nc 
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

tim_range=[201401,201809];lat_range=[5, 34];dep_range=[0, 2000]

monstr = '201401'
# vn = 'temperature'
vn = 'salinity'
key_str = '137section_observation.nc'
# pnstr = '137section-' + vn + '-' + monstr 
pnstr = '137section-' + vn + '-win' 
path = './' 

nowdir = os.path.abspath(path)
for item in os.listdir(nowdir):
   fn = []
   if key_str in item and '.png' not in item:
        fn = path + item

        print('reading netcdf '+fn+' -'+vn)

        ncdata = nc.Dataset(fn)
        var0 = ncdata.variables[vn][:]
        # lon0 = ncdata.variables['lon'][:]
        lat0 = ncdata.variables['latitude'][:]
        depth0=ncdata.variables['depth_axis'][:]
        time0= ncdata.variables['time'][:]

        lon_loc = [0, 0];lat_loc = [0, 0];dep_loc = [0, 0];tim_loc = [0, 0]
        # lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
        # lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
        lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
        lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
        dep_loc[0] = np.argmin(abs(depth0-dep_range[0])) 
        dep_loc[1] = np.argmin(abs(depth0-dep_range[1])) 
        tim_loc[0] = np.argmin(abs(time0-tim_range[0])) 
        tim_loc[1] = np.argmin(abs(time0-tim_range[1])) + 1 
        # lon = lon0[lon_loc[0]:lon_loc[1]]
        lat = lat0[lat_loc[0]:lat_loc[1]+1]
        depth = depth0[dep_loc[0]:dep_loc[1]+1]
        time = time0[tim_loc[0]:tim_loc[1]+1]
        print(time[0])

        var = np.nanmean(var0[[2, 4, 6, 8, 10],lat_loc[0]:lat_loc[1]+1, dep_loc[0]:dep_loc[1]+1], axis=0)
        # var = np.nanmean(var0[[3, 5, 7, 9, 11],lat_loc[0]:lat_loc[1]+1, dep_loc[0]:dep_loc[1]+1], axis=0)

        lat_tick = np.arange(lat[-1], lat[0], 5);lat_tick = lat_tick[::-1]

        if vn == 'depth':
            bar_limit0 = 0
            bar_limit1 = 2000
        elif vn == 'temperature':
            bar_limit0 = 0
            bar_limit1 = 30
        elif vn == 'salinity':
            bar_limit0 = 33.6 
            bar_limit1 = 35

        nlevel = np.linspace(bar_limit0, bar_limit1, 16)

        print('contourf '+fn)

        fig, ax = plt.subplots(1, 1, figsize=(6, 1.9))

        cf = ax.contourf(lat,depth,var.T, levels=nlevel, cmap='RdBu_r', extend='both')

        ax.set_ylim(2.5,depth[-1])
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_xticks(lat_tick)

        cbar = plt.colorbar(cf, shrink=0.8)
        cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
        cbar.set_label(r"psu")

        xticks = ax.get_xticks()
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{int(tick)}°N" for tick in xticks])
        ax.set_xlim(lat_range)
        ax.set_ylabel('depth(m)')


        pn = 'figures/'+pnstr+'.png'
        plt.savefig(pn, dpi=600)
        plt.close()

