import numpy as np
import netCDF4 as nc 
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.dates as mdates
from datetime import datetime, timedelta

lon_range=[137];lat_range=[0,35];lev_range=[0,1000]

# pnstr = 'ssh137E-day'
# pnstr = 'ssh137E-force'
pnstr = 'ssh137E-all'

# key_strs = ['z0-2014-2018-era.nc','z0-2014-2018-glorys.nc','z0-2014-2018-AVISO.nc']
key_strs = ['z0-2014-2018-era-glorys.nc','z0-2014-2018-jra-glorys.nc','z0-2014-2018-jra-licom3.nc',
            'z0-2014-2018-glorys.nc','z0-2014-2018-AVISO.nc']
path = "./"
nowdir = os.path.abspath(path)
number_day = 0
regionalmean = np.zeros((1200,5))
weighted_averages = np.zeros((55))
meansst = np.zeros((1200,2))
for key_str in key_strs:
    for item in sorted(os.listdir(nowdir)):
        fn = []
        if key_str in item:
            fn = path + item
            if number_day < 2.5 :
                vn = 'z0'
            elif number_day == 3:
                vn = 'ssh'
            elif number_day == 4:
                vn = 'adt'

            print('reading netcdf '+fn+' -'+vn)

            ncdata = nc.Dataset(fn)
            var0 = ncdata.variables[vn][:].filled(np.nan)
            lon0 = ncdata.variables['lon'][:]
            lat0 = ncdata.variables['lat'][:]

            lon_loc = [0, 0];lat_loc = [0, 0];lev_loc = [0, 0]
            lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
            # lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
            lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
            lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
            # lev_loc[0] = np.argmin(abs(depth0-lev_range[0])) 
            # lev_loc[1] = np.argmin(abs(depth0-lev_range[1])) 
            lon = lon0[lon_loc[0]:lon_loc[1]]
            lat = lat0[lat_loc[0]:lat_loc[1]];lat1=lat
            # depth = depth0[lev_loc[0]:lev_loc[1]]

            if number_day < 2.5 :
                var = np.squeeze(var0[0, 0, lat_loc[0]:lat_loc[1], lon_loc[0]])
            else :
                var = np.squeeze(var0[0, lat_loc[0]:lat_loc[1], lon_loc[0]])

            nlat0 = lat.shape
            nlat = nlat0[0]
            regionalmean[0:nlat,number_day] = var

            number_day = number_day + 1

fig, axs = plt.subplots(1, 1, figsize=(6, 3), sharex=True)
colors = ['blue','red','green','blue','green','gray']
lnstr = ['-', '-', '-', '--', '--']
# axs.plot(lat, regionalmean[0:nlat,0], label='RLICOM')
# axs.plot(lat, regionalmean[0:nlat,1], label='GLORYS')
# axs.plot(lat, regionalmean[0:nlat,2], label='AVISO')
axs.plot(lat, regionalmean[0:nlat,0], linestyle=lnstr[0], label='Reg-EG', color=colors[0])
axs.plot(lat, regionalmean[0:nlat,1], linestyle=lnstr[1], label='Reg-JG', color=colors[1])
axs.plot(lat, regionalmean[0:nlat,2], linestyle=lnstr[2], label='Reg-JL', color=colors[2])
axs.plot(lat, regionalmean[0:nlat,3], linestyle=lnstr[3], label='GLORYS', color=colors[3])
axs.plot(lat, regionalmean[0:nlat,4], linestyle=lnstr[4], label='AVISO', color=colors[4])
# axs.set_title('', loc='left')
axs.set_ylabel('SSH (m)')
axs.set_xlabel('Latitude')
axs.legend(loc='upper left',ncol=2,frameon=False)
xticks = axs.get_xticks()
axs.set_xticks(xticks)
axs.set_xticklabels([f"{int(tick)}°N" for tick in xticks])

plt.xlim(lat[-1]-3, lat[0])
# plt.ylim(0, 2)
plt.ylim(-0.1, 2.1)

pn = 'figures/' + pnstr + '.png'
plt.savefig(pn, dpi=600)
plt.close()
