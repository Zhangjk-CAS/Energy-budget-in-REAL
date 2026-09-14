import numpy as np
import netCDF4 as nc 
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

lon_range=[137];lat_range=[5,34]

vn = 'ss' 
pnstr = '137section-' + vn + '-win-glorys'
key_strs = ['12.1.2-5year-ss.nc'] 
# key_strs = ['woa-ss2015-2022.nc'] 
# key_strs = ['argo-ss-2014-2018-10km.nc']
# pathstr = ["../ave-TH/eg/"] 
pathstr = ["../exe/boundary-data-glorys-daily/ave/"] 
# pathstr = ["../../observation/woa/"] 
# pathstr = ["./"] 
# fig, axs = plt.subplots(3, 1, figsize=(6, 7), sharex=True)
number_fig = 0 
if vn in ['ss', 's_an', 'PSAL']:
    bar_limit0 = 33.6
    bar_limit1 = 35 
    bar_limitdiff0 = -0.5
    bar_limitdiff1 =  0.5
elif vn in ['tt', 't_an', 'TEMP']:
    bar_limit0 = 0
    bar_limit1 = 30
    bar_limitdiff0 = -1
    bar_limitdiff1 =  1
elif vn == 'uu':
    bar_limit0 = -0.5
    bar_limit1 =  0.5
    bar_limitdiff0 = -0.5
    bar_limitdiff1 =  0.5
elif vn == 'vv':
    bar_limit0 = -0.5
    bar_limit1 =  0.5
    bar_limitdiff0 = -0.5
    bar_limitdiff1 =  0.5

for key_str in key_strs:
    var_z = np.zeros((56,2))
    var_compare = np.zeros((3, 56, 290))
    number_var = 0
    for path in pathstr:
        nowdir = os.path.abspath(path)
        for item in sorted(os.listdir(nowdir)):
            # fn = []
            if key_str in item and '.png' not in item:
                fn = path + item

                print('reading netcdf '+fn+' -'+vn)

                dz = np.zeros((55))
                ncdata = nc.Dataset(fn)
                if number_var == 1:
                    var0 = ncdata.variables[vn][:,1:,:,:].filled(np.nan)
                    lon0 = ncdata.variables['lon'][:]
                    lat0 = ncdata.variables['lat'][:]
                    depth= ncdata.variables['depth'][1:]
                else:
                    var0 = ncdata.variables[vn][:].filled(np.nan)
                    lon0 = ncdata.variables['lon'][:]
                    lat0 = ncdata.variables['lat'][:]
                    depth= abs(ncdata.variables['depth'][:])

                lon_loc = [0, 0];lat_loc = [0, 0]
                lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
                # lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
                lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) 
                lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
                lon = lon0[lon_loc[0]]
                lat = lat0[lat_loc[0]:lat_loc[1]]

                var = (np.squeeze(var0[0,:, lat_loc[0]:lat_loc[1], lon_loc[0]])) 
                print(var.shape)

                # dz = dz.reshape(55, 1)
                var_compare[number_var, :, :] = var

                number_var = number_var + 1

    var_compare[2, :, :] = var_compare[0, :, :] - var_compare[1, :, :]

nlevel = np.linspace(bar_limit0, bar_limit1, 16) 
lat_tick = np.arange(lat[-1], lat[0], 5);lat_tick = lat_tick[::-1] 

fig, ax = plt.subplots(1, 1, figsize=(6, 1.9))

cf = ax.contourf(lat,depth,var_compare[0, :, :], levels=nlevel, cmap='RdBu_r', extend='both')
ax.pcolormesh(lat, depth, np.isnan(var_compare[0,:,:]), cmap=mpl.colors.ListedColormap(['none','gray']))

ax.set_ylim(2.5,depth[-1])
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_xticks(lat_tick)

cbar = plt.colorbar(cf, shrink=0.8)
cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
# cbar.set_label(r"℃")

xticks = ax.get_xticks()
ax.set_xticks(xticks)
ax.set_xticklabels([f"{int(tick)}°N" for tick in xticks])
ax.set_xlim(lat_range)
ax.set_ylabel('depth(m)')

pn = 'figures/'+pnstr+'.png'
plt.savefig(pn, dpi=600)
plt.close()

# # ---------------------------------------------------------------------------------
# fig, ax = plt.subplots(1, 1, figsize=(6, 1.9))

# cf = ax.contourf(lat,depth,var_compare[1, :, :], levels=nlevel, cmap='RdBu_r', extend='both')
# ax.pcolormesh(lat, depth, np.isnan(var_compare[1,:,:]), cmap=mpl.colors.ListedColormap(['none','gray']))

# ax.set_ylim(2.5,depth[-1])
# ax.invert_yaxis()
# ax.set_yscale('log')
# ax.set_xticks(lat_tick)

# cbar = plt.colorbar(cf, shrink=0.8)
# cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))

# xticks = ax.get_xticks()
# ax.set_xticks(xticks)
# ax.set_xticklabels([f"{int(tick)}°N" for tick in xticks])
# ax.set_xlim(lat_range)

# pn = 'figures/'+pnstr+'-glorys.png'
# plt.savefig(pn, dpi=600)
# plt.close()