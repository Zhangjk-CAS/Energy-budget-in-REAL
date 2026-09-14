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

# lon_range=[95, 150];lat_range=[-16,47]
lon_range=[101, 147];lat_range=[-14,44]
# lon_range=[100, 126];lat_range=[0,25]

# pnstr = 'windwork-mke-jl'
# key_str1 = 'windwork-mke-jra-licom3-2.nc' # str1 for ej and jl ; str2 for jg 
# key_str2 = 'windwork-mke-jra-glorys-2.nc'
pnstr = 'sla-eg-aviso'
key_str1 = 'slastd-eg.nc' # str1 for ej and jl ; str2 for jg 
key_str2 = 'slastd-aviso-10km.nc'
vn = 'STDSLA'
pathstr = ["./data/", "./data/"]
tn1 = 'a) Model'
tn2 = 'b) GLORYS'
tn3 = 'c) difference'
number_var = 0
for path in pathstr:
    nowdir = os.path.abspath(path)
    if number_var == 0 :
        key_str = key_str1
    elif number_var == 1 :
        key_str = key_str2
    for item in sorted(os.listdir(nowdir)):
        # fn = []
        if key_str in item and '.png' not in item:
            fn = path + item
            
            lev = 1

            print('reading netcdf '+fn+' -'+vn)

            ncdata = nc.Dataset(fn)
            var0 = ncdata.variables[vn][:]
            lon0 = ncdata.variables['lon'][:]
            lat0 = ncdata.variables['lat'][:]

            lon_loc = [0, 0];lat_loc = [0, 0]
            lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
            lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
            lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
            lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
            lon = lon0[lon_loc[0]:lon_loc[1]]
            lat = lat0[lat_loc[0]:lat_loc[1]]
            nlon = len(lon);nlat = len(lat)

            R = 6371000 
            dlon = np.deg2rad(np.diff(lon)[0]) 
            dlat =-np.deg2rad(np.diff(lat)[0])
            lon2d, lat2d = np.meshgrid(lon, lat)
            dy = R * dlat
            dx = abs(R * np.cos(np.deg2rad(lat2d)) * dlon)
            areas = dx * dy 
            # if vn == 'ss':
            #     var0 = var0*1000+35
            # elif vn == 'vv':
            #     var0 = -var0

            var = np.squeeze(var0[lev-1,lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
            var[var.mask] = np.nan

            areas[np.isnan(var)] = np.nan 
            e_total = np.nansum(var*areas)
            print(f'The total wind work is {e_total/1e9}')

            if np.sum(var) < 1e-20:
                var = np.squeeze(var0[lev, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])

            var = var[np.newaxis, :, :]
            var[var==0] = np.nan 

            lon_tick = np.arange(115, lon[-1], 15)
            lat_tick = np.arange(40, -15, -10);lat_tick = lat_tick[::-1]
            # lon_tick = np.arange(105, lon[-1], 10)
            # lat_tick = np.arange(25, 1, -5);lat_tick = lat_tick[::-1]

            bar_limit0 = 0
            bar_limit1 = 0.6

            nlevel = np.linspace(bar_limit0, bar_limit1, 64)

            if number_var == 0:
                var_fig = var
            else:
                var_fig = np.concatenate([var_fig, var], axis=0)

            number_var = number_var + 1

##################################### R-LICOM ######################################################################

mse = np.nanmean((var_fig[0,:,:].flatten() - var_fig[1,:,:].flatten()) ** 2)
rmse = np.sqrt(mse)
print(rmse)

print('contourf nc1')

# fig = plt.figure()
fig = plt.figure(figsize=(4, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='gray')
ax.add_feature(cfeature.LAND, facecolor='gray')

cf = ax.contourf(lon,lat,var_fig[0,:,:]-var_fig[1,:,:]+var_fig[1,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap='Spectral_r', extend='both')

ax.set_xticks(lon_tick)
ax.set_yticks(lat_tick)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
# ax.tick_params(axis='both', labelsize=15)

cbar = plt.colorbar(cf, orientation="vertical", shrink=0.6)
cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
cbar.set_label('m')

# plt.title(tn1, loc='left', pad=2, fontsize=18)
plt.tight_layout()
pn = 'figures/'+ pnstr + '-1.png'
plt.savefig(pn, dpi=600)
plt.close()

##################################### mercator ######################################################################

print('contourf nc2')

# fig = plt.figure()
fig = plt.figure(figsize=(4, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='gray')
ax.add_feature(cfeature.LAND, facecolor='gray')

cf = ax.contourf(lon,lat,var_fig[1,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap='Spectral_r', extend='both')

ax.set_xticks(lon_tick)
ax.set_yticks(lat_tick)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
# ax.tick_params(axis='both', labelsize=15)

cbar = plt.colorbar(cf, orientation="vertical", shrink=0.6)
cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
cbar.set_label('m')

# plt.title(tn2, loc='left', pad=2, fontsize=18)
plt.tight_layout()
pn = 'figures/'+ pnstr + '-2.png'
plt.savefig(pn, dpi=600)
plt.close()

##################################### v2503 - mercator ######################################################################
bar_limit0 =-0.4
bar_limit1 = 0.4

nlevel = np.linspace(bar_limit0, bar_limit1, 64)

print('contourf 1 - 2')

# fig = plt.figure()
fig = plt.figure(figsize=(4, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='gray')
ax.add_feature(cfeature.LAND, facecolor='gray')

cf = ax.contourf(lon,lat,var_fig[0,:,:]-var_fig[1,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap='RdBu_r', extend='both')

ax.set_xticks(lon_tick)
ax.set_yticks(lat_tick)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
# ax.tick_params(axis='both', labelsize=15)

cbar = plt.colorbar(cf, orientation="vertical", shrink=0.6)
cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
if vn == 'ss': 
    cbar.set_label('psu')
    plt.title(f'RMSE={rmse:.3f}psu', loc='right', pad=2, fontsize=12)
elif vn == 'tt':
    cbar.set_label('℃')
    plt.title(f'RMSE={rmse:.2f}℃', loc='right', pad=2, fontsize=12)
elif vn == 'STDSLA' or vn == 'vv':
    # cbar.set_label(r'$m^2/s^2$')
    # plt.title(f'RMSE={rmse:.2f}'+r'$m^2/s^2$', loc='right', pad=2, fontsize=12)
    cbar.set_label('m')
    plt.title(f'RMSE={rmse:.3f}m', loc='right', pad=2, fontsize=12)
# plt.title(tn3, loc='left', pad=2, fontsize=18)
plt.tight_layout()
pn = 'figures/'+ pnstr + '-1m2.png'
plt.savefig(pn, dpi=600)
plt.close()
