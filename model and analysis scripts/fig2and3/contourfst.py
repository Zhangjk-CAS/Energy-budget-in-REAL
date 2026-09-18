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
lon_range=[100, 148];lat_range=[-15,45]
level = 2.5

pnstr = 'sss-eg-woa-win'
key_str1 = '12.1.2-5year-ss.nc'
key_str2 = '12.1.2-ss-woa.nc'
# key_str1 = 'ss-2014-2018.nc'
# key_str2 = 'ss-2014-2018.nc'
# key_str2 = 'oisst.2014-2018.10km'
# key_str2 = 'woa-tt2015-2022.nc'
# key_str2 = '12.1.2-5year-tt.nc'
# tn1 = 'a) GLORYS'
# tn2 = 'b) Model'
# tn3 = 'c) difference'
vn1 = 'ss'
vn2 = 's_an'
# pathstr = ["../ave-TH/eg/", "../exe/boundary-data-glorys-daily/ave/"]
pathstr = ["../ave-TH/eg/", "../ave-TH/woa/"]

number_var = 0
for path in pathstr:
    nowdir = os.path.abspath(path)
    if number_var == 0 :
        key_str = key_str1
        vn = vn1
    elif number_var == 1 :
        key_str = key_str2
        vn = vn2
    for item in sorted(os.listdir(nowdir)):
        # fn = []
        if key_str in item and '.png' not in item:
            fn = path + item
            
            lev = 1

            print('reading netcdf '+fn+' -'+vn)

            ncdata = nc.Dataset(fn)
            var0 = ncdata.variables[vn][:]
            # lon0 = ncdata.variables['lon'][:]
            # lat0 = ncdata.variables['lat'][:]

            ncdata_ind = nc.Dataset('ind_final_eas.nc')
            # ind0 = ncdata_ind.variables['ind'][:]
            lon0 = ncdata.variables['lon'][:]
            lat0 = ncdata.variables['lat'][:]
            if number_var == 0:
                dep0 = abs(ncdata.variables['lev1'][:])
            elif number_var == 1:
                dep0 = abs(ncdata.variables['depth'][:])

            lon_loc = [0, 0];lat_loc = [0, 0]
            lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
            lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
            lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
            lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
            lon = lon0[lon_loc[0]:lon_loc[1]]
            lat = lat0[lat_loc[0]:lat_loc[1]]
            dep_loc =np.argmin(abs(dep0-level)) 
            print(dep_loc)

            # if vn == 'ss':
            #     var0 = var0*1000+35
            # elif vn == 'vv':
            #     var0 = -var0
            
            # ind = np.squeeze(ind0[lev-1,lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) + 0.0
            # ind[ind==0] = np.nan 
            if number_var == 0:
                var = np.squeeze(var0[0, dep_loc, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]).filled(np.nan)
            else:
                var = np.squeeze(var0[0, dep_loc, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]).filled(np.nan)
                
            var = var[np.newaxis, :, :]

            lon_tick = np.arange(115, lon[-1], 15)
            lat_tick = np.arange(40, -15, -10);lat_tick = lat_tick[::-1]

            if vn1 == 'tt':
                bar_limit0 = 10
                bar_limit1 = 30
            elif vn1 == 'ss':
                bar_limit0 = 31
                bar_limit1 = 35

            nlevel = np.linspace(bar_limit0, bar_limit1, 16)

            if number_var == 0:
                var_fig = var
            else:
                var_fig = np.concatenate([var_fig, var], axis=0)

            number_var = number_var + 1

##################################### R-LICOM ######################################################################

def corr2d_ignore_nan(a, b):

    # 将数组展平
    a_flat = a.flatten()
    b_flat = b.flatten()
    
    # 找到非NaN的索引
    mask = ~(np.isnan(a_flat) | np.isnan(b_flat))
    a_valid = a_flat[mask]
    b_valid = b_flat[mask]
    
    # 计算相关系数
    if len(a_valid) < 2:
        return np.nan
    
    correlation = np.corrcoef(a_valid, b_valid)[0, 1]
    return correlation

corr = corr2d_ignore_nan(var_fig[0,:,:], var_fig[1,:,:])
# print(f"\n相关系数: {corr:.4f}")

mse = np.nanmean((var_fig[0,:,:].flatten() - var_fig[1,:,:].flatten()) ** 2)
rmse = np.sqrt(mse)
print(corr)
print(rmse)

print('contourf data1')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='gray')
ax.add_feature(cfeature.LAND, facecolor='gray')

cf = ax.contourf(lon,lat,var_fig[0,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap='RdBu_r', extend='both')

ax.set_xticks(lon_tick)
ax.set_yticks(lat_tick)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.tick_params(axis='both', labelsize=15)

cbar = plt.colorbar(cf, shrink=0.8)
cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))

# plt.title(tn1, loc='left', pad=2, fontsize=18)

if vn1 == 'ss': 
    cbar.set_label('psu')
elif vn1 == 'tt':
    cbar.set_label('℃')
elif vn1 == 'uu' or vn1 == 'vv':
    cbar.set_label('m/s')

pn = 'figures/obs/'+ pnstr + '-1.png'
plt.savefig(pn, dpi=600)
plt.close()

##################################### mercator ######################################################################

print('contourf data2')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='gray')
ax.add_feature(cfeature.LAND, facecolor='gray')

cf = ax.contourf(lon,lat,var_fig[1,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap='RdBu_r', extend='both')

ax.set_xticks(lon_tick)
ax.set_yticks(lat_tick)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.tick_params(axis='both', labelsize=15)

cbar = plt.colorbar(cf, shrink=0.8)
cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))

# plt.title(tn2, loc='left', pad=2, fontsize=18)

if vn == 'ss': 
    cbar.set_label('psu')
elif vn == 'tt':
    cbar.set_label('℃')
elif vn == 'uu' or vn == 'vv':
    cbar.set_label('m/s')

pn = 'figures/obs/'+ pnstr + '-2.png'
plt.savefig(pn, dpi=600)
plt.close()

##################################### v2503 - mercator ######################################################################
if vn1 == 'tt':
    bar_limit0 = -2
    bar_limit1 = 2
elif vn1 == 'ss':
    bar_limit0 = -1
    bar_limit1 = 1

nlevel = np.linspace(bar_limit0, bar_limit1, 16)

print('contourf 1-2')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='gray')
ax.add_feature(cfeature.LAND, facecolor='gray')

cf = ax.contourf(lon,lat,var_fig[0,:,:]-var_fig[1,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap='RdBu_r', extend='both')

ax.set_xticks(lon_tick)
ax.set_yticks(lat_tick)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.tick_params(axis='both', labelsize=15)

cbar = plt.colorbar(cf, shrink=0.8)
cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))

if vn1 == 'ss': 
    cbar.set_label('psu')
    plt.title(f'RMSE={rmse:.3f}psu', loc='right', pad=2, fontsize=12)
elif vn1 == 'tt':
    cbar.set_label('℃')
    plt.title(f'RMSE={rmse:.2f}℃', loc='right', pad=2, fontsize=12)
elif vn1 == 'uu' or vn == 'vv':
    cbar.set_label(r'$m^2/s^2$')
    plt.title(f'RMSE={rmse:.2f}'+r'$m^2/s^2$', loc='right', pad=2, fontsize=12)

pn = 'figures/obs/'+ pnstr + '-1m2.png'
plt.savefig(pn, dpi=600)
plt.close()
