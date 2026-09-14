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
import cmaps 

# lon_range=[95, 150];lat_range=[-16,47]
lon_range=[100, 148];lat_range=[-15,45];dep_range=[0,6000]
# lon_range=[100, 126];lat_range=[0,25]

pnstr = 'ape'
vn = 'eke'
path = './data/'

number_var = 0
for item in ['ape-era-glorys-nwp.nc', 'ape-jra-glorys-nwp.nc', 'ape-jra-licom3-nwp.nc', 'ape-glorys-nwp.nc', 'ape-licom3Globalnwp.nc']:

    fn = path + item

    print('reading netcdf '+fn+' -'+vn)

    ncdata = nc.Dataset(fn)
    if number_var == 3:
        var0 = ncdata.variables[vn][1:,:,:].filled(np.nan)
        dep0 = abs(ncdata.variables['depth'][1:])
    else :
        var0 = ncdata.variables[vn][:].filled(np.nan)
        dep0 = abs(ncdata.variables['depth'][:])
    lon0 = ncdata.variables['lon'][:]
    lat0 = ncdata.variables['lat'][:]
    var0 = np.array(var0)
    print(np.nanmean(var0))

    lon_loc = [0, 0];lat_loc = [0, 0];dep_loc = [0, 0]
    lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
    lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
    lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
    lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
    dep_loc[0] = np.argmin(abs(dep0-dep_range[0]))
    dep_loc[1] = np.argmin(abs(dep0-dep_range[1])) + 1
    lon = lon0[lon_loc[0]:lon_loc[1]]
    lat = lat0[lat_loc[0]:lat_loc[1]]
    dep = dep0[dep_loc[0]:dep_loc[1]]
    nlon = len(lon) 
    nlat = len(lat) 
    ndep = len(dep)

    if number_var == 0:
        varz = np.zeros((5,ndep))

    dep3d = np.broadcast_to(dep.reshape(-1, 1, 1), (ndep,nlat,nlon))
    dz_mid = (dep3d[1:, :, :] + dep3d[:-1, :, :]) / 2
    dz = np.zeros_like(dep3d)
    dz[0, :, :] = dep[0] / 2
    dz[1:-1, :, :] = dz_mid[1:, :, :] - dz_mid[:-1, :, :]
    dz[-1, :, :] = (dep[-1] - dz_mid[-1, :, :]) * 2
    
    var_nan = np.squeeze(var0[dep_loc[0]:dep_loc[1], lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
    var1 = var_nan #/ (2191-365) 
    # print(dz)
    var = np.nansum(var1*dz, axis=0)/1e5
    var[var == 0] = np.nan 

    lon_tick = np.arange(115, lon[-1], 15)
    lat_tick = np.arange(40, -15, -10);lat_tick = lat_tick[::-1]

    bar_limit0 = 0.0
    bar_limit1 = 5.0 

    nlevel = np.linspace(bar_limit0, bar_limit1, 33)

    R = 6371000 
    dlon = np.deg2rad(np.diff(lon)[0]) 
    dlat =-np.deg2rad(np.diff(lat)[0])
    lon2d, lat2d = np.meshgrid(lon, lat)
    dy = R * dlat
    dx2 = abs(R * np.cos(np.deg2rad(lat2d)) * dlon)
    dx = np.repeat(dx2[np.newaxis, :, :], ndep, axis=0)
    areas = dx * dy 
    dep3d = np.broadcast_to(dep.reshape(-1, 1, 1), (ndep,nlat,nlon))
    dz_mid = (dep3d[1:, :, :] + dep3d[:-1, :, :]) / 2
    dz = np.zeros_like(dep3d)
    dz[0, :, :] = dep[0] / 2
    dz[1:-1, :, :] = dz_mid[1:, :, :] - dz_mid[:-1, :, :]
    dz[-1, :, :] = (dep[-1] - dz_mid[-1, :, :]) * 2

    varz[number_var, :] = np.nansum(var1*areas, axis=(1,2))

    number_var = number_var + 1

##################################### R-LICOM ######################################################################

    print('contourf nc1')
    cmap = cmaps.WhiteBlueGreenYellowRed
    # fig = plt.figure()
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
    ax.add_feature(cfeature.OCEAN, facecolor='gray')
    ax.add_feature(cfeature.LAND, facecolor='gray')
    print(np.nanmean(var))
    cf = ax.contourf(lon,lat,var, transform=ccrs.PlateCarree(), levels=nlevel, cmap=cmap, extend='both')
    ax.set_extent(
        [lon_range[0], lon_range[1],
        lat_range[0], lat_range[1]],
        crs=ccrs.PlateCarree(),
    )

    ax.set_xticks(lon_tick)
    ax.set_yticks(lat_tick)
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    # ax.tick_params(axis='both', labelsize=15)

    cbar = plt.colorbar(cf, orientation="vertical", shrink=0.7)
    cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
    cbar.set_label('×10⁵J/m²')
    # plt.title(tn1, loc='left', pad=2, fontsize=18)
    plt.tight_layout()
    pn = 'figures/'+ pnstr + f'-{number_var}.png'
    plt.savefig(pn, dpi=600)
    plt.close()
    print(pn)

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(3, 4), sharey=True, sharex=True)
plt.yscale('log')
plt.xscale('log')
axs.plot(varz[0,:]/1e15, dep, color='blue', alpha=0.8, label='Reg-EG',linestyle='-')
axs.plot(varz[1,:]/1e15, dep, color='red', alpha=0.8, label='Reg-JG',linestyle='-')
axs.plot(varz[2,:]/1e15, dep, color='green', alpha=0.8, label='Reg-JL',linestyle='-')
axs.plot(varz[3,:]/1e15, dep, color='blue', alpha=0.8, label='Reg-JG',linestyle='--')
axs.plot(varz[4,:]/1e15, dep, color='green', alpha=0.8, label='Reg-JL',linestyle='--')

axs.set_ylim(1000, 25)
axs.set_xlim(1e-1, 20)
# axs.legend(loc='lower right',frameon=False)
axs.set_ylabel('depth(m)')
axs.set_xlabel('PJ/m')
# axs[2, j].set_xlabel('×10$^1$$^5$ J')
print('finish ploting ' + fn)
# axs[nn].axvspan(days[0], days[365], facecolor='none', alpha=0.6, edgecolor='black', linestyle='--')

plt.tight_layout()
pn = 'figures/'+pnstr+'-depth.png'
plt.savefig(pn, dpi=600)
plt.close()
print(pn)


# import numpy as np
# import netCDF4 as nc 
# import cartopy.crs as ccrs
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# import cartopy.feature as cfeature
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import os
# import matplotlib.cm as cm
# import cmaps 
# cmap = cmaps.WhiteBlueGreenYellowRed

# # lon_range=[95, 150];lat_range=[-16,47]
# lon_range=[100, 148];lat_range=[-15,45];dep_range=[0,500]
# # lon_range=[100, 126];lat_range=[0,25]

# pnstr = 'ape-era-glorys'
# key_str1 = 'ape-era-glorys-nwp.nc' # str1 for ej and jl ; str2 for jg 
# key_str2 = 'ape-era-glorys-nwp.nc'
# vn = 'eke'
# pathstr = ["./data/", "./data/"]
# tn1 = 'a) Model'
# tn2 = 'b) GLORYS'
# tn3 = 'c) difference'
# number_var = 0
# dz = np.zeros((55))
# for path in pathstr:
#     nowdir = os.path.abspath(path)
#     if number_var == 0 :
#         key_str = key_str1
#     elif number_var == 1 :
#         key_str = key_str2
#     for item in sorted(os.listdir(nowdir)):
#         # fn = []
#         if key_str in item and '.png' not in item:
#             fn = path + item

#             print('reading netcdf '+fn+' -'+vn)

#             ncdata = nc.Dataset(fn)
#             var0 = ncdata.variables[vn][:]
#             lon0 = ncdata.variables['lon'][:]
#             lat0 = ncdata.variables['lat'][:]
#             dep0 = abs(ncdata.variables['depth'][:])

#             ncdata_ind = nc.Dataset('ind_final_eas.nc')
#             ind0 = ncdata_ind.variables['ind'][:]

#             lon_loc = [0, 0];lat_loc = [0, 0];dep_loc = [0, 0]
#             lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
#             lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
#             lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
#             lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
#             dep_loc[0] = np.argmin(abs(dep0-dep_range[0]))
#             dep_loc[1] = np.argmin(abs(dep0-dep_range[1])) + 1
#             lon = lon0[lon_loc[0]:lon_loc[1]]
#             lat = lat0[lat_loc[0]:lat_loc[1]]
#             dep = dep0[dep_loc[0]:dep_loc[1]]
#             nlon = len(lon) 
#             nlat = len(lat) 
#             ndep = len(dep)
#             # print(dep)
#             dep3d = np.broadcast_to(dep.reshape(-1, 1, 1), (ndep,nlat,nlon))
#             dz_mid = (dep3d[1:, :, :] + dep3d[:-1, :, :]) / 2
#             dz = np.zeros_like(dep3d)
#             dz[0, :, :] = dep[0] / 2
#             dz[1:-1, :, :] = dz_mid[1:, :, :] - dz_mid[:-1, :, :]
#             dz[-1, :, :] = (dep[-1] - dz_mid[-1, :, :]) * 2
            
#             if number_var == 1:
#                 var_nan = np.squeeze(var0[dep_loc[0]:dep_loc[1],lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
#                 var1 = var_nan#np.nan_to_num(var_nan, nan=0.0) #* ind
#             else:
#                 var_nan = np.squeeze(var0[dep_loc[0]:dep_loc[1], lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
#                 var1 = var_nan#np.nan_to_num(var_nan, nan=0.0) #* ind
#             # var1[var1.mask] = 0
#             # var1[var1==0] = np.nan 
#             var = np.nansum(var1*dz, axis=0)/1e5

#             var = var[np.newaxis, :, :]

#             # lon_tick = np.arange(115, lon[-1], 15)
#             # lat_tick = np.arange(40, -15, -10);lat_tick = lat_tick[::-1]
#             lon_tick = np.arange(115, lon[-1], 15)
#             lat_tick = np.arange(40, -15, -10);lat_tick = lat_tick[::-1]

#             bar_limit0 = 0.0
#             bar_limit1 = 3.0 

#             nlevel = np.linspace(bar_limit0, bar_limit1, 16)

#             if number_var == 0:
#                 var_fig = var
#             else:
#                 var_fig = np.concatenate([var_fig, var], axis=0)

#             var_fig[var_fig==0] = np.nan 
#             number_var = number_var + 1

# ##################################### R-LICOM ######################################################################

# mse = np.nanmean((var_fig[0,:,:].flatten() - var_fig[1,:,:].flatten()) ** 2)
# rmse = np.sqrt(mse)
# # print(rmse)

# print('contourf nc1')

# # fig = plt.figure()
# fig = plt.figure(figsize=(6, 4))
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
# ax.add_feature(cfeature.OCEAN, facecolor='gray')
# ax.add_feature(cfeature.LAND, facecolor='gray')

# cf = ax.contourf(lon,lat,var_fig[0,:,:]-var_fig[1,:,:]+var_fig[1,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap=cmap, extend='both')

# ax.set_xticks(lon_tick)
# ax.set_yticks(lat_tick)
# ax.xaxis.set_major_formatter(LongitudeFormatter())
# ax.yaxis.set_major_formatter(LatitudeFormatter())
# # ax.tick_params(axis='both', labelsize=15)

# cbar = plt.colorbar(cf, orientation="vertical", shrink=0.7)
# cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
# cbar.set_label('×10⁵J/m²')
# # plt.title(tn1, loc='left', pad=2, fontsize=18)
# plt.tight_layout()
# pn = 'figures/'+ pnstr + '-1.png'
# plt.savefig(pn, dpi=600)
# plt.close()
# print(pn)
# ##################################### mercator ######################################################################

# print('contourf nc2')

# # fig = plt.figure()
# fig = plt.figure(figsize=(6, 4))
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
# ax.add_feature(cfeature.OCEAN, facecolor='gray')
# ax.add_feature(cfeature.LAND, facecolor='gray')

# cf = ax.contourf(lon,lat,var_fig[1,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap=cmap, extend='both')

# ax.set_xticks(lon_tick)
# ax.set_yticks(lat_tick)
# ax.xaxis.set_major_formatter(LongitudeFormatter())
# ax.yaxis.set_major_formatter(LatitudeFormatter())
# # ax.tick_params(axis='both', labelsize=15)

# cbar = plt.colorbar(cf, orientation="vertical", shrink=0.7)
# cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
# cbar.set_label('×10⁵J/m²')
# # plt.title(tn2, loc='left', pad=2, fontsize=18)
# plt.tight_layout()
# pn = 'figures/'+ pnstr + '-2.png'
# plt.savefig(pn, dpi=600)
# plt.close()
# print(pn)
# ##################################### v2503 - mercator ######################################################################
# if vn == 'tt':
#     bar_limit0 = -2
#     bar_limit1 = 2
# elif vn == 'ss':
#     bar_limit0 = -1
#     bar_limit1 = 1

# bar_limit0 =-0.1
# bar_limit1 = 0.1

# nlevel = np.linspace(bar_limit0, bar_limit1, 64)

# print('contourf 1 - 2')

# # fig = plt.figure()
# fig = plt.figure(figsize=(4, 3))
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5)
# ax.add_feature(cfeature.OCEAN, facecolor='gray')
# ax.add_feature(cfeature.LAND, facecolor='gray')

# cf = ax.contourf(lon,lat,var_fig[0,:,:]-var_fig[1,:,:], transform=ccrs.PlateCarree(), levels=nlevel, cmap=cmap, extend='both')

# ax.set_xticks(lon_tick)
# ax.set_yticks(lat_tick)
# ax.xaxis.set_major_formatter(LongitudeFormatter())
# ax.yaxis.set_major_formatter(LatitudeFormatter())
# # ax.tick_params(axis='both', labelsize=15)

# cbar = plt.colorbar(cf, orientation="vertical", shrink=0.7)
# cbar.set_ticks(np.arange(bar_limit0, bar_limit1+(bar_limit1-bar_limit0)*0.25, (bar_limit1-bar_limit0)/2))
# if vn == 'ss': 
#     cbar.set_label('psu')
#     plt.title(f'RMSE={rmse:.3f}psu', loc='right', pad=2, fontsize=12)
# elif vn == 'tt':
#     cbar.set_label('℃')
#     plt.title(f'RMSE={rmse:.2f}℃', loc='right', pad=2, fontsize=12)
# elif vn == 'eke' or vn == 'vv':
#     # cbar.set_label(r'$m^2/s^2$')
#     # plt.title(f'RMSE={rmse:.2f}'+r'$m^2/s^2$', loc='right', pad=2, fontsize=12)
#     cbar.set_label('m²/s²')
#     plt.title(f'RMSE={rmse:.3f}m²/s²', loc='right', pad=2, fontsize=12)
# # plt.title(tn3, loc='left', pad=2, fontsize=18)
# plt.tight_layout()
# pn = 'figures/'+ pnstr + '-1m2.png'
# plt.savefig(pn, dpi=600)
# plt.close()
# print(pn)
