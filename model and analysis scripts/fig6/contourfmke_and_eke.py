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

pnstr = 'mke'
vn = 'eke'
path = './data/'

number_var = 0
for item in ['mke-era-glorys-nwp.nc','mke-jra-glorys-nwp.nc','mke-jra-licom3-nwp.nc'
             ,'mke-glorys-nwp.nc','data-licom3high/mke-licom3Globalnwp.nc']:
# for item in ['eke-era-glorys-nwp.nc','eke-jra-glorys-nwp.nc','eke-jra-licom3-nwp.nc'
#              ,'eke-glorys-nwp.nc','data-licom3high/eke-licom3Globalnwp.nc']:
# for item in ['eape-era-glorys-nwp.nc','eape-jra-glorys-nwp.nc','eape-jra-licom3-nwp.nc'
#              ,'eape-glorys-nwp.nc','data-licom3high/eape-licom3Globalnwp.nc']:

    fn = path + item

    print('reading netcdf '+fn+' -'+vn)

    if number_var == 3:
        ncdata = nc.Dataset(fn)
        var0 = ncdata.variables[vn][1:,:,:].filled(np.nan)
        lon0 = ncdata.variables['lon'][:]
        lat0 = ncdata.variables['lat'][:]
    else:
        ncdata = nc.Dataset(fn)
        var0 = ncdata.variables[vn][:,:,:].filled(np.nan)
        lon0 = ncdata.variables['lon'][:]
        lat0 = ncdata.variables['lat'][:]
    
    if number_var == 0:
        dep0 = abs(ncdata.variables['depth'][:])

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
    if pnstr == 'mke':
        var1 = var_nan 
    else:
        var1 = var_nan / (2191-365) 

    var = np.nansum(var1*dz, axis=0)/1e5
    var[var == 0] = np.nan 

    lon_tick = np.arange(115, lon[-1], 15)
    lat_tick = np.arange(40, -15, -10);lat_tick = lat_tick[::-1]

    bar_limit0 = 0.0
    bar_limit1 = 1.0 

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
    print(varz[number_var, :])

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

    cf = ax.contourf(lon,lat,var, transform=ccrs.PlateCarree(), levels=nlevel, cmap=cmap, extend='both')

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

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(3, 4), sharey=True, sharex=True)
plt.yscale('log')
plt.xscale('log')
axs.plot(varz[0,:]/1e15, dep, color='blue', alpha=0.8, label='Reg-EG',linestyle='-')
axs.plot(varz[1,:]/1e15, dep, color='red', alpha=0.8, label='Reg-JG',linestyle='-')
axs.plot(varz[2,:]/1e15, dep, color='green', alpha=0.8, label='Reg-JL',linestyle='-')
axs.plot(varz[3,:]/1e15, dep, color='blue', alpha=0.8, label='GLORYS',linestyle='--')
axs.plot(varz[4,:]/1e15, dep, color='green', alpha=0.8, label='LICOM3',linestyle='--')

axs.set_ylim(5000, 25)
axs.set_xlim(1e-5, 10)
# axs.legend(loc='lower right',frameon=False)
axs.set_ylabel('depth(m)')
axs.set_xlabel('PW/m')
# axs[2, j].set_xlabel('×10$^1$$^5$ J')
print('finish ploting ' + fn)
# axs[nn].axvspan(days[0], days[365], facecolor='none', alpha=0.6, edgecolor='black', linestyle='--')

plt.tight_layout()
pn = 'figures/'+pnstr+'-depth.png'
plt.savefig(pn, dpi=600)
plt.close()