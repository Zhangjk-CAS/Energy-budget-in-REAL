import numpy as np
import netCDF4 as nc 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import convolve1d

# colors = ['black', 'blue', 'red', 'green', 'blue', 'green']
# linest = ['-', '-', '-', '-', '--', '--']
# legends= ['dailyBC', 'Reg-EG', 'Reg-JG', 'Reg-JL', 'GLORYS', 'LICOM3']
colors = ['blue', 'red', 'green', 'blue', 'green']
linest = ['-', '-', '-', '--', '--']
legends= ['Reg-EG', 'Reg-JG', 'Reg-JL', 'GLORYS', 'LICOM3']
# for forcestr in ['mke-eg-dailyBC', 'mke-eg', 'mke-jg', 'mke-jl', 'mke-glorys']: 
# for term in ['ape', 'eape', 'umke', 'ueke', 'rhoaaa', 'uaemke', 'up', 'uapa']:
for term in ['umke']:
    nn = 0
    # for forcestr in ['eg-dailyBC', 'eg', 'jg', 'jl', 'glorys', 'licom3']: 
    fig, axs = plt.subplots(1, 1, figsize=(5, 2), sharey=True)
    for forcestr in ['eg', 'jg', 'jl', 'glorys', 'licom3']: 
    # for forcestr in ['eg-dailyBC']: 

        pnstr = f'edgework-{term}-with-distance-TH' #f'edgework-{term}-distance'
        fn = f'data/TH/edge{term}-{forcestr}-all.nc'
        print(f'reading {fn}')

        lon_range=[120, 149];lat_range=[-15,45]

        ntime = 2191 - 365 
        ncpu = nc.Dataset(fn)

        varpu0 =-ncpu.variables['varpu'][:].filled(np.nan)
        varpv0 = ncpu.variables['varpv'][:].filled(np.nan)

        lon0 = ncpu.variables['lon'][:]
        lat0 = ncpu.variables['lat'][:]
        if nn == 0:
            depth= abs(ncpu.variables['depth'][:])

        lon_loc = [0, 0];lat_loc = [0, 0]
        lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
        lon_loc[1] = np.argmin(abs(lon0-lon_range[1])) + 1
        lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) 
        lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) + 1
        lon = lon0[lon_loc[0]:lon_loc[1]]
        lat = lat0[lat_loc[0]:lat_loc[1]]
        nlon = len(lon);nlat = len(lat);ndep = len(depth)

        R = 6371000 
        dlon = np.deg2rad(np.diff(lon)[0]) 
        dlat =-np.deg2rad(np.diff(lat)[0])

        # Calculate grid distance 
        dy = R * dlat
        dx = abs(R * np.cos(np.deg2rad(lat)) * dlon)
        dz_mid = (depth[1:] + depth[:-1]) / 2
        dz = np.zeros_like(depth)
        dz[0] = depth[0] / 2
        dz[1:-1] = dz_mid[1:] - dz_mid[:-1]
        dz[-1] = (depth[-1] - dz_mid[-1]) * 2

        varpu = varpu0[:,lat_loc[0]:lat_loc[1],lon_loc[0]:lon_loc[1]] 
        varpv = varpv0[:,lat_loc[0]:lat_loc[1],lon_loc[0]:lon_loc[1]] 
        # print(np.nanmean(varpu,axis=(1,2)))

        vare_itg = 1e-9*np.nansum(varpu * dz[:, np.newaxis, np.newaxis] * dy, axis=(0,1));vare_itg[vare_itg==0] = np.nan 

        e_loc = np.argmin(abs(lon0-148))
        vare = varpu0[:,lat_loc[0]:lat_loc[1],e_loc]
        vare_total = np.nansum(vare * dz[:, np.newaxis] * dy)
        print(np.nanmean(vare))
        print(vare.shape)
        print(lon0[e_loc])
        print(f'energy input at east boundary is {vare_total/1e9}')

        print(varpu.shape)
        # print(vare_itg)

        # if term == 'up' or term == 'uapa':
        #     axs2 = axs.twinx()
        #     axs2.plot(lon, vare_itg, color=colors[nn], linestyle=[nn], alpha=0.8)
        #     axs2.set_ylim([-300,300])
        # else :
        #     axs.plot(lon, vare_itg, color=colors[nn], alpha=0.1)
        #     axs.set_ylim([-50,50])

        axs.plot(lon, vare_itg, color=colors[nn], linestyle=linest[nn], label=legends[nn], alpha=0.8)
        # axs.set_ylim([-20,200])
        axs.set_ylim([10,-40])
        # axs.set_ylim([-10,40])
        xticks = np.arange(120, 151, 5)  # 0, 20, 40, ..., 360
        xtick_labels = [f'{int(x)}°E' for x in xticks]
        axs.set_xticks(xticks)
        axs.set_xticklabels(xtick_labels)
        axs.set_ylabel('GW')
        axs.legend(fontsize=8, frameon=False, ncol=2)
        nn = nn + 1 

        # plt.tight_layout()
    pn = 'figures/' + pnstr + '.png'
    plt.savefig(pn, dpi=600)
    plt.close()
