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
import warnings

warnings.filterwarnings("ignore", message="Mean of empty slice")

fig, axs = plt.subplots(1, 2, figsize=(5, 5), sharey=True)

terms = [('edge3pwork-mke-era-glorys', 
          'edge3pwork-mke-jra-glorys', 
          'edge3pwork-mke-jra-licom3', 
          'edgepwork-mke-glorys', 
          'data-licom3high/edgepwork-mke-licom3Global'),
              
         ('edge3pwork-era-glorys', 
          'edge3pwork-jra-glorys', 
          'edge3pwork-jra-licom3', 
          'edgepwork-eke-glorys', 
          'data-licom3high/edgepwork-licom3Global')] # for mke 

colors = ['blue', 'red', 'green', 'blue', 'green'] 
tn = ['(a)', '(b)', '(c)', '(d)'] 
ln = ['Reg-EG','Reg-JG', 'Reg-JL', 'GLORYS', 'LICOM3']
linesty = ['-', '-', '-', '--', '--']

for nnterm in range(2):
    nnterm1 = nnterm // 2
    nnterm2 = nnterm % 2
    for nndata in range(5): 

        budget = terms[nnterm][nndata]
        pnstr = f'{budget}'
        fn1 = f'data/{budget}.nc'

        lon_range=[100, 148];lat_range=[-15,45]

        ntime = 2191 - 365 
        ncpu = nc.Dataset(fn1)
        varpu_e =-ncpu.variables['deudx_e'][:].filled(np.nan)
        varpu_w = ncpu.variables['deudx_w'][:].filled(np.nan)
        varpu_s = ncpu.variables['devdy_s'][:].filled(np.nan)
        varpu_n =-ncpu.variables['devdy_n'][:].filled(np.nan)

        lon0 = ncpu.variables['lon'][:]
        lat0 = ncpu.variables['lat'][:]
        # depth= abs(ncpu.variables['depth'][:])
        depth = [-2.50300002, -7.52099991, -12.5748997, -17.6883011, -22.8841515, 
                    -28.1846504, -33.6110001, -39.1831512, -44.9196472, -50.8374481, 
                    -56.9518509, -63.2761993, -69.8218002, -76.5978546, -83.6112976, 
                    -90.8668518, -98.3668518, -106.111298, -114.097855, -122.3218, 
                    -130.776199, -139.451843, -148.337463, -157.419647, -166.683151, 
                    -176.110992, -185.684647, -195.384155, -205.188293, -215.07489, 
                    -225.020996, -235.002991, -245.519348, -258.625793, -277.897034, 
                    -307.315247, -350.673523, -411.505432, -493.018921, -598.036621, 
                    -728.943481, -887.64209, -1075.51746, -1293.41052, -1541.604, 
                    -1819.81592, -2127.20459, -2462.38452, -2823.45361, -3208.02686, 
                    -3613.28247, -4036.01074, -4472.67969, -4919.49561, -5372.47607] 
        if nndata == 3: 
            depth = [0, -2.50300002, -7.52099991, -12.5748997, -17.6883011, -22.8841515, 
                    -28.1846504, -33.6110001, -39.1831512, -44.9196472, -50.8374481, 
                    -56.9518509, -63.2761993, -69.8218002, -76.5978546, -83.6112976, 
                    -90.8668518, -98.3668518, -106.111298, -114.097855, -122.3218, 
                    -130.776199, -139.451843, -148.337463, -157.419647, -166.683151, 
                    -176.110992, -185.684647, -195.384155, -205.188293, -215.07489, 
                    -225.020996, -235.002991, -245.519348, -258.625793, -277.897034, 
                    -307.315247, -350.673523, -411.505432, -493.018921, -598.036621, 
                    -728.943481, -887.64209, -1075.51746, -1293.41052, -1541.604, 
                    -1819.81592, -2127.20459, -2462.38452, -2823.45361, -3208.02686, 
                    -3613.28247, -4036.01074, -4472.67969, -4919.49561] 
        
        lev_n = 55    
        depth = np.array(depth)  
        depth = -depth

        lon_loc = [0, 0];lat_loc = [0, 0]
        lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
        lon_loc[1] = np.argmin(abs(lon0-lon_range[1])) + 1
        lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) 
        lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) + 1
        lon = lon0[lon_loc[0]:lon_loc[1]]
        lat = lat0[lat_loc[0]:lat_loc[1]]
        nlon = len(lon);nlat = len(lat);ndep = len(depth)

        if nnterm in [0,2]:
            vare_pu = (varpu_e[0, :lev_n, lat_loc[0]:lat_loc[1]]) 
            varw_pu = (varpu_w[0, :lev_n, lat_loc[0]:lat_loc[1]]) 
            vars_pu = (varpu_s[0, :lev_n, lon_loc[0]:lon_loc[1]]) 
            varn_pu = (varpu_n[0, :lev_n, lon_loc[0]:lon_loc[1]]) 
        elif nnterm in [1,3]:
            vare_pu = np.nanmean(varpu_e[0:ntime, :lev_n, lat_loc[0]:lat_loc[1]], axis=0) 
            varw_pu = np.nanmean(varpu_w[0:ntime, :lev_n, lat_loc[0]:lat_loc[1]], axis=0) 
            vars_pu = np.nanmean(varpu_s[0:ntime, :lev_n, lon_loc[0]:lon_loc[1]], axis=0) 
            varn_pu = np.nanmean(varpu_n[0:ntime, :lev_n, lon_loc[0]:lon_loc[1]], axis=0) 

        vare = vare_pu 
        varw = varw_pu 
        vars = vars_pu 
        varn = varn_pu 
        kernel_y = np.ones(30) / 30
        vare_sm = convolve1d(vare, weights=kernel_y, axis=1, mode='nearest')
        varw_sm = convolve1d(varw, weights=kernel_y, axis=1, mode='nearest')
        vars_sm = convolve1d(vars, weights=kernel_y, axis=1, mode='nearest')
        varn_sm = convolve1d(varn, weights=kernel_y, axis=1, mode='nearest')

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

        vare_itg = 1e-5*np.nansum(vare * dz[:, np.newaxis], axis=0);vare_itg[vare_itg==0] = np.nan 
        varw_itg = 1e-5*np.nansum(varw * dz[:, np.newaxis], axis=0);varw_itg[varw_itg==0] = np.nan 
        vars_itg = 1e-5*np.nansum(vars * dz[:, np.newaxis], axis=0);vars_itg[vars_itg==0] = np.nan 
        varn_itg = 1e-5*np.nansum(varn * dz[:, np.newaxis], axis=0);varn_itg[varn_itg==0] = np.nan 

        vare_total = np.nansum(vare * dz[:, np.newaxis] * dy)
        varw_total = np.nansum(varw * dz[:, np.newaxis] * dy)
        vars_total = np.nansum(vars * dz[:, np.newaxis] * dx[-1])
        varn_total = np.nansum(varn * dz[:, np.newaxis] * dx[ 0])

        vare_z = np.nansum(vare * dy, axis=1)
        varw_z = np.nansum(varw * dy, axis=1)
        vars_z = np.nansum(vars * dx[-1], axis=1)
        varn_z = np.nansum(varn * dx[ 0], axis=1)
        var_z = vare_z + varw_z + vars_z + varn_z 
        print(var_z.shape)

        print('-------------------------')
        print(f'{budget}:')
        print([vare_total/1e9, varw_total/1e9, vars_total/1e9, varn_total/1e9])
        print(vare_total/1e9+ varw_total/1e9+ vars_total/1e9+ varn_total/1e9)

        # axs[nnterm].plot(lat, vare_itg, color=colors[nndata], alpha=0.2, 
        #                  linewidth=1.5, linestyle=linesty[nndata])

        # kernel_y = np.ones(30) / 30
        # vare_itg = convolve1d(vare_itg, weights=kernel_y, axis=0, mode='nearest')
        if nndata < 3:
            axs[nnterm].plot(var_z/1e9, depth, color=colors[nndata], linewidth=0.9, linestyle='-')
        else:
            axs[nnterm].plot(var_z/1e9, depth, color=colors[nndata], linewidth=0.9, linestyle='--')

        # axs[nnterm].plot(lat, vare_itg, color=colors[nndata], linewidth=0.9, linestyle='--')

        # if nndata == 0:
        #     axs[nnterm].set_title(tn[nnterm], loc='left', pad=0.1)


axs[0].set_ylim(1000, 2.5)
# axs[1].set_ylim(1000, 2.5)
# axs[1].set_ylim(1000, 0)
# axs[0].set_xticks([])
# axs[0,1].set_xticks([])
# axs[1].set_xticks(range(-10,45,10))
# axs[1,1].set_xticks(range(-15,46,15))
# axs[0].set_xlim(lat_range)
# axs[1].set_xlim(lat_range)
# axs[0].set_xlim(lat_range)
# axs[1].set_xlim(lat_range)

# latticklabel = ['10°S',  '0°N', '10°N', '20°N', '30°N', '40°N']
# axs[1].set_xticklabels(latticklabel)
# axs[1,1].set_xticklabels(latticklabel)
axs[0].set_ylabel('depth(m)')
axs[0].set_xlabel('GW/m')
# axs[1].set_ylabel('×10⁵W/m')

pnstr2 = f'edgework-itg-pwork-dep'
pn = 'figures/' + pnstr2 + '.png'
plt.savefig(pn, dpi=600)
plt.close()