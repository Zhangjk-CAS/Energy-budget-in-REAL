import numpy as np
import netCDF4 as nc
import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt


lon_range=[100, 148];lat_range=[-15,45];dep_range=[0,1000]

# pathstr = ["../output-FaD-era/", "../output-jra-glorys/", "../output-FaD-jra/", 
#            "../../mercatordata/mercator_compare/RegionalData/annual-ave/", "./", "../../observation/woa/"]
# pathstr = ["../output-era-glorys/fort.22/TH/", "../output-jra-glorys/", "../output-FaD-jra/", 
#            "../../mercatordata/mercator_compare/RegionalData/annual-ave/", "../../observation/woa/"]
# pathstr = ["../output-FaD-era/", 
#            "../../mercatordata/mercator_compare/RegionalData/annual-ave/", "./data/data-licom3high/", "../../observation/woa/"]
pathstr = ["../ave-TH/woa/","../ave-TH/glorys/", "../ave-TH/eg/", "../ave-TH/jg/", "../ave-TH/jl/st/"
           ]

pnstr = 'vertical-ss-sum-minusWOA'

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(2.6, 4),)

# labelstr = ['Reg-EG','Reg-JG','Reg-JL''GLORYS','Argo','WOA']
labelstr = ['WOA','GLORYS','Reg-EG','Reg-JG','Reg-JL']
# labelstr = ['Reg-EG','GLORYS','WOA']
colors = ['green','blue','blue','red','green']
# colors = ['blue','blue','gray']
lnstr = ['--', '--','-', '-', '-']

number_var = 0
for path in pathstr:
    number_day = 0 
    if pnstr[0:11] == 'vertical-ss':
        # key_str = ['ss-2014-2018', 'ss-2014-2018', 'ss-2014-2018', 
        #            'ss-2014-2018', 'woa-ss2015-2022.nc'][number_var]
        # key_str = ['12.1.2-ss-woa.nc','12.1.2-5year-ss.nc','12.1.2-5year-ss.nc', '12.1.2-5year-ss.nc', '12.1.2-5year-ss.nc'
        #             ][number_var]
        key_str = ['6.7.8-ss-woa.nc','6.7.8-5year-ss.nc','6.7.8-5year-ss.nc', '6.7.8-5year-ss.nc', '6.7.8-5year-ss.nc'
                    ][number_var]
        # key_str = ['6.7.8-5year-ss.nc', '6.7.8-5year-ss.nc', '6.7.8-5year-ss.nc', 
        #            '6.7.8-5year-ss.nc', '6.7.8-ss-woa.nc'][number_var]
        vn = ['s_an', 'ss', 'ss', 'ss', 'ss'][number_var] 
        xn = 'psu' 
    elif pnstr[0:11] == 'vertical-tt':
        # key_str = ['tt-2014-2018', 'tt-2014-2018', 'tt-2014-2018', 
        #            'tt-2014-2018', 'woa-tt2015-2022.nc'][number_var]
        # key_str = ['12.1.2-5year-tt.nc', '12.1.2-5year-tt.nc', '12.1.2-5year-tt.nc', 
        #            '12.1.2-5year-tt.nc', '12.1.2-tt-woa.nc'][number_var]
        # key_str = ['6.7.8-tt-woa.nc','6.7.8-5year-tt.nc', '6.7.8-5year-tt.nc', '6.7.8-5year-tt.nc', '6.7.8-5year-tt.nc' 
        #            ][number_var]
        key_str = ['12.1.2-tt-woa.nc','12.1.2-5year-tt.nc', '12.1.2-5year-tt.nc', '12.1.2-5year-tt.nc', '12.1.2-5year-tt.nc' 
                   ][number_var]
        vn = ['t_an','tt', 'tt', 'tt', 'tt'][number_var] 
        
        xn = '℃'

    nowdir = os.path.abspath(path)
    for item in sorted(os.listdir(nowdir)):

        if key_str in item and '.png' not in item:
            fn = path + item

            print('reading netcdf '+fn+' -'+vn)

            ncdata = nc.Dataset(fn)
            if number_var < 1.5:
                var0 = ncdata.variables[vn][:].filled(np.nan)
                lon0 = ncdata.variables['lon'][:]
                lat0 = ncdata.variables['lat'][:]
                depth= ncdata.variables['depth'][:]

            # elif number_var == 1:
            #     var0 = ncdata.variables[vn][:,1:,:,:].filled(np.nan)
            #     lon0 = ncdata.variables['lon'][:]
            #     lat0 = ncdata.variables['lat'][:]
            #     depth= ncdata.variables['depth'][1:]
            elif number_var == 4:
                var0 = ncdata.variables[vn][:].filled(np.nan)
                lon0 = ncdata.variables['lon'][:]
                lat0 = ncdata.variables['lat'][:]
                depth= abs(ncdata.variables['lev'][:])
            else :
                var0 = ncdata.variables[vn][:].filled(np.nan)
                lon0 = ncdata.variables['lon'][:]
                lat0 = ncdata.variables['lat'][:]
                depth= abs(ncdata.variables['lev1'][:])

            lon_loc = [0, 0];lat_loc = [0, 0];dep_loc = [0, 0]
            lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
            lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))
            lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
            lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) 
            # if number_var == 2 :
            #     lat_loc[0] = np.argmin(abs(lat0-lat_range[0])) 
            #     lat_loc[1] = np.argmin(abs(lat0-lat_range[1])) 
            dep_loc[0] = np.argmin(abs(depth-dep_range[0]))
            dep_loc[1] = np.argmin(abs(depth-dep_range[1]))
            lon = lon0[lon_loc[0]:lon_loc[1]]
            lat = lat0[lat_loc[0]:lat_loc[1]]
            dep = depth[dep_loc[0]:dep_loc[1]+1]
            ndep = dep.shape[0]
            print(lat.shape)
            print(var0.shape)

            weighted_averages = np.zeros((ndep))
            dz = np.zeros((ndep))

            dz[0] = dep[0] - 0
            dz[1:] = np.diff(dep)
            R = 6371000  # Earth radius in meters
            dz = abs(dz)
            R = 6371000 
            dlon = np.deg2rad(np.diff(lon)[0]) 
            dlat = -np.deg2rad(np.diff(lat)[0])

            # Create 2D arrays of latitude and longitude
            lon2d, lat2d = np.meshgrid(lon, lat)

            # Calculate grid cell areas
            dy = R * dlat
            dx = R * np.cos(np.deg2rad(lat2d)) * dlon
            areas = dx * dy  # 2D array of cell areas

            ind = np.squeeze(var0[0, :, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])

            for lev in range(0,ndep):

                var = np.squeeze(var0[0,lev, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
                areas = areas 

                ind[lev, :, :] = var

                weighted_sum = np.nansum(var * areas)
                total_area = np.nansum(areas[~np.isnan(var)])
                weighted_averages[lev] = weighted_sum / total_area 

            print(weighted_averages)
            
            if number_var == 0:
                ref0 = weighted_averages.copy()
                depref = dep 
                
            ref = np.interp(dep, depref, ref0)
            weighted_averages = ref - weighted_averages 

    axs.plot(weighted_averages, dep, alpha=0.8, label=labelstr[number_var], 
             color=colors[number_var], linestyle=lnstr[number_var]) 
    axs.legend(loc='upper left',frameon=False)
    # axs.legend(loc='lower left',frameon=False)
    axs.set_yscale('log')
    # axs.set_xscale('log')
    axs.set_ylim(1000,2.5)
    # axs.invert_yaxis()
    axs.set_xlabel(xn)
    axs.set_ylabel('depth(m)')
    
    number_var = number_var + 1 

pn = 'figures/' + pnstr + '-bias.png'
plt.savefig(pn, dpi=600)
plt.close()
print('finish:'+ pn)