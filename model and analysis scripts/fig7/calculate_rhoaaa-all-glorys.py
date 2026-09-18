import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import os

# lon_range=[95.3, 149.8];lat_range=[-15.7,47.8]
# lon_range=[100, 148];lat_range=[-15,45]
lon_range=[100, 150];lat_range=[-16,48]
# lon_range=[105, 125];lat_range=[0, 25]
# lon_range=[105, 140];lat_range=[0,40]

# key_strs = ['uu2014','uu2015','uu2016','uu2017','uu2018']
key_strs = ['uu-2014','uu-2015','uu-2016','uu-2017','uu-2018']
vn = 'uu'
vn2= 'vv'
vn3= 'rho'

# pnstr = f'edgepwork-eke-jra-licom3'
# # path = "../output-FaD-era/fort.22/"
# # path = "../output-jra-glorys/fort.22/"
# path = "../output-FaD-jra/fort.22/"
# # path = "../../mercatordata/mercator_compare/RegionalData/"

#pnstrs = [f'edgerhoaaa-jl']
#paths = ["../output-jra-licom3/fort.22/"]
pnstrs = [f'TH/edgerhoaaa-glorys-all']
paths = ["../exe/boundary-data-glorys-daily/"]

for nforce in range(1):
    pnstr = pnstrs[nforce]
    path = paths[nforce]

    nowdir = os.path.abspath(path)
    number_day = 0 
    regionalmean = np.zeros((2400,1))
    weighted_averages = np.zeros((55))
    meansst = np.zeros((2400,1))
    dz = np.zeros((55))

    # dsind = xr.open_dataset('ind_final_eas_uv.nc', decode_times=False)
    # ind0 = dsind['ind'].values[:,:,:] + 0.0 
    # ind0[ind0<0.5] = np.nan 

    # dsdep = xr.open_dataset('./optimized_ocean_bathymetry_10_RL.nc', decode_times=False)
    # oceandepth0 = dsdep['depth'].values

    # dsxm = xr.open_dataset(path + 'tx-2014-2018.nc', decode_times=False)
    # dsym = xr.open_dataset(path + 'ty-2014-2018.nc', decode_times=False)
    dsum = xr.open_dataset(path + 'ave/uu-2014-2018.nc', decode_times=False)
    dsvm = xr.open_dataset(path + 'ave/vv-2014-2018.nc', decode_times=False)
    dsrm = xr.open_dataset(path + 'ave/rho-2014-2018.nc', decode_times=False)
    # xm = dsxm['tx'].values[0,0,:,:]
    # ym = dsym['ty'].values[0,0,:,:]
    um = dsum['uu'].values[0,1:,:,:]
    vm = dsvm['vv'].values[0,1:,:,:]
    rhom = dsrm['rho'].values[0,1:,:,:]
    lon0 = dsum['lon'].values
    lat0 = dsum['lat'].values
    depth = abs(dsum['depth'].values[1:]) 

    g = 9.8 
    rho0 = 1025 

    for item in sorted(os.listdir(nowdir)):
        fn = []
        if any(key in item for key in key_strs):
            fn = path + item
            fn2= fn.replace(vn,vn2)
            fn3= fn.replace(vn,vn3)

            print('reading netcdf '+fn+' -'+vn)
            print(fn3)

            ds = xr.open_dataset(fn, decode_times=False)
            var0 = ds[vn].values[0,1:,:,:] - um 
            
            ds2 = xr.open_dataset(fn2, decode_times=False)
            var2 = ds2[vn2].values[0,1:,:,:] -vm  

            ds3 = xr.open_dataset(fn3, decode_times=False)
            var3 = ds3[vn3].values[0,1:,:,:]  

            lon_loc = [0, 0];lat_loc = [0, 0]
            lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
            lon_loc[1] = np.argmin(abs(lon0-lon_range[1])) + 1
            lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
            lat_loc[1] = np.argmin(abs(lat0-lat_range[0])) + 1 
            lon = lon0[lon_loc[0]:lon_loc[1]]
            lat = lat0[lat_loc[0]:lat_loc[1]]
            nlon = len(lon) 
            nlat = len(lat) 
            ndep = len(depth)
            if number_day == 0 :
                deudx_e = np.zeros((2400,ndep,nlat))
                deudx_w = np.zeros((2400,ndep,nlat))
                devdy_s = np.zeros((2400,ndep,nlon))
                devdy_n = np.zeros((2400,ndep,nlon))

            R = 6371000 
            dlon = np.deg2rad(np.diff(lon)[0]) 
            dlat = -np.deg2rad(np.diff(lat)[0])

            lon2d, lat2d = np.meshgrid(lon, lat)

            # Calculate grid cell areas
            dy = R * dlat
            dx2 = abs(R * np.cos(np.deg2rad(lat2d)) * dlon)
            dx = np.repeat(dx2[np.newaxis, :, :], ndep, axis=0)
            areas = dx * dy  
            dz_mid = (depth[1:] + depth[:-1]) / 2
            dz = np.zeros_like(depth)
            dz[0] = depth[0] / 2
            dz[1:-1] = dz_mid[1:] - dz_mid[:-1]
            dz[-1] = (depth[-1] - dz_mid[-1]) * 2

            # varx = np.squeeze(var3[lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
            # vary = np.squeeze(var4[lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]])
            varu = np.squeeze(var0[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            varv = np.squeeze(var2[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            varrho = np.squeeze(var3[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 
            rhomin = np.squeeze(rhom[:, lat_loc[0]:lat_loc[1], lon_loc[0]:lon_loc[1]]) 

            areas = areas * varu/varu  

            # get reference N^2 
            # rhoref = np.nansum(areas * rhomin, axis=(1,2))/np.nansum(areas, axis=(1,2))
            fntxt = f'rhoref-glorys.txt'
            print(fntxt)
            with open(fntxt, 'r') as f:
                rhoref = np.array(eval(f.read()))
            rhoa = rhomin - rhoref[:, np.newaxis, np.newaxis] 
            rhoaa= varrho - rhomin 
            N2 = g / rho0 * np.gradient(np.squeeze(rhoref), axis=0) / dz 
            cc = g**2 / (rho0 * N2[:, np.newaxis, np.newaxis]) 

            # varpu = - cc * varu * rhoa * rhoaa 
            # varpv = - cc * varv * rhoa * rhoaa 
            if number_day == 0:
                varpu = np.zeros_like(varu)
                varpv = np.zeros_like(varu)
            varpu += - cc * varu * rhoa * rhoaa 
            varpv += - cc * varv * rhoa * rhoaa 
            print(np.nanmean(varpu,axis=(1,2)))
            print(varpu.shape)

            # deudx_e[number_day, :, :] = (varpu[:, :,-1] - varpu[:, :,-2]) / ((dx[:, :,-1]+dx[:, :,-2])/2) 
            # deudx_w[number_day, :, :] = (varpu[:, :, 1] - varpu[:, :, 0]) / ((dx[:, :, 1]+dx[:, :, 0])/2) 
            # devdy_s[number_day, :, :] = (varpv[:, -1,:] - varpv[:, -2,:]) / dy 
            # devdy_n[number_day, :, :] = (varpv[:,  1,:] - varpv[:,  0,:]) / dy 
            # deudx_e[number_day, :, :] = (varpu[:, :,-1]) 
            # deudx_w[number_day, :, :] = (varpu[:, :, 0]) 
            # devdy_s[number_day, :, :] = (varpv[:, -1,:]) 
            # devdy_n[number_day, :, :] = (varpv[:,  0,:]) 
            # print(varpu.shape)

            print([np.nansum(deudx_e[number_day, 0, :]),np.nansum(deudx_w[number_day, 0, :]),
                np.nansum(devdy_s[number_day, 0, :]),np.nansum(devdy_n[number_day, 0, :])])
            number_day = number_day + 1

            ds.close()
            ds2.close()
            ds3.close()

    varpu = varpu / number_day 
    varpv = varpv / number_day
    print(np.nanmean(varpu))

    def save_to_netcdf(filename, varpu, varpv, lon, lat, depth, time_index=None):

        if time_index is None:
            ds = xr.Dataset(
                {
                    "varpu": (("depth", "lat", "lon"), varpu),
                    "varpv": (("depth", "lat", "lon"), varpv),
                },
                coords={
                    "depth": depth,
                    "lat": lat,
                    "lon": lon,
                },
            )
        else:
            ds = xr.Dataset(
                {
                    "varpu": (("time", "depth", "lat", "lon"), varpu[np.newaxis, ...]),
                    "varpv": (("time", "depth", "lat", "lon"), varpv[np.newaxis, ...]),
                },
                coords={
                    "time": [time_index],
                    "depth": depth,
                    "lat": lat,
                    "lon": lon,
                },
            )

        ds.to_netcdf(filename)
        ds.close()


    fnot = 'data/' + pnstr + '.nc'
    save_to_netcdf(fnot, varpu, varpv, lon, lat, depth)
