import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# === 1. 读取数据 ===
ds = xr.open_dataset('optimized_ocean_bathymetry_10_GL.nc')
depth0 = ds['depth']   # (lat, lon)

lon0 = ds['lon'].values
lat0 = ds['lat'].values
print('1')
lon_range = [95, 150];lat_range = [-16, 48]
lon_loc = [0, 0];lat_loc = [0, 0]
lon_loc[0] = np.argmin(abs(lon0-lon_range[0]))
lon_loc[1] = np.argmin(abs(lon0-lon_range[1]))+1
lat_loc[0] = np.argmin(abs(lat0-lat_range[1])) #reversed latitude
lat_loc[1] = np.argmin(abs(lat0-lat_range[0]))+1 
lon = lon0[lon_loc[0]-200:lon_loc[1]+200]
lat = lat0[lat_loc[0]-200:lat_loc[1]+200]

n = 4
depth0[lat_loc[0]:lat_loc[1], lon_loc[0]-n:lon_loc[0]+n] = np.nan 
depth0[lat_loc[0]:lat_loc[1], lon_loc[1]-n:lon_loc[1]+n] = np.nan 
depth0[lat_loc[0]-n:lat_loc[0]+n, lon_loc[0]:lon_loc[1]] = np.nan 
depth0[lat_loc[1]-n:lat_loc[1]+n, lon_loc[0]:lon_loc[1]] = np.nan 
depth = depth0[lat_loc[0]-200:lat_loc[1]+200, lon_loc[0]-200:lon_loc[1]+200]
print('2')
# === 2. 构造陆地和海洋掩膜 ===
ocean = np.where(depth < -200, depth, np.nan)
land  = np.where(depth >= -200, depth, np.nan)

# === 3. 画图 ===
fig = plt.figure(figsize=(6, 4))
proj = ccrs.PlateCarree()

ax = plt.axes(projection=proj)
ax.set_extent([80, 170, -30, 65], crs=proj)

# === 4. 海洋（蓝色）===
cmap_ocean = plt.cm.GnBu_r
levels_ocean = np.linspace(-5600, -200, 32)

im_ocean = ax.contourf(
    lon, lat, ocean,
    levels=levels_ocean,
    cmap=cmap_ocean,
    transform=proj,
    extend='both'
)
print('3')
# === 5. 陆地（绿色）===
cmap_land = plt.cm.YlGn
levels_land = np.linspace(-200, 0, 32)

im_land = ax.contourf(
    lon, lat, land,
    levels=levels_land,
    cmap=cmap_land,
    transform=proj,
    extend='both'
)
print('4')
# === 6. 红色等值线（海岸线/0 m）===
# ax.contour(
#     lon, lat, depth,
#     levels=[-200],
#     colors='red',
#     linewidths=2.0,
#     linestyles='-',
#     transform=proj
# )
print('5')
# === 7. 添加地理元素 ===
ax.coastlines(resolution='10m', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='gray')
ax.add_feature(cfeature.OCEAN, facecolor='black')

# # === 8. 两个 colorbar ===

# # 海洋 colorbar
# cbar_ocean = plt.colorbar(
#     im_ocean,
#     ax=ax,
#     orientation='vertical',
#     shrink=0.8,
#     # fraction=0.046,
#     pad=0.06, 
#     ticks=[-5600, -2900, -200, ]
# )
# cbar_ocean.set_label('depth (m)')

# # 陆地 colorbar（单独放）
# cbar_land = plt.colorbar(
#     im_land,
#     ax=ax,
#     orientation='vertical',
#     shrink=0.8,
#     # fraction=0.046,
#     pad=0.04,
#     ticks=[-200, -100, -0]
# )

# ax.plot(
#     [137, 137],        # 经度固定 137°E
#     [5, 35],         # 纬度范围（与你图一致）
#     linestyle='-',    # 虚线
#     color='r',         # 颜色
#     linewidth=1.5,
#     transform=ccrs.PlateCarree()
# )
# ax.text(
#     142, 25,
#     '137°E',
#     ha='center',
#     va='bottom',
#     fontsize=8,   
#     color='red',   
#     transform=ccrs.PlateCarree()
# )

# # taiwan
# ax.plot(
#     [118, 120],      
#     [25, 23.5],        
#     linestyle='-',    
#     color='r',         
#     linewidth=1.5,
#     transform=ccrs.PlateCarree()
# )
# ax.text(
#     119, 25,
#     'Taiwan Strait',
#     ha='center',
#     va='bottom',
#     fontsize=8,   
#     color='red',   
#     transform=ccrs.PlateCarree()
# )
# # Luzon
# ax.plot(
#     [121, 121],     
#     [17, 22],        
#     linestyle='-',  
#     color='r',        
#     linewidth=1.5,
#     transform=ccrs.PlateCarree()
# )
# ax.text(
#     111, 18.5,
#     'Luzon Strait',
#     ha='center',
#     va='bottom',
#     fontsize=8,   
#     color='red',   
#     transform=ccrs.PlateCarree()
# )
# # Mindoro
# ax.plot(
#     [117, 121],        
#     [6, 13],        
#     linestyle='-',    
#     color='r',         
#     linewidth=1.5,
#     transform=ccrs.PlateCarree()
# )
# ax.text(
#     108, 9,
#     'Mindoro Strait',
#     ha='center',
#     va='bottom',
#     fontsize=8,   
#     color='red',   
#     transform=ccrs.PlateCarree()
# )
# # Karimata
# ax.plot(
#     [104, 111],     
#     [2, 2],      
#     linestyle='-',   
#     color='r',      
#     linewidth=1.5,
#     transform=ccrs.PlateCarree()
# )
# ax.text(
#     108, 3,
#     'Karimata Strait',
#     ha='center',
#     va='bottom',
#     fontsize=8,   
#     color='red',   
#     transform=ccrs.PlateCarree()
# )
# print('6')
# === 9. 坐标轴 ===
# ax.set_xticks(np.arange(100, 151, 15), crs=proj)
# ax.set_yticks(np.arange(-10, 46, 10), crs=proj)

# lon_tick = np.arange(80, lon[-1], 15)
# lat_tick = np.arange(65, -30, -10);lat_tick = lat_tick[::-1]
lon_tick = np.arange(80, lon[-1], 15)
lat_tick = np.arange(65, -30, -10);lat_tick = lat_tick[::-1]
ax.tick_params(labelsize=10)
ax.set_xticks(lon_tick)
ax.set_yticks(lat_tick)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
print('7')
plt.tight_layout()
pn = 'figures/opt.png'
plt.savefig(pn, dpi=600)
plt.close()