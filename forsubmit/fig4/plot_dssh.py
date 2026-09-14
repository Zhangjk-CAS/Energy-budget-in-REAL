import os
import re
import csv
from datetime import datetime, timedelta

import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


source_names = ['Reg-EG', 'Reg-JG', 'Reg-JL', 'GLORYS', 'AVISO']

path_list = [
    '../output-FaD-era/',
    '../output-jra-glorys/',
    '../output-FaD-jra/',
    '../../mercatordata/mercator_compare/RegionalData/',
    '../../observation/AVISO/',
]


key_strs = ['z0', 'z0', 'z0', 'z0', '_phy_l4_']
vn_list = ['ssh', 'ssh', 'ssh', 'ssh', 'ssh']

lon_range = [137.0]

lat_ranges_necc = [
    [2.0, 7.5], [2.0, 7.5], [2.0, 7.5], [2.0, 7.5], [2.0, 7.5],
]
lat_ranges_kuro = [
    [29.0, 31.0], [29.0, 34.6], [29.0, 31.0], [31.0, 32.0], [31.0, 32.0],
]

start_date = datetime(2013, 1, 1)
end_date = datetime(2018, 12, 31)
analysis_start_date = datetime(2014, 1, 1)

output_data_dir = 'data'
output_figure_dir = 'figures'
figure_name = 'transport-ncckuro-errorbar-all.png'
os.makedirs(output_data_dir, exist_ok=True)
os.makedirs(output_figure_dir, exist_ok=True)

text_names = [
    'ncc-dz-era.txt',
    'ncc-dz-jra-glorys.txt',
    'ncc-dz-jra-licom3.txt',
    'ncc-dz-glorys.txt',
    'ncc-dz-aviso.txt',
    'kuro-dz-era-glorys.txt',
    'kuro-dz-jra-glorys.txt',
    'kuro-dz-jra-licom3.txt',
    'kuro-dz-glorys0124.txt',
    'kuro-dz-aviso0124.txt',
]


number_days = (end_date - start_date).days + 1
date_list = [start_date + timedelta(days=n) for n in range(number_days)]
date_to_index = {day.date(): n for n, day in enumerate(date_list)}

dz = np.full((number_days, 10), np.nan)
already_read = np.zeros((number_days, 5), dtype=bool)


for source_index in range(5):
    path = path_list[source_index]
    if not os.path.isdir(path):
        raise FileNotFoundError('Input directory not found for %s: %s' %
                                (source_names[source_index], path))

    nowdir = os.path.abspath(path)
    files = sorted(os.listdir(nowdir))
    matched_files = 0

    for item in files:
        if key_strs[source_index] not in item:
            continue
        if not item.lower().endswith(('.nc', '.nc4')):
            continue

        fn = os.path.join(nowdir, item)
        print('reading netcdf ' + fn + ' -' + vn_list[source_index])
        matched_files += 1

        with nc.Dataset(fn) as ncdata:
            vn = vn_list[source_index]
            if vn not in ncdata.variables:
                if source_names[source_index] == 'AVISO' and 'adt' in ncdata.variables:
                    vn = 'adt'
                    print('  AVISO has no ssh variable; using adt instead')
                else:
                    raise KeyError('%s does not contain variable %s' % (fn, vn))

            var0 = ncdata.variables[vn]  
            if var0.ndim != 3:
                raise ValueError('%s: %s is not a 3-D (time, lat, lon) variable' % (fn, vn))

            lon_name = 'lon' if 'lon' in ncdata.variables else 'longitude'
            lat_name = 'lat' if 'lat' in ncdata.variables else 'latitude'
            lon0 = ncdata.variables[lon_name][:]
            lat0 = ncdata.variables[lat_name][:]
            if lon0.ndim != 1 or lat0.ndim != 1:
                raise ValueError('%s: longitude and latitude must be 1-D coordinates' % fn)

            lon_diff = (np.asarray(lon0) - lon_range[0] + 180.0) % 360.0 - 180.0
            lon_loc = int(np.argmin(np.abs(lon_diff)))

            lat_range = lat_ranges_necc[source_index]
            necc_south = int(np.argmin(abs(lat0 - lat_range[0])))
            necc_north = int(np.argmin(abs(lat0 - lat_range[1])))

            lat_range = lat_ranges_kuro[source_index]
            kuro_south = int(np.argmin(abs(lat0 - lat_range[0])))
            kuro_north = int(np.argmin(abs(lat0 - lat_range[1])))

            time_name = None
            for dim_name in var0.dimensions:
                if 'time' in dim_name.lower() and dim_name in ncdata.variables:
                    time_name = dim_name
                    break

            if time_name is not None and hasattr(ncdata.variables[time_name], 'units'):
                time_var = ncdata.variables[time_name]
                raw_dates = nc.num2date(time_var[:], time_var.units,
                                       calendar=getattr(time_var, 'calendar', 'standard'))
                file_dates = [datetime(int(t.year), int(t.month), int(t.day)).date()
                              for t in raw_dates]
                if len(file_dates) != var0.shape[0]:
                    raise ValueError('%s: time and SSH dimensions have different lengths' % fn)
            else:
                match = re.search(r'(?<!\d)(20\d{2})[-_]?([01]\d)[-_]?([0-3]\d)(?!\d)', item)
                if match is None or var0.shape[0] != 1:
                    raise ValueError('%s: cannot determine dates; use CF time or daily filenames containing YYYYMMDD' % fn)
                try:
                    file_dates = [datetime(*map(int, match.groups())).date()]
                except ValueError as exc:
                    raise ValueError('%s: invalid date in filename' % fn) from exc

            for time_index, file_date in enumerate(file_dates):
                if file_date not in date_to_index:
                    continue
                day_index = date_to_index[file_date]
                if already_read[day_index, source_index]:
                    raise ValueError('%s has a duplicate daily record for %s' %
                                     (source_names[source_index], file_date))
                already_read[day_index, source_index] = True

                necc_south_ssh = float(np.ma.filled(
                    var0[time_index, necc_south, lon_loc], np.nan))
                necc_north_ssh = float(np.ma.filled(
                    var0[time_index, necc_north, lon_loc], np.nan))
                kuro_south_ssh = float(np.ma.filled(
                    var0[time_index, kuro_south, lon_loc], np.nan))
                kuro_north_ssh = float(np.ma.filled(
                    var0[time_index, kuro_north, lon_loc], np.nan))

                dz[day_index, source_index] = necc_south_ssh - necc_north_ssh
                dz[day_index, source_index + 5] = kuro_south_ssh - kuro_north_ssh

    if matched_files == 0:
        raise FileNotFoundError('%s: no NetCDF files in %s matching %s' %
                                (source_names[source_index], nowdir,
                                 key_strs[source_index]))

    valid_necc = np.count_nonzero(np.isfinite(dz[:, source_index]))
    valid_kuro = np.count_nonzero(np.isfinite(dz[:, source_index + 5]))
    print('%s: NECC %d/%d valid days, Kuroshio %d/%d valid days' %
          (source_names[source_index], valid_necc, number_days,
           valid_kuro, number_days))


for nn in range(10):
    fnot = os.path.join(output_data_dir, text_names[nn])
    np.savetxt(fnot, dz[:, nn], fmt='%.6f')
    print('saved ' + fnot)

csv_path = os.path.join(output_data_dir, 'dssh_daily_2013_2018.csv')
with open(csv_path, 'w', newline='', encoding='utf-8') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['date'] + [name[:-4] for name in text_names])
    for day_index, day in enumerate(date_list):
        writer.writerow([day.strftime('%Y-%m-%d')] +
                        ['%.6f' % value for value in dz[day_index, :]])
print('saved ' + csv_path)


analysis_loc = [nn for nn, day in enumerate(date_list)
                if day >= analysis_start_date]
dz_analysis = dz[analysis_loc, :]

means = np.full(10, np.nan)
stds = np.full(10, np.nan)
corr = np.full(10, np.nan)

for nn in range(10):
    valid = np.isfinite(dz_analysis[:, nn])
    if not np.any(valid):
        raise ValueError('%s has no valid data during the analysis period' % text_names[nn])
    means[nn] = np.mean(dz_analysis[valid, nn])
    stds[nn] = np.std(dz_analysis[valid, nn])

    reference_index = 4 if nn < 5 else 9
    paired = valid & np.isfinite(dz_analysis[:, reference_index])
    if np.count_nonzero(paired) >= 2:
        a = dz_analysis[paired, nn]
        b = dz_analysis[paired, reference_index]
        if np.std(a) > 0 and np.std(b) > 0:
            corr[nn] = np.corrcoef(a, b)[0, 1]

    print('%-27s mean=%8.4f m  std=%8.4f m  corr=%6.3f  valid days=%d' %
          (text_names[nn], means[nn], stds[nn], corr[nn], np.count_nonzero(valid)))


x = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11]
x_labels = [3, 9]
labels = ['North Equatorial Countercurrent', 'Kuroshio']
colors = ['blue', 'red', 'green', 'blue', 'green',
          'blue', 'red', 'green', 'blue', 'green']
lnstr = ['o', 'o', 'o', '<', '<',
         'o', 'o', 'o', '<', '<']
labelstr = ['Reg-EG', 'Reg-JG', 'Reg-JL', 'GLORYS', 'AVISO']

fig, ax = plt.subplots(1, 1, figsize=(6, 3))
for nn in range(10):
    if nn < 5:
        ax.errorbar(x[nn], means[nn], yerr=stds[nn], fmt=lnstr[nn],
                    color=colors[nn], label=labelstr[nn],
                    capsize=4, elinewidth=2, markersize=5)
    else:
        ax.errorbar(x[nn], means[nn], yerr=stds[nn], fmt=lnstr[nn],
                    color=colors[nn], capsize=4,
                    elinewidth=2, markersize=5)

ax.grid(True, linestyle='--', alpha=0.4)
ax.set_ylabel('dSSH (m)', fontsize=10)
ax.set_xticks(x_labels)
ax.set_xticklabels(labels, ha='center', fontsize=11)
ax.set_xlim(-0.5, 12)
ax.legend(loc='upper left', ncol=2, frameon=False)

ax2 = ax.twinx()
ax2.plot(x[0:5], corr[0:5], marker='o', linestyle='--')
ax2.plot(x[5:10], corr[5:10], marker='o', linestyle='--')
ax2.set_ylim(-0.04, 1.84)
ax2.set_yticks([0, 0.3, 0.6, 0.9])
ax2.set_ylabel('correlation')

pn = os.path.join(output_figure_dir, figure_name)
plt.savefig(pn, dpi=600)
plt.close(fig)
print('saved ' + pn)
