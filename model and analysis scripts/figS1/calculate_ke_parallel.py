#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate regional volume-integrated total kinetic energy for three experiments.

For every time file:

    KE = integral[0.5 * rho0 * (u**2 + v**2) dV]    [J]

The three experiments are processed sequentially. Within each experiment,
each worker process calculates one matching uu/vv file pair.
"""

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import xarray as xr


# ============================================================
# User settings
# ============================================================
LON_RANGE = [100, 148]
LAT_RANGE = [-15, 45]
DEP_RANGE = [0, -6000]

UPPER_DEPTH = 1000.0
RHO0 = 1025.0                 # kg m-3
EARTH_RADIUS = 6371000.0      # m
BAD_VALUE_ABS = 1.0e30

EXPERIMENT_DIRS = {
    # "eg": "../output-era-glorys/",
    # "jg": "../output-jra-glorys/",
    "jl": "../output-jra-licom3/",
}

YEARS = range(2013, 2019)
OUTPUT_DIR = "data"

# One process reads one uu/vv file pair. Each process therefore uses about
MAX_WORKERS = 600


def calculate_layer_thickness(depth, upper_depth=1000.0):
    """Return full layer thickness and thickness within 0-upper_depth."""
    z = np.abs(np.asarray(depth, dtype=np.float64))

    if z.ndim != 1 or z.size == 0:
        raise ValueError("Depth coordinate must be a non-empty 1D array")

    order = np.argsort(z)
    z_sorted = z[order]
    interfaces = np.empty(z_sorted.size + 1, dtype=np.float64)
    interfaces[0] = 0.0

    if z_sorted.size > 1:
        interfaces[1:-1] = 0.5 * (z_sorted[:-1] + z_sorted[1:])
        interfaces[-1] = (
            z_sorted[-1] + 0.5 * (z_sorted[-1] - z_sorted[-2])
        )
    else:
        interfaces[-1] = max(2.0 * z_sorted[0], upper_depth)

    thickness_sorted = np.diff(interfaces)
    upper_thickness_sorted = np.maximum(
        np.minimum(interfaces[1:], upper_depth)
        - np.maximum(interfaces[:-1], 0.0),
        0.0,
    )

    thickness = np.empty_like(thickness_sorted)
    upper_thickness = np.empty_like(upper_thickness_sorted)
    thickness[order] = thickness_sorted
    upper_thickness[order] = upper_thickness_sorted

    return thickness, upper_thickness


def nearest_index(values, target):
    return int(np.argmin(np.abs(values - target)))


def index_slice(values, value1, value2, include_end=False):
    index1 = nearest_index(values, value1)
    index2 = nearest_index(values, value2)
    start = min(index1, index2)
    end = max(index1, index2) + int(include_end)
    return slice(start, end)


def calculate_grid_area(lon, lat):
    """Calculate grid-cell areas for a regular longitude-latitude grid."""
    if lon.size < 2 or lat.size < 2:
        raise ValueError("At least two longitude and latitude points are needed")

    dlon = np.abs(np.deg2rad(np.nanmedian(np.diff(lon))))
    dlat = np.abs(np.deg2rad(np.nanmedian(np.diff(lat))))
    _, lat2d = np.meshgrid(lon, lat)

    dx = EARTH_RADIUS * np.cos(np.deg2rad(lat2d)) * dlon
    dy = EARTH_RADIUS * dlat
    return dx * dy


def read_regional_variable(ds, variable, dep_slice, lat_slice, lon_slice):
    """Read one time record over the selected lev/lat/lon region."""
    data_array = ds[variable]

    lev_dim = ds["lev"].dims[0]
    lat_dim = ds["lat"].dims[0]
    lon_dim = ds["lon"].dims[0]

    indexers = {
        lev_dim: dep_slice,
        lat_dim: lat_slice,
        lon_dim: lon_slice,
    }

    # Select the first record of any remaining dimension, normally time.
    for dim in data_array.dims:
        if dim not in indexers:
            indexers[dim] = 0

    data = data_array.isel(indexers).transpose(lev_dim, lat_dim, lon_dim)
    return np.asarray(data.values, dtype=np.float64)


def process_one_file_pair(uu_path, vv_path):
    """Calculate integrated KE for one matching uu/vv file pair."""
    filename = os.path.basename(uu_path)

    try:
        with xr.open_dataset(
            uu_path,
            decode_times=False,
            mask_and_scale=True,
            cache=False,
        ) as uds, xr.open_dataset(
            vv_path,
            decode_times=False,
            mask_and_scale=True,
            cache=False,
        ) as vds:
            lon0 = np.asarray(uds["lon"].values)
            lat0 = np.asarray(uds["lat"].values)
            depth0 = np.asarray(uds["lev"].values)

            lon_slice = index_slice(lon0, LON_RANGE[0], LON_RANGE[1])
            lat_slice = index_slice(lat0, LAT_RANGE[0], LAT_RANGE[1])
            dep_slice = index_slice(
                depth0, DEP_RANGE[0], DEP_RANGE[1], include_end=True
            )

            lon = lon0[lon_slice]
            lat = lat0[lat_slice]
            depth = depth0[dep_slice]

            uu = read_regional_variable(
                uds, "uu", dep_slice, lat_slice, lon_slice
            )
            vv = read_regional_variable(
                vds, "vv", dep_slice, lat_slice, lon_slice
            )

        expected_shape = (depth.size, lat.size, lon.size)
        if uu.shape != expected_shape or vv.shape != expected_shape:
            raise ValueError(
                f"Unexpected shapes: uu={uu.shape}, vv={vv.shape}, "
                f"expected={expected_shape}"
            )

        uu[np.abs(uu) >= BAD_VALUE_ABS] = np.nan
        vv[np.abs(vv) >= BAD_VALUE_ABS] = np.nan

        areas = calculate_grid_area(lon, lat)
        dz, dz_upper1000 = calculate_layer_thickness(depth, UPPER_DEPTH)

        full_ke_j = 0.0
        upper1000_ke_j = 0.0
        full_volume_m3 = 0.0
        upper1000_volume_m3 = 0.0

        for k in range(depth.size):
            valid = np.isfinite(uu[k]) & np.isfinite(vv[k])
            if not np.any(valid):
                continue

            # Kinetic-energy density, J m-3.
            ke_density = 0.5 * RHO0 * (
                uu[k, valid] ** 2 + vv[k, valid] ** 2
            )
            valid_area = areas[valid]

            full_cell_volume = valid_area * dz[k]
            full_ke_j += np.sum(ke_density * full_cell_volume)
            full_volume_m3 += np.sum(full_cell_volume)

            if dz_upper1000[k] > 0.0:
                upper_cell_volume = valid_area * dz_upper1000[k]
                upper1000_ke_j += np.sum(ke_density * upper_cell_volume)
                upper1000_volume_m3 += np.sum(upper_cell_volume)

        full_mean_ke_j_m3 = (
            full_ke_j / full_volume_m3 if full_volume_m3 > 0.0 else np.nan
        )
        upper1000_mean_ke_j_m3 = (
            upper1000_ke_j / upper1000_volume_m3
            if upper1000_volume_m3 > 0.0
            else np.nan
        )

        return {
            "filename": filename,
            "success": True,
            "full_ke_j": full_ke_j,
            "upper1000_ke_j": upper1000_ke_j,
            "full_mean_ke_j_m3": full_mean_ke_j_m3,
            "upper1000_mean_ke_j_m3": upper1000_mean_ke_j_m3,
            "error": "",
        }

    except Exception as error:
        return {
            "filename": filename,
            "success": False,
            "full_ke_j": np.nan,
            "upper1000_ke_j": np.nan,
            "full_mean_ke_j_m3": np.nan,
            "upper1000_mean_ke_j_m3": np.nan,
            "error": repr(error),
        }


def find_file_pairs(input_dir):
    """Find uu files for 2013-2018 and pair them with corresponding vv files."""
    filenames = os.listdir(input_dir)
    valid_year_keys = [key for year in YEARS for key in (f"daym-uu{year}", f"daym-uu-{year}")]

    uu_files = sorted(
        name
        for name in filenames
        if any(key in name for key in valid_year_keys)
    )

    pairs = []
    missing_vv = []

    for uu_name in uu_files:
        vv_name = uu_name.replace("uu", "vv", 1)
        uu_path = os.path.join(input_dir, uu_name)
        vv_path = os.path.join(input_dir, vv_name)

        if os.path.isfile(vv_path):
            pairs.append((uu_path, vv_path))
        else:
            missing_vv.append((uu_name, vv_name))

    return pairs, missing_vv


def save_results(results, pnstr):
    """Save ordered numeric time series and a traceable combined table."""
    output_fields = {
        f"regional-integrated-total-ke-{pnstr}.txt": "full_ke_j",
        f"regional-integrated-total-ke-upper1000-{pnstr}.txt": "upper1000_ke_j",
        f"regional-mean-total-ke-density-{pnstr}.txt": "full_mean_ke_j_m3",
        f"regional-mean-total-ke-density-upper1000-{pnstr}.txt": (
            "upper1000_mean_ke_j_m3"
        ),
    }

    for output_name, field in output_fields.items():
        values = np.asarray([result[field] for result in results])
        np.savetxt(os.path.join(OUTPUT_DIR, output_name), values, fmt="%.10e")

    table_path = os.path.join(OUTPUT_DIR, f"regional-total-ke-{pnstr}.dat")
    with open(table_path, "w", encoding="utf-8") as output_file:
        output_file.write(
            "filename full_ke_J upper1000_ke_J "
            "full_mean_ke_J_m3 upper1000_mean_ke_J_m3\n"
        )
        for result in results:
            output_file.write(
                f"{result['filename']} "
                f"{result['full_ke_j']:.10e} "
                f"{result['upper1000_ke_j']:.10e} "
                f"{result['full_mean_ke_j_m3']:.10e} "
                f"{result['upper1000_mean_ke_j_m3']:.10e}\n"
            )

    return table_path


def run_one_experiment(pnstr, input_dir):
    input_dir = os.path.abspath(input_dir)

    print("\n" + "=" * 72)
    print(f"Starting experiment: {pnstr}")
    print(f"Input directory: {input_dir}")

    if not os.path.isdir(input_dir):
        print(f"SKIPPED: directory does not exist: {input_dir}")
        return False

    pairs, missing_vv = find_file_pairs(input_dir)

    for uu_name, vv_name in missing_vv:
        print(f"Missing pair: {uu_name} -> {vv_name}")

    if not pairs:
        print("SKIPPED: no complete uu/vv file pairs were found")
        return False

    worker_number = min(MAX_WORKERS, len(pairs))
    print(f"Number of complete file pairs: {len(pairs)}")
    print(f"Number of processes: {worker_number}")

    results_by_name = {}
    mp_context = mp.get_context("spawn")

    with ProcessPoolExecutor(
        max_workers=worker_number,
        mp_context=mp_context,
    ) as executor:
        future_to_pair = {
            executor.submit(process_one_file_pair, uu_path, vv_path): (
                uu_path,
                vv_path,
            )
            for uu_path, vv_path in pairs
        }

        for finished, future in enumerate(as_completed(future_to_pair), 1):
            uu_path, _ = future_to_pair[future]
            filename = os.path.basename(uu_path)

            try:
                result = future.result()
            except Exception as error:
                result = {
                    "filename": filename,
                    "success": False,
                    "full_ke_j": np.nan,
                    "upper1000_ke_j": np.nan,
                    "full_mean_ke_j_m3": np.nan,
                    "upper1000_mean_ke_j_m3": np.nan,
                    "error": repr(error),
                }

            results_by_name[filename] = result

            if result["success"]:
                print(
                    f"[{finished:4d}/{len(pairs):4d}] {filename}: "
                    f"KE={result['full_ke_j']:.6e} J, "
                    f"upper1000={result['upper1000_ke_j']:.6e} J"
                )
            else:
                print(
                    f"[{finished:4d}/{len(pairs):4d}] FAILED "
                    f"{filename}: {result['error']}"
                )

    ordered_results = [
        results_by_name[os.path.basename(uu_path)] for uu_path, _ in pairs
    ]
    table_path = save_results(ordered_results, pnstr)
    failed_number = sum(not result["success"] for result in ordered_results)

    print(f"Completed experiment: {pnstr}")
    print(f"Failed files: {failed_number}")
    print(f"Combined output: {table_path}")

    return failed_number == 0 and not missing_vv


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    experiment_status = []

    # Experiments run sequentially; file pairs within each experiment run in parallel.
    for pnstr, input_dir in EXPERIMENT_DIRS.items():
        success = run_one_experiment(pnstr, input_dir)
        experiment_status.append((pnstr, success))

    print("\n" + "=" * 72)
    print("All requested experiments finished:")
    for pnstr, success in experiment_status:
        status = "OK" if success else "FAILED/SKIPPED/INCOMPLETE"
        print(f"  {pnstr}: {status}")


if __name__ == "__main__":
    mp.freeze_support()
    main()