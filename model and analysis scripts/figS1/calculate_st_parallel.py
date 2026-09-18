#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate regional means for three experiments and two variables.

The six cases (tt/ss x eg/jg/jl) are processed sequentially. Within each
case, every worker process calculates one NetCDF file.
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
DEP_RANGE = [0, 6000]

UPPER_DEPTH = 1000.0
EARTH_RADIUS = 6371000.0
BAD_VALUE_ABS = 1.0e30

VARIABLES = ["tt"]

# EXPERIMENT_DIRS = {
#     "eg": "../output-era-glorys/",
#     "jg": "../output-jra-glorys/",
#     "jl": "../output-jra-licom3/",
# }
EXPERIMENT_DIRS = {
    "eg": "../licomreal-TH/ts-eg/",
#     "jg": "../licomreal-TH/ts-jg/",
#     "jl": "../licomreal-TH/ts-jl/",
}

YEARS = range(2014, 2019)
OUTPUT_DIR = "data"

# Each process calculates one file. Reduce this value if memory is insufficient.
MAX_WORKERS = 60


# ============================================================
# Utility functions
# ============================================================
def calculate_layer_thickness(depth, upper_depth=1000.0):
    """Calculate full layer thickness and thickness within 0-upper_depth."""
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
    """Calculate grid-cell area for a regular lon-lat grid."""
    if lon.size < 2 or lat.size < 2:
        raise ValueError("At least two longitude and latitude points are needed")

    dlon = np.abs(np.deg2rad(np.nanmedian(np.diff(lon))))
    dlat = np.abs(np.deg2rad(np.nanmedian(np.diff(lat))))
    _, lat2d = np.meshgrid(lon, lat)

    dx = EARTH_RADIUS * np.cos(np.deg2rad(lat2d)) * dlon
    dy = EARTH_RADIUS * dlat
    return dx * dy


def get_variable_subset(data_array, dep_slice, lat_slice, lon_slice):
    """Read only the requested region, assuming time/lev/lat/lon order."""
    dims = data_array.dims

    if data_array.ndim != 4:
        raise ValueError(
            f"Variable {data_array.name} must be 4D, but dimensions are {dims}"
        )

    subset = data_array.isel(
        {
            dims[0]: 0,
            dims[1]: dep_slice,
            dims[2]: lat_slice,
            dims[3]: lon_slice,
        }
    )
    return np.asarray(subset.values, dtype=np.float64)


# ============================================================
# Worker: one process calculates one file
# ============================================================
def process_one_file(filepath, vn):
    """Calculate surface, full-depth and upper-1000 m means for one file."""
    filename = os.path.basename(filepath)

    try:
        with xr.open_dataset(
            filepath,
            decode_times=False,
            mask_and_scale=True,
            cache=False,
        ) as ds:
            lon0 = np.asarray(ds["lon"].values)
            lat0 = np.asarray(ds["lat"].values)
            depth0 = np.asarray(ds["lev1"].values)

            lon_slice = index_slice(lon0, LON_RANGE[0], LON_RANGE[1])
            lat_slice = index_slice(lat0, LAT_RANGE[0], LAT_RANGE[1])
            dep_slice = index_slice(
                depth0, DEP_RANGE[0], DEP_RANGE[1], include_end=True
            )

            lon = lon0[lon_slice]
            lat = lat0[lat_slice]
            depth = depth0[dep_slice]
            data = get_variable_subset(
                ds[vn], dep_slice, lat_slice, lon_slice
            )

        expected_shape = (depth.size, lat.size, lon.size)
        if data.shape != expected_shape:
            raise ValueError(
                f"Unexpected subset shape {data.shape}; expected {expected_shape}"
            )

        data[np.abs(data) >= BAD_VALUE_ABS] = np.nan
        areas = calculate_grid_area(lon, lat)
        dz, dz_upper1000 = calculate_layer_thickness(depth, UPPER_DEPTH)

        surface = data[0]
        surface_valid = np.isfinite(surface)
        surface_mean = (
            np.sum(surface[surface_valid] * areas[surface_valid])
            / np.sum(areas[surface_valid])
            if np.any(surface_valid)
            else np.nan
        )

        full_value_sum = 0.0
        full_volume_sum = 0.0
        upper_value_sum = 0.0
        upper_volume_sum = 0.0

        for k, var in enumerate(data):
            valid = np.isfinite(var)
            if not np.any(valid):
                continue

            valid_area = areas[valid]
            full_volume = valid_area * dz[k]
            full_value_sum += np.sum(var[valid] * full_volume)
            full_volume_sum += np.sum(full_volume)

            if dz_upper1000[k] > 0.0:
                upper_volume = valid_area * dz_upper1000[k]
                upper_value_sum += np.sum(var[valid] * upper_volume)
                upper_volume_sum += np.sum(upper_volume)

        full_mean = (
            full_value_sum / full_volume_sum
            if full_volume_sum > 0.0
            else np.nan
        )
        upper1000_mean = (
            upper_value_sum / upper_volume_sum
            if upper_volume_sum > 0.0
            else np.nan
        )

        return {
            "filename": filename,
            "success": True,
            "surface_mean": surface_mean,
            "full_mean": full_mean,
            "upper1000_mean": upper1000_mean,
            "error": "",
        }

    except Exception as error:
        return {
            "filename": filename,
            "success": False,
            "surface_mean": np.nan,
            "full_mean": np.nan,
            "upper1000_mean": np.nan,
            "error": repr(error),
        }


# ============================================================
# Run one variable for one experiment
# ============================================================
def run_one_case(vn, pnstr, input_dir):
    input_dir = os.path.abspath(input_dir)
    key_strs = [f"{vn}-{year}" for year in YEARS]

    print("\n" + "=" * 70)
    print(f"Starting case: variable={vn}, experiment={pnstr}")
    print(f"Input directory: {input_dir}")

    if not os.path.isdir(input_dir):
        print(f"SKIPPED: input directory does not exist: {input_dir}")
        return False

    files = sorted(
        os.path.join(input_dir, item)
        for item in os.listdir(input_dir)
        if any(key in item for key in key_strs)
    )

    if not files:
        print(f"SKIPPED: no {vn} files were found")
        return False

    worker_number = min(MAX_WORKERS, len(files))
    print(f"Number of files: {len(files)}")
    print(f"Number of processes: {worker_number}")

    results_by_name = {}
    mp_context = mp.get_context("spawn")

    with ProcessPoolExecutor(
        max_workers=worker_number,
        mp_context=mp_context,
    ) as executor:
        future_to_file = {
            executor.submit(process_one_file, filepath, vn): filepath
            for filepath in files
        }

        for finished, future in enumerate(as_completed(future_to_file), 1):
            filepath = future_to_file[future]
            filename = os.path.basename(filepath)

            try:
                result = future.result()
            except Exception as error:
                result = {
                    "filename": filename,
                    "success": False,
                    "surface_mean": np.nan,
                    "full_mean": np.nan,
                    "upper1000_mean": np.nan,
                    "error": repr(error),
                }

            results_by_name[filename] = result

            if result["success"]:
                print(
                    f"[{finished:4d}/{len(files):4d}] {filename}: "
                    f"surface={result['surface_mean']:.6f}, "
                    f"full={result['full_mean']:.6f}, "
                    f"upper1000={result['upper1000_mean']:.6f}"
                )
            else:
                print(
                    f"[{finished:4d}/{len(files):4d}] FAILED "
                    f"{filename}: {result['error']}"
                )

    ordered_results = [
        results_by_name[os.path.basename(filepath)] for filepath in files
    ]
    failed_results = [r for r in ordered_results if not r["success"]]

    surface_means = np.asarray(
        [r["surface_mean"] for r in ordered_results], dtype=np.float64
    )
    full_means = np.asarray(
        [r["full_mean"] for r in ordered_results], dtype=np.float64
    )
    upper1000_means = np.asarray(
        [r["upper1000_mean"] for r in ordered_results], dtype=np.float64
    )

    if vn == "ss":
        output_names = (
            f"average-sss-{pnstr}-TH.txt",
            f"average-ss-{pnstr}-TH.txt",
            f"average-ss-upper1000-{pnstr}-TH.txt",
        )
    else:
        output_names = (
            f"average-sst-{pnstr}-TH.txt",
            f"average-st-{pnstr}-TH.txt",
            f"average-st-upper1000-{pnstr}-TH.txt",
        )

    output_paths = [os.path.join(OUTPUT_DIR, name) for name in output_names]
    np.savetxt(output_paths[0], surface_means, fmt="%.6f")
    np.savetxt(output_paths[1], full_means, fmt="%.6f")
    np.savetxt(output_paths[2], upper1000_means, fmt="%.6f")

    print(f"Completed case: variable={vn}, experiment={pnstr}")
    print(f"Failed files: {len(failed_results)}")
    for output_path in output_paths:
        print(f"Saved: {output_path}")

    return len(failed_results) == 0


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    case_status = []

    # Six cases are run sequentially; files within each case run in parallel.
    for vn in VARIABLES:
        for pnstr, input_dir in EXPERIMENT_DIRS.items():
            success = run_one_case(vn, pnstr, input_dir)
            case_status.append((vn, pnstr, success))

    print("\n" + "=" * 70)
    print("All requested cases finished:")
    for vn, pnstr, success in case_status:
        status = "OK" if success else "FAILED/SKIPPED"
        print(f"  {pnstr} {vn}: {status}")


if __name__ == "__main__":
    mp.freeze_support()
    main()