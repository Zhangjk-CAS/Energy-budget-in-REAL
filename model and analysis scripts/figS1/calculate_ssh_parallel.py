#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate area-weighted regional mean sea surface height (SSH).

The experiments are processed sequentially. Within each experiment, every
worker process calculates one SSH NetCDF file. Both 2D (lat, lon) and arrays
with additional singleton dimensions such as (time, lat, lon) are supported.
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

EARTH_RADIUS = 6371000.0
BAD_VALUE_ABS = 1.0e30

# The program searches these variable names in order. Modify this tuple if the
# SSH variable in the files has another name.
SSH_VARIABLE_CANDIDATES = ("ssh", "z0", "h0")

# Supported filename forms include ssh2014..., ssh-2014..., z02014... and
# z0-2014.... A prefix can also occur after another prefix, e.g. daym-z0-2014.
SSH_FILE_PREFIXES = ("ssh", "z0")

EXPERIMENT_DIRS = {
    "eg": "../output-era-glorys/",
    "jg": "../output-jra-glorys/",
    "jl": "../output-jra-licom3/",
}

YEARS = range(2013, 2019)
OUTPUT_DIR = "data"

# SSH is only two-dimensional, so its per-process memory use is much smaller
# than that of the temperature/salinity volume-mean calculation.
MAX_WORKERS = 60


def nearest_index(values, target):
    return int(np.argmin(np.abs(values - target)))


def index_slice(values, value1, value2, include_end=False):
    """Return an index slice regardless of coordinate storage direction."""
    index1 = nearest_index(values, value1)
    index2 = nearest_index(values, value2)
    start = min(index1, index2)
    end = max(index1, index2) + int(include_end)
    return slice(start, end)


def calculate_grid_area(lon, lat):
    """Calculate grid-cell areas for a regular longitude-latitude grid."""
    if lon.ndim != 1 or lat.ndim != 1:
        raise ValueError("This program requires one-dimensional lon and lat")
    if lon.size < 2 or lat.size < 2:
        raise ValueError("At least two longitude and latitude points are needed")

    dlon = np.abs(np.deg2rad(np.nanmedian(np.diff(lon))))
    dlat = np.abs(np.deg2rad(np.nanmedian(np.diff(lat))))
    _, lat2d = np.meshgrid(lon, lat)

    dx = EARTH_RADIUS * np.cos(np.deg2rad(lat2d)) * dlon
    dy = EARTH_RADIUS * dlat
    return dx * dy


def find_ssh_variable(ds):
    """Return the first configured SSH variable found in a dataset."""
    for variable in SSH_VARIABLE_CANDIDATES:
        if variable in ds.variables:
            return variable

    raise KeyError(
        "None of the SSH variables "
        f"{SSH_VARIABLE_CANDIDATES} was found. Available data variables: "
        f"{list(ds.data_vars)}"
    )


def read_ssh_subset(ds, variable, lat_slice, lon_slice):
    """Read one SSH record over the selected horizontal region."""
    data_array = ds[variable]
    lat_dim = ds["lat"].dims[0]
    lon_dim = ds["lon"].dims[0]

    if lat_dim not in data_array.dims or lon_dim not in data_array.dims:
        raise ValueError(
            f"Variable {variable} dimensions {data_array.dims} do not contain "
            f"latitude/longitude dimensions {lat_dim}/{lon_dim}"
        )

    indexers = {
        lat_dim: lat_slice,
        lon_dim: lon_slice,
    }

    # Select the first record from time or other extra dimensions.
    for dim in data_array.dims:
        if dim not in indexers:
            indexers[dim] = 0

    subset = data_array.isel(indexers).transpose(lat_dim, lon_dim)
    return np.asarray(subset.values, dtype=np.float64)


def process_one_file(filepath):
    """Calculate the area-weighted regional mean SSH for one file."""
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

            lon_slice = index_slice(lon0, LON_RANGE[0], LON_RANGE[1])
            lat_slice = index_slice(lat0, LAT_RANGE[0], LAT_RANGE[1])

            lon = lon0[lon_slice]
            lat = lat0[lat_slice]
            variable = find_ssh_variable(ds)
            ssh = read_ssh_subset(ds, variable, lat_slice, lon_slice)

        expected_shape = (lat.size, lon.size)
        if ssh.shape != expected_shape:
            raise ValueError(
                f"Unexpected SSH shape {ssh.shape}; expected {expected_shape}"
            )

        ssh[np.abs(ssh) >= BAD_VALUE_ABS] = np.nan
        areas = calculate_grid_area(lon, lat)
        valid = np.isfinite(ssh)

        if not np.any(valid):
            regional_mean_ssh = np.nan
        else:
            regional_mean_ssh = (
                np.sum(ssh[valid] * areas[valid])
                / np.sum(areas[valid])
            )

        return {
            "filename": filename,
            "variable": variable,
            "success": True,
            "regional_mean_ssh": regional_mean_ssh,
            "error": "",
        }

    except Exception as error:
        return {
            "filename": filename,
            "variable": "",
            "success": False,
            "regional_mean_ssh": np.nan,
            "error": repr(error),
        }


def is_matching_ssh_file(entry):
    """Check file type and supported SSH/year filename patterns."""
    if not entry.is_file():
        return False

    name = entry.name
    keys = [
        key
        for year in YEARS
        for prefix in SSH_FILE_PREFIXES
        for key in (f"daym-{prefix}{year}", f"daym-{prefix}-{year}")
    ]
    return any(key in name for key in keys)


def run_one_experiment(pnstr, input_dir):
    """Run all SSH files for one experiment and save ordered results."""
    input_dir = os.path.abspath(input_dir)

    print("\n" + "=" * 70)
    print(f"Starting SSH calculation: experiment={pnstr}")
    print(f"Input directory: {input_dir}")

    if not os.path.isdir(input_dir):
        print(f"SKIPPED: input directory does not exist: {input_dir}")
        return False

    files = sorted(
        os.path.join(input_dir, entry.name)
        for entry in os.scandir(input_dir)
        if is_matching_ssh_file(entry)
    )

    if not files:
        print("SKIPPED: no matching SSH files were found")
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
            executor.submit(process_one_file, filepath): filepath
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
                    "variable": "",
                    "success": False,
                    "regional_mean_ssh": np.nan,
                    "error": repr(error),
                }

            results_by_name[filename] = result

            if result["success"]:
                print(
                    f"[{finished:4d}/{len(files):4d}] {filename}: "
                    f"{result['variable']}={result['regional_mean_ssh']:.8f} m"
                )
            else:
                print(
                    f"[{finished:4d}/{len(files):4d}] FAILED "
                    f"{filename}: {result['error']}"
                )

    ordered_results = [
        results_by_name[os.path.basename(filepath)] for filepath in files
    ]
    failed_results = [result for result in ordered_results if not result["success"]]
    means = np.asarray(
        [result["regional_mean_ssh"] for result in ordered_results],
        dtype=np.float64,
    )

    output_txt = os.path.join(OUTPUT_DIR, f"average-ssh-{pnstr}.txt")
    np.savetxt(output_txt, means, fmt="%.10e")

    output_dat = os.path.join(OUTPUT_DIR, f"average-ssh-{pnstr}.dat")
    with open(output_dat, "w", encoding="utf-8") as output_file:
        output_file.write("filename variable regional_mean_ssh_m\n")
        for result in ordered_results:
            output_file.write(
                f"{result['filename']} {result['variable']} "
                f"{result['regional_mean_ssh']:.10e}\n"
            )

    print(f"Completed experiment: {pnstr}")
    print(f"Failed files: {len(failed_results)}")
    print(f"Saved: {output_txt}")
    print(f"Saved: {output_dat}")

    return len(failed_results) == 0


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    experiment_status = []

    for pnstr, input_dir in EXPERIMENT_DIRS.items():
        success = run_one_experiment(pnstr, input_dir)
        experiment_status.append((pnstr, success))

    print("\n" + "=" * 70)
    print("All requested SSH calculations finished:")
    for pnstr, success in experiment_status:
        status = "OK" if success else "FAILED/SKIPPED"
        print(f"  {pnstr}: {status}")


if __name__ == "__main__":
    mp.freeze_support()
    main()