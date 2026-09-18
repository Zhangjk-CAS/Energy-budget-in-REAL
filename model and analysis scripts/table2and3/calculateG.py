#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
import math
import sys
import numpy as np
import netCDF4 as nc

from gsw_alternative import gsw_rho, gsw_SA_from_SP, gsw_CT_from_t


# ============================================================
# 1. USER SETTINGS
# ============================================================

START_YEAR = 2014
END_YEAR = 2018

# ------------------------------------------------------------
# Calculation region
# ------------------------------------------------------------
# Grid-point centers inside these bounds are retained, including boundaries.
LON_RANGE = (100.0, 148.0)
LAT_RANGE = (-15.0, 45.0)

# All files were shown in one directory. Change these if needed.
DAILY_ROOT = Path("/lustre/gmmcfs/zhangjk/licomreal/output-jra-licom3/fort.22")
MEAN_ROOT = Path("/lustre/gmmcfs/zhangjk/licomreal/output-jra-licom3/fort.22/ave")

OUTPUT_NC = Path(
    f"G-potdens-rhoN0-{START_YEAR}-{END_YEAR}-JL-"
    f"{LON_RANGE[0]:g}_{LON_RANGE[1]:g}E-"
    f"{abs(LAT_RANGE[0]):g}S_{LAT_RANGE[1]:g}N.nc"
)
OUTPUT_DAT = Path(
    f"G-potdens-rhoN0-{START_YEAR}-{END_YEAR}-JL-"
    f"{LON_RANGE[0]:g}_{LON_RANGE[1]:g}E-"
    f"{abs(LAT_RANGE[0]):g}S_{LAT_RANGE[1]:g}N.dat"
)

# File/variable names.
VAR_NAMES = {
    "fr": "fr",
    "lw": "lw",
    "sw": "sw",
    "sh": "sh",
    "lh": "lh",
    "rho": "rho",
    "tx": "tx",
    "ty": "ty",
    "uu": "uu",
    "vv": "vv",
    "tt": "tt",
    "ss": "ss",
    "z0": "ssh",
}

# Daily and mean filename templates.
DAILY_TEMPLATE = "{var}-{year:04d}-{month:02d}-{day:02d}.nc"
MEAN_TEMPLATE = "{var}-" + f"{START_YEAR}-{END_YEAR}" + ".nc"

# 365-day model calendar.
DAYS_IN_MONTH = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Missing files:
# False -> stop immediately if any required daily file is absent.
# True  -> skip the whole day and print it.
SKIP_MISSING_DAYS = False

# ------------------------------------------------------------
# Heat-flux sign convention
# ------------------------------------------------------------
# H must be positive INTO THE OCEAN.
#
# The ncdump metadata only give W/m^2, not positive direction.
# Default below means "use the signs already stored in the files".
#
# If your archive stores:
#   sw > 0 downward into ocean
#   lw, sh, lh > 0 upward out of ocean
# then change to:
#   {"sw": +1.0, "lw": -1.0, "sh": -1.0, "lh": -1.0}
#
HEAT_FLUX_SIGNS = {
    "sw": +1.0,
    "lw": +1.0,
    "sh": +1.0,
    "lh": +1.0,
}


# ------------------------------------------------------------
# Heat-flux diagnostics
# ------------------------------------------------------------
# These checks do NOT change the heat flux used in the LEC calculation.
# They only print metadata/statistics and compare two common conventions:
#
#   candidate A: H_all_stored = sw + lw + sh + lh
#       Appropriate if every component is already stored with its physical sign
#       relative to the ocean (e.g. upward losses are negative).
#
#   candidate B: H_sw_minus_outgoing = sw - lw - sh - lh
#       Appropriate if sw is positive downward while lw/sh/lh are stored as
#       positive upward loss magnitudes.
#
# The configured heat flux actually used by the calculation is still:
#       H = sum(HEAT_FLUX_SIGNS[name] * component[name])
#
PRINT_HEAT_FILE_METADATA = True
PRINT_HEAT_FIRST_DAY_DIAGNOSTICS = True
PRINT_HEAT_PERIOD_DIAGNOSTICS = True
PRINT_HEAT_G_DIAGNOSTICS = True

# The diagnostics below do not alter the LEC result.  They decompose the heat
# contribution into SW/LW/SH/LH, print the spatial rho*-J relationship, and
# compare alternative rho_s choices in J=H/(rho_s cp).

# If >= this fraction of the area-time samples have the same sign, print a
# heuristic sign clue.  This is diagnostic only; metadata/source code remains
# the authoritative way to determine the sign convention.
HEAT_SIGN_CLUE_FRACTION = 0.80

# ------------------------------------------------------------
# GSW / potential-density thermodynamic settings
# ------------------------------------------------------------
# In-situ pressure is needed only to convert SP,t -> SA,CT.  This preserves
# the convention in the user-supplied density-generation script.
# Verify the pressure unit expected by your custom gsw_alternative module.
GSW_PRESSURE_PER_M = 0.1

# Potential density reference pressure.  0 gives potential density referenced
# to the sea surface, equivalent to sigma0 + 1000 kg m-3.
POTENTIAL_DENSITY_REF_PRESSURE = 0.0

# Include daily z0 in the in-situ pressure used for SP,t -> SA,CT conversion.
# It does NOT alter the fixed potential-density reference pressure above.
USE_Z0_IN_PRESSURE = True

# Finite-difference increments for d(rho_pot)/dT and d(rho_pot)/dS.
DT_GSW = 1.0e-3
DS_GSW = 1.0e-3

# ------------------------------------------------------------
# Mean fields used in covariance decomposition
# ------------------------------------------------------------
# True:
#   use xx-2014-2018.nc for mean tx,ty,u,v.
# False:
#   use wind/velocity means reconstructed from the same daily samples.
#
# Density covariance terms ALWAYS use potential density reconstructed from the
# same daily tt/ss samples.
USE_PROVIDED_MEANS = True

# ------------------------------------------------------------
# Physical/numerical constants
# ------------------------------------------------------------
G_GRAV = 9.8                 # m s-2
CP_SEAWATER = 4000.0         # J kg-1 K-1, following von Storch et al.
RHO_FRESHWATER = 1000.0      # kg m-3
R_EARTH = 6_371_000.0        # m

BAD_ABS = 1.0e30
FILL_VALUE = 1.0e35

# Physical validity range used for potential and precomputed actual density.
# Values outside this range are treated as invalid/land for the n0 diagnostic.
RHO_VALID_MIN = 900.0
RHO_VALID_MAX = 1100.0

# Print the full n0(z) profile before the daily calculation starts.
PRINT_N0_PROFILE = True

# Print progress every N successfully processed days.
PRINT_INTERVAL = 30

# NetCDF compression.
DEFLATE_LEVEL = 1


# ============================================================
# 2. FILE / ARRAY UTILITIES
# ============================================================

def daily_path(var: str, year: int, month: int, day: int) -> Path:
    return DAILY_ROOT / DAILY_TEMPLATE.format(
        var=var, year=year, month=month, day=day
    )


def mean_path(var: str) -> Path:
    return MEAN_ROOT / MEAN_TEMPLATE.format(var=var)


def iter_365_dates():
    for year in range(START_YEAR, END_YEAR + 1):
        for month, ndays in enumerate(DAYS_IN_MONTH, start=1):
            for day in range(1, ndays + 1):
                yield year, month, day


def clean_array(a) -> np.ndarray:
    """Convert masked/missing values to NaN and return float64."""
    if np.ma.isMaskedArray(a):
        a = a.filled(np.nan)
    a = np.asarray(a, dtype=np.float64)
    a[~np.isfinite(a)] = np.nan
    a[np.abs(a) >= BAD_ABS] = np.nan
    return a


def _surface_index(var) -> tuple:
    """
    Construct an index that selects time=0 and vertical level=0,
    while retaining horizontal dimensions.
    """
    idx = []
    for dim in var.dimensions:
        d = dim.lower()
        if d.startswith("time"):
            idx.append(0)
        elif d.startswith("lev") or d in {"depth", "z"}:
            idx.append(0)
        else:
            idx.append(slice(None))
    return tuple(idx)


def read_surface_file(
    path: Path,
    varname: str,
    j_slice: slice | None = None,
    i_slice: slice | None = None,
) -> np.ndarray:
    with nc.Dataset(path, "r") as ds:
        if varname not in ds.variables:
            raise KeyError(
                f"Variable '{varname}' not found in {path}. "
                f"Available: {list(ds.variables.keys())}"
            )
        var = ds.variables[varname]
        a = clean_array(var[_surface_index(var)])
        a = np.squeeze(a)
        if a.ndim != 2:
            raise ValueError(
                f"Expected a 2-D surface field from {path}:{varname}, "
                f"but got shape {a.shape}, dimensions={var.dimensions}"
            )

        if j_slice is not None or i_slice is not None:
            js = j_slice if j_slice is not None else slice(None)
            is_ = i_slice if i_slice is not None else slice(None)
            a = a[js, is_]

        return a


def read_lon_lat(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with nc.Dataset(path, "r") as ds:
        lon = clean_array(ds.variables["lon"][:]).squeeze()
        lat = clean_array(ds.variables["lat"][:]).squeeze()
    if lon.ndim != 1 or lat.ndim != 1:
        raise ValueError("This program expects 1-D lon and lat coordinates.")
    return lon, lat


def region_slices_from_coords(
    lon: np.ndarray,
    lat: np.ndarray,
) -> tuple[slice, slice]:
    """Return contiguous (j_slice, i_slice) for the requested lon/lat region."""
    lon_min, lon_max = sorted(LON_RANGE)
    lat_min, lat_max = sorted(LAT_RANGE)

    ii = np.where((lon >= lon_min) & (lon <= lon_max))[0]
    jj = np.where((lat >= lat_min) & (lat <= lat_max))[0]

    if ii.size == 0:
        raise ValueError(
            f"No longitude grid points found in requested range {LON_RANGE}. "
            f"Available lon range is {np.nanmin(lon)} .. {np.nanmax(lon)}."
        )
    if jj.size == 0:
        raise ValueError(
            f"No latitude grid points found in requested range {LAT_RANGE}. "
            f"Available lat range is {np.nanmin(lat)} .. {np.nanmax(lat)}."
        )

    if not np.all(np.diff(ii) == 1):
        raise ValueError("Selected longitude indices are not contiguous.")
    if not np.all(np.diff(jj) == 1):
        raise ValueError("Selected latitude indices are not contiguous.")

    return slice(jj[0], jj[-1] + 1), slice(ii[0], ii[-1] + 1)


def read_lev(path: Path, varname: str = "rho") -> np.ndarray:
    with nc.Dataset(path, "r") as ds:
        var = ds.variables[varname]
        lev_dim = None
        for dim in var.dimensions:
            d = dim.lower()
            if d.startswith("lev") or d in {"depth", "z"}:
                lev_dim = dim
                break
        if lev_dim is None:
            raise ValueError(f"No vertical dimension found in {path}:{varname}")
        if lev_dim not in ds.variables:
            raise ValueError(
                f"Vertical dimension '{lev_dim}' has no coordinate variable in {path}"
            )
        lev = clean_array(ds.variables[lev_dim][:]).squeeze()
    return lev


def read_level_from_open_dataset(
    ds,
    varname: str,
    k: int,
    j_slice: slice | None = None,
    i_slice: slice | None = None,
) -> np.ndarray:
    var = ds.variables[varname]
    idx = []
    for dim in var.dimensions:
        d = dim.lower()
        if d.startswith("time"):
            idx.append(0)
        elif d.startswith("lev") or d in {"depth", "z"}:
            idx.append(k)
        else:
            idx.append(slice(None))
    a = clean_array(var[tuple(idx)])
    a = np.squeeze(a)
    if a.ndim != 2:
        raise ValueError(
            f"Expected 2-D level from {varname}, got {a.shape}, "
            f"dimensions={var.dimensions}"
        )

    if j_slice is not None or i_slice is not None:
        js = j_slice if j_slice is not None else slice(None)
        is_ = i_slice if i_slice is not None else slice(None)
        a = a[js, is_]

    return a


def check_shape(name: str, a: np.ndarray, expected: tuple[int, int]):
    if a.shape != expected:
        raise ValueError(
            f"{name}: shape {a.shape} does not match expected {expected}"
        )


# ============================================================
# 3. GRID / AREA
# ============================================================

def centers_to_edges(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 1 or x.size < 2:
        raise ValueError("Coordinate axis must be 1-D with at least 2 points.")
    edges = np.empty(x.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return edges


def spherical_cell_area(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """
    Exact spherical quadrilateral area for a rectilinear lon-lat grid.
    Returned shape is (lat, lon).
    """
    lon_e = np.deg2rad(centers_to_edges(lon))
    lat_e_deg = centers_to_edges(lat)
    lat_e_deg = np.clip(lat_e_deg, -90.0, 90.0)
    lat_e = np.deg2rad(lat_e_deg)

    dlon = np.abs(np.diff(lon_e))
    dsin = np.abs(np.sin(lat_e[1:]) - np.sin(lat_e[:-1]))

    return (R_EARTH ** 2) * dsin[:, None] * dlon[None, :]


def area_weighted_mean_2d(field: np.ndarray, area: np.ndarray) -> float:
    valid = np.isfinite(field) & np.isfinite(area) & (area > 0.0)
    if not np.any(valid):
        return np.nan
    return float(np.sum(field[valid] * area[valid]) / np.sum(area[valid]))


# ============================================================
# 4. EOS80 SURFACE DENSITY AND DERIVATIVES
# ============================================================

def rho_unesco_surface(T: np.ndarray, S: np.ndarray) -> np.ndarray:
    """
    UNESCO 1983 / EOS-80 density at atmospheric pressure.

    T : degC
    S : psu
    rho : kg m-3

    This is used only to estimate d(rho)/dT and d(rho)/dS at the
    mean surface state.  It does NOT replace the supplied LICOM rho.
    """
    T = np.asarray(T, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)

    rho_w = (
        999.842594
        + 6.793952e-2 * T
        - 9.095290e-3 * T**2
        + 1.001685e-4 * T**3
        - 1.120083e-6 * T**4
        + 6.536332e-9 * T**5
    )

    A = (
        0.824493
        - 4.0899e-3 * T
        + 7.6438e-5 * T**2
        - 8.2467e-7 * T**3
        + 5.3875e-9 * T**4
    )
    B = (
        -5.72466e-3
        + 1.0227e-4 * T
        - 1.6546e-6 * T**2
    )
    C = 4.8314e-4

    # Avoid invalid S**1.5 over land/bad points.
    S_nonneg = np.where(np.isfinite(S) & (S >= 0.0), S, np.nan)

    return rho_w + A * S_nonneg + B * S_nonneg**1.5 + C * S_nonneg**2


def eos80_surface_derivatives(
    T: np.ndarray,
    S: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return:
        alpha0 = d(rho)/dT [kg m-3 K-1]
        beta0  = d(rho)/dS [kg m-3 psu-1]

    Note that alpha0 here follows the paper's derivative notation and is
    normally NEGATIVE in warm seawater. It is not the conventional
    positive thermal-expansion coefficient alpha = -(1/rho)d(rho)/dT.
    """
    drho_dT = (
        rho_unesco_surface(T + DT_EOS, S)
        - rho_unesco_surface(T - DT_EOS, S)
    ) / (2.0 * DT_EOS)

    # Use one-sided difference near S=0 if necessary.
    S_minus = S - DS_EOS
    centered = S_minus >= 0.0

    drho_dS = np.full_like(S, np.nan, dtype=np.float64)
    if np.any(centered):
        drho_dS[centered] = (
            rho_unesco_surface(T[centered], S[centered] + DS_EOS)
            - rho_unesco_surface(T[centered], S[centered] - DS_EOS)
        ) / (2.0 * DS_EOS)

    one_sided = np.isfinite(S) & ~centered
    if np.any(one_sided):
        drho_dS[one_sided] = (
            rho_unesco_surface(T[one_sided], S[one_sided] + DS_EOS)
            - rho_unesco_surface(T[one_sided], S[one_sided])
        ) / DS_EOS

    return drho_dT, drho_dS



# ============================================================
# GSW THERMODYNAMICS AND POTENTIAL DENSITY (rho terms; not n0)
# ============================================================

def _clean_gsw_result(a) -> np.ndarray:
    """Convert a GSW/custom-GSW result to finite float64/NaN."""
    return clean_array(a)


def pressure_at_level(
    lev_m: float,
    shape2d: tuple[int, int],
    z0: np.ndarray | None = None,
) -> np.ndarray:
    """
    Pressure field using the convention in the user's density script:
        p = GSW_PRESSURE_PER_M * (abs(lev) + z0)

    This is intentionally configurable because the custom gsw_alternative
    pressure unit must be checked by the user.
    """
    depth = abs(float(lev_m))
    if z0 is None:
        return np.full(shape2d, GSW_PRESSURE_PER_M * depth, dtype=np.float64)

    return GSW_PRESSURE_PER_M * (depth + z0)


def gsw_state_from_SP_t(
    SP: np.ndarray,
    t: np.ndarray,
    p: np.ndarray,
    lon2d: np.ndarray,
    lat2d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert Practical Salinity + in-situ temperature to SA, CT, and
    POTENTIAL DENSITY referenced to POTENTIAL_DENSITY_REF_PRESSURE.

        SA      = gsw_SA_from_SP(SP, p, lon, lat)
        CT      = gsw_CT_from_t(SA, t, p)
        rho_pot = gsw_rho(SA, CT, p_ref)

    For p_ref=0 this is the full potential density corresponding to sigma0+1000.
    The actual pressure p is used only in the SP/t -> SA/CT conversion.
    """
    SA = _clean_gsw_result(gsw_SA_from_SP(SP, p, lon2d, lat2d))
    CT = _clean_gsw_result(gsw_CT_from_t(SA, t, p))
    p_ref = np.full_like(p, POTENTIAL_DENSITY_REF_PRESSURE, dtype=np.float64)
    rho_pot = _clean_gsw_result(gsw_rho(SA, CT, p_ref))
    return SA, CT, rho_pot


def potential_density_gradient_two_levels(
    rho_pot0: np.ndarray,
    rho_pot1: np.ndarray,
    z0_m: float,
    z1_m: float,
) -> np.ndarray:
    """Local vertical gradient of potential density between two levels."""
    out = np.full(rho_pot0.shape, np.nan, dtype=np.float64)
    valid = valid_density_mask(rho_pot0) & valid_density_mask(rho_pot1)
    dz = float(z1_m - z0_m)
    if abs(dz) < 1.0e-12:
        raise RuntimeError("Zero vertical separation in potential-density gradient.")
    out[valid] = (rho_pot1[valid] - rho_pot0[valid]) / dz
    return out


def gsw_surface_derivatives_from_SP_t(
    SP: np.ndarray,
    t: np.ndarray,
    p: np.ndarray,
    lon2d: np.ndarray,
    lat2d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute derivatives of POTENTIAL DENSITY with respect to the stored
    surface variables at fixed actual pressure:

        alpha0 = d(rho_pot)/dT
        beta0  = d(rho_pot)/dSP

    Every perturbed state is first converted to SA/CT at the actual pressure,
    then evaluated at the fixed potential-density reference pressure.
    """
    p_ref = np.full_like(p, POTENTIAL_DENSITY_REF_PRESSURE, dtype=np.float64)

    # d rho_pot / d t
    SA0 = _clean_gsw_result(gsw_SA_from_SP(SP, p, lon2d, lat2d))
    CT_p = _clean_gsw_result(gsw_CT_from_t(SA0, t + DT_GSW, p))
    CT_m = _clean_gsw_result(gsw_CT_from_t(SA0, t - DT_GSW, p))
    rho_tp = _clean_gsw_result(gsw_rho(SA0, CT_p, p_ref))
    rho_tm = _clean_gsw_result(gsw_rho(SA0, CT_m, p_ref))
    alpha0 = (rho_tp - rho_tm) / (2.0 * DT_GSW)

    # d rho_pot / d SP
    SPp = SP + DS_GSW
    SPm = SP - DS_GSW
    SA_p = _clean_gsw_result(gsw_SA_from_SP(SPp, p, lon2d, lat2d))
    SA_m = _clean_gsw_result(gsw_SA_from_SP(SPm, p, lon2d, lat2d))
    CT_sp = _clean_gsw_result(gsw_CT_from_t(SA_p, t, p))
    CT_sm = _clean_gsw_result(gsw_CT_from_t(SA_m, t, p))
    rho_sp = _clean_gsw_result(gsw_rho(SA_p, CT_sp, p_ref))
    rho_sm = _clean_gsw_result(gsw_rho(SA_m, CT_sm, p_ref))
    beta0 = (rho_sp - rho_sm) / (2.0 * DS_GSW)

    bad = ~np.isfinite(SP) | ~np.isfinite(t) | ~np.isfinite(p)
    alpha0[bad] = np.nan
    beta0[bad] = np.nan
    return alpha0, beta0



def build_reference_potential_density_profile_from_mean_ts(
    ss_mean_file: Path,
    tt_mean_file: Path,
    area: np.ndarray,
    j_slice: slice,
    i_slice: slice,
    lon2d: np.ndarray,
    lat2d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build only the potential-density reference profile used in rho*.

    This function deliberately does NOT calculate n0.  n0 is obtained
    separately from the original precomputed rho files.
    """
    lev = read_lev(ss_mean_file, VAR_NAMES["ss"])
    nz = lev.size
    rho_ref_pot = np.full(nz, np.nan, dtype=np.float64)

    with nc.Dataset(ss_mean_file, "r") as dss, nc.Dataset(tt_mean_file, "r") as dst:
        for k in range(nz):
            SP = read_level_from_open_dataset(
                dss, VAR_NAMES["ss"], k, j_slice, i_slice
            )
            t = read_level_from_open_dataset(
                dst, VAR_NAMES["tt"], k, j_slice, i_slice
            )
            p = pressure_at_level(lev[k], area.shape)
            _, _, rho_pot = gsw_state_from_SP_t(
                SP, t, p, lon2d, lat2d
            )
            valid = valid_density_mask(rho_pot)
            rho_ref_pot[k] = area_weighted_mean_with_mask(
                rho_pot, area, valid
            )

    if np.any(~np.isfinite(rho_ref_pot)):
        bad = np.where(~np.isfinite(rho_ref_pot))[0]
        raise RuntimeError(
            f"Potential-density rho_ref invalid at levels: {bad.tolist()}"
        )

    return lev, rho_ref_pot


def actual_density_gradient_two_levels(
    rho0: np.ndarray,
    rho1: np.ndarray,
    z0_m: float,
    z1_m: float,
) -> np.ndarray:
    """Local d(rho_actual)/dz between the first two precomputed rho levels."""
    out = np.full(rho0.shape, np.nan, dtype=np.float64)
    valid = valid_density_mask(rho0) & valid_density_mask(rho1)
    dz = float(z1_m - z0_m)
    if abs(dz) < 1.0e-12:
        raise RuntimeError("Zero vertical separation in actual-density n0.")
    out[valid] = (rho1[valid] - rho0[valid]) / dz
    return out


def UNUSED_build_reference_potential_density_and_n0_from_mean_ts(
    ss_mean_file: Path,
    tt_mean_file: Path,
    area: np.ndarray,
    j_slice: slice,
    i_slice: slice,
    lon2d: np.ndarray,
    lat2d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build rho_ref(z) and n0(z) from 2014-2018 mean T/S using POTENTIAL DENSITY.

    rho_ref(z): horizontal area mean of potential density at each level.

    n0(z): first compute local d(rho_pot)/dz between adjacent model levels,
    horizontally area average those pair gradients, then interpolate the pair
    values from layer mid-depths to model-level depths.

    The surface n0 used by G is later replaced by the time mean of DAILY
    potential-density gradients between levels 0 and 1.
    """
    lev = read_lev(ss_mean_file, VAR_NAMES["ss"])
    zcoord = -np.abs(lev.astype(np.float64))
    dep = np.abs(lev.astype(np.float64))
    nz = lev.size
    if nz < 2:
        raise RuntimeError("At least two T/S levels are required.")

    rho_ref = np.full(nz, np.nan, dtype=np.float64)
    pair_n0 = np.full(nz - 1, np.nan, dtype=np.float64)
    pair_count = np.zeros(nz - 1, dtype=np.int64)

    with nc.Dataset(ss_mean_file, "r") as dss, nc.Dataset(tt_mean_file, "r") as dst:
        SP_prev = read_level_from_open_dataset(
            dss, VAR_NAMES["ss"], 0, j_slice, i_slice
        )
        t_prev = read_level_from_open_dataset(
            dst, VAR_NAMES["tt"], 0, j_slice, i_slice
        )
        p_prev = pressure_at_level(lev[0], area.shape)
        SA_prev, CT_prev, rho_prev = gsw_state_from_SP_t(
            SP_prev, t_prev, p_prev, lon2d, lat2d
        )
        valid_prev = valid_density_mask(rho_prev)
        rho_ref[0] = area_weighted_mean_with_mask(rho_prev, area, valid_prev)

        for k in range(1, nz):
            SP_curr = read_level_from_open_dataset(
                dss, VAR_NAMES["ss"], k, j_slice, i_slice
            )
            t_curr = read_level_from_open_dataset(
                dst, VAR_NAMES["tt"], k, j_slice, i_slice
            )
            p_curr = pressure_at_level(lev[k], area.shape)
            SA_curr, CT_curr, rho_curr = gsw_state_from_SP_t(
                SP_curr, t_curr, p_curr, lon2d, lat2d
            )

            valid_curr = valid_density_mask(rho_curr)
            rho_ref[k] = area_weighted_mean_with_mask(rho_curr, area, valid_curr)

            grad = potential_density_gradient_two_levels(
                rho_prev, rho_curr, zcoord[k-1], zcoord[k]
            )
            valid_grad = np.isfinite(grad)
            pair_n0[k-1] = area_weighted_mean_with_mask(grad, area, valid_grad)
            pair_count[k-1] = int(np.count_nonzero(valid_grad))

            SP_prev = SP_curr
            t_prev = t_curr
            p_prev = p_curr
            SA_prev = SA_curr
            CT_prev = CT_curr
            rho_prev = rho_curr

    if np.any(~np.isfinite(rho_ref)):
        bad = np.where(~np.isfinite(rho_ref))[0]
        raise RuntimeError(f"Potential-density rho_ref invalid at levels: {bad.tolist()}")
    if np.any(~np.isfinite(pair_n0)):
        bad = np.where(~np.isfinite(pair_n0))[0]
        raise RuntimeError(f"Potential-density n0 invalid at level pairs: {bad.tolist()}")

    dep_mid = 0.5 * (dep[:-1] + dep[1:])
    n0 = np.interp(dep, dep_mid, pair_n0, left=pair_n0[0], right=pair_n0[-1])

    n0_valid_count = np.empty(nz, dtype=np.int64)
    n0_valid_count[0] = pair_count[0]
    n0_valid_count[-1] = pair_count[-1]
    if nz > 2:
        n0_valid_count[1:-1] = np.minimum(pair_count[:-1], pair_count[1:])

    return lev, rho_ref, n0, n0_valid_count


def print_potential_density_surface_diagnostics(
    lon: np.ndarray,
    lat: np.ndarray,
    lev: np.ndarray,
    rho0_mean: np.ndarray,
    rho1_mean: np.ndarray,
    n0_surface_local_mean: np.ndarray,
    n0_surface: float,
    n0_profile_mean_ts: np.ndarray,
    rho_ref: np.ndarray,
    n0_valid_count: np.ndarray,
    area: np.ndarray,
):
    """Print diagnostics for the potential-density reference state and n0."""
    zcoord = -np.abs(lev.astype(np.float64))
    valid = (
        valid_density_mask(rho0_mean)
        & valid_density_mask(rho1_mean)
        & np.isfinite(n0_surface_local_mean)
    )
    if not np.any(valid):
        raise RuntimeError("No valid surface potential-density gradient points.")

    dz = zcoord[1] - zcoord[0]
    direct_grad = np.full(area.shape, np.nan, dtype=np.float64)
    direct_grad[valid] = (rho1_mean[valid] - rho0_mean[valid]) / dz

    direct_area_mean = area_weighted_mean_with_mask(direct_grad, area, valid)
    dailymean_area_mean = area_weighted_mean_with_mask(
        n0_surface_local_mean, area, valid
    )

    nneg = int(np.count_nonzero(valid & (n0_surface_local_mean < 0.0)))
    npos = int(np.count_nonzero(valid & (n0_surface_local_mean > 0.0)))
    nzero = int(np.count_nonzero(valid & (n0_surface_local_mean == 0.0)))
    ntot = nneg + npos + nzero

    print("\n" + "=" * 104)
    print("n0 CHECK: daily T/S -> GSW potential density -> d(rho_pot)/dz")
    print("=" * 104)
    print(
        f"potential-density reference pressure = "
        f"{POTENTIAL_DENSITY_REF_PRESSURE:g} (gsw_alternative pressure units)"
    )
    print(f"level 0 = {lev[0]:.8f} m, level 1 = {lev[1]:.8f} m, dz={dz:.8f} m")
    print(f"valid surface points = {ntot}")
    print(
        "area-mean gradient from mean rho_pot levels = "
        f"{direct_area_mean:.10e} kg m-4"
    )
    print(
        "area-mean of time-mean daily gradients      = "
        f"{dailymean_area_mean:.10e} kg m-4"
    )
    print(f"surface n0 used in G                       = {n0_surface:.10e} kg m-4")
    print(
        "potential-density n0 signs: "
        f"negative={nneg} ({100.0*nneg/ntot:.2f}%), "
        f"positive={npos} ({100.0*npos/ntot:.2f}%), "
        f"zero={nzero} ({100.0*nzero/ntot:.2f}%)"
    )

    tmp = np.where(valid, n0_surface_local_mean, np.nan)
    jmin, imin = np.unravel_index(np.nanargmin(tmp), tmp.shape)
    jmax, imax = np.unravel_index(np.nanargmax(tmp), tmp.shape)
    print(
        "most negative potential-density n0: "
        f"lon={lon[imin]:.6f}, lat={lat[jmin]:.6f}, n0={tmp[jmin,imin]:.10e}"
    )
    print(
        "most positive potential-density n0: "
        f"lon={lon[imax]:.6f}, lat={lat[jmax]:.6f}, n0={tmp[jmax,imax]:.10e}"
    )

    if PRINT_N0_PROFILE:
        print("\nFull diagnostic potential-density profile from 2014-2018 mean T/S:")
        print(
            f"{'k':>3s} {'lev(m)':>12s} {'rho_ref_pot':>15s} "
            f"{'n0=d(rho_pot)/dz':>20s} {'valid_pts':>12s}"
        )
        print("-" * 78)
        for k in range(lev.size):
            print(
                f"{k:3d} {lev[k]:12.5f} {rho_ref[k]:15.8f} "
                f"{n0_profile_mean_ts[k]:20.8e} {int(n0_valid_count[k]):12d}"
            )
    print("=" * 104 + "\n")


# ============================================================
# 5. SAFE ACCUMULATION
# ============================================================

def add_sum(sum_array, count_array, value, valid):
    sum_array[valid] += value[valid]
    count_array[valid] += 1


def safe_mean(sum_array, count_array):
    out = np.full(sum_array.shape, np.nan, dtype=np.float64)
    valid = count_array > 0
    out[valid] = sum_array[valid] / count_array[valid]
    return out


def finite_common(*arrays):
    valid = np.ones(arrays[0].shape, dtype=bool)
    for a in arrays:
        valid &= np.isfinite(a)
    return valid



# ============================================================
# HEAT-FLUX DIAGNOSTICS
# ============================================================

def read_variable_metadata(path: Path, varname: str) -> dict:
    """Read useful NetCDF metadata for one variable."""
    with nc.Dataset(path, "r") as ds:
        if varname not in ds.variables:
            raise KeyError(f"Variable '{varname}' not found in {path}")
        var = ds.variables[varname]

        meta = {
            "path": str(path),
            "variable": varname,
            "dimensions": tuple(var.dimensions),
            "dtype": str(var.dtype),
        }

        for att in (
            "units",
            "long_name",
            "standard_name",
            "positive",
            "comment",
            "description",
            "sign_convention",
            "coordinates",
        ):
            if att in var.ncattrs():
                meta[att] = getattr(var, att)

        return meta


def print_heat_file_metadata(paths: dict):
    """Print metadata that may reveal the stored heat-flux sign convention."""
    print("\n" + "=" * 96)
    print("HEAT-FLUX FILE METADATA CHECK")
    print("=" * 96)
    print("Required by this LEC program: net H must be POSITIVE INTO THE OCEAN.")
    print("Metadata such as 'positive', 'long_name', 'comment', or 'standard_name'")
    print("may identify whether each stored component is upward/downward positive.")
    print("-" * 96)

    for key in ("sw", "lw", "sh", "lh"):
        meta = read_variable_metadata(paths[key], VAR_NAMES[key])
        print(f"{key.upper():>3s}: {meta['path']}")
        print(
            f"     variable={meta['variable']}, "
            f"dtype={meta['dtype']}, dimensions={meta['dimensions']}"
        )

        found = False
        for att in (
            "units",
            "long_name",
            "standard_name",
            "positive",
            "sign_convention",
            "comment",
            "description",
        ):
            if att in meta:
                print(f"     {att} = {meta[att]}")
                found = True

        if not found:
            print("     No useful sign-convention metadata found.")

    print("-" * 96)
    print(
        "Configured expression used in the LEC calculation: H = "
        + " ".join(
            f"{HEAT_FLUX_SIGNS[k]:+g}*{k}"
            for k in ("sw", "lw", "sh", "lh")
        )
    )
    print("=" * 96 + "\n")


def heat_weighted_stats(
    field: np.ndarray,
    area: np.ndarray,
    ocean_mask: np.ndarray,
) -> dict:
    """
    Return spatial diagnostics over valid ocean points.

    Positive/negative fractions are AREA-weighted, not simple grid-point counts.
    """
    valid = (
        ocean_mask
        & np.isfinite(field)
        & np.isfinite(area)
        & (area > 0.0)
    )

    if not np.any(valid):
        return {
            "n": 0,
            "min": np.nan,
            "mean": np.nan,
            "max": np.nan,
            "pos_frac": np.nan,
            "neg_frac": np.nan,
            "zero_frac": np.nan,
        }

    aa = area[valid]
    ff = field[valid]
    atot = float(np.sum(aa))

    pos = ff > 0.0
    neg = ff < 0.0
    zero = ff == 0.0

    return {
        "n": int(np.count_nonzero(valid)),
        "min": float(np.min(ff)),
        "mean": float(np.sum(ff * aa) / atot),
        "max": float(np.max(ff)),
        "pos_frac": float(np.sum(aa[pos]) / atot),
        "neg_frac": float(np.sum(aa[neg]) / atot),
        "zero_frac": float(np.sum(aa[zero]) / atot),
    }


def heat_sign_clue(name: str, stats: dict) -> str:
    """
    Return a HEURISTIC sign clue.

    This is intentionally conservative: it does not decide the convention.
    """
    if stats["n"] == 0:
        return "no valid ocean values"

    p = stats["pos_frac"]
    n = stats["neg_frac"]
    threshold = HEAT_SIGN_CLUE_FRACTION

    if name == "sw":
        if p >= threshold:
            return (
                "SW is predominantly positive; this is compatible with "
                "downward-positive shortwave input."
            )
        if n >= threshold:
            return (
                "SW is predominantly negative; check whether downward input is "
                "stored negative."
            )
        return "SW has mixed signs; sign cannot be inferred from values alone."

    # LW/SH/LH are commonly ocean heat losses, but can locally reverse.
    if n >= threshold:
        return (
            f"{name.upper()} is predominantly negative; this is compatible with "
            "an already-signed convention where ocean heat loss is negative."
        )
    if p >= threshold:
        return (
            f"{name.upper()} is predominantly positive; this is compatible with "
            "an upward-positive loss magnitude and may require a minus sign in H."
        )
    return (
        f"{name.upper()} has mixed signs; sign cannot be inferred from values alone."
    )


def print_one_day_heat_diagnostics(
    year: int,
    month: int,
    day: int,
    sw: np.ndarray,
    lw: np.ndarray,
    sh: np.ndarray,
    lh: np.ndarray,
    H_configured: np.ndarray,
    area: np.ndarray,
    ocean_mask: np.ndarray,
):
    """Print detailed first-day numerical checks for heat-flux convention."""
    H_all_stored = sw + lw + sh + lh
    H_sw_minus_outgoing = sw - lw - sh - lh

    print("\n" + "=" * 108)
    print(
        f"HEAT-FLUX NUMERICAL CHECK FOR {year:04d}-{month:02d}-{day:02d} "
        "(regional ocean points)"
    )
    print("=" * 108)
    print(
        f"{'field':<24s} {'min':>14s} {'area mean':>14s} {'max':>14s} "
        f"{'positive%':>11s} {'negative%':>11s}"
    )
    print("-" * 108)

    arrays = {
        "sw (stored)": sw,
        "lw (stored)": lw,
        "sh (stored)": sh,
        "lh (stored)": lh,
        "H configured": H_configured,
        "H all stored signs": H_all_stored,
        "H sw-lw-sh-lh": H_sw_minus_outgoing,
    }

    stats_map = {}
    for name, arr in arrays.items():
        stats = heat_weighted_stats(arr, area, ocean_mask)
        stats_map[name] = stats
        print(
            f"{name:<24s} "
            f"{stats['min']:14.6e} {stats['mean']:14.6e} {stats['max']:14.6e} "
            f"{100.0*stats['pos_frac']:10.2f}% "
            f"{100.0*stats['neg_frac']:10.2f}%"
        )

    print("-" * 108)
    print("Heuristic component sign clues (NOT a substitute for model documentation):")
    for key in ("sw", "lw", "sh", "lh"):
        stats = stats_map[f"{key} (stored)"]
        print(f"  {key.upper():>3s}: {heat_sign_clue(key, stats)}")

    print("\nTwo common candidate conventions:")
    print("  A) H_all_stored       = sw + lw + sh + lh")
    print("  B) H_sw_minus_outgoing = sw - lw - sh - lh")
    print(
        "  Actual program setting: H = "
        + " ".join(
            f"{HEAT_FLUX_SIGNS[k]:+g}*{k}"
            for k in ("sw", "lw", "sh", "lh")
        )
    )

    # Warn when the configured sign visibly conflicts with a strong numerical clue.
    warnings = []
    for key in ("lw", "sh", "lh"):
        stats = stats_map[f"{key} (stored)"]
        if (
            stats["pos_frac"] >= HEAT_SIGN_CLUE_FRACTION
            and HEAT_FLUX_SIGNS[key] > 0.0
        ):
            warnings.append(
                f"{key.upper()} is mostly positive but is ADDED by the current "
                "HEAT_FLUX_SIGNS. If it is an upward-positive ocean heat loss, "
                "its configured sign should be -1."
            )
        if (
            stats["neg_frac"] >= HEAT_SIGN_CLUE_FRACTION
            and HEAT_FLUX_SIGNS[key] < 0.0
        ):
            warnings.append(
                f"{key.upper()} is mostly negative but is SUBTRACTED by the current "
                "HEAT_FLUX_SIGNS. If negative already denotes ocean heat loss, "
                "its configured sign should be +1."
            )

    if warnings:
        print("\n*** HEAT SIGN WARNINGS (heuristic) ***")
        for msg in warnings:
            print("  WARNING:", msg)
    else:
        print("\nNo strong numerical conflict with the configured signs was detected.")

    print("=" * 108 + "\n")


def init_heat_period_diagnostics() -> dict:
    """Initialize area-time accumulators for whole-period heat diagnostics."""
    names = (
        "sw",
        "lw",
        "sh",
        "lh",
        "H_configured",
        "H_all_stored",
        "H_sw_minus_outgoing",
    )
    return {
        name: {
            "weighted_sum": 0.0,
            "area_sum": 0.0,
            "positive_area": 0.0,
            "negative_area": 0.0,
            "zero_area": 0.0,
            "min": np.inf,
            "max": -np.inf,
            "samples": 0,
        }
        for name in names
    }


def update_heat_period_diagnostic(
    diag: dict,
    name: str,
    field: np.ndarray,
    area: np.ndarray,
    ocean_mask: np.ndarray,
):
    """Accumulate one day's area-weighted heat statistics."""
    valid = (
        ocean_mask
        & np.isfinite(field)
        & np.isfinite(area)
        & (area > 0.0)
    )
    if not np.any(valid):
        return

    ff = field[valid]
    aa = area[valid]
    rec = diag[name]

    rec["weighted_sum"] += float(np.sum(ff * aa))
    rec["area_sum"] += float(np.sum(aa))
    rec["positive_area"] += float(np.sum(aa[ff > 0.0]))
    rec["negative_area"] += float(np.sum(aa[ff < 0.0]))
    rec["zero_area"] += float(np.sum(aa[ff == 0.0]))
    rec["min"] = min(rec["min"], float(np.min(ff)))
    rec["max"] = max(rec["max"], float(np.max(ff)))
    rec["samples"] += int(ff.size)


def heat_period_stats(rec: dict) -> dict:
    """Convert whole-period accumulator to printable statistics, including sample count."""
    if rec["area_sum"] <= 0.0:
        return {
            "n": 0,
            "mean": np.nan,
            "min": np.nan,
            "max": np.nan,
            "pos_frac": np.nan,
            "neg_frac": np.nan,
            "zero_frac": np.nan,
        }

    return {
        "n": rec["samples"],
        "mean": rec["weighted_sum"] / rec["area_sum"],
        "min": rec["min"],
        "max": rec["max"],
        "pos_frac": rec["positive_area"] / rec["area_sum"],
        "neg_frac": rec["negative_area"] / rec["area_sum"],
        "zero_frac": rec["zero_area"] / rec["area_sum"],
    }


def print_heat_period_diagnostics(diag: dict, valid_days: int):
    """Print the whole-period heat-component and convention comparison."""
    print("\n" + "=" * 108)
    print(
        f"HEAT-FLUX WHOLE-PERIOD CHECK: {valid_days} valid days "
        "(area-time weighted)"
    )
    print("=" * 108)
    print(
        f"{'field':<24s} {'min':>14s} {'mean':>14s} {'max':>14s} "
        f"{'positive%':>11s} {'negative%':>11s}"
    )
    print("-" * 108)

    labels = (
        ("sw", "sw (stored)"),
        ("lw", "lw (stored)"),
        ("sh", "sh (stored)"),
        ("lh", "lh (stored)"),
        ("H_configured", "H configured"),
        ("H_all_stored", "H all stored signs"),
        ("H_sw_minus_outgoing", "H sw-lw-sh-lh"),
    )

    stats_map = {}
    for key, label in labels:
        stats = heat_period_stats(diag[key])
        stats_map[key] = stats
        print(
            f"{label:<24s} "
            f"{stats['min']:14.6e} {stats['mean']:14.6e} {stats['max']:14.6e} "
            f"{100.0*stats['pos_frac']:10.2f}% "
            f"{100.0*stats['neg_frac']:10.2f}%"
        )

    print("-" * 108)
    print("Whole-period heuristic sign clues:")
    for key in ("sw", "lw", "sh", "lh"):
        print(f"  {key.upper():>3s}: {heat_sign_clue(key, stats_map[key])}")

    print("\nInterpretation guide:")
    print(
        "  * If LW/SH/LH are predominantly NEGATIVE, they are likely already "
        "signed as ocean heat losses; sw+lw+sh+lh is then a plausible net H."
    )
    print(
        "  * If LW/SH/LH are predominantly POSITIVE and represent upward losses, "
        "sw-lw-sh-lh is the plausible net H."
    )
    print(
        "  * This numerical check is not definitive. The LICOM source code or "
        "forcing/output variable definition is authoritative."
    )
    print(
        "  * The LEC code requires H > 0 INTO THE OCEAN. Current configured H: "
        + " ".join(
            f"{HEAT_FLUX_SIGNS[k]:+g}*{k}"
            for k in ("sw", "lw", "sh", "lh")
        )
    )
    print("=" * 108 + "\n")



def area_integral_masked(field, area, mask=None) -> float:
    valid = np.isfinite(field) & np.isfinite(area) & (area > 0.0)
    if mask is not None:
        valid &= mask
    if not np.any(valid):
        return np.nan
    return float(np.sum(field[valid] * area[valid]))


def area_mean_masked(field, area, mask=None) -> float:
    valid = np.isfinite(field) & np.isfinite(area) & (area > 0.0)
    if mask is not None:
        valid &= mask
    if not np.any(valid):
        return np.nan
    aa = area[valid]
    return float(np.sum(field[valid] * aa) / np.sum(aa))


def area_covariance(a, b, area, mask=None) -> float:
    valid = (
        np.isfinite(a) & np.isfinite(b)
        & np.isfinite(area) & (area > 0.0)
    )
    if mask is not None:
        valid &= mask
    if not np.any(valid):
        return np.nan

    aa = area[valid]
    av = a[valid]
    bv = b[valid]
    asum = np.sum(aa)
    am = np.sum(av * aa) / asum
    bm = np.sum(bv * aa) / asum
    return float(np.sum((av-am)*(bv-bm)*aa) / asum)


def area_correlation(a, b, area, mask=None) -> float:
    valid = (
        np.isfinite(a) & np.isfinite(b)
        & np.isfinite(area) & (area > 0.0)
    )
    if mask is not None:
        valid &= mask
    if not np.any(valid):
        return np.nan

    aa = area[valid]
    av = a[valid]
    bv = b[valid]
    asum = np.sum(aa)
    am = np.sum(av * aa) / asum
    bm = np.sum(bv * aa) / asum
    da = av-am
    db = bv-bm
    cov = np.sum(da*db*aa) / asum
    va = np.sum(da*da*aa) / asum
    vb = np.sum(db*db*aa) / asum
    if va <= 0.0 or vb <= 0.0:
        return np.nan
    return float(cov / np.sqrt(va*vb))


def print_heat_G_diagnostics(
    area,
    rho_star_mean,
    rho_bar_file,
    alpha0,
    n0_surface,
    J_bar,
    H_bar,
    component_J_bar,
    component_rhoJ_bar,
    rho_bar_heat,
):
    """
    Diagnose WHY G(Pm)_heat has its sign.

    This prints:
      1) rho*-J spatial covariance/correlation,
      2) SW/LW/SH/LH contributions to G(Pm)_heat and G(Pe)_heat,
      3) four physical quadrants (light/heavy x heating/cooling),
      4) sensitivity of G(Pm)_heat to rho_s in J = H/(rho_s cp).

    Nothing here changes the calculation.
    """
    factor_t = -G_GRAV * alpha0 / n0_surface
    common = finite_common(
        area, rho_star_mean, rho_bar_file, alpha0, J_bar, H_bar
    )

    print("\n" + "="*112)
    print("DETAILED G(Pm)_heat / G(Pe)_heat DIAGNOSTICS")
    print("="*112)
    print("Paper sign convention used here:")
    print("  alpha0 = d(rho)/dT < 0")
    print("  n0 < 0 for stable stratification")
    print("  H > 0 and J > 0 mean heating of the ocean")
    print("  factor_t = -g*alpha0/n0 is therefore normally NEGATIVE")
    print()

    print("Regional mean-field checks:")
    print(
        f"  area-mean rho*              = "
        f"{area_mean_masked(rho_star_mean, area, common): .8e} kg m-3"
    )
    print(
        f"  area-mean H                 = "
        f"{area_mean_masked(H_bar, area, common): .8e} W m-2"
    )
    print(
        f"  area-mean J                 = "
        f"{area_mean_masked(J_bar, area, common): .8e} K m s-1"
    )
    print(
        f"  area covariance(J, rho*)    = "
        f"{area_covariance(J_bar, rho_star_mean, area, common): .8e}"
    )
    print(
        f"  area correlation(J, rho*)   = "
        f"{area_correlation(J_bar, rho_star_mean, area, common): .6f}"
    )
    print(
        f"  integral[J*rho* dA]         = "
        f"{area_integral_masked(J_bar*rho_star_mean, area, common): .8e}"
    )
    print(
        f"  integral[alpha*J*rho* dA]   = "
        f"{area_integral_masked(alpha0*J_bar*rho_star_mean, area, common): .8e}"
    )

    # Component decomposition. Linearity means these should sum to total heat G.
    print("\nHeat-component decomposition:")
    print(
        f"{'component':<10s} {'G(Pm)_heat TW':>18s} "
        f"{'G(Pe)_heat TW':>18s} {'mean J':>16s} {'corr(J,rho*)':>16s}"
    )
    print("-"*84)

    gpm_comp_sum = 0.0
    gpe_comp_sum = 0.0
    for key in ("sw", "lw", "sh", "lh"):
        jb = component_J_bar[key]
        rjb = component_rhoJ_bar[key]
        cov_rhoj = rjb - rho_bar_heat * jb
        gpm_k = factor_t * jb * rho_star_mean
        gpe_k = factor_t * cov_rhoj

        gpm_w = area_integral_masked(gpm_k, area)
        gpe_w = area_integral_masked(gpe_k, area)
        gpm_comp_sum += gpm_w
        gpe_comp_sum += gpe_w

        print(
            f"{key.upper():<10s} {gpm_w/1e12:18.8f} "
            f"{gpe_w/1e12:18.8f} "
            f"{area_mean_masked(jb,area):16.8e} "
            f"{area_correlation(jb,rho_star_mean,area):16.6f}"
        )

    gpm_total = factor_t * J_bar * rho_star_mean
    cov_total = component_rhoJ_bar["total"] - rho_bar_heat * J_bar
    gpe_total = factor_t * cov_total
    gpm_total_w = area_integral_masked(gpm_total, area)
    gpe_total_w = area_integral_masked(gpe_total, area)

    print("-"*84)
    print(
        f"{'SUM4':<10s} {gpm_comp_sum/1e12:18.8f} "
        f"{gpe_comp_sum/1e12:18.8f}"
    )
    print(
        f"{'TOTAL':<10s} {gpm_total_w/1e12:18.8f} "
        f"{gpe_total_w/1e12:18.8f}"
    )
    print(
        f"  closure G(Pm)_heat SUM4-TOTAL = "
        f"{(gpm_comp_sum-gpm_total_w)/1e12:.6e} TW"
    )
    print(
        f"  closure G(Pe)_heat SUM4-TOTAL = "
        f"{(gpe_comp_sum-gpe_total_w)/1e12:.6e} TW"
    )

    # Quadrant decomposition.
    # With factor_t < 0:
    #   light+heating, dense+cooling => positive generation
    #   dense+heating, light+cooling => negative generation
    print("\nG(Pm)_heat physical-quadrant decomposition:")
    quadrants = (
        ("light + heating",  (rho_star_mean < 0.0) & (J_bar > 0.0)),
        ("dense + cooling",  (rho_star_mean > 0.0) & (J_bar < 0.0)),
        ("dense + heating",  (rho_star_mean > 0.0) & (J_bar > 0.0)),
        ("light + cooling",  (rho_star_mean < 0.0) & (J_bar < 0.0)),
    )
    for label, mask in quadrants:
        w = area_integral_masked(gpm_total, area, mask)
        a = area_integral_masked(np.ones_like(area), area, mask)
        atot = area_integral_masked(np.ones_like(area), area, common)
        frac = 100.0*a/atot if np.isfinite(a) and atot > 0 else np.nan
        print(f"  {label:<18s}: {w/1e12: .8f} TW   area={frac:6.2f}%")

    # Test whether use of instantaneous rho_s in J can explain the sign.
    # This comparison is ONLY for G(Pm), because an offline mean field is
    # insufficient to reconstruct G(Pe) for the alternative definitions.
    J_const1025 = H_bar / (1025.0 * CP_SEAWATER)
    J_meanrho = H_bar / (rho_bar_file * CP_SEAWATER)
    gpm_const1025 = factor_t * J_const1025 * rho_star_mean
    gpm_meanrho = factor_t * J_meanrho * rho_star_mean

    print("\nSensitivity of G(Pm)_heat to rho_s used in J=H/(rho_s cp):")
    print(
        f"  current mean[J_daily=H/(rho_daily cp)] : "
        f"{gpm_total_w/1e12: .8f} TW"
    )
    print(
        f"  J=mean(H)/(mean(rho) cp)                : "
        f"{area_integral_masked(gpm_meanrho,area)/1e12: .8f} TW"
    )
    print(
        f"  J=mean(H)/(1025 cp)                     : "
        f"{area_integral_masked(gpm_const1025,area)/1e12: .8f} TW"
    )

    print("\nInterpretation:")
    print("  * If the four-component closure is near zero, SW/LW/SH/LH summation is internally consistent.")
    print("  * The component table shows exactly which heat-flux term drives the sign.")
    print("  * The quadrant table shows whether negative G(Pm)_heat comes from dense-water heating")
    print("    or light-water cooling.")
    print("  * If all three rho_s choices have the same sign and similar magnitude, rho_s is not the cause.")
    print("="*112 + "\n")


# ============================================================
# REFERENCE n0 PROFILE FROM PRECOMPUTED ACTUAL DENSITY
# ============================================================

def valid_density_mask(rho: np.ndarray) -> np.ndarray:
    """Physical/finite mask used for density-gradient calculations."""
    return (
        np.isfinite(rho)
        & (rho >= RHO_VALID_MIN)
        & (rho <= RHO_VALID_MAX)
    )


def area_weighted_mean_with_mask(
    field: np.ndarray,
    area: np.ndarray,
    valid: np.ndarray,
) -> float:
    valid = valid & np.isfinite(field) & np.isfinite(area) & (area > 0.0)
    if not np.any(valid):
        return np.nan
    return float(np.sum(field[valid] * area[valid]) / np.sum(area[valid]))


def three_point_derivative_nonuniform(
    f0: np.ndarray,
    f1: np.ndarray,
    f2: np.ndarray,
    x0: float,
    x1: float,
    x2: float,
) -> np.ndarray:
    """Second-order derivative at x1 on a nonuniform three-point stencil."""
    w0 = (x1 - x2) / ((x0 - x1) * (x0 - x2))
    w1 = (2.0 * x1 - x0 - x2) / ((x1 - x0) * (x1 - x2))
    w2 = (x1 - x0) / ((x2 - x0) * (x2 - x1))
    return w0 * f0 + w1 * f1 + w2 * f2


def build_reference_density_and_n0_from_actual_density(
    rho_mean_file: Path,
    rho_var: str,
    area: np.ndarray,
    j_slice: slice,
    i_slice: slice,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """
    Construct rho_ref(z) and n0(z), but compute n0 from ACTUAL density.

    rho_ref(z):
        Horizontal area-weighted mean density at each model level.  This is
        still needed for rho* = rho - rho_ref in the APE generation terms.

    n0(z):
        First calculate local d(rho)/dz from the actual 2014-2018 mean
        density field using common wet points, then horizontally area-average
        that local gradient.  Thus the derivative is taken BEFORE horizontal
        averaging.

    Surface level:
        n0_local = [rho(level 1)-rho(level 0)] / [z1-z0]
        using only points wet/valid in BOTH levels.

    Interior levels:
        nonuniform 3-point derivative using levels k-1,k,k+1 and only points
        valid in all three levels.

    Bottom level:
        first-order backward difference using the last two levels.

    z is positive upward, so stable stratification should normally produce
    negative n0.
    """
    lev = read_lev(rho_mean_file, rho_var)
    zcoord = -np.abs(lev.astype(np.float64))
    nz = lev.size

    if nz < 2:
        raise RuntimeError("At least two density levels are required for n0.")

    rho_ref = np.full(nz, np.nan, dtype=np.float64)
    n0 = np.full(nz, np.nan, dtype=np.float64)
    n0_valid_count = np.zeros(nz, dtype=np.int64)

    # --------------------------------------------------------
    # First pass: rho_ref(z) from the actual mean density field.
    # --------------------------------------------------------
    with nc.Dataset(rho_mean_file, "r") as ds:
        for k in range(nz):
            rho_k = read_level_from_open_dataset(ds, rho_var, k, j_slice, i_slice)
            check_shape(f"rho mean level {k}", rho_k, area.shape)
            valid_k = valid_density_mask(rho_k)
            rho_ref[k] = area_weighted_mean_with_mask(
                rho_k, area, valid_k
            )

    if np.any(~np.isfinite(rho_ref)):
        bad = np.where(~np.isfinite(rho_ref))[0]
        raise RuntimeError(
            f"rho_ref contains invalid levels: {bad.tolist()}"
        )

    # --------------------------------------------------------
    # Second pass: local d(rho)/dz first, horizontal mean second.
    # --------------------------------------------------------
    with nc.Dataset(rho_mean_file, "r") as ds:
        rho0 = read_level_from_open_dataset(ds, rho_var, 0, j_slice, i_slice)
        rho1 = read_level_from_open_dataset(ds, rho_var, 1, j_slice, i_slice)
        check_shape("rho mean level 0", rho0, area.shape)
        check_shape("rho mean level 1", rho1, area.shape)

        # Surface: first-order derivative from the ACTUAL first two levels.
        valid01 = valid_density_mask(rho0) & valid_density_mask(rho1)
        n0_surface_local = np.full(area.shape, np.nan, dtype=np.float64)
        dz01 = zcoord[1] - zcoord[0]
        n0_surface_local[valid01] = (
            rho1[valid01] - rho0[valid01]
        ) / dz01
        n0[0] = area_weighted_mean_with_mask(
            n0_surface_local, area, valid01
        )
        n0_valid_count[0] = int(np.count_nonzero(valid01))

        # Interior levels: local nonuniform 3-point derivative.
        rho_prev = rho0
        rho_curr = rho1

        for k in range(1, nz - 1):
            rho_next = read_level_from_open_dataset(ds, rho_var, k + 1, j_slice, i_slice)
            check_shape(f"rho mean level {k+1}", rho_next, area.shape)

            valid = (
                valid_density_mask(rho_prev)
                & valid_density_mask(rho_curr)
                & valid_density_mask(rho_next)
            )

            grad_local = three_point_derivative_nonuniform(
                rho_prev,
                rho_curr,
                rho_next,
                zcoord[k - 1],
                zcoord[k],
                zcoord[k + 1],
            )
            grad_local[~valid] = np.nan

            n0[k] = area_weighted_mean_with_mask(
                grad_local, area, valid
            )
            n0_valid_count[k] = int(np.count_nonzero(valid))

            rho_prev = rho_curr
            rho_curr = rho_next

        # Bottom: first-order backward difference on common wet points.
        valid_last = valid_density_mask(rho_prev) & valid_density_mask(rho_curr)
        grad_last = np.full(area.shape, np.nan, dtype=np.float64)
        dz_last = zcoord[-1] - zcoord[-2]
        grad_last[valid_last] = (
            rho_curr[valid_last] - rho_prev[valid_last]
        ) / dz_last
        n0[-1] = area_weighted_mean_with_mask(
            grad_last, area, valid_last
        )
        n0_valid_count[-1] = int(np.count_nonzero(valid_last))

    return (
        lev,
        rho_ref,
        n0,
        n0_valid_count,
        rho0,
        rho1,
        n0_surface_local,
    )


def print_n0_diagnostics(
    lon: np.ndarray,
    lat: np.ndarray,
    lev: np.ndarray,
    rho_ref: np.ndarray,
    n0: np.ndarray,
    n0_valid_count: np.ndarray,
    rho0: np.ndarray,
    rho1: np.ndarray,
    n0_surface_local: np.ndarray,
    area: np.ndarray,
):
    """Print detailed diagnostics for checking the actual-density n0."""
    zcoord = -np.abs(lev.astype(np.float64))
    valid01 = (
        valid_density_mask(rho0)
        & valid_density_mask(rho1)
        & np.isfinite(n0_surface_local)
    )

    if not np.any(valid01):
        raise RuntimeError("No common valid wet points in the first two rho levels.")

    drho01 = np.full(area.shape, np.nan, dtype=np.float64)
    drho01[valid01] = rho1[valid01] - rho0[valid01]

    nneg = int(np.count_nonzero(valid01 & (n0_surface_local < 0.0)))
    npos = int(np.count_nonzero(valid01 & (n0_surface_local > 0.0)))
    nzero = int(np.count_nonzero(valid01 & (n0_surface_local == 0.0)))
    ntot = nneg + npos + nzero

    # Comparison with the old "differentiate independently averaged rho_ref" method.
    old_style_n0_01 = (
        (rho_ref[1] - rho_ref[0]) / (zcoord[1] - zcoord[0])
    )

    print("\n" + "=" * 92)
    print("n0 CHECK: actual-density gradient first, horizontal average second")
    print("=" * 92)
    print(f"level 0: lev={lev[0]:.8f} m, z={zcoord[0]:.8f} m")
    print(f"level 1: lev={lev[1]:.8f} m, z={zcoord[1]:.8f} m")
    print(f"dz(1-0) = {zcoord[1]-zcoord[0]:.8f} m")
    print(f"common wet/valid points = {ntot}")
    print(
        "rho(level 0) [kg m-3]: "
        f"min={np.nanmin(rho0[valid01]):.8f}, "
        f"mean={np.nanmean(rho0[valid01]):.8f}, "
        f"max={np.nanmax(rho0[valid01]):.8f}"
    )
    print(
        "rho(level 1) [kg m-3]: "
        f"min={np.nanmin(rho1[valid01]):.8f}, "
        f"mean={np.nanmean(rho1[valid01]):.8f}, "
        f"max={np.nanmax(rho1[valid01]):.8f}"
    )
    print(
        "delta rho (level1-level0) [kg m-3]: "
        f"min={np.nanmin(drho01[valid01]):.8e}, "
        f"mean={np.nanmean(drho01[valid01]):.8e}, "
        f"max={np.nanmax(drho01[valid01]):.8e}"
    )
    print(
        "local n0(surface) [kg m-4]: "
        f"min={np.nanmin(n0_surface_local[valid01]):.8e}, "
        f"mean={np.nanmean(n0_surface_local[valid01]):.8e}, "
        f"max={np.nanmax(n0_surface_local[valid01]):.8e}"
    )
    print(
        "area-weighted local n0(surface) = "
        f"{area_weighted_mean_with_mask(n0_surface_local, area, valid01):.10e} kg m-4"
    )
    print(
        "old-style [rho_ref(1)-rho_ref(0)]/dz = "
        f"{old_style_n0_01:.10e} kg m-4"
    )
    print(
        "surface local n0 signs: "
        f"negative={nneg} ({100.0*nneg/ntot:.2f}%), "
        f"positive={npos} ({100.0*npos/ntot:.2f}%), "
        f"zero={nzero} ({100.0*nzero/ntot:.2f}%)"
    )

    # Print the locations of the most negative and most positive local n0 values.
    tmp_min = np.where(valid01, n0_surface_local, np.nan)
    tmp_max = tmp_min
    jmin, imin = np.unravel_index(np.nanargmin(tmp_min), tmp_min.shape)
    jmax, imax = np.unravel_index(np.nanargmax(tmp_max), tmp_max.shape)
    print(
        "most negative local n0: "
        f"lon={lon[imin]:.6f}, lat={lat[jmin]:.6f}, "
        f"rho0={rho0[jmin,imin]:.8f}, rho1={rho1[jmin,imin]:.8f}, "
        f"n0={n0_surface_local[jmin,imin]:.10e}"
    )
    print(
        "most positive local n0: "
        f"lon={lon[imax]:.6f}, lat={lat[jmax]:.6f}, "
        f"rho0={rho0[jmax,imax]:.8f}, rho1={rho1[jmax,imax]:.8f}, "
        f"n0={n0_surface_local[jmax,imax]:.10e}"
    )

    if PRINT_N0_PROFILE:
        print("\nFull horizontally averaged n0 profile from local actual-density gradients:")
        print(
            f"{'k':>3s} {'lev(m)':>12s} {'z(m)':>12s} "
            f"{'rho_ref':>15s} {'n0(kg m-4)':>16s} {'valid_pts':>12s}"
        )
        print("-" * 78)
        for k in range(lev.size):
            print(
                f"{k:3d} {lev[k]:12.5f} {zcoord[k]:12.5f} "
                f"{rho_ref[k]:15.8f} {n0[k]:16.8e} "
                f"{int(n0_valid_count[k]):12d}"
            )
    print("=" * 92 + "\n")


# ============================================================
# 7. DIAGNOSTICS
# ============================================================

def rms_difference(a: np.ndarray, b: np.ndarray) -> float:
    valid = np.isfinite(a) & np.isfinite(b)
    if not np.any(valid):
        return np.nan
    d = a[valid] - b[valid]
    return float(np.sqrt(np.mean(d * d)))


def max_abs_difference(a: np.ndarray, b: np.ndarray) -> float:
    valid = np.isfinite(a) & np.isfinite(b)
    if not np.any(valid):
        return np.nan
    return float(np.max(np.abs(a[valid] - b[valid])))


def integrate_w(field: np.ndarray, area: np.ndarray) -> float:
    valid = np.isfinite(field) & np.isfinite(area)
    if not np.any(valid):
        return np.nan
    return float(np.sum(field[valid] * area[valid]))


# ============================================================
# 8. NETCDF OUTPUT
# ============================================================

def write_output(
    lon,
    lat,
    lev,
    area,
    rho_ref,
    n0,
    n0_surface_local,
    rho_surface_level0,
    rho_actual_level0,
    rho_actual_level1,
    alpha0,
    beta0,
    rho_star_mean,
    j_mean,
    gs_mean,
    fields,
    totals,
    valid_day_count,
):
    OUTPUT_NC.parent.mkdir(parents=True, exist_ok=True)

    with nc.Dataset(OUTPUT_NC, "w", format="NETCDF4") as ds:
        ds.createDimension("lon", lon.size)
        ds.createDimension("lat", lat.size)
        ds.createDimension("lev", lev.size)

        vlon = ds.createVariable("lon", "f8", ("lon",))
        vlat = ds.createVariable("lat", "f8", ("lat",))
        vlev = ds.createVariable("lev", "f8", ("lev",))

        vlon[:] = lon
        vlat[:] = lat
        vlev[:] = lev

        vlon.units = "degrees_east"
        vlat.units = "degrees_north"
        vlev.units = "m"
        vlev.long_name = "input vertical coordinate"

        varea = ds.createVariable(
            "cell_area", "f8", ("lat", "lon"),
            zlib=True, complevel=DEFLATE_LEVEL
        )
        varea[:] = area
        varea.units = "m2"

        vrref = ds.createVariable("rho_ref", "f8", ("lev",))
        vn0 = ds.createVariable("n0", "f8", ("lev",))
        vrref[:] = rho_ref
        vn0[:] = n0
        vrref.units = "kg m-3"
        vrref.long_name = "horizontal area-weighted mean potential density"
        vn0.units = "kg m-4"
        vn0.long_name = (
            "actual-density gradient from precomputed rho files; "
            "surface value is the 2014-2018 mean of daily level-0/1 gradients"
        )

        auxiliary = {
            "n0_surface_local": (
                n0_surface_local, "kg m-4",
                "2014-2018 mean local d(rho_actual)/dz from daily rho files, levels 0 and 1"
            ),
            "rho_surface_level0": (
                rho_surface_level0, "kg m-3",
                "2014-2018 mean GSW potential density reconstructed from daily T/S at level 0"
            ),
            "rho_actual_level0": (
                rho_actual_level0, "kg m-3",
                "2014-2018 mean precomputed actual density from daily rho files at level 0"
            ),
            "rho_actual_level1": (
                rho_actual_level1, "kg m-3",
                "2014-2018 mean precomputed actual density from daily rho files at level 1"
            ),
            "alpha0": (
                alpha0, "kg m-3 K-1",
                "surface d(rho_pot)/dT used in APE generation"
            ),
            "beta0": (
                beta0, "kg m-3 psu-1",
                "surface d(rho_pot)/dS used in APE generation"
            ),
            "rho_star_mean": (
                rho_star_mean, "kg m-3",
                "mean surface potential density minus surface potential-density rho_ref"
            ),
            "J_mean": (
                j_mean, "K m s-1",
                "mean surface temperature flux H/(rho_pot*cp)"
            ),
            "Gs_mean": (
                gs_mean, "psu m s-1",
                "mean surface salinity flux Sbar*fr/rho_fw"
            ),
        }

        for name, (data, units, long_name) in auxiliary.items():
            v = ds.createVariable(
                name, "f4", ("lat", "lon"),
                fill_value=FILL_VALUE,
                zlib=True, complevel=DEFLATE_LEVEL
            )
            v[:] = np.ma.masked_invalid(data.astype(np.float32))
            v.units = units
            v.long_name = long_name

        long_names = {
            "gkm": "generation density of mean kinetic energy",
            "gke": "generation density of eddy kinetic energy from daily-resolved covariance",
            "gpm_heat": "mean APE generation density due to heat flux",
            "gpm_fw": "mean APE generation density due to freshwater flux",
            "gpm": "total mean APE generation density",
            "gpe_heat": "eddy APE generation density due to daily-resolved heat-flux covariance",
            "gpe_fw": "eddy APE generation density due to daily-resolved freshwater-flux covariance",
            "gpe": "total eddy APE generation density from daily-resolved covariance",
        }

        for name, data in fields.items():
            v = ds.createVariable(
                name, "f4", ("lat", "lon"),
                fill_value=FILL_VALUE,
                zlib=True, complevel=DEFLATE_LEVEL
            )
            v[:] = np.ma.masked_invalid(data.astype(np.float32))
            v.units = "W m-2"
            v.long_name = long_names[name]

        # Scalar integrated values as global attributes.
        for key, value in totals.items():
            ds.setncattr(f"{key}_W", float(value))
            ds.setncattr(f"{key}_TW", float(value) / 1.0e12)

        ds.title = (
            f"Daily-resolved Lorenz energy generation terms, "
            f"{START_YEAR}-{END_YEAR}"
        )
        ds.method = (
            "Surface generation formulas following von Storch et al. (2012); "
            "eddy covariances evaluated from daily fields."
        )
        ds.calendar = "365days"
        ds.valid_days = int(valid_day_count)
        ds.region_lon_min = float(LON_RANGE[0])
        ds.region_lon_max = float(LON_RANGE[1])
        ds.region_lat_min = float(LAT_RANGE[0])
        ds.region_lat_max = float(LAT_RANGE[1])
        ds.fr_units_used = "kg m-2 s-1"
        ds.fr_definition_used = "E - P_rain - P_snow - runoff"
        ds.freshwater_conversion = (
            f"Fw=fr/{RHO_FRESHWATER:g}; Gs=Sbar_surface*Fw"
        )
        ds.heat_flux_positive_direction = "into ocean"
        ds.heat_flux_expression = " + ".join(
            f"({HEAT_FLUX_SIGNS[k]:+g})*{k}"
            for k in ("sw", "lw", "sh", "lh")
        )
        ds.cp_seawater = CP_SEAWATER
        ds.rho_freshwater = RHO_FRESHWATER
        ds.gravity = G_GRAV
        ds.earth_radius = R_EARTH
        ds.ape_coefficient_mode = "gsw_alternative finite differences"
        ds.gsw_pressure_per_m = float(GSW_PRESSURE_PER_M)
        ds.gsw_pressure_note = (
            "Pressure convention follows user-supplied density script; verify "
            "custom gsw_alternative pressure units."
        )
        ds.potential_density_reference_pressure = float(POTENTIAL_DENSITY_REF_PRESSURE)
        ds.n0_method = (
            "n0 read from precomputed daily rho files: "
            "[rho_actual(level1)-rho_actual(level0)]/[z1-z0]; "
            "time mean local gradient then horizontal area mean"
        )
        ds.density_method = (
            "daily potential density reconstructed from ss/tt using "
            "gsw_SA_from_SP, gsw_CT_from_t, then gsw_rho at fixed reference pressure"
        )
        ds.density_type = "potential density"
        ds.daily_sampling_warning = (
            "G(Ke) and G(Pe) omit covariance at periods shorter than one day."
        )


def write_summary(totals, valid_days, n0_surface):
    OUTPUT_DAT.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append(
        f"Lorenz Energy Cycle generation, {START_YEAR}-{END_YEAR}, daily-resolved"
    )
    lines.append(f"Valid days: {valid_days}")
    lines.append(
        f"Region: lon=[{LON_RANGE[0]:g}, {LON_RANGE[1]:g}], "
        f"lat=[{LAT_RANGE[0]:g}, {LAT_RANGE[1]:g}]"
    )
    lines.append(
        f"n0 at surface (time+area mean d(rho_actual)/dz from daily rho files): "
        f"{n0_surface:.10e} kg m-4"
    )
    lines.append("")
    lines.append(f"{'term':<16s} {'W':>20s} {'TW':>16s}")
    lines.append("-" * 56)

    order = [
        "GKm",
        "GKe",
        "GPm_heat",
        "GPm_fw",
        "GPm",
        "GPe_heat",
        "GPe_fw",
        "GPe",
    ]

    for key in order:
        w = totals[key]
        lines.append(f"{key:<16s} {w:20.10e} {w/1.0e12:16.8f}")

    text = "\n".join(lines) + "\n"
    OUTPUT_DAT.write_text(text, encoding="utf-8")
    print("\n" + text)


# ============================================================
# 9. MAIN CALCULATION
# ============================================================

def main():
    print("=" * 80)
    print("Daily-resolved LEC generation: potential density for rho terms, rho files for n0")
    print(f"Period: {START_YEAR}-{END_YEAR}   calendar: 365days")
    print(f"Daily root: {DAILY_ROOT}")
    print(f"Mean root : {MEAN_ROOT}")
    print(
        "Thermodynamics: daily ss/tt -> potential density for rho terms; "
        "n0 from precomputed rho-YYYY-MM-DD.nc"
    )
    print(
        f"Pressure convention: p = {GSW_PRESSURE_PER_M:g} * (depth + z0) "
        f"when z0 is available"
    )
    print("=" * 80)

    # --------------------------------------------------------
    # 9.1 Grid and supplied mean fields
    # --------------------------------------------------------
    ss_mean_file = mean_path("ss")
    tt_mean_file = mean_path("tt")
    if not ss_mean_file.exists():
        raise FileNotFoundError(ss_mean_file)
    if not tt_mean_file.exists():
        raise FileNotFoundError(tt_mean_file)

    lon_full, lat_full = read_lon_lat(ss_mean_file)
    j_slice, i_slice = region_slices_from_coords(lon_full, lat_full)

    area_full = spherical_cell_area(lon_full, lat_full)
    lon = lon_full[i_slice]
    lat = lat_full[j_slice]
    area = area_full[j_slice, i_slice]
    shape2d = (lat.size, lon.size)
    lon2d, lat2d = np.meshgrid(lon, lat)

    print(
        f"Requested region: lon={LON_RANGE[0]:g}..{LON_RANGE[1]:g}, "
        f"lat={LAT_RANGE[0]:g}..{LAT_RANGE[1]:g}"
    )
    print(
        f"Selected grid centers: lon={lon[0]:.6f}..{lon[-1]:.6f}, "
        f"lat={lat[0]:.6f}..{lat[-1]:.6f}"
    )
    print(f"Regional grid: nx={lon.size}, ny={lat.size}")
    print(f"Regional rectangular area sum: {np.sum(area):.6e} m2")

    mean_fields = {}
    for key in ("tx", "ty", "uu", "vv"):
        path = mean_path(key)
        if not path.exists():
            raise FileNotFoundError(f"Required mean file not found: {path}")
        mean_fields[key] = read_surface_file(
            path, VAR_NAMES[key], j_slice, i_slice
        )
        check_shape(f"mean {key}", mean_fields[key], shape2d)

    ss_bar_file = read_surface_file(
        ss_mean_file, VAR_NAMES["ss"], j_slice, i_slice
    )
    tt_bar_file = read_surface_file(
        tt_mean_file, VAR_NAMES["tt"], j_slice, i_slice
    )
    check_shape("mean ss", ss_bar_file, shape2d)
    check_shape("mean tt", tt_bar_file, shape2d)

    # Potential-density reference profile used in rho*.
    # No n0 is calculated from GSW in this path.
    lev, rho_ref = build_reference_potential_density_profile_from_mean_ts(
        ss_mean_file=ss_mean_file,
        tt_mean_file=tt_mean_file,
        area=area,
        j_slice=j_slice,
        i_slice=i_slice,
        lon2d=lon2d,
        lat2d=lat2d,
    )
    zcoord = -np.abs(lev.astype(np.float64))

    # Full-depth n0 diagnostic profile from the original precomputed density.
    # The surface value used in G is later replaced by the time mean of DAILY
    # rho-YYYY-MM-DD.nc level-0/1 gradients.
    rho_mean_file = mean_path("rho")
    if not rho_mean_file.exists():
        raise FileNotFoundError(
            f"Required mean rho file for the diagnostic n0 profile not found: {rho_mean_file}"
        )

    (
        lev_rho,
        rho_ref_actual_for_n0,
        n0_profile_actual,
        n0_valid_count,
        _rho_actual_mean0_file,
        _rho_actual_mean1_file,
        _n0_surface_local_meanfile,
    ) = build_reference_density_and_n0_from_actual_density(
        rho_mean_file=rho_mean_file,
        rho_var=VAR_NAMES["rho"],
        area=area,
        j_slice=j_slice,
        i_slice=i_slice,
    )

    if lev_rho.shape != lev.shape or not np.allclose(lev_rho, lev, equal_nan=True):
        raise ValueError("rho and T/S vertical coordinates do not match.")

    # Mean z0 is optional; only affects the very small surface pressure correction.
    z0_mean_path = mean_path("z0")
    if USE_Z0_IN_PRESSURE and z0_mean_path.exists():
        z0_bar = read_surface_file(
            z0_mean_path, VAR_NAMES["z0"], j_slice, i_slice
        )
        check_shape("mean z0", z0_bar, shape2d)
    else:
        z0_bar = np.zeros(shape2d, dtype=np.float64)

    p_surface_mean = pressure_at_level(lev[0], shape2d, z0_bar)

    # Potential-density-consistent alpha0/beta0 from supplied 2014-2018 mean T/S.
    alpha0, beta0 = gsw_surface_derivatives_from_SP_t(
        ss_bar_file, tt_bar_file, p_surface_mean, lon2d, lat2d
    )

    print(
        "alpha0=d(rho_pot)/dT finite range: "
        f"{np.nanmin(alpha0):.6e} .. {np.nanmax(alpha0):.6e} kg m-3 K-1"
    )
    print(
        "beta0=d(rho_pot)/dS finite range : "
        f"{np.nanmin(beta0):.6e} .. {np.nanmax(beta0):.6e} kg m-3 psu-1"
    )

    # --------------------------------------------------------
    # 9.2 Accumulators
    # --------------------------------------------------------
    z = np.zeros(shape2d, dtype=np.float64)
    iz = np.zeros(shape2d, dtype=np.int32)

    # Wind means/products.
    sum_tx = z.copy()
    sum_ty = z.copy()
    sum_u = z.copy()
    sum_v = z.copy()
    sum_txu = z.copy()
    sum_tyv = z.copy()
    count_wind = iz.copy()

    # Reconstructed surface T/S/potential-density means.
    sum_ss0 = z.copy()
    sum_tt0 = z.copy()
    sum_rho0 = z.copy()
    count_ts0 = iz.copy()

    # Precomputed actual density used ONLY for n0.
    sum_rho_actual0 = z.copy()
    sum_rho_actual1 = z.copy()
    count_rho_actual0 = iz.copy()
    count_rho_actual1 = iz.copy()

    # Daily actual-density-gradient mean field used for n0.
    sum_n0_actual = z.copy()
    count_n0_actual = iz.copy()

    # Heat-flux/rho statistics.
    sum_rho_heat = z.copy()
    sum_J = z.copy()
    sum_rhoJ = z.copy()
    count_heat = iz.copy()

    sum_H = z.copy()
    sum_J_component = {key: z.copy() for key in ("sw", "lw", "sh", "lh")}
    sum_rhoJ_component = {key: z.copy() for key in ("sw", "lw", "sh", "lh")}

    # Freshwater/salinity-flux/rho statistics.
    sum_rho_fw = z.copy()
    sum_Gs = z.copy()
    sum_rhoGs = z.copy()
    count_fw = iz.copy()

    total_calendar_days = sum(DAYS_IN_MONTH) * (END_YEAR - START_YEAR + 1)
    valid_days = 0
    heat_period_diag = init_heat_period_diagnostics()
    first_checked = False
    warned_missing_z0 = False

    required_daily = (
        "fr", "lw", "sw", "sh", "lh",
        "tx", "ty", "uu", "vv", "ss", "tt", "rho"
    )

    # --------------------------------------------------------
    # 9.3 Daily loop
    # --------------------------------------------------------
    for iday, (year, month, day) in enumerate(iter_365_dates(), start=1):
        paths = {
            key: daily_path(key, year, month, day)
            for key in required_daily
        }

        missing = [p for p in paths.values() if not p.exists()]
        if missing:
            msg = (
                f"Missing files on {year:04d}-{month:02d}-{day:02d}: "
                + ", ".join(str(p) for p in missing)
            )
            if SKIP_MISSING_DAYS:
                print("SKIP:", msg)
                continue
            raise FileNotFoundError(msg)

        if valid_days == 0 and PRINT_HEAT_FILE_METADATA:
            print_heat_file_metadata(paths)

        # Surface forcing/dynamic fields.
        fr = read_surface_file(paths["fr"], VAR_NAMES["fr"], j_slice, i_slice)
        lw = read_surface_file(paths["lw"], VAR_NAMES["lw"], j_slice, i_slice)
        sw = read_surface_file(paths["sw"], VAR_NAMES["sw"], j_slice, i_slice)
        sh = read_surface_file(paths["sh"], VAR_NAMES["sh"], j_slice, i_slice)
        lh = read_surface_file(paths["lh"], VAR_NAMES["lh"], j_slice, i_slice)
        tx = read_surface_file(paths["tx"], VAR_NAMES["tx"], j_slice, i_slice)
        ty = read_surface_file(paths["ty"], VAR_NAMES["ty"], j_slice, i_slice)
        uu = read_surface_file(paths["uu"], VAR_NAMES["uu"], j_slice, i_slice)
        vv = read_surface_file(paths["vv"], VAR_NAMES["vv"], j_slice, i_slice)

        # Surface T/S is used to reconstruct POTENTIAL DENSITY for rho*, rho'J'
        # and rho'Gs'.  The first two PRECOMPUTED rho levels are read separately
        # and are used ONLY to calculate n0.
        with nc.Dataset(paths["ss"], "r") as dss,              nc.Dataset(paths["tt"], "r") as dst,              nc.Dataset(paths["rho"], "r") as dsr:
            ss0 = read_level_from_open_dataset(
                dss, VAR_NAMES["ss"], 0, j_slice, i_slice
            )
            tt0 = read_level_from_open_dataset(
                dst, VAR_NAMES["tt"], 0, j_slice, i_slice
            )
            rho_actual0 = read_level_from_open_dataset(
                dsr, VAR_NAMES["rho"], 0, j_slice, i_slice
            )
            rho_actual1 = read_level_from_open_dataset(
                dsr, VAR_NAMES["rho"], 1, j_slice, i_slice
            )

        for name, a in (
            ("fr", fr), ("lw", lw), ("sw", sw), ("sh", sh), ("lh", lh),
            ("tx", tx), ("ty", ty), ("uu", uu), ("vv", vv),
            ("ss0", ss0), ("tt0", tt0),
            ("rho_actual0", rho_actual0), ("rho_actual1", rho_actual1),
        ):
            check_shape(name, a, shape2d)

        # Optional daily z0, following the supplied density-generation script.
        z0_daily = None
        if USE_Z0_IN_PRESSURE:
            z0_path = daily_path("z0", year, month, day)
            if z0_path.exists():
                z0_daily = read_surface_file(
                    z0_path, VAR_NAMES["z0"], j_slice, i_slice
                )
                check_shape("z0", z0_daily, shape2d)
            else:
                if not warned_missing_z0:
                    print(
                        "WARNING: z0 daily file not found; pressure will use "
                        "z0=0. This warning is printed only once."
                    )
                    warned_missing_z0 = True

        p0 = pressure_at_level(lev[0], shape2d, z0_daily)

        # Potential density used in the APE density/covariance terms.
        SA0, CT0, rho0 = gsw_state_from_SP_t(
            ss0, tt0, p0, lon2d, lat2d
        )

        # n0 is calculated ONLY from the original precomputed actual-density files.
        n0_actual_daily = actual_density_gradient_two_levels(
            rho_actual0, rho_actual1, zcoord[0], zcoord[1]
        )

        if not first_checked:
            lon_d_full, lat_d_full = read_lon_lat(paths["tx"])
            if (
                lon_d_full.shape != lon_full.shape
                or lat_d_full.shape != lat_full.shape
                or not np.allclose(lon_d_full, lon_full, equal_nan=True)
                or not np.allclose(lat_d_full, lat_full, equal_nan=True)
            ):
                raise ValueError("Daily-grid lon/lat do not match the mean T/S grid.")
            first_checked = True

        # Total heat flux H > 0 into ocean.
        H = (
            HEAT_FLUX_SIGNS["sw"] * sw
            + HEAT_FLUX_SIGNS["lw"] * lw
            + HEAT_FLUX_SIGNS["sh"] * sh
            + HEAT_FLUX_SIGNS["lh"] * lh
        )
        H_all_stored = sw + lw + sh + lh
        H_sw_minus_outgoing = sw - lw - sh - lh

        heat_ocean_mask = valid_density_mask(rho0)

        for diag_name, diag_field in (
            ("sw", sw),
            ("lw", lw),
            ("sh", sh),
            ("lh", lh),
            ("H_configured", H),
            ("H_all_stored", H_all_stored),
            ("H_sw_minus_outgoing", H_sw_minus_outgoing),
        ):
            update_heat_period_diagnostic(
                heat_period_diag, diag_name, diag_field, area, heat_ocean_mask
            )

        # Temperature flux J = H/(rho_s cp).
        J = np.full(shape2d, np.nan, dtype=np.float64)
        valid_j = np.isfinite(H) & valid_density_mask(rho0)
        J[valid_j] = H[valid_j] / (
            rho0[valid_j] * CP_SEAWATER
        )

        heat_components = {"sw": sw, "lw": lw, "sh": sh, "lh": lh}
        J_components = {}
        for key, q in heat_components.items():
            jq = np.full(shape2d, np.nan, dtype=np.float64)
            valid_q = np.isfinite(q) & valid_density_mask(rho0)
            jq[valid_q] = (
                HEAT_FLUX_SIGNS[key] * q[valid_q]
                / (rho0[valid_q] * CP_SEAWATER)
            )
            J_components[key] = jq

        # Freshwater salinity flux uses the time-mean uppermost salinity,
        # exactly as in the von Storch/Koldunov formulation.
        Fw = fr / RHO_FRESHWATER
        Gs = ss_bar_file * Fw

        # ---- Wind statistics ----
        valid_w = finite_common(tx, ty, uu, vv)
        sum_tx[valid_w] += tx[valid_w]
        sum_ty[valid_w] += ty[valid_w]
        sum_u[valid_w] += uu[valid_w]
        sum_v[valid_w] += vv[valid_w]
        sum_txu[valid_w] += tx[valid_w] * uu[valid_w]
        sum_tyv[valid_w] += ty[valid_w] * vv[valid_w]
        count_wind[valid_w] += 1

        # ---- T/S and reconstructed density statistics ----
        valid_ts0 = finite_common(ss0, tt0, rho0)
        sum_ss0[valid_ts0] += ss0[valid_ts0]
        sum_tt0[valid_ts0] += tt0[valid_ts0]
        sum_rho0[valid_ts0] += rho0[valid_ts0]
        count_ts0[valid_ts0] += 1

        valid_ra0 = valid_density_mask(rho_actual0)
        sum_rho_actual0[valid_ra0] += rho_actual0[valid_ra0]
        count_rho_actual0[valid_ra0] += 1

        valid_ra1 = valid_density_mask(rho_actual1)
        sum_rho_actual1[valid_ra1] += rho_actual1[valid_ra1]
        count_rho_actual1[valid_ra1] += 1

        valid_n0 = np.isfinite(n0_actual_daily)
        sum_n0_actual[valid_n0] += n0_actual_daily[valid_n0]
        count_n0_actual[valid_n0] += 1

        # ---- Heat/rho statistics ----
        valid_h = finite_common(rho0, J)
        sum_rho_heat[valid_h] += rho0[valid_h]
        sum_J[valid_h] += J[valid_h]
        sum_rhoJ[valid_h] += rho0[valid_h] * J[valid_h]
        count_heat[valid_h] += 1

        sum_H[valid_h] += H[valid_h]
        for key in ("sw", "lw", "sh", "lh"):
            jq = J_components[key]
            valid_q = finite_common(rho0, jq)
            sum_J_component[key][valid_q] += jq[valid_q]
            sum_rhoJ_component[key][valid_q] += rho0[valid_q] * jq[valid_q]

        # ---- Freshwater/rho statistics ----
        valid_f = finite_common(rho0, Gs)
        sum_rho_fw[valid_f] += rho0[valid_f]
        sum_Gs[valid_f] += Gs[valid_f]
        sum_rhoGs[valid_f] += rho0[valid_f] * Gs[valid_f]
        count_fw[valid_f] += 1

        valid_days += 1

        if (
            valid_days == 1
            or valid_days % PRINT_INTERVAL == 0
            or iday == total_calendar_days
        ):
            print(
                f"[{iday:4d}/{total_calendar_days}] "
                f"{year:04d}-{month:02d}-{day:02d}  "
                f"valid_days={valid_days}"
            )

            if valid_days == 1:
                if PRINT_HEAT_FIRST_DAY_DIAGNOSTICS:
                    print_one_day_heat_diagnostics(
                        year=year,
                        month=month,
                        day=day,
                        sw=sw,
                        lw=lw,
                        sh=sh,
                        lh=lh,
                        H_configured=H,
                        area=area,
                        ocean_mask=heat_ocean_mask,
                    )

                fr_stats = heat_weighted_stats(fr, area, heat_ocean_mask)
                print(
                    "First-day freshwater flux fr [kg m-2 s-1]: "
                    f"min={fr_stats['min']:.6e}, "
                    f"area_mean={fr_stats['mean']:.6e}, "
                    f"max={fr_stats['max']:.6e}"
                )

                print(
                    "First-day reconstructed GSW surface potential density: "
                    f"min={np.nanmin(rho0):.6f}, "
                    f"mean={np.nanmean(rho0):.6f}, "
                    f"max={np.nanmax(rho0):.6f} kg m-3"
                )
                print(
                    "First-day n0 from precomputed rho files: "
                    f"min={np.nanmin(n0_actual_daily):.6e}, "
                    f"mean={np.nanmean(n0_actual_daily):.6e}, "
                    f"max={np.nanmax(n0_actual_daily):.6e} kg m-4"
                )

    if valid_days == 0:
        raise RuntimeError("No daily files were processed.")

    print(f"Processed {valid_days} valid days.")

    if PRINT_HEAT_PERIOD_DIAGNOSTICS:
        print_heat_period_diagnostics(heat_period_diag, valid_days)

    # --------------------------------------------------------
    # 9.4 Means and n0 from precomputed rho files
    # --------------------------------------------------------
    tx_bar_daily = safe_mean(sum_tx, count_wind)
    ty_bar_daily = safe_mean(sum_ty, count_wind)
    u_bar_daily = safe_mean(sum_u, count_wind)
    v_bar_daily = safe_mean(sum_v, count_wind)

    ss_bar_daily = safe_mean(sum_ss0, count_ts0)
    tt_bar_daily = safe_mean(sum_tt0, count_ts0)

    rho_bar_heat_daily = safe_mean(sum_rho_heat, count_heat)
    rho_bar_fw_daily = safe_mean(sum_rho_fw, count_fw)

    rho_actual0_bar_daily = safe_mean(
        sum_rho_actual0, count_rho_actual0
    )
    rho_actual1_bar_daily = safe_mean(
        sum_rho_actual1, count_rho_actual1
    )

    n0_surface_local_mean = safe_mean(
        sum_n0_actual, count_n0_actual
    )
    valid_n0_mean = np.isfinite(n0_surface_local_mean)
    n0_surface = area_weighted_mean_with_mask(
        n0_surface_local_mean, area, valid_n0_mean
    )

    if not np.isfinite(n0_surface) or abs(n0_surface) < 1.0e-12:
        raise RuntimeError(
            f"Invalid/too-small surface n0 from daily rho files = {n0_surface}."
        )

    # For rho* use the potential-density field reconstructed from daily T/S.
    # The regional reference density is its horizontal area mean.
    rho_ref_surface = area_weighted_mean_with_mask(
        rho_bar_heat_daily,
        area,
        valid_density_mask(rho_bar_heat_daily),
    )
    rho_ref[0] = rho_ref_surface

    # n0 profile comes from the PRECOMPUTED actual-density mean file.
    # The surface value is replaced by the time mean of DAILY rho-file gradients,
    # which is the value actually used in G.
    n0 = n0_profile_actual.copy()
    n0[0] = n0_surface

    rho_star_mean = rho_bar_heat_daily - rho_ref_surface

    print_n0_diagnostics(
        lon=lon,
        lat=lat,
        lev=lev,
        rho_ref=rho_ref_actual_for_n0,
        n0=n0,
        n0_valid_count=n0_valid_count,
        rho0=rho_actual0_bar_daily,
        rho1=rho_actual1_bar_daily,
        n0_surface_local=n0_surface_local_mean,
        area=area,
    )

    print(
        f"Surface rho_ref from reconstructed daily GSW potential density = "
        f"{rho_ref_surface:.8f} kg m-3"
    )
    print(
        "Surface n0 used in G, from daily precomputed rho files = "
        f"{n0_surface:.8e} kg m-4"
    )

    print("\nMean-field consistency: daily reconstruction vs supplied mean files")
    for name, daily_bar, supplied_bar in (
        ("tx", tx_bar_daily, mean_fields["tx"]),
        ("ty", ty_bar_daily, mean_fields["ty"]),
        ("uu", u_bar_daily, mean_fields["uu"]),
        ("vv", v_bar_daily, mean_fields["vv"]),
        ("ss", ss_bar_daily, ss_bar_file),
        ("tt", tt_bar_daily, tt_bar_file),
    ):
        print(
            f"  {name:<15s} "
            f"RMS diff={rms_difference(daily_bar, supplied_bar):.6e}, "
            f"max abs diff={max_abs_difference(daily_bar, supplied_bar):.6e}"
        )

    # Wind/velocity Reynolds means can use supplied means.
    if USE_PROVIDED_MEANS:
        tx_bar = mean_fields["tx"]
        ty_bar = mean_fields["ty"]
        u_bar = mean_fields["uu"]
        v_bar = mean_fields["vv"]
        print("Wind/velocity covariance decomposition uses supplied 2014-2018 means.")
    else:
        tx_bar = tx_bar_daily
        ty_bar = ty_bar_daily
        u_bar = u_bar_daily
        v_bar = v_bar_daily
        print("Wind/velocity covariance decomposition uses daily reconstructed means.")

    # Density covariance MUST use the same daily reconstructed potential-density samples.
    rho_bar_heat = rho_bar_heat_daily
    rho_bar_fw = rho_bar_fw_daily
    print(
        "Density covariance decomposition uses GSW potential density reconstructed "
        "from the same daily tt/ss samples."
    )

    # Product means.
    txu_bar = safe_mean(sum_txu, count_wind)
    tyv_bar = safe_mean(sum_tyv, count_wind)

    J_bar = safe_mean(sum_J, count_heat)
    rhoJ_bar = safe_mean(sum_rhoJ, count_heat)
    H_bar = safe_mean(sum_H, count_heat)

    component_J_bar = {}
    component_rhoJ_bar = {}
    for key in ("sw", "lw", "sh", "lh"):
        component_J_bar[key] = safe_mean(sum_J_component[key], count_heat)
        component_rhoJ_bar[key] = safe_mean(sum_rhoJ_component[key], count_heat)
    component_rhoJ_bar["total"] = rhoJ_bar

    Gs_bar = safe_mean(sum_Gs, count_fw)
    rhoGs_bar = safe_mean(sum_rhoGs, count_fw)

    # --------------------------------------------------------
    # 9.5 G(Km), G(Ke)
    # --------------------------------------------------------
    gkm = tx_bar * u_bar + ty_bar * v_bar
    cov_txu = txu_bar - tx_bar * u_bar
    cov_tyv = tyv_bar - ty_bar * v_bar
    gke = cov_txu + cov_tyv

    # --------------------------------------------------------
    # 9.6 G(Pm), G(Pe)
    # --------------------------------------------------------
    factor_t = -G_GRAV * alpha0 / n0_surface
    factor_s = -G_GRAV * beta0 / n0_surface

    gpm_heat = factor_t * J_bar * rho_star_mean
    gpm_fw = factor_s * Gs_bar * rho_star_mean
    gpm = gpm_heat + gpm_fw

    cov_rhoJ = rhoJ_bar - rho_bar_heat * J_bar
    cov_rhoGs = rhoGs_bar - rho_bar_fw * Gs_bar

    gpe_heat = factor_t * cov_rhoJ
    gpe_fw = factor_s * cov_rhoGs
    gpe = gpe_heat + gpe_fw

    if PRINT_HEAT_G_DIAGNOSTICS:
        print_heat_G_diagnostics(
            area=area,
            rho_star_mean=rho_star_mean,
            rho_bar_file=rho_bar_heat_daily,
            alpha0=alpha0,
            n0_surface=n0_surface,
            J_bar=J_bar,
            H_bar=H_bar,
            component_J_bar=component_J_bar,
            component_rhoJ_bar=component_rhoJ_bar,
            rho_bar_heat=rho_bar_heat,
        )

    km_mask = finite_common(tx_bar, ty_bar, u_bar, v_bar)
    ape_mask = finite_common(
        rho_bar_heat_daily, alpha0, beta0, rho_star_mean
    )

    gkm[~km_mask] = np.nan
    gke[~km_mask] = np.nan
    for a in (gpm_heat, gpm_fw, gpm, gpe_heat, gpe_fw, gpe):
        a[~ape_mask] = np.nan

    fields = {
        "gkm": gkm,
        "gke": gke,
        "gpm_heat": gpm_heat,
        "gpm_fw": gpm_fw,
        "gpm": gpm,
        "gpe_heat": gpe_heat,
        "gpe_fw": gpe_fw,
        "gpe": gpe,
    }

    totals = {
        "GKm": integrate_w(gkm, area),
        "GKe": integrate_w(gke, area),
        "GPm_heat": integrate_w(gpm_heat, area),
        "GPm_fw": integrate_w(gpm_fw, area),
        "GPm": integrate_w(gpm, area),
        "GPe_heat": integrate_w(gpe_heat, area),
        "GPe_fw": integrate_w(gpe_fw, area),
        "GPe": integrate_w(gpe, area),
    }

    # --------------------------------------------------------
    # 9.7 Output
    # --------------------------------------------------------
    write_output(
        lon=lon,
        lat=lat,
        lev=lev,
        area=area,
        rho_ref=rho_ref,
        n0=n0,
        n0_surface_local=n0_surface_local_mean,
        rho_surface_level0=rho_bar_heat_daily,
        rho_actual_level0=rho_actual0_bar_daily,
        rho_actual_level1=rho_actual1_bar_daily,
        alpha0=alpha0,
        beta0=beta0,
        rho_star_mean=rho_star_mean,
        j_mean=J_bar,
        gs_mean=Gs_bar,
        fields=fields,
        totals=totals,
        valid_day_count=valid_days,
    )

    write_summary(
        totals=totals,
        valid_days=valid_days,
        n0_surface=n0_surface,
    )

    print(f"NetCDF output: {OUTPUT_NC}")
    print(f"Text output  : {OUTPUT_DAT}")
    print("\nIMPORTANT:")
    print("  * n0 is calculated from the original precomputed daily rho files.")
    print("  * rho_pot/rho* are still reconstructed from daily tt/ss with gsw_alternative.")
    print("  * Only n0 changed; the potential-density APE covariance calculation is retained.")
    print("  * G(Ke) and G(Pe) are DAILY-resolved covariances and omit sub-daily covariance.")


if __name__ == "__main__":
    main()
