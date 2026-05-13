"""Test interpvert via GriddedProfile.interp3d_to."""

import os
from time import time
from pathlib import Path
import traceback

import numpy as np
import xarray as xr

from src.interp.wrapper import GriddedProfile

verbose = False
aero_clim_fpath = Path(os.path.join("data/cams_clim", "aerosol_cams_climatology_4D_test.nc"))

# Keep this test reasonably small while preserving realistic dimensions.
dclim = xr.open_dataset(aero_clim_fpath).isel(epoch=[0, 1], month=[0, 1, 2])

print("Loaded dataset:")
print(dclim.__str__())

gp_dclim = GriddedProfile(dclim, profile_coord="pressure", lev_dim="lev")

# 1D target pressure levels in Pa: 100, 200, ..., 102000
plevs = np.logspace(2, np.log10(102001), 30, dtype=np.float64)


# 1D tests with different level-dimension names.
for lev_name in ["plev", "lev", "level"]:
    print(f"Testing 1D vertical interpolation with dim '{lev_name}'...")
    start = time()
    try:
        ptgt_1d = xr.DataArray(
            plevs,
            dims=(lev_name,),
            coords={lev_name: plevs},
            name="ptgt",
        )
        gp_interp = gp_dclim.interp3d_to(
            ptgt=ptgt_1d,
            tgt_coords=None,
            lev_dim_tgt=lev_name,
            verbose=verbose,
        )
        print(gp_interp.profile.array.__str__(), flush=True)
    except Exception as exc:
        print(f"1D interpolation for dim '{lev_name}' failed with error: {exc}")
        traceback.print_exc()
    print(
        f"1D interpolation with dim '{lev_name}' took {time() - start:.2f} seconds"
    )


print("Testing 3D vertical interpolation on rectangular grid (lat, lon, lev)...")
start = time()
try:
    lat_rect = np.arange(90.0, -90.0, -5.0)
    lon_rect = np.arange(0.0, 360.0, 5.0)

    # Same vertical profile at each horizontal point.
    ptgt_rect_vals = np.broadcast_to(
        plevs[np.newaxis, np.newaxis, :],
        (lat_rect.size, lon_rect.size, plevs.size),
    ).copy()

    ptgt_rect = xr.DataArray(
        data=ptgt_rect_vals,
        dims=("lat", "lon", "lev"),
        coords={
            "lat": ("lat", lat_rect),
            "lon": ("lon", lon_rect),
            "lev": ("lev", plevs),
        },
        name="ptgt",
    )

    gp_interp_rect = gp_dclim.interp3d_to(
        ptgt=ptgt_rect,
        tgt_coords=ptgt_rect.coords,
        lev_dim_tgt="lev",
        verbose=verbose,
    )
    print(gp_interp_rect.profile.array.__str__(), flush=True)
except Exception as exc:
    print(f"3D rectangular-grid interpolation failed with error: {exc}")
    traceback.print_exc()
print("3D rectangular-grid interpolation took {:.2f} seconds".format(time() - start))


print("Testing 3D vertical interpolation on reduced grid (col, lev)...")
start = time()
try:
    nredpts = 30
    lats_red = np.linspace(-90.0, 90.0, nredpts)
    reduced_pts = (np.cos(np.deg2rad(lats_red)) * nredpts * 2).clip(min=4).astype(int)
    ncol = int(np.sum(reduced_pts))

    # Same vertical profile at each reduced-grid column.
    ptgt_red_vals = np.broadcast_to(
        plevs[np.newaxis, :],
        (ncol, plevs.size),
    ).copy()

    ds_ptgt_red = xr.Dataset(
        data_vars={"ptgt" : (("col","lev"), ptgt_red_vals)},
        coords={
            "col": ("col", np.arange(ncol, dtype=np.int64)),
            "lev": ("lev", plevs),
            "lat": ("lat", lats_red),
            "lon": ("lat", np.zeros_like(lats_red)),
            "reduced_points": ("lat", reduced_pts),
        }
    )

    gp_interp_red = gp_dclim.interp3d_to(
        ptgt=ds_ptgt_red["ptgt"],
        tgt_coords=ds_ptgt_red.coords,
        tgt_datadims=["col"],
        lev_dim_tgt="lev",
        verbose=verbose,
    )
    print(gp_interp_red.profile.array.__str__(), flush=True)
except Exception as exc:
    print(f"3D reduced-grid interpolation failed with error: {exc}")
    traceback.print_exc()
print("3D reduced-grid interpolation took {:.2f} seconds".format(time() - start))
