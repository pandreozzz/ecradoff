"""Test interp2d"""

import os
import sys

from time import time

from pathlib import Path
import traceback
import xarray as xr
import numpy as np

from src.interp.wrapper import GriddedArray

verbose=False
aero_clim_fpath = Path(os.path.join("data/cams_clim", "aerosol_cams_climatology_4D_test.nc"))

dclim = xr.open_dataset(aero_clim_fpath).isel(epoch=[0,1], month=[0,1,2])

print("Loaded dataset:")
print(dclim.__str__())
ga_dclim = GriddedArray(dclim)

dgrid = 0.25
tgt_grid1 = xr.Dataset(
    coords={
        "lon": ("lon", np.arange(0, 360, dgrid)),
        "lat": ("lat", np.arange(90, -90, -dgrid))
    }
).coords

nredpts = 180
lats = np.linspace(-90, 90, nredpts)
reduced_pts = (np.cos(np.deg2rad(lats))*nredpts*2).clip(min=4).astype(int)
tgt_grid2 = xr.Dataset(
    coords={
        "lon": ("lat", np.zeros_like(lats)),
        "lat": ("lat", lats),
        "reduced_points": ("lat", reduced_pts)
    }
).coords


nunstruct = 1000
lon = np.random.uniform(0, 360, nunstruct)
lat = np.random.uniform(-90, 90, nunstruct)
tgt_grid3 = xr.Dataset(
    coords={
        "lon": ("col", lon),
        "lat": ("col", lat)
    }
).coords

nredpts = 60
lats = np.linspace(-90, 90, nredpts)
reduced_pts = (np.cos(np.deg2rad(lats))*nredpts*2).clip(min=4).astype(int)
tgt_grid4 = xr.Dataset(
    coords={
        "lon": ("lat", np.zeros_like(lats)),
        "lat": ("lat", lats),
        "reduced_points": ("lat", reduced_pts)
    }
).coords


nunstruct = 1300
lon = np.random.uniform(0, 360, nunstruct)
lat = np.random.uniform(-90, 90, nunstruct)
tgt_grid5 = xr.Dataset(
    coords={
        "lon": ("col", lon),
        "lat": ("col", lat)
    }
).coords


print("Testing regular grid interpolation...")
start = time()
try:
    ga_interp1 = ga_dclim.interp2d_to(tgt_grid1, verbose=verbose)
    print(ga_interp1.array.__str__(), flush=True)
except Exception as e:
    print(f"Regular grid interpolation failed with error: {e}")
    traceback.print_exc()
print("Regular grid interpolation took {:.2f} seconds".format(time()-start))

print("Testing reduced grid interpolation...")
start = time()
try:
    ga_interp2 = ga_dclim.interp2d_to(tgt_grid2, verbose=verbose)
    print(ga_interp2.array.__str__(), flush=True)
except Exception as e:  
    print(f"Reduced grid interpolation failed with error: {e}")
    traceback.print_exc()
    ga_interp2 = None
print("Reduced grid interpolation took {:.2f} seconds".format(time()-start))

print("Testing unstructured grid interpolation...")
start = time()
try:
    ga_interp3 = ga_dclim.interp2d_to(tgt_grid3, verbose=verbose)
    print(ga_interp3.array.__str__(), flush=True)
except Exception as e:
    print(f"Unstructured grid interpolation failed with error: {e}")
    traceback.print_exc()
    ga_interp3 = None
print("Unstructured grid interpolation took {:.2f} seconds".format(time()-start))

if ga_interp2 is not None:
    start = time()
    print("Testing reduced to regular grid interpolation")
    try:
        ga_interp4 = ga_interp2.interp2d_to(tgt_grid1, verbose=verbose)
        print(ga_interp4.array.__str__(), flush=True)
    except Exception as e:
        print(f"Reduced to regular grid interpolation failed with error: {e}")
        traceback.print_exc()
    print("Reduced to regular grid interpolation took {:.2f} seconds".format(time()-start))

    start = time()
    print("Testing reduced to reduced grid interpolation")
    try:
        ga_interp4 = ga_interp2.interp2d_to(tgt_grid4, verbose=verbose)
        print(ga_interp4.array.__str__(), flush=True)
    except Exception as e:
        print(f"Reduced to reduced grid interpolation failed with error: {e}")
        traceback.print_exc()
    print("Reduced to reduced grid interpolation took {:.2f} seconds".format(time()-start))

    start = time()
    print("Testing reduced to unstructured grid interpolation")
    try:
        ga_interp5 = ga_interp2.interp2d_to(tgt_grid5, verbose=verbose)
        print(ga_interp5.array.__str__(), flush=True)
    except Exception as e:
        print(f"Reduced to unstructured grid interpolation failed with error: {e}")
        traceback.print_exc()
    print("Reduced to unstructured grid interpolation took {:.2f} seconds".format(time()-start))

# if ga_interp3 is not None:
#     print("Testing unstructured to regular grid interpolation")
#     try:
#         ga_interp6 = ga_interp3.interp2d_to(tgt_grid1, verbose=verbose)
#         print(ga_interp6.array.__str__(), flush=True)
#     except Exception as e:
#         print(f"Unstructured to regular grid interpolation failed with error: {e}")
#         traceback.print_exc()
    
#     print("Testing unstructured to reduced grid interpolation")
#     try:
#         ga_interp7 = ga_interp3.interp2d_to(tgt_grid2, verbose=verbose)
#         print(ga_interp7.array.__str__(), flush=True)
#     except Exception as e:
#         print(f"Unstructured to reduced grid interpolation failed with error: {e}")
#         traceback.print_exc()
    
#     print("Testing unstructured to unstructured grid interpolation")
#     try:
#         ga_interp8 = ga_interp3.interp2d_to(tgt_grid5, verbose=verbose)
#         print(ga_interp8.array.__str__(), flush=True)
#     except Exception as e:
#         print(f"Unstructured to unstructured grid interpolation failed with error: {e}")
#         traceback.print_exc()