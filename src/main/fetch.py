
import os

from typing import List, Optional, Union
from pathlib import Path

import xarray as xr
import numpy as np

from ..interp.wrapper import XRCoords
from .config import HOR_DIM_1D, TIME_DIM, LEV_DIM, HLEV_DIM, USE_DASK

from ..aerodefs.optics import OpticsData

from .types import FLOAT_DTYPE

RENAMEDIC = {
    "longitude" : "lon",
    "latitude" : "lat",
    "level" : LEV_DIM,
    "plev" : LEV_DIM,
    "rgrid" : HOR_DIM_1D, # For cdo nc4
}

#ReducedLonLat = namedtuple("redlonlat", ["lon", "lat"])
def gen_reduced_lonlat(ds : XRCoords,
                       dim_1d : str = HOR_DIM_1D,
                       rpoint_coord : str = "reduced_points"):
    from ..interp.wrapper import get_arraygriddesc

    newlon = np.array([])
    newlat = np.array([])
    for rpn, lat in zip(ds[rpoint_coord].values, ds["lat"].values):
        newlon = np.concat([newlon, np.arange(0, 360, 360/rpn)[0:rpn]], axis=0)
        newlat = np.concat([newlat, [lat]*rpn], axis=0)

    grid_coords = xr.Dataset(
        coords = {
            "lon": ([dim_1d], newlon),
            "lat": ([dim_1d], newlat)
        }
    ).coords

    return get_arraygriddesc(grid_coords=grid_coords,
                             coordims=[dim_1d])

def cast_floats(ds : xr.Dataset) -> xr.Dataset:
    """Casts all float variables in the dataset to FLOAT_DTYPE."""
    
    caster = {}


    for name, da in ds.data_vars.items():
        if np.issubdtype(da.dtype, np.floating):
            caster[name] = FLOAT_DTYPE
        else:
            caster[name] = da.dtype
    
    return ds.astype(caster, copy=False)

def _full_rename(ds : xr.Dataset) -> xr.Dataset:
    """Renames variables, coordinates, and dimensions in the dataset according to RENAMEDIC."""
    rename_dict = {k:v for k, v in RENAMEDIC.items()
                   if k in ds.data_vars or k in ds.coords or k in ds.dims}
    return ds.rename(rename_dict)

def model_fields(model_files: Union[str, List[str]],
                 time_sel: Union[None, List[np.datetime64], List[int]] = None):
    from .config import TIME_DIM, USE_DASK
    from .types import NS_DTTYPE
    from ..tools.ifs import populate_mlfields, compute_sp_from_lnsp

    from glob import glob
    if not isinstance(model_files, list):
        model_files = [model_files]
    all_files = []
    for path in model_files:
        all_files.extend(glob(os.path.realpath(path)))
    all_files = list(set(all_files))

    if all_files == []:
        raise ValueError("No files found for the given output_path(s).")

    opends_kwargs = {"combine": "by_coords", "parallel":False}
    # Are they zarr archives?
    if all_files[0].endswith("zarr"):
        print("Opening archives as zarr")
        opends_kwargs["engine"] = "zarr"
    else:
        print("Opening files as netCDF")
        opends_kwargs["engine"] = "netcdf4"

    model_fields = xr.open_mfdataset(all_files, **opends_kwargs) #type: ignore

    if time_sel is None or len(time_sel) == 0:
        print("Using all model times")
        times = model_fields[TIME_DIM].values.astype(NS_DTTYPE)
    else:
        if isinstance(time_sel[0], int):
            print("Using model times at input indices:")
            times = model_fields[TIME_DIM].isel({TIME_DIM: time_sel}).values.astype(NS_DTTYPE)
        else:
            print("Using model times at input datetimes:")
            times = time_sel
    
    # Sanity check
    assert isinstance(times[0], np.datetime64)

    model_fields.sel({TIME_DIM: times}, method="nearest").assign_coords({TIME_DIM: times})

    # Rename to standard names if needed
    model_fields = _full_rename(model_fields)

    # if "lon" not in model_fields:
    #     model_fields = gen_reduced_lon(model_fields)
    if "reduced_points" in model_fields:
        model_fields["reduced_points"] = xr.DataArray(
            data=model_fields["reduced_points"].values,
            dims=["lat"]
            )
        
    # Populate ML and HL fields
    if "sp" not in model_fields and "lnsp" in model_fields:
            print("Computing sp from lnsp")
            model_fields["sp"] = compute_sp_from_lnsp(model_fields["lnsp"])
    populate_mlfields(model_fields, lev_dim=LEV_DIM, half_lev_dim=HLEV_DIM)

    print(f"Variables in model_fields: {', '.join([str(v) for v in model_fields.variables])}")

    return cast_floats(model_fields)

def aerosol_clim(model_data : xr.Dataset,
                 model_pres_var : Optional[str] = None) -> OpticsData:

    from ..aerodefs.optics import get_optical_aerosol
    from .config import get_cams_clim_path, CONFIGDICT

    from .types import NS_DTTYPE

    from ..tools.cams_clim import interpolate_monthly_clim, PRES_VAR

    from ..interp.wrapper import GriddedProfile

    # Data collection
    cams_dset = get_optical_aerosol(
        xr.open_dataset(get_cams_clim_path()),
        CONFIGDICT["aeosol_optics_version"]
        )
    cams_dset.data = _full_rename(cast_floats(cams_dset.data))
    
    # Time interpolation
    if TIME_DIM in model_data.coords:
        model_times = model_data.coords[TIME_DIM]
        # Use unique dates for interpolation
        model_dates = xr.DataArray(
            data=np.unique(model_times.dt.date.astype(NS_DTTYPE)),
            dims=TIME_DIM)

        cams_intp_dset = interpolate_monthly_clim(cams_dset.data, model_dates)
        assert isinstance(cams_intp_dset, xr.Dataset)
        cams_dset.data = cams_intp_dset.sel({TIME_DIM: model_times}, method="nearest").assign_coords({TIME_DIM: model_times})

    if model_pres_var is None:
        model_pres_var = PRES_VAR
    
    if model_pres_var not in model_data.data_vars:
        raise ValueError(f"Model pressure variable '{model_pres_var}' not found in model data variables: {list(model_data.data_vars.keys())}")

    cams_dset.data = xr.Dataset(
        data_vars={
            "aerosol_mmr": xr.concat(
                [cams_dset.data.data_vars[var].expand_dims(aero_type=[v])
                 for v,var in enumerate(cams_dset.optics_var)], dim="aero_type"),
                 PRES_VAR: cams_dset.data[PRES_VAR]
                 },
        coords=cams_dset.data.coords,
        attrs=cams_dset.data.attrs
    )
    
    # 3D interpolation to model grid and vertical levels
    gp_clim = GriddedProfile(cams_dset.data, profile_coord=PRES_VAR, lev_dim=LEV_DIM)
    ptgt = model_data[model_pres_var]
    if USE_DASK:
        ptgt = ptgt.compute(scheduler="threads")
    else:
        ptgt = ptgt.load()
    
    cams_dset.data = gp_clim.interp3d_to(
        ptgt=ptgt,
        tgt_coords=model_data.coords,
        lev_dim_tgt=LEV_DIM,
        verbose=False
    )
    
    return cams_dset

def ghg_data(model_times : Union[None, xr.DataArray] = None) -> xr.Dataset:
    from .config import get_ghg_path
    from .types import NS_TDTYPE, NS_DTTYPE

    def preprocess_ghg_dset(ds):
        """Sets the time dimension"""
        ds = ds.sel(time=slice(1768, 2262))
        years = np.floor(ds["time"].values).astype(int)
        fracs = ds["time"].values - years

        ds = ds.drop_vars("time").rename(time=TIME_DIM)
        one_year = np.timedelta64(1, 'Y').astype(NS_TDTYPE)
        ds = ds.assign_coords({
            TIME_DIM: xr.DataArray(
            data=[np.datetime64(f"{y:4d}-01-01").astype(NS_DTTYPE)+one_year*f
                  for y, f in zip(years, fracs)],
            dims=[TIME_DIM,])})
        return ds

    # Time should be from fcdate
    ghg_data = preprocess_ghg_dset(xr.open_dataset(
        get_ghg_path(),
        decode_times=False))
    
    if model_times is None:
        print("No model times provided, using all GHG dataset times")
        return ghg_data
    
    model_tmin, model_tmax = model_times.min(), model_times.max()
    ghg_tmin, ghg_tmax = ghg_data[TIME_DIM].min(), ghg_data[TIME_DIM].max()
    if model_tmin < ghg_tmin or model_tmax > ghg_tmax:
        raise ValueError(f"Model times are out of bounds of GHG dataset times. " +\
                         f"Model time range: {model_tmin} to {model_tmax}. " +\
                         f"GHG dataset time range: {ghg_tmin} to {ghg_tmax}.")
    
    ghg_data = ghg_data.interp(**{
        TIME_DIM: model_times, "method":"linear"
        }).astype(FLOAT_DTYPE).load()

    return ghg_data