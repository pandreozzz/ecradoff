"""Handle CAMS climatology structure"""

from typing import Union

import xarray as xr
import numpy as np


PRES_VAR = "pressure"
HLDELPRES_VAR = "half_level_delta_pressure"
HLPRES_VAR = "half_level_pressure"

# -----------------------------------------------------------------------------
# Monthly climatology interpolation
# -----------------------------------------------------------------------------
def interpolate_monthly_clim(
    dset: Union[xr.Dataset, xr.DataArray],
    dates: xr.DataArray,
) -> Union[xr.Dataset, xr.DataArray]:
    """Linearly interpolate a monthly climatology to arbitrary `dates`.

    The input `dset` **must** have a `month` coordinate/dimension.
    dates must have 'time' coordinate and dimension, but can also have other coordinates
    along the time dimension. If those coordinates are also a dimension in dset, then these are considered
    aligned to the time interpolation and the result will be contracted along those dimensions as well.

    For example, if dset has dimensions (month, lat) and dates has coordinates (time(time), lat(time)), 
    then the output will be the interpolation along time of the corresponding lat,
    so will have only the time dimension.

    Parameters
    ----------
    dset : xr.Dataset or xr.DataArray
        Monthly climatology with a `month` dimension.
    dates : xr.DataArray
        Target dates (datetime64).

    Returns
    -------
    xr.Dataset or xr.DataArray
        The same kind as input, interpolated to `dates` along a new/used `time` coordinate.
    """
    if len(dates.time) == 0:
        collapse_output = True
        dates = dates.isel(time=[0])
    else:
        collapse_output = False

    # Previous and following "anchor" mid-months around the target dates
    prev_month = (dates.values - np.timedelta64(14, "D")).astype("datetime64[M]") +\
          np.timedelta64(14, "D")
    foll_month = (prev_month + np.timedelta64(18, "D")).astype("datetime64[M]") +\
          np.timedelta64(14, "D")

    monthdelta = foll_month - prev_month
    thisdelta = dates - prev_month.astype("datetime64[ns]")
    timeweight_m = thisdelta / monthdelta  # in [0, 1]

    months_coords = {
        **{"time": dates},
        **{cname: coord
           for cname, coord in dates.coords.items()
           if cname != "month"},
    }
    
    intmonths_bot = xr.DataArray(
        data=prev_month.astype("datetime64[ns]"),
        dims=["time"],
        coords=months_coords
        ).dt.month
    intmonths_top = xr.DataArray(
        data=foll_month.astype("datetime64[ns]"),
        dims=["time"],
        coords=months_coords
        ).dt.month

    add_sel_kwargs = {str(cname): coord
                      for cname, coord in dates.coords.items()
                      if cname in dset.dims}

    # Take 
    if "epoch" in dset.coords:
        epochs = dset.epoch.values
        max_epoch, min_epoch = epochs.max(), epochs.min()
        target_years = dates.dt.year.values
        lower_idx  = np.clip(np.searchsorted(epochs, target_years, side="right") - 1, 0, len(epochs)-1)
        upper_idx  = np.clip(np.searchsorted(epochs, target_years, side="left"), 0, len(epochs)-1)
        epoch_lower = epochs[lower_idx]
        epoch_upper = epochs[upper_idx]

        span = epoch_upper - epoch_lower
        timeweight_e = np.divide(
            target_years - epoch_lower,
            span,
            out=np.zeros_like(span, dtype=float),
            where=~np.isclose(span, 0.0, atol=1.e-2)
            )
        needed_epochs = np.unique(np.concat([epoch_lower, epoch_upper]))
        add_sel_kwargs["epoch"] = needed_epochs
    else:
        timeweight_e = []
        epoch_lower = []
        epoch_upper = []
    
    lower = dset.sel(month=intmonths_bot, **add_sel_kwargs).drop_vars("month")
    upper = dset.sel(month=intmonths_top, **add_sel_kwargs).drop_vars("month")

    dset_intp_m = (1 - timeweight_m) * lower + timeweight_m * upper
    if "epoch" in dset_intp_m.coords:
        timeweight_e = xr.DataArray(data=timeweight_e, dims=["time"])
        epoch_lower = xr.DataArray(data=epoch_lower, dims=["time"])
        epoch_upper = xr.DataArray(data=epoch_upper, dims=["time"])

        dset_intp = (1-timeweight_e) * dset_intp_m.sel(epoch=epoch_lower) +\
            timeweight_e * dset_intp_m.sel(epoch=epoch_upper)
    else:
        dset_intp = dset_intp_m

    if collapse_output:
        dset_intp = dset_intp.isel(time=0, drop=False)

    return dset_intp