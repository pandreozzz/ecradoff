"""Utilities to stack/unstack xarray fields and LUT structures for fast NumPy access."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union, NamedTuple, Literal

import numpy as np
import xarray as xr

from ..interp.wrapper import ArrayGridDesc


class StackTools(NamedTuple):
    """Helper information to reorder and reshape arrays for interpolation."""
    src_dim_order: List[str]
    dst_dim_order: List[str]
    out_dim_order: List[str]
    src_stackshape: Tuple[int, ...]
    dst_stackshape: Tuple[int, ...]
    out_coords: Dict[str, Any]
    out_shape: Tuple[int, ...]

def tools_to_stack_2dgrids(src_arr : xr.DataArray,
                           srcgrid : ArrayGridDesc,
                           tgtgrid : ArrayGridDesc,
                           lat_row : bool = True
                           ):
    """
    lat_row : bool (default True)
        whether latitude is major dimension for flattening (lat, lon) grid
        otherwise longitude major (lon, lat)
    """
    from ..interp.grids import GridType 
    import xarray as xr

    def get_ordered_intp_dims(grid : ArrayGridDesc, lat_row : bool,
                              coord_dims : bool = False,) -> List[str]:
        these_dims = grid.datadims if not coord_dims else grid.coordims

        if grid.gtyp == GridType.REGULAR:
            assert len(these_dims) == 2, "2 data dimensions needed for regular grids!"
            lat_dim_name = None
            lon_dim_name = None
            for dn in these_dims:
                if dn.lower().startswith("lat"):
                    lat_dim_name = dn
                if dn.lower().startswith("lon"):
                    lon_dim_name = dn
            assert lat_dim_name is not None, "could not find latitude dim but grid is regular!"
            assert lon_dim_name is not None, "could not find longitude dim but grid is regular!"
            if lat_row:
                return [lat_dim_name, lon_dim_name]
            else:
                return [lon_dim_name, lat_dim_name]
        else:
            assert len(these_dims) == 1, "1 data dimension needed for non-regular grids!"
            return these_dims
    def get_ordered_datadims(grid : ArrayGridDesc, lat_row : bool) -> List[str]:
        return get_ordered_intp_dims(grid, lat_row, coord_dims=False)
    def get_ordered_coorddims(grid : ArrayGridDesc, lat_row : bool) -> List[str]:
        return get_ordered_intp_dims(grid, lat_row, coord_dims=True)
    
    # Interpolation dimension of data fields
    intp_dims_src = get_ordered_datadims(srcgrid, lat_row)
    # The rest comes after
    nonintp_dims_src : List[str] = [str(d) for d in src_arr.dims if d not in intp_dims_src]
    nonintp_ndims_src = [len(src_arr[src_dim]) for src_dim in nonintp_dims_src]

    intp_ndims_src = [len(src_arr[src_dim]) for src_dim in intp_dims_src]

    src_dim_order = nonintp_dims_src+intp_dims_src
    
    src_stackshape = [int(max(np.prod(nonintp_ndims_src),1))]+\
    [int(np.prod(intp_ndims_src))]

    latattrs = {"units" : "degrees_north", "long_name":"Latitude",
    "standard_name":"latitude", "axis":"X"}
    lonattrs = {"units" : "degrees_east", "long_name":"Longitude",
    "standard_name":"longitude", "axis":"Y"}

    # Interpolation dimensions for data in target
    intp_dims_tgt = get_ordered_datadims(tgtgrid, lat_row)
    # Dimensions for coords in target (for reduced, they are different than for data!)
    # lat_row=true so that they are always (lat, lon) for regular
    coordims_tgt = get_ordered_coorddims(tgtgrid, lat_row=True)
    if (tgtgrid.gtyp == GridType.REGULAR):
        assert len(intp_dims_tgt) == 2
        assert len(coordims_tgt) == 2

        intp_coords_tgt = {
            "lat" : xr.DataArray(data=tgtgrid.lat, dims=coordims_tgt[0], attrs=latattrs),
            "lon" : xr.DataArray(data=tgtgrid.lon, dims=coordims_tgt[1], attrs=lonattrs),
        }
        intp_ndims_tgt = [len(getattr(tgtgrid, dim)) for dim in intp_dims_tgt]
        
    else:
        assert len(coordims_tgt) == 1
        assert len(intp_dims_tgt) == 1
        intp_ndims_tgt = [tgtgrid.npts]
        # Reduced
        if (tgtgrid.gtyp == GridType.REDUCED):
            assert tgtgrid.reduced_pts is not None
            intp_coords_tgt = {
            "lat" : xr.DataArray(data=tgtgrid.lat, dims=coordims_tgt, attrs=latattrs),
            "reduced_points" : xr.DataArray(data=tgtgrid.reduced_pts, dims=coordims_tgt),
        }
        # Unstructured 
        elif (tgtgrid.gtyp == GridType.UNSTRUC):
            intp_coords_tgt = {
            "lat" : xr.DataArray(data=tgtgrid.lat, dims=coordims_tgt, attrs=latattrs),
            "lon" : xr.DataArray(data=tgtgrid.lon, dims=coordims_tgt, attrs=lonattrs),
        }
        else:
            raise ValueError(f"Unrecognised target grid type {tgtgrid.gtyp}!")

    out_dim_order = nonintp_dims_src+intp_dims_tgt
    out_shape = tuple(nonintp_ndims_src+intp_ndims_tgt)
    out_coords = {
        **{dim: src_arr.coords[dim] for dim in nonintp_dims_src},
        **intp_coords_tgt
    }
    
    return StackTools(
        src_dim_order = src_dim_order,
        dst_dim_order = [],
        out_dim_order = out_dim_order,
        src_stackshape = tuple(src_stackshape),
        dst_stackshape = tuple(),
        out_coords = out_coords,
        out_shape = out_shape
    )

def tools_to_stack_xarrays(
    src_arr: Union[xr.Dataset, xr.DataArray],
    dst_arr: Union[xr.Dataset, xr.DataArray],
    intp_dim_name: Optional[str] = None,
) -> StackTools:
    """
    Compute all tools to reorder/reshape arrays with NumPy for efficient interpolation.

    Parameters
    ----------
    src_arr : xr.Dataset
        Source array(s).
    dst_arr : xr.Dataset
        Destination array(s).
    intp_dim_name : str
        Interpolation dimension name (e.g., vertical level). If None, ignored.

    Returns
    -------
    StackTools
        All metadata and shapes needed to flatten/reshape consistently.
    """
    nonintp_dims_src = [str(d) for d in src_arr.dims if d != intp_dim_name]
    nonintp_dims_dst = [str(d) for d in dst_arr.dims if d != intp_dim_name]

    common_dim_names = list(set(nonintp_dims_src) & set(nonintp_dims_dst))
    unique_dim_names = list(set(nonintp_dims_src) ^ set(nonintp_dims_dst))
    onlysrc_dim_names = [d for d in nonintp_dims_src if d in unique_dim_names]
    onlydst_dim_names = [d for d in nonintp_dims_dst if d in unique_dim_names]

    src_dim_order = common_dim_names + onlysrc_dim_names
    dst_dim_order = common_dim_names + onlydst_dim_names
    if intp_dim_name:
        src_dim_order = src_dim_order + [intp_dim_name]
        dst_dim_order = dst_dim_order + [intp_dim_name]

    com_ndims = [len(src_arr[com_dim]) for com_dim in common_dim_names]
    for com_dim in common_dim_names:
        if len(src_arr[com_dim]) != len(dst_arr[com_dim]):
            raise ValueError(
                f"Non-interpolation common dimension '{com_dim}' mismatch length in source ({len(src_arr[com_dim])}) and destination ({len(dst_arr[com_dim])})"
            )
    src_ndims = [len(src_arr[src_dim]) for src_dim in onlysrc_dim_names]
    src_intp_ndims = [len(src_arr[intp_dim_name])] if intp_dim_name in src_arr.dims else []
    dst_ndims = [len(dst_arr[dst_dim]) for dst_dim in onlydst_dim_names]
    dst_intp_ndims = [len(dst_arr[intp_dim_name])] if intp_dim_name in dst_arr.dims else []

    def _prod_or_one(ndims: List[int]) -> int:
        return int(np.array(ndims).prod()) if ndims != [] else 1

    src_stackshape = tuple(_prod_or_one(ndims) for ndims in [com_ndims, src_ndims, src_intp_ndims])
    dst_stackshape = tuple(_prod_or_one(ndims) for ndims in [com_ndims, dst_ndims, dst_intp_ndims])

    out_shape: Optional[Tuple[int, ...]] = tuple(com_ndims + src_ndims + dst_ndims + dst_intp_ndims)
    out_dim_order: Optional[List[str]]
    if len(out_shape) == 0:
        out_shape = None
        out_dim_order = None
    else:
        out_dim_order = common_dim_names + onlysrc_dim_names + onlydst_dim_names
        if intp_dim_name:
            out_dim_order = out_dim_order + [intp_dim_name]

    out_coords: Dict[str, Any] = {
        **{str(com_dim): src_arr.coords[com_dim] for com_dim in common_dim_names},
        **{str(src_dim): src_arr.coords[src_dim] for src_dim in onlysrc_dim_names},
        **{str(dst_dim): dst_arr.coords[dst_dim] for dst_dim in onlydst_dim_names},
    }
    if intp_dim_name:
        out_coords[intp_dim_name] = dst_arr.coords[intp_dim_name]

    arglist = [src_dim_order, dst_dim_order, out_dim_order] + [
        src_stackshape,
        dst_stackshape,
        out_coords,
        out_shape,
    ]
    return StackTools(*arglist)


class StackAeroTuple(NamedTuple):
    """Flattened aerosol dataset and related metadata."""
    data: np.ndarray
    w_data: Optional[np.ndarray]
    varlist: List[str]
    w_varlist: List[str]
    dimorder: Tuple[str, ...]
    fldshape: Tuple[int, ...]
    coords: Dict[str, Any]
    c_order: bool


def get_stacked_aero(
    fld: xr.Dataset,
    this_lutaero: List[Any],
    include_w: Optional[List[str]] = None,
    include_rainratio: bool = False,
    flat_order: Union[Literal['C'], Literal['F']] = 'C',
) -> StackAeroTuple:
    """
    Stack aerosol fields to a 2D array (n_points, n_vars) with deterministic dim order.

    Parameters
    ----------
    fld : xr.Dataset
        Dataset with aerosol variables and coordinates.
    this_lutaero : list
        List of LUT aerosol specs (objects with `.name` attribute).
    include_w : list[str], optional
        Extra variables to include (e.g., ["w_mean", "w_prime"]).
    include_rainratio : bool
        If True, also include "rainratio".
    flat_order : {"C", "F"}
        Flattening order.

    Returns
    -------
    StackAeroTuple
        Flattened arrays plus metadata for reshaping.
    """
    include_w = include_w or []

    varlist = [v.name for v in this_lutaero if v.name in list(fld.data_vars)]
    if not varlist:
        raise ValueError(
            "No LUT variables found in field. "
            f"LUT specs: {this_lutaero}; available vars: {list(fld.data_vars)}"
        )

    if include_rainratio:
        varlist = varlist + ["rainratio"]

    dimorder = tuple(fld[varlist[0]].dims)
    fldshape = tuple(len(fld[dim]) for dim in dimorder)
    print(f"Enforcing dimension order on aero fields: {dimorder} {fldshape}")

    # Stack main variables
    fld_ordered = np.ascontiguousarray(
        np.concatenate(
            [
                fld[var].transpose(*dimorder).values.flatten(order=flat_order)[:, None]
                for var in varlist
            ],
            axis=-1,
        )
    )

    # Stack w-variables (optional)
    w_fld_ordered = (
        np.ascontiguousarray(
            np.concatenate(
                [
                    fld[var].transpose(*dimorder).values.flatten(order=flat_order)[:, None]
                    for var in include_w
                ],
                axis=-1,
            )
        )
        if include_w
        else None
    )

    coords = {d: fld.coords[d] for d in dimorder if d in fld.coords}
    return StackAeroTuple(
        fld_ordered,
        w_fld_ordered,
        varlist,
        include_w,
        dimorder,
        fldshape,
        coords,
        flat_order == "C",
    )


class StackLutTuple(NamedTuple):
    """Flattened LUT maps and bins with metadata."""
    lut_maps: np.ndarray
    maplist: List[str]
    maptypes: List[str]
    lut_bins: xr.Dataset
    c_order: bool


def get_stacked_lut(
    pyrcel_lut: xr.Dataset,
    wspeed_type: int,
    lutkin: bool,
    include_mass_act_lut: bool = False,
    flat_order: Union[Literal['C'], Literal['F']] = 'C',
) -> StackLutTuple:
    """
    Stack LUT maps (num_act[/mass_act] per species) into a 2D array and collect bin axes.

    Parameters
    ----------
    pyrcel_lut : xr.Dataset
        LUT dataset with maps and bin coordinates.
    wspeed_type : int
        Vertical velocity parameterization (0..4).
    lutkin : bool
        If True, use kinetically limited variants.
    include_mass_act_lut : bool
        If True, also stack mass activation maps.
    flat_order : {"C", "F"}
        Flattening order.

    Returns
    -------
    StackLutTuple
        Flattened maps, their names, and the associated bin coordinates.
    """
    num_act_var = "num_act_kn" if lutkin else "num_act"
    mass_act_var = "mass_act_kn" if lutkin else "mass_act"

    maptypes = [num_act_var]
    if include_mass_act_lut:
        maptypes.append(mass_act_var)

    maplist = [
        f"aero{spn}_{maptype}" for maptype in maptypes for spn in pyrcel_lut.spec_num.values
    ]

    lut_maps = np.ascontiguousarray(
        np.concatenate([pyrcel_lut[mapname].values.flatten(order=flat_order)[:, None]
                        for mapname in maplist],
                       axis=-1)
    )

    lut_bins_vars: List[str] = [f"aero{spn}_nccn" for spn in pyrcel_lut.spec_num.values if f"aero{spn}" in pyrcel_lut.dims]
    if wspeed_type in [1, 2]:
        lut_bins_vars = ["w_mean"] + lut_bins_vars
    elif wspeed_type in [3, 4]:
        lut_bins_vars = ["w_mean", "w_prime"] + lut_bins_vars

    return StackLutTuple(lut_maps, maplist, maptypes, pyrcel_lut[lut_bins_vars], flat_order == "C")
