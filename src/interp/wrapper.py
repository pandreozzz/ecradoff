"""Tools to support interpolation"""
from __future__ import annotations

from typing import Optional, List, Union
import numpy as np
import xarray as xr
from xarray.core.coordinates import DatasetCoordinates, DataArrayCoordinates

XRCoords = Union[DatasetCoordinates, DataArrayCoordinates]

from .grids import GridDesc, GridType


class ArrayGridDesc(GridDesc):
    """Grid descriptor with dimensions"""
    def __init__(self, gtyp : GridType,
                 lons : np.ndarray, lats : np.ndarray,
                 datadims : List[str],
                 coordims : List[str],
                 reduced_pts : Optional[np.ndarray] = None
                 ):
        super().__init__(gtyp, lons, lats, reduced_pts)
        self.datadims : List[str] = datadims
        self.coordims : List[str] = coordims

    def is_equal_to(self, other : ArrayGridDesc,
                    atol : float = 1.e-4) -> bool:
        """Check if grid descriptors are equal"""
        def equal_arr(a : np.ndarray, b : np.ndarray,
                      atol : float) -> bool:
            if a.shape != b.shape:
                return False
            return np.allclose(a, b, atol=atol)
        if self.gtyp != other.gtyp:
            return False
        if not equal_arr(self.lon, other.lon, atol=atol):
            return False
        if not equal_arr(self.lat, other.lat, atol=atol):
            return False
        if self.gtyp == GridType.REDUCED:
            assert self.reduced_pts is not None and other.reduced_pts is not None
            if not equal_arr(self.reduced_pts,
                             other.reduced_pts, atol=atol):
                return False
        return True

    def get_lat(self):
        if self.gtyp == GridType.REDUCED:
            assert self.reduced_pts is not None
            lats = np.repeat(self.lat, self.reduced_pts)
        else:
            lats = self.lat
        return lats

    @property
    def xa_lat(self):

        lats = self.get_lat()

        # Dimension name
        if self.gtyp == GridType.REGULAR:
            dims = [d for d in self.datadims
                    if "lat" in d.lower()]
        else:
            dims = self.datadims

        if len(dims) == 0:
            raise ValueError("Could not find latitude dimension in coordims!")
        elif len(dims) > 1:
            raise ValueError("Multiple latitude dimensions found in coordims!")
        else:
            dim = dims[0]

        if self.gtyp == GridType.REGULAR:
            coords = {dim: lats}
        else:
            coords = {}

        # Attributes
        lat_attrs = {
            "long_name": "latitude",
            "units": "degrees_north"
        }
        return xr.DataArray(data=lats, dims=dim,
                            coords=coords, attrs=lat_attrs)
    @property
    def xa_lon(self):

        lons = self.get_lon()

        # Dimension name
        # Dimension name
        if self.gtyp == GridType.REGULAR:
            dims = [d for d in self.datadims
                    if "lon" in d.lower()]
        else:
            dims = self.datadims

        if len(dims) == 0:
            raise ValueError("Could not find latitude dimension in coordims!")
        elif len(dims) > 1:
            raise ValueError("Multiple latitude dimensions found in coordims!")
        else:
            dim = dims[0]

        if self.gtyp == GridType.REGULAR:
            coords = {dim: lons}
        else:
            coords = {}

        # Attributes
        lon_attrs = {
            "long_name": "longitude",
            "units": "degrees_east"
        }
        return xr.DataArray(data=lons, dims=dim, coords=coords,
                            attrs=lon_attrs)


def get_arraygriddesc(grid_coords : XRCoords,
                      datadims : Optional[List[str]] = None,
                      coordims : Optional[List[str]] = None,
                      ) -> ArrayGridDesc:
    """Get grid type from xarray dataset"""
    from .grids import GridType

    if "reduced_points" in grid_coords:
        # Assume lats and reduced_pts use the "lat" dimension
        # Lons are the initial longitude points for each latitude
        lats = grid_coords["lat"].values
        lons = np.zeros_like(lats)
        reduced_pts = grid_coords["reduced_points"].values
        assert len(lats) == len(reduced_pts)
        gtyp = GridType.REDUCED

        datadims = ["col"] if datadims is None else datadims
        coordims = ["lat"] if coordims is None else coordims

    else:
        lats = grid_coords["lat"].values
        lons = grid_coords["lon"].values
        reduced_pts = None

        # if datadims is None, defaults.
        # if grid_coords is nonambiguous, datadims is needed,
        # otherwise assume regular
        if datadims is None:
            if coordims is not None:
                datadims = coordims
            else:
                if "lon" in grid_coords.dims and "lat" in grid_coords.dims:
                    datadims = ["lat", "lon"] #["lon", "lat"]
                elif "col" in grid_coords.dims:
                    datadims = ["col"]
                else:
                    raise ValueError("Could not determine grid type from xarray dataset!")
        else:
            if coordims is not None:
                assert datadims == coordims, "datadims and coordims must be the same for non-reduced grids!"

        if len(datadims) == 2 and "lon" in datadims and "lat" in datadims:
            # Canonical internal order for regular grids is (lat, lon).
            # This keeps flatten/reshape consistent with Fortran (y,x) indexing
            # even if callers provide datadims as (lon, lat).
            lat_dim_name = next(d for d in datadims if d.lower().startswith("lat"))
            lon_dim_name = next(d for d in datadims if d.lower().startswith("lon"))
            datadims = [lat_dim_name, lon_dim_name]
            coordims = datadims
            gtyp = GridType.REGULAR
        # Unstructured dim name is flexible
        elif len(datadims) == 1:
            # Coordims same as datadims for unstructured grids
            coordims = datadims
            gtyp = GridType.UNSTRUC
        else:
            raise ValueError(f"Could not determine grid type from xarray dataset and provided dims {datadims}!")

    return ArrayGridDesc(
        gtyp=gtyp,
        lons=lons, lats=lats,
        datadims=datadims,
        coordims=coordims,
        reduced_pts=reduced_pts
        )


class GriddedArray():
    """Dataset and grids"""
    def __init__(self, array : xr.Dataset,
                 datadims : Optional[List[str]] = None,
                 ):
        """Wrapper for Dataset to use GriddedData tools"""

        self.array : xr.Dataset = array
        self.grid : ArrayGridDesc = get_arraygriddesc(array.coords, datadims=datadims)

    @property
    def ordered_array(self) -> xr.Dataset:
        """Returns data view ensuring that dims are consistent with grid description"""
        if self.grid.gtyp == GridType.REGULAR:
            selcoords = {"lon": self.grid.lon,
                         "lat": self.grid.lat}
        elif self.grid.gtyp == GridType.REDUCED:
            assert self.grid.reduced_pts is not None
            selcoords = {self.grid.coordims[0]: self.grid.lat}
        elif self.grid.gtyp == GridType.UNSTRUC:
            selcoords = {self.grid.coordims[0]: np.arange(self.grid.npts)}
        else:
            raise ValueError(f"Unrecognised grid type {self.grid.gtyp}!")

        return self.array.sel(selcoords)

    def interp2d_to(self, tgt_grid : XRCoords,
                    tgt_datadims : Optional[List[str]] = None,
                    verbose : bool = False,
                    lat_row : bool = True,
                    **interp_kwargs) -> GriddedArray:
        """Interpolate field to target grid and return interpolated GriddedArray.
        interp_kwargs passed directly to interp2d
        lat_row : bool (default True)
            For rectangular grids, what is the shape before flattening if true then (lat, lon) otherwise (lon, lat) - Since the parallelisation in the interp2d library is done along the first dimension, this gives some control for cases with a very different number of lat and lon dimension sizes.

        """
        from ..tools.stack import tools_to_stack_2dgrids
        from .grids import GriddedData

        tgt_grid_desc : ArrayGridDesc = get_arraygriddesc(tgt_grid,
                                                          datadims=tgt_datadims)

        if self.grid.is_equal_to(tgt_grid_desc):
            if verbose:
                print("Source and target grids are the same, skipping interpolation.", flush=True)
            return self

        data_vars = {}
        xtra_coords = {}
        for var in self.array.data_vars:
            if verbose:
                print(f"Interpolating variable {var}...", flush=True)
            # Be aware that if the interpolation coordinates are
            # not in the coordinates vars, but appears in data_vars
            # then this gets broken
            if not np.all([dim in self.array[var].dims
                           for dim in self.grid.datadims]):
                # copy non-gridded variables
                data_vars[var] = self.array[var]
                continue
            this_array = self.ordered_array[var]
            stacktools = tools_to_stack_2dgrids(
                src_arr=this_array,
                srcgrid=self.grid,
                tgtgrid=tgt_grid_desc,
                lat_row=lat_row
                )
            out_xrda_dim_order = [d for d in self.array.dims
                                  if d in stacktools.out_dim_order]

            # Only keep keys that are in the output dimension order
            out_coords = {k:v for k,v in stacktools.out_coords.items() if k in stacktools.out_dim_order}
            # All the rest of the coordinates
            xtra_coords = {**xtra_coords,
                           **{k:v for k,v in stacktools.out_coords.items()
                              if (k not in xtra_coords) and (k not in out_coords)}}

            data_vars[var] = xr.DataArray(
                data=GriddedData(data=this_array.transpose(
                *stacktools.src_dim_order
                ).values.reshape(
                    stacktools.src_stackshape
                    ),
                    grid=self.grid).interp2d_to(
                        tgt_grid=tgt_grid_desc,
                        **{
                            **{"lat_row": lat_row},
                            **interp_kwargs
                        }).data.reshape(stacktools.out_shape),
                dims=stacktools.out_dim_order,
                coords=out_coords
            ).transpose(..., *out_xrda_dim_order)

        return GriddedArray(
            array=xr.Dataset(
                data_vars=data_vars
                ).assign_coords(**xtra_coords),
            datadims=tgt_grid_desc.datadims
        )



class GriddedProfile():
    """Dataset and grids"""
    def __init__(self, profile : Union[xr.Dataset, GriddedArray],
                 profile_coord : str = "pressure",
                 lev_dim : Optional[str] = "lev",
                 ):
        """Wrapper for Dataset to use GriddedData tools
        profile_data must contain profile_coord as a data variable
        """

        if isinstance(profile, GriddedArray):
            self.profile : GriddedArray = profile
        else:
            self.profile : GriddedArray = GriddedArray(profile)

        if lev_dim is not None:
            if lev_dim not in self.profile.array.dims:
                raise ValueError(f"Specified lev_dim {lev_dim} not found in xarray dataset dimensions!")
        else:
            # Try to infer lev_dim
            for candidate in ["lev", "plev", "level"]:
                if candidate in self.profile.array.dims:
                    lev_dim = candidate
                    break
            if lev_dim is None:
                raise ValueError("Could not determine vertical level dimension from xarray dataset!")
        self.lev_dim = lev_dim

        if profile_coord not in self.profile.array.data_vars:
            raise ValueError(f"Specified profile_coord {profile_coord} not found in xarray dataset data variables!")
        self.profile_coord : str = profile_coord

        if self.lev_dim not in self.profile.array[self.profile_coord].dims:
            raise ValueError(f"Specified lev_dim {self.lev_dim} not found in profile_coord {self.profile_coord} dimensions!")

    def interp3d_to(self, ptgt : xr.DataArray,
                    tgt_coords : Optional[XRCoords] = None,
                    tgt_datadims : Optional[List[str]] = None,
                    lev_dim_tgt : Optional[str] = None,
                    verbose : bool = False,
                    out_chunks : Optional[dict] = None,
                    **interp_kwargs) -> GriddedProfile:
        """Interpolate profile to target grid and to target pressure"""

        # Check if ptgt and tgt_coords are consistent and check if tgt_coords has grid info. If not, just do vertical interpolation.
        if tgt_coords is not None:
            try:
                _ = get_arraygriddesc(tgt_coords, datadims=tgt_datadims)
                tgt_has_grid_info = True
            except:
                tgt_has_grid_info = False
        else:
            tgt_has_grid_info = False

        if tgt_has_grid_info:
            if verbose:
                print("Interpolating horizontally to target grid...", flush=True)
            assert tgt_coords is not None
            profile_interp2d = self.interp2d_to(
                tgt_grid=tgt_coords,
                tgt_datadims=tgt_datadims,
                verbose=verbose,
                **interp_kwargs
            )
        else:
            profile_interp2d = self
        if verbose:
            print("Interpolating vertically to target pressure levels...", flush=True)
        profile_interp3d = profile_interp2d.interpvert_to(
            ptgt=ptgt,
            lev_dim_tgt=lev_dim_tgt,
            vert_kwargs=interp_kwargs,
            interp_kwargs=interp_kwargs,
            out_chunks=out_chunks,
        )
        return GriddedProfile(
            profile=profile_interp3d.profile,
            profile_coord=self.profile_coord,
            lev_dim=profile_interp3d.lev_dim
        )

    def interp2d_to(self, tgt_grid : XRCoords,
                    tgt_datadims : Optional[List[str]] = None,
                    verbose : bool = False,
                    **interp_kwargs) -> GriddedProfile:
        """Interpolate profile to target grid and return interpolated GriddedProfile.
        interp_kwargs passed directly to interp2d
        """
        return GriddedProfile(
            profile=self.profile.interp2d_to(
                tgt_grid=tgt_grid,
                tgt_datadims=tgt_datadims,
                verbose=verbose,
                **interp_kwargs
            ),
            profile_coord=self.profile_coord,
            lev_dim=self.lev_dim
        )

    def interpvert_to(self, ptgt : xr.DataArray,
                      lev_dim_tgt : Optional[str] = None,
                      vert_kwargs : dict = {},
                      interp_kwargs : dict = {},
                      out_chunks : Optional[dict] = None
                      ) -> GriddedProfile:
        """Interpolate to target pressure levels."""
        from .verticals import ProfileData, InterpWeights
        from ..tools.stack import tools_to_stack_xarrays

        def get_tmp_dim_name(base_name : str, existing_dims : List[str]) -> str:
            """Get a temporary dimension name that doesn't conflict with existing dimensions"""
            if base_name not in existing_dims:
                return base_name
            else:
                i = 1
                while f"{base_name}_{i}" in existing_dims:
                    i += 1
                return f"{base_name}_{i}"

        if lev_dim_tgt is None:
            lev_dim_tgt = self.lev_dim
            comm_lev_dim = self.lev_dim
        else:
            if lev_dim_tgt != self.lev_dim:
                comm_lev_dim = get_tmp_dim_name(
                    "lev",
                    [str(d) for d in self.profile.array.dims] +\
                            [str(d) for d in ptgt.dims])
            else:
                comm_lev_dim = self.lev_dim

        if lev_dim_tgt not in ptgt.dims:
            raise ValueError(f"Specified lev_dim_tgt {lev_dim_tgt} not found in target xarray dataset dimensions!")

        rename_dim_src = {self.lev_dim: comm_lev_dim} if comm_lev_dim != self.lev_dim else {}
        rename_dim_tgt = {lev_dim_tgt: comm_lev_dim} if comm_lev_dim != lev_dim_tgt else {}

        xa_psrc = self.profile.array[self.profile_coord].rename(rename_dim_src)
        xa_ptgt = ptgt.rename(rename_dim_tgt)

        tmp_stacktools = tools_to_stack_xarrays(
            src_arr=xa_psrc,
            dst_arr=xa_ptgt,
            intp_dim_name=comm_lev_dim
        )
        out_xrda_dim_order = [d for d in xa_ptgt.dims
                              if d in tmp_stacktools.out_dim_order]

        pd_psrc = ProfileData(
            data = None,
            pres = xa_psrc.transpose(*tmp_stacktools.src_dim_order).values.reshape(tmp_stacktools.src_stackshape),
            nlevs = len(xa_psrc[comm_lev_dim])
        )
        pd_ptgt = ProfileData(
            data = None,
            pres = xa_ptgt.transpose(*tmp_stacktools.dst_dim_order).values.reshape(tmp_stacktools.dst_stackshape),
            nlevs = len(xa_ptgt[comm_lev_dim])
        )

        interp_weights = pd_psrc.get_weights_to(
            pd_ptgt,
            **vert_kwargs
        )
        # Generate DataArrays
        # Variable might have extra dimensions: in this case, weights are broadcast.
        xa_tgtidxs = xr.DataArray(interp_weights.tgtidxs.reshape(tmp_stacktools.out_shape), dims=tmp_stacktools.out_dim_order, coords=tmp_stacktools.out_coords)
        xa_weights = xr.DataArray(interp_weights.weights.reshape(tmp_stacktools.out_shape), dims=tmp_stacktools.out_dim_order, coords=tmp_stacktools.out_coords)

        del tmp_stacktools, interp_weights

        data_vars = {}
        xtra_coords = {}

        if self.profile.grid.gtyp == GridType.REDUCED:
            for nondim_coord in ["lat", "reduced_points"]:
                if nondim_coord in self.profile.array.coords:
                    xtra_coords[nondim_coord] = self.profile.array.coords[nondim_coord]
                else:
                    xtra_coords[nondim_coord] = xr.DataArray(data=(
                        "lat",
                        getattr(self.profile.grid, nondim_coord))
                        )

        for var in self.profile.array.data_vars:
            # No need to interpolate the profile coordinate!
            if var == self.profile_coord:
                continue
            # Variable is not a profile
            if self.lev_dim not in self.profile.array[var].dims:
                data_vars[var] = self.profile.array[var]
                continue

            xa_fld_src = self.profile.array[var].rename(rename_dim_src)

            tmp_stacktools = tools_to_stack_xarrays(
                src_arr=xa_fld_src,
                dst_arr=xa_tgtidxs,
                intp_dim_name=comm_lev_dim
            )
            out_coords = {k:v for k,v in tmp_stacktools.out_coords.items() if k in tmp_stacktools.out_dim_order}
            xtra_coords = {**xtra_coords,
                            **{k:v for k,v in tmp_stacktools.out_coords.items()
                                if (k not in xtra_coords) and (k not in out_coords)}}

            this_interp_weights = InterpWeights(
                tgtidxs=xa_tgtidxs.transpose(*tmp_stacktools.dst_dim_order).values.reshape(tmp_stacktools.dst_stackshape),
                weights=xa_weights.transpose(*tmp_stacktools.dst_dim_order).values.reshape(tmp_stacktools.dst_stackshape)
            )

            pd_fld_out = ProfileData(
                data = xa_fld_src.transpose(*tmp_stacktools.src_dim_order).values.reshape(tmp_stacktools.src_stackshape),
                pres = None,
                nlevs = pd_psrc.nlevs
            ).interp_fld(this_interp_weights, **interp_kwargs)

            data_vars[var] = xr.DataArray(
                data=pd_fld_out.reshape(tmp_stacktools.out_shape),
                dims=tmp_stacktools.out_dim_order,
                coords=out_coords
            ).transpose(..., *out_xrda_dim_order)
            if out_chunks is not None:
                data_vars[var] = data_vars[var].chunk(out_chunks)
            del tmp_stacktools, this_interp_weights, pd_fld_out

        xds_out = xr.Dataset(
            data_vars=data_vars
            ).assign_coords(**xtra_coords)
        if rename_dim_tgt != {}:
            xds_out = xds_out.rename({v: k for k,v in rename_dim_tgt.items()})

        xds_out[self.profile_coord] = xa_ptgt.rename({comm_lev_dim: lev_dim_tgt})

        return GriddedProfile(
            profile=xds_out,
            profile_coord=self.profile_coord,
            lev_dim=lev_dim_tgt
        )

