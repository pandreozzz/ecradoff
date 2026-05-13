def interpolate_monthly_aerosols(dset: xr.Dataset or xr.DataArray,
                                 dates: xr.DataArray or np.ndarray) -> xr.Dataset or xr.DataArray:
    """
    Interpolation of aerosol monthly climatologies
    """

    prev_month     = (dates.values - np.timedelta64(14,'D')).astype("datetime64[M]")+ np.timedelta64(14,'D')
    foll_month     = (prev_month + np.timedelta64(18,'D')).astype("datetime64[M]") + np.timedelta64(14,'D')

    monthdelta     = foll_month - prev_month
    thisdelta      = dates  - prev_month.astype("datetime64[ns]")

    timeweight     = thisdelta/monthdelta # in [0,1]

    intmonths_bot  = xr.DataArray(data=prev_month.astype("datetime64[ns]"), coords={"time":dates}).dt.month
    intmonths_top  = xr.DataArray(data=foll_month.astype("datetime64[ns]"), coords={"time":dates}).dt.month

    return  (1-timeweight) * dset.sel(month=intmonths_bot).drop_vars("month") +\
        timeweight*dset.sel(month=intmonths_top).drop_vars("month")

def complete_lon_periodic(dset: xr.Dataset or xr.DataArray, method="linear") -> xr.Dataset or xr.DataArray:
    """
    Completes eventually missing 0. and 360. longitude values on periodic domain
    only supported method is linear. If either 0 or 360 are present,
    then those are copied to fill the missing values.
    """
    this_dset = dset.copy()

    minlon = this_dset.lon.min().values
    maxlon = this_dset.lon.max().values

    appendmax = maxlon < 360
    appendmin = minlon > 0


    # Missing 360
    if appendmin or appendmax:
        wei = minlon/(360-maxlon+minlon)
        borderslice = wei*this_dset.sel(lon=maxlon, drop=True) + (1-wei)*this_dset.sel(lon=minlon, drop=True)

    # Missing 0
    if appendmax:
        maxslice = borderslice.assign_coords(lon=360)
        this_dset = xr.concat(
            [this_dset, maxslice], dim="lon")

    if appendmin:
        minslice = borderslice.assign_coords(lon=0)
        this_dset = xr.concat(
            [minslice, this_dset], dim="lon")

    this_dset = this_dset.sortby("lon")

    return this_dset

def complete_lat_boundaries(dset: xr.Dataset or xr.DataArray) -> xr.Dataset or xr.DataArray:
    """
    Completes latitudes 90N and 90S boundary values
    """

    this_dset = dset.copy()

    minlat = this_dset.lat.min().values
    maxlat = this_dset.lat.max().values

    appendmax = maxlat < 90
    appendmin = minlat > -90

    if appendmin:
        minslice = this_dset.sel(lat=minlat).assign_coords(lat=-90)
        this_dset = xr.concat(
            [minslice, this_dset], dim="lat"
        )
    if appendmax:
        maxslice = this_dset.sel(lat=maxlat).assign_coords(lat=90)
        this_dset = xr.concat(
            [this_dset, maxslice], dim="lat"
        )

    this_dset = this_dset.sortby("lat", ascending=False)

    return this_dset


def interpolate_3d_aerosols(aerosol_fields : xr.DataArray, model_pres : xr.DataArray,
global_domain : bool = True, reduced_pts_src = None, reduced_pts_tgt = None, flat_dim : str = "col"):
    import interp_2d_iface as i2d
    import fvertintp_iface as fvint
    import stack_tools as stack

    if reduced_pts_src is not None:
        lats_src = reduced_pts_src.lat.values
        lons_src = np.zeros_like(lats_src)
        redpts_src = reduced_pts_src.values
        gtyp_src = 2
    else:
        lats_src = aerosol_fields.lat.values
        lons_src = aerosol_fields.lon.values
        redpts_src = None
        gtyp_src = 3 if "col" in aerosol_fields.dims else 1

    if reduced_pts_tgt is not None:
        lats_tgt = reduced_pts_tgt.lat.values
        lons_tgt = np.zeros_like(lats_tgt)
        redpts_tgt = reduced_pts_tgt.values
        gtyp_tgt = 2
    else:
        lats_tgt = model_pres.lat.values
        lons_tgt = model_pres.lon.values
        redpts_tgt = None
        gtyp_tgt = 3 if "col" in model_pres.dims else 1
    
        
    src2dgrid = i2d.GridDef("srcgrid", gtyp=gtyp_src, lons=lons_src, lats=lats_src, reduced_pts=redpts_src)
    tgt2dgrid = i2d.GridDef("tgtgrid", gtyp=gtyp_tgt, lons=lons_tgt, lats=lats_tgt, reduced_pts=redpts_tgt)

    aerosol_hintp = xr.Dataset()
    xtra_coords = {}
    for var in ["pressure", "aerosol_mmr"]:
        tmp_stacktools = stack.tools_to_stack_2dgrids(aerosol_fields[var], srcgrid=src2dgrid,
                                                      tgtgrid=tgt2dgrid, flat_dim=flat_dim)
        out_coords = {k:v for k,v in tmp_stacktools.out_coords.items() if k in tmp_stacktools.out_dim_order}
        xtra_coords = {**xtra_coords,
                       **{k:v for k,v in tmp_stacktools.out_coords.items()
                          if (k not in xtra_coords) and (k not in out_coords)}}
        aerosol_hintp[var] = xr.DataArray(
            data=i2d.interp_2d(
                aerosol_fields[var].transpose(*tmp_stacktools.src_dim_order).values.reshape(tmp_stacktools.src_stackshape),
                srcgrid=src2dgrid, tgtgrid=tgt2dgrid
            ).reshape(tmp_stacktools.out_shape),
            dims=tmp_stacktools.out_dim_order,
            coords=out_coords,
        )
        del tmp_stacktools, out_coords
    for k,v in xtra_coords.items():
        aerosol_hintp[k] = v
    del xtra_coords

    # Reorder aerosol fields
    tmp_src = aerosol_hintp[PDIM]
    tmp_dst = model_pres

    tmp_stacktools = stack.tools_to_stack_xarrays(
        src_arr=tmp_src, dst_arr=tmp_dst,
        intp_dim_name="lev")

    tgtlevs, weights = fvint.interp(
        psrc=tmp_src.transpose(*tmp_stacktools.src_dim_order).values.reshape(tmp_stacktools.src_stackshape),
        ptgt=tmp_dst.transpose(*tmp_stacktools.dst_dim_order).values.reshape(tmp_stacktools.dst_stackshape)
    )
    tgtlevs = xr.DataArray(data=tgtlevs.reshape(tmp_stacktools.out_shape),
                            dims=tmp_stacktools.out_dim_order, coords=tmp_stacktools.out_coords)
    weights = xr.DataArray(data=weights.reshape(tmp_stacktools.out_shape),
                            dims=tmp_stacktools.out_dim_order, coords=tmp_stacktools.out_coords)

    del tmp_stacktools, tmp_src, tmp_dst

    aero_mmr = aerosol_hintp["aerosol_mmr"]

    tmp_stacktools = stack.tools_to_stack_xarrays(
        src_arr=aero_mmr, dst_arr=tgtlevs,
        intp_dim_name="lev")

    print("Vertically interpolating aerosol fields...")
    aero_mmr_vintp = xr.DataArray(
        data=fvint.interp_fld(
            fsrc=aero_mmr.transpose(*tmp_stacktools.src_dim_order).values.reshape(tmp_stacktools.src_stackshape),
            tgtlevs=tgtlevs.transpose(*tmp_stacktools.dst_dim_order).values.reshape(tmp_stacktools.dst_stackshape),
            weights=weights.transpose(*tmp_stacktools.dst_dim_order).values.reshape(tmp_stacktools.dst_stackshape)
            ).reshape(tmp_stacktools.out_shape),
        dims=tmp_stacktools.out_dim_order, coords=tmp_stacktools.out_coords
    )

    print("Interpolation done")
    return aero_mmr_vintp