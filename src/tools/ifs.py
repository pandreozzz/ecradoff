"""Numpy version of IFS subroutines"""

import logging
from typing import Optional, Union

import numpy as np
import xarray as xr

from ..main.types import FLOAT_DTYPE

# Magic constants
RNAVO = 6.0221367E+23
RKBOL = 1.380658E-23
RMD   = 28.9644
R     = RNAVO*RKBOL
RD    = 1000*R/RMD

GCNST = 9.81 # m/s2
CPAIR = 1004.5 # J/Kg/K

def pick_level(da : xr.DataArray, levsel : xr.DataArray, lev_dim : str = "lev") -> xr.DataArray:
    """Pick values from da at levsel. Uses numpy instead of isel
    """

    # Expected dims after removing lev
    target_dims = tuple(d for d in da.dims if d != lev_dim)

    # Check that levsel dims are subset of target dims
    missing_dims = set(levsel.dims) - set(target_dims)
    if missing_dims:
        raise ValueError(
            f"levsel has extra dims {missing_dims} not in da without '{lev_dim}'. "
            f"da dims: {da.dims}, levsel dims: {levsel.dims}"
        )

    da_data = da.data
    da_out = da.isel({lev_dim: 0}, drop=True)

    # Broadcast levsel to match da_out
    try:
        levsel_broad = levsel.broadcast_like(da_out)
    except Exception as e:
        raise ValueError(
            f"Cannot broadcast levsel to da without '{lev_dim}'. "
            f"da_out dims: {da_out.dims}, levsel dims: {levsel.dims}"
        ) from e

    # Final sanity check: shapes must match exactly
    if levsel_broad.shape != tuple(da_out.sizes[d] for d in da_out.dims):
        raise ValueError(
            f"Broadcasted levsel shape mismatch.\n"
            f"Expected shape: {tuple(da_out.sizes[d] for d in da_out.dims)}\n"
            f"Got shape: {levsel_broad.shape}\n"
            f"dims: {da_out.dims}"
        )

    if hasattr(da_data, "compute"):
        logging.warning("Computing da inside pick_level. consider computing beforehand.")
        da_data = da_data.compute()
    lev_axis = da.get_axis_num(lev_dim)
    lev_idx = levsel_broad.data

    if hasattr(lev_idx, "compute"):
        logging.warning("Computing lev_idx inside pick_level. consider computing beforehand.")
        lev_idx = lev_idx.compute()
    lev_idx = np.asarray(lev_idx, dtype=int)

    da_out.data = np.take_along_axis(
        da_data,
        np.expand_dims(lev_idx, axis=lev_axis),
        axis=lev_axis
    ).squeeze(axis=lev_axis)

    return da_out



def compute_ndlut(model_fields : xr.Dataset,
                  aerosol_mmr : xr.DataArray,
                  lut_dset : Optional[xr.Dataset] = None,
                  lut_recipes : Optional[xr.Dataset] = None,
                  deardorff_factor : float = 0.6,
                  w_mean_required : bool = True,
                  w_std_required : bool = True,
                  nlev_below : Optional[int] = None,
                  lwp_frac_threshold : float = 0.98,
                  lev_nan : int = 0
                  ) -> xr.DataArray:
    """
    Compute Nd from aerosol mmr fields using the NDLUT.
    """

    from ..main.config import get_ndlut_path, get_ndlut_rec_path
    import dask

    if lut_dset is None:
        lut_path = get_ndlut_path()
        print(f"Using NDLUT from {lut_path}")
        lut_dset = xr.open_dataset(lut_path)
    if lut_recipes is None:
        lut_rec_path = get_ndlut_rec_path()
        print(f"Using NDLUT recipe from {lut_rec_path}")
        lut_recipes = xr.open_dataset(lut_rec_path)

    dz = model_fields["p_half"].diff(dim="half_lev").rename(half_lev="lev").assign_coords(lev=model_fields["lev"])/GCNST
    rho_air = (model_fields["p"] / (model_fields["t"] * RD))
    lwc = (model_fields["clwc"] * rho_air * dz)
    blh = model_fields["blh"]

    # Cloud base!
    repr_cb_lev = get_cloud_level(
        lwc,
        nlev_below=0,
        lwp_frac_threshold=0.98,
        lev_nan=lev_nan
        )

    # Find BLH model level
    hei = dz.isel(lev=slice(None,None,-1)).cumsum(dim="lev").isel(lev=slice(None,None,-1))
    blh_lev = (hei < blh).argmax(dim="lev") #type: ignore
    blh_lev = blh_lev.where(blh > 0, other=repr_cb_lev) #type: ignore

    # Should do with BLH, use cloud bottom
    if nlev_below is None:
        repr_aero_lev = None
    else:
        repr_aero_lev = get_cloud_level(
            lwc,
            nlev_below=nlev_below,
            lwp_frac_threshold=lwp_frac_threshold,
            lev_nan=lev_nan
        )

    # Compute here
    rho_air, lwc, blh, repr_cb_lev,\
        blh_lev, repr_aero_lev = dask.compute(rho_air, lwc, blh, repr_cb_lev,
                                              blh_lev, repr_aero_lev)

    # We extract aerosols at one level relative to the cloud profile
    if repr_aero_lev is not None:
        aerarr = (aerosol_mmr*rho_air).isel(lev=repr_aero_lev).drop_vars("lev", errors="ignore")
        # aerarr = pick_level(
        #     da=aerosol_mmr*rho_air,
        #     levsel=repr_aero_lev,
        #     lev_dim="lev"
        # ).drop_vars("lev", errors="ignore")
    else:
        aerarr = (aerosol_mmr*rho_air)
    aerarr = aerarr.compute()

    # Integer indexing should guarantee return of views
    aero_dict = {
        aero: aerarr.isel(aero_type=a)
        for a, aero in enumerate(aerarr["aero_type"].values)
    }

    w_mean = pick_level(
        da=model_fields["w"]/(-1.0*rho_air*GCNST),
        levsel=blh_lev,
        lev_dim="lev"
    ).drop_vars("lev", errors="ignore").assign_attrs(units="m/s")

    # TEMPORARY 3D aero fields but this should be done on the fly by the interpolation
    if repr_aero_lev is None:
        raise NotImplementedError("3D ND not implemented yet")
        #w_mean = w_mean.broadcast_like(model_fields["w"])

    w_std = get_wstar(
        blh=model_fields["blh"],
        tsurf=model_fields["2t"],
        heat_flx=model_fields["ishf"]
    )*deardorff_factor

    # TEMPORARY 3D aero fields but this should be done on the fly by the interpolation
    if repr_aero_lev is None:
        raise NotImplementedError("3D ND not implemented yet")
        # w_std = w_std.broadcast_like(model_fields["w"])

    w_mean, w_std = dask.compute(w_mean, w_std)

    from .ndlut import NdLUT, LUTAerosol
    ThisLUTAero = LUTAerosol(
    aerosol_mcon_fields=aero_dict,
    ndlut=NdLUT(
        lut_dset=lut_dset,
        lut_recipes=lut_recipes
        )
    )

    return ThisLUTAero.compute_nd(w_mean=w_mean, w_std=w_std)

def get_wstar(blh : Union[xr.DataArray, np.ndarray],
              tsurf : Union[xr.DataArray, np.ndarray],
              heat_flx : Union[xr.DataArray, np.ndarray],
              q_flx : Union[None, xr.DataArray, np.ndarray] = None
               ) -> xr.DataArray:
    """
     Calculate convective velocity scale (w*) from boundary layer parameters.

    Args:
        blh: Boundary layer height
        tsurf: Surface temperature
        q_flx: Specific humidity flux (Kg/m2/s)
        heat_flx: Heat flux (J/m2/s)

    Returns:
        xr.DataArray: Convective velocity scale (w*)
        assumes rho_air = 1 Kg/m3
    """
    rho_air = 1. # Kg/m3

    if q_flx is not None:
        lat_vap = 2.46*1.e6 #J/Kg
        wst1 = (GCNST/tsurf*blh*(-q_flx * lat_vap / (CPAIR*rho_air)))
    else:
        wst1 = 0.
    wst2 = (GCNST/tsurf*blh*(-heat_flx / (CPAIR*rho_air)))
    if isinstance(wst2, xr.DataArray):
        return (wst1+wst2).clip(min=0)**(1/3) #type: ignore
    else:
        assert not isinstance(wst1, xr.DataArray)
        assert not isinstance(wst2, xr.DataArray)
        return np.clip(wst1+wst2, a_min=0.)**(1/3) #type: ignore


def get_cloud_level(lwc : xr.DataArray, nlev_below : int = 0,
                    lwp_frac_threshold : float = 0.98,
                    lev_nan : float = 0,
                    lev_dim : str = "lev"
                    ) -> xr.DataArray:
    """Get a cloud level index inside the cloud profile
    at required fraction of liquid water path (LWP) from cloud top.
    If no cloud is present, returns lev_nan.
    Could be cloud top (lwp_frac_threshold << 1)
    or cloud base (lwp_frac_threshold ~ 1)
    or any level in between.
    Does not discern multiple cloud layers.

    Args
    ------
    lwc:
        cloud liquid water content profile (normally clwc * rho_air)
    nlev_below:
        number of offset level from threshold level
        (positive means go below, negative means go above)
    lwp_frac_threshold:
        fraction of LWP to identify the "cloud level"
    lev_nan:
        value for no clouds (0 is default, nan do not exist for integers)

    """

    lwp_dn = lwc.cumsum(dim=lev_dim)
    total = lwp_dn.isel({lev_dim: -1}, drop=True)

    lwp_frac = lwp_dn / total

    cond = lwp_frac > lwp_frac_threshold
    mask_2d = (total > 1e-4) & cond.any(dim=lev_dim)
    return xr.where(
        mask_2d,
        (cond.argmax(dim=lev_dim) + nlev_below).clip( #type: ignore
            max=len(lwp_dn[lev_dim])-1
            ),
        lev_nan
    )


def compute_ice_reff_ifs(dset, lats = None, default_re : float = 10.,
                         max_diameter : float = 155,
                         min_diameter : Optional[float] = None):
    """
    Reproduced computations of the ice effective radius as in IFS (RADIP = 3)
    dset xarray.Dataset : must contain the following fields:
                          * cc         -> grid-point cloud cover
                          * ciwc, cswc -> cloud ice/snow water content
                          * lsm        -> land fraction
                          * p          -> full-level pressure
                          * t          -> full-level temperature

    Returns:
    --------
    xarray.DataArray "re_ice" containing ice effective radius in micrometers.
    """

    rre2de : float = 0.64952

    if lats is None:
        if "latitude" in dset.coords:
            lats = dset.latitude
        elif "lat" in dset.coords:
            lats = dset.lat
        else:
            raise ValueError(f"Error in compute_ice_reff_ifs: could not finde latitude coords in dset!")

    rtt = 273.15

    if min_diameter is None:
        min_diameter = 20 + 40*np.cos(np.deg2rad(lats))

    # Do computations where clouds are
    mask        = np.logical_and(dset["cc"] >= 0.001, (dset["ciwc"]+dset["cswc"]) > 0)

    t_safe = dset["t"].clip(min=1.0) # avoid negative temperatures for the computation of reff
    inv_cc_safe = 1/dset["cc"].clip(min=1e-6) # avoid division by zero for the computation of reff

    # Compute liquid and rain water contents
    air_density_gm3 = (1000*dset["p"]/(t_safe*RD)).rename("density")
    iwc_incloud_gm3 = air_density_gm3*(dset["ciwc"]+dset["cswc"])*inv_cc_safe

    temp_celsius    = dset["t"] - rtt

    aiwc = 25.8966 * iwc_incloud_gm3**0.2214
    biwc = 0.7957  * iwc_incloud_gm3**0.2535


    diameter_um = (1.2351 + 0.0105 * temp_celsius) * (aiwc + biwc*(dset["t"] - 83.15))

    diameter_um = diameter_um.clip(min_diameter, max_diameter)
    reff = diameter_um * rre2de

    return xr.where(mask, reff, default_re)

def compute_ccn_ifs(ws :xr.DataArray, lsm : xr.DataArray):
    """
    Ws : absolute 10m wind speed
    lsm : land-sea mask
    """
    landmask = lsm > 0.5
    #seamask  = np.logical_not(landmask)
    wind_lt15 = ws < 15
    #wind_lt30 = ws < 30

    a_par = xr.where(wind_lt15, 0.16, 0.13)
    b_par = xr.where(wind_lt15, 1.45, 1.89)
    qa = np.minimum(np.exp(a_par*ws+b_par), 327)

    c_par = xr.where(landmask, 2.21, 1.2)
    d_par = xr.where(landmask, 0.3,  0.5)
    na = 10**(c_par + d_par*np.log10(qa))

    nd = xr.where(landmask,
                  -2.10e-4*na**2 + 0.568*na - 27.9,
                  -1.15e-3*na**2 + 0.963*na + 5.30
                 )

    return nd

def compute_liquid_reff(dset, nd_fields, wood_correction : bool = True,
                            use_rwc : bool = False,
                            min_reff : float = 4., max_reff : float = 30., cle_reff : float = 2.,
                            min_nd: float = 1., max_nd: float = 3000.,
                            spectr_disp_land : float = 0.69, spectr_disp_sea : float = 0.77):
    """
    Reproduced computation of the liquid effective radius as in IFS (RADLP = 2)
    dset xarray.Dataset : must contain the following fields:
                          * cc         -> grid-point cloud cover
                          * clwc, crwc -> cloud liquid/rain water content
                          * lsm        -> land fraction
                          * p          -> full-level pressure
                          * t          -> full-level temperature

    Returns:
    --------
    xarray.DataArray "re_liquid" containing liquid effective radius in micrometers.
    """

    # Global IFS variables STILL TO FILL
    repscw = 1.e-12 #sec. epsilon for abs. amount in laplace transform
    #replog = 1.e-12 #sec. epsilon for cloud liquid water path

    eps_cc = 1.e-6
    eps_lwc = repscw
    eps_cdnc = 1.e-6

    def _to_numpy(x, target_dims):
        xd = x.transpose(*target_dims).data
        if hasattr(xd, "compute"):
            logging.warning("Computing array inside compute_liquid_reff. consider computing beforehand.")
            return xd.compute()
        else:
            return xd
    ref_da = dset["t"]
    target_dims = ref_da.dims
    t = _to_numpy(dset["t"], target_dims)
    p = _to_numpy(dset["p"], target_dims)
    cc = _to_numpy(dset["cc"], target_dims)
    clwc = _to_numpy(dset["clwc"], target_dims)
    crwc = _to_numpy(dset["crwc"], target_dims)
    lsm = _to_numpy(dset["lsm"].broadcast_like(ref_da), target_dims)

    inv_cc_safe = 1/cc.clip(min=eps_cc)


    mask  = np.logical_and(cc >= 0.001, (clwc+crwc) > 0)

    # Spectral dispersion (land vs. sea) (ZSPECTRAL_DISPERSION)
    #inv_spectr_disp = 1/xr.where(dset["lsm"] > 0.5, spectr_disp_land, spectr_disp_sea)
    inv_spectr_disp = 1/(spectr_disp_sea + (lsm > 0.5) * (spectr_disp_land - spectr_disp_sea))

    ratio  = np.cbrt(0.222*inv_spectr_disp)

    # Compute liquid and rain water contents
    air_density_gm3 = 1000 * (p/(t*RD))#.rename("density")
    # In-cloud mean water contents found by dividing by cloud
    # fraction
    lwc_gm3     = air_density_gm3 * clwc * inv_cc_safe
    rwc_gm3     = air_density_gm3 * crwc * inv_cc_safe

    # Where to get cdcn_fields from
    pot_cdnc = np.clip(
        _to_numpy(
            nd_fields.broadcast_like(ref_da),
            target_dims
        ),
        min_nd, max_nd)

    if wood_correction:
        # Wood's (2000, eq. 19) adjustment to Martin et al's
        # parametrization (ZWOOD_FACTOR)
        rain_ratio  = np.clip(rwc_gm3/lwc_gm3, eps_lwc, None)
        wood_factor = 1 + (lwc_gm3 > repscw) * ((1 + ratio)**(0.666)/\
                                (1+0.2*ratio*rain_ratio) - 1)
    else:
        wood_factor = 1.

    # g m-3 and cm-3 units cancel out with density of water
    # 10^6/(1000*1000); need a factor of 10^6 to convert to
    # microns and cubed root is factor of 100 which appears in
    # equation below
    if use_rwc:
        lwc_gm3 = rwc_gm3
    reff = np.cbrt((0.2387*inv_spectr_disp * lwc_gm3/pot_cdnc))*100*wood_factor

    #replog_mask = re_cubed > replog
    #mask &= re_cubed > replog
    reff = np.clip(reff, min_reff, max_reff)
    reff = cle_reff + mask * (reff - cle_reff)
    reff_out = nd_fields.broadcast_like(ref_da).transpose(*target_dims).copy(data=reff).rename("re_liquid")

    # fill non-cloud areas with clear values and return fields in micrometers
    return reff_out

def compute_sp_from_lnsp(lsnp : xr.DataArray) -> xr.DataArray:
    """Compute surface pressure from its logarithm."""
    new_attrs = lsnp.attrs.copy()
    new_attrs["long_name"] = "Surface pressure"
    new_attrs["units"] = "Pa"

    sp = np.exp(lsnp).assign_attrs(new_attrs).rename("sp") # type: ignore
    if "lev_2" in sp.dims:
        sp = sp.squeeze("lev_2", drop=True)
    return sp

def populate_mlfields(ds : xr.Dataset, keep_p_half : bool = True,
                      lev_dim : str = "lev", half_lev_dim : str = "half_lev"
                      ) -> None:
    """
    Compute 3D meteorological fields from hybrid coordinates.

    Calculates pressure, half-level pressure, air density, and layer pressure
    difference from IFS hybrid coordinate coefficients and surface pressure.

    Args:
        ds: xarray Dataset - Must contain sp, hybm, hyam, hybi, hyai, t
        keep_p_half: should p_half (needed for computing p at full levels) be kept?
    Returns:
        None (modifies dataset in place)

    Side Effects:
        Adds variables to dataset: p, p_half, rho_air, dp
    """

    # check variables needed for computation are present
    required_vars = ["sp", "hybm", "hyam", "hybi", "hyai", "t"]
    missing_vars = []
    for var in required_vars:
        if var not in ds.data_vars:
            missing_vars.append(var)
    assert missing_vars == [], f"Dataset is missing required variables for populate_mlfields: {missing_vars}"

    def _as_float_dtype(da):
        return da.astype(FLOAT_DTYPE, copy=False)

    sp   = _as_float_dtype(ds["sp"])
    t    = _as_float_dtype(ds["t"])
    hybm = _as_float_dtype(ds["hybm"])
    hyam = _as_float_dtype(ds["hyam"])
    hybi = _as_float_dtype(ds["hybi"])
    hyai = _as_float_dtype(ds["hyai"])

    lev_vals = ds[lev_dim].values.astype(int)
    minlev = lev_vals.min().item()
    maxlev = lev_vals.max().item()

    fullplevslice = slice(minlev-1, maxlev)
    halfplevslice = slice(minlev-1, maxlev+1)

    ds["p"] = (sp*hybm.isel(nhym=fullplevslice) +
              hyam.isel(nhym=fullplevslice)).rename(nhym=lev_dim).transpose(*t.dims)

    ds["p_half"] = (sp*hybi.isel(nhyi=halfplevslice) +\
                   hyai.isel(nhyi=halfplevslice) \
                   ).rename(nhyi=half_lev_dim)
    ds = ds.assign_coords(
                   half_lev=np.arange(minlev, maxlev+2, dtype=FLOAT_DTYPE)
                   )

    ds["rho_air"] = ds["p"] / (FLOAT_DTYPE(RD) * t)

    ds["dp"] = ds["p_half"].diff(dim=half_lev_dim).rename({half_lev_dim: lev_dim})

    if not keep_p_half:
        del ds["p_half"]

def compute_hl_temperature(
    p_half: xr.DataArray,
    p_full: xr.DataArray,
    t_full: xr.DataArray,
    skt_sfc: xr.DataArray,
    lev_dim: str = "lev",
    half_lev_dim: str = "half_lev",
):
    """
    Compute half-level temperature using a vectorized formulation with no in-place writes.
    This avoids repeated label-based selections and xarray assignment overhead.
    """

    # Put all 3D fields in a common order with level last for simpler slicing/alignment.
    non_lev_dims = [d for d in p_full.dims if d != lev_dim]
    full_dims = non_lev_dims + [lev_dim]

    p_full = p_full.transpose(*full_dims)
    t_full = t_full.transpose(*full_dims)
    p_half_lev = p_half.rename({half_lev_dim: lev_dim}).transpose(*full_dims)
    skt_sfc = skt_sfc.transpose(*non_lev_dims)

    n_lev = p_full.sizes[lev_dim]
    inner_coord = p_full[lev_dim].isel({lev_dim: slice(0, n_lev - 1)})

    def _as_inner(da: xr.DataArray) -> xr.DataArray:
        return da.assign_coords({lev_dim: inner_coord})

    # All levels but bottom and top half levels
    # k / k+1 pairs on full levels and corresponding half-level k+1/2.
    p_k = _as_inner(p_full.isel({lev_dim: slice(0, n_lev - 1)}))
    p_k1 = _as_inner(p_full.isel({lev_dim: slice(1, n_lev)}))
    t_k = _as_inner(t_full.isel({lev_dim: slice(0, n_lev - 1)}))
    t_k1 = _as_inner(t_full.isel({lev_dim: slice(1, n_lev)}))
    p_h = _as_inner(p_half_lev.isel({lev_dim: slice(1, n_lev)}))

    # Inner half levels (Eq. 2.52) (half levels 0 and 137 are excluded)
    t_inner = (t_k * p_k * (p_k1 - p_h) + t_k1 * p_k1 * (p_h - p_k)) / (p_h * (p_k1 - p_k))

    # Top half level: equal to first full level.
    first_full_lev = p_half_lev[lev_dim].values[0]
    t_top = t_full.isel(
        {lev_dim: [first_full_lev]}
        ).assign_coords({lev_dim: [first_full_lev-1]})

    # Bottom half level: extrapolation using logs of pressure.
    logp_h_n = np.log(p_half_lev.isel({lev_dim: -1}))
    logp_h_nm1 = np.log(p_half_lev.isel({lev_dim: -2}))
    logp_f_n = np.log(p_full.isel({lev_dim: -1}))

    t_f_n = t_full.isel({lev_dim: -1})
    t_h_nm1 = t_inner.isel({lev_dim: -1})

    t_bottom = 0.5 * (
        skt_sfc + t_f_n + (t_f_n - t_h_nm1) / (logp_f_n - logp_h_nm1) * (logp_h_n - logp_f_n)
    )
    t_bottom = t_bottom.expand_dims({lev_dim: [p_half_lev[lev_dim].values[-1]]})

    out = xr.concat([t_top, t_inner, t_bottom], dim=lev_dim)
    out = out.assign_coords({lev_dim: p_half_lev[lev_dim]})
    out = out.rename({lev_dim: half_lev_dim})
    out = out.transpose(*p_half.dims)

    return out.rename("temperature_hl")
