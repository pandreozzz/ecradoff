"""Numpy version of IFS subroutines"""

from typing import Optional

import numpy as np
import xarray as xr

from ..main.types import FLOAT_DTYPE

# Magic constants
RNAVO = 6.0221367E+23
RKBOL = 1.380658E-23
RMD   = 28.9644
R     = RNAVO*RKBOL
RD    = 1000*R/RMD

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
    cc_safe = dset["cc"].clip(min=1e-6) # avoid division by zero for the computation of reff

    # Compute liquid and rain water contents
    air_density_gm3 = (1000*dset["p"]/(t_safe*RD)).rename("density")
    iwc_incloud_gm3 = air_density_gm3*(dset["ciwc"]+dset["cswc"])/cc_safe

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

def compute_liquid_reff_nrdalp3(dset, cdnc_fields, min_cdnc : float = 1., max_cdnc : float = 3000.,
                               min_reff : float = 4., max_reff : float = 30., cle_reff : float = 2.,
                               spectr_disp_land :float = 0.73, spectr_disp_sea : float = 0.73):
    """
    Reproduced computation of the liquid effective radius from LUT (RADLP = 3)
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

    mask     = np.logical_and(dset["cc"] >= 0.001, (dset["clwc"]+dset["crwc"]) > 0)

    # Spectral dispersion (land vs. sea) (ZSPECTRAL_DISPERSION)
    spectr_disp = xr.where(dset["lsm"] > 0.5, spectr_disp_land, spectr_disp_sea)

    # Compute liquid and rain water contents
    air_density_gm3 = 1000 * (dset["p"]/(dset["t"]*RD)).rename("density")
    # In-cloud mean water contents found by dividing by cloud
    # fraction
    lwc_gm3     = air_density_gm3 * dset["clwc"] / dset["cc"]

    reff = 100*(0.2387/spectr_disp * lwc_gm3/cdnc_fields.clip(min_cdnc, max_cdnc))**0.333
    reff = reff.clip(min_reff, max_reff).rename("re_liquid")

    return xr.where(mask, reff, cle_reff)

def compute_liquid_reff_ifs(dset, ccn_fields = None, wood_correction : bool = True,
                            min_reff : float = 4., max_reff : float = 30., cle_reff : float = 2.,
                            min_ccn: float = 1., max_ccn: float = 3000.,
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

    cc_safe = dset["cc"].clip(min=eps_cc)

    mask     = np.logical_and(dset["cc"] >= 0.001, (dset["clwc"]+dset["crwc"]) > 0)

    # Spectral dispersion (land vs. sea) (ZSPECTRAL_DISPERSION)
    spectr_disp = xr.where(dset["lsm"] > 0.5, spectr_disp_land, spectr_disp_sea)

    ratio       = (0.222/spectr_disp)**0.333

    # Compute liquid and rain water contents
    air_density_gm3 = 1000 * (dset["p"]/(dset["t"]*RD)).rename("density")
    # In-cloud mean water contents found by dividing by cloud
    # fraction
    lwc_gm3     = air_density_gm3 * dset["clwc"] / cc_safe
    rwc_gm3     = air_density_gm3 * dset["crwc"] / cc_safe

    # Where to get cdcn_fields from
    if ccn_fields is None:
        pot_cdnc = 150
    else:
        pot_cdnc = np.clip(ccn_fields, min_ccn, max_ccn)

    if wood_correction:
        # Wood's (2000, eq. 19) adjustment to Martin et al's
        # parametrization (ZWOOD_FACTOR)
        wood_mask   = lwc_gm3 > repscw
        rain_ratio  = rwc_gm3/lwc_gm3.clip(min=eps_lwc)
        wood_factor = xr.where(wood_mask,
                               (1 + rain_ratio)**(0.666)/(1+0.2*ratio*rain_ratio),
                               1)
    else:
        wood_factor = 1.

    # g m-3 and cm-3 units cancel out with density of water
    # 10^6/(1000*1000); need a factor of 10^6 to convert to
    # microns and cubed root is factor of 100 which appears in
    # equation below
    re_cubed = (3 * (lwc_gm3 + rwc_gm3))/(4*np.pi*pot_cdnc*spectr_disp)

    #replog_mask = re_cubed > replog
    #mask &= re_cubed > replog
    reff        = xr.where(mask,
                           wood_factor*100*re_cubed**0.333,
                           min_reff).clip(min_reff, max_reff).rename("re_liquid")

    # fill non-cloud areas with clear values and return fields in micrometers
    return xr.where(mask, reff, cle_reff)

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
              hyam.isel(nhym=fullplevslice)).rename(nhym=lev_dim)

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