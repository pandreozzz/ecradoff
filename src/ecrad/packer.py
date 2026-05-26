"""Pack the ecrad input dataset"""

from typing import Optional

import xarray as xr

from ..interp.wrapper import ArrayGridDesc
from ..interp.grids import GridDesc, GridType
from ..main.types import FLOAT_DTYPE
from ..main.config import ADDITIONALVARS, GHG_VARS, HOR_DIM_1D, LEV_DIM, HLEV_DIM

def gen_ecrad_dset(model_fields : xr.Dataset, solar_irradiance : float,
                   cos_sza : xr.DataArray, aerosol_mmr: xr.DataArray,
                   ghg_data: xr.Dataset, grid_desc : ArrayGridDesc,
                   cos_satza : Optional[xr.DataArray] = None,
                   sun_azimuth : Optional[xr.DataArray] = None,
                   sat_azimuth : Optional[xr.DataArray] = None,
                   re_liquid_m : Optional[xr.DataArray] = None,
                   iseed : Optional[xr.DataArray] = None
                  ) -> xr.Dataset:
    import numpy as np
    import cf_xarray as cfxr

    from ..tools import ifs as ifst
    
    # Shortcuts
    cc = model_fields["cc"]
    clwc = model_fields["clwc"]
    ciwc = model_fields["ciwc"]
    cswc = model_fields["cswc"]
    q = model_fields["q"]
    o3 = model_fields["o3"]

    p_full = model_fields["p"]
    t_full = model_fields["t"]

    p_half = model_fields["p_half"]
    t_half = getattr(model_fields, "t_half", None)

    # Surface
    skt = model_fields["skt"]
    lsm = model_fields["lsm"]
    fal = model_fields["fal"]
    u10 = model_fields["10u"]
    v10 = model_fields["10v"]

    for var in GHG_VARS:
        if var not in ghg_data.data_vars:
            raise ValueError(f"GHG variable '{var}' not found in ghg_data dataset: {list(ghg_data.data_vars.keys())}!")

    # Take a required 2D field as mold
    if iseed is None:
        iseed = xr.DataArray(
            data=np.random.randint(
                low=0,
                high=np.iinfo(np.int32).max,
                size=skt.shape),
            dims=skt.dims, coords=skt.coords,
            name="iseed"
        ).astype(np.int32)

    # ecRad ice content is sum of ice and snow water content
    q_ice = ciwc + cswc

    # Compute ml t
    if t_half is None:
        t_half = ifst.compute_hl_temperature(
            p_half=p_half, p_full=p_full,
            t_full=t_full, skt_sfc=skt,
            lev_dim=LEV_DIM, half_lev_dim=HLEV_DIM
            ).compute()

    # Compute liquid and ice effective radius in meters (requires p in model_fields)
    # if lut_recipes and lut_dset:
    #     pass
        # print("re_liquid from LUT!")
        # llut = True
        # import lut_tools as lutt
        # from ifs_tools import RD
        # dens_full = p_full/(RD*model_fields["t"].astype(np.float32))
        # aero_mcon_for_lut = (dens_full*aerosol_mmr).sel(lev=[NRADLLEV,])
        # with_extra_fields = "bin_num" in ADDITIONALVARS or "ccn_num" in ADDITIONALVARS or "ccn_act" in ADDITIONALVARS
        # lut_species_mcon = lutt.compute_lut_species_from_ifs_species(lut_recipes, aero_mcon_for_lut)
        # ccn_fields = lutt.compute_cdnc(lut_species_mcon, lut_dset,
        #                               with_extra_fields=with_extra_fields).squeeze(drop=True)

        # re_liquid = ifst.compute_liquid_reff_ifslut(dset=model_fields, cdnc_fields=ccn_fields["cdnc"]).compute()*1.e-6
        # #re_liquid=xr.ones_like(p_full)
    # else:
    #     llut = False
    if re_liquid_m is None:
        re_liquid_m = ifst.compute_liquid_reff(
            dset=model_fields,
            nd_fields=ifst.compute_ccn_ifs(
                ws=np.sqrt(u10**2+v10**2), #type: ignore
                lsm=lsm
            ),
            min_reff=4, max_reff=30,
            min_nd=10, max_nd=3000,
            spectr_disp_land=0.69,
            spectr_disp_sea=0.77
        ).compute()*1.e-6

    # Ice particle size
    re_ice_m = ifst.compute_ice_reff_ifs(
        dset=model_fields,
        lats=grid_desc.xa_lat
        ).compute()*1.e-6

    data_vars = {
        "solar_irradiance"       : FLOAT_DTYPE(solar_irradiance),
        "skin_temperature"       : skt,
        "cos_solar_zenith_angle" : cos_sza,
        "sw_albedo"              : fal, # use spectrally constant albedo
        "lw_emissivity"          : xr.full_like(fal, 0.99),
        "iseed"                  : iseed,
        "pressure_hl"            : p_half,
        "temperature_hl"         : t_half,
        "h2o_mmr"                : q,
        "o3_mmr"                 : o3,
        "o2_vmr"                 : FLOAT_DTYPE(0.20944),
        #"co_vmr"                 : xr.DataArray(np.float32(1.e-6)),
        "hcfc22_vmr"             : FLOAT_DTYPE(0.), #xr.DataArray(np.float32(240*1.e-12)),
        "ccl4_vmr"               : FLOAT_DTYPE(0.), #xr.DataArray(np.float32(77*1.e-12)),
        "no2_vmr"                : FLOAT_DTYPE(1.e-8),
        "aerosol_mmr"            : aerosol_mmr,
        "q_liquid"               : clwc,
        "q_ice"                  : q_ice,
        "re_liquid"              : re_liquid_m,
        "re_ice"                 : re_ice_m,
        "cloud_fraction"         : cc,
        **{var: FLOAT_DTYPE(ghg_data[var])
           for var in GHG_VARS}
    }

    if cos_satza is not None:
        data_vars["cos_sensor_zenith_angle"] = cos_satza
    if sun_azimuth is not None:
        data_vars["solar_azimuth_angle"] = sun_azimuth
    if sat_azimuth is not None:
        data_vars["sensor_azimuth_angle"] = sat_azimuth

    for additional_var in ADDITIONALVARS:
        if additional_var in model_fields.variables:
            data_vars = {**data_vars,**{additional_var : model_fields[additional_var]}}
            continue
        # match additional_var:
        #     case "lre":
        #         data_vars["lre"] = re_liquid
        #     case "ire":
        #         data_vars["ire"] = re_ice
        #     case "cdnc":
        #         data_vars["cdnc"] = ccn_fields["cdnc"] if llut else ccn_fields
        #     case "bin_num" if llut:
        #         data_vars["bin_num"] = xr.concat([ccn_fields[f"aero{i}_bin"].assign_coords(lutspec=i) for i in range(1,5)], dim="lutspec")
        #     case "ccn_num" if llut:
        #         data_vars["ccn_num"] = xr.concat([ccn_fields[f"aero{i}_ccn"].assign_coords(lutspec=i) for i in range(1,5)], dim="lutspec")
        #     case "ccn_act" if llut:
        #         data_vars["ccn_act"] = xr.concat([ccn_fields[f"aero{i}_act"].assign_coords(lutspec=i) for i in range(1,5)], dim="lutspec")
        #     case _:
        #         print(f"Warning in packer: Ignoring {additional_var} from ADDITIONALVARS")

    xtr_coords = {}
    if grid_desc.gtyp == GridType.REDUCED:
        red_coordims = grid_desc.coordims
        xtr_coords["lat"] = (red_coordims, grid_desc.lat)
        xtr_coords["reduced_points"] = (red_coordims, grid_desc.reduced_pts)
    ecrad_dset = xr.Dataset(
        data_vars=data_vars,
    ).assign_attrs({
         "Peculiarities" : "None"
        }
    ).assign_coords(**xtr_coords)
    if "reduced_points" in model_fields:
        ecrad_dset = ecrad_dset.assign_coords(lat=model_fields.lat,
                                              reduced_points=model_fields.reduced_points)
    if grid_desc.gtyp == GridType.REGULAR:
        #ecrad_dset = ecrad_dset.transpose("lon", "lat", ..., "half_level", "lev")
        ecrad_dset = cfxr.encode_multi_index_as_compress( # type: ignore
            ecrad_dset.stack({HOR_DIM_1D: grid_desc.coordims}),
            HOR_DIM_1D
            )

    # THIS is the order the ecRad expects.
    ecrad_dset = ecrad_dset.transpose(HOR_DIM_1D, ..., HLEV_DIM, LEV_DIM)

    for var_to_drop in ("time", "fcdate", "levaux"):
        if var_to_drop in ecrad_dset:
            ecrad_dset = ecrad_dset.drop_vars(var_to_drop)

    return ecrad_dset
