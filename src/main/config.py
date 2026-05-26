import os
from pathlib import Path

from typing import Any, Dict, Optional

USE_DASK = True
USE_DP = False

# ---------------------------------------------------------------------
# Paths (computed from script location)
# ---------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
ECRADOFF_DIR = Path(os.path.join(SCRIPT_DIR, "../../")).resolve()


# Main directories
DATA_DIR = os.path.join(ECRADOFF_DIR, "data")
SRC_DIR = os.path.join(ECRADOFF_DIR, "src")
NMLIST_DIR = os.path.join(ECRADOFF_DIR, "namelists")

# Shared Libraries
SHAREDLIBS_DIR = os.path.join(ECRADOFF_DIR, "locals/lib")
INTERP2D_LIB = os.path.join(SHAREDLIBS_DIR, "interp2d.so")
INTERPVERT_LIB = os.path.join(SHAREDLIBS_DIR, "interpvert.so")
INTERPLUT_LIB = os.path.join(SHAREDLIBS_DIR, "interplut.so")

# Data
CAMS_CLIM_DIR = os.path.join(DATA_DIR, "cams_clim")
ERA5_DATA_DIR = os.path.join(DATA_DIR, "era5")
GHG_SERIES_DIR = os.path.join(DATA_DIR, "ghg")
NDLUT_DIR = os.path.join(DATA_DIR, "ndlut")

# ecRad
ECRAD_DIR = os.path.join(ECRADOFF_DIR, "ecrad")
ECRAD_DATA_DIR = os.path.join(ECRAD_DIR, "data")

# General info
# Used by parser for now
SUPPORTED_GRIDS = ["rectangular", "reduced", "unstructured"]

CONFIGDICT_DEF = {
    "ghg_file": "greenhouse_gas_timeseries_CMIP6_SSP370_CFC11equiv_47r1.nc",
    "cams_clim_file": "aerosol_cams_climatology_49r2_1951-2019_4D.nc",
    "aeosol_optics_version": "49r1",
    "tot_solar_irr_file": "total_solar_irradiance_CMIP6_49r1.nc",
    "nml_template": "config_flotsam_template.nam",
    "ndlut_file": "pyrcel_lut_flat_clim_ad3_20260407.nc",
    "ndlut_rec_file": "pyrcel_lut_flat_clim_ad3_recipe_20260407.nc",
    "nd_from_aerosols" : True,
}

# ---------------------------------------------------------------------
# Global constants
# ---------------------------------------------------------------------

TIME_DIM = "time"
HOR_DIM_1D = "col"
LEV_DIM = "lev"
HLEV_DIM = "half_lev"


# ---------------------------------------------------------------------
# ecRad settings (Flotsam)
# ---------------------------------------------------------------------

ADDITIONALVARS = ["tcc", "cdnc", "bin_num", "ccn_act", "ccn_num", "mu0_spreader",
                  "cos_sensor_zenith_angle", "cos_solar_zenith_angle",
                  "sensor_azimuth_angle", "solar_azimuth_angle", "sea_fraction"]

SENSOR_GEOSTATIONARY : bool = True
SENSOR_LAT_DEG_N : Optional[float] = 0.0
SENSOR_LON_DEG_E : Optional[float] = 9.5

GHG_VARS = ["co2_vmr", "ch4_vmr", "ch4_vmr", "n2o_vmr", "cfc11_vmr", "cfc12_vmr"]

# ---------------------------------------------------------------------
# CONFIGURATION

CONFIGDICT : Dict[str, Any] = {}

def get_ndlut_path() -> str:
    """Get the path to the NDLUT file based on the current configuration."""
    if "ndlut_file" not in CONFIGDICT:
        raise ValueError("ndlut_file not set in CONFIGDICT. Call digest_config first.")
    return os.path.join(NDLUT_DIR, CONFIGDICT["ndlut_file"])

def get_ndlut_rec_path() -> str:
    """Get the path to the NDLUT recipe file based on the current configuration."""
    if "ndlut_rec_file" not in CONFIGDICT:
        raise ValueError("ndlut_rec_file not set in CONFIGDICT. Call digest_config first.")
    return os.path.join(NDLUT_DIR, CONFIGDICT["ndlut_rec_file"])

def get_nml_template_path() -> str:
    """Get the path to the namelist template file based on the current configuration."""
    if "nml_template" not in CONFIGDICT:
        raise ValueError("nml_template not set in CONFIGDICT. Call digest_config first.")
    return os.path.join(NMLIST_DIR, CONFIGDICT["nml_template"])

def get_cams_clim_path() -> str:
    """Get the path to the CAMS climatology file based on the current configuration."""
    if "cams_clim_file" not in CONFIGDICT:
        raise ValueError("cams_clim_file not set in CONFIGDICT. Call digest_config first.")
    return os.path.join(CAMS_CLIM_DIR, CONFIGDICT["cams_clim_file"])

def get_ghg_path() -> str:
    """Get the path to the GHG timeseries file based on the current configuration."""
    if "ghg_file" not in CONFIGDICT:
        raise ValueError("ghg_file not set in CONFIGDICT. Call digest_config first.")
    return os.path.join(GHG_SERIES_DIR, CONFIGDICT["ghg_file"])

def get_tsirr_path() -> str:
    """Get the path to the total solar irradiance file based on the current configuration."""
    if "tot_solar_irr_file" not in CONFIGDICT:
        raise ValueError("tot_solar_irr_file not set in CONFIGDICT. Call digest_config first.")
    return os.path.join(ECRAD_DATA_DIR, CONFIGDICT["tot_solar_irr_file"])

def clear_config() -> None:
    """
    Clears configuration. Further usage needs a new call to digest_config.
    """
    CONFIGDICT = {}

def digest_config(config_path: Optional[str] = None) -> None:
    """
    Load and validate configuration from a JSON file and merge into CONFIGDICT.

    Parameters
    ----------
    config_path : str
        Path to the JSON configuration file.

    Raises
    ------
    ValueError
        If a config key in the file is not present in the default CONFIGDICT.
    """
    import json

    # Read configurations
    if config_path is not None:
        with open(config_path, "r", encoding="utf-8") as config_file:
            cfg_in: Dict[str, Any] = json.load(config_file)
    else:
        cfg_in = {}

    # Validate keys: anything not in defaults (and not prefixed with "other_") is an error
    for key in cfg_in:
        if not key.startswith("other_") and key not in CONFIGDICT_DEF:
            raise ValueError(
                f"config key {key} unknown. Verify spelling errors."
            )

    # Merge (excluding "other_*" keys which are ignored by design)
    # and log CONFIGDICT population
    for key, val in CONFIGDICT_DEF.items():
        if key.startswith("other_"):
            continue
        if key in cfg_in:
            CONFIGDICT[key] = cfg_in[key]
            print(f"Set {key:>20} to {CONFIGDICT[key]}")
        else:
            print(f"Using default value for {key}: {val}")
            CONFIGDICT[key] = val


    print(f"Successfully read configuration from {config_path}")