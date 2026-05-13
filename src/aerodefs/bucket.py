from dataclasses import dataclass
import numpy as np
import xarray as xr

#############################
# CAMS SPECIES NAMES
#############################
@dataclass(frozen=True)
class ClimSpecies:
    longname: str
    shortname: str
    spectype: str
    spechydro: bool
    specbin: int

AEROCAMSBUCKET = {
    "Sea_Salt_bin1"                  : ClimSpecies("Sea_Salt_bin1",
                                                   "SS1", "SS", True, 1),
    "Sea_Salt_bin2"                  : ClimSpecies("Sea_Salt_bin2",
                                                   "SS2", "SS", True, 2),
    "Sea_Salt_bin3"                  : ClimSpecies("Sea_Salt_bin3",
                                                   "SS3", "SS", True, 3),
    "Mineral_Dust_bin1"              : ClimSpecies("Mineral_Dust_bin1",
                                                   "DD1", "DD", False, 1),
    "Mineral_Dust_bin2"              : ClimSpecies("Mineral_Dust_bin2",
                                                   "DD2", "DD", False, 2),
    "Mineral_Dust_bin3"              : ClimSpecies("Mineral_Dust_bin3",
                                                   "DD3", "DD", False, 3),
    "Organic_Matter_hydrophilic"     : ClimSpecies("Organic_Matter_hydrophilic",
                                                   "OMH", "OM", True, 0),
    "Organic_Matter_hydrophobic"     : ClimSpecies("Organic_Matter_hydrophobic",
                                                   "OMN", "OM", False, 0),
    "Black_Carbon_hydrophilic"       : ClimSpecies("Black_Carbon_hydrophilic",
                                                   "BCH", "BC", True, 0),
    "Black_Carbon_hydrophobic"       : ClimSpecies("Black_Carbon_hydrophobic",
                                                   "BCN", "BC", False, 0),
    "Sulfates"                       : ClimSpecies("Sulfates",
                                                   "SU", "SU", True, 0),
    "Nitrate_fine"                   : ClimSpecies("Nitrate_fine",
                                                   "NI1", "NI", True, 1),
    "Nitrate_coarse"                 : ClimSpecies("Nitrate_coarse",
                                                   "NI2", "NI", True, 2),
    "Ammonium"                       : ClimSpecies("Ammonium",
                                                   "AM", "AM", True, 0),
    "Biogenic_Secondary_Organic"     : ClimSpecies("Biogenic_Secondary_Organic",
                                                   "BSO", "OB", True, 0),
    "Anthropogenic_Secondary_Organic": ClimSpecies("Anthropogenic_Secondary_Organic",
                                                   "ASO", "OA", True, 0),
    "Stratospheric_Sulfate"          : ClimSpecies("Stratospheric_Sulfate",
                                                   "SSU", "SSU", False, 0),
}

def get_aero_longname(aerotype, aerohydro, aerobin=None):
    """
        Fetches aerosol long name from properties.
        Returns None if not found.
        If aerobin is not specified and aerotype+aerohydro
        match an aerosol species, the first bin in the bucket is returned
    """
    for long_name,aero in AEROCAMSBUCKET.items():
        if aerotype == aero.spectype and aerohydro == aero.spechydro:
            if aerobin is None or aerobin == aero.specbin:
                return long_name
    return None