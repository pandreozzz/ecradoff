from dataclasses import dataclass
import numpy as np
import xarray as xr

#############################
# CAMS SPECIES NAMES
#############################

@dataclass(frozen=True)
class AeroPSD:
    pass

@dataclass(frozen=True)
class AeroDesc:
    longname: str
    shortname: str
    spectype: str
    spechydro: bool
    specbin: int
    psd : AeroPSD = AeroPSD()

AEROCAMSBUCKET = {
    "Sea_Salt_bin1"                  : AeroDesc("Sea_Salt_bin1",
                                                   "SS1", "SS", True, 1),
    "Sea_Salt_bin2"                  : AeroDesc("Sea_Salt_bin2",
                                                   "SS2", "SS", True, 2),
    "Sea_Salt_bin3"                  : AeroDesc("Sea_Salt_bin3",
                                                   "SS3", "SS", True, 3),
    "Mineral_Dust_bin1"              : AeroDesc("Mineral_Dust_bin1",
                                                   "DD1", "DD", False, 1),
    "Mineral_Dust_bin2"              : AeroDesc("Mineral_Dust_bin2",
                                                   "DD2", "DD", False, 2),
    "Mineral_Dust_bin3"              : AeroDesc("Mineral_Dust_bin3",
                                                   "DD3", "DD", False, 3),
    "Organic_Matter_hydrophilic"     : AeroDesc("Organic_Matter_hydrophilic",
                                                   "OMH", "OM", True, 0),
    "Organic_Matter_hydrophobic"     : AeroDesc("Organic_Matter_hydrophobic",
                                                   "OMN", "OM", False, 0),
    "Black_Carbon_hydrophilic"       : AeroDesc("Black_Carbon_hydrophilic",
                                                   "BCH", "BC", True, 0),
    "Black_Carbon_hydrophobic"       : AeroDesc("Black_Carbon_hydrophobic",
                                                   "BCN", "BC", False, 0),
    "Sulfates"                       : AeroDesc("Sulfates",
                                                   "SU", "SU", True, 0),
    "Nitrate_fine"                   : AeroDesc("Nitrate_fine",
                                                   "NI1", "NI", True, 1),
    "Nitrate_coarse"                 : AeroDesc("Nitrate_coarse",
                                                   "NI2", "NI", True, 2),
    "Ammonium"                       : AeroDesc("Ammonium",
                                                   "AM", "AM", True, 0),
    "Biogenic_Secondary_Organic"     : AeroDesc("Biogenic_Secondary_Organic",
                                                   "SOB", "OB", True, 0),
    "Anthropogenic_Secondary_Organic": AeroDesc("Anthropogenic_Secondary_Organic",
                                                   "SOA", "OA", True, 0),
    "Stratospheric_Sulfate"          : AeroDesc("Stratospheric_Sulfate",
                                                   "SSU", "SSU", False, 0),
}

def get_longname(aerotype, aerohydro, aerobin=None):
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

def get_longname_from_short(shortname):
    for long_name,aero in AEROCAMSBUCKET.items():
        if shortname == aero.shortname:
            return long_name
    return None