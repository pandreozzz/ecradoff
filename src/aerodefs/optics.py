
from dataclasses import dataclass

from typing import Dict, Any, Tuple

import xarray as xr

PHILIC_DESCR = \
"""
1: Sea salt, bin 1, 0.03-0.5 micron, OPAC
2: Sea salt, bin 2, 0.50-5.0 micron, OPAC
3: Sea salt, bin 3, 5.0-20.0 micron, OPAC
4: Hydrophilic organic matter, OPAC
5: Ammonium sulfate (for sulfate), GACP Lacis et al https://gacp.giss.nasa.gov/data_sets/
6: Secondary organic aerosol - biogenic, Moise et al 2015
7: Secondary organic aerosol - anthropogenic, Moise et al 2015
8: Fine mode Ammonium sulfate (for ammonia), GACP Lacis et al https://gacp.giss.nasa.gov/data_sets/
9: Fine mode Nitrate, GLOMAP
10: Coarse mode Nitrate, GLOMAP
11: Hydrophilic organic matter, Brown et al 2018
12: Sulfate, GACP Lacis et al https://gacp.giss.nasa.gov/data_sets/
13: Sulfate, GACP Lacis et al https://gacp.giss.nasa.gov/data_sets/ with modified size distribution
14: Desert dust, bin 1, 0.03-0.55 micron, Composite-Philic Non-Sphere-Scaling-Kandler (Balkanski et 2007 , Di Baggio 2017, Ryder et al 2019)
15: Desert dust, bin 2, 0.55-0.90 micron, Composite-Philic Non-Sphere-Scaling-Kandler (Balkanski el 2007 , Di Baggio 2017, Ryder et al 2019)
16: Desert dust, bin 3, 0.90-20.0 micron, Composite-Philic Non-Sphere-Scaling-Kandler (Balkanski el 2007 , Di Baggio 2017, Ryder et al 2019)
17: Desert dust, bin 1, 0.03-0.55 micron, Composite-Philic (Balkanski et 2007 , Di Baggio 2017, Ryder et al 2019)
18: Desert dust, bin 2, 0.55-0.90 micron, Composite-Philic (Balkanski el 2007 , Di Baggio 2017, Ryder et al 2019)
19: Desert dust, bin 3, 0.90-20.0 micron, Composite-Philic (Balkanski el 2007 , Di Baggio 2017, Ryder et al 2019)
"""
PHOBIC_DESCR = \
"""
1: Desert dust, bin 1, 0.03-0.55 micron, (SW) Dubovik et al. 2002 (LW) Fouquart et al. 1987
2: Desert dust, bin 2, 0.55-0.90 micron, (SW) Dubovik et al. 2002 (LW) Fouquart et al. 1987
3: Desert dust, bin 3, 0.90-20.0 micron, (SW) Dubovik et al. 2002 (LW) Fouquart et al. 1987
4: Desert dust, bin 1, 0.03-0.55 micron, Fouquart et al 1987
5: Desert dust, bin 2, 0.55-0.90 micron, Fouquart et al 1987
6: Desert dust, bin 3, 0.90-20.0 micron, Fouquart et al 1987
7: Desert dust, bin 1, 0.03-0.55 micron, Woodward 2001, Table 2
8: Desert dust, bin 2, 0.55-0.90 micron, Woodward 2001, Table 2
9: Desert dust, bin 3, 0.90-20.0 micron, Woodward 2001, Table 2
10: Hydrophobic organic matter, OPAC (hydrophilic at RH=20%)
11: Black carbon, OPAC
12: Black carbon, Bond and Bergstrom 2006
13: Black carbon, Stier et al 2007
14: Stratospheric sulfate (hydrophilic ammonium sulfate at RH 20%-30%)
15: Desert dust, bin 1, 0.03-0.55 micron, Composite (Balkanski et 2007 , Di Baggio 2017, Ryder et al 2019)
16: Desert dust, bin 2, 0.55-0.90 micron, Composite (Balkanski el 2007 , Di Baggio 2017, Ryder et al 2019)
17: Desert dust, bin 3, 0.90-20.0 micron, Composite (Balkanski el 2007 , Di Baggio 2017, Ryder et al 2019)
18: Hydrophobic organic matter, Brown et al 2018 (hydrophilic at RH=20%)
19: Black carbon, Williams 2007
"""


@dataclass(frozen=True)
class OpticIndex():
    index: int
    hydrophilic: bool

OPTICS_AERO_MAP : Dict[str, Dict[str, OpticIndex]]= {}
###
## Prognostic 43r3 and Bozzo climatology
###
OPTICS_AERO_MAP["43r3"] = {
    "Sea_Salt_bin1"                   : OpticIndex(1, True),
    "Sea_Salt_bin2"                   : OpticIndex(2, True),
    "Sea_Salt_bin3"                   : OpticIndex(3, True),
    "Mineral_Dust_bin1"               : OpticIndex(7, False), # Composite-Phobic
    "Mineral_Dust_bin2"               : OpticIndex(8, False), # Composite-Phobic
    "Mineral_Dust_bin3"               : OpticIndex(9, False), # Composite-Phobic
    "Organic_Matter_hydrophilic"      : OpticIndex(4, True),
    "Organic_Matter_hydrophobic"      : OpticIndex(4, False),
    "Black_Carbon_hydrophilic"        : OpticIndex(11, False),
    "Black_Carbon_hydrophobic"        : OpticIndex(11, False),
    "Sulfates"                        : OpticIndex(5, True)
}

###
## Prognostic 48r1
###
OPTICS_AERO_MAP["48r1"] = {
    **OPTICS_AERO_MAP["43r3"].copy(),
    **{
    "Nitrate_fine"                    : OpticIndex(9, True),
    "Nitrate_coarse"                  : OpticIndex(10, True),
    "Ammonium"                        : OpticIndex(8, True),  
    "Biogenic_Secondary_Organic"      : OpticIndex(6, True),  
    "Anthropogenic_Secondary_Organic" : OpticIndex(7, True),
    "Stratospheric_Sulfate"           : OpticIndex(14, False),  
    }
}
# Composite phobic dust
OPTICS_AERO_MAP["48r1"]["Mineral Dust_bin1"] = OpticIndex(15, False)
OPTICS_AERO_MAP["48r1"]["Mineral Dust_bin2"] = OpticIndex(16, False)
OPTICS_AERO_MAP["48r1"]["Mineral Dust_bin3"] = OpticIndex(17, False)
# Brown OM
OPTICS_AERO_MAP["48r1"]["Organic_Matter_hydrophilic"] = OpticIndex(11, True)
OPTICS_AERO_MAP["48r1"]["Organic_Matter_hydrophobic"] = OpticIndex(10, False)

###
## IFS-COMPO 48r1-based 4D climatology (Tim's) deployed in IFS 49R2
## has inconsistent sulfates and uses the new PSD
###
OPTICS_AERO_MAP["48r1_4dclim"]  = OPTICS_AERO_MAP["48r1"].copy()
OPTICS_AERO_MAP["48r1_4dclim"]["Sulfates"] = OpticIndex(13, True)

###
## Prognostic 49r1
###
OPTICS_AERO_MAP["49r1"] = OPTICS_AERO_MAP["48r1"].copy()

# New PSD for Sulfates
OPTICS_AERO_MAP["48r1_4dclim"]["Sulfates"] = OpticIndex(13, True)

# Hydrophilic Dust
OPTICS_AERO_MAP["49r1"]["Mineral Dust_bin1"] = OpticIndex(14, True)
OPTICS_AERO_MAP["49r1"]["Mineral Dust_bin2"] = OpticIndex(15, True)
OPTICS_AERO_MAP["49r1"]["Mineral Dust_bin3"] = OpticIndex(16, True)

# Bond BC
OPTICS_AERO_MAP["49r1"]["Black_Carbon_hydrophobic"] = OpticIndex(12, False)
OPTICS_AERO_MAP["49r1"]["Black_Carbon_hydrophilic"] = OpticIndex(12, False)


@dataclass
class OpticsData:
    data : Any
    optics_var: Tuple[str]
    optics_idx: Tuple[int]

    @property
    def optics_map(self) -> Dict[str, int]:
        assert len(self.optics_var) == len(self.optics_idx), "optics_var and optics_idx must have the same length"
        return {var: idx for var, idx in zip(self.optics_var, self.optics_idx)}

def get_optical_aerosol(cams_dset: xr.Dataset,
                        aero_version: str,
                        verbose: bool = False) -> OpticsData:
    '''
    generates dataset with radiatively-active species and maps to optical properties
    '''
    from .bucket import AEROCAMSBUCKET

    aero_map = []
    aero_typ = []
    var_to_drop = []

    aeroptics = OPTICS_AERO_MAP[aero_version]

    non_aero_vars = []
    optics_map = {}
    optics_var = []
    optics_idx = []
    for var in cams_dset.data_vars.keys():
        # Ignore non-aerosol variables
        if var not in AEROCAMSBUCKET:
            non_aero_vars.append(var)
            continue
        if var not in aeroptics:
            if verbose:
                print(f"Warning: Aerosol {var} in dataset not used in optics version {aero_version}. " +\
                      f"Available aerosols for optics version {aero_version}: {list(aeroptics.keys())}.")
            var_to_drop.append(var)
            continue
        this_optics = aeroptics[var]

        optics_var.append(var)
        optics_idx.append(-this_optics.index if this_optics.hydrophilic else this_optics.index)
    
    if optics_var == []:
        raise ValueError(f"No aerosols in dataset match optics version {aero_version}. " +\
                         f"Available aerosols: {list(cams_dset.data_vars.keys())}. " +\
                         f"Expected aerosols for optics version {aero_version}: {list(aeroptics.keys())}.")
    var_sel = non_aero_vars+optics_var

    return OpticsData(cams_dset[var_sel], tuple(optics_var), tuple(optics_idx))