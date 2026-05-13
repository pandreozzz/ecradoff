"""Manage parsing of command line arguments"""

from typing import List, Union
import numpy as np


def get_parser():
    """Get the argument parser."""

    import argparse
    from .config import SUPPORTED_GRIDS

    parser = argparse.ArgumentParser(prog='Ecrad Offline', description="tbd",
                                     epilog="Something informative")
    parser.add_argument("-n", "--exp-name",
                        type=str, required=True,
                        help="Name of the experiment")

    parser.add_argument("-i", "--model-files",
                        type=str, nargs="+", required=True,
                        help="The IFS output fields to use for offline computations. " +\
                        "Typically 1 ml and 1 sfc file. Only times present in both " +\
                        "datasets are loaded."
                       )
    parser.add_argument("-t", "--times",
                        type=str, nargs="+", default=["all",],
                        help="Model time to use (can be repeated). "+\
                        "Format is YYYY-MM-HH:THH:MM or index starting from 0." +\
                        "Note that ecRad input at each timestep is stored in a separate file." +\
                        "Default is \"all\" - for all model timesteps"
                       )
    parser.add_argument("-g", "--grid",
                        type=str,
                        choices=SUPPORTED_GRIDS+["auto"],
                        default="auto",
                        help="By default grid deduced from model input."
                       )
    parser.add_argument("-a", "--aerosol-version",
                        type=int, default=3,
                        choices=[3, 4, 5],
                        help="Version of aerosol fields to use." +\
                        "v3: CY43R3-CY49R1 Bozzo et al. 2020 3D climatology" +\
                        "v4: CY49R2 4D climatology" +\
                        "v5: CY49R1 prognostic - hydrophilic dust and new PSD for sulfates"
                       )
    parser.add_argument("-lld", "--liquid-lut-dset",
                        type=str, default=None,
                        help="LUT dset to compute cdnc from aerosol fields")

    parser.add_argument("-llr", "--liquid-lut-recipes",
                        type=str, default=None,
                        help="LUT recipes to compute cdnc from aerosol fields")

    parser.add_argument("-zz", "--zamu0-cosine_sz_angle",
                        action="store_true",
                        help="If the cosine of solar zenith angle should go to zero " +\
                        "or to a finite value (zamu0 true, as in ecrad calls from IFS)"
                        )
    parser.add_argument("-si", "--spreader-interval",
                        type=str, default="0s",
                        help="Define a spreader (requires zamu0 true) to simulate " +\
                        "time interpolation of fluxes. Interval is e.g. `15m` for 15 minutes"
                        )
    return parser

def parse_times(time_strs: List[str]) -> Union[None, List[np.datetime64], List[int]]:
    """Parse time strings into a list of datetime objects or indices."""
    from datetime import datetime
    parsed_times = []
    for tstr in time_strs:
        if tstr.lower() == "all":
            return None  # Use None to indicate all times
        try:
            # Try parsing as datetime
            parsed_time = datetime.strptime(tstr, "%Y-%m-%dT%H:%M")
            parsed_times.append(parsed_time)
        except ValueError:
            try:
                # Try parsing as integer index
                parsed_time = int(tstr)
                parsed_times.append(parsed_time)
            except ValueError:
                raise ValueError(f"Invalid time format: {tstr}. " +\
                                 "Expected YYYY-MM-DDTHH:MM or integer index.")
            
    # Sanity check for parsed times
    list_of_ints = False
    mixed_time_fmt_errstr = "Mixed time formats are not allowed. " +\
                            "All times must be either datetime strings or integer indices."
    for p, parst in enumerate(parsed_times):
        if isinstance(parst, int):
            if (p == 0) or list_of_ints:
                list_of_ints = True
            else:
                raise ValueError(mixed_time_fmt_errstr)
                
        else:
            if list_of_ints:
                raise ValueError(mixed_time_fmt_errstr)
        
    return parsed_times
