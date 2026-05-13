"""The ecradoff driver - launch this and cross your fingers!"""

def main():
    """Driver function for the ecradoff"""
    from .parser import get_parser, parse_times
    import fetch
    from .config import digest_config, TIME_DIM
    
    parser = get_parser()
    args = parser.parse_args()

    # Eventually handle json config in the future
    digest_config()

    # Fetch model fields
    model_fields = fetch.model_fields(
        args.model_files,
        time_sel=parse_times(args.times))
    
    aerosol_fields = fetch.aerosol_clim(model_fields)

    lut_dset = None
    lut_recipes = None

    ghg_data = fetch.ghg_data(model_fields[TIME_DIM])