import os
import src.main.fetch as fetch
from src.main.config import digest_config
digest_config()
model_files = "./data/era5/era5_ecrad_1956-10-28_N48_*_zarr"

model_fields = fetch.model_fields(model_files=model_files)

aerosol_fields = fetch.aerosol_clim(model_data=model_fields, model_pres_var="p")