# Installation

ecRadOff is configured with several submodules dependencies, that allow version pinpointing to ensure reproducibility and compatibility. In addition to those, some native Fortran libraries are required and are compiled during the installation process.

## Cloning the repository
Recursing submodules is essential to ensure that all dependencies are available:
```bash
$ git clone --recurse-submodules git@github.com:pandreozzz/ecradoff.git
$ cd ecradoff
```
The `main` branch is aligned with one of the latest stable releases of ecRad (see <https://github.com/ecmwf-ifs/ecrad/tags> for available tags) and corresponds to the most recent available offline version of the ECMWF atmospheric radiation scheme.

### FLOTSAM
A branch of the ecRadOff repository is dedicated to the FLOTSAM version of ecRad. This is used for the simulation of satellite spectral reflectances and radiances and requires additional dependencies that are configured as submodules as well and include FLOTSAM, the "Forward-Lobe Two-Stream Radiance Model" as well as the ADEPT-2 library.

If you want to switch to the FLOTSAM branch:
```bash
$ git checkout flotsam
```


### Submodule recursion

After checking out the branch, make sure all submodules are updated:
```bash
$ git submodule update --init --recursive
```

Alternatively, to avoid manually updating submodules each time, you can optionally set the local config to automatically recurse submodules:
```bash
$ git config submodule.recurse true
```

## Dependencies and python installation

It is recommended to setup a Python virtual environment before proceeding:

```bash
$ python3 -m venv .venv
$ source .venv/bin/activate
```

At this point, proceed with the installation. This will handle all Python dependencies, but some system dependencies may still need to be installed manually. Check out the requirements for the required system dependencies in the `docs/user_guide/system_dependencies.md` file.

### Installation

To install ecRadOff and its Python dependencies:
```bash
$ pip install -v -e .
```
This will call the `setup.py` script, compile and install all submodule dependencies as well as the native ecRadOff libraries. Once the Python dependencies are installed, the `setup.py` can also be invoked directly to recompile and reinstall all required libraries. A possible use case is the usage of debugging flags or switching float precision settings, for example:

```bash
$ python setup.py clean
$ ECRAD_SP=1 python setup.py build_ext
```

For a list of available build options:
```bash
$ python setup.py help
```

#### FLOTSAM
You will need a BLAS/LAPACK library and pass the environment variable `BLAS_LIB` e.g.:

```bash
$ BLAS_LIB=openblas pip install -v -e .
```

To switch between main and flotsam branches:
```bash
$ python setup.py clean
$ git checkout flotsam
$ git submodule update --init --recursive --force
$ BLAS_LIB=openblas pip install -v -e .
```

## Required support fields
Some additional support physical fields are required by ecRad to function. Several lookup tables of aerosol, cloud hydrometeors and files containing solar irradiance values are included in the ecRad package under `ecrad/data`; advanced users can also decide to use their own versions of these files, that need one edit of the `get_tsirr_path()` function inside `src/main/config.py` as well as of the template namelist for the ecRad setup. Soon this will be made more straightforward by switching to JSON configuration files.

Among the fields not included in the ecRad repository are a greenhouse gas climatology and an aerosol climatology and need to be provided by the user. These are NetCDF files and are store under the `data` directory.

### Greenhouse gas climatology
To be placed in the `data/ghg` directory, the filename must be specified accordingly into the `src/main/config.py` file, within the `CONFIGDICT_DEF` dictionary under the key `ghg_file`.

```python
CONFIGDICT_DEF = {
    ...
    "ghg_file": "your_ghg_file.nc",
    ...
}
```

To download a recent publicly available greenhouse gas climatology provided by ECMWF, you can use:
```bash
$ wget "https://confluence.ecmwf.int/download/attachments/70945505/greenhouse_gas_timeseries_CMIP6_SSP370_CFC11equiv_47r1.nc?version=1&modificationDate=1606823141142&api=v2&download=true" -O data/ghg/greenhouse_gas_timeseries_CMIP6_SSP370_CFC11equiv_47r1.nc
```

### Aerosol climatology

Similarly to the greenhouse gas climatology, the aerosol climatology must be placed in the `data/aerosol` directory, and the filename must be specified accordingly into the `src/main/config.py` file, within the `CONFIGDICT_DEF` dictionary under the key `aerosol_file`.

```python
CONFIGDICT_DEF = {
    ...
    "cams_clim_file": "your_aerosol_file.nc",
    ...
}
```
Since IFS cycle CY47R3, the aerosol climatology consists of profiles of aerosol species represented by the IFS-COMPO system used by ECMWF within the Copernicus Atmosphere Monitoring Service (CAMS). For more details, check out the confluence web page on [Aerosol-radiation interactions in the IFS](https://confluence.ecmwf.int/spaces/ECRAD/pages/473845531/Aerosol-radiation+interactions+in+the+IFS).

The monthly climatology introduced by Bozzo et al. (2020) [^1] can be downloaded from:
```bash
$ wget "https://aux.ecmwf.int/ecpds/home/radiation_public/aerosol_climatology/aerosol_cams_climatology_43r3_v2_3D.nc" -O data/cams_clim
```

From the IFS CY49R2 uses a monthly climatology with an additional `epoch` dimension, produced by Stockdale et al. (2023) [^2] that can be downloaded from:
```bash
$ wget "https://aux.ecmwf.int/ecpds/home/radiation_public/aerosol_climatology/aerosol_cams_climatology_49r2_1951-2019_4D.nc" -O data/cams_clim
```

[^1] : Bozzo, A., Benedetti, A., Flemming, J., Kipling, Z., and Rémy, S.: An aerosol climatology for global models based on the tropospheric aerosol scheme in the Integrated Forecasting System of ECMWF, Geosci. Model Dev., 13, 1007–1034, [https://doi.org/10.5194/gmd-13-1007-2020], 2020.

[^2] : Stockdale, T., Senan, R., Hogan, R., Kipling, Z., and Flemming, J.: A new time-varying tropospheric aerosol climatology for the IFS, ECMWF Tech. Rep., [https://www.ecmwf.int/en/elibrary/81558-new-time-varying-tropospheric-aerosol-climatology-ifs], 2024