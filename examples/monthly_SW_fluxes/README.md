# Offline calculation of monthly fluxes

This example performs an offline calculation of monthly radiative fluxes using meteorological fields from ERA5.

Download the example support fields, that include data on single levels (`sfc`) and on model levels (`ml`) from [this public repository](https://sites.ecmwf.int/repository/ecradoff_examples_v1/). The file size is ~1.4GB, therefore you might need to download to a storage location if ecRadOff is on a limited-storage filesystem:
```bash
$ cd <storage_location>
$ wget https://sites.ecmwf.int/repository/ecradoff_examples_v1/era5_ecrad_1984_04_1p0_eur_3hourly.tar
$ tar -xvf era5_ecrad_1984_04_1p0_eur_3hourly.tar && rm era5_ecrad_1984_04_1p0_eur_3hourly.tar
$ cd <ecradoff_location>/examples/monthly_SW_fluxes
$ ln -s <storage_location>/era5_ecrad_1984-04_1p0_eur_3hourly_*_zarr era5/
```