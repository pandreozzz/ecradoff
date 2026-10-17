# System dependencies

Python requirements are handled via a `pyproject.toml` file and are automatically installed via `pip` (see {doc}`installation`). In addition to these, ecRadOff requires a few system-level dependencies such as compilers and libraries needed to build some package components such as native libraries, ecRad, FLOTSAM and Adept-2. These last packages are all included as Git submodules and are handled automatically by the {doc}`installation` process. Depending on your usage case, you might want to check the [ecRad repository](https://github.com/ecmwf/ecrad), the [FLOTSAM repository](https://github.com/rjhogan/flotsam), and the [Adept-2 repository](https://github.com/rjhogan/Adept-2).

A complete source build requires:
* GNU Make
* a Fortran compiler
* netCDF-C with netCDF-4 support
* netCDF-Fortran
* OpenMP support

### FLOTSAM
For the FLOTSAM configuration, additional requirements are:
* GCC and G++ compilers
* Autoconf/Automake
* a BLAS/LAPACK implementation (e.g. OpenBLAS)

