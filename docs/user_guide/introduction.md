
<img src="../_static/logo_ecradoff.svg"
     alt="ecRadOff logo"
     width="120"
     height="120">

# Introduction
ecRadOff is a Python package meant to assist researchers and developers in working offline with the ECWMF radiation scheme, ecRad. While all radiation calculations are performed by ecRad, ecRadOff constitutes a bridge between typical model data (e.g. from ERA5) and the input required by ecRad.

While ecRad is designed to be used within the ECMWF IFS model, and its offline version exposes an intuitive interface requiring as input a Fortran namelist and a netCDF input file, some preprocessing is necessary to prepare diagnostics required for the radiation calculations, that are typically not available from model output. Notably, these relate to aerosol fields and clouds properties. In fact, most of the ECMWF IFS model configurations use a climatological representation of aerosols, which are therefore not stored as output, while the current single-moment cloud scheme of the IFS requires a dagnostic representation of cloud microphysical properties, i.e. liquid and ice particle number concentrations. For a scientific description of the production of such diagnostics within ecRadOff, please refer to the {doc}`../methods/aerosol` and {doc}`../methods/clouds` sections.

This package aims to streamline the preparation of these inputs, providing a flexible set of tools to generate the necessary fields from available model data, ensuring compatibility with ecRad's requirements as well as reproducibility of the calculations, therefore facilitating offline radiation experiments.
