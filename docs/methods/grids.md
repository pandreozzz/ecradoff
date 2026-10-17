# Grid interpolation interface

## Supported grids
Provides xArray/NumPy interpolation interfaces for the following grid types:
* Rectangular (regular) grids
* Reduced Gaussian grids
* Unstructured grids


## Supported interpolations
Supported interpolations:
* Rectangular → rectangular
* Rectangular ↔ reduced Gaussian
* Any → unstructured

Interpolation from unstructured grids not supported.

## Lookup-table accessor
The lookup-table (LUT) accessor is also included in the interpolation module.