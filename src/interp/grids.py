"""Handle grid descriptors for interpolation.
Full numpy interfaces. xArray interfaces are in wrapper.py
* Rectangular (regular) grids
* Reduced Gaussian grids
* Unstructured grids
---
Supported bilinear interpolations:
* Rectangular to rectangular
* Rectangular to/from reduced Gaussian
* Any to unstructured
* None from unstructured
"""
from __future__ import annotations

from typing import Optional, Union
from enum import Enum

import numpy as np

class GridType(Enum):
    """Grid types supported for interpolation."""
    REGULAR = 1
    REDUCED = 2
    UNSTRUC = 3

class GridDesc:
    def __init__(self, gtyp : GridType,
                 lon : np.ndarray, lat : np.ndarray,
                 reduced_pts : Optional[np.ndarray] = None,
                ):
        """Initialise grid descriptor and validate inputs."""
        
        self.gtyp : GridType = gtyp
        self.lon : np.ndarray = lon
        self.lat : np.ndarray = lat
        self.reduced_pts : Union[None, np.ndarray] = reduced_pts 

        if gtyp == GridType.REGULAR:
            n_points = len(lon)*len(lat)

        elif gtyp == GridType.REDUCED:
            if self.reduced_pts is None:
                raise ValueError("reduced_pts must be specified for reduced grids!")
        
            if (len(self.lat) != len(self.lon)) or \
                (len(self.reduced_pts) != len(self.lat)):
                raise ValueError(
                    f"For grid type reduced "+\
                    f"lon ({len(self.lon)}), lat ({len(self.lat)}), reduced_pts ({len(self.reduced_pts)}) "+\
                    "must have the same length!")

            n_points = int(self.reduced_pts.sum())

        elif gtyp == GridType.UNSTRUC:

            if (len(lat) != len(lon)):
                raise ValueError(
                    f"For grid type unstructured "+\
                    f"lon ({len(lon)}) and lat ({len(lat)}) "+\
                    "must have the same length!")
            n_points = len(lon)
        
        self.npts : int = n_points
    
    def get_lat(self):
        """Get lats on data grid points"""
        if self.gtyp == GridType.REDUCED:
            assert self.reduced_pts is not None
            return np.repeat(self.lat, self.reduced_pts)
        else:
            return self.lat
        
    def get_lon(self):
        """Get lons on data grid points"""
        if self.gtyp == GridType.REDUCED:
            assert self.reduced_pts is not None
            return np.concatenate([np.arange(0, 360, 360/rpn)[0:rpn]
                                   for rpn in self.reduced_pts], axis=0)
        else:
            return self.lon


class GriddedData:
    """Data and grids"""
    def __init__(self, data : np.ndarray, grid : GridDesc):
        """Data must have last dimension of size grid.npts"""
        self.data : np.ndarray = data
        self.grid : GridDesc = grid

        if len(self.data.shape) < 1:
            raise ValueError(f"Data must have at least one dimension, got shape {self.data.shape}!")

        # Check grid is consistent with data shape
        if grid.npts != data.shape[-1]:
            raise ValueError(f"Grid points ({grid.npts}) do not match last dim of data ({data.shape})")
        
    def interp2d_to(self, tgt_grid : GridDesc, **interp_kwargs) -> GriddedData:
        """Interpolate field to target grid and return interpolated GriddedField.
        interp_kwargs passed directly to interp2d
        """
        from .interface import interp_2d

        data_shape = self.data.shape
        if len(data_shape) == 1:
            inp_shape = (1,)
            out_shape = (tgt_grid.npts,)
        else:
            inp_shape = (np.prod(np.asarray(data_shape[:-1])),)
            out_shape = data_shape[:-1]+(tgt_grid.npts,)
        inp_shape = inp_shape+(self.grid.npts,)

        interp_data = interp_2d(self.data.reshape(inp_shape), self.grid, tgt_grid, **interp_kwargs).reshape(out_shape)
        
        return GriddedData(interp_data, tgt_grid)
