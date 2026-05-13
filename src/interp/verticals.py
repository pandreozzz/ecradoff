"""Handle vertical interpolation in pressure coordinates.
Full numpy interfaces. xArray interfaces are in wrapper.py
Supported linear interpolation between any of 
* Model levels (sigma levels, pressures are any D)
* Pressure levels (1D, actually subcategory of Model levels)
"""

from __future__ import annotations

from typing import Optional
from enum import Enum

import numpy as np

class InterpWeights:
    """Interpolation weights and indices"""
    def __init__(self, tgtidxs : np.ndarray, weights : np.ndarray,):
        self.tgtidxs : np.ndarray = tgtidxs
        self.weights : np.ndarray = weights

class ProfileData:
    """Data and vertical levels"""

    def __init__(self,
                 data : Optional[np.ndarray],
                 pres : Optional[np.ndarray],
                 nlevs : Optional[int] = None):

        if nlevs is None:
            if data is not None:
                nlevs = data.shape[-1]
            elif pres is not None:
                nlevs = pres.shape[-1]
            else:
                raise ValueError("At least one of data or pres must be specified to infer nlevs!")

        if data is not None and data.shape[-1] != nlevs:
            raise ValueError(f"Last dimension of data ({data.shape[-1]}) must match nlevs ({nlevs})!")
        if pres is not None and pres.shape[-1] != nlevs:
            raise ValueError(f"Last dimension of pres ({pres.shape[-1]}) must match nlevs ({nlevs})!")
        
        self.data = data
        self.pres = pres

        assert nlevs is not None # Ensured above
        self.nlevs : int = nlevs
    
    def get_weights_to(self, ptgt : ProfileData,
                       **kwargs) -> InterpWeights:
        """Interpolate to target pressure levels."""
        from .interface import interp_vert

        assert self.pres is not None
        assert ptgt.pres is not None

        tgtidxs, weights = interp_vert(self.pres, ptgt.pres, **kwargs)
        return InterpWeights(tgtidxs=tgtidxs, weights=weights)
    
    def interp_fld(self, interp_weights : InterpWeights,
                **kwargs) -> np.ndarray:
        """Interpolate fields according to weights."""
        from .interface import interp_fld

        assert self.data is not None

        fdst = interp_fld(self.data,
                          interp_weights.tgtidxs,
                          interp_weights.weights,
                          **kwargs)
        return fdst

