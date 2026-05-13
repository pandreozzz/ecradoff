"""Shared types"""
import numpy as np
import ctypes as ct

from .config import USE_DP

FLOAT_DTYPE = np.float64 if USE_DP else np.float32

DA_DTTYPE = "datetime64[D]"
NS_DTTYPE = "datetime64[ns]"
NS_TDTYPE = "timedelta64[ns]"


CTYPES_TO_NUMPY = {
    ct.c_double: np.float64,
    ct.c_float: np.float32,
    ct.c_int32: np.int32,
}

def ensure_dtype_and_contiguous(arr : np.ndarray, dtype) -> np.ndarray:
    """Array  of the specified dtype and C-contiguous."""
    if isinstance(dtype, type) and dtype in CTYPES_TO_NUMPY:
        dtype = CTYPES_TO_NUMPY[dtype]
    dtype = np.dtype(dtype)

    if arr.dtype != dtype:
        arr = arr.astype(dtype, copy=False)
    if not arr.flags['C_CONTIGUOUS']:
        arr = np.ascontiguousarray(arr)
    return arr