import ctypes as ct
import numpy as np

from typing import Optional, List

from ..main.types import ensure_dtype_and_contiguous
from .grids import GridDesc, GridType
from ..main.config import INTERP2D_LIB, INTERPVERT_LIB, INTERPLUT_LIB, USE_DP

# Types
np_real = np.float64 if USE_DP else np.float32
ct_real = ct.c_double if USE_DP else ct.c_float
np_int = np.int32
ct_int = ct.c_int32
ct_bool = ct.c_bool

c_real_ptr = np.ctypeslib.ndpointer(np_real)
c_int_ptr = np.ctypeslib.ndpointer(np_int)

# 2D interpolation
f_lib_2d = ct.CDLL(INTERP2D_LIB)
# Set types for all interpolation interfaces
head_2dtypes = [c_real_ptr]*6+[ct_int]*7
tail_2dtypes = [ct_int]+[ct_real]+[ct_bool]
# Interface type declarations
f_lib_2d.interp_2d_rec2rec.argtypes = head_2dtypes+tail_2dtypes
f_lib_2d.interp_2d_red2rec.argtypes = head_2dtypes+[c_int_ptr]+tail_2dtypes
f_lib_2d.interp_2d_rec2red.argtypes = head_2dtypes+[c_int_ptr]+[ct_int]+[ct_real]
f_lib_2d.interp_2d_red2red.argtypes = head_2dtypes+[c_int_ptr]*2+tail_2dtypes
f_lib_2d.interp_2d_rec2uns.argtypes = head_2dtypes+tail_2dtypes
f_lib_2d.interp_2d_red2uns.argtypes = head_2dtypes+[c_int_ptr]+tail_2dtypes

# Vertical interpolation
f_lib_vert = ct.CDLL(INTERPVERT_LIB)
f_lib_vert.interp.argtypes = [c_real_ptr,]*2+[c_int_ptr, c_real_ptr]+[ct_int,]*6
f_lib_vert.interp_fld.argtypes = [c_real_ptr,]*2+[c_int_ptr, c_real_ptr]+[ct_int,]*6

# LUT interpolation
f_lib_lut = ct.CDLL(INTERPLUT_LIB)
f_lib_lut.getvals.argtypes = [c_real_ptr,] + [ct_int,]*6 + [c_real_ptr,] *2 + [c_int_ptr, c_real_ptr] + [c_real_ptr,] * 2 + [c_int_ptr,] + [ct_int, ct_bool]
# f_lib_lut.getvals_nowsp.argtypes =  [c_real_ptr,] + [ct_int,]*6 + [c_real_ptr,]*2 + [c_int_ptr, c_real_ptr] + [ct_int, ct_bool]


def interp_2d(fsrc, srcgrid : GridDesc, tgtgrid : GridDesc,
              lat_row : bool = True,
              chunk_size_max : Optional[int] = None,
              abs_tolerance : float = 1.e-3):
    """ 2D bilinear interpolation from srcgrid to tgtgrid

    Args:
        fsrc: Source fields, shape (nflds, nxysrc) and nxysrc must match srcgrid.npts
        srcgrid: Source grid descriptor
        tgtgrid: Target grid descriptor
        lat_row: is (lat,lon) (true) or (lon,lat) (false) for the flattened fsrc fields?
        chunk_size_max: Maximum chunk size for interpolation (very relevant for unstructured targets)
        abs_tolerance: Absolute tolerance for interpolation (passed to C library)
    """

    nflds, nxysrc = fsrc.shape
    if (nxysrc != srcgrid.npts):
        raise ValueError(f"Incompatible number of gridpoints for fsrc ({nxysrc}) and grid specs ({srcgrid.npts})!")

    outshape = (nflds, tgtgrid.npts)

    # Enforce contiguity
    fdst = ensure_dtype_and_contiguous(np.empty(outshape, dtype=np_real), np_real)

    # Chunk size
    if chunk_size_max is None:
        if tgtgrid.gtyp == GridType.REGULAR:
            chunk_size_max = int(np.ceil(len(tgtgrid.lat)/128))
        elif tgtgrid.gtyp == GridType.REDUCED:
            chunk_size_max = 1
        elif tgtgrid.gtyp == GridType.UNSTRUC:
            chunk_size_max = int(np.ceil(tgtgrid.npts/512))
        else:
            raise ValueError(f"Unknown grid type {tgtgrid.gtyp}!")


    # Args shared by all interpolation interfaces
    head_args = [
        ensure_dtype_and_contiguous(fsrc, np_real), fdst,]


    # ! lat_row - meaning
    # ! Whether the dimension ordering is (Fortran writing, column-major)
    # ! x and y are assigned according to Fortran column-major convention:
    # ! (x,y)=(lon,lat): lat_row=true or (x,y)=(lat,lon): lat_row=false
    # ! for the flattened fields on regular grid. For reduced
    # ! and unstructured grids, lat_row is ignored

    # ysrc, xsrc, ydst, xdst
    if lat_row: # y->lat, x->lon
        head_args += [
        ensure_dtype_and_contiguous(srcgrid.lat, np_real),
        ensure_dtype_and_contiguous(srcgrid.lon, np_real),
        ensure_dtype_and_contiguous(tgtgrid.lat, np_real),
        ensure_dtype_and_contiguous(tgtgrid.lon, np_real)
        ]
    else: # y->lon, x->lat
        head_args += [
        ensure_dtype_and_contiguous(srcgrid.lon, np_real),
        ensure_dtype_and_contiguous(srcgrid.lat, np_real),
        ensure_dtype_and_contiguous(tgtgrid.lon, np_real),
        ensure_dtype_and_contiguous(tgtgrid.lat, np_real)
        ]

    if lat_row: # y->lat, x->lon
        head_args += [ct_int(nflds),
        ct_int(len(srcgrid.lat)), ct_int(len(srcgrid.lon)), ct_int(srcgrid.npts),
        ct_int(len(tgtgrid.lat)), ct_int(len(tgtgrid.lon)), ct_int(tgtgrid.npts)
        ]
    else: # y->lon, x->lat
        head_args += [ct_int(nflds),
        ct_int(len(srcgrid.lon)), ct_int(len(srcgrid.lat)), ct_int(srcgrid.npts),
        ct_int(len(tgtgrid.lon)), ct_int(len(tgtgrid.lat)), ct_int(tgtgrid.npts)
        ]

    tail_args = [
        ct_int(chunk_size_max), ct_real(abs_tolerance), ct_bool(lat_row)
    ]

    if (tgtgrid.gtyp == GridType.REGULAR):
        if (srcgrid.gtyp == GridType.REGULAR):
            f_lib_2d.interp_2d_rec2rec(
                *head_args,
                *tail_args
                )
        if (srcgrid.gtyp == GridType.REDUCED):
            assert srcgrid.reduced_pts is not None
            f_lib_2d.interp_2d_red2rec(
                *head_args,
                ensure_dtype_and_contiguous(srcgrid.reduced_pts, np_int),
                *tail_args
                )
    elif (tgtgrid.gtyp == GridType.REDUCED):
        assert tgtgrid.reduced_pts is not None
        if (srcgrid.gtyp == GridType.REGULAR):
            f_lib_2d.interp_2d_rec2red(
                *head_args,
                ensure_dtype_and_contiguous(tgtgrid.reduced_pts, np_int),
                *tail_args
                )
        if (srcgrid.gtyp == GridType.REDUCED):
            assert srcgrid.reduced_pts is not None
            f_lib_2d.interp_2d_red2red(
                *head_args,
                ensure_dtype_and_contiguous(srcgrid.reduced_pts, np_int),
                ensure_dtype_and_contiguous(tgtgrid.reduced_pts, np_int),
                *tail_args
                )
    elif tgtgrid.gtyp == GridType.UNSTRUC:
        if (srcgrid.gtyp == GridType.REGULAR):
            f_lib_2d.interp_2d_rec2uns(
                *head_args,
                *tail_args
                )
        if (srcgrid.gtyp == GridType.REDUCED):
            assert srcgrid.reduced_pts is not None
            f_lib_2d.interp_2d_red2uns(
                *head_args,
                ensure_dtype_and_contiguous(srcgrid.reduced_pts, np_int),
                *tail_args
                )

    return fdst

def interp_vert(psrc : np.ndarray, ptgt : np.ndarray,
                chunk_size_max : int = 1000):
    """
        Interpolates psrc to ptgt
        psrc(ncom,nsrc,nlevsrc)
        ptgt(ncom,ntgt,nlevtgt)

        returns
        tgtlevs(ncom,nsrc,ntgt,nlevtgt), weights(ncom,nsrc,ntgt,nlevtgt)
    """

    ncom,ntgt,nlevtgt  = ptgt.shape
    ncom2,nsrc,nlevsrc = psrc.shape

    if (ncom != ncom2):
        raise ValueError(
                "Different common dimension size between source and target!"+\
                f"({ncom2},{nsrc}) and ({ncom},{nsrc})")

    outshape = (ncom,nsrc,ntgt,nlevtgt)

    tgtlevs = np.empty(outshape, dtype=np_int)
    weights = np.empty(outshape, dtype=np_real)

    # Enforce contiguity
    tgtlevs = ensure_dtype_and_contiguous(tgtlevs, np_int)
    weights = ensure_dtype_and_contiguous(weights, np_real)

    f_lib_vert.interp(
        ensure_dtype_and_contiguous(psrc, np_real),
        ensure_dtype_and_contiguous(ptgt, np_real),
        tgtlevs, weights,
        ct_int(ncom), ct_int(nsrc), ct_int(ntgt),
        ct_int(nlevsrc), ct_int(nlevtgt), ct_int(chunk_size_max)
        )
    del psrc, ptgt

    return tgtlevs, weights

def interp_fld(fsrc : np.ndarray,
               tgtlevs : np.ndarray, weights : np.ndarray,
               chunk_size_max : int = 1000):
    """
        Interpolates fields according to weights
        fsrc(ncom,nsrc,nlevsrc)
        tgtlevs(ncom,ntgt,nlevtgt)
        weights(ncom,ntgt,nlevtgt)

        returns
        fdst(ncom,nsrc,ntgt,nlevtgt)
    """

    ncom,nsrc,nlevsrc = fsrc.shape

    ncom2,ntgt,nlevtgt   = tgtlevs.shape
    ncom3,ntgt2,nlevtgt2 = weights.shape

    if (ncom != ncom2) or (ncom2 != ncom3) or (ntgt != ntgt2) or (nlevtgt != nlevtgt2):
        raise ValueError(
                "Some dimensions are incompatible!!"+\
                f"fsrc(ncom={ncom},nsrc={nsrc},nlevsrc{nlevsrc})"+\
                f"tgtlevs(ncom={ncom2},ntgt={ntgt},nlevtgt{nlevtgt})"+\
                f"tgtlevs(ncom={ncom3},ntgt={ntgt2},nlevtgt{nlevtgt2})")

    fdstshape = (ncom,nsrc,ntgt,nlevtgt)

    # Enforce contiguity
    fdst = ensure_dtype_and_contiguous(np.empty(fdstshape, dtype=np_real), np_real)
    tgtlevs_cont = np.ascontiguousarray(tgtlevs)
    weights_cont = np.ascontiguousarray(weights)

    f_lib_vert.interp_fld(
            ensure_dtype_and_contiguous(fsrc, np_real),
            fdst,
            tgtlevs_cont.astype(np_int),
            weights_cont.astype(np_real),
            ct_int(ncom), ct_int(nsrc), ct_int(ntgt),
            ct_int(nlevsrc), ct_int(nlevtgt), ct_int(chunk_size_max)
            )
    del tgtlevs_cont, weights_cont

    return ensure_dtype_and_contiguous(fdst, np_real)

def get_lutvals(lut_stack : np.ndarray,
                species_fields : np.ndarray,
                lut_species_bins : List[np.ndarray],
                wspeeds_fields : Optional[np.ndarray] = None,
                lut_wspeeds_bins : Optional[List[np.ndarray]] = None,
                c_order : bool = True,
                chunksize : int = 10000):
    """Interpolate the LUT stack to fields values.

    lut_stack : np.ndarray (map_size, nmaps)
        The LUT stack, where each column contains the LUT map
    _fields : np.ndarray (nvals, nfields=nwspeed | nspecies)
        The field values to interpolate the LUT stack to.
        wspeeds_fields is optional

    Returns:
    val_out_data : np.ndarray (nvals, nmaps)
        The interpolated LUT values for each field value in _fields
    """

    map_size, nmaps = lut_stack.shape
    nvals, nspec = species_fields.shape
    if len(lut_species_bins) != nspec:
        raise ValueError(f"Got {len(lut_species_bins)} species bins, but expected {nspec}!")

    if lut_wspeeds_bins is None:
        lut_wspeeds_bins = []

    if wspeeds_fields is None:
        wspeeds_fields = np.empty((nvals, 0), dtype=np_real)

    if len(wspeeds_fields) > 0:
        nvals2, nwspeed = wspeeds_fields.shape
        if nvals2 != nvals:
            raise ValueError(f"Got wspeeds with nvals = {nvals2}, but expected nvals = {nvals}!")
        if len(lut_wspeeds_bins) != nwspeed:
            raise ValueError(f"Got {len(lut_wspeeds_bins)} wspeed bins, but expected {nwspeed}!")
        nwspeedbins = ensure_dtype_and_contiguous(np.array([len(wsp) for wsp in lut_wspeeds_bins], dtype=np_int), np_int)
    else:
        nwspeed = 0
        nwspeedbins = np.array([], dtype=np_int)

    nspecbins = ensure_dtype_and_contiguous(
        np.array([len(spec) for spec in lut_species_bins], np_int),
        np_int)

    max_lut_bins = int(max(nspecbins.max(), nwspeedbins.max(initial=0)))

    lut_species_bins_matrix = ensure_dtype_and_contiguous(np.full((nspec, max_lut_bins), np.nan, dtype=np_real), np_real)

    for i, spec in enumerate(lut_species_bins):
        lut_species_bins_matrix[i, :len(spec)] = spec.astype(np_real, copy=False)

    lut_wspeeds_bins_matrix = ensure_dtype_and_contiguous(np.full((nwspeed, max_lut_bins), np.nan, dtype=np_real), np_real) if nwspeed > 0 else np.empty((0,0), dtype=np_real)
    for i, wsp in enumerate(lut_wspeeds_bins):
        lut_wspeeds_bins_matrix[i, :len(wsp)] = wsp.astype(np_real, copy=False)

    val_out_data = np.ascontiguousarray(np.empty((nvals, nmaps), dtype=np_real))


    f_lib_lut.getvals(
        ensure_dtype_and_contiguous(lut_stack, np_real),
        ct_int(nmaps), ct_int(map_size),
        ct_int(nwspeed), ct_int(nspec), ct_int(nvals), ct_int(max_lut_bins),
        ensure_dtype_and_contiguous(species_fields, np_real),
        ensure_dtype_and_contiguous(lut_species_bins_matrix, np_real),
        ensure_dtype_and_contiguous(nspecbins, np_int),
        val_out_data,
        ensure_dtype_and_contiguous(wspeeds_fields, np_real),
        ensure_dtype_and_contiguous(lut_wspeeds_bins_matrix, np_real),
        ensure_dtype_and_contiguous(nwspeedbins, np_int),
        ct_int(chunksize), ct_bool(c_order)
    )

    return val_out_data
