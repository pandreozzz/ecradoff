"""Module to handle LUT representations and interpolation"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Union
from enum import Enum

from ..main.types import FLOAT_DTYPE

import numpy as np
import xarray as xr
import dask.array as da
from xarray.core.coordinates import DatasetCoordinates, DataArrayCoordinates

XrCoords = Union[DatasetCoordinates, DataArrayCoordinates]

class LUTStruct(Enum):
    """LUT structure types"""

    EXPLODED = 0
    C_FLAT = 1
    F_FLAT = 2

class ActivWType(Enum):
    """LUT types"""

    NONE = 0
    W_MEAN = 1
    W_MEAN_STD = 2
    W_STD = 3


class NdLUT:

    """Class to represent a multi-dimensional LUT and perform interpolation"""
    def __init__(self, lut_dset : xr.Dataset,
                 lut_recipes : xr.Dataset,
                 lut_map_kinds : List[str] = ["num_act"]):
        """Setup the LUT from the given path"""

        # Once finalized, the LUT should not be modified anymore
        self.frozen = False

        #self.lut = lut_dset
        # LUT type and structure
        self._xr_lut = lut_dset.copy()

        renamedic = self.set_activw(self._xr_lut.coords)
        if renamedic != {}:
            self._xr_lut = self._xr_lut.rename(renamedic)
        self.set_lut_structure(self._xr_lut)

        # Aero indices and naming
        self.lut_aeros = [f"aero{spnum}" for spnum in self._xr_lut.spec_num.values]
        # Maps longnames to/from aeroX names
        self.aeronamemap = {
            lut_aero: str(spname)
            for lut_aero, spname in zip(self.lut_aeros, self._xr_lut.name.values)}
        self.pyrcnamemap = {v: k for k, v in self.aeronamemap.items()}

        for lut_map_kind in lut_map_kinds:
            for aero in self.lut_aeros:
                if f"{aero}_{lut_map_kind}" not in self._xr_lut.data_vars:
                    raise ValueError(f"Expected LUT map {aero}_{lut_map_kind} not found!")
        self.lut_maps = lut_map_kinds

        self.lut_aero_nccn_bins = {}
        for lut_aero in self.lut_aeros:
            self.lut_aero_nccn_bins[lut_aero] = self._xr_lut[f"{lut_aero}_nccn"].values

        self.lut_wspeeds_bins = {}
        if self.activw in [ActivWType.W_MEAN_STD, ActivWType.W_STD]:
            self.lut_wspeeds_bins["w_std"] = self._xr_lut["w_std"].values
        if self.activw in [ActivWType.W_MEAN_STD, ActivWType.W_MEAN]:
            self.lut_wspeeds_bins["w_mean"] = self._xr_lut["w_mean"].values

        self.recipes = self._unpack_recipes(lut_recipes)

        self.aerospecs = {}
        self.num_mas_ratio = self._infer_num_mas_ratio(self._xr_lut, lut_recipes)


    def _check_ccn_names(self, lut_recipes : xr.Dataset):
        """Check that the CCN names in the LUT recipes match those in the LUT dataset"""
        assert not self.frozen

        if "ccn_names" not in lut_recipes:
            raise ValueError("No ccn_names found in LUT recipes!")
        if "name" not in self._xr_lut:
            raise ValueError("No name coordinate found in LUT dataset!")
        lut_ccn_names = set(self._xr_lut["name"].values)
        rec_ccn_names = set(lut_recipes["ccn_names"].values)
        if not rec_ccn_names.issubset(lut_ccn_names):
            raise ValueError(f"CCN names in LUT recipes {rec_ccn_names} do not match those in LUT dataset {lut_ccn_names}!")

    @staticmethod
    def _infer_num_mas_ratio(lut_dset : xr.Dataset,
                             lut_recipes : xr.Dataset) -> Dict:
        """Get number-to-mass ratios for the aerosols"""

        from ..main.types import FLOAT_DTYPE

        if "nums_per_mass" in lut_recipes:
            logging.info("[NDLUT] Number-to-mass ratios from values found in LUT recipes")
            return {
                f"aero{spidx+1}": FLOAT_DTYPE(
                    lut_recipes["nums_per_mass"].values[spidx].item()
                )
                for spidx in range(len(lut_dset.spec_num))
            }
        # Infer from PSD definitions if available
        # Check overrides in recipe file
        # Use defaults in LUT file
        else:
            pass
        raise ValueError("Could not find number-to-mass ratios for aerosols in LUT recipes!")

    def set_activw(self, lutcoords : XrCoords):
        """Determine the type of w-activation in the LUT based on its dimensions
        bit of a poor man detection now, might need upgrade at some point
        """
        assert not self.frozen

        renamedic = {}
        w_mean_found = False
        w_prime_found = False
        for d in lutcoords:
            if d == "w_mean":
                w_mean_found = True
            # considered equivalent atm
            elif d == "w_prime" or d == "w_std":
                if w_prime_found:
                    raise ValueError("More than one w_prime or w_std found in LUT dataset!")
                if d == "w_prime":
                    renamedic["w_prime"] = "w_std"
                w_prime_found = True
            # Ignore
            else:
                pass

        if w_mean_found and w_prime_found:
            self.activw = ActivWType.W_MEAN_STD
        elif w_mean_found:
            self.activw = ActivWType.W_MEAN
        elif w_prime_found:
            self.activw = ActivWType.W_STD
        else:
            self.activw = ActivWType.NONE

        return renamedic

    def set_lut_structure(self, lut_dset : xr.Dataset):
        """Determine the structure of the LUT dataset"""
        assert not self.frozen

        # Clearly flat, C or F order?
        if "c_order" in lut_dset.variables:
            if lut_dset["c_order"].values.item():
                self._xr_struct = LUTStruct.C_FLAT
            else:
                self._xr_struct = LUTStruct.F_FLAT
        else:
            self._xr_struct = LUTStruct.EXPLODED

        expected_dims = None
        for var in lut_dset.data_vars:
            if not str(var).startswith("aero"):
                continue

            var_dims = tuple([(str(d), len(lut_dset[d])) for d in lut_dset[var].dims])
            if self._xr_struct != LUTStruct.EXPLODED and len(var_dims) != 1:
                raise ValueError(f"Expected 1D data variable {var} but has dims {var_dims}!")

            if expected_dims is None:
                expected_dims = var_dims
            elif var_dims != expected_dims:
                raise ValueError("All flat data variables should have same dimensions "
                                 f"{expected_dims} but {var} has dims {var_dims}!")

        if self._xr_struct == LUTStruct.EXPLODED:
            assert expected_dims is not None
            self.expl_xrlut_dims = expected_dims
            self.flat_xrlut_dims = (("flatidx",np.asarray([dd[1] for dd in self.expl_xrlut_dims]).prod().item()),)
        else:
            self.flat_xrlut_dims = expected_dims
            expl_xrlut_dims = [(f"aero{spnum}", len(lut_dset[f"aero{spnum}_nccn"]))
                               for spnum in lut_dset.spec_num.values]
            # THIS ORDER IS IMPORTANT (assumes w_mean comes before w_std)
            if self.activw in [ActivWType.W_MEAN_STD, ActivWType.W_STD]:
                expl_xrlut_dims = [("w_std", len(lut_dset["w_std"]))] + expl_xrlut_dims
            if self.activw in [ActivWType.W_MEAN, ActivWType.W_MEAN_STD]:
                expl_xrlut_dims = [("w_mean", len(lut_dset["w_mean"])),] + expl_xrlut_dims
            self.expl_xrlut_dims = tuple(expl_xrlut_dims)


    def _unpack_recipes(self, lut_recipes : xr.Dataset) -> Dict:
        """Unpack the LUT recipes from either adataset or a dict (JSON file)"""
        from ..main.types import FLOAT_DTYPE
        assert not self.frozen
        """Unpacks recipes and checks aero-CCN name consistency with LUT maps"""

        # # It is already a dict, just enforce type
        # if isinstance(lut_recipes, dict):
        #     return {k : {vk : FLOAT_DTYPE(vv)
        #                  for vk, vv in v.items()}
        #             for k, v in lut_recipes.items()}

        # it is xr.Dataset, must be unpacked
        from ..aerodefs.bucket import get_longname_from_short

        recipes = {}
        num_ccn_species = len(lut_recipes["ccn_species"])
        for spidx in range(num_ccn_species):
            spn = spidx + 1
            this_ccn_key = f"aero{spn}"
            this_ccn_nam = lut_recipes["ccn_names"].values[spidx]

            # Check CCN name consistency
            # Stripping is necessary for netCDF classic format
            if self.aeronamemap[this_ccn_key].strip() != this_ccn_nam.strip():
                raise ValueError(f"CCN name mismatch for {this_ccn_key}: LUT dataset has {self.aeronamemap[this_ccn_key]}, but LUT recipes has {this_ccn_nam}!")

            logging.debug(f"[NDLUT] Reading recipe for {this_ccn_key}")
            this_recipe = {}
            for ifs_idx, ifs_shtnam in enumerate(lut_recipes["ifs_species"].values):
                this_ingredient_qty = lut_recipes["ccn_recipes"].isel(
                    ccn_species=spidx,
                    ifs_species=ifs_idx).values.item()
                if this_ingredient_qty > 0:
                    ifs_lngnam = get_longname_from_short(ifs_shtnam)
                    this_recipe[ifs_lngnam] = FLOAT_DTYPE(
                        this_ingredient_qty
                        )
                else:
                    logging.debug(f"[NDLUT] Skipping {ifs_shtnam} as 0-contributor to {this_ccn_key}")
            recipes[this_ccn_key] = this_recipe

        return recipes

    # Getters for dimensions and shapes
    def get_expl_dims(self):
        """Get the exploded lut dimensions"""
        if not self.frozen:
            assert hasattr(self, "expl_xrlut_dims")
            assert self.expl_xrlut_dims is not None
            return tuple(dim[0] for dim in self.expl_xrlut_dims)
        else:
            return self.expl_lut_dims
    def get_expl_shape(self):
        """Get the exploded lut shape"""
        if not self.frozen:
            assert hasattr(self, "expl_xrlut_dims")
            assert self.expl_xrlut_dims is not None
            return tuple(dim[1] for dim in self.expl_xrlut_dims)
        else:
            return self.expl_lut_shape
    def get_flat_dim(self):
        """Get the flat lut dimensions"""
        if not self.frozen:
            assert hasattr(self, "flat_xrlut_dims")
            assert self.flat_xrlut_dims is not None
            assert len(self.flat_xrlut_dims) == 1
            return self.flat_xrlut_dims[0][0]
        else:
            assert len(self.flat_lut_dims) == 1
            return self.flat_lut_dims[0]
    def get_flat_len(self):
        """Get the flat lut shape"""

        if not self.frozen:
            assert hasattr(self, "flat_xrlut_dims")
            assert self.flat_xrlut_dims is not None
            assert len(self.flat_xrlut_dims) == 1
            return self.flat_xrlut_dims[0][1]
        else:
            assert len(self.flat_lut_shape) == 1
            return self.flat_lut_shape[0]
    def get_lut_species_bins(self):
        """Get the LUT species bins in the order of the exploded dimensions"""
        return [self.lut_aero_nccn_bins[dim]
                for dim in self.get_expl_dims()
                if dim in self.lut_aeros]
    def get_lut_wspeeds_bins(self):
        """Get the LUT wspeeds bins in the order of the exploded dimensions"""
        return [self.lut_wspeeds_bins[dim]
                for dim in self.get_expl_dims()
                if dim in self.lut_wspeeds_bins]

    def _get_xr_lut_exploded(self) -> xr.Dataset:
        """Converts the given xr.Dataset to the target structure"""
        if self._xr_struct == LUTStruct.EXPLODED:
            return self._xr_lut

        # Assign multiindex based on the flat dimension
        flat_dim_name = self.get_flat_dim()
        return self._xr_lut.set_index({flat_dim_name: self.get_expl_dims()}).unstack(flat_dim_name)

    def interp_lut(self, interp_points : xr.Dataset):
        assert self.frozen
        pass

    def finalize(self, represented_aeros : List[str],
                 lut_sel : List[str] = ["num_act"]) -> List[str]:
        assert not self.frozen
        # extract subset of LUT
        # convert_lut_structure to flat

        # Something like
        # # Which ones are needed for the LUT recipes?
        # needed_aeros = set()
        # for ccn_name, rec in ndlut.recipes.items():
        #     for aero, qty in rec.items():
        #         if qty > 0:
        #             needed_aeros.add(aero)

        lut_zero_slicer = {}
        for ccn_name, rec in self.recipes.items():
            aeros_to_drop = []
            for aero, qty in rec.items():
                if qty > 0 and aero not in represented_aeros:
                    logging.info(f"Removing {aero} ingredient of ccn {ccn_name}")
                    aeros_to_drop.append(aero)
            if aeros_to_drop != []:
                for aero in aeros_to_drop:
                    rec.pop(aero)
            if all(qty == 0 for qty in rec.values()):
                print(f"Deactivating ccn {ccn_name}")
                lut_zero_slicer[ccn_name] = 0
        for ccn_name in lut_zero_slicer.keys():
            self.recipes.pop(ccn_name)
            self.lut_aeros.remove(ccn_name)

        # Use original flat structure otherwise C order is default
        ini_struct = self._xr_struct
        if ini_struct != LUTStruct.EXPLODED:
            final_struct = self._xr_struct
        else:
            final_struct = LUTStruct.C_FLAT

        # Do we need to slice?
        expl_lut_dims = self.get_expl_dims()
        expl_lut_shape = self.get_expl_shape()
        if lut_zero_slicer != {}:
            ini_struct = LUTStruct.EXPLODED
            xr_lut = self._get_xr_lut_exploded().isel(**lut_zero_slicer, drop=True)
            dimfilter = [dim not in lut_zero_slicer.keys() for dim in expl_lut_dims]
            expl_lut_dims = tuple([dim for dim, keep in zip(expl_lut_dims, dimfilter) if keep])
            expl_lut_shape = tuple([size for size, keep in zip(expl_lut_shape, dimfilter) if keep])
        else:
            xr_lut = self._xr_lut

        self.lut_maps = tuple(
            [f"{aero}_{lut_map_kind}" for aero in self.lut_aeros for lut_map_kind in self.lut_maps]
        )
        self.expl_lut_dims = expl_lut_dims
        self.expl_lut_shape = expl_lut_shape
        self.flat_lut_dims = ("flatidx",)
        self.flat_lut_shape = (np.asarray(expl_lut_shape).prod().item(),)
        logging.info(f"Using LUT maps: {', '.join(self.lut_maps)}")
        # Convert to flat if needed
        if ini_struct != final_struct:
            flat_order = "C" if final_struct == LUTStruct.C_FLAT else "F"
            self.lut = np.ascontiguousarray(
                np.concatenate(
                    [xr_lut[lut_map].transpose(*expl_lut_dims).values.flatten(order=flat_order)[:, None]
                    for lut_map in self.lut_maps],
                axis=-1
                )
            )
        else:
            self.lut = np.ascontiguousarray(
                np.concatenate(
                    [xr_lut[lut_map].values[:, None] for lut_map in self.lut_maps],
                    axis=-1
                )
            )
        self.lut_struct = final_struct

        self.frozen = True
        return represented_aeros

class LUTAerosol:
    """Handle aerosols to use for Nd calculation"""
    import xarray as xr

    def __init__(self,
                 aerosol_mcon_fields : Union[xr.Dataset, Dict[str, xr.DataArray]],
                 ndlut : NdLUT, flat_order : str = "C"):

        from ..aerodefs.bucket import AEROCAMSBUCKET, AeroDesc

        # What aero_vars are needed
        if isinstance(aerosol_mcon_fields, xr.Dataset):
            aero_avail_vars = [str(v) for v in aerosol_mcon_fields.data_vars]
        else:
            aero_avail_vars = list(aerosol_mcon_fields.keys())
        self.ndlut = ndlut
        self.aero_vars = self.ndlut.finalize(
            [str(var) for var in aero_avail_vars
             if str(var) in AEROCAMSBUCKET.keys()])
        self.aero_fields = {var: aerosol_mcon_fields[var] for var in self.aero_vars}
        assert self.aero_fields != {}

        # dimorder and shape
        self.dimorder = tuple(self.aero_fields[self.aero_vars[0]].dims)
        self.fldshape = tuple([len(self.aero_fields[self.aero_vars[0]][dim]) for dim in self.dimorder])

        logging.debug("[NDLUT] Considering the following aerosol variables: "+\
                      ", ".join(self.aero_vars))

    def _get_ccn_ncon(self) -> np.ndarray:
        assert self.ndlut.frozen

        flat_order = None
        if self.ndlut.lut_struct == LUTStruct.C_FLAT:
            flat_order = "C"
        elif self.ndlut.lut_struct == LUTStruct.F_FLAT:
            flat_order = "F"

        assert flat_order is not None

        # out = np.empty((np.prod(self.fldshape), len(self.ndlut.lut_aeros)), dtype=FLOAT_DTYPE)
        # for j, lut_aero in enumerate(self.ndlut.lut_aeros):
        #     out[:, j] = self._get_ccn_ncon_species(lut_aero).transpose(*self.dimorder).values.flatten(order=flat_order)
        # return out

        # Compute once and ensure contiguous
        ccn_out = xr.concat(
                [self._get_ccn_ncon_species(lut_aero) for lut_aero in self.ndlut.lut_aeros],
                dim=xr.IndexVariable("lut_aero", self.ndlut.lut_aeros)
            ).transpose(
                *self.dimorder, "lut_aero"
            ).data

        if hasattr(ccn_out, "compute"):
            ccn_out = da.reshape(
                ccn_out,
                (int(np.prod(self.fldshape)), len(self.ndlut.lut_aeros))
            ).compute()
        else:
            ccn_out = ccn_out.reshape(
                (int(np.prod(self.fldshape)), len(self.ndlut.lut_aeros)),
                order=flat_order
            )
        return np.asarray(ccn_out, dtype=FLOAT_DTYPE, order=flat_order)




    def _get_ccn_ncon_species(self, lut_aero : str):
        """Compute the CCN species by applying recipe"""
        if lut_aero not in self.ndlut.recipes:
            raise ValueError(f"CCN {lut_aero} not found in LUT recipes!")
        else:
            this_ccn = xr.zeros_like(self.aero_fields[self.aero_vars[0]])

        this_recipe = self.ndlut.recipes[lut_aero]
        for aero, qty in this_recipe.items():
            if qty > 0:
                this_ccn += self.aero_fields[aero] * qty
        return this_ccn * self.ndlut.num_mas_ratio[lut_aero]


    def compute_nd(self,
                   w_mean : Optional[xr.DataArray] = None,
                   w_std : Optional[xr.DataArray] = None,
                   map_kind : str = "num_act"):
        """w_mean and w_std if given must have the same dimensions as the aerosol fields"""
        from ..interp.interface import get_lutvals

        required_lut_maps = [f"{aero}_{map_kind}" for aero in self.ndlut.lut_aeros]
        found_idxs = []
        for lut_map in required_lut_maps:
            found_idx = np.where(np.array(self.ndlut.lut_maps) == lut_map)[0]
            if len(found_idx) == 0:
                raise ValueError(f"Required {lut_map} not in LUT!")
            elif len(found_idx) > 1:
                raise ValueError(f"Multiple matches for {lut_map} in LUT!")
            found_idxs.append(found_idx[0])
        lut_stack = np.ascontiguousarray(
            np.concatenate(
                [self.ndlut.lut[:,idx][:,None]
                 for idx in found_idxs],
                axis=-1
            )
        )

        flat_order = None
        if self.ndlut.lut_struct == LUTStruct.C_FLAT:
            flat_order = "C"
        elif self.ndlut.lut_struct == LUTStruct.F_FLAT:
            flat_order = "F"


        lut_wspeeds_bins = self.ndlut.get_lut_wspeeds_bins()
        wspeeds_fields = []
        for w_fld in [w_mean, w_std]:
            if w_fld is not None:

                    #w_fld.transpose(*self.dimorder).values.flatten(order=flat_order)[:, None]
                this_reshaped = w_fld.transpose(*self.dimorder).data
                if hasattr(this_reshaped, "compute"):
                    this_reshaped = da.reshape(
                        this_reshaped,
                        this_reshaped.size,
                        ).compute()
                else:
                    this_reshaped = this_reshaped.reshape(
                        -1,
                        order=flat_order
                    )

                wspeeds_fields.append(
                    np.asarray(this_reshaped, dtype=FLOAT_DTYPE, order=flat_order)[:,None]
                    )

        if wspeeds_fields != [] or self.ndlut.get_lut_wspeeds_bins() != []:
            if len(wspeeds_fields) != len(lut_wspeeds_bins):
                raise ValueError(f"Got {len(wspeeds_fields)} wspeed fields, but LUT expects {len(lut_wspeeds_bins)}!")
            w_present = []
            for w_probe, w_fld in zip(["w_mean", "w_std"], [w_mean, w_std]):
                if w_fld is not None:
                    w_present.append(w_probe)
            nw_present = len(w_present)
            lut_wspdim = self.ndlut.get_expl_dims()[:nw_present]
            if lut_wspdim == tuple(w_present):
                do_transpose_w = False
            elif lut_wspdim == tuple(reversed(w_present)):
                do_transpose_w = True
            else:
                raise ValueError(f"Got wspeed fields {w_present} but LUT expects {lut_wspdim}!")
            wspeeds_fields = np.concatenate(wspeeds_fields, axis=-1)
            if do_transpose_w:
                wspeeds_fields = wspeeds_fields[:, ::-1]
        else:
            wspeeds_fields = np.array([[]])

        species_ccn_ncon_fields = self._get_ccn_ncon()

        out_vals = get_lutvals(
            lut_stack = lut_stack,
            species_fields = species_ccn_ncon_fields,
            lut_species_bins = self.ndlut.get_lut_species_bins(),
            wspeeds_fields = wspeeds_fields,
            lut_wspeeds_bins = self.ndlut.get_lut_wspeeds_bins(),
            chunksize=10000
        )

        return xr.DataArray(
            data=np.sum(out_vals*species_ccn_ncon_fields, axis=1).reshape(self.fldshape),
            dims=self.dimorder,
            coords={dim: self.aero_fields[self.aero_vars[0]].coords[dim] for dim in self.dimorder}
        )

