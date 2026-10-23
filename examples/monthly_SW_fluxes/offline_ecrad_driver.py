"""Driver for ecrad fluxes computations"""

from asyncio import subprocess
import os
import sys

from typing import Tuple

import numpy as np
import xarray as xr

import f90nml

AERO_OPTICS_PATH = "ecradoff/data/optics_lut/aerosol_ifs_49R1_20230725.nc"

def compute_ndlut_reff(model_fields, aerosol_fields):
    """ Computes the liquid droplet effective radius based on aerosol activation diagnostics.

    Parameters
    ----------
    model_fields : xarray.Dataset
        Model variables
    aerosol_fields : xarray.Dataset
        Aerosol variables

    Returns
    -------
    nd_aero, re_liquid_m : xarray.DataArray
        Nd (Cloud liquid droplet number concentration) and liquid effective radius (in meters)
    """
    model_fields_nd = model_fields[var_ndlut].isel(time=tstep).compute()

    nd_aero = compute_ndlut(
        model_fields=model_fields_nd,
        aerosol_mmr=aerosol_fields.data.profile.array["aerosol_mmr"].isel(time=tstep),
        lut_dset=lut_dset, lut_recipes=lut_recipes,
        nlev_below=3, lwp_frac_threshold=0.98, lev_nan=0
        )

    # with cProfile.Profile() as pr:
    model_fields_reff = xr.merge(
        [model_fields_nd[var_common],
        model_fields[var_reff_missing].isel(time=tstep).compute()]
        )

    del model_fields_nd
    re_liquid_m = compute_liquid_reff(
        dset=model_fields_reff,
        nd_fields=nd_aero,
        wood_correction=False,
        use_rwc = False,
    )*1.e-6

    return nd_aero, re_liquid_m

def process_nml(nml_template : f90nml.Namelist,
                aerosol_indexes : Tuple[int]) -> f90nml.Namelist:
    """Processing of the Fortran namelist template for ecRad

    Parameters
    ----------
    nml_template : f90nml.Namelist
        Template namelist
    aerosol_indexes : Tuple[int]
        Indexes of aerosol types to use in the namelist

    Returns
    -------
    f90nml.Namelist
        Processed namelist
    """
    import copy
    # Accessing a specific section and variablie
    var_to_pop = {
        "radiation_driver": [
            "solar_latitude",
            "solar_longitude",
            "sensor_latitude",
            "sensor_longitude",
            "solar_irradiance_override",
        ],
        "radiation" : [
            "sw_albedo_wavelength_bound",
            "i_sw_albedo_index",
            "lw_emiss_wavelength_bound",
            "i_lw_emiss_index",
        ]
    }
    # Meteosat-11 (9.5E, 0.0N)
    var_to_update = {
        "radiation_driver": {
            "iverbose" : 1,
            "do_write_double_precision": False,
            #"nblocksize" : 512,
            # "sensor_latitude": SENSOR_LAT_DEG_N,
            # "sensor_longitude": SENSOR_LON_DEG_E,
        },
        "radiation" : {
            "aerosol_optics_override_file_name": os.path.abspath(AERO_OPTICS_PATH),
            "do_save_spectral_flux": False,
            "iverbose" : 1,
            #"n_cloudy_subcolumns_sw" : 2,
            "n_aerosol_types" : len(aerosol_indexes),
            "i_aerosol_type_map" : list(aerosol_indexes)
        }
    }
    nml_out = copy.deepcopy(nml_template)
    for section, vars in var_to_pop.items():
        for var in vars:
            if var in nml_template[section]:
                print(f"Popping {var}")
                nml_out[section].pop(var)
    for section, vars in var_to_update.items():
        for var, value in vars.items():
            nml_out[section][var] = value
    return nml_out

def run_cmd(cmd):
    """Run the ecRad command"""
    import subprocess

    proc = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    # Read both streams line by line
    assert proc.stdout is not None and proc.stderr is not None, "stdout or stderr is None"
    while True:
        out = proc.stdout.readline()
        err = proc.stderr.readline()

        if out:
            print(out, end="")
        if err:
            print(err, end="", file=sys.stderr)

        # Exit condition
        if out == "" and err == "" and proc.poll() is not None:
            break

    returncode = proc.poll()

    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, cmd)

    return returncode


if __name__ == "__main__":
    # Dask threading
    import dask
    num_dask_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
    dask.config.set(scheduler="threads", num_workers=num_dask_threads)
    client = None

    # Poor man argument parsing
    this_month=sys.argv[1]
    this_day=int(sys.argv[2]) # "1"
    this_grid=sys.argv[3] #"1p0"

    this_exp : str = sys.argv[4] #"v1" # Passed from caller
    timeshift_y : int = int(sys.argv[5]) #30 # Passed from caller
    timefreq_h : int = int(sys.argv[6]) #3 # Passed from caller

    out_dir : str = sys.argv[7] # Passed from caller
#    out_dir = f"/home/papa/Documents/longrange_aerclo_cy49r2b/offline_ecrad_fluxes/{this_exp}_{this_grid}"
    os.makedirs(out_dir, exist_ok=True)

    ini_hour_fc=12

    # File naming convention for the model files
    this_fname = f"era5_ecrad_{this_month}_{this_grid}_eur_{timefreq_h}hourly_*_zarr"


    import ecradoff.src.main.fetch as fetch
    from ecradoff.src.main.config import digest_config, INTERP2D_LIB
    print(INTERP2D_LIB)

    # Digest config (initialise ecRadOff)
    digest_config()

    model_file = os.path.join("/home/papa/Documents/longrange_aerclo_cy49r2b/era5", this_fname)

    ini_dt = np.datetime64(f"{this_month}-{this_day:02d}T{ini_hour_fc:02d}:00:00")
    end_dt = ini_dt + np.timedelta64(1,'D') - np.timedelta64(1,'s')

    model_fields = fetch.model_fields(model_files=model_file).sel(time=slice(ini_dt, end_dt))

    model_fields["p"] = model_fields["p"].compute()

    model_surf_vars_to_compute = ["skt", "lsm", "fal", "10u", "10v", "blh", "2t", "ishf", "tcc"]
    print("Loading surface variables", flush=True)
    for var in model_surf_vars_to_compute:
        if var in model_fields:
            print(var, flush=True)
            model_fields[var] = model_fields[var].compute()

    this_model_data = model_fields[["p"]].assign_coords(model_fields.coords).compute()
    print("Aerosol climatology interpolation", flush=True)
    aerosol_fields = fetch.aerosol_clim(model_data=this_model_data, model_pres_var="p")

    print(f"Aerosol climatology shifted by {timeshift_y:+d} years", flush=True)
    shifted_time = model_fields.time + xr.DataArray(np.timedelta64(timeshift_y,'Y'))
    model_fields_shifted = model_fields.assign_coords(time=shifted_time)
    this_model_data_shifted = this_model_data.assign_coords(time=shifted_time)
    aerosol_fields_shifted = fetch.aerosol_clim(model_data=this_model_data_shifted, model_pres_var="p")
    aerosol_fields_shifted.data.profile.array = aerosol_fields_shifted.data.profile.array.assign_coords(time=model_fields.time)

    from ecradoff.src.tools.ifs import compute_ndlut, compute_liquid_reff
    from ecradoff.src.main.config import get_ndlut_path, get_ndlut_rec_path
    lut_dset = xr.open_dataset(get_ndlut_path())
    lut_recipes = xr.open_dataset(get_ndlut_rec_path())

    from tqdm import tqdm
    # with dask.config.set(scheduler="threads"):
    var_ndlut = ["p", "p_half", "t",
                "w", "clwc",
                "blh", "2t", "ishf"
                ]
    var_reff = ["p", "t", "clwc", "crwc", "cc", "lsm"]
    var_common = set(var_ndlut) & set(var_reff)
    var_reff_missing = set(var_reff) - set(var_ndlut)


    # Pre-load all vars needed by compute_ndlut_reff for all timesteps at once
    nd_vars_to_preload = list((set(var_ndlut) | set(var_reff_missing)) - set(model_surf_vars_to_compute) - {"p"})
    print(f"Pre-loading Nd loop vars: {nd_vars_to_preload}", flush=True)
    for v in nd_vars_to_preload:
        if v in model_fields:
            model_fields[v] = model_fields[v].compute()
        if v in model_fields_shifted:
            model_fields_shifted[v] = model_fields_shifted[v].compute()

    print("Calculating Nd from aerosols and liquid clouds reff")
    nd_aero_list = []
    reff_list = []
    nd_aero_list_shifted = []
    reff_list_shifted = []
    for this_model_fields, this_aerosol_fields, this_nd_aero_list, this_reff_list in zip(
        [model_fields, model_fields_shifted],
        [aerosol_fields, aerosol_fields_shifted],
        [nd_aero_list, nd_aero_list_shifted],
        [reff_list, reff_list_shifted]
        ):
        for tstep in tqdm(range(len(this_model_fields.time))):

            nd_aero, re_liquid_m = compute_ndlut_reff(
                model_fields=this_model_fields,
                aerosol_fields=this_aerosol_fields)

            this_nd_aero_list.append(nd_aero)
            this_reff_list.append(re_liquid_m)

            del nd_aero, re_liquid_m

    nd_aero = xr.concat(nd_aero_list, dim="time").assign_coords(time=model_fields.time)
    re_liquid_m = xr.concat(reff_list, dim="time").assign_coords(time=model_fields.time)
    nd_aero_shifted = xr.concat(nd_aero_list_shifted, dim="time").assign_coords(time=model_fields.time)
    re_liquid_m_shifted = xr.concat(reff_list_shifted, dim="time").assign_coords(time=model_fields.time)

    # Diagnose CCN numbers for ifsdef
    from ecradoff.src.tools.ifs import compute_ccn_ifs
    nd_ifsdef = compute_ccn_ifs(
        ws=np.sqrt(model_fields["10u"]**2 + model_fields["10v"]**2),
        lsm=model_fields["lsm"],
    ).compute()

    # GHG data
    ghg_data = fetch.ghg_data(model_times=model_fields["time"])

    from ecradoff.src.main.config import get_tsirr_path
    from ecradoff.src.tools.eot import Irradiance, cos_satza_deg, sat_azimuthrad_deg
    alltimes = model_fields.time
    ThisIrr = Irradiance(date=alltimes,
                        solar_irr_path=get_tsirr_path())
    map_lons = aerosol_fields.data.profile.grid.xa_lon
    map_lats = aerosol_fields.data.profile.grid.xa_lat
    cos_sza = ThisIrr.mu0_cos_sza_deg(
        phi=map_lats,
        lam=map_lons, zamu0=True,
        night_to_zero=False)

    solar_irr = ThisIrr.solar_irr
    sun_azimuth = None
    cos_satza = None
    sat_azimuth = None

    # Spreader
    delta_freq = np.timedelta64(10, "m")

    all_steps = alltimes.diff(dim="time").values
    if any(all_steps != all_steps[0]):
        raise ValueError("Time steps are not consistent!")
    spreaders = []
    for tstep in alltimes.values:
        timerange = np.arange(
            tstep - all_steps[0]/2,
            tstep + all_steps[0]/2,
            delta_freq)
        timerange = xr.DataArray(data=timerange, dims=["time"], coords={"time": timerange})
        this_spreader = (Irradiance(timerange).mu0_cos_sza_deg(
            phi=map_lats,
            lam=map_lons, zamu0=False,
            night_to_zero=True)/cos_sza.sel(time=tstep, drop=True)).clip(max=1).mean(dim="time").expand_dims(time=[tstep])
        spreaders.append(this_spreader)
    spreader = xr.concat(spreaders, dim="time")


    from ecradoff.src.main.config import get_nml_template_path
    with open(get_nml_template_path().replace(".nml", ".nam"), "r") as fnml:
        nml_template = f90nml.read(fnml)

    from ecradoff.src.ecrad.packer import gen_ecrad_dset

    tmp_root = os.environ.get("ECRADOFF_TMP", None) or os.environ.get("TMPDIR", None)
    import time
    tstamp = time.strftime("%Y%m%d%H%M%S")

    import tempfile
    def prepare_ecrad_rundir(tmp_root : str) -> tempfile.TemporaryDirectory:

        tmp_ecrad = tempfile.TemporaryDirectory(dir=tmp_root)

        if not os.path.exists(os.path.join(tmp_ecrad.name, "data")):
            # create symlink to ecrad/data
            from ecradoff.src.main.config import ECRAD_DATA_DIR
            os.symlink(ECRAD_DATA_DIR, os.path.join(tmp_ecrad.name, "data"))
        if not os.path.exists(os.path.join(tmp_ecrad.name, "ecrad")):
            # create symlink to ecrad
            from ecradoff.src.main.config import ECRAD_DIR
            os.symlink(os.path.join(ECRAD_DIR, "bin/ecrad"), os.path.join(tmp_ecrad.name, "ecrad"))
        return tmp_ecrad

    assert(isinstance(cos_sza, xr.DataArray))
    assert(isinstance(solar_irr, xr.DataArray))

    from ecradoff.src.main.fetch import cast_floats

    labels = ["ndlut", "ndlutshift", "ndonlyshift", "ifsdef"]
    re_liquid_m_dict = {
        "ndlut" : re_liquid_m,
        "ndlutshift": re_liquid_m_shifted,
        "ndonlyshift": re_liquid_m_shifted,
        "ifsdef" : None
    }
    nd_aero_m_dict = {
        "ndlut" : nd_aero,
        "ndlutshift" : nd_aero_shifted,
        "ndonlyshift" : nd_aero_shifted,
        "ifsdef" : nd_ifsdef
    }
    aero_mmr_dict = {
        "ndlut" : aerosol_fields.data.profile.array["aerosol_mmr"],
        "ndlutshift" : aerosol_fields_shifted.data.profile.array["aerosol_mmr"],
        "ndonlyshift" : aerosol_fields.data.profile.array["aerosol_mmr"],
        "ifsdef" : aerosol_fields.data.profile.array["aerosol_mmr"]
    }

    model_ml_vars_to_compute = ["cc", "clwc", "crwc", "cswc",
                                "q", "o3", "p", "t", "p_half", "t_half"]

    import cf_xarray as cfxr
    from concurrent.futures import ThreadPoolExecutor, as_completed
    groups = {}
    cmd_returns = {}
    for lab in labels:
        groups[lab] = []
        cmd_returns[lab] = []

    n_parallel_ecrad = len(labels)
    omp_threads = max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", 8)) // n_parallel_ecrad)
    print(f"Running {n_parallel_ecrad} ecrad labels in parallel, OMP_NUM_THREADS={omp_threads}", flush=True)

    for tstep in range(len(model_fields.time)):
        tmp_ecrad = prepare_ecrad_rundir(tmp_root=tmp_root)
        print(tmp_ecrad.name)

        print(f"Loading vars for time step {tstep}", flush=True)
        this_model_fields = model_fields.isel(time=tstep)
        vars_to_load = [v for v in model_ml_vars_to_compute if v in this_model_fields]
        print(vars_to_load, flush=True)
        computed_block = this_model_fields[vars_to_load].compute()
        this_model_fields = this_model_fields.assign({v: computed_block[v] for v in vars_to_load})
        del computed_block
        print("Done loading model fields", flush=True)

        if cos_satza is not None:
            this_cos_satza = cos_satza.isel(time=tstep) if "time" in cos_satza.dims else  cos_satza
            this_cos_satza = this_cos_satza.compute()
        else:
            this_cos_satza = None

        if sat_azimuth is not None:
            this_sat_azimuth = sat_azimuth.isel(time=tstep) if "time" in sat_azimuth.dims else  sat_azimuth
            this_sat_azimuth = this_sat_azimuth.compute()
        else:
            this_sat_azimuth = None

        if sun_azimuth is not None:
            this_sun_azimuth = sun_azimuth.isel(time=tstep) if "time" in sun_azimuth.dims else  sun_azimuth
            this_sun_azimuth = this_sun_azimuth.compute()
        else:
            this_sun_azimuth = None

        # --- Sequential: prepare all ecrad inputs (not thread-safe due to numba/HDF5) ---
        ecrad_dsets = {}
        ecrad_cmds = {}
        ecrad_path = os.path.join(tmp_ecrad.name, "ecrad")
        for lab in labels:
            print(f"Generating ecrad_dset for {lab} at time step {tstep}", flush=True)

            ecrad_dset = gen_ecrad_dset(model_fields=this_model_fields,
                                        cos_sza=cos_sza.isel(time=tstep),
                                        solar_irradiance=solar_irr.isel(time=tstep).values.item(),
                                        aerosol_mmr=aero_mmr_dict[lab].isel(time=tstep),
                                        ghg_data=ghg_data.isel(time=tstep),
                                        grid_desc=aerosol_fields.data.profile.grid,
                                        cos_satza=this_cos_satza,
                                        sun_azimuth=this_sun_azimuth,
                                        sat_azimuth=this_sat_azimuth,
                                        re_liquid_m=re_liquid_m_dict[lab].isel(time=tstep) \
                                            if re_liquid_m_dict[lab] is not None else None
                                        )
            ecrad_dset = ecrad_dset.assign_coords(
                lon_full=aerosol_fields.data.profile.grid.xa_lon,
                lat_full=aerosol_fields.data.profile.grid.xa_lat,
                time=model_fields.time.isel(time=tstep)
                )

            print(f"Calculating ecrad_dset for {lab}", flush=True)
            ecrad_dset = cast_floats(ecrad_dset).compute()
            ecrad_dsets[lab] = ecrad_dset

            infile = os.path.join(tmp_ecrad.name, f"ecrad_input_{lab}_{tstep:04d}.nc")
            outfile = os.path.join(tmp_ecrad.name, f"ecrad_output_{lab}_{tstep:04d}.nc")
            namfile = os.path.join(tmp_ecrad.name, f"ecrad_flotsam_{lab}_{tstep:04d}.nam")

            ecrad_dset.to_netcdf(infile, format="NETCDF4_CLASSIC")

            nml_out = process_nml(nml_template, aerosol_fields.optics_idx)
            with open(namfile, "w") as fnml:
                f90nml.write(nml_out, fnml)

            groups[lab].append((infile, namfile, outfile))

            subsubdir = os.path.join(tmp_ecrad.name, f"sub_{lab}", "sub2")
            os.makedirs(subsubdir, exist_ok=True)
            cmd = f"cd {subsubdir} && OMP_NUM_THREADS={omp_threads} {ecrad_path} {namfile} {infile} {outfile}"
            ecrad_cmds[lab] = (cmd, infile, outfile, namfile)
            print(cmd, flush=True)

        # --- Parallel: run ecrad subprocesses (subprocesses are fully isolated) ---
        def run_ecrad_subprocess(lab):
            cmd, infile, outfile, namfile = ecrad_cmds[lab]
            try:
                rc = run_cmd(cmd)
            except Exception:
                print(f"Command failed: {cmd}")
                return lab, -1
            return lab, rc

        with ThreadPoolExecutor(max_workers=n_parallel_ecrad) as pool:
            futures = {pool.submit(run_ecrad_subprocess, lab): lab for lab in labels}
            for f in as_completed(futures):
                lab, rc = f.result()
                cmd_returns[lab].append(rc)

        # --- Sequential: process outputs ---
        inp_diag_vars = ["tcc"]
        out_sfc_sw_vars = ["flux_dn_sw", "flux_dn_sw_clear"]
        for lab in labels:
            cmd, infile, outfile, namfile = ecrad_cmds[lab]
            ecrad_dset = ecrad_dsets[lab]
            if not os.path.exists(outfile):
                continue

            out_ds = xr.open_dataset(outfile).rename(column="col")
            out_ds = out_ds[[v for v in out_sfc_sw_vars if v in out_ds]].isel(half_level=-1)
            inp_ds = ecrad_dset
            for var in inp_diag_vars:
                out_ds[var] = inp_ds[var]

            col_decoded = cfxr.decode_compress_to_multi_index(inp_ds).col
            out_ds = out_ds.assign_coords({"col": col_decoded}).unstack("col")
            out_ds = out_ds.assign_coords(time=inp_ds.time).expand_dims("time")
            for v in out_sfc_sw_vars:
                if v in out_ds:
                    out_ds[v] = out_ds[v]*spreader.isel(time=tstep)

            if nd_aero_m_dict[lab] is not None:
                nd_aero_diag = nd_aero_m_dict[lab]
                if "time" in nd_aero_diag.dims:
                    nd_aero_diag = nd_aero_diag.isel(time=tstep)
                # Keep diagnostics dimensionally consistent with other output variables.
                if "time" not in nd_aero_diag.dims:
                    nd_aero_diag = nd_aero_diag.expand_dims(time=[inp_ds.time.values])
                out_ds["nd_aero"] = nd_aero_diag.astype(np.float32)

            tstep_str = str(out_ds.time.values[0])[:13]
            out_ds.to_netcdf(
                os.path.join(out_dir, f"ecrad_output_{lab}_{tstep_str}.nc"),
                unlimited_dims=["time"]
                )

            os.remove(infile)
            os.remove(outfile)
            os.remove(namfile)

        tmp_ecrad.cleanup()
        del this_model_fields






