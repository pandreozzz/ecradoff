"""Setup for ecradoff: ecrad + interpolation utilities
Global variables:
ECRADOFF_DP: whether to use double precision for ecradoff. Default is single - no real reason to go for double
PROFILE: which profile to use (takes Makefile.PROFILE from ecrad/)
SKIP_ECRAD_BUILD: if set to 1, skip building ecrad
SKIP_LIBS_BUILD: if set to 1, skip building libs for ECRADOFF
ECRAD_SP: precision for ecrad. default is double. With ecrad flotsam sp build fails because the gauss quadrature (cgqf) code has real(kind=8) hardcoded
"""


from setuptools import setup
from setuptools.command.build_ext import build_ext
import subprocess
import os
import sys

from typing import Optional, List, Union

# Global switches
# do not touch!
USE_FLOTSAM = False

# All absolute paths

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

LOCALS_DIR = os.path.join(
    ROOT_DIR,
    "locals"
    )

ADEPT2_SRC_DIR = os.path.join(
    ROOT_DIR,
    "Adept-2"
)
FLOTSAM_SRC_DIR = os.path.join(
    ROOT_DIR,
    "flotsam"
)
ECRAD_SRC_DIR = os.path.join(
    ROOT_DIR,
    "ecrad"
)

ECRADOFF_SRC_DIR = os.path.join(
    ROOT_DIR,
    "src/libs"
)

def prGreen(s): print("\033[92m {}\033[00m".format(s))

class HelperExt(build_ext):
    """Just print help message about environment variables"""

    def run(self):
        help_message = """
        Environment variables that can be set to control the build:

        NUM_THREADS: number of threads to use for build (arg of -j), 0 to use all available, default 1
        PROFILE: which profile to use (takes Makefile.PROFILE from ecrad/)
        ECRAD_SP: precision for ecrad. default is double. With ecrad flotsam sp build fails because the gauss quadrature (cgqf) code has real(kind=8) hardcoded
        ECRADOFF_DP: whether to use double precision for ecradoff. Default is single - no real reason to go for double
        SKIP_ECRAD_BUILD: if set to 1, skip building ecrad
        SKIP_ECRADOFF_LIBS_BUILD: if set to 1, skip building ecradoff shared libs (in libs/src)
        FLOTSAM_DEBUG: if set to 1, build flotsam with debug symbols and no optimisation. Default is 0 (no debug, but also no optimisation).
        FLOTSAM_THREAD_SAFE: if set to 1, build flotsam with storage thread-safe flags. Default is 1 (thread safe).
        FLOTSAM_MINIMAL_CHECKS: if set to 1, build flotsam with minimal checks (no alias or dimension checking). Default is 0 (full checks).
        FLOTSAM_NAN_INIT: if set to 1, build flotsam with NaN initialization. Default is 0 (no NaN initialization).
        FC: Fortran compiler to use for building ecradoff libs. If not set, uses the default compiler in the Makefiles.
        BLAS_LIB or OPENBLAS_LIB: BLAS/LAPACK library to use for building Adept2 and Flotsam (e.g. "openblas", "mkl", "atlas"). Required for FLOTSAM builds.

        """
        print(help_message)

def clean_all_src_dirs():
    """Call make clean everywgere"""
    for src_dir in [ECRAD_SRC_DIR, FLOTSAM_SRC_DIR, ADEPT2_SRC_DIR, ECRADOFF_SRC_DIR]:
        for cleanop in ['clean', 'distclean']:
            if os.path.exists(src_dir):
                res_clean = subprocess.run(['make', '-C', src_dir, cleanop], capture_output=True, text=True)
                if res_clean.returncode != 0 and cleanop == 'clean': # Not all have distclean
                    print(f"Warning: make {cleanop} failed in {src_dir} with error:\n{res_clean.stderr}")
        else:
            print(f"Cleaned {src_dir}.")

class CleanBuildExt(build_ext):
    def run(self):
        clean_all_src_dirs()
        super().run()

class BuildExt(build_ext):
    """Custom build that compiles ecrad and libs/src Fortran extensions"""

    @staticmethod
    def _write_cmd_logs(log_dir: str, stdout: str, stderr: str, label : str = "build"):
        """Always persist command output for post-mortem debugging."""
        out_file = os.path.join(log_dir, f"{label}_out.log")
        err_file = os.path.join(log_dir, f"{label}_err.log")
        with open(out_file, "w") as f:
            f.write(stdout)
        with open(err_file, "w") as f:
            f.write(stderr)

        print(f"stdout and stderr put into {out_file} and {err_file}.")

    def _easy_logged_run(self, log_dir : str, cmd : Union[str, List[str]],
                         label : str, message : str, shell : bool = False):
        """Just wrap cmd run and logging"""
        if type(cmd) == str:
            cmd_str = cmd
            if not shell:
                raise ValueError("If cmd is a string, shell must be True.")
        else:
            cmd_str = ' '.join(cmd)
        import subprocess
        print(f"Running {message}:\n\t{cmd_str}\n")
        result = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
        self._write_cmd_logs(log_dir, result.stdout, result.stderr, label=label)
        if result.returncode != 0:
            print(result.stderr)
            raise RuntimeError(f"{message} failed with return code {result.returncode}")
        else:
            print(prGreen(f"[SUCCESS] {message} completed."))

    def run(self):

        # Multithreading
        # ecRad Flotsam only works single-thread
        num_threads = int(os.environ.get('NUM_THREADS', '1'))

        # Collect environment variables
        profile = os.environ.get('PROFILE', "")
        # Default SP
        ecradoff_dp = os.environ.get('ECRADOFF_DP', '0')
        # Default DP
        l_ecrad_sp = os.environ.get('ECRAD_SP','0') == '1'

        # Flotsam debug, minimal checks and nan init
        flotsam_debug = os.environ.get('FLOTSAM_DEBUG', '0') == '1'
        flotsam_thread_safe = os.environ.get('FLOTSAM_THREAD_SAFE', '1') == '1'
        flotsam_minimal_checks = os.environ.get('FLOTSAM_MINIMAL_CHECKS', '0') == '1'
        flotsam_nan_init = os.environ.get('FLOTSAM_NAN_INIT', '0') == '1'

        # Fortran compiler and FLAGS
        fc = os.environ.get('FC', "")
        # fcflags = os.environ.get('FFLAGS', "")

        # Is a blas / lapack api available?
        # The user can set BLAS_LIB to control what to use
        # Otherwise it looks for openblas
        blas_lib = os.environ.get(
            'BLAS_LIB',
            os.environ.get(
                'OPENBLAS_LIB',
                None
            )
        )

        # skip_ecrad_build?
        skip_ecrad_build = os.environ.get('SKIP_ECRAD_BUILD', '0') == '1'
        skip_ecradoff_libs_build = os.environ.get('SKIP_ECRADOFF_LIBS_BUILD', '0') == '1'

        # Do we build ecrad (and dependencies)
        if not skip_ecrad_build:

            # -------------------
            # Required system dependencies (depend on the actual version of ecrad)
            # -------------------

            # NetCDF4 with Fortran support required
            cmd_nc4 = ['nc-config', '--has-nc4']
            result_nc4 = subprocess.run(cmd_nc4, capture_output=True, text=True)
            if result_nc4.returncode != 0 or "yes" not in result_nc4.stdout.lower():
                raise RuntimeError("netCDF4 Fortran support is required to build ecrad. Please install netCDF4 with Fortran support and ensure nc-config is in your PATH.")

            if USE_FLOTSAM:
                print("This is a FLOTSAM installation")
                # Did we find a BLAS/LAPACK library?
                if blas_lib is None:
                    raise ValueError("BLAS_LIB or OPENBLAS_LIB environment required to build Adept2")

                # GCC required for flotsam
                cmd_gcc_test = ["gcc", '--version']
                result_gcc_test = subprocess.run(cmd_gcc_test, capture_output=True, text=True)
                if result_gcc_test.returncode != 0:
                    raise RuntimeError(f"GCC is required to build flotsam. {result_gcc_test.stderr}")

                adept2_cpp_flags = self._get_adept2_cppflags(thread_safe=flotsam_thread_safe,
                                                             minimal_checks=flotsam_minimal_checks,
                                                             nan_init=flotsam_nan_init)
                self._build_adept2(ADEPT2_SRC_DIR,
                                   blas_lapack_lib=blas_lib,
                                   adept_cpp_flags=adept2_cpp_flags,
                                   prefix=LOCALS_DIR, debug=flotsam_debug,
                                   num_threads=num_threads
                                   )
                self._build_flotsam(FLOTSAM_SRC_DIR,
                                    cc="gcc", adept2_dir=LOCALS_DIR,
                                    adept_cpp_flags=adept2_cpp_flags,
                                    prefix=LOCALS_DIR, debug=flotsam_debug,
                                    num_threads=num_threads
                                    )
                flotsam_dir = LOCALS_DIR
            else:
                flotsam_dir = None

            self._build_ecrad(ECRAD_SRC_DIR,
                              profile=profile,
                              sp_switch=l_ecrad_sp,
                              flotsam_dir=flotsam_dir)

        if not skip_ecradoff_libs_build:
            self._build_ecradoff_libs(ROOT_DIR, fc=fc,
                             profile=profile,
                             ecradoff_dp=ecradoff_dp)

    def _build_ecradoff_libs(self, root : str, fc : str = "",
                             profile : str = "",
                             ecradoff_dp : str = "0"):
        """Build libs/src Fortran sources into shared libraries using the project Makefile"""

        cmd = ['make', '-C', root]
        cmd.append(f"ECRADOFF_DP={ecradoff_dp}")
        # if fcflags != "":
        #     cmd.append(f"FCFLAGS={fcflags}")
        if profile != "":
            cmd.append(f"PROFILE={profile}")
        elif fc != "":
            cmd.append(f"FC={fc}")

        self._easy_logged_run(root, cmd, label="libs_build", message="ECRADOFF libs build")

    @staticmethod
    def _get_adept2_cppflags(thread_safe : bool,
                             minimal_checks : bool,
                             nan_init : bool) -> str:

        cppflags = ""

        if thread_safe:
            cppflags = "-DADEPT_STORAGE_THREAD_SAFE"

        if minimal_checks:
            print("Building flotsam with minimal checks")
            cppflags += " -DADEPT_NO_ALIAS_CHECKING -DADEPT_NO_DIMENSION_CHECKING"
        else:
            print("Building flotsam with full checks")
            cppflags += " -DADEPT_BOUNDS_CHECKING"

        if nan_init:
            print("Building flotsam with NaN initialization")
            cppflags += " -DADEPT_INIT_REAL_SNAN"

        return cppflags

    def _build_adept2(self, adept2_path,
                      blas_lapack_lib : Optional[str] = None,
                      adept_cpp_flags : str = "",
                      prefix : str = "", skip_checks : bool = True,
                      debug : bool = False,
                      num_threads : int = 1):
        """Build adept2 library using its Makefile. BLAS/LAPACK can be specified via blas_lapack argument (e.g. "openblas", "mkl", "atlas")"""


        if not os.path.exists(os.path.join(adept2_path, "Makefile")):
            # First autoreconf
            cmd_ar = ["autoreconf",'-i', '-f', adept2_path]
            self._easy_logged_run(adept2_path, cmd_ar,
                                label="autoreconf", message="Adept2 autoreconf")


            # Configure
            cmd_cf = f"cd {os.path.join(adept2_path)} && ./configure"
            if debug:
                cxxflags = "-g -O0"
            else:
                cxxflags = "-O0"
            cxxflags = cxxflags + "-Wall -march=native"

            cmd_cf = cmd_cf + f" CXXFLAGS='{cxxflags}'"

            if adept_cpp_flags != "":
                cmd_cf = cmd_cf + f" CPPFLAGS='{adept_cpp_flags}'"
            if prefix != "":
                cmd_cf = cmd_cf + f" --prefix={prefix}"
            if blas_lapack_lib is not None:
                cmd_cf =  cmd_cf + f" --with-blas='{blas_lapack_lib}'"

            cmd_cf = cmd_cf + " && cd -"

            self._easy_logged_run(adept2_path, cmd_cf, shell=True, label="configure", message="Adept2 configure")
        else:
            print("[Info] Adept2 Makefile already exists, skipping autoreconf and configure.")

        # Build
        cmd = ['make', '-C', adept2_path]
        if num_threads > 0:
            cmd.append(f"-j{num_threads}")
        else:
            cmd.append("-j")
        self._easy_logged_run(adept2_path, cmd, label="build", message="Adept2 build")

        # Check
        if not skip_checks:
            cmd = ['make', '-C', adept2_path, 'check']
            self._easy_logged_run(adept2_path, cmd, label="check", message="Adept2 check")

        # Install
        cmd_install = ['make', '-C', adept2_path, 'install']
        self._easy_logged_run(adept2_path, cmd_install, label="install", message="Adept2 install")


    def _build_flotsam(self, flotsam_path, adept2_dir : str,
                       adept_cpp_flags : str = "",
                       prefix : str = "", cc : str = "gcc",
                       cxxflags : str = "g++",
                       debug : bool = False,
                       num_threads : int = 1
                       ):
        """Build flotsam library using its Makefile. CC can be specified via cc argument. adept2_dir should point to the built Adept2 library (e.g. "adept2/build/")"""

        # Check adept2_dir exists (NO OPTION TO USE SYSTEM ADEPT2 INSTALLATION)
        if not os.path.exists(adept2_dir):
            raise ValueError(f"Adept2 directory {adept2_dir} does not exist.")
        adept2_lib = os.path.join(adept2_dir, "lib64")
        if not os.path.exists(adept2_lib):
            adept2_lib = os.path.join(adept2_dir, "lib")
            if not os.path.exists(adept2_lib):
                raise ValueError(f"Adept2 lib or lib64 not found in {adept2_dir}.")
        adept2_include = os.path.join(adept2_dir, "include")
        if not os.path.exists(adept2_include):
            raise ValueError(f"Adept2 include not found in {adept2_dir}.")

        if not os.path.exists(os.path.join(flotsam_path, "Makefile")):
            # Autoreconf
            cmd_ar = ['autoreconf', '-i', '-f', flotsam_path]
            self._easy_logged_run(flotsam_path, cmd_ar, label="autoreconf", message="Flotsam autoreconf")

            # configure

            # CC, CXX, CPPFLAGS and CXXFLAGS
            # TODO: make these environment-dependent
            cxx = "g++"
            cc = "gcc"

            # Start writing the flags
            if debug:
                # debug flag no optimisation
                cxxflags = "-Wall -g -O0"
            else:
                # no debug, but no optimisation either
                cxxflags = "-Wall -O0"
            cxxflags = cxxflags+" -march=native"
            cppflags = f"-I{adept2_include}"
            if adept_cpp_flags != "":
                cppflags = cppflags + " " + adept_cpp_flags

            # LD flags
            ldflags = f"-L{adept2_lib} -Wl,-rpath,{adept2_lib}"
            cmd_cf = f"cd {os.path.join(flotsam_path)} && " +\
                    "./configure " +\
                    f"CC='{cc}' CXX='{cxx}' CXXFLAGS='{cxxflags}' CPPFLAGS='{cppflags}' LDFLAGS='{ldflags}'"
            # where to install Flotsam
            if prefix != "":
                cmd_cf = cmd_cf +f" --prefix={prefix}"
            cmd_cf = cmd_cf + " && cd -"
            self._easy_logged_run(flotsam_path, cmd_cf, shell=True, label="configure", message="Flotsam configure")
        else:
            print("[Info] Flotsam Makefile already exists, skipping autoreconf and configure.")

        cmd = ['make', '-C', flotsam_path]
        if num_threads > 0:
            cmd.append(f"-j{num_threads}")
        else:
            cmd.append("-j")
        self._easy_logged_run(flotsam_path, cmd, label="build", message="Flotsam build")

    def _build_ecrad(self, ecrad_path, profile : str = "",
                     sp_switch : bool = False,
                     flotsam_dir : Optional[str] = None):
        """Build ecrad using specified compiler configuration"""


        cmd = ['make', '-C', ecrad_path]

        # ecrad does not check the value, but just if set
        if sp_switch:
            cmd.append(f"SINGLE_PRECISION=1")

        if flotsam_dir is not None:
            cmd.append(f"FLOTSAM_DIR={flotsam_dir}")

        # Try compiler-specific Makefile_include
        if profile != "":
            makefile_include = os.path.join(ecrad_path, f'Makefile_include.{profile}')
            if os.path.exists(makefile_include):
                cmd.append(f"PROFILE={profile}")
            else:
                available = [f.split("Makefile_include.")[1] for f in os.listdir(ecrad_path) if f.startswith("Makefile_include.")]
                print(f"Warning: Makefile_include.{profile} not found")
                print(f"Available profiles: {', '.join(available)}")
                sys.exit(1)

        print(f"Building ecrad:\n\t{' '.join(cmd)}\n")

        self._easy_logged_run(ecrad_path, cmd, label="build", message="ecRad build")

setup(
    name='ecradoff',
    version='0.1.0',
    description='ECRad radiation scheme with interpolation utilities',
    packages=['src'],
    cmdclass={
        'build_ext': BuildExt,
        'clean': CleanBuildExt,
        'help' : HelperExt
        },
    python_requires='>=3.8',
)
