#!/bin/sh
# Help for build options

help_message="""
    Environment variables that can be set to control the build:

    NUM_THREADS:
        number of threads to use for build (arg of -j), 0 to use all available, default 1

    PROFILE:
        which profile to use (takes Makefile.PROFILE from ecrad/)

    ECRAD_SP:
        precision for ecrad. default is double. With ecrad flotsam sp build fails because the gauss quadrature (cgqf) code has real(kind=8) hardcoded

    ECRADOFF_DP:
        whether to use double precision for ecradoff. Default is single - no real reason to go for double


    SKIP_ECRAD_BUILD:
        if set to 1, skip building ecrad

    SKIP_ECRADOFF_LIBS_BUILD:
        if set to 1, skip building ecradoff shared libs (in libs/src)
    FLOTSAM_DEBUG:
        if set to 1, build flotsam with debug symbols and no optimisation. Default is 0 (no debug, but also no optimisation).
    FLOTSAM_THREAD_SAFE:
        if set to 1, build flotsam with storage thread-safe flags. Default is 1 (thread safe).
    FLOTSAM_MINIMAL_CHECKS:
        if set to 1, build flotsam with minimal checks (no alias or dimension checking). Default is 0 (full checks).
    FLOTSAM_NAN_INIT:
        if set to 1, build flotsam with NaN initialization. Default is 0 (no NaN initialization).
    FC:
        Fortran compiler to use for building ecradoff libs. If not set, uses the default compiler in the Makefiles.
    BLAS_LIB or OPENBLAS_LIB:
        BLAS/LAPACK library to use for building Adept2 and Flotsam (e.g. "openblas", "mkl", "atlas"). Required for FLOTSAM builds.
    CLEAN_ALL:
        if set to 1, clean all source dirs before build. Default is 0.
"""
echo "$help_message"

