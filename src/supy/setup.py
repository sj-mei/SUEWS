import os
from signal import raise_signal
from setuptools import setup
import json
from pathlib import Path

# write version info using git commit
import subprocess
import warnings
import re

from setuptools import Distribution
from numpy.distutils.core import Extension, setup
import platform
import glob
import os
from pathlib import Path
import subprocess
import shutil
# from nonstopf2py import f2py

# ISRELEASED = True
# # if a release, use strict requirement for supy-driver; otehrwise, use a loose requirement
# DRIVER_REQ = "supy_driver==2021a15" if ISRELEASED else "supy_driver"
# FULLVERSION += '.dev'

# pipe = None
# p_fn_ver = Path("./supy/supy_version.json")

# # force remove the version info file
# flag_dirty = False

# for cmd in ["git", "/usr/bin/git", "git.cmd"]:

#     try:
#         pipe = subprocess.Popen(
#             [cmd, "describe", "--tags", "--match", "2[0-9]*", "--dirty=-dirty"],
#             stdout=subprocess.PIPE,
#         )
#         (sout, serr) = pipe.communicate()
#         # parse version info from git
#         list_str_ver = sout.decode("utf-8").strip().split("-")

#         if list_str_ver[-1].lower() == "dirty":
#             flag_dirty = True
#             # remove the "dirty" part from version info list
#             list_str_ver = list_str_ver[:-1]

#         ver_main = list_str_ver[0]
#         print("ver_main", ver_main)
#         if len(list_str_ver) > 1:
#             ver_post = list_str_ver[1]
#             ver_git_commit = list_str_ver[2]
#         else:
#             ver_post = ""
#             ver_git_commit = ""

#         # save version info to json file
#         p_fn_ver.unlink(missing_ok=True)
#         with open(p_fn_ver, "w") as f:
#             json.dump(
#                 {
#                     "version": ver_main + ("" if ISRELEASED else ".dev"),
#                     "iter": ver_post,
#                     "git_commit": ver_git_commit
#                     + (
#                         "-dirty"
#                         if (flag_dirty and len(ver_git_commit) > 0)
#                         else (
#                             "dirty" if (flag_dirty and len(ver_git_commit) == 0) else ""
#                         )
#                     ),
#                 },
#                 f,
#             )
#         if pipe.returncode == 0:
#             print(f"in {cmd}, git version info saved to", p_fn_ver)
#             break
#     except Exception as e:
#         pass

# if pipe is None or pipe.returncode != 0:
#     # no git, or not in git dir

#     if p_fn_ver.exists():
#         warnings.warn(
#             f"WARNING: Couldn't get git revision, using existing {p_fn_ver.as_posix()}"
#         )
#         write_version = False
#     else:
#         warnings.warn(
#             "WARNING: Couldn't get git revision, using generic " "version string"
#         )
# else:
#     # have git, in git dir, but may have used a shallow clone (travis)
#     rev = sout.strip()
#     rev = rev.decode("ascii")

# if p_fn_ver.exists():
#     with open(p_fn_ver, "r") as f:
#         dict_ver = json.load(f)
#         ver_main = dict_ver["version"]
#         ver_post = dict_ver["iter"]
#         ver_git_commit = dict_ver["git_commit"]

#     # print(dict_ver)
#     __version__ = "-".join(filter(None, [ver_main, ver_post, ver_git_commit]))
#     # raise ValueError(f"version info found: {__version__}")
# else:
#     __version__ = "0.0.0"
#     raise ValueError("version info not found")


# wrap OS-specific `SUEWS_driver` libs
sysname = platform.system()
lib_basename = "supy_driver"
if sysname == "Windows":
    lib_name = lib_basename + ".pyd"
elif sysname == "Darwin":
    lib_name = lib_basename + ".so"
elif sysname == "Linux":
    lib_name = lib_basename + ".so"

# change compiler settings
if sysname == "Windows":
    pfn = Path.cwd() / "setup.cfg"
    try:
        shutil.copyfile("win-setup.cfg", pfn)
    except:
        pass

# load SUEWS Fortran source files
dir_f95 = "../suews/src"
path_src = Path(dir_f95)
path_mod = (path_src.parent / "mod").resolve()
path_lib = (path_src.parent / "lib").resolve()
path_target_f95 = [
    (path_src / f)
    for f in [
        "suews_ctrl_const.f95",
        "suews_ctrl_error.f95",
        "suews_util_meteo.f95",
        "suews_phys_waterdist.f95",
        "suews_phys_narp.f95",
        "suews_phys_atmmoiststab.f95",
        "suews_phys_resist.f95",
        "suews_phys_evap.f95",
        "suews_phys_snow.f95",
        "suews_phys_dailystate.f95",
        "suews_phys_lumps.f95",
        "suews_phys_anemsn.f95",
        "suews_phys_rslprof.f95",
        "suews_phys_biogenco2.f95",
        "suews_phys_ohm.f95",
        "suews_phys_estm.f95",
        "suews_ctrl_driver.f95",
    ]
]
# all_f95 = glob.glob(os.path.join(dir_f95, "*.f95"))
path_all_f95 = [f for f in path_src.glob("*.f95")]
path_exclude_f95 = [
    (path_src / f)
    for f in [
        "suews_c_wrapper.f95",
        "suews_ctrl_sumin.f95",
        "suews_program.f95",
        "suews_ctrl_init.f95",
        "suews_ctrl_calculations.f95",
        "suews_ctrl_translate.f95",
    ]
]
path_other_f95 = list(set(path_all_f95) - set(path_target_f95) - set(path_exclude_f95))
fn_other_obj = [str(f).replace(".f95", ".o") for f in path_other_f95]
if sysname == "Windows":
    fn_other_obj.append(os.path.join(dir_f95, "strptime.o"))

src_f95 = path_target_f95 + path_other_f95


def readme():
    try:
        with open("../README.md", encoding="utf-8") as f:
            return f.read()
    except:
        return f"SuPy package"

class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True

    def is_pure(self):
        return False

ext_modules = [
    Extension(
        "supy.supy_driver.suews_driver",
        [str(p) for p in path_target_f95],
        extra_compile_args=[
            "-D_POSIX_C_SOURCE=200809L",
            "-fbracket-depth=1024" if sysname == "Darwin" else '-Wall', # for clang on MacOS
        ],
        extra_f90_compile_args=["-cpp", f"-I{str(path_mod)}"],
        f2py_options=[
            # '--quiet',
            # "--verbose",
            # "--debug-capi",  # this is for debugging data types
            # '--f2cmap="f2py_f2cmap"',
            # ('-DF2PY_REPORT_ATEXIT' if sysname == 'Linux' else ''),
        ],
        extra_objects=fn_other_obj,
        # "-v" under Linux is necessary because it can avoid the blank variable issue
        # ref: https://github.com/metomi/fcm/issues/220
        extra_link_args=["-v" if sysname == "Linux" else "-static"]
        + [f"-L{str(path_lib)}", "-lspartacus"],
    )
]

setup(
    name="supy",
    # version=__version__,
    description="the SUEWS model that speaks python",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/UMEP-Dev/SuPy",
    author=", ".join(
        [
            "Dr Ting Sun",
            "Dr Hamidreza Omidvar",
            "Prof Sue Grimmond",
        ]
    ),
    author_email=", ".join(
        [
            "ting.sun@ucl.ac.uk",
            "h.omidvar@reading.ac.uk",
            "c.s.grimmond@reading.ac.uk",
        ]
    ),
    license="GPL-V3.0",
    packages=["supy"],
    package_data={
        "supy": [
            "sample_run/*",
            "sample_run/Input/*",
            "*.json",
            "util/*",
            "cmd/*",
        ]
    },
    distclass=BinaryDistribution,
    ext_modules=ext_modules,
    install_requires=[
        "pandas< 1.5; python_version <= '3.9'", # to fix scipy dependency issue in UMEP under QGIS3 wtih python 3.9
        "pandas; python_version > '3.9'",
        "matplotlib",
        "chardet",
        "scipy",
        "dask",  # needs dask for parallel tasks
        "f90nml",  # utility for namelist files
        "seaborn",  # stat plotting
        "atmosp",  # my own `atmosp` module forked from `atmos-python`
        "cdsapi",  # ERA5 data
        "xarray",  # utility for high-dimensional datasets
        "multiprocess",  # a better multiprocessing library
        "click",  # cmd tool
        "lmfit",  # optimiser
        "numdifftools",  # required by `lmfit` for uncertainty estimation
        "pvlib",  # TMY-related solar radiation calculations
        "platypus-opt==1.0.4",  # a multi-objective optimiser
        # DRIVER_REQ,  # a separate f2py-based driver
    ],
    extras_require={
        "hdf": [
            "tables",  # for dumping in hdf5
        ]
    },
    entry_points={
        #   command line tools
        "console_scripts": [
            "suews-run=supy.cmd.SUEWS:SUEWS",
            "suews-convert=supy.cmd.table_converter:convert_table_cmd",
        ]
    },
    include_package_data=True,
    python_requires="~=3.7",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
    ],
    zip_safe=False,
)
