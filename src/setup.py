import os
from signal import raise_signal
from setuptools import setup
import json
from pathlib import Path

# write version info using git commit
import subprocess
import warnings
import re

ISRELEASED = True
# if a release, use strict requirement for supy-driver; otehrwise, use a loose requirement
DRIVER_REQ = "supy_driver==2021a15" if ISRELEASED else "supy_driver"
# FULLVERSION += '.dev'

pipe = None
p_fn_ver = Path("./supy/supy_version.json")

# force remove the version info file
flag_dirty = False

for cmd in ["git", "/usr/bin/git", "git.cmd"]:

    try:
        pipe = subprocess.Popen(
            [cmd, "describe", "--always", "--match", "2[0-9]*", "--dirty=-dirty"],
            stdout=subprocess.PIPE,
        )
        (sout, serr) = pipe.communicate()
        # parse version info from git
        list_str_ver = sout.decode("utf-8").strip().split("-")

        if list_str_ver[-1].lower() == "dirty":
            flag_dirty = True
            # remove the "dirty" part from version info list
            list_str_ver = list_str_ver[:-1]

        ver_main = list_str_ver[0]
        print("ver_main", ver_main)
        if len(list_str_ver) > 1:
            ver_post = list_str_ver[1]
            ver_git_commit = list_str_ver[2]
        else:
            ver_post = ""
            ver_git_commit = ""

        # save version info to json file
        p_fn_ver.unlink(missing_ok=True)
        with open(p_fn_ver, "w") as f:
            json.dump(
                {
                    "version": ver_main + ("" if ISRELEASED else ".dev"),
                    "iter": ver_post,
                    "git_commit": ver_git_commit + ("-dirty" if flag_dirty else ""),
                },
                f,
            )
        if pipe.returncode == 0:
            print(f"in {cmd}, git version info saved to", p_fn_ver)
            break
    except Exception as e:
        pass

if pipe is None or pipe.returncode != 0:
    # no git, or not in git dir

    if p_fn_ver.exists():
        warnings.warn(
            f"WARNING: Couldn't get git revision, using existing {p_fn_ver.as_posix()}"
        )
        write_version = False
    else:
        warnings.warn(
            "WARNING: Couldn't get git revision, using generic " "version string"
        )
else:
    # have git, in git dir, but may have used a shallow clone (travis)
    rev = sout.strip()
    rev = rev.decode("ascii")

if p_fn_ver.exists():
    with open(p_fn_ver, "r") as f:
        dict_ver = json.load(f)
        ver_main = dict_ver["version"]
        ver_post = dict_ver["iter"]
        ver_git_commit = dict_ver["git_commit"]

    # print(dict_ver)
    __version__ = f"{ver_main}-{ver_post}-{ver_git_commit}".strip()
    # raise ValueError(f"version info found: {__version__}")
else:
    __version__ = "0.0.0"
    raise ValueError("version info not found")


def readme():
    try:
        with open("../README.md", encoding="utf-8") as f:
            return f.read()
    except:
        return f"SuPy package"


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
    # distclass=BinaryDistribution,
    ext_modules=[],
    install_requires=[
        "pandas>=1.3",
        "matplotlib",
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
        DRIVER_REQ,  # a separate f2py-based driver
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
