# %%
from shutil import rmtree
import Test_SUEWS as ts
import f90nml
import os
import numpy as np
import unittest
from pathlib import Path
from tempfile import gettempdir, TemporaryDirectory
from shutil import copyfile, copytree
# %%
fn_nml = "BTS_config.nml"
####################################
# load basic configurations
nml = f90nml.read(fn_nml)
cfg_file = nml["file"]

# load path
path_dir_exe = Path(cfg_file["dir_exe"]).resolve()
path_dir_src = path_dir_exe.parent / "src"

# identify programme name from SUEWS source code
path_suews_const = path_dir_src / "suews_ctrl_const.f95"
str_const = path_suews_const.read_text()
for ln in str_const.split("\n"):
    if "progname" in ln:
        name_exe = ln.split("progname = ")[-1].replace("'", "")
        break

path_baserun = Path(cfg_file["dir_baserun"]).resolve() / name_exe.replace("SUEWS_V", "")
path_exe=path_baserun/name_exe
# %%
# make longterm met as input
fn_met='test_2004_data_60.txt'
fn_met_long='.'.join(['test_2004_data_60.txt','long'])
fn_met_short='.'.join(['test_2004_data_60.txt','short'])

Path(path_baserun/'Input'/fn_met).unlink()
(path_baserun/'Input'/fn_met).symlink_to(fn_met_long)

# %%
# save current path for later use
dir_sys = os.getcwd()

# change to path for simulation
os.chdir(path_exe.parent)

# run simulation
# suppress output info
os.system("./" + name_exe + " &>/dev/null")


# change back to previous path
os.chdir(dir_sys)

# link back to short-term met
Path(path_baserun/'Input'/fn_met).unlink()
(path_baserun/'Input'/fn_met).symlink_to(fn_met_short)


# %%
