from numpy.distutils.mingw32ccompiler import Mingw32CCompiler as mgwin32C

from numpy.distutils.mingw32ccompiler import find_python_dll,generate_def
import sys,os
from pathlib import Path
from time import sleep


print("Fixing def file for Windows")


dll_file = find_python_dll()
p_dll=Path(dll_file)
print('Here is the pythonlib.dll:', dll_file)

# generate symbol list from this library
def_name = "python%d%d.def" % tuple(sys.version_info[:2])
def_file = os.path.join(sys.prefix, 'libs', def_name)

# fix path for def file
Path(def_file).parent.mkdir(parents=True, exist_ok=True)
print('OK to create? ', Path(def_file).parent.exists())

# t_sleep_sec=1
# print(f'sleeping for {t_sleep_sec} seconds')
# sleep(t_sleep_sec)


generate_def(dll_file, def_file)