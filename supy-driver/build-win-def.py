from numpy.distutils.mingw32ccompiler import Mingw32CCompiler as mgwin32C

from numpy.distutils.mingw32ccompiler import find_python_dll,generate_def
import sys,os
from pathlib import Path
from time import sleep


print("Fixing def file for Windows")


# mgwin32=mgwin32C()
# mgwin32.build_import_library()

dll_file = find_python_dll()
p_dll=Path(dll_file)
print('here is the dll:', dll_file)
# log.info('Building import library (arch=AMD64): "%s" (from %s)' %
#             (out_file, dll_file))

# generate symbol list from this library
def_name = "python%d%d.def" % tuple(sys.version_info[:2])
def_file = os.path.join(sys.prefix, 'libs', def_name)
Path(def_file).parent.mkdir(parents=True, exist_ok=True)
print('ok to create? ', Path(def_file).parent.exists())
# print('ok to create? ', Path(def_file).parent.exists())

# print('sleeping for 999 seconds')
# sleep(999)

generate_def(dll_file, def_file)