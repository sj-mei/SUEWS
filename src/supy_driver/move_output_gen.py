print("""
===============================================================================
This is a Python script to move generated files to proper location for meson build.
It is called by the meson build system.

""")

from pathlib import Path
import os
import sys

# receive two arguments from the command line
OUTPUT = sys.argv[1:-1]
OUTDIR = sys.argv[-1]


list_fn_out = [Path(f).name for f in OUTPUT]
print("output files:", list_fn_out)


list_p_out = [Path.cwd() / f for f in list_fn_out]
for f in list_p_out:
    if not f.exists():
        # raise FileNotFoundError(f)
        # generate an empty file
        f.touch()
        print(f, "does not exist; but I've created an empty file.")
    else:
        print(f, "exists")

p_outdir = Path(OUTDIR)
print("Output directory:", p_outdir)
# mv the generated files to the source directory
for f in list_p_out:
    f.rename(p_outdir / f.name)
