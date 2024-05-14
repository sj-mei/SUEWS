
print("Hello, world!")
print("This is a Python script to move generated files to proper location for meson build.")

from pathlib import Path
import os
import sys

# receive two arguments from the command line
OUTPUT = sys.argv[1:-1]
OUTDIR = sys.argv[-1]


list_fn_out = [Path(f).name for f in OUTPUT]
print("output files:", list_fn_out)

list_p_out = [Path.cwd()/f for f in list_fn_out]
for f in list_p_out:
  if not f.exists():
    raise FileNotFoundError(f)
  else:
    print(f, "exists")

p_outdir = Path(OUTDIR)
print("Output directory:", p_outdir)

# # mv the generated files to the source directory
for f in list_p_out:
  f.rename(p_outdir / f.name)