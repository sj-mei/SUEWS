# This is the `meson.build` file for building `supy_driver`, which relies on SUEWS - a Fortran-dominant project with some C code.

project('supy_driver', 'fortran', 'c')

fortran_compiler = meson.get_compiler('fortran')
if fortran_compiler.get_id() == 'gcc'
  f90flags = ['-fPIC', '-fconvert=big-endian']
  fortran_program = find_program('gfortran')
elif fortran_compiler.get_id() == 'intel'
  f90flags = ['-fpscomp logicals', '-fPIC']
  fortran_program = find_program('ifort')
endif

# flags for converting f95 to fpp
fpp_flags = ['-E', '-x', 'f95-cpp-input', '-fPIC']

add_global_arguments(f90flags, language: 'fortran')

# lib_suews = meson.get_compiler('fortran').find_library('suews', dirs: '../suews/lib')
include_dir_suews = include_directories('../suews/mod')

dir_src_suews = '../suews/src'

src_f95 = [
  'suews_ctrl_const.f95',
  'suews_util_stringmod.f95',
  'suews_ctrl_input.f95',
  'suews_ctrl_output.f95',
  'suews_util_time.f95',
  'suews_util_meteo.f95',
  'suews_phys_waterdist.f95',
  'suews_phys_narp.f95',
  'suews_phys_beers.f95',
  'suews_phys_atmmoiststab.f95',
  'suews_phys_resist.f95',
  'suews_phys_evap.f95',
  'suews_phys_snow.f95',
  'suews_phys_dailystate.f95',
  'suews_phys_lumps.f95',
  'suews_phys_anemsn.f95',
  'suews_phys_rslprof.f95',
  'suews_phys_biogenco2.f95',
  'suews_phys_ohm.f95',
  'suews_phys_estm.f95',
  'suews_phys_spartacus.f95',
  'suews_ctrl_driver.f95',
]

src_f95_files = []
foreach f : src_f95
  src_f95_files += join_paths(dir_src_suews, f)
endforeach

src_f95_files = files(src_f95_files)

src_to_wrap = [
  'suews_ctrl_driver',
  'suews_phys_narp',
  'suews_phys_atmmoiststab',
  'suews_phys_waterdist',
  'suews_phys_resist',
  'suews_phys_evap',
  'suews_phys_dailystate',
  'suews_phys_anemsn',
  'suews_phys_rslprof',
  'suews_util_meteo',
  'suews_phys_spartacus',
]

src_wrap_f95 = []
foreach f : src_to_wrap
  src_wrap_f95 += f + '.f95'
endforeach

src_wrap_objects = []
foreach f : src_to_wrap
  src_wrap_objects += f + '.o'
endforeach

src_wrap_fpp = []
foreach f : src_to_wrap
  src_wrap_fpp += f + '.fpp'
endforeach

# convert f95 to fpp files
f95_to_fpp = generator(
  fortran_program,
  output: '@BASENAME@.fpp',
  arguments: [fpp_flags, '@INPUT@', '-o', '@OUTPUT@'],
)

# Use the generator to process each source file
gen_fpp_files = f95_to_fpp.process(src_f95_files)

# # Compile the f95 files
# f95_to_o = fortran_compiler.compile(f95_to_fpp)

# invoke make in the suews directory to build the objects
custom_target(
  'suews_objects',
  command: ['make', '-C', '../../suews/', '&&', 'touch', 'dummy.flag'],
  build_by_default: true,
  output: ['dummy.flag'],
)

# build the static library for later use
static_library(
  'src',
  src_f95_files,
  include_directories: include_dir_suews,
)

# build the driver
prog_f90wrap = find_program('f90wrap')
file_kind_map = files('kind_map')
supy_driver = custom_target(
  'supy_driver',
  command: [
    prog_f90wrap,
    '-m', 'supy_driver',
    '@INPUT@',
    '-k', file_kind_map,
    '-v',
    '--skip', 'error_hint',
  ],
  input: gen_fpp_files,
  build_by_default: true,
  output: 'supy_driver',
)

# build the driver
prog_f2py_f90wrap = find_program('f2py-f90wrap')
py3 = find_program('python3')
# Create a new custom target that depends on 'supy_driver'
f90wrap_files_target = custom_target(
  'f90wrap_files',
  command: [
    py3,
    '-c', 'import glob; print("\\n".join(glob.glob("builddir/f90wrap*.f90")))',
  ],
  # depend_files: files_f90wrap, # If applicable, specify files to check for changes
  depends: [supy_driver], # Ensure this runs after supy_driver
  output: 'f90wrap_files_list.txt',
  build_by_default: true,
)
# files_f90wrap = run_command(
#     py3,
#     '-c', 'import glob; print("\\n".join(glob.glob("builddir/f90wrap*.f90")))',
#     check: true,
#   ).stdout().strip().split('\n')

# custom_target(
#   '_supy_driver',
#   command: [
#     prog_f2py_f90wrap,
#     '--f90flags=\"-cpp -I../suews/mod\"',
#     '-c',
#     '-m', '_supy_driver',
#     '-L.',
#     '-lsrc', '@INPUT@',
#     '-L../suews/lib',
#     '-lsuewsdriver',
#     '-lspartacus',

#   ],
#   input: files_f90wrap,
#   build_by_default: true,
#   output: '_supy_driver',
# )