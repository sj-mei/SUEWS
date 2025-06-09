project(
    'supy',
    'fortran',
    'c',
    version: run_command('get_ver_git.py',check: true).stdout().strip(),
    meson_version: '>=0.64.0',
)

if host_machine.system() == 'windows' and host_machine.cpu_family() == 'x86'
  add_project_arguments('-D__MINGW_USE_VC_SETJMP=1', language: ['c', 'cpp'])
endif

# find the python installation
py = import('python').find_installation(pure: false)
dep_py = py.dependency()

subdir('src/supy_driver')

subdir('src/supy')
