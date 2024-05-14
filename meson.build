project('supy', 'fortran', 'c',)


# find the python installation
py = import('python').find_installation(pure: false)
dep_py = py.dependency()


subdir('src/supy_driver')

subdir('src/supy')

