# Meson Build Fix: NumPy Absolute Paths

## Issue Description

When running `make dev` with Meson 1.8.2, the build fails with:

```
ERROR: Tried to form an absolute path to a dir in the source tree.
```

This occurs at line 249 in `src/supy_driver/meson.build` when trying to use `include_directories()` with NumPy's include path.

## Root Cause

- Meson's `include_directories()` function doesn't accept absolute paths (long-standing restriction since ~2017)
- `numpy.get_include()` returns an absolute path to NumPy's headers
- This restriction is enforced more strictly on Linux/macOS than Windows

## Solution

Replace `include_directories()` with compiler arguments (`-I` flags) for external absolute paths:

### Changes Made

1. **Remove invalid include_directories calls**:
   ```meson
   # OLD (causes error):
   inc_np = include_directories(incdir_numpy)
   inc_f2py = include_directories(incdir_f2py)
   ```

2. **Create compiler arguments instead**:
   ```meson
   # NEW:
   numpy_inc_args = ['-I' + incdir_numpy, '-I' + incdir_f2py]
   ```

3. **Update static library to use c_args**:
   ```meson
   lib_fortranobject = static_library(
     '_fortranobject',
     fortranobject_c,
     dependencies: dep_py,
     c_args: numpy_inc_args,  # Changed from include_directories
   )
   ```

4. **Update dependency declaration**:
   ```meson
   dep_fortranobject = declare_dependency(
     link_with: lib_fortranobject,
     compile_args: numpy_inc_args,  # Changed from include_directories
   )
   ```

5. **Update extension module**:
   ```meson
   ext_supy_driver = py.extension_module(
     '_supy_driver',
     sources: [...],
     include_directories: [inc_suews_mod],  # Keep relative paths
     c_args: c_args_windows + numpy_inc_args,  # Add numpy args to c_args
     ...
   )
   ```

## Why This Works

- Using `compile_args` with `-I` flags is the standard pattern for external absolute paths
- This approach is used by major projects like SciPy, scikit-learn
- Meson properly passes these flags to the compiler
- Maintains compatibility across platforms

## Testing

After applying this fix:
- `make dev` completes successfully
- `make test` passes all tests
- Build works on macOS with Meson 1.8.2

## Future Considerations

1. **Native NumPy dependency support**: Work is ongoing to support `dependency('numpy')` in Meson, which would eliminate this workaround
2. **Cross-platform testing**: Ensure this fix works on Windows and Linux
3. **Version compatibility**: Test with different NumPy and Meson versions

## References

- [Meson issue #1535: Absolute-path requirement for include_directories() breaks existing projects](https://github.com/mesonbuild/meson/issues/1535)
- [Meson discussion #12244: Thoughts on include_directories absolute paths](https://github.com/mesonbuild/meson/discussions/12244)
- [SciPy issue #16312: Meson complains about an absolute include path](https://github.com/scipy/scipy/issues/16312)

## Creating a Separate PR

This fix should be submitted as a separate PR to the main branch:

```bash
# From main branch:
git checkout master
git pull origin master
git checkout -b fix/meson-numpy-absolute-paths

# Apply the changes to src/supy_driver/meson.build
# Then:
git add src/supy_driver/meson.build
git commit -m "fix: handle NumPy absolute paths in Meson build"
git push origin fix/meson-numpy-absolute-paths
```

PR title: `fix: handle NumPy absolute paths in Meson build`