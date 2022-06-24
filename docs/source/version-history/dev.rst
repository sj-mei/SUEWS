
.. warning:: Information here is ONLY for developers.


Version 2021a (in development)
----------------------------------------------------

- **Improvement**

  1. Added a new `RoughLenMomMethod` (``4``) to calculate roughness and displacement height as a function of plan area index and effective height of roughness elements following the ensemble mean of Fig 1a in :cite:`GO99`
  2. Coupled `SPARCATUS <https://github.com/Urban-Meteorology-Reading/spartacus-surface>`_ into SUEWS for detailed modelling of radiation balance.
  3. Added a new option `DiagMethod` in `RunControl` to control the output of radiation balance.


- **Changes**

  1. TO ADD


- **Fix**

  #. fixed a bug in radiation scheme: observed incoming longwave radiation cannot be used.

- **Known issues**

  #. Wind direction is not currently downscaled so non -999 values will cause an error.
