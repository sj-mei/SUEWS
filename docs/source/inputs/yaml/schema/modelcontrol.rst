Modelcontrol
============

**Parameters:**

.. option:: tstep <int>

   Time step in seconds for model calculations

   :Default: ``300``

.. option:: forcing_file <RefValue[str] | str>

   Path to meteorological forcing data file. The forcing file contains time-series meteorological data that drives SUEWS simulations. For detailed information about required variables, file format, and data preparation guidelines, see :ref:`met_input`.

   :Default: ``'forcing.txt'``

.. option:: kdownzen <RefValue[int] | int | NoneType>

   Use zenithal correction for downward shortwave radiation

   :Default: Not specified

.. option:: output_file <str>

   Path to model output file. SUEWS generates multiple output files with time-series results and diagnostic information. For detailed information about output file formats, variables, and interpretation, see :ref:`output_files`.

   :Default: ``'output.txt'``

.. option:: diagnose <int>

   Level of diagnostic output (0=none, 1=basic, 2=detailed)

   :Default: ``0``

.. option:: start_time <str (Optional)>

   Start time of model run. If None use forcing data bounds.

   :Default: Not specified

.. option:: end_time <str (Optional)>

   End time of model run. If None use forcing data bounds.

   :Default: Not specified

.. option:: ref <Reference (Optional)>

   :Default: Not specified

   For ``ref``, if using the Reference structure, see :doc:`reference` for details.
