.. _yaml_input:

YAML Configuration Format
=========================

The YAML configuration format is the recommended method for providing inputs to SUEWS. It uses a single, structured `config_suews.yml` file to define all model parameters, making simulations easier to manage and reproduce.

Overview
--------

A SUEWS YAML configuration file is organized into two main sections:

1. :ref:`model <model>`: Contains global settings that control the simulation, such as physics options, time stepping, and file paths.
2. :ref:`sites <site>`: A list of one or more sites to be simulated. Each site has its own set of properties, initial conditions, and land cover characteristics.

Here's a minimal example of the YAML structure:

.. code-block:: yaml

   # SUEWS configuration file
   model:
     control:
       tstep: 3600
       forcing_file: "forcing.txt"
     physics:
       net_radiation_method: 3
   
   sites:
     - name: "London_KCL"
       latitude: 51.5115
       longitude: -0.1160
       land_cover:
         paved: 0.38
         bldgs: 0.37
         grass: 0.14

For a complete working example, please refer to the `sample configuration file <https://github.com/UMEP-dev/SUEWS/blob/master/src/supy/sample_run/sample_config.yml>`_ provided with SuPy.


Data Files
----------

In addition to the YAML configuration file, SUEWS works with input and output data files:

**Input Data:**
- **Forcing data**: Meteorological time-series data file specified by :ref:`model.control.forcing_file <modelcontrol>`

**Output Data:**
- **Model results**: Time-series output files specified by :ref:`model.control.output_file <modelcontrol>`

For detailed information about:

- **Input data format and variables**: see :ref:`met_input`
- **Output file formats and variables**: see :ref:`output_files`

Configuration Builder (Experimental)
--------------------------------------

.. note::

   **Interactive Configuration Builder**: An experimental web-based tool for creating YAML configuration files with guided forms.

   📝 **Access the builder**: `SUEWS Configuration Builder <../../_static/index.html>`__

   .. warning::

      **Beta Status**: This configuration builder is currently in development and has not been fully tested. Please verify generated configurations carefully before use.

   Features:

   - Interactive forms for parameter input
   - Basic validation for parameter values
   - YAML file export functionality
   - Integrated parameter documentation

   **Alternative**: You can also create YAML files manually using the schema documentation below or by adapting the sample configuration files provided with SuPy.

Data Model Schema
-----------------

The following pages provide a detailed reference for every component of the YAML data model. Each page corresponds to a specific parameter group and details its available keys, expected data types, units, and default values.

**Key Features of the Schema Documentation:**

- **Clear Method Descriptions**: All model physics methods now include explicit explanations of what each numeric option represents (e.g., 0=MOST, 1=RST, 2=VARIABLE for diagnostic methods)
- **Comprehensive Parameter Documentation**: Every parameter includes its purpose, units, default values, and constraints
- **Hierarchical Organization**: The schema is organized hierarchically, with top-level components (model, site) linking to detailed sub-components
- **Cross-References**: Related parameters and methods are cross-referenced for easier navigation

.. toctree::
   :maxdepth: 2
   :caption: Schema Reference

   schema/model
   schema/site