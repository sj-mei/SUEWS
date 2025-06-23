.. _api_simulation:

SUEWSSimulation API
===================

The SUEWSSimulation class provides an object-oriented interface for running SUEWS simulations with a modern, intuitive workflow.

Overview
--------

The :class:`supy.SUEWSSimulation` class offers:

- YAML-based configuration management
- Chainable method design for clear workflows
- Built-in validation and error handling
- Multiple output formats (CSV, Excel, NetCDF, Pickle)
- Integration with existing SuPy infrastructure

Quick Start
-----------

.. code-block:: python

    from supy import SUEWSSimulation
    
    # Create simulation from YAML configuration
    sim = SUEWSSimulation.from_yaml('config.yml')
    
    # Setup forcing data and run
    sim.setup_forcing('forcing_data.txt')
    sim.run()
    
    # Access results
    results = sim.get_results()
    
    # Save outputs
    sim.save('outputs.csv')

Class Reference
---------------

.. currentmodule:: supy

.. autoclass:: SUEWSSimulation
    :members:
    :undoc-members:
    :show-inheritance:

    .. automethod:: __init__
    .. automethod:: from_yaml
    .. automethod:: setup_forcing
    .. automethod:: run
    .. automethod:: get_results
    .. automethod:: save
    .. automethod:: summary
    .. automethod:: see
    .. automethod:: quick_plot
    .. automethod:: clone
    .. automethod:: reset
    .. automethod:: validate

Examples
--------

Basic Simulation
~~~~~~~~~~~~~~~~

.. code-block:: python

    from supy import SUEWSSimulation
    from pathlib import Path
    
    # Load configuration
    sim = SUEWSSimulation.from_yaml('test/benchmark1/benchmark1.yml')
    
    # Setup forcing data
    forcing_file = Path('test/benchmark1/forcing/Kc1_2011_data_5.txt')
    sim.setup_forcing(forcing_file)
    
    # Run simulation
    sim.run()
    
    # Check results
    print(sim.summary())
    
    # Access specific variables
    results = sim.get_results()
    qh = results[('SUEWS', 'QH')]  # Sensible heat flux
    qe = results[('SUEWS', 'QE')]  # Latent heat flux

Parameter Overrides
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Override specific parameters at initialisation
    sim = SUEWSSimulation.from_yaml(
        'config.yml',
        tstep=300,  # 5-minute timestep
        debug_mode=True,
        output_dir='/path/to/outputs'
    )

Multiple Output Formats
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Save results in different formats
    sim.save('results.csv')           # CSV format
    sim.save('results.xlsx')          # Excel format
    sim.save('results.nc')            # NetCDF format
    sim.save('results.pkl')           # Pickle format
    
    # Custom save options
    sim.save('results.csv', 
             variables=['QH', 'QE', 'QS'],  # Specific variables
             resample='1h')                  # Hourly averages

Validation and Debugging
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Validate setup before running
    validation = sim.validate()
    if validation['status'] != 'valid':
        print("Warnings:", validation['warnings'])
        print("Errors:", validation['errors'])
    
    # Enable debug logging
    sim = SUEWSSimulation(config, debug_mode=True)
    
    # Check simulation state
    print(f"Simulation complete: {sim.is_complete}")
    print(f"Results available: {sim.has_results}")

Integration with Existing SuPy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import supy as sp
    
    # Traditional SuPy workflow
    df_state, df_forcing = sp.load_SampleData()
    df_output, df_state_final = sp.run_supy(df_forcing, df_state)
    
    # New object-oriented workflow
    sim = SUEWSSimulation(df_state)
    sim.setup_forcing(df_forcing)
    sim.run()
    
    # Results are compatible
    df_output_new = sim.get_results()

Data Structures
---------------

The SUEWSSimulation class works with standard SuPy data structures:

- **Configuration**: YAML files or SUEWSConfig objects from :mod:`supy.data_model`
- **Forcing Data**: Pandas DataFrame with datetime index and meteorological variables
- **Results**: Multi-level DataFrame with (grid, datetime) index and (group, variable) columns

See Also
--------

- :ref:`data_model` - Configuration data model reference
- :ref:`df_forcing` - Forcing data format specification
- :ref:`df_output` - Output data format specification
- :doc:`/inputs/yaml/index` - YAML configuration guide