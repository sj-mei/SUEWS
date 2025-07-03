.. _config_builder_guide:

SUEWS Configuration Builder User Guide
======================================

The SUEWS Configuration Builder is an interactive web-based tool that helps you create valid YAML configuration files for SUEWS simulations. This guide explains how to use the builder effectively.

.. warning::

   **Beta Status**: The configuration builder is currently in beta development. While it has been tested with common configurations, please verify generated files carefully before use in production simulations.

Accessing the Builder
---------------------

The configuration builder is available in two ways:

1. **Online**: `Access the builder <../../_static/index.html>`_ directly from the documentation
2. **Local**: Open ``docs/source/_static/index.html`` in your web browser from your SUEWS installation

Getting Started
---------------

Basic Workflow
~~~~~~~~~~~~~~

1. **Start with Basic Mode**: The builder opens in Basic Mode by default, showing only essential parameters
2. **Navigate sections**: Use the sidebar to navigate between configuration sections
3. **Fill in parameters**: Enter values for required parameters (marked with asterisks)
4. **Preview YAML**: See your configuration in real-time in the preview panel
5. **Export**: Save your configuration as a YAML file

Interface Overview
~~~~~~~~~~~~~~~~~~

The builder interface consists of:

- **Header**: Contains action buttons (New, Import, Export, Validate)
- **Sidebar**: Navigation between configuration sections
- **Main Form**: Parameter input fields organised by category
- **Preview Panel**: Real-time YAML preview of your configuration

Key Features
------------

Basic vs Advanced Mode
~~~~~~~~~~~~~~~~~~~~~~

- **Basic Mode** (default): Shows only commonly-used parameters
- **Advanced Mode**: Reveals all parameters including advanced physics options

Toggle between modes using the switch in the control bar. Your entered values are preserved when switching modes.

Parameter Search
~~~~~~~~~~~~~~~~

Use the search box to quickly find parameters:

1. Type any part of a parameter name or description
2. The interface highlights matching parameters
3. Non-matching sections are temporarily hidden
4. Clear the search to show all parameters again

Real-time Validation
~~~~~~~~~~~~~~~~~~~~

The builder validates your input as you type:

- **Red borders**: Invalid values (e.g., text in numeric fields)
- **Required markers**: Asterisks (*) indicate mandatory parameters
- **Tooltips**: Hover over field labels for parameter descriptions

Working with Arrays
~~~~~~~~~~~~~~~~~~~

The builder handles various array types:

**Simple Arrays** (e.g., surface fractions):
   - Click "Add Item" to add elements
   - Click "Remove" to delete elements
   - Enter values directly in the input fields

**Vertical Layer Arrays** (e.g., building_frac, veg_frac):
   - Automatically synchronised with the ``nlayer`` parameter
   - Cannot manually add/remove items
   - All layers displayed as inline inputs

**Complex Arrays** (e.g., sites):
   - Each item expands to show nested parameters
   - Use "Add Site" to create new site configurations
   - Collapse/expand items for easier navigation

Import and Export
-----------------

Importing Configurations
~~~~~~~~~~~~~~~~~~~~~~~~

1. Click the **Import** button
2. Paste your existing YAML configuration
3. Click **Import Configuration**
4. The builder populates all recognised parameters

.. note::

   The import feature preserves all valid parameters from your YAML, even those not visible in Basic Mode.

Exporting Configurations
~~~~~~~~~~~~~~~~~~~~~~~~

1. Click **Export YAML** button
2. The generated YAML opens in a modal
3. Click **Copy to Clipboard** or **Download**
4. Save with a ``.yml`` extension

Common Tasks
------------

Creating a Minimal Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Start in Basic Mode
2. Fill in **Model** section:
   
   - Set timestep (e.g., 3600 for hourly)
   - Specify forcing file path
   - Set output file path

3. Add at least one **Site**:
   
   - Click "Add Site" 
   - Enter site name
   - Set latitude and longitude
   - Define surface fractions (must sum to 1.0)

4. Export your configuration

Adding Multiple Sites
~~~~~~~~~~~~~~~~~~~~~

1. Navigate to the **Sites** section
2. Click **Add Site** for each location
3. Configure each site independently
4. Sites can share model settings but have unique properties

Configuring Advanced Physics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Switch to **Advanced Mode**
2. Navigate to **Model → Physics** section
3. Select physics options:
   
   - Storage heat method
   - Radiation schemes
   - Stability methods
   
4. Additional parameters appear based on selections

Tips and Best Practices
-----------------------

Surface Fractions
~~~~~~~~~~~~~~~~~

- Must sum to exactly 1.0
- Use the preview to verify totals
- Common urban fractions:
  
  - Buildings: 0.3-0.5
  - Paved: 0.2-0.4
  - Vegetation: 0.1-0.3
  - Water: 0.0-0.1

Parameter Units
~~~~~~~~~~~~~~~

- Always check parameter units in tooltips
- Common units:
  
  - Heights: metres [m]
  - Temperature: Celsius [°C] or Kelvin [K]
  - Radiation: W/m²
  - Time: seconds [s]

Validation Errors
~~~~~~~~~~~~~~~~~

If validation fails:

1. Check the browser console for detailed errors
2. Ensure all required fields are filled
3. Verify surface fractions sum to 1.0
4. Check that file paths are properly formatted

Troubleshooting
---------------

Page Not Loading
~~~~~~~~~~~~~~~~

- Ensure JavaScript is enabled in your browser
- Try a different browser (Chrome, Firefox, Safari recommended)
- Check browser console for error messages

Import Not Working
~~~~~~~~~~~~~~~~~~

- Verify YAML syntax is valid
- Ensure parameter names match current schema
- Try importing smaller sections individually

Array Display Issues
~~~~~~~~~~~~~~~~~~~~

- Vertical layer arrays auto-sync with ``nlayer``
- Clear browser cache if display seems corrupted
- Refresh page to reset the interface

Known Limitations
-----------------

- Some complex nested structures may require manual editing
- Browser-specific display differences may occur
- Large configurations may impact performance
- Not all SUEWS parameters are exposed in the interface

Getting Help
------------

- Check parameter descriptions via tooltips
- Refer to :ref:`yaml_input` for detailed parameter documentation
- Consult example configurations in the SUEWS repository
- Report issues on the `SUEWS GitHub page <https://github.com/UMEP-dev/SUEWS/issues>`_

Future Improvements
-------------------

Planned enhancements include:

- Configuration templates for common scenarios
- Batch site configuration import
- Enhanced validation with physics-based rules
- Configuration comparison tools
- Direct integration with SUEWS runner