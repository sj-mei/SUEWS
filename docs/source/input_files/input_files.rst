.. _input_files:

Input files
===========

.. note::

   Since 2025, SUEWS adopts a new input file format, which is a YAML file, please refer to :ref:`yaml_input` for more details. To accomdate larger compatibility of input files, we sill keep the following sections operational as of May 2025 but will deprecate the table-based input files in the future (end of 2025 estimated; content may still be kept for reference).




.. toctree::
   :maxdepth: 1
   :numbered:

   yaml_input
   RunControl/RunControl
   SUEWS_SiteInfo/SUEWS_SiteInfo
   Initial_Conditions/Initial_Conditions
   met_input
   CBL_input/CBL_input
   ESTM_input/ESTM_input
   SS_input/SS_input
   input_converter



SUEWS allows you to input a large number of parameters to describe the characteristics of your site.

.. warning::

   You should not assume that the example values provided in sample files are appropriate.

YAML-based input file format
----------------------------

For the new YAML-based input file format, please refer to :ref:`yaml_input`.




Conventional SUEWS input files (prior to 2025)
----------------------------------------------

Values marked with 'MD' are examples of recommended values (see the suggested references to help decide how appropriate these are for your site/model domain);
values marked with 'MU' need to be set (i.e. changed from the example) for your site/model domain.

