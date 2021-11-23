.. _index_page:

SUEWS: Surface Urban Energy and Water Balance Scheme
----------------------------------------------------

.. image:: http://readthedocs.org/projects/suews/badge/?version=latest
    :target: https://suews.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


How to get SUEWS?
=================

- **Latest release:**

  The **latest formal** release of SUEWS is `new_latest` and can be downloaded via `our Zenodo repository`_ (a sample input dataset is included in the release archive).

- **Previous releases:**

  Previous releases can be downloaded via `our GitHub page`_.



How to use SUEWS?
=================

- **For existing users:**

  Overview of changes in this version, see :ref:`new_latest`.
  If these changes impact your existing simulations, please see appropriate parts of the manual. It may be necessary to adapt some of your input files for for the current version.

  .. tip::

      A helper python script, `SUEWS table converter <input_converter>`, is provided to help facilitate the conversion of input files between different SUEWS versions.

  Additionally, the manuals for previous versions can be accessed in respective sections under `version_history`.


- **For new users:**


  Before performing SUEWS simulations, new users should read the overview :ref:`introduction`, then follow the steps in `Preparing_to_run_the_model` to prepare `input files <input_files>` for SUEWS.

  Note there are tutorials learning about running SUEWS available :ref:`the tutorial. <tutorials_suews>`


How to get help in using SUEWS?
==================================

Please let us know in the `UMEP Community`_.
The developers and other users are willing to help you :)


How has SUEWS been used?
==================================

The scientific details and application examples of SUEWS can be found in `Recent_publications`.

.. _cite_suews:


How to cite SUEWS?
==================================

Please go to `our Zenodo repository`_ for a proper citation of SUEWS.

.. tip::

    Visit the repositories below for different citation styles.




How to support SUEWS?
==================================

#. `Cite SUEWS <cite_suews>` appropriately in your work.
#. Contribute to the `development <Development_Suggestions_Support>`.
#. Report issues via the `GitHub page <https://github.com/UMEP-dev/SUEWS/issues/new/choose>`_.
#. Provide `suggestions and feedback <Development_Suggestions_Support>`.


.. _our GitHub page: https://github.com/UMEP-dev/SUEWS
.. _our Zenodo repository: `Zenodo page`_
.. _this form: `dowload form`_
.. _dowload form: https://github.com/UMEP-dev/SUEWS/releases

.. |doi_software| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3533450.svg
    :target: http://doi.org/10.5281/zenodo.3533450


.. toctree::
   :caption: For users
   :maxdepth: 2
   :numbered:
   :hidden:

   introduction
   prepare-to-run-the-model
   input_files/input_files
   output_files/output_files
   troubleshooting
   related_softwares
   sub-tutorials/tutorials
   benchmark/benchmark_report
   notation


.. toctree::
   :caption: For developers/contributors
   :maxdepth: 2
   :numbered:
   :hidden:

   contributing/contributing
   api
   acknowledgement
   version-history/version-history
   GitHub repository <https://github.com/UMEP-dev/SUEWS>


.. toctree::
   :caption: Scientific background
   :maxdepth: 2
   :numbered:
   :hidden:

   parameterisations-and-sub-models
   related_publications
   references



.. toctree::
   :caption: Community
   :maxdepth: 2
   :numbered:
   :hidden:

   GitHub discussion <https://github.com/UMEP-dev/UMEP/discussions>
   SUEWS mailing list <https://www.lists.rdg.ac.uk/mailman/listinfo/met-suews>

