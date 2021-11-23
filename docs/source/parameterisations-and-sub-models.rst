.. _physics_schemes:

Parameterisations and sub-models within SUEWS
=============================================

Net all-wave radiation, Q\*
---------------------------

There are several options for modelling or using observed radiation
components depending on the data available. As a minimum, SUEWS requires
incoming shortwave radiation to be provided.

#. Observed net all-wave radiation can be provided as input instead of
   being calculated by the model.
#. Observed incoming shortwave and incoming longwave components can be
   provided as input, instead of incoming longwave being calculated by
   the model.
#. Other data can be provided as input, such as cloud fraction (see
   options in `RunControl.nml`).
#. **NARP** (Net All-wave Radiation Parameterization) :cite:`O03,L11` scheme calculates outgoing
   shortwave and incoming and outgoing longwave radiation components
   based on incoming shortwave radiation, temperature, relative humidity
   and surface characteristics (albedo, emissivity).



Anthropogenic heat flux, Q\ :sub:`F`
------------------------------------


#. Two simple anthropogenic heat flux sub-models exist within SUEWS:

   -  :cite:t:`J11` approach, based on heating and cooling degree days and population density (allows distinction between weekdays and weekends).
   -  :cite:t:`L11` approach, based on a linear piece-wise relation with air temperature.

#. Pre-calculated values can be supplied with the meteorological forcing data, either derived from knowledge of the study site, or obtained from other models, for example:

   -  **LUCY** :cite:`A11,L13`. A new version has been now included in UMEP. To distinguish it is referred to as `LQF`_
   -  **GreaterQF** :cite:`I11`. A new version has been now included in UMEP. To distinguish it is referred to as `GQF`_

Storage heat flux, ΔQ\ :sub:`S`
-------------------------------

#. Three sub-models are available to estimate the storage heat flux:

   -  **OHM** (Objective Hysteresis Model) :cite:`G91,GO99,GO02`. Storage heat heat flux is calculated using empirically-fitted relations with net all-wave radiation and the rate of change in net all-wave radiation.
   -  **AnOHM** (Analytical Objective Hysteresis Model) :cite:`S17`. OHM approach using analytically-derived coefficients. |NotRecmd|
   -  **ESTM** (Element Surface Temperature Method) :cite:`O05`. Heat transfer through urban facets (roof, wall, road, interior) is calculated from surface temperature measurements and knowledge of material properties. |NotRecmd|

#. Alternatively, 'observed' storage heat flux can be supplied with the meteorological forcing data.

Turbulent heat fluxes, Q\ :sub:`H` and Q\ :sub:`E`
--------------------------------------------------

#. **LUMPS** (Local-scale Urban Meteorological Parameterization Scheme) :cite:`GO02` provides a simple means of estimating sensible and latent heat fluxes based on the proportion of vegetation in the study area.

#. **SUEWS** adopts a more biophysical approach to calculate the latent heat flux; the sensible heat flux is then calculated as the residual of the energy balance.
   The initial estimate of stability is based on the LUMPS calculations of sensible and latent heat flux.
   Future versions will have alternative sensible heat and storage heat flux options.

Sensible and latent heat fluxes from both LUMPS and SUEWS are provided in the `output_files`.
Whether the turbulent heat fluxes are calculated using LUMPS or SUEWS can have a major impact on the results.
For SUEWS, an appropriate surface conductance parameterisation is also critical :cite:`J11` :cite:`W16`.
For more details see `Differences_between_SUEWS_LUMPS_and_FRAISE` .

Water balance
-------------

The running water balance at each time step is based on the urban water balance model of :cite:t:`G86` and urban evaporation-interception scheme of :cite:t:`GO91`.

-  Precipitation is a required variable in the meteorological forcing file.
-  Irrigation can be modelled :cite:`J11` or observed values can be provided if data are available.
-  Drainage equations and coefficients to use must be specified in the input files.
-  Soil moisture can be calculated by the model.
-  Runoff is permitted:

   -  between surface types within each model grid
   -  between model grids (|NotAvail|)
   -  to deep soil
   -  to pipes.

Snowmelt
--------

The snowmelt model is described in :cite:t:`J14`.
Changes since v2016a:
1) previously all surface states could freeze in 1-h time step, now the freezing surface state is
calculated similarly as melt water and can freeze within the snow pack.
2) Snowmelt-related coefficients have also slightly changed (see `SUEWS_Snow.txt`).

Convective boundary layer
-------------------------

A convective boundary layer (CBL) slab model :cite:`CG01` calculates the CBL height, temperature and humidity during daytime :cite:`O15`.

.. SOLWEIG is fully removed since 2019a

.. Thermal comfort
.. ---------------

.. **SOLWEIG** (Solar and longwave environmental irradiance geometry model,
.. Lindberg et al. 2008 :cite:`F08`, Lindberg and Grimmond 2011 :cite:`FG11`) is a 2D
.. radiation model to estimate mean radiant temperature.

.. .. figure:: /assets/img/Bluews_2.jpg
..     :alt:  Overview of scales. Source: Onomura et al. (2015) :cite:`O15`

..     Overview of scales. Source: Onomura et al. (2015) :cite:`O15`




.. _LQF: http://umep-docs.readthedocs.io/en/latest/OtherManuals/LQF_Manual.html
.. _GQF: http://umep-docs.readthedocs.io/en/latest/OtherManuals/GQF_Manual.html

.. _rsl_mod:

Wind, Temperature and Humidity Profiles in the Roughness Sublayer
----------------------------------------------------------------------------
A dignostic RSL scheme for calculating the wind, temperature and humidity profiles in the roughness sublayer is implemented in 2020a following :cite:t:`HF07, HF08` and :cite:t:`T19`.
An recent application of this RSL scheme can be found in :cite:t:`T21`.

The diagnostic profiles are outputed in 30 uneven levels between the ground and forcing height, which are divided into two groups:

- One group of levels are evenly distributed within the urban canopy layer characterised by mean height of roughness elements (e.g. buildings, trees, etc.) :math:`z_H`, which determines the number of layers within urban canopy :math:`n_{can}`:

.. math::
   :nowrap:

   \[
         n_{can} =
   \begin{cases}
      3 & \text{if } z_H \leq \text{2 m} \\
      10 & \text{if } \text{2 m} \lt z_H \leq \text{10 m} \\
      15 & \text{if } z_H \gt \text{10 m} \\

   \end{cases}
   \]

- The other levels are evenly distributed between the urban canopy layer top and forcing height.


.. note::

   All the diagnostic profiles (wind speed, temperature and humidity) are calculated
   from the forcing data down into the canopy.
   Therefore it is assumed that the forcing temperature and humidity
   are above the blending height.



Common near-surface diagnostics:

   -  T2: air temperature at 2 m agl
   -  Q2: air specific humidity at 2 m agl
   -  RH2: air relative humidity at 2 m agl
   -  U10: wind speed at 10 m agl

are calculated by the `RSL scheme <rsl_mod>` by interpolating RSL profile results to the corresonding diagnostic heights.
