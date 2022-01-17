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

are calculated by the `RSL scheme <rsl_mod>` by interpolating RSL profile results to the corresponding diagnostic heights.


SPARTACUS-Surface (SS)
----------------------

.. warning:: This module is highly experimental and not yet fully tested: description here is not yet complete, either. Please refer to the original `SPARTACUS-Surface page <https://github.com/ecmwf/spartacus-surface>`_ for more details, which may differ from the coupled version in SUEWS described below due to possibly different implementations.


.. note:: Future Work

   -  New SUEWS input table containing SPARTACUS profiles

   -  Add check for consistency of SUEWS and SS surface fractions

   -  Include snow

Introduction to SS
******************

The `SPARTACUS-Surface module <https://github.com/ecmwf/spartacus-surface>`_ computes the 3D interaction of shortwave and longwave radiation with complex surface canopies, including vegetated and urban canopies (with or without vegetation).


.. _SPARTACUS-Surface:
.. figure:: /assets/img/SUEWS002.jpg
	:alt: Multi-layer structure of SS

	Multi-layer structure (horizontal dashed lines) used in SS to characterise differences in the canopy (Cyan building, Green – vegetation). Source: `SPARTACUS-Surface GH page`_

It uses a multi-layer description of the canopy (:numref:`SPARTACUS-Surface`), with a statistical description of the horizontal distribution of trees and buildings.
Assumptions include:

-  Trees are randomly distributed.

-  Wall-to-wall separation distances follow an exponential probability distribution.

-  From a statistical representation of separation distances one can determine the probabilities of light being intercepted by trees, walls and the ground.

In the tree canopy (i.e. between buildings) there are two or three regions (based on user choice) (:numref:`schematic_tree_canopy`): clear-air and either one vegetated region or two vegetated regions of equal fractional cover but different extinction coefficient.
Assumptions include:

-  The rate of exchange of radiation between the clear and vegetated parts of a layer are assumed to be proportional to the length of the interface between them.

-  Likewise for the rate of interception of radiation by building walls.


.. _schematic_tree_canopy:
.. figure:: /assets/img/SUEWS003.jpg
   :alt: Areas between trees

   Areas between trees. Source: `SPARTACUS-Surface GH page`_

.. _SPARTACUS-Surface GH page: https://github.com/ecmwf/spartacus-surface


Each time light is intercepted it can undergo diffuse or specular reflection, be absorbed or be transmitted (as diffuse radiation).
The probabilities for buildings and the ground are determined by albedos and emissivities, and for trees are determined by extinction coefficients and single scattering albedos.

SUEWS-SS Implementation
************************

-  Maximum of 15 vertical layers

-  Building and tree fractions, building and tree dimensions, building albedo and emissivity, and diffuse versus specular reflection, can be treated as vertically heterogenous or uniform with height depending on parameter choices

-  As tree fraction increases towards 1 it is assumed that the tree crown merges when calculating tree perimeters.

-  Representing horizontal heterogeneity in the tree crowns is optional. When represented it is assumed that heterogeneity in leaf area index is between the core and periphery of the tree, not between trees.

-  When calculating building perimeters it is assumed that buildings do not touch (analogous to crown shyness) as building fraction increases towards 1.

-  Vegetation extinction coefficients (calculated from leaf area index, LAI) are assumed to be the same in all vegetated layers.

.. margin::

  .. [#estm_coupling] Confirming the ESTM coupling will allow this to be modified



-  Building facet and ground temperatures are equal to SUEWS TSfc_C (i.e.surface temperature) [#estm_coupling]_


.. margin::

  .. [#rsl_layers] This is from the RSL model (XX) and varies with layer

-  Leaf temperatures are equal to SUEWS temp_C (i.e. air temperature within the canopy) [#rsl_layers]_


-  Ground albedo and emissivity are an area weighted average of SUEWS paved, grass, bare soil and water values

-  Inputs from SUEWS: sfr, zenith_deg, TSfc_C, avKdn, ldown, temp_c, alb_next, emis, LAI_id

-  SS specific input parameters: read in from SUEWS_SPARTACUS.nml

-  Outputs used by SUEWS: alb_spc, emis_spc, lw_emission_spc

-  Although the radiation is calculated in multiple vertical layers within SS it is only the upwelling top-of-canopy fluxes:

   -  alb_spc*avKdn,

      -  (1emis_spc)*ldown,

      -  lw_emission_spc)

   -  that are used by SUEWS.

.. margin::

  .. [#ss_output] this will be updated but requires other updates first as of December 2021


- Output variables (including multi-layer ones) are in SUEWS-SS output file `SSss_YYYY_SPARTACUS.txt`. [#ss_output]_



RSL and SS Canopy Representation Comparison
*******************************************


-  The RSL has 30 levels but when the average building height is <2 m, < 12 m and > 12 m there are 3, 10 and 15 evenly spaced layers in the canopy.
-  The remaining levels are evenly spaced up to the forcing level (:numref:`SUEWS-RSL`).
-  The buildings are assumed to be uniform height.


.. _SUEWS-RSL:
.. figure:: /assets/img/SUEWS004.png
   :alt: SUEWS-RSL

   SUEWS-RSL module assumes the RSL has 30 layers that are spread between the canopy and within the atmosphere above


A maximum of 15 layers are used by SS (:numref:`vertial_layers_SS-RSL`), with the top of the highest layer at the tallest building height.
The layer heights are user defined and there is no limit on maximum building height.
The buildings are allowed to vary in height.


.. _vertial_layers_SS:
.. figure:: /assets/img/SUEWS005.png
   :alt: Vertical layers used by SS

   Vertical layers used by SS

.. .. |SUEWS005|

How to use SUEWS-SS
*******************

To run SUEWS-SS the SS specific files that need to be modified are:

- `RunControl.nml`

- `SUEWS_SPARTACUS.nml`

Non-SS specific SUEWS input file parameters also need to have appropriate values.
For example, LAI, albedos and emissivities are used by SUEWS-SS as explained in `more_SS_details`.



Inputs
^^^^^^

RunControl.nml
###################

SS is used within SUEWS when RunControl.nml parameter:

-  NetRadiationMethod is greater than 1000.

   -  Sensible values are 1001, 1002 or 1003.

   -  The remainder after division by 1000 is 1, 2 or 3 then the Ldown method of NetRadiationMethod 1, 2 or 3 is used, respectively.

SUEWS_SPARTACUS.nml
###################

-  File used to specify the SS model options when coupled to SUEWS.

-  **Configuration / Default Value**

-  SUEWS \* - Layer 1 SS parameter should equal the SUEWS parameter specified

.. list-table::
   :widths: 14 14 14 14 14 14 14
   :header-rows: 0


   * - **Parameter**
     - * *Type**
     - **De fault**
     - **Use**
     - * *Descri ption**
     - **Com ments**
     - **SUEWS \**

   * - nlayers
     - Integer
     - 1
     - MU
     - Number of vertical layers in the canopy
     - Max value is 15
     -

   * - use_sw_direct_albedo
     - Boolean
     - False – default

       True
     - MD
     - Specify ground and roof albedos sep arately for direct solar ra diation
     - Only con sidered when sw_dn_ direct_ frac>0.

       False reco mmended
       - unless user knows the diffuse and direct albedo.

       False (d efault)
       - Same ground albedo for diffuse and direct solar ra diation

       True - Specify ground albedo sep arately for direct solar ra diation
     -

   * - n_vegetation_region_urban
     - Integer
     - 1.
       -

       default

       2
     - MD
     - Number of regions used to describe vegetation
     - 1 might be okay de pending on the level of a ccuracy needed.

       See Hogan et al.
       (2018) – details of SPARTA CUS-Veg etation for more infor mation.

       1 (d efault) – hetero geneity of veg etation not con sidered

       2 – hetero geneity of veg etation con sidered
     -

   * - n_stream_sw_urban
     - Integer
     - 4
     - MD
     - SW diffuse streams per hem isphere (note: this is the number of qua drature points so a value of 4 corr esponds to an ‘8 -stream s cheme’)
     - 4 is recom mended.
       At large comput ational cost small impro vements in a ccuracy can be made by inc reasing from 4 (Hogan 2019b).

       Comp arisons with DART (S tretton et al.
       1)    found that 8 may be better (i.e.
       16 streams total)
     -

   * - n_stream_lw_urban
     - Integer
     - 4
     - MD
     - LW streams per hem isphere (note: this is the number of qua drature points so a value of 4 corr esponds to an ‘8 -stream s cheme’)
     - 4 is recom mended.
       At large comput ational cost small impro vements in a ccuracy can be made by inc reasing from 4 (Hogan 2019b).

       Eva luation against DART has used N = 8 (i.e 16 streams total).
     -

   * - sw_dn_direct_frac
     - Float
     - 0.45
     - MD
     - F raction of down welling sh ortwave ra diation that is direct
     - 0.45 is based on Berri zbeitia et al.
       (2020) ( Belgium and Berlin annual av erage), but could be im proved.
     -

   * - air_ext_sw
     - Float
     - 0.0
     - MD
     - Sh ortwave w aveleng th-inde pendent air ext inction coef ficient (m\ :su p:`-1`) (i.e.
       number of r adiance e-f oldings per metre)
     - (Hogan 2019b).

       Rea sonable approx imation (p ersonal commun ication Robin Hogan).
     -

   * - air_ssa_sw
     - Float
     - 0.95
     - MD
     - Sh ortwave single sca ttering albedo of air
     - 0.95 (Hogan 2019b)

       Air _ext_sw is not used if a ir_ext_ sw=0.0.
     -

   * - veg_ssa_sw
     - Float
     - 0.46
     - MD
     - Sh ortwave single sca ttering albedo of leaves
     - Br oadband sh ortwave veg etation SSA values range between 0.41 and 0.52 for RAMI5 Jä rvselja birch stand forest trees (see section 5.3).
       0.46 is the default value but users can choose their own value, for example by using the most appr opriate RAMI5 Jä rvselja birch stand forest tree
       type in section 5.3.
     -

   * - air_ext_lw
     - Float
     - 0.
     - MD
     - L ongwave w aveleng th-inde pendent air ext inction coef ficient (m\ :su p:`-1`) (i.e.
       number of r adiance e-f oldings per metre)
     - is a bad approx imation (Hogan 2019b)

       better represe ntation r equires several band tr eatment which is not in SS yet
     -

   * - air_ssa_lw
     - Float
     - 0.
     - MD
     - L ongwave single sca ttering albedo of air
     - is from Hogan 2019b.

       air _ssa_lw is not used when a ir_ext_ lw=0.0.
     -

   * - veg_ssa_lw
     - Float
     - 0.06
     - MD
     - L ongwave single sca ttering albedo of leaves
     - Should be es timated using a veg etation type in h ttps:// speclib .jpl.na sa.gov/ library (see section 5.3 for de tails).

       Refl ectance values are:

       ~0.04 for *Acer Pensylv anicum*

       ~0.02 for
       * Quercus Robur*

       ~0.04 for *Betula Lenta*.

       SSA ~2*refl ectance so 0.06 is chosen as the d efault.
     -

   * - veg_fsd
     - Float
     - 0.7
     - MD
     - Fra ctional s tandard de viation of the veg etation exti nction.
     - 0.7 has been used in SS for the RAMI-V radiat ion-veg etation inte rcompar ison\ : sup:`5` and the value should be updated based on the fi ndings.

       Det ermines ext inction coef ficient in the inner and outer layers of the tree crown when n_veget ation_r egion_u rban=2.
     - -

   * - veg_contact_fraction
     - Float
     - 0.
     - MD
     - F raction of veg etation edge in contact with b uilding walls
     - Change from 0.
       if d etailed canopy g eometry data are ava ilable.
     - -

   * - ground_albedo_dir_mult_fact
     - Float
     - 1.
     - MD
     - Ratio of the direct and diffuse albedo of the ground
     - Can from 1.
       if d etailed kn owledge of the direct and diffuse albedo is ava ilable.
     - -

   * - height
     - Float array (dim: nl ayer+1)
     - < 16
     - MU
     - Height of the layer int erfaces
     - -
     - -

   * - building_frac
     - Float array (dim: nlayer)
     - -
     - MU
     - B uilding plan area density
     - -
     - F r_Bldgs

   * - veg_frac
     - Float array (dim: nlayer)
     - -
     - MU
     - Tree plan area density
     - -
     - F r_EveTr + F r_DecTr

   * - building_scale
     - Float array (dim: nlayer)
     - -
     - MU
     - Building hor izontal scale (m)
     - Effective b uilding d iameter (See Fig.
       5, S tretton et al.
       (in prep))

       e stimate from ins pecting bu ildings (e.g. GIS data)

       used with buildi ng_frac to ca lculate the average b uilding pe rimeter length fo llowing Eq.
       8 of S partacu s_surfa ce_docu mentati on.pdf.
     - -

   * - veg_scale
     - Float array (dim: nlayer)
     - -
     - MU
     - Veg etation hor izontal scale (m)
     - Veg etation scale.

       e stimate from veg etation data (e.g.
       Google street view)

       used with veg_f raction to ca lculate the average veg etation pe rimeter length fo llowing Eq.
       2 of Hogan et al.
       (2018)
     - -

   * - roof_albedo
     - Float array (dim: nlayer)
     - -
     - MU
     - Roof albedo
     - If unknown
       - values can be found in SUE WS_NonV eg.txt.
     - -

   * - wall_albedo
     - Float array (dim: nlayer)
     - -
     - MU
     - Wall albedo
     - If unknown
       - values can be found in SUE WS_NonVeg.txt.
     - -

   * - roof_emissivity
     - Float array (dim: nlayer)
     - -
     - MU
     - Roof emi ssivity
     - If unknown
       - values can be found in SUE WS_NonV eg.txt.
     - -

   * - wall_emissivity
     - Float array (dim: nlayer)
     - -
     - MU
     - Wall emi ssivity
     - If unknown
       - values can be found in SUE WS_NonV eg.txt.
     - -

   * - roof_albedo_dir_mult_fact
     - Float array (dim: nlayer)
     - 1., 1., etc.
     - MD
     - Ratio of the direct and diffuse albedo of the roof
     - updated from 1.
       if know direct and diffuse albedo values
     - -

   * - wall_specular_frac
     - Float array (dim: nlayer)
     - 0., 0., etc
     - MD
     - F raction of wall ref lection that is s pecular
     - Updated from 0.
       if s pecular ref lection is known.
     - -



Outputs
^^^^^^^^^^^^

SSss_YYYY_SPARTACUS.txt
#######################

-  SSss - FileCode in RunControl.nml

-  YYYY - year the output file contains

-  For 12-71, 75-134 # represents each layer in nlayer

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 0


   * - Columns
     - Name
     - Description
     - SS Name

   * - 1
     - Year
     - Year [YYYY]
     - -

   * - 2
     - DOY
     - Day of year [DOY]
     - -

   * - 3
     - Hour
     - Hour [HH]
     - -

   * - 4
     - Min
     - Minute [MM]
     - -

   * - 5
     - Dectime
     - Decimal time [-]
     - -

   * - 6
     - alb
     - Albedo at top-of-canopy.

       Average of diffuse and direct albedos weighted by the amount of diffuse and direct shortwave radiation.
     - bc_out%sw_albedo and b c_out%sw_albedo_dir

   * - 7
     - emis
     - Emissivity at top-of-canopy
     - b c_out%lw_emissivity

   * - 8
     - Lemission
     - Longwave upward emission at top-of-canopy [W m\ :sup:`-2`]
     - bc_out%lw_emission

   * - 9
     - Lup
     - Longwave upward ( emission+reflected) at top-of-canopy [W m\ :sup:`-2`]
     - -

   * - 10
     - Kup
     - Shortwave upward (reflected) at top-of-canopy [W m\ :sup:`-2`]
     - -

   * - 11
     - Qn
     - Net all-wave radiation at top-of-canopy [W m\ :sup:`-2`]
     - -

   * - 12-26
     - LCAAbs#
     - Longwave absorption rate in clear-air part of layer [W m\ :sup:`-2`]
     - lw _flux%clear_air_abs

   * - 27-41
     - LWallNet#
     - Net longwave flux into walls [W m\ :sup:`-2`]
     - lw_flux%wall_net

   * - 42-56
     - LRfNet#
     - Net longwave flux into roofs [W m\ :sup:`-2`]
     - lw_flux%roof_net

   * - 57-71
     - LRfIn#
     - Longwave flux into walls [W m-2]
     - lw_flux%roof_in

   * - 72
     - LTopNet
     - Top-of-canopy net longwave flux [W m-2]
     - lw_flux%top_net

   * - 73
     - LGrndNet
     - Net longwave flux into the ground [W m-2]
     - lw_flux%ground_net

   * - 74
     - LTopDn
     - Top-of-canopy downwelling longwave flux [W m-2]
     - lw_flux%top_dn

   * - 75-89
     - KCAAbs#
     - Shortwave absorption rate in clear-air part of layer [W m-2]
     - sw _flux%clear_air_abs

   * - 90-104
     - KWallNet#
     - Net shortwave flux into walls [W m-2]
     - sw_flux%wall_net

   * - 105-119
     - KRfNet#
     - Net shortwave flux into roofs [W m-2]
     - sw_flux%roof_net

   * - 120-134
     - KRfIn#
     - Shortwave flux into walls [W m-2]
     - sw_flux%roof_in

   * - 135
     - KTopDnDir
     - Direct shortwave flux into roofs [W m-2]
     - sw_flux%top_dn_dir

   * - 136
     - KTopNet
     - Top-of-canopy net shortwave flux [W m-2]
     - sw_flux%top_net

   * - 137
     - KGrndDnDir
     - Direct downwelling shortwave flux into the ground [W m-2]
     - sw _flux%ground_dn_dir

   * - 138
     - KGrndNet
     - Net shortwave flux into the ground [W m-2]
     - sw_flux%ground_net

.. _more_SS_details:

More background information
***************************

Vegetation single scattering albedo (SSA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **shortwave** broadband SSA is equal to the sum of the broadband reflectance and transmittance :cite:`Yang2020Sep`.
Given reflectance and transmittance spectra the SSA is calculated according to

.. math:: \text{SSA} = \ \frac{\int_{\sim 400\ \text{nm}}^{\sim 2200\ \text{nm}}{R \times S}\text{dλ}}{\int_{\sim 400\ \text{nm}}^{\sim 2200\ \text{nm}}S\text{dλ}} + \frac{\int_{\sim 400\ \text{nm}}^{\sim 2200\ \text{nm}}{T \times S}\text{dλ}}{\int_{\sim 400\ \text{nm}}^{\sim 2200\ \text{nm}}S\text{dλ}}

where :math:`R` leaf reflectance spectrum, :math:`T` leaf transmittance spectrum and :math:`S` clear-sky surface spectrum.

The integrals are performed between 400 nm and 2200 nm because this is the spectral range that RAMI5\ :sup:`5` Järvselja birch stand forest spectra are available.
This is a reasonable approximation since it is where the majority of incoming SW energy resides (as seen from the clear-sky surface spectrum in Fig. 6).

Users can use the default value of 0.46, from RAMI5 Järvselja birch stand forest tree types or calculate their own SSA (:numref:`rami5`).
There are more tree R and T profiles `here <https://rami-benchmark.jrc.ec.europa.eu/_www/phase_descr.php?strPhase=RAMI5>`__\ :sup:`5`,




.. _rami5:
.. figure:: /assets/img/SUEWS006.png
	:alt: Overview of SUEWS

	RAMI5\ :sup:`5` data used to calculate R, T, and SSA, and R, T, and SSA values: (a) top-of-atmosphere incoming solar flux and clear-sky surface spectrum :cite:`Hogan2020Dec` (b) RAMI5 R and T spectra, and (c) calculated broadband R, T, and SSA values.


The **longwave** broadband SSA could be calculated in the same way but with the integral over the thermal infra-red (8-14 𝜇m), S replaced with the Plank function at Earth surface temperature, and R + T for the thermal infra-red.
The approximation R + T = 2R can be made.
R for different materials is available at https://speclib.jpl.nasa.gov/library. The peak in the thermal infra-red is ~10 𝜇m.
Based on inspection of R profiles for several tree species SSA=0.06 is the default value.

Building albedo and emissivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use broadband values in Table C.1 of :cite:t:`Kotthaus2014Aug`.
Full spectra can be found in the `spectral library documentation <http://micromet.reading.ac.uk/spectral-library/>`__.

Ground albedo and emissivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In SUEWS-SS this is calculated as::

   (𝛼(1)*sfr(PavSurf)+𝛼(5)*sfr(GrassSurf)+𝛼(6)*sfr(BSoilSurf)+𝛼(7)*sfr(WaterSurf))/ (sfr(PavSurf) + sfr(GrassSurf) + sfr(BSoilSurf) + sfr(WaterSurf))

where 𝛼 is either the ground albedo or emissivity.

𝛼 values for the surfaces should be set by specifying surface codes in `SUEWS_SiteSelect.txt`.
Codes should correspond to existing appropriate surfaces in `SUEWS_NonVeg.txt` and `SUEWS_NonVeg.txt`.
Alternatively, new surfaces can be made in `SUEWS_NonVeg.txt` and `SUEWS_NonVeg.txt` with 𝛼 values obtained for example from the spectral library.

Consistency of SUEWS and SS parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SUEWS building and tree (evergreen+deciduous) fractions in `SUEWS_SiteSelect.txt` should be consistent with the `SUEWS_SPARTACUS.nml` `building_frac` and `veg_frac` of the lowest model layer.

Leaf area index (LAI)
^^^^^^^^^^^^^^^^^^^^^^^^

The total vertically integrated LAI provided by SUEWS is used in SS to determine the LAI and vegetation extinction coefficient in each layer.
Surface codes in `SUEWS_SiteSelect.txt` should correspond to appropriate LAI values in `SUEWS_veg.txt`.



.. References
.. ******************
.. .. un-cited:
.. Berrizbeitia SE, EJ Gago, T Muneer 2020: Empirical Models for the Estimation of Solar Sky-Diffuse Radiation.
.. A Review and Experimental Analysis.
.. *Energies*, 13, 701

.. Hogan RJ, T Quaife, R Braghiere 2018: Fast matrix treatment of 3-D radiative transfer in vegetation canopies: SPARTACUS-Vegetation 1.1.
.. *Geoscientific Model Development* 11.1, 339-350.

.. Hogan RJ 2019a: An exponential model of urban geometry for use in radiative transfer applications.
.. *Boundary-Layer Meteorology* 170.3, 357-372.

.. Hogan RJ 2019b: Flexible treatment of radiative transfer in complex urban canopies for use in weather and climate models.
.. *Boundary-Layer Meteorology* 173.1, 53-78.

.. Hogan RJ and M Matricardi 2020: Evaluating and improving the treatment of gases in radiation schemes: the Correlated K-Distribution Model
.. Intercomparison Project (CKDMIP).
.. *Geoscientific* *Model Development* 13, 6501–6521.
