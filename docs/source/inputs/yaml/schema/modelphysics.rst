Modelphysics
============

Model physics configuration options.

Key method interactions:
- diagmethod: Determines HOW near-surface values (2m temp, 10m wind) are calculated from forcing data
- stabilitymethod: Provides stability correction functions used BY diagmethod calculations  
- localclimatemethod: Uses the near-surface values FROM diagmethod to modify vegetation processes
- gsmodel: Stomatal conductance model that may be influenced by localclimatemethod adjustments

**Parameters:**

.. option:: netradiationmethod <RefValue[NetRadiationMethod] | NetRadiationMethod>

   Method used to calculate net all-wave radiation (Q*). Options include observed values, modelled with various longwave parameterisations, and SPARTACUS-Surface integration

   :Unit: dimensionless
   :Default: ``3``

.. option:: emissionsmethod <RefValue[EmissionsMethod] | EmissionsMethod>

   Method used to calculate anthropogenic heat flux (QF) and CO2 emissions. Options include observed values, Loridan et al. (2011) SAHP, Järvi et al. (2011) SAHP_2, and Järvi et al. (2019) methods

   :Unit: dimensionless
   :Default: ``2``

.. option:: storageheatmethod <RefValue[StorageHeatMethod] | StorageHeatMethod>

   Method used to calculate storage heat flux (ΔQS). Options include observed values, Objective Hysteresis Model (OHM), AnOHM, Element Surface Temperature Method (ESTM), and extended ESTM

   :Unit: dimensionless
   :Default: ``1``

.. option:: ohmincqf <RefValue[OhmIncQf] | OhmIncQf>

   Whether to include anthropogenic heat flux (QF) in OHM storage heat calculations. 0: use Q* only, 1: use Q*+QF

   :Unit: dimensionless
   :Default: ``0``

.. option:: roughlenmommethod <RefValue[RoughnessMethod] | RoughnessMethod>

   Method used to calculate momentum roughness length (z0). Options include fixed values, variable based on vegetation, MacDonald (1998), and Grimmond & Oke (1999) methods

   :Unit: dimensionless
   :Default: ``2``

.. option:: roughlenheatmethod <RefValue[RoughnessMethod] | RoughnessMethod>

   Method used to calculate heat roughness length (z0h). Options include fixed values, variable based on vegetation, MacDonald (1998), and Grimmond & Oke (1999) methods

   :Unit: dimensionless
   :Default: ``2``

.. option:: stabilitymethod <RefValue[StabilityMethod] | StabilityMethod>

   Method used for atmospheric stability correction functions. Options include Dyer (1974)/Högström (1988), Campbell & Norman (1998), and Businger et al. (1971) formulations

   :Unit: dimensionless
   :Default: ``3``

.. option:: smdmethod <RefValue[SMDMethod] | SMDMethod>

   Method used to calculate soil moisture deficit (SMD). Options include modelled using parameters, or observed volumetric/gravimetric soil moisture from forcing file

   :Unit: dimensionless
   :Default: ``0``

.. option:: waterusemethod <RefValue[WaterUseMethod] | WaterUseMethod>

   Method used to calculate external water use for irrigation. Options include modelled using parameters or observed values from forcing file

   :Unit: dimensionless
   :Default: ``0``

.. option:: diagmethod <RefValue[DiagMethod] | DiagMethod>

   Method for calculating near-surface meteorological diagnostics (2m temperature, 2m humidity, 10m wind speed). Options: 0 (MOST) = Monin-Obukhov Similarity Theory for homogeneous surfaces; 1 (RST) = Roughness Sublayer Theory for heterogeneous urban surfaces; 2 (VARIABLE) = Automatic selection based on surface morphology (plan area index, frontal area index, and roughness element heights)

   :Unit: dimensionless
   :Default: ``2``

.. option:: faimethod <RefValue[FAIMethod] | FAIMethod>

   Method used to calculate frontal area index (FAI). Options include fixed values or variable based on vegetation state

   :Unit: dimensionless
   :Default: ``1``

.. option:: localclimatemethod <RefValue[LocalClimateMethod] | LocalClimateMethod>

   Method for incorporating urban microclimate feedbacks on vegetation and evapotranspiration. Options: 0 (NONE) = No local climate adjustments, use forcing file meteorology directly; 1 (BASIC) = Simple adjustments for urban temperature effects on leaf area index and growing degree days; 2 (DETAILED) = Comprehensive feedbacks including moisture stress, urban CO2 dome effects, and modified phenology cycles

   :Unit: dimensionless
   :Default: ``0``

.. option:: gsmodel <RefValue[GSModel]>

   Stomatal conductance parameterisation method for vegetation surfaces. Options: 1 (JARVI) = Original parameterisation (Järvi et al. 2011) based on environmental controls; 2 (WARD) = Updated parameterisation (Ward et al. 2016) with improved temperature and VPD responses

   :Unit: dimensionless
   :Default: ``GSModel.WARD``

.. option:: snowuse <RefValue[SnowUse] | SnowUse>

   Whether to include snow calculations in the model. 0: snow calculations disabled, 1: snow calculations enabled

   :Unit: dimensionless
   :Default: ``0``

.. option:: stebbsmethod <RefValue[StebbsMethod] | StebbsMethod>

   Method used for STEBBS (Surface Temperature Energy Balance Based Scheme) calculations. Options include none, default parameters, or user-provided parameters

   :Unit: dimensionless
   :Default: ``0``

.. option:: ref <Reference (Optional)>

   :Default: Not specified

   For ``ref``, if using the Reference structure, see :doc:`reference` for details.
