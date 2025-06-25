Archetypeproperties
===================

**Parameters:**

.. option:: BuildingType <str>

   :Default: ``'SampleType'``

.. option:: BuildingName <str>

   :Default: ``'SampleBuilding'``

.. option:: BuildingCount <RefValue[int] | int>

   Number of buildings of this archetype [-]

   :Unit: dimensionless
   :Default: ``1``

.. option:: Occupants <RefValue[int] | int>

   Number of occupants present in building [-]

   :Unit: dimensionless
   :Default: ``1``

.. option:: stebbs_Height <RefValue[float] | float>

   Building height [m]

   :Unit: m
   :Default: ``10.0``

.. option:: FootprintArea <RefValue[float] | float>

   Building footprint area [m2]

   :Unit: m^2
   :Default: ``64.0``

.. option:: WallExternalArea <RefValue[float] | float>

   External wall area (including window area) [m2]

   :Unit: m^2
   :Default: ``80.0``

.. option:: RatioInternalVolume <RefValue[float] | float>

   Ratio of internal mass volume to total building volume [-]

   :Unit: dimensionless
   :Default: ``0.01``

.. option:: WWR <RefValue[float] | float>

   window to wall ratio [-]

   :Unit: dimensionless
   :Default: ``0.2``

.. option:: WallThickness <RefValue[float] | float>

   Thickness of external wall and roof (weighted) [m]

   :Unit: m
   :Default: ``20.0``

.. option:: WallEffectiveConductivity <RefValue[float] | float>

   Effective thermal conductivity of walls and roofs (weighted) [W m-1 K-1]

   :Unit: W m^-1 K^-1
   :Default: ``60.0``

.. option:: WallDensity <RefValue[float] | float>

   Effective density of the walls and roof (weighted) [kg m-3]

   :Unit: kg m^-3
   :Default: ``1600.0``

.. option:: WallCp <RefValue[float] | float>

   Effective specific heat capacity of walls and roof (weighted) [J kg-1 K-1]

   :Unit: J kg^-1 K^-1
   :Default: ``850.0``

.. option:: Wallx1 <RefValue[float] | float>

   Weighting factor for heat capacity of walls and roof [-]

   :Unit: dimensionless
   :Default: ``1.0``

.. option:: WallExternalEmissivity <RefValue[float] | float>

   Emissivity of the external surface of walls and roof [-]

   :Unit: dimensionless
   :Default: ``0.9``

.. option:: WallInternalEmissivity <RefValue[float] | float>

   Emissivity of the internal surface of walls and roof [-]

   :Unit: dimensionless
   :Default: ``0.9``

.. option:: WallTransmissivity <RefValue[float] | float>

   Transmissivity of walls and roof [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: WallAbsorbtivity <RefValue[float] | float>

   Absorbtivity of walls and roof [-]

   :Unit: dimensionless
   :Default: ``0.8``

.. option:: WallReflectivity <RefValue[float] | float>

   Reflectivity of the external surface of walls and roof [-]

   :Unit: dimensionless
   :Default: ``0.2``

.. option:: FloorThickness <RefValue[float] | float>

   Thickness of ground floor [m]

   :Unit: m
   :Default: ``0.2``

.. option:: GroundFloorEffectiveConductivity <RefValue[float] | float>

   Effective thermal conductivity of ground floor [W m-1 K-1]

   :Unit: W m^-1 K^-1
   :Default: ``0.15``

.. option:: GroundFloorDensity <RefValue[float] | float>

   Density of the ground floor [kg m-3]

   :Unit: kg m^-3
   :Default: ``500.0``

.. option:: GroundFloorCp <RefValue[float] | float>

   Effective specific heat capacity of the ground floor [J kg-1 K-1]

   :Unit: J kg^-1 K^-1
   :Default: ``1500.0``

.. option:: WindowThickness <RefValue[float] | float>

   Window thickness [m]

   :Unit: m
   :Default: ``0.015``

.. option:: WindowEffectiveConductivity <RefValue[float] | float>

   Effective thermal conductivity of windows [W m-1 K-1]

   :Unit: W m^-1 K^-1
   :Default: ``1.0``

.. option:: WindowDensity <RefValue[float] | float>

   Effective density of the windows [kg m-3]

   :Unit: kg m^-3
   :Default: ``2500.0``

.. option:: WindowCp <RefValue[float] | float>

   Effective specific heat capacity of windows [J kg-1 K-1]

   :Unit: J kg^-1 K^-1
   :Default: ``840.0``

.. option:: WindowExternalEmissivity <RefValue[float] | float>

   Emissivity of the external surface of windows [-]

   :Unit: dimensionless
   :Default: ``0.9``

.. option:: WindowInternalEmissivity <RefValue[float] | float>

   Emissivity of the internal surface of windows [-]

   :Unit: dimensionless
   :Default: ``0.9``

.. option:: WindowTransmissivity <RefValue[float] | float>

   Transmissivity of windows [-]

   :Unit: dimensionless
   :Default: ``0.9``

.. option:: WindowAbsorbtivity <RefValue[float] | float>

   Absorbtivity of windows [-]

   :Unit: dimensionless
   :Default: ``0.01``

.. option:: WindowReflectivity <RefValue[float] | float>

   Reflectivity of the external surface of windows [-]

   :Unit: dimensionless
   :Default: ``0.09``

.. option:: InternalMassDensity <RefValue[float] | float>

   Effective density of the internal mass [kg m-3]

   :Unit: kg m^-3
   :Default: ``0.0``

.. option:: InternalMassCp <RefValue[float] | float>

   Specific heat capacity of internal mass [J kg-1 K-1]

   :Unit: J kg^-1 K^-1
   :Default: ``0.0``

.. option:: InternalMassEmissivity <RefValue[float] | float>

   Emissivity of internal mass [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: MaxHeatingPower <RefValue[float] | float>

   Maximum power demand of heating system [W]

   :Unit: W
   :Default: ``0.0``

.. option:: WaterTankWaterVolume <RefValue[float] | float>

   Volume of water in hot water tank [m3]

   :Unit: m^3
   :Default: ``0.0``

.. option:: MaximumHotWaterHeatingPower <RefValue[float] | float>

   Maximum power demand of water heating system [W]

   :Unit: W
   :Default: ``0.0``

.. option:: HeatingSetpointTemperature <RefValue[float] | float>

   Heating setpoint temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: CoolingSetpointTemperature <RefValue[float] | float>

   Cooling setpoint temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: ref <Reference (Optional)>

   :Default: Not specified

   For ``ref``, if using the Reference structure, see :doc:`reference` for details.
