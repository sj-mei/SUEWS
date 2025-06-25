Stebbsproperties
================

**Parameters:**

.. option:: WallInternalConvectionCoefficient <RefValue[float] | float>

   Internal convection coefficient of walls and roof [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: InternalMassConvectionCoefficient <RefValue[float] | float>

   Convection coefficient of internal mass [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: FloorInternalConvectionCoefficient <RefValue[float] | float>

   Internal convection coefficient of ground floor [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: WindowInternalConvectionCoefficient <RefValue[float] | float>

   Internal convection coefficient of windows [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: WallExternalConvectionCoefficient <RefValue[float] | float>

   Initial external convection coefficient of walls and roof [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: WindowExternalConvectionCoefficient <RefValue[float] | float>

   Initial external convection coefficient of windows [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: GroundDepth <RefValue[float] | float>

   Depth of external ground (deep soil) [m]

   :Unit: m
   :Default: ``0.0``

.. option:: ExternalGroundConductivity <RefValue[float] | float>

   External ground thermal conductivity

   :Unit: W m^-1 K^-1
   :Default: ``0.0``

.. option:: IndoorAirDensity <RefValue[float] | float>

   Density of indoor air [kg m-3]

   :Unit: kg m^-3
   :Default: ``0.0``

.. option:: IndoorAirCp <RefValue[float] | float>

   Specific heat capacity of indoor air [J kg-1 K-1]

   :Unit: J kg^-1 K^-1
   :Default: ``0.0``

.. option:: WallBuildingViewFactor <RefValue[float] | float>

   Building view factor of external walls [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: WallGroundViewFactor <RefValue[float] | float>

   Ground view factor of external walls [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: WallSkyViewFactor <RefValue[float] | float>

   Sky view factor of external walls [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: MetabolicRate <RefValue[float] | float>

   Metabolic rate of building occupants [W]

   :Unit: W
   :Default: ``0.0``

.. option:: LatentSensibleRatio <RefValue[float] | float>

   Latent-to-sensible ratio of metabolic energy release of occupants [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: ApplianceRating <RefValue[float] | float>

   Power demand of single appliance [W]

   :Unit: W
   :Default: ``0.0``

.. option:: TotalNumberofAppliances <RefValue[float] | float>

   Number of appliances present in building [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: ApplianceUsageFactor <RefValue[float] | float>

   Number of appliances in use [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: HeatingSystemEfficiency <RefValue[float] | float>

   Efficiency of space heating system [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: MaxCoolingPower <RefValue[float] | float>

   Maximum power demand of cooling system [W]

   :Unit: W
   :Default: ``0.0``

.. option:: CoolingSystemCOP <RefValue[float] | float>

   Coefficient of performance of cooling system [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: VentilationRate <RefValue[float] | float>

   Ventilation rate (air changes per hour, ACH) [h-1]

   :Unit: h^-1
   :Default: ``0.0``

.. option:: IndoorAirStartTemperature <RefValue[float] | float>

   Initial indoor air temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: IndoorMassStartTemperature <RefValue[float] | float>

   Initial indoor mass temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: WallIndoorSurfaceTemperature <RefValue[float] | float>

   Initial wall/roof indoor surface temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: WallOutdoorSurfaceTemperature <RefValue[float] | float>

   Initial wall/roof outdoor surface temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: WindowIndoorSurfaceTemperature <RefValue[float] | float>

   Initial window indoor surface temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: WindowOutdoorSurfaceTemperature <RefValue[float] | float>

   Initial window outdoor surface temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: GroundFloorIndoorSurfaceTemperature <RefValue[float] | float>

   Initial ground floor indoor surface temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: GroundFloorOutdoorSurfaceTemperature <RefValue[float] | float>

   Initial ground floor outdoor surface temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: WaterTankTemperature <RefValue[float] | float>

   Initial water temperature in hot water tank [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: InternalWallWaterTankTemperature <RefValue[float] | float>

   Initial hot water tank internal wall temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: ExternalWallWaterTankTemperature <RefValue[float] | float>

   Initial hot water tank external wall temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: WaterTankWallThickness <RefValue[float] | float>

   Hot water tank wall thickness [m]

   :Unit: m
   :Default: ``0.0``

.. option:: MainsWaterTemperature <RefValue[float] | float>

   Temperature of water coming into the water tank [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: WaterTankSurfaceArea <RefValue[float] | float>

   Surface area of hot water tank cylinder [m2]

   :Unit: m^2
   :Default: ``0.0``

.. option:: HotWaterHeatingSetpointTemperature <RefValue[float] | float>

   Water tank setpoint temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: HotWaterTankWallEmissivity <RefValue[float] | float>

   Effective external wall emissivity of the hot water tank [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: DomesticHotWaterTemperatureInUseInBuilding <RefValue[float] | float>

   Initial water temperature of water held in use in building [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: InternalWallDHWVesselTemperature <RefValue[float] | float>

   Initial hot water vessel internal wall temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: ExternalWallDHWVesselTemperature <RefValue[float] | float>

   Initial hot water vessel external wall temperature [degC]

   :Unit: degC
   :Default: ``0.0``

.. option:: DHWVesselWallThickness <RefValue[float] | float>

   Hot water vessel wall thickness [m]

   :Unit: m
   :Default: ``0.0``

.. option:: DHWWaterVolume <RefValue[float] | float>

   Volume of water held in use in building [m3]

   :Unit: m^3
   :Default: ``0.0``

.. option:: DHWSurfaceArea <RefValue[float] | float>

   Surface area of hot water in vessels in building [m2]

   :Unit: m^2
   :Default: ``0.0``

.. option:: DHWVesselEmissivity <RefValue[float] | float>

   NEEDS CHECKED! NOT USED (assumed same as DHWVesselWallEmissivity) [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: HotWaterFlowRate <RefValue[float] | float>

   Hot water flow rate from tank to vessel [m3 s-1]

   :Unit: m^3 s^-1
   :Default: ``0.0``

.. option:: DHWDrainFlowRate <RefValue[float] | float>

   Flow rate of hot water held in building to drain [m3 s-1]

   :Unit: m^3 s^-1
   :Default: ``0.0``

.. option:: DHWSpecificHeatCapacity <RefValue[float] | float>

   Specific heat capacity of hot water [J kg-1 K-1]

   :Unit: J kg^-1 K^-1
   :Default: ``0.0``

.. option:: HotWaterTankSpecificHeatCapacity <RefValue[float] | float>

   Specific heat capacity of hot water tank wal [J kg-1 K-1]

   :Unit: J kg^-1 K^-1
   :Default: ``0.0``

.. option:: DHWVesselSpecificHeatCapacity <RefValue[float] | float>

   Specific heat capacity of vessels containing hot water in use in buildings [J kg-1 K-1]

   :Unit: J kg^-1 K^-1
   :Default: ``0.0``

.. option:: DHWDensity <RefValue[float] | float>

   Density of hot water in use [kg m-3]

   :Unit: kg m^-3
   :Default: ``0.0``

.. option:: HotWaterTankWallDensity <RefValue[float] | float>

   Density of hot water tank wall [kg m-3]

   :Unit: kg m^-3
   :Default: ``0.0``

.. option:: DHWVesselDensity <RefValue[float] | float>

   Density of vessels containing hot water in use [kg m-3]

   :Unit: kg m^-3
   :Default: ``0.0``

.. option:: HotWaterTankBuildingWallViewFactor <RefValue[float] | float>

   Water tank/vessel internal building wall/roof view factor [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: HotWaterTankInternalMassViewFactor <RefValue[float] | float>

   Water tank/vessel building internal mass view factor [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: HotWaterTankWallConductivity <RefValue[float] | float>

   Effective wall conductivity of the hot water tank [W m-1 K-1]

   :Unit: W m^-1 K^-1
   :Default: ``0.0``

.. option:: HotWaterTankInternalWallConvectionCoefficient <RefValue[float] | float>

   Effective internal wall convection coefficient of the hot water tank [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: HotWaterTankExternalWallConvectionCoefficient <RefValue[float] | float>

   Effective external wall convection coefficient of the hot water tank [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: DHWVesselWallConductivity <RefValue[float] | float>

   Effective wall conductivity of the hot water tank [W m-1 K-1]

   :Unit: W m^-1 K^-1
   :Default: ``0.0``

.. option:: DHWVesselInternalWallConvectionCoefficient <RefValue[float] | float>

   Effective internal wall convection coefficient of the vessels holding hot water in use in building [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: DHWVesselExternalWallConvectionCoefficient <RefValue[float] | float>

   Effective external wall convection coefficient of the vessels holding hot water in use in building [W m-2 K-1]

   :Unit: W m^-2 K^-1
   :Default: ``0.0``

.. option:: DHWVesselWallEmissivity <RefValue[float] | float>

   Effective external wall emissivity of hot water being used within building [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: HotWaterHeatingEfficiency <RefValue[float] | float>

   Efficiency of hot water system [-]

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: MinimumVolumeOfDHWinUse <RefValue[float] | float>

   Minimum volume of hot water in use [m3]

   :Unit: m^3
   :Default: ``0.0``

.. option:: ref <Reference (Optional)>

   :Default: Not specified

   For ``ref``, if using the Reference structure, see :doc:`reference` for details.
