Initialstatebsoil
=================

**Parameters:**

.. option:: state <RefValue[float] | float>

   Initial water state of the surface

   :Unit: mm
   :Default: ``0.0``

.. option:: soilstore <RefValue[float] | float>

   Initial soil store (essential for QE)

   :Unit: mm
   :Default: ``150.0``

.. option:: snowfrac <RefValue[float] | float | NoneType>

   Snow fraction

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: snowpack <RefValue[float] | float | NoneType>

   Snow pack

   :Unit: mm
   :Default: ``0.0``

.. option:: icefrac <RefValue[float] | float | NoneType>

   Ice fraction

   :Unit: dimensionless
   :Default: ``0.0``

.. option:: snowwater <RefValue[float] | float | NoneType>

   Snow water

   :Unit: mm
   :Default: ``0.0``

.. option:: snowdens <RefValue[float] | float | NoneType>

   Snow density

   :Unit: kg m^-3
   :Default: ``0.0``

.. option:: temperature <RefValue[List[float]] | List of float>

   Initial temperature for each thermal layer

   :Unit: degC
   :Default: ``[15.0, 15.0, 15.0, 15.0, 15.0]``

.. option:: tsfc <RefValue[float] | float | NoneType>

   Initial exterior surface temperature

   :Unit: degC
   :Default: ``15.0``

.. option:: tin <RefValue[float] | float | NoneType>

   Initial interior surface temperature

   :Unit: degC
   :Default: ``20.0``

.. option:: ref <Reference (Optional)>

   :Default: Not specified

   For ``ref``, if using the Reference structure, see :doc:`reference` for details.
