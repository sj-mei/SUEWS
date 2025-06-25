Storagedrainparams
==================

**Parameters:**

.. option:: store_min <RefValue[float] | float>

   Minimum water storage capacity

   :Unit: mm
   :Default: ``0.0``

.. option:: store_max <RefValue[float] | float>

   Maximum water storage capacity

   :Unit: mm
   :Default: ``10.0``

.. option:: store_cap <RefValue[float] | float>

   Current water storage capacity - the actual storage capacity available for surface water retention. This represents the depth of water that can be stored on or in the surface before drainage begins. For paved surfaces, this might represent depression storage; for vegetated surfaces, it includes canopy interception storage.

   :Unit: mm
   :Default: ``10.0``

.. option:: drain_eq <RefValue[int] | int>

   Drainage equation selection (0: linear, 1: exponential)

   :Unit: dimensionless
   :Default: ``0``

.. option:: drain_coef_1 <RefValue[float] | float>

   Drainage coefficient 1 (rate parameter)

   :Unit: mm h^-1
   :Default: ``0.013``

.. option:: drain_coef_2 <RefValue[float] | float>

   Drainage coefficient 2 (shape parameter)

   :Unit: dimensionless
   :Default: ``1.71``

.. option:: ref <Reference (Optional)>

   :Default: Not specified

   For ``ref``, if using the Reference structure, see :doc:`reference` for details.
