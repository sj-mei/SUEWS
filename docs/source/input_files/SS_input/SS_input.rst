SUEWS-SPARTACUS (SS) input files
--------------------------------



To run SUEWS-SS the SS specific files that need to be modified are:

- `RunControl.nml`

- `SUEWS_SPARTACUS.nml`

Non-SS specific SUEWS input file parameters also need to have appropriate values.
For example, LAI, albedos and emissivities are used by SUEWS-SS as explained in `more_SS_details`.


RunControl.nml
~~~~~~~~~~~~~~~
See `NetRadiationMethod` (Sensible values are 1001, 1002 or 1003) in `RunControl.nml` parameter.


SUEWS_SPARTACUS.nml
~~~~~~~~~~~~~~~~~~~~


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

