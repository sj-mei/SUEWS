# Migration Guide: Default Value Changes in SUEWS Data Models

## Overview

As of version 2025.6.2, SUEWS has updated its data model to remove many default values that could lead to incorrect model runs when using partially complete YAML configuration files. This change addresses issue [#428](https://github.com/UMEP-dev/SUEWS/issues/428).

## Why This Change?

Previously, many physical parameters had arbitrary default values (e.g., 0.0 for CO2 emissions, 15.0m for tree heights). When users provided incomplete YAML files, these defaults would be silently used, potentially leading to unrealistic model outputs without any warning.

Now, these parameters default to `None`, which means:
1. Missing parameters will trigger validation warnings
2. Users are explicitly notified about missing values
3. The model still runs (backward compatibility) but with clear warnings

## What Changed?

### Phase 1: Critical Zero Defaults
Parameters that previously defaulted to 0.0 but are critical for accurate simulations:

- **CO2 Emissions**: `co2pointsource`, `ef_umolco2perj`, `enef_v_jkm`, `frfossilfuel_heat`, `frfossilfuel_nonheat`, `maxfcmetab`, `maxqfmetab`, `minfcmetab`, `minqfmetab`, `trafficunits`
- **Irrigation Timing**: `ie_start`, `ie_end`, `internalwateruse_h`, `h_maintain`  
- **OHM Coefficients**: `a1`, `a2`, `a3`
- **Daylight Savings**: `startdls`, `enddls`

### Phase 2: Arbitrary Physical Defaults
Parameters that had non-zero defaults representing physical properties:

- **Conductance** (11 parameters): `g_max`, `g_k`, `g_q_base`, `g_q_shape`, `g_t`, `g_sm`, `kmax`, `s1`, `s2`, `tl`, `th`
- **Drainage** (4 parameters): `store_max`, `store_cap`, `drain_coef_1`, `drain_coef_2`
- **Thermal Layers** (3 parameters): `dz`, `k`, `rho_cp`
- **Vegetation** (9 parameters): `beta_bioco2`, `alpha_bioco2`, `resp_a`, `resp_b`, `theta_bioco2`, `baset`, `gddfull`, `basete`, `sddfull`, `laimax`
- **Tree Properties** (4 parameters): `evetreeh`, `dectreeh`, `faievetree`, `faidectree`
- **Snow** (3 parameters): `preciplimit`, `snowdensmin`, `snowdensmax`
- **Soil** (3 parameters): `soildepth`, `soilstorecap`, `sathydraulicconduct`
- **Buildings** (2 parameters): `bldgh`, `faibldg`

### Phase 3: Validation Warnings
New validation warnings are issued when critical parameters are missing:

```python
UserWarning: Missing critical CO2 emission parameters which may affect model accuracy:
  - co2pointsource (CO2 point source emission factor)
  - ef_umolco2perj (CO2 emission factor per unit of fuel energy)
  - frfossilfuel_heat (Fraction of heating energy from fossil fuels)
  - frfossilfuel_nonheat (Fraction of non-heating energy from fossil fuels)
Consider providing values for these parameters in your configuration.
```

## Migration Steps

### 1. Review Your YAML Files

Check if your YAML configuration files explicitly set all the parameters your model needs. Parameters that were previously omitted (relying on defaults) now need to be explicitly specified.

### 2. Run Your Model and Check Warnings

When you run SUEWS with your existing configuration, you'll see warnings for any missing critical parameters:

```bash
python -m supy.run your_config.yaml
# Look for UserWarning messages about missing parameters
```

**Note**: Warnings appear during model initialization, not during the run. They are Python UserWarnings that can be captured or filtered using Python's warnings module if needed.

### 3. Update Your Configuration

Add the missing parameters to your YAML file. For example:

```yaml
# Before (relying on defaults)
site:
  - name: "My Site"
    properties:
      land_cover:
        bldgs:
          sfr: 0.3  # 30% building fraction

# After (explicit values)
site:
  - name: "My Site"  
    properties:
      land_cover:
        bldgs:
          sfr: 0.3
          bldgh: 15.0      # Building height in meters
          faibldg: 0.35    # Frontal area index
```

### 4. Use Previous Default Values (If Appropriate)

If the previous default values were appropriate for your use case, you can explicitly set them in your configuration. Here are the previous defaults:

| Parameter | Previous Default | Unit | Description |
|-----------|-----------------|------|-------------|
| **CO2 Emissions** |
| co2pointsource | 0.0 | kg m⁻² s⁻¹ | CO2 point source emission |
| ef_umolco2perj | 0.0 | μmol J⁻¹ | CO2 emission factor |
| frfossilfuel_heat | 0.0 | - | Fossil fuel fraction for heating |
| frfossilfuel_nonheat | 0.0 | - | Fossil fuel fraction for non-heating |
| **Conductance** |
| g_max | 30.0 | mm s⁻¹ | Maximum surface conductance |
| g_k | 0.2 | - | Solar radiation parameter |
| g_sm | 0.5 | - | Soil moisture parameter |
| s1 | 0.05 | - | Lower soil moisture threshold |
| s2 | 200.0 | mm | Soil moisture parameter |
| **Vegetation** |
| beta_bioco2 | 0.6 | - | Biogenic CO2 coefficient |
| alpha_bioco2 | 0.8 | - | Biogenic CO2 coefficient |
| resp_a | 1.0 | μmol m⁻² s⁻¹ | Respiration coefficient |
| resp_b | 1.1 | - | Respiration coefficient |
| **Trees** |
| evetreeh | 15.0 | m | Evergreen tree height |
| dectreeh | 15.0 | m | Deciduous tree height |
| faievetree | 0.1 | - | Evergreen tree FAI |
| faidectree | 0.1 | - | Deciduous tree FAI |
| **Buildings** |
| bldgh | 10.0 | m | Building height |
| faibldg | 0.3 | - | Building FAI |
| **Soil** |
| soildepth | 150.0 | mm | Soil depth |
| soilstorecap | 150.0 | mm | Soil storage capacity |
| sathydraulicconduct | 0.0001 | mm s⁻¹ | Saturated hydraulic conductivity |

## Best Practices

1. **Always specify critical parameters explicitly** - Don't rely on defaults for physical parameters
2. **Use validation warnings as a checklist** - The warnings indicate which parameters might be important for your simulation
3. **Document your parameter choices** - Include comments in your YAML explaining why certain values were chosen
4. **Test with known configurations** - Verify your model outputs against previous results

## Additional Resources

For detailed parameter descriptions and guidance:
- [SUEWS Documentation](https://suews.readthedocs.io/)
- [Site Parameters Guide](https://suews.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo.html)
- [Model Physics Parameters](https://suews.readthedocs.io/en/latest/input_files/RunControl/RunControl.html)
- [Surface Properties Reference](https://suews.readthedocs.io/en/latest/notation.html)

For specific parameter groups:
- [Building Parameters](https://suews.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo.html#suews-surface-characteristics)
- [Vegetation Parameters](https://suews.readthedocs.io/en/latest/input_files/SUEWS_Veg.html)
- [Conductance Parameters](https://suews.readthedocs.io/en/latest/parameterisations.html#conductances)
- [Thermal Properties](https://suews.readthedocs.io/en/latest/recent-publications.html#surface-temperature-modelling)

## Backward Compatibility

The model maintains backward compatibility by:
1. Converting `None` values to appropriate defaults in the DataFrame state (for FORTRAN compatibility)
2. Issuing warnings rather than errors for missing parameters
3. Continuing to run with default values when parameters are missing

However, we strongly recommend updating your configurations to explicitly specify all required parameters.

## Example: Complete Building Configuration

```yaml
site:
  - name: "Urban Site"
    properties:
      land_cover:
        bldgs:
          sfr: 0.25               # 25% building fraction
          bldgh: 20.0            # 20m average building height
          faibldg: 0.4           # Frontal area index
          emis: 0.95             # Emissivity
          alb: 0.15              # Albedo
          # Storage and drainage
          storedrainprm:
            store_max: 10.0
            store_cap: 5.0
            drain_coef_1: 0.1
            drain_coef_2: 0.05
          # Thermal properties
          thermal_layers:
            dz: [0.01, 0.05, 0.1, 0.2, 0.5]
            k: [1.2, 1.0, 0.8, 0.6, 0.5]
            rho_cp: [2.0e6, 1.8e6, 1.6e6, 1.4e6, 1.2e6]
```

## Getting Help

If you encounter issues or need assistance:
1. Check the validation warnings for guidance on missing parameters
2. Refer to the SUEWS documentation for parameter descriptions
3. Open an issue on GitHub if you believe a default should be retained