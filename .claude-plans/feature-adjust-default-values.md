# Feature: Adjust Default Values in Pydantic Models

## Context
GitHub issue #428 reports that generated default values are causing issues in model runs with partially complete YAML files. The suggestion is to remove defaults so missing values are raised and bad "fallback" values are not used.

## Progress Tracking
- [x] Analyse current default values in pydantic models
- [x] Identify problematic default values
- [x] Phase 1: Update zero defaults for critical parameters (CO2, irrigation, OHM)
- [x] Debug runtime issue: Resolved by rebuilding with `make dev`
- [x] Update dependent code (to_df_state methods handle None → 0.0 conversion)
- [x] Test with partial YAML files - works correctly
- [x] Phase 2: Replace arbitrary physical defaults with None
  - [x] Conductance parameters (11 fields)
  - [x] Drainage coefficients (4 fields)
  - [x] Thermal layer parameters (3 fields)
  - [x] Vegetation parameters (5 fields: beta_bioco2, alpha_bioco2, resp_a, resp_b, theta_bioco2)
  - [x] Tree properties (4 fields: evetreeh, dectreeh, faievetree, faidectree)
  - [x] Snow parameters (3 fields: preciplimit, snowdensmin, snowdensmax)
  - [x] Soil parameters (3 fields: soildepth, soilstorecap, sathydraulicconduct)
  - [x] Building parameters (2 fields: bldgh, faibldg)
  - [x] Additional parameters (h_maintain, startdls, enddls)
- [x] Phase 3: Add validation/warnings
  - [x] CO2 emission parameters validation
  - [x] Conductance parameters validation
  - [x] Building parameters validation (checks if fraction > 5%)
  - [x] Vegetation parameters validation
  - [x] Thermal layer parameters validation
  - [x] Test validation warnings
- [x] Update documentation
  - [x] Created comprehensive migration guide (MIGRATION_GUIDE_DEFAULT_VALUES.md)
  - [x] Documented all changed parameters with previous defaults
  - [x] Provided best practices and examples

## Key Decisions
1. **Using None vs removing defaults**: Use None as default to make fields optional but detectable
2. **Phased approach**: Start with most critical parameters (CO2, irrigation) to test impact
3. **Validation strategy**: Add model validators to check parameter consistency
4. **Documentation**: Create migration guide for users

## Resolution Summary

Successfully implemented Phase 1 changes:
- Changed critical zero defaults to None for CO2, irrigation, and OHM parameters
- Updated to_df_state methods to handle None → 0.0 conversion for DataFrame compatibility
- Resolved runtime issue by rebuilding with `make dev`
- Verified with partial YAML tests that missing parameters now use None instead of 0.0

This addresses the core issue in #428 where arbitrary default values were causing problems in model runs.

## Implementation Notes

### Categories of Defaults to Address:

1. **Critical Zero Defaults (Phase 1)**:
   - All CO2 emission parameters (co2pointsource, ef_umolco2perj, etc.)
   - Irrigation timing parameters (ie_start, ie_end)
   - OHM coefficients (a1, a2, a3)
   - Traffic units

2. **Arbitrary Physical Values (Phase 2)**:
   - Conductance parameters (g_max, g_k, etc.)
   - Drainage coefficients
   - Surface properties

3. **Reasonable Defaults to Keep**:
   - Initial temperatures (15.0°C)
   - Some profile factories
   - Model control parameters

### Implications of Using None:
- Fields become optional in Pydantic validation
- Need to handle None checks in code using these values
- Better error messages when values are missing
- Can add validators to ensure consistency

## Files to Modify
- `src/supy/data_model/human_activity.py` - CO2 and irrigation parameters
- `src/supy/data_model/ohm.py` - OHM coefficients
- `src/supy/data_model/surface.py` - Surface conductance parameters
- `src/supy/data_model/hydro.py` - Drainage parameters
- `src/supy/data_model/state.py` - Initial state parameters
- Any code that directly accesses these fields