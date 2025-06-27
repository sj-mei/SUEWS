# Feature: Adjust Default Values in Pydantic Models

## Context
GitHub issue #428 reports that generated default values are causing issues in model runs with partially complete YAML files. The suggestion is to remove defaults so missing values are raised and bad "fallback" values are not used.

## Progress Tracking
- [x] Analyse current default values in pydantic models
- [x] Identify problematic default values
- [ ] Phase 1: Remove zero defaults for critical parameters (CO2, irrigation)
- [ ] Phase 2: Replace arbitrary physical defaults with None
- [ ] Phase 3: Add validation/warnings
- [ ] Update dependent code
- [ ] Test with partial YAML files
- [ ] Update documentation

## Key Decisions
1. **Using None vs removing defaults**: Use None as default to make fields optional but detectable
2. **Phased approach**: Start with most critical parameters (CO2, irrigation) to test impact
3. **Validation strategy**: Add model validators to check parameter consistency
4. **Documentation**: Create migration guide for users

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