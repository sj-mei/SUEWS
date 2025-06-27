# Feature: RSL Physics Fixes

## Context
This branch addresses issues with the Roughness Sublayer (RSL) physics scheme in SUEWS, particularly for short buildings and high roughness (z0) cases. The RSL scheme is critical for accurate turbulent flux calculations in urban environments.

## GitHub Issues to Address
- **#419**: Issue using MOST in RSL scheme for short buildings (RSL) - PRIMARY ISSUE
- **#338**: RSL issue - high z0 for low FAI cases (RSL, WIP)
- **#349**: Couple STEBBS to the RSL scheme (RSL, STEBBS) - CLOSED but verify

## Progress Tracking
- [ ] Fix MOST implementation for short buildings (#419)
  - [ ] Identify the issue with Monin-Obukhov Similarity Theory
  - [ ] Implement proper scaling for short building cases
  - [ ] Add validation tests for various building heights
  - [ ] Document limitations and valid ranges
- [ ] Address high z0 issues for low FAI (#338)
  - [ ] Investigate roughness length calculations
  - [ ] Fix scaling issues with low frontal area index
  - [ ] Implement bounds checking and warnings
  - [ ] Add diagnostic output for debugging
- [ ] Verify RSL-STEBBS coupling (#349)
  - [ ] Check integration is working correctly
  - [ ] Test heat flux consistency
  - [ ] Ensure proper variable exchange
- [ ] General RSL improvements
  - [ ] Review and update RSL documentation
  - [ ] Add more comprehensive test cases
  - [ ] Improve numerical stability
  - [ ] Optimise performance

## Key Decisions
- Maintain physical realism over empirical fits
- Ensure smooth transitions between RSL and inertial sublayer
- Provide clear warnings for out-of-range conditions
- Keep backward compatibility where possible

## Implementation Notes
- RSL issues often manifest in extreme urban geometries
- Short buildings (< 10m) show particular problems
- High z0/low FAI combinations need special handling
- Consider implementing alternative RSL formulations
- May need to adjust stability functions

## Files to Modify
- `src/suews/src/suews_phys_rslprof.f95` - Main RSL implementation
- `src/suews/src/suews_phys_atmmoiststab.f95` - Stability functions
- `src/suews/src/suews_phys_resist.f95` - Resistance calculations
- `src/suews/src/suews_ctrl_const.f95` - Physical constants
- `test/test_rsl_physics.py` - New test suite (create)

## Physics Background
The RSL scheme accounts for the complex flow patterns within and just above the urban canopy. Key issues include:
1. Breakdown of MOST assumptions near surfaces
2. Complex scaling with building height
3. Interaction with stability corrections
4. Proper roughness length formulations

## Testing Strategy
1. Create test cases for various building heights (2m to 200m)
2. Test extreme FAI values (0.01 to 5.0)
3. Verify flux conservation
4. Compare with field measurements where available
5. Test stability across all atmospheric conditions

## Validation Data
- Use urban flux tower data for validation
- Compare with LES results for idealised cases
- Ensure consistency with published RSL studies