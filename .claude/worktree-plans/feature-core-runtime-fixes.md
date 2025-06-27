# Feature: Core Runtime Fixes

## Context
This branch addresses critical runtime bugs and stability issues in SUEWS/SuPy core functionality. Focus on fixing bugs that affect model execution and data processing reliability.

## GitHub Issues to Address
- **#391**: Test reinitialisation of SuPy (bug, STEBBS)
- **#406**: LAI and SMD not correct using sample data (pre-processing, user question)
- **#335**: Water balance issue (water balance, WIP, P0) - CLOSED but verify fix
- **#354**: Issue with STEBBS water system when in>out - CLOSED but verify fix

## Progress Tracking
- [ ] Investigate SuPy reinitialisation issues (#391)
  - [ ] Create test case for reinitialisation
  - [ ] Fix any state persistence issues
  - [ ] Ensure STEBBS properly resets
- [ ] Fix LAI and SMD sample data issues (#406)
  - [ ] Verify sample data processing
  - [ ] Fix any data type conversion issues
  - [ ] Update documentation if needed
- [ ] Verify water balance fixes are complete
  - [ ] Run comprehensive water balance tests
  - [ ] Check edge cases (in>out scenarios)
- [ ] General runtime stability improvements
  - [ ] Fix any remaining FlexibleRefValue access issues
  - [ ] Ensure proper error handling and recovery

## Key Decisions
- Prioritise fixes that affect multiple users
- Maintain backward compatibility where possible
- Add comprehensive tests for each fix
- Document any breaking changes clearly

## Implementation Notes
- Many runtime issues stem from the pydantic migration
- Pay attention to data type conversions and validation
- STEBBS water system needs particular attention
- Consider performance implications of fixes

## Files to Modify
- `src/supy/_run.py` - Core runtime logic
- `src/supy/_check.py` - Validation and preprocessing
- `src/supy/data_model/` - Data model fixes
- `src/suews/src/suews_phys_stebbs.f95` - STEBBS water system
- `test/test_supy.py` - Add regression tests

## Testing Strategy
- Add unit tests for each bug fix
- Create integration tests for reinitialisation
- Verify against benchmark data
- Test with various input configurations