# Feature: Adjust Default Values

## Context
This branch focuses on removing or adjusting default values in the SUEWS configuration system to improve usability and prevent common configuration errors. The goal is to make the model more explicit about required inputs while providing sensible defaults only where appropriate.

## GitHub Issues to Address
- **#428**: Removing or adjusting default values (WIP, pydantic) - PRIMARY ISSUE
- **#412**: Update sample data with the benchmark one (pre-processing)

## Progress Tracking
- [ ] Audit all default values in the data model
  - [ ] Identify which defaults are sensible vs problematic
  - [ ] Document rationale for each decision
  - [ ] Create list of breaking changes
- [ ] Remove problematic defaults (#428)
  - [ ] Remove defaults that mask missing required inputs
  - [ ] Keep only universally applicable defaults
  - [ ] Update validation to catch missing values
- [ ] Update sample configurations (#412)
  - [ ] Replace sample data with benchmark configurations
  - [ ] Ensure all required fields are explicitly set
  - [ ] Add comments explaining each parameter
- [ ] Update documentation
  - [ ] Document all removed defaults
  - [ ] Provide migration guide for existing users
  - [ ] Update tutorials with explicit values

## Key Decisions
- Default values should only exist for truly optional parameters
- Physical parameters should generally not have defaults
- Site-specific parameters must always be explicit
- Maintain a clear distinction between optional and required fields

## Implementation Notes
- This is a breaking change that will require a major version bump
- Need to coordinate with documentation updates
- Consider providing a migration tool or script
- Ensure error messages clearly indicate missing required fields

## Files to Modify
- `src/supy/data_model/*.py` - Remove defaults from field definitions
- `test/benchmark1/benchmark1.yml` - Update benchmark config
- `src/supy/sample_run/sample_config.yml` - Update sample config
- `docs/source/inputs/yaml/*.rst` - Update documentation
- `src/supy/_check.py` - Enhance validation messages

## Migration Strategy
1. Create comprehensive list of removed defaults
2. Provide clear error messages for missing values
3. Update all example configurations
4. Create migration guide with before/after examples
5. Consider automated migration script

## Testing Approach
- Test that missing required fields produce clear errors
- Verify benchmark still runs correctly
- Test migration of existing configurations
- Ensure backward compatibility mode if needed