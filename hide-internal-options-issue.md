# GitHub Issue: Hide options for developers or internal use to users

## Description
Following the discussion about forcing file corrections and other advanced options, we need to hide certain options that are meant for developer/internal use only from the standard user interface while keeping them available in the code.

## Context
Some SUEWS options are currently visible to all users but are intended for advanced/internal use cases only. These include experimental features, deprecated methods, and debugging options. We want to keep these options in the code but hide them from the public-facing documentation to reduce complexity and avoid confusion for standard users.

## Identified Internal Options

### Physics Method Options
- **NetRadiationMethod**: Options 11-13, 100-300, 1001-1003 (not recommended/experimental)
- **StorageHeatMethod**: Options 3-4 (ANOHM, ESTM - not recommended)
- **StabilityMethod**: Options 0-1, 2, 4 (reserved/not recommended)
- **EmissionsMethod**: Options 3, 5 (updated/experimental variants)
- **RoughnessMethod**: Option 5 (vague alternative method)
- **FAIMethod**: Option 0 (undocumented ZERO option)

### Control Parameters
- **kdownzen**: Zenithal correction for shortwave radiation
- **diagnose**: Diagnostic output level for debugging

## Tasks
- [ ] Add `_internal` attribute to internal enum values in `src/supy/data_model/model.py`
- [ ] Add `internal_only: true` to json_schema_extra for internal fields
- [ ] Update `docs/generate_datamodel_rst.py` to filter out internal options
- [ ] Add `--include-internal` flag for generating developer documentation
- [ ] Test that internal options remain functional in code
- [ ] Regenerate user documentation without internal options
- [ ] Update developer documentation to explain how to access internal options

## Implementation Details

1. **Mark internal options in data model**
   - Add metadata to enum values and fields
   - No changes to actual functionality

2. **Update documentation generator**
   - Filter internal options from descriptions
   - Skip internal fields in user documentation
   - Provide flag for complete documentation

3. **Testing**
   - Verify YAML configs with internal options still work
   - Check documentation correctly hides/shows options
   - Ensure backwards compatibility

## Acceptance Criteria
- Standard users only see recommended, well-documented options in documentation
- Internal options remain fully functional when used directly in YAML configs
- Developer documentation available with all options when needed
- No breaking changes to existing configurations
- Clear separation between user-facing and developer options

## Benefits
- Cleaner, less confusing documentation for users
- Reduced support burden from questions about experimental features
- Maintained flexibility for developers and power users
- Easy path to promote internal options to public when ready

**Assignee**: @sunt05

**Labels**: documentation, enhancement, user-experience