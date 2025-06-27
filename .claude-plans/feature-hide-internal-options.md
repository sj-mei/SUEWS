# Feature: Hide Internal/Developer Options from User Documentation

## Context
Following discussion about forcing file corrections and other advanced options, we need to hide certain options that are meant for developer/internal use only from the standard user documentation while keeping them functional in the code.

## Progress Tracking
- [x] Analyze all model options to identify internal/developer-only ones
- [x] Create implementation plan
- [x] Draft GitHub issue
- [ ] Add metadata markers to internal options in data model
- [ ] Update documentation generator to filter internal options
- [ ] Test that internal options remain functional
- [ ] Generate updated user documentation
- [ ] Create developer documentation section

## Key Decisions
- Use metadata approach rather than code changes to maintain backwards compatibility
- Add `_internal` attribute to enum values that should be hidden
- Add `internal_only` flag to field json_schema_extra
- Keep all options functional, only hide from documentation

## Identified Internal Options

### NetRadiationMethod
- Options 11-13: Surface temperature variants (not recommended)
- Options 100-300: Zenith angle variants with NARP output (not recommended)
- Options 1001-1003: SPARTACUS-Surface variants (experimental)

### StorageHeatMethod
- Option 3: ANOHM (not recommended)
- Option 4: ESTM (not recommended)

### StabilityMethod
- Options 0-1: Reserved/NOT_USED
- Option 2: HOEGSTROM (not recommended)
- Option 4: BUSINGER_HOEGSTROM (not recommended)

### EmissionsMethod
- Option 3: L11_UPDATED (internal variant)
- Option 5: J19_UPDATED (experimental with CO2)

### RoughnessMethod
- Option 5: ALTERNATIVE (vague/undocumented)

### FAIMethod
- Option 0: ZERO (marked as "Not documented")

### ModelControl Fields
- `kdownzen`: Zenithal correction for shortwave radiation
- `diagnose`: Diagnostic output level for debugging

## Implementation Notes

### Step 1: Update model.py
Add metadata to enum values. For example:
```python
class NetRadiationMethod(Enum):
    LDOWN_SURFACE = 11  # Add: _internal = True
```

Add to field definitions:
```python
kdownzen: Optional[FlexibleRefValue(int)] = Field(
    default=None,
    description="Use zenithal correction for downward shortwave radiation",
    json_schema_extra={"internal_only": True}
)
```

### Step 2: Update generate_datamodel_rst.py
1. Check for `internal_only` flag in field metadata
2. Filter enum options marked with `_internal`
3. Add `--include-internal` flag for developer docs

### Step 3: Testing
- Create test YAML with internal options
- Verify it still works
- Check documentation generation filters correctly

## Files to Modify
- `src/supy/data_model/model.py` - Add internal markers
- `docs/generate_datamodel_rst.py` - Add filtering logic
- `docs/Makefile` - Add developer docs target (optional)

## GitHub Issue
Draft available in: `hide-internal-options-issue.md`

## Notes
- This approach maintains full backwards compatibility
- No changes to SUEWS Fortran code required
- Documentation can be regenerated at any time
- Easy to promote internal options to public later