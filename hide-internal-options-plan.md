# Implementation Plan: Hide Internal/Developer Options from User Documentation

## Overview
This document outlines the plan to hide certain SUEWS options that are meant for internal/developer use from the standard user documentation while keeping them functional in the code.

## Identified Internal Options

### NetRadiationMethod
- **Options 11-13** (Surface temperature variants) - marked as "not recommended"
  - `LDOWN_SURFACE` (11)
  - `LDOWN_CLOUD_SURFACE` (12)
  - `LDOWN_AIR_SURFACE` (13)
- **Options 100-300** (Zenith angle correction variants) - marked as "not recommended"
  - `LDOWN_ZENITH` (100)
  - `LDOWN_CLOUD_ZENITH` (200)
  - `LDOWN_AIR_ZENITH` (300)
- **Options 1001-1003** (SPARTACUS-Surface variants) - marked as "experimental"
  - `LDOWN_SS_OBSERVED` (1001)
  - `LDOWN_SS_CLOUD` (1002)
  - `LDOWN_SS_AIR` (1003)

### StorageHeatMethod
- **Option 3** (`ANOHM`) - marked as "not recommended"
- **Option 4** (`ESTM`) - marked as "not recommended"

### StabilityMethod
- **Options 0-1** - Reserved/not used
  - `NOT_USED` (0)
  - `NOT_USED2` (1)
- **Option 2** (`HOEGSTROM`) - marked as "not recommended"
- **Option 4** (`BUSINGER_HOEGSTROM`) - marked as "not recommended"

### EmissionsMethod
- **Option 3** (`L11_UPDATED`) - Modified variant for internal use
- **Option 5** (`J19_UPDATED`) - Updated variant with CO2 emissions (experimental)

### RoughnessMethod
- **Option 5** (`ALTERNATIVE`) - Vague "Alternative method" without clear documentation

### FAIMethod
- **Option 0** (`ZERO`) - Marked as "Not documented" in code comments

### ModelControl Fields
- **`kdownzen`** - Zenithal correction for downward shortwave radiation (advanced forcing correction)
- **`diagnose`** - Level of diagnostic output (debugging/development feature)

## Implementation Strategy

### 1. Add Metadata to Mark Internal Options

#### For Enum Values
Add an `_internal` attribute to enum values that should be hidden:

```python
class NetRadiationMethod(Enum):
    '''Method for calculating net all-wave radiation (Q*)...'''
    OBSERVED = 0
    LDOWN_OBSERVED = 1
    LDOWN_CLOUD = 2
    LDOWN_AIR = 3
    LDOWN_SURFACE = 11  # _internal = True
    LDOWN_CLOUD_SURFACE = 12  # _internal = True
    LDOWN_AIR_SURFACE = 13  # _internal = True
    LDOWN_ZENITH = 100  # _internal = True
    LDOWN_CLOUD_ZENITH = 200  # _internal = True
    LDOWN_AIR_ZENITH = 300  # _internal = True
    LDOWN_SS_OBSERVED = 1001  # _internal = True
    LDOWN_SS_CLOUD = 1002  # _internal = True
    LDOWN_SS_AIR = 1003  # _internal = True
```

#### For Field Definitions
Add `internal_only` flag to `json_schema_extra`:

```python
kdownzen: Optional[FlexibleRefValue(int)] = Field(
    default=None,
    description="Use zenithal correction for downward shortwave radiation",
    json_schema_extra={"internal_only": True}
)

diagnose: int = Field(
    default=0,
    description="Level of diagnostic output (0=none, 1=basic, 2=detailed)",
    json_schema_extra={"internal_only": True}
)
```

### 2. Update Documentation Generator

Modify `docs/generate_datamodel_rst.py` to:

1. **Check for internal flags on fields:**
```python
# In generate_rst_for_model() function
for field_name, field_info in model_class.model_fields.items():
    # Check if field is marked as internal
    if isinstance(field_info.json_schema_extra, dict):
        if field_info.json_schema_extra.get("internal_only", False):
            continue  # Skip internal fields in user docs
```

2. **Filter out internal enum options from descriptions:**
```python
def filter_internal_options(enum_class, description):
    """Remove internal options from enum descriptions"""
    lines = description.split('\n')
    filtered_lines = []
    
    for line in lines:
        # Check if line describes an internal option
        if any(internal_val in line for internal_val in get_internal_values(enum_class)):
            continue
        filtered_lines.append(line)
    
    return '\n'.join(filtered_lines)

def get_internal_values(enum_class):
    """Get list of internal enum values"""
    internal_values = []
    for member in enum_class:
        if hasattr(member, '_internal') and member._internal:
            internal_values.append(str(member.value))
    return internal_values
```

3. **Update option parsing to exclude internal options:**
```python
# In parse_method_options() function
def parse_method_options(description: str, enum_class=None) -> tuple[str, list[str]]:
    """Parse method options from description string, filtering internal ones."""
    # ... existing code ...
    
    if enum_class:
        internal_values = get_internal_values(enum_class)
        options = [opt for opt in options if not any(val in opt for val in internal_values)]
    
    return main_desc, options
```

### 3. Create Developer Documentation

Add a command-line flag to generate developer documentation that includes internal options:

```python
# In main() function
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--include-internal', action='store_true', 
                    help='Include internal/developer options in documentation')
args = parser.parse_args()

# Pass flag to generation functions
generate_rst_for_model(model_class, output_dir, processed_models, 
                      all_supy_models, include_internal=args.include_internal)
```

### 4. Testing Plan

1. **Verify internal options still work:**
   - Create test YAML configs using internal options
   - Ensure they parse and run correctly
   - Verify enum integer values remain unchanged

2. **Check documentation generation:**
   - Run standard generation: internal options should be hidden
   - Run with `--include-internal`: all options should appear
   - Verify cross-references still work

3. **Backwards compatibility:**
   - Test existing YAML configs that use internal options
   - Ensure no breaking changes for current users

## Benefits

1. **Cleaner user documentation** - Users see only recommended, well-documented options
2. **Reduced confusion** - No exposure to experimental or deprecated features
3. **Maintained functionality** - All options remain available for those who need them
4. **Future flexibility** - Easy to promote internal options to public when ready

## Migration Path

1. Phase 1: Add metadata markers (no visible changes)
2. Phase 2: Update documentation generator
3. Phase 3: Generate new documentation
4. Phase 4: Add developer documentation section if needed

## Notes

- This approach maintains full backwards compatibility
- No changes to the actual SUEWS model code required
- Documentation can be regenerated at any time
- Internal options can be promoted to public by removing the `_internal` flag