# Important Update: New Validation Approach for Issue #400

Hi Silvia,

I've made significant changes to the parameter validation approach that affects your work on issue #400. Please read this carefully before continuing with the conditional validation implementation.

## What Changed

### Old Approach (Component-Level Validation)
- Validation warnings were generated at component creation time (e.g., when creating `CO2Params`, `BldgsProperties`)
- This caused spurious warnings during model import and default object creation
- Users complained about false positive warnings (issue #444)

### New Approach (Top-Down Validation)
- **All validation now happens at the SUEWSConfig level**
- No warnings during component creation
- Validation only runs when explicitly called or when loading from YAML
- Clear, consolidated validation summary for users

## Key Changes You Need to Know

### 1. Validation Method Location
```python
# NEW: Validation is now a method on SUEWSConfig
class SUEWSConfig(BaseModel):
    def validate_parameter_completeness(self) -> ValidationResult:
        """Run complete parameter validation."""
        # This is where ALL validation happens now
```

### 2. No More Component Warnings
```python
# OLD WAY (Don't do this anymore):
class CO2Params(BaseModel):
    @model_validator(mode='after')
    def validate_params(self):
        warnings.warn("Missing CO2 parameters...")  # NO!

# NEW WAY:
# Components have no validation warnings
# All validation logic moves to SUEWSConfig.validate_parameter_completeness()
```

### 3. Validation Results Structure
```python
@dataclass
class ValidationResult:
    """Results from parameter validation."""
    issues: List[ValidationIssue]
    total_issues: int
    sites_affected: int
    issue_types: List[str]
```

## How to Update Your Conditional Validation Work

### 1. Move All Validation Logic to Top Level

Instead of implementing validation in individual classes, add your conditional validation logic to `SUEWSConfig.validate_parameter_completeness()`:

```python
def validate_parameter_completeness(self) -> ValidationResult:
    issues = []
    
    # Existing parameter completeness checks...
    
    # ADD YOUR CONDITIONAL VALIDATION HERE:
    # Example structure:
    for site in self.sites:
        # Check model physics settings
        if self.model.physics.rslmethod == RSLMethod.VARIABLE:
            # Validate variable roughness parameters
            issues.extend(self._validate_variable_roughness(site))
        
        if self.model.physics.netradiationmethod >= 1000:
            # Validate SPARTACUS parameters
            issues.extend(self._validate_spartacus(site))
    
    return ValidationResult(issues=issues, ...)
```

### 2. Use the Existing Validation Infrastructure

The validation system now provides:
- Structured issue reporting (`ValidationIssue` class)
- Automatic YAML annotation generation
- Consolidated user feedback

### 3. Example Implementation Pattern

```python
def _validate_spartacus(self, site) -> List[ValidationIssue]:
    """Validate SPARTACUS-specific parameters."""
    issues = []
    
    # Check required SPARTACUS parameters
    if site.properties.spartacus_params is None:
        issues.append(ValidationIssue(
            level="ERROR",
            site_id=site.site_id,
            param="spartacus_params",
            message="SPARTACUS parameters required when NetRadiationMethod >= 1000",
            fix="Add spartacus_params configuration",
            example={"wall_specular_frac": 0.1, "ground_albedo": 0.2}
        ))
    
    return issues
```

## Benefits of This Approach

1. **No spurious warnings** - Users only see validation when they need it
2. **Better user experience** - Clear, actionable feedback
3. **Easier testing** - All validation in one place
4. **Consistent behavior** - Same validation for all entry points

## Migration Steps for Your Work

1. **Don't add validators to individual model classes**
2. **Move all conditional logic to `SUEWSConfig.validate_parameter_completeness()`**
3. **Use `ValidationIssue` for reporting problems**
4. **Test at the config level, not component level**

## Example Test Pattern

```python
def test_spartacus_validation():
    """Test SPARTACUS conditional validation."""
    config = SUEWSConfig(
        model={"physics": {"netradiationmethod": 1000}},
        sites=[{...}]  # Missing SPARTACUS params
    )
    
    result = config.validate_parameter_completeness()
    
    # Check for SPARTACUS issues
    spartacus_issues = [i for i in result.issues if "SPARTACUS" in i.message]
    assert len(spartacus_issues) > 0
```

## Files to Check

- `src/supy/data_model/core.py` - See `validate_parameter_completeness()` method
- `src/supy/data_model/validation_controller.py` - Existing conditional validation
- `test/test_parameter_validation_simplified.py` - New test patterns

## Questions?

Feel free to ask if anything is unclear. The key point is: **all validation now happens at the top level**, not in individual components.

Best,
Ting