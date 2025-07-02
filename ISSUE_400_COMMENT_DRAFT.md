# Draft Comment for Issue #400

Hi @dayantur,

I've just completed some changes to the parameter validation system (see #444) that affect how we implement validation. Since you're working on the conditional validation, I wanted to share some suggestions on how to adapt your work to the new approach.

Would this approach work for what you're implementing? Happy to discuss any questions or concerns you might have about adapting to these changes.

The key point to remember: **all validation now happens at the top level**, not in individual components. This gives us much better control and eliminates the spurious warnings users were experiencing.

Please see details below.

<details>
<summary>Click to expand technical migration details</summary>

## What Changed

The validation system has been restructured to eliminate spurious warnings. All validation now happens at the top level (`SUEWSConfig`) rather than in individual components. This gives us better control over when validation runs and provides clearer feedback to users.

## Migration Suggestions

### 1. Important: No More Component-Level Validators

The key change is that we no longer use validators in individual model classes:

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

### 2. Move All Validation Logic to SUEWSConfig

Instead of adding validators to individual model classes, I'd suggest implementing all conditional validation in the `validate_parameter_completeness()` method:

```python
# In src/supy/data_model/core.py
class SUEWSConfig(BaseModel):
    def validate_parameter_completeness(self) -> ValidationResult:
        issues = []

        # Existing parameter completeness checks...

        # Add conditional validation here
        issues.extend(self._validate_conditional_parameters())

        return ValidationResult(issues=issues, ...)

    def _validate_conditional_parameters(self) -> List[ValidationIssue]:
        """Validate parameters based on active physics methods."""
        issues = []

        for site in self.sites:
            # Example: Variable roughness validation
            if self.model.physics.roughnessmethod == RoughnessMethod.VARIABLE:
                if not site.properties.variable_roughness_params:
                    issues.append(ValidationIssue(
                        level="ERROR",
                        site_id=site.site_id,
                        param="variable_roughness_params",
                        message="Variable roughness parameters required when RoughnessMethod=VARIABLE",
                        fix="Add variable_roughness_params configuration",
                        example={...}  # Example configuration
                    ))

        return issues
```

### 3. Use the ValidationIssue Structure

The new system uses a structured approach for reporting issues:

```python
@dataclass
class ValidationIssue:
    level: str  # "ERROR" or "WARNING"
    site_id: str
    param: str
    message: str
    fix: str
    example: dict
```

This allows the system to generate helpful annotated YAML files for users.

### 4. Update Tests

The tests in `test_conditional_validation.py` have been temporarily disabled (marked with `@pytest.mark.skip`). You might want to:

1. Remove the skip decorator
2. Update tests to check validation at the config level:

```python
def test_spartacus_validation():
    """Test SPARTACUS conditional validation."""
    config = SUEWSConfig(
        model={"physics": {"netradiationmethod": 1000}},
        sites=[{...}]  # Missing SPARTACUS params
    )

    # Validation happens automatically on from_yaml()
    # or can be called explicitly:
    result = config.validate_parameter_completeness()

    # Check for expected issues
    spartacus_issues = [i for i in result.issues if "SPARTACUS" in i.message]
    assert len(spartacus_issues) > 0
```

### 5. Benefits of This Approach

- No spurious warnings during model creation
- All validation logic in one place
- Easier to test and maintain
- Users get consolidated feedback with clear fixes
- Automatic YAML annotation generation

## Example Implementation Pattern

Here's how you might structure the conditional validation for different physics methods:

```python
def _validate_conditional_parameters(self) -> List[ValidationIssue]:
    issues = []

    # Group validations by physics method
    if self._needs_variable_roughness_validation():
        issues.extend(self._validate_variable_roughness())

    if self._needs_spartacus_validation():
        issues.extend(self._validate_spartacus())

    if self._needs_estm_validation():
        issues.extend(self._validate_estm_thermal())

    return issues

def _needs_spartacus_validation(self) -> bool:
    return self.model.physics.netradiationmethod >= 1000
```

## Key Files to Check

- `src/supy/data_model/core.py` - See `validate_parameter_completeness()` method
- `src/supy/data_model/validation_controller.py` - Existing conditional validation (to be migrated)
- `test/test_validation_topdown.py` - New test patterns and examples
- `test/test_conditional_validation.py` - Your tests (currently skipped, need updating)

## Migration Steps Summary

1. **Don't add validators to individual model classes**
2. **Move all conditional logic to `SUEWSConfig.validate_parameter_completeness()`**
3. **Use `ValidationIssue` for reporting problems**
4. **Test at the config level, not component level**



</details>