"""
Test module for validators migrated from Model class to SUEWSConfig.

This module contains tests for validators that have been migrated from 
individual data model classes (Model, ModelPhysics, etc.) to the top-level
SUEWSConfig class for better centralized validation.

Each migrated validator should have comprehensive tests covering:
- Valid configurations that should pass
- Invalid configurations that should raise ValidationError
- Edge cases and boundary conditions
- Integration with existing SUEWSConfig validation
"""

import pytest
from pydantic import ValidationError

from supy.data_model import SUEWSConfig
from supy.data_model.model import OutputConfig


class TestMigratedValidators:
    """Test cases for validators migrated to SUEWSConfig."""

    def test_validate_model_output_config_valid_freq(self):
        """Test that valid output frequency configurations pass validation."""
        # Test case: freq is multiple of tstep
        config = SUEWSConfig(
            name="test_valid_output_config",
            description="Test valid output configuration",
            model={
                "control": {
                    "tstep": 300,
                    "output_file": {
                        "format": "parquet",
                        "freq": 3600  # 3600 % 300 = 0, should pass
                    }
                }
            },
            sites=[{
                "name": "test_site",
                "gridiv": 1,
                "properties": {
                    "lat": {"value": 51.5},
                    "lng": {"value": -0.1},
                    "alt": {"value": 10}
                }
            }]
        )
        # Should not raise any validation errors
        assert config.model.control.tstep == 300
        assert config.model.control.output_file.freq == 3600

    def test_validate_model_output_config_invalid_freq(self):
        """Test that invalid output frequency configurations raise ValidationError."""
        # Test case: freq is not multiple of tstep
        with pytest.raises(ValidationError, match="Output frequency.*multiple of timestep"):
            SUEWSConfig(
                name="test_invalid_output_config", 
                description="Test invalid output configuration",
                model={
                    "control": {
                        "tstep": 300,
                        "output_file": {
                            "format": "parquet",
                            "freq": 3601  # 3601 % 300 = 1, should fail
                        }
                    }
                },
                sites=[{
                    "name": "test_site",
                    "gridiv": 1,
                    "properties": {
                        "lat": {"value": 51.5},
                        "lng": {"value": -0.1},
                        "alt": {"value": 10}
                    }
                }]
            )

    def test_validate_model_output_config_none_freq(self):
        """Test that output config with None freq passes validation."""
        config = SUEWSConfig(
            name="test_none_freq",
            description="Test output config with None freq",
            model={
                "control": {
                    "tstep": 300,
                    "output_file": {
                        "format": "parquet",
                        "freq": None  # Should pass validation
                    }
                }
            },
            sites=[{
                "name": "test_site",
                "gridiv": 1,
                "properties": {
                    "lat": {"value": 51.5},
                    "lng": {"value": -0.1},
                    "alt": {"value": 10}
                }
            }]
        )
        assert config.model.control.output_file.freq is None

    def test_validate_model_output_config_deprecated_string(self):
        """Test that deprecated string output_file raises deprecation warning."""
        import warnings
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            config = SUEWSConfig(
                name="test_deprecated_string",
                description="Test deprecated string output_file",
                model={
                    "control": {
                        "tstep": 300,
                        "output_file": "custom_output.txt"  # Deprecated format
                    }
                },
                sites=[{
                    "name": "test_site",
                    "gridiv": 1,
                    "properties": {
                        "lat": {"value": 51.5},
                        "lng": {"value": -0.1},
                        "alt": {"value": 10}
                    }
                }]
            )
            
            # Should have issued deprecation warning
            deprecation_warnings = [warning for warning in w if issubclass(warning.category, DeprecationWarning)]
            assert len(deprecation_warnings) > 0
            assert "deprecated and was never used" in str(deprecation_warnings[0].message)

    def test_validate_model_output_config_default_string_no_warning(self):
        """Test that default 'output.txt' doesn't raise deprecation warning."""
        import warnings
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            config = SUEWSConfig(
                name="test_default_string",
                description="Test default string output_file",
                model={
                    "control": {
                        "tstep": 300,
                        "output_file": "output.txt"  # Default - should not warn
                    }
                },
                sites=[{
                    "name": "test_site",
                    "gridiv": 1,
                    "properties": {
                        "lat": {"value": 51.5},
                        "lng": {"value": -0.1},
                        "alt": {"value": 10}
                    }
                }]
            )
            
            # Should not have issued deprecation warning for default value
            deprecation_warnings = [warning for warning in w if issubclass(warning.category, DeprecationWarning)]
            output_file_warnings = [w for w in deprecation_warnings if "output_file" in str(w.message)]
            assert len(output_file_warnings) == 0

    def test_validate_model_output_config_multiple_tstep_values(self):
        """Test output frequency validation with different tstep values."""
        test_cases = [
            (300, 600, True),   # 600 % 300 = 0 ✓
            (300, 900, True),   # 900 % 300 = 0 ✓
            (300, 1800, True),  # 1800 % 300 = 0 ✓
            (300, 3600, True),  # 3600 % 300 = 0 ✓
            (300, 301, False),  # 301 % 300 = 1 ✗
            (600, 1200, True),  # 1200 % 600 = 0 ✓
            (600, 1201, False), # 1201 % 600 = 1 ✗
        ]
        
        for tstep, freq, should_pass in test_cases:
            config_data = {
                "name": f"test_tstep_{tstep}_freq_{freq}",
                "description": f"Test tstep={tstep}, freq={freq}",
                "model": {
                    "control": {
                        "tstep": tstep,
                        "output_file": {
                            "format": "parquet",
                            "freq": freq
                        }
                    }
                },
                "sites": [{
                    "name": "test_site",
                    "gridiv": 1,
                    "properties": {
                        "lat": {"value": 51.5},
                        "lng": {"value": -0.1},
                        "alt": {"value": 10}
                    }
                }]
            }
            
            if should_pass:
                config = SUEWSConfig(**config_data)
                assert config.model.control.tstep == tstep
                assert config.model.control.output_file.freq == freq
            else:
                with pytest.raises(ValidationError, match="Output frequency.*multiple of timestep"):
                    SUEWSConfig(**config_data)


class TestValidatorMigrationFramework:
    """Test the validator migration framework itself."""
    
    def test_migrated_validator_runs_after_base_validation(self):
        """Test that migrated validators run after basic model validation."""
        # This test ensures the migration framework works correctly
        # by confirming that model_validator(mode="after") runs after
        # the base model is fully constructed
        
        # Create minimal valid config
        config = SUEWSConfig(
            name="framework_test",
            description="Test validator migration framework",
            model={
                "control": {"tstep": 300},
            },
            sites=[{
                "name": "test_site",
                "gridiv": 1,
                "properties": {
                    "lat": {"value": 51.5},
                    "lng": {"value": -0.1},
                    "alt": {"value": 10}
                }
            }]
        )
        
        # Confirm the config was created successfully
        assert config.name == "framework_test"
        assert config.model.control.tstep == 300
        
    def test_multiple_migrated_validators_coexist(self):
        """Test that multiple migrated validators can coexist in SUEWSConfig."""
        # This will become relevant as we migrate more validators
        # For now, test that the existing migrated validators work together
        
        config = SUEWSConfig(
            name="multiple_validators_test",
            description="Test multiple migrated validators",
            model={
                "control": {
                    "tstep": 300,
                    "output_file": {
                        "format": "parquet",
                        "freq": 3600  # Valid freq
                    }
                },
            },
            sites=[{
                "name": "test_site",
                "gridiv": 1, 
                "properties": {
                    "lat": {"value": 51.5},
                    "lng": {"value": -0.1},
                    "alt": {"value": 10}
                }
            }]
        )
        
        # Both parameter completeness and output config validation should pass
        assert config.model.control.tstep == 300
        assert config.model.control.output_file.freq == 3600


class TestRadiationMethodValidator:
    """Test cases for the validate_radiation_method validator migrated to SUEWSConfig."""

    def test_validate_radiation_method_compatible_config(self):
        """Test that compatible radiation method configurations pass validation."""
        # Test case: netradiationmethod=3 with any forcing file (should pass)
        config = SUEWSConfig(
            name="test_compatible_radiation_method",
            description="Test compatible radiation method configuration",
            model={
                "control": {
                    "forcing_file": "forcing.txt"
                },
                "physics": {
                    "netradiationmethod": 3  # LDOWN_AIR method
                }
            },
            sites=[{
                "name": "test_site",
                "gridiv": 1,
                "properties": {
                    "lat": {"value": 51.5},
                    "lng": {"value": -0.1},
                    "alt": {"value": 10}
                }
            }]
        )
        # Should not raise any validation errors
        assert config.model.physics.netradiationmethod.value == 3

    def test_validate_radiation_method_incompatible_config(self):
        """Test that incompatible radiation method configurations raise ValidationError."""
        # Test case: netradiationmethod=1 with forcing.txt (should fail)
        with pytest.raises(ValidationError, match="NetRadiationMethod is set to 1.*observed Ldown"):
            SUEWSConfig(
                name="test_incompatible_radiation_method",
                description="Test incompatible radiation method configuration",
                model={
                    "control": {
                        "forcing_file": "forcing.txt"
                    },
                    "physics": {
                        "netradiationmethod": 1  # LDOWN_OBSERVED method
                    }
                },
                sites=[{
                    "name": "test_site",
                    "gridiv": 1,
                    "properties": {
                        "lat": {"value": 51.5},
                        "lng": {"value": -0.1},
                        "alt": {"value": 10}
                    }
                }]
            )

    def test_validate_radiation_method_with_refvalue_wrapper(self):
        """Test radiation method validation with RefValue wrapper objects."""
        # Test case: RefValue wrapped values should work correctly
        from supy.data_model.type import RefValue
        
        config = SUEWSConfig(
            name="test_refvalue_radiation_method",
            description="Test radiation method with RefValue wrappers",
            model={
                "control": {
                    "forcing_file": RefValue("custom_forcing.txt")
                },
                "physics": {
                    "netradiationmethod": RefValue(1)  # LDOWN_OBSERVED with custom file
                }
            },
            sites=[{
                "name": "test_site",
                "gridiv": 1,
                "properties": {
                    "lat": {"value": 51.5},
                    "lng": {"value": -0.1},
                    "alt": {"value": 10}
                }
            }]
        )
        # Should pass because forcing file is not "forcing.txt"
        # Note: RefValue(1) gets converted to NetRadiationMethod.LDOWN_OBSERVED enum
        method_val = config.model.physics.netradiationmethod.value
        assert int(method_val) == 1
        assert config.model.control.forcing_file.value == "custom_forcing.txt"

    def test_validate_radiation_method_with_enum_values(self):
        """Test radiation method validation with Enum values."""
        from supy.data_model.model import NetRadiationMethod
        
        # Test case: Direct Enum value should be handled correctly
        config = SUEWSConfig(
            name="test_enum_radiation_method",
            description="Test radiation method with Enum values",
            model={
                "control": {
                    "forcing_file": "custom_forcing.txt"
                },
                "physics": {
                    "netradiationmethod": NetRadiationMethod.LDOWN_OBSERVED  # Enum value = 1
                }
            },
            sites=[{
                "name": "test_site",
                "gridiv": 1,
                "properties": {
                    "lat": {"value": 51.5},
                    "lng": {"value": -0.1},
                    "alt": {"value": 10}
                }
            }]
        )
        # Should pass because forcing file is not "forcing.txt"
        method_val = config.model.physics.netradiationmethod.value
        assert int(method_val) == 1

    def test_validate_radiation_method_edge_case_forcing_files(self):
        """Test radiation method validation with various forcing file configurations."""
        test_cases = [
            ("forcing.txt", 1, False),    # Should fail: method=1 + forcing.txt
            ("forcing.txt", 3, True),     # Should pass: method=3 + any file
            ("custom.txt", 1, True),      # Should pass: method=1 + custom file
            ("observed_ldown.txt", 1, True),  # Should pass: method=1 + custom file
        ]
        
        for forcing_file, method, should_pass in test_cases:
            config_data = {
                "name": f"test_forcing_{forcing_file.replace('.', '_')}_method_{method}",
                "description": f"Test forcing_file={forcing_file}, method={method}",
                "model": {
                    "control": {
                        "forcing_file": forcing_file
                    },
                    "physics": {
                        "netradiationmethod": method
                    }
                },
                "sites": [{
                    "name": "test_site",
                    "gridiv": 1,
                    "properties": {
                        "lat": {"value": 51.5},
                        "lng": {"value": -0.1},
                        "alt": {"value": 10}
                    }
                }]
            }
            
            if should_pass:
                config = SUEWSConfig(**config_data)
                method_val = config.model.physics.netradiationmethod.value
                assert int(method_val) == method
                forcing_file_val = getattr(config.model.control.forcing_file, 'value', config.model.control.forcing_file)
                assert forcing_file_val == forcing_file
            else:
                with pytest.raises(ValidationError, match="NetRadiationMethod is set to 1.*observed Ldown"):
                    SUEWSConfig(**config_data)

    def test_validate_radiation_method_mixed_types(self):
        """Test radiation method validation with mixed RefValue and direct value types."""
        from supy.data_model.type import RefValue
        from supy.data_model.model import NetRadiationMethod
        
        test_cases = [
            # (forcing_file_value, netradiationmethod_value, should_pass)
            (RefValue("forcing.txt"), 1, False),           # RefValue + int: should fail
            ("forcing.txt", RefValue(1), False),           # str + RefValue: should fail  
            (RefValue("custom.txt"), RefValue(1), True),   # RefValue + RefValue: should pass
            ("custom.txt", NetRadiationMethod.LDOWN_OBSERVED, True),  # str + Enum: should pass
        ]
        
        for i, (forcing_file_val, method_val, should_pass) in enumerate(test_cases):
            config_data = {
                "name": f"test_mixed_types_{i}",
                "description": f"Test mixed types case {i}",
                "model": {
                    "control": {
                        "forcing_file": forcing_file_val
                    },
                    "physics": {
                        "netradiationmethod": method_val
                    }
                },
                "sites": [{
                    "name": "test_site",
                    "gridiv": 1,
                    "properties": {
                        "lat": {"value": 51.5},
                        "lng": {"value": -0.1},
                        "alt": {"value": 10}
                    }
                }]
            }
            
            if should_pass:
                config = SUEWSConfig(**config_data)
                # Verify the values are properly accessed
                # Note: All values get converted to Enum, so we expect the .value to be 1
                method_val = config.model.physics.netradiationmethod.value
                assert int(method_val) == 1
            else:
                with pytest.raises(ValidationError, match="NetRadiationMethod is set to 1.*observed Ldown"):
                    SUEWSConfig(**config_data)