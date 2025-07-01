"""Test that spurious warnings are not generated during normal YAML loading."""

import pytest
import warnings
import yaml
from supy.data_model import SUEWSConfig
from supy.data_model.surface import BldgsProperties
from supy.data_model.type import RefValue


def test_direct_instantiation_warnings():
    """Test that direct instantiation without parameters generates warnings."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        # Create BldgsProperties with significant fraction but no params
        bldgs = BldgsProperties(sfr=RefValue(0.3))
        
        # Should have warnings
        building_warnings = [
            warn for warn in w 
            if "Missing critical building parameters" in str(warn.message)
        ]
        
        assert len(building_warnings) >= 1, "Should warn when creating objects without critical params"
        
        # Check the warning mentions the specific parameters
        warning_msg = str(building_warnings[0].message)
        assert "bldgh" in warning_msg
        assert "faibldg" in warning_msg


def test_no_warnings_with_all_params():
    """Test that no warnings are generated when all parameters are provided."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        # Create BldgsProperties with all required params
        bldgs = BldgsProperties(
            sfr=RefValue(0.3),
            bldgh=RefValue(20.0),
            faibldg=RefValue(0.5)
        )
        
        # Should have no building warnings
        building_warnings = [
            warn for warn in w 
            if "Missing critical building parameters" in str(warn.message)
        ]
        
        assert len(building_warnings) == 0, "Should not warn when all params provided"


def test_no_spurious_warnings_during_import():
    """Test that importing the module doesn't generate spurious warnings."""
    # This test verifies that our fix for default=[Site()] works
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        # Re-import to check for warnings
        import importlib
        import supy.data_model.core
        importlib.reload(supy.data_model.core)
        
        # Check for any validation warnings during import
        validation_warnings = [
            warn for warn in w
            if any(keyword in str(warn.message) for keyword in [
                'Missing critical',
                'parameters which may affect'
            ])
        ]
        
        assert len(validation_warnings) == 0, f"Import generated {len(validation_warnings)} validation warnings"


def test_validation_warnings_shown_when_appropriate():
    """Test that validation warnings are shown for user errors but not internal construction."""
    # This is a focused test of our decorator behavior
    
    # Test 1: User creates empty object - should warn
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        from supy.data_model.human_activity import CO2Params
        
        co2 = CO2Params()
        
        co2_warnings = [
            warn for warn in w
            if "Missing critical CO2 emission parameters" in str(warn.message)
        ]
        
        assert len(co2_warnings) >= 1, "Should warn when user creates object without params"
    
    # Test 2: Object created with values - should not warn
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        co2 = CO2Params(
            co2pointsource=0.1,
            ef_umolco2perj=0.2,
            frfossilfuel_heat=0.8,
            frfossilfuel_nonheat=0.7
        )
        
        co2_warnings = [
            warn for warn in w
            if "Missing critical CO2 emission parameters" in str(warn.message)
        ]
        
        assert len(co2_warnings) == 0, "Should not warn when params provided"