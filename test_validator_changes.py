#!/usr/bin/env python
"""
Quick test to verify our validator improvements work correctly.
"""

import sys
from src.supy.data_model.core import SUEWSConfig, _unwrap_value
from src.supy.data_model.type import RefValue
from pydantic import ValidationError

def test_unwrap_helper():
    """Test the _unwrap_value helper function."""
    print("Testing _unwrap_value helper...")
    
    # Test raw value
    assert _unwrap_value(42) == 42
    print("✓ Raw value unwrapping works")
    
    # Test RefValue
    ref = RefValue(value=100)
    assert _unwrap_value(ref) == 100
    print("✓ RefValue unwrapping works")
    
    return True

def test_albedo_bounds():
    """Test albedo bounds validation."""
    print("\nTesting albedo bounds validation...")
    
    # Test invalid albedo below 0
    config_data = {
        "sites": [{
            "gridiv": 1,
            "properties": {
                "land_cover": {
                    "grass": {
                        "sfr": 1.0,
                        "alb_min": -0.1,  # Invalid
                        "alb_max": 0.3
                    }
                }
            }
        }]
    }
    
    try:
        SUEWSConfig(**config_data)
        print("✗ Failed to catch negative albedo")
        return False
    except ValidationError as e:
        if "alb_min (-0.1) must be in range [0, 1]" in str(e):
            print("✓ Negative albedo validation works")
        else:
            print(f"✗ Wrong error message: {e}")
            return False
    
    # Test albedo equality (should fail with strict inequality)
    config_data["sites"][0]["properties"]["land_cover"]["grass"]["alb_min"] = 0.3
    config_data["sites"][0]["properties"]["land_cover"]["grass"]["alb_max"] = 0.3
    
    try:
        SUEWSConfig(**config_data)
        print("✗ Failed to catch albedo equality")
        return False
    except ValidationError as e:
        if "alb_min (0.3) must be less than alb_max (0.3)" in str(e):
            print("✓ Albedo strict inequality validation works")
        else:
            print(f"✗ Wrong error message: {e}")
            return False
    
    return True

def test_output_frequency():
    """Test output frequency validation."""
    print("\nTesting output frequency validation...")
    
    # Test negative frequency
    config_data = {
        "model": {
            "control": {
                "tstep": 300,
                "output_file": {
                    "format": "txt",
                    "freq": -3600  # Invalid
                }
            }
        },
        "sites": [{"gridiv": 1}]
    }
    
    try:
        SUEWSConfig(**config_data)
        print("✗ Failed to catch negative frequency")
        return False
    except ValidationError as e:
        if "Output frequency must be positive, got -3600s" in str(e):
            print("✓ Negative frequency validation works")
        else:
            print(f"✗ Wrong error message: {e}")
            return False
    
    return True

def main():
    """Run all tests."""
    print("Testing validator improvements...\n")
    
    tests = [
        test_unwrap_helper,
        test_albedo_bounds,
        test_output_frequency
    ]
    
    all_passed = True
    for test in tests:
        try:
            if not test():
                all_passed = False
        except Exception as e:
            print(f"✗ Test {test.__name__} failed with exception: {e}")
            all_passed = False
    
    print("\n" + "="*50)
    if all_passed:
        print("✓ All validator improvement tests passed!")
    else:
        print("✗ Some tests failed!")
    
    return all_passed

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)