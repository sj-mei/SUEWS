"""Test that spurious warnings are not generated during normal YAML loading."""

import pytest
import warnings
import yaml
from supy.data_model import SUEWSConfig
from supy.data_model.surface import BldgsProperties
from supy.data_model.type import RefValue


def test_no_spurious_warnings_yaml_loading():
    """Test that loading a properly configured YAML doesn't generate spurious warnings."""
    # YAML with all required parameters properly placed
    yaml_content = """
site:
  name: "Test Site"
  location:
    latitude: 51.5
    longitude: -0.1
  time_zone: "Europe/London"
  properties:
    land_cover:
      - surface_type: "paved"
        fraction: 0.3
      - surface_type: "bldgs"
        fraction: 0.4
        parameters:
          sfr: 0.3
          bldgh: 20.0
          faibldg: 0.5
      - surface_type: "grass"
        fraction: 0.3
    conductance:
      g_max: 30.0
      g_k: 0.4
      g_q_base: 0.001
      g_q_shape: 0.5
      g_t: 3.5
      g_sm: 2.0
      kmax: 1200.0
      s1: 5.56
      s2: 0.005
    vegetation_params:
      gdd_id: 0
      sdd_id: 0
      porosity_id: 0.2
      lai:
        baset: 5.0
        basete: 10.0
        gddfull: 300.0
        laimax: 5.5
        laimin: 1.0
        sddfull: 50.0
      ie_a: 0.8
      ie_m: 1.0
    anthropogenic_emissions:
      heat:
        qf0_beu:
          weekday: [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]
      co2:
        co2pointsource: 0.1
        ef_umolco2perj: 0.2
        frfossilfuel_heat: 0.8
        frfossilfuel_nonheat: 0.7
  surfaces:
    paved:
      thermal_layers:
        dz: [0.1, 0.2, 0.3, 0.4, 0.5]
        k: [1.5, 1.5, 1.5, 1.5, 1.5]
        rho_cp: [2e6, 2e6, 2e6, 2e6, 2e6]
    bldgs:
      thermal_layers:
        dz: [0.05, 0.1, 0.15, 0.2, 0.25]
        k: [1.0, 1.0, 1.0, 1.0, 1.0]
        rho_cp: [1.5e6, 1.5e6, 1.5e6, 1.5e6, 1.5e6]
    grass:
      beta_bioco2: 0.001
      alpha_bioco2: 0.002
      resp_a: 0.1
      resp_b: 0.2
"""
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        # Load the configuration
        config_dict = yaml.safe_load(yaml_content)
        config = SUEWSConfig(**config_dict)
        
        # Check for spurious warnings
        spurious_keywords = [
            'Missing critical building parameters',
            'Missing critical conductance parameters',
            'Missing critical vegetation parameters',
            'Missing critical CO2 emission parameters',
            'Missing critical thermal layer parameters'
        ]
        
        spurious_warnings = []
        for warning in w:
            msg = str(warning.message)
            for keyword in spurious_keywords:
                if keyword in msg:
                    spurious_warnings.append(msg)
                    break
        
        # Should have no spurious warnings
        assert len(spurious_warnings) == 0, f"Found {len(spurious_warnings)} spurious warnings: {spurious_warnings}"


def test_legitimate_warnings_still_shown():
    """Test that legitimate warnings are still shown when parameters are missing."""
    # YAML missing critical parameters
    yaml_content = """
site:
  name: "Test Site"
  location:
    latitude: 51.5
    longitude: -0.1
  time_zone: "Europe/London"
  properties:
    land_cover:
      - surface_type: "bldgs"
        fraction: 0.4
        parameters:
          sfr: 0.3
          # Missing bldgh and faibldg - should warn
      - surface_type: "grass"
        fraction: 0.6
"""
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        # Load the configuration
        config_dict = yaml.safe_load(yaml_content)
        config = SUEWSConfig(**config_dict)
        
        # Should have warnings about missing building parameters
        building_warnings = [
            warn for warn in w 
            if "Missing critical building parameters" in str(warn.message)
        ]
        
        assert len(building_warnings) >= 1, "Should warn about missing building parameters"
        
        # Check the warning mentions the specific parameters
        warning_msg = str(building_warnings[0].message)
        assert "bldgh" in warning_msg
        assert "faibldg" in warning_msg


def test_direct_instantiation_warnings():
    """Test that direct instantiation without parameters still generates warnings."""
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