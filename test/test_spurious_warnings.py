"""Test suite for spurious validation warnings."""

import warnings
import pytest

from supy.data_model import SUEWSConfig
from supy.data_model.surface import BldgsProperties
from supy.data_model.human_activity import CO2Params
from supy.data_model.site import Conductance


class TestSpuriousWarnings:
    """Test that spurious warnings are suppressed while legitimate warnings are shown."""
    
    def test_import_no_warnings(self):
        """Test that importing SUEWSConfig doesn't generate warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # Force reimport to test import-time behavior
            import importlib
            import supy.data_model
            importlib.reload(supy.data_model)
            
        # Check no warnings were generated
        assert len(w) == 0, f"Import generated {len(w)} unexpected warnings: {[str(x.message)[:80] for x in w[:3]]}"
    
    def test_co2_params_missing_parameters(self):
        """Test that CO2Params without parameters doesn't generate warnings (trade-off)."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            co2 = CO2Params()
            
        # With top-level validation, direct object creation won't warn
        co2_warnings = [
            warn for warn in w
            if "Missing critical CO2 emission parameters" in str(warn.message)
        ]
        assert len(co2_warnings) == 0, f"Expected 0 CO2 warnings (trade-off), got {len(co2_warnings)}"
    
    def test_co2_params_complete_no_warnings(self):
        """Test that CO2Params with all parameters doesn't generate warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            co2 = CO2Params(
                co2pointsource=0.0,
                ef_umolco2perj=1.159,
                frfossilfuel_heat=0.7,
                frfossilfuel_nonheat=0.7
            )
            
        # Should have no warnings
        co2_warnings = [
            warn for warn in w
            if "Missing critical CO2 emission parameters" in str(warn.message)
        ]
        assert len(co2_warnings) == 0, "Should not warn with complete CO2 parameters"
    
    def test_conductance_missing_parameters(self):
        """Test that Conductance without parameters doesn't generate warnings (trade-off)."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cond = Conductance()
            
        # With top-level validation, direct object creation won't warn
        cond_warnings = [
            warn for warn in w
            if "Missing critical surface conductance parameters" in str(warn.message)
        ]
        assert len(cond_warnings) == 0, f"Expected 0 conductance warnings (trade-off), got {len(cond_warnings)}"
    
    def test_conductance_complete_no_warnings(self):
        """Test that Conductance with all parameters doesn't generate warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cond = Conductance(
                g_max=3.5,
                g_k=200.0,
                g_q_base=0.13,
                g_q_shape=0.7,
                g_t=30.0,
                g_sm=0.05,
                kmax=1200.0,
                s1=5.56,
                s2=0.0
            )
            
        # Should have no warnings
        cond_warnings = [
            warn for warn in w
            if "Missing critical surface conductance parameters" in str(warn.message)
        ]
        assert len(cond_warnings) == 0, "Should not warn with complete conductance parameters"
    
    def test_building_missing_parameters(self):
        """Test that BldgsProperties without critical params doesn't generate warnings (trade-off)."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # Building with significant fraction but missing critical params
            bldgs = BldgsProperties(sfr=0.3)
            
        # With top-level validation, direct object creation won't warn
        building_warnings = [
            warn for warn in w
            if "Missing critical building parameters" in str(warn.message)
        ]
        thermal_warnings = [
            warn for warn in w
            if "Missing critical thermal layer parameters" in str(warn.message)
        ]
        
        assert len(building_warnings) == 0, f"Expected 0 building warnings (trade-off), got {len(building_warnings)}"
        assert len(thermal_warnings) == 0, f"Expected 0 thermal warnings (trade-off), got {len(thermal_warnings)}"
    
    def test_building_complete_no_warnings(self):
        """Test that BldgsProperties with all parameters doesn't generate warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bldgs = BldgsProperties(
                sfr=0.3,
                faibldg=0.5,
                bldgh=20.0,
                thermal_layers={
                    "dz": [0.2, 0.1, 0.1, 0.5, 1.6],
                    "k": [1.2, 1.1, 1.1, 1.5, 1.6],
                    "rho_cp": [1.2e6, 1.1e6, 1.1e6, 1.5e6, 1.6e6]
                }
            )
            
        # Should have no warnings
        building_warnings = [
            warn for warn in w
            if "Missing critical building parameters" in str(warn.message)
        ]
        thermal_warnings = [
            warn for warn in w
            if "Missing critical thermal layer parameters" in str(warn.message)
        ]
        
        assert len(building_warnings) == 0, "Should not warn with complete building parameters"
        assert len(thermal_warnings) == 0, "Should not warn with complete thermal parameters"
    
    def test_minimal_config_creation(self):
        """Test creating a minimal SUEWSConfig doesn't generate spurious warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            # Create minimal config
            config = SUEWSConfig(
                sites=[{
                    "name": "test",
                    "gridiv": 1,
                    "properties": {
                        "lat": {"value": 51.5},
                        "lng": {"value": -0.1},
                        "alt": {"value": 10.0},
                        "timezone": {"value": 0}
                    }
                }]
            )
        
        # Filter out legitimate warnings (missing params) and check for spurious ones
        # Spurious warnings would come from internal validation during import/default creation
        spurious_warnings = [
            warn for warn in w
            if "lambda" in str(warn.filename) or 
               "default_factory" in str(warn.message) or
               warn.lineno == -1  # Internal pydantic warnings often have -1 lineno
        ]
        
        assert len(spurious_warnings) == 0, f"Found {len(spurious_warnings)} spurious warnings"