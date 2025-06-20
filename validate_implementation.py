#!/usr/bin/env python3
"""
Final validation script for the SUEWSSimulation implementation.

This script validates all key functionality using the standalone implementation
and generates a comprehensive report.
"""

import pandas as pd
import numpy as np
from datetime import datetime
import tempfile
import yaml
from pathlib import Path
import sys

# Import our standalone implementation
from suews_simulation_standalone import SUEWSSimulation


class ValidationTester:
    """Comprehensive validation testing for SUEWSSimulation."""
    
    def __init__(self):
        self.passed_tests = 0
        self.failed_tests = 0
        self.test_results = []
    
    def test(self, test_name: str, test_func, *args, **kwargs):
        """Run a test and record the result."""
        try:
            print(f"  🔍 {test_name}...")
            result = test_func(*args, **kwargs)
            if result:
                self.passed_tests += 1
                self.test_results.append((test_name, True, None))
                print(f"    ✅ PASSED")
            else:
                self.failed_tests += 1
                self.test_results.append((test_name, False, "Test returned False"))
                print(f"    ❌ FAILED")
            return result
        except Exception as e:
            self.failed_tests += 1
            self.test_results.append((test_name, False, str(e)))
            print(f"    ❌ FAILED: {e}")
            return False
    
    def report(self):
        """Generate validation report."""
        total = self.passed_tests + self.failed_tests
        success_rate = (self.passed_tests / total * 100) if total > 0 else 0
        
        print(f"\n{'='*70}")
        print(f"VALIDATION REPORT")
        print(f"{'='*70}")
        print(f"Total tests: {total}")
        print(f"Passed: {self.passed_tests}")
        print(f"Failed: {self.failed_tests}")
        print(f"Success rate: {success_rate:.1f}%")
        
        if self.failed_tests > 0:
            print(f"\nFailed tests:")
            for test_name, passed, error in self.test_results:
                if not passed:
                    print(f"  ❌ {test_name}: {error}")
        
        return success_rate >= 95


def create_test_config():
    """Create a test configuration."""
    return {
        'name': 'validation_test',
        'model': {'control': {'tstep': 300}},
        'sites': [{'name': 'test', 'gridiv': 1, 
                  'properties': {'lat': {'value': 51.5}, 'lng': {'value': -0.1}}}]
    }


def create_test_forcing(n_hours=48):
    """Create test forcing data."""
    dates = pd.date_range('2012-06-01', periods=n_hours, freq='h')
    return pd.DataFrame({
        'Tair': np.random.normal(15, 5, n_hours),
        'RH': np.random.normal(70, 15, n_hours),
        'U': np.random.normal(3, 1, n_hours),
        'pres': np.random.normal(1013, 10, n_hours),
        'rain': np.random.exponential(0.1, n_hours),
        'kdown': np.random.normal(200, 100, n_hours).clip(0)
    }, index=dates)


# Test functions
def test_basic_instantiation():
    """Test basic class instantiation."""
    config = create_test_config()
    sim = SUEWSSimulation(config)
    return isinstance(sim, SUEWSSimulation) and sim.config is not None


def test_from_yaml_classmethod():
    """Test from_yaml class method."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(create_test_config(), f)
        yaml_path = f.name
    
    try:
        sim = SUEWSSimulation.from_yaml(yaml_path)
        return isinstance(sim, SUEWSSimulation)
    finally:
        Path(yaml_path).unlink()


def test_forcing_setup():
    """Test forcing data setup."""
    sim = SUEWSSimulation(create_test_config())
    forcing = create_test_forcing()
    sim.setup_forcing(forcing)
    return sim.forcing is not None and len(sim.forcing) > 0


def test_validation():
    """Test configuration validation."""
    sim = SUEWSSimulation(create_test_config())
    forcing = create_test_forcing()
    sim.setup_forcing(forcing)
    validation = sim.validate()
    return validation['status'] in ['valid', 'valid_with_warnings']


def test_simulation_run():
    """Test simulation execution."""
    sim = SUEWSSimulation(create_test_config())
    forcing = create_test_forcing()
    sim.setup_forcing(forcing)
    results = sim.run()
    return sim.is_complete and results is not None and len(results) > 0


def test_result_access():
    """Test result access methods."""
    sim = SUEWSSimulation(create_test_config())
    forcing = create_test_forcing()
    sim.setup_forcing(forcing)
    sim.run()
    
    # Test various access methods
    all_results = sim.get_results()
    filtered_results = sim.get_results(['QH', 'QE'])
    summary = sim.summary()
    
    return (all_results is not None and 
            filtered_results is not None and 
            summary is not None and
            'run_metadata' in summary)


def test_plotting_capability():
    """Test plotting functionality."""
    try:
        import matplotlib.pyplot as plt
        sim = SUEWSSimulation(create_test_config())
        forcing = create_test_forcing()
        sim.setup_forcing(forcing)
        sim.run()
        sim.quick_plot(['QH', 'QE'])
        plt.close('all')  # Clean up
        return True
    except ImportError:
        return True  # Skip if matplotlib not available
    except Exception:
        return False


def test_save_functionality():
    """Test save functionality."""
    sim = SUEWSSimulation(create_test_config())
    forcing = create_test_forcing()
    sim.setup_forcing(forcing)
    sim.run()
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Test CSV save
        csv_path = sim.save(Path(temp_dir) / "test.csv", format='csv')
        csv_exists = csv_path.exists() and csv_path.stat().st_size > 0
        
        # Test Excel save
        excel_path = sim.save(Path(temp_dir) / "test.xlsx", format='excel')
        excel_exists = excel_path.exists() and excel_path.stat().st_size > 0
        
        return csv_exists and excel_exists


def test_clone_and_reset():
    """Test clone and reset functionality."""
    sim = SUEWSSimulation(create_test_config())
    forcing = create_test_forcing()
    sim.setup_forcing(forcing)
    sim.run()
    
    # Test clone
    clone = sim.clone()
    clone_valid = (isinstance(clone, SUEWSSimulation) and 
                   not clone.is_complete and 
                   clone.forcing is not None)
    
    # Test reset
    original_complete = sim.is_complete
    sim.reset()
    reset_valid = original_complete and not sim.is_complete
    
    return clone_valid and reset_valid


def test_error_handling():
    """Test error handling."""
    try:
        # Test invalid config file
        try:
            SUEWSSimulation("nonexistent.yaml")
            return False  # Should have raised error
        except FileNotFoundError:
            pass  # Expected
        
        # Test run without forcing
        try:
            sim = SUEWSSimulation(create_test_config())
            sim.run()
            return False  # Should have raised error
        except (RuntimeError, ValueError):
            pass  # Expected
        
        # Test results before run
        try:
            sim = SUEWSSimulation(create_test_config())
            sim.summary()
            return False  # Should have raised error
        except RuntimeError:
            pass  # Expected
        
        return True
    except Exception:
        return False


def test_dataframe_integration():
    """Test pandas DataFrame integration."""
    sim = SUEWSSimulation(create_test_config())
    forcing = create_test_forcing()
    sim.setup_forcing(forcing)
    results = sim.run()
    
    # Test DataFrame properties
    is_dataframe = isinstance(results, pd.DataFrame)
    has_multiindex = results.columns.nlevels == 2
    has_datetime_index = pd.api.types.is_datetime64_any_dtype(results.index)
    
    # Test DataFrame operations
    mean_values = results.mean()
    subset = results.head(10)
    
    return (is_dataframe and has_multiindex and has_datetime_index and
            mean_values is not None and len(subset) == 10)


def test_properties_and_metadata():
    """Test properties and metadata."""
    sim = SUEWSSimulation(create_test_config())
    forcing = create_test_forcing()
    sim.setup_forcing(forcing)
    
    # Test properties before run
    config_valid = sim.config is not None
    forcing_valid = sim.forcing is not None
    not_complete = not sim.is_complete
    no_results = sim.results is None
    
    # Test properties after run
    sim.run()
    now_complete = sim.is_complete
    has_results = sim.results is not None
    has_metadata = len(sim.metadata) > 0
    
    return (config_valid and forcing_valid and not_complete and no_results and
            now_complete and has_results and has_metadata)


def main():
    """Run comprehensive validation."""
    print("🔬 SUEWSSimulation Comprehensive Validation")
    print("=" * 70)
    
    tester = ValidationTester()
    
    # Core functionality tests
    print("\n📋 Core Functionality Tests")
    tester.test("Basic instantiation", test_basic_instantiation)
    tester.test("YAML class method", test_from_yaml_classmethod)
    tester.test("Forcing data setup", test_forcing_setup)
    tester.test("Configuration validation", test_validation)
    tester.test("Simulation execution", test_simulation_run)
    
    # Advanced functionality tests  
    print("\n🔬 Advanced Functionality Tests")
    tester.test("Result access methods", test_result_access)
    tester.test("Plotting capability", test_plotting_capability)
    tester.test("Save functionality", test_save_functionality)
    tester.test("Clone and reset", test_clone_and_reset)
    tester.test("Error handling", test_error_handling)
    
    # Integration tests
    print("\n🔗 Integration Tests")
    tester.test("DataFrame integration", test_dataframe_integration)
    tester.test("Properties and metadata", test_properties_and_metadata)
    
    # Generate final report
    success = tester.report()
    
    if success:
        print(f"\n🎉 VALIDATION SUCCESSFUL! 🎉")
        print("The SUEWSSimulation class meets all requirements.")
        return True
    else:
        print(f"\n❌ VALIDATION FAILED")
        print("Some tests did not pass. Check the report above.")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)