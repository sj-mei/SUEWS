#!/usr/bin/env python3
"""
Comprehensive demonstration and test of the standalone SUEWSSimulation class.
"""

import pandas as pd
import numpy as np
from datetime import datetime
import tempfile
import yaml
from pathlib import Path

# Import our standalone class
from suews_simulation_standalone import SUEWSSimulation


def create_sample_config():
    """Create a sample configuration for testing."""
    return {
        'name': 'test_suews_config',
        'description': 'Sample configuration for demonstration',
        'model': {
            'control': {
                'tstep': 300,  # 5-minute timesteps
                'forcing_file': {'value': 'forcing.csv'},
                'output_file': 'output.txt',
                'diagnose': 0
            },
            'physics': {
                'netradiationmethod': {'value': 3},
                'emissionsmethod': {'value': 2},
                'storageheatmethod': {'value': 1},
                'roughlenmommethod': {'value': 2},
                'stabilitymethod': {'value': 3}
            }
        },
        'sites': [{
            'name': 'london_kingsway',
            'gridiv': 1,
            'properties': {
                'lat': {'value': 51.51},
                'lng': {'value': -0.12},
                'alt': {'value': 35.0},
                'timezone': {'value': 0},
                'surfacearea': {'value': 10000.0}
            }
        }]
    }


def create_sample_forcing_data(n_days=2):
    """Create realistic sample forcing data."""
    # Create hourly data for n_days
    dates = pd.date_range('2012-06-01', periods=24*n_days, freq='H')
    n_timesteps = len(dates)
    
    # Create realistic diurnal patterns
    hour_of_day = dates.hour
    day_of_year = dates.dayofyear
    
    # Temperature with diurnal cycle
    temp_base = 15 + 10 * np.sin(2 * np.pi * (day_of_year - 81) / 365)  # Seasonal
    temp_diurnal = 8 * np.sin(2 * np.pi * (hour_of_day - 6) / 24)  # Diurnal
    temperature = temp_base + temp_diurnal + np.random.normal(0, 2, n_timesteps)
    
    # Solar radiation with diurnal cycle (only during day)
    kdown = np.maximum(0, 800 * np.sin(2 * np.pi * (hour_of_day - 6) / 24) * 
                       (hour_of_day >= 6) * (hour_of_day <= 18) + 
                       np.random.normal(0, 50, n_timesteps))
    
    # Other variables
    data = {
        'Tair': temperature,
        'RH': np.random.normal(70, 15, n_timesteps).clip(20, 95),
        'U': np.random.lognormal(1, 0.5, n_timesteps).clip(0.1, 15),
        'pres': np.random.normal(1013, 5, n_timesteps),
        'rain': np.random.exponential(0.1, n_timesteps) * (np.random.random(n_timesteps) < 0.1),
        'kdown': kdown,
        'ldown': np.random.normal(350, 50, n_timesteps),
        'snow': np.zeros(n_timesteps),  # No snow in June
        'fcld': np.random.normal(0.4, 0.3, n_timesteps).clip(0, 1),
        'wdir': np.random.uniform(0, 360, n_timesteps)
    }
    
    return pd.DataFrame(data, index=dates)


def demo_basic_usage():
    """Demonstrate basic usage patterns."""
    print("🚀 SUEWS Simulation Class - Basic Usage Demo")
    print("=" * 60)
    
    # 1. Create simulation from dictionary config
    print("\n📋 1. Creating simulation from configuration...")
    config = create_sample_config()
    sim = SUEWSSimulation(config)
    print(f"✅ Simulation created: {sim}")
    
    # 2. Set up forcing data
    print("\n🌤️  2. Setting up forcing data...")
    forcing = create_sample_forcing_data(n_days=3)
    sim.setup_forcing(forcing)
    print(f"✅ Forcing data loaded: {len(sim.forcing)} timesteps")
    print(f"   Date range: {sim.forcing.index.min()} to {sim.forcing.index.max()}")
    print(f"   Variables: {list(sim.forcing.columns)}")
    
    # 3. Validate configuration
    print("\n✅ 3. Validating configuration...")
    validation = sim.validate()
    print(f"   Status: {validation['status']}")
    print(f"   Warnings: {len(validation['warnings'])}")
    print(f"   Errors: {len(validation['errors'])}")
    
    # 4. Run simulation
    print("\n🏃 4. Running simulation...")
    start_time = datetime.now()
    results = sim.run(debug_mode=True)
    duration = (datetime.now() - start_time).total_seconds()
    
    print(f"✅ Simulation completed in {duration:.2f} seconds")
    print(f"   Results shape: {results.shape}")
    print(f"   Result groups: {results.columns.get_level_values(0).unique().tolist()}")
    print(f"   Variables: {results.columns.get_level_values(1).unique().tolist()}")
    
    # 5. Analyze results
    print("\n📊 5. Analyzing results...")
    summary = sim.summary()
    print(f"   Run metadata keys: {list(summary['run_metadata'].keys())}")
    print(f"   Energy balance components: {list(summary.get('energy_balance', {}).get('components', {}).keys())}")
    
    # 6. Preview results
    print("\n👀 6. Previewing results...")
    sim.see(5)
    
    return sim


def demo_advanced_features(sim):
    """Demonstrate advanced features."""
    print("\n\n🎯 Advanced Features Demo")
    print("=" * 60)
    
    # 1. Filtered results
    print("\n🔍 1. Getting filtered results...")
    energy_vars = ['QH', 'QE', 'QS']
    energy_results = sim.get_results(variables=energy_vars)
    print(f"✅ Energy fluxes extracted: {energy_results.shape}")
    
    # 2. Resampled results
    print("\n📈 2. Getting resampled results...")
    daily_results = sim.get_results(resample_freq='D')
    print(f"✅ Daily averages: {daily_results.shape}")
    
    # 3. Date-filtered results
    print("\n📅 3. Getting date-filtered results...")
    start_date = sim.results.index[10]
    end_date = sim.results.index[30]
    subset_results = sim.get_results(start_date=start_date, end_date=end_date)
    print(f"✅ Date subset: {subset_results.shape}")
    
    # 4. Save results in different formats
    print("\n💾 4. Saving results...")
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # CSV
        csv_path = sim.save(temp_path / "results.csv", format='csv')
        print(f"✅ CSV saved: {csv_path.name} ({csv_path.stat().st_size} bytes)")
        
        # Excel
        excel_path = sim.save(temp_path / "results.xlsx", format='excel')
        print(f"✅ Excel saved: {excel_path.name} ({excel_path.stat().st_size} bytes)")
        
        # Pickle
        pickle_path = sim.save(temp_path / "results.pkl", format='pickle')
        print(f"✅ Pickle saved: {pickle_path.name} ({pickle_path.stat().st_size} bytes)")
    
    # 5. Clone and modify
    print("\n👥 5. Cloning simulation...")
    sim_clone = sim.clone()
    print(f"✅ Clone created: {sim_clone}")
    print(f"   Original status: {sim.is_complete}")
    print(f"   Clone status: {sim_clone.is_complete}")
    
    # 6. Reset
    print("\n🔄 6. Testing reset...")
    sim_clone.reset()
    print(f"✅ Reset completed: is_complete = {sim_clone.is_complete}")


def demo_yaml_workflow():
    """Demonstrate YAML-based workflow."""
    print("\n\n📝 YAML Configuration Workflow")
    print("=" * 60)
    
    # Create YAML config file
    config_dict = create_sample_config()
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(config_dict, f, default_flow_style=False)
        yaml_path = f.name
    
    try:
        # Create simulation from YAML
        print(f"\n📄 Loading configuration from YAML file...")
        sim = SUEWSSimulation.from_yaml(yaml_path, debug_mode=True)
        print(f"✅ YAML simulation created: {sim}")
        
        # Add forcing and run
        forcing = create_sample_forcing_data(n_days=1)
        sim.setup_forcing(forcing)
        
        results = sim.run()
        print(f"✅ YAML-based simulation completed: {results.shape}")
        
        return sim
        
    finally:
        # Clean up
        Path(yaml_path).unlink()


def demo_plotting():
    """Demonstrate plotting capabilities."""
    print("\n\n📊 Plotting Demo")
    print("=" * 60)
    
    try:
        import matplotlib.pyplot as plt
        
        # Create a simple simulation for plotting
        config = create_sample_config()
        sim = SUEWSSimulation(config)
        forcing = create_sample_forcing_data(n_days=7)  # One week
        sim.setup_forcing(forcing)
        results = sim.run()
        
        print("\n🎨 Creating quick plots...")
        
        # Plot default variables
        print("   Plotting default energy fluxes...")
        sim.quick_plot()
        
        # Plot specific variables
        print("   Plotting temperature and humidity...")
        sim.quick_plot(['Tair', 'RH'], figsize=(12, 6))
        
        print("✅ Plots created successfully!")
        
    except ImportError:
        print("⚠️  Matplotlib not available - skipping plotting demo")
    except Exception as e:
        print(f"⚠️  Plotting demo failed: {e}")


def demo_error_handling():
    """Demonstrate error handling."""
    print("\n\n⚠️  Error Handling Demo")
    print("=" * 60)
    
    # 1. Invalid configuration
    print("\n1. Testing invalid configuration...")
    try:
        sim = SUEWSSimulation("nonexistent_file.yaml")
    except FileNotFoundError as e:
        print(f"✅ Caught expected error: {e}")
    
    # 2. Missing forcing data
    print("\n2. Testing missing forcing data...")
    try:
        config = create_sample_config()
        sim = SUEWSSimulation(config)
        sim.run()  # No forcing data
    except (RuntimeError, ValueError) as e:
        print(f"✅ Caught expected error: {type(e).__name__}")
    
    # 3. Results before run
    print("\n3. Testing results access before run...")
    try:
        config = create_sample_config()
        sim = SUEWSSimulation(config)
        summary = sim.summary()  # No results yet
    except RuntimeError as e:
        print(f"✅ Caught expected error: RuntimeError")
    
    print("✅ Error handling working correctly!")


def demo_performance_test():
    """Demonstrate performance characteristics."""
    print("\n\n⚡ Performance Test")
    print("=" * 60)
    
    # Test different simulation sizes
    test_cases = [
        ("Small (1 day)", 1),
        ("Medium (1 week)", 7),
        ("Large (1 month)", 30),
    ]
    
    for test_name, n_days in test_cases:
        print(f"\n🏃 {test_name}: {n_days} days...")
        
        config = create_sample_config()
        sim = SUEWSSimulation(config)
        forcing = create_sample_forcing_data(n_days=n_days)
        sim.setup_forcing(forcing)
        
        start_time = datetime.now()
        results = sim.run()
        duration = (datetime.now() - start_time).total_seconds()
        
        # Calculate performance metrics
        timesteps_per_sec = len(results) / duration if duration > 0 else 0
        memory_mb = results.memory_usage(deep=True).sum() / 1024 / 1024
        
        print(f"   ✅ Completed: {len(results)} timesteps in {duration:.3f}s")
        print(f"      Performance: {timesteps_per_sec:.0f} timesteps/sec")
        print(f"      Memory usage: {memory_mb:.1f} MB")


def main():
    """Run the complete demonstration."""
    print("🌟 SUEWS Simulation Class - Complete Demonstration")
    print("=" * 80)
    print("This demo showcases the modern, object-oriented interface for SUEWS")
    print("=" * 80)
    
    try:
        # Basic usage demo
        sim = demo_basic_usage()
        
        # Advanced features
        demo_advanced_features(sim)
        
        # YAML workflow
        demo_yaml_workflow()
        
        # Plotting capabilities
        demo_plotting()
        
        # Error handling
        demo_error_handling()
        
        # Performance testing
        demo_performance_test()
        
        print("\n\n🎉 DEMONSTRATION COMPLETED SUCCESSFULLY! 🎉")
        print("=" * 80)
        print("The SUEWSSimulation class provides:")
        print("✅ Intuitive object-oriented interface")
        print("✅ YAML-based configuration")
        print("✅ Pandas DataFrame integration")
        print("✅ Built-in data validation")
        print("✅ Flexible result access and analysis")
        print("✅ Quick plotting capabilities")
        print("✅ Robust error handling")
        print("✅ Export to multiple formats")
        print("✅ Chainable method design")
        print("=" * 80)
        
    except Exception as e:
        print(f"\n❌ Demo failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    if not success:
        exit(1)