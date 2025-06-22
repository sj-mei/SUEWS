"""
Comprehensive test suite for SUEWSSimulation class.

Tests cover:
- Basic functionality
- Configuration management  
- Simulation execution
- Result handling
- Error conditions
- Integration with existing SuPy infrastructure

Author: Claude Code Integration
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import tempfile
import shutil
from unittest.mock import Mock, patch, MagicMock

# Import the class under test
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import supy as sp
from supy.suews_sim import SUEWSSimulation
from supy.data_model import SUEWSConfig


class TestSUEWSSimulationInit:
    """Test SUEWSSimulation initialisation and configuration."""
    
    def test_init_from_yaml_config(self, sample_yaml_config):
        """Test initialisation from YAML configuration file."""
        sim = SUEWSSimulation(sample_yaml_config)
        
        assert sim.config is not None
        assert sim.is_complete is False
        assert sim.results is None
        assert isinstance(sim.metadata, dict)
    
    def test_init_from_suews_config_object(self, sample_suews_config):
        """Test initialisation from SUEWSConfig object."""
        sim = SUEWSSimulation(sample_suews_config)
        
        assert sim.config is sample_suews_config
        assert sim.is_complete is False
    
    def test_init_from_dict_config(self, sample_config_dict):
        """Test initialisation from dictionary configuration."""
        sim = SUEWSSimulation(sample_config_dict)
        
        assert sim.config is not None
        assert sim.is_complete is False
    
    def test_init_with_forcing_file(self, sample_yaml_config, sample_forcing_data):
        """Test initialisation with forcing file provided."""
        sim = SUEWSSimulation(sample_yaml_config, forcing_file=sample_forcing_data)
        
        assert sim.forcing is not None
        assert len(sim.forcing) > 0
    
    def test_init_with_parameter_overrides(self, sample_yaml_config):
        """Test initialisation with parameter overrides."""
        sim = SUEWSSimulation(
            sample_yaml_config,
            tstep=600,
            debug_mode=True,
            output_dir="/tmp/test"
        )
        
        assert sim.config is not None
    
    def test_init_with_invalid_config_path(self):
        """Test error handling for invalid configuration path."""
        with pytest.raises(FileNotFoundError):
            SUEWSSimulation("nonexistent_config.yaml")
    
    def test_init_with_invalid_config_type(self):
        """Test error handling for invalid configuration type."""
        with pytest.raises(ValueError):
            SUEWSSimulation(123)  # Invalid type
    
    def test_from_yaml_classmethod(self, sample_yaml_config):
        """Test from_yaml class method."""
        sim = SUEWSSimulation.from_yaml(sample_yaml_config)
        
        assert isinstance(sim, SUEWSSimulation)
        assert sim.config is not None


class TestSUEWSSimulationForcing:
    """Test forcing data management."""
    
    def test_setup_forcing_from_dataframe(self, sample_suews_config, sample_forcing_dataframe):
        """Test setting up forcing from DataFrame."""
        sim = SUEWSSimulation(sample_suews_config)
        sim.setup_forcing(sample_forcing_dataframe)
        
        assert sim.forcing is not None
        assert len(sim.forcing) == len(sample_forcing_dataframe)
        pd.testing.assert_frame_equal(sim.forcing, sample_forcing_dataframe)
    
    def test_setup_forcing_from_file(self, sample_suews_config, sample_forcing_file):
        """Test setting up forcing from file."""
        sim = SUEWSSimulation(sample_suews_config)
        sim.setup_forcing(sample_forcing_file)
        
        assert sim.forcing is not None
        assert len(sim.forcing) > 0
    
    def test_setup_forcing_invalid_file(self, sample_suews_config):
        """Test error handling for invalid forcing file."""
        sim = SUEWSSimulation(sample_suews_config)
        
        with pytest.raises(FileNotFoundError):
            sim.setup_forcing("nonexistent_forcing.csv")
    
    def test_setup_forcing_invalid_type(self, sample_suews_config):
        """Test error handling for invalid forcing data type."""
        sim = SUEWSSimulation(sample_suews_config)
        
        with pytest.raises(ValueError):
            sim.setup_forcing(123)  # Invalid type
    
    def test_validate_forcing_missing_columns(self, sample_suews_config):
        """Test validation with missing required columns."""
        sim = SUEWSSimulation(sample_suews_config)
        
        # Create forcing with missing columns
        incomplete_forcing = pd.DataFrame({
            'Tair': [20.0, 21.0],
            'RH': [60.0, 65.0]
            # Missing other required columns
        }, index=pd.date_range('2012-01-01', periods=2, freq='h'))
        
        # Should warn but not fail
        sim.setup_forcing(incomplete_forcing)
        assert sim.forcing is not None


class TestSUEWSSimulationValidation:
    """Test simulation validation functionality."""
    
    def test_validate_complete_setup(self, configured_simulation):
        """Test validation of complete simulation setup."""
        validation = configured_simulation.validate()
        
        assert validation['status'] in ['valid', 'valid_with_warnings']
        assert isinstance(validation['warnings'], list)
        assert isinstance(validation['errors'], list)
    
    def test_validate_missing_config(self):
        """Test validation with missing configuration."""
        sim = SUEWSSimulation.__new__(SUEWSSimulation)  # Create without __init__
        sim._config = None
        sim._df_state_init = None
        sim._df_forcing = None
        
        validation = sim.validate()
        
        assert validation['status'] == 'invalid'
        assert 'No configuration loaded' in validation['errors']
    
    def test_validate_missing_forcing(self, sample_suews_config):
        """Test validation with missing forcing data."""
        sim = SUEWSSimulation(sample_suews_config)
        
        validation = sim.validate()
        
        assert validation['status'] == 'invalid'
        assert 'forcing' in str(validation['errors']).lower()


class TestSUEWSSimulationExecution:
    """Test simulation execution methods."""
    
    @patch('supy._run.run_supy_ser')
    def test_run_basic_execution(self, mock_run_supy_ser, configured_simulation):
        """Test basic simulation execution."""
        # Mock successful simulation
        mock_output = pd.DataFrame(
            np.random.randn(48, 5),
            columns=pd.MultiIndex.from_product([['SUEWS'], ['QH', 'QE', 'QS', 'Tair', 'RH']]),
            index=pd.date_range('2012-01-01', periods=48, freq='30min')
        )
        mock_state = pd.DataFrame({'dummy': [1]})
        mock_run_supy_ser.return_value = (mock_output, mock_state, None, None)
        
        results = configured_simulation.run()
        
        assert configured_simulation.is_complete is True
        assert configured_simulation.results is not None
        assert len(results) == 48
        mock_run_supy_ser.assert_called_once()
    
    @patch('supy._run.run_supy_par')
    def test_run_parallel_execution(self, mock_run_supy_par, configured_simulation):
        """Test parallel simulation execution."""
        # Mock successful parallel simulation
        mock_output = pd.DataFrame(
            np.random.randn(48, 5),
            columns=pd.MultiIndex.from_product([['SUEWS'], ['QH', 'QE', 'QS', 'Tair', 'RH']]),
            index=pd.date_range('2012-01-01', periods=48, freq='30min')
        )
        mock_state = pd.DataFrame({'dummy': [1]})
        mock_run_supy_par.return_value = (mock_output, mock_state)
        
        # Add multi-grid setup
        configured_simulation._df_state_init = pd.DataFrame({'dummy': [1, 2]})  # 2 grids
        
        results = configured_simulation.run(parallel=True)
        
        assert configured_simulation.is_complete is True
        mock_run_supy_par.assert_called_once()
    
    def test_run_with_date_range(self, configured_simulation):
        """Test simulation with specific date range."""
        with patch('supy._run.run_supy_ser') as mock_run:
            mock_run.return_value = (pd.DataFrame(), pd.DataFrame(), None, None)
            
            # Use dates that exist in our sample forcing data (48 timesteps from 2012-01-01)
            start_date = datetime(2012, 1, 1, 1, 0)  # Start from second timestep
            end_date = datetime(2012, 1, 1, 2, 0)    # End at third timestep
            
            configured_simulation.run(start_date=start_date, end_date=end_date)
            
            # Verify filtering was applied
            call_args = mock_run.call_args[0]
            forcing_used = call_args[0]
            assert forcing_used.index.min() >= start_date
            assert forcing_used.index.max() <= end_date
    
    def test_run_with_options(self, configured_simulation):
        """Test simulation with various options."""
        with patch('supy._run.run_supy_ser') as mock_run:
            mock_run.return_value = (pd.DataFrame(), pd.DataFrame(), None, None)
            
            configured_simulation.run(
                save_state=True,
                debug_mode=True,
                chunk_day=1800
            )
            
            # Verify options were passed
            call_args = mock_run.call_args
            assert call_args[0][2] is True  # save_state
            assert call_args[0][3] == 1800  # chunk_day
            assert call_args[0][4] is True  # debug_mode
    
    def test_run_validation_failure(self, sample_suews_config):
        """Test run with validation failure."""
        sim = SUEWSSimulation(sample_suews_config)
        # Don't setup forcing to cause validation failure
        
        with pytest.raises(RuntimeError, match="validation failed"):
            sim.run()
    
    def test_run_simulation_failure(self, configured_simulation):
        """Test handling of simulation execution failure."""
        with patch('supy._run.run_supy_ser') as mock_run:
            mock_run.side_effect = Exception("Simulation kernel error")
            
            with pytest.raises(RuntimeError, match="Simulation execution failed"):
                configured_simulation.run()
    
    def test_run_metadata_tracking(self, configured_simulation):
        """Test run metadata is properly tracked."""
        with patch('supy._run.run_supy_ser') as mock_run:
            mock_run.return_value = (pd.DataFrame(), pd.DataFrame(), None, None)
            
            start_time = datetime.now()
            configured_simulation.run(debug_mode=True, chunk_day=1000)
            
            metadata = configured_simulation.metadata
            assert 'start_time' in metadata
            assert 'end_time' in metadata
            assert 'duration' in metadata
            assert metadata['debug_mode'] is True
            assert metadata['chunk_day'] == 1000
            assert metadata['start_time'] >= start_time


class TestSUEWSSimulationResults:
    """Test result handling and analysis."""
    
    def test_get_results_all_variables(self, simulation_with_results):
        """Test getting all results."""
        results = simulation_with_results.get_results()
        
        assert isinstance(results, pd.DataFrame)
        assert len(results) > 0
        assert results.columns.nlevels == 2  # Multi-index columns
    
    def test_get_results_specific_variables(self, simulation_with_results):
        """Test getting specific variables."""
        variables = ['QH', 'QE']
        results = simulation_with_results.get_results(variables=variables)
        
        # Check that only requested variables are present
        result_vars = results.columns.get_level_values(1).unique()
        assert all(var in result_vars for var in variables)
    
    def test_get_results_with_resampling(self, simulation_with_results):
        """Test getting results with resampling."""
        results = simulation_with_results.get_results(resample_freq='H')
        
        assert isinstance(results, pd.DataFrame)
        # Should have fewer rows due to resampling
        original_length = len(simulation_with_results.results)
        assert len(results) <= original_length
    
    def test_get_results_date_filtering(self, simulation_with_results):
        """Test getting results with date filtering."""
        start_date = simulation_with_results.results.index[10]
        end_date = simulation_with_results.results.index[20]
        
        results = simulation_with_results.get_results(
            start_date=start_date,
            end_date=end_date
        )
        
        assert len(results) <= 11  # 10 to 20 inclusive
        assert results.index.min() >= start_date
        assert results.index.max() <= end_date
    
    def test_get_results_before_run(self, configured_simulation):
        """Test error when getting results before running simulation."""
        with pytest.raises(RuntimeError, match="No simulation results available"):
            configured_simulation.get_results()
    
    def test_summary_generation(self, simulation_with_results):
        """Test summary generation."""
        summary = simulation_with_results.summary()
        
        assert isinstance(summary, dict)
        assert 'run_metadata' in summary
        assert 'statistics' in summary
        assert 'data_coverage' in summary
        assert 'n_timesteps' in summary['data_coverage']
    
    def test_summary_before_run(self, configured_simulation):
        """Test error when generating summary before running simulation."""
        with pytest.raises(RuntimeError, match="No simulation results available"):
            configured_simulation.summary()
    
    def test_see_method(self, simulation_with_results, capsys):
        """Test see method for result preview."""
        simulation_with_results.see(n_rows=5)
        
        captured = capsys.readouterr()
        assert "SUEWS Simulation Results Preview" in captured.out
        assert "Total timesteps:" in captured.out
    
    def test_see_before_run(self, configured_simulation, capsys):
        """Test see method before running simulation."""
        configured_simulation.see()
        
        captured = capsys.readouterr()
        assert "No simulation results available" in captured.out


class TestSUEWSSimulationPlotting:
    """Test plotting and visualisation functionality."""
    
    @patch('pandas.Series.plot')
    @patch('matplotlib.pyplot.tight_layout')
    @patch('matplotlib.pyplot.show')
    @patch('matplotlib.pyplot.subplots')
    def test_quick_plot_default_variables(self, mock_subplots, mock_show, mock_tight_layout, mock_plot, simulation_with_results):
        """Test quick plot with default variables."""
        # Mock matplotlib
        mock_fig = Mock()
        mock_axes = [Mock() for _ in range(4)]
        mock_subplots.return_value = (mock_fig, mock_axes)
        mock_plot.return_value = Mock()
        
        simulation_with_results.quick_plot()
        
        mock_subplots.assert_called_once()
        mock_show.assert_called_once()
    
    @patch('pandas.Series.plot')
    @patch('matplotlib.pyplot.tight_layout')
    @patch('matplotlib.pyplot.show')
    @patch('matplotlib.pyplot.subplots')
    def test_quick_plot_specific_variables(self, mock_subplots, mock_show, mock_tight_layout, mock_plot, simulation_with_results):
        """Test quick plot with specific variables."""
        mock_fig = Mock()
        mock_axes = [Mock(), Mock()]
        mock_subplots.return_value = (mock_fig, mock_axes)
        mock_plot.return_value = Mock()
        
        simulation_with_results.quick_plot(['QH', 'QE'])
        
        mock_subplots.assert_called_once()
    
    def test_quick_plot_before_run(self, configured_simulation, capsys):
        """Test quick plot before running simulation."""
        configured_simulation.quick_plot()
        
        captured = capsys.readouterr()
        assert "No simulation results available" in captured.out


class TestSUEWSSimulationSaveLoad:
    """Test save/load functionality."""
    
    def test_save_csv(self, simulation_with_results, tmp_path):
        """Test saving results to CSV."""
        output_path = tmp_path / "results.csv"
        
        saved_path = simulation_with_results.save(output_path, format='csv')
        
        assert saved_path.exists()
        assert saved_path == output_path
        
        # Verify file can be read back
        loaded_df = pd.read_csv(saved_path, index_col=0, header=[0, 1])
        assert len(loaded_df) > 0
    
    def test_save_excel(self, simulation_with_results, tmp_path):
        """Test saving results to Excel."""
        output_path = tmp_path / "results.xlsx"
        
        saved_path = simulation_with_results.save(output_path, format='excel')
        
        assert saved_path.exists()
    
    def test_save_pickle(self, simulation_with_results, tmp_path):
        """Test saving results to pickle."""
        output_path = tmp_path / "results.pkl"
        
        saved_path = simulation_with_results.save(output_path, format='pickle')
        
        assert saved_path.exists()
        
        # Verify file can be read back
        loaded_df = pd.read_pickle(saved_path)
        pd.testing.assert_frame_equal(loaded_df, simulation_with_results.results)
    
    @patch('xarray.Dataset.from_dataframe')
    def test_save_netcdf(self, mock_from_dataframe, simulation_with_results, tmp_path):
        """Test saving results to netCDF."""
        output_path = tmp_path / "results.nc"
        
        # Mock xarray functionality
        mock_dataset = Mock()
        mock_from_dataframe.return_value = mock_dataset
        
        saved_path = simulation_with_results.save(output_path, format='netcdf')
        
        assert saved_path == output_path
        mock_dataset.to_netcdf.assert_called_once()
    
    def test_save_unsupported_format(self, simulation_with_results, tmp_path):
        """Test error handling for unsupported format."""
        output_path = tmp_path / "results.txt"
        
        with pytest.raises(ValueError, match="Unsupported format"):
            simulation_with_results.save(output_path, format='txt')
    
    def test_save_before_run(self, configured_simulation, tmp_path):
        """Test error when saving before running simulation."""
        output_path = tmp_path / "results.csv"
        
        with pytest.raises(RuntimeError, match="No simulation results available"):
            configured_simulation.save(output_path)


class TestSUEWSSimulationUtility:
    """Test utility methods."""
    
    def test_clone_simulation(self, configured_simulation):
        """Test cloning simulation."""
        clone = configured_simulation.clone()
        
        assert isinstance(clone, SUEWSSimulation)
        assert clone is not configured_simulation
        assert clone.config is configured_simulation.config
        assert clone.is_complete is False
    
    def test_reset_simulation(self, simulation_with_results):
        """Test resetting simulation."""
        # Verify initial state
        assert simulation_with_results.is_complete is True
        assert simulation_with_results.results is not None
        
        # Reset
        simulation_with_results.reset()
        
        # Verify reset state
        assert simulation_with_results.is_complete is False
        assert simulation_with_results.results is None
        assert len(simulation_with_results.metadata) == 0
    
    def test_property_access(self, simulation_with_results):
        """Test property access."""
        assert simulation_with_results.config is not None
        assert simulation_with_results.results is not None
        assert simulation_with_results.forcing is not None
        assert simulation_with_results.is_complete is True
        assert isinstance(simulation_with_results.metadata, dict)
    
    def test_string_representation(self, configured_simulation):
        """Test string representation."""
        repr_str = repr(configured_simulation)
        
        assert "SUEWSSimulation" in repr_str
        assert "grids=" in repr_str
        assert "status=" in repr_str


class TestSUEWSSimulationEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_empty_forcing_data(self, sample_suews_config):
        """Test handling of empty forcing data."""
        sim = SUEWSSimulation(sample_suews_config)
        
        empty_forcing = pd.DataFrame(columns=['Tair', 'RH', 'U'])
        sim.setup_forcing(empty_forcing)
        
        with pytest.raises(ValueError, match="No forcing data available"):
            sim.run()
    
    def test_invalid_date_range(self, configured_simulation):
        """Test handling of invalid date range."""
        future_start = datetime(2050, 1, 1)
        future_end = datetime(2050, 12, 31)
        
        with pytest.raises(ValueError, match="No forcing data available"):
            configured_simulation.run(start_date=future_start, end_date=future_end)
    
    def test_missing_required_variables(self, simulation_with_results):
        """Test handling of missing variables in get_results."""
        results = simulation_with_results.get_results(['nonexistent_var'])
        
        # Should return empty DataFrame or handle gracefully
        assert isinstance(results, pd.DataFrame)


# Test Fixtures

# Get SuPy module directory
supy_module_dir = Path(sp.__file__).parent

@pytest.fixture
def sample_yaml_config():
    """Use the actual sample YAML configuration file from SuPy."""
    return supy_module_dir / "sample_run" / "sample_config.yml"


@pytest.fixture
def sample_supy_data():
    """Load real sample data from SuPy."""
    df_state_init, df_forcing = sp.load_SampleData()
    return df_state_init, df_forcing


@pytest.fixture
def sample_suews_config():
    """Create a sample SUEWSConfig object using real SuPy initialization."""
    from supy.data_model import init_config_from_yaml
    config_path = supy_module_dir / "sample_run" / "sample_config.yml"
    return init_config_from_yaml(config_path)


@pytest.fixture
def sample_config_dict():
    """Create a sample configuration dictionary."""
    return {
        'name': 'test config',
        'model': {'control': {'tstep': 300}},
        'sites': [{'name': 'test site', 'gridiv': 1}]
    }


@pytest.fixture
def sample_forcing_dataframe(sample_supy_data):
    """Get sample forcing data from SuPy."""
    _, df_forcing = sample_supy_data
    # Return first 48 timesteps for quick testing
    return df_forcing.iloc[:48]


@pytest.fixture
def sample_forcing_file(sample_forcing_dataframe, tmp_path):
    """Create sample forcing data file."""
    forcing_file = tmp_path / "forcing.csv"
    sample_forcing_dataframe.to_csv(forcing_file)
    return forcing_file


@pytest.fixture
def sample_forcing_data(sample_forcing_dataframe):
    """Alias for sample_forcing_dataframe."""
    return sample_forcing_dataframe


@pytest.fixture
def configured_simulation(sample_yaml_config, sample_forcing_dataframe):
    """Create a configured simulation ready for testing."""
    sim = SUEWSSimulation(sample_yaml_config)
    sim.setup_forcing(sample_forcing_dataframe)
    return sim


@pytest.fixture
def simulation_with_results(configured_simulation):
    """Create a simulation with mock results."""
    # Create mock results matching SuPy output structure
    dates = pd.date_range('2012-01-01', periods=48, freq='5min')
    
    # Create multi-level columns matching SuPy output structure
    # SuPy uses (group, variable) structure where group is like 'SUEWS'
    columns = pd.MultiIndex.from_product([
        ['SUEWS'],  # group
        ['QH', 'QE', 'QS', 'Tair', 'RH']  # variables
    ], names=['group', 'var'])
    
    # Create multi-level index with grid level for resample_output compatibility
    index = pd.MultiIndex.from_product([
        [1],  # grid
        dates  # datetime
    ], names=['grid', 'datetime'])
    
    # Create DataFrame with proper structure
    data = np.random.randn(48, 5)
    results = pd.DataFrame(
        data,
        columns=columns,
        index=index
    )
    
    # Mock the simulation as completed
    configured_simulation._df_output = results
    configured_simulation._run_completed = True
    configured_simulation._run_metadata = {
        'start_time': datetime.now(),
        'end_time': datetime.now(),
        'duration': 10.5,
        'n_timesteps': 48,
        'n_grids': 1
    }
    
    return configured_simulation


if __name__ == "__main__":
    pytest.main([__file__])