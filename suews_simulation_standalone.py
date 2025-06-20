"""
SUEWS Simulation Class - Standalone Version

Modern, object-oriented interface for SUEWS urban climate model simulations.
This is a standalone version that doesn't depend on the existing SuPy infrastructure
during development and testing.

Author: Claude Code Integration
"""

from pathlib import Path
from typing import Union, Optional, Dict, Any, List, Tuple
import warnings
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import tempfile
import yaml


class SUEWSSimulation:
    """
    Modern SUEWS simulation class providing an intuitive interface for urban climate modelling.
    
    This class is designed to wrap the existing SuPy infrastructure with enhanced usability features:
    - Simple YAML-based configuration
    - Intelligent parameter overriding
    - Built-in result management and analysis
    - Chainable method design
    - Enhanced error handling
    
    Examples
    --------
    Basic usage:
    >>> sim = SUEWSSimulation.from_yaml("config.yaml")
    >>> results = sim.run()
    >>> sim.quick_plot(['QH', 'QE'])
    >>> sim.save("output.csv")
    
    Advanced usage with parameter overrides:
    >>> sim = SUEWSSimulation("config.yaml", 
    ...                       forcing_file="custom_forcing.txt",
    ...                       model_params={"tstep": 600})
    >>> results = sim.run(debug_mode=True, chunk_day=1800)
    >>> summary = sim.summary()
    """
    
    def __init__(
        self,
        config: Union[str, Path, Dict[str, Any]],
        forcing_file: Optional[Union[str, Path, pd.DataFrame]] = None,
        **kwargs
    ):
        """
        Initialise SUEWS simulation.
        
        Parameters
        ----------
        config : str, Path, or dict
            Configuration source:
            - Path to YAML configuration file
            - Dictionary with configuration parameters
        forcing_file : str, Path, or DataFrame, optional
            Meteorological forcing data source.
            If not provided, uses forcing_file from config.
        **kwargs
            Additional parameters to override configuration values.
        """
        self._config = None
        self._df_state_init = None
        self._df_forcing = None
        self._df_output = None
        self._df_state_final = None
        self._run_completed = False
        self._run_metadata = {}
        
        # Load and validate configuration
        self._load_config(config, **kwargs)
        
        # Set up forcing data if provided
        if forcing_file is not None:
            self.setup_forcing(forcing_file)
            
        self._log("SUEWSSimulation initialised successfully")
    
    def _log(self, message: str, level: str = "info"):
        """Simple logging method."""
        print(f"[{level.upper()}] {message}")
    
    def _load_config(self, config: Union[str, Path, Dict], **kwargs):
        """Load and validate configuration from various sources."""
        try:
            if isinstance(config, dict):
                self._config = config.copy()
            elif isinstance(config, (str, Path)):
                config_path = Path(config)
                if not config_path.exists():
                    raise FileNotFoundError(f"Configuration file not found: {config_path}")
                
                with open(config_path, 'r') as f:
                    self._config = yaml.safe_load(f)
            else:
                raise ValueError(f"Unsupported configuration type: {type(config)}")
            
            # Apply parameter overrides
            if kwargs:
                self._apply_parameter_overrides(**kwargs)
                
            # Convert to DataFrame state format
            self._df_state_init = self._config_to_dataframe_state()
            
        except Exception as e:
            self._log(f"Failed to load configuration: {e}", "error")
            raise
    
    def _apply_parameter_overrides(self, **kwargs):
        """Apply parameter overrides to configuration."""
        self._log(f"Applied {len(kwargs)} parameter overrides")
        # In a full implementation, this would update nested configuration parameters
        for key, value in kwargs.items():
            self._config[key] = value
    
    def _config_to_dataframe_state(self) -> pd.DataFrame:
        """Convert configuration to DataFrame state format."""
        self._log("Converting configuration to DataFrame state format")
        
        # Create a minimal valid state representation
        # In real implementation, this would use the SUEWS data model
        state_data = {
            'grid_id': [1],
            'lat': [self._config.get('sites', [{}])[0].get('properties', {}).get('lat', {}).get('value', 51.5)],
            'lon': [self._config.get('sites', [{}])[0].get('properties', {}).get('lng', {}).get('value', -0.1)],
            'tstep': [self._config.get('model', {}).get('control', {}).get('tstep', 300)]
        }
        
        return pd.DataFrame(state_data)
    
    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path], **kwargs) -> 'SUEWSSimulation':
        """
        Create SUEWSSimulation from YAML configuration file.
        
        Parameters
        ----------
        yaml_path : str or Path
            Path to YAML configuration file
        **kwargs
            Additional parameters to override configuration values
            
        Returns
        -------
        SUEWSSimulation
            Configured simulation instance
        """
        return cls(yaml_path, **kwargs)
    
    def setup_forcing(self, forcing_data: Union[str, Path, pd.DataFrame]):
        """
        Load and validate meteorological forcing data.
        
        Parameters
        ----------
        forcing_data : str, Path, or DataFrame
            Forcing data source
        """
        try:
            if isinstance(forcing_data, pd.DataFrame):
                self._df_forcing = forcing_data.copy()
            elif isinstance(forcing_data, (str, Path)):
                forcing_path = Path(forcing_data)
                if not forcing_path.exists():
                    raise FileNotFoundError(f"Forcing file not found: {forcing_path}")
                self._df_forcing = self._load_forcing_file(forcing_path)
            else:
                raise ValueError(f"Unsupported forcing data type: {type(forcing_data)}")
            
            # Validate forcing data
            self._validate_forcing()
            
            self._log(f"Loaded forcing data: {len(self._df_forcing)} timesteps")
            
        except Exception as e:
            self._log(f"Failed to setup forcing data: {e}", "error")
            raise
    
    def _load_forcing_file(self, forcing_path: Path) -> pd.DataFrame:
        """Load forcing data from file."""
        if forcing_path.suffix.lower() == '.csv':
            df = pd.read_csv(forcing_path, index_col=0, parse_dates=True)
        else:
            df = pd.read_csv(forcing_path, index_col=0, parse_dates=True)  # Simplified
        
        return df
    
    def _validate_forcing(self):
        """Validate forcing data format and content."""
        if self._df_forcing is None:
            raise ValueError("No forcing data loaded")
        
        # Check for required columns
        required_cols = ['Tair', 'RH', 'U', 'pres', 'rain', 'kdown']
        missing_cols = [col for col in required_cols if col not in self._df_forcing.columns]
        
        if missing_cols:
            self._log(f"Missing forcing columns: {missing_cols}", "warning")
    
    def validate(self) -> Dict[str, Any]:
        """
        Validate simulation configuration and data.
        
        Returns
        -------
        dict
            Validation report with status and any warnings/errors
        """
        validation_report = {
            'status': 'valid',
            'warnings': [],
            'errors': []
        }
        
        try:
            # Check configuration
            if self._config is None:
                validation_report['errors'].append("No configuration loaded")
            
            # Check state data
            if self._df_state_init is None:
                validation_report['errors'].append("No initial state data")
            
            # Check forcing data
            if self._df_forcing is None:
                validation_report['warnings'].append("No forcing data loaded")
            
            # Update status based on errors
            if validation_report['errors']:
                validation_report['status'] = 'invalid'
            elif validation_report['warnings']:
                validation_report['status'] = 'valid_with_warnings'
                
        except Exception as e:
            validation_report['status'] = 'error'
            validation_report['errors'].append(str(e))
        
        return validation_report
    
    def run(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        save_state: bool = False,
        chunk_day: int = 3660,
        debug_mode: bool = False,
        parallel: bool = False,
        **kwargs
    ) -> pd.DataFrame:
        """
        Execute SUEWS simulation.
        
        Parameters
        ----------
        start_date : datetime, optional
            Simulation start date
        end_date : datetime, optional
            Simulation end date
        save_state : bool, default False
            Save model states at each timestep
        chunk_day : int, default 3660
            Chunk size in days for memory management
        debug_mode : bool, default False
            Enable debug output
        parallel : bool, default False
            Use parallel execution
        **kwargs
            Additional parameters
            
        Returns
        -------
        pd.DataFrame
            Simulation results
        """
        # Validate before running
        validation = self.validate()
        if validation['status'] == 'invalid':
            raise RuntimeError(f"Simulation validation failed: {validation['errors']}")
        
        if validation['warnings']:
            for warning in validation['warnings']:
                warnings.warn(warning)
        
        try:
            # Prepare forcing data for date range
            df_forcing = self._prepare_forcing_for_run(start_date, end_date)
            
            # Record run metadata
            self._run_metadata = {
                'start_time': datetime.now(),
                'start_date': start_date,
                'end_date': end_date,
                'save_state': save_state,
                'chunk_day': chunk_day,
                'debug_mode': debug_mode,
                'parallel': parallel,
                'n_timesteps': len(df_forcing),
                'n_grids': len(self._df_state_init)
            }
            
            self._log(f"Starting SUEWS simulation: {self._run_metadata['n_timesteps']} timesteps, "
                     f"{self._run_metadata['n_grids']} grids")
            
            # Create mock results for demonstration
            # In real implementation, this would call the SuPy simulation functions
            self._df_output = self._create_mock_results(df_forcing)
            self._df_state_final = self._df_state_init.copy()
            
            # Update metadata
            self._run_metadata['end_time'] = datetime.now()
            self._run_metadata['duration'] = (
                self._run_metadata['end_time'] - self._run_metadata['start_time']
            ).total_seconds()
            
            self._run_completed = True
            
            self._log(f"SUEWS simulation completed successfully in "
                     f"{self._run_metadata['duration']:.2f} seconds")
            
            return self._df_output
            
        except Exception as e:
            self._log(f"SUEWS simulation failed: {e}", "error")
            raise RuntimeError(f"Simulation execution failed: {e}") from e
    
    def _prepare_forcing_for_run(
        self, 
        start_date: Optional[datetime], 
        end_date: Optional[datetime]
    ) -> pd.DataFrame:
        """Prepare forcing data for simulation run."""
        if self._df_forcing is None:
            raise ValueError("No forcing data available. Use setup_forcing() first.")
        
        df_forcing = self._df_forcing.copy()
        
        # Filter by date range if provided
        if start_date is not None:
            df_forcing = df_forcing[df_forcing.index >= start_date]
        if end_date is not None:
            df_forcing = df_forcing[df_forcing.index <= end_date]
        
        if len(df_forcing) == 0:
            raise ValueError("No forcing data available for specified date range")
        
        return df_forcing
    
    def _create_mock_results(self, df_forcing: pd.DataFrame) -> pd.DataFrame:
        """Create mock simulation results for testing."""
        n_timesteps = len(df_forcing)
        
        # Create realistic-looking SUEWS output variables
        data = {
            ('SUEWS', 'QH'): np.random.normal(50, 20, n_timesteps),  # Sensible heat flux
            ('SUEWS', 'QE'): np.random.normal(30, 15, n_timesteps),  # Latent heat flux  
            ('SUEWS', 'QS'): np.random.normal(-10, 30, n_timesteps), # Storage heat flux
            ('SUEWS', 'Tair'): np.random.normal(15, 8, n_timesteps), # Air temperature
            ('SUEWS', 'RH'): np.random.normal(70, 20, n_timesteps).clip(0, 100), # Relative humidity
            ('SUEWS', 'U10'): np.random.normal(3, 1.5, n_timesteps).clip(0), # Wind speed
            ('SUEWS', 'QN'): np.random.normal(70, 40, n_timesteps),  # Net radiation
            ('daily', 'Tmax'): np.random.normal(20, 10, n_timesteps), # Daily max temp
            ('daily', 'Tmin'): np.random.normal(8, 8, n_timesteps),   # Daily min temp
        }
        
        # Create multi-index columns
        columns = pd.MultiIndex.from_tuples(data.keys(), names=['group', 'variable'])
        
        # Create DataFrame with forcing data index
        df_results = pd.DataFrame(data, index=df_forcing.index, columns=columns)
        
        return df_results
    
    def get_results(
        self,
        variables: Optional[List[str]] = None,
        resample_freq: Optional[str] = None,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None
    ) -> pd.DataFrame:
        """
        Retrieve simulation results with optional filtering and resampling.
        
        Parameters
        ----------
        variables : list of str, optional
            Specific variables to return
        resample_freq : str, optional
            Pandas frequency string for resampling
        start_date : datetime, optional
            Start date for result filtering
        end_date : datetime, optional
            End date for result filtering
            
        Returns
        -------
        pd.DataFrame
            Filtered and optionally resampled results
        """
        if not self._run_completed:
            raise RuntimeError("No simulation results available. Run simulation first.")
        
        df_results = self._df_output.copy()
        
        # Filter by date range
        if start_date is not None:
            df_results = df_results[df_results.index >= start_date]
        if end_date is not None:
            df_results = df_results[df_results.index <= end_date]
        
        # Filter by variables
        if variables is not None:
            available_vars = df_results.columns.get_level_values(1).unique()
            missing_vars = [v for v in variables if v not in available_vars]
            if missing_vars:
                warnings.warn(f"Variables not found: {missing_vars}")
            
            # Select available variables
            valid_vars = [v for v in variables if v in available_vars]
            if valid_vars:
                df_results = df_results.loc[:, (slice(None), valid_vars)]
        
        # Resample if requested
        if resample_freq is not None:
            # Simple resampling - in real implementation would use SuPy's resample_output
            df_results = df_results.resample(resample_freq).mean()
        
        return df_results
    
    def summary(self) -> Dict[str, Any]:
        """
        Generate statistical summary of simulation results.
        
        Returns
        -------
        dict
            Statistical summary including means, extremes, and metadata
        """
        if not self._run_completed:
            raise RuntimeError("No simulation results available. Run simulation first.")
        
        # Basic statistics
        stats = self._df_output.describe()
        
        # Energy balance components
        energy_vars = ['QH', 'QE', 'QS']
        available_energy_vars = [v for v in energy_vars if ('SUEWS', v) in self._df_output.columns]
        
        energy_balance = {}
        if len(available_energy_vars) >= 2:
            energy_sums = {var: self._df_output[('SUEWS', var)].sum() for var in available_energy_vars}
            energy_balance['components'] = energy_sums
            energy_balance['total'] = sum(energy_sums.values())
        
        summary = {
            'run_metadata': self._run_metadata,
            'statistics': stats.to_dict(),
            'energy_balance': energy_balance,
            'data_coverage': {
                'n_timesteps': len(self._df_output),
                'start_date': self._df_output.index.min(),
                'end_date': self._df_output.index.max(),
                'missing_values': self._df_output.isnull().sum().to_dict()
            }
        }
        
        return summary
    
    def see(self, n_rows: int = 10) -> None:
        """
        Display preview of simulation results.
        
        Parameters
        ----------
        n_rows : int, default 10
            Number of rows to display
        """
        if not self._run_completed:
            print("No simulation results available. Run simulation first.")
            return
        
        print("SUEWS Simulation Results Preview")
        print("=" * 40)
        print(f"Total timesteps: {len(self._df_output)}")
        print(f"Date range: {self._df_output.index.min()} to {self._df_output.index.max()}")
        print(f"Variables: {self._df_output.columns.nlevels} groups, {len(self._df_output.columns)} total")
        print(f"\nFirst {min(n_rows, len(self._df_output))} rows:")
        print(self._df_output.head(n_rows))
    
    def quick_plot(
        self,
        variables: Optional[List[str]] = None,
        **plot_kwargs
    ):
        """
        Create quick visualisation of results.
        
        Parameters
        ----------
        variables : list of str, optional
            Variables to plot
        **plot_kwargs
            Additional arguments passed to plot function
        """
        if not self._run_completed:
            print("No simulation results available. Run simulation first.")
            return
        
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Matplotlib not available for plotting")
            return
        
        # Default variables for plotting
        if variables is None:
            variables = ['QH', 'QE', 'QS', 'Tair']
        
        # Get available variables
        all_vars = self._df_output.columns.get_level_values(1).unique()
        plot_vars = [v for v in variables if v in all_vars]
        
        if not plot_vars:
            print(f"None of the requested variables {variables} are available.")
            print(f"Available variables: {list(all_vars)}")
            return
        
        # Create subplots
        n_vars = len(plot_vars)
        fig, axes = plt.subplots(n_vars, 1, figsize=plot_kwargs.get('figsize', (12, 3*n_vars)))
        
        if n_vars == 1:
            axes = [axes]
        
        for i, var in enumerate(plot_vars):
            # Find the variable in the multi-index columns
            var_data = None
            for col in self._df_output.columns:
                if col[1] == var:
                    var_data = self._df_output[col]
                    break
            
            if var_data is not None:
                var_data.plot(ax=axes[i], title=f'{var}')
                axes[i].set_ylabel(var)
        
        plt.tight_layout()
        plt.show()
    
    def save(
        self,
        output_path: Union[str, Path],
        format: str = 'csv',
        **save_kwargs
    ) -> Path:
        """
        Save simulation results to file.
        
        Parameters
        ----------
        output_path : str or Path
            Output file path
        format : str, default 'csv'
            Output format: 'csv', 'excel', 'pickle', 'netcdf'
        **save_kwargs
            Additional arguments passed to save function
            
        Returns
        -------
        Path
            Path to saved file
        """
        if not self._run_completed:
            raise RuntimeError("No simulation results available. Run simulation first.")
        
        output_path = Path(output_path)
        
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if format.lower() == 'csv':
            self._df_output.to_csv(output_path, **save_kwargs)
        elif format.lower() == 'excel':
            self._df_output.to_excel(output_path, **save_kwargs)
        elif format.lower() == 'pickle':
            self._df_output.to_pickle(output_path, **save_kwargs)
        elif format.lower() == 'netcdf':
            # Convert to xarray and save as netCDF
            try:
                import xarray as xr
                ds = xr.Dataset.from_dataframe(self._df_output)
                ds.to_netcdf(output_path, **save_kwargs)
            except ImportError:
                raise ImportError("xarray required for netCDF export")
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        self._log(f"Results saved to: {output_path}")
        return output_path
    
    def clone(self) -> 'SUEWSSimulation':
        """
        Create a copy of the simulation for parameter studies.
        
        Returns
        -------
        SUEWSSimulation
            New simulation instance with same configuration
        """
        new_sim = SUEWSSimulation(self._config)
        if self._df_forcing is not None:
            new_sim._df_forcing = self._df_forcing.copy()
        return new_sim
    
    def reset(self):
        """
        Reset simulation to initial state, clearing results.
        """
        self._df_output = None
        self._df_state_final = None
        self._run_completed = False
        self._run_metadata = {}
        self._log("Simulation reset to initial state")
    
    @property
    def config(self) -> Dict[str, Any]:
        """Access to simulation configuration."""
        return self._config
    
    @property
    def results(self) -> Optional[pd.DataFrame]:
        """Access to simulation results."""
        return self._df_output
    
    @property
    def forcing(self) -> Optional[pd.DataFrame]:
        """Access to forcing data."""
        return self._df_forcing
    
    @property
    def is_complete(self) -> bool:
        """Check if simulation has been completed."""
        return self._run_completed
    
    @property
    def metadata(self) -> Dict[str, Any]:
        """Access to run metadata."""
        return self._run_metadata
    
    def __repr__(self) -> str:
        """String representation of simulation."""
        status = "completed" if self._run_completed else "not run"
        n_grids = len(self._df_state_init) if self._df_state_init is not None else 0
        forcing_info = f"{len(self._df_forcing)} timesteps" if self._df_forcing is not None else "no forcing"
        
        return (f"SUEWSSimulation(grids={n_grids}, forcing={forcing_info}, "
                f"status={status})")