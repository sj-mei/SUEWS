"""
Test supy post-processing functionality.

This module tests all post-processing functions in supy, including
output resampling, aggregation, and data manipulation utilities.
"""

import tempfile
import warnings
from pathlib import Path
from unittest import TestCase

import numpy as np
import pandas as pd
import pytest

import supy as sp
from supy._post import resample_output, dict_var_aggm


class TestResampleOutput(TestCase):
    """Test output resampling functionality."""
    
    def setUp(self):
        """Set up test environment."""
        warnings.simplefilter("ignore", category=ImportWarning)
        
        # Create sample output data for testing
        self.df_state_init, self.df_forcing = sp.load_SampleData()
        
        # Run a short simulation to get output
        df_forcing_short = self.df_forcing.iloc[:288*7]  # One week
        self.df_output, self.df_state = sp.run_supy(
            df_forcing_short, 
            self.df_state_init,
            check_input=False
        )
    
    def test_resample_output_hourly(self):
        """Test resampling output to hourly frequency."""
        print("\n========================================")
        print("Testing resample_output to hourly...")
        
        # Resample to hourly
        df_hourly = resample_output(self.df_output, freq="60T")
        
        # Validation
        self.assertIsInstance(df_hourly, pd.DataFrame)
        self.assertEqual(df_hourly.index.freq.n, 60)  # 60 minutes
        
        # Check that we have fewer rows (aggregated)
        self.assertLess(len(df_hourly), len(self.df_output))
        
        # Check that columns are preserved
        original_groups = self.df_output.columns.get_level_values('group').unique()
        resampled_groups = df_hourly.columns.get_level_values('group').unique()
        for group in original_groups:
            if group != 'DailyState':  # DailyState handled differently
                self.assertIn(group, resampled_groups)
        
        print(f"✓ Resampled from {len(self.df_output)} to {len(df_hourly)} hourly records")
    
    def test_resample_output_daily(self):
        """Test resampling output to daily frequency."""
        print("\n========================================")
        print("Testing resample_output to daily...")
        
        # Resample to daily
        df_daily = resample_output(self.df_output, freq="D")
        
        # Validation
        self.assertIsInstance(df_daily, pd.DataFrame)
        self.assertEqual(len(df_daily), 7)  # 7 days
        
        # Check aggregation methods
        # Energy fluxes should be averaged
        qn_hourly = self.df_output.SUEWS['QN'].resample('D').mean()
        qn_daily = df_daily.SUEWS['QN']
        pd.testing.assert_series_equal(qn_hourly, qn_daily, check_names=False)
        
        # Rain should be summed
        rain_hourly = self.df_output.SUEWS['Rain'].resample('D').sum()
        rain_daily = df_daily.SUEWS['Rain']
        pd.testing.assert_series_equal(rain_hourly, rain_daily, check_names=False)
        
        print(f"✓ Resampled to {len(df_daily)} daily records with correct aggregation")
    
    def test_resample_output_custom_freq(self):
        """Test resampling with custom frequency."""
        print("\n========================================")
        print("Testing resample_output with custom frequency...")
        
        # Resample to 30 minutes
        df_30min = resample_output(self.df_output, freq="30T")
        
        # Should have twice as many records as hourly
        df_hourly = resample_output(self.df_output, freq="60T")
        self.assertAlmostEqual(len(df_30min), len(df_hourly) * 2, delta=2)
        
        print(f"✓ Custom frequency resampling works correctly")
    
    def test_resample_output_custom_aggm(self):
        """Test resampling with custom aggregation methods."""
        print("\n========================================")
        print("Testing resample_output with custom aggregation...")
        
        # Create custom aggregation dictionary
        custom_aggm = {
            'QN': 'max',  # Maximum instead of mean
            'Rain': 'mean',  # Mean instead of sum
        }
        
        # Resample with custom aggregation
        df_custom = resample_output(self.df_output, freq="D", dict_aggm=custom_aggm)
        
        # Verify custom aggregation
        qn_max = self.df_output.SUEWS['QN'].resample('D').max()
        qn_custom = df_custom.SUEWS['QN']
        pd.testing.assert_series_equal(qn_max, qn_custom, check_names=False)
        
        print("✓ Custom aggregation methods work correctly")
    
    def test_resample_preserves_multiindex(self):
        """Test that resampling preserves MultiIndex column structure."""
        print("\n========================================")
        print("Testing MultiIndex preservation in resampling...")
        
        # Resample
        df_resampled = resample_output(self.df_output, freq="60T")
        
        # Check MultiIndex structure
        self.assertIsInstance(df_resampled.columns, pd.MultiIndex)
        self.assertEqual(df_resampled.columns.nlevels, 2)
        self.assertIn('group', df_resampled.columns.names)
        self.assertIn('var', df_resampled.columns.names)
        
        print("✓ MultiIndex structure preserved correctly")


class TestAggregationMethods(TestCase):
    """Test aggregation method definitions."""
    
    def test_dict_var_aggm_completeness(self):
        """Test that aggregation dictionary covers common variables."""
        print("\n========================================")
        print("Testing aggregation method dictionary...")
        
        # Check essential variables have aggregation methods
        essential_vars = [
            'QN', 'QF', 'QS', 'QE', 'QH',  # Energy fluxes
            'Rain', 'Evap', 'RO',  # Water fluxes
            'Tair', 'RH', 'Pres', 'U',  # Meteorological
            'Kdown', 'Kup', 'Ldown', 'Lup',  # Radiation
        ]
        
        for var in essential_vars:
            self.assertIn(var, dict_var_aggm, f"Missing aggregation method for {var}")
        
        # Check aggregation methods are valid
        valid_methods = ['mean', 'sum', 'max', 'min', 'first', 'last']
        for var, method in dict_var_aggm.items():
            self.assertIn(method, valid_methods, f"Invalid method '{method}' for {var}")
        
        print(f"✓ Aggregation dictionary contains {len(dict_var_aggm)} variables")
    
    def test_aggregation_method_logic(self):
        """Test that aggregation methods make physical sense."""
        print("\n========================================")
        print("Testing aggregation method logic...")
        
        # Cumulative variables should be summed
        cumulative_vars = ['Rain', 'Irr', 'Evap', 'RO', 'ROSoil', 'ROPipe', 'ROWater', 'ROPav', 'ROVeg']
        for var in cumulative_vars:
            if var in dict_var_aggm:
                self.assertEqual(dict_var_aggm[var], 'sum', f"{var} should be summed")
        
        # Instantaneous variables should be averaged
        instant_vars = ['Tair', 'RH', 'Pres', 'U', 'QN', 'QF', 'QS', 'QE', 'QH']
        for var in instant_vars:
            if var in dict_var_aggm:
                self.assertEqual(dict_var_aggm[var], 'mean', f"{var} should be averaged")
        
        print("✓ Aggregation methods follow physical logic")


class TestPostProcessingUtilities(TestCase):
    """Test other post-processing utilities."""
    
    def setUp(self):
        """Set up test environment."""
        warnings.simplefilter("ignore", category=ImportWarning)
        
        # Run a minimal simulation
        df_state_init, df_forcing = sp.load_SampleData()
        df_forcing_short = df_forcing.iloc[:288*2]  # Two days
        self.df_output, self.df_state = sp.run_supy(
            df_forcing_short,
            df_state_init,
            check_input=False
        )
    
    def test_output_groups(self):
        """Test output group structure."""
        print("\n========================================")
        print("Testing output group structure...")
        
        # Get unique groups
        groups = self.df_output.columns.get_level_values('group').unique()
        
        # Check essential groups exist
        essential_groups = ['SUEWS', 'DailyState']
        for group in essential_groups:
            self.assertIn(group, groups, f"Missing output group: {group}")
        
        # Access groups using attribute notation
        self.assertTrue(hasattr(self.df_output, 'SUEWS'))
        self.assertIsInstance(self.df_output.SUEWS, pd.DataFrame)
        
        print(f"✓ Output contains {len(groups)} groups")
    
    def test_dailystate_extraction(self):
        """Test DailyState data extraction."""
        print("\n========================================")
        print("Testing DailyState extraction...")
        
        # Extract DailyState using xs()
        df_daily = self.df_output.xs('DailyState', level='group', axis=1)
        
        # Check that we have daily data
        mask_has_data = df_daily.notna().any(axis=1)
        n_days = mask_has_data.sum()
        
        # Should have approximately 2 days of data
        self.assertGreaterEqual(n_days, 1)
        self.assertLessEqual(n_days, 3)
        
        # Check some DailyState variables exist
        daily_vars = ['HDD1_h', 'HDD2_c']
        for var in daily_vars:
            if var in df_daily.columns:
                self.assertTrue(df_daily[var].notna().any())
        
        print(f"✓ Extracted {n_days} days of DailyState data")
    
    def test_energy_balance(self):
        """Test energy balance calculations from output."""
        print("\n========================================")
        print("Testing energy balance calculations...")
        
        # Calculate energy balance residual
        df_suews = self.df_output.SUEWS
        residual = df_suews['QN'] - (df_suews['QH'] + df_suews['QE'] + df_suews['QS'])
        
        # Check that residual is small (within 5% of QN)
        relative_error = (residual.abs() / df_suews['QN'].abs().replace(0, np.nan)).mean()
        self.assertLess(relative_error, 0.05, "Energy balance residual too large")
        
        print(f"✓ Energy balance residual: {relative_error*100:.2f}% of QN")
    
    def test_water_balance(self):
        """Test water balance calculations from output."""
        print("\n========================================")
        print("Testing water balance calculations...")
        
        df_suews = self.df_output.SUEWS
        
        # Water balance: Rain + Irr = Evap + RO + ΔStorage
        water_in = df_suews['Rain'] + df_suews['Irr']
        water_out = df_suews['Evap'] + df_suews['RO']
        storage_change = df_suews['State'].diff()
        
        # Calculate balance (ignoring first timestep due to diff)
        balance = water_in[1:] - water_out[1:] - storage_change[1:]
        
        # Check balance is close to zero
        self.assertLess(balance.abs().mean(), 0.01, "Water balance error too large")
        
        print(f"✓ Water balance error: {balance.abs().mean():.6f} mm")


class TestMultiGridPostProcessing(TestCase):
    """Test post-processing for multi-grid simulations."""
    
    def setUp(self):
        """Set up test environment."""
        warnings.simplefilter("ignore", category=ImportWarning)
        
        # Create multi-grid simulation
        df_state_single, df_forcing = sp.load_SampleData()
        
        # Duplicate for 3 grids
        n_grids = 3
        df_state_multi = pd.concat([df_state_single for _ in range(n_grids)])
        df_state_multi.index = pd.RangeIndex(n_grids, name='grid')
        
        # Run short simulation
        df_forcing_short = df_forcing.iloc[:288]  # One day
        self.df_output, self.df_state = sp.run_supy(
            df_forcing_short,
            df_state_multi,
            check_input=False
        )
        self.n_grids = n_grids
    
    def test_multigrid_resample(self):
        """Test resampling multi-grid output."""
        print("\n========================================")
        print("Testing multi-grid resampling...")
        
        # Resample to hourly
        df_hourly = resample_output(self.df_output, freq="60T")
        
        # Check that grid structure is preserved
        self.assertEqual(df_hourly.index.nlevels, 2)  # datetime and grid
        self.assertIn('grid', df_hourly.index.names)
        
        # Check we have data for all grids
        grids = df_hourly.index.get_level_values('grid').unique()
        self.assertEqual(len(grids), self.n_grids)
        
        print(f"✓ Resampled {self.n_grids} grids correctly")
    
    def test_grid_aggregation(self):
        """Test aggregating across grids."""
        print("\n========================================")
        print("Testing grid aggregation...")
        
        # Get mean values across all grids
        df_suews = self.df_output.SUEWS
        
        # Group by time and calculate mean across grids
        df_mean = df_suews.groupby(level='datetime').mean()
        
        # Validate
        self.assertEqual(len(df_mean), 288)  # One day of 5-min data
        self.assertFalse(df_mean.index.names[0])  # Single index, not MultiIndex
        
        # Check that mean is reasonable
        for grid in range(self.n_grids):
            grid_data = df_suews.xs(grid, level='grid')
            # Mean should be between min and max of individual grids
            self.assertTrue((df_mean >= grid_data.min()).all().all())
            self.assertTrue((df_mean <= grid_data.max()).all().all())
        
        print("✓ Grid aggregation works correctly")


class TestErrorHandling(TestCase):
    """Test error handling in post-processing."""
    
    def test_invalid_frequency(self):
        """Test handling of invalid resampling frequency."""
        print("\n========================================")
        print("Testing error handling for invalid frequency...")
        
        # Create minimal output
        df_state, df_forcing = sp.load_SampleData()
        df_output, _ = sp.run_supy(df_forcing.iloc[:288], df_state)
        
        # Test invalid frequency
        with self.assertRaises((ValueError, KeyError)):
            resample_output(df_output, freq="invalid")
        
        print("✓ Error handling works correctly")
    
    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        print("\n========================================")
        print("Testing error handling for empty DataFrame...")
        
        # Create empty DataFrame with correct structure
        df_empty = pd.DataFrame(
            columns=pd.MultiIndex.from_tuples([('SUEWS', 'QN'), ('SUEWS', 'QH')]),
            index=pd.DatetimeIndex([])
        )
        
        # Should handle gracefully
        df_result = resample_output(df_empty, freq="60T")
        self.assertTrue(df_result.empty)
        
        print("✓ Empty DataFrame handled correctly")


if __name__ == "__main__":
    import unittest
    unittest.main()