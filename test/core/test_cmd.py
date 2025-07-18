"""
Test command line interface tools.

This module tests the command line tools in supy.cmd,
including SUEWS simulation runner and table converter.
"""

import unittest
from unittest.mock import patch, MagicMock, mock_open, call
import warnings
from pathlib import Path
import tempfile
import shutil

import pandas as pd
import numpy as np
import click
from click.testing import CliRunner
import pytest

from supy.cmd.SUEWS import SUEWS
from supy.cmd.table_converter import convert_table_cmd

# Try to import to_yaml, but make it optional since it's not exposed in __init__.py
try:
    from supy.cmd.to_yaml import to_yaml
    HAS_TO_YAML = True
except ImportError:
    HAS_TO_YAML = False


class TestSUEWSCommand(unittest.TestCase):
    """Test SUEWS command line interface."""

    def setUp(self):
        """Set up test environment."""
        warnings.simplefilter("ignore", category=ImportWarning)
        self.runner = CliRunner()
        self.temp_dir = tempfile.mkdtemp()
        
        # Create a mock RunControl.nml file
        self.runcontrol_path = Path(self.temp_dir) / "RunControl.nml"
        self.runcontrol_content = """
&RunControl
    FileCode = 'Test'
    FileInputPath = './'
    FileOutputPath = './Output/'
    ResolutionFilesIn = 3600
    ResolutionFilesOut = 3600
    MultipleMetFiles = 0
/
"""
        self.runcontrol_path.write_text(self.runcontrol_content)

    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch('supy.cmd.SUEWS.init_supy')
    @patch('supy.cmd.SUEWS.load_forcing_grid')
    @patch('supy.cmd.SUEWS.run_supy')
    @patch('supy.cmd.SUEWS.save_supy')
    @patch('supy.cmd.SUEWS.load_SUEWS_nml_simple')
    def test_suews_single_grid_success(self, mock_load_nml, mock_save, mock_run, 
                                      mock_load_forcing, mock_init):
        """Test successful SUEWS run with single grid."""
        print("\n========================================")
        print("Testing SUEWS command with single grid...")

        # Set up mocks
        mock_df_state = pd.DataFrame({'test': [1]}, index=[1])  # Single grid
        mock_init.return_value = mock_df_state
        
        mock_load_nml.return_value = MagicMock(runcontrol=MagicMock(multiplemetfiles=0))
        
        mock_forcing = pd.DataFrame({
            'Tair': np.random.uniform(15, 25, 24),
            'RH': np.random.uniform(40, 80, 24)
        }, index=pd.date_range('2023-01-01', periods=24, freq='h'))
        mock_load_forcing.return_value = mock_forcing
        
        mock_output = pd.DataFrame({'QH': np.random.uniform(0, 200, 24)})
        mock_state_final = pd.DataFrame({'test': [2]})
        mock_run.return_value = (mock_output, mock_state_final)
        
        mock_save.return_value = ['output1.txt', 'output2.txt']

        # Run command
        result = self.runner.invoke(SUEWS, ['--path_runcontrol', str(self.runcontrol_path)])

        # Verify
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Running SUEWS via SuPy", result.output)
        self.assertIn("1 grids detected", result.output)
        self.assertIn("Same forcing conditions will be used for all grids", result.output)
        self.assertIn("SUEWS run successfully done!", result.output)

        # Check function calls
        mock_init.assert_called_once()
        mock_run.assert_called_once()
        mock_save.assert_called_once()

        print("✓ Single grid SUEWS run successful")

    @patch('supy.cmd.SUEWS.init_supy')
    @patch('supy.cmd.SUEWS.load_forcing_grid')
    @patch('supy.cmd.SUEWS.run_supy')
    @patch('supy.cmd.SUEWS.save_supy')
    @patch('supy.cmd.SUEWS.load_SUEWS_nml_simple')
    def test_suews_multi_grid_success(self, mock_load_nml, mock_save, mock_run, 
                                     mock_load_forcing, mock_init):
        """Test successful SUEWS run with multiple grids."""
        print("\n========================================")
        print("Testing SUEWS command with multiple grids...")

        # Set up mocks for multiple grids
        mock_df_state = pd.DataFrame({'test': [1, 2, 3]}, index=[1, 2, 3])  # 3 grids
        mock_init.return_value = mock_df_state
        
        mock_load_nml.return_value = MagicMock(runcontrol=MagicMock(multiplemetfiles=1))
        
        # Different forcing for each grid
        mock_forcings = []
        for i in range(3):
            forcing = pd.DataFrame({
                'Tair': np.random.uniform(15+i, 25+i, 24),
                'RH': np.random.uniform(40, 80, 24)
            }, index=pd.date_range(f'2023-01-0{i+1}', periods=24, freq='h'))
            mock_forcings.append(forcing)
        
        mock_load_forcing.side_effect = mock_forcings * 2  # Called twice per grid
        
        # Mock outputs for each grid
        mock_outputs = []
        mock_states = []
        for i in range(3):
            output = pd.DataFrame({'QH': np.random.uniform(0, 200, 24)})
            state = pd.DataFrame({'test': [i+10]})
            mock_outputs.append(output)
            mock_states.append(state)
        
        # Mock the parallel computation
        with patch('supy.cmd.SUEWS.db') as mock_db:
            mock_bag = MagicMock()
            mock_bag.map.return_value.compute.return_value = list(zip(mock_outputs, mock_states))
            mock_db.from_sequence.return_value = mock_bag
            
            mock_save.return_value = ['output1.txt', 'output2.txt']

            # Run command
            result = self.runner.invoke(SUEWS, ['--path_runcontrol', str(self.runcontrol_path)])

        # Verify
        self.assertEqual(result.exit_code, 0)
        self.assertIn("3 grids detected", result.output)
        self.assertIn("Grid-specific forcing conditions will be used", result.output)
        self.assertIn("SUEWS run successfully done!", result.output)

        print("✓ Multi-grid SUEWS run successful")

    def test_suews_missing_runcontrol(self):
        """Test SUEWS command with missing RunControl file."""
        print("\n========================================")
        print("Testing SUEWS command with missing RunControl...")

        # Run with non-existent file
        result = self.runner.invoke(SUEWS, ['--path_runcontrol', '/nonexistent/path.nml'])

        # Should fail
        self.assertNotEqual(result.exit_code, 0)

        print("✓ Missing RunControl properly handled")

    @patch('supy.cmd.SUEWS.init_supy')
    def test_suews_kernel_error(self, mock_init):
        """Test SUEWS command with kernel error."""
        print("\n========================================")
        print("Testing SUEWS command with kernel error...")

        # Set up mock to raise error
        mock_init.side_effect = RuntimeError("Kernel error")

        # Run command
        result = self.runner.invoke(SUEWS, ['--path_runcontrol', str(self.runcontrol_path)])

        # The SUEWS function calls sys.exit() on any error, which returns exit code 0
        # We can only check that it exits early without completing successfully
        self.assertNotIn("SUEWS run successfully done!", result.output)

        print("✓ Kernel error handled gracefully")


class TestTableConverterCommand(unittest.TestCase):
    """Test table converter command line interface."""

    def setUp(self):
        """Set up test environment."""
        warnings.simplefilter("ignore", category=ImportWarning)
        self.runner = CliRunner()
        self.temp_dir = tempfile.mkdtemp()
        self.input_dir = Path(self.temp_dir) / "input"
        self.output_dir = Path(self.temp_dir) / "output"
        self.input_dir.mkdir()
        
        # Create mock RunControl.nml
        (self.input_dir / "RunControl.nml").write_text("Mock RunControl")

    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch('supy.cmd.table_converter.convert_table')
    def test_convert_table_success(self, mock_convert):
        """Test successful table conversion."""
        print("\n========================================")
        print("Testing table converter command...")

        # Run command
        result = self.runner.invoke(convert_table_cmd, [
            '--from', '2019a',
            '--to', '2020a',
            '--input', str(self.input_dir),
            '--output', str(self.output_dir)
        ])

        # Verify
        self.assertEqual(result.exit_code, 0)
        mock_convert.assert_called_once_with(
            str(self.input_dir),
            str(self.output_dir),
            '2019a',
            '2020a'
        )

        print("✓ Table conversion successful")

    def test_convert_table_missing_input(self):
        """Test table converter with missing input directory."""
        print("\n========================================")
        print("Testing table converter with missing input...")

        # Run with non-existent directory
        result = self.runner.invoke(convert_table_cmd, [
            '--from', '2019a',
            '--to', '2020a',
            '--input', '/nonexistent/path',
            '--output', str(self.output_dir)
        ])

        # Should fail with missing directory
        self.assertNotEqual(result.exit_code, 0)

        print("✓ Missing input directory handled")

    def test_convert_table_missing_required_options(self):
        """Test table converter with missing required options."""
        print("\n========================================")
        print("Testing table converter with missing options...")

        # Run without required options
        result = self.runner.invoke(convert_table_cmd, [
            '--input', str(self.input_dir)
        ])

        # Should fail
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("Missing option", result.output)

        print("✓ Missing required options handled")


@pytest.mark.skipif(not HAS_TO_YAML, reason="to_yaml module not available")
class TestToYamlCommand(unittest.TestCase):
    """Test YAML converter command line interface."""

    def setUp(self):
        """Set up test environment."""
        warnings.simplefilter("ignore", category=ImportWarning)
        self.runner = CliRunner()
        self.temp_dir = tempfile.mkdtemp()
        self.input_dir = Path(self.temp_dir) / "input"
        self.input_dir.mkdir()
        
        # Create mock RunControl.nml
        (self.input_dir / "RunControl.nml").write_text("Mock RunControl")

    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch('supy.cmd.to_yaml.load_InitialCond_grid_df')
    @patch('supy.cmd.to_yaml.SUEWSConfig')
    def test_to_yaml_without_conversion(self, mock_config_class, mock_load):
        """Test YAML conversion without version conversion."""
        print("\n========================================")
        print("Testing to_yaml command without version conversion...")

        # Set up mocks
        mock_df_state = pd.DataFrame({'test': [1]})
        mock_load.return_value = mock_df_state
        
        mock_config = MagicMock()
        mock_config_class.from_df_state.return_value = mock_config

        # Run command
        output_file = self.temp_dir / "output.yaml"
        result = self.runner.invoke(to_yaml, [
            '--input-dir', str(self.input_dir),
            '--output-file', str(output_file)
        ])

        # Verify
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Successfully converted", result.output)
        
        mock_load.assert_called_once()
        mock_config_class.from_df_state.assert_called_once_with(mock_df_state)
        mock_config.to_yaml.assert_called_once()

        print("✓ YAML conversion without version conversion successful")

    @patch('supy.cmd.to_yaml.convert_table')
    @patch('supy.cmd.to_yaml.load_InitialCond_grid_df')
    @patch('supy.cmd.to_yaml.SUEWSConfig')
    def test_to_yaml_with_conversion(self, mock_config_class, mock_load, mock_convert):
        """Test YAML conversion with version conversion."""
        print("\n========================================")
        print("Testing to_yaml command with version conversion...")

        # Set up mocks
        mock_df_state = pd.DataFrame({'test': [1]})
        mock_load.return_value = mock_df_state
        
        mock_config = MagicMock()
        mock_config_class.from_df_state.return_value = mock_config

        # Run command with version conversion
        output_file = self.temp_dir / "output.yaml"
        result = self.runner.invoke(to_yaml, [
            '--input-dir', str(self.input_dir),
            '--output-file', str(output_file),
            '--from-ver', '2019a'
        ])

        # Verify
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Converting tables from version 2019a", result.output)
        self.assertIn("Successfully converted", result.output)
        
        # Should call convert_table
        mock_convert.assert_called_once()
        mock_config.to_yaml.assert_called_once()

        print("✓ YAML conversion with version conversion successful")

    def test_to_yaml_missing_runcontrol(self):
        """Test YAML conversion with missing RunControl.nml."""
        print("\n========================================")
        print("Testing to_yaml with missing RunControl...")

        # Remove RunControl.nml
        (self.input_dir / "RunControl.nml").unlink()

        # Run command
        output_file = self.temp_dir / "output.yaml"
        result = self.runner.invoke(to_yaml, [
            '--input-dir', str(self.input_dir),
            '--output-file', str(output_file)
        ])

        # Should fail
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("RunControl.nml not found", result.output)

        print("✓ Missing RunControl.nml handled properly")


class TestCommandLineHelp(unittest.TestCase):
    """Test command line help messages."""

    def setUp(self):
        """Set up test environment."""
        self.runner = CliRunner()

    def test_suews_help(self):
        """Test SUEWS command help."""
        print("\n========================================")
        print("Testing SUEWS command help...")

        result = self.runner.invoke(SUEWS, ['--help'])
        
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Run SUEWS simulation", result.output)
        self.assertIn("--path_runcontrol", result.output)

        print("✓ SUEWS help displayed correctly")

    def test_convert_table_help(self):
        """Test table converter command help."""
        print("\n========================================")
        print("Testing table converter help...")

        result = self.runner.invoke(convert_table_cmd, ['--help'])
        
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Convert SUEWS input tables", result.output)
        self.assertIn("--from", result.output)
        self.assertIn("--to", result.output)

        print("✓ Table converter help displayed correctly")

    @pytest.mark.skipif(not HAS_TO_YAML, reason="to_yaml module not available")
    def test_to_yaml_help(self):
        """Test YAML converter command help."""
        print("\n========================================")
        print("Testing to_yaml help...")

        result = self.runner.invoke(to_yaml, ['--help'])
        
        self.assertEqual(result.exit_code, 0)
        self.assertIn("Convert SUEWS table-based input to a YAML configuration", result.output)
        self.assertIn("--input-dir", result.output)
        self.assertIn("--output-file", result.output)

        print("✓ to_yaml help displayed correctly")


if __name__ == "__main__":
    unittest.main()