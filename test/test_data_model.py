import unittest
import pandas as pd
import supy as sp
from pathlib import Path

class TestSUEWSConfig(unittest.TestCase):
    def test_config_conversion_cycle(self):
        """Test if SUEWS configuration can be correctly converted between YAML and DataFrame formats."""
        print("\n========================================")
        print("Testing YAML-DataFrame-YAML conversion cycle for SUEWS configuration...")

        # Get path to sample config from supy package
        path_sample_config = Path(sp.__file__).parent / "sample_run" / "sample_config.yml"

        # Load sample config from YAML
        config_orig = sp.data_model.SUEWSConfig.from_yaml(path_sample_config)

        # Convert to DataFrame
        df_state = config_orig.to_df_state()

        # Convert back to config object
        config_reconst = sp.data_model.SUEWSConfig.from_df_state(df_state)

        # Compare the two configs
        # We'll compare a few key attributes as a basic test
        self.assertEqual(config_orig.name, config_reconst.name)
        self.assertEqual(config_orig.description, config_reconst.description)
        self.assertEqual(config_orig.model.control.tstep, config_reconst.model.control.tstep)
        self.assertEqual(
            config_orig.model.physics.netradiationmethod.value,
            config_reconst.model.physics.netradiationmethod.value
        )
        self.assertEqual(
            config_orig.site[0].properties.lat.value,
            config_reconst.site[0].properties.lat.value
        )

        # Test if DataFrame conversion preserves structure
        df_state_2 = config_reconst.to_df_state()
        pd.testing.assert_frame_equal(df_state, df_state_2)

if __name__ == '__main__':
    unittest.main()

