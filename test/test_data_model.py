import unittest
import pandas as pd
import supy as sp
from pathlib import Path
import numpy as np
from supy.data_model import (
    SUEWSConfig,
    Model,
    Site,
    SiteProperties,
    InitialStates,
    AnthropogenicEmissions,
    AnthropogenicHeat,
    CO2Params,
    ValueWithDOI,
)


class TestSUEWSConfig(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.path_sample_config = Path(sp.__file__).parent / "sample_run" / "sample_config.yml"
        self.config = SUEWSConfig.from_yaml(self.path_sample_config)

    def test_config_conversion_cycle(self):
        """Test if SUEWS configuration can be correctly converted between YAML and DataFrame formats."""
        print("\n========================================")
        print("Testing YAML-DataFrame-YAML conversion cycle for SUEWS configuration...")

        # Convert to DataFrame
        df_state = self.config.to_df_state()

        # Convert back to config object
        config_reconst = SUEWSConfig.from_df_state(df_state)

        # Compare the two configs
        # We'll compare a few key attributes as a basic test
        self.assertEqual(self.config.name, config_reconst.name)
        self.assertEqual(self.config.description, config_reconst.description)
        self.assertEqual(self.config.model.control.tstep, config_reconst.model.control.tstep)
        self.assertEqual(
            self.config.model.physics.netradiationmethod.value,
            config_reconst.model.physics.netradiationmethod.value
        )
        self.assertEqual(
            self.config.site[0].properties.lat.value,
            config_reconst.site[0].properties.lat.value
        )

        # Test if DataFrame conversion preserves structure
        df_state_2 = config_reconst.to_df_state()
        pd.testing.assert_frame_equal(df_state, df_state_2)

    def test_anthropogenic_emissions(self):
        """Test anthropogenic emissions data model."""
        # Create test data
        emissions = AnthropogenicEmissions(
            startdls=ValueWithDOI(80.0),
            enddls=ValueWithDOI(300.0),
            heat=AnthropogenicHeat(
                popdensnighttime=100.0,
            ),
            co2=CO2Params(
                co2pointsource=ValueWithDOI(0.5),
                ef_umolco2perj=ValueWithDOI(0.1),
            )
        )

        # Convert to DataFrame
        df_state = emissions.to_df_state(grid_id=0)

        # Convert back to object
        emissions_reconst = AnthropogenicEmissions.from_df_state(df_state, grid_id=0)

        # Test key values
        self.assertEqual(emissions.startdls.value, emissions_reconst.startdls.value)
        self.assertEqual(emissions.enddls.value, emissions_reconst.enddls.value)
        self.assertEqual(emissions.heat.popdensnighttime, emissions_reconst.heat.popdensnighttime)
        self.assertEqual(emissions.co2.co2pointsource.value, emissions_reconst.co2.co2pointsource.value)
        self.assertEqual(emissions.co2.ef_umolco2perj.value, emissions_reconst.co2.ef_umolco2perj.value)

    def test_model_physics_validation(self):
        """Test model physics validation rules."""
        model = Model()

        # Test storageheatmethod and ohmincqf validation
        with self.assertRaises(ValueError):
            model.physics.storageheatmethod = ValueWithDOI(1)
            model.physics.ohmincqf = ValueWithDOI(1)
            model.physics.model_validate(model.physics)

        with self.assertRaises(ValueError):
            model.physics.storageheatmethod = ValueWithDOI(2)
            model.physics.ohmincqf = ValueWithDOI(0)
            model.physics.model_validate(model.physics)

    def test_site_properties(self):
        """Test site properties data model."""
        site = Site()

        # Test latitude bounds
        with self.assertRaises(ValueError):
            site.properties.lat = ValueWithDOI(91.0)  # Invalid latitude
            site.properties.model_validate(site.properties)

        with self.assertRaises(ValueError):
            site.properties.lat = ValueWithDOI(-91.0)  # Invalid latitude
            site.properties.model_validate(site.properties)

        # Test valid latitude
        site.properties.lat = ValueWithDOI(51.5)  # London's latitude
        site.properties.model_validate(site.properties)

    def test_initial_states(self):
        """Test initial states data model."""
        states = InitialStates()

        # Test snow albedo bounds
        with self.assertRaises(ValueError):
            states.snowalb = ValueWithDOI(1.5)  # Invalid albedo > 1
            states.model_validate(states)

        with self.assertRaises(ValueError):
            states.snowalb = ValueWithDOI(-0.1)  # Invalid albedo < 0
            states.model_validate(states)

        # Test valid snow albedo
        states.snowalb = ValueWithDOI(0.8)
        states.model_validate(states)

    def test_multi_site_config(self):
        """Test configuration with multiple sites."""
        config = SUEWSConfig(
            name="Multi-site test",
            description="Test configuration with multiple sites",
            site=[
                Site(gridiv=0),
                Site(gridiv=1),
                Site(gridiv=2),
            ]
        )

        # Convert to DataFrame
        df_state = config.to_df_state()

        # Check if all sites are present
        self.assertEqual(len(df_state.index), 3)
        self.assertTrue(all(idx in df_state.index for idx in [0, 1, 2]))

        # Convert back and check if sites are preserved
        config_reconst = SUEWSConfig.from_df_state(df_state)
        self.assertEqual(len(config_reconst.site), 3)
        self.assertEqual([site.gridiv for site in config_reconst.site], [0, 1, 2])

if __name__ == '__main__':
    unittest.main()

