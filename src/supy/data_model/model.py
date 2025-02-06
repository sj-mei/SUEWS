
import yaml
from typing import Optional
import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator
import pandas as pd
from enum import Enum

from .type import ValueWithDOI, Reference
from .type import init_df_state


class EmissionsMethod(Enum):
    # just a demo to show how to use Enum for emissionsmethod
    NO_EMISSIONS = 0
    CO2_ONLY = 1
    CO2_AND_ENERGY = 2
    CO2_AND_ENERGY_AND_VOC = 3
    CO2_AND_ENERGY_AND_VOC_AND_NOX = 4
    CO2_AND_ENERGY_AND_VOC_AND_NOX_AND_SO2 = 45
    CO2_AND_ENERGY_AND_VOC_AND_NOX_AND_SO2_AND_PM = 111

    def __int__(self):
        """Representation showing just the value"""
        return self.value

    def __repr__(self):
        """Representation showing the name and value"""
        return str(self.value)


def yaml_equivalent_of_default(dumper, data):
    return dumper.represent_scalar("tag:yaml.org,2002:int", str(data.value))


yaml.add_representer(EmissionsMethod, yaml_equivalent_of_default)


class ModelPhysics(BaseModel):
    netradiationmethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(3), description="Method used to calculate net radiation"
    )
    emissionsmethod: ValueWithDOI[EmissionsMethod] = Field(
        default=ValueWithDOI(EmissionsMethod.CO2_AND_ENERGY),
        description="Method used to calculate anthropogenic emissions",
    )
    storageheatmethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(1),
        description="Method used to calculate storage heat flux",
    )
    ohmincqf: ValueWithDOI[int] = Field(
        default=ValueWithDOI(0),
        description="Include anthropogenic heat in OHM calculations (1) or not (0)",
    )
    roughlenmommethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(2),
        description="Method used to calculate momentum roughness length",
    )
    roughlenheatmethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(2),
        description="Method used to calculate heat roughness length",
    )
    stabilitymethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(2),
        description="Method used for atmospheric stability calculation",
    )
    smdmethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(1),
        description="Method used to calculate soil moisture deficit",
    )
    waterusemethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(1), description="Method used to calculate water use"
    )
    diagmethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(1), description="Method used for model diagnostics"
    )
    faimethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(1),
        description="Method used to calculate frontal area index",
    )
    localclimatemethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(0),
        description="Method used for local climate zone calculations",
    )
    snowuse: ValueWithDOI[int] = Field(
        default=ValueWithDOI(0),
        description="Include snow calculations (1) or not (0)",
        enum=[0, 1],
    )
    stebbsmethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(0), description="Method used for stebbs calculations"
    )

    ref: Optional[Reference] = None

    @model_validator(mode="after")
    def check_all(self) -> "ModelPhysics":
        """
        Collect and aggregate all validation checks for ModelPhysics,
        so multiple errors (if any) can be raised together
        """
        errors = []
        # storageheatmethod check
        if self.storageheatmethod == 1 and self.ohmincqf != 0:
            errors.append(
                f"\nStorageHeatMethod is set to {self.storageheatmethod} and OhmIncQf is set to {self.ohmincqf}.\n"
                f"You should switch to OhmIncQf=0.\n"
            )
        elif self.storageheatmethod == 2 and self.ohmincqf != 1:
            errors.append(
                f"\nStorageHeatMethod is set to {self.storageheatmethod} and OhmIncQf is set to {self.ohmincqf}.\n"
                f"You should switch to OhmIncQf=1.\n"
            )

        # snowusemethod check
        if self.snowuse == 1:
            errors.append(
                f"\nSnowUse is set to {self.snowuse}.\n"
                f"There are no checks implemented for this case (snow calculations included in the run).\n"
                f"You should switch to SnowUse=0.\n"
            )

        # emissionsmethod check
        # if self.emissionsmethod == 45:
        #     errors.append(
        #         f"\nEmissionsMethod is set to {self.emissionsmethod}.\n"
        #         f"There are no checks implemented for this case (CO2 calculations included in the run).\n"
        #         f"You should switch to EmissionsMethod=0, 1, 2, 3, or 4.\n"
        #     )

        if errors:
            raise ValueError("\n".join(errors))

        return self

    # We then need to set to 0 (or None) all the CO2-related parameters or rules
    # in the code and return them accordingly in the yml file.

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert model physics properties to DataFrame state format."""
        df_state = init_df_state(grid_id)

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = "0"
            if (col_name, idx_str) not in df_state.columns:
                # df_state[(col_name, idx_str)] = np.nan
                df_state[(col_name, idx_str)] = None
            df_state.at[grid_id, (col_name, idx_str)] = int(value.value)

        list_attr = [
            "netradiationmethod",
            "emissionsmethod",
            "storageheatmethod",
            "ohmincqf",
            "roughlenmommethod",
            "roughlenheatmethod",
            "stabilitymethod",
            "smdmethod",
            "waterusemethod",
            "diagmethod",
            "faimethod",
            "localclimatemethod",
            "snowuse",
            "stebbsmethod",
        ]
        for attr in list_attr:
            set_df_value(attr, getattr(self, attr))
        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "ModelPhysics":
        """
        Reconstruct ModelPhysics from a DataFrame state format.

        Args:
            df: DataFrame containing model physics properties
            grid_id: Grid ID for the DataFrame index

        Returns:
            ModelPhysics: Instance of ModelPhysics
        """

        properties = {}

        list_attr = [
            "netradiationmethod",
            "emissionsmethod",
            "storageheatmethod",
            "ohmincqf",
            "roughlenmommethod",
            "roughlenheatmethod",
            "stabilitymethod",
            "smdmethod",
            "waterusemethod",
            "diagmethod",
            "faimethod",
            "localclimatemethod",
            "snowuse",
            "stebbsmethod",
        ]

        for attr in list_attr:
            try:
                properties[attr] = ValueWithDOI(int(df.loc[grid_id, (attr, "0")]))
            except KeyError:
                raise ValueError(f"Missing attribute '{attr}' in the DataFrame")

        return cls(**properties)


class ModelControl(BaseModel):
    tstep: int = Field(
        default=300, description="Time step in seconds for model calculations"
    )
    forcing_file: ValueWithDOI[str] = Field(
        default=ValueWithDOI("forcing.txt"),
        description="Path to meteorological forcing data file",
    )
    kdownzen: Optional[ValueWithDOI[int]] = Field(
        default=None,
        description="Use zenithal correction for downward shortwave radiation",
    )
    output_file: str = Field(
        default="output.txt", description="Path to model output file"
    )
    # daylightsaving_method: int
    diagnose: int = Field(
        default=0,
        description="Level of diagnostic output (0=none, 1=basic, 2=detailed)",
    )

    ref: Optional[Reference] = None

    @field_validator("tstep", "diagnose", mode="after")
    def validate_int_float(cls, v):
        if isinstance(v, (np.float64, np.float32)):
            return float(v)
        elif isinstance(v, (np.int64, np.int32)):
            return int(v)
        return v

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert model control properties to DataFrame state format."""
        df_state = init_df_state(grid_id)

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = "0"
            if (col_name, idx_str) not in df_state.columns:
                # df_state[(col_name, idx_str)] = np.nan
                df_state[(col_name, idx_str)] = None
            df_state.at[grid_id, (col_name, idx_str)] = value

        list_attr = ["tstep", "diagnose"]
        for attr in list_attr:
            set_df_value(attr, getattr(self, attr))
        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "ModelControl":
        """Reconstruct model control properties from DataFrame state format."""
        instance = cls()
        for attr in ["tstep", "diagnose"]:
            value = df.loc[grid_id, (attr, "0")]
            setattr(
                instance,
                attr,
                int(value) if isinstance(value, (np.int64, np.int32)) else value,
            )
        return instance


class Model(BaseModel):
    control: ModelControl = Field(
        default_factory=ModelControl,
        description="Model control parameters including timestep, output options, etc.",
    )
    physics: ModelPhysics = Field(
        default_factory=ModelPhysics,
        description="Model physics parameters including surface properties, coefficients, etc.",
    )

    @model_validator(mode="after")
    def validate_radiation_method(self) -> "Model":
        if (
            self.physics.netradiationmethod.value == 1
            and self.control.forcing_file.value == "forcing.txt"
        ):
            raise ValueError(
                "NetRadiationMethod is set to 1 (using observed Ldown). "
                "The sample forcing file lacks observed Ldown. Use netradiation = 3 for sample forcing. "
                "If not using sample forcing, ensure that the forcing file contains Ldown and rename from forcing.txt."
                # TODO: This is a temporary solution. We need to provide a better way to catch this.
            )
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert model to DataFrame state format"""
        df_state = init_df_state(grid_id)
        df_control = self.control.to_df_state(grid_id)
        df_physics = self.physics.to_df_state(grid_id)
        df_state = pd.concat([df_state, df_control, df_physics], axis=1)
        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "Model":
        """Reconstruct Model from DataFrame state format."""
        # Extract control and physics parameters
        control = ModelControl.from_df_state(df, grid_id)
        physics = ModelPhysics.from_df_state(df, grid_id)

        # Create an instance using the extracted parameters
        return cls(control=control, physics=physics)
