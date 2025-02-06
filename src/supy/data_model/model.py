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





class NetRadiationMethod(Enum):
    OBSERVED_LDOWN = 1  # Using observed Ldown
    OBSERVED_LDOWN_LOUT = 2  # Using observed Ldown and Lout
    MODELLED = 3  # Using modelled radiation components

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class StorageHeatMethod(Enum):
    OHM_WITHOUT_QF = 1  # OHM without anthropogenic heat
    OHM_WITH_QF = 2  # OHM with anthropogenic heat

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)

class OhmIncQf(Enum):
    INCLUDE = 1  # Include anthropogenic heat in OHM calculations
    EXCLUDE = 0  # Exclude anthropogenic heat in OHM calculations

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class RoughnessMethod(Enum):
    FIXED = 1  # Fixed roughness length
    VARIABLE = 2  # Variable roughness length based on vegetation state

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class StabilityMethod(Enum):
    NEUTRAL = 1  # Neutral stability
    VARIABLE = 2  # Variable stability based on atmospheric conditions

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class SMDMethod(Enum):
    BASIC = 1  # Basic soil moisture deficit calculation
    ADVANCED = 2  # Advanced soil moisture deficit calculation with more processes

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class WaterUseMethod(Enum):
    NONE = 0  # No water use calculation
    BASIC = 1  # Basic water use calculation
    ADVANCED = 2  # Advanced water use calculation with irrigation

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class DiagMethod(Enum):
    NONE = 0  # No diagnostics
    BASIC = 1  # Basic diagnostics
    DETAILED = 2  # Detailed diagnostics

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class FAIMethod(Enum):
    FIXED = 1  # Fixed frontal area index
    VARIABLE = 2  # Variable frontal area index based on vegetation state

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class LocalClimateMethod(Enum):
    NONE = 0  # No local climate zone calculations
    BASIC = 1  # Basic local climate zone calculations
    DETAILED = 2  # Detailed local climate zone calculations

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


class StebbsMethod(Enum):
    NONE = 0  # No STEBBS calculations
    BASIC = 1  # Basic STEBBS calculations
    DETAILED = 2  # Detailed STEBBS calculations

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)

class SnowUse(Enum):
    ENABLED = 1  # Using snow calculations
    DISABLED = 0  # Not using snow calculations

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)


def yaml_equivalent_of_default(dumper, data):
    """Convert enum values to YAML scalar integers.

    This function is used as a YAML representer for enum classes. It converts the enum value
    to a string and represents it as a YAML integer scalar.

    Args:
        dumper: The YAML dumper instance
        data: The enum value to be converted

    Returns:
        A YAML scalar node containing the integer value of the enum
    """
    return dumper.represent_scalar("tag:yaml.org,2002:int", str(data.value))

# Register YAML representers for all enums
for enum_class in [
    NetRadiationMethod,
    EmissionsMethod,
    StorageHeatMethod,
    RoughnessMethod,
    StabilityMethod,
    SMDMethod,
    WaterUseMethod,
    DiagMethod,
    FAIMethod,
    LocalClimateMethod,
    StebbsMethod,
    SnowUse,
    OhmIncQf,
]:
    yaml.add_representer(enum_class, yaml_equivalent_of_default)


class ModelPhysics(BaseModel):
    netradiationmethod: ValueWithDOI[NetRadiationMethod] = Field(
        default=ValueWithDOI(NetRadiationMethod.MODELLED),
        description="Method used to calculate net radiation",
    )
    emissionsmethod: ValueWithDOI[EmissionsMethod] = Field(
        default=ValueWithDOI(EmissionsMethod.CO2_AND_ENERGY),
        description="Method used to calculate anthropogenic emissions",
    )
    storageheatmethod: ValueWithDOI[StorageHeatMethod] = Field(
        default=ValueWithDOI(StorageHeatMethod.OHM_WITHOUT_QF),
        description="Method used to calculate storage heat flux",
    )
    ohmincqf: ValueWithDOI[OhmIncQf] = Field(
        default=ValueWithDOI(OhmIncQf.EXCLUDE),
        description="Include anthropogenic heat in OHM calculations (1) or not (0)",
    )
    roughlenmommethod: ValueWithDOI[RoughnessMethod] = Field(
        default=ValueWithDOI(RoughnessMethod.VARIABLE),
        description="Method used to calculate momentum roughness length",
    )
    roughlenheatmethod: ValueWithDOI[RoughnessMethod] = Field(
        default=ValueWithDOI(RoughnessMethod.VARIABLE),
        description="Method used to calculate heat roughness length",
    )
    stabilitymethod: ValueWithDOI[StabilityMethod] = Field(
        default=ValueWithDOI(StabilityMethod.VARIABLE),
        description="Method used for atmospheric stability calculation",
    )
    smdmethod: ValueWithDOI[SMDMethod] = Field(
        default=ValueWithDOI(SMDMethod.BASIC),
        description="Method used to calculate soil moisture deficit",
    )
    waterusemethod: ValueWithDOI[WaterUseMethod] = Field(
        default=ValueWithDOI(WaterUseMethod.BASIC),
        description="Method used to calculate water use",
    )
    diagmethod: ValueWithDOI[DiagMethod] = Field(
        default=ValueWithDOI(DiagMethod.BASIC),
        description="Method used for model diagnostics",
    )
    faimethod: ValueWithDOI[FAIMethod] = Field(
        default=ValueWithDOI(FAIMethod.FIXED),
        description="Method used to calculate frontal area index",
    )
    localclimatemethod: ValueWithDOI[LocalClimateMethod] = Field(
        default=ValueWithDOI(LocalClimateMethod.NONE),
        description="Method used for local climate zone calculations",
    )
    snowuse: ValueWithDOI[SnowUse] = Field(
        default=ValueWithDOI(SnowUse.DISABLED),
        description="Include snow calculations (1) or not (0)",
    )
    stebbsmethod: ValueWithDOI[StebbsMethod] = Field(
        default=ValueWithDOI(StebbsMethod.NONE),
        description="Method used for stebbs calculations",
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
