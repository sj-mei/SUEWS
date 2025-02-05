"""Physics-related model settings and enumerations."""

from enum import Enum
from typing import Optional
import pandas as pd
from pydantic import BaseModel, Field, model_validator

from .base import ValueWithDOI, Reference
from .surface.base import init_df_state


class EmissionsMethod(Enum):
    """Methods for calculating anthropogenic emissions."""
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


class ModelPhysics(BaseModel):
    """Physics scheme settings for SUEWS model."""

    netradiationmethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(3),
        description="Method used to calculate net radiation",
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
        default=ValueWithDOI(1),
        description="Method used to calculate water use",
    )
    diagmethod: ValueWithDOI[int] = Field(
        default=ValueWithDOI(1),
        description="Method used for model diagnostics",
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
        default=ValueWithDOI(0),
        description="Method used for STEBBS calculations",
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
        if self.storageheatmethod.value == 1 and self.ohmincqf.value != 0:
            errors.append(
                f"\nStorageHeatMethod is set to {self.storageheatmethod.value} and OhmIncQf is set to {self.ohmincqf.value}.\n"
                f"You should switch to OhmIncQf=0.\n"
            )
        elif self.storageheatmethod.value == 2 and self.ohmincqf.value != 1:
            errors.append(
                f"\nStorageHeatMethod is set to {self.storageheatmethod.value} and OhmIncQf is set to {self.ohmincqf.value}.\n"
                f"You should switch to OhmIncQf=1.\n"
            )

        # snowusemethod check
        if self.snowuse.value == 1:
            errors.append(
                f"\nSnowUse is set to {self.snowuse.value}.\n"
                f"There are no checks implemented for this case (snow calculations included in the run).\n"
                f"You should switch to SnowUse=0.\n"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert model physics properties to DataFrame state format."""
        df_state = init_df_state(grid_id)

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = "0"
            if (col_name, idx_str) not in df_state.columns:
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
                value = int(df.loc[grid_id, (attr, "0")])
                if attr == "emissionsmethod":
                    value = EmissionsMethod(value)
                properties[attr] = ValueWithDOI(value)
            except KeyError:
                raise ValueError(f"Missing attribute '{attr}' in the DataFrame")

        return cls(**properties)
