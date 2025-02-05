"""Core model configuration and control classes."""

from typing import Dict, List, Optional, Union, Literal
import numpy as np
import pandas as pd
import yaml
from pydantic import (
    BaseModel,
    Field,
    model_validator,
    field_validator,
)

from .base import ValueWithDOI, Reference
from .physics import ModelPhysics
from .surface.base import init_df_state


class ModelControl(BaseModel):
    """Model control parameters including timestep, output options, etc."""
    tstep: int = Field(
        default=300,
        description="Time step in seconds for model calculations",
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
        default="output.txt",
        description="Path to model output file",
    )
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
            setattr(instance, attr, int(value) if isinstance(value, (np.int64, np.int32)) else value)
        return instance


class Model(BaseModel):
    """Complete model configuration including control and physics settings."""
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
        if self.physics.netradiationmethod.value == 1 and self.control.forcing_file.value == "forcing.txt":
            raise ValueError(
                "NetRadiationMethod is set to 1 (using observed Ldown). "
                "The sample forcing file lacks observed Ldown. Use netradiation = 3 for sample forcing. "
                "If not using sample forcing, ensure that the forcing file contains Ldown and rename from forcing.txt."
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
