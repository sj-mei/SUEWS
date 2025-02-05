"""Water surface models for SUEWS."""

from typing import Dict, List, Optional, Union, Literal
import pandas as pd
from pydantic import BaseModel, Field, model_validator, PrivateAttr

from ..base import ValueWithDOI, Reference
from .base import SurfaceType, init_df_state
from .urban import NonVegetatedSurfaceProperties


class WaterProperties(NonVegetatedSurfaceProperties):
    """Properties for water surfaces."""
    _surface_type: Literal[SurfaceType.WATER] = SurfaceType.WATER
    waterdepth: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.5),
        description="Water depth [m]",
    )
    watervolume: ValueWithDOI[float] = Field(
        default=ValueWithDOI(1000.0),
        description="Water volume [m3]",
    )
    waterdensity: ValueWithDOI[float] = Field(
        default=ValueWithDOI(1000.0),
        description="Water density [kg m-3]",
    )
    wateralbedo: ValueWithDOI[float] = Field(
        ge=0, le=1,
        default=ValueWithDOI(0.1),
        description="Water surface albedo",
    )
    wateremis: ValueWithDOI[float] = Field(
        ge=0, le=1,
        default=ValueWithDOI(0.95),
        description="Water surface emissivity",
    )
    waterroughness: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.001),
        description="Water surface roughness length [m]",
    )
    waterheatcapacity: ValueWithDOI[float] = Field(
        default=ValueWithDOI(4.18e6),
        description="Water heat capacity [J m-3 K-1]",
    )
    watertransmissivity: ValueWithDOI[float] = Field(
        ge=0, le=1,
        default=ValueWithDOI(0.2),
        description="Water transmissivity",
    )
    waterextinction: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.3),
        description="Water extinction coefficient [m-1]",
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert water surface properties to DataFrame state format."""
        df_state = super().to_df_state(grid_id)

        # Get surface index for water (6)
        surf_idx = 6

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = None
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Add water-specific properties
        for attr in [
            "waterdepth",
            "watervolume",
            "waterdensity",
            "wateralbedo",
            "wateremis",
            "waterroughness",
            "waterheatcapacity",
            "watertransmissivity",
            "waterextinction",
        ]:
            set_df_value(attr, getattr(self, attr).value)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "WaterProperties":
        """Reconstruct water surface properties from DataFrame state format."""
        surf_idx = 6
        instance = super().from_df_state(df, grid_id, surf_idx)

        # Helper function to get values from DataFrame
        def get_df_value(col_name: str) -> float:
            return df.loc[grid_id, (col_name, f"({surf_idx},)")]

        # Set water-specific properties
        for attr in [
            "waterdepth",
            "watervolume",
            "waterdensity",
            "wateralbedo",
            "wateremis",
            "waterroughness",
            "waterheatcapacity",
            "watertransmissivity",
            "waterextinction",
        ]:
            try:
                value = get_df_value(attr)
                setattr(instance, attr, ValueWithDOI(value))
            except KeyError:
                continue

        return instance
