"""State management for SUEWS model."""

from typing import Dict, List, Optional, Union, Literal
import pandas as pd
from pydantic import BaseModel, Field, model_validator, PrivateAttr

from .base import ValueWithDOI, Reference
from .surface.base import SurfaceType, init_df_state


class SurfaceState(BaseModel):
    """Surface state variables for each surface type."""
    surface_type: SurfaceType = Field(
        description="Type of surface this state belongs to",
    )
    state: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="Surface water state [mm]",
    )
    gdd: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="Growing degree days [°C]",
    )
    sdd: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="Senescence degree days [°C]",
    )
    lai: ValueWithDOI[float] = Field(
        default=ValueWithDOI(1.0),
        description="Leaf area index [m² m⁻²]",
    )
    albedo: ValueWithDOI[float] = Field(
        ge=0, le=1,
        default=ValueWithDOI(0.2),
        description="Surface albedo",
    )
    porosity: Optional[ValueWithDOI[float]] = Field(
        default=None,
        description="Porosity (deciduous trees only)",
    )
    soilstore: ValueWithDOI[float] = Field(
        default=ValueWithDOI(150.0),
        description="Soil moisture store [mm]",
    )
    snowpack: Optional[ValueWithDOI[float]] = Field(
        default=None,
        description="Snow water equivalent [mm]",
    )
    snowfrac: Optional[ValueWithDOI[float]] = Field(
        default=None,
        description="Snow covered fraction",
    )
    snowdens: Optional[ValueWithDOI[float]] = Field(
        default=None,
        description="Snow density [kg m⁻³]",
    )
    snowdepth: Optional[ValueWithDOI[float]] = Field(
        default=None,
        description="Snow depth [mm]",
    )
    snowpack_old: Optional[ValueWithDOI[float]] = Field(
        default=None,
        description="Previous timestep snow water equivalent [mm]",
    )
    snowdens_old: Optional[ValueWithDOI[float]] = Field(
        default=None,
        description="Previous timestep snow density [kg m⁻³]",
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert surface state to DataFrame format."""
        df_state = init_df_state(grid_id)

        # Map surface types to indices
        surf_idx = {
            SurfaceType.PAVED: 0,
            SurfaceType.BLDGS: 1,
            SurfaceType.EVETR: 2,
            SurfaceType.DECTR: 3,
            SurfaceType.GRASS: 4,
            SurfaceType.BSOIL: 5,
            SurfaceType.WATER: 6,
        }[self.surface_type]

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx},)"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = None
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Set basic state variables
        for attr in [
            "state",
            "gdd",
            "sdd",
            "lai",
            "albedo",
            "soilstore",
        ]:
            value = getattr(self, attr)
            if value is not None:
                set_df_value(attr, value.value)

        # Set porosity for deciduous trees
        if self.surface_type == SurfaceType.DECTR and self.porosity is not None:
            set_df_value("porosity", self.porosity.value)

        # Set snow-related variables if present
        if any(getattr(self, attr) is not None for attr in [
            "snowpack",
            "snowfrac",
            "snowdens",
            "snowdepth",
            "snowpack_old",
            "snowdens_old",
        ]):
            for attr in [
                "snowpack",
                "snowfrac",
                "snowdens",
                "snowdepth",
                "snowpack_old",
                "snowdens_old",
            ]:
                value = getattr(self, attr)
                if value is not None:
                    set_df_value(attr, value.value)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surface_type: SurfaceType) -> "SurfaceState":
        """Reconstruct surface state from DataFrame format."""
        # Map surface types to indices
        surf_idx = {
            SurfaceType.PAVED: 0,
            SurfaceType.BLDGS: 1,
            SurfaceType.EVETR: 2,
            SurfaceType.DECTR: 3,
            SurfaceType.GRASS: 4,
            SurfaceType.BSOIL: 5,
            SurfaceType.WATER: 6,
        }[surface_type]

        # Helper function to get values from DataFrame
        def get_df_value(col_name: str) -> float:
            return df.loc[grid_id, (col_name, f"({surf_idx},)")]

        # Initialize with surface type
        params = {"surface_type": surface_type}

        # Get basic state variables
        for attr in [
            "state",
            "gdd",
            "sdd",
            "lai",
            "albedo",
            "soilstore",
        ]:
            try:
                params[attr] = ValueWithDOI(get_df_value(attr))
            except KeyError:
                continue

        # Get porosity for deciduous trees
        if surface_type == SurfaceType.DECTR:
            try:
                params["porosity"] = ValueWithDOI(get_df_value("porosity"))
            except KeyError:
                pass

        # Get snow-related variables if present
        for attr in [
            "snowpack",
            "snowfrac",
            "snowdens",
            "snowdepth",
            "snowpack_old",
            "snowdens_old",
        ]:
            try:
                params[attr] = ValueWithDOI(get_df_value(attr))
            except KeyError:
                continue

        return cls(**params)


class MetState(BaseModel):
    """Meteorological state variables."""
    temp_c: ValueWithDOI[float] = Field(
        default=ValueWithDOI(20.0),
        description="Air temperature [°C]",
    )
    pressure: ValueWithDOI[float] = Field(
        default=ValueWithDOI(101300.0),
        description="Air pressure [Pa]",
    )
    rh: ValueWithDOI[float] = Field(
        ge=0, le=100,
        default=ValueWithDOI(50.0),
        description="Relative humidity [%]",
    )
    vpd: ValueWithDOI[float] = Field(
        default=ValueWithDOI(1000.0),
        description="Vapor pressure deficit [Pa]",
    )
    wind_speed: ValueWithDOI[float] = Field(
        default=ValueWithDOI(3.0),
        description="Wind speed [m s⁻¹]",
    )
    wind_dir: ValueWithDOI[float] = Field(
        ge=0, le=360,
        default=ValueWithDOI(180.0),
        description="Wind direction [degrees]",
    )
    kdown: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="Incoming shortwave radiation [W m⁻²]",
    )
    kup: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="Outgoing shortwave radiation [W m⁻²]",
    )
    ldown: ValueWithDOI[float] = Field(
        default=ValueWithDOI(340.0),
        description="Incoming longwave radiation [W m⁻²]",
    )
    lup: ValueWithDOI[float] = Field(
        default=ValueWithDOI(400.0),
        description="Outgoing longwave radiation [W m⁻²]",
    )
    rain: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.0),
        description="Rainfall rate [mm h⁻¹]",
    )
    snow: Optional[ValueWithDOI[float]] = Field(
        default=None,
        description="Snowfall rate [mm h⁻¹]",
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert meteorological state to DataFrame format."""
        df_state = init_df_state(grid_id)

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = "0"
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = None
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Set all meteorological variables
        for attr in [
            "temp_c",
            "pressure",
            "rh",
            "vpd",
            "wind_speed",
            "wind_dir",
            "kdown",
            "kup",
            "ldown",
            "lup",
            "rain",
        ]:
            value = getattr(self, attr)
            if value is not None:
                set_df_value(attr, value.value)

        # Set snow if present
        if self.snow is not None:
            set_df_value("snow", self.snow.value)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "MetState":
        """Reconstruct meteorological state from DataFrame format."""
        # Helper function to get values from DataFrame
        def get_df_value(col_name: str) -> float:
            return df.loc[grid_id, (col_name, "0")]

        # Initialize parameters dictionary
        params = {}

        # Get all meteorological variables
        for attr in [
            "temp_c",
            "pressure",
            "rh",
            "vpd",
            "wind_speed",
            "wind_dir",
            "kdown",
            "kup",
            "ldown",
            "lup",
            "rain",
            "snow",
        ]:
            try:
                params[attr] = ValueWithDOI(get_df_value(attr))
            except KeyError:
                continue

        return cls(**params)


class ModelState(BaseModel):
    """Complete model state including surface and meteorological states."""
    surfaces: Dict[SurfaceType, SurfaceState] = Field(
        description="State variables for each surface type",
    )
    meteorology: MetState = Field(
        default_factory=MetState,
        description="Meteorological state variables",
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert complete model state to DataFrame format."""
        # Initialize with meteorological state
        df_state = self.meteorology.to_df_state(grid_id)

        # Add each surface state
        for surface_state in self.surfaces.values():
            df_surface = surface_state.to_df_state(grid_id)
            df_state = pd.concat([df_state, df_surface], axis=1)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "ModelState":
        """Reconstruct complete model state from DataFrame format."""
        # Get meteorological state
        met_state = MetState.from_df_state(df, grid_id)

        # Get surface states
        surfaces = {}
        for surface_type in SurfaceType:
            surface_state = SurfaceState.from_df_state(df, grid_id, surface_type)
            surfaces[surface_type] = surface_state

        return cls(surfaces=surfaces, meteorology=met_state)
