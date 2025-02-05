"""Vegetation surface models for SUEWS."""

from typing import Dict, List, Optional, Union, Literal
import pandas as pd
from pydantic import BaseModel, Field, model_validator, PrivateAttr

from ..base import ValueWithDOI, Reference
from .base import SurfaceType, init_df_state
from .urban import WaterDistribution, StorageDrainParams, NonVegetatedSurfaceProperties


class LAIPowerCoefficients(BaseModel):
    """Power coefficients for LAI calculations."""
    growth_lai: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.1),
        description="Power coefficient for LAI in growth equation",
    )
    growth_gdd: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.1),
        description="Power coefficient for GDD in growth equation",
    )
    senescence_lai: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.1),
        description="Power coefficient for LAI in senescence equation",
    )
    senescence_sdd: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.1),
        description="Power coefficient for SDD in senescence equation",
    )

    ref: Optional[Reference] = None

    def to_df_state(self, grid_id: int, veg_idx: int) -> pd.DataFrame:
        """Convert LAI power coefficients to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            veg_idx: Vegetation index (0: EVETR, 1: DECTR, 2: GRASS)

        Returns:
            pd.DataFrame: DataFrame containing LAI power coefficients
        """
        df_state = init_df_state(grid_id)

        # Set power coefficients in order
        for i, value in enumerate(
            [
                self.growth_lai,
                self.growth_gdd,
                self.senescence_lai,
                self.senescence_sdd,
            ]
        ):
            df_state.loc[grid_id, ("laipower", f"({i}, {veg_idx})")] = value.value

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, veg_idx: int) -> "LAIPowerCoefficients":
        """Reconstruct LAIPowerCoefficients from DataFrame state format.

        Args:
            df: DataFrame containing LAI power coefficients
            grid_id: Grid ID for the DataFrame index
            veg_idx: Vegetation index (0: EVETR, 1: DECTR, 2: GRASS)

        Returns:
            LAIPowerCoefficients: Instance of LAIPowerCoefficients
        """
        # Map each coefficient to its corresponding index
        coefficients = [
            ValueWithDOI(df.loc[grid_id, ("laipower", f"(0, {veg_idx})")]),
            ValueWithDOI(df.loc[grid_id, ("laipower", f"(1, {veg_idx})")]),
            ValueWithDOI(df.loc[grid_id, ("laipower", f"(2, {veg_idx})")]),
            ValueWithDOI(df.loc[grid_id, ("laipower", f"(3, {veg_idx})")]),
        ]

        # Return the instance with coefficients
        return cls(
            growth_lai=coefficients[0],
            growth_gdd=coefficients[1],
            senescence_lai=coefficients[2],
            senescence_sdd=coefficients[3],
        )


class LAIParams(BaseModel):
    """Parameters for LAI calculations."""
    baset: ValueWithDOI[float] = Field(
        default=ValueWithDOI(10.0),
        description="Base Temperature for initiating growing degree days (GDD) for leaf growth [degC]",
    )
    gddfull: ValueWithDOI[float] = Field(
        default=ValueWithDOI(100.0),
        description="Growing degree days (GDD) needed for full capacity of LAI [degC]",
    )
    basete: ValueWithDOI[float] = Field(
        default=ValueWithDOI(10.0),
        description="Base temperature for initiating senescence degree days (SDD) for leaf off [degC]",
    )
    sddfull: ValueWithDOI[float] = Field(
        default=ValueWithDOI(100.0),
        description="Senescence degree days (SDD) needed to initiate leaf off [degC]",
    )
    laimin: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.1),
        description="Leaf-off wintertime value [m2 m-2]",
    )
    laimax: ValueWithDOI[float] = Field(
        default=ValueWithDOI(10.0),
        description="Full leaf-on summertime value [m2 m-2]",
    )
    laipower: LAIPowerCoefficients = Field(
        default_factory=LAIPowerCoefficients,
        description="LAI calculation power parameters for growth and senescence",
    )
    laitype: ValueWithDOI[int] = Field(
        default=ValueWithDOI(0),
        description="LAI calculation choice (0: original, 1: new high latitude)",
    )

    ref: Optional[Reference] = None

    @model_validator(mode="after")
    def validate_lai_ranges(self) -> "LAIParams":
        if self.laimin.value > self.laimax.value:
            raise ValueError(
                f"laimin ({self.laimin.value}) must be less than or equal to laimax ({self.laimax.value})."
            )
        if self.baset.value > self.gddfull.value:
            raise ValueError(
                f"baset {self.baset.value} must be less than gddfull ({self.gddfull.value})."
            )
        return self

    def to_df_state(self, grid_id: int, surf_idx: int) -> pd.DataFrame:
        """Convert LAI parameters to DataFrame state format.

        Args:
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index for vegetation (2: EVETR, 3: DECTR, 4: GRASS)

        Returns:
            pd.DataFrame: DataFrame containing LAI parameters
        """
        df_state = init_df_state(grid_id)

        # Adjust index for vegetation surfaces (surface index - 2)
        veg_idx = surf_idx - 2

        # Set basic LAI parameters
        for attr in ["baset", "gddfull", "basete", "sddfull", "laimin", "laimax", "laitype"]:
            df_state.loc[grid_id, (attr, f"({veg_idx},)")] = getattr(self, attr).value

        # Add LAI power coefficients
        df_power = self.laipower.to_df_state(grid_id, veg_idx)
        df_state = pd.concat([df_state, df_power], axis=1)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "LAIParams":
        """Reconstruct LAIParams from DataFrame state format.

        Args:
            df: DataFrame containing LAI parameters
            grid_id: Grid ID for the DataFrame index
            surf_idx: Surface index for vegetation (2: EVETR, 3: DECTR, 4: GRASS)

        Returns:
            LAIParams: Instance of LAIParams
        """
        # Adjust index for vegetation surfaces (surface index - 2)
        veg_idx = surf_idx - 2

        # Extract basic LAI parameters
        params = {}
        for attr in ["baset", "gddfull", "basete", "sddfull", "laimin", "laimax", "laitype"]:
            params[attr] = ValueWithDOI(df.loc[grid_id, (attr, f"({veg_idx},)")])

        # Extract LAI power coefficients
        laipower = LAIPowerCoefficients.from_df_state(df, grid_id, veg_idx)

        return cls(**params, laipower=laipower)


class VegetatedSurfaceProperties(NonVegetatedSurfaceProperties):
    """Base properties for vegetated surfaces."""
    alb_min: ValueWithDOI[float] = Field(
        ge=0, le=1,
        description="Minimum albedo",
        default=ValueWithDOI(0.2),
    )
    alb_max: ValueWithDOI[float] = Field(
        ge=0, le=1,
        description="Maximum albedo",
        default=ValueWithDOI(0.3),
    )
    beta_bioco2: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.6),
        description="Biogenic CO2 exchange coefficient",
    )
    beta_enh_bioco2: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.7),
        description="Enhanced biogenic CO2 exchange coefficient",
    )
    alpha_bioco2: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.8),
        description="Biogenic CO2 exchange coefficient",
    )
    alpha_enh_bioco2: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.9),
        description="Enhanced biogenic CO2 exchange coefficient",
    )
    resp_a: ValueWithDOI[float] = Field(
        default=ValueWithDOI(1.0),
        description="Respiration coefficient",
    )
    resp_b: ValueWithDOI[float] = Field(
        default=ValueWithDOI(1.1),
        description="Respiration coefficient",
    )
    theta_bioco2: ValueWithDOI[float] = Field(
        default=ValueWithDOI(1.2),
        description="Biogenic CO2 exchange coefficient",
    )
    maxconductance: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.5),
        description="Maximum surface conductance",
    )
    min_res_bioco2: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.1),
        description="Minimum respiratory biogenic CO2",
    )
    lai: LAIParams = Field(
        default_factory=LAIParams,
        description="Leaf area index parameters",
    )
    ie_a: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.5),
        description="Irrigation efficiency coefficient-automatic",
    )
    ie_m: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.6),
        description="Irrigation efficiency coefficient-manual",
    )

    ref: Optional[Reference] = None

    @model_validator(mode="after")
    def validate_albedo_range(self) -> "VegetatedSurfaceProperties":
        if self.alb_min.value > self.alb_max.value:
            raise ValueError(
                f"alb_min ({self.alb_min.value}) must be less than or equal to alb_max ({self.alb_max.value})."
            )
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert vegetated surface properties to DataFrame state format."""
        # Get base properties
        df_state = super().to_df_state(grid_id)

        # Get surface index
        surf_idx = self.get_surface_index()

        # Helper function to set values in DataFrame
        def set_df_value(col_name: str, value: float):
            idx_str = f"({surf_idx-2},)"  # Adjust for vegetation index
            if (col_name, idx_str) not in df_state.columns:
                df_state[(col_name, idx_str)] = None
            df_state.loc[grid_id, (col_name, idx_str)] = value

        # Add vegetation-specific properties
        for attr in [
            "beta_bioco2",
            "beta_enh_bioco2",
            "alpha_bioco2",
            "alpha_enh_bioco2",
            "resp_a",
            "resp_b",
            "theta_bioco2",
            "maxconductance",
            "min_res_bioco2",
            "ie_a",
            "ie_m",
        ]:
            set_df_value(attr, getattr(self, attr).value)

        # Add LAI parameters
        df_lai = self.lai.to_df_state(grid_id, surf_idx)
        df_state = pd.concat([df_state, df_lai], axis=1)

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int, surf_idx: int) -> "VegetatedSurfaceProperties":
        """Reconstruct vegetated surface properties from DataFrame state format."""
        instance = super().from_df_state(df, grid_id, surf_idx)

        # Helper function to get values from DataFrame
        def get_df_value(col_name: str) -> float:
            return df.loc[grid_id, (col_name, f"({surf_idx-2},)")]

        # Set vegetation-specific properties
        for attr in [
            "beta_bioco2",
            "beta_enh_bioco2",
            "alpha_bioco2",
            "alpha_enh_bioco2",
            "resp_a",
            "resp_b",
            "theta_bioco2",
            "maxconductance",
            "min_res_bioco2",
            "ie_a",
            "ie_m",
        ]:
            try:
                value = get_df_value(attr)
                setattr(instance, attr, ValueWithDOI(value))
            except KeyError:
                continue

        # Set LAI parameters
        instance.lai = LAIParams.from_df_state(df, grid_id, surf_idx)

        return instance


class EvetrProperties(VegetatedSurfaceProperties):
    """Properties for evergreen tree surfaces."""
    _surface_type: Literal[SurfaceType.EVETR] = SurfaceType.EVETR
    faievetree: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.1),
        description="Frontal area index of evergreen trees",
    )
    evetreeh: ValueWithDOI[float] = Field(
        default=ValueWithDOI(15.0),
        description="Evergreen tree height",
    )
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.EVETR),
        description="Water distribution for evergreen trees",
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert evergreen tree properties to DataFrame state format."""
        df_state = super().to_df_state(grid_id)

        # Add evergreen-specific properties
        df_state.loc[grid_id, ("faievetree", "0")] = self.faievetree.value
        df_state.loc[grid_id, ("evetreeh", "0")] = self.evetreeh.value
        df_state.loc[grid_id, ("albevetr_id", "0")] = self.alb_min.value

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "EvetrProperties":
        """Reconstruct evergreen tree properties from DataFrame state format."""
        surf_idx = 2
        instance = super().from_df_state(df, grid_id, surf_idx)

        # Set evergreen-specific properties
        instance.faievetree = ValueWithDOI(df.loc[grid_id, ("faievetree", "0")])
        instance.evetreeh = ValueWithDOI(df.loc[grid_id, ("evetreeh", "0")])
        instance.alb_min = ValueWithDOI(df.loc[grid_id, ("albevetr_id", "0")])

        return instance


class DectrProperties(VegetatedSurfaceProperties):
    """Properties for deciduous tree surfaces."""
    _surface_type: Literal[SurfaceType.DECTR] = SurfaceType.DECTR
    faidectree: ValueWithDOI[float] = Field(
        default=ValueWithDOI(0.1),
        description="Frontal area index of deciduous trees",
    )
    dectreeh: ValueWithDOI[float] = Field(
        default=ValueWithDOI(15.0),
        description="Deciduous tree height",
    )
    pormin_dec: ValueWithDOI[float] = Field(
        ge=0.1, le=0.9,
        default=ValueWithDOI(0.2),
        description="Minimum porosity",
    )
    pormax_dec: ValueWithDOI[float] = Field(
        ge=0.1, le=0.9,
        default=ValueWithDOI(0.6),
        description="Maximum porosity",
    )
    capmax_dec: ValueWithDOI[float] = Field(
        default=ValueWithDOI(100.0),
        description="Maximum capacity",
    )
    capmin_dec: ValueWithDOI[float] = Field(
        default=ValueWithDOI(10.0),
        description="Minimum capacity",
    )
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.DECTR),
        description="Water distribution for deciduous trees",
    )

    @model_validator(mode="after")
    def validate_porosity_range(self) -> "DectrProperties":
        if self.pormin_dec.value >= self.pormax_dec.value:
            raise ValueError(
                f"pormin_dec ({self.pormin_dec.value}) must be less than pormax_dec ({self.pormax_dec.value})."
            )
        return self

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert deciduous tree properties to DataFrame state format."""
        df_state = super().to_df_state(grid_id)

        # Add deciduous-specific properties
        df_state.loc[grid_id, ("faidectree", "0")] = self.faidectree.value
        df_state.loc[grid_id, ("dectreeh", "0")] = self.dectreeh.value
        df_state.loc[grid_id, ("pormin_dec", "0")] = self.pormin_dec.value
        df_state.loc[grid_id, ("pormax_dec", "0")] = self.pormax_dec.value
        df_state.loc[grid_id, ("capmax_dec", "0")] = self.capmax_dec.value
        df_state.loc[grid_id, ("capmin_dec", "0")] = self.capmin_dec.value
        df_state.loc[grid_id, ("albdectr_id", "0")] = self.alb_min.value

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "DectrProperties":
        """Reconstruct deciduous tree properties from DataFrame state format."""
        surf_idx = 3
        instance = super().from_df_state(df, grid_id, surf_idx)

        # Set deciduous-specific properties
        instance.faidectree = ValueWithDOI(df.loc[grid_id, ("faidectree", "0")])
        instance.dectreeh = ValueWithDOI(df.loc[grid_id, ("dectreeh", "0")])
        instance.pormin_dec = ValueWithDOI(df.loc[grid_id, ("pormin_dec", "0")])
        instance.pormax_dec = ValueWithDOI(df.loc[grid_id, ("pormax_dec", "0")])
        instance.capmax_dec = ValueWithDOI(df.loc[grid_id, ("capmax_dec", "0")])
        instance.capmin_dec = ValueWithDOI(df.loc[grid_id, ("capmin_dec", "0")])
        instance.alb_min = ValueWithDOI(df.loc[grid_id, ("albdectr_id", "0")])

        return instance


class GrassProperties(VegetatedSurfaceProperties):
    """Properties for grass surfaces."""
    _surface_type: Literal[SurfaceType.GRASS] = SurfaceType.GRASS
    waterdist: WaterDistribution = Field(
        default_factory=lambda: WaterDistribution(SurfaceType.GRASS),
        description="Water distribution for grass",
    )

    def to_df_state(self, grid_id: int) -> pd.DataFrame:
        """Convert grass properties to DataFrame state format."""
        df_state = super().to_df_state(grid_id)

        # Add grass-specific properties
        df_state.loc[grid_id, ("albgrass_id", "0")] = self.alb_min.value

        return df_state

    @classmethod
    def from_df_state(cls, df: pd.DataFrame, grid_id: int) -> "GrassProperties":
        """Reconstruct grass properties from DataFrame state format."""
        surf_idx = 4
        instance = super().from_df_state(df, grid_id, surf_idx)

        # Set grass-specific properties
        instance.alb_min = ValueWithDOI(df.loc[grid_id, ("albgrass_id", "0")])

        return instance
